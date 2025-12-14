# -*- coding: utf-8 -*-
"""
HADDOCK3 Integration Module
最小可用封装：
- 检测 haddock3 是否可用（优先 CLI，其次 Python 包）
- 生成一个最小配置（topoaa -> rigidbody），执行 `haddock3 <config>`
- 收集 rigidbody 阶段生成的若干 PDB，合并为多模型 PDB 以便一次性载入 PyMOL

注意：当前版本未将 rsite/lsite 转换为 HADDOCK 约束（.tbl）。
后续可通过 haddock3-cfg / haddock-restraints 扩展。
"""
from __future__ import annotations

import os
import sys
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def _validate_haddock_path(path: str) -> bool:
    """验证路径是否符合 HADDOCK3 的字符要求: [a-zA-Z0-9._-/\\]"""
    import re
    return bool(re.match(r'^[a-zA-Z0-9._/\\-]+$', path))


def _which(exe: str) -> Optional[str]:
    # 更稳健的 CLI 查找：PATH、CONDA_PREFIX、常见用户路径、Homebrew
    path = shutil.which(exe)
    if path:
        return path
    candidates = []
    # 1) Conda 前缀
    conda_prefix = os.environ.get('CONDA_PREFIX')
    if conda_prefix:
        candidates.append(os.path.join(conda_prefix, 'bin', exe))
    # 2) 常见路径
    home = os.path.expanduser('~')
    candidates += [
        os.path.join(home, 'bin', exe),
        os.path.join(home, '.local', 'bin', exe),
        os.path.join(home, 'local', 'bin', exe),
        f'/usr/local/bin/{exe}',
    ]
    # 3) Homebrew on macOS
    if sys.platform == 'darwin':
        candidates.append(f'/opt/homebrew/bin/{exe}')
    for p in candidates:
        if p and os.path.exists(p) and os.access(p, os.X_OK):
            return p
    return None


def check_haddock3_available() -> Dict:
    """
    检测 HADDOCK3 是否可用。

    返回:
        dict: {
            'available': bool,
            'via': 'cli' | 'module' | None,
            'version': Optional[str],
            'detail': Optional[str]
        }
    """
    # 1) 优先检查 CLI
    cli = _which('haddock3')
    if cli:
        ver = None
        try:
            out = subprocess.run([cli, '-v'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=5)
            ver = (out.stdout or '').strip().splitlines()[-1]
        except Exception:
            pass
        return {'available': True, 'via': 'cli', 'version': ver, 'detail': cli}

    # 2) 退回 Python 包探测（import 名称保持兼容性）
    try:
        import importlib.util as _ilus  # noqa
        spec = _ilus.find_spec('haddock') or _ilus.find_spec('haddock3')
        if spec is not None:
            # 读取版本（尽力而为）
            ver = None
            try:
                import importlib.metadata as _ilm  # py3.8+
                for pkg in ('haddock3', 'haddock'):
                    try:
                        ver = _ilm.version(pkg)
                        if ver:
                            break
                    except Exception:
                        continue
            except Exception:
                pass
            return {'available': True, 'via': 'module', 'version': ver, 'detail': str(spec.origin)}
    except Exception:
        pass

    return {'available': False, 'via': None, 'version': None, 'detail': None}


class Haddock3Runner:
    def __init__(self, use_cli: Optional[bool] = None, ncores: Optional[int] = None):
        """
        初始化 Runner。
        - use_cli: 强制使用 CLI；若为 None 则按可用性自动选择（优先 CLI）。
        - ncores: 全局并行核数（缺省不写入，走默认）。
        """
        info = check_haddock3_available()
        if not info['available']:
            raise RuntimeError('HADDOCK3 未检测到，请先 `pip install haddock3` 或将其加入 PATH。')
        if use_cli is None:
            self.use_cli = (info['via'] == 'cli')
        else:
            self.use_cli = bool(use_cli)
        self.ncores = ncores
        self._cli_path = info['detail'] if info['via'] == 'cli' else _which('haddock3')

    @staticmethod
    def _parse_resid_ranges(txt: str) -> List[Tuple[str, int, int]]:
        """
        解析形如 "195:A,203-206:A,108:B" 的字符串，返回 [(chain, start, end), ...]
        """
        if not txt:
            return []
        import re
        txt = txt.replace(';', ',').replace('\n', ',').replace('\t', ',')
        toks = [t.strip() for t in txt.split(',') if t.strip()]
        out: List[Tuple[str, int, int]] = []
        for t in toks:
            m = re.match(r"^(\d+)(?:-(\d+))?:(\w)$", t)
            if not m:
                # 放宽：若匹配不到，尝试无链信息（默认链 A）
                m2 = re.match(r"^(\d+)(?:-(\d+))?$", t)
                if m2:
                    s = int(m2.group(1)); e = int(m2.group(2) or m2.group(1)); out.append(('A', s, e))
                continue
            s = int(m.group(1)); e = int(m.group(2) or m.group(1)); ch = m.group(3)
            if s > e:
                s, e = e, s
            out.append((ch, s, e))
        return out

    @staticmethod
    def _which_restraints_cli() -> Optional[str]:
        for exe in ('haddock3-restraints', 'haddock-restraints'):
            p = _which(exe)
            if p:
                return p
        return None

    def _write_act_pass_files(self, out_dir: Path, receptor_pdb: str, ligand_pdb: str,
                              rsite: Optional[str], lsite: Optional[str], expand_passive: bool) -> Tuple[Path, Path]:
        """
        生成两个 .act-pass 文件，格式为：
        第一行：active residues（空格分隔）
        第二行（可选）：passive residues（空格分隔）
        若 expand_passive 为 True 且可用 CLI，则自动计算被动残基；否则仅写 active 行。
        """
        r_ranges = self._parse_resid_ranges(rsite or '')
        l_ranges = self._parse_resid_ranges(lsite or '')
        def _expand_ranges(ranges: List[Tuple[str, int, int]]) -> Tuple[str, List[int]]:
            if not ranges:
                return ('A', [])
            chain = ranges[0][0]
            res = []
            for ch, s, e in ranges:
                chain = ch  # 简化：假设同一链
                res.extend(list(range(s, e+1)))
            return (chain, sorted(set(res)))
        r_chain, r_act = _expand_ranges(r_ranges)
        l_chain, l_act = _expand_ranges(l_ranges)
        if not r_act or not l_act:
            raise ValueError('限制口袋模式需要同时提供受体与配体的残基列表（rsite/lsite）。')
        # 写 active 行
        ap1 = out_dir / 'receptor.act-pass'
        ap2 = out_dir / 'ligand.act-pass'
        ap1.write_text(' '.join(str(x) for x in r_act) + '\n')
        ap2.write_text(' '.join(str(x) for x in l_act) + '\n')
        # 需要被动残基时，尽量调用 CLI 计算
        if expand_passive:
            cli = self._which_restraints_cli()
            if cli:
                try:
                    # receptor
                    cmd1 = [cli, 'passive_from_active', receptor_pdb, '--active-residues', *[str(x) for x in r_act]]
                    p1 = subprocess.run(cmd1, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
                    r_pass = p1.stdout.strip().split()
                    if r_pass:
                        with open(ap1, 'a') as f:
                            f.write(' '.join(r_pass) + '\n')
                    # ligand
                    cmd2 = [cli, 'passive_from_active', ligand_pdb, '--active-residues', *[str(x) for x in l_act]]
                    p2 = subprocess.run(cmd2, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
                    l_pass = p2.stdout.strip().split()
                    if l_pass:
                        with open(ap2, 'a') as f:
                            f.write(' '.join(l_pass) + '\n')
                except Exception:
                    # 忽略被动扩展失败，保留 active 行
                    pass
        return ap1, ap2

    def _prepare_restraints(self, out_dir: Path, receptor_pdb: str, ligand_pdb: str,
                            rsite: Optional[str], lsite: Optional[str], mode: str,
                            expand_passive: bool) -> Optional[Path]:
        """按模式准备约束文件，返回 ambig.tbl 路径或 None。"""
        if mode == 'air_from_residues':
            ap1, ap2 = self._write_act_pass_files(out_dir, receptor_pdb, ligand_pdb, rsite, lsite, expand_passive)
            # 优先使用官方 CLI 生成 AMBIG
            cli = self._which_restraints_cli()
            if cli:
                # 获取默认链标识（从 rsite/lsite 解析）
                r_ranges = self._parse_resid_ranges(rsite or '')
                l_ranges = self._parse_resid_ranges(lsite or '')
                r_chain = r_ranges[0][0] if r_ranges else 'A'
                l_chain = l_ranges[0][0] if l_ranges else 'B'
                ambig = out_dir / 'ambig.tbl'
                cmd = [cli, 'active_passive_to_ambig', str(ap1), str(ap2), '--segid-one', r_chain, '--segid-two', l_chain]
                try:
                    with open(ambig, 'w') as fout:
                        subprocess.run(cmd, stdout=fout, stderr=subprocess.STDOUT, text=True, check=True)
                    return ambig
                except Exception:
                    # 回退到简易 AMBIG 生成
                    return self._write_simple_ambig(out_dir, rsite, lsite)
            else:
                return self._write_simple_ambig(out_dir, rsite, lsite)
        else:
            # 盲对接模式：不写 ambig，依靠 cmrest/ranair/surfrest
            return None

    @staticmethod
    def _write_simple_ambig(out_dir: Path, rsite: Optional[str], lsite: Optional[str]) -> Path:
        """
        简易 AMBIG 生成：使用 active 残基全集做两两配对，写 CNS 约束。
        注意：这是兜底方案，质量不及官方生成器，推荐用户安装 haddock3-restraints。
        """
        r = Haddock3Runner._parse_resid_ranges(rsite or '')
        l = Haddock3Runner._parse_resid_ranges(lsite or '')
        if not r or not l:
            raise ValueError('无法生成 AMBIG：残基列表为空。')
        r_chain = r[0][0] if r else 'A'
        l_chain = l[0][0] if l else 'B'
        r_list: List[int] = []
        l_list: List[int] = []
        for ch, s, e in r:
            r_list.extend(range(s, e+1))
        for ch, s, e in l:
            l_list.extend(range(s, e+1))
        r_list = sorted(set(r_list))
        l_list = sorted(set(l_list))
        ambig = out_dir / 'ambig.tbl'
        with open(ambig, 'w') as f:
            for i in r_list:
                for j in l_list:
                    # 采用 CA-CA 简化约束，目标上限 8Å（示例），CNS 三元为：target lower-upper
                    f.write(f"assign (segid {r_chain} and resi {i} and name CA) (segid {l_chain} and resi {j} and name CA) 8.0 0.0 0.0\n")
        return ambig

    @staticmethod
    def _make_config_text(run_dir: str, receptor_pdb: str, ligand_pdb: str, ncores: Optional[int], rigidbody_sampling: int,
                          mode: str,
                          ambig_fname: Optional[str]) -> str:
        """
        生成最小可用配置。遵循 HADDOCK3 文档格式：
        - 全局: run_dir, molecules, (可选) ncores
        - 步骤: [topoaa], [rigidbody]（指定 sampling）
        参考: docs “Workflow configuration files”。
        """
        lines: List[str] = []
        lines.append(f'run_dir = "{run_dir}"')
        lines.append('')
        lines.append('molecules = [')
        lines.append(f'    "{receptor_pdb}",')
        lines.append(f'    "{ligand_pdb}"')
        lines.append(']')
        if ncores and ncores > 0:
            lines.append('')
            lines.append(f'ncores = {int(ncores)}')
        # steps
        lines.append('')
        lines.append('[topoaa]')
        # 保守使用默认参数，避免在不同体系上过拟合
        lines.append('')
        lines.append('[rigidbody]')
        lines.append(f'sampling = {int(rigidbody_sampling)}')
        # 约束/模式选择：至少提供一种约束，避免 guardrail 报错
        if ambig_fname:
            lines.append(f'ambig_fname = "{ambig_fname}"')
        if mode == 'blind_cm':
            lines.append('cmrest = true')
        elif mode == 'blind_ranair':
            lines.append('ranair = true')
        elif mode == 'blind_surf':
            lines.append('surfrest = true')
        # 可在后续扩展: randremoval, npart, flexref, caprieval 等
        return "\n".join(lines) + "\n"

    @staticmethod
    def _gather_pdbs(run_dir: Path, max_models: int = 10) -> List[Path]:
        """在 run 目录下查找 rigidbody 产物 PDB，返回最多 max_models 个路径。"""
        candidates: List[Path] = []
        for root, _dirs, files in os.walk(run_dir):
            for fn in files:
                if fn.lower().endswith('.pdb'):
                    p = Path(root) / fn
                    # 过滤常见的中间文件（策略保守，仅排除明显非模型）
                    if 'rigidbody' in str(p.parent):
                        candidates.append(p)
        # 排序策略：按文件大小&文件名自然序综合排序，尽量把体积更大的模型靠前
        candidates.sort(key=lambda p: (-(p.stat().st_size), str(p)))
        return candidates[:max_models]

    @staticmethod
    def _merge_models(pdb_paths: List[Path], out_path: Path) -> None:
        """将多个单模型 PDB 合并为多模型 PDB（MODEL/ENDMDL 包裹）。"""
        with open(out_path, 'w') as fout:
            for i, p in enumerate(pdb_paths, 1):
                fout.write(f"MODEL     {i}\n")
                with open(p, 'r') as fin:
                    for line in fin:
                        # 粗略过滤 END, TER，避免重复终止记录
                        if line.startswith(('END', 'TER')):
                            continue
                        fout.write(line)
                fout.write("ENDMDL\n")
        # 追加标准 END 记录，兼容部分解析器
        with open(out_path, 'a') as fout:
            fout.write('END\n')

    def run_docking(
        self,
        receptor_pdb: str,
        ligand_pdb: str,
        output_dir: Optional[str] = None,
        rsite: Optional[str] = None,
        lsite: Optional[str] = None,
        n_models: int = 10,
        rigidbody_sampling: int = 60,
        mode: str = 'blind_ranair',
        expand_passive: bool = True,
    ) -> Dict:
        """
        运行 HADDOCK3 最小对接流程。

        参数:
            receptor_pdb: 受体 PDB 路径
            ligand_pdb: 配体 PDB 路径
            output_dir: 输出目录（默认当前目录下 haddock3_results）
            rsite/lsite: 站位残基（限制口袋模式下用于生成 AIR）
            n_models: 合并输出的模型数量（默认 10，均衡）
            rigidbody_sampling: 刚体采样数量（默认 60，均衡；越大越久）
            mode: 'blind_ranair'|'blind_cm'|'blind_surf'|'air_from_residues'
            expand_passive: 在基于残基 AIR 模式下，是否自动扩展被动残基（6.5Å）
        返回:
            dict: {'success': bool, 'output_dir': str, 'models_pdb': str, 'run_dir': str}
        """
        receptor_pdb = os.path.abspath(receptor_pdb)
        ligand_pdb = os.path.abspath(ligand_pdb)
        if not os.path.exists(receptor_pdb) or not os.path.exists(ligand_pdb):
            return {'success': False, 'error': '输入 PDB 文件不存在。'}

        # 智能选择输出目录:
        # 1. 优先使用用户指定的 output_dir
        # 2. 其次在受体文件所在目录创建输出(结果和输入放一起)
        # 3. 兜底:使用当前目录,如果不可写则用用户主目录
        # 4. HADDOCK3 路径限制:仅允许 [a-zA-Z0-9._-/\\] 字符
        if output_dir:
            out_dir = Path(output_dir)
        else:
            # 尝试在受体文件所在目录创建输出
            receptor_dir = Path(receptor_pdb).parent
            if os.access(receptor_dir, os.W_OK):
                out_dir = receptor_dir / 'haddock3_results'
            else:
                # 受体目录不可写,尝试当前目录
                cwd = Path(os.getcwd())
                if cwd != Path('/') and os.access(cwd, os.W_OK):
                    out_dir = cwd / 'haddock3_results'
                else:
                    # 兜底:使用用户主目录
                    out_dir = Path.home() / 'haddock3_results'
        
        # 验证路径是否符合 HADDOCK3 要求
        path_fallback_warning = None
        out_dir_str = str(out_dir.resolve())
        if not _validate_haddock_path(out_dir_str):
            # 路径包含非法字符(如中文),回退到临时目录
            original_path = out_dir_str
            import tempfile
            tmp_base = Path(tempfile.gettempdir()) / 'haddock3_results'
            # 使用时间戳确保唯一性
            import time
            timestamp = int(time.time())
            out_dir = tmp_base / f'run_{timestamp}'
            out_dir_str = str(out_dir.resolve())
            # 再次验证(临时目录应该总是安全的)
            if not _validate_haddock_path(out_dir_str):
                return {'success': False, 'error': f'无法找到符合 HADDOCK3 路径要求的输出目录。HADDOCK3 仅支持 ASCII 字符路径。'}
            path_fallback_warning = f'路径包含非 ASCII 字符，已自动回退到临时目录: {out_dir_str}'
        
        out_dir.mkdir(parents=True, exist_ok=True)
        run_dir = out_dir / 'run1'
        # 配置文件路径
        cfg_path = out_dir / 'config.toml'

        # 约束文件（可选）
        ambig_path: Optional[Path] = None
        try:
            ambig_path = self._prepare_restraints(out_dir, receptor_pdb, ligand_pdb, rsite, lsite, mode, expand_passive)
        except Exception as e:
            return {'success': False, 'error': f'生成约束失败: {e}'}

        # 若是盲对接，经验上需要增大采样
        if mode.startswith('blind_') and rigidbody_sampling < 200:
            rigidbody_sampling = 300

        # 生成配置
        cfg_text = self._make_config_text(
            run_dir=str(run_dir),
            receptor_pdb=receptor_pdb,
            ligand_pdb=ligand_pdb,
            ncores=self.ncores,
            rigidbody_sampling=rigidbody_sampling,
            mode=mode,
            ambig_fname=(str(ambig_path) if ambig_path else None),
        )
        cfg_path.write_text(cfg_text)

        # 执行
        try:
            if self.use_cli:
                exe = self._cli_path or 'haddock3'
                cmd = [exe, str(cfg_path)]
                proc = subprocess.run(
                    cmd, cwd=str(out_dir), text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True
                )
                log = proc.stdout
            else:
                # 兜底：通过 Python 包入口执行（与 CLI 等价，仍以子进程形式，避免长时阻塞当前进程上下文）
                exe = sys.executable
                code = (
                    'import sys; from haddock.clis import haddock3 as _h; '
                    'sys.argv=["haddock3", sys.argv[1]]; _h.main()'
                )
                cmd = [exe, '-c', code, str(cfg_path)]
                proc = subprocess.run(
                    cmd, cwd=str(out_dir), text=True,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True
                )
                log = proc.stdout
        except subprocess.CalledProcessError as e:
            return {'success': False, 'error': f'HADDOCK3 运行失败: {e.stdout}'}
        except Exception as e:
            return {'success': False, 'error': f'HADDOCK3 执行异常: {e}'}

        # 收集模型
        pdbs = self._gather_pdbs(run_dir, max_models=max(1, int(n_models)))
        if not pdbs:
            return {'success': False, 'error': '未在 rigidbody 结果中找到 PDB 模型。'}

        merged = out_dir / f'haddock3_top{len(pdbs)}.pdb'
        try:
            self._merge_models(pdbs, merged)
        except Exception as e:
            return {'success': False, 'error': f'合并模型失败: {e}'}

        result = {
            'success': True,
            'output_dir': str(out_dir),
            'run_dir': str(run_dir),
            'models_pdb': str(merged),
            'log_tail': (log.splitlines()[-50:] if isinstance(log, str) else None)
        }
        if path_fallback_warning:
            result['warning'] = path_fallback_warning
        return result
