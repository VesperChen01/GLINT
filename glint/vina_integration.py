# -*- coding: utf-8 -*-
"""
vina_integration.py
GLINT - AutoDock Vina 完整集成Module
"""

from __future__ import print_function
import os
import sys
import subprocess
import tempfile
import shutil
import glob
from pathlib import Path
from typing import List, Dict, Any, Optional, Callable
from pymol import cmd

# Version标记 - 用于Confirm代码是否被正确Load
_VINA_MODULE_VERSION = "2026-01-03-v2"
print(f"[vina_integration] ModuleVersion: {_VINA_MODULE_VERSION}")


def find_vina_executable():
    """自动Find Vina 可执行File"""
    vina = shutil.which('vina') or shutil.which('vina.exe')
    if vina:
        return vina

    if 'CONDA_PREFIX' in os.environ:
        conda_vina = os.path.join(os.environ['CONDA_PREFIX'], 'bin', 'vina')
        if os.path.exists(conda_vina):
            return conda_vina

    # 3. 遍历常见 conda 安装Path回退（与 find_obabel_executable 保持一致）
    home = os.path.expanduser('~')
    for env_name in ["glint", "base"]:
        for base in [
            f"{home}/miniconda3",
            f"{home}/anaconda3",
            f"{home}/opt/miniconda3",
            "/opt/miniconda3",
            "/opt/anaconda3",
            "/opt/homebrew/Caskroom/miniconda/base",
            "/usr/local/Caskroom/miniconda/base",
        ]:
            vina_path = f"{base}/envs/{env_name}/bin/vina"
            if os.path.exists(vina_path):
                return vina_path

    # 4. 其他常见Path
    user_paths = [
        os.path.join(home, 'bin', 'vina'),
        os.path.join(home, '.local', 'bin', 'vina'),
        '/usr/local/bin/vina',
        '/opt/homebrew/bin/vina',
    ]

    for path in user_paths:
        if '*' in path:
            matches = glob.glob(path)
            if matches:
                return matches[0]
        elif os.path.exists(path):
            return path
    return None


def find_obabel_executable():
    """Find Open Babel 可执行File"""
    if os.environ.get("OBABEL_BINARY"):
        return os.environ.get("OBABEL_BINARY")

    obabel = shutil.which("obabel") or shutil.which("obabel.exe")
    if obabel:
        return obabel

    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        conda_obabel = os.path.join(conda_prefix, "bin", "obabel")
        if os.path.exists(conda_obabel):
            return conda_obabel

    home = os.path.expanduser("~")
    for env_name in ["glint", "base"]:
        for base in [f"{home}/miniconda3", f"{home}/anaconda3", "/opt/homebrew/Caskroom/miniconda/base"]:
            obabel_path = f"{base}/envs/{env_name}/bin/obabel"
            if os.path.exists(obabel_path):
                return obabel_path
    return None


def find_mgltools_scripts():
    """Find MGLTools 脚本Path"""
    possible_paths = ['/usr/local/MGLTools', '/opt/mgltools', os.path.join(os.path.expanduser('~'), 'MGLTools')]
    if 'MGLTOOLS_HOME' in os.environ:
        possible_paths.insert(0, os.environ['MGLTOOLS_HOME'])

    for base_path in possible_paths:
        prepare_ligand = os.path.join(base_path, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_ligand4.py')
        prepare_receptor = os.path.join(base_path, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_receptor4.py')
        if os.path.exists(prepare_ligand) and os.path.exists(prepare_receptor):
            return {'prepare_ligand': prepare_ligand, 'prepare_receptor': prepare_receptor,
                    'python': os.path.join(base_path, 'bin', 'pythonsh')}
    return None


def check_vina_available():
    """检查 Vina 是否可用"""
    vina = find_vina_executable()
    if vina:
        print(f"[vina_integration] ✅ Found Vina: {vina}")
        return True
    print("[vina_integration] ⚠️ Vina not found. Install: conda install -c conda-forge vina")
    return False


def export_to_pdbqt(obj_name, selection, output_pdbqt, is_receptor=True):
    """标准 PDBQT Export

    对于受体：using -xr Parameters生成刚性受体格式（无 ROOT/ENDROOT Label）
    对于配体：using -h ParametersAdd氢原子
    """
    print(f"[vina_integration] export_to_pdbqt: is_receptor={is_receptor}, output={output_pdbqt}")

    temp_pdb = output_pdbqt.replace('.pdbqt', '_temp.pdb')

    # 防御性检查：避免无效/空 selection 导致 cmd.save Selector-Error 或导出空文件
    try:
        atom_count = cmd.count_atoms(selection)
    except Exception as e:
        raise ValueError(f"[vina_integration] Invalid PyMOL selection: '{selection}'. Error: {e}")
    if atom_count <= 0:
        raise ValueError(
            f"[vina_integration] Selection matches 0 atoms: '{selection}'. Check obj_name/remove_selection."
        )

    cmd.save(temp_pdb, selection)

    if not os.path.exists(temp_pdb):
        print(f"[vina_integration] ❌ 临时 PDB FileCreateFailed")
        return False

    obabel = find_obabel_executable()
    if obabel:
        try:
            if is_receptor:
                # 受体：using -xr 生成刚性受体格式，不Package含 ROOT Label
                cmd_args = [obabel, temp_pdb, '-O', output_pdbqt, '-xr']
                print(f"[vina_integration] using -xr Parameters（刚性受体格式）")
            else:
                # 配体：Add氢原子
                cmd_args = [obabel, temp_pdb, '-O', output_pdbqt, '-h']
                print(f"[vina_integration] using -h Parameters（配体格式）")

            result = subprocess.run(cmd_args, capture_output=True, text=True, timeout=60)
            if result.returncode == 0 and os.path.exists(output_pdbqt):
                os.unlink(temp_pdb)
                print(f"[vina_integration] ✅ PDBQT ExportSuccess")
                return True
            else:
                print(f"[vina_integration] ❌ obabel Failed: {result.stderr}")
        except Exception as e:
            print(f"[vina_integration] ❌ 异常: {e}")

    if os.path.exists(temp_pdb):
        os.unlink(temp_pdb)
    return False


def convert_ligand_to_pdbqt(ligand_file: str, output_pdbqt: str) -> bool:
    """将配体File转换为 PDBQT（不通过 PyMOL）"""
    # 如果源File和目标File相同，直接ReturnSuccess
    if os.path.abspath(ligand_file) == os.path.abspath(output_pdbqt):
        return True
    obabel = find_obabel_executable()
    if not obabel:
        return False
    try:
        result = subprocess.run([obabel, ligand_file, '-O', output_pdbqt, '-h'],
                               capture_output=True, text=True, timeout=60)
        return result.returncode == 0 and os.path.exists(output_pdbqt)
    except Exception:
        return False


def generate_vina_config(receptor_pdbqt: str, ligand_pdbqt: str, box_params: Dict[str, float],
                        output_pdbqt: str, config_path: Optional[str] = None,
                        exhaustiveness: int = 8, num_modes: int = 9) -> str:
    """生成 Vina ConfigurationFile"""
    if config_path is None:
        fd, config_path = tempfile.mkstemp(suffix='_vina_config.txt', text=True)
        os.close(fd)

    config_content = f"""receptor = {os.path.abspath(receptor_pdbqt)}
ligand = {os.path.abspath(ligand_pdbqt)}
out = {os.path.abspath(output_pdbqt)}
center_x = {box_params['center_x']:.3f}
center_y = {box_params['center_y']:.3f}
center_z = {box_params['center_z']:.3f}
size_x = {box_params['size_x']:.1f}
size_y = {box_params['size_y']:.1f}
size_z = {box_params['size_z']:.1f}
exhaustiveness = {exhaustiveness}
num_modes = {num_modes}
"""
    with open(config_path, 'w') as f:
        f.write(config_content)
    return config_path


def run_vina_docking(receptor_pdbqt: str, ligand_pdbqt: str, box_params: Dict[str, float],
                     output_dir: str, ligand_name: str, exhaustiveness: int = 8,
                     num_modes: int = 9) -> Dict[str, Any]:
    """运行单个 Vina 对接"""
    vina_bin = find_vina_executable()
    if not vina_bin:
        return {'success': False, 'error': 'Vina not found'}

    output_pdbqt = os.path.join(output_dir, f"{ligand_name}_out.pdbqt")
    log_file = os.path.join(output_dir, f"{ligand_name}.log")
    config_path = os.path.join(output_dir, f"{ligand_name}_config.txt")

    generate_vina_config(receptor_pdbqt, ligand_pdbqt, box_params, output_pdbqt,
                        config_path, exhaustiveness, num_modes)

    try:
        print(f"[vina_integration] 运行 Vina: {vina_bin}")
        # 新版 Vina 不支持 --log Parameters，输出到 stdout
        result = subprocess.run([vina_bin, '--config', config_path],
                      capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            error_msg = result.stderr.strip() if result.stderr else 'Vina execution failed'
            print(f"[vina_integration] ❌ Vina Failed: {error_msg}")
            return {'success': False, 'error': error_msg, 'ligand': ligand_name}

        # Save输出到日志File
        with open(log_file, 'w') as f:
            f.write(result.stdout)

    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Timeout', 'ligand': ligand_name}
    except Exception as e:
        return {'success': False, 'error': str(e), 'ligand': ligand_name}

    # 从 stdout 或日志File中解析Results
    affinity = None
    output_text = result.stdout if result.stdout else ""
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            output_text = f.read()

    for line in output_text.split('\n'):
        if line.strip().startswith('1 '):
            parts = line.split()
            try:
                affinity = float(parts[1])
                break
            except (ValueError, IndexError):
                pass

    if affinity is None:
        return {'success': False, 'error': 'Parse failed', 'ligand': ligand_name}

    print(f"[vina_integration] ✅ 对接Success: {affinity:.2f} kcal/mol")
    return {'success': True, 'affinity': affinity, 'output_pdbqt': output_pdbqt,
            'log_file': log_file, 'ligand': ligand_name}



def _normalize_remove_selection(obj_name: str, remove_selection: Optional[str]) -> Optional[str]:
    """将用户输入的 remove_selection 规范为 PyMOL 可识别的 selection 表达式。

    兼容场景：
    - 已存在的 selection 名称 / object 名称：直接使用
    - 链名简写：'A' 或 'EF2' -> 'chain A' / 'chain EF2'
    - 残基名简写：'LIG' -> 'resn LIG'
    - 复合选择：'chain A or resn LIG' -> 原样使用

    返回值为可直接放入 "not (<sel>)" 的 selection 字符串；若输入为空则返回 None。
    """
    if not remove_selection:
        return None

    sel = str(remove_selection).strip()
    if not sel:
        return None

    def _try_count(sel_expr: str) -> Optional[int]:
        """对 selection 进行安全计数：解析异常返回 None。"""
        try:
            cnt = cmd.count_atoms(sel_expr)
            return int(cnt) if cnt is not None else 0
        except Exception:
            return None

    # 复合表达式通常包含空白/括号/逻辑操作符：原样保留，但需要验证可解析且命中>0，否则忽略
    sel_lower = sel.lower()
    if (
        any(ch in sel for ch in (" ", "\t", "(", ")"))
        or " and " in sel_lower
        or " or " in sel_lower
        or sel_lower.startswith("not ")
    ):
        cnt = _try_count(f"({obj_name}) and ({sel})")
        if cnt is not None and cnt > 0:
            return sel
        return None

    # 向后兼容：如果用户传入的是已存在 selection / object 名称，则直接使用（但仍需验证命中>0）
    try:
        existing_sels = set(cmd.get_names("selections"))
    except Exception:
        existing_sels = set()

    try:
        existing_objs = set(cmd.get_names("objects"))
    except Exception:
        existing_objs = set()

    if sel in existing_sels or sel in existing_objs:
        cnt = _try_count(f"({obj_name}) and ({sel})")
        if cnt is not None and cnt > 0:
            return sel
        return None

    token = sel

    # 尝试用实际 atom count 判断用户想表达 chain / resn / segi / resi
    candidates: List[str] = []

    # 纯数字或范围更像 resi
    token_no_dash = token.replace("-", "")
    if token_no_dash.isdigit():
        candidates.append(f"resi {token}")

    candidates.extend([
        f"chain {token}",
        f"segi {token}",
        f"resn {token}",
    ])

    # 按候选优先顺序逐一验证：只有 count>0 才返回；非法 selection 则跳过
    for cand in candidates:
        cnt = _try_count(f"({obj_name}) and ({cand})")
        if cnt is not None and cnt > 0:
            return cand

    # 所有候选均无匹配或解析失败：忽略 remove_selection
    return None


def manual_box_docking(obj_name: str, ligand_file: str, box_params: Dict[str, float],
                       output_dir: Optional[str] = None, exhaustiveness: int = 8,
                       num_modes: int = 9, remove_selection: Optional[str] = None) -> Dict[str, Any]:
    """using自定义对接盒Parameters进行 Vina 对接"""
    if not check_vina_available():
        return {'success': False, 'error': 'Vina not available'}

    try:
        objs = cmd.get_names("objects")
    except AttributeError:
        objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []

    if obj_name not in objs:
        return {'success': False, 'error': f"Object '{obj_name}' not found"}

    if output_dir is None:
        output_dir = f"{obj_name}_docking"
    os.makedirs(output_dir, exist_ok=True)

    # Export受体 PDBQT
    receptor_pdbqt = os.path.join(output_dir, f"{obj_name}_receptor.pdbqt")
    normalized_remove_sel = _normalize_remove_selection(obj_name, remove_selection)
    receptor_sel = f"({obj_name})"
    if normalized_remove_sel and str(normalized_remove_sel).strip():
        receptor_sel += f" and not ({normalized_remove_sel})"

    try:
        if not export_to_pdbqt(obj_name, receptor_sel, receptor_pdbqt, is_receptor=True):
            return {'success': False, 'error': 'Receptor export failed'}
    except ValueError as e:
        return {'success': False, 'error': str(e)}

    # Export配体 PDBQT
    ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
    ligand_obj = "temp_ligand"
    cmd.load(ligand_file, ligand_obj)
    try:
        if not export_to_pdbqt(ligand_obj, ligand_obj, ligand_pdbqt, is_receptor=False):
            cmd.delete(ligand_obj)
            return {'success': False, 'error': 'Ligand export failed'}
    except ValueError as e:
        cmd.delete(ligand_obj)
        return {'success': False, 'error': str(e)}
    cmd.delete(ligand_obj)

    ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
    result = run_vina_docking(receptor_pdbqt, ligand_pdbqt, box_params, output_dir,
                              ligand_name, exhaustiveness, num_modes)

    if result['success']:
        return {'success': True, 'output_dir': output_dir, 'results': [result]}
    return result


# ==================== 批量对接功能 ====================
def get_ligand_files(ligand_path: str) -> List[str]:
    """获取配体File列表（支持单File或File夹）"""
    SUPPORTED_EXTS = {'.mol2', '.sdf', '.pdbqt', '.mol', '.pdb'}

    if os.path.isfile(ligand_path):
        return [ligand_path]
    elif os.path.isdir(ligand_path):
        files = []
        for f in os.listdir(ligand_path):
            if os.path.splitext(f)[1].lower() in SUPPORTED_EXTS:
                files.append(os.path.join(ligand_path, f))
        return sorted(files)
    return []


def batch_docking(obj_name: str, ligand_path: str, box_params: Dict[str, float],
                  output_dir: Optional[str] = None, exhaustiveness: int = 8,
                  num_modes: int = 9, remove_selection: Optional[str] = None,
                  progress_callback: Optional[Callable[[int, int, str], None]] = None) -> Dict[str, Any]:
    """
    批量对接：支持单个配体File或配体File夹

    Parameters:
        obj_name: 受体的 PyMOL 对象名
        ligand_path: 配体File或File夹Path
        box_params: 对接盒Parameters
        output_dir: 输出Directory
        exhaustiveness: Vina Search精degrees
        num_modes: 输出构象Count
        remove_selection: 要Remove的共晶配体 selection
        progress_callback: 进degrees回调Function (current, total, ligand_name)
    """
    print(f"[vina_integration] batch_docking Start，ModuleVersion: {_VINA_MODULE_VERSION}")

    if not check_vina_available():
        return {'success': False, 'error': 'Vina not available'}

    ligand_files = get_ligand_files(ligand_path)
    if not ligand_files:
        return {'success': False, 'error': 'No ligand files found'}

    try:
        objs = cmd.get_names("objects")
    except AttributeError:
        objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []

    if obj_name not in objs:
        return {'success': False, 'error': f"Object '{obj_name}' not found"}

    if output_dir is None:
        output_dir = f"{obj_name}_batch_docking"
    os.makedirs(output_dir, exist_ok=True)

    # Export受体 PDBQT（只做一次）
    receptor_pdbqt = os.path.join(output_dir, f"{obj_name}_receptor.pdbqt")
    normalized_remove_sel = _normalize_remove_selection(obj_name, remove_selection)
    receptor_sel = f"({obj_name})"
    if normalized_remove_sel and str(normalized_remove_sel).strip():
        receptor_sel += f" and not ({normalized_remove_sel})"

    print(f"[vina_integration] Export受体: {receptor_pdbqt}")
    try:
        if not export_to_pdbqt(obj_name, receptor_sel, receptor_pdbqt, is_receptor=True):
            return {'success': False, 'error': 'Receptor export failed'}
    except ValueError as e:
        return {'success': False, 'error': str(e)}

    results = []
    total = len(ligand_files)

    for i, lig_file in enumerate(ligand_files):
        ligand_name = os.path.splitext(os.path.basename(lig_file))[0]

        if progress_callback:
            progress_callback(i + 1, total, ligand_name)

        # 转换配体为 PDBQT
        ligand_pdbqt = os.path.join(output_dir, f"{ligand_name}.pdbqt")
        # 如果已经是 pdbqt 格式
        if lig_file.lower().endswith(".pdbqt"):
            # 检查源File和目标File是否相同，避免 shutil.copy 报错
            if os.path.abspath(lig_file) != os.path.abspath(ligand_pdbqt):
                shutil.copy(lig_file, ligand_pdbqt)
            # 如果相同，直接using原FilePath
            else:
                ligand_pdbqt = lig_file
        elif not convert_ligand_to_pdbqt(lig_file, ligand_pdbqt):
            # 尝试通过 PyMOL 转换
            ligand_obj = f"temp_lig_{i}"
            try:
                cmd.load(lig_file, ligand_obj)
                if not export_to_pdbqt(ligand_obj, ligand_obj, ligand_pdbqt, is_receptor=False):
                    cmd.delete(ligand_obj)
                    results.append({'success': False, 'error': 'Conversion failed', 'ligand': ligand_name})
                    continue
                cmd.delete(ligand_obj)
            except ValueError as e:
                try:
                    cmd.delete(ligand_obj)
                except Exception:
                    pass
                results.append({'success': False, 'error': str(e), 'ligand': ligand_name})
                continue
            except Exception as e:
                results.append({'success': False, 'error': str(e), 'ligand': ligand_name})
                continue

        # 运行对接
        result = run_vina_docking(receptor_pdbqt, ligand_pdbqt, box_params, output_dir,
                                  ligand_name, exhaustiveness, num_modes)
        results.append(result)

    # 按亲和力Sort
    successful = [r for r in results if r.get('success')]
    successful.sort(key=lambda x: x.get('affinity', 0))

    # 生成汇总 CSV
    csv_path = os.path.join(output_dir, "docking_results.csv")
    with open(csv_path, 'w') as f:
        f.write("Ligand,Affinity (kcal/mol),Output File\n")
        for r in successful:
            f.write(f"{r['ligand']},{r['affinity']:.2f},{r.get('output_pdbqt', '')}\n")

    return {
        'success': True,
        'output_dir': output_dir,
        'results': results,
        'csv_path': csv_path,
        'total': total,
        'successful': len(successful),
        'failed': total - len(successful)
    }


if __name__ == "__main__":
    print("[vina_integration] This is a PyMOL plugin module")