# -*- coding: utf-8 -*-
"""
GlueTK 蛋白突变分析模块
支持多种突变方法和 ΔΔG 计算策略
"""

import os
import sys
import tempfile
import subprocess
from typing import List, Tuple, Dict, Optional, Any
import csv
import math

try:
    from pymol import cmd
    from .highlight_residues import set_b_factors, color_by_b
except ImportError:
    cmd = None
    # Mock for testing
    def set_b_factors(*args): pass
    def color_by_b(*args): pass


# ==================== 氨基酸转换表 ====================

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items()}


# ==================== 工具检测 ====================\n
def _which(exe: str) -> Optional[str]:
    path = os.environ.get("FOLDX") if exe.lower() == "foldx" else None
    if path and os.path.isfile(path) and os.access(path, os.X_OK):
        return path
    for p in os.environ.get("PATH", "").split(os.pathsep):
        cand = os.path.join(p, exe)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
        # macOS app bundle
        if exe.lower() == "foldx" and os.path.isdir(os.path.join(p, "FoldX.app")):
            app_bin = os.path.join(p, "FoldX.app", "Contents", "MacOS", "FoldX")
            if os.path.isfile(app_bin) and os.access(app_bin, os.X_OK):
                return app_bin
    return None

def _detect_foldx() -> Optional[str]:
    """检测 FoldX 可执行文件"""
    return _which("foldx")
        return None


def _detect_pyrosetta() -> bool:
    """检测 PyRosetta 是否可用"""
    try:
        import pyrosetta
        return True
    except ImportError:
        return False


# ==================== 突变字符串解析 ====================

def _parse_mutation_string(mut_str: str) -> Tuple[str, str, str]:
    """
    解析突变字符串
    
    支持格式：
    - "A:123:ALA" 或 "A:123:A"
    - "A123A"
    - "A:123A"
    
    返回：(chain, resi, target_aa_3letter)
    """
    mut_str = mut_str.strip()
    
    # 格式 1: A:123:ALA 或 A:123:A
    if mut_str.count(':') == 2:
        parts = mut_str.split(':')
        chain = parts[0].strip()
        resi = parts[1].strip()
        target = parts[2].strip().upper()
        
        # 转换为三字母代码
        if len(target) == 1:
            target = AA_1TO3.get(target, target)
        
        return (chain, resi, target)
    
    # 格式 2: A:123A
    elif ':' in mut_str:
        parts = mut_str.split(':')
        chain = parts[0].strip()
        rest = parts[1].strip()
        
        # 提取残基号和目标氨基酸
        import re
        match = re.match(r'(\d+)([A-Z])', rest)
        if match:
            resi = match.group(1)
            target = match.group(2)
            target = AA_1TO3.get(target, target)
            return (chain, resi, target)
    
    # 格式 3: A123A
    else:
        import re
        match = re.match(r'([A-Z])(\d+)([A-Z])', mut_str)
        if match:
            chain = match.group(1)
            resi = match.group(2)
            target = match.group(3)
            target = AA_1TO3.get(target, target)
            return (chain, resi, target)
    
    raise ValueError(f"无法解析突变字符串: {mut_str}")


# ==================== PyMOL 突变实现 ====================

def _pymol_mutate(obj_name: str, chain: str, resi: str, target_aa: str) -> bool:
    """
    使用 PyMOL 执行单点突变
    
    参数：
        obj_name: PyMOL 对象名
        chain: 链 ID
        resi: 残基号
        target_aa: 目标氨基酸（三字母代码）
    
    返回：
        True 如果成功
    """
    if not cmd:
        print("❌ PyMOL 不可用")
        return False
    
    try:
        # 确保目标氨基酸是三字母代码
        if len(target_aa) == 1:
            target_aa = AA_1TO3.get(target_aa.upper(), target_aa)
        
        target_aa = target_aa.upper()
        
        # 选择要突变的残基
        selection = f"{obj_name} and chain {chain} and resi {resi}"
        
        # 检查残基是否存在
        if cmd.count_atoms(selection) == 0:
            print(f"❌ 未找到残基: chain {chain} resi {resi}")
            return False
        
        # 使用 PyMOL 的 wizard 进行突变
        # 注意：这是一个简化版本，实际可能需要更复杂的处理
        cmd.wizard("mutagenesis")
        cmd.get_wizard().set_mode(target_aa)
        cmd.get_wizard().do_select(selection)
        cmd.get_wizard().apply()
        cmd.set_wizard()
        
        print(f"✅ 突变成功: chain {chain} resi {resi} → {target_aa}")
        return True
        
    except Exception as e:
        print(f"❌ PyMOL 突变失败: {e}")
        import traceback
        traceback.print_exc()
        return False


# ==================== 能量最小化 ====================

def minimize_energy(obj_name: str, selection: str = "all", cycles: int = 100, 
                   method: str = "pymol") -> Optional[str]:
    """
    能量最小化
    
    参数：
        obj_name: PyMOL 对象名
        selection: 选择表达式
        cycles: 迭代次数
        method: 'pymol' 或 'rosetta'
    
    返回：
        最小化后的对象名，失败返回 None
    """
    if not cmd:
        print("❌ PyMOL 不可用")
        return None
    
    if method == "pymol":
        try:
            # 创建副本
            min_obj = f"{obj_name}_minimized"
            cmd.create(min_obj, obj_name)
            
            # 使用 PyMOL sculpting
            cmd.sculpt_activate(min_obj)
            cmd.sculpt_iterate(min_obj, cycles=cycles)
            cmd.sculpt_deactivate(min_obj)
            
            print(f"✅ PyMOL 能量最小化完成: {cycles} 次迭代")
            return min_obj
            
        except Exception as e:
            print(f"❌ 能量最小化失败: {e}")
            return None
    
    elif method == "rosetta":
        print("⚠️ Rosetta 能量最小化尚未实现")
        return None
    
    else:
        print(f"❌ 未知的最小化方法: {method}")
        return None


# ==================== 主要接口函数 ====================

def perform_mutation(obj_name: str, mutations: List[Tuple[str, str, str]], 
                    method: str = "pymol") -> Optional[str]:
    """
    执行突变
    
    参数：
        obj_name: PyMOL 对象名
        mutations: 突变列表 [(chain, resi, target_aa), ...]
        method: 'pymol', 'foldx', 或 'pyrosetta'
    
    返回：
        突变后的对象名，失败返回 None
    
    方法说明：
        - pymol: 使用 PyMOL 内置 mutagenesis wizard（快速，无需外部工具）
        - foldx: 使用 FoldX BuildModel（更准确，需要安装 FoldX）
        - pyrosetta: 使用 PyRosetta（最准确，跨平台，需要许可证）
    """
    if not cmd:
        print("❌ PyMOL 不可用")
        return None
    
    # 创建突变对象
    mut_obj = f"{obj_name}_mut"
    cmd.create(mut_obj, obj_name)
    
    if method == "pymol":
        success_count = 0
        for chain, resi, target_aa in mutations:
            if _pymol_mutate(mut_obj, chain, resi, target_aa):
                success_count += 1
        
        if success_count == 0:
            print("❌ 所有突变均失败")
            cmd.delete(mut_obj)
            return None
        
        print(f"✅ 完成 {success_count}/{len(mutations)} 个突变")
        return mut_obj
    
    elif method == "foldx":
        print("⚠️ FoldX 突变尚未实现")
        print("💡 提示：FoldX 集成将在后续版本中添加")
        cmd.delete(mut_obj)
        return None
    
    elif method == "pyrosetta":
        if not _detect_pyrosetta():
            print("❌ PyRosetta 未安装")
            print("💡 安装方法：")
            print("   1. 申请许可证（学术免费）: https://www.pyrosetta.org/")
            print("   2. 使用 Conda 安装（推荐）:")
            print("      conda install -c https://USERNAME:PASSWORD@conda.graylab.jhu.edu pyrosetta")
            print("   3. 或使用 pip 安装:")
            print("      pip install pyrosetta-*.whl")
            cmd.delete(mut_obj)
            return None
        
        print("⚠️ PyRosetta 突变尚未实现")
        print("💡 这将在后续版本中添加")
        cmd.delete(mut_obj)
        return None
    
    else:
        print(f"❌ 未知的突变方法: {method}")
        print("💡 可用方法: pymol, foldx, pyrosetta")
        cmd.delete(mut_obj)
        return None


def calculate_mutation_ddg(wt_obj: str, mut_obj: str, partner_sel: Optional[str] = None,
                          method: str = "auto") -> Optional[Dict[str, Any]]:
    """
    计算突变的 ΔΔG
    
    参数：
        wt_obj: 野生型对象名
        mut_obj: 突变型对象名
        partner_sel: 结合伴侣选择（可选）
        method: 'auto', 'foldx', 'pyrosetta'
    
    返回：
        结果字典 {'ddg': float, 'method': str, 'details': dict}
    
    方法说明：
        - auto: 自动选择最佳可用方法（FoldX > PyRosetta）
        - foldx: FoldX BuildModel（推荐，中等准确性，跨平台）
        - pyrosetta: PyRosetta 计算（高准确性，需要许可证）
    """
    if not cmd:
        print("❌ PyMOL 不可用")
        return None
    
    # 自动选择方法
    if method == "auto":
        if _detect_foldx():
            method = "foldx"
        elif _detect_pyrosetta():
            method = "pyrosetta"
        else:
            print("❌ 未检测到 FoldX 或 PyRosetta")
            print("� 请安装其中一个：")
            print("   - FoldX: https://foldxsuite.crg.eu/")
            print("   - PyRosetta: https://www.pyrosetta.org/")
            return None
        
        print(f"🔍 自动选择方法: {method}")
    
    if method == "foldx":
        if not _detect_foldx():
            print("❌ FoldX 未检测到")
            print("💡 安装方法：")
            print("   1. 下载 FoldX: https://foldxsuite.crg.eu/")
            print("   2. 设置环境变量: export FOLDX=/path/to/foldx")
            return None
        
        print("⚠️ FoldX ΔΔG 计算尚未实现")
        print("💡 提示：可以使用现有的 ddg_heatmap 命令")
        print("   示例: ddg_heatmap('CRBN_sel', 'POI_sel', method='foldx')")
        return None
    
    elif method == "pyrosetta":
        if not _detect_pyrosetta():
            print("❌ PyRosetta 未安装")
            print("💡 安装方法：")
            print("   1. 申请许可证（学术免费）: https://www.pyrosetta.org/")
            print("   2. 使用 Conda 安装（推荐）:")
            print("      conda install -c https://USERNAME:PASSWORD@conda.graylab.jhu.edu pyrosetta")
            print("   3. 或使用 pip 安装:")
            print("      pip install pyrosetta-*.whl")
            print("   注意：Windows 用户现已支持原生安装（无需 WSL）")
            return None
        
        print("⚠️ PyRosetta ΔΔG 计算尚未实现")
        print("💡 这将在后续版本中添加")
        return None
    
    else:
        print(f"❌ 未知的方法: {method}")
        print("💡 可用方法: auto, foldx, pyrosetta")
        return None


def analyze_mutation_effects(obj_name: str, mutations: List[Tuple[str, str, str]],
                            partner_sel: Optional[str] = None,
                            output_csv: Optional[str] = None,
                            method: str = "auto") -> Optional[Dict[str, Any]]:
    """
    完整的突变分析流程
    
    参数：
        obj_name: PyMOL 对象名
        mutations: 突变列表
        partner_sel: 结合伴侣选择
        output_csv: 输出 CSV 路径
        method: ΔΔG 计算方法
    
    返回：
        分析结果字典
    """
    if not cmd:
        print("❌ PyMOL 不可用")
        return None
    
    print(f"\n{'='*60}")
    print(f"🧬 开始突变分析: {obj_name}")
    print(f"{'='*60}\n")
    
    # 1. 执行突变
    print("步骤 1/3: 执行突变...")
    mut_obj = perform_mutation(obj_name, mutations, method="pymol")
    if not mut_obj:
        return None
    
    # 2. 能量最小化
    print("\n步骤 2/3: 能量最小化...")
    min_obj = minimize_energy(mut_obj, cycles=100)
    if not min_obj:
        min_obj = mut_obj  # 如果最小化失败，使用未最小化的对象
    
    # 3. 计算 ΔΔG
    print("\n步骤 3/3: 计算 ΔΔG...")
    ddg_result = calculate_mutation_ddg(obj_name, min_obj, partner_sel, method)
    
    # 汇总结果
    results = {
        'wt_obj': obj_name,
        'mut_obj': mut_obj,
        'min_obj': min_obj,
        'mutations': mutations,
        'ddg_result': ddg_result,
    }
    
    # 导出 CSV
    if output_csv and ddg_result:
        try:
            import csv
            with open(output_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Chain', 'Residue', 'Target_AA', 'ddG', 'Method'])
                for (chain, resi, target_aa) in mutations:
                    writer.writerow([chain, resi, target_aa, 
                                   ddg_result.get('ddg', 'N/A'),
                                   ddg_result.get('method', 'N/A')])
            print(f"\\n✅ 结果已保存到: {output_csv}")
        except Exception as e:
            print(f"⚠️ CSV 导出失败: {e}")
    
    print(f"\\n{'='*60}")
    print(f"✅ 突变分析完成")
    print(f"{'='*60}\\n")
    
    return results


# ==================== ΔΔG Heatmap (Legacy/Quick Mode) ====================

def _prep_complex_tmp(crbn_sel: str, poi_sel: str) -> Tuple[str, Dict[Tuple[str,str,str], Tuple[str,str,str]]]:
    # Export CRBN+POI complex to PDB with consistent chain/resi mapping
    tmp = tempfile.mkdtemp(prefix="pymol_crbn_")
    pdb_path = os.path.join(tmp, "complex.pdb")
    # Create a temporary object combining selections
    tmp_obj = "__crbn_tmp_complex__"
    cmd.delete(tmp_obj)
    cmd.create(tmp_obj, f"({crbn_sel}) or ({poi_sel})")
    cmd.sort(tmp_obj)
    cmd.save(pdb_path, tmp_obj)
    # Build atom-to-resi mapping for later results
    mapping = {}
    model = cmd.get_model(poi_sel)
    for a in model.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        mapping[key] = key
    return pdb_path, mapping

def _foldx_alanine_scan(foldx_bin: str, pdb_path: str, poi_sel: str) -> Dict[Tuple[str,str,str], float]:
    # Interface residues on POI
    poi_iface = cmd.get_model(f"byres ({poi_sel} within 5.0 of not {poi_sel})")
    # Mutation list format for FoldX: "<PDB>; <chain><resi><icode><WT><Mut>"
    muts = []
    seen = set()
    for a in poi_iface.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        if key in seen:
            continue
        seen.add(key)
        wt1 = AA_3TO1.get(a.resn.upper(), 'X')
        if wt1 == "X":
            continue
        code = (a.chain or "A") + a.resi + (getattr(a, 'q', a.icode) or " ") + wt1 + "A"
        muts.append(code)
    if not muts:
        return {}

    tmpdir = os.path.dirname(pdb_path)
    mut_file = os.path.join(tmpdir, "mutations.txt")
    with open(mut_file, "w") as f:
        for m in muts:
            f.write(m + "\\n")

    # Run FoldX BuildModel
    cmd_line = [foldx_bin, "--command=BuildModel", f"--pdb={os.path.basename(pdb_path)}", f"--pdb-dir={tmpdir}", f"--mutant-file={mut_file}"]
    try:
        subprocess.run(cmd_line, cwd=tmpdir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception:
        return {}

    # Parse DifferencesBetweenMutantAndWildType_fxout.csv if exists
    out_csv = os.path.join(tmpdir, "DifferencesBetweenMutantAndWildType_fxout.csv")
    ddg = {}
    if os.path.isfile(out_csv):
        with open(out_csv, newline='') as f:
            reader = csv.DictReader(f, delimiter=';')
            for row in reader:
                # Name like: <pdb>;<mutation> e.g., complex.pdb;A123A->A
                mut = row.get('Mutation', '') or row.get('mutation', '')
                # Extract chain/resi
                if len(mut) >= 5:
                    chain = mut[0]
                    resi = mut[1:4].strip()
                    icode = ""
                    key = (chain, resi, icode)
                    try:
                        ddg[key] = float(row.get('Total Energy', row.get('Total energy', '0')))
                    except Exception:
                        pass
    return ddg

def _asa_ddg_proxy(poi_sel: str, complex_sel: str, scale: float = 0.025) -> Dict[Tuple[str,str,str], float]:
    # Use PyMOL get_area per residue in monomer vs complex; ΔASA scaled to ddG
    poi = f"({poi_sel})"
    # Build list of residues
    model = cmd.get_model(poi)
    residues = []
    seen = set()
    for a in model.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        if key not in seen:
            seen.add(key)
            residues.append(key)

    ddg = {}
    # Ensure surface areas are computed with dot settings
    prev_dot_solvent = cmd.get("dot_solvent")
    prev_dot_density = cmd.get("dot_density")
    cmd.set("dot_solvent", 1)
    cmd.set("dot_density", 3)

    for (ch, resi, icode) in residues:
        sel_res = f"{poi} and chain {ch} and resi {resi}"
        asa_complex = cmd.get_area(sel_res, load_b=0, state=1)
        # To get monomer ASA, duplicate POI only into temp object and measure
        tmp_obj = "__poi_tmp__"
        cmd.delete(tmp_obj)
        cmd.create(tmp_obj, sel_res)
        asa_monomer = cmd.get_area(tmp_obj, load_b=0, state=1)
        cmd.delete(tmp_obj)
        dASA = max(0.0, asa_monomer - asa_complex)
        ddg[(ch, resi, icode)] = -scale * dASA  # burial stabilizes binding (negative)

    # Restore settings
    cmd.set("dot_solvent", prev_dot_solvent)
    cmd.set("dot_density", prev_dot_density)
    return ddg

def ddg_heatmap(CRBN_sel: str, POI_sel: str, method: str = "auto", name: str = "ddg", scale: float = 0.025):
    """
    Color POI by ΔΔG. method: auto|foldx|asa. Stores value to b-factor and colors by spectrum.
    """
    poi = f"({POI_sel})"
    crbn = f"({CRBN_sel})"
    # Try FoldX if auto
    ddg = {}
    use_foldx = False
    if method in ("auto", "foldx"):
        fx = _detect_foldx()
        if fx:
            use_foldx = True
            print(f"[ddg_heatmap] ✅ Using FoldX at: {fx}")
            pdb_path, _ = _prep_complex_tmp(crbn, poi)
            ddg = _foldx_alanine_scan(fx, pdb_path, poi)
        elif method == "foldx":
            # User explicitly requested FoldX but it's not available
            print("[ddg_heatmap] ❌ FoldX not found in PATH or $FOLDX")
            print("[ddg_heatmap] 💡 Download FoldX from: https://foldxsuite.crg.eu/")
            if method == "foldx":  # Don't fallback if explicitly requested
                return
    
    if (not ddg) and method in ("auto", "asa"):
        # ASA proxy fallback
        print("[ddg_heatmap] ⚠️ Using ASA-based proxy (ΔΔG approximation, not publication quality)")
        print("[ddg_heatmap] 💡 For accurate ΔΔG values, install FoldX: https://foldxsuite.crg.eu/")
        ddg = _asa_ddg_proxy(poi, f"{crbn} or {poi}", scale=scale)

    if not ddg:
        cmd.feedback("pop", "all", "actions")
        print("[ddg_heatmap] ❌ No ΔΔG values computed. Check selections.")
        return

    set_b_factors(poi, ddg)
    color_by_b(poi, palette="blue_white_red", ramp_name=f"{name}_ramp")
    cmd.show("surface", poi)

# ==================== PyMOL 命令注册 ====================\n
if cmd:
    cmd.extend("perform_mutation", perform_mutation)
    cmd.extend("minimize_energy", minimize_energy)
    cmd.extend("calculate_mutation_ddg", calculate_mutation_ddg)
    cmd.extend("analyze_mutation_effects", analyze_mutation_effects)
    cmd.extend("ddg_heatmap", ddg_heatmap)
    cmd.extend("analyze_mutation_effects", analyze_mutation_effects)
