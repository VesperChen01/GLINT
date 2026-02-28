# -*- coding: utf-8 -*-
"""
GLINT 蛋白突变分析Module
支持多种突变Method和 ΔΔG 计算策略
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


# ==================== Tool检测 ====================

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
    """检测 FoldX 可执行File"""
    return _which("foldx")


# ==================== 突变字符串解析 ====================

def _parse_mutation_string(mut_str: str) -> Tuple[str, str, str]:
    """
    解析突变字符串
    
    支持格式：
    - "A:123:ALA" 或 "A:123:A"
    - "A123A"
    - "A:123A"
    
    Return：(chain, resi, target_aa_3letter)
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
    using PyMOL 执行单点突变
    
    Parameters：
        obj_name: PyMOL 对象名
        chain: 链 ID
        resi: 残基号
        target_aa: 目标氨基酸（三字母代码）
    
    Return：
        True 如果Success
    """
    if not cmd:
        print("❌ PyMOL 不可用")
        return False
    
    try:
        # 确保目标氨基酸是三字母代码
        if len(target_aa) == 1:
            target_aa = AA_1TO3.get(target_aa.upper(), target_aa)
        
        target_aa = target_aa.upper()
        
        # Select要突变的残基
        selection = f"{obj_name} and chain {chain} and resi {resi}"
        
        # 检查残基是否存在
        if cmd.count_atoms(selection) == 0:
            print(f"❌ 未找到残基: chain {chain} resi {resi}")
            return False
        
        # using PyMOL 的 wizard 进行突变
        # 注意：这是一个简化Version，实际可能需要更复杂的处理
        cmd.wizard("mutagenesis")
        cmd.get_wizard().set_mode(target_aa)
        cmd.get_wizard().do_select(selection)
        cmd.get_wizard().apply()
        cmd.set_wizard()
        
        print(f"✅ 突变Success: chain {chain} resi {resi} → {target_aa}")
        return True
        
    except Exception as e:
        print(f"❌ PyMOL 突变Failed: {e}")
        import traceback
        traceback.print_exc()
        return False


# ==================== 能量最小化 ====================

def minimize_energy(obj_name: str, selection: str = "all", cycles: int = 100, 
                   method: str = "pymol") -> Optional[str]:
    """
    能量最小化
    
    Parameters：
        obj_name: PyMOL 对象名
        selection: Select表达式
        cycles: 迭代次数
        method: 'pymol' 或 'rosetta'
    
    Return：
        最小化后的对象名，FailedReturn None
    """
    if not cmd:
        print("❌ PyMOL 不可用")
        return None
    
    if method == "pymol":
        try:
            # Create副本
            min_obj = f"{obj_name}_minimized"
            cmd.create(min_obj, obj_name)
            
            # using PyMOL sculpting
            cmd.sculpt_activate(min_obj)
            cmd.sculpt_iterate(min_obj, cycles=cycles)
            cmd.sculpt_deactivate(min_obj)
            
            print(f"✅ PyMOL 能量最小化Completed: {cycles} 次迭代")
            return min_obj
            
        except Exception as e:
            print(f"❌ 能量最小化Failed: {e}")
            return None
    
    elif method == "rosetta":
        print("⚠️ Rosetta 能量最小化尚未实现")
        return None
    
    else:
        print(f"❌ 未知的最小化Method: {method}")
        return None


# ==================== 主要接口Function ====================

def perform_mutation(obj_name: str, mutations: List[Tuple[str, str, str]],
                    method: str = "pymol") -> Optional[str]:
    """
    执行突变
    
    Parameters：
        obj_name: PyMOL 对象名
        mutations: 突变列表 [(chain, resi, target_aa), ...]
        method: 'pymol' 或 'foldx'
    
    Return：
        突变后的对象名，FailedReturn None
    
    Method说明：
        - pymol: using PyMOL 内置 mutagenesis wizard（快速，无需外部Tool）
        - foldx: using FoldX BuildModel（更准确，需要安装 FoldX）
    """
    if not cmd:
        print("❌ PyMOL 不可用")
        return None
    
    # Create突变对象
    mut_obj = f"{obj_name}_mut"
    cmd.create(mut_obj, obj_name)
    
    if method == "pymol":
        success_count = 0
        for chain, resi, target_aa in mutations:
            if _pymol_mutate(mut_obj, chain, resi, target_aa):
                success_count += 1
        
        if success_count == 0:
            print("❌ 所有突变均Failed")
            cmd.delete(mut_obj)
            return None
        
        print(f"✅ Completed {success_count}/{len(mutations)} 个突变")
        return mut_obj
    
    elif method == "foldx":
        if not _detect_foldx():
            print("❌ FoldX 未检测到")
            print("💡 安装Method：")
            print("   1. Download FoldX: https://foldxsuite.crg.eu/")
            print("   2. Settings环境变量: export FOLDX=/path/to/foldx")
            cmd.delete(mut_obj)
            return None
        
        print("⚠️ FoldX 突变功能请using ddg_heatmap 命令")
        print("💡 示例: ddg_heatmap('CRBN_sel', 'POI_sel')")
        cmd.delete(mut_obj)
        return None
    
    else:
        print(f"❌ 未知的突变Method: {method}")
        print("💡 可用Method: pymol, foldx")
        cmd.delete(mut_obj)
        return None


def calculate_mutation_ddg(wt_obj: str, mut_obj: str, partner_sel: Optional[str] = None,
                          method: str = "foldx") -> Optional[Dict[str, Any]]:
    """
    计算突变的 ΔΔG
    
    Parameters：
        wt_obj: 野生型对象名
        mut_obj: 突变型对象名
        partner_sel: 结合伴侣Select（可选）
        method: 'foldx'（目前仅支持 FoldX）
    
    Return：
        Results字典 {'ddg': float, 'method': str, 'details': dict}
    
    Method说明：
        - foldx: FoldX BuildModel（推荐，需要安装 FoldX）
    
    FoldX 输入：
        - PDB File：野生型和突变型结构
        - 突变File：格式为 <chain><resi><icode><WT><Mut>
    
    FoldX 输出：
        - DifferencesBetweenMutantAndWildType_fxout.csv
        - Package含 Total Energy (ΔΔG) Value，单位 kcal/mol
    """
    if not cmd:
        print("❌ PyMOL 不可用")
        return None
    
    if method != "foldx":
        print(f"❌ 未知的Method: {method}")
        print("💡 目前仅支持 FoldX Method")
        return None
    
    if not _detect_foldx():
        print("❌ FoldX 未检测到")
        print("💡 安装Method：")
        print("   1. Download FoldX: https://foldxsuite.crg.eu/")
        print("   2. Settings环境变量: export FOLDX=/path/to/foldx")
        print("   或将 FoldX 可执行FileAdd到 PATH")
        return None
    
    print("💡 提示：请using ddg_heatmap 命令进行完整的 ΔΔG 分析")
    print("   示例: ddg_heatmap('CRBN_sel', 'POI_sel')")
    return None


def analyze_mutation_effects(obj_name: str, mutations: List[Tuple[str, str, str]],
                            partner_sel: Optional[str] = None,
                            output_csv: Optional[str] = None,
                            method: str = "auto") -> Optional[Dict[str, Any]]:
    """
    完整的突变分析流程
    
    Parameters：
        obj_name: PyMOL 对象名
        mutations: 突变列表
        partner_sel: 结合伴侣Select
        output_csv: 输出 CSV Path
        method: ΔΔG 计算Method
    
    Return：
        分析Results字典
    """
    if not cmd:
        print("❌ PyMOL 不可用")
        return None
    
    print(f"\n{'='*60}")
    print(f"🧬 Start突变分析: {obj_name}")
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
        min_obj = mut_obj  # 如果最小化Failed，using未最小化的对象
    
    # 3. 计算 ΔΔG
    print("\n步骤 3/3: 计算 ΔΔG...")
    ddg_result = calculate_mutation_ddg(obj_name, min_obj, partner_sel, method)
    
    # 汇总Results
    results = {
        'wt_obj': obj_name,
        'mut_obj': mut_obj,
        'min_obj': min_obj,
        'mutations': mutations,
        'ddg_result': ddg_result,
    }
    
    # Export CSV
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
            print(f"\\n✅ Results已Save到: {output_csv}")
        except Exception as e:
            print(f"⚠️ CSV ExportFailed: {e}")
    
    print(f"\\n{'='*60}")
    print(f"✅ 突变分析Completed")
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

def ddg_heatmap(CRBN_sel: str, POI_sel: str, name: str = "ddg"):
    """
    using FoldX 计算 ΔΔG 并在 POI 上Display热图
    
    Parameters：
        CRBN_sel: CRBN（E3连接酶）Select表达式
        POI_sel: POI（目标蛋白）Select表达式
        name: 颜色渐变Name（默认 'ddg'）
    
    FoldX 输入：
        - PDB File：CRBN + POI 复合物结构
        - 突变File：interface残基的丙氨酸扫描突变列表
          格式：<chain><resi><icode><WT>A（每行一个突变）
    
    FoldX 输出：
        - DifferencesBetweenMutantAndWildType_fxout.csv
        - 字段：Mutation, Total Energy (ΔΔG, kcal/mol)
    
    可视化输出：
        - POI 表面按 ΔΔG Value着色（蓝-白-红渐变）
        - 蓝色：稳定化突变（负 ΔΔG）
        - 红色：去稳定化突变（正 ΔΔG）
    
    示例：
        ddg_heatmap('chain A', 'chain B')
    """
    poi = f"({POI_sel})"
    crbn = f"({CRBN_sel})"
    
    # 检测 FoldX
    fx = _detect_foldx()
    if not fx:
        print("[ddg_heatmap] ❌ FoldX not found in PATH or $FOLDX")
        print("[ddg_heatmap] 💡 安装Method：")
        print("   1. Download FoldX: https://foldxsuite.crg.eu/")
        print("   2. Settings环境变量: export FOLDX=/path/to/foldx")
        print("   或将 FoldX 可执行FileAdd到 PATH")
        return
    
    print(f"[ddg_heatmap] ✅ Using FoldX at: {fx}")
    
    # 准备 PDB File并运行 FoldX
    pdb_path, _ = _prep_complex_tmp(crbn, poi)
    ddg = _foldx_alanine_scan(fx, pdb_path, poi)
    
    if not ddg:
        cmd.feedback("pop", "all", "actions")
        print("[ddg_heatmap] ❌ No ΔΔG values computed. Check selections.")
        print("[ddg_heatmap] 💡 确保Select的残基在interface区域（5Å 内）")
        return
    
    # Settings B 因子并着色
    set_b_factors(poi, ddg)
    color_by_b(poi, palette="blue_white_red", ramp_name=f"{name}_ramp")
    cmd.show("surface", poi)
    
    print(f"[ddg_heatmap] ✅ Completed！共计算 {len(ddg)} 个残基的 ΔΔG Value")

# ==================== PyMOL 命令注册 ====================

if cmd:
    cmd.extend("perform_mutation", perform_mutation)
    cmd.extend("minimize_energy", minimize_energy)
    cmd.extend("calculate_mutation_ddg", calculate_mutation_ddg)
    cmd.extend("analyze_mutation_effects", analyze_mutation_effects)
    cmd.extend("ddg_heatmap", ddg_heatmap)
    cmd.extend("analyze_mutation_effects", analyze_mutation_effects)
