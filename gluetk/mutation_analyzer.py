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

try:
    from pymol import cmd
except ImportError:
    cmd = None


# ==================== 氨基酸转换表 ====================

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items()}


# ==================== 工具检测 ====================

def _detect_foldx() -> Optional[str]:
    """检测 FoldX 可执行文件"""
    # 复用 pymol_crbn_tools.py 中的检测逻辑
    try:
        from .pymol_crbn_tools import _detect_foldx as detect
        return detect()
    except ImportError:
        # 简化版检测
        foldx_path = os.environ.get("FOLDX")
        if foldx_path and os.path.isfile(foldx_path):
            return foldx_path
        
        # 检查 PATH
        for path in os.environ.get("PATH", "").split(os.pathsep):
            candidate = os.path.join(path, "foldx")
            if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                return candidate
        
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
            print(f"\n✅ 结果已保存到: {output_csv}")
        except Exception as e:
            print(f"⚠️ CSV 导出失败: {e}")
    
    print(f"\n{'='*60}")
    print(f"✅ 突变分析完成")
    print(f"{'='*60}\n")
    
    return results


# ==================== PyMOL 命令注册 ====================

if cmd:
    cmd.extend("perform_mutation", perform_mutation)
    cmd.extend("minimize_energy", minimize_energy)
    cmd.extend("calculate_mutation_ddg", calculate_mutation_ddg)
    cmd.extend("analyze_mutation_effects", analyze_mutation_effects)
