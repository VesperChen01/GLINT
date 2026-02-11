# -*- coding: utf-8 -*-
"""
interaction_2d_plot.py
蛋白质-配体 2D 相互作用图生成器 (Discovery Studio 风格)

特性：
- 基于 3D 坐标的精确原子/环映射
- 使用 RDKit CoordGen 实现先进的 2D 布局
- 支持大环分子的优雅展示
- 统一的配色方案 (Schrödinger 风格)
- 精确的 Pi-相互作用环中心连接
- 与 interaction_analyzer.py 的 3D 可视化完全一致的配色和类型定义
"""

import os
import csv
import tempfile
import math
import numpy as np
import re
from PIL import Image
import matplotlib
matplotlib.use('Agg')  # 非交互式后端，避免与 PyMOL Qt 事件循环死锁
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    # 尝试导入 rdDetermineBonds 模块 (RDKit 2022.09+)，用于从 3D 坐标推断键级
    try:
        from rdkit.Chem import rdDetermineBonds
        HAS_DETERMINE_BONDS = True
    except ImportError:
        HAS_DETERMINE_BONDS = False
except ImportError:
    print("[GLINT] RDKit required for 2D plotting.")
    HAS_DETERMINE_BONDS = False

# Import unified color scheme
try:
    from .color_scheme import (
        INTERACTION_COLORS_HEX,
        INTERACTION_LINE_STYLES,
        RESIDUE_COLORS_HEX,
    )
except ImportError:
    from color_scheme import (
        INTERACTION_COLORS_HEX,
        INTERACTION_LINE_STYLES,
        RESIDUE_COLORS_HEX,
    )

# ============================================================================
# 残基类型分类 (与 3D 视图一致)
# ============================================================================
RESIDUE_TYPES = {
    'ALA': 'hydrophobic', 'VAL': 'hydrophobic', 'LEU': 'hydrophobic',
    'ILE': 'hydrophobic', 'MET': 'hydrophobic', 'PHE': 'hydrophobic',
    'TRP': 'hydrophobic', 'PRO': 'hydrophobic',
    'SER': 'polar', 'THR': 'polar', 'ASN': 'polar', 'GLN': 'polar',
    'TYR': 'polar', 'CYS': 'polar',
    'ASP': 'negative', 'GLU': 'negative',
    'LYS': 'positive', 'ARG': 'positive', 'HIS': 'positive',
    'GLY': 'nonpolar',
}

# Discovery Studio 风格颜色方案 (残基气泡) - 使用统一配色
DS_RESIDUE_STYLE = {
    'hydrophobic': {
        'facecolor': RESIDUE_COLORS_HEX['hydrophobic']['face'],
        'edgecolor': RESIDUE_COLORS_HEX['hydrophobic']['edge'],
        'textcolor': RESIDUE_COLORS_HEX['hydrophobic']['text']
    },
    'nonpolar': {
        'facecolor': RESIDUE_COLORS_HEX['nonpolar']['face'],
        'edgecolor': RESIDUE_COLORS_HEX['nonpolar']['edge'],
        'textcolor': RESIDUE_COLORS_HEX['nonpolar']['text']
    },
    'polar': {
        'facecolor': RESIDUE_COLORS_HEX['polar']['face'],
        'edgecolor': RESIDUE_COLORS_HEX['polar']['edge'],
        'textcolor': RESIDUE_COLORS_HEX['polar']['text']
    },
    'negative': {
        'facecolor': RESIDUE_COLORS_HEX['negative']['face'],
        'edgecolor': RESIDUE_COLORS_HEX['negative']['edge'],
        'textcolor': RESIDUE_COLORS_HEX['negative']['text']
    },
    'positive': {
        'facecolor': RESIDUE_COLORS_HEX['positive']['face'],
        'edgecolor': RESIDUE_COLORS_HEX['positive']['edge'],
        'textcolor': RESIDUE_COLORS_HEX['positive']['text']
    },
}

# ============================================================================
# 统一配色方案 (使用 color_scheme.py 中的定义)
# Schrödinger 风格 - Hex 颜色值
# ============================================================================
UNIFIED_INTERACTION_COLORS = INTERACTION_COLORS_HEX

# 相互作用连线样式 (使用统一的 INTERACTION_LINE_STYLES)
INTERACTION_LINE_STYLE = INTERACTION_LINE_STYLES

# ============================================================================
# Interaction type mapping (unified English names, compatible with 3D view)
# ============================================================================
INTERACTION_TYPE_MAP = {
    # Primary English names -> internal type
    'Hydrogen Bond': 'hbond',
    'Salt Bridge': 'salt',
    'Pi-Pi Stacking': 'pipi',
    'Pi-Cation': 'pication',
    'Hydrophobic': 'hydrophobic',
    'Halogen Bond': 'halogen',
    'Metal Coordination': 'metal',
    'Water Bridge': 'water',
    'Weak Polar Contact': 'other',
    'Disulfide Bond': 'other',
    'van der Waals': 'other',
    # Legacy Chinese mappings (backward compatibility)
    '氢键': 'hbond',
    '盐桥': 'salt',
    'π–π 堆积': 'pipi',
    'π–阳离子相互作用': 'pication',
    '疏水相互作用': 'hydrophobic',
    '卤素键': 'halogen',
    '金属配位': 'metal',
    '水桥': 'water',
    '疏水接触': 'hydrophobic',
    'π-π堆积': 'pipi',
    '阳离子-π': 'pication',
    '弱极性接触': 'other',
    'polar contact': 'other',
    # Lowercase English -> internal type
    'hydrogen bond': 'hbond',
    'hbond': 'hbond',
    'h-bond': 'hbond',
    'salt bridge': 'salt',
    'saltbridge': 'salt',
    'pi-pi': 'pipi',
    'pi-pi stacking': 'pipi',
    'stacking': 'pipi',
    'pi-cation': 'pication',
    'cation-pi': 'pication',
    'hydrophobic': 'hydrophobic',
    'alkyl': 'hydrophobic',
    'halogen': 'halogen',
    'halogen bond': 'halogen',
    'metal': 'metal',
    'metal coordination': 'metal',
    'metalcoord': 'metal',
    'water bridge': 'water',
    'waterbridge': 'water',
}

def get_residue_style(resname):
    """
    获取残基的显示样式
    
    参数:
        resname: 残基名称 (3字母代码)
    
    返回:
        (style_dict, res_type): 样式字典和残基类型
    """
    resname = resname.upper()[:3]
    res_type = RESIDUE_TYPES.get(resname, 'polar')
    return DS_RESIDUE_STYLE[res_type], res_type


def normalize_interaction_type(itype):
    """
    将相互作用类型标准化为内部类型名
    
    支持中英文混合输入，与 3D 视图 (interaction_analyzer.py) 完全一致
    
    参数:
        itype: 相互作用类型字符串 (中文或英文)
    
    返回:
        str: 标准化的内部类型名 ('hbond', 'salt', 'pipi', 等)
    """
    if not itype:
        return 'other'
    
    itype_lower = itype.lower().strip()
    
    # 首先尝试精确匹配
    if itype_lower in INTERACTION_TYPE_MAP:
        return INTERACTION_TYPE_MAP[itype_lower]
    
    # Then try partial matching (by priority order)
    # Hydrogen bond
    if 'hydrogen' in itype_lower or 'hbond' in itype_lower or 'h-bond' in itype_lower or '氢键' in itype:
        return 'hbond'
    # Salt bridge
    if 'salt' in itype_lower or '盐桥' in itype:
        return 'salt'
    # Pi-Pi stacking
    if 'pi-pi' in itype_lower or 'pipi' in itype_lower or 'stacking' in itype_lower or 'π–π' in itype or 'π-π' in itype:
        return 'pipi'
    # Pi-Cation
    if 'pi-cation' in itype_lower or 'pication' in itype_lower or 'cation-pi' in itype_lower or 'π–阳离子' in itype or 'π-阳离子' in itype:
        return 'pication'
    # Metal coordination
    if 'metal' in itype_lower or '金属' in itype:
        return 'metal'
    # Water bridge
    if 'water' in itype_lower or '水桥' in itype:
        return 'water'
    # Halogen bond
    if 'halogen' in itype_lower or '卤素' in itype:
        return 'halogen'
    # Hydrophobic (check last, as many type names may contain related words)
    if 'hydrophobic' in itype_lower or 'alkyl' in itype_lower or '疏水' in itype:
        return 'hydrophobic'
    
    return 'other'


def get_interaction_line_style(itype):
    """
    获取相互作用的连线样式
    
    与 3D 视图 (interaction_analyzer.py visualize_protein_ligand_3d) 使用相同的配色方案
    
    参数:
        itype: 相互作用类型字符串 (中文或英文)
    
    返回:
        dict: 包含 color, linewidth, linestyle, label 的样式字典
    """
    normalized_type = normalize_interaction_type(itype)
    return INTERACTION_LINE_STYLE.get(normalized_type, INTERACTION_LINE_STYLE['other'])

def parse_residue_label(prot_res):
    parts = prot_res.strip().split()
    if len(parts) >= 2:
        return parts[0][:3].upper(), parts[1]
    m = re.match(r'([A-Za-z]+)(\d+)', prot_res)
    if m: return m.group(1).upper(), m.group(2)
    return prot_res[:3].upper(), ""

def clean_pdb_block(pdb_block, ignore_connect=False):
    """
    纯文本方式清理 PDB 内容：
    1. 保留 CONECT 记录 (关键！用于正确的键连接) - 除非 ignore_connect=True
    2. 移除所有氢原子行 (基于原子名称判断)
    """
    lines = pdb_block.split('\n')
    new_lines = []
    h_atom_serials = set()  # 记录氢原子的序列号，用于清理 CONECT
    
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # PDB 原子名称在 12-16 列 (0-indexed: 12,13,14,15)
            # 严格判定：原子符号在 76-78 列 (如果有)，或者基于名称推断
            # 这里继续沿用名称判断，但增强逻辑
            atom_name = line[12:16].strip()
            element = line[76:78].strip().upper()
            
            is_h = False
            # 1. 优先检查元素符号
            if element == 'H':
                is_h = True
            # 2. 如果没有元素符号，检查名称
            elif not element:
                # 常见氢原子命名模式
                if atom_name.startswith("H") and not (len(atom_name) > 1 and atom_name[1].islower() and atom_name[:2] not in ["He", "Hf", "Hg", "Ho", "Hs"]):
                    # 排除 He, Hf, Hg, Ho, Hs 等元素，其他以 H 开头的通常是氢
                    # 注意：PDB原子名通常是大写。如果有小写可能是特殊情况。
                    # 简单规则：如果是 H 开头，且不是 Hg, Hf 等已知元素的缩写
                    is_h = True
                elif atom_name[0].isdigit() and len(atom_name) > 1 and atom_name[1] == 'H':
                    # 如 1H, 2H
                    is_h = True
            
            if is_h:
                try:
                    serial = int(line[6:11].strip())
                    h_atom_serials.add(serial)
                except (ValueError, TypeError):  # int() 转换可能失败
                    pass
                continue
        new_lines.append(line)
    
    # 彻底忽略 CONECT 记录（如果请求）
    if ignore_connect:
        final_lines = [line for line in new_lines if not line.startswith("CONECT")]
        return '\n'.join(final_lines)

    # 清理 CONECT 记录中的氢原子引用
    if h_atom_serials:
        final_lines = []
        for line in new_lines:
            if line.startswith("CONECT"):
                parts = line.split()
                if len(parts) > 1:
                    try:
                        main_atom = int(parts[1])
                        if main_atom in h_atom_serials:
                            continue  # 跳过以氢原子为主的 CONECT
                        new_parts = ["CONECT", str(main_atom)]
                        for p in parts[2:]:
                            try:
                                if int(p) not in h_atom_serials:
                                    new_parts.append(p)
                            except (ValueError, TypeError):  # int() 转换可能失败
                                pass
                        if len(new_parts) > 2:
                            final_lines.append(" ".join(new_parts))
                    except (ValueError, IndexError):  # CONECT 记录解析可能失败
                        final_lines.append(line)
                else:
                    final_lines.append(line)
            else:
                final_lines.append(line)
        return '\n'.join(final_lines)
    
    return '\n'.join(new_lines)


def try_load_mol_from_smiles(pdb_file, ligand_resname):
    """
    尝试从 PDB 文件中提取配体的 SMILES 并用 RDKit 加载
    这是处理复杂配体（如多肽）的备用方案
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        # 尝试使用 Open Babel 转换 (如果可用)
        import subprocess
        import tempfile
        import shutil
        
        # 1) 在当前进程环境中定位 obabel
        #    优先使用环境变量 OBABEL_BINARY，其次用 shutil.which('obabel')
        obabel_bin = os.environ.get("OBABEL_BINARY") or shutil.which("obabel")
        if not obabel_bin:
            print("[2D Diagram] SMILES fallback skipped: 'obabel' binary not found in PyMOL PATH.")
            print("             如果已安装，请在 PyMOL 启动环境中设置 OBABEL_BINARY=/full/path/to/obabel")
            return None
        
        # 创建临时文件
        with tempfile.NamedTemporaryFile(suffix='.smi', delete=False) as tmp:
            tmp_smi = tmp.name
        
        try:
            # 使用 obabel 将 PDB 转为 SMILES
            cmd_args = [obabel_bin, pdb_file, "-O", tmp_smi, "-osmi"]
            print(f"[2D Diagram] Trying Open Babel: {' '.join(cmd_args)}")
            result = subprocess.run(
                cmd_args,
                capture_output=True,
                text=True,
                timeout=20,
            )
            if result.returncode != 0:
                print(f"[2D Diagram] Open Babel failed (code={result.returncode}): {result.stderr.strip()}")
            elif os.path.exists(tmp_smi):
                with open(tmp_smi, 'r') as f:
                    smiles_line = f.readline().strip()
                    if smiles_line:
                        smiles = smiles_line.split()[0]
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            print(f"[2D Diagram] Loaded via SMILES (Open Babel): {smiles[:80]}...")
                            # 为后续 2D 布局生成一个合理的 2D 构象
                            try:
                                AllChem.Compute2DCoords(mol)
                            except Exception:
                                pass
                            return mol
        except subprocess.TimeoutExpired:
            print("[2D Diagram] Open Babel timeout while generating SMILES.")
        except FileNotFoundError:
            # 理论上不会走到这里，因为前面已经用 which 检查过
            print("[2D Diagram] Open Babel executable not found at runtime.")
        finally:
            if os.path.exists(tmp_smi):
                os.remove(tmp_smi)
    except Exception as e:
        print(f"[2D Diagram] SMILES fallback failed: {e}")

    return None


def _parse_element_from_atom_name(atom_name):
    """
    从 PDB 原子名解析元素类型。
    
    PDB 原子名规则：
    - 常见单字母元素 + 数字：N1 → N, O2 → O, C3 → C, S1 → S
    - 双字母元素 + 数字：CL1 → Cl, BR1 → Br, FE1 → Fe, ZN → Zn
    - 纯元素名：N → N, O → O
    
    Returns:
        元素符号（首字母大写），如果解析失败返回 None
    """
    if not atom_name:
        return None
    
    name = atom_name.strip()
    
    # 已知的双字母元素（PDB 中常见）
    two_letter_elements = {
        'CL', 'BR', 'FE', 'ZN', 'MG', 'MN', 'CO', 'CU', 'NI', 
        'SE', 'SI', 'NA', 'CA',
    }
    
    # 先尝试匹配双字母元素（大写形式）
    prefix2 = name[:2].upper()
    if prefix2 in two_letter_elements:
        return prefix2[0] + prefix2[1].lower()  # CL → Cl
    
    # 单字母元素 + 可选数字/后缀
    first_char = name[0].upper()
    if first_char in ('C', 'N', 'O', 'S', 'P', 'F', 'H', 'I', 'B', 'K'):
        # 确认第二个字符不是小写字母（否则可能是双字母元素）
        if len(name) == 1 or not name[1].isalpha() or name[1].isupper():
            return first_char
    
    # 最后尝试：用正则提取前导字母部分
    match = re.match(r'^([A-Za-z]{1,2})', name)
    if match:
        letters = match.group(1)
        if len(letters) == 2 and letters.upper() in two_letter_elements:
            return letters[0].upper() + letters[1].lower()
        return letters[0].upper()
    
    return None


def _transfer_pdb_info(mol_src, mol_dst):
    """
    尝试将原始分子的 PDB 原子名和坐标信息传递到 SMILES 重建的分子。
    
    策略：
    1. 如果原子数一致，按元素类型顺序逐一映射 PDBResidueInfo
    2. 如果原子数不匹配，至少在分子上挂载一个 _pdb_atom_map 属性
       供后续元素回退匹配使用
    """
    if mol_src is None or mol_dst is None:
        return
    
    # 收集原始分子的 PDB 原子名 → 元素映射
    pdb_atom_map = {}  # {atom_name: element_symbol}
    src_pdb_atoms = []  # [(idx, atom_name, element, x, y, z)]
    
    src_has_conf = mol_src.GetNumConformers() > 0
    src_conf = mol_src.GetConformer() if src_has_conf else None
    
    for atom in mol_src.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        elem = atom.GetSymbol()
        aname = ""
        if pdb_info:
            aname = pdb_info.GetName().strip()
            pdb_atom_map[aname] = elem
        
        pos = (0, 0, 0)
        if src_conf:
            p = src_conf.GetAtomPosition(atom.GetIdx())
            pos = (p.x, p.y, p.z)
        
        src_pdb_atoms.append((atom.GetIdx(), aname, elem, pos[0], pos[1], pos[2]))
    
    # 在目标分子上挂载映射信息（用于第三层回退匹配）
    mol_dst._pdb_atom_map = pdb_atom_map
    mol_dst._src_pdb_atoms = src_pdb_atoms
    
    # 尝试按元素对应关系复制 PDBResidueInfo
    n_src = mol_src.GetNumAtoms()
    n_dst = mol_dst.GetNumAtoms()
    
    if n_src == n_dst:
        # 原子数一致 → 按元素匹配尝试逐一复制
        # 按元素分组建立映射
        src_by_elem = {}
        for atom in mol_src.GetAtoms():
            elem = atom.GetSymbol()
            src_by_elem.setdefault(elem, []).append(atom.GetIdx())
        
        dst_by_elem = {}
        for atom in mol_dst.GetAtoms():
            elem = atom.GetSymbol()
            dst_by_elem.setdefault(elem, []).append(atom.GetIdx())
        
        # 检查每个元素的数量是否一致
        elem_match = all(
            len(src_by_elem.get(e, [])) == len(dst_by_elem.get(e, []))
            for e in set(list(src_by_elem.keys()) + list(dst_by_elem.keys()))
        )
        
        if elem_match:
            # 按元素分组、按顺序对应复制 PDB 信息
            copied = 0
            for elem in src_by_elem:
                src_indices = src_by_elem[elem]
                dst_indices = dst_by_elem.get(elem, [])
                for si, di in zip(src_indices, dst_indices):
                    src_atom = mol_src.GetAtomWithIdx(si)
                    dst_atom = mol_dst.GetAtomWithIdx(di)
                    pdb_info = src_atom.GetPDBResidueInfo()
                    if pdb_info:
                        # 深拷贝 PDB 信息到目标原子
                        new_info = Chem.AtomPDBResidueInfo()
                        new_info.SetName(pdb_info.GetName())
                        new_info.SetResidueName(pdb_info.GetResidueName())
                        new_info.SetResidueNumber(pdb_info.GetResidueNumber())
                        new_info.SetChainId(pdb_info.GetChainId())
                        new_info.SetIsHeteroAtom(pdb_info.GetIsHeteroAtom())
                        dst_atom.SetMonomerInfo(new_info)
                        copied += 1
            
            if copied > 0:
                print(f"[2D Diagram] ✅ 成功复制 {copied} 个原子的 PDB 信息到 SMILES 分子")
        else:
            print(f"[2D Diagram] ⚠️ 原子元素组成不匹配，无法复制 PDB 信息（src={n_src}, dst={n_dst}）")
    else:
        print(f"[2D Diagram] ⚠️ 原子数不匹配（src={n_src}, dst={n_dst}），跳过 PDB 信息复制，保留 _pdb_atom_map 供回退匹配")


def determine_bond_orders_from_3d(mol, temp_pdb_path=None, ligand_resname=None):
    """
    使用多种策略从 3D 坐标推断键级（单键/双键/芳香键）

    策略优先级：
    1. Open Babel SMILES 转换（最可靠）
    2. rdDetermineBonds 推断（RDKit 2022.09+）
    3. 返回原始分子（回退）

    Args:
        mol: RDKit 分子对象（必须有 3D 构象）
        temp_pdb_path: 配体 PDB 文件路径（用于 Open Babel 转换）
        ligand_resname: 配体残基名称（用于日志）

    Returns:
        修改后的分子对象（如果推断失败则返回原始分子）
    """
    if mol is None:
        print("[2D Diagram] ⚠️ mol 为 None，跳过键级推断")
        return mol

    # ============ 策略 1: Open Babel SMILES 转换 ============
    # 这是最可靠的方法，因为 Open Babel 可以正确解析 PDB 并生成包含键级的 SMILES
    if temp_pdb_path and os.path.exists(temp_pdb_path):
        print(f"[2D Diagram] 🔧 尝试使用 Open Babel 推断键级...")
        print(f"[2D Diagram] PDB 文件: {temp_pdb_path}")

        import subprocess
        import tempfile
        import shutil

        obabel_bin = os.environ.get("OBABEL_BINARY") or shutil.which("obabel")
        print(f"[2D Diagram] obabel 路径: {obabel_bin}")

        if obabel_bin:
            tmp_smi = None
            try:
                with tempfile.NamedTemporaryFile(suffix='.smi', delete=False, mode='w') as tmp:
                    tmp_smi = tmp.name

                cmd_args = [obabel_bin, temp_pdb_path, "-O", tmp_smi, "-osmi"]
                print(f"[2D Diagram] 执行命令: {' '.join(cmd_args)}")

                # 使用 Popen 以获得更好的控制，避免卡住
                proc = subprocess.Popen(
                    cmd_args,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    # 关键：设置 stdin 为 DEVNULL，防止 obabel 等待输入
                    stdin=subprocess.DEVNULL
                )

                try:
                    stdout, stderr = proc.communicate(timeout=30)
                    print(f"[2D Diagram] obabel 返回码: {proc.returncode}")
                    if stderr:
                        print(f"[2D Diagram] obabel stderr: {stderr[:200]}")
                    if stdout:
                        print(f"[2D Diagram] obabel stdout: {stdout[:200]}")
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.communicate()
                    print("[2D Diagram] Open Babel 超时 (30s)")
                    if tmp_smi and os.path.exists(tmp_smi):
                        os.remove(tmp_smi)
                    # 继续尝试其他策略
                else:
                    if proc.returncode == 0 and tmp_smi and os.path.exists(tmp_smi):
                        with open(tmp_smi, 'r') as f:
                            smiles_line = f.readline().strip()

                        print(f"[2D Diagram] SMILES: {smiles_line[:100] if smiles_line else 'empty'}")

                        if smiles_line:
                            smiles = smiles_line.split()[0]
                            mol_from_smiles = Chem.MolFromSmiles(smiles)

                            if mol_from_smiles:
                                # 使用 SMILES 分子的键级，但保留原始分子的 3D 坐标
                                # 通过 AssignBondOrdersFromTemplate 实现
                                try:
                                    from rdkit.Chem import AllChem
                                    # 将原始分子作为 3D 模板
                                    mol_with_orders = AllChem.AssignBondOrdersFromTemplate(mol_from_smiles, mol)

                                    # 统计键类型
                                    bond_types = {}
                                    for bond in mol_with_orders.GetBonds():
                                        bt = str(bond.GetBondType())
                                        bond_types[bt] = bond_types.get(bt, 0) + 1

                                    print(f"[2D Diagram] ✅ Open Babel + AssignBondOrdersFromTemplate 成功: {bond_types}")

                                    if tmp_smi and os.path.exists(tmp_smi):
                                        os.remove(tmp_smi)

                                    return mol_with_orders
                                except Exception as e:
                                    print(f"[2D Diagram] AssignBondOrdersFromTemplate 失败: {e}")
                                    # 回退：直接使用 SMILES 分子（会丢失 3D 坐标，但键级正确）
                                    # 🔧 修复：尝试将原始 mol 的 PDB 原子名和坐标映射到新分子
                                    try:
                                        AllChem.Compute2DCoords(mol_from_smiles)

                                        # 尝试从原始分子复制 PDB 原子名到 SMILES 分子
                                        # 这让后续的名称匹配策略仍然可用
                                        _transfer_pdb_info(mol, mol_from_smiles)

                                        bond_types = {}
                                        for bond in mol_from_smiles.GetBonds():
                                            bt = str(bond.GetBondType())
                                            bond_types[bt] = bond_types.get(bt, 0) + 1

                                        print(f"[2D Diagram] ✅ 使用 Open Babel SMILES 分子（键级正确）: {bond_types}")

                                        if tmp_smi and os.path.exists(tmp_smi):
                                            os.remove(tmp_smi)

                                        return mol_from_smiles
                                    except Exception:  # RDKit/文件操作可能失败
                                        pass

                    if tmp_smi and os.path.exists(tmp_smi):
                        os.remove(tmp_smi)

            except Exception as e:
                print(f"[2D Diagram] Open Babel 失败: {e}")
                import traceback
                traceback.print_exc()
                if tmp_smi and os.path.exists(tmp_smi):
                    try:
                        os.remove(tmp_smi)
                    except OSError:  # 临时文件删除可能失败
                        pass
        else:
            print("[2D Diagram] Open Babel 未找到，跳过此策略")

    # ============ 策略 2: rdDetermineBonds ============
    if HAS_DETERMINE_BONDS and mol.GetNumConformers() > 0:
        print(f"[2D Diagram] 🔧 尝试使用 rdDetermineBonds 推断键级...")

        charges_to_try = [0, -1, -2, 1, 2]

        for charge in charges_to_try:
            try:
                rw_mol = Chem.RWMol(mol)
                rdDetermineBonds.DetermineBonds(rw_mol, charge=charge)
                result_mol = rw_mol.GetMol()

                # 检查结果是否合理（不应该有三键，除非配体确实有炔基）
                bond_types = {}
                for bond in result_mol.GetBonds():
                    bt = str(bond.GetBondType())
                    bond_types[bt] = bond_types.get(bt, 0) + 1

                # 如果有三键，可能是推断错误，跳过这个电荷值
                if bond_types.get('TRIPLE', 0) > 0:
                    print(f"[2D Diagram] ⚠️ charge={charge} 产生了 {bond_types.get('TRIPLE', 0)} 个三键，可能不正确")
                    continue

                print(f"[2D Diagram] ✅ rdDetermineBonds 成功 (charge={charge}): {bond_types}")
                return result_mol

            except Exception as e:
                if "does not match input" in str(e) or "valence" in str(e).lower():
                    continue
                print(f"[2D Diagram] rdDetermineBonds 失败 (charge={charge}): {e}")

    # ============ 策略 3: 回退 ============
    print("[2D Diagram] ⚠️ 所有键级推断策略失败，使用原始分子")
    return mol


def rebuild_bonds_by_distance(mol):
    """
    基于原子间距离重建分子键连接
    用于修复 RDKit 自动推断键连接失败的情况
    
    使用标准共价键半径来判断原子间是否应该成键
    """
    from rdkit import Chem
    
    # 标准共价键半径 (Å)
    COVALENT_RADII = {
        6: 0.77,   # C
        7: 0.75,   # N
        8: 0.73,   # O
        9: 0.71,   # F
        15: 1.07,  # P
        16: 1.05,  # S
        17: 0.99,  # Cl
        35: 1.14,  # Br
        53: 1.33,  # I
        1: 0.31,   # H
    }
    DEFAULT_RADIUS = 0.77
    
    if mol.GetNumConformers() == 0:
        print("[2D Diagram] No conformer available for distance-based bond rebuild")
        return mol
    
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()
    
    # 创建新的可编辑分子
    emol = Chem.RWMol(Chem.Mol())
    
    # 复制原子
    atom_map = {}
    for i in range(num_atoms):
        atom = mol.GetAtomWithIdx(i)
        new_idx = emol.AddAtom(Chem.Atom(atom.GetAtomicNum()))
        atom_map[i] = new_idx
    
    # 基于距离添加键
    # 为了避免多肽折叠时不同片段之间“跨链连线成网”，我们增加一个
    # 以 PDB 原子序号为基础的局部窗口约束：只允许在序号相差较小的
    # 原子之间尝试成键（典型情况下这些原子在同一残基或相邻残基）。
    tolerance = 0.4  # 容差因子
    
    # 预先提取 PDB 原子序号（如果存在），否则退回到 RDKit 索引
    pdb_serials = []
    for i in range(num_atoms):
        atom = mol.GetAtomWithIdx(i)
        info = atom.GetPDBResidueInfo()
        if info is not None:
            try:
                pdb_serials.append(int(info.GetSerialNumber()))
            except Exception:
                pdb_serials.append(i)
        else:
            pdb_serials.append(i)
    
    for i in range(num_atoms):
        pos_i = conf.GetAtomPosition(i)
        r_i = COVALENT_RADII.get(mol.GetAtomWithIdx(i).GetAtomicNum(), DEFAULT_RADIUS)
        
        for j in range(i + 1, num_atoms):
            # 关键约束：只在原子序号相差较小的局部范围内尝试成键，
            # 避免由于多肽折叠导致的远程近距离原子误连。
            if abs(pdb_serials[i] - pdb_serials[j]) > 6:
                continue
            
            pos_j = conf.GetAtomPosition(j)
            r_j = COVALENT_RADII.get(mol.GetAtomWithIdx(j).GetAtomicNum(), DEFAULT_RADIUS)
            
            # 计算距离
            dist = ((pos_i.x - pos_j.x)**2 + (pos_i.y - pos_j.y)**2 + (pos_i.z - pos_j.z)**2)**0.5
            
            # 判断是否应该成键
            max_bond_dist = (r_i + r_j) * (1 + tolerance)
            if dist <= max_bond_dist:
                try:
                    emol.AddBond(atom_map[i], atom_map[j], Chem.BondType.SINGLE)
                except Exception:  # RDKit 添加键可能失败
                    pass
    
    # 添加构象
    new_conf = Chem.Conformer(emol.GetNumAtoms())
    for i in range(num_atoms):
        pos = conf.GetAtomPosition(i)
        new_conf.SetAtomPosition(atom_map[i], pos)
    emol.AddConformer(new_conf, assignId=True)
    
    result_mol = emol.GetMol()
    
    # 尝试推断键级
    try:
        Chem.SanitizeMol(result_mol, Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                        Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                        Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
                        Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION |
                        Chem.SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
    except Exception:  # RDKit SanitizeMol 可能失败
        pass
    
    print(f"[2D Diagram] Rebuilt molecule: {result_mol.GetNumAtoms()} atoms, {result_mol.GetNumBonds()} bonds")
    return result_mol


def check_bond_sanity(mol):
    """
    检查分子键连接是否合理
    返回 True 如果键连接看起来正常，False 如果可能有问题
    """
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    
    if num_atoms == 0:
        return False
    
    # 正常有机分子的键数应该接近原子数
    # 典型范围: bonds ≈ atoms * 0.8 ~ 1.3 (对于线性分子接近1.0，有环的分子稍高)
    # 大分子（多肽等）如果出现 bonds/atoms 明显 > 1.2 往往就是“蜘蛛网”式错误连接
    bond_ratio = num_bonds / num_atoms if num_atoms > 0 else 0
    
    print(f"[2D Diagram] Bond sanity check: {num_bonds} bonds / {num_atoms} atoms = {bond_ratio:.2f}")
    
    # 分级阈值：对大分子更严格
    if num_atoms >= 150:
        ratio_threshold = 1.20
    elif num_atoms >= 80:
        ratio_threshold = 1.30
    else:
        ratio_threshold = 1.40
    
    if bond_ratio > ratio_threshold:
        print(f"[2D Diagram] ⚠️ Abnormal bond ratio: {bond_ratio:.2f} (bonds={num_bonds}, atoms={num_atoms}, threshold={ratio_threshold:.2f})")
        return False
    
    # 检查是否有原子连接了太多键 (正常最多 4 个，特殊情况如 S, P 可能有 5-6 个)
    max_degree = 0
    high_degree_count = 0
    for atom in mol.GetAtoms():
        degree = atom.GetDegree()
        if degree > max_degree:
            max_degree = degree
        if degree > 4:
            high_degree_count += 1
    
    # 如果有超过 10% 的原子连接度 > 4，说明有问题
    if high_degree_count > num_atoms * 0.1:
        print(f"[2D Diagram] ⚠️ Too many high-degree atoms: {high_degree_count}/{num_atoms}")
        return False
    
    if max_degree > 6:
        print(f"[2D Diagram] ⚠️ Abnormal max atom degree: {max_degree}")
        return False
    
    return True

def generate_2d_interaction_diagram(csv_path, ligand_resname, pdb_file=None, obj_name=None,
                                     output_path=None, width=1600, height=1200, dpi=150,
                                     min_confidence=0.8):
    print(f"[2D Diagram] ========================================")
    print(f"[2D Diagram] CSV Path: {csv_path}")
    print(f"[2D Diagram] Ligand: {ligand_resname}")
    print(f"[2D Diagram] ========================================")

    # 局部工具函数：从 RDKit 分子中提取 3D 坐标和环中心
    # 必须在函数开始时定义，以便在所有代码路径中都可用
    def extract_3d_info(mol):
        coords = {}
        rings = []
        if mol and mol.GetNumConformers() > 0:
            conf = mol.GetConformer()
            for atom in mol.GetAtoms():
                idx = atom.GetIdx()
                pos = conf.GetAtomPosition(idx)
                coords[idx] = (pos.x, pos.y, pos.z)
            
            try:
                for ri, ring in enumerate(mol.GetRingInfo().AtomRings()):
                    pts = [conf.GetAtomPosition(aid) for aid in ring]
                    cen = (
                        sum(p.x for p in pts) / len(pts),
                        sum(p.y for p in pts) / len(pts),
                        sum(p.z for p in pts) / len(pts),
                    )
                    rings.append({"id": ri, "coord": cen, "indices": list(ring)})
            except Exception:
                pass
        return coords, rings

    # 1. 提取结构并清理氢原子
    mol_draw = None
    saved_temp_pdb = None  # 初始化变量，避免后续引用错误
    
    if obj_name:
        try:
            from pymol import cmd
            temp_pdb = tempfile.mktemp(suffix=".pdb")
            space = {'c': [], 'r': []}
            cmd.iterate(f"first ({obj_name} and resn {ligand_resname})", "c.append(chain); r.append(resi)", space=space)
            sel = f"{obj_name} and resn {ligand_resname}"
            if space['c']: sel += f" and chain {space['c'][0]} and resi {space['r'][0]}"
            print(f"[2D Diagram] Extracting: {sel}")
            cmd.save(temp_pdb, sel)
            


            # 保存 temp_pdb 路径供后续使用
            saved_temp_pdb = temp_pdb
            
            # 策略调整：优先使用标准加载 (自动去氢 + 圣化)
            try:
                mol_draw = Chem.MolFromPDBFile(temp_pdb, removeHs=True, sanitize=True)
                if mol_draw:
                    print(f"[2D Diagram] Standard load: {mol_draw.GetNumAtoms()} atoms, {mol_draw.GetNumBonds()} bonds")

                    # 🔧 关键修复：推断正确的键级（双键、芳香键等）
                    # PDB 格式不存储键级，需要通过 Open Babel 或 rdDetermineBonds 推断
                    mol_draw = determine_bond_orders_from_3d(mol_draw, temp_pdb, ligand_resname)

            except Exception as e:
                print(f"[2D Diagram] Standard load exception: {e}")
                mol_draw = None
                
            # 如果标准加载失败，使用终极文本清洗加载
            if mol_draw is None:
                print(f"[2D Diagram] Standard load failed, using Text-Based Cleaning...")
                try:
                    with open(temp_pdb, 'r') as f:
                        raw_pdb = f.read()
                    
                    
                    # 关键修复：不保留 CONECT 记录，让 RDKit 基于距离推断键
                    # 因为 PyMOL 导出的 CONECT 记录可能不完整或有问题
                    clean_pdb = clean_pdb_block(raw_pdb)
                    
                    # 尝试方案1：使用清洗后的 PDB（保留 CONECT）
                    mol_draw = Chem.MolFromPDBBlock(clean_pdb, removeHs=True, sanitize=False)
                    
                    if mol_draw:
                        print(f"[2D Diagram] Text cleaning load (with CONECT): {mol_draw.GetNumAtoms()} atoms, {mol_draw.GetNumBonds()} bonds")
                        # 重新计算拓扑
                        mol_draw.UpdatePropertyCache(strict=False)
                        # 尝试圣化以获得键级 (如果可能)
                        try:
                            Chem.SanitizeMol(mol_draw, Chem.SanitizeFlags.SANITIZE_FINDRADICALS|
                                           Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|
                                           Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|
                                           Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|
                                           Chem.SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
                        except Exception: pass  # RDKit SanitizeMol 可能失败

                        # 🔧 推断键级
                        mol_draw = determine_bond_orders_from_3d(mol_draw, temp_pdb, ligand_resname)
                    else:
                        print(f"[2D Diagram] Text cleaning load failed, mol_draw is None")
                    
                    # 蜘蛛网检测与修复：如果加载出的分子键太疯狂，优先尝试“丢弃 CONECT 重新加载”
                    # 真正的“无键 3D 投影”只在后面的统一修复逻辑里触发，避免过早丢失所有键信息。
                    if mol_draw and not check_bond_sanity(mol_draw):
                        print(f"[2D Diagram] ⚠️ Detected corrupted CONECT records (Spiderweb). Retrying without CONECT...")
                        clean_pdb_no_conect = clean_pdb_block(raw_pdb, ignore_connect=True)
                        mol_retry = Chem.MolFromPDBBlock(clean_pdb_no_conect, removeHs=True, sanitize=False)
                        if mol_retry:
                            mol_draw = mol_retry
                            print(f"[2D Diagram] ✅ Spiderweb fix applied: Used structure without CONECT records.")
                            mol_draw.UpdatePropertyCache(strict=False)
                            try:
                                Chem.SanitizeMol(mol_draw, catchErrors=True)
                            except Exception:
                                pass
                            # 再次检查仅用于日志，是否依然异常由后续统一逻辑处理
                            check_bond_sanity(mol_draw)
                        else:
                            print(f"[2D Diagram] Retry without CONECT failed; keeping original topology for later fixes.")

                except Exception as e:
                    print(f"[2D Diagram] Text cleaning failed: {e}")
                    import traceback
                    traceback.print_exc()

            # 不要删除 temp_pdb，后面可能需要用于 SMILES 转换
            # if os.path.exists(temp_pdb): os.remove(temp_pdb)
        except Exception as e:
            print(f"[2D Diagram] PyMOL Error: {e}")
    
    elif pdb_file and os.path.exists(pdb_file):
        print(f"[2D Diagram] Loading from PDB file: {pdb_file}")
        try:
            # 首先从 PDB 文件中提取配体部分
            with open(pdb_file, 'r') as f:
                raw_pdb = f.read()
            
            # 提取配体行 (HETATM 或 ATOM 中匹配 ligand_resname 的行)
            ligand_lines = []
            conect_lines = []
            ligand_serials = set()
            
            for line in raw_pdb.split('\n'):
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    # 残基名称在 17-20 列 (0-indexed: 17,18,19)
                    resname = line[17:20].strip()
                    if resname == ligand_resname:
                        ligand_lines.append(line)
                        try:
                            serial = int(line[6:11].strip())
                            ligand_serials.add(serial)
                        except (ValueError, TypeError):  # int() 转换可能失败
                            pass
                elif line.startswith("CONECT"):
                    conect_lines.append(line)
            
            if not ligand_lines:
                print(f"[2D Diagram] No ligand '{ligand_resname}' found in PDB file")
                mol_draw = None
            else:
                # 过滤 CONECT 记录，只保留配体原子之间的连接
                filtered_conect = []
                for line in conect_lines:
                    parts = line.split()
                    if len(parts) > 1:
                        try:
                            main_atom = int(parts[1])
                            if main_atom in ligand_serials:
                                new_parts = ["CONECT", str(main_atom)]
                                for p in parts[2:]:
                                    try:
                                        if int(p) in ligand_serials:
                                            new_parts.append(p)
                                    except (ValueError, TypeError):  # int() 转换可能失败
                                        pass
                                if len(new_parts) > 2:
                                    filtered_conect.append(" ".join(new_parts))
                        except (ValueError, IndexError):  # CONECT 记录解析可能失败
                            pass
                
                # 构建配体 PDB 块
                ligand_pdb = '\n'.join(ligand_lines + filtered_conect + ['END'])
                print(f"[2D Diagram] Extracted ligand: {len(ligand_lines)} atoms, {len(filtered_conect)} CONECT records")
                
                # 保存临时文件供后续 SMILES 转换使用
                temp_ligand_pdb = tempfile.mktemp(suffix=".pdb")
                with open(temp_ligand_pdb, 'w') as f:
                    f.write(ligand_pdb)
                saved_temp_pdb = temp_ligand_pdb
                
                # 清理氢原子
                clean_pdb = clean_pdb_block(ligand_pdb)
                
                # 尝试加载
                mol_draw = Chem.MolFromPDBBlock(clean_pdb, removeHs=True, sanitize=True)
                if not mol_draw:
                    # 备用方案：不圣化
                    mol_draw = Chem.MolFromPDBBlock(clean_pdb, removeHs=True, sanitize=False)
                    if mol_draw:
                        try:
                            Chem.SanitizeMol(mol_draw, catchErrors=True)
                        except Exception: pass  # RDKit SanitizeMol 可能失败

                if mol_draw:
                    print(f"[2D Diagram] Loaded ligand: {mol_draw.GetNumAtoms()} atoms, {mol_draw.GetNumBonds()} bonds")
                    # 🔧 推断键级
                    mol_draw = determine_bond_orders_from_3d(mol_draw, saved_temp_pdb, ligand_resname)
        except Exception as e:
            print(f"[2D Diagram] PDB Load Error: {e}")
            import traceback
            traceback.print_exc()

    if not mol_draw:
        print(f"[2D Diagram] Failed to load molecule.")
        return None

    # 检测大型多肽配体并给出警告
    num_atoms = mol_draw.GetNumAtoms()
    is_large_peptide = num_atoms > 50
    if is_large_peptide:
        print(f"[2D Diagram] ⚠️ 检测到大型多肽配体 ({num_atoms} 个原子)")
        print(f"[2D Diagram] ⚠️ 对于大型多肽，2D 相互作用图可能布局不佳")
        print(f"[2D Diagram] ⚠️ 建议使用 3D 可视化查看相互作用")

    # 再次确认清理 (三层保险：处理可能的残留)
    try:
        # 强制移除所有氢原子 - 使用更彻底的方法
        mol_draw = Chem.RemoveAllHs(mol_draw)  # 🔧 使用 RemoveAllHs 替代 RemoveHs，更彻底

        # 再次检查原子序数，手动移除任何残留的氢原子
        mw = Chem.RWMol(mol_draw)
        atoms_to_remove = [i for i in range(mw.GetNumAtoms()) if mw.GetAtomWithIdx(i).GetAtomicNum() <= 1]
        if atoms_to_remove:
            print(f"[2D Diagram] 🔧 手动移除 {len(atoms_to_remove)} 个残留氢原子")
            atoms_to_remove.sort(reverse=True)
            for i in atoms_to_remove:
                mw.RemoveAtom(i)
        mol_draw = mw.GetMol()

        # 🔧 再次尝试 RemoveAllHs，确保完全清除
        try:
            mol_draw = Chem.RemoveAllHs(mol_draw)
        except Exception:  # RDKit RemoveAllHs 可能失败
            pass

    except Exception as e:
        print(f"[2D Diagram] H Removal Error: {e}")
    
    # 2. 提前提取 3D 坐标 (关键！必须在投影到 2D 之前完成)
    atom_coords_3d = {}
    ring_centroids_3d = []

    # 初始提取
    atom_coords_3d, ring_centroids_3d = extract_3d_info(mol_draw)
    
    # 检查键连接是否合理，如果有问题则使用 Open Babel SMILES 重建
    # 但对于大型多肽配体（>50 原子），跳过 SMILES 重建，因为会丢失 3D 坐标导致布局变差
    bond_sanity_ok = check_bond_sanity(mol_draw)
    
    if not bond_sanity_ok and not is_large_peptide:
        print(f"[2D Diagram] Attempting to fix bond connectivity using Open Babel...")
        
        # 使用 Open Babel SMILES 重建分子 (主要方案)
        pdb_for_smiles = (
            saved_temp_pdb
            if 'saved_temp_pdb' in dir() and saved_temp_pdb and os.path.exists(saved_temp_pdb)
            else pdb_file
        )
        mol_from_smiles = try_load_mol_from_smiles(pdb_for_smiles, ligand_resname)
        if mol_from_smiles and check_bond_sanity(mol_from_smiles):
            mol_draw = mol_from_smiles
            # SMILES 重建的分子已经有 2D 坐标，更新 3D 坐标信息
            atom_coords_3d, ring_centroids_3d = extract_3d_info(mol_draw)
            print(f"[2D Diagram] ✅ Fixed using Open Babel SMILES reconstruction")
        else:
            # 如果 SMILES 方法失败，尝试基于距离重建键
            if not check_bond_sanity(mol_draw):
                num_atoms_for_rebuild = mol_draw.GetNumAtoms()
                print(f"[2D Diagram] Trying distance-based bond inference on {num_atoms_for_rebuild} atoms...")
                mol_rebuilt = rebuild_bonds_by_distance(mol_draw)
                if check_bond_sanity(mol_rebuilt):
                    mol_draw = mol_rebuilt
                    atom_coords_3d, ring_centroids_3d = extract_3d_info(mol_draw)
                    print(f"[2D Diagram] ✅ Fixed using distance-based bond rebuild")
                else:
                    print(f"[2D Diagram] ⚠️ Bond connectivity issues persist, proceeding with current structure")
    elif not bond_sanity_ok and is_large_peptide:
        print(f"[2D Diagram] ⚠️ 跳过 SMILES 重建（大型多肽配体），保留原始 3D 坐标以获得更好的 2D 布局")
    
    # 清理临时文件
    if 'saved_temp_pdb' in dir() and saved_temp_pdb and os.path.exists(saved_temp_pdb):
        os.remove(saved_temp_pdb)
    
    # 3. 解析相互作用并匹配 (与 3D 视图一致的置信度过滤)
    interactions = []
    interaction_types_found = set()  # 记录找到的相互作用类型，用于动态图例
    # 🔧 统计未匹配交互，用于最终报告
    unmatched_interactions = []  # [(lig_atom_name, itype)]
    total_csv_interactions = 0   # CSV 中满足置信度的总交互数
    try:
        with open(csv_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            for row in reader:
                itype = row.get("Interaction", "").strip()
                conf_str = row.get("Confidence", "1.0")
                try:
                    conf = float(conf_str)
                except (ValueError, TypeError):
                    conf = 1.0  # 没有置信度字段时默认视为 1.0
                if conf < min_confidence:
                    continue
                
                total_csv_interactions += 1
                lig_atom_name = row.get("Ligand_Atom", "").strip()
                prot_res = row.get("Protein_Residue", "").strip()
                lx, ly, lz = row.get("Ligand_Atom_X"), row.get("Ligand_Atom_Y"), row.get("Ligand_Atom_Z")
                
                target_idx, target_type = None, "atom"
                
                # 方案1：通过坐标匹配（如果 CSV 中有坐标列）
                if lx and ly and lz:
                    try:
                        lx, ly, lz = float(lx), float(ly), float(lz)
                    except (ValueError, TypeError):
                        lx, ly, lz = None, None, None
                    
                    if lx is not None:
                        # 环匹配 (π-π 堆积, π-阳离子)
                        if "Ring" in lig_atom_name or "Cation" in lig_atom_name:
                            best_rid, best_rdist = None, float('inf')
                            for r_data in ring_centroids_3d:
                                rc = r_data["coord"]
                                d = ((rc[0]-lx)**2 + (rc[1]-ly)**2 + (rc[2]-lz)**2)**0.5
                                if d < best_rdist: best_rdist, best_rid = d, r_data["id"]
                            if best_rid is not None and best_rdist < 0.6:
                                target_idx, target_type = best_rid, "ring"
                        
                        # 原子匹配 (如果环没匹配上或不是环)
                        # 放宽阈值：从 0.3 Å 改为 1.5 Å，因为 RDKit 加载后原子可能有轻微偏移
                        if target_idx is None:
                            best_aid, best_adist = None, float('inf')
                            for aid, ac in atom_coords_3d.items():
                                d = ((ac[0]-lx)**2 + (ac[1]-ly)**2 + (ac[2]-lz)**2)**0.5
                                if d < best_adist: best_adist, best_aid = d, aid
                            # 使用更宽松的阈值（1.5 Å）以确保匹配成功
                            if best_aid is not None and best_adist < 1.5:
                                target_idx, target_type = best_aid, "atom"
                                print(f"[2D DEBUG] Coordinate matched: {lig_atom_name} -> atom {target_idx} (dist={best_adist:.3f}Å)")

                # 方案2：通过原子名称匹配（如果 CSV 中没有坐标列）
                if target_idx is None and lig_atom_name:
                    # 特殊处理：Ring(...) 或 ring 表示环中心
                    if lig_atom_name.lower() == "ring" or lig_atom_name.lower().startswith("ring("):
                        # 对于 π-π 堆积，选择第一个环
                        if ring_centroids_3d:
                            target_idx = ring_centroids_3d[0]["id"]
                            target_type = "ring"
                            print(f"[2D DEBUG] Matched Ring to ring centroid {target_idx}")
                    # 特殊处理：Cation(...) 表示阳离子
                    elif lig_atom_name.lower().startswith("cation("):
                        # 对于 π-阳离子相互作用，尝试匹配环中心
                        if ring_centroids_3d:
                            target_idx = ring_centroids_3d[0]["id"]
                            target_type = "ring"
                            print(f"[2D DEBUG] Matched Cation to ring centroid {target_idx}")
                    else:
                        # 通过 PDB 原子名称匹配 - 改进版：找到所有匹配的原子，选择离 3D 坐标最近的
                        # 这样即使配体有多个同名原子（如多个 O），也能正确区分
                        matching_atoms = []
                        for aid in range(mol_draw.GetNumAtoms()):
                            atom = mol_draw.GetAtomWithIdx(aid)
                            pdb_info = atom.GetPDBResidueInfo()
                            if pdb_info:
                                atom_name = pdb_info.GetName().strip()
                                if atom_name == lig_atom_name:
                                    matching_atoms.append(aid)

                        if len(matching_atoms) == 1:
                            # 只有一个匹配，直接使用
                            target_idx = matching_atoms[0]
                            target_type = "atom"
                            print(f"[2D DEBUG] Name matched (unique): {lig_atom_name} -> atom {target_idx}")
                        elif len(matching_atoms) > 1:
                            # 多个同名原子，无法区分（没有坐标信息）
                            # 警告用户并使用第一个
                            print(f"[2D DEBUG] ⚠️ Multiple atoms named '{lig_atom_name}' found ({len(matching_atoms)} total)")
                            print(f"[2D DEBUG] ⚠️ Cannot distinguish without coordinates - using first match (atom {matching_atoms[0]})")
                            print(f"[2D DEBUG] ⚠️ For accurate mapping, re-run interaction analysis to generate CSV with coordinates")
                            target_idx = matching_atoms[0]
                            target_type = "atom"
                
                # 🔧 方案3：元素+名称模糊匹配（第三层回退）
                # 当 SMILES 重建导致分子丢失 PDB 信息和 3D 坐标时，
                # 通过解析 CSV 原子名中的元素类型，在 mol_draw 中查找相同元素的原子
                if target_idx is None and lig_atom_name:
                    # 跳过 Ring/Cation 等特殊名称（已在方案2中处理过）
                    if not (lig_atom_name.lower().startswith("ring") or 
                            lig_atom_name.lower().startswith("cation")):
                        parsed_elem = _parse_element_from_atom_name(lig_atom_name)
                        if parsed_elem:
                            # 在 mol_draw 中找所有相同元素的原子
                            same_elem_atoms = []
                            for aid in range(mol_draw.GetNumAtoms()):
                                atom = mol_draw.GetAtomWithIdx(aid)
                                if atom.GetSymbol().upper() == parsed_elem.upper():
                                    same_elem_atoms.append(aid)
                            
                            if len(same_elem_atoms) == 1:
                                # 唯一同元素原子 → 直接匹配
                                target_idx = same_elem_atoms[0]
                                target_type = "atom"
                                print(f"[2D DEBUG] 元素回退匹配: {lig_atom_name} -> atom {target_idx} (element={parsed_elem}, unique)")
                            elif len(same_elem_atoms) > 1:
                                # 多个同元素原子 → 尝试用 3D 坐标选最近的
                                matched_by_3d = False
                                if lx is not None and ly is not None and lz is not None and atom_coords_3d:
                                    # CSV 有坐标且分子有 3D 信息 → 选最近的同元素原子
                                    best_aid, best_adist = None, float('inf')
                                    for aid in same_elem_atoms:
                                        ac = atom_coords_3d.get(aid)
                                        if ac:
                                            d = ((ac[0]-lx)**2 + (ac[1]-ly)**2 + (ac[2]-lz)**2)**0.5
                                            if d < best_adist:
                                                best_adist, best_aid = d, aid
                                    if best_aid is not None:
                                        target_idx = best_aid
                                        target_type = "atom"
                                        matched_by_3d = True
                                        print(f"[2D DEBUG] 元素回退匹配: {lig_atom_name} -> atom {target_idx} (element={parsed_elem}, 3D nearest, dist={best_adist:.3f}Å)")
                                
                                # 尝试用 _pdb_atom_map（从原始分子传递的映射）辅助匹配
                                if not matched_by_3d and hasattr(mol_draw, '_src_pdb_atoms'):
                                    src_atoms = mol_draw._src_pdb_atoms
                                    # 找到原始分子中与 lig_atom_name 同名的原子索引
                                    for s_idx, s_name, s_elem, sx, sy, sz in src_atoms:
                                        if s_name == lig_atom_name and s_idx < len(same_elem_atoms):
                                            # 使用原始索引在同元素列表中找对应位置
                                            # 按元素在 mol_draw 中的出现顺序映射
                                            src_same_elem = [si for si, sn, se, *_ in src_atoms if se.upper() == parsed_elem.upper()]
                                            if s_idx in src_same_elem:
                                                pos_in_group = src_same_elem.index(s_idx)
                                                if pos_in_group < len(same_elem_atoms):
                                                    target_idx = same_elem_atoms[pos_in_group]
                                                    target_type = "atom"
                                                    matched_by_3d = True
                                                    print(f"[2D DEBUG] 元素回退匹配: {lig_atom_name} -> atom {target_idx} (element={parsed_elem}, PDB映射)")
                                            break
                                
                                if not matched_by_3d:
                                    # 无法进一步区分 → 使用第一个同元素原子
                                    target_idx = same_elem_atoms[0]
                                    target_type = "atom"
                                    print(f"[2D DEBUG] 元素回退匹配: {lig_atom_name} -> atom {target_idx} (element={parsed_elem}, first of {len(same_elem_atoms)})")
                            else:
                                print(f"[2D DEBUG] 元素回退匹配失败: {lig_atom_name} (element={parsed_elem}) - 分子中无此元素原子")
                
                # 🔧 未匹配交互的调试日志
                if target_idx is None:
                    print(f"[2D DEBUG] ⚠️ 未能匹配: {lig_atom_name} ({itype}) - 该相互作用将不会显示在2D图中")
                    unmatched_interactions.append((lig_atom_name, itype))

                if target_idx is not None:
                    # 标准化相互作用类型 (与 3D 视图一致)
                    normalized_itype = normalize_interaction_type(itype)

                    # 🔧 特殊处理：Pi-Pi 和 Pi-Cation 应该连接到环中心
                    # 如果当前匹配的是原子，检查该原子是否属于某个环
                    if normalized_itype in ['pipi', 'pication'] and target_type == "atom":
                        atom_idx = target_idx
                        # 检查该原子是否在已知环中
                        best_rid = None
                        for r_data in ring_centroids_3d:
                            if atom_idx in r_data["indices"]:
                                best_rid = r_data["id"]
                                break

                        if best_rid is not None:
                            target_idx = best_rid
                            target_type = "ring"
                            print(f"   -> Upgraded {lig_atom_name} to Ring {target_idx} Center (pipi/pication)")
                        else:
                            # 原子不在任何环中，但可能环检测失败
                            # 尝试找离这个原子最近的环中心
                            if ring_centroids_3d:
                                atom_coord = atom_coords_3d.get(atom_idx)
                                if atom_coord:
                                    best_rid, best_dist = None, float('inf')
                                    for r_data in ring_centroids_3d:
                                        rc = r_data["coord"]
                                        d = ((rc[0]-atom_coord[0])**2 + (rc[1]-atom_coord[1])**2 + (rc[2]-atom_coord[2])**2)**0.5
                                        if d < best_dist:
                                            best_dist, best_rid = d, r_data["id"]
                                    # 如果距离小于 3 Å，认为是同一个环
                                    if best_rid is not None and best_dist < 3.0:
                                        target_idx = best_rid
                                        target_type = "ring"
                                        print(f"   -> Upgraded {lig_atom_name} to nearest Ring {target_idx} (dist={best_dist:.2f}Å)")

                    interaction_types_found.add(normalized_itype)
                    interactions.append({
                        "res": prot_res,
                        "idx": target_idx,
                        "type": target_type,
                        "itype": itype,
                        "normalized_itype": normalized_itype,
                        "confidence": conf
                    })
                    print(f"[2D DEBUG] Mapped {lig_atom_name} to {target_type} {target_idx} ({normalized_itype})")
    except Exception as e:
        print(f"[2D Diagram] CSV Error: {e}")
        import traceback
        traceback.print_exc()

    # 🔧 输出未匹配交互的汇总统计
    if unmatched_interactions:
        # 按交互类型分组统计
        unmatched_by_type = {}
        for name, itype in unmatched_interactions:
            norm_type = normalize_interaction_type(itype)
            unmatched_by_type[norm_type] = unmatched_by_type.get(norm_type, 0) + 1
        type_summary = ", ".join(f"{count}个{t}" for t, count in unmatched_by_type.items())
        print(f"[2D Diagram] ⚠️ {len(unmatched_interactions)}/{total_csv_interactions} 个相互作用未能匹配到2D图（{type_summary}）")
    else:
        print(f"[2D Diagram] ✅ 全部 {total_csv_interactions} 个相互作用已成功匹配")

    # 4. 生成 2D 布局 (使用 RDKit CoordGen，带超时保护)
    import threading as _threading
    
    def _compute_2d_layout(mol, result_holder):
        """在超时保护下计算 2D 布局坐标。"""
        try:
            rdDepictor.SetPreferCoordGen(True)
            try:
                from rdkit.Chem import rdCoordGen
                rdCoordGen.AddCoords(mol)
                result_holder['ok'] = True
                print("[2D Diagram] Using rdCoordGen.AddCoords for layout")
            except Exception:
                rdDepictor.Compute2DCoords(mol)
                result_holder['ok'] = True
        except Exception as e:
            result_holder['error'] = str(e)
    
    layout_result = {'ok': False, 'error': None}
    layout_thread = _threading.Thread(target=_compute_2d_layout, args=(mol_draw, layout_result), daemon=True)
    layout_thread.start()
    layout_thread.join(timeout=30)  # 2D 布局最多 30 秒
    
    if layout_thread.is_alive() or not layout_result['ok']:
        if layout_thread.is_alive():
            print("[2D Diagram] ⚠️ CoordGen 2D 布局超时 (30s)，尝试 fallback...")
        elif layout_result['error']:
            print(f"[2D Diagram] Layout error (CoordGen failed): {layout_result['error']}")
        
        # Fallback：使用更简单的 AllChem.Compute2DCoords
        try:
            AllChem.Compute2DCoords(mol_draw)
            print("[2D Diagram] Fallback: AllChem.Compute2DCoords succeeded")
        except Exception as e2:
            print(f"[2D Diagram] Layout error (fallback failed): {e2}")
            return None
    
    # 动态画布
    conf_2d = mol_draw.GetConformer()
    pts = [conf_2d.GetAtomPosition(i) for i in range(mol_draw.GetNumAtoms())]
    phys_span = max([p.x for p in pts]) - min([p.x for p in pts]) if pts else 10.0
    dw = int(max(1200, phys_span * 45.0 + 800))
    dh = int(dw * 0.75)
    
    # 限制画布最大像素，防止内存爆炸
    # 对于大型多肽配体，使用更小的画布
    if is_large_peptide:
        max_width = 2000
        max_height = 1500
        print(f"[2D Diagram] 大型多肽配体：限制画布大小为 {max_width}x{max_height}")
    else:
        max_width = 3000
        max_height = 2250
    
    dw = min(dw, max_width)
    dh = min(dh, max_height)
    
    try:
        drawer = rdMolDraw2D.MolDraw2DCairo(dw, dh)
        opts = drawer.drawOptions()
        opts.padding = 0.20  # 适度的 padding
        opts.prepareMolsBeforeDrawing = True
        opts.bondLineWidth = 3.0
        opts.addStereoAnnotation = False  # 不显示 R/S 标签
        opts.addAtomIndices = False
        opts.includeAtomTags = False
        opts.explicitMethyl = False
        # opts.addHBonds = False # Not supported in some RDKit versions
        opts.minFontSize = 10
        # 🔧 关键修复：调整字体大小和原子标签显示
        opts.maxFontSize = 14  # 限制最大字体大小
        opts.annotationFontScale = 0.8  # 减小注释字体比例
        # 🔧 尝试设置 addChiralHs（某些 RDKit 版本可能不支持）
        try:
            opts.addChiralHs = False  # 不显示手性中心的氢原子
        except AttributeError:
            pass  # 如果属性不存在，忽略
        # 不使用 fixedBondLength，让 RDKit 自动计算更好的布局

        # 🔧 不显示原子电荷和氢原子计数
        opts.noAtomLabels = False  # 仍然显示原子符号

        # 🔧 关键修复：直接修改 mol_draw 的原子属性，而不是创建副本
        # 这样可以确保绘图和坐标提取使用的是同一个分子对象
        for atom in mol_draw.GetAtoms():
            # 清除形式电荷（不显示 +/-）
            atom.SetFormalCharge(0)
            # 清除显式氢原子计数，让 RDKit 只显示原子符号
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
            # 清除手性标签（不显示 R/S）
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        # 🔧 直接使用修改后的 mol_draw 绘制
        drawer.DrawMolecule(mol_draw)
        drawer.FinishDrawing()
    except Exception as e:
        print(f"[2D Diagram] Drawing error: {e}")
        return None

    # 提取渲染后的像素坐标 (用于连线)
    # 🔧 现在 drawer 绑定的就是 mol_draw，坐标索引一致
    px_atoms = {i: (drawer.GetDrawCoords(i).x, drawer.GetDrawCoords(i).y) for i in range(mol_draw.GetNumAtoms())}
    px_rings = {}
    for ri, ring in enumerate(mol_draw.GetRingInfo().AtomRings()):
        r_pts = [px_atoms[aid] for aid in ring if aid in px_atoms]
        if r_pts:
            px_rings[ri] = (sum(p[0] for p in r_pts)/len(r_pts), sum(p[1] for p in r_pts)/len(r_pts))

    # 5. Matplotlib 组装 (使用 io.BytesIO 优化内存)
    import io
    
    # 对于大型多肽配体，降低 DPI 以减少内存使用
    actual_dpi = 100 if is_large_peptide else dpi
    if is_large_peptide:
        print(f"[2D Diagram] 大型多肽配体：降低 DPI 为 {actual_dpi} 以减少内存使用")
    
    fig = plt.figure(figsize=(dw/100, dh/100), dpi=actual_dpi)
    ax = fig.add_subplot(111)
    ax.set_facecolor('white')
    
    try:
        with Image.open(io.BytesIO(drawer.GetDrawingText())) as lig_img:
            ax.imshow(lig_img, origin='upper', extent=[0, dw, dh, 0])
    except Exception as e:
        print(f"[2D Diagram] Image assembly error: {e}")
        plt.close(fig)
        return None
    
    ax.set_xlim(-50, dw + 50); ax.set_ylim(dh + 50, -50); ax.axis('off')

    # 计算配体中心
    c_x = sum(p[0] for p in px_atoms.values()) / len(px_atoms)
    c_y = sum(p[1] for p in px_atoms.values()) / len(px_atoms)
    
    # 绘制气泡与连线
    res_map = {}
    for inter in interactions:
        res = inter["res"]
        if res not in res_map: res_map[res] = []
        res_map[res].append(inter)
    
    # ========================================================================
    # Discovery Studio layout algorithm v5
    # Convex hull boundary (radial scan) + ray intersection placement
    # + clamped angle redistribution + multi-layer radius + smart anchor
    # v5 improvements (over v4):
    #   1. Larger standoff: min_standoff = max_boundary_radius + badge_radius + 80
    #   2. Larger badges (30px) and font (8pt) for better readability
    #   3. Wider angle gap: 30° min (≤12 residues), 20° min (more)
    #   4. Larger min_badge_dist (badge_radius*2 + 30) for clearer separation
    #   5. Center crossing threshold 0.7 for better avoidance routing
    #   6. MAX_ANGLE_DISPLACEMENT ±60° (was ±55°)
    #   7. Multi-layer threshold: full min_angle_gap (was 0.6x), multiplier 1.4 (was 1.45)
    #   8. Smart anchor threshold 0.6 (was 0.65) for more aggressive nearest-atom rerouting
    # ========================================================================
    placed_badges = []  # list of {"x": x, "y": y, "r": radius, "res": res_name}
    badge_radius = 30.0  # v5: enlarged to 30px for better label readability
    padding = 16.0  # extra padding between badges
    used_residue_types = set()

    sorted_res = sorted(res_map.items())

    # --- Step 1: Compute ligand convex-hull boundary (radial scan, no scipy) ---
    # Scan 360 angles (1° each) to find the farthest atom/ring center, forming boundary profile
    NUM_SCAN_ANGLES = 360
    boundary_radii = [0.0] * NUM_SCAN_ANGLES

    # Scan all atom pixel coordinates
    for _aid, (ax_p, ay_p) in px_atoms.items():
        dx_a = ax_p - c_x
        dy_a = ay_p - c_y
        r_a = math.hypot(dx_a, dy_a)
        if r_a < 1.0:
            continue
        angle_deg = int(math.degrees(math.atan2(dy_a, dx_a))) % 360
        # ±18° sector expansion (covers atom symbol display width, prevents gaps)
        for delta in range(-18, 19):
            idx = (angle_deg + delta) % NUM_SCAN_ANGLES
            if r_a > boundary_radii[idx]:
                boundary_radii[idx] = r_a

    # Scan all ring centers (rings are part of the molecular skeleton)
    for _rid, (rx_p, ry_p) in px_rings.items():
        dx_r = rx_p - c_x
        dy_r = ry_p - c_y
        r_r = math.hypot(dx_r, dy_r)
        if r_r < 1.0:
            continue
        angle_deg = int(math.degrees(math.atan2(dy_r, dx_r))) % 360
        for delta in range(-18, 19):
            idx = (angle_deg + delta) % NUM_SCAN_ANGLES
            if r_r > boundary_radii[idx]:
                boundary_radii[idx] = r_r

    # Smooth boundary (moving average to remove jaggedness)
    smoothed = list(boundary_radii)
    SMOOTH_WINDOW = 12
    for i in range(NUM_SCAN_ANGLES):
        total = 0.0
        for j in range(-SMOOTH_WINDOW, SMOOTH_WINDOW + 1):
            total += boundary_radii[(i + j) % NUM_SCAN_ANGLES]
        smoothed[i] = total / (2 * SMOOTH_WINDOW + 1)
    boundary_radii = smoothed

    # Global maximum boundary radius
    max_boundary_radius = max(boundary_radii) if boundary_radii else 100.0
    # v5: global minimum standoff = max_boundary_radius + badge_radius + 80
    # Ensures all badges stay well outside the molecular structure
    min_standoff_global = max_boundary_radius + badge_radius + 80

    def get_boundary_at_angle(angle_rad):
        """Get molecular boundary radius at specified angle (with floor value)."""
        deg = int(math.degrees(angle_rad)) % NUM_SCAN_ANGLES
        # Floor: 30% of max boundary radius, prevents zero boundary in sparse directions
        return max(boundary_radii[deg], max_boundary_radius * 0.3)

    def ray_boundary_distance(angle_rad):
        """Ray-boundary intersection: returns badge center distance along given angle."""
        bnd = get_boundary_at_angle(angle_rad)
        # v5: badge center = boundary radius + badge radius + safety margin (55px)
        local_standoff = bnd + badge_radius + 55
        # Must not be below global minimum standoff
        return max(local_standoff, min_standoff_global)

    print(f"\n[2D Layout] Boundary scan complete: max_boundary_radius={max_boundary_radius:.1f}px, "
          f"min_standoff_global={min_standoff_global:.1f}px")

    # --- Step 2: Pre-compute target angle and interaction points for each residue ---
    residue_layout_data = []  # [{res, inters, target_x, target_y, angle_rad}]

    for i, (res, inters) in enumerate(sorted_res):
        # Compute weighted average position of all interaction points for this residue
        tx_sum, ty_sum = 0.0, 0.0
        count = 0
        for inter in inters:
            tid, ttype = inter["idx"], inter["type"]
            if ttype == "ring" and tid in px_rings:
                tx, ty = px_rings[tid]
            elif tid in px_atoms:
                tx, ty = px_atoms[tid]
            else:
                continue
            tx_sum += tx
            ty_sum += ty
            count += 1

        if count == 0:
            print(f"[2D Layout] WARNING: residue {res} has no valid interaction points, skipping")
            continue

        target_x, target_y = tx_sum / count, ty_sum / count

        # Angle from molecular center to target point
        vx, vy = target_x - c_x, target_y - c_y
        if math.hypot(vx, vy) < 1.0:
            # Target nearly at center, use evenly-distributed default angle
            angle_rad = math.radians(i * (360.0 / max(len(sorted_res), 1)))
        else:
            angle_rad = math.atan2(vy, vx)

        residue_layout_data.append({
            "res": res,
            "inters": inters,
            "target_x": target_x,
            "target_y": target_y,
            "angle_rad": angle_rad,
            "natural_angle": angle_rad,  # v4: remember original angle for displacement limit
        })

    # --- Step 3: Angle redistribution with displacement limit (v4) ---
    # Sort by original angle
    residue_layout_data.sort(key=lambda d: d["angle_rad"])
    n_residues = len(residue_layout_data)

    # v4: Maximum angle displacement from natural position (radians)
    MAX_ANGLE_DISPLACEMENT = math.radians(60.0)  # ±60° max displacement

    if n_residues > 0:
        # v5: dynamic minimum angle gap — at least 30°, degrades no lower than 20°
        if n_residues <= 12:
            min_angle_gap_deg = 30.0  # 12 or fewer: at least 30°
        else:
            ideal_gap = 360.0 / n_residues
            min_angle_gap_deg = max(20.0, ideal_gap * 0.85)
        min_angle_gap = math.radians(min_angle_gap_deg)

        print(f"[2D Layout] n_residues={n_residues}, min_angle_gap={min_angle_gap_deg:.1f}°, "
              f"max_displacement={math.degrees(MAX_ANGLE_DISPLACEMENT):.0f}°")

        # Iterative angle redistribution with displacement clamp
        for _pass in range(30):
            adjusted = False
            for j in range(n_residues):
                next_j = (j + 1) % n_residues
                a1 = residue_layout_data[j]["angle_rad"]
                a2 = residue_layout_data[next_j]["angle_rad"]

                # Angle difference (considering 2π wrap)
                delta = (a2 - a1) % (2 * math.pi)
                if delta < 0:
                    delta += 2 * math.pi

                if 0 < delta < min_angle_gap:
                    push = (min_angle_gap - delta) / 2.0 * 0.7

                    # v4: Clamp displacement from natural angle
                    new_a1 = residue_layout_data[j]["angle_rad"] - push
                    nat1 = residue_layout_data[j]["natural_angle"]
                    disp1 = (new_a1 - nat1 + math.pi) % (2 * math.pi) - math.pi
                    if abs(disp1) > MAX_ANGLE_DISPLACEMENT:
                        new_a1 = nat1 + MAX_ANGLE_DISPLACEMENT * (1 if disp1 > 0 else -1)

                    new_a2 = residue_layout_data[next_j]["angle_rad"] + push
                    nat2 = residue_layout_data[next_j]["natural_angle"]
                    disp2 = (new_a2 - nat2 + math.pi) % (2 * math.pi) - math.pi
                    if abs(disp2) > MAX_ANGLE_DISPLACEMENT:
                        new_a2 = nat2 + MAX_ANGLE_DISPLACEMENT * (1 if disp2 > 0 else -1)

                    residue_layout_data[j]["angle_rad"] = new_a1
                    residue_layout_data[next_j]["angle_rad"] = new_a2
                    adjusted = True

            if not adjusted:
                break

        # v4: Mark residues that are still too close (need second radius layer)
        for j in range(n_residues):
            next_j = (j + 1) % n_residues
            a1 = residue_layout_data[j]["angle_rad"]
            a2 = residue_layout_data[next_j]["angle_rad"]
            delta = (a2 - a1) % (2 * math.pi)
            if delta < 0:
                delta += 2 * math.pi
            # If still too close after clamped redistribution, mark one for outer layer
            if 0 < delta < min_angle_gap:
                # Push the one with larger displacement to outer layer
                disp_j = abs((residue_layout_data[j]["angle_rad"] - residue_layout_data[j]["natural_angle"] + math.pi) % (2 * math.pi) - math.pi)
                disp_next = abs((residue_layout_data[next_j]["angle_rad"] - residue_layout_data[next_j]["natural_angle"] + math.pi) % (2 * math.pi) - math.pi)
                if disp_j > disp_next:
                    residue_layout_data[j]["outer_layer"] = True
                else:
                    residue_layout_data[next_j]["outer_layer"] = True

    # --- Step 4: Ray-boundary intersection placement (v4: multi-layer support) ---
    for data in residue_layout_data:
        angle = data["angle_rad"]
        ux = math.cos(angle)
        uy = math.sin(angle)

        # Ray-boundary intersection: get standoff distance for this direction
        standoff = ray_boundary_distance(angle)
        # v4: If marked for outer layer, push further out
        if data.get("outer_layer", False):
            standoff *= 1.4
            print(f"[2D Layout] {data['res']}: pushed to outer layer (standoff x1.40)")

        bx = c_x + ux * standoff
        by = c_y + uy * standoff

        data["bx"] = bx
        data["by"] = by
        data["ux"] = ux
        data["uy"] = uy

    # --- Step 5: Force-directed refinement (constraint: badges cannot enter convex hull, tangential or outward only) ---
    min_badge_dist = badge_radius * 2 + 30  # v5: minimum distance between badges, prevents label overlap
    min_atom_clearance = badge_radius + 50  # v5: minimum distance between badge and atom

    for iteration in range(100):
        any_moved = False

        for j, data in enumerate(residue_layout_data):
            bx, by = data["bx"], data["by"]
            # Accumulated force
            force_x, force_y = 0.0, 0.0

            # (a) Badge-badge repulsion
            for k, other in enumerate(residue_layout_data):
                if j == k:
                    continue
                obx, oby = other["bx"], other["by"]
                dx = bx - obx
                dy = by - oby
                d = math.hypot(dx, dy)

                if d < min_badge_dist:
                    if d < 1.0:
                        # Nearly coincident: use index to create a slight directional offset
                        dx = math.cos(j * 0.7 + k * 0.3)
                        dy = math.sin(j * 0.7 + k * 0.3)
                        d = 1.0
                    overlap = min_badge_dist - d
                    force_x += (dx / d) * overlap * 0.5
                    force_y += (dy / d) * overlap * 0.5

            # (b) Badge-atom repulsion
            for _aid, (ax_p, ay_p) in px_atoms.items():
                dx = bx - ax_p
                dy = by - ay_p
                d = math.hypot(dx, dy)

                if d < min_atom_clearance:
                    if d < 1.0:
                        # Coincident: push outward along radial direction
                        dx = bx - c_x
                        dy = by - c_y
                        d = math.hypot(dx, dy)
                        if d < 1.0:
                            dx, dy, d = 1.0, 0.0, 1.0
                    overlap = min_atom_clearance - d
                    force_x += (dx / d) * overlap * 1.0
                    force_y += (dy / d) * overlap * 1.0

            # Decompose accumulated force into radial and tangential components
            if abs(force_x) > 0.1 or abs(force_y) > 0.1:
                # Radial direction: from molecular center toward current badge
                rx = bx - c_x
                ry = by - c_y
                r_dist = math.hypot(rx, ry)
                if r_dist < 1.0:
                    r_dist = 1.0
                    rx, ry = 1.0, 0.0
                rad_ux, rad_uy = rx / r_dist, ry / r_dist

                # Tangential direction: perpendicular to radial
                tan_ux, tan_uy = -rad_uy, rad_ux

                # Radial component (outward only, inward push forbidden)
                radial_force = force_x * rad_ux + force_y * rad_uy
                if radial_force < 0:
                    radial_force = 0.0  # Forbid inward push

                # Tangential component (fully preserved, allows circumferential sliding)
                tangential_force = force_x * tan_ux + force_y * tan_uy

                bx += rad_ux * radial_force + tan_ux * tangential_force
                by += rad_uy * radial_force + tan_uy * tangential_force
                any_moved = True

            # (c) Hard constraint: badge center must not enter convex hull
            # v5: distance from center must not be less than max(boundary + badge_radius + 30, global min standoff)
            dist_from_center = math.hypot(bx - c_x, by - c_y)
            cur_angle = math.atan2(by - c_y, bx - c_x)
            min_dist = max(
                get_boundary_at_angle(cur_angle) + badge_radius + 30,
                min_standoff_global
            )
            if dist_from_center < min_dist:
                if dist_from_center < 1.0:
                    dist_from_center = 1.0
                scale = min_dist / dist_from_center
                bx = c_x + (bx - c_x) * scale
                by = c_y + (by - c_y) * scale
                any_moved = True

            data["bx"] = bx
            data["by"] = by

        if not any_moved:
            print(f"[2D Layout] Force-directed converged after {iteration + 1} iterations")
            break
    else:
        print(f"[2D Layout] Force-directed reached max iterations (100)")

    # --- Step 6: Build placed_badges and log results ---
    res_to_layout = {}
    for data in residue_layout_data:
        bx, by = data["bx"], data["by"]
        res = data["res"]
        placed_badges.append({"x": bx, "y": by, "r": badge_radius, "res": res})
        res_to_layout[res] = data

        print(f"[2D Layout] {res}: angle={math.degrees(data['angle_rad']):.1f}deg, "
              f"pos=({bx:.1f}, {by:.1f}), dist_from_center={math.hypot(bx-c_x, by-c_y):.1f}px")

    # Draw all residue badges and interaction lines
    for data in residue_layout_data:
        res = data["res"]
        inters = data["inters"]
        bx, by = data["bx"], data["by"]

        # Draw badge
        style, res_type = get_residue_style(res)
        used_residue_types.add(res_type)
        ax.add_patch(Circle((bx, by), badge_radius, facecolor=style['facecolor'], edgecolor=style['edgecolor'], lw=1.5, zorder=10))
        r_name, r_num = parse_residue_label(res)
        # v5: font 8pt, matching 30px badge
        ax.text(bx, by, f"{r_name}\n{r_num}", ha='center', va='center', fontweight='bold', fontsize=8, color=style['textcolor'], zorder=11)

        # Draw interaction lines (one per interaction)
        # v4: Smart anchor point selection - connect to nearest atom when target is too far
        mol_diameter = 2 * max_boundary_radius if max_boundary_radius > 0 else 400.0

        for inter in inters:
            tid, ttype = inter["idx"], inter["type"]
            if ttype == "ring" and tid in px_rings:
                tx, ty = px_rings[tid]
            elif tid in px_atoms:
                tx, ty = px_atoms[tid]
            else:
                continue

            normalized_itype = inter.get("normalized_itype", normalize_interaction_type(inter["itype"]))

            if normalized_itype == 'hydrophobic':
                # Hydrophobic interactions: no line drawn (Discovery Studio convention)
                # Badge color already indicates hydrophobic residue type
                continue

            # v4: Smart anchor - if target atom is too far from badge, connect to nearest atom instead
            badge_to_target = math.hypot(tx - bx, ty - by)
            far_threshold = mol_diameter * 0.6  # threshold for "too far"

            if badge_to_target > far_threshold and normalized_itype not in ['pipi', 'pication']:
                # Find the nearest ligand atom to the badge
                best_aid = None
                best_dist = float('inf')
                for aid, (ax_p, ay_p) in px_atoms.items():
                    d = math.hypot(ax_p - bx, ay_p - by)
                    if d < best_dist:
                        best_dist = d
                        best_aid = aid
                        tx, ty = ax_p, ay_p
                if best_aid is not None:
                    print(f"[2D Layout] {res}: rerouted line from atom {tid} to nearest atom {best_aid} "
                          f"(dist {badge_to_target:.0f} > {far_threshold:.0f})")

            ls = INTERACTION_LINE_STYLE.get(normalized_itype, INTERACTION_LINE_STYLE['other'])
            interaction_types_found.add(normalized_itype)

            # Compute clipped start point: line originates from badge edge
            angle_to_target = math.atan2(ty - by, tx - bx)
            sx = bx + badge_radius * math.cos(angle_to_target)
            sy = by + badge_radius * math.sin(angle_to_target)

            # Key fix: line endpoint is retracted from target atom (tx,ty) toward line start (sx,sy)
            # Not from badge center direction, but from actual line start direction
            atom_margin = 18  # approximate atom symbol radius (pixels)

            # Compute actual distance and direction from start to target
            line_dx = tx - sx
            line_dy = ty - sy
            line_dist = math.hypot(line_dx, line_dy)

            if line_dist > atom_margin + 5:
                # Retract from target toward start by atom_margin distance
                # Unit vector: from target pointing toward start
                ux = -line_dx / line_dist
                uy = -line_dy / line_dist
                # Endpoint = target + unit_vector * margin
                tx_end = tx + ux * atom_margin
                ty_end = ty + uy * atom_margin
            else:
                tx_end, ty_end = tx, ty

            # Check if the line would pass through other badges or atoms; if so, route around
            blocking_obstacles = []

            # Check if line passes through other badges
            for other_badge in placed_badges:
                if other_badge["res"] == res:
                    continue  # Skip self
                obx, oby, obr = other_badge["x"], other_badge["y"], other_badge["r"]

                # Line segment (sx, sy) -> (tx_end, ty_end)
                line_len = math.hypot(tx_end - sx, ty_end - sy)
                if line_len < 1.0:
                    continue

                # Compute shortest distance from point to line segment
                t_param = max(0, min(1, ((obx - sx) * (tx_end - sx) + (oby - sy) * (ty_end - sy)) / (line_len * line_len)))
                closest_x = sx + t_param * (tx_end - sx)
                closest_y = sy + t_param * (ty_end - sy)
                dist_to_line = math.hypot(obx - closest_x, oby - closest_y)

                # Check if line passes through badge (with margin)
                if dist_to_line < obr + 15:
                    blocking_obstacles.append({
                        "x": obx, "y": oby, "r": obr,
                        "t": t_param,  # position on line segment
                        "dist": dist_to_line,
                        "type": "badge"
                    })

            # Check if line passes through other atoms in the molecule
            atom_display_radius = 20  # atom symbol display radius (pixels)
            for atom_idx, (ax_pos, ay_pos) in px_atoms.items():
                # Skip the target atom itself
                if atom_idx == tid:
                    continue

                line_len = math.hypot(tx_end - sx, ty_end - sy)
                if line_len < 1.0:
                    continue

                # Compute shortest distance from atom to line segment
                t_param = max(0, min(1, ((ax_pos - sx) * (tx_end - sx) + (ay_pos - sy) * (ty_end - sy)) / (line_len * line_len)))

                # Only check middle portion of line (avoid false positives near endpoints)
                if t_param < 0.1 or t_param > 0.9:
                    continue

                closest_x = sx + t_param * (tx_end - sx)
                closest_y = sy + t_param * (ty_end - sy)
                dist_to_line = math.hypot(ax_pos - closest_x, ay_pos - closest_y)

                # Check if line passes through atom symbol
                if dist_to_line < atom_display_radius:
                    blocking_obstacles.append({
                        "x": ax_pos, "y": ay_pos, "r": atom_display_radius,
                        "t": t_param,
                        "dist": dist_to_line,
                        "type": "atom"
                    })

            # v3: Check if line crosses through molecular center region
            # If midpoint of line is too close to center, the line crosses the molecule
            mid_line_x = (sx + tx_end) / 2
            mid_line_y = (sy + ty_end) / 2
            mid_to_center = math.hypot(mid_line_x - c_x, mid_line_y - c_y)
            center_threshold = max_boundary_radius * 0.7  # v5: molecular center region threshold
            if mid_to_center < center_threshold and not blocking_obstacles:
                # Line crosses molecular center, add virtual obstacle to force routing
                blocking_obstacles.append({
                    "x": c_x, "y": c_y, "r": max_boundary_radius * 0.5,
                    "t": 0.5,
                    "dist": mid_to_center,
                    "type": "center"
                })

            if blocking_obstacles:
                # Need to route around obstacles
                # Strategy: use quadratic Bezier curve with control point on the outside of the obstacle
                # Find the obstacle closest to the line midpoint
                blocking_obstacles.sort(key=lambda b: abs(b["t"] - 0.5))
                main_blocker = blocking_obstacles[0]

                # Compute routing direction: perpendicular to line direction
                line_dx, line_dy = tx_end - sx, ty_end - sy
                line_len = math.hypot(line_dx, line_dy)
                perp_x, perp_y = -line_dy / line_len, line_dx / line_len

                # Determine routing side (prefer direction away from molecular center)
                # v3: prioritize perpendicular direction farther from molecular center
                blocker_to_closest_x = main_blocker["x"] - (sx + main_blocker["t"] * (tx_end - sx))
                blocker_to_closest_y = main_blocker["y"] - (sy + main_blocker["t"] * (ty_end - sy))

                # Determine which perpendicular direction is farther from molecular center
                mid_x = (sx + tx_end) / 2
                mid_y = (sy + ty_end) / 2
                # Direction 1: (perp_x, perp_y)
                test1_x = mid_x + perp_x * 10
                test1_y = mid_y + perp_y * 10
                dist1_to_center = math.hypot(test1_x - c_x, test1_y - c_y)
                # Direction 2: (-perp_x, -perp_y)
                test2_x = mid_x - perp_x * 10
                test2_y = mid_y - perp_y * 10
                dist2_to_center = math.hypot(test2_x - c_x, test2_y - c_y)

                # Choose the direction farther from molecular center
                if dist1_to_center < dist2_to_center:
                    perp_x, perp_y = -perp_x, -perp_y

                # Control point: offset from line midpoint along perpendicular direction
                mid_x = (sx + tx_end) / 2
                mid_y = (sy + ty_end) / 2
                offset_dist = main_blocker["r"] + 50  # v5: routing offset distance
                ctrl_x = mid_x + perp_x * offset_dist
                ctrl_y = mid_y + perp_y * offset_dist

                # Draw quadratic Bezier curve
                # matplotlib supports Path/PathPatch, but for simplicity we approximate with polyline

                # Generate curve points
                curve_points = []
                for t in [i / 20.0 for i in range(21)]:
                    # Quadratic Bezier curve formula
                    px = (1 - t) ** 2 * sx + 2 * (1 - t) * t * ctrl_x + t ** 2 * tx_end
                    py = (1 - t) ** 2 * sy + 2 * (1 - t) * t * ctrl_y + t ** 2 * ty_end
                    curve_points.append((px, py))

                # Draw curve
                xs, ys = zip(*curve_points)
                ax.plot(xs, ys, color=ls['color'], lw=ls['linewidth'] + 0.5,
                       ls=ls['linestyle'], zorder=2, alpha=0.9)
            else:
                # No obstacles, draw straight line (using retracted endpoint)
                ax.plot([sx, tx_end], [sy, ty_end], color=ls['color'], lw=ls['linewidth'] + 0.5,
                       ls=ls['linestyle'], zorder=2, alpha=0.9)

    # Legend: Discovery Studio style
    # 1) Interaction type lines (Hydrogen Bond, Pi-Pi, Pi-Cation, etc.)
    # 2) Residue type badges (Hydrophobic / Nonpolar / Polar / Negative / Positive)
    legend_handles = []
    
    # Interaction type legend (in priority order, all 7 interaction types)
    interaction_legend_order = ['hbond', 'salt', 'pipi', 'pication', 'halogen', 'metal', 'water']
    interaction_labels = {
        'hbond': 'Hydrogen Bond',
        'salt': 'Salt Bridge',
        'pipi': 'Pi-Pi Stacking',
        'pication': 'Pi-Cation',
        'halogen': 'Halogen Bond',
        'metal': 'Metal Coordination',
        'water': 'Water Bridge',
    }
    
    for itype in interaction_legend_order:
        if itype in interaction_types_found:
            s = INTERACTION_LINE_STYLE.get(itype, INTERACTION_LINE_STYLE['other'])
            legend_handles.append(
                Line2D([0], [0], color=s['color'], lw=2, ls=s.get('linestyle', '-'), label=interaction_labels.get(itype, itype))
            )
    
    # Residue type legend
    residue_labels = {
        'hydrophobic': 'Hydrophobic',
        'nonpolar': 'Nonpolar',
        'polar': 'Polar',
        'negative': 'Negative',
        'positive': 'Positive',
    }
    # Fixed display order for consistent legend
    residue_order = ['hydrophobic', 'nonpolar', 'polar', 'negative', 'positive']
    for r_type in residue_order:
        if r_type in used_residue_types:
            style = DS_RESIDUE_STYLE[r_type]
            legend_handles.append(
                mpatches.Circle(
                    (0, 0),
                    radius=6,
                    facecolor=style['facecolor'],
                    edgecolor=style['edgecolor'],
                    linewidth=1.5,
                    label=residue_labels[r_type],
                )
            )
    
    if legend_handles:
        ax.legend(
            handles=legend_handles,
            loc='upper left',
            frameon=True,
            fontsize=8,
            fancybox=True,
            framealpha=0.9,
            edgecolor='gray',
        )

    if not output_path:
        output_path = os.path.join(os.path.expanduser("~"), "Desktop", f"{ligand_resname}_2d.png")
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    print(f"[2D Diagram] ✅ Saved to {output_path}")
    return output_path