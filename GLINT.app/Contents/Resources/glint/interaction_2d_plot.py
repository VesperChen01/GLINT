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
# 相互作用类型映射 (中英文统一，与 3D 视图一致)
# ============================================================================
INTERACTION_TYPE_MAP = {
    # 中文 -> 内部类型
    '氢键': 'hbond',
    '盐桥': 'salt',
    'π–π 堆积': 'pipi',
    'π–阳离子相互作用': 'pication',
    '疏水相互作用': 'hydrophobic',
    '卤素键': 'halogen',
    '金属配位': 'metal',
    '水桥': 'water',
    # 英文 -> 内部类型
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
    
    # 然后尝试部分匹配 (按优先级顺序)
    # 氢键
    if '氢键' in itype or 'hbond' in itype_lower or 'hydrogen' in itype_lower or 'h-bond' in itype_lower:
        return 'hbond'
    # 盐桥
    if '盐桥' in itype or 'salt' in itype_lower:
        return 'salt'
    # π-π 堆积
    if 'π–π' in itype or 'π-π' in itype or 'pi-pi' in itype_lower or 'pipi' in itype_lower or 'stacking' in itype_lower:
        return 'pipi'
    # π-阳离子
    if 'π–阳离子' in itype or 'π-阳离子' in itype or 'pi-cation' in itype_lower or 'pication' in itype_lower or 'cation-pi' in itype_lower:
        return 'pication'
    # 金属配位
    if '金属' in itype or 'metal' in itype_lower:
        return 'metal'
    # 水桥
    if '水桥' in itype or 'water' in itype_lower:
        return 'water'
    # 卤素键
    if '卤素' in itype or 'halogen' in itype_lower:
        return 'halogen'
    # 疏水相互作用 (最后检查，因为很多类型名称可能包含相关词)
    if '疏水' in itype or 'hydrophobic' in itype_lower or 'alkyl' in itype_lower:
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
                except:
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
                            except:
                                pass
                        if len(new_parts) > 2:
                            final_lines.append(" ".join(new_parts))
                    except:
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
                                    try:
                                        AllChem.Compute2DCoords(mol_from_smiles)

                                        bond_types = {}
                                        for bond in mol_from_smiles.GetBonds():
                                            bt = str(bond.GetBondType())
                                            bond_types[bt] = bond_types.get(bt, 0) + 1

                                        print(f"[2D Diagram] ✅ 使用 Open Babel SMILES 分子（键级正确）: {bond_types}")

                                        if tmp_smi and os.path.exists(tmp_smi):
                                            os.remove(tmp_smi)

                                        return mol_from_smiles
                                    except:
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
                    except:
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
                except:
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
    except:
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
                        except: pass

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
                        except:
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
                                    except:
                                        pass
                                if len(new_parts) > 2:
                                    filtered_conect.append(" ".join(new_parts))
                        except:
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
                        except: pass

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
        except:
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

    # 4. 生成 2D 布局 (使用 RDKit CoordGen)
    try:
        # 确保使用 CoordGen，而不是基于当前 3D 构象的小幅调整
        rdDepictor.SetPreferCoordGen(True)
        # 🔧 使用 coordGen 时尝试更宽松的参数
        try:
            from rdkit.Chem import rdCoordGen
            # 设置 CoordGen 参数以获得更好的布局
            rdCoordGen.AddCoords(mol_draw)
            print("[2D Diagram] Using rdCoordGen.AddCoords for layout")
        except Exception:
            rdDepictor.Compute2DCoords(mol_draw)
    except Exception as e:
        print(f"[2D Diagram] Layout error (CoordGen failed, falling back to RDKit 2D): {e}")
        try:
            AllChem.Compute2DCoords(mol_draw)
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
    
    # 智能布局算法 (基于相互作用点的局部扩展)
    placed_badges = [] # list of {"x": x, "y": y, "r": radius, "res": res_name}
    badge_radius = 22.0 # 再次缩小气泡半径（从 28 减小到 22）
    padding = 12.0  # 减小气泡之间的间距（从 15 减小到 12）
    used_residue_types = set()  # 用于构建 Discovery Studio 风格的残基图例
    
    sorted_res = sorted(res_map.items())

    for i, (res, inters) in enumerate(sorted_res):
        # 1. 计算该残基所有相互作用点的平均位置 (目标中心)
        tx_sum, ty_sum = 0, 0
        count = 0
        interaction_points = []  # 调试用
        for inter in inters:
            tid, ttype = inter["idx"], inter["type"]
            if ttype == "ring" and tid in px_rings:
                tx, ty = px_rings[tid]
                interaction_points.append((tid, ttype, tx, ty))
            elif tid in px_atoms:
                tx, ty = px_atoms[tid]
                interaction_points.append((tid, ttype, tx, ty))
            else:
                # 🔧 修复：如果找不到原子索引，尝试跳过而不是忽略整个残基
                print(f"  ⚠️ 警告：残基 {res} 的相互作用 {ttype} 索引 {tid} 未找到对应原子")
                continue
            tx_sum += tx
            ty_sum += ty
            count += 1

        if count == 0:
            print(f"  ⚠️ 警告：残基 {res} 没有有效的相互作用点，跳过")
            continue

        target_x, target_y = tx_sum / count, ty_sum / count

        # 调试输出：显示相互作用点
        print(f"\n残基 {res} 的相互作用点:")
        for tid, ttype, tx, ty in interaction_points:
            print(f"  - {ttype} {tid}: ({tx:.1f}, {ty:.1f})")
        print(f"  配体中心: ({c_x:.1f}, {c_y:.1f}), 目标点: ({target_x:.1f}, {target_y:.1f})")

        # 2. 计算从配体中心指向目标中心的向量
        vx, vy = target_x - c_x, target_y - c_y
        dist = math.hypot(vx, vy)

        # 计算配体的边界框（用于后续的方向判断）
        min_y_mol = min(p[1] for p in px_atoms.values())
        max_y_mol = max(p[1] for p in px_atoms.values())
        min_x_mol = min(p[0] for p in px_atoms.values())
        max_x_mol = max(p[0] for p in px_atoms.values())
        mol_width = max(max_x_mol - min_x_mol, 100)
        mol_height = max(max_y_mol - min_y_mol, 100)

        # 计算配体的最大半径（从中心到最远原子）
        max_radius = 0
        for atom_idx, (atom_x, atom_y) in px_atoms.items():
            r = math.hypot(atom_x - c_x, atom_y - c_y)
            if r > max_radius:
                max_radius = r

        # 🔧 新策略：气泡放在目标原子的「外侧」（远离分子中心的方向）
        # 这样连线就会从分子边缘延伸出去，而不是穿过分子

        # 打印调试信息
        mol_center_y = (min_y_mol + max_y_mol) / 2
        print(f"  分子边界: y=[{min_y_mol:.1f}, {max_y_mol:.1f}], x=[{min_x_mol:.1f}, {max_x_mol:.1f}]")
        print(f"  目标位置: ({target_x:.1f}, {target_y:.1f}), 分子中心: ({c_x:.1f}, {c_y:.1f})")

        # 增加边距，让气泡远离分子边界
        margin = 50  # 50 像素边距

        # 🔧 核心改变：总是沿着「目标原子相对于分子中心」的方向放置气泡
        # 这样气泡就会在目标原子的外侧，连线不会穿过分子
        if dist < 1.0:
            # 距离太近，选择默认方向
            default_angles = [270, 180, 0, 90, 225, 315, 45, 135]
            angle_idx = len(placed_badges) % len(default_angles)
            angle_deg = default_angles[angle_idx]
            angle_rad = math.radians(angle_deg)
            ux, uy = math.cos(angle_rad), math.sin(angle_rad)
            print(f"  方向：默认角度 {angle_deg}° (距离太近)")
        else:
            # 气泡方向 = 从分子中心指向目标原子的方向（向外延伸）
            ux, uy = vx / dist, vy / dist
            print(f"  方向：外侧延伸 ({ux:.2f}, {uy:.2f})")

        # 不使用边界固定位置，使用 standoff 方式
        use_boundary_placement = False
        bx_start, by_start = 0, 0  # 占位符

        # 🔧 关键修改：缩短 standoff 距离，让气泡更靠近目标原子
        # 这样连线就不会太长，不容易穿过分子
        # standoff = 从目标原子位置 + 一点点距离（而不是从分子中心 + 最大半径）

        # 计算目标点到分子中心的距离
        target_dist_from_center = math.hypot(target_x - c_x, target_y - c_y)

        # 气泡位置 = 目标点位置 + 向外延伸一小段距离
        # 这样连线从气泡到目标原子只有很短的距离
        extra_offset = badge_radius + 40  # 气泡半径 + 40 像素间距
        standoff = target_dist_from_center + extra_offset

        # ⚠️ 关键改进：扇形分布
        # 如果多个残基与同一个原子相互作用，需要在角度上分散它们
        # 检查已放置的气泡中，有多少个与相同的目标点相互作用
        angle_offset = 0.0
        similar_count = 0

        for placed in placed_badges:
            # 计算已放置气泡的目标点（反推）
            placed_vx = placed["x"] - c_x
            placed_vy = placed["y"] - c_y
            placed_dist = math.hypot(placed_vx, placed_vy)

            if placed_dist > 1.0:
                placed_ux = placed_vx / placed_dist
                placed_uy = placed_vy / placed_dist

                # 计算方向相似度（点积）
                dot = ux * placed_ux + uy * placed_uy

                # 如果方向非常接近（cos > 0.85，即角度 < 32°）
                # 说明它们可能与同一个原子或非常接近的原子相互作用
                if dot > 0.85:
                    similar_count += 1

        # 根据相似气泡的数量，计算角度偏移
        if similar_count > 0:
            # 每个相似气泡增加 45° 的偏移（从 30° 增加到 45°）
            angle_offset = similar_count * (math.pi / 4.0)  # 45° = π/4

            # 应用角度偏移（旋转方向向量）
            cos_offset = math.cos(angle_offset)
            sin_offset = math.sin(angle_offset)
            ux_new = ux * cos_offset - uy * sin_offset
            uy_new = ux * sin_offset + uy * cos_offset
            ux, uy = ux_new, uy_new

            print(f"  ⚠️ 检测到 {similar_count} 个相似方向的气泡，应用 {math.degrees(angle_offset):.1f}° 角度偏移")

        # 计算气泡位置
        if use_boundary_placement:
            # 使用边界计算的位置
            bx = bx_start
            by = by_start
        else:
            # 使用 standoff 计算位置
            bx = c_x + ux * standoff
            by = c_y + uy * standoff

        # 调试输出
        print(f"  配体中心: ({c_x:.1f}, {c_y:.1f})")
        print(f"  目标点平均: ({target_x:.1f}, {target_y:.1f})")
        print(f"  方向向量: ({ux:.2f}, {uy:.2f})")
        print(f"  max_radius={max_radius:.1f}, standoff={standoff:.1f}")
        print(f"  初始气泡位置: ({bx:.1f}, {by:.1f})")
        print(f"  气泡到中心距离: {math.hypot(bx - c_x, by - c_y):.1f}")
        
        # 4. 碰撞检测与解决 (增强的迭代推挤)
        # 增加迭代次数，确保充分分离
        for iteration in range(50):  # 从 30 增加到 50
            moved = False
            # 5. 增强的碰撞检测：同时避开其他气泡和配体原子
            # 检查与已放置气泡的重叠
            for badge in placed_badges:
                dx = bx - badge["x"]
                dy = by - badge["y"]
                d = math.hypot(dx, dy)
                min_dist = badge_radius * 2 + padding + 10  # 额外增加 10 像素间距

                if d < min_dist:
                    # 发生重叠，推开 (气泡之间互斥)
                    if d < 1.0: dx, dy, d = 1.0, 0.0, 1.0
                    overlap = min_dist - d
                    # 增加推挤力度，确保充分分离
                    push_factor = 1.5 if iteration < 15 else 1.0  # 前期强力推挤，后期微调
                    push_x = (dx / d) * overlap * push_factor
                    push_y = (dy / d) * overlap * push_factor
                    bx += push_x
                    by += push_y
                    moved = True

            # 检查与配体原子的重叠 (关键修复！)
            # 遍历所有配体原子，确保气泡不覆盖它们
            # 为了效率，可以只检查凸包或者采样点，但原子数不多，直接遍历即可
            min_atom_dist = badge_radius + 40  # 气泡半径 + 安全边距（从 20 增加到 40）
            
            
            for aid, (atom_x, atom_y) in px_atoms.items():
                dx = bx - atom_x
                dy = by - atom_y
                d = math.hypot(dx, dy)
                
                if d < min_atom_dist:
                    # 严重重叠！必须强力推开
                    # 推开方向：从配体中心向气泡方向 (ux, uy)
                    # 或者从原子向气泡方向
                    if d < 1.0: 
                        # 如果正好重合，沿径向推
                        push_dir_x, push_dir_y = ux, uy
                    else:
                        push_dir_x, push_dir_y = dx/d, dy/d
                    
                    overlap = min_atom_dist - d
                    bx += push_dir_x * overlap * 1.5 # 推挤系数从 1.2 增加到 1.5
                    by += push_dir_y * overlap * 1.5
                    moved = True
            
            # 额外的几何约束：确保气泡不"陷入"配体内部 concave 区域并遮挡内部原子
            # 检查气泡中心到配体中心的距离，不能小于某个阈值(动态)
            # 或者简单地：气泡必须比它所作用的原子更靠外
            
            # 投影检查 (气泡必须位于目标原子沿径向的外侧)
            vec_c_target = (target_x - c_x, target_y - c_y) # 中心->目标原子
            vec_c_badge = (bx - c_x, by - c_y) # 中心->气泡
            
            # 计算气泡在 (中心->目标) 方向上的投影长度
            # 如果气泡跑到了目标原子的内侧 (投影长度 < 目标投影长度)，这是不合理的
            proj_len = (vec_c_badge[0]*ux + vec_c_badge[1]*uy)
            target_len = (vec_c_target[0]*ux + vec_c_target[1]*uy)
            
            min_proj = target_len + badge_radius * 0.8 # 至少要在目标原子外侧
            if proj_len < min_proj:
                diff = min_proj - proj_len
                bx += ux * diff
                by += uy * diff
                moved = True

            if not moved:
                break
        
        placed_badges.append({"x": bx, "y": by, "r": badge_radius, "res": res})

        # 绘制
        style, res_type = get_residue_style(res)
        used_residue_types.add(res_type)
        ax.add_patch(Circle((bx, by), badge_radius, facecolor=style['facecolor'], edgecolor=style['edgecolor'], lw=1.2, zorder=10))
        r_name, r_num = parse_residue_label(res)
        # 字体再小一点以适应更小的气泡
        ax.text(bx, by, f"{r_name}\n{r_num}", ha='center', va='center', fontweight='bold', fontsize=6.5, color=style['textcolor'], zorder=11)

        # 绘制相互作用连线（每个相互作用都画一条线）
        # 关键改进：连线绕开其他气泡（使用贝塞尔曲线或折线）
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
                # 疏水相互作用不画连线
                continue

            ls = INTERACTION_LINE_STYLE.get(normalized_itype, INTERACTION_LINE_STYLE['other'])
            interaction_types_found.add(normalized_itype)

            # 计算截断点：从气泡边缘发出
            angle_to_target = math.atan2(ty - by, tx - bx)
            sx = bx + badge_radius * math.cos(angle_to_target)
            sy = by + badge_radius * math.sin(angle_to_target)

            # 🔧 关键修复：连线终点从目标原子 (tx,ty) 往线起点 (sx,sy) 方向缩进
            # 不是从气泡中心方向，而是从实际连线起点方向
            atom_margin = 18  # 原子符号的大约半径（像素）

            # 计算从起点到目标的实际距离和方向
            line_dx = tx - sx
            line_dy = ty - sy
            line_dist = math.hypot(line_dx, line_dy)

            if line_dist > atom_margin + 5:
                # 从目标点往起点方向缩进 atom_margin 距离
                # 单位向量：从目标指向起点
                ux = -line_dx / line_dist
                uy = -line_dy / line_dist
                # 终点 = 目标点 + 单位向量 * margin
                tx_end = tx + ux * atom_margin
                ty_end = ty + uy * atom_margin
            else:
                tx_end, ty_end = tx, ty

            # 检查连线是否会穿过其他气泡或原子，如果是则添加中间控制点绕开
            blocking_obstacles = []

            # 🔧 检查是否穿过其他气泡
            for other_badge in placed_badges:
                if other_badge["res"] == res:
                    continue  # 跳过自己
                obx, oby, obr = other_badge["x"], other_badge["y"], other_badge["r"]

                # 线段 (sx, sy) -> (tx_end, ty_end)
                line_len = math.hypot(tx_end - sx, ty_end - sy)
                if line_len < 1.0:
                    continue

                # 计算点到线段的最短距离
                t_param = max(0, min(1, ((obx - sx) * (tx_end - sx) + (oby - sy) * (ty_end - sy)) / (line_len * line_len)))
                closest_x = sx + t_param * (tx_end - sx)
                closest_y = sy + t_param * (ty_end - sy)
                dist_to_line = math.hypot(obx - closest_x, oby - closest_y)

                # 检查是否穿过气泡（包含一定边距）
                if dist_to_line < obr + 15:
                    blocking_obstacles.append({
                        "x": obx, "y": oby, "r": obr,
                        "t": t_param,  # 在线段上的位置
                        "dist": dist_to_line,
                        "type": "badge"
                    })

            # 🔧 新增：检查是否穿过分子中的其他原子
            atom_display_radius = 20  # 原子符号的显示半径（像素）
            for atom_idx, (ax_pos, ay_pos) in px_atoms.items():
                # 跳过目标原子本身
                if atom_idx == tid:
                    continue

                line_len = math.hypot(tx_end - sx, ty_end - sy)
                if line_len < 1.0:
                    continue

                # 计算原子到线段的最短距离
                t_param = max(0, min(1, ((ax_pos - sx) * (tx_end - sx) + (ay_pos - sy) * (ty_end - sy)) / (line_len * line_len)))

                # 只检查线段中间部分（避免起点和终点附近的误判）
                if t_param < 0.1 or t_param > 0.9:
                    continue

                closest_x = sx + t_param * (tx_end - sx)
                closest_y = sy + t_param * (ty_end - sy)
                dist_to_line = math.hypot(ax_pos - closest_x, ay_pos - closest_y)

                # 检查是否穿过原子符号
                if dist_to_line < atom_display_radius:
                    blocking_obstacles.append({
                        "x": ax_pos, "y": ay_pos, "r": atom_display_radius,
                        "t": t_param,
                        "dist": dist_to_line,
                        "type": "atom"
                    })

            if blocking_obstacles:
                # 需要绕开障碍物
                # 策略：使用二次贝塞尔曲线，控制点在障碍物外侧
                # 找出最靠近线段中点的障碍物
                blocking_obstacles.sort(key=lambda b: abs(b["t"] - 0.5))
                main_blocker = blocking_obstacles[0]

                # 计算绕开方向：垂直于线段方向
                line_dx, line_dy = tx_end - sx, ty_end - sy
                line_len = math.hypot(line_dx, line_dy)
                perp_x, perp_y = -line_dy / line_len, line_dx / line_len

                # 确定绕开方向（选择远离其他障碍物的方向）
                # 简单策略：使用气泡中心到线段最近点的方向的反方向
                blocker_to_closest_x = main_blocker["x"] - (sx + main_blocker["t"] * (tx_end - sx))
                blocker_to_closest_y = main_blocker["y"] - (sy + main_blocker["t"] * (ty_end - sy))

                # 选择垂直方向中远离障碍物的那个
                if blocker_to_closest_x * perp_x + blocker_to_closest_y * perp_y > 0:
                    perp_x, perp_y = -perp_x, -perp_y

                # 控制点：在线段中点沿垂直方向偏移
                mid_x = (sx + tx_end) / 2
                mid_y = (sy + ty_end) / 2
                offset_dist = main_blocker["r"] + 40  # 绕开距离
                ctrl_x = mid_x + perp_x * offset_dist
                ctrl_y = mid_y + perp_y * offset_dist

                # 使用二次贝塞尔曲线绘制
                # matplotlib 可以用 Path 和 PathPatch，但为简化，用多段折线近似

                # 生成曲线点
                curve_points = []
                for t in [i / 20.0 for i in range(21)]:
                    # 二次贝塞尔曲线公式
                    px = (1 - t) ** 2 * sx + 2 * (1 - t) * t * ctrl_x + t ** 2 * tx_end
                    py = (1 - t) ** 2 * sy + 2 * (1 - t) * t * ctrl_y + t ** 2 * ty_end
                    curve_points.append((px, py))

                # 绘制曲线
                xs, ys = zip(*curve_points)
                ax.plot(xs, ys, color=ls['color'], lw=ls['linewidth'] + 0.5,
                       ls=ls['linestyle'], zorder=2, alpha=0.9)
            else:
                # 无障碍物，直接画直线（使用缩进后的终点）
                ax.plot([sx, tx_end], [sy, ty_end], color=ls['color'], lw=ls['linewidth'] + 0.5,
                       ls=ls['linestyle'], zorder=2, alpha=0.9)

    # 图例：Discovery Studio 风格
    # 1) 相互作用类型连线 (氢键、Pi-Pi、Pi-阳离子等)
    # 2) 各类残基类型气泡 (Hydrophobic / Nonpolar / Polar / Negative / Positive)
    legend_handles = []
    
    # 相互作用类型图例（按优先级顺序）
    interaction_legend_order = ['hbond', 'pipi', 'pication', 'salt']
    interaction_labels = {
        'hbond': 'Hydrogen Bond',
        'pipi': 'Pi-Pi Stacking',
        'pication': 'Pi-Cation',
        'salt': 'Salt Bridge',
    }
    
    for itype in interaction_legend_order:
        if itype in interaction_types_found:
            s = INTERACTION_LINE_STYLE.get(itype, INTERACTION_LINE_STYLE['other'])
            legend_handles.append(
                Line2D([0], [0], color=s['color'], lw=2, ls=s.get('linestyle', '-'), label=interaction_labels.get(itype, itype))
            )
    
    # 残基类型图例
    residue_labels = {
        'hydrophobic': 'Hydrophobic',
        'nonpolar': 'Nonpolar',
        'polar': 'Polar',
        'negative': 'Negative',
        'positive': 'Positive',
    }
    # 按固定顺序展示，使图例稳定
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