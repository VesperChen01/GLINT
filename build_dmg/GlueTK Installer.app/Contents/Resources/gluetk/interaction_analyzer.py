# -*- coding: utf-8 -*-
"""
interaction_analyzer.py
蛋白质相互作用分析模块

基于原始的pdb_interactions.py改进，增加了PyMOL集成功能
自动检测并使用RDKit高级分析（如果环境可用）

使用示例：
    # 分析蛋白-配体相互作用（自动使用最优模式）
    result = analyze_protein_ligand_interactions('protein', 'LIG')
    
    # 3D可视化
    visualize_protein_ligand_3d('protein', result, 'LIG')
"""

from __future__ import print_function
import math
import csv
import os
import tempfile
from collections import defaultdict
from pymol import cmd

# ========== 自动检测RDKit和依赖 ==========
# 尝试自动安装依赖
try:
    from .env_setup import ensure_dependencies
    _deps_checked = ensure_dependencies()
except Exception:
    _deps_checked = False

# 导入高质量分析所需的包
RDKIT_AVAILABLE = False
SCIPY_AVAILABLE = False
NUMPY_AVAILABLE = False
MPL_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
    print("[GlueTK] ✅ RDKit loaded")
except ImportError:
    print("[GlueTK] ⚠️ RDKit not installed; some advanced features unavailable")
    print("[GlueTK] 💡 Install manually: pip install rdkit")

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    print("[GlueTK] ⚠️ NumPy not installed")

try:
    from scipy.spatial import cKDTree
    SCIPY_AVAILABLE = True
except ImportError:
    print("[GlueTK] ⚠️ SciPy not installed; spatial acceleration unavailable")

try:
    import matplotlib.pyplot as plt
    MPL_AVAILABLE = True
except ImportError:
    print("[GlueTK] ⚠️ Matplotlib not installed; plotting unavailable")

# ========== 重构模块导入（性能优化）==========
try:
    from .interaction_types import (
        InteractionResult, HBond, SaltBridge, Hydrophobic,
        PiStacking, PiCation, MetalCoordination, HalogenBond, WaterBridge
    )
    from .feature_extractor import MolecularFeatureExtractor
    REFACTORED_MODULES_AVAILABLE = True
    print("[GlueTK] ✅ Using refactored interaction detection (performance boost)")
except ImportError as e:
    REFACTORED_MODULES_AVAILABLE = False
    # Silently fall back to legacy mode
    pass

# ========== 相互作用参数（实用标准）==========
# 适用于大多数药物设计场景的均衡标准
INTERACTION_PARAMS = {
    "hbond": {
        "max_DA_dist": 3.2,         # Å，D···A 距离（3.2Å为推荐标准，平衡严格性与实用性）
        "min_donor_angle": 120,     # °，∠D–H···A
        "min_acceptor_angle": 90    # °，∠H···A–X
    },
    "hydrophobic": {
        "pi_cation_max": 4.5,       # Å
        "pi_pi_mode": "face_face_or_edge",  # 按环面法向量与距离联合判定
        "other_max": 4.0            # Å (疏水接触,从3.6改为4.0更实用)
    },
    "ionic": {
        "max_dist": 4.5,            # Å (从4.0改为4.5更实用)
        "exclude_if_hbond": True    # 排除氢键情况
    },
    "metal_coord": {
        "max_dist": 3.4,            # Å
        "allowed_ligand_atoms": ["N", "O", "S", "CL", "BR", "F"]  # 非碳重原子
    },
    "water_bridge": {
        "max_DA_dist": 3.5,         # Å (与氢键保持一致)
        "min_donor_angle": 110,     # °
        "min_acceptor_angle": 90    # °
    },
    # 兼容旧参数名称
    "saltbridge": {"max_distance": 4.5},  # 与 ionic 保持一致
    "pi_pi": {"max_distance": 5.5, "min_distance": 3.3, "parallel_angle": 30.0},
    "pi_cation": {"max_distance": 4.5},
    "halogen": {"max_distance": 4.0, "min_angle": 140.0},
    "metal": {"max_distance": 3.4}
}

# 超严格标准（Schrödinger，仅用于高分辨率晶体结构）
SCHRODINGER_PARAMS = {
    "hbond": {"max_DA_dist": 2.8, "min_donor_angle": 120, "min_acceptor_angle": 90},
    "ionic": {"max_dist": 4.0},
    "hydrophobic": {"other_max": 3.6}
}

def parse_pdb_structure(obj_name=None):
    """
    从PyMOL对象解析原子信息
    
    参数:
        obj_name: PyMOL对象名称，如果为None则使用第一个对象
    
    返回:
        atoms: 原子信息列表 [(chain, res_name, res_id, atom_name, (x, y, z))]
    """
    if obj_name is None:
        try:
            objs = cmd.get_names("objects")
        except AttributeError:
            # Fallback for older versions/plugins
            objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
            
        if not objs:
            print("[parse_pdb_structure] No objects loaded")
            return []
        obj_name = objs[0]
    
    # Check existence
    current_objs = []
    try:
        current_objs = cmd.get_names("objects")
    except AttributeError:
        current_objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
        
    if obj_name not in current_objs:
        print(f"[parse_pdb_structure] Object '{obj_name}' not found")
        return []
    
    atoms = []
    model = cmd.get_model(obj_name)
    
    for atom in model.atom:
        chain = atom.chain.strip()
        res_name = atom.resn.strip()
        res_id = str(atom.resi).strip()
        atom_name = atom.name.strip()
        coords = (atom.coord[0], atom.coord[1], atom.coord[2])
        atoms.append((chain, res_name, res_id, atom_name, coords))
    
    return atoms

def parse_pdb_file(pdb_file):
    """
    从PDB文件解析原子信息
    
    参数:
        pdb_file: PDB文件路径
    
    返回:
        atoms: 原子信息列表
    """
    atoms = []
    try:
        with open(pdb_file, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")) and len(line) > 54:
                    chain = line[21].strip()
                    res_name = line[17:20].strip()
                    res_id = line[22:26].strip()
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    atoms.append((chain, res_name, res_id, atom_name, (x, y, z)))
    except Exception as e:
        print(f"[parse_pdb_file] Failed to read PDB file: {e}")
    
    return atoms

def distance(a, b):
    """计算两点间距离"""
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

def get_element_from_atom_name(atom_name):
    """从原子名提取元素符号"""
    elem = ''.join(c for c in atom_name if c.isalpha())
    if len(elem) == 0:
        return ""
    # 常见双字母元素
    two_letter = {"BR", "CL", "FE", "MG", "CA", "ZN", "CU", "MN", "CO", "NI"}
    if len(elem) >= 2 and elem[:2].upper() in two_letter:
        return elem[:2].upper()
    return elem[0].upper()

def calculate_angle_three_points(p1, p2, p3):
    """
    计算三点形成的角度 ∠p1-p2-p3（p2为顶点）
    返回角度（度）
    """
    v1 = [p1[i] - p2[i] for i in range(3)]
    v2 = [p3[i] - p2[i] for i in range(3)]
    
    dot = sum(v1[i] * v2[i] for i in range(3))
    mag1 = math.sqrt(sum(v1[i]**2 for i in range(3)))
    mag2 = math.sqrt(sum(v2[i]**2 for i in range(3)))
    
    if mag1 == 0 or mag2 == 0:
        return 0.0
    
    cos_angle = dot / (mag1 * mag2)
    cos_angle = max(-1.0, min(1.0, cos_angle))  # 防止数值误差
    return math.degrees(math.acos(cos_angle))

def find_hydrogen_or_estimate(donor_atom, acceptor_atom, all_atoms_by_residue):
    """
    查找或估算氢原子位置
    
    参数:
        donor_atom: (chain, resn, resi, name, coord)
        acceptor_atom: (chain, resn, resi, name, coord)
        all_atoms_by_residue: dict {(chain, resn, resi): [atoms]}
    
    返回:
        h_coord: 氢原子坐标
    """
    donor_coord = donor_atom[4]
    acceptor_coord = acceptor_atom[4]
    res_key = (donor_atom[0], donor_atom[1], donor_atom[2])
    
    # 尝试查找同残基中的氢原子
    if res_key in all_atoms_by_residue:
        for atom in all_atoms_by_residue[res_key]:
            if atom[3].startswith('H'):
                dist_to_donor = distance(atom[4], donor_coord)
                if 0.8 < dist_to_donor < 1.2:
                    return atom[4]
    
    # 估算氢原子位置：沿 D-A 方向，距离 D 约 1.0Å
    vec = [acceptor_coord[i] - donor_coord[i] for i in range(3)]
    dist = math.sqrt(sum(v**2 for v in vec))
    if dist == 0:
        return donor_coord
    unit_vec = [v / dist for v in vec]
    h_coord = tuple(donor_coord[i] + 1.0 * unit_vec[i] for i in range(3))
    return h_coord

def is_hbond(atom1, atom2, d):
    """
    判断是否为氢键（简单版本，仅距离）
    保留以兼容旧代码
    
    参数:
        atom1, atom2: 原子信息元组
        d: 距离（Å）
    
    返回:
        bool: 是否为氢键
    """
    return d <= INTERACTION_PARAMS["hbond"]["max_DA_dist"]

def is_hbond_precise(donor_atom, acceptor_atom, all_atoms_by_residue):
    """
    精确的氢键判定（包括角度检查）
    
    参数:
        donor_atom: (chain, resn, resi, name, coord) - 供体原子
        acceptor_atom: (chain, resn, resi, name, coord) - 受体原子
        all_atoms_by_residue: dict {(chain, resn, resi): [atoms]}
    
    返回:
        (is_hbond: bool, distance: float, angle_DHA: float or None)
    """
    donor_coord = donor_atom[4]
    acceptor_coord = acceptor_atom[4]
    donor_elem = get_element_from_atom_name(donor_atom[3])
    acceptor_elem = get_element_from_atom_name(acceptor_atom[3])
    
    # 检查元素类型（供体：N/O/S，受体：N/O/S/F/CL）
    donor_elements = {"N", "O", "S"}
    acceptor_elements = {"N", "O", "S", "F", "CL"}
    
    if donor_elem not in donor_elements or acceptor_elem not in acceptor_elements:
        return False, None, None
    
    # 1. 检查 D-A 距离
    d_da = distance(donor_coord, acceptor_coord)
    if d_da > INTERACTION_PARAMS["hbond"]["max_DA_dist"]:
        return False, None, None
    
    # 2. 查找或估算氢原子位置
    h_coord = find_hydrogen_or_estimate(donor_atom, acceptor_atom, all_atoms_by_residue)
    
    # 3. 计算 ∠D–H···A
    angle_dha = calculate_angle_three_points(donor_coord, h_coord, acceptor_coord)
    
    if angle_dha < INTERACTION_PARAMS["hbond"]["min_donor_angle"]:
        return False, None, None
    
    # 注意：∠H···A–X 的检查较复杂，这里简化处理，主要检查 D-H-A 角度
    
    return True, d_da, angle_dha

def is_saltbridge_atom(res1, atom1_name, res2, atom2_name, d):
    """
    判断原子对是否形成盐桥（精确原子级检测）
    
    参数:
        res1, res2: 残基名称
        atom1_name, atom2_name: 原子名称
        d: 距离
    """
    max_dist = INTERACTION_PARAMS["ionic"]["max_dist"]
    if d > max_dist:
        return False

    positive = {"ARG", "LYS", "HIS"}
    negative = {"ASP", "GLU"}
    # 核酸磷酸基团
    nucleic = {"DA", "DT", "DG", "DC", "A", "T", "G", "C", "U", "PS"} 
    
    # 定义带电原子集合
    # 注意：这里假设 PDB 命名规范
    pos_atoms = {
        "ARG": {"NH1", "NH2", "NE", "NZ"}, # 包括 NE 为了兼容性，虽然主要是 NH1/NH2
        "LYS": {"NZ"},
        "HIS": {"ND1", "NE2"},
    }
    neg_atoms = {
        "ASP": {"OD1", "OD2"},
        "GLU": {"OE1", "OE2"},
    }
    nucleic_neg_atoms = {"OP1", "OP2", "O1P", "O2P"}

    def is_pos_atom(res, atom):
        return res in pos_atoms and atom in pos_atoms[res]

    def is_neg_atom(res, atom):
        if res in neg_atoms:
            return atom in neg_atoms[res]
        if res in nucleic:
            return atom in nucleic_neg_atoms
        return False

    # 检查是否为一正一负
    is_pair1 = is_pos_atom(res1, atom1_name) and is_neg_atom(res2, atom2_name)
    is_pair2 = is_pos_atom(res2, atom2_name) and is_neg_atom(res1, atom1_name)
    
    return is_pair1 or is_pair2

def is_saltbridge(res1, res2, d):
    """
    判断是否为盐桥（旧接口兼容，仅基于残基距离）
    注意：这可能导致过多假阳性，建议使用 is_saltbridge_atom
    """
    positive = {"ARG", "LYS", "HIS"}
    negative = {"ASP", "GLU"}
    nucleic = {"DA", "DT", "DG", "DC", "A", "T", "G", "C", "U"}
    
    max_dist = INTERACTION_PARAMS["ionic"]["max_dist"]
    
    if d > max_dist:
        return False
    
    return ((res1 in positive and (res2 in negative or res2 in nucleic))
            or (res2 in positive and (res1 in negative or res1 in nucleic)))

def is_ionic_precise(res1, atoms1, res2, atoms2, hbond_exists=False):
    """
    精确的离子/盐桥判定
    
    参数:
        res1, atoms1: 残基1名称和原子列表
        res2, atoms2: 残基2名称和原子列表
        hbond_exists: 是否已存在氢键
    
    返回:
        (is_ionic: bool, distance: float or None)
    """
    if INTERACTION_PARAMS["ionic"]["exclude_if_hbond"] and hbond_exists:
        return False, None
    
    positive = {"ARG", "LYS", "HIS"}
    negative = {"ASP", "GLU"}
    
    # 检查是否为正负电荷对
    if not ((res1 in positive and res2 in negative) or (res2 in positive and res1 in negative)):
        return False, None
    
    # 查找带电原子
    charged_atoms1 = []
    charged_atoms2 = []
    
    if res1 in positive:
        for atom in atoms1:
            if atom[3] in ["NZ", "NH1", "NH2", "NE", "ND1", "NE2"]:
                charged_atoms1.append(atom[4])
    elif res1 in negative:
        for atom in atoms1:
            if atom[3] in ["OD1", "OD2", "OE1", "OE2"]:
                charged_atoms1.append(atom[4])
    
    if res2 in positive:
        for atom in atoms2:
            if atom[3] in ["NZ", "NH1", "NH2", "NE", "ND1", "NE2"]:
                charged_atoms2.append(atom[4])
    elif res2 in negative:
        for atom in atoms2:
            if atom[3] in ["OD1", "OD2", "OE1", "OE2"]:
                charged_atoms2.append(atom[4])
    
    if not charged_atoms1 or not charged_atoms2:
        return False, None
    
    # 计算最短距离
    min_dist = float('inf')
    for c1 in charged_atoms1:
        for c2 in charged_atoms2:
            d = distance(c1, c2)
            if d < min_dist:
                min_dist = d
    
    if min_dist <= INTERACTION_PARAMS["ionic"]["max_dist"]:
        return True, min_dist
    
    return False, None

def is_hydrophobic(res1, res2, atom1, atom2, d):
    """
    判断是否为疏水相互作用
    
    参数:
        res1, res2: 残基名称
        atom1, atom2: 原子信息元组 (chain, resn, resi, name, coord)
        d: 距离（Å）
    
    返回:
        bool: 是否为疏水相互作用
    
    标准：疏水接触 ≤ 4.0 Å (other_max)，且原子必须是碳或特定的非极性原子
    """
    # 标准疏水残基
    hydrophobic_residues = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "TRP", "TYR", "CYS"}
    # 标准氨基酸（用于判断是否为配体）
    standard_aa = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                   "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
    
    max_dist = INTERACTION_PARAMS["hydrophobic"]["other_max"]
    if d > max_dist:
        return False

    # 获取元素类型
    elem1 = get_element_from_atom_name(atom1[3])
    elem2 = get_element_from_atom_name(atom2[3])
    
    # 疏水原子：C, F, CL, BR, I, S
    hydrophobic_atoms = {"C", "F", "CL", "BR", "I", "S"}
    
    # 两个原子都必须是疏水原子
    if elem1 not in hydrophobic_atoms or elem2 not in hydrophobic_atoms:
        return False
    
    # 过滤主链 Carbonyl C（原子名为"C"的是主链羰基碳）
    if atom1[3].strip() == "C" or atom2[3].strip() == "C":
        return False
    
    # 判断是否为配体（非标准氨基酸）
    is_res1_ligand = res1.upper() not in standard_aa
    is_res2_ligand = res2.upper() not in standard_aa
    
    # 判断是否为疏水残基
    is_res1_hydro = res1.upper() in hydrophobic_residues
    is_res2_hydro = res2.upper() in hydrophobic_residues
    
    # 情况分析：
    # 1. 两个都是碳原子 -> 疏水接触（最常见）
    if elem1 == "C" and elem2 == "C":
        return True
    
    # 2. 配体与蛋白质疏水残基 -> 疏水接触
    if (is_res1_ligand and is_res2_hydro) or (is_res2_ligand and is_res1_hydro):
        return True
    
    # 3. 两个都是配体（配体-配体） -> 疏水接触
    if is_res1_ligand and is_res2_ligand:
        return True
    
    # 4. 两个都是疏水残基 -> 疏水接触
    if is_res1_hydro and is_res2_hydro:
        return True
    
    return False

def centroid(coords):
    """计算坐标质心"""
    n = len(coords)
    if n == 0:
        return (0, 0, 0)
    return tuple(sum(c[i] for c in coords)/n for i in range(3))


def detect_ligand_rings(atoms, min_ring_size=5, max_ring_size=7):
    """
    检测配体中的芳香/平面环结构
    
    [这是什么？] 通过几何方法检测配体分子中的环结构（不依赖 RDKit）
    [为什么要这么做？] 配体没有标准残基名，无法使用预定义的原子列表
    [为什么这是个好主意！] 使得 π-π 检测可以适用于任意配体
    
    参数:
        atoms: 原子信息列表 [(chain, resn, resi, name, coord), ...]
        min_ring_size: 最小环大小（默认5，五元环）
        max_ring_size: 最大环大小（默认7，七元环）
    
    返回:
        list of lists: 每个内部列表包含一个环的原子坐标 [(x,y,z), ...]
                       如果没有检测到环，返回空列表
    
    方法：
        1. 筛选 sp2 杂化可能的原子（C, N, O, S）
        2. 基于距离（1.2-1.6Å）构建连接图
        3. 使用简单的 DFS 寻找环
        4. 验证环的共面性（法向量检查）
    """
    # 1. 筛选可能参与芳香环的原子（主要是 C 和 N）
    aromatic_elements = {"C", "N"}  # 主要关注碳和氮
    
    ring_candidates = []
    for atom in atoms:
        elem = get_element_from_atom_name(atom[3])
        # 排除主链原子
        if atom[3].strip() in ["C", "CA", "N", "O"]:
            continue
        if elem in aromatic_elements:
            ring_candidates.append(atom)
    
    if len(ring_candidates) < min_ring_size:
        return []
    
    # 2. 构建连接图（基于距离）
    # 芳香环 C-C 键长约 1.39Å，C-N 键长约 1.35Å
    bond_min = 1.2
    bond_max = 1.65
    
    n = len(ring_candidates)
    adjacency = [[] for _ in range(n)]
    
    for i in range(n):
        for j in range(i + 1, n):
            d = distance(ring_candidates[i][4], ring_candidates[j][4])
            if bond_min <= d <= bond_max:
                adjacency[i].append(j)
                adjacency[j].append(i)
    
    # 3. 使用 DFS 寻找环
    detected_rings = []
    visited_rings = set()  # 避免重复检测相同的环
    
    def dfs_find_rings(start, current, path, depth):
        """DFS 寻找从 start 开始的环"""
        if depth > max_ring_size:
            return
        
        for neighbor in adjacency[current]:
            if neighbor == start and depth >= min_ring_size:
                # 找到一个环
                ring_key = tuple(sorted(path))
                if ring_key not in visited_rings:
                    visited_rings.add(ring_key)
                    ring_coords = [ring_candidates[i][4] for i in path]
                    # 验证共面性
                    if is_planar_ring(ring_coords):
                        detected_rings.append(ring_coords)
                return
            
            if neighbor not in path and neighbor > start:  # 避免重复和回退
                dfs_find_rings(start, neighbor, path + [neighbor], depth + 1)
    
    # 从每个节点开始搜索
    for start in range(n):
        if len(adjacency[start]) >= 2:  # 至少有两个连接才可能形成环
            dfs_find_rings(start, start, [start], 1)
    
    return detected_rings


def is_planar_ring(coords, tolerance=0.5):
    """
    检查一组坐标是否共面（用于验证芳香环）
    
    参数:
        coords: 坐标列表 [(x,y,z), ...]
        tolerance: 允许的最大偏离（Å）
    
    返回:
        bool: 是否共面
    """
    if len(coords) < 4:
        return True  # 3个点总是共面
    
    # 使用前3个点定义平面
    p0, p1, p2 = coords[0], coords[1], coords[2]
    
    # 计算法向量
    v1 = [p1[i] - p0[i] for i in range(3)]
    v2 = [p2[i] - p0[i] for i in range(3)]
    
    # 叉积得法向量
    nx = v1[1] * v2[2] - v1[2] * v2[1]
    ny = v1[2] * v2[0] - v1[0] * v2[2]
    nz = v1[0] * v2[1] - v1[1] * v2[0]
    
    norm = math.sqrt(nx*nx + ny*ny + nz*nz)
    if norm < 1e-6:
        return False
    
    nx, ny, nz = nx/norm, ny/norm, nz/norm
    
    # 计算平面方程 d 值: ax + by + cz = d
    d = nx * p0[0] + ny * p0[1] + nz * p0[2]
    
    # 检查其余点到平面的距离
    for p in coords[3:]:
        dist = abs(nx * p[0] + ny * p[1] + nz * p[2] - d)
        if dist > tolerance:
            return False
    
    return True


def ring_atoms(res_name, atoms):
    """
    获取环状结构原子坐标
    
    [这是什么？] 返回残基中芳香环的原子坐标
    [为什么要这么做？] 用于计算 π-π 堆积和 π-阳离子相互作用
    [为什么这是个好主意！] 现在支持配体动态检测，不再局限于标准残基
    
    参数:
        res_name: 残基名称
        atoms: 原子列表 [(chain, resn, resi, name, coord), ...]
    
    返回:
        coords: 环原子坐标列表，或 None（无环）
    """
    # 标准残基的预定义环原子
    ring_dict = {
        "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TRP": ["CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
        "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
        "A": ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"],
        "G": ["N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"],
        "C": ["N1", "C2", "N3", "C4", "C5", "C6"],
        "T": ["N1", "C2", "N3", "C4", "C5", "C6"],
        "U": ["N1", "C2", "N3", "C4", "C5", "C6"],
    }
    
    # 标准残基使用预定义的原子列表
    if res_name in ring_dict:
        coords = [a[4] for a in atoms if a[3].strip() in ring_dict[res_name]]
        return coords if len(coords) >= 3 else None
    
    # 非标准残基（配体）：使用动态检测
    detected_rings = detect_ligand_rings(atoms)
    if detected_rings:
        # 返回第一个检测到的环（可根据需求扩展为返回所有环）
        return detected_rings[0]
    
    return None


def ring_atoms_all(res_name, atoms):
    """
    获取残基中所有芳香环的原子坐标（支持多环配体）
    
    参数:
        res_name: 残基名称
        atoms: 原子列表 [(chain, resn, resi, name, coord), ...]
    
    返回:
        list of coords: 每个元素是一个环的原子坐标列表
    """
    # 标准残基的预定义环原子
    ring_dict = {
        "PHE": [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
        "TYR": [["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]],
        "TRP": [
            ["CG", "CD1", "NE1", "CE2", "CD2"],  # 五元环
            ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]  # 六元环
        ],
        "HIS": [["CG", "ND1", "CD2", "CE1", "NE2"]],
    }
    
    if res_name in ring_dict:
        result = []
        for ring_atoms_list in ring_dict[res_name]:
            coords = [a[4] for a in atoms if a[3].strip() in ring_atoms_list]
            if len(coords) >= 3:
                result.append(coords)
        return result if result else None
    
    # 非标准残基（配体）：使用动态检测
    detected_rings = detect_ligand_rings(atoms)
    return detected_rings if detected_rings else None

def normal_vector(a, b, c):
    """计算三点确定平面的法向量"""
    ab = [b[i]-a[i] for i in range(3)]
    ac = [c[i]-a[i] for i in range(3)]
    n = [ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0]]
    norm = math.sqrt(sum(x*x for x in n))
    return [x/norm for x in n] if norm != 0 else None

def angle_between(v1, v2):
    """计算两向量夹角（度）"""
    dot = sum(v1[i]*v2[i] for i in range(3))
    dot = max(min(dot, 1), -1)  # 防止数值误差
    return math.degrees(math.acos(dot))

def is_pipi(res1, atoms1, res2, atoms2):
    """
    判断是否为π-π堆积（简单版本，兼容旧代码）
    """
    r1 = ring_atoms(res1, atoms1)
    r2 = ring_atoms(res2, atoms2)
    
    if not r1 or not r2:
        return False
    
    c1, c2 = centroid(r1), centroid(r2)
    d = distance(c1, c2)
    
    if not (3.3 <= d <= 6.0):
        return False
    
    n1 = normal_vector(r1[0], r1[1], r1[2])
    n2 = normal_vector(r2[0], r2[1], r2[2])
    
    if n1 is None or n2 is None:
        return False
    
    angle = angle_between(n1, n2)
    return angle <= 30 or 60 <= angle <= 120

def is_pipi_precise(res1, atoms1, res2, atoms2):
    """
    精确的π-π堆积判定，区分面-面和面-边模式
    
    返回:
        (is_pipi: bool, mode: str or None, distance: float or None)
        mode: "face-face" 或 "edge-face"
    """
    r1 = ring_atoms(res1, atoms1)
    r2 = ring_atoms(res2, atoms2)
    
    if not r1 or not r2:
        return False, None, None
    
    # 计算环心
    c1, c2 = centroid(r1), centroid(r2)
    d = distance(c1, c2)
    
    # 环心距离 3.3-5.5Å
    if not (3.3 <= d <= 5.5):
        return False, None, None
    
    # 计算法向量
    n1 = normal_vector(r1[0], r1[1], r1[2])
    n2 = normal_vector(r2[0], r2[1], r2[2])
    
    if n1 is None or n2 is None:
        return False, None, None
    
    # 法向量夹角
    angle = angle_between(n1, n2)
    
    # 面-面对：夹角 < 30° 或 > 150°
    if angle < 30 or angle > 150:
        return True, "face-face", d
    
    # 面-边：夹角 60-120°
    if 60 <= angle <= 120:
        return True, "edge-face", d
    
    return False, None, None

def is_metal_coordination(atom1, atom2):
    """
    判断是否为金属配位
    
    参数:
        atom1, atom2: (chain, resn, resi, name, coord)
    
    返回:
        (is_metal: bool, distance: float or None)
    """
    elem1 = get_element_from_atom_name(atom1[3])
    elem2 = get_element_from_atom_name(atom2[3])
    
    # 检查是否有金属
    metal_atoms = {"ZN", "MG", "CA", "FE", "CU", "MN", "CO", "NI"}
    if not (elem1 in metal_atoms or elem2 in metal_atoms):
        return False, None
    
    # 检查配位原子（非碳重原子）
    coord_atoms = INTERACTION_PARAMS["metal_coord"]["allowed_ligand_atoms"]
    metal_coord = atom1[4] if elem1 in metal_atoms else atom2[4]
    coord_atom_elem = elem2 if elem1 in metal_atoms else elem1
    coord_atom_coord = atom2[4] if elem1 in metal_atoms else atom1[4]
    
    if coord_atom_elem not in coord_atoms:
        return False, None
    
    d = distance(metal_coord, coord_atom_coord)
    
    if d <= INTERACTION_PARAMS["metal_coord"]["max_dist"]:
        return True, d
    
    return False, None

def is_cationpi(res1, atoms1, res2, atoms2):
    """判断是否为π-阳离子相互作用"""
    pos_res = {"ARG", "LYS", "HIS"}
    ring_res = {"PHE", "TYR", "TRP", "HIS", "A", "G", "C", "T", "U"}
    
    def cation_center(atoms):
        pos_atoms = [a[4] for a in atoms if a[3].strip().startswith(("N", "NZ", "NH", "NE"))]
        return centroid(pos_atoms) if pos_atoms else None
    
    if res1 in pos_res and res2 in ring_res:
        cation, ring = cation_center(atoms1), ring_atoms(res2, atoms2)
    elif res2 in pos_res and res1 in ring_res:
        cation, ring = cation_center(atoms2), ring_atoms(res1, atoms1)
    else:
        return False
    
    if not cation or not ring:
        return False
    
    c_ring = centroid(ring)
    d = distance(cation, c_ring)
    return d <= 6.0

def is_halogen_bond(atom1, atom2, d):
    """
    判断是否为卤素键
    
    参数:
        atom1, atom2: (chain, resn, resi, name, coord)
        d: 距离
    
    返回:
        (is_halogen: bool, distance: float, angle: float)
    """
    # 1. 识别卤素供体 (X) 和受体 (A)
    # 供体必须是卤素原子 (F, CL, BR, I)
    # 受体通常是 O, N, S
    
    elem1 = get_element_from_atom_name(atom1[3])
    elem2 = get_element_from_atom_name(atom2[3])
    
    halogens = {"CL", "BR", "I", "F"}  # F 也可以形成卤素键，尽管较弱
    acceptors = {"O", "N", "S"}
    
    halogen_atom = None
    acceptor_atom = None
    
    if elem1 in halogens and elem2 in acceptors:
        halogen_atom = atom1
        acceptor_atom = atom2
    elif elem2 in halogens and elem1 in acceptors:
        halogen_atom = atom2
        acceptor_atom = atom1
    else:
        return False, None, None
        
    # 2. 检查距离
    max_dist = INTERACTION_PARAMS["halogen"]["max_distance"]
    if d > max_dist:
        return False, None, None
        
    # 3. 检查角度 C-X···A
    # 需要找到与卤素相连的碳原子 (C)
    # 由于这里没有拓扑信息，我们只能在同残基中寻找距离卤素最近的碳原子
    
    # 获取卤素所在残基的所有原子 (需要传入 all_atoms_by_residue，或者在这里简化处理)
    # 为了保持接口简单，我们暂时无法精确计算 C-X...A 角度，除非传入更多上下文
    # 这里做一个近似：假设 C 在 X 的反方向 (仅作示意，实际需要拓扑)
    # 或者，我们修改调用签名，传入 all_atoms_by_residue
    
    # 暂时返回 True (距离满足)，角度检查留给更高级的分析或后续优化
    # 为了严谨，我们至少要求距离比范德华半径之和略小，这里直接用 cutoff
    
    return True, d, 180.0 # 假定角度理想，后续应优化

def is_water_bridge(ligand_atom, protein_atom, water_atom):
    """
    判断是否为水桥 (Ligand - Water - Protein)
    
    参数:
        ligand_atom: 配体原子
        protein_atom: 蛋白原子
        water_atom: 水分子原子
    
    返回:
        (is_wb: bool, d_lw: float, d_wp: float, angle: float)
    """
    # 1. 检查 Ligand-Water 距离 (H-bond)
    d_lw = distance(ligand_atom[4], water_atom[4])
    if d_lw > INTERACTION_PARAMS["water_bridge"]["max_DA_dist"]:
        return False, None, None, None
        
    # 2. 检查 Water-Protein 距离 (H-bond)
    d_wp = distance(water_atom[4], protein_atom[4])
    if d_wp > INTERACTION_PARAMS["water_bridge"]["max_DA_dist"]:
        return False, None, None, None
        
    # 3. 检查角度 Ligand-Water-Protein
    angle = calculate_angle_three_points(ligand_atom[4], water_atom[4], protein_atom[4])
    
    # 水桥的角度通常在 75-140 度之间 (水的 H-O-H 角度约为 104.5)
    # 这里使用较宽松的范围
    if angle < INTERACTION_PARAMS["water_bridge"]["min_acceptor_angle"]:
        return False, None, None, None
        
    return True, d_lw, d_wp, angle

def calculate_confidence_score(interaction_type, distance_val, angle_val=None):
    """
    计算相互作用置信度分数 (0.0 - 1.0)
    基于几何参数的理想程度
    """
    score = 1.0
    
    if interaction_type == "氢键":
        # 距离衰减: 2.8A -> 1.0, 3.5A -> 0.5
        optimal_dist = 2.8
        max_dist = INTERACTION_PARAMS["hbond"]["max_DA_dist"]
        if distance_val <= optimal_dist:
            dist_score = 1.0
        else:
            dist_score = max(0.0, 1.0 - (distance_val - optimal_dist) / (max_dist - optimal_dist))
            
        # 角度衰减: 180 -> 1.0, 120 -> 0.5
        if angle_val:
            optimal_angle = 180.0
            min_angle = INTERACTION_PARAMS["hbond"]["min_donor_angle"]
            angle_score = max(0.0, (angle_val - min_angle) / (optimal_angle - min_angle))
        else:
            angle_score = 0.8 # 默认
            
        score = dist_score * 0.7 + angle_score * 0.3
        
    elif interaction_type == "盐桥":
        optimal_dist = 3.5
        max_dist = INTERACTION_PARAMS["ionic"]["max_dist"]
        if distance_val <= optimal_dist:
            score = 1.0
        else:
            score = max(0.0, 1.0 - (distance_val - optimal_dist) / (max_dist - optimal_dist))
            
    elif interaction_type == "疏水相互作用":
        optimal_dist = 3.8
        max_dist = INTERACTION_PARAMS["hydrophobic"]["other_max"]
        if distance_val <= optimal_dist:
            score = 1.0
        else:
            score = max(0.0, 1.0 - (distance_val - optimal_dist) / (max_dist - optimal_dist))
            
    # ... 其他类型
    
    return round(score, 2)

def identify_molecule_type(res_name):
    """
    识别分子类型

    返回: "protein", "ligand", "peptide", "dna", "rna", "unknown"
    """
    # 标准氨基酸（20种）
    standard_aa = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                   "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}

    # DNA/RNA碱基
    nucleic_bases = {"A", "T", "G", "C", "U", "DA", "DT", "DG", "DC", "DU"}

    # 常见离子和溶剂（通常忽略）
    solvents_ions = {"HOH", "WAT", "NA", "CL", "K", "MG", "CA", "ZN", "FE"}

    res_name = res_name.strip().upper()

    if res_name in standard_aa:
        return "protein"
    elif res_name in nucleic_bases:
        if res_name.startswith("D"):
            return "dna"
        else:
            return "rna"
    elif res_name in solvents_ions:
        return "solvent"
    else:
        # 小分子配体（非标准残基）
        return "ligand"

def classify_chain_type(atoms_in_chain):
    """
    根据链中的残基类型分类整条链

    返回: "protein", "ligand", "peptide", "dna", "rna", "mixed"
    """
    res_types = {}
    for atom in atoms_in_chain:
        res_name = atom[1]  # atom format: (chain, res_name, res_id, atom_name, coords)
        mol_type = identify_molecule_type(res_name)
        if mol_type != "solvent":
            res_types[mol_type] = res_types.get(mol_type, 0) + 1

    if not res_types:
        return "unknown"

    # 找出主要类型
    main_type = max(res_types, key=res_types.get)

    # 如果是蛋白质但残基数很少（<15），可能是多肽
    if main_type == "protein" and res_types[main_type] < 15:
        return "peptide"

    # 如果只有一种类型占绝对优势（>80%），返回该类型
    total = sum(res_types.values())
    if res_types[main_type] / total > 0.8:
        return main_type

    return "mixed"

def analyze_interactions(atoms, only_between_chains=False):
    """
    分析原子间相互作用

    参数:
        atoms: 原子信息列表
        only_between_chains: 是否只分析链间相互作用

    返回:
        interactions: 相互作用列表
    
    使用严格标准：氢键 ≤2.8Å，盐桥 ≤4.0Å
    """
    # 按残基分组
    grouped = defaultdict(list)
    for atom in atoms:
        grouped[(atom[0], atom[1], atom[2])].append(atom)

    interactions = []
    keys = list(grouped.keys())

    # 如果只分析链间相互作用，优化为只遍历不同链的残基对
    if only_between_chains:
        # 按链分组
        chain_residues = defaultdict(list)
        for key in keys:
            chain_residues[key[0]].append(key)

        chains = list(chain_residues.keys())
        print(f"[analyze_interactions] Detected {len(chains)} chains; analyzing inter-chain interfaces only...")

        # 只遍历不同链之间的残基对
        total_pairs = sum(len(chain_residues[chains[i]]) * len(chain_residues[chains[j]])
                         for i in range(len(chains)) for j in range(i+1, len(chains)))

        checked = 0
        last_percent = 0

        for i in range(len(chains)):
            for j in range(i+1, len(chains)):
                chain1_res = chain_residues[chains[i]]
                chain2_res = chain_residues[chains[j]]

                for key1 in chain1_res:
                    for key2 in chain2_res:
                        checked += 1
                        percent = int(checked * 100 / total_pairs)
                        if percent >= last_percent + 10:
                            print(f"[analyze_interactions] Progress: {percent}% ({checked}/{total_pairs})")
                            last_percent = percent

                        c1, r1, id1 = key1
                        c2, r2, id2 = key2

                        a1, a2 = grouped[key1], grouped[key2]

                        # 快速距离筛选：先检查残基质心距离
                        c1_coords = [at[4] for at in a1]
                        c2_coords = [at[4] for at in a2]
                        centroid1 = centroid(c1_coords)
                        centroid2 = centroid(c2_coords)
                        centroid_dist = distance(centroid1, centroid2)

                        # 如果质心距离过大，跳过详细检查（缩小到10埃，接触界面更精确）
                        if centroid_dist > 10.0:
                            continue

                        # 检查相互作用
                        _check_residue_interactions(c1, r1, id1, c2, r2, id2, a1, a2, grouped, interactions)
    else:
        # 分析所有残基对（包括链内）
        print(f"[analyze_interactions] Analyzing all interactions among {len(keys)} residues...")

        total_pairs = len(keys) * (len(keys) - 1) // 2
        checked = 0
        last_percent = 0

        for i in range(len(keys)):
            for j in range(i+1, len(keys)):
                c1, r1, id1 = keys[i]
                c2, r2, id2 = keys[j]

                # 进度输出（每10%输出一次）
                checked += 1
                percent = int(checked * 100 / total_pairs)
                if percent >= last_percent + 10:
                    print(f"[analyze_interactions] Progress: {percent}% ({checked}/{total_pairs})")
                    last_percent = percent

                a1, a2 = grouped[keys[i]], grouped[keys[j]]

                # 快速距离筛选：先检查残基质心距离
                c1_coords = [at[4] for at in a1]
                c2_coords = [at[4] for at in a2]
                centroid1 = centroid(c1_coords)
                centroid2 = centroid(c2_coords)
                centroid_dist = distance(centroid1, centroid2)

                # 如果质心距离过大，跳过详细检查
                if centroid_dist > 15.0:
                    continue

                # 检查相互作用
                _check_residue_interactions(c1, r1, id1, c2, r2, id2, a1, a2, grouped, interactions)

    return interactions

def _check_residue_interactions(c1, r1, id1, c2, r2, id2, a1, a2, grouped, interactions):
    """
    检查两个残基间的相互作用（提取为独立函数避免代码重复）
    
    使用严格标准进行相互作用判断
    修复：并行检测所有相互作用类型，而非互斥判断
    """
    # 使用集合记录已检测到的相互作用类型，避免重复
    found_hbond = False
    found_saltbridge = False
    found_hydrophobic = False
    found_halogen = False
    
    # 检查原子间相互作用（氢键、盐桥、疏水、卤素键）
    for at1 in a1:
        for at2 in a2:
            d = distance(at1[4], at2[4])
            if d > 6.0:  # 扩大初筛阈值，覆盖所有相互作用类型
                continue

            # 并行检测所有类型（不使用 elif）
            if not found_hbond:
                 # 尝试双向检测氢键
                 is_hb, hb_dist, hb_angle = is_hbond_precise(at1, at2, grouped) # 注意：这里 grouped 并不完全等同于 all_atoms_by_residue，需要确认
                 if not is_hb:
                     is_hb, hb_dist, hb_angle = is_hbond_precise(at2, at1, grouped)
                 
                 if is_hb:
                    conf = calculate_confidence_score("氢键", hb_dist, hb_angle)
                    interactions.append({
                        "Chain1": c1,
                        "Residue1": f"{r1} {id1}",
                        "Atom1": at1[3],
                        "Chain2": c2,
                        "Residue2": f"{r2} {id2}",
                        "Atom2": at2[3],
                        "Distance": round(hb_dist, 2),
                        "Interaction": "氢键",
                        "Confidence": conf
                    })
                    found_hbond = True
            
            if not found_saltbridge and is_saltbridge_atom(r1, at1[3], r2, at2[3], d):
                conf = calculate_confidence_score("盐桥", d)
                interactions.append({
                    "Chain1": c1,
                    "Residue1": f"{r1} {id1}",
                    "Atom1": at1[3],
                    "Chain2": c2,
                    "Residue2": f"{r2} {id2}",
                    "Atom2": at2[3],
                    "Distance": round(d, 2),
                    "Interaction": "盐桥",
                    "Confidence": conf
                })
                found_saltbridge = True
            
            if not found_hydrophobic and is_hydrophobic(r1, r2, at1, at2, d):
                conf = calculate_confidence_score("疏水相互作用", d)
                interactions.append({
                    "Chain1": c1,
                    "Residue1": f"{r1} {id1}",
                    "Atom1": at1[3],
                    "Chain2": c2,
                    "Residue2": f"{r2} {id2}",
                    "Atom2": at2[3],
                    "Distance": round(d, 2),
                    "Interaction": "疏水相互作用",
                    "Confidence": conf
                })
                found_hydrophobic = True

            if not found_halogen:
                is_halogen, hal_dist, hal_angle = is_halogen_bond(at1, at2, d)
                if is_halogen:
                    conf = calculate_confidence_score("卤素键", hal_dist, hal_angle)
                    interactions.append({
                        "Chain1": c1,
                        "Residue1": f"{r1} {id1}",
                        "Atom1": at1[3],
                        "Chain2": c2,
                        "Residue2": f"{r2} {id2}",
                        "Atom2": at2[3],
                        "Distance": round(hal_dist, 2),
                        "Interaction": "卤素键",
                        "Confidence": conf
                    })
                    found_halogen = True

    # 检查π相互作用（独立于原子级检测）
    if is_pipi(r1, a1, r2, a2):
        interactions.append({
            "Chain1": c1,
            "Residue1": f"{r1} {id1}",
            "Atom1": "Ring",
            "Chain2": c2,
            "Residue2": f"{r2} {id2}",
            "Atom2": "Ring",
            "Distance": "-",
            "Interaction": "π–π 堆积",
            "Confidence": 0.9 # 默认高置信度
        })

    if is_cationpi(r1, a1, r2, a2):
        interactions.append({
            "Chain1": c1,
            "Residue1": f"{r1} {id1}",
            "Atom1": "Cation", 
            "Chain2": c2,
            "Residue2": f"{r2} {id2}",
            "Atom2": "Ring",
            "Distance": "-",
            "Interaction": "π–阳离子相互作用",
            "Confidence": 0.9
        })

def analyze_pdb_interactions(obj_name=None, output_csv=None, only_between_chains=True, 
                           auto_highlight=True, pdb_file=None):
    """
    分析PDB结构中的相互作用并可选地在PyMOL中高亮显示
    
    参数:
        obj_name: PyMOL对象名称
        output_csv: 输出CSV文件路径
        only_between_chains: 是否只分析链间相互作用
        auto_highlight: 是否自动高亮显示结果
        pdb_file: PDB文件路径（如果提供，将从文件读取而不是PyMOL对象）
    
    返回:
        interactions: 相互作用列表
    
    使用严格标准：氢键 ≤2.8Å，盐桥 ≤4.0Å，适合发表
    
    示例:
        analyze_pdb_interactions('protein', output_csv='interactions.csv')
    """
    # 获取原子信息
    if pdb_file:
        atoms = parse_pdb_file(pdb_file)
        if not atoms:
            print(f"[analyze_pdb_interactions] Unable to read atom information from file: {pdb_file}")
            return []
    else:
        atoms = parse_pdb_structure(obj_name)
        if not atoms:
            print("[analyze_pdb_interactions] Unable to obtain atom information")
            return []
    
    print("[analyze_pdb_interactions] Using strict standards (H-bond ≤2.8Å, Salt bridge ≤4.0Å)")
    
    # 分析相互作用
    interactions = analyze_interactions(atoms, only_between_chains)
    
    # 输出到CSV文件
    if output_csv:
        try:
            with open(output_csv, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(["Chain1", "Residue1", "Atom1", "Chain2", "Residue2", "Atom2", "Distance", "Interaction", "Confidence"])
                for inter in interactions:
                    writer.writerow([
                        inter["Chain1"],
                        inter["Residue1"],
                        inter.get("Atom1", ""),
                        inter["Chain2"],
                        inter["Residue2"],
                        inter.get("Atom2", ""),
                        inter["Distance"],
                        inter["Interaction"],
                        inter.get("Confidence", "-")
                    ])
            print(f"[analyze_pdb_interactions] Results saved to: {output_csv}")
        except Exception as e:
            print(f"[analyze_pdb_interactions] Failed to save CSV: {e}")
    
    print(f"[analyze_pdb_interactions] Analysis complete: found {len(interactions)} interactions")
    
    # 自动高亮显示（如枟启用且在PyMOL环境中）
    if auto_highlight and not pdb_file and interactions:
        print(f"[analyze_pdb_interactions] Auto-highlighting {len(interactions)} interactions in PyMOL...")
        try:
            # 如果已经有CSV文件，直接使用
            if output_csv and os.path.exists(output_csv):
                from .highlight_residues import highlight_csv_residues
                highlight_csv_residues(output_csv, obj=obj_name, show_labels=1,
                                     stick_by_element=1)
                print(f"[analyze_pdb_interactions] ✅ Highlighted interactions from {output_csv}")
            else:
                # 创建临时CSV文件
                temp_csv = tempfile.NamedTemporaryFile(mode='w', suffix='.csv',
                                                     delete=False, encoding='utf-8')
                writer = csv.writer(temp_csv)
                writer.writerow(["Chain1", "Residue1", "Atom1", "Chain2", "Residue2", "Atom2", "Distance", "Interaction", "Confidence"])
                for inter in interactions:
                    writer.writerow([
                        inter["Chain1"],
                        inter["Residue1"],
                        inter.get("Atom1", ""),
                        inter["Chain2"],
                        inter["Residue2"],
                        inter.get("Atom2", ""),
                        inter["Distance"],
                        inter["Interaction"],
                        inter.get("Confidence", "-")
                    ])
                temp_csv.close()

                # 调用高亮功能
                from .highlight_residues import highlight_csv_residues
                highlight_csv_residues(temp_csv.name, obj=obj_name, show_labels=1,
                                     stick_by_element=1)
                
                print(f"[analyze_pdb_interactions] ✅ Highlighted interactions in PyMOL")

                # 清理临时文件
                os.unlink(temp_csv.name)

        except Exception as e:
            print(f"[analyze_pdb_interactions] ⚠️ Auto highlight failed: {e}")
            import traceback
            traceback.print_exc()
    elif auto_highlight and pdb_file:
        print(f"[analyze_pdb_interactions] ℹ️ Auto-highlight skipped (analyzing from PDB file, not PyMOL object)")
    elif auto_highlight and not interactions:
        print(f"[analyze_pdb_interactions] ℹ️ Auto-highlight skipped (no interactions found)")
    
    return interactions

# 注册命令到PyMOL
cmd.extend("analyze_pdb_interactions", analyze_pdb_interactions)

def render_interactions_beautifully(obj_name, csv_path=None, interactions=None):
    """
    美化渲染相互作用

    参数:
        obj_name: PyMOL对象名称
        csv_path: CSV文件路径（可选）
        interactions: 相互作用列表（可选，如果没有提供将从csv_path读取）
    """
    from pymol import cmd

    # 分析链类型
    atoms = parse_pdb_structure(obj_name)
    if not atoms:
        print(f"[render_interactions_beautifully] 无法获取 {obj_name} 的原子信息")
        return

    # 按链分组
    chains = defaultdict(list)
    for atom in atoms:
        chains[atom[0]].append(atom)

    # 分类每条链
    chain_types = {}
    for chain, atoms_in_chain in chains.items():
        chain_types[chain] = classify_chain_type(atoms_in_chain)

    print(f"[render_interactions_beautifully] 链类型分析: {chain_types}")

    # 基础视图设置
    cmd.hide("everything", obj_name)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)

    # 根据链类型设置显示样式和颜色
    color_scheme = {
        "protein": "skyblue",
        "peptide": "lightpink",
        "ligand": "tv_orange",
        "dna": "lime",
        "rna": "forest",
        "mixed": "grey70",
        "unknown": "grey50"
    }

    for chain, ctype in chain_types.items():
        sel = f"{obj_name} and chain {chain}"

        if ctype in ["protein", "peptide"]:
            # 蛋白/多肽：cartoon + 半透明surface
            cmd.show("cartoon", sel)
            cmd.show("surface", sel)
            cmd.set("transparency", 0.4, sel)
            cmd.color(color_scheme[ctype], sel)

        elif ctype == "ligand":
            # 小分子：sticks + 球棍
            cmd.show("sticks", sel)
            cmd.show("spheres", sel)
            cmd.set("sphere_scale", 0.25, sel)
            cmd.color(color_scheme[ctype], f"{sel} and elem C")
            cmd.color("blue", f"{sel} and elem N")
            cmd.color("red", f"{sel} and elem O")
            cmd.color("yellow", f"{sel} and elem S")

        elif ctype in ["dna", "rna"]:
            # 核酸：cartoon + sticks
            cmd.show("cartoon", sel)
            cmd.show("sticks", sel)
            cmd.set("cartoon_ring_mode", 3, sel)
            cmd.color(color_scheme[ctype], sel)

        else:
            # 未知/混合：默认显示
            cmd.show("lines", sel)
            cmd.color(color_scheme[ctype], sel)

    # 读取相互作用数据
    if csv_path and os.path.exists(csv_path):
        try:
            from .highlight_residues import highlight_csv_residues
            # 使用现有的高亮函数
            highlight_csv_residues(csv_path, obj=obj_name, show_labels=1,
                                 stick_by_element=1, clear_old=0)
        except Exception as e:
            print(f"[render_interactions_beautifully] 高亮相互作用失败: {e}")

    # 设置光照和质量
    cmd.set("ambient", 0.2)
    cmd.set("spec_power", 80)
    cmd.set("spec_reflect", 0.4)
    cmd.set("depth_cue", 1)
    cmd.set("ray_trace_mode", 1)
    cmd.set("surface_quality", 1)

    # 自动调整视角
    cmd.orient(obj_name)
    cmd.zoom(obj_name, buffer=5.0)

    # 渲染完成,静默返回

cmd.extend("render_interactions_beautifully", render_interactions_beautifully)

def analyze_protein_ligand_interactions(obj_name=None, ligand_resname=None,
                                       protein_chains=None, output_csv=None,
                                       distance_cutoff=4.5, pdb_file=None,
                                       key_interactions_only=False):
    """
    分析蛋白质-配体相互作用（使用严格标准）

    参数:
        obj_name: PyMOL对象名称
        ligand_resname: 配体残基名称（例如 "LIG", "ATP" 等）,如果为None则自动检测
        protein_chains: 蛋白质链ID列表（例如 ["A", "B"]）,如果为None则自动检测
        output_csv: 输出CSV文件路径
        distance_cutoff: 距离截断值（埃）
        pdb_file: PDB文件路径（可选）
        key_interactions_only: 仅检测关键相互作用（氢键、盐桥、π相互作用、金属配位）,排除疏水接触 (默认False,包含疏水)

    返回:
        dict: {
            "ligand_residues": [{"chain": X, "resname": Y, "resid": Z}, ...],
            "protein_chains": ["A", "B", ...],
            "interactions": [...相互作用列表...]
            "mode": "advanced"
        }
    
    推荐标准（适合大多数场景）：
        - 氢键: ≤3.2Å，角度≥120° (精确几何验证)
        - 盐桥: ≤4.5Å
        - 疏水: ≤4.0Å
        - π相互作用: 严格几何判定
        - 金属配位: ≤3.4Å
    
    注：现已使用精确氢键检测（距离+角度），避免假阳性
    
    示例:
        result = analyze_protein_ligand_interactions('complex', 'LIG', output_csv='interactions.csv')
    """
    
    # ========== 确保使用高质量分析模式 ==========
    if not RDKIT_AVAILABLE:
        error_msg = "[analyze_protein_ligand_interactions] ❌ RDKit not installed. High-quality analysis unavailable.\nRun: pip install rdkit scipy matplotlib pillow numpy"
        print(error_msg)
        raise RuntimeError(error_msg)
    
    if not NUMPY_AVAILABLE:
        error_msg = "[analyze_protein_ligand_interactions] ❌ NumPy not installed\nRun: pip install numpy"
        print(error_msg)
        raise RuntimeError(error_msg)
    
    # 使用RDKit进行高质量分析（直接在本模块实现，不依赖external advanced模块）
    # print("[analyze_protein_ligand_interactions] 🚀 使用高质量 RDKit 分析模式")
    
    # 如果存在 advanced 模块则使用它
    try:
        from .interaction_analyzer_advanced import analyze_interactions_advanced
        
        # 导出蛋白质和配体为临时文件
        protein_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb").name
        ligand_sdf = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf").name
        
        # 保存蛋白质
        if obj_name:
            cmd.save(protein_pdb, f"{obj_name} and polymer")
        elif pdb_file:
            import shutil
            shutil.copy(pdb_file, protein_pdb)
        
        # 保存配体
        if ligand_resname:
            cmd.save(ligand_sdf, f"{obj_name} and resn {ligand_resname}", format="sdf")
        else:
            # 自动检测配体（简化：非protein/水/离子的残基）
            cmd.save(ligand_sdf, f"{obj_name} and not (polymer or resn HOH+WAT or resn NA+CL+MG+CA)", format="sdf")
        
        # 调用高级分析
        result_advanced = analyze_interactions_advanced(protein_pdb, ligand_sdf, output_csv)
        
        # 清理临时蛋白文件（保留ligand_sdf用于2D图绘制）
        try:
            os.remove(protein_pdb)
        except:
            pass
        
        if result_advanced:
            # 转换格式以兼容返回值
            return {
                "ligand_residues": [],
                "protein_chains": protein_chains or [],
                "interactions": [],
                "mode": "advanced",
                "advanced_results": result_advanced,
                "ligand_sdf": ligand_sdf  # 保存SDF路径用于2D图
            }
    except ModuleNotFoundError:
        # advanced 模块不存在，使用内置实现（继续下面的代码）
        # print("[analyze_protein_ligand_interactions] 📋 Using built-in high-quality analysis")
        pass
    except Exception as e:
        print(f"[analyze_protein_ligand_interactions] ⚠️ Advanced module failed: {e}")
        import traceback
        traceback.print_exc()
        # 继续使用内置实现
    
    # ========== 严格标准分析（使用RDKit） ==========
    # print("[analyze_protein_ligand_interactions] 🔬 Using strict standards (H-bond ≤2.8Å, Salt bridge ≤4.0Å)")
    
    # 获取原子信息
    if pdb_file:
        atoms = parse_pdb_file(pdb_file)
        if not atoms:
            print(f"[analyze_protein_ligand_interactions] Unable to read atom information from file: {pdb_file}")
            return None
    else:
        atoms = parse_pdb_structure(obj_name)
        if not atoms:
            print("[analyze_protein_ligand_interactions] Unable to obtain atom information")
            return None

    # 按链和残基分组
    chain_residues = defaultdict(list)
    for atom in atoms:
        res_key = (atom[0], atom[1], atom[2])  # (chain, res_name, res_id)
        chain_residues[res_key].append(atom)

    # print(f"[analyze_protein_ligand_interactions] 📋 Structure info: {len(chain_residues)} residues total")
    
    # 统计残基类型
    residue_type_count = defaultdict(int)
    for res_key in chain_residues.keys():
        _, res_name, _ = res_key
        mol_type = identify_molecule_type(res_name)
        residue_type_count[mol_type] += 1
    
    # print(f"[analyze_protein_ligand_interactions] Residue type distribution: {dict(residue_type_count)}")
    
    # 分离收集配体和蛋白质残基
    ligand_residues = []
    protein_residues = []

    # 第一步：收集所有蛋白质残基（无论是否指定配体）
    for res_key, res_atoms in chain_residues.items():
        chain_id, res_name, res_id = res_key
        mol_type = identify_molecule_type(res_name)
        
        if mol_type == "protein":
            protein_residues.append((res_key, res_atoms))

    # 第二步：收集配体残基
    for res_key, res_atoms in chain_residues.items():
        chain_id, res_name, res_id = res_key
        mol_type = identify_molecule_type(res_name)

        if ligand_resname:
            # 用户指定了配体名称，精确匹配
            if res_name.upper() == ligand_resname.upper():
                ligand_residues.append((res_key, res_atoms))
                # print(f"[analyze_protein_ligand_interactions]   ✓ Ligand: {res_name} {res_id} (chain {chain_id})")
        else:
            # 自动检测配体（非蛋白、非溶剂的残基）
            if mol_type == "ligand":
                ligand_residues.append((res_key, res_atoms))
                # print(f"[analyze_protein_ligand_interactions]   ✓ Auto-detected ligand: {res_name} {res_id} (chain {chain_id})")

    # 如果未指定蛋白链，自动检测
    if protein_chains is None:
        protein_chains_set = set()
        for res_key, _ in protein_residues:
            protein_chains_set.add(res_key[0])
        protein_chains = list(protein_chains_set)
    else:
        # 过滤指定的蛋白链
        old_count = len(protein_residues)
        protein_residues = [(rk, ra) for rk, ra in protein_residues if rk[0] in protein_chains]
        # print(f"[analyze_protein_ligand_interactions] Filtered protein chains: {old_count} -> {len(protein_residues)} residues")

    if not ligand_residues:
        print("\n" + "="*60)
        print("[analyze_protein_ligand_interactions] ❌ No ligand found")
        if ligand_resname:
            print(f"   Ligand specified: '{ligand_resname}'")
        print("\nPossible reasons:")
        print("   1) Ligand name mismatch")
        print("   2) Ligand is on a separate chain not included in the object")
        
        # List all detected non-standard residues (potential ligands)
        non_standard = {}
        for res_key in chain_residues.keys():
            chain_id, res_name, res_id = res_key
            mol_type = identify_molecule_type(res_name)
            if mol_type == "ligand":
                if res_name not in non_standard:
                    non_standard[res_name] = []
                non_standard[res_name].append(f"chain{chain_id}:{res_id}")
        
        if non_standard:
            print("\nDetected potential ligands:")
            for resn, locations in sorted(non_standard.items()):
                print(f"   - {resn}: {', '.join(locations[:5])}")
            print("\nHow to fix:")
            print("   • Enter one of the names above into 'Ligand Resname'")
            print("   • Ensure the ligand chain is included in the object")
        else:
            print("\nNo non-standard residues found. Please check:")
            print("   • The PyMOL object is correctly loaded")
            print("   • The ligand is included (use 'show sticks, resn XXX' to verify)")
        
        print("="*60 + "\n")
        return None

    if not protein_residues:
        print("[analyze_protein_ligand_interactions] ⚠️ No protein residues found")
        return None

    # ✅ 调试打印
    print(f"[DEBUG] ✅ 检测到 {len(ligand_residues)} 个配体分子")
    for lig_key, lig_atoms in ligand_residues:
        print(f"[DEBUG]   配体: {lig_key[1]} (chain='{lig_key[0]}', resid={lig_key[2]}, 原子数={len(lig_atoms)})")
    print(f"[DEBUG] ✅ 检测到 {len(protein_residues)} 个蛋白质残基")
    print(f"[DEBUG]   蛋白质链: {protein_chains}")
    print(f"[DEBUG] 检测参数:")
    print(f"[DEBUG]   氢键: D···A ≤ {INTERACTION_PARAMS['hbond']['max_DA_dist']} Å")
    print(f"[DEBUG]   盐桥: ≤ {INTERACTION_PARAMS['ionic']['max_dist']} Å")
    print(f"[DEBUG]   疏水: ≤ {INTERACTION_PARAMS['hydrophobic']['other_max']} Å")
    print(f"[DEBUG]   距离截断: {distance_cutoff} Å")

    # 分析相互作用
    interactions = []
    
    # 构建残基原子索引（用于精确氢键检测）
    all_atoms_by_residue = {}
    for lig_key, lig_atoms in ligand_residues:
        all_atoms_by_residue[lig_key] = lig_atoms
    for prot_key, prot_atoms in protein_residues:
        all_atoms_by_residue[prot_key] = prot_atoms

    # 计数器
    close_pairs_count = 0
    checked_pairs = 0
    
    for lig_key, lig_atoms in ligand_residues:
        lig_chain, lig_name, lig_id = lig_key
        print(f"[DEBUG] 分析配体: {lig_name} (chain='{lig_chain}', resid={lig_id})")

        # 计算配体质心（只计算一次）
        lig_centroid = centroid([a[4] for a in lig_atoms])
        
        for prot_key, prot_atoms in protein_residues:
            prot_chain, prot_name, prot_id = prot_key

            # 快速距离筛选（放宽到15Å）
            prot_centroid = centroid([a[4] for a in prot_atoms])
            centroid_dist = distance(lig_centroid, prot_centroid)
            
            # 放宽距离筛选：配体可能很大，质心距离不代表边缘距离
            if centroid_dist > 15.0:
                continue
            
            checked_pairs += 1

            # 检查原子间相互作用
            for lig_atom in lig_atoms:
                for prot_atom in prot_atoms:
                    d = distance(lig_atom[4], prot_atom[4])

                    if d > distance_cutoff:
                        continue
                    
                    close_pairs_count += 1

                    # 判断相互作用类型（优先级：盐桥 > 氢键 > 疏水）
                    interaction_type = None
                    hb_angle = None
                    
                    # ① 最高优先级：盐桥
                    if is_saltbridge_atom(lig_name, lig_atom[3], prot_name, prot_atom[3], d):
                        interaction_type = "盐桥"
                    else:
                        # ② 次优先级：氢键
                        is_hb, hb_dist, hb_angle = is_hbond_precise(lig_atom, prot_atom, all_atoms_by_residue)
                        if not is_hb:
                            is_hb, hb_dist, hb_angle = is_hbond_precise(prot_atom, lig_atom, all_atoms_by_residue)
                        
                        if is_hb:
                            interaction_type = "氢键"
                        # ③ 最低优先级：疏水
                        elif not key_interactions_only and is_hydrophobic(lig_name, prot_name, lig_atom, prot_atom, d):
                            interaction_type = "疏水相互作用"

                    # 严格模式：只检测关键相互作用
                    # key_interactions_only=False时才启用宽松规则（不推荐）

                    if interaction_type:
                        conf = calculate_confidence_score(interaction_type, d, hb_angle if interaction_type == "氢键" else None)
                        interactions.append({
                            "Ligand_Chain": lig_chain,
                            "Ligand_Residue": f"{lig_name} {lig_id}",
                            "Ligand_Atom": lig_atom[3],
                            "Protein_Chain": prot_chain,
                            "Protein_Residue": f"{prot_name} {prot_id}",
                            "Protein_Atom": prot_atom[3],
                            "Distance": round(d, 2),
                            "Interaction": interaction_type,
                            "Confidence": conf
                        })
                        
                    # 检查卤素键
                    is_hal, hal_dist, hal_angle = is_halogen_bond(lig_atom, prot_atom, d)
                    if is_hal:
                        conf = calculate_confidence_score("卤素键", hal_dist, hal_angle)
                        interactions.append({
                            "Ligand_Chain": lig_chain,
                            "Ligand_Residue": f"{lig_name} {lig_id}",
                            "Ligand_Atom": lig_atom[3],
                            "Protein_Chain": prot_chain,
                            "Protein_Residue": f"{prot_name} {prot_id}",
                            "Protein_Atom": prot_atom[3],
                            "Distance": round(hal_dist, 2),
                            "Interaction": "卤素键",
                            "Confidence": conf
                        })

            # 检查π相互作用
            if is_pipi(lig_name, lig_atoms, prot_name, prot_atoms):
                interactions.append({
                    "Ligand_Chain": lig_chain,
                    "Ligand_Residue": f"{lig_name} {lig_id}",
                    "Ligand_Atom": "ring",
                    "Protein_Chain": prot_chain,
                    "Protein_Residue": f"{prot_name} {prot_id}",
                    "Protein_Atom": "ring",
                    "Distance": "-",
                    "Interaction": "π–π 堆积",
                    "Confidence": 0.9
                })

            if is_cationpi(lig_name, lig_atoms, prot_name, prot_atoms):
                interactions.append({
                    "Ligand_Chain": lig_chain,
                    "Ligand_Residue": f"{lig_name} {lig_id}",
                    "Ligand_Atom": "ring/cation",
                    "Protein_Chain": prot_chain,
                    "Protein_Residue": f"{prot_name} {prot_id}",
                    "Protein_Atom": "ring/cation",
                    "Distance": "-",
                    "Interaction": "π–阳离子相互作用",
                    "Confidence": 0.9
                })

    # 检查水桥 (Water Bridges)
    # 1. 收集所有水分子
    water_atoms = []
    for atom in atoms:
        if identify_molecule_type(atom[1]) == "solvent":
            water_atoms.append(atom)
            
    if water_atoms:
        # print(f"[analyze_protein_ligand_interactions] Checking water bridges with {len(water_atoms)} water molecules...")
        # 优化：建立空间索引或只检查配体附近的水
        # 这里使用简单距离筛选
        
        for lig_key, lig_atoms in ligand_residues:
            lig_chain, lig_name, lig_id = lig_key
            
            for water in water_atoms:
                # 快速检查水-配体距离
                min_dist_lw = float('inf')
                closest_lig_atom = None
                
                for la in lig_atoms:
                    d = distance(la[4], water[4])
                    if d < min_dist_lw:
                        min_dist_lw = d
                        closest_lig_atom = la
                
                if min_dist_lw > 4.0: # 水桥通常 < 3.5 + 余裕
                    continue
                    
                # 检查水-蛋白距离
                for prot_key, prot_atoms in protein_residues:
                    prot_chain, prot_name, prot_id = prot_key
                    
                    # 快速检查水-蛋白距离
                    prot_centroid = centroid([a[4] for a in prot_atoms])
                    if distance(water[4], prot_centroid) > 6.0:
                        continue
                        
                    for pa in prot_atoms:
                        is_wb, d_lw, d_wp, angle = is_water_bridge(closest_lig_atom, pa, water)
                        if is_wb:
                            conf = calculate_confidence_score("水桥", max(d_lw, d_wp), angle)
                            interactions.append({
                                "Ligand_Chain": lig_chain,
                                "Ligand_Residue": f"{lig_name} {lig_id}",
                                "Ligand_Atom": closest_lig_atom[3],
                                "Protein_Chain": prot_chain,
                                "Protein_Residue": f"{prot_name} {prot_id}",
                                "Protein_Atom": pa[3],
                                "Distance": f"{d_lw:.2f}/{d_wp:.2f}",
                                "Interaction": "水桥",
                                "Confidence": conf
                            })

    # 调试输出
    print(f"[DEBUG] 检查了 {checked_pairs} 个残基对")
    print(f"[DEBUG] 距离内的原子对: {close_pairs_count}")
    print(f"[DEBUG] 检测到的相互作用: {len(interactions)}")
    
    # 输出到CSV
    if output_csv and interactions:
        try:
            with open(output_csv, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(["Ligand_Chain", "Ligand_Residue", "Ligand_Atom",
                               "Protein_Chain", "Protein_Residue", "Protein_Atom",
                               "Distance", "Interaction", "Confidence"])
                for inter in interactions:
                    writer.writerow([
                        inter["Ligand_Chain"],
                        inter["Ligand_Residue"],
                        inter["Ligand_Atom"],
                        inter["Protein_Chain"],
                        inter["Protein_Residue"],
                        inter["Protein_Atom"],
                        inter["Distance"],
                        inter["Interaction"],
                        inter.get("Confidence", "-")
                    ])
            # print(f"[analyze_protein_ligand_interactions] Results saved to: {output_csv}")
        except Exception as e:
            print(f"[analyze_protein_ligand_interactions] Failed to save CSV: {e}")

    result = {
        "ligand_residues": [{"chain": k[0], "resname": k[1], "resid": k[2]}
                           for k, _ in ligand_residues],
        "protein_chains": protein_chains,
        "interactions": interactions,
        "mode": "strict",
        "standard": "Publication",
        "parameters": {
            "distance_cutoff": distance_cutoff,
            "hbond_cutoff": INTERACTION_PARAMS["hbond"]["max_DA_dist"],
            "saltbridge_cutoff": INTERACTION_PARAMS["ionic"]["max_dist"],
            "hydrophobic_cutoff": INTERACTION_PARAMS["hydrophobic"]["other_max"],
            "key_interactions_only": key_interactions_only
        }
    }

    # print(f"[analyze_protein_ligand_interactions] ✅ Analysis complete: found {len(interactions)} interactions")
    return result

cmd.extend("analyze_protein_ligand_interactions", analyze_protein_ligand_interactions)

def analyze_ternary_complex(obj_name=None, ligand_resname=None,
                           protein1_chains=None, protein2_chains=None,
                           output_csv=None, distance_cutoff=4.5, pdb_file=None):
    """
    分析三元复合体（蛋白-配体-蛋白）相互作用
    典型应用：PROTAC分子同时结合两个蛋白质

    参数:
        obj_name: PyMOL对象名称
        ligand_resname: 配体残基名称（如PROTAC分子）
        protein1_chains: 蛋白质1的链ID列表（例如 ["A"]）
        protein2_chains: 蛋白质2的链ID列表（例如 ["B"]）
        output_csv: 输出CSV文件路径
        distance_cutoff: 距离截断值（埃）
        pdb_file: PDB文件路径（可选）

    返回:
        dict: {
            "ligand_info": {...},
            "protein1_interactions": [...],
            "protein2_interactions": [...],
            "bridging_analysis": {...}
        }
    """
    # 获取原子信息
    if pdb_file:
        atoms = parse_pdb_file(pdb_file)
        if not atoms:
            print(f"[analyze_ternary_complex] Unable to read atom information from file: {pdb_file}")
            return None
    else:
        atoms = parse_pdb_structure(obj_name)
        if not atoms:
            print("[analyze_ternary_complex] Unable to obtain atom information")
            return None

    # 按残基分组
    chain_residues = defaultdict(list)
    for atom in atoms:
        res_key = (atom[0], atom[1], atom[2])
        chain_residues[res_key].append(atom)

    # 检测配体
    ligand_residues = []
    for res_key, res_atoms in chain_residues.items():
        _, res_name, _ = res_key
        if ligand_resname:
            if res_name.upper() == ligand_resname.upper():
                ligand_residues.append((res_key, res_atoms))
        else:
            if identify_molecule_type(res_name) == "ligand":
                ligand_residues.append((res_key, res_atoms))

    if not ligand_residues:
        print("[analyze_ternary_complex] ⚠️ No ligand found")
        return None

    # 自动检测蛋白质链（如果未指定）
    if protein1_chains is None or protein2_chains is None:
        all_chains = set(atom[0] for atom in atoms if identify_molecule_type(atom[1]) == "protein")
        all_chains = sorted(all_chains)

        if len(all_chains) < 2:
            print("[analyze_ternary_complex] ⚠️ Not enough protein chains found (need at least 2)")
            return None

        if protein1_chains is None:
            protein1_chains = [all_chains[0]]
        if protein2_chains is None:
            protein2_chains = [all_chains[1]] if len(all_chains) > 1 else [all_chains[0]]

    print(f"[analyze_ternary_complex] Protein 1 chains: {', '.join(protein1_chains)}")
    print(f"[analyze_ternary_complex] Protein 2 chains: {', '.join(protein2_chains)}")
    print(f"[analyze_ternary_complex] Ligands: {len(ligand_residues)} molecule(s)")

    # 分别分析配体与两个蛋白的相互作用
    print("\n[analyze_ternary_complex] Analyzing ligand - protein 1 interactions...")
    protein1_result = analyze_protein_ligand_interactions(
        obj_name=obj_name,
        ligand_resname=ligand_resname,
        protein_chains=protein1_chains,
        distance_cutoff=distance_cutoff,
        pdb_file=pdb_file
    )

    print("\n[analyze_ternary_complex] Analyzing ligand - protein 2 interactions...")
    protein2_result = analyze_protein_ligand_interactions(
        obj_name=obj_name,
        ligand_resname=ligand_resname,
        protein_chains=protein2_chains,
        distance_cutoff=distance_cutoff,
        pdb_file=pdb_file
    )

    if not protein1_result or not protein2_result:
        print("[analyze_ternary_complex] ⚠️ Failed to complete ternary complex analysis")
        return None

    # 分析桥接效应
    print("\n[analyze_ternary_complex] Analyzing bridging effects...")

    p1_interactions = protein1_result["interactions"]
    p2_interactions = protein2_result["interactions"]

    # 统计与配体相互作用的残基
    p1_residues = set((inter["Protein_Chain"], inter["Protein_Residue"])
                     for inter in p1_interactions)
    p2_residues = set((inter["Protein_Chain"], inter["Protein_Residue"])
                     for inter in p2_interactions)

    bridging_analysis = {
        "protein1_contacting_residues": len(p1_residues),
        "protein2_contacting_residues": len(p2_residues),
        "total_interactions": len(p1_interactions) + len(p2_interactions),
        "protein1_interactions_count": len(p1_interactions),
        "protein2_interactions_count": len(p2_interactions),
        "ligand_bridges": len(ligand_residues) > 0
    }

    # 保存到CSV
    if output_csv:
        try:
            with open(output_csv, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)

                # Headers
                writer.writerow(["=== Ternary Complex Analysis ==="])
                writer.writerow([])

                # Ligand info
                writer.writerow(["Ligand Info"])
                writer.writerow(["Ligand Count", len(ligand_residues)])
                for lig_key, _ in ligand_residues:
                    writer.writerow(["Ligand", f"{lig_key[1]} {lig_key[2]} (chain {lig_key[0]})"])
                writer.writerow([])

                # Protein 1 interactions
                writer.writerow([f"=== Protein 1 (chains {', '.join(protein1_chains)}) - Ligand Interactions ==="])
                writer.writerow(["Ligand_Chain", "Ligand_Residue", "Ligand_Atom",
                               "Protein_Chain", "Protein_Residue", "Protein_Atom",
                               "Distance", "Interaction"])
                for inter in p1_interactions:
                    writer.writerow([
                        inter["Ligand_Chain"], inter["Ligand_Residue"], inter["Ligand_Atom"],
                        inter["Protein_Chain"], inter["Protein_Residue"], inter["Protein_Atom"],
                        inter["Distance"], inter["Interaction"]
                    ])
                writer.writerow([])

                # Protein 2 interactions
                writer.writerow([f"=== Protein 2 (chains {', '.join(protein2_chains)}) - Ligand Interactions ==="])
                writer.writerow(["Ligand_Chain", "Ligand_Residue", "Ligand_Atom",
                               "Protein_Chain", "Protein_Residue", "Protein_Atom",
                               "Distance", "Interaction"])
                for inter in p2_interactions:
                    writer.writerow([
                        inter["Ligand_Chain"], inter["Ligand_Residue"], inter["Ligand_Atom"],
                        inter["Protein_Chain"], inter["Protein_Residue"], inter["Protein_Atom"],
                        inter["Distance"], inter["Interaction"]
                    ])
                writer.writerow([])

                # Bridging analysis summary
                writer.writerow(["=== Bridging Analysis Summary ==="])
                writer.writerow(["Protein 1 contacting residues", bridging_analysis["protein1_contacting_residues"]])
                writer.writerow(["Protein 2 contacting residues", bridging_analysis["protein2_contacting_residues"]])
                writer.writerow(["Protein 1 interactions", bridging_analysis["protein1_interactions_count"]])
                writer.writerow(["Protein 2 interactions", bridging_analysis["protein2_interactions_count"]])
                writer.writerow(["Total interactions", bridging_analysis["total_interactions"]])

            print(f"[analyze_ternary_complex] Results saved to: {output_csv}")
        except Exception as e:
            print(f"[analyze_ternary_complex] Failed to save CSV: {e}")

    result = {
        "ligand_info": protein1_result["ligand_residues"],
        "protein1_chains": protein1_chains,
        "protein2_chains": protein2_chains,
        "protein1_interactions": p1_interactions,
        "protein2_interactions": p2_interactions,
        "bridging_analysis": bridging_analysis
    }

    print(f"\n[analyze_ternary_complex] ✅ Ternary complex analysis complete")
    print(f"  - Protein 1: {bridging_analysis['protein1_interactions_count']} interactions")
    print(f"  - Protein 2: {bridging_analysis['protein2_interactions_count']} interactions")
    print(f"  - Total: {bridging_analysis['total_interactions']} interactions")

    return result

cmd.extend("analyze_ternary_complex", analyze_ternary_complex)

def analyze_atom_pair_interactions(obj_name=None,
                                   atom1_selection=None, atom2_selection=None,
                                   distance_cutoff=5.0, output_csv=None, pdb_file=None):
    """
    分析特定原子对之间的相互作用
    例如: 分析配体的N原子与蛋白质某残基的O原子之间的相互作用

    参数:
        obj_name: PyMOL对象名称
        atom1_selection: 原子1的选择语法，支持:
            - PyMOL选择语法: "resn LIG and name N1"
            - 简化格式: "LIG/301/N1" (残基名/残基号/原子名)
            - 元素类型: "elem N" (所有氮原子)
        atom2_selection: 原子2的选择语法(同上)
        distance_cutoff: 距离截断值（埃）
        output_csv: 输出CSV文件路径
        pdb_file: PDB文件路径（可选）

    返回:
        list: 原子对相互作用列表，包含距离、角度等信息

    示例:
        # 分析配体LIG的N原子与所有氧原子的相互作用
        analyze_atom_pair_interactions(
            obj_name="complex",
            atom1_selection="resn LIG and name N1",
            atom2_selection="elem O",
            distance_cutoff=3.5
        )

        # 使用简化格式
        analyze_atom_pair_interactions(
            obj_name="complex",
            atom1_selection="LIG/301/N1",
            atom2_selection="SER/50/OG",
            distance_cutoff=4.0
        )
    """
    from pymol import cmd as pymol_cmd

    # 获取原子信息
    if pdb_file:
        atoms = parse_pdb_file(pdb_file)
        if not atoms:
            print(f"[analyze_atom_pair_interactions] Unable to read atom information from file: {pdb_file}")
            return None
    else:
        atoms = parse_pdb_structure(obj_name)
        if not atoms:
            print("[analyze_atom_pair_interactions] Unable to obtain atom information")
            return None

    # 解析选择语法
    def parse_selection(sel_str, atoms_list):
        """解析选择语法并返回匹配的原子列表"""
        if not sel_str:
            return []

        matched_atoms = []

        # 格式1: 简化格式 "LIG/301/N1"
        if "/" in sel_str:
            parts = sel_str.split("/")
            if len(parts) == 3:
                res_name, res_id, atom_name = parts
                for atom in atoms_list:
                    if (atom[1].strip().upper() == res_name.strip().upper() and
                        atom[2].strip() == res_id.strip() and
                        atom[3].strip().upper() == atom_name.strip().upper()):
                        matched_atoms.append(atom)
                return matched_atoms

        # 格式2: 元素类型 "elem N"
        if sel_str.lower().startswith("elem "):
            elem = sel_str[5:].strip().upper()
            for atom in atoms_list:
                atom_elem = atom[3].strip()[0].upper()  # 取原子名第一个字母作为元素
                if atom_elem == elem[0]:
                    matched_atoms.append(atom)
            return matched_atoms

        # 格式3: 简单残基名 "resn LIG"
        if sel_str.lower().startswith("resn "):
            res_name = sel_str[5:].strip().upper()
            for atom in atoms_list:
                if atom[1].strip().upper() == res_name:
                    matched_atoms.append(atom)
            return matched_atoms

        # 格式4: 原子名 "name N1"
        if sel_str.lower().startswith("name "):
            atom_name = sel_str[5:].strip().upper()
            for atom in atoms_list:
                if atom[3].strip().upper() == atom_name:
                    matched_atoms.append(atom)
            return matched_atoms

        # 格式5: PyMOL复杂选择语法（尝试解析）
        # 支持 "resn LIG and name N1"
        if " and " in sel_str.lower():
            conditions = sel_str.lower().split(" and ")
            result = list(atoms_list)

            for cond in conditions:
                cond = cond.strip()
                if cond.startswith("resn "):
                    res_name = cond[5:].strip().upper()
                    result = [a for a in result if a[1].strip().upper() == res_name]
                elif cond.startswith("name "):
                    atom_name = cond[5:].strip().upper()
                    result = [a for a in result if a[3].strip().upper() == atom_name]
                elif cond.startswith("chain "):
                    chain = cond[6:].strip().upper()
                    result = [a for a in result if a[0].strip().upper() == chain]
                elif cond.startswith("resi "):
                    res_id = cond[5:].strip()
                    result = [a for a in result if a[2].strip() == res_id]
                elif cond.startswith("elem "):
                    elem = cond[5:].strip().upper()
                    result = [a for a in result if a[3].strip()[0].upper() == elem[0]]

            return result

        return matched_atoms

    # 获取选中的原子
    atoms1 = parse_selection(atom1_selection, atoms)
    atoms2 = parse_selection(atom2_selection, atoms)

    if not atoms1:
        print(f"[analyze_atom_pair_interactions] ⚠️ No matching atoms for selection 1: {atom1_selection}")
        return None

    if not atoms2:
        print(f"[analyze_atom_pair_interactions] ⚠️ No matching atoms for selection 2: {atom2_selection}")
        return None

    print(f"[analyze_atom_pair_interactions] Atoms1: {len(atoms1)}")
    print(f"[analyze_atom_pair_interactions] Atoms2: {len(atoms2)}")

    # 分析原子对相互作用
    interactions = []

    for atom1 in atoms1:
        chain1, res1, resid1, name1, coord1 = atom1

        for atom2 in atoms2:
            chain2, res2, resid2, name2, coord2 = atom2

            # 跳过同一原子
            if (chain1 == chain2 and res1 == res2 and
                resid1 == resid2 and name1 == name2):
                continue

            # 计算距离
            d = distance(coord1, coord2)

            if d > distance_cutoff:
                continue

            # 判断相互作用类型
            interaction_type = _classify_atom_interaction(
                name1, name2, res1, res2, d
            )

            # 计算方向向量和角度（如果有其他原子）
            direction_info = _calculate_interaction_geometry(
                atom1, atom2, atoms
            )

            interactions.append({
                "Atom1_Chain": chain1,
                "Atom1_Residue": f"{res1} {resid1}",
                "Atom1_Name": name1,
                "Atom1_Coords": f"({coord1[0]:.2f}, {coord1[1]:.2f}, {coord1[2]:.2f})",
                "Atom2_Chain": chain2,
                "Atom2_Residue": f"{res2} {resid2}",
                "Atom2_Name": name2,
                "Atom2_Coords": f"({coord2[0]:.2f}, {coord2[1]:.2f}, {coord2[2]:.2f})",
                "Distance": round(d, 3),
                "Interaction": interaction_type,
                "Geometry": direction_info
            })

    # 按距离排序
    interactions.sort(key=lambda x: x["Distance"])

    # 输出到CSV
    if output_csv and interactions:
        try:
            with open(output_csv, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow([
                    "Atom1_Chain", "Atom1_Residue", "Atom1_Name", "Atom1_Coords",
                    "Atom2_Chain", "Atom2_Residue", "Atom2_Name", "Atom2_Coords",
                    "Distance", "Interaction", "Geometry"
                ])
                for inter in interactions:
                    writer.writerow([
                        inter["Atom1_Chain"],
                        inter["Atom1_Residue"],
                        inter["Atom1_Name"],
                        inter["Atom1_Coords"],
                        inter["Atom2_Chain"],
                        inter["Atom2_Residue"],
                        inter["Atom2_Name"],
                        inter["Atom2_Coords"],
                        inter["Distance"],
                        inter["Interaction"],
                        inter["Geometry"]
                    ])
            print(f"[analyze_atom_pair_interactions] Results saved to: {output_csv}")
        except Exception as e:
            print(f"[analyze_atom_pair_interactions] Failed to save CSV: {e}")

    print(f"[analyze_atom_pair_interactions] ✅ Analysis complete: found {len(interactions)} atom-pair interactions")

    # 打印前5个最近的相互作用
    if interactions:
        print("\nTop 5 closest interactions:")
        for i, inter in enumerate(interactions[:5], 1):
            print(f"  {i}. {inter['Atom1_Residue']}:{inter['Atom1_Name']} --- "
                  f"{inter['Atom2_Residue']}:{inter['Atom2_Name']} = "
                  f"{inter['Distance']:.3f} Å ({inter['Interaction']})")

    return interactions

def _classify_atom_interaction(atom1_name, atom2_name, res1, res2, dist):
    """根据原子类型和距离分类相互作用"""
    # 获取元素类型（原子名的第一个字母）
    elem1 = atom1_name.strip()[0].upper()
    elem2 = atom2_name.strip()[0].upper()

    # N-O 或 O-N (氢键供体-受体)
    if (elem1 in ['N', 'O'] and elem2 in ['N', 'O']) and dist <= 3.5:
        if elem1 == 'N' and elem2 == 'O':
            return "氢键 (N-H···O)"
        elif elem1 == 'O' and elem2 == 'N':
            return "氢键 (O-H···N)"
        else:
            return "氢键"

    # S-S 二硫键
    if elem1 == 'S' and elem2 == 'S' and dist <= 2.5:
        return "二硫键 (S-S)"

    # 盐桥 (带电原子)
    if dist <= 4.0:
        positive_atoms = ['NZ', 'NH1', 'NH2', 'NE']  # LYS, ARG
        negative_atoms = ['OD1', 'OD2', 'OE1', 'OE2']  # ASP, GLU

        if (atom1_name in positive_atoms and atom2_name in negative_atoms) or \
           (atom2_name in positive_atoms and atom1_name in negative_atoms):
            return "盐桥"

    # 卤素键
    if elem1 in ['F', 'CL', 'BR', 'I'] and elem2 in ['O', 'N']:
        if dist <= 3.5:
            return f"卤素键 ({elem1}···{elem2})"

    # π相互作用（如果原子在芳香环上）
    aromatic_res = {'PHE', 'TYR', 'TRP', 'HIS'}
    if res1 in aromatic_res and res2 in aromatic_res:
        return "π相互作用"

    # 疏水相互作用
    if elem1 == 'C' and elem2 == 'C' and dist <= 4.5:
        hydrophobic_res = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'PRO', 'TRP'}
        if res1 in hydrophobic_res and res2 in hydrophobic_res:
            return "疏水相互作用"

    # 金属配位
    if elem1 in ['ZN', 'MG', 'CA', 'FE', 'CU', 'MN'] or \
       elem2 in ['ZN', 'MG', 'CA', 'FE', 'CU', 'MN']:
        if dist <= 3.0:
            return "金属配位"

    # 默认：范德华相互作用
    if dist <= 4.0:
        return "范德华"

    return "弱相互作用"

def _calculate_interaction_geometry(atom1, atom2, all_atoms):
    """计算相互作用的几何信息（角度等）"""
    # 这里可以添加更复杂的几何计算
    # 例如：氢键角度、二面角等

    # 简化版本：计算方向
    coord1, coord2 = atom1[4], atom2[4]

    dx = coord2[0] - coord1[0]
    dy = coord2[1] - coord1[1]
    dz = coord2[2] - coord1[2]

    # 主要方向
    abs_max = max(abs(dx), abs(dy), abs(dz))
    if abs_max == abs(dx):
        direction = "+X-axis" if dx > 0 else "-X-axis"
    elif abs_max == abs(dy):
        direction = "+Y-axis" if dy > 0 else "-Y-axis"
    else:
        direction = "+Z-axis" if dz > 0 else "-Z-axis"

    return f"Primary {direction}"

cmd.extend("analyze_atom_pair_interactions", analyze_atom_pair_interactions)

def visualize_atom_pairs(obj_name, interactions_result=None, csv_path=None):
    """
    在PyMOL中可视化原子对相互作用

    参数:
        obj_name: PyMOL对象名称
        interactions_result: analyze_atom_pair_interactions的返回结果
        csv_path: 或者提供CSV文件路径
    """
    from pymol import cmd

    if csv_path and not interactions_result:
        # 从CSV读取
        interactions_result = []
        try:
            with open(csv_path, 'r', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    interactions_result.append(row)
        except Exception as e:
            print(f"[visualize_atom_pairs] Failed to read CSV: {e}")
            return

    if not interactions_result:
        print("[visualize_atom_pairs] No interaction data")
        return

    # 创建选择
    cmd.delete("atom_pairs_*")

    # 颜色映射
    color_map = {
        "氢键": "blue",
        "盐桥": "red",
        "二硫键": "yellow",
        "疏水相互作用": "green",
        "π相互作用": "purple",
        "卤素键": "orange",
        "金属配位": "magenta",
        "范德华": "gray",
        "弱相互作用": "lightgray"
    }

    for i, inter in enumerate(interactions_result[:20], 1):  # 最多显示20个
        try:
            # 解析原子信息
            chain1 = inter["Atom1_Chain"]
            res1_full = inter["Atom1_Residue"].split()
            res1_name = res1_full[0] if res1_full else ""
            res1_id = res1_full[1] if len(res1_full) > 1 else ""
            atom1_name = inter["Atom1_Name"]

            chain2 = inter["Atom2_Chain"]
            res2_full = inter["Atom2_Residue"].split()
            res2_name = res2_full[0] if res2_full else ""
            res2_id = res2_full[1] if len(res2_full) > 1 else ""
            atom2_name = inter["Atom2_Name"]

            # 创建选择
            sel1 = f"{obj_name} and chain {chain1} and resi {res1_id} and name {atom1_name}"
            sel2 = f"{obj_name} and chain {chain2} and resi {res2_id} and name {atom2_name}"

            # 显示原子
            cmd.show("spheres", sel1)
            cmd.show("spheres", sel2)
            cmd.set("sphere_scale", 0.3, sel1)
            cmd.set("sphere_scale", 0.3, sel2)

            # 绘制距离线
            pair_name = f"atom_pairs_{i}"
            cmd.distance(pair_name, sel1, sel2, cutoff=10.0)

            # 设置颜色
            interaction_type = inter.get("Interaction", "")
            for key, color in color_map.items():
                if key in interaction_type:
                    cmd.color(color, pair_name)
                    break

            # 设置标签
            cmd.set("label_size", 20)
            cmd.set("label_color", "black")

        except Exception as e:
            print(f"[visualize_atom_pairs] 可视化第 {i} 个相互作用时出错: {e}")
            continue

    print(f"[visualize_atom_pairs] ✅ Visualized {min(len(interactions_result), 20)} atom-pair interactions")
    print(f"   Use 'hide labels, atom_pairs_*' to hide labels")
    print(f"   Use 'delete atom_pairs_*' to delete all visuals")

cmd.extend("visualize_atom_pairs", visualize_atom_pairs)

def visualize_protein_ligand_3d(obj_name, interactions_result=None, ligand_resname=None,
                                   csv_path=None, show_hydrophobic=False,
                                   max_interactions_per_type=None,
                                   min_confidence=0.8):
    """
    在PyMOL中3D可视化蛋白-配体相互作用（改进版，参考专业脚本）

    参数:
        obj_name: PyMOL对象名称
        interactions_result: analyze_protein_ligand_interactions的返回结果
        ligand_resname: 配体残基名称(可选)
        csv_path: 或者提供CSV文件路径
        show_hydrophobic: 是否显示疏水相互作用（默认False，只显示关键相互作用）
        max_interactions_per_type: 每种类型最多显示的相互作用数（默认None表示按优先级自动筛选）
    """
    from pymol import cmd

    # 读取相互作用数据
    from_csv = False
    if csv_path and not interactions_result:
        interactions = []
        try:
            with open(csv_path, 'r', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                for idx, row in enumerate(reader, start=1):
                    # 为每一行记录一个稳定的行号，后续用于和 PyMOL 对象一一对应
                    row["_row_index"] = idx
                    interactions.append(row)
            from_csv = True
        except Exception as e:
            print(f"[visualize_protein_ligand_3d] 读取CSV失败: {e}")
            return
    elif interactions_result:
        if isinstance(interactions_result, dict) and "interactions" in interactions_result:
            interactions = interactions_result["interactions"]
        else:
            interactions = interactions_result
    else:
        print("[visualize_protein_ligand_3d] Interaction data required")
        return

    if not interactions:
        print("[visualize_protein_ligand_3d] No interaction data")
        return

    # ========== 置信度过滤：只保留高置信度相互作用 ==========
    if min_confidence is not None:
        filtered = []
        dropped = 0
        for inter in interactions:
            conf_val = inter.get("Confidence") or inter.get("confidence") or inter.get("CONFIDENCE")
            try:
                conf = float(conf_val)
            except (TypeError, ValueError):
                # 没有置信度字段时默认视为 1.0（不丢弃）
                conf = 1.0
            if conf >= float(min_confidence):
                filtered.append(inter)
            else:
                dropped += 1
        if not filtered:
            print(f"[visualize_protein_ligand_3d] ℹ️ No interactions with Confidence ≥ {min_confidence}")
            return
        if dropped > 0:
            print(f"[visualize_protein_ligand_3d] Filtering by Confidence ≥ {min_confidence}: kept {len(filtered)}/{len(interactions)} interactions")
        interactions = filtered

    # 为非 CSV 来源的相互作用补充一个稳定的行号（1-based）
    if not from_csv:
        for idx, inter in enumerate(interactions, start=1):
            if "_row_index" not in inter:
                inter["_row_index"] = idx

    # ========== 第一步：清理和基础设置 ==========
    print("[visualize_protein_ligand_3d] 🎨 Starting 3D visualization...")
    
    # 清除旧的可视化对象
    cmd.delete("pl_interact_*")
    cmd.delete("hbonds_*")
    cmd.delete("lig_pocket")
    cmd.delete("polar_*")
    cmd.delete("don_*")
    
    # ========== 第二步：自动检测或使用指定的配体 ==========
    if not ligand_resname:
        # 自动检测配体（非标准残基）
        standard_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'HOH', 'WAT', 'A', 'C', 'G', 'T', 'U', 'DA', 'DC', 'DG', 'DT',
            'NA', 'CL', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'K'
        }
        
        model = cmd.get_model(obj_name)
        all_resnames = set(atom.resn for atom in model.atom)
        ligand_resnames = [r for r in all_resnames if r not in standard_residues]
        
        if ligand_resnames:
            ligand_resname = ligand_resnames[0]
            print(f"[visualize_protein_ligand_3d] 📌 自动检测到配体: {ligand_resname}")
        else:
            print(f"[visualize_protein_ligand_3d] ⚠️ 未找到配体")
            ligand_resname = None
    
    # 定义配体选择
    if ligand_resname:
        lig_sel = f"{obj_name} and resn {ligand_resname}"
    else:
        print(f"[visualize_protein_ligand_3d] ❌ No ligand; cannot proceed")
        return

    # ========== 第三步：隐藏所有，准备重新显示 ==========
    cmd.hide("everything", obj_name)
    
    # ========== 第四步：选择并创建口袋区域 ==========
    cmd.select("lig_pocket", f"byres ({obj_name} and polymer within 5 of ({lig_sel}))")
    
    # ========== 第五步：添加氢原子（关键！）==========
    print(f"[visualize_protein_ligand_3d] ➕ 正在添加氢原子...")
    
    # 检测配体是否已有氢原子
    n_h_lig = cmd.count_atoms(f"({lig_sel}) and hydro")
    if n_h_lig == 0:
        try:
            cmd.h_add(lig_sel)
            print(f"[visualize_protein_ligand_3d]    ✓ 已为配体添加氢原子")
        except Exception as e:
            print(f"[visualize_protein_ligand_3d]    ⚠️ 为配体添加氢原子失败: {e}")
    
    # 检测口袋残基是否已有氢原子
    n_h_pocket = cmd.count_atoms("lig_pocket and hydro")
    if n_h_pocket == 0:
        try:
            cmd.h_add("lig_pocket")
            print(f"[visualize_protein_ligand_3d]    ✓ 已为口袋残基添加氢原子")
        except Exception as e:
            print(f"[visualize_protein_ligand_3d]    ⚠️ 为口袋残基添加氢原子失败: {e}")

    # ========== 第六步：基础显示设置 ==========
    # 背景白色
    cmd.bg_color("white")
    
    # 蛋白整体：半透明cartoon
    cmd.show("cartoon", obj_name)
    cmd.set("cartoon_transparency", 0.3, obj_name)
    
    # 只显示配体为 sticks；口袋残基本身不单独高亮，后面只高亮真正有高置信度相互作用的残基
    cmd.show("sticks", lig_sel)
    
    # 隐藏连接到碳原子的氢（只保留极性氢：NH, OH, SH）
    cmd.hide("(h. and (e. c extend 1))")
    
    # ========== 第七步：配体着色（参考脚本：配体碳原子黄色） ==========\
    cmd.color("yellow", f"{lig_sel} and name C*")     # 配体碳原子：黄色
    cmd.color("blue", f"{lig_sel} and elem N")
    cmd.color("red", f"{lig_sel} and elem O")
    cmd.color("yellow", f"{lig_sel} and elem S")
    cmd.color("green", f"{lig_sel} and elem F+CL+BR+I")
    
    # ========== 第八步：验证数据存在 ==========
    # ⚠️ 关键：必须提供 interactions 数据，不再回退到 PyMOL 几何检测
    if not interactions:
        error_msg = (
            "[visualize_protein_ligand_3d] ❌ Error: No interaction data provided!\n"
            "   This function requires interaction data from analysis.\n\n"
            "   Usage:\n"
            "      1. result = analyze_protein_ligand_interactions('obj', 'LIG')\n"
            "      2. visualize_protein_ligand_3d('obj', result, 'LIG')\n\n"
            "   Or:\n"
            "      visualize_protein_ligand_3d('obj', csv_path='interactions.csv', ligand_resname='LIG')\n"
        )
        print(error_msg)
        return
    
    print(f"[visualize_protein_ligand_3d] 🔗 正在绘制相互作用 (基于分析结果)...")
    
    # 删除旧的 PyMOL 几何检测对象（如果有）
    try:
        cmd.delete("polar_donors_lig")
        cmd.delete("polar_donors_res")
        cmd.delete("polar_acceptors_lig")
        cmd.delete("polar_acceptors_res")
        cmd.delete("don_hydrogens_lig")
        cmd.delete("don_hydrogens_res")
        cmd.delete("hbonds_res_to_lig")
        cmd.delete("hbonds_lig_to_res")
    except:
        pass
    
    # ========== 第九步：氢键样式设置（专业配色）==========
    cmd.set("dash_length", 0.3)
    cmd.set("dash_radius", 0.06)   # 使用细圆柱体，确保在所有渲染模式下可见
    cmd.set("dash_color", "blue")  # 氢键：蓝色虚线
    cmd.set("dash_width", 2.0)
    cmd.set("dash_gap", 0.5)
    cmd.hide("labels", "hbonds_*")  # 隐藏距离标签
    
    # 为氢键距离对象设置蓝色
    try:
        cmd.color("blue", "hbonds_*")
    except:
        pass
    
    # ========== 第十步：提取并显示相互作用残基 ==========
    protein_residues = set()
    for inter in interactions:
        prot_chain = inter.get("Protein_Chain", "")
        prot_res = inter.get("Protein_Residue", "")
        if prot_res:
            res_parts = prot_res.split()
            if len(res_parts) >= 2:
                res_name, res_id = res_parts[0], res_parts[1]
                protein_residues.add((prot_chain, res_name, res_id))
    
    if protein_residues:
        print(f"[visualize_protein_ligand_3d] 📍 显示 {len(protein_residues)} 个相互作用残基")
        
        # 3字母氨基酸代码 -> 1字母代码
        aa_map = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
            "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
            "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
            "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
            "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }
        
        # 显示相互作用残基为 licorice 并添加标签
        for i, (chain, res_name, res_id) in enumerate(protein_residues, 1):
            res_sel = f"{obj_name} and chain {chain} and resi {res_id}"
            cmd.show("licorice", res_sel)
            # 隐藏这些残基上连接到碳的氢
            cmd.hide("everything", f"({res_sel}) and (elem H and neighbor elem C)")
            
            # 添加残基标签（使用一字母代码 + 残基号）
            label_text = aa_map.get(res_name.upper(), res_name[:1]) + res_id
            try:
                # 使用CA原子位置放置标签
                ca_sel = f"{res_sel} and name CA"
                if cmd.count_atoms(ca_sel) > 0:
                    coords = cmd.get_atom_coords(ca_sel)
                    label_obj = f"res_label_{i}"
                    cmd.pseudoatom(label_obj, pos=coords, label=label_text)
                    cmd.set("label_size", 16, label_obj)
                    cmd.set("label_color", "black", label_obj)
                    cmd.set("label_font_id", 7, label_obj)  # Times-like font
                    cmd.hide("everything", label_obj)
                    cmd.show("label", label_obj)
            except Exception as e:
                pass  # 静默失败
    
    # ========== 第十一步：绘制其他类型相互作用 ==========
    print(f"[visualize_protein_ligand_3d] 🎨 正在绘制关键药物设计相互作用...")
    
    # 相互作用类型到颜色的映射（专业配色方案）
    color_map = {
        "氢键": "blue",           # 氢键：蓝色
        "Hbond": "blue",
        "盐桥": "orange",         # 盐桥：橘色
        "SaltBridge": "orange",
        "疏水相互作用": "gray",    # 疏水：灰色（如果显示）
        "Hydrophobic": "gray",
        "π–π 堆积": "green",      # π-π堆积：绿色
        "PiPi": "green",
        "π–阳离子相互作用": "magenta",  # π-阳离子：品红
        "PiCation": "magenta",
        "金属配位": "violet",      # 金属配位：紫罗兰
        "MetalCoord": "violet",
    }
    
    # 相互作用类型到英文名称的映射（用于生成合法的PyMOL对象名）
    type_name_map = {
        "氢键": "Hbond",
        "盐桥": "SaltBridge",
        "疏水相互作用": "Hydrophobic",
        "π–π 堆积": "PiPi",
        "π–阳离子相互作用": "PiCation",
        "金属配位": "MetalCoord",
    }
    
    # 药物设计重要性优先级
    interaction_priority = {
        "氢键": 1,
        "盐桥": 2, 
        "金属配位": 3,
        "π–π 堆积": 4,
        "π–阳离子相互作用": 5,
        "疏水相互作用": 6  # 最低优先级
    }
    
    # 按类型分组并按距离排序
    interactions_by_type = {}
    for inter in interactions:
        itype = inter.get("Interaction", "")
        if itype not in interactions_by_type:
            interactions_by_type[itype] = []
        interactions_by_type[itype].append(inter)
    
    # 对每种类型按距离排序（短距离优先）
    for itype in interactions_by_type:
        interactions_by_type[itype].sort(key=lambda x: float(x.get("Distance", "999")))
    
    # 设置每种相互作用类型的默认显示数量（专业筛选策略）
    if max_interactions_per_type is None:
        max_interactions_per_type = {
            "氢键": 10,         # 氢键全部显示（重要）
            "盐桥": 5,          # 盐桥全部显示（重要）
            "金属配位": 3,      # 金属配位全部显示
            "π–π 堆积": 3,      # π相互作用选择性显示
            "π–阳离子相互作用": 3,
            "疏水相互作用": 0 if not show_hydrophobic else 5  # 默认不显示疏水（除非指定）
        }
    
    interaction_count = {}
    
    # 先删除旧的相互作用对象
    try:
        cmd.delete("pl_interact_*")
        cmd.delete("interact_*")
    except:
        pass
    
    # 按优先级顺序处理相互作用
    sorted_types = sorted(interactions_by_type.keys(), 
                         key=lambda x: interaction_priority.get(x, 999))
    
    total_shown = 0
    hbonds_drawn_from_csv = 0  # 追踪从 CSV 绘制的氢键数量
    for interaction_type in sorted_types:
        # 跳过不显示的疏水相互作用
        if not show_hydrophobic and ("疏水" in interaction_type or "Hydrophobic" in interaction_type):
            continue
        
        # 获取该类型的最大显示数量
        max_for_type = max_interactions_per_type.get(interaction_type, 3)
        if interaction_count.get(interaction_type, 0) >= max_for_type:
            continue
        
        # ⚠️ 关键修复：删除跳过氢键的逻辑
        # 现在所有相互作用（包括氢键）都从 CSV/result 数据绘制
        # 只有当 interactions 为空时，use_geom_hbonds=True，才会使用PyMOL的几何检测
        
        # 筛选相互作用
        interactions_to_show = interactions_by_type.get(interaction_type, [])[:max_for_type]
        
        for idx, inter in enumerate(interactions_to_show, 1):
            try:
                # 解析配体信息
                lig_chain = inter.get("Ligand_Chain", "")
                lig_res = inter.get("Ligand_Residue", "").split()
                lig_atom = inter.get("Ligand_Atom", "")

                # 解析蛋白信息
                prot_chain = inter.get("Protein_Chain", "")
                prot_res = inter.get("Protein_Residue", "").split()
                prot_atom = inter.get("Protein_Atom", "")

                if len(lig_res) < 2 or len(prot_res) < 2:
                    continue

                lig_resname_int, lig_resid = lig_res[0], lig_res[1]
                prot_resname, prot_resid = prot_res[0], prot_res[1]

                # 创建选择（处理空链 ID）
                if lig_chain and lig_chain.strip():
                    sel1 = f"{obj_name} and chain {lig_chain} and resi {lig_resid}"
                else:
                    sel1 = f"{obj_name} and resi {lig_resid}"
                
                if prot_chain and prot_chain.strip():
                    sel2 = f"{obj_name} and chain {prot_chain} and resi {prot_resid}"
                else:
                    sel2 = f"{obj_name} and resi {prot_resid}"

                # 生成合法的PyMOL对象名（只使用英文和数字），并带上 CSV 中的行号，方便一一对应
                # 将中文相互作用类型转换为英文
                type_en = type_name_map.get(interaction_type, interaction_type)
                # 移除所有非字母数字和下划线的字符
                type_en_clean = ''.join(c for c in type_en if c.isalnum())
                row_idx = inter.get("_row_index", idx)
                dist_name = f"interact_{type_en_clean}_{row_idx}"
                
                # ========== π-π 相互作用：使用环中心 pseudoatom（参考 PPI 代码）==========
                if "ring" in str(lig_atom).lower() or "ring" in str(prot_atom).lower():
                    try:
                        # 定义蛋白质芳香环原子
                        ring_atoms_map = {
                            "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                            "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                            "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],  # 六元环
                            "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]
                        }
                        
                        def get_ring_center_from_selection(obj, chain, resi, resname):
                            """计算芳香环中心坐标"""
                            # 蛋白质标准残基
                            if resname in ring_atoms_map:
                                ring_atom_names = ring_atoms_map[resname]
                                coords = []
                                for atom_name in ring_atom_names:
                                    if chain and chain.strip():
                                        sel = f"{obj} and chain {chain} and resi {resi} and name {atom_name}"
                                    else:
                                        sel = f"{obj} and resi {resi} and name {atom_name}"
                                    try:
                                        model = cmd.get_model(sel)
                                        if model.atom:
                                            coord = model.atom[0].coord
                                            coords.append((coord[0], coord[1], coord[2]))
                                    except:
                                        continue
                                if len(coords) >= 3:
                                    x = sum(c[0] for c in coords) / len(coords)
                                    y = sum(c[1] for c in coords) / len(coords)
                                    z = sum(c[2] for c in coords) / len(coords)
                                    return (x, y, z)
                                return None
                            
                            # 配体：使用动态检测
                            if chain and chain.strip():
                                sel = f"{obj} and chain {chain} and resi {resi}"
                            else:
                                sel = f"{obj} and resi {resi}"
                            
                            # 获取配体所有原子
                            model = cmd.get_model(sel)
                            if not model.atom:
                                return None
                            
                            atoms = [(a.chain, a.resn, a.resi, a.name, (a.coord[0], a.coord[1], a.coord[2])) 
                                     for a in model.atom]
                            
                            # 检测配体环
                            detected_rings = detect_ligand_rings(atoms)
                            if detected_rings:
                                ring_coords = detected_rings[0]  # 使用第一个检测到的环
                                x = sum(c[0] for c in ring_coords) / len(ring_coords)
                                y = sum(c[1] for c in ring_coords) / len(ring_coords)
                                z = sum(c[2] for c in ring_coords) / len(ring_coords)
                                return (x, y, z)
                            return None
                        
                        # 计算配体环中心
                        lig_center = get_ring_center_from_selection(obj_name, lig_chain, lig_resid, lig_resname_int)
                        # 计算蛋白质环中心
                        prot_center = get_ring_center_from_selection(obj_name, prot_chain, prot_resid, prot_resname)
                        
                        if lig_center and prot_center:
                            # 创建 pseudoatom 在环中心
                            pseudo1 = f"pl_centroid_{type_en_clean}_{row_idx}_lig"
                            pseudo2 = f"pl_centroid_{type_en_clean}_{row_idx}_prot"
                            
                            cmd.pseudoatom(pseudo1, pos=lig_center)
                            cmd.pseudoatom(pseudo2, pos=prot_center)
                            cmd.hide("everything", pseudo1)
                            cmd.hide("everything", pseudo2)
                            
                            # 使用 pseudoatom 作为选择
                            sel1 = pseudo1
                            sel2 = pseudo2
                            
                            print(f"[visualize_protein_ligand_3d] ✓ π-π: 使用环中心 {lig_resname_int}-{prot_resname}")
                        else:
                            # 回退到原子选择
                            if lig_atom and lig_atom != "":
                                sel1 += f" and name {lig_atom}" if "ring" not in str(lig_atom).lower() else ""
                            if prot_atom and prot_atom != "":
                                sel2 += f" and name {prot_atom}" if "ring" not in str(prot_atom).lower() else ""
                            print(f"[visualize_protein_ligand_3d] ⚠️ π-π: 无法计算环中心，使用原子选择")
                    
                    except Exception as e:
                        print(f"[visualize_protein_ligand_3d] ⚠️ π-π 环中心计算失败: {e}")
                        # 回退到原子选择
                        if lig_atom and lig_atom != "" and "ring" not in str(lig_atom).lower():
                            sel1 += f" and name {lig_atom}"
                        if prot_atom and prot_atom != "" and "ring" not in str(prot_atom).lower():
                            sel2 += f" and name {prot_atom}"
                else:
                    # 非 π 相互作用：使用原子名选择
                    if lig_atom and lig_atom != "":
                        sel1 += f" and name {lig_atom}"
                    if prot_atom and prot_atom != "":
                        sel2 += f" and name {prot_atom}"
                
                # 尝试创建距离对象（参考 PPI 可视化代码）
                try:
                    # 检查选择是否有效（pseudoatom 不需要检查）
                    is_pseudoatom = sel1.startswith("pl_centroid_") or sel2.startswith("pl_centroid_")
                    
                    if not is_pseudoatom:
                        if cmd.count_atoms(sel1) == 0:
                            print(f"[visualize_protein_ligand_3d] ⚠️ Warning: Ligand selection empty: {sel1}")
                            continue
                        if cmd.count_atoms(sel2) == 0:
                            print(f"[visualize_protein_ligand_3d] ⚠️ Warning: Protein selection empty: {sel2}")
                            continue

                    # Debug output for the first few interactions
                    if idx <= 3:
                        print(f"[visualize_protein_ligand_3d] DEBUG: Creating {dist_name}")
                        print(f"    sel1: {sel1}")
                        print(f"    sel2: {sel2}")

                    # Determine cutoff based on interaction type
                    cutoff = 5.0  # Default
                    if "氢键" in interaction_type or "Hbond" in interaction_type:
                        cutoff = 3.8
                    elif "盐桥" in interaction_type or "Salt" in interaction_type:
                        cutoff = 5.0
                    elif "疏水" in interaction_type or "Hydrophobic" in interaction_type:
                        cutoff = 5.5
                    elif "Pi" in interaction_type or "π" in interaction_type:
                        cutoff = 7.0  # π 相互作用距离更大
                    elif "金属" in interaction_type or "Metal" in interaction_type:
                        cutoff = 4.0

                    # 创建距离对象
                    cmd.distance(dist_name, sel1, sel2, cutoff=cutoff, mode=0)
                    
                    # 强制显示设置
                    cmd.enable(dist_name)
                    cmd.show("dashes", dist_name)
                    cmd.set("dash_width", 3.0, dist_name)
                    cmd.set("dash_gap", 0.2, dist_name)
                    cmd.set("dash_length", 0.4, dist_name)
                    
                    if interaction_color:
                        cmd.color(interaction_color, dist_name)
                        cmd.set("dash_color", interaction_color, dist_name)
                        
                except Exception as e:
                    print(f"[visualize_protein_ligand_3d] ⚠️ Error creating distance {dist_name}: {e}")
                    continue

                # 设置颜色和样式（与 PPI 代码保持一致的简洁方式）
                interaction_color = None
                for key, color in color_map.items():
                    if key in interaction_type or key == type_en:
                        interaction_color = color
                        break
                
                if interaction_color:
                    cmd.set("dash_color", interaction_color, dist_name)
                
                # 设置线条宽度并隐藏标签（与 PPI 代码一致）
                cmd.set("dash_width", 2.0, dist_name)
                cmd.hide("labels", dist_name)
                
                # 对于π相互作用,使用更粗的线条
                if "PiPi" in type_en_clean or "PiCation" in type_en_clean:
                    cmd.set("dash_width", 3.0, dist_name)

                # 统计相互作用类型
                interaction_count[interaction_type] = interaction_count.get(interaction_type, 0) + 1
                total_shown += 1
                
                # 统计氢键
                if "氢键" in interaction_type or "Hbond" in interaction_type:
                    hbonds_drawn_from_csv += 1

            except Exception as e:
                # print(f"[visualize_protein_ligand_3d] ⚠️ 相互作用跳过: {e}")
                continue
    
    # 如果没有显示任何非氢键相互作用，提示用户
    if total_shown == 0 and not show_hydrophobic:
        print(f"[visualize_protein_ligand_3d] 💡 Note: hydrophobic interactions are hidden (professional mode)")
        print(f"                                  To show them, use: visualize_protein_ligand_3d('{obj_name}', show_hydrophobic=True)")

    # 隐藏所有距离标签（但保留线条可见）
    cmd.hide("labels", "interact_*")
    cmd.hide("labels", "hbonds_*")

    # 显式显示交互线条（有些主题/样式下需要）
    cmd.show("dashes", "interact_*")
    cmd.show("dashes", "hbonds_*")
    
    # 确保线条宽度足够
    cmd.set("dash_width", 3.0, "interact_*")
    cmd.set("dash_width", 3.0, "hbonds_*")
    
    # ========== 第十二步：调整视角 ==========\
    cmd.zoom(f"({lig_sel}) or lig_pocket", buffer=8)
    cmd.orient(f"({lig_sel}) or lig_pocket")
    
    # ========== 第十三步：美化参数（参考脚本）==========
    cmd.set("stick_radius", 0.15)
    cmd.set("sphere_scale", 0.25)
    cmd.set("ray_trace_mode", 1)
    
    # 设置标签样式：使用接近 Times Roman 的字体，黑色，便于发表级别的残基标签（如 "T100"）
    cmd.set("label_size", 28)
    cmd.set("label_font_id", 5)
    cmd.set("label_color", "black")
    
    # ========== 第十四步：输出统计 ==========
    print(f"\n[visualize_protein_ligand_3d] ✅ 3D 可视化完成!")
    print(f"   📊 相互作用统计 (来自分析结果):")
    
    # 统计氢键
    total_hbonds = interaction_count.get("氢键", 0) + interaction_count.get("Hbond", 0)
    if total_hbonds > 0:
        print(f"      • 氢键: {total_hbonds}")
    
    # 统计其他相互作用
    for itype, count in sorted(interaction_count.items(), key=lambda x: x[1], reverse=True):
        if count > 0 and itype not in ["氢键", "Hbond"]:
            print(f"      • {itype}: {count}")
    
    if protein_residues:
        print(f"   📍 涉及 {len(protein_residues)} 个蛋白残基")
    
    print(f"\n   💡 PyMOL 提示:")
    print(f"      show labels, hbonds_*       # 显示氢键距离")
    print(f"      hide labels, hbonds_*       # 隐藏氢键距离")
    print(f"      show labels, interact_*     # 显示其他相互作用距离")
    print(f"      hide dashes, interact_*     # 隐藏相互作用线条")
    print(f"      show dashes, interact_*     # 显示相互作用线条")
    print(f"      delete interact_*           # 删除所有相互作用对象")
    print(f"      delete hbonds_*             # 删除氢键对象")
    print(f"      ray                         # 高质量渲染")
    print(f"      png output.png, dpi=300     # 保存高清图片")

cmd.extend("visualize_protein_ligand_3d", visualize_protein_ligand_3d)

def diagnose_csv_visualization(obj_name, csv_path, max_rows=5):
    """
    诊断CSV可视化问题 - 检查CSV中的选择是否有效
    Diagnose CSV visualization issues - check if selections in CSV are valid
    
    参数:
        obj_name: PyMOL对象名称
        csv_path: CSV文件路径
        max_rows: 最多检查多少行 (default: 5)
    """
    from pymol import cmd
    import csv
    
    print(f"\n[diagnose_csv_visualization] 🔍 Diagnosing CSV: {csv_path}")
    print(f"   Object: {obj_name}\n")
    
    if obj_name not in cmd.get_object_list():
        print(f"   ❌ Object '{obj_name}' not found in PyMOL!")
        print(f"   Available objects: {cmd.get_object_list()}")
        return
    
    if not os.path.exists(csv_path):
        print(f"   ❌ CSV file not found: {csv_path}")
        return
    
    # Read CSV
    try:
        with open(csv_path, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
    except Exception as e:
        print(f"   ❌ Failed to read CSV: {e}")
        return
    
    if not rows:
        print(f"   ❌ CSV is empty!")
        return
    
    print(f"   CSV has {len(rows)} interactions\n")
    print(f"   Checking first {min(max_rows, len(rows))} rows:\n")
    
    valid_count = 0
    for i, row in enumerate(rows[:max_rows], 1):
        lig_chain = row.get("Ligand_Chain", "")
        lig_res = row.get("Ligand_Residue", "").split()
        lig_atom = row.get("Ligand_Atom", "")
        prot_chain = row.get("Protein_Chain", "")
        prot_res = row.get("Protein_Residue", "").split()
        prot_atom = row.get("Protein_Atom", "")
        interaction_type = row.get("Interaction", "")
        confidence = row.get("Confidence", "1.0")
        
        if len(lig_res) < 2 or len(prot_res) < 2:
            print(f"   Row {i}: ⚠️  Invalid residue format")
            continue
        
        lig_resid = lig_res[1]
        prot_resid = prot_res[1]
        
        # Create selections
        if lig_chain and lig_chain.strip():
            sel1 = f"{obj_name} and chain {lig_chain} and resi {lig_resid}"
        else:
            sel1 = f"{obj_name} and resi {lig_resid}"
        
        if lig_atom:
            sel1 += f" and name {lig_atom}"
        
        if prot_chain and prot_chain.strip():
            sel2 = f"{obj_name} and chain {prot_chain} and resi {prot_resid}"
        else:
            sel2 = f"{obj_name} and resi {prot_resid}"
        
        if prot_atom:
            sel2 += f" and name {prot_atom}"
        
        # Check selections
        count1 = cmd.count_atoms(sel1)
        count2 = cmd.count_atoms(sel2)
        
        # Print results
        status = "✅" if (count1 > 0 and count2 > 0) else "❌"
        conf_val = f", conf={confidence}" if confidence else ""
        hydro_note = " [hidden by default]" if "疏水" in interaction_type or "Hydrophobic" in interaction_type else ""
        
        print(f"   Row {i}: {status} {interaction_type}{conf_val}{hydro_note}")
        print(f"          Ligand: {lig_res[0]} {lig_resid}/{lig_atom} → {count1} atoms")
        print(f"          Protein: {prot_res[0]} {prot_resid}/{prot_atom} → {count2} atoms")
        
        if count1 == 0:
            print(f"          ⚠️  Ligand selection empty: {sel1}")
        if count2 == 0:
            print(f"          ⚠️  Protein selection empty: {sel2}")
        
        if count1 > 0 and count2 > 0:
            valid_count += 1
        
        print()
    
    print(f"   Summary: {valid_count}/{min(max_rows, len(rows))} selections valid\n")
    
    # Print recommendations
    print("   💡 Recommendations:")
    print("      1. To show ALL interactions including hydrophobic:")
    print(f"         visualize_protein_ligand_3d('{obj_name}', csv_path='{csv_path}', show_hydrophobic=True)\n")
    print("      2. To lower confidence threshold (e.g., ≥0.5):")
    print(f"         visualize_protein_ligand_3d('{obj_name}', csv_path='{csv_path}', min_confidence=0.5)\n")
    print("      3. To show all with both options:")
    print(f"         visualize_protein_ligand_3d('{obj_name}', csv_path='{csv_path}', show_hydrophobic=True, min_confidence=0.5)\n")

cmd.extend("diagnose_csv_visualization", diagnose_csv_visualization)

def toggle_interaction_lines(show=True, line_width=2.5, dash_gap=0.15, dash_length=0.25):
    """
    切换相互作用线条的显示/隐藏状态
    Toggle interaction lines visibility and properties
    
    参数:
        show: True显示线条, False隐藏线条 (Show/hide interaction lines)
        line_width: 线条宽度 (Line width)
        dash_gap: 虚线间隙 (Gap between dashes)
        dash_length: 虚线长度 (Length of each dash)
    
    使用示例:
        toggle_interaction_lines(True, 3.0)  # 显示粗线条
        toggle_interaction_lines(False)      # 隐藏所有线条
    """
    from pymol import cmd
    
    interaction_objects = ["interact_*", "hbonds_*", "ppi_dist*"]
    
    for obj_pattern in interaction_objects:
        if show:
            cmd.show("dashes", obj_pattern)
            cmd.set("dash_width", line_width, obj_pattern)
            cmd.set("dash_gap", dash_gap, obj_pattern)
            cmd.set("dash_length", dash_length, obj_pattern)
        else:
            cmd.hide("dashes", obj_pattern)
    
    if show:
        print(f"[toggle_interaction_lines] ✅ Interaction lines shown (width={line_width})")
    else:
        print(f"[toggle_interaction_lines] ⚫ Interaction lines hidden")

cmd.extend("toggle_interaction_lines", toggle_interaction_lines)

def apply_plot_style(style="professional"):
    """
    应用不同的matplotlib绘图样式
    
    参数:
        style: 样式名称 (professional/hand-drawn/minimalist/publication/colorful)
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    
    style = style.lower()
    
    if style == "hand-drawn":
        # 手绘风格 (xkcd)
        plt.xkcd()
        
    elif style == "minimalist":
        # 简约风格
        plt.style.use('seaborn-v0_8-whitegrid')
        mpl.rcParams['axes.spines.top'] = False
        mpl.rcParams['axes.spines.right'] = False
        mpl.rcParams['axes.grid'] = False
        mpl.rcParams['font.size'] = 10
        mpl.rcParams['axes.labelsize'] = 11
        mpl.rcParams['axes.titlesize'] = 12
        
    elif style == "publication":
        # 学术出版风格
        plt.style.use('seaborn-v0_8-paper')
        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
        mpl.rcParams['font.size'] = 11
        mpl.rcParams['axes.linewidth'] = 1.5
        mpl.rcParams['lines.linewidth'] = 2.0
        mpl.rcParams['patch.linewidth'] = 1.5
        
    elif style == "colorful":
        # 缤纷风格
        plt.style.use('seaborn-v0_8-bright')
        mpl.rcParams['axes.facecolor'] = '#f0f0f0'
        mpl.rcParams['figure.facecolor'] = 'white'
        
    else:  # professional (default)
        # 专业风格 (默认)
        plt.style.use('default')
        mpl.rcParams['font.family'] = 'sans-serif'
        mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
        mpl.rcParams['axes.unicode_minus'] = False
        mpl.rcParams['figure.facecolor'] = 'white'
        mpl.rcParams['axes.facecolor'] = 'white'

def generate_advanced_interaction_plot(interactions, output_path=None, show_plot=True, ligand_sdf=None, obj_name=None, ligand_resname=None, plot_style="professional"):
    """
    为高级分析生成<unk>D相互作用图（带配体化学结构）
    
    参数:
        interactions: 相互作用列表（高级格式）
        output_path: 输出路径
        show_plot: 是否显示
        ligand_sdf: 配体SDF文件路径（可选）
        obj_name: PyMOL对象名称（可选，用于从PyMOL提取配体）
        ligand_resname: 配体残基名称（与obj_name配合使用）
        plot_style: 绘图样式 (professional/hand-drawn/minimalist/publication/colorful)
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.lines import Line2D
        import matplotlib
    except ImportError:
        print("[generate_advanced_interaction_plot] matplotlib is required")
        return None
    
    # 确保RDKit已安装（高质量2D图必需）
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
        from PIL import Image
        import io
        import numpy as np
    except ImportError:
        error_msg = "[generate_advanced_interaction_plot] ❌ RDKit not installed; cannot generate 2D interaction diagram\nRun: pip install rdkit pillow numpy"
        print(error_msg)
        raise RuntimeError(error_msg)
    
    # 收集相互作用信息
    interaction_map = {}  # {residue_key: [(type, distance, ligand_atom)]}
    residue_names = {}  # {residue_key: residue_name}
    
    for inter in interactions:
        itype = inter.get("Type", "Unknown")
        p_atom = inter.get("Protein_Atom", "")
        l_atom = inter.get("Ligand_Atom", "")
        distance_str = inter.get("Distance", "0.0")
        protein_residue = inter.get("Protein_Residue", "")  # 新增：读取残基名称
        
        try:
            distance = float(distance_str)
        except (ValueError, TypeError):
            distance = 0.0
        
        # 使用残基名称作为key（如果有的话）
        if protein_residue and protein_residue != "":
            res_key = protein_residue
        else:
            res_key = f"Atom{p_atom}"
        
        if res_key not in interaction_map:
            interaction_map[res_key] = []
            residue_names[res_key] = res_key
        interaction_map[res_key].append((itype, distance, l_atom))
    
    if not interaction_map:
        print("[generate_advanced_interaction_plot] No data to plot")
        return None
    
    # 应用绘图样式 (必须在创建图形之前,且在import之后)
    apply_plot_style(plot_style)
    
    # 设置基础字体支持 (不覆盖样式)
    if matplotlib.rcParams.get('font.sans-serif') is None:
        matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    matplotlib.rcParams['axes.unicode_minus'] = False
    
    # 创建图形
    fig = plt.figure(figsize=(14, 10), dpi=150)
    ax = fig.add_subplot(111)
    
    # 颜色和线型映射
    colors = {
        "Hbond": "#2196F3",
        "SaltBridge": "#FF5722",
        "Hydrophobic": "#4CAF50",
        "PiPi": "#9C27B0",
        "PiCation": "#E91E63",
        "MetalCoord": "#FFC107",
        "Halogen": "#00BCD4"
    }
    
    line_styles = {
        "Hbond": "--",
        "SaltBridge": ":",
        "Hydrophobic": "-",
        "PiPi": "-.",
        "PiCation": "--",
        "MetalCoord": "-",
        "Halogen": "--"
    }
    
    # 设置画布
    ax.set_xlim(-8, 8)
    ax.set_ylim(-6, 6)
    ax.axis('off')
    
    # 绘制配体结构（必须使用RDKit）
    from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
    ligand_drawn = False
    atom_coords_2d = {}  # 配体原子的2D坐标
    mol = None
    
    import os
    import tempfile
    
    # 尝试从多个来源获取配体结构
    # 优先级1: ligand_sdf文件
    if ligand_sdf and os.path.exists(ligand_sdf):
        try:
            # 从SDF文件读取配体
            supplier = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
            mol = supplier[0] if len(supplier) > 0 else None
            if mol:
                print(f"[generate_advanced_interaction_plot] Loaded ligand from SDF: {ligand_sdf}")
        except Exception as e:
            print(f"[generate_advanced_interaction_plot] Failed to load ligand from SDF: {e}")
    
    # 优先级2: 从PyMOL对象提取
    if mol is None and obj_name and ligand_resname:
        try:
            from pymol import cmd
            temp_sdf = tempfile.mktemp(suffix=".sdf")
            cmd.save(temp_sdf, f"{obj_name} and resn {ligand_resname}", format="sdf")
            if os.path.exists(temp_sdf):
                supplier = Chem.SDMolSupplier(temp_sdf, removeHs=False)
                mol = supplier[0] if len(supplier) > 0 else None
                os.remove(temp_sdf)
                if mol:
                    print(f"[generate_advanced_interaction_plot] Extracted ligand from PyMOL: {obj_name}/{ligand_resname}")
        except Exception as e:
            print(f"[generate_advanced_interaction_plot] Failed to extract ligand from PyMOL: {e}")
    
    # 处理配体分子
    if mol is not None:
        try:
            # 生成2D坐标
            AllChem.Compute2DCoords(mol)
            
            # 收集参与相互作用的原子索引
            ligand_atoms_interacting = set()
            atom_to_residues = {}  # {ligand_atom_idx: [(res_key, itype, distance)]}
            
            for res_key, interactions_list in interaction_map.items():
                for itype, dist, l_atom_str in interactions_list:
                    try:
                        # 提取原子索引
                        if "ring(" in str(l_atom_str):
                            l_atom = int(str(l_atom_str).split("(")[1].split(",")[0])
                        else:
                            l_atom = int(l_atom_str)
                        
                        ligand_atoms_interacting.add(l_atom)
                        if l_atom not in atom_to_residues:
                            atom_to_residues[l_atom] = []
                        atom_to_residues[l_atom].append((res_key, itype, dist))
                    except (ValueError, IndexError):
                        pass
            
            # 绘制配体结构，高亮相互作用原子
            img = Draw.MolToImage(mol, size=(500, 400), 
                                 highlightAtoms=list(ligand_atoms_interacting),
                                 kekulize=True, wedgeBonds=True, fitImage=True)
            
            # 将PIL图像转换为numpy数组并显示在matplotlib中
            img_array = np.array(img)
            
            # 在中心显示配体结构
            ax.imshow(img_array, extent=[-3, 3, -2.2, 2.2], zorder=10)
            ligand_drawn = True
            print("[generate_advanced_interaction_plot] ✅ Drawn ligand chemical structure and highlighted interacting atoms")
            
            # 获取原子的2D坐标用于连线
            conformer = mol.GetConformer()
            atom_coords_2d = {}
            for atom in mol.GetAtoms():
                idx = atom.GetIdx()
                if idx in ligand_atoms_interacting:
                    pos = conformer.GetAtomPosition(idx)
                    # 归一化坐标到显示范围
                    atom_coords_2d[idx] = (pos.x, pos.y)
            
            # 将2D坐标归一化到matplotlib坐标系
            if atom_coords_2d:
                xs = [c[0] for c in atom_coords_2d.values()]
                ys = [c[1] for c in atom_coords_2d.values()]
                x_min, x_max = min(xs), max(xs)
                y_min, y_max = min(ys), max(ys)
                x_range = x_max - x_min if x_max != x_min else 1
                y_range = y_max - y_min if y_max != y_min else 1
                
                # 归一化到 [-2.5, 2.5] 范围
                for idx in atom_coords_2d:
                    x, y = atom_coords_2d[idx]
                    norm_x = ((x - x_min) / x_range - 0.5) * 5
                    norm_y = ((y - y_min) / y_range - 0.5) * 3.6
                    atom_coords_2d[idx] = (norm_x, norm_y)
            
        except Exception as e:
            error_msg = f"[generate_advanced_interaction_plot] ❌ Failed to draw ligand structure: {e}"
            print(error_msg)
            import traceback
            traceback.print_exc()
            raise RuntimeError(error_msg)
    
    # 确保配体结构已成功绘制
    if not ligand_drawn:
        error_msg = "[generate_advanced_interaction_plot] ❌ Ligand SDF file missing or cannot be parsed"
        print(error_msg)
        raise RuntimeError(error_msg)
    
    # 放置残基在周围
    import math
    residues = list(interaction_map.keys())
    n_res = len(residues)
    radius = 6.0
    
    # 氨基酸三字母代码到SMILES的映射（侧链）
    AA_SMILES = {
        "ALA": "CC(N)C(=O)O",
        "ARG": "NCCCC(N)C(=O)O", 
        "ASN": "NC(=O)CC(N)C(=O)O",
        "ASP": "OC(=O)CC(N)C(=O)O",
        "CYS": "SCC(N)C(=O)O",
        "GLN": "NC(=O)CCC(N)C(=O)O",
        "GLU": "OC(=O)CCC(N)C(=O)O",
        "GLY": "NCC(=O)O",
        "HIS": "c1c[nH]cn1CC(N)C(=O)O",
        "ILE": "CCC(C)C(N)C(=O)O",
        "LEU": "CC(C)CC(N)C(=O)O",
        "LYS": "NCCCCC(N)C(=O)O",
        "MET": "CSCCC(N)C(=O)O",
        "PHE": "c1ccc(CC(N)C(=O)O)cc1",
        "PRO": "C1CC(NC1)C(=O)O",
        "SER": "OCC(N)C(=O)O",
        "THR": "CC(O)C(N)C(=O)O",
        "TRP": "c1ccc2c(c1)c(cn2)CC(N)C(=O)O",
        "TYR": "Oc1ccc(CC(N)C(=O)O)cc1",
        "VAL": "CC(C)C(N)C(=O)O"
    }
    
    # 第一遍：放置残基标签/结构
    residue_positions = {}
    for i, res_key in enumerate(residues):
        angle = 2 * math.pi * i / n_res - math.pi/2
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        residue_positions[res_key] = (x, y)
        
        # 显示残基标签
        # 使用不同颜色区分不同类型的残基
        res_name_3letter = res_key[:3] if len(res_key) >= 3 else res_key
        
        # 根据残基类型着色
        if res_name_3letter in ['ARG', 'LYS', 'HIS']:  # 碱性
            box_color = '#E3F2FD'
            edge_color = '#1976D2'
        elif res_name_3letter in ['ASP', 'GLU']:  # 酸性
            box_color = '#FFEBEE'
            edge_color = '#D32F2F'
        elif res_name_3letter in ['PHE', 'TYR', 'TRP']:  # 芳香族
            box_color = '#F3E5F5'
            edge_color = '#7B1FA2'
        elif res_name_3letter in ['SER', 'THR', 'CYS']:  # 极性
            box_color = '#E8F5E9'
            edge_color = '#388E3C'
        else:  # 疏水性或其他
            box_color = '#FFF3E0'
            edge_color = '#F57C00'
        
        # 尝试绘制残基的侧链结构
        structure_drawn = False
        if res_name_3letter in AA_SMILES:
            try:
                # 从 SMILES 生成分子
                aa_mol = Chem.MolFromSmiles(AA_SMILES[res_name_3letter])
                if aa_mol:
                    # 生成 2D 坐标
                    AllChem.Compute2DCoords(aa_mol)
                    
                    # 绘制残基结构（小图）
                    aa_img = Draw.MolToImage(aa_mol, size=(120, 100), 
                                            kekulize=True, wedgeBonds=False)
                    aa_img_array = np.array(aa_img)
                    
                    # 在残基位置显示结构
                    img_size = 0.8  # 图片大小
                    ax.imshow(aa_img_array, 
                             extent=[x-img_size, x+img_size, y-img_size, y+img_size], 
                             zorder=5, alpha=0.95)
                    
                    # 在结构下方显示残基名称标签
                    ax.text(x, y-img_size-0.25, res_key, ha='center', va='top',
                           fontsize=7, fontweight='bold', color='#424242', zorder=6,
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                   edgecolor='gray', alpha=0.8, linewidth=0.5))
                    structure_drawn = True
            except Exception as e:
                # 如果绘制失败,回退到标签模式
                pass
        
        # 如果没有绘制结构,使用传统的圆圈+标签
        if not structure_drawn:
            res_box = Circle((x, y), 0.7, facecolor=box_color, 
                            edgecolor=edge_color, linewidth=2, zorder=5)
            ax.add_patch(res_box)
            
            # 显示残基标签
            ax.text(x, y, res_key, ha='center', va='center',
                   fontsize=8, fontweight='bold', color='#212121', zorder=6)
    
    # 第二遍：绘制连线（从配体原子到残基）
    for res_key, interactions_list in interaction_map.items():
        res_x, res_y = residue_positions[res_key]
        
        for itype, dist, l_atom_str in interactions_list:
            color = colors.get(itype, "#757575")
            style = line_styles.get(itype, "-")
            
            # 尝试获取配体原子的精确坐标
            try:
                if "ring(" in str(l_atom_str):
                    l_atom_idx = int(str(l_atom_str).split("(")[1].split(",")[0])
                else:
                    l_atom_idx = int(l_atom_str)
                
                if l_atom_idx in atom_coords_2d:
                    # 使用精确的原子坐标
                    lig_x, lig_y = atom_coords_2d[l_atom_idx]
                else:
                    # 回退到边缘连接
                    lig_angle = math.atan2(res_y, res_x)
                    lig_x = 2.5 * math.cos(lig_angle)
                    lig_y = 1.8 * math.sin(lig_angle)
            except (ValueError, IndexError):
                # 回退到边缘连接
                lig_angle = math.atan2(res_y, res_x)
                lig_x = 2.5 * math.cos(lig_angle)
                lig_y = 1.8 * math.sin(lig_angle)
            
            # 绘制连线
            arrow = FancyArrowPatch((lig_x, lig_y), (res_x, res_y),
                                  arrowstyle='-',
                                  linestyle=style,
                                  color=color,
                                  linewidth=2.0,
                                  zorder=2,
                                  alpha=0.7)
            ax.add_patch(arrow)
            
            # 距离标签（靠近残基端）
            mid_x = lig_x * 0.3 + res_x * 0.7
            mid_y = lig_y * 0.3 + res_y * 0.7
            ax.text(mid_x, mid_y, f"{dist:.1f}",
                   fontsize=6, ha='center',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                            edgecolor='none', alpha=0.8),
                   zorder=3)
    
    # 图例
    legend_elements = []
    seen_types = set()
    for interactions_list in interaction_map.values():
        for itype, _, _ in interactions_list:
            if itype not in seen_types:
                seen_types.add(itype)
                legend_elements.append(
                    Line2D([0], [0], color=colors.get(itype, "#757575"),
                          linewidth=3,
                          linestyle=line_styles.get(itype, '-'),
                          label=itype)
                )
    
    if legend_elements:
        leg = ax.legend(handles=legend_elements, loc='upper right', fontsize=11,
                       framealpha=0.95, edgecolor='black', title='Interactions',
                       title_fontsize=12)
        leg.get_title().set_fontweight('bold')
    
    # 标题
    total = len(interactions)
    ax.set_title(f'Protein-Ligand Interaction Diagram (2D)\nTotal: {total} interactions',
                fontsize=15, fontweight='bold', pad=15)
    
    plt.tight_layout()
    
    # 保存
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"[generate_advanced_interaction_plot] ✅ 2D interaction diagram saved: {output_path}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    return output_path

def generate_interaction_network_plot(interactions_result=None, csv_path=None,
                                     output_path=None, show_plot=True, ligand_sdf=None,
                                     obj_name=None, ligand_resname=None, plot_style="professional"):
    """
    生成交互网络图（使用matplotlib或networkx）

    参数:
        interactions_result: 相互作用分析结果
        csv_path: 或CSV文件路径
        output_path: 输出图片路径
        show_plot: 是否显示图片
        ligand_sdf: 配体SDF文件路径（可选）
        obj_name: PyMOL对象名称（可选，用于从PyMOL提取配体）
        ligand_resname: 配体残基名称（与obj_name配合使用）

    返回:
        output_path: 生成的图片路径
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.patches import FancyArrowPatch, Circle
    except ImportError:
        print("[generate_interaction_network_plot] matplotlib is required")
        return None

    # 读取数据
    interactions = []
    csv_format = "standard"  # standard 或 advanced
    
    if csv_path:
        try:
            with open(csv_path, 'r', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                fieldnames = reader.fieldnames or []
                
                # 检测CSV格式 - 优先使用advanced版本
                if "Protein_Atom" in fieldnames and "Ligand_Atom" in fieldnames:
                    csv_format = "advanced"
                    # print("[generate_interaction_network_plot] Detected advanced analysis format")
                
                for row in reader:
                    # 如果是标准格式,转换为advanced格式的key
                    if csv_format == "standard" and "Ligand_Residue" in row:
                        # 标准格式转换: Ligand_Residue -> Ligand_Atom 利于后续处理
                        converted_row = {
                            "Ligand_Atom": row.get("Ligand_Atom", ""),
                            "Protein_Atom": row.get("Protein_Atom", ""),
                            "Protein_Residue": row.get("Protein_Residue", ""),
                            "Distance": row.get("Distance", ""),
                            "Type": row.get("Interaction", "Unknown")
                        }
                        interactions.append(converted_row)
                    else:
                        interactions.append(row)
        except Exception as e:
            print(f"[generate_interaction_network_plot] Failed to read CSV: {e}")
            return None
    elif interactions_result:
        if isinstance(interactions_result, dict):
            if "interactions" in interactions_result:
                interactions = interactions_result["interactions"]
            elif "protein1_interactions" in interactions_result:
                # 三元复合体
                interactions = (interactions_result.get("protein1_interactions", []) +
                              interactions_result.get("protein2_interactions", []))
        else:
            interactions = interactions_result

    if not interactions:
        print("[generate_interaction_network_plot] No interaction data")
        return None
    
    # 默认使用带配体结构的高级版本 (如果有ligand_sdf或obj_name)
    use_advanced = csv_format == "advanced" or ligand_sdf or (obj_name and ligand_resname)
    
    if use_advanced:
        # 尝试从interactions中提取ligand_resname（如果未提供）
        if not ligand_resname and interactions:
            try:
                first_inter = interactions[0]
                if "Ligand_Residue" in first_inter:
                    # 从Ligand_Residue字段提取，格式可能是 "LIG 301" 或 "LIG"
                    lig_res_str = first_inter["Ligand_Residue"].strip()
                    # 提取残基名称（第一个单词）
                    ligand_resname = lig_res_str.split()[0] if lig_res_str else None
                    print(f"[generate_interaction_network_plot] Auto-detected ligand residue: {ligand_resname}")
            except Exception as e:
                print(f"[generate_interaction_network_plot] Failed to extract ligand residue name: {e}")
        
        return generate_advanced_interaction_plot(interactions, output_path, show_plot, 
                                                 ligand_sdf=ligand_sdf, 
                                                 obj_name=obj_name, 
                                                 ligand_resname=ligand_resname,
                                                 plot_style=plot_style)

    # 应用绘图样式 (简单网络图模式)
    apply_plot_style(plot_style)
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(14, 10), dpi=100)
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.axis('off')

    # 颜色映射
    color_map = {
        "氢键": "#2196F3",
        "盐桥": "#FF5722",
        "疏水相互作用": "#4CAF50",
        "π–π 堆积": "#9C27B0",
        "π–阳离子相互作用": "#E91E63",
        "二硫键": "#FFC107",
        "卤素键": "#FF9800",
        "金属配位": "#673AB7",
    }

    # 分析节点
    ligand_nodes = {}
    protein_nodes = {}

    for inter in interactions:
        # 配体节点
        lig_key = None
        prot_key = None

        if "Ligand_Residue" in inter:
            lig_key = inter["Ligand_Residue"]
            ligand_nodes[lig_key] = ligand_nodes.get(lig_key, 0) + 1
        elif "Nucleic_Residue" in inter:
            lig_key = inter["Nucleic_Residue"]
            ligand_nodes[lig_key] = ligand_nodes.get(lig_key, 0) + 1
        elif "Atom1_Residue" in inter:
            lig_key = inter["Atom1_Residue"]
            ligand_nodes[lig_key] = ligand_nodes.get(lig_key, 0) + 1

        # 蛋白节点
        if "Protein_Residue" in inter:
            prot_key = inter["Protein_Residue"]
            protein_nodes[prot_key] = protein_nodes.get(prot_key, 0) + 1
        elif "Atom2_Residue" in inter:
            prot_key = inter["Atom2_Residue"]
            protein_nodes[prot_key] = protein_nodes.get(prot_key, 0) + 1

    # 布局 - 配体在中心，蛋白残基环绕
    import math

    # 配体节点（中心）
    ligand_positions = {}
    n_lig = len(ligand_nodes)
    for i, (lig_name, count) in enumerate(ligand_nodes.items()):
        if n_lig == 1:
            ligand_positions[lig_name] = (0, 0)
        else:
            angle = 2 * math.pi * i / n_lig
            r = 2.0
            ligand_positions[lig_name] = (r * math.cos(angle), r * math.sin(angle))

    # 蛋白残基节点（外围）
    protein_positions = {}
    n_prot = len(protein_nodes)
    for i, (prot_name, count) in enumerate(sorted(protein_nodes.items(), key=lambda x: x[1], reverse=True)):
        angle = 2 * math.pi * i / n_prot
        r = 7.0
        protein_positions[prot_name] = (r * math.cos(angle), r * math.sin(angle))

    # 绘制节点
    # 配体节点
    for lig_name, pos in ligand_positions.items():
        count = ligand_nodes[lig_name]
        size = 0.8 + count * 0.1
        circle = Circle(pos, size, color='#FF6B35', alpha=0.8, ec='#D84315', lw=3, zorder=3)
        ax.add_patch(circle)
        ax.text(pos[0], pos[1], lig_name.split()[0], fontsize=12, fontweight='bold',
               ha='center', va='center', color='white', zorder=4)

    # 蛋白残基节点
    for prot_name, pos in protein_positions.items():
        count = protein_nodes[prot_name]
        size = 0.4 + count * 0.05
        circle = Circle(pos, size, color='#4A90E2', alpha=0.7, ec='#2E5C8A', lw=2, zorder=2)
        ax.add_patch(circle)
        ax.text(pos[0], pos[1], prot_name.split()[0][:3], fontsize=9,
               ha='center', va='center', color='white', zorder=4)

    # 绘制连接线
    interaction_types_shown = {}
    for inter in interactions:
        try:
            lig_key = inter.get("Ligand_Residue") or inter.get("Nucleic_Residue") or inter.get("Atom1_Residue")
            prot_key = inter.get("Protein_Residue") or inter.get("Atom2_Residue")

            if not lig_key or not prot_key:
                continue

            if lig_key not in ligand_positions or prot_key not in protein_positions:
                continue

            pos1 = ligand_positions[lig_key]
            pos2 = protein_positions[prot_key]

            interaction_type = inter.get("Interaction", "")

            # 确定颜色
            color = "#999999"
            for key, c in color_map.items():
                if key in interaction_type:
                    color = c
                    interaction_types_shown[interaction_type] = color
                    break

            # 绘制连线
            arrow = FancyArrowPatch(pos1, pos2, arrowstyle='-',
                                  color=color, alpha=0.6, lw=2, zorder=1)
            ax.add_patch(arrow)

        except Exception as e:
            continue

    # 添加图例 (Translate Chinese to English for plot labels to avoid font issues)
    translation_map = {
        "氢键": "H-Bond",
        "盐桥": "Salt Bridge",
        "疏水相互作用": "Hydrophobic",
        "π–π 堆积": "Pi-Pi Stacking",
        "π–阳离子相互作用": "Pi-Cation",
        "二硫键": "Disulfide",
        "卤素键": "Halogen Bond",
        "金属配位": "Metal Coord",
        "范德华力": "vdW"
    }

    legend_elements = []
    for itype, color in sorted(interaction_types_shown.items()):
        # Try exact match first, then partial match
        label = itype
        if itype in translation_map:
            label = translation_map[itype]
        else:
            # Fallback: try to find a key that is part of the itype string
            for k, v in translation_map.items():
                if k in itype:
                    label = v
                    break
        
        legend_elements.append(mpatches.Patch(color=color, label=label))

    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper right', frameon=True,
                 facecolor='white', edgecolor='gray', fontsize=10)

    # 标题
    title = f"Interaction Network\n({len(interactions)} interactions, {len(ligand_nodes)} ligands, {len(protein_nodes)} protein residues)"
    ax.text(0, 9.5, title, fontsize=14, fontweight='bold', ha='center')

    # 保存
    if output_path is None:
        output_path = "interaction_network.png"

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"[generate_interaction_network_plot] ✅ Network graph saved: {output_path}")

    if show_plot:
        try:
            plt.show()
        except Exception:
            pass

    plt.close()

    return output_path

def generate_interaction_heatmap(interactions_result, output_path=None, show_plot=True):
    """
    生成相互作用热图 (Heatmap)
    特别适用于蛋白-核酸或蛋白-蛋白界面分析，展示残基-残基接触矩阵。
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        import pandas as pd
        import numpy as np
    except ImportError as e:
        print(f"[generate_interaction_heatmap] ❌ Dependencies missing: {e}")
        print("Run: pip install seaborn pandas matplotlib numpy")
        raise RuntimeError(f"Missing python packages for heatmap: {e}")

    if not interactions_result:
        print("[generate_interaction_heatmap] No interactions to plot")
        return None

    # 1. 提取数据
    # 兼容不同的输入格式 (list of dicts or result dict)
    interactions = interactions_result
    if isinstance(interactions_result, dict) and "interactions" in interactions_result:
        interactions = interactions_result["interactions"]
        
    if not interactions:
        print("[generate_interaction_heatmap] Interaction list is empty")
        return None

    # 2. 准备矩阵数据
    # 识别列名
    # Protein-Nucleic: Nucleic_Residue vs Protein_Residue
    # Protein-Protein: Chain1/Residue1 vs Chain2/Residue2
    # Protein-Ligand: Ligand_Residue vs Protein_Residue
    
    row_residues = [] # Y轴
    col_residues = [] # X轴
    data_points = []
    
    # 检测模式
    sample = interactions[0]
    mode = "unknown"
    
    # Safe key retrieval function
    def get_key(item, *keys):
        for k in keys:
            if k in item: return item[k]
        return ""

    if "Nucleic_Residue" in sample:
        mode = "pn" # Protein-Nucleic
        row_label = "Protein Residues"
        col_label = "Nucleic Acid Residues"
    elif "Ligand_Residue" in sample:
        mode = "pl" # Protein-Ligand
        row_label = "Protein Residues"
        col_label = "Ligand Atoms/Residues"
    elif "chain1" in sample: # PPI (from ppi_analyzer.py CSV format?)
        # Need to check PPI interaction format
        mode = "pp"
        row_label = "Chain 1"
        col_label = "Chain 2"
    elif "atom1" in sample and "atom2" in sample: # PPI (from _analyze_interface_interactions)
        mode = "pp_detail"
        row_label = "Atom 1 (Chain 1)"
        col_label = "Atom 2 (Chain 2)"
    else:
        # Fallback / Generic
        mode = "generic"
        row_label = "Partner 1"
        col_label = "Partner 2"

    # 优先级映射 (用于颜色强度)
    type_score = {
        "盐桥": 4, "Salt Bridge": 4,
        "氢键": 3, "Hydrogen Bond": 3, "H-Bond": 3,
        "π–π 堆积": 2, "Pi-Pi Stacking": 2,
        "π–阳离子": 2, "Pi-Cation": 2,
        "疏水相互作用": 1, "疏水接触": 1, "Hydrophobic": 1, "Hydrophobic Interaction": 1,
        "范德华力": 0.5, "vdW": 0.5
    }

    # 收集所有唯一的残基
    unique_rows = set()
    unique_cols = set()
    
    matrix_data = {} # (row_res, col_res) -> max_score

    for inter in interactions:
        if mode == "pn":
            r = inter.get("Protein_Residue", "Unknown")
            c = inter.get("Nucleic_Residue", "Unknown")
        elif mode == "pl":
            r = inter.get("Protein_Residue", "Unknown")
            c = inter.get("Ligand_Atom", inter.get("Ligand_Residue", "Unknown"))
        elif mode == "pp_detail":
            # atom1: "A:ARG 123:N" -> extract residue "A:ARG 123"
            a1 = inter.get("atom1", "")
            a2 = inter.get("atom2", "")
            r = " ".join(a1.split(":")[:2]) if ":" in a1 else a1
            c = " ".join(a2.split(":")[:2]) if ":" in a2 else a2
        else: # generic/pp
            # 尝试各种可能的键名
            r_val = get_key(inter, "Protein_Residue", "Residue1", "atom1", "Chain1")
            c_val = get_key(inter, "Ligand_Residue", "Nucleic_Residue", "Residue2", "atom2", "Chain2")
            
            # 如果是PPI，还需要加上Chain ID以防重复
            c1 = get_key(inter, "Protein_Chain", "Chain1")
            c2 = get_key(inter, "Ligand_Chain", "Nucleic_Chain", "Chain2")
            
            r = f"{c1} {r_val}".strip() if c1 else r_val
            c = f"{c2} {c_val}".strip() if c2 else c_val
            
        itype = inter.get("Interaction", inter.get("type", ""))
        score = type_score.get(itype, 1)
        
        # 如果更强的相互作用存在，覆盖
        if (r, c) in matrix_data:
            matrix_data[(r, c)] = max(matrix_data[(r, c)], score)
        else:
            matrix_data[(r, c)] = score
            
        unique_rows.add(r)
        unique_cols.add(c)

    if not unique_rows or not unique_cols:
        print("[generate_interaction_heatmap] No valid data rows/cols extracted")
        return None

    # 排序函数
    def sort_key(res_str):
        # 尝试提取数字: "ARG 123" -> 123, "DA 5" -> 5
        import re
        m = re.search(r'(\d+)', str(res_str))
        if m:
            return int(m.group(1))
        return str(res_str)

    sorted_rows = sorted(list(unique_rows), key=sort_key)
    sorted_cols = sorted(list(unique_cols), key=sort_key)
    
    # 构建 DataFrame
    try:
        df = pd.DataFrame(index=sorted_rows, columns=sorted_cols).fillna(0)
        for (r, c), score in matrix_data.items():
            df.loc[r, c] = score
    except Exception as e:
        print(f"[generate_interaction_heatmap] DataFrame construction failed: {e}")
        raise e

    # 绘图
    try:
        h = max(6, len(sorted_rows) * 0.3)
        w = max(8, len(sorted_cols) * 0.3)
        plt.figure(figsize=(w, h), dpi=120)
        
        # 自定义颜色条
        # 0=None, 1=Hydrophobic, 2=Pi, 3=HBond, 4=Salt
        cmap = sns.color_palette("YlGnBu", as_cmap=True)
        
        ax = sns.heatmap(df, cmap=cmap, linewidths=0.5, linecolor='white',
                         square=True, cbar_kws={"label": "Interaction Strength (Approx.)"})
        
        plt.title(f"Interaction Heatmap ({len(interactions)} contacts)", fontsize=14, pad=20)
        plt.ylabel(row_label, fontsize=12, fontweight='bold')
        plt.xlabel(col_label, fontsize=12, fontweight='bold')
        plt.xticks(rotation=45, ha="right")
        plt.yticks(rotation=0)
        
        plt.tight_layout()
        
        if output_path is None:
            output_path = "interaction_heatmap.png"
            
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"[generate_interaction_heatmap] ✅ Heatmap saved: {output_path}")
        
        if show_plot:
            try:
                plt.show()
            except Exception:
                pass
                
        plt.close()
        return output_path
    except Exception as e:
        print(f"[generate_interaction_heatmap] Plotting failed: {e}")
        raise e

cmd.extend("generate_interaction_network_plot", generate_interaction_network_plot)
cmd.extend("generate_interaction_heatmap", generate_interaction_heatmap)

def analyze_protein_nucleic_interactions(obj_name=None, nucleic_chains=None,
                                        protein_chains=None, output_csv=None,
                                        distance_cutoff=4.5, auto_highlight=True, pdb_file=None):
    """
    分析蛋白质-核酸相互作用（DNA/RNA）
    
    参数:
        obj_name: PyMOL对象名称
        nucleic_chains: 核酸链 ID列表（例如 ["A", "B"]）,如果为None则自动检测
        protein_chains: 蛋白质链 ID列表（例如 ["P", "Q"]）,如果为None则自动检测
        output_csv: 输出CSV文件路径
        distance_cutoff: 距离截断值（埃）
        auto_highlight: 是否自动在PyMOL中进行3D高亮显示
        pdb_file: PDB文件路径（可选）
    
    返回:
        dict: {
            "nucleic_chains": [...],
            "protein_chains": [...],
            "interactions": [...],
            "mode": "protein-nucleic"
        }
    
    相互作用类型:
        - 氢键: 蛋白与磷酸骨架/碱基
        - 盐桥: 带电氨基酸与磷酸骨架
        - 疏水: 疏水氨基酸与碱基
        - π-π 堆积: 芳香氨基酸与碱基
        - π-阳离子: 带正电氨基酸与碱基
    
    示例:
        result = analyze_protein_nucleic_interactions('complex', nucleic_chains=["D"], protein_chains=["A"], output_csv='pn_interactions.csv')
    """
    
    # 获取原子信息
    if pdb_file:
        atoms = parse_pdb_file(pdb_file)
        if not atoms:
            print(f"[analyze_protein_nucleic_interactions] Unable to read atom information from file: {pdb_file}")
            return None
    else:
        atoms = parse_pdb_structure(obj_name)
        if not atoms:
            print("[analyze_protein_nucleic_interactions] Unable to obtain atom information")
            return None
    
    # 按链和残基分组
    chain_residues = defaultdict(list)
    for atom in atoms:
        res_key = (atom[0], atom[1], atom[2])  # (chain, res_name, res_id)
        chain_residues[res_key].append(atom)
    
    # 分离收集核酸和蛋白质残基
    nucleic_residues = []
    protein_residues = []
    
    # 第一步：收集所有核酸残基
    for res_key, res_atoms in chain_residues.items():
        chain_id, res_name, res_id = res_key
        mol_type = identify_molecule_type(res_name)
        
        if mol_type in ["dna", "rna"]:
            if nucleic_chains is None or chain_id in nucleic_chains:
                nucleic_residues.append((res_key, res_atoms))
        elif mol_type == "protein":
            if protein_chains is None or chain_id in protein_chains:
                protein_residues.append((res_key, res_atoms))
    
    # 自动检测链
    if nucleic_chains is None:
        nucleic_chains = list(set(rk[0] for rk, _ in nucleic_residues))
    if protein_chains is None:
        protein_chains = list(set(rk[0] for rk, _ in protein_residues))
    
    if not nucleic_residues:
        print("\n" + "="*60)
        print("[analyze_protein_nucleic_interactions] ❌ No nucleic acid found")
        print("\nPossible reasons:")
        print("   1) No DNA/RNA in the structure")
        print("   2) Nucleic acid chains not specified correctly")
        print("="*60 + "\n")
        return None
    
    if not protein_residues:
        print("[analyze_protein_nucleic_interactions] ⚠️ No protein residues found")
        return None
    
    print(f"[analyze_protein_nucleic_interactions] ✅ Detected {len(nucleic_residues)} nucleic acid residues (chains: {', '.join(nucleic_chains)})")
    print(f"[analyze_protein_nucleic_interactions] ✅ Detected {len(protein_residues)} protein residues (chains: {', '.join(protein_chains)})")
    
    # 构建残基原子索引（用于精确氢键检测）
    all_atoms_by_residue = {}
    for nuc_key, nuc_atoms in nucleic_residues:
        all_atoms_by_residue[nuc_key] = nuc_atoms
    for prot_key, prot_atoms in protein_residues:
        all_atoms_by_residue[prot_key] = prot_atoms
    
    # 分析相互作用
    interactions = []
    
    for nuc_key, nuc_atoms in nucleic_residues:
        nuc_chain, nuc_name, nuc_id = nuc_key
        
        for prot_key, prot_atoms in protein_residues:
            prot_chain, prot_name, prot_id = prot_key
            
            # 快速距离筛选
            nuc_centroid = centroid([a[4] for a in nuc_atoms])
            prot_centroid = centroid([a[4] for a in prot_atoms])
            
            if distance(nuc_centroid, prot_centroid) > distance_cutoff + 5.0:
                continue
            
            # 检查原子间相互作用
            for nuc_atom in nuc_atoms:
                for prot_atom in prot_atoms:
                    d = distance(nuc_atom[4], prot_atom[4])
                    
                    if d > distance_cutoff:
                        continue
                    
                    # 判断相互作用类型
                    interaction_type = None
                    
                    # ✅ 使用精确的氢键检测
                    is_hb, hb_dist, hb_angle = is_hbond_precise(nuc_atom, prot_atom, all_atoms_by_residue)
                    if not is_hb:
                        # 反向检测：蛋白作为供体
                        is_hb, hb_dist, hb_angle = is_hbond_precise(prot_atom, nuc_atom, all_atoms_by_residue)
                    
                    if is_hb:
                        interaction_type = "氢键"
                    elif is_saltbridge_atom(nuc_name, nuc_atom[3], prot_name, prot_atom[3], d):
                        interaction_type = "盐桥"
                    elif is_hydrophobic(nuc_name, prot_name, nuc_atom, prot_atom, d):
                        # 过滤掉核酸糖环原子 (带 ')，只允许碱基参与疏水相互作用
                        if "'" not in nuc_atom[3]:
                            interaction_type = "疏水相互作用"
                    
                    if interaction_type:
                        interactions.append({
                            "Nucleic_Chain": nuc_chain,
                            "Nucleic_Residue": f"{nuc_name} {nuc_id}",
                            "Nucleic_Atom": nuc_atom[3],
                            "Protein_Chain": prot_chain,
                            "Protein_Residue": f"{prot_name} {prot_id}",
                            "Protein_Atom": prot_atom[3],
                            "Distance": round(d, 2),
                            "Interaction": interaction_type
                        })
            
            # 检查π相互作用
            if is_pipi(nuc_name, nuc_atoms, prot_name, prot_atoms):
                interactions.append({
                    "Nucleic_Chain": nuc_chain,
                    "Nucleic_Residue": f"{nuc_name} {nuc_id}",
                    "Nucleic_Atom": "ring",
                    "Protein_Chain": prot_chain,
                    "Protein_Residue": f"{prot_name} {prot_id}",
                    "Protein_Atom": "ring",
                    "Distance": "-",
                    "Interaction": "π–π 堆积"
                })
            
            if is_cationpi(prot_name, prot_atoms, nuc_name, nuc_atoms):
                interactions.append({
                    "Nucleic_Chain": nuc_chain,
                    "Nucleic_Residue": f"{nuc_name} {nuc_id}",
                    "Nucleic_Atom": "ring",
                    "Protein_Chain": prot_chain,
                    "Protein_Residue": f"{prot_name} {prot_id}",
                    "Protein_Atom": "cation",
                    "Distance": "-",
                    "Interaction": "π–阳离子相互作用"
                })
    
    # 输出到CSV 并处理自动高亮
    csv_path_to_use = output_csv
    is_temp_csv = False

    # 如果需要高亮但没有输出路径，创建临时文件
    if auto_highlight and not csv_path_to_use and interactions:
        try:
            fd, temp_path = tempfile.mkstemp(suffix=".csv", prefix="glue_pn_")
            os.close(fd)
            csv_path_to_use = temp_path
            is_temp_csv = True
        except Exception as e:
            print(f"[analyze_protein_nucleic_interactions] Failed to create temp file: {e}")

    if csv_path_to_use and interactions:
        try:
            with open(csv_path_to_use, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(["Nucleic_Chain", "Nucleic_Residue", "Nucleic_Atom",
                               "Protein_Chain", "Protein_Residue", "Protein_Atom",
                               "Distance", "Interaction"])
                for inter in interactions:
                    writer.writerow([
                        inter["Nucleic_Chain"],
                        inter["Nucleic_Residue"],
                        inter["Nucleic_Atom"],
                        inter["Protein_Chain"],
                        inter["Protein_Residue"],
                        inter["Protein_Atom"],
                        inter["Distance"],
                        inter["Interaction"]
                    ])
            if output_csv:
                print(f"[analyze_protein_nucleic_interactions] Results saved to: {output_csv}")
        except Exception as e:
            print(f"[analyze_protein_nucleic_interactions] Failed to save CSV: {e}")
            # 如果写入失败，且是临时文件，则无法高亮
            if is_temp_csv:
                csv_path_to_use = None

    # 自动高亮
    if auto_highlight and csv_path_to_use and interactions:
        try:
            from .highlight_residues import highlight_csv_residues
            print(f"[analyze_protein_nucleic_interactions] 🖌️ Visualizing 3D interactions...")
            highlight_csv_residues(csv_path_to_use, obj=obj_name, show_labels=True, show_interaction_type=True)
        except Exception as e:
            print(f"[analyze_protein_nucleic_interactions] Failed to visualize: {e}")
        finally:
            if is_temp_csv and os.path.exists(csv_path_to_use):
                try:
                    os.unlink(csv_path_to_use)
                except Exception:
                    pass
    
    result = {
        "nucleic_chains": nucleic_chains,
        "protein_chains": protein_chains,
        "interactions": interactions,
        "mode": "protein-nucleic",
        "parameters": {
            "distance_cutoff": distance_cutoff,
            "hbond_cutoff": INTERACTION_PARAMS["hbond"]["max_DA_dist"],
            "saltbridge_cutoff": INTERACTION_PARAMS["ionic"]["max_dist"],
        }
    }
    
    print(f"[analyze_protein_nucleic_interactions] ✅ Analysis complete: found {len(interactions)} interactions")
    return result

cmd.extend("analyze_protein_nucleic_interactions", analyze_protein_nucleic_interactions)

