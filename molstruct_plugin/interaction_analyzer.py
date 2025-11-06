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
RDKIT_AVAILABLE = False
SCIPY_AVAILABLE = False

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import numpy as np
    RDKIT_AVAILABLE = True
    print("[MolStruct] ✅ RDKit 可用，启用高级分析模式")
except ImportError:
    print("[MolStruct] ℹ️ RDKit 未安装，使用基础分析模式")
    print("[MolStruct] 💡 安装提示: pip install rdkit scipy")

try:
    from scipy.spatial import cKDTree
    SCIPY_AVAILABLE = True
except ImportError:
    pass

# ========== 相互作用参数（精确判定标准）==========
INTERACTION_PARAMS = {
    "hbond": {
        "max_DA_dist": 2.5,         # Å，D···A 距离
        "min_donor_angle": 120,     # °，∠D–H···A
        "min_acceptor_angle": 90    # °，∠H···A–X
    },
    "hydrophobic": {
        "pi_cation_max": 4.5,       # Å
        "pi_pi_mode": "face_face_or_edge",  # 按环面法向量与距离联合判定
        "other_max": 3.6            # Å
    },
    "ionic": {
        "max_dist": 3.7,            # Å
        "exclude_if_hbond": True    # 排除氢键情况
    },
    "metal_coord": {
        "max_dist": 3.4,            # Å
        "allowed_ligand_atoms": ["N", "O", "S", "CL", "BR", "F"]  # 非碳重原子
    },
    "water_bridge": {
        "max_DA_dist": 2.8,         # Å
        "min_donor_angle": 110,     # °
        "min_acceptor_angle": 90    # °
    },
    # 保留旧参数以兼容其他函数
    "saltbridge": {"max_distance": 3.7},
    "pi_pi": {"max_distance": 5.5, "min_distance": 3.3, "parallel_angle": 30.0},
    "pi_cation": {"max_distance": 4.5},
    "halogen": {"max_distance": 4.0, "min_angle": 140.0},
    "metal": {"max_distance": 3.4}
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
        objs = cmd.get_object_list()
        if not objs:
            print("[parse_pdb_structure] 没有加载的对象")
            return []
        obj_name = objs[0]
    
    if obj_name not in cmd.get_object_list():
        print(f"[parse_pdb_structure] 对象 '{obj_name}' 未找到")
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
        print(f"[parse_pdb_file] 读取PDB文件失败: {e}")
    
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

def is_saltbridge(res1, res2, d):
    """
    判断是否为盐桥（简单版本，兼容旧代码）
    
    参数:
        res1, res2: 残基名称
        d: 距离（Å）
    
    返回:
        bool: 是否为盐桥
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
        atom1, atom2: 原子名称
        d: 距离（Å）
    
    返回:
        bool: 是否为疏水相互作用
    
    标准：疏水接触 ≤ 4.5 Å
    """
    hydrophobic = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "TRP", "TYR", "CYS"}
    max_dist = INTERACTION_PARAMS["hydrophobic"]["max_distance"]
    return res1 in hydrophobic and res2 in hydrophobic and d <= max_dist

def centroid(coords):
    """计算坐标质心"""
    n = len(coords)
    if n == 0:
        return (0, 0, 0)
    return tuple(sum(c[i] for c in coords)/n for i in range(3))

def ring_atoms(res_name, atoms):
    """获取环状结构原子坐标"""
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
    
    if res_name not in ring_dict:
        return None
    
    coords = [a[4] for a in atoms if a[3].strip() in ring_dict[res_name]]
    return coords if len(coords) >= 3 else None

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
        print(f"[analyze_interactions] 检测到 {len(chains)} 条链，只分析链间接触界面...")

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
                            print(f"[analyze_interactions] 进度: {percent}% ({checked}/{total_pairs})")
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
        print(f"[analyze_interactions] 分析 {len(keys)} 个残基的所有相互作用...")

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
                    print(f"[analyze_interactions] 进度: {percent}% ({checked}/{total_pairs})")
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
    """
    # 检查原子间相互作用
    interaction_type = None
    for at1 in a1:
        for at2 in a2:
            d = distance(at1[4], at2[4])
            if d > 4.5:  # 距离阈值
                continue

            if is_hbond(at1, at2, d):
                interaction_type = "氢键"
            elif is_saltbridge(r1, r2, d):
                interaction_type = "盐桥"
            elif is_hydrophobic(r1, r2, at1, at2, d):
                interaction_type = "疏水相互作用"

            if interaction_type:
                interactions.append({
                    "Chain1": c1,
                    "Residue1": f"{r1} {id1}",
                    "Chain2": c2,
                    "Residue2": f"{r2} {id2}",
                    "Distance": round(d, 2),
                    "Interaction": interaction_type
                })
                break  # 找到一个相互作用就跳出内层循环
        if interaction_type:
            break  # 找到一个相互作用就跳出外层循环

    # 检查π相互作用
    if is_pipi(r1, a1, r2, a2):
        interactions.append({
            "Chain1": c1,
            "Residue1": f"{r1} {id1}",
            "Chain2": c2,
            "Residue2": f"{r2} {id2}",
            "Distance": "-",
            "Interaction": "π–π 堆积"
        })

    if is_cationpi(r1, a1, r2, a2):
        interactions.append({
            "Chain1": c1,
            "Residue1": f"{r1} {id1}",
            "Chain2": c2,
            "Residue2": f"{r2} {id2}",
            "Distance": "-",
            "Interaction": "π–阳离子相互作用"
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
            print(f"[analyze_pdb_interactions] 无法从文件读取原子信息: {pdb_file}")
            return []
    else:
        atoms = parse_pdb_structure(obj_name)
        if not atoms:
            print("[analyze_pdb_interactions] 无法获取原子信息")
            return []
    
    print("[analyze_pdb_interactions] 使用严格标准（氢键 ≤2.8Å）")
    
    # 分析相互作用
    interactions = analyze_interactions(atoms, only_between_chains)
    
    # 输出到CSV文件
    if output_csv:
        try:
            with open(output_csv, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(["Chain1", "Residue1", "Chain2", "Residue2", "Distance", "Interaction"])
                for inter in interactions:
                    writer.writerow([
                        inter["Chain1"],
                        inter["Residue1"],
                        inter["Chain2"],
                        inter["Residue2"],
                        inter["Distance"],
                        inter["Interaction"]
                    ])
            print(f"[analyze_pdb_interactions] 结果已保存到: {output_csv}")
        except Exception as e:
            print(f"[analyze_pdb_interactions] 保存CSV文件失败: {e}")
    
    # 自动高亮显示（如果启用且在PyMOL环境中）
    if auto_highlight and not pdb_file and interactions:
        try:
            # 创建临时CSV文件
            temp_csv = tempfile.NamedTemporaryFile(mode='w', suffix='.csv',
                                                 delete=False, encoding='utf-8')
            writer = csv.writer(temp_csv)
            writer.writerow(["Chain1", "Residue1", "Chain2", "Residue2", "Distance", "Interaction"])
            for inter in interactions:
                writer.writerow([
                    inter["Chain1"],
                    inter["Residue1"],
                    inter["Chain2"],
                    inter["Residue2"],
                    inter["Distance"],
                    inter["Interaction"]
                ])
            temp_csv.close()

            # 调用高亮功能
            from .highlight_residues import highlight_csv_residues
            highlight_csv_residues(temp_csv.name, obj=obj_name, show_labels=1,
                                 stick_by_element=1)

            # 清理临时文件
            os.unlink(temp_csv.name)

        except Exception as e:
            print(f"[analyze_pdb_interactions] 自动高亮失败: {e}")
    
    print(f"[analyze_pdb_interactions] 分析完成: 发现 {len(interactions)} 个相互作用")
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

    print(f"[render_interactions_beautifully] 渲染完成")

cmd.extend("render_interactions_beautifully", render_interactions_beautifully)

def analyze_protein_ligand_interactions(obj_name=None, ligand_resname=None,
                                       protein_chains=None, output_csv=None,
                                       distance_cutoff=4.5, pdb_file=None):
    """
    分析蛋白质-配体相互作用（自动使用RDKit高级模式，如果可用）

    参数:
        obj_name: PyMOL对象名称
        ligand_resname: 配体残基名称（例如 "LIG", "ATP" 等）,如果为None则自动检测
        protein_chains: 蛋白质链ID列表（例如 ["A", "B"]）,如果为None则自动检测
        output_csv: 输出CSV文件路径
        distance_cutoff: 距离截断值（埃）
        pdb_file: PDB文件路径（可选）

    返回:
        dict: {
            "ligand_residues": [{"chain": X, "resname": Y, "resid": Z}, ...],
            "protein_chains": ["A", "B", ...],
            "interactions": [...相互作用列表...]
            "mode": "advanced" or "basic"
        }
    
    使用严格标准：氢键 ≤2.8Å，盐桥 ≤4.0Å，适合发表
    
    示例:
        result = analyze_protein_ligand_interactions('protein', 'LIG', output_csv='interactions.csv')
    """
    
    # ========== 自动检测并使用RDKit高级模式 ==========
    if RDKIT_AVAILABLE and SCIPY_AVAILABLE:
        try:
            print("[analyze_protein_ligand_interactions] 🚀 使用 RDKit 高级分析模式")
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
        except Exception as e:
            print(f"[analyze_protein_ligand_interactions] RDKit模式失败，降级到基础模式: {e}")
            import traceback
            traceback.print_exc()
    
    # ========== 基础分析模式 ==========
    print("[analyze_protein_ligand_interactions] 📊 使用基础分析模式")
    
    # 获取原子信息
    if pdb_file:
        atoms = parse_pdb_file(pdb_file)
        if not atoms:
            print(f"[analyze_protein_ligand_interactions] 无法从文件读取原子信息: {pdb_file}")
            return None
    else:
        atoms = parse_pdb_structure(obj_name)
        if not atoms:
            print("[analyze_protein_ligand_interactions] 无法获取原子信息")
            return None

    # 按链和残基分组
    chain_residues = defaultdict(list)
    for atom in atoms:
        res_key = (atom[0], atom[1], atom[2])  # (chain, res_name, res_id)
        chain_residues[res_key].append(atom)

    print(f"[analyze_protein_ligand_interactions] 📋 结构信息：共 {len(chain_residues)} 个残基")
    
    # 统计残基类型
    residue_type_count = defaultdict(int)
    for res_key in chain_residues.keys():
        _, res_name, _ = res_key
        mol_type = identify_molecule_type(res_name)
        residue_type_count[mol_type] += 1
    
    print(f"[analyze_protein_ligand_interactions] 残基类型分布：{dict(residue_type_count)}")
    
    # 自动检测配体（如果未指定）
    ligand_residues = []
    protein_residues = []

    for res_key, res_atoms in chain_residues.items():
        chain_id, res_name, res_id = res_key
        mol_type = identify_molecule_type(res_name)

        if ligand_resname:
            # 用户指定了配体
            if res_name.upper() == ligand_resname.upper():
                ligand_residues.append((res_key, res_atoms))
                print(f"[analyze_protein_ligand_interactions]   ✓ 配体: {res_name} {res_id} (链 {chain_id})")
        else:
            # 自动检测
            if mol_type == "ligand":
                ligand_residues.append((res_key, res_atoms))
                print(f"[analyze_protein_ligand_interactions]   ✓ 自动检测配体: {res_name} {res_id} (链 {chain_id})")
            elif mol_type == "protein":
                protein_residues.append((res_key, res_atoms))

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
        print(f"[analyze_protein_ligand_interactions] 过滤蛋白链: {old_count} -> {len(protein_residues)} 个残基")

    if not ligand_residues:
        print("[analyze_protein_ligand_interactions] ⚠️ 未找到配体分子")
        print(f"[analyze_protein_ligand_interactions]    提示：检查配体残基名是否正确，或结构中是否包含非标准残基")
        # 列出所有非标准残基
        non_standard = []
        for res_key in chain_residues.keys():
            _, res_name, res_id = res_key
            if identify_molecule_type(res_name) not in ("protein", "solvent"):
                non_standard.append(f"{res_name}:{res_id}")
        if non_standard:
            print(f"[analyze_protein_ligand_interactions]    发现的非标准残基: {', '.join(non_standard[:10])}")
        return None

    if not protein_residues:
        print("[analyze_protein_ligand_interactions] ⚠️ 未找到蛋白质残基")
        return None

    print(f"[analyze_protein_ligand_interactions] ✓ 检测到 {len(ligand_residues)} 个配体分子")
    print(f"[analyze_protein_ligand_interactions] ✓ 检测到 {len(protein_residues)} 个蛋白质残基（链: {', '.join(protein_chains)}）")
    print("[analyze_protein_ligand_interactions] 🔬 使用精确判定标准")
    print(f"  - 氢键: D···A ≤ {INTERACTION_PARAMS['hbond']['max_DA_dist']} Å, ∠D-H-A ≥ {INTERACTION_PARAMS['hbond']['min_donor_angle']}°")
    print(f"  - 离子/盐桥: ≤ {INTERACTION_PARAMS['ionic']['max_dist']} Å")
    print(f"  - π-阳离子: ≤ {INTERACTION_PARAMS['hydrophobic']['pi_cation_max']} Å")
    print(f"  - 金属配位: ≤ {INTERACTION_PARAMS['metal_coord']['max_dist']} Å")
    print(f"  - 疏水: ≤ {INTERACTION_PARAMS['hydrophobic']['other_max']} Å")

    # 分析相互作用
    interactions = []

    for lig_key, lig_atoms in ligand_residues:
        lig_chain, lig_name, lig_id = lig_key

        print(f"[analyze_protein_ligand_interactions] 分析配体: {lig_name} {lig_id} (链 {lig_chain})")

        for prot_key, prot_atoms in protein_residues:
            prot_chain, prot_name, prot_id = prot_key

            # 快速距离筛选
            lig_centroid = centroid([a[4] for a in lig_atoms])
            prot_centroid = centroid([a[4] for a in prot_atoms])

            if distance(lig_centroid, prot_centroid) > distance_cutoff + 5.0:
                continue

            # 检查原子间相互作用
            for lig_atom in lig_atoms:
                for prot_atom in prot_atoms:
                    d = distance(lig_atom[4], prot_atom[4])

                    if d > distance_cutoff:
                        continue

                    # 判断相互作用类型
                    interaction_type = None

                    if is_hbond(lig_atom, prot_atom, d):
                        interaction_type = "氢键"
                    elif is_saltbridge(lig_name, prot_name, d):
                        interaction_type = "盐桥"
                    elif is_hydrophobic(lig_name, prot_name, lig_atom, prot_atom, d):
                        interaction_type = "疏水相互作用"

                    # —— 兼容性增强（蛋白-配体专用，放宽与更合理的规则）——
                    if interaction_type is None:
                        try:
                            elem_l = lig_atom[3].strip()[0].upper()
                            elem_p = prot_atom[3].strip()[0].upper()
                        except Exception:
                            elem_l, elem_p = "", ""

                        # 1) 更宽松的氢键近邻：N/O 之间 ≤ 3.5 Å（无角度，基础模式）
                        if interaction_type is None and (elem_l in ("N", "O") and elem_p in ("N", "O")) and d <= 3.5:
                            interaction_type = "氢键"

                        # 2) 疏水相互作用：C-C 之间 ≤ 4.5 Å（不要求配体残基名）
                        if interaction_type is None and elem_l == "C" and elem_p == "C" and d <= 4.5:
                            interaction_type = "疏水相互作用"

                        # 3) 简化盐桥：蛋白正电侧链原子(NZ/NH1/NH2/NE) 与 配体氧原子(O*) ≤ 4.0 Å
                        if interaction_type is None and d <= 4.0:
                            pos_atoms = {"NZ", "NH1", "NH2", "NE"}
                            lig_is_o = elem_l == "O"
                            prot_is_pos = prot_atom[3].strip().upper() in pos_atoms
                            if prot_is_pos and lig_is_o:
                                interaction_type = "盐桥"

                    if interaction_type:
                        interactions.append({
                            "Ligand_Chain": lig_chain,
                            "Ligand_Residue": f"{lig_name} {lig_id}",
                            "Ligand_Atom": lig_atom[3],
                            "Protein_Chain": prot_chain,
                            "Protein_Residue": f"{prot_name} {prot_id}",
                            "Protein_Atom": prot_atom[3],
                            "Distance": round(d, 2),
                            "Interaction": interaction_type
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
                    "Interaction": "π–π 堆积"
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
                    "Interaction": "π–阳离子相互作用"
                })

    # 输出到CSV
    if output_csv and interactions:
        try:
            with open(output_csv, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(["Ligand_Chain", "Ligand_Residue", "Ligand_Atom",
                               "Protein_Chain", "Protein_Residue", "Protein_Atom",
                               "Distance", "Interaction"])
                for inter in interactions:
                    writer.writerow([
                        inter["Ligand_Chain"],
                        inter["Ligand_Residue"],
                        inter["Ligand_Atom"],
                        inter["Protein_Chain"],
                        inter["Protein_Residue"],
                        inter["Protein_Atom"],
                        inter["Distance"],
                        inter["Interaction"]
                    ])
            print(f"[analyze_protein_ligand_interactions] 结果已保存到: {output_csv}")
        except Exception as e:
            print(f"[analyze_protein_ligand_interactions] 保存CSV文件失败: {e}")

    result = {
        "ligand_residues": [{"chain": k[0], "resname": k[1], "resid": k[2]}
                           for k, _ in ligand_residues],
        "protein_chains": protein_chains,
        "interactions": interactions,
        "mode": "basic",
        "parameters": {
            "distance_cutoff": distance_cutoff,
            "hbond_cutoff": INTERACTION_PARAMS["hbond"]["max_distance"],
            "saltbridge_cutoff": INTERACTION_PARAMS["saltbridge"]["max_distance"],
            "hydrophobic_cutoff": INTERACTION_PARAMS["hydrophobic"]["max_distance"]
        }
    }

    print(f"[analyze_protein_ligand_interactions] ✅ 分析完成: 发现 {len(interactions)} 个相互作用")
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
            print(f"[analyze_ternary_complex] 无法从文件读取原子信息: {pdb_file}")
            return None
    else:
        atoms = parse_pdb_structure(obj_name)
        if not atoms:
            print("[analyze_ternary_complex] 无法获取原子信息")
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
        print("[analyze_ternary_complex] ⚠️ 未找到配体分子")
        return None

    # 自动检测蛋白质链（如果未指定）
    if protein1_chains is None or protein2_chains is None:
        all_chains = set(atom[0] for atom in atoms if identify_molecule_type(atom[1]) == "protein")
        all_chains = sorted(all_chains)

        if len(all_chains) < 2:
            print("[analyze_ternary_complex] ⚠️ 未找到足够的蛋白质链（需要至少2条）")
            return None

        if protein1_chains is None:
            protein1_chains = [all_chains[0]]
        if protein2_chains is None:
            protein2_chains = [all_chains[1]] if len(all_chains) > 1 else [all_chains[0]]

    print(f"[analyze_ternary_complex] 蛋白质1链: {', '.join(protein1_chains)}")
    print(f"[analyze_ternary_complex] 蛋白质2链: {', '.join(protein2_chains)}")
    print(f"[analyze_ternary_complex] 配体: {len(ligand_residues)} 个分子")

    # 分别分析配体与两个蛋白的相互作用
    print("\n[analyze_ternary_complex] 分析配体 - 蛋白质1相互作用...")
    protein1_result = analyze_protein_ligand_interactions(
        obj_name=obj_name,
        ligand_resname=ligand_resname,
        protein_chains=protein1_chains,
        distance_cutoff=distance_cutoff,
        pdb_file=pdb_file
    )

    print("\n[analyze_ternary_complex] 分析配体 - 蛋白质2相互作用...")
    protein2_result = analyze_protein_ligand_interactions(
        obj_name=obj_name,
        ligand_resname=ligand_resname,
        protein_chains=protein2_chains,
        distance_cutoff=distance_cutoff,
        pdb_file=pdb_file
    )

    if not protein1_result or not protein2_result:
        print("[analyze_ternary_complex] ⚠️ 未能完成三元复合体分析")
        return None

    # 分析桥接效应
    print("\n[analyze_ternary_complex] 分析桥接效应...")

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

                # 写入标题
                writer.writerow(["=== 三元复合体分析结果 ==="])
                writer.writerow([])

                # 配体信息
                writer.writerow(["配体信息"])
                writer.writerow(["配体数量", len(ligand_residues)])
                for lig_key, _ in ligand_residues:
                    writer.writerow(["配体", f"{lig_key[1]} {lig_key[2]} (链 {lig_key[0]})"])
                writer.writerow([])

                # 蛋白质1相互作用
                writer.writerow([f"=== 蛋白质1（链 {', '.join(protein1_chains)}）- 配体相互作用 ==="])
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

                # 蛋白质2相互作用
                writer.writerow([f"=== 蛋白质2（链 {', '.join(protein2_chains)}）- 配体相互作用 ==="])
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

                # 桥接分析统计
                writer.writerow(["=== 桥接分析统计 ==="])
                writer.writerow(["蛋白质1接触残基数", bridging_analysis["protein1_contacting_residues"]])
                writer.writerow(["蛋白质2接触残基数", bridging_analysis["protein2_contacting_residues"]])
                writer.writerow(["蛋白质1相互作用数", bridging_analysis["protein1_interactions_count"]])
                writer.writerow(["蛋白质2相互作用数", bridging_analysis["protein2_interactions_count"]])
                writer.writerow(["总相互作用数", bridging_analysis["total_interactions"]])

            print(f"[analyze_ternary_complex] 结果已保存到: {output_csv}")
        except Exception as e:
            print(f"[analyze_ternary_complex] 保存CSV文件失败: {e}")

    result = {
        "ligand_info": protein1_result["ligand_residues"],
        "protein1_chains": protein1_chains,
        "protein2_chains": protein2_chains,
        "protein1_interactions": p1_interactions,
        "protein2_interactions": p2_interactions,
        "bridging_analysis": bridging_analysis
    }

    print(f"\n[analyze_ternary_complex] ✅ 三元复合体分析完成")
    print(f"  - 蛋白质1: {bridging_analysis['protein1_interactions_count']} 个相互作用")
    print(f"  - 蛋白质2: {bridging_analysis['protein2_interactions_count']} 个相互作用")
    print(f"  - 总计: {bridging_analysis['total_interactions']} 个相互作用")

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
            print(f"[analyze_atom_pair_interactions] 无法从文件读取原子信息: {pdb_file}")
            return None
    else:
        atoms = parse_pdb_structure(obj_name)
        if not atoms:
            print("[analyze_atom_pair_interactions] 无法获取原子信息")
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
        print(f"[analyze_atom_pair_interactions] ⚠️ 未找到匹配的原子1: {atom1_selection}")
        return None

    if not atoms2:
        print(f"[analyze_atom_pair_interactions] ⚠️ 未找到匹配的原子2: {atom2_selection}")
        return None

    print(f"[analyze_atom_pair_interactions] 原子1: {len(atoms1)} 个")
    print(f"[analyze_atom_pair_interactions] 原子2: {len(atoms2)} 个")

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
            print(f"[analyze_atom_pair_interactions] 结果已保存到: {output_csv}")
        except Exception as e:
            print(f"[analyze_atom_pair_interactions] 保存CSV文件失败: {e}")

    print(f"[analyze_atom_pair_interactions] ✅ 分析完成: 发现 {len(interactions)} 个原子对相互作用")

    # 打印前5个最近的相互作用
    if interactions:
        print("\n前5个最近的相互作用:")
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
        direction = "X轴" if dx > 0 else "-X轴"
    elif abs_max == abs(dy):
        direction = "Y轴" if dy > 0 else "-Y轴"
    else:
        direction = "Z轴" if dz > 0 else "-Z轴"

    return f"主向{direction}"

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
            print(f"[visualize_atom_pairs] 读取CSV失败: {e}")
            return

    if not interactions_result:
        print("[visualize_atom_pairs] 没有相互作用数据")
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
            cmd.distance(pair_name, sel1, sel2)

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

    print(f"[visualize_atom_pairs] ✅ 已可视化 {min(len(interactions_result), 20)} 个原子对相互作用")
    print(f"   使用 'hide labels, atom_pairs_*' 隐藏标签")
    print(f"   使用 'delete atom_pairs_*' 删除所有可视化")

cmd.extend("visualize_atom_pairs", visualize_atom_pairs)

def visualize_protein_ligand_3d(obj_name, interactions_result=None, ligand_resname=None, csv_path=None):
    """
    在PyMOL中3D可视化蛋白-配体相互作用（改进版，参考专业脚本）

    参数:
        obj_name: PyMOL对象名称
        interactions_result: analyze_protein_ligand_interactions的返回结果
        ligand_resname: 配体残基名称(可选)
        csv_path: 或者提供CSV文件路径
    """
    from pymol import cmd

    # 读取相互作用数据
    if csv_path and not interactions_result:
        interactions = []
        try:
            with open(csv_path, 'r', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    interactions.append(row)
        except Exception as e:
            print(f"[visualize_protein_ligand_3d] 读取CSV失败: {e}")
            return
    elif interactions_result:
        if isinstance(interactions_result, dict) and "interactions" in interactions_result:
            interactions = interactions_result["interactions"]
        else:
            interactions = interactions_result
    else:
        print("[visualize_protein_ligand_3d] 需要提供相互作用数据")
        return

    if not interactions:
        print("[visualize_protein_ligand_3d] 没有相互作用数据")
        return

    # ========== 第一步：清理和基础设置 ==========
    print("[visualize_protein_ligand_3d] 🎨 开始3D可视化...")
    
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
        print(f"[visualize_protein_ligand_3d] ❌ 无配体，无法继续")
        return

    # ========== 第三步：隐藏所有，准备重新显示 ==========
    cmd.hide("everything", obj_name)
    
    # ========== 第四步：选择并创建口袋区域 ==========
    cmd.select("lig_pocket", f"byres ({obj_name} and polymer within 5 of ({lig_sel}))")
    
    # ========== 第五步：添加氢原子（关键！）==========
    print(f"[visualize_protein_ligand_3d] ➕ 添加氢原子...")
    
    # 检测配体是否已有氢原子
    n_h_lig = cmd.count_atoms(f"({lig_sel}) and hydro")
    if n_h_lig == 0:
        try:
            cmd.h_add(lig_sel)
            print(f"[visualize_protein_ligand_3d]    ✓ 配体氢原子已添加")
        except Exception as e:
            print(f"[visualize_protein_ligand_3d]    ⚠️ 配体氢原子添加失败: {e}")
    
    # 检测口袋残基是否已有氢原子
    n_h_pocket = cmd.count_atoms("lig_pocket and hydro")
    if n_h_pocket == 0:
        try:
            cmd.h_add("lig_pocket")
            print(f"[visualize_protein_ligand_3d]    ✓ 口袋残基氢原子已添加")
        except Exception as e:
            print(f"[visualize_protein_ligand_3d]    ⚠️ 口袋残基氢原子添加失败: {e}")

    # ========== 第六步：基础显示设置 ==========
    # 背景白色
    cmd.bg_color("white")
    
    # 蛋白整体：半透明cartoon
    cmd.show("cartoon", obj_name)
    cmd.set("cartoon_transparency", 0.3, obj_name)
    
    # 显示配体和口袋残基为sticks
    cmd.show("sticks", lig_sel)
    cmd.show("sticks", "lig_pocket")
    
    # 隐藏连接到碳原子的氢（只保留极性氢：NH, OH, SH）
    cmd.hide("(h. and (e. c extend 1))")
    
    # ========== 第七步：配体着色（参考脚本：配体碳原子黄色） ==========
    cmd.color("yellow", f"{lig_sel} and name C*")     # 配体碳原子：黄色
    cmd.color("blue", f"{lig_sel} and elem N")
    cmd.color("red", f"{lig_sel} and elem O")
    cmd.color("yellow", f"{lig_sel} and elem S")
    cmd.color("green", f"{lig_sel} and elem F+CL+BR+I")
    
    # 口袋残基碳原子：青色（参考脚本）
    cmd.color("cyan", "lig_pocket and name C*")
    cmd.color("blue", "lig_pocket and elem N")
    cmd.color("red", "lig_pocket and elem O")
    cmd.color("yellow", "lig_pocket and elem S")
    
    # ========== 第八步：检测并绘制氢键（参考脚本逻辑）==========
    print(f"[visualize_protein_ligand_3d] 🔗 检测氢键...")
    
    try:
        # 定义极性供体（N/O 且带氢）
        cmd.select("polar_donors_lig", f"{lig_sel} and elem n,o and (neighbor hydro)")
        cmd.select("polar_donors_res", "lig_pocket and elem n,o and (neighbor hydro)")
        
        # 定义极性受体（O 或者 N且不带氢）
        cmd.select("polar_acceptors_lig", f"{lig_sel} and (elem o or (elem n and not (neighbor hydro)))")
        cmd.select("polar_acceptors_res", "lig_pocket and (elem o or (elem n and not (neighbor hydro)))")
        
        # 定义氢原子（连接到供体的）
        cmd.select("don_hydrogens_lig", "hydro and (neighbor polar_donors_lig)")
        cmd.select("don_hydrogens_res", "hydro and (neighbor polar_donors_res)")
        
        # 绘制氢键（双向）
        # 蛋白供体 H -> 配体受体
        n_hbonds1 = 0
        try:
            cmd.distance("hbonds_res_to_lig", "don_hydrogens_res", "polar_acceptors_lig", 3.2)
            n_hbonds1 = cmd.count_atoms("hbonds_res_to_lig") // 2
        except Exception as e:
            print(f"[visualize_protein_ligand_3d]    ℹ️ 蛋白→配体氢键: {e}")
        
        # 配体供体 H -> 蛋白受体
        n_hbonds2 = 0
        try:
            cmd.distance("hbonds_lig_to_res", "don_hydrogens_lig", "polar_acceptors_res", 3.2)
            n_hbonds2 = cmd.count_atoms("hbonds_lig_to_res") // 2
        except Exception as e:
            print(f"[visualize_protein_ligand_3d]    ℹ️ 配体→蛋白氢键: {e}")
        
        total_hbonds = n_hbonds1 + n_hbonds2
        if total_hbonds > 0:
            print(f"[visualize_protein_ligand_3d]    ✅ 检测到 {total_hbonds} 个氢键")
        else:
            print(f"[visualize_protein_ligand_3d]    ℹ️ 未检测到氢键（可能需要调整参数）")
            
    except Exception as e:
        print(f"[visualize_protein_ligand_3d]    ⚠️ 氢键检测失败: {e}")
        import traceback
        traceback.print_exc()
    
    # ========== 第九步：氢键样式设置（参考脚本）==========
    cmd.set("dash_length", 0.3)
    cmd.set("dash_radius", 0.08)
    cmd.set("dash_color", "yellow")  # 黄色虚线（参考脚本）
    cmd.set("dash_width", 2.0)
    cmd.set("dash_gap", 0.5)
    cmd.hide("labels", "hbonds_*")  # 隐藏距离标签
    
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
        
        # 显示相互作用残基为 licorice（参考脚本）
        for i, (chain, res_name, res_id) in enumerate(protein_residues, 1):
            res_sel = f"{obj_name} and chain {chain} and resi {res_id}"
            cmd.show("licorice", res_sel)
            # 隐藏这些残基上连接到碳的氢
            cmd.hide("(elem H and neighbor elem C)", res_sel)
    
    # ========== 第十一步：绘制其他类型相互作用 ==========
    print(f"[visualize_protein_ligand_3d] 🎨 绘制其他相互作用...")
    
    color_map = {
        "盐桥": "red",
        "Salt": "red",
        "疏水相互作用": "green",
        "Hydrophobic": "green",
        "π–π 堆积": "purple",
        "PiPi": "purple",
        "π–阳离子相互作用": "magenta",
        "PiCation": "magenta",
    }
    
    interaction_count = {}
    
    for i, inter in enumerate(interactions[:50], 1):  # 最多显示50个
        interaction_type = inter.get("Interaction", "")
        
        # 跳过氢键（已经用专门的方法绘制了）
        if "氢键" in interaction_type or "Hbond" in interaction_type:
            continue
        
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

            # 创建选择
            sel1 = f"{obj_name} and chain {lig_chain} and resi {lig_resid}"
            sel2 = f"{obj_name} and chain {prot_chain} and resi {prot_resid}"

            if lig_atom and lig_atom not in ("ring", "ring/cation"):
                sel1 += f" and name {lig_atom}"
            if prot_atom and prot_atom not in ("ring", "ring/cation"):
                sel2 += f" and name {prot_atom}"

            # 绘制距离
            dist_name = f"pl_interact_{interaction_type.replace(' ', '_').replace('–', '_')}_{i}"
            cmd.distance(dist_name, sel1, sel2)

            # 设置颜色
            for key, color in color_map.items():
                if key in interaction_type:
                    cmd.color(color, dist_name)
                    break

            # 统计相互作用类型
            interaction_count[interaction_type] = interaction_count.get(interaction_type, 0) + 1

        except Exception as e:
            # print(f"[visualize_protein_ligand_3d] ⚠️ 第 {i} 个相互作用跳过: {e}")
            continue

    # 隐藏其他距离标签
    cmd.hide("labels", "pl_interact_*")
    
    # ========== 第十二步：调整视角 ==========
    cmd.zoom(f"({lig_sel}) or lig_pocket", buffer=8)
    cmd.orient(f"({lig_sel}) or lig_pocket")
    
    # ========== 第十三步：美化参数（参考脚本）==========
    cmd.set("stick_radius", 0.15)
    cmd.set("sphere_scale", 0.25)
    cmd.set("ray_trace_mode", 1)
    
    # 设置标签样式
    cmd.set("label_size", 28)
    cmd.set("label_font_id", 8)
    cmd.set("label_color", "grey")
    
    # ========== 第十四步：输出统计 ==========
    print(f"\n[visualize_protein_ligand_3d] ✅ 3D可视化完成！")
    print(f"   📊 相互作用统计:")
    
    # 统计氢键
    total_hbonds = 0
    try:
        n1 = cmd.count_atoms("hbonds_res_to_lig")
        n2 = cmd.count_atoms("hbonds_lig_to_res")
        total_hbonds = (n1 + n2) // 2
        if total_hbonds > 0:
            print(f"      • 氢键: {total_hbonds}")
    except:
        pass
    
    # 统计其他相互作用
    for itype, count in sorted(interaction_count.items(), key=lambda x: x[1], reverse=True):
        if count > 0:
            print(f"      • {itype}: {count}")
    
    if protein_residues:
        print(f"   📍 涉及 {len(protein_residues)} 个蛋白残基")
    
    print(f"\n   💡 PyMOL 命令提示:")
    print(f"      show labels, hbonds_*       # 显示氢键距离")
    print(f"      hide labels, hbonds_*       # 隐藏氢键距离")
    print(f"      show labels, pl_interact_*  # 显示其他相互作用距离")
    print(f"      delete pl_interact_*        # 删除所有相互作用可视化")
    print(f"      delete hbonds_*             # 删除氢键可视化")
    print(f"      ray                         # 生成高质量渲染")
    print(f"      png output.png, dpi=300     # 保存高分辨率图片")
    print(f"      set cartoon_transparency, 0.5  # 调整cartoon透明度")

cmd.extend("visualize_protein_ligand_3d", visualize_protein_ligand_3d)

def generate_advanced_interaction_plot(interactions, output_path=None, show_plot=True, ligand_sdf=None):
    """
    为高级分析生成2D相互作用图（带配体化学结构）
    
    参数:
        interactions: 相互作用列表（高级格式）
        output_path: 输出路径
        show_plot: 是否显示
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.lines import Line2D
        import matplotlib
        matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
        matplotlib.rcParams['axes.unicode_minus'] = False
    except ImportError:
        print("[generate_advanced_interaction_plot] 需要安装matplotlib")
        return None
    
    # 尝试导入RDKit用于绘制配体结构
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
        from PIL import Image
        import io
        import numpy as np
        RDKIT_AVAILABLE = True
    except ImportError:
        RDKIT_AVAILABLE = False
        print("[generate_advanced_interaction_plot] RDKit未安装，使用简化版本")
    
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
        print("[generate_advanced_interaction_plot] 没有可绘制的数据")
        return None
    
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
    
    # 尝试绘制配体结构
    from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
    ligand_drawn = False
    atom_coords_2d = {}  # 配体原子的2D坐标
    
    import os
    if RDKIT_AVAILABLE and ligand_sdf and os.path.exists(ligand_sdf):
        try:
            # 从SDF文件读取配体
            supplier = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
            mol = supplier[0] if len(supplier) > 0 else None
            
            if mol is not None:
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
                
                print("[generate_advanced_interaction_plot] ✅ 已绘制配体化学结构，高亮相互作用原子")
        except Exception as e:
            print(f"[generate_advanced_interaction_plot] 绘制配体结构失败: {e}")
            import traceback
            traceback.print_exc()
    
    # 如果无法绘制真实结构，使用占位符
    if not ligand_drawn:
        ligand_box = FancyBboxPatch((-2, -1.5), 4, 3,
                                    boxstyle="round,pad=0.1",
                                    facecolor='#FFE082',
                                    edgecolor='black',
                                    linewidth=3,
                                    zorder=10)
        ax.add_patch(ligand_box)
        ax.text(0, 0, 'LIGAND\nStructure', ha='center', va='center',
               fontsize=14, fontweight='bold', zorder=11)
        print("[generate_advanced_interaction_plot] 使用配体占位符")
    
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
        print(f"[generate_advanced_interaction_plot] ✅ 2D相互作用图已保存: {output_path}")
    
    if show_plot:
        plt.show()
    else:
        plt.close()
    
    return output_path

def generate_interaction_network_plot(interactions_result=None, csv_path=None,
                                     output_path=None, show_plot=True, ligand_sdf=None):
    """
    生成交互网络图（使用matplotlib或networkx）

    参数:
        interactions_result: 相互作用分析结果
        csv_path: 或CSV文件路径
        output_path: 输出图片路径
        show_plot: 是否显示图片

    返回:
        output_path: 生成的图片路径
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.patches import FancyArrowPatch, Circle
    except ImportError:
        print("[generate_interaction_network_plot] 需要安装matplotlib")
        return None

    # 读取数据
    interactions = []
    csv_format = "standard"  # standard 或 advanced
    
    if csv_path:
        try:
            with open(csv_path, 'r', encoding='utf-8-sig') as f:
                reader = csv.DictReader(f)
                fieldnames = reader.fieldnames or []
                
                # 检测CSV格式
                if "Protein_Atom" in fieldnames and "Ligand_Atom" in fieldnames:
                    csv_format = "advanced"
                    print("[generate_interaction_network_plot] 检测到高级分析格式")
                
                for row in reader:
                    interactions.append(row)
        except Exception as e:
            print(f"[generate_interaction_network_plot] 读取CSV失败: {e}")
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
        print("[generate_interaction_network_plot] 没有相互作用数据")
        return None
    
    # 如果是高级分析格式，生成2D相互作用图（带配体结构）
    if csv_format == "advanced":
        return generate_advanced_interaction_plot(interactions, output_path, show_plot, ligand_sdf)

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
            lig_key = inter.get("Ligand_Residue") or inter.get("Atom1_Residue")
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

    # 添加图例
    legend_elements = []
    for itype, color in sorted(interaction_types_shown.items()):
        legend_elements.append(mpatches.Patch(color=color, label=itype))

    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper right', frameon=True,
                 facecolor='white', edgecolor='gray', fontsize=10)

    # 标题
    title = f"相互作用网络图\n({len(interactions)} 个相互作用, {len(ligand_nodes)} 个配体, {len(protein_nodes)} 个蛋白残基)"
    ax.text(0, 9.5, title, fontsize=14, fontweight='bold', ha='center')

    # 保存
    if output_path is None:
        output_path = "interaction_network.png"

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"[generate_interaction_network_plot] ✅ 网络图已生成: {output_path}")

    if show_plot:
        try:
            plt.show()
        except Exception:
            pass

    plt.close()

    return output_path

cmd.extend("generate_interaction_network_plot", generate_interaction_network_plot)

