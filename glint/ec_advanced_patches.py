# -*- coding: utf-8 -*-
"""
ec_advanced_patches.py
高级 EC 分析补丁Module - 解决 σ-hole、孤对电子、桥联水等方向性静电问题

主要功能：
1. SigmaHoleGenerator - 为卤素（Cl/Br/I）Add σ-hole 虚拟正电点
2. LonePairGenerator - 为羰基/胺Add孤对电子虚拟点
3. BridgingWaterFilter - 筛选结构性桥联水
4. EnhancedChargeCalculator - 整合所有补丁的增强电荷Calculator

跨平台兼容：
- Windows / macOS / Linux
- 纯 Python 实现，无外部二进制Dependencies
- 可选Dependencies：RDKit（推荐）、OpenBabel（备选）

Author: GLINT Team
"""

from __future__ import print_function
import os
import sys
import math
from typing import List, Dict, Tuple, Optional, Any, Union
from collections import defaultdict

# NumPy - 必需
try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False
    print("[ec_advanced_patches] ⚠️ NumPy not installed")

# RDKit - 推荐
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# SciPy - 可选（用于 KDTree 加速）
try:
    from scipy.spatial import cKDTree
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False


# ========== 常量定义 ==========

# σ-hole Parameters（基于文献Value）
# 参考: Clark et al., J. Mol. Model. 2007, 13, 291-296
SIGMA_HOLE_PARAMS = {
    'Cl': {
        'distance': 1.20,      # Å，虚拟点到卤素的距离
        'charge': +0.08,       # e，虚拟点电荷
        'halogen_adjust': -0.08,  # 相应减少卤素电荷以保持总电荷
        'vdw_radius': 0.3,     # 虚拟点的 vdW 半径
    },
    'Br': {
        'distance': 1.30,
        'charge': +0.10,
        'halogen_adjust': -0.10,
        'vdw_radius': 0.3,
    },
    'I': {
        'distance': 1.40,
        'charge': +0.12,
        'halogen_adjust': -0.12,
        'vdw_radius': 0.3,
    },
}

# 孤对电子Parameters（基于 TIP5P 水模型思想）
LONE_PAIR_PARAMS = {
    'carbonyl_O': {
        'distance': 0.70,      # Å，LP 到 O 的距离
        'charge': -0.20,       # e，每个 LP 的电荷
        'angle': 120.0,        # degrees，两个 LP 之间的angle
        'oxygen_adjust': +0.40,  # 相应增加 O 的电荷
    },
    'ether_O': {
        'distance': 0.60,
        'charge': -0.15,
        'angle': 109.5,
        'oxygen_adjust': +0.30,
    },
    'amine_N': {
        'distance': 0.60,
        'charge': -0.15,
        'angle': 109.5,
        'nitrogen_adjust': +0.15,
    },
}

# 桥联水筛选Parameters
BRIDGING_WATER_PARAMS = {
    'hbond_distance_max': 3.5,    # Å，氢Key最大距离
    'hbond_angle_min': 120.0,     # degrees，氢Key最小角degrees
    'min_protein_contacts': 2,    # 最少蛋白contact count
    'ligand_distance_max': 4.0,   # Å，到配体的最大距离（桥水）
}


class SigmaHoleGenerator:
    """
    σ-hole 虚拟点生成器
    
    在卤素（Cl/Br/I）的 C-X Key延长方向Add虚拟正电点，
    模拟卤素Key的方向性静电效应。
    
    原理：
    - 卤素原子在 C-X Key方向有电子密degrees缺失（σ-hole）
    - 这导致该方向呈现正电势，可与 Lewis 碱形成卤素Key
    - 标准点电荷模型无法表达这种各向异性
    
    usingMethod：
        generator = SigmaHoleGenerator()
        virtual_points = generator.generate(mol)
        # virtual_points Package含坐标、电荷、关联Atom information
    """
    
    def __init__(self, params: Dict = None):
        """
        Args:
            params: 自定义 σ-hole Parameters，默认using SIGMA_HOLE_PARAMS
        """
        self.params = params or SIGMA_HOLE_PARAMS
    
    def generate(self, mol: 'Chem.Mol', conformer_id: int = 0) -> Dict[str, Any]:
        """
        为分子生成 σ-hole 虚拟点
        
        Args:
            mol: RDKit 分子对象
            conformer_id: 构象 ID
            
        Returns:
            字典Package含：
            - coords: 虚拟点坐标 (N, 3)
            - charges: 虚拟点电荷 (N,)
            - radii: 虚拟点半径 (N,)
            - parent_atoms: 关联的卤素原子索引
            - halogen_charge_adjustments: 卤素电荷调整Value
        """
        if not RDKIT_AVAILABLE:
            print("[SigmaHoleGenerator] ⚠️ RDKit required")
            return self._empty_result()
        
        if not NUMPY_AVAILABLE:
            print("[SigmaHoleGenerator] ⚠️ NumPy required")
            return self._empty_result()
        
        conf = mol.GetConformer(conformer_id)
        
        virtual_coords = []
        virtual_charges = []
        virtual_radii = []
        parent_atoms = []
        halogen_adjustments = {}
        
        # 遍历所有Key，找 C-X Key
        for bond in mol.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            
            # 检查是否为 C-X Key（X = Cl/Br/I）
            carbon, halogen = None, None
            
            if atom1.GetSymbol() == 'C' and atom2.GetSymbol() in self.params:
                carbon, halogen = atom1, atom2
            elif atom2.GetSymbol() == 'C' and atom1.GetSymbol() in self.params:
                carbon, halogen = atom2, atom1
            
            if carbon is None or halogen is None:
                continue
            
            halogen_symbol = halogen.GetSymbol()
            params = self.params[halogen_symbol]
            
            # 获取原子坐标
            c_pos = np.array(conf.GetAtomPosition(carbon.GetIdx()))
            x_pos = np.array(conf.GetAtomPosition(halogen.GetIdx()))
            
            # 计算 C→X 方向向量
            direction = x_pos - c_pos
            direction = direction / np.linalg.norm(direction)
            
            # 在 X 延长方向放置虚拟点
            virtual_pos = x_pos + direction * params['distance']
            
            virtual_coords.append(virtual_pos)
            virtual_charges.append(params['charge'])
            virtual_radii.append(params['vdw_radius'])
            parent_atoms.append(halogen.GetIdx())
            
            # 记录卤素电荷调整
            halogen_idx = halogen.GetIdx()
            if halogen_idx not in halogen_adjustments:
                halogen_adjustments[halogen_idx] = 0.0
            halogen_adjustments[halogen_idx] += params['halogen_adjust']
        
        result = {
            'coords': np.array(virtual_coords) if virtual_coords else np.zeros((0, 3)),
            'charges': np.array(virtual_charges) if virtual_charges else np.zeros(0),
            'radii': np.array(virtual_radii) if virtual_radii else np.zeros(0),
            'parent_atoms': parent_atoms,
            'halogen_charge_adjustments': halogen_adjustments,
            'n_points': len(virtual_coords),
        }
        
        if result['n_points'] > 0:
            print(f"[SigmaHoleGenerator] Generated {result['n_points']} σ-hole virtual points")
            for idx, (coord, charge) in enumerate(zip(virtual_coords, virtual_charges)):
                parent_idx = parent_atoms[idx]
                parent_symbol = mol.GetAtomWithIdx(parent_idx).GetSymbol()
                print(f"  - {parent_symbol}[{parent_idx}]: pos={coord}, charge={charge:+.3f}")
        
        return result
    
    def generate_from_coords(self, coords: np.ndarray, elements: List[str],
                            bonds: List[Tuple[int, int]]) -> Dict[str, Any]:
        """
        从原始坐标生成 σ-hole 虚拟点（不Dependencies RDKit）
        
        Args:
            coords: 原子坐标 (N, 3)
            elements: 元素符号列表
            bonds: Key连接列表 [(i, j), ...]
            
        Returns:
            同 generate() Method
        """
        if not NUMPY_AVAILABLE:
            return self._empty_result()
        
        virtual_coords = []
        virtual_charges = []
        virtual_radii = []
        parent_atoms = []
        halogen_adjustments = {}
        
        for i, j in bonds:
            elem_i, elem_j = elements[i].upper(), elements[j].upper()
            
            carbon_idx, halogen_idx = None, None
            
            if elem_i == 'C' and elem_j in self.params:
                carbon_idx, halogen_idx = i, j
            elif elem_j == 'C' and elem_i in self.params:
                carbon_idx, halogen_idx = j, i
            
            if carbon_idx is None:
                continue
            
            halogen_symbol = elements[halogen_idx].upper()
            params = self.params[halogen_symbol]
            
            c_pos = coords[carbon_idx]
            x_pos = coords[halogen_idx]
            
            direction = x_pos - c_pos
            direction = direction / np.linalg.norm(direction)
            
            virtual_pos = x_pos + direction * params['distance']
            
            virtual_coords.append(virtual_pos)
            virtual_charges.append(params['charge'])
            virtual_radii.append(params['vdw_radius'])
            parent_atoms.append(halogen_idx)
            
            if halogen_idx not in halogen_adjustments:
                halogen_adjustments[halogen_idx] = 0.0
            halogen_adjustments[halogen_idx] += params['halogen_adjust']
        
        return {
            'coords': np.array(virtual_coords) if virtual_coords else np.zeros((0, 3)),
            'charges': np.array(virtual_charges) if virtual_charges else np.zeros(0),
            'radii': np.array(virtual_radii) if virtual_radii else np.zeros(0),
            'parent_atoms': parent_atoms,
            'halogen_charge_adjustments': halogen_adjustments,
            'n_points': len(virtual_coords),
        }
    
    def _empty_result(self) -> Dict[str, Any]:
        """Return空Results"""
        return {
            'coords': np.zeros((0, 3)) if NUMPY_AVAILABLE else [],
            'charges': np.zeros(0) if NUMPY_AVAILABLE else [],
            'radii': np.zeros(0) if NUMPY_AVAILABLE else [],
            'parent_atoms': [],
            'halogen_charge_adjustments': {},
            'n_points': 0,
        }


class LonePairGenerator:
    """
    孤对电子虚拟点生成器
    
    为羰基氧、醚氧、胺氮等Add孤对电子虚拟点，
    模拟氢Key受体的方向性。
    
    原理：
    - 羰基 C=O 的氧有两个孤对电子，位于 C=O Key两侧
    - 标准点电荷只给 O 一个负电荷，无法表达方向性
    - Add LP 虚拟点可改善氢Key方向性预测
    
    usingMethod：
        generator = LonePairGenerator()
        virtual_points = generator.generate(mol)
    """
    
    def __init__(self, params: Dict = None):
        """
        Args:
            params: 自定义孤对电子Parameters
        """
        self.params = params or LONE_PAIR_PARAMS
    
    def generate(self, mol: 'Chem.Mol', conformer_id: int = 0) -> Dict[str, Any]:
        """
        为分子生成孤对电子虚拟点
        
        Args:
            mol: RDKit 分子对象
            conformer_id: 构象 ID
            
        Returns:
            字典Package含虚拟点Information
        """
        if not RDKIT_AVAILABLE:
            print("[LonePairGenerator] ⚠️ RDKit required")
            return self._empty_result()
        
        if not NUMPY_AVAILABLE:
            print("[LonePairGenerator] ⚠️ NumPy required")
            return self._empty_result()
        
        conf = mol.GetConformer(conformer_id)
        
        virtual_coords = []
        virtual_charges = []
        virtual_radii = []
        parent_atoms = []
        atom_adjustments = {}
        
        # 遍历所有原子，找羰基 O、醚 O、胺 N
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            idx = atom.GetIdx()
            
            if symbol == 'O':
                # 检查是否为羰基 O (C=O)
                is_carbonyl = False
                is_ether = False
                carbon_neighbor = None
                
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C':
                        bond = mol.GetBondBetweenAtoms(idx, neighbor.GetIdx())
                        if bond.GetBondType() == Chem.BondType.DOUBLE:
                            is_carbonyl = True
                            carbon_neighbor = neighbor
                            break
                        else:
                            is_ether = True
                            carbon_neighbor = neighbor
                
                if is_carbonyl and carbon_neighbor is not None:
                    lp_points = self._generate_carbonyl_lp(
                        conf, atom, carbon_neighbor, self.params['carbonyl_O']
                    )
                    virtual_coords.extend(lp_points['coords'])
                    virtual_charges.extend(lp_points['charges'])
                    virtual_radii.extend(lp_points['radii'])
                    parent_atoms.extend([idx] * len(lp_points['coords']))
                    atom_adjustments[idx] = self.params['carbonyl_O']['oxygen_adjust']
                
                elif is_ether and carbon_neighbor is not None:
                    # 醚氧的孤对电子（简化处理）
                    pass  # 可选实现
            
            elif symbol == 'N':
                # 检查是否为胺 N（sp3 杂化，有孤对）
                # 简化：只处理伯胺/仲胺
                if atom.GetHybridization() == Chem.HybridizationType.SP3:
                    # 可选实现
                    pass
        
        result = {
            'coords': np.array(virtual_coords) if virtual_coords else np.zeros((0, 3)),
            'charges': np.array(virtual_charges) if virtual_charges else np.zeros(0),
            'radii': np.array(virtual_radii) if virtual_radii else np.zeros(0),
            'parent_atoms': parent_atoms,
            'atom_charge_adjustments': atom_adjustments,
            'n_points': len(virtual_coords),
        }
        
        if result['n_points'] > 0:
            print(f"[LonePairGenerator] Generated {result['n_points']} lone pair virtual points")
        
        return result
    
    def _generate_carbonyl_lp(self, conf: 'Chem.Conformer', oxygen: 'Chem.Atom',
                              carbon: 'Chem.Atom', params: Dict) -> Dict[str, Any]:
        """
        为羰基氧生成两个孤对电子虚拟点
        
        LP1 和 LP2 位于 C=O Key两侧，与 C=O Key成 120° 角
        """
        o_pos = np.array(conf.GetAtomPosition(oxygen.GetIdx()))
        c_pos = np.array(conf.GetAtomPosition(carbon.GetIdx()))
        
        # C→O 方向
        co_vec = o_pos - c_pos
        co_vec = co_vec / np.linalg.norm(co_vec)
        
        # 找一个垂直于 C=O 的向量
        # using羰基碳的另一个邻居来确定平面
        perp_vec = None
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetIdx() != oxygen.GetIdx():
                n_pos = np.array(conf.GetAtomPosition(neighbor.GetIdx()))
                cn_vec = n_pos - c_pos
                cn_vec = cn_vec / np.linalg.norm(cn_vec)
                perp_vec = np.cross(co_vec, cn_vec)
                perp_vec = perp_vec / np.linalg.norm(perp_vec)
                break
        
        if perp_vec is None:
            # 如果找不到邻居，using任意垂直向量
            if abs(co_vec[0]) < 0.9:
                perp_vec = np.cross(co_vec, np.array([1, 0, 0]))
            else:
                perp_vec = np.cross(co_vec, np.array([0, 1, 0]))
            perp_vec = perp_vec / np.linalg.norm(perp_vec)
        
        # 计算 LP 方向（与 C=O 成 120° 角）
        angle_rad = math.radians(params['angle'] / 2)  # 每个 LP 与 C=O 反方向的angle
        
        # LP 基础方向（指向 C 的方向）
        lp_base = -co_vec
        
        # 旋转得到两个 LP 方向
        lp1_dir = self._rotate_vector(lp_base, perp_vec, angle_rad)
        lp2_dir = self._rotate_vector(lp_base, perp_vec, -angle_rad)
        
        # 计算 LP 位置
        lp1_pos = o_pos + lp1_dir * params['distance']
        lp2_pos = o_pos + lp2_dir * params['distance']
        
        return {
            'coords': [lp1_pos, lp2_pos],
            'charges': [params['charge'], params['charge']],
            'radii': [0.3, 0.3],
        }
    
    def _rotate_vector(self, vec: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
        """using Rodrigues 旋转公式绕轴旋转向量"""
        axis = axis / np.linalg.norm(axis)
        cos_a = math.cos(angle)
        sin_a = math.sin(angle)
        
        rotated = (vec * cos_a + 
                   np.cross(axis, vec) * sin_a + 
                   axis * np.dot(axis, vec) * (1 - cos_a))
        
        return rotated / np.linalg.norm(rotated)
    
    def _empty_result(self) -> Dict[str, Any]:
        """Return空Results"""
        return {
            'coords': np.zeros((0, 3)) if NUMPY_AVAILABLE else [],
            'charges': np.zeros(0) if NUMPY_AVAILABLE else [],
            'radii': np.zeros(0) if NUMPY_AVAILABLE else [],
            'parent_atoms': [],
            'atom_charge_adjustments': {},
            'n_points': 0,
        }


class BridgingWaterFilter:
    """
    桥联水筛选器
    
    从晶体结构中筛选出结构性重要的水分子：
    1. 与蛋白形成多个氢Key的水
    2. 同时连接蛋白和配体的桥水
    
    原理：
    - 结构性水分子对结合亲和力有重要贡献
    - 简单Delete所有水会丢失这些Information
    - 保留关Key水分子可改善 EC 分析准确性
    
    usingMethod：
        filter = BridgingWaterFilter()
        waters_to_keep = filter.filter_pdb('complex.pdb', 'LIG')
    """
    
    def __init__(self, params: Dict = None):
        """
        Args:
            params: 自定义筛选Parameters
        """
        self.params = params or BRIDGING_WATER_PARAMS
    
    def filter_pdb(self, pdb_file: str, ligand_resname: str = None) -> Dict[str, Any]:
        """
        从 PDB File筛选桥联水
        
        Args:
            pdb_file: PDB FilePath
            ligand_resname: 配体残基名（用于识别桥水）
            
        Returns:
            字典Package含：
            - waters_to_keep: 应保留的水分子列表 [(chain, resnum), ...]
            - bridging_waters: 桥水列表
            - protein_bound_waters: 蛋白结合水列表
            - statistics: 统计Information
        """
        if not NUMPY_AVAILABLE:
            print("[BridgingWaterFilter] ⚠️ NumPy required")
            return self._empty_result()
        
        # 解析 PDB File
        protein_atoms = []
        ligand_atoms = []
        water_atoms = []
        
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if not line.startswith(('ATOM', 'HETATM')):
                        continue
                    
                    atom_name = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21:22].strip() or 'A'
                    resnum = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                    
                    atom_info = {
                        'name': atom_name,
                        'resname': resname,
                        'chain': chain,
                        'resnum': resnum,
                        'coords': np.array([x, y, z]),
                        'element': element,
                    }
                    
                    if resname in ('HOH', 'WAT', 'H2O', 'TIP', 'TIP3', 'SPC'):
                        water_atoms.append(atom_info)
                    elif ligand_resname and resname == ligand_resname:
                        ligand_atoms.append(atom_info)
                    elif resname in self._standard_aa():
                        protein_atoms.append(atom_info)
        
        except Exception as e:
            print(f"[BridgingWaterFilter] Error reading PDB: {e}")
            return self._empty_result()
        
        if not water_atoms:
            print("[BridgingWaterFilter] No water molecules found")
            return self._empty_result()
        
        # 获取蛋白和配体的重原子坐标（用于距离计算）
        protein_coords = np.array([a['coords'] for a in protein_atoms 
                                   if a['element'] not in ('H', 'D')])
        ligand_coords = np.array([a['coords'] for a in ligand_atoms 
                                  if a['element'] not in ('H', 'D')]) if ligand_atoms else None
        
        # 构建 KDTree 加速距离查询
        if SCIPY_AVAILABLE and len(protein_coords) > 0:
            protein_tree = cKDTree(protein_coords)
        else:
            protein_tree = None
        
        if SCIPY_AVAILABLE and ligand_coords is not None and len(ligand_coords) > 0:
            ligand_tree = cKDTree(ligand_coords)
        else:
            ligand_tree = None
        
        # 分析每个水分子
        waters_to_keep = []
        bridging_waters = []
        protein_bound_waters = []
        
        # 按残基分组水分子
        water_residues = defaultdict(list)
        for atom in water_atoms:
            key = (atom['chain'], atom['resnum'])
            water_residues[key].append(atom)
        
        for (chain, resnum), atoms in water_residues.items():
            # 获取水氧原子坐标
            o_atom = None
            for atom in atoms:
                if atom['element'] == 'O' or atom['name'] == 'O':
                    o_atom = atom
                    break
            
            if o_atom is None:
                continue
            
            water_coord = o_atom['coords']
            
            # 计算到蛋白的距离
            protein_contacts = 0
            if protein_tree is not None:
                distances, _ = protein_tree.query(water_coord, k=min(10, len(protein_coords)))
                if isinstance(distances, float):
                    distances = [distances]
                protein_contacts = sum(1 for d in distances 
                                       if d <= self.params['hbond_distance_max'])
            elif len(protein_coords) > 0:
                distances = np.linalg.norm(protein_coords - water_coord, axis=1)
                protein_contacts = np.sum(distances <= self.params['hbond_distance_max'])
            
            # 计算到配体的距离
            ligand_contact = False
            if ligand_tree is not None:
                dist, _ = ligand_tree.query(water_coord)
                ligand_contact = dist <= self.params['ligand_distance_max']
            elif ligand_coords is not None and len(ligand_coords) > 0:
                distances = np.linalg.norm(ligand_coords - water_coord, axis=1)
                ligand_contact = np.min(distances) <= self.params['ligand_distance_max']
            
            # 判断是否保留
            is_bridging = protein_contacts >= 1 and ligand_contact
            is_protein_bound = protein_contacts >= self.params['min_protein_contacts']
            
            if is_bridging:
                bridging_waters.append((chain, resnum))
                waters_to_keep.append((chain, resnum))
            elif is_protein_bound:
                protein_bound_waters.append((chain, resnum))
                waters_to_keep.append((chain, resnum))
        
        result = {
            'waters_to_keep': waters_to_keep,
            'bridging_waters': bridging_waters,
            'protein_bound_waters': protein_bound_waters,
            'statistics': {
                'total_waters': len(water_residues),
                'kept_waters': len(waters_to_keep),
                'bridging_count': len(bridging_waters),
                'protein_bound_count': len(protein_bound_waters),
            }
        }
        
        print(f"[BridgingWaterFilter] Water analysis:")
        print(f"  Total waters: {result['statistics']['total_waters']}")
        print(f"  Bridging waters: {result['statistics']['bridging_count']}")
        print(f"  Protein-bound waters: {result['statistics']['protein_bound_count']}")
        print(f"  Waters to keep: {result['statistics']['kept_waters']}")
        
        return result
    
    def write_filtered_pdb(self, input_pdb: str, output_pdb: str,
                          waters_to_keep: List[Tuple[str, int]],
                          keep_all_protein: bool = True,
                          keep_ligand: bool = True,
                          ligand_resname: str = None) -> bool:
        """
        写入筛选后的 PDB File
        
        Args:
            input_pdb: 输入 PDB File
            output_pdb: 输出 PDB File
            waters_to_keep: 要保留的水分子列表
            keep_all_protein: 是否保留所有蛋白原子
            keep_ligand: 是否保留配体
            ligand_resname: 配体残基名
            
        Returns:
            是否Success
        """
        waters_set = set(waters_to_keep)
        
        try:
            with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
                for line in f_in:
                    if not line.startswith(('ATOM', 'HETATM')):
                        f_out.write(line)
                        continue
                    
                    resname = line[17:20].strip()
                    chain = line[21:22].strip() or 'A'
                    resnum = int(line[22:26])
                    
                    # 水分子
                    if resname in ('HOH', 'WAT', 'H2O', 'TIP', 'TIP3', 'SPC'):
                        if (chain, resnum) in waters_set:
                            f_out.write(line)
                    # 蛋白
                    elif resname in self._standard_aa():
                        if keep_all_protein:
                            f_out.write(line)
                    # 配体
                    elif ligand_resname and resname == ligand_resname:
                        if keep_ligand:
                            f_out.write(line)
                    # 其他（离子等）
                    else:
                        f_out.write(line)
            
            print(f"[BridgingWaterFilter] Wrote filtered PDB: {output_pdb}")
            return True
            
        except Exception as e:
            print(f"[BridgingWaterFilter] Error writing PDB: {e}")
            return False
    
    def _standard_aa(self) -> set:
        """标准氨基酸残基名"""
        return {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
            'HIE', 'HID', 'HIP', 'CYX',  # 常见变体
        }
    
    def _empty_result(self) -> Dict[str, Any]:
        """Return空Results"""
        return {
            'waters_to_keep': [],
            'bridging_waters': [],
            'protein_bound_waters': [],
            'statistics': {
                'total_waters': 0,
                'kept_waters': 0,
                'bridging_count': 0,
                'protein_bound_count': 0,
            }
        }


class EnhancedChargeCalculator:
    """
    增强电荷Calculator
    
    整合所有高级补丁，提供统一的电荷计算接口：
    1. 基础 Gasteiger 电荷
    2. σ-hole 虚拟点
    3. 孤对电子虚拟点
    
    usingMethod：
        calc = EnhancedChargeCalculator(use_sigma_holes=True, use_lone_pairs=True)
        result = calc.calculate(mol)
        # result Package含原子电荷 + 虚拟点电荷
    """
    
    def __init__(self,
                 use_sigma_holes: bool = True,
                 use_lone_pairs: bool = False,  # 默认Close，因为实现较复杂
                 dielectric: float = 4.0,
                 sigma_hole_params: Dict = None,
                 lone_pair_params: Dict = None):
        """
        Args:
            use_sigma_holes: 是否Add σ-hole 虚拟点
            use_lone_pairs: 是否Add孤对电子虚拟点
            dielectric: 介电常数
            sigma_hole_params: 自定义 σ-hole Parameters
            lone_pair_params: 自定义孤对电子Parameters
        """
        self.use_sigma_holes = use_sigma_holes
        self.use_lone_pairs = use_lone_pairs
        self.dielectric = dielectric
        
        self.sigma_hole_gen = SigmaHoleGenerator(sigma_hole_params) if use_sigma_holes else None
        self.lone_pair_gen = LonePairGenerator(lone_pair_params) if use_lone_pairs else None
    
    def calculate(self, mol: 'Chem.Mol', conformer_id: int = 0) -> Dict[str, Any]:
        """
        计算增强电荷（原子 + 虚拟点）
        
        Args:
            mol: RDKit 分子对象
            conformer_id: 构象 ID
            
        Returns:
            字典Package含：
            - atom_coords: 原子坐标 (N, 3)
            - atom_charges: 原子电荷 (N,)
            - atom_radii: 原子半径 (N,)
            - virtual_coords: 虚拟点坐标 (M, 3)
            - virtual_charges: 虚拟点电荷 (M,)
            - virtual_radii: 虚拟点半径 (M,)
            - all_coords: 所有点坐标 (N+M, 3)
            - all_charges: 所有点电荷 (N+M,)
            - all_radii: 所有点半径 (N+M,)
        """
        if not RDKIT_AVAILABLE:
            print("[EnhancedChargeCalculator] ⚠️ RDKit required")
            return self._empty_result()
        
        if not NUMPY_AVAILABLE:
            print("[EnhancedChargeCalculator] ⚠️ NumPy required")
            return self._empty_result()
        
        from rdkit.Chem import AllChem
        
        # 计算 Gasteiger 电荷
        AllChem.ComputeGasteigerCharges(mol)
        
        conf = mol.GetConformer(conformer_id)
        
        atom_coords = []
        atom_charges = []
        atom_radii = []
        
        # VDW 半径
        vdw_radii = {
            'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
            'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98,
        }
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            atom_coords.append([pos.x, pos.y, pos.z])
            
            charge = atom.GetDoubleProp('_GasteigerCharge')
            if np.isnan(charge):
                charge = 0.0
            atom_charges.append(charge)
            
            symbol = atom.GetSymbol()
            radius = vdw_radii.get(symbol, 1.70)
            atom_radii.append(radius)
        
        atom_coords = np.array(atom_coords)
        atom_charges = np.array(atom_charges)
        atom_radii = np.array(atom_radii)
        
        # 收集虚拟点
        virtual_coords = []
        virtual_charges = []
        virtual_radii = []
        
        # Add σ-hole 虚拟点
        if self.use_sigma_holes and self.sigma_hole_gen:
            sigma_result = self.sigma_hole_gen.generate(mol, conformer_id)
            
            if sigma_result['n_points'] > 0:
                virtual_coords.append(sigma_result['coords'])
                virtual_charges.append(sigma_result['charges'])
                virtual_radii.append(sigma_result['radii'])
                
                # 调整卤素原子电荷
                for atom_idx, adjustment in sigma_result['halogen_charge_adjustments'].items():
                    atom_charges[atom_idx] += adjustment
        
        # Add孤对电子虚拟点
        if self.use_lone_pairs and self.lone_pair_gen:
            lp_result = self.lone_pair_gen.generate(mol, conformer_id)
            
            if lp_result['n_points'] > 0:
                virtual_coords.append(lp_result['coords'])
                virtual_charges.append(lp_result['charges'])
                virtual_radii.append(lp_result['radii'])
                
                # 调整原子电荷
                for atom_idx, adjustment in lp_result['atom_charge_adjustments'].items():
                    atom_charges[atom_idx] += adjustment
        
        # 合并虚拟点
        if virtual_coords:
            virtual_coords = np.vstack(virtual_coords)
            virtual_charges = np.concatenate(virtual_charges)
            virtual_radii = np.concatenate(virtual_radii)
        else:
            virtual_coords = np.zeros((0, 3))
            virtual_charges = np.zeros(0)
            virtual_radii = np.zeros(0)
        
        # 合并所有点
        all_coords = np.vstack([atom_coords, virtual_coords]) if len(virtual_coords) > 0 else atom_coords
        all_charges = np.concatenate([atom_charges, virtual_charges]) if len(virtual_charges) > 0 else atom_charges
        all_radii = np.concatenate([atom_radii, virtual_radii]) if len(virtual_radii) > 0 else atom_radii
        
        result = {
            'atom_coords': atom_coords,
            'atom_charges': atom_charges,
            'atom_radii': atom_radii,
            'virtual_coords': virtual_coords,
            'virtual_charges': virtual_charges,
            'virtual_radii': virtual_radii,
            'all_coords': all_coords,
            'all_charges': all_charges,
            'all_radii': all_radii,
            'n_atoms': len(atom_coords),
            'n_virtual': len(virtual_coords),
            'total_charge': float(np.sum(all_charges)),
        }
        
        print(f"[EnhancedChargeCalculator] Calculated charges:")
        print(f"  Atoms: {result['n_atoms']}")
        print(f"  Virtual points: {result['n_virtual']}")
        print(f"  Total charge: {result['total_charge']:.4f}")
        
        return result
    
    def calculate_potential(self, mol: 'Chem.Mol', points: np.ndarray,
                           conformer_id: int = 0) -> np.ndarray:
        """
        计算增强静电势（Package含虚拟点贡献）
        
        Args:
            mol: RDKit 分子对象
            points: 评估点坐标 (M, 3)
            conformer_id: 构象 ID
            
        Returns:
            静电势Value (M,)
        """
        charge_result = self.calculate(mol, conformer_id)
        
        all_coords = charge_result['all_coords']
        all_charges = charge_result['all_charges']
        
        # 计算 Coulomb 势能
        potentials = np.zeros(len(points))
        
        for i, point in enumerate(points):
            distances = np.linalg.norm(all_coords - point, axis=1)
            distances = np.maximum(distances, 0.1)  # 避免除零
            potentials[i] = np.sum(all_charges / (self.dielectric * distances))
        
        # 转换为 kT/e 单位
        potentials *= 332.0637 / 0.593
        
        return potentials
    
    def _empty_result(self) -> Dict[str, Any]:
        """Return空Results"""
        empty_array = np.zeros((0, 3)) if NUMPY_AVAILABLE else []
        empty_1d = np.zeros(0) if NUMPY_AVAILABLE else []
        
        return {
            'atom_coords': empty_array,
            'atom_charges': empty_1d,
            'atom_radii': empty_1d,
            'virtual_coords': empty_array,
            'virtual_charges': empty_1d,
            'virtual_radii': empty_1d,
            'all_coords': empty_array,
            'all_charges': empty_1d,
            'all_radii': empty_1d,
            'n_atoms': 0,
            'n_virtual': 0,
            'total_charge': 0.0,
        }


# ========== PQR FileTool ==========

def write_enhanced_pqr(output_file: str,
                       atom_coords: np.ndarray,
                       atom_charges: np.ndarray,
                       atom_radii: np.ndarray,
                       atom_names: List[str],
                       atom_elements: List[str],
                       virtual_coords: np.ndarray = None,
                       virtual_charges: np.ndarray = None,
                       virtual_radii: np.ndarray = None,
                       virtual_names: List[str] = None) -> bool:
    """
    写入增强 PQR File（Package含虚拟点）
    
    Args:
        output_file: 输出FilePath
        atom_coords: 原子坐标
        atom_charges: 原子电荷
        atom_radii: 原子半径
        atom_names: 原子Name
        atom_elements: 元素符号
        virtual_coords: 虚拟点坐标
        virtual_charges: 虚拟点电荷
        virtual_radii: 虚拟点半径
        virtual_names: 虚拟点Name
        
    Returns:
        是否Success
    """
    try:
        with open(output_file, 'w') as f:
            f.write("REMARK Enhanced PQR file with virtual points\n")
            f.write("REMARK Generated by GLINT ec_advanced_patches\n")
            
            atom_idx = 1
            
            # 写入原子
            for i in range(len(atom_coords)):
                name = atom_names[i] if i < len(atom_names) else 'X'
                elem = atom_elements[i] if i < len(atom_elements) else 'X'
                coord = atom_coords[i]
                charge = atom_charges[i]
                radius = atom_radii[i]
                
                line = (f"ATOM  {atom_idx:5d} {name:4s} LIG A   1    "
                       f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                       f"{charge:8.4f}{radius:7.4f}          {elem:2s}\n")
                f.write(line)
                atom_idx += 1
            
            # 写入虚拟点
            if virtual_coords is not None and len(virtual_coords) > 0:
                for i in range(len(virtual_coords)):
                    name = virtual_names[i] if virtual_names and i < len(virtual_names) else 'VP'
                    coord = virtual_coords[i]
                    charge = virtual_charges[i] if virtual_charges is not None else 0.0
                    radius = virtual_radii[i] if virtual_radii is not None else 0.3
                    
                    line = (f"HETATM{atom_idx:5d} {name:4s} VPT A   2    "
                           f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                           f"{charge:8.4f}{radius:7.4f}          X\n")
                    f.write(line)
                    atom_idx += 1
            
            f.write("END\n")
        
        print(f"[write_enhanced_pqr] Wrote {atom_idx - 1} entries to {output_file}")
        return True
        
    except Exception as e:
        print(f"[write_enhanced_pqr] Error: {e}")
        return False


# ========== 便捷Function ==========

def apply_all_patches(mol: 'Chem.Mol',
                      use_sigma_holes: bool = True,
                      use_lone_pairs: bool = False,
                      conformer_id: int = 0) -> Dict[str, Any]:
    """
    Apply所有高级补丁的便捷Function
    
    Args:
        mol: RDKit 分子对象
        use_sigma_holes: 是否Add σ-hole
        use_lone_pairs: 是否Add孤对电子
        conformer_id: 构象 ID
        
    Returns:
        增强电荷计算Results
    """
    calc = EnhancedChargeCalculator(
        use_sigma_holes=use_sigma_holes,
        use_lone_pairs=use_lone_pairs
    )
    return calc.calculate(mol, conformer_id)


def filter_and_prepare_complex(pdb_file: str,
                               output_dir: str,
                               ligand_resname: str = None,
                               keep_bridging_waters: bool = True) -> Dict[str, str]:
    """
    筛选桥联水并准备复合物结构
    
    Args:
        pdb_file: 输入 PDB File
        output_dir: 输出Directory
        ligand_resname: 配体残基名
        keep_bridging_waters: 是否保留桥联水
        
    Returns:
        输出FilePath字典
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    result = {
        'protein_pdb': os.path.join(output_dir, 'protein.pdb'),
        'ligand_pdb': os.path.join(output_dir, 'ligand.pdb'),
        'complex_with_waters_pdb': os.path.join(output_dir, 'complex_with_waters.pdb'),
    }
    
    if keep_bridging_waters:
        water_filter = BridgingWaterFilter()
        water_result = water_filter.filter_pdb(pdb_file, ligand_resname)
        
        water_filter.write_filtered_pdb(
            pdb_file,
            result['complex_with_waters_pdb'],
            water_result['waters_to_keep'],
            keep_all_protein=True,
            keep_ligand=True,
            ligand_resname=ligand_resname
        )
        
        result['bridging_waters'] = water_result['bridging_waters']
        result['water_statistics'] = water_result['statistics']
    
    return result


# ========== ModuleInformation ==========

def print_module_info():
    """PrintModuleInformation"""
    print("="*60)
    print("GLINT EC Advanced Patches Module")
    print("="*60)
    print("\nFeatures:")
    print("  ✓ σ-hole virtual points for halogens (Cl/Br/I)")
    print("  ✓ Lone pair virtual points for carbonyl O")
    print("  ✓ Bridging water filter")
    print("  ✓ Enhanced charge calculator")
    print("\nDependencies:")
    print(f"  NumPy: {'✅' if NUMPY_AVAILABLE else '❌'}")
    print(f"  RDKit: {'✅' if RDKIT_AVAILABLE else '❌'}")
    print(f"  SciPy: {'✅ (optional)' if SCIPY_AVAILABLE else '❌ (optional)'}")
    print("\nUsage:")
    print("  from glint.ec_advanced_patches import EnhancedChargeCalculator")
    print("  calc = EnhancedChargeCalculator(use_sigma_holes=True)")
    print("  result = calc.calculate(mol)")
    print("="*60)


if __name__ == '__main__':
    print_module_info()