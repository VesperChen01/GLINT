# -*- coding: utf-8 -*-
"""
feature_extractor.py
分子特征提取器

基于 PLIP 设计理念，一次遍历提取所有特征并缓存
减少重复计算，提升性能
"""

from __future__ import annotations
import math
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict


# ============================================================================
# 辅助函数（从 interaction_analyzer.py 复用）
# ============================================================================

def get_element_from_atom_name(atom_name: str) -> str:
    """从原子名提取元素符号"""
    elem = ''.join(c for c in atom_name if c.isalpha())
    if len(elem) == 0:
        return ""
    # 常见双字母元素
    two_letter = {"BR", "CL", "FE", "MG", "CA", "ZN", "CU", "MN", "CO", "NI"}
    if len(elem) >= 2 and elem[:2].upper() in two_letter:
        return elem[:2].upper()
    return elem[0].upper()


def distance(coord1: Tuple[float, float, float], 
             coord2: Tuple[float, float, float]) -> float:
    """计算两点间距离"""
    return math.sqrt(sum((coord1[i] - coord2[i])**2 for i in range(3)))


def centroid(coords: List[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    """计算坐标质心"""
    n = len(coords)
    if n == 0:
        return (0.0, 0.0, 0.0)
    return tuple(sum(c[i] for c in coords)/n for i in range(3))


def normal_vector(a: Tuple[float, float, float], 
                  b: Tuple[float, float, float], 
                  c: Tuple[float, float, float]) -> Optional[List[float]]:
    """计算三点确定平面的法向量"""
    ab = [b[i]-a[i] for i in range(3)]
    ac = [c[i]-a[i] for i in range(3)]
    n = [ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0]]
    norm = math.sqrt(sum(x*x for x in n))
    return [x/norm for x in n] if norm != 0 else None


# ============================================================================
# 分子特征提取器
# ============================================================================

class MolecularFeatureExtractor:
    """
    分子特征提取器（参考 PLIP）
    
    一次遍历提取所有化学特征并缓存，避免重复计算
    
    支持的特征：
    - 氢键供体/受体
    - 疏水原子
    - 芳香环
    - 带电基团（正/负）
    - 金属离子
    - 卤素原子
    
    Attributes:
        atoms (List[Tuple]): 原子列表，格式 (chain, resn, resi, atom_name, coords)
    """
    
    def __init__(self, atoms: List[Tuple]):
        """
        初始化特征提取器
        
        Args:
            atoms: 原子列表 [(chain, resn, resi, atom_name, (x,y,z)), ...]
        """
        self.atoms = atoms
        
        # 缓存标志
        self._hbond_donors: Optional[List[Tuple]] = None
        self._hbond_acceptors: Optional[List[Tuple]] = None
        self._hydrophobic_atoms: Optional[List[Tuple]] = None
        self._aromatic_rings: Optional[List[Dict]] = None
        self._charged_groups: Optional[List[Dict]] = None
        self._metal_ions: Optional[List[Tuple]] = None
        self._halogen_atoms: Optional[List[Tuple]] = None
        
        # 残基分组缓存
        self._residues: Optional[Dict[Tuple, List[Tuple]]] = None
    
    @property
    def residues(self) -> Dict[Tuple[str, str, str], List[Tuple]]:
        """按残基分组的原子（缓存）"""
        if self._residues is None:
            self._residues = self._group_by_residue()
        return self._residues
    
    @property
    def hbond_donors(self) -> List[Tuple]:
        """氢键供体原子（缓存）"""
        if self._hbond_donors is None:
            self._hbond_donors = self._extract_hbond_donors()
        return self._hbond_donors
    
    @property
    def hbond_acceptors(self) -> List[Tuple]:
        """氢键受体原子（缓存）"""
        if self._hbond_acceptors is None:
            self._hbond_acceptors = self._extract_hbond_acceptors()
        return self._hbond_acceptors
    
    @property
    def hydrophobic_atoms(self) -> List[Tuple]:
        """疏水原子（缓存）"""
        if self._hydrophobic_atoms is None:
            self._hydrophobic_atoms = self._extract_hydrophobic_atoms()
        return self._hydrophobic_atoms
    
    @property
    def aromatic_rings(self) -> List[Dict]:
        """芳香环（缓存）"""
        if self._aromatic_rings is None:
            self._aromatic_rings = self._extract_aromatic_rings()
        return self._aromatic_rings
    
    @property
    def charged_groups(self) -> List[Dict]:
        """带电基团（缓存）"""
        if self._charged_groups is None:
            self._charged_groups = self._extract_charged_groups()
        return self._charged_groups
    
    @property
    def metal_ions(self) -> List[Tuple]:
        """金属离子（缓存）"""
        if self._metal_ions is None:
            self._metal_ions = self._extract_metal_ions()
        return self._metal_ions
    
    @property
    def halogen_atoms(self) -> List[Tuple]:
        """卤素原子（缓存）"""
        if self._halogen_atoms is None:
            self._halogen_atoms = self._extract_halogen_atoms()
        return self._halogen_atoms
    
    # ========================================================================
    # 内部提取方法
    # ========================================================================
    
    def _group_by_residue(self) -> Dict[Tuple[str, str, str], List[Tuple]]:
        """按残基分组原子"""
        residues = defaultdict(list)
        for atom in self.atoms:
            res_key = (atom[0], atom[1], atom[2])  # (chain, resn, resi)
            residues[res_key].append(atom)
        return dict(residues)
    
    def _extract_hbond_donors(self) -> List[Tuple]:
        """
        提取氢键供体
        
        标准：N, O, S 原子（可能连接氢）
        """
        donors = []
        donor_elements = {'N', 'O', 'S'}
        
        for atom in self.atoms:
            elem = get_element_from_atom_name(atom[3])
            if elem in donor_elements:
                donors.append(atom)
        
        return donors
    
    def _extract_hbond_acceptors(self) -> List[Tuple]:
        """
        提取氢键受体
        
        标准：N, O, S, F, Cl 原子（有孤对电子）
        """
        acceptors = []
        acceptor_elements = {'N', 'O', 'S', 'F', 'CL'}
        
        for atom in self.atoms:
            elem = get_element_from_atom_name(atom[3])
            if elem in acceptor_elements:
                acceptors.append(atom)
        
        return acceptors
    
    def _extract_hydrophobic_atoms(self) -> List[Tuple]:
        """
        提取疏水原子
        
        标准：疏水氨基酸的碳原子
        """
        hydrophobic_residues = {
            'ALA', 'VAL', 'LEU', 'ILE', 'MET', 
            'PHE', 'TRP', 'PRO', 'TYR', 'CYS'
        }
        
        hydrophobic = []
        for atom in self.atoms:
            resn = atom[1]
            elem = get_element_from_atom_name(atom[3])
            
            # 疏水残基的碳原子
            if resn in hydrophobic_residues and elem == 'C':
                # 排除主链碳
                atom_name = atom[3].strip().upper()
                if atom_name not in ['C', 'CA', 'N', 'O']:
                    hydrophobic.append(atom)
        
        return hydrophobic
    
    def _extract_aromatic_rings(self) -> List[Dict]:
        """
        提取芳香环
        
        返回格式：
        [{
            'residue': (chain, resn, resi),
            'coords': [(x,y,z), ...],
            'center': (x,y,z),
            'normal': [nx, ny, nz]
        }, ...]
        """
        # 芳香环原子定义
        ring_definitions = {
            'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
            'TRP': ['CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
            'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
            # 核酸碱基
            'A': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9'],
            'G': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9'],
            'C': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
            'T': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
            'U': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
            'DA': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9'],
            'DG': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N7', 'C8', 'N9'],
            'DC': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
            'DT': ['N1', 'C2', 'N3', 'C4', 'C5', 'C6'],
        }
        
        rings = []
        
        for res_key, res_atoms in self.residues.items():
            resn = res_key[1]
            
            if resn not in ring_definitions:
                continue
            
            # 提取环原子坐标
            ring_atom_names = ring_definitions[resn]
            ring_coords = []
            
            for atom in res_atoms:
                atom_name = atom[3].strip().upper()
                if atom_name in ring_atom_names:
                    ring_coords.append(atom[4])
            
            # 至少需要3个原子才能定义平面
            if len(ring_coords) >= 3:
                ring_center = centroid(ring_coords)
                ring_normal = normal_vector(ring_coords[0], ring_coords[1], ring_coords[2])
                
                rings.append({
                    'residue': res_key,
                    'coords': ring_coords,
                    'center': ring_center,
                    'normal': ring_normal
                })
        
        return rings
    
    def _extract_charged_groups(self) -> List[Dict]:
        """
        提取带电基团
        
        基于生理 pH (~7.4) 下的 pKa 值判断残基离子化状态：
        - ARG (pKa ~12.5), LYS (pKa ~10.5): 始终带正电
        - HIP/HSP（双质子化 HIS）: 带正电
        - HIS (pKa ~6.0): 中性，不计入正电基团
        - ASP (pKa ~3.7), GLU (pKa ~4.1): 始终带负电
        
        返回格式：
        [{'residue': (chain, resn, resi), 'charge': '+' or '-',
          'center': (x,y,z), 'atoms': [(x,y,z), ...]}, ...]
        """
        # 与 interaction_analyzer.py 中的常量定义保持一致
        positive_residues = {'ARG', 'LYS'}
        # 双质子化组氨酸（HIP=AMBER, HSP=CHARMM）在 PDB 中明确带正电
        protonated_his_residues = {'HIP', 'HSP'}
        # 合并所有正电残基
        all_positive = positive_residues | protonated_his_residues
        negative_residues = {'ASP', 'GLU'}
        
        # 带电原子定义（PDB 命名规范）
        positive_atoms_arg_lys = {'NZ', 'NH1', 'NH2', 'NE'}
        positive_atoms_his = {'ND1', 'NE2'}  # HIP/HSP 咪唑环 N 原子
        negative_atoms = {'OD1', 'OD2', 'OE1', 'OE2'}
        
        charged = []
        
        for res_key, res_atoms in self.residues.items():
            resn = res_key[1]
            
            if resn in all_positive:
                # 正电荷中心
                # 根据残基类型选择正确的带电原子集合
                if resn in protonated_his_residues:
                    target_atoms = positive_atoms_his
                else:
                    target_atoms = positive_atoms_arg_lys
                
                charge_coords = []
                for atom in res_atoms:
                    if atom[3].strip().upper() in target_atoms:
                        charge_coords.append(atom[4])
                
                if charge_coords:
                    charged.append({
                        'residue': res_key,
                        'charge': '+',
                        'center': centroid(charge_coords),
                        'atoms': charge_coords
                    })
            
            elif resn in negative_residues:
                # 负电荷中心
                charge_coords = []
                for atom in res_atoms:
                    if atom[3].strip().upper() in negative_atoms:
                        charge_coords.append(atom[4])
                
                if charge_coords:
                    charged.append({
                        'residue': res_key,
                        'charge': '-',
                        'center': centroid(charge_coords),
                        'atoms': charge_coords
                    })
        
        return charged
    
    def _extract_metal_ions(self) -> List[Tuple]:
        """
        提取金属离子
        
        标准金属：Zn, Mg, Ca, Fe, Cu, Mn, Co, Ni
        """
        metal_elements = {'ZN', 'MG', 'CA', 'FE', 'CU', 'MN', 'CO', 'NI'}
        
        metals = []
        for atom in self.atoms:
            elem = get_element_from_atom_name(atom[3])
            if elem in metal_elements:
                metals.append(atom)
        
        return metals
    
    def _extract_halogen_atoms(self) -> List[Tuple]:
        """
        提取卤素原子
        
        卤素：F, Cl, Br, I
        """
        halogen_elements = {'F', 'CL', 'BR', 'I'}
        
        halogens = []
        for atom in self.atoms:
            elem = get_element_from_atom_name(atom[3])
            if elem in halogen_elements:
                halogens.append(atom)
        
        return halogens
    
    # ========================================================================
    # 空间索引（可选，需要 scipy）
    # ========================================================================
    
    def build_spatial_index(self):
        """
        构建空间索引（KD-Tree）
        
        需要 scipy.spatial.cKDTree
        """
        try:
            from scipy.spatial import cKDTree
            
            coords = [atom[4] for atom in self.atoms]
            self.kdtree = cKDTree(coords)
            self.use_kdtree = True
            
            print("[FeatureExtractor] Spatial index built (KD-Tree)")
        except ImportError:
            self.kdtree = None
            self.use_kdtree = False
            print("[FeatureExtractor] scipy not available, spatial index disabled")
    
    def find_neighbors(self, query_coord: Tuple[float, float, float], 
                       max_distance: float) -> List[Tuple]:
        """
        查找邻近原子（使用空间索引加速）
        
        Args:
            query_coord: 查询坐标
            max_distance: 最大距离
        
        Returns:
            List[Tuple]: 邻近原子列表
        """
        if hasattr(self, 'use_kdtree') and self.use_kdtree:
            # 使用 KD-Tree
            indices = self.kdtree.query_ball_point(query_coord, max_distance)
            return [self.atoms[i] for i in indices]
        else:
            # 暴力搜索
            neighbors = []
            for atom in self.atoms:
                if distance(query_coord, atom[4]) <= max_distance:
                    neighbors.append(atom)
            return neighbors
    
    # ========================================================================
    # 实用方法
    # ========================================================================
    
    def get_residue_atoms(self, chain: str, resn: str, resi: str) -> List[Tuple]:
        """
        获取指定残基的所有原子
        
        Args:
            chain: 链ID
            resn: 残基名
            resi: 残基编号
        
        Returns:
            List[Tuple]: 原子列表
        """
        res_key = (chain, resn, resi)
        return self.residues.get(res_key, [])
    
    def get_atom_by_name(self, chain: str, resn: str, resi: str, 
                         atom_name: str) -> Optional[Tuple]:
        """
        按名称查找特定原子
        
        Args:
            chain: 链ID
            resn: 残基名
            resi: 残基编号
            atom_name: 原子名
        
        Returns:
            Tuple or None: 原子信息
        """
        res_atoms = self.get_residue_atoms(chain, resn, resi)
        atom_name_upper = atom_name.strip().upper()
        
        for atom in res_atoms:
            if atom[3].strip().upper() == atom_name_upper:
                return atom
        
        return None
    
    def summary(self) -> Dict[str, int]:
        """
        返回特征统计摘要
        
        Returns:
            dict: 各类特征数量
        """
        return {
            'total_atoms': len(self.atoms),
            'total_residues': len(self.residues),
            'hbond_donors': len(self.hbond_donors),
            'hbond_acceptors': len(self.hbond_acceptors),
            'hydrophobic_atoms': len(self.hydrophobic_atoms),
            'aromatic_rings': len(self.aromatic_rings),
            'charged_groups': len(self.charged_groups),
            'metal_ions': len(self.metal_ions),
            'halogen_atoms': len(self.halogen_atoms),
        }
    
    def __str__(self) -> str:
        summary = self.summary()
        return (f"MolecularFeatureExtractor("
                f"atoms={summary['total_atoms']}, "
                f"residues={summary['total_residues']}, "
                f"rings={summary['aromatic_rings']}, "
                f"charged={summary['charged_groups']})")
    
    def __repr__(self) -> str:
        return self.__str__()
