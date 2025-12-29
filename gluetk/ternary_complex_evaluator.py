#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
三元复合物评估模块 (Ternary Complex Evaluator)

整合三类计算性质：
1. Interface Module - BSA/SASA表面积计算、接触数统计
2. Ligand Module - 小分子理化参数（MW/LogP/TPSA等）
3. Ternary Geometry Module - 三元复合物重心、几何特征

依赖：BioPython, RDKit, numpy
"""

import numpy as np
import logging
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from pathlib import Path

logger = logging.getLogger(__name__)

# ==================== 数据结构 ====================

@dataclass
class AtomInfo:
    """原子信息"""
    atom_id: int
    atom_name: str
    element: str
    coords: np.ndarray
    residue_name: str
    residue_id: int
    chain_id: str
    is_hetatm: bool = False

@dataclass
class ChainInfo:
    """链信息"""
    chain_id: str
    chain_type: str  # "protein" or "ligand"
    atoms: List[AtomInfo] = field(default_factory=list)
    
    @property
    def center_of_mass(self) -> np.ndarray:
        """计算质心"""
        if not self.atoms:
            return np.zeros(3)
        coords = np.array([a.coords for a in self.atoms])
        return coords.mean(axis=0)

@dataclass
class InterfaceFeatures:
    """界面特征"""
    bsa_total: float = 0.0           # 总掩埋表面积
    bsa_mg_e3: float = 0.0           # MG-E3界面BSA
    bsa_mg_poi: float = 0.0          # MG-POI界面BSA
    bsa_e3_poi: float = 0.0          # E3-POI界面BSA
    contact_count_45: int = 0        # 4.5Å接触数
    contact_count_50: int = 0        # 5.0Å接触数
    min_inter_chain_dist: float = 0.0

@dataclass
class LigandFeatures:
    """配体特征"""
    molecular_weight: float = 0.0
    logp: float = 0.0
    tpsa: float = 0.0
    rotatable_bonds: int = 0
    hbd_count: int = 0
    hba_count: int = 0
    fsp3: float = 0.0
    num_rings: int = 0
    net_charge_ph74: float = 0.0

@dataclass
class GeometryFeatures:
    """几何特征"""
    cog_e3: np.ndarray = field(default_factory=lambda: np.zeros(3))
    cog_poi: np.ndarray = field(default_factory=lambda: np.zeros(3))
    cog_mg: np.ndarray = field(default_factory=lambda: np.zeros(3))
    cog_shift: float = 0.0           # 重心偏移（最重要特征）
    angle_deg: float = 0.0           # 关键向量夹角
    dist_e3_poi: float = 0.0
    dist_e3_mg: float = 0.0
    dist_poi_mg: float = 0.0

@dataclass
class TernaryComplexFeatures:
    """三元复合物综合特征"""
    interface: InterfaceFeatures = field(default_factory=InterfaceFeatures)
    ligand: LigandFeatures = field(default_factory=LigandFeatures)
    geometry: GeometryFeatures = field(default_factory=GeometryFeatures)
    cooperativity_energy: Optional[float] = None
    hook_risk_score: Optional[float] = None
    duality_index: Optional[float] = None

# ==================== SASA计算器 ====================

class SASACalculator:
    """SASA计算器 - Shrake-Rupley算法"""
    
    VDW_RADII = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52,
        'F': 1.47, 'P': 1.80, 'S': 1.80, 'CL': 1.75
    }
    DEFAULT_RADIUS = 1.70
    PROBE_RADIUS = 1.40

    def __init__(self, n_points: int = 100):
        self.n_points = n_points
        self.sphere_points = self._generate_sphere_points(n_points)

    def _generate_sphere_points(self, n: int) -> np.ndarray:
        """Golden Spiral均匀分布点"""
        indices = np.arange(0, n, dtype=float) + 0.5
        phi = np.arccos(1 - 2*indices/n)
        theta = np.pi * (1 + 5**0.5) * indices
        return np.column_stack((
            np.cos(theta) * np.sin(phi),
            np.sin(theta) * np.sin(phi),
            np.cos(phi)
        ))

    def calculate_sasa(self, atoms: List[AtomInfo]) -> float:
        """计算原子列表的总SASA"""
        if not atoms:
            return 0.0
        
        coords = np.array([a.coords for a in atoms])
        radii = np.array([self.VDW_RADII.get(a.element, self.DEFAULT_RADIUS) for a in atoms])
        effective_radii = radii + self.PROBE_RADIUS
        
        total_area = 0.0
        for i in range(len(atoms)):
            current_radius = effective_radii[i]
            test_points = coords[i] + self.sphere_points * current_radius
            
            # 检查邻居遮挡
            accessible = 0
            for point in test_points:
                dists = np.linalg.norm(coords - point, axis=1)
                dists[i] = np.inf  # 排除自身
                if np.all(dists > effective_radii):
                    accessible += 1
            
            total_area += 4 * np.pi * current_radius**2 * (accessible / self.n_points)
        
        return total_area

    def calculate_bsa(self, atoms1: List[AtomInfo], atoms2: List[AtomInfo]) -> float:
        """计算埋藏表面积"""
        sasa1 = self.calculate_sasa(atoms1)
        sasa2 = self.calculate_sasa(atoms2)
        sasa_complex = self.calculate_sasa(atoms1 + atoms2)
        return (sasa1 + sasa2 - sasa_complex) / 2.0

# ==================== 结构解析器 ====================

class StructureParser:
    """PDB结构解析器"""
    
    STANDARD_RESIDUES = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    SOLVENT = {'HOH', 'WAT', 'NA', 'CL', 'MG', 'CA', 'ZN'}

    def __init__(self):
        self.structure = None
        self._biopython_available = False
        try:
            from Bio.PDB import PDBParser
            self.pdb_parser = PDBParser(QUIET=True)
            self._biopython_available = True
        except ImportError:
            logger.warning("BioPython不可用，使用简化解析")

    def parse_structure(self, pdb_path: str) -> Dict[str, ChainInfo]:
        """解析PDB文件"""
        if self._biopython_available:
            return self._parse_with_biopython(pdb_path)
        return self._parse_simple(pdb_path)

    def _parse_with_biopython(self, pdb_path: str) -> Dict[str, ChainInfo]:
        """使用BioPython解析"""
        self.structure = self.pdb_parser.get_structure("struct", pdb_path)
        chains = {}
        
        for model in self.structure:
            for chain in model:
                atoms = []
                for residue in chain:
                    if residue.get_resname() in self.SOLVENT:
                        continue
                    for atom in residue:
                        atoms.append(AtomInfo(
                            atom_id=atom.get_serial_number(),
                            atom_name=atom.get_name(),
                            element=atom.element,
                            coords=atom.get_coord(),
                            residue_name=residue.get_resname(),
                            residue_id=residue.id[1],
                            chain_id=chain.id,
                            is_hetatm=(residue.id[0] != ' ')
                        ))
                
                # 判断链类型
                standard_count = sum(1 for a in atoms if a.residue_name in self.STANDARD_RESIDUES)
                chain_type = "protein" if len(atoms) > 0 and standard_count / len(atoms) > 0.5 else "ligand"
                
                chains[chain.id] = ChainInfo(chain_id=chain.id, chain_type=chain_type, atoms=atoms)
        
        return chains

    def _parse_simple(self, pdb_path: str) -> Dict[str, ChainInfo]:
        """简化PDB解析"""
        chains = {}
        with open(pdb_path, 'r') as f:
            for line in f:
                if not line.startswith(('ATOM', 'HETATM')):
                    continue
                
                chain_id = line[21]
                if chain_id not in chains:
                    chains[chain_id] = ChainInfo(chain_id=chain_id, chain_type="unknown", atoms=[])
                
                atom = AtomInfo(
                    atom_id=int(line[6:11]),
                    atom_name=line[12:16].strip(),
                    element=line[76:78].strip() or line[12:14].strip()[0],
                    coords=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                    residue_name=line[17:20].strip(),
                    residue_id=int(line[22:26]),
                    chain_id=chain_id,
                    is_hetatm=line.startswith('HETATM')
                )
                chains[chain_id].atoms.append(atom)
        
        # 判断链类型
        for chain in chains.values():
            standard = sum(1 for a in chain.atoms if a.residue_name in self.STANDARD_RESIDUES)
            chain.chain_type = "protein" if len(chain.atoms) > 0 and standard / len(chain.atoms) > 0.5 else "ligand"
        
        return chains

# ==================== 小分子计算器 ====================

class LigandCalculator:
    """小分子性质计算器"""
    
    def __init__(self):
        self._rdkit_available = False
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, AllChem
            self._rdkit_available = True
            self.Chem = Chem
            self.Descriptors = Descriptors
            self.AllChem = AllChem
        except ImportError:
            logger.warning("RDKit不可用，小分子计算功能受限")

    def calculate(self, smiles: Optional[str] = None, mol: Any = None) -> LigandFeatures:
        """计算小分子性质"""
        features = LigandFeatures()
        
        if not self._rdkit_available:
            return features
        
        if mol is None and smiles:
            mol = self.Chem.MolFromSmiles(smiles)
        
        if mol is None:
            return features
        
        mol = self.Chem.AddHs(mol)
        
        features.molecular_weight = self.Descriptors.MolWt(mol)
        features.logp = self.Descriptors.MolLogP(mol)
        features.tpsa = self.Descriptors.TPSA(mol)
        features.rotatable_bonds = self.Descriptors.NumRotatableBonds(mol)
        features.hbd_count = self.Descriptors.NumHDonors(mol)
        features.hba_count = self.Descriptors.NumHAcceptors(mol)
        features.num_rings = self.Descriptors.RingCount(mol)
        
        try:
            from rdkit.Chem import rdMolDescriptors
            features.fsp3 = rdMolDescriptors.CalcFractionCsp3(mol)
        except:
            pass
        
        return features

# ==================== 三元复合物评估器 ====================

class TernaryComplexEvaluator:
    """三元复合物评估器 - 整合Interface、Ligand、Geometry计算"""
    
    def __init__(self):
        self.parser = StructureParser()
        self.sasa_calc = SASACalculator(n_points=100)
        self.ligand_calc = LigandCalculator()

    def evaluate(
        self,
        pdb_path: str,
        e3_chain: str,
        poi_chain: str,
        ligand_chain: str,
        ligand_smiles: Optional[str] = None
    ) -> TernaryComplexFeatures:
        """
        评估三元复合物
        
        Args:
            pdb_path: PDB文件路径
            e3_chain: E3链ID
            poi_chain: POI链ID
            ligand_chain: 配体链ID
            ligand_smiles: 配体SMILES（可选，用于计算小分子性质）
        
        Returns:
            TernaryComplexFeatures: 综合特征
        """
        features = TernaryComplexFeatures()
        
        # 1. 解析结构
        chains = self.parser.parse_structure(pdb_path)
        
        e3 = chains.get(e3_chain)
        poi = chains.get(poi_chain)
        mg = chains.get(ligand_chain)
        
        if not all([e3, poi, mg]):
            missing = [c for c, v in [('E3', e3), ('POI', poi), ('MG', mg)] if not v]
            logger.error(f"缺少链: {missing}")
            return features
        
        # 2. 计算界面特征
        features.interface = self._calculate_interface(e3, poi, mg)
        
        # 3. 计算配体特征
        if ligand_smiles:
            features.ligand = self.ligand_calc.calculate(smiles=ligand_smiles)
        
        # 4. 计算几何特征
        features.geometry = self._calculate_geometry(e3, poi, mg)
        
        # 5. 计算协同性指标
        features.cooperativity_energy = self._estimate_cooperativity(features)
        features.hook_risk_score = self._estimate_hook_risk(features)
        features.duality_index = self._calculate_duality(features)
        
        return features

    def _calculate_interface(self, e3: ChainInfo, poi: ChainInfo, mg: ChainInfo) -> InterfaceFeatures:
        """计算界面特征"""
        interface = InterfaceFeatures()
        
        # BSA计算
        interface.bsa_mg_e3 = self.sasa_calc.calculate_bsa(mg.atoms, e3.atoms)
        interface.bsa_mg_poi = self.sasa_calc.calculate_bsa(mg.atoms, poi.atoms)
        interface.bsa_e3_poi = self.sasa_calc.calculate_bsa(e3.atoms, poi.atoms)
        interface.bsa_total = interface.bsa_mg_e3 + interface.bsa_mg_poi + interface.bsa_e3_poi
        
        # 接触数统计
        interface.contact_count_45 = self._count_contacts(mg.atoms, e3.atoms + poi.atoms, 4.5)
        interface.contact_count_50 = self._count_contacts(mg.atoms, e3.atoms + poi.atoms, 5.0)
        
        # 最小链间距离
        interface.min_inter_chain_dist = self._min_distance(e3.atoms, poi.atoms)
        
        return interface

    def _calculate_geometry(self, e3: ChainInfo, poi: ChainInfo, mg: ChainInfo) -> GeometryFeatures:
        """计算几何特征"""
        geom = GeometryFeatures()
        
        # 质心坐标
        geom.cog_e3 = e3.center_of_mass
        geom.cog_poi = poi.center_of_mass
        geom.cog_mg = mg.center_of_mass
        
        # 质心间距离
        geom.dist_e3_poi = np.linalg.norm(geom.cog_e3 - geom.cog_poi)
        geom.dist_e3_mg = np.linalg.norm(geom.cog_e3 - geom.cog_mg)
        geom.dist_poi_mg = np.linalg.norm(geom.cog_poi - geom.cog_mg)
        
        # 重心偏移（MG到E3-POI连线的距离）
        e3_poi_vec = geom.cog_poi - geom.cog_e3
        e3_mg_vec = geom.cog_mg - geom.cog_e3
        
        if np.linalg.norm(e3_poi_vec) > 0:
            proj = np.dot(e3_mg_vec, e3_poi_vec) / np.dot(e3_poi_vec, e3_poi_vec) * e3_poi_vec
            geom.cog_shift = np.linalg.norm(e3_mg_vec - proj)
        
        # 向量夹角
        if np.linalg.norm(e3_mg_vec) > 0 and np.linalg.norm(e3_poi_vec) > 0:
            cos_angle = np.dot(e3_mg_vec, e3_poi_vec) / (np.linalg.norm(e3_mg_vec) * np.linalg.norm(e3_poi_vec))
            geom.angle_deg = np.degrees(np.arccos(np.clip(cos_angle, -1, 1)))
        
        return geom

    def _count_contacts(self, atoms1: List[AtomInfo], atoms2: List[AtomInfo], cutoff: float) -> int:
        """统计接触数"""
        count = 0
        for a1 in atoms1:
            for a2 in atoms2:
                if np.linalg.norm(a1.coords - a2.coords) < cutoff:
                    count += 1
        return count

    def _min_distance(self, atoms1: List[AtomInfo], atoms2: List[AtomInfo]) -> float:
        """计算最小距离"""
        if not atoms1 or not atoms2:
            return float('inf')
        
        min_dist = float('inf')
        for a1 in atoms1:
            for a2 in atoms2:
                d = np.linalg.norm(a1.coords - a2.coords)
                if d < min_dist:
                    min_dist = d
        return min_dist

    def _estimate_cooperativity(self, features: TernaryComplexFeatures) -> float:
        """估算协同性能量"""
        # 简化估算：基于BSA和几何特征
        # 实际应用中应使用Rosetta或MM/GBSA
        bsa_score = -0.01 * features.interface.bsa_total
        geom_score = -0.1 * (1.0 / (1.0 + features.geometry.cog_shift))
        return bsa_score + geom_score

    def _estimate_hook_risk(self, features: TernaryComplexFeatures) -> float:
        """估算Hook效应风险"""
        # Hook_risk = max(BSA_E3, BSA_POI) / BSA_total
        bsa_e3 = features.interface.bsa_mg_e3
        bsa_poi = features.interface.bsa_mg_poi
        bsa_total = features.interface.bsa_total
        
        if bsa_total > 0:
            return max(bsa_e3, bsa_poi) / bsa_total
        return 0.0

    def _calculate_duality(self, features: TernaryComplexFeatures) -> float:
        """计算双面性指数"""
        bsa_e3 = features.interface.bsa_mg_e3
        bsa_poi = features.interface.bsa_mg_poi
        
        if max(bsa_e3, bsa_poi) > 0:
            return min(bsa_e3, bsa_poi) / max(bsa_e3, bsa_poi)
        return 0.0

    def to_dict(self, features: TernaryComplexFeatures) -> Dict[str, Any]:
        """转换为字典格式"""
        return {
            # Interface
            'bsa_total': features.interface.bsa_total,
            'bsa_mg_e3': features.interface.bsa_mg_e3,
            'bsa_mg_poi': features.interface.bsa_mg_poi,
            'bsa_e3_poi': features.interface.bsa_e3_poi,
            'contact_count_45': features.interface.contact_count_45,
            'contact_count_50': features.interface.contact_count_50,
            'min_inter_chain_dist': features.interface.min_inter_chain_dist,
            
            # Ligand
            'ligand_mw': features.ligand.molecular_weight,
            'ligand_logp': features.ligand.logp,
            'ligand_tpsa': features.ligand.tpsa,
            'ligand_rotatable_bonds': features.ligand.rotatable_bonds,
            'ligand_hbd': features.ligand.hbd_count,
            'ligand_hba': features.ligand.hba_count,
            'ligand_fsp3': features.ligand.fsp3,
            'ligand_rings': features.ligand.num_rings,
            
            # Geometry
            'cog_e3_x': features.geometry.cog_e3[0],
            'cog_e3_y': features.geometry.cog_e3[1],
            'cog_e3_z': features.geometry.cog_e3[2],
            'cog_poi_x': features.geometry.cog_poi[0],
            'cog_poi_y': features.geometry.cog_poi[1],
            'cog_poi_z': features.geometry.cog_poi[2],
            'cog_mg_x': features.geometry.cog_mg[0],
            'cog_mg_y': features.geometry.cog_mg[1],
            'cog_mg_z': features.geometry.cog_mg[2],
            'geom_cog_shift': features.geometry.cog_shift,
            'geom_angle_deg': features.geometry.angle_deg,
            'dist_e3_poi': features.geometry.dist_e3_poi,
            'dist_e3_mg': features.geometry.dist_e3_mg,
            'dist_poi_mg': features.geometry.dist_poi_mg,
            
            # Derived
            'cooperativity_energy': features.cooperativity_energy,
            'hook_risk_score': features.hook_risk_score,
            'duality_index': features.duality_index,
        }


# ==================== 便捷函数 ====================

def evaluate_ternary_complex(
    pdb_path: str,
    e3_chain: str,
    poi_chain: str,
    ligand_chain: str,
    ligand_smiles: Optional[str] = None
) -> Dict[str, Any]:
    """
    评估三元复合物的便捷函数
    
    Args:
        pdb_path: PDB文件路径
        e3_chain: E3链ID
        poi_chain: POI链ID
        ligand_chain: 配体链ID
        ligand_smiles: 配体SMILES（可选）
    
    Returns:
        Dict: 特征字典
    
    Example:
        >>> features = evaluate_ternary_complex(
        ...     "complex.pdb",
        ...     e3_chain="A",
        ...     poi_chain="B", 
        ...     ligand_chain="C",
        ...     ligand_smiles="CC(=O)Nc1ccc(O)cc1"
        ... )
        >>> print(f"BSA Total: {features['bsa_total']:.1f} Å²")
        >>> print(f"COG Shift: {features['geom_cog_shift']:.2f} Å")
    """
    evaluator = TernaryComplexEvaluator()
    features = evaluator.evaluate(pdb_path, e3_chain, poi_chain, ligand_chain, ligand_smiles)
    return evaluator.to_dict(features)


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 5:
        print("Usage: python ternary_complex_evaluator.py <pdb_path> <e3_chain> <poi_chain> <ligand_chain> [smiles]")
        sys.exit(1)
    
    pdb_path = sys.argv[1]
    e3_chain = sys.argv[2]
    poi_chain = sys.argv[3]
    ligand_chain = sys.argv[4]
    smiles = sys.argv[5] if len(sys.argv) > 5 else None
    
    features = evaluate_ternary_complex(pdb_path, e3_chain, poi_chain, ligand_chain, smiles)
    
    print("\n=== 三元复合物评估结果 ===\n")
    print("【界面特征】")
    print(f"  总BSA: {features['bsa_total']:.1f} Å²")
    print(f"  MG-E3 BSA: {features['bsa_mg_e3']:.1f} Å²")
    print(f"  MG-POI BSA: {features['bsa_mg_poi']:.1f} Å²")
    print(f"  接触数(4.5Å): {features['contact_count_45']}")
    
    print("\n【几何特征】")
    print(f"  重心偏移: {features['geom_cog_shift']:.2f} Å")
    print(f"  向量夹角: {features['geom_angle_deg']:.1f}°")
    print(f"  E3-POI距离: {features['dist_e3_poi']:.1f} Å")
    
    print("\n【协同性指标】")
    print(f"  协同能: {features['cooperativity_energy']:.2f} kcal/mol")
    print(f"  Hook风险: {features['hook_risk_score']:.2f}")
    print(f"  双面性指数: {features['duality_index']:.2f}")
    
    if smiles:
        print("\n【配体特征】")
        print(f"  分子量: {features['ligand_mw']:.1f} Da")
        print(f"  LogP: {features['ligand_logp']:.2f}")
        print(f"  TPSA: {features['ligand_tpsa']:.1f} Å²")