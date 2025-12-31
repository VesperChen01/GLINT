#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
三元复合物评估模块 (Ternary Complex Evaluator)

整合三类计算性质：
1. Interface Module - BSA/SASA表面积计算、接触数统计
2. Ligand Module - 小分子理化参数（MW/LogP/TPSA等）
3. Ternary Geometry Module - 三元复合物重心、几何特征

依赖：BioPython, RDKit, numpy, PyMOL
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

# ==================== BSA计算器（使用PyMOL） ====================

class BSACalculator:
    """BSA计算器 - 使用PyMOL的cmd.get_area()"""
    
    def __init__(self, obj_name: str = None, ligand_resn: str = None):
        self.obj_name = obj_name
        self.ligand_resn = ligand_resn
        self._pymol_available = False
        try:
            from pymol import cmd
            self.cmd = cmd
            self._pymol_available = True
        except ImportError:
            logger.warning("PyMOL不可用，BSA计算功能受限")
    
    def calculate_bsa_chain_chain(self, chain1: str, chain2: str) -> float:
        """
        计算两条蛋白链之间的埋藏表面积 (Buried Surface Area)
        
        BSA = (SA_chain1 + SA_chain2 - SA_complex) / 2
        
        参数:
            chain1: 第一条链ID
            chain2: 第二条链ID
        
        返回:
            float: 埋藏表面积 (Å²)
        """
        if not self._pymol_available or not self.obj_name:
            return 0.0
        
        import uuid
        suffix = str(uuid.uuid4())[:8]
        
        obj_c1 = f"temp_c1_{suffix}"
        obj_c2 = f"temp_c2_{suffix}"
        obj_complex = f"temp_complex_{suffix}"
        
        try:
            solvent_sel = "resn HOH+WAT+NA+CL+MG+CA+ZN"
            # 排除配体残基，只计算蛋白质部分
            ligand_exclude = f" and not resn {self.ligand_resn}" if self.ligand_resn else ""
            sel_base = f"{self.obj_name} and not ({solvent_sel}){ligand_exclude}"
            
            self.cmd.create(obj_c1, f"{sel_base} and chain {chain1}")
            area_chain1 = self.cmd.get_area(obj_c1)
            
            self.cmd.create(obj_c2, f"{sel_base} and chain {chain2}")
            area_chain2 = self.cmd.get_area(obj_c2)
            
            self.cmd.create(obj_complex, f"{sel_base} and (chain {chain1} or chain {chain2})")
            area_complex = self.cmd.get_area(obj_complex)
            
            bsa = (area_chain1 + area_chain2 - area_complex) / 2.0
            
            self.cmd.delete(obj_c1)
            self.cmd.delete(obj_c2)
            self.cmd.delete(obj_complex)
            
            return max(0.0, bsa)
        
        except Exception as e:
            logger.error(f"BSA计算失败: {e}")
            try:
                self.cmd.delete(obj_c1)
                self.cmd.delete(obj_c2)
                self.cmd.delete(obj_complex)
            except:
                pass
            return 0.0
    
    def calculate_bsa_ligand_chain(self, chain: str) -> float:
        """
        计算配体与蛋白链之间的埋藏表面积 (Buried Surface Area)
        
        使用标准 BSA 公式：BSA = (SA_ligand + SA_chain - SA_complex) / 2
        
        由于 PyMOL 的 get_area() 对 HETATM 原子可能返回 0，
        我们尝试使用 FreeSASA 库作为备选方案。
        
        参数:
            chain: 蛋白链ID
        
        返回:
            float: 埋藏表面积 (Å²)
        """
        if not self._pymol_available or not self.obj_name or not self.ligand_resn:
            print(f"[BSA Debug] Skipping: pymol={self._pymol_available}, obj={self.obj_name}, resn={self.ligand_resn}")
            return 0.0
        
        import uuid
        suffix = str(uuid.uuid4())[:8]
        
        obj_lig = f"temp_lig_{suffix}"
        obj_chain = f"temp_chain_{suffix}"
        obj_complex = f"temp_complex_{suffix}"
        
        try:
            solvent_sel = "resn HOH+WAT+NA+CL+MG+CA+ZN"
            sel_base = f"{self.obj_name} and not ({solvent_sel})"
            
            # 选择配体（按 residue name）
            lig_sel = f"{sel_base} and resn {self.ligand_resn}"
            # 选择蛋白链（排除配体）
            chain_sel = f"{sel_base} and chain {chain} and not resn {self.ligand_resn}"
            
            print(f"[BSA Debug] === Calculating BSA for chain {chain} ===")
            
            # 检查选择是否有原子
            lig_count = self.cmd.count_atoms(lig_sel)
            chain_count = self.cmd.count_atoms(chain_sel)
            print(f"[BSA Debug] Ligand atoms: {lig_count}, Chain {chain} atoms: {chain_count}")
            
            if lig_count == 0 or chain_count == 0:
                return 0.0
            
            # 保存当前设置
            old_dot_solvent = self.cmd.get("dot_solvent")
            old_dot_density = self.cmd.get("dot_density")
            
            # 设置表面积计算参数
            self.cmd.set("dot_solvent", 1)
            self.cmd.set("dot_density", 3)
            
            # 尝试使用 FreeSASA 计算配体面积
            area_lig = 0.0
            try:
                import freesasa
                # 导出配体到临时 PDB 文件
                import tempfile
                import os
                fd, temp_pdb = tempfile.mkstemp(suffix=".pdb")
                os.close(fd)
                self.cmd.save(temp_pdb, lig_sel)
                
                # 使用 FreeSASA 计算
                structure = freesasa.Structure(temp_pdb)
                result = freesasa.calc(structure)
                area_lig = result.totalArea()
                print(f"[BSA Debug] Ligand area (FreeSASA): {area_lig:.1f}")
                
                os.remove(temp_pdb)
            except ImportError:
                print(f"[BSA Debug] FreeSASA not available, using PyMOL")
                # 使用 PyMOL 计算
                self.cmd.create(obj_lig, lig_sel)
                area_lig = self.cmd.get_area(obj_lig, state=1)
                print(f"[BSA Debug] Ligand area (PyMOL): {area_lig:.1f}")
                self.cmd.delete(obj_lig)
            except Exception as e:
                print(f"[BSA Debug] FreeSASA error: {e}, using PyMOL")
                self.cmd.create(obj_lig, lig_sel)
                area_lig = self.cmd.get_area(obj_lig, state=1)
                print(f"[BSA Debug] Ligand area (PyMOL): {area_lig:.1f}")
                self.cmd.delete(obj_lig)
            
            # 如果配体面积仍为 0，使用基于接触原子的估算
            if area_lig == 0:
                print(f"[BSA Debug] Ligand area is 0, using contact-based estimation")
                # 估算配体面积：每个重原子约 15-20 Å²
                area_lig = lig_count * 17.0  # 平均每个原子 17 Å²
                print(f"[BSA Debug] Estimated ligand area: {lig_count} atoms × 17 Å² = {area_lig:.1f} Å²")
            
            # 计算蛋白链面积
            self.cmd.create(obj_chain, chain_sel)
            area_chain = self.cmd.get_area(obj_chain, state=1)
            print(f"[BSA Debug] Chain {chain} area: {area_chain:.1f}")
            self.cmd.delete(obj_chain)
            
            # 计算复合物面积
            complex_sel = f"({lig_sel}) or ({chain_sel})"
            self.cmd.create(obj_complex, complex_sel)
            area_complex = self.cmd.get_area(obj_complex, state=1)
            
            # 如果复合物面积等于蛋白链面积（说明配体没有被计入），需要调整
            if abs(area_complex - area_chain) < 1.0:
                print(f"[BSA Debug] Complex area equals chain area, adjusting...")
                # 复合物面积 = 蛋白链面积 + 配体面积 - 2×BSA
                # 我们需要估算 BSA，使用接触原子方法
                contact_dist = 4.5
                lig_contact_sel = f"({lig_sel}) within {contact_dist} of ({chain_sel})"
                chain_contact_sel = f"({chain_sel}) within {contact_dist} of ({lig_sel})"
                lig_contact_count = self.cmd.count_atoms(lig_contact_sel)
                chain_contact_count = self.cmd.count_atoms(chain_contact_sel)
                
                # 估算 BSA：每对接触原子贡献约 15 Å²
                contact_pairs = min(lig_contact_count, chain_contact_count)
                bsa = contact_pairs * 15.0
                print(f"[BSA Debug] Contact-based BSA: {contact_pairs} pairs × 15 Å² = {bsa:.1f} Å²")
                
                self.cmd.delete(obj_complex)
                self.cmd.set("dot_solvent", old_dot_solvent)
                self.cmd.set("dot_density", old_dot_density)
                return bsa
            
            print(f"[BSA Debug] Complex area: {area_complex:.1f}")
            self.cmd.delete(obj_complex)
            
            # 恢复设置
            self.cmd.set("dot_solvent", old_dot_solvent)
            self.cmd.set("dot_density", old_dot_density)
            
            # 计算 BSA
            bsa = (area_lig + area_chain - area_complex) / 2.0
            print(f"[BSA Debug] BSA = ({area_lig:.1f} + {area_chain:.1f} - {area_complex:.1f}) / 2 = {bsa:.1f} Å²")
            
            return max(0.0, bsa)
        
        except Exception as e:
            logger.error(f"配体-蛋白BSA计算失败: {e}")
            print(f"[BSA Debug] Error: {e}")
            import traceback
            traceback.print_exc()
            try:
                self.cmd.delete(obj_lig)
                self.cmd.delete(obj_chain)
                self.cmd.delete(obj_complex)
            except:
                pass
            return 0.0

    def calculate_bsa_ternary(self, e3_chain: str, poi_chain: str) -> Dict[str, float]:
        """
        使用三元复合物公式计算 BSA
        
        公式：BSA_MG = (SA_E3-MG + SA_POI-MG) - SA_E3-POI-MG
        
        其中：
        - SA_E3-MG = E3 和 MG 组成的二元复合物的表面积
        - SA_POI-MG = POI 和 MG 组成的二元复合物的表面积
        - SA_E3-POI-MG = 三元复合物的总表面积
        
        参数:
            e3_chain: E3 链 ID
            poi_chain: POI 链 ID
        
        返回:
            Dict: 包含 bsa_mg_e3, bsa_mg_poi, bsa_e3_poi, bsa_total 的字典
        """
        result = {
            'bsa_mg_e3': 0.0,
            'bsa_mg_poi': 0.0,
            'bsa_e3_poi': 0.0,
            'bsa_total': 0.0,
            # 中间计算值
            'sa_e3_mg': 0.0,      # E3-MG 二元复合物表面积
            'sa_poi_mg': 0.0,     # POI-MG 二元复合物表面积
            'sa_ternary': 0.0,    # 三元复合物表面积
            'sa_e3': 0.0,         # E3 单独表面积
            'sa_poi': 0.0,        # POI 单独表面积
            'sa_e3_poi': 0.0,     # E3-POI 二元复合物表面积
            'contact_mg_e3': 0,   # MG-E3 接触原子数
            'contact_mg_poi': 0,  # MG-POI 接触原子数
        }
        
        if not self._pymol_available or not self.obj_name or not self.ligand_resn:
            print(f"[BSA Debug] Skipping ternary BSA: pymol={self._pymol_available}, obj={self.obj_name}, resn={self.ligand_resn}")
            return result
        
        import uuid
        suffix = str(uuid.uuid4())[:8]
        
        try:
            solvent_sel = "resn HOH+WAT+NA+CL+MG+CA+ZN"
            sel_base = f"{self.obj_name} and not ({solvent_sel})"
            
            # 选择各组分
            lig_sel = f"{sel_base} and resn {self.ligand_resn}"
            e3_sel = f"{sel_base} and chain {e3_chain} and not resn {self.ligand_resn}"
            poi_sel = f"{sel_base} and chain {poi_chain} and not resn {self.ligand_resn}"
            
            print(f"[BSA Debug] === Calculating Ternary BSA ===")
            
            # 检查选择是否有原子
            lig_count = self.cmd.count_atoms(lig_sel)
            e3_count = self.cmd.count_atoms(e3_sel)
            poi_count = self.cmd.count_atoms(poi_sel)
            print(f"[BSA Debug] Atoms: Ligand={lig_count}, E3={e3_count}, POI={poi_count}")
            
            if lig_count == 0 or e3_count == 0 or poi_count == 0:
                return result
            
            # 保存当前设置
            old_dot_solvent = self.cmd.get("dot_solvent")
            old_dot_density = self.cmd.get("dot_density")
            
            # 设置表面积计算参数
            self.cmd.set("dot_solvent", 1)
            self.cmd.set("dot_density", 3)
            
            # 创建临时对象
            obj_e3_mg = f"temp_e3_mg_{suffix}"
            obj_poi_mg = f"temp_poi_mg_{suffix}"
            obj_ternary = f"temp_ternary_{suffix}"
            obj_e3_poi = f"temp_e3_poi_{suffix}"
            
            # 计算 E3-MG 二元复合物面积
            self.cmd.create(obj_e3_mg, f"({e3_sel}) or ({lig_sel})")
            area_e3_mg = self.cmd.get_area(obj_e3_mg, state=1)
            print(f"[BSA Debug] SA(E3-MG): {area_e3_mg:.1f}")
            
            # 计算 POI-MG 二元复合物面积
            self.cmd.create(obj_poi_mg, f"({poi_sel}) or ({lig_sel})")
            area_poi_mg = self.cmd.get_area(obj_poi_mg, state=1)
            print(f"[BSA Debug] SA(POI-MG): {area_poi_mg:.1f}")
            
            # 计算三元复合物面积
            self.cmd.create(obj_ternary, f"({e3_sel}) or ({poi_sel}) or ({lig_sel})")
            area_ternary = self.cmd.get_area(obj_ternary, state=1)
            print(f"[BSA Debug] SA(E3-POI-MG): {area_ternary:.1f}")
            
            # 计算 E3-POI 二元复合物面积（不含配体）
            self.cmd.create(obj_e3_poi, f"({e3_sel}) or ({poi_sel})")
            area_e3_poi = self.cmd.get_area(obj_e3_poi, state=1)
            print(f"[BSA Debug] SA(E3-POI): {area_e3_poi:.1f}")
            
            # 计算各组分单独的面积
            obj_e3 = f"temp_e3_{suffix}"
            obj_poi = f"temp_poi_{suffix}"
            
            self.cmd.create(obj_e3, e3_sel)
            area_e3 = self.cmd.get_area(obj_e3, state=1)
            print(f"[BSA Debug] SA(E3): {area_e3:.1f}")
            
            self.cmd.create(obj_poi, poi_sel)
            area_poi = self.cmd.get_area(obj_poi, state=1)
            print(f"[BSA Debug] SA(POI): {area_poi:.1f}")
            
            # 清理临时对象
            self.cmd.delete(obj_e3_mg)
            self.cmd.delete(obj_poi_mg)
            self.cmd.delete(obj_ternary)
            self.cmd.delete(obj_e3_poi)
            self.cmd.delete(obj_e3)
            self.cmd.delete(obj_poi)
            
            # 恢复设置
            self.cmd.set("dot_solvent", old_dot_solvent)
            self.cmd.set("dot_density", old_dot_density)
            
            # 使用三元复合物公式计算 BSA
            # BSA_MG_total = (SA_E3-MG + SA_POI-MG) - SA_E3-POI-MG
            bsa_mg_total = (area_e3_mg + area_poi_mg) - area_ternary
            print(f"[BSA Debug] BSA_MG_total = ({area_e3_mg:.1f} + {area_poi_mg:.1f}) - {area_ternary:.1f} = {bsa_mg_total:.1f}")
            
            # 计算 E3-POI BSA（标准公式）
            bsa_e3_poi = (area_e3 + area_poi - area_e3_poi) / 2.0
            print(f"[BSA Debug] BSA_E3-POI = ({area_e3:.1f} + {area_poi:.1f} - {area_e3_poi:.1f}) / 2 = {bsa_e3_poi:.1f}")
            
            # 分配 MG 的 BSA 到 E3 和 POI
            # 使用接触原子数来估算比例
            contact_dist = 4.5
            lig_e3_contact = self.cmd.count_atoms(f"({lig_sel}) within {contact_dist} of ({e3_sel})")
            lig_poi_contact = self.cmd.count_atoms(f"({lig_sel}) within {contact_dist} of ({poi_sel})")
            total_contact = lig_e3_contact + lig_poi_contact
            
            print(f"[BSA Debug] Contact atoms: MG-E3={lig_e3_contact}, MG-POI={lig_poi_contact}")
            
            if total_contact > 0:
                ratio_e3 = lig_e3_contact / total_contact
                ratio_poi = lig_poi_contact / total_contact
            else:
                ratio_e3 = 0.5
                ratio_poi = 0.5
            
            bsa_mg_e3 = bsa_mg_total * ratio_e3
            bsa_mg_poi = bsa_mg_total * ratio_poi
            
            print(f"[BSA Debug] BSA_MG-E3 = {bsa_mg_total:.1f} × {ratio_e3:.2f} = {bsa_mg_e3:.1f}")
            print(f"[BSA Debug] BSA_MG-POI = {bsa_mg_total:.1f} × {ratio_poi:.2f} = {bsa_mg_poi:.1f}")
            
            result['bsa_mg_e3'] = max(0.0, bsa_mg_e3)
            result['bsa_mg_poi'] = max(0.0, bsa_mg_poi)
            result['bsa_e3_poi'] = max(0.0, bsa_e3_poi)
            result['bsa_total'] = result['bsa_mg_e3'] + result['bsa_mg_poi'] + result['bsa_e3_poi']
            
            # 保存中间计算值
            result['sa_e3_mg'] = area_e3_mg
            result['sa_poi_mg'] = area_poi_mg
            result['sa_ternary'] = area_ternary
            result['sa_e3'] = area_e3
            result['sa_poi'] = area_poi
            result['sa_e3_poi'] = area_e3_poi
            result['contact_mg_e3'] = lig_e3_contact
            result['contact_mg_poi'] = lig_poi_contact
            
            print(f"[BSA Debug] Total BSA = {result['bsa_total']:.1f}")
            
            return result
        
        except Exception as e:
            logger.error(f"三元复合物BSA计算失败: {e}")
            print(f"[BSA Debug] Error: {e}")
            import traceback
            traceback.print_exc()
            return result

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
        self.ligand_calc = LigandCalculator()
        self._obj_name = None  # PyMOL对象名

    def evaluate(
        self,
        pdb_path: str,
        e3_chain: str,
        poi_chain: str,
        ligand_resn: str,
        ligand_smiles: Optional[str] = None,
        obj_name: Optional[str] = None
    ) -> TernaryComplexFeatures:
        """
        评估三元复合物
        
        Args:
            pdb_path: PDB文件路径
            e3_chain: E3链ID
            poi_chain: POI链ID
            ligand_resn: 配体残基名称（如 UNL, LIG, MOL）
            ligand_smiles: 配体SMILES（可选，用于计算小分子性质）
            obj_name: PyMOL对象名（可选，用于BSA计算）
        
        Returns:
            TernaryComplexFeatures: 综合特征
        """
        features = TernaryComplexFeatures()
        
        # 1. 解析结构
        chains = self.parser.parse_structure(pdb_path)
        
        e3 = chains.get(e3_chain)
        poi = chains.get(poi_chain)
        
        # 按 residue name 查找配体（而不是 chain ID）
        mg = self._find_ligand_by_resn(chains, ligand_resn)
        
        if not all([e3, poi, mg]):
            missing = []
            if not e3:
                missing.append(f'E3 (chain {e3_chain})')
            if not poi:
                missing.append(f'POI (chain {poi_chain})')
            if not mg:
                missing.append(f'Ligand (resn {ligand_resn})')
            logger.error(f"缺少组分: {missing}")
            print(f"缺少组分: {missing}")
            return features
        
        # 2. 计算界面特征（使用PyMOL）
        self._obj_name = obj_name
        self._ligand_resn = ligand_resn
        if not obj_name:
            # 尝试加载PDB到PyMOL
            try:
                from pymol import cmd
                import os
                obj_name = os.path.splitext(os.path.basename(pdb_path))[0]
                cmd.load(pdb_path, obj_name)
                self._obj_name = obj_name
            except Exception as e:
                logger.warning(f"无法加载PDB到PyMOL: {e}")
        
        # 获取配体所在的链ID（用于BSA计算）
        mg_chain = mg.chain_id if mg else ""
        features.interface, self._intermediate_values = self._calculate_interface(e3, poi, mg, e3_chain, poi_chain, mg_chain)
        
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

    def _find_ligand_by_resn(self, chains: Dict[str, ChainInfo], ligand_resn: str) -> Optional[ChainInfo]:
        """
        按残基名称查找配体
        
        Args:
            chains: 所有链的字典
            ligand_resn: 配体残基名称（如 UNL, LIG, MOL）
        
        Returns:
            ChainInfo: 包含配体原子的虚拟链信息，如果未找到则返回 None
        """
        ligand_atoms = []
        ligand_chain_id = None
        
        for chain_id, chain in chains.items():
            for atom in chain.atoms:
                if atom.residue_name.upper() == ligand_resn.upper():
                    ligand_atoms.append(atom)
                    if ligand_chain_id is None:
                        ligand_chain_id = chain_id
        
        if ligand_atoms:
            # 创建一个虚拟的 ChainInfo 来存储配体原子
            return ChainInfo(
                chain_id=ligand_chain_id or "LIG",
                chain_type="ligand",
                atoms=ligand_atoms
            )
        
        return None

    def _calculate_interface(self, e3: ChainInfo, poi: ChainInfo, mg: ChainInfo,
                            e3_chain: str, poi_chain: str, mg_chain: str) -> Tuple[InterfaceFeatures, Dict[str, float]]:
        """
        计算界面特征 - 使用三元复合物 BSA 公式
        
        BSA_MG = (SA_E3-MG + SA_POI-MG) - SA_E3-POI-MG
        
        返回:
            Tuple[InterfaceFeatures, Dict]: 界面特征和中间计算值
        """
        interface = InterfaceFeatures()
        
        # 使用PyMOL计算BSA，传入配体残基名称
        bsa_calc = BSACalculator(self._obj_name, self._ligand_resn)
        
        # 使用三元复合物公式计算 BSA
        bsa_results = bsa_calc.calculate_bsa_ternary(e3_chain, poi_chain)
        
        interface.bsa_mg_e3 = bsa_results.get('bsa_mg_e3', 0.0)
        interface.bsa_mg_poi = bsa_results.get('bsa_mg_poi', 0.0)
        interface.bsa_e3_poi = bsa_results.get('bsa_e3_poi', 0.0)
        interface.bsa_total = bsa_results.get('bsa_total', 0.0)
        
        # 接触数统计（仍使用原子坐标）
        interface.contact_count_45 = self._count_contacts(mg.atoms, e3.atoms + poi.atoms, 4.5)
        interface.contact_count_50 = self._count_contacts(mg.atoms, e3.atoms + poi.atoms, 5.0)
        
        # 最小链间距离
        interface.min_inter_chain_dist = self._min_distance(e3.atoms, poi.atoms)
        
        # 返回中间计算值
        intermediate = {
            'sa_e3_mg': bsa_results.get('sa_e3_mg', 0.0),
            'sa_poi_mg': bsa_results.get('sa_poi_mg', 0.0),
            'sa_ternary': bsa_results.get('sa_ternary', 0.0),
            'sa_e3': bsa_results.get('sa_e3', 0.0),
            'sa_poi': bsa_results.get('sa_poi', 0.0),
            'sa_e3_poi': bsa_results.get('sa_e3_poi', 0.0),
            'contact_mg_e3': bsa_results.get('contact_mg_e3', 0),
            'contact_mg_poi': bsa_results.get('contact_mg_poi', 0),
        }
        
        return interface, intermediate

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
        # 获取中间计算值
        intermediate = getattr(self, '_intermediate_values', {})
        
        return {
            # Interface
            'bsa_total': features.interface.bsa_total,
            'bsa_mg_e3': features.interface.bsa_mg_e3,
            'bsa_mg_poi': features.interface.bsa_mg_poi,
            'bsa_e3_poi': features.interface.bsa_e3_poi,
            'contact_count_45': features.interface.contact_count_45,
            'contact_count_50': features.interface.contact_count_50,
            'min_inter_chain_dist': features.interface.min_inter_chain_dist,
            
            # 中间计算值（用于显示计算过程）
            'sa_e3_mg': intermediate.get('sa_e3_mg', 0.0),
            'sa_poi_mg': intermediate.get('sa_poi_mg', 0.0),
            'sa_ternary': intermediate.get('sa_ternary', 0.0),
            'sa_e3': intermediate.get('sa_e3', 0.0),
            'sa_poi': intermediate.get('sa_poi', 0.0),
            'sa_e3_poi_complex': intermediate.get('sa_e3_poi', 0.0),
            'contact_mg_e3': intermediate.get('contact_mg_e3', 0),
            'contact_mg_poi': intermediate.get('contact_mg_poi', 0),
            
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
    ligand_smiles: Optional[str] = None,
    obj_name: Optional[str] = None
) -> Dict[str, Any]:
    """
    评估三元复合物的便捷函数
    
    Args:
        pdb_path: PDB文件路径
        e3_chain: E3链ID
        poi_chain: POI链ID
        ligand_chain: 配体链ID
        ligand_smiles: 配体SMILES（可选）
        obj_name: PyMOL对象名（可选）
    
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
    features = evaluator.evaluate(pdb_path, e3_chain, poi_chain, ligand_chain, ligand_smiles, obj_name)
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