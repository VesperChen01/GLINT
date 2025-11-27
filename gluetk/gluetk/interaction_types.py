# -*- coding: utf-8 -*-
"""
interaction_types.py
相互作用数据结构定义

基于 PLIP 设计理念，使用 namedtuple 定义各类相互作用
"""

from __future__ import annotations
from collections import namedtuple
from typing import List, Dict, Any, Optional, Tuple


# ============================================================================
# 相互作用 namedtuple 定义（参考 PLIP）
# ============================================================================

# 氢键
HBond = namedtuple('HBond', [
    'donor_chain', 'donor_res', 'donor_resi', 'donor_atom',
    'acceptor_chain', 'acceptor_res', 'acceptor_resi', 'acceptor_atom',
    'distance', 'angle_dha', 'donor_type', 'acceptor_type'
])

# 疏水相互作用
Hydrophobic = namedtuple('Hydrophobic', [
    'chain1', 'res1', 'resi1', 'atom1',
    'chain2', 'res2', 'resi2', 'atom2',
    'distance', 'contact_type'  # 'aliphatic' or 'aromatic'
])

# 盐桥/离子相互作用
SaltBridge = namedtuple('SaltBridge', [
    'positive_chain', 'positive_res', 'positive_resi',
    'negative_chain', 'negative_res', 'negative_resi',
    'distance', 'center_positive', 'center_negative'
])

# π-π堆积
PiStacking = namedtuple('PiStacking', [
    'chain1', 'res1', 'resi1', 'ring1_center',
    'chain2', 'res2', 'resi2', 'ring2_center',
    'distance', 'angle', 'offset', 'stacking_type'  # 'face-face', 'edge-face', 'T-shaped'
])

# π-阳离子
PiCation = namedtuple('PiCation', [
    'ring_chain', 'ring_res', 'ring_resi', 'ring_center',
    'cation_chain', 'cation_res', 'cation_resi', 'cation_center',
    'distance', 'offset'
])

# 金属配位
MetalCoordination = namedtuple('MetalCoordination', [
    'metal_chain', 'metal_res', 'metal_resi', 'metal_atom', 'metal_type',
    'ligand_chain', 'ligand_res', 'ligand_resi', 'ligand_atom', 'ligand_type',
    'distance', 'coordination_num'
])

# 卤素键
HalogenBond = namedtuple('HalogenBond', [
    'donor_chain', 'donor_res', 'donor_resi', 'donor_atom', 'halogen_type',
    'acceptor_chain', 'acceptor_res', 'acceptor_resi', 'acceptor_atom',
    'distance', 'angle', 'acceptor_type'
])

# 水桥
WaterBridge = namedtuple('WaterBridge', [
    'donor_chain', 'donor_res', 'donor_resi', 'donor_atom',
    'acceptor_chain', 'acceptor_res', 'acceptor_resi', 'acceptor_atom',
    'water_chain', 'water_resi', 'water_coord',
    'distance_dw', 'distance_wa', 'angle'
])


# ============================================================================
# 相互作用结果类
# ============================================================================

class InteractionResult:
    """
    相互作用检测结果容器
    
    参考 PLIP 设计，提供统一的结果管理接口
    
    Attributes:
        interaction_type (str): 相互作用类型 'PP'/'PL'/'PN'/'LL'
        hydrogen_bonds (List[HBond]): 氢键列表
        hydrophobic (List[Hydrophobic]): 疏水相互作用列表
        salt_bridges (List[SaltBridge]): 盐桥列表
        pi_stackings (List[PiStacking]): π-π堆积列表
        pi_cations (List[PiCation]): π-阳离子列表
        metal_coordinations (List[MetalCoordination]): 金属配位列表
        halogen_bonds (List[HalogenBond]): 卤素键列表
        water_bridges (List[WaterBridge]): 水桥列表
    """
    
    def __init__(self, interaction_type: str = 'UNKNOWN'):
        self.interaction_type = interaction_type
        self.hydrogen_bonds: List[HBond] = []
        self.hydrophobic: List[Hydrophobic] = []
        self.salt_bridges: List[SaltBridge] = []
        self.pi_stackings: List[PiStacking] = []
        self.pi_cations: List[PiCation] = []
        self.metal_coordinations: List[MetalCoordination] = []
        self.halogen_bonds: List[HalogenBond] = []
        self.water_bridges: List[WaterBridge] = []
        
        # 元数据
        self.metadata: Dict[str, Any] = {}
    
    def to_dict(self) -> Dict[str, Any]:
        """
        转换为字典格式（向后兼容旧API）
        
        Returns:
            dict: 包含所有相互作用的字典
        """
        return {
            'interaction_type': self.interaction_type,
            'hydrogen_bonds': self._convert_to_old_format(self.hydrogen_bonds, 'hbond'),
            'hydrophobic': self._convert_to_old_format(self.hydrophobic, 'hydrophobic'),
            'salt_bridges': self._convert_to_old_format(self.salt_bridges, 'saltbridge'),
            'pi_stacking': self._convert_to_old_format(self.pi_stackings, 'pi_stacking'),
            'pi_cation': self._convert_to_old_format(self.pi_cations, 'pi_cation'),
            'metal_complexation': self._convert_to_old_format(self.metal_coordinations, 'metal'),
            'halogen_bonds': self._convert_to_old_format(self.halogen_bonds, 'halogen'),
            'water_bridges': self._convert_to_old_format(self.water_bridges, 'water_bridge'),
        }
    
    def _convert_to_old_format(self, interactions: List, interaction_type: str) -> List[Dict]:
        """
        将 namedtuple 转换为旧的 dict 格式
        
        Args:
            interactions: namedtuple 列表
            interaction_type: 相互作用类型
        
        Returns:
            List[Dict]: 字典列表
        """
        if not interactions:
            return []
        
        return [self._namedtuple_to_dict(item) for item in interactions]
    
    def _namedtuple_to_dict(self, item) -> Dict:
        """将 namedtuple 转换为 dict"""
        if hasattr(item, '_asdict'):
            return item._asdict()
        return dict(item)
    
    def summary(self) -> Dict[str, Any]:
        """
        返回统计摘要
        
        Returns:
            dict: 包含各类相互作用数量的字典
        """
        return {
            'interaction_type': self.interaction_type,
            'total_interactions': self.total_count(),
            'hbond_count': len(self.hydrogen_bonds),
            'hydrophobic_count': len(self.hydrophobic),
            'saltbridge_count': len(self.salt_bridges),
            'pi_stacking_count': len(self.pi_stackings),
            'pi_cation_count': len(self.pi_cations),
            'metal_coord_count': len(self.metal_coordinations),
            'halogen_bond_count': len(self.halogen_bonds),
            'water_bridge_count': len(self.water_bridges),
        }
    
    def total_count(self) -> int:
        """返回所有相互作用总数"""
        return sum([
            len(self.hydrogen_bonds),
            len(self.hydrophobic),
            len(self.salt_bridges),
            len(self.pi_stackings),
            len(self.pi_cations),
            len(self.metal_coordinations),
            len(self.halogen_bonds),
            len(self.water_bridges),
        ])
    
    def get_all_interactions(self) -> List[Tuple[str, Any]]:
        """
        获取所有相互作用（带类型标签）
        
        Returns:
            List[Tuple[str, namedtuple]]: (类型名, 相互作用对象) 列表
        """
        all_interactions = []
        
        for hb in self.hydrogen_bonds:
            all_interactions.append(('hydrogen_bond', hb))
        for hydro in self.hydrophobic:
            all_interactions.append(('hydrophobic', hydro))
        for sb in self.salt_bridges:
            all_interactions.append(('salt_bridge', sb))
        for ps in self.pi_stackings:
            all_interactions.append(('pi_stacking', ps))
        for pc in self.pi_cations:
            all_interactions.append(('pi_cation', pc))
        for mc in self.metal_coordinations:
            all_interactions.append(('metal_coordination', mc))
        for hb in self.halogen_bonds:
            all_interactions.append(('halogen_bond', hb))
        for wb in self.water_bridges:
            all_interactions.append(('water_bridge', wb))
        
        return all_interactions
    
    def filter_by_distance(self, max_distance: float) -> 'InteractionResult':
        """
        按距离过滤相互作用
        
        Args:
            max_distance: 最大距离（Å）
        
        Returns:
            InteractionResult: 新的结果对象
        """
        filtered = InteractionResult(self.interaction_type)
        
        filtered.hydrogen_bonds = [hb for hb in self.hydrogen_bonds if hb.distance <= max_distance]
        filtered.hydrophobic = [h for h in self.hydrophobic if h.distance <= max_distance]
        filtered.salt_bridges = [sb for sb in self.salt_bridges if sb.distance <= max_distance]
        filtered.pi_stackings = [ps for ps in self.pi_stackings if ps.distance <= max_distance]
        filtered.pi_cations = [pc for pc in self.pi_cations if pc.distance <= max_distance]
        filtered.metal_coordinations = [mc for mc in self.metal_coordinations if mc.distance <= max_distance]
        filtered.halogen_bonds = [hb for hb in self.halogen_bonds if hb.distance <= max_distance]
        
        return filtered
    
    def filter_by_residue(self, residue_key: Tuple[str, str, str]) -> 'InteractionResult':
        """
        过滤包含特定残基的相互作用
        
        Args:
            residue_key: (chain, resname, resi)
        
        Returns:
            InteractionResult: 新的结果对象
        """
        filtered = InteractionResult(self.interaction_type)
        chain, resname, resi = residue_key
        
        # 氢键
        for hb in self.hydrogen_bonds:
            if ((hb.donor_chain == chain and hb.donor_res == resname and hb.donor_resi == resi) or
                (hb.acceptor_chain == chain and hb.acceptor_res == resname and hb.acceptor_resi == resi)):
                filtered.hydrogen_bonds.append(hb)
        
        # 疏水
        for h in self.hydrophobic:
            if ((h.chain1 == chain and h.res1 == resname and h.resi1 == resi) or
                (h.chain2 == chain and h.res2 == resname and h.resi2 == resi)):
                filtered.hydrophobic.append(h)
        
        # ... 其他类型类似
        
        return filtered
    
    def __str__(self) -> str:
        """字符串表示"""
        return (f"InteractionResult(type={self.interaction_type}, "
                f"total={self.total_count()}, "
                f"H-bonds={len(self.hydrogen_bonds)}, "
                f"hydrophobic={len(self.hydrophobic)}, "
                f"salt-bridges={len(self.salt_bridges)})")
    
    def __repr__(self) -> str:
        return self.__str__()


# ============================================================================
# 辅助函数
# ============================================================================

def format_residue_id(chain: str, resname: str, resi: str) -> str:
    """
    格式化残基ID
    
    Args:
        chain: 链ID
        resname: 残基名
        resi: 残基编号
    
    Returns:
        str: 格式化的残基ID，如 "A:ARG:123"
    """
    return f"{chain}:{resname}:{resi}"


def parse_residue_id(residue_id: str) -> Tuple[str, str, str]:
    """
    解析残基ID字符串
    
    Args:
        residue_id: 格式化的残基ID，如 "A:ARG:123"
    
    Returns:
        Tuple[str, str, str]: (chain, resname, resi)
    """
    parts = residue_id.split(':')
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    return '', '', ''


def interaction_to_csv_row(interaction_type: str, interaction) -> Dict[str, str]:
    """
    将相互作用转换为CSV行格式
    
    Args:
        interaction_type: 相互作用类型
        interaction: 相互作用 namedtuple
    
    Returns:
        Dict[str, str]: CSV行字典
    """
    if interaction_type == 'hydrogen_bond':
        return {
            'Chain1': interaction.donor_chain,
            'Residue1': f"{interaction.donor_res} {interaction.donor_resi}",
            'Atom1': interaction.donor_atom,
            'Chain2': interaction.acceptor_chain,
            'Residue2': f"{interaction.acceptor_res} {interaction.acceptor_resi}",
            'Atom2': interaction.acceptor_atom,
            'Distance': f"{interaction.distance:.2f}",
            'Angle': f"{interaction.angle_dha:.1f}" if interaction.angle_dha else '',
            'Interaction': 'Hydrogen Bond',
            'Confidence': f"{getattr(interaction, 'confidence', 1.0):.2f}"
        }
    
    elif interaction_type == 'salt_bridge':
        return {
            'Chain1': interaction.positive_chain,
            'Residue1': f"{interaction.positive_res} {interaction.positive_resi}",
            'Chain2': interaction.negative_chain,
            'Residue2': f"{interaction.negative_res} {interaction.negative_resi}",
            'Distance': f"{interaction.distance:.2f}",
            'Interaction': 'Salt Bridge',
            'Confidence': f"{getattr(interaction, 'confidence', 1.0):.2f}"
        }
    
    elif interaction_type == 'hydrophobic':
        return {
            'Chain1': interaction.chain1,
            'Residue1': f"{interaction.res1} {interaction.resi1}",
            'Atom1': interaction.atom1,
            'Chain2': interaction.chain2,
            'Residue2': f"{interaction.res2} {interaction.resi2}",
            'Atom2': interaction.atom2,
            'Distance': f"{interaction.distance:.2f}",
            'Interaction': f"Hydrophobic ({interaction.contact_type})",
            'Confidence': f"{getattr(interaction, 'confidence', 1.0):.2f}"
        }
        
    elif interaction_type == 'halogen_bond':
        return {
            'Chain1': interaction.donor_chain,
            'Residue1': f"{interaction.donor_res} {interaction.donor_resi}",
            'Atom1': interaction.donor_atom,
            'Chain2': interaction.acceptor_chain,
            'Residue2': f"{interaction.acceptor_res} {interaction.acceptor_resi}",
            'Atom2': interaction.acceptor_atom,
            'Distance': f"{interaction.distance:.2f}",
            'Angle': f"{interaction.angle:.1f}",
            'Interaction': 'Halogen Bond',
            'Confidence': f"{getattr(interaction, 'confidence', 1.0):.2f}"
        }

    elif interaction_type == 'water_bridge':
        return {
            'Chain1': interaction.donor_chain,
            'Residue1': f"{interaction.donor_res} {interaction.donor_resi}",
            'Atom1': interaction.donor_atom,
            'Chain2': interaction.acceptor_chain,
            'Residue2': f"{interaction.acceptor_res} {interaction.acceptor_resi}",
            'Atom2': interaction.acceptor_atom,
            'Distance': f"{interaction.distance_dw:.2f}/{interaction.distance_wa:.2f}",
            'Angle': f"{interaction.angle:.1f}",
            'Interaction': 'Water Bridge',
            'Confidence': f"{getattr(interaction, 'confidence', 1.0):.2f}"
        }
    
    return {}
