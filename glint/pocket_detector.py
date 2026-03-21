# -*- coding: utf-8 -*-
"""
GLINT Pocket Detector
======================
纯 Python 实现的蛋白口袋检测与性质分析Module

功能：
- 基于网格法的口袋检测
- 几何性质：体积、表面积、深degrees、开口Size、球形degrees
- 化学性质：疏水性、极性、电荷、芳香性、氢Key供受体
- 可成药性评分
- 与 PPI interface、相互作用、静电势联动

Dependencies：NumPy, SciPy (现有Dependencies)
"""

import os
import sys
import numpy as np
from scipy import ndimage
from scipy.spatial.distance import cdist
import tempfile
import csv

try:
    from pymol import cmd
except ImportError:
    cmd = None

# 氨基酸性质Category
HYDROPHOBIC_RESIDUES = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO', 'GLY'}
POLAR_RESIDUES = {'SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN'}
POSITIVE_RESIDUES = {'LYS', 'ARG', 'HIS'}
NEGATIVE_RESIDUES = {'ASP', 'GLU'}
AROMATIC_RESIDUES = {'PHE', 'TRP', 'TYR', 'HIS'}

# VDW 半径 (Å)
VDW_RADII = {
    'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
    'P': 1.80, 'H': 1.20, 'F': 1.47, 'CL': 1.75,
    'BR': 1.85, 'I': 1.98, 'MG': 1.73, 'CA': 2.31,
    'FE': 2.00, 'ZN': 1.39, 'CU': 1.40, 'MN': 1.61
}

# 氢Key供体/受体原子Type
HBOND_DONOR_ATOMS = {'N', 'O'}  # 简化：带H的N/O
HBOND_ACCEPTOR_ATOMS = {'O', 'N', 'S'}


# 从Package级别Import统一的中文检测Function，避免重复定义
from . import _zh


def _info(cn, en):
    """双语Information输出"""
    print(cn if _zh() else en)


class PocketDetector:
    """口袋检测器（严格标准Parameters）"""
    
    def __init__(self, grid_spacing=0.5, probe_radius=1.4, 
                 min_volume=30.0, min_depth=2.5, max_solvent_access=0.2):
        """
        Parameters：
            grid_spacing: 网格间距 (Å)
            probe_radius: 探针半径 (Å, 模拟水分子)
            min_volume: 最小口袋体积 (ų)
            min_depth: 最小埋藏深degrees (Å)
            max_solvent_access: 最大溶剂可及degrees (0=完全埋藏, 1=完全暴露)
        """
        self.grid_spacing = grid_spacing
        self.probe_radius = probe_radius
        self.min_volume = min_volume
        self.min_depth = min_depth
        self.max_solvent_access = max_solvent_access
        
    def detect_pockets(self, obj_name=None, pdb_file=None, selection='all'):
        """
        检测口袋
        
        Return：[{
            'id': int,
            'volume': float,
            'surface_area': float,
            'depth': float,
            'mouth_size': float,
            'sphericity': float,
            'center': (x, y, z),
            'grid_coords': [(i,j,k), ...],
            'residues': [{'chain': str, 'resn': str, 'resi': str}, ...]
        }, ...]
        """
        # 获取原子坐标
        atoms = self._get_atoms(obj_name, pdb_file, selection)
        if len(atoms) == 0:
            _info("Error：未找到原子", "Error: No atoms found")
            return []
        
        _info(f"检测口袋：{len(atoms)} 个原子", 
              f"Detecting pockets: {len(atoms)} atoms")
        
        # 构建网格
        grid, origin, atom_coords = self._build_grid(atoms)
        
        # 标记蛋白占据的格点
        occupied = self._mark_occupied(grid, atom_coords, origin)
        
        # 标记外部溶剂可及区域
        solvent = self._mark_solvent(occupied)
        
        # 找出口袋（未被占据且未暴露的连通区域）
        cavity = (~occupied) & (~solvent)
        labeled, num_pockets = ndimage.label(cavity)
        
        _info(f"初步检测到 {num_pockets} 个口袋", 
              f"Initially found {num_pockets} pockets")
        
        # 分析每个口袋
        pockets = []
        for pocket_id in range(1, num_pockets + 1):
            pocket_mask = (labeled == pocket_id)
            pocket_info = self._analyze_pocket(
                pocket_id, pocket_mask, grid, origin, 
                occupied, solvent, atoms
            )
            
            # Filter
            if (pocket_info['volume'] >= self.min_volume and 
                pocket_info['depth'] >= self.min_depth and
                pocket_info['solvent_access'] <= self.max_solvent_access):
                pockets.append(pocket_info)
        
        _info(f"Filter后保留 {len(pockets)} 个口袋", 
              f"Retained {len(pockets)} pockets after filtering")
        
        return pockets
    
    def _get_atoms(self, obj_name, pdb_file, selection):
        """获取Atom information"""
        atoms = []
        
        if obj_name and cmd:
            # 从 PyMOL 对象获取
            model = cmd.get_model(f"{obj_name} and {selection}")
            for atom in model.atom:
                atoms.append({
                    'coord': np.array([atom.coord[0], atom.coord[1], atom.coord[2]]),
                    'element': atom.symbol.upper(),
                    'chain': atom.chain.strip(),
                    'resn': atom.resn.strip(),
                    'resi': str(atom.resi).strip() + atom.segi.strip(),
                    'name': atom.name.strip()
                })
        elif pdb_file:
            # 从 PDB File读取
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        try:
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            element = line[76:78].strip().upper()
                            if not element:
                                element = line[12:16].strip()[0].upper()
                            atoms.append({
                                'coord': np.array([x, y, z]),
                                'element': element,
                                'chain': line[21:22].strip(),
                                'resn': line[17:20].strip(),
                                'resi': line[22:27].strip(),
                                'name': line[12:16].strip()
                            })
                        except (ValueError, IndexError):  # PDB 行解析可能Failed
                            continue
        
        return atoms
    
    def _build_grid(self, atoms):
        """构建网格"""
        coords = np.array([a['coord'] for a in atoms])
        
        # 边界框 + padding
        padding = self.probe_radius + 5.0
        min_coords = coords.min(axis=0) - padding
        max_coords = coords.max(axis=0) + padding
        
        origin = min_coords
        dimensions = max_coords - min_coords
        grid_size = (np.ceil(dimensions / self.grid_spacing) + 1).astype(int)
        
        grid = {
            'size': grid_size,
            'spacing': self.grid_spacing,
            'origin': origin,
            'dimensions': dimensions
        }
        
        return grid, origin, coords
    
    def _mark_occupied(self, grid, atom_coords, origin):
        """标记被蛋白占据的格点"""
        occupied = np.zeros(grid['size'], dtype=bool)
        
        # 为每个原子标记其 VDW + probe 范围内的格点
        for coord in atom_coords:
            # 计算原子影响的格点范围
            radius = 2.0 + self.probe_radius  # 保守估计，using通用半径
            grid_coord = ((coord - origin) / self.grid_spacing).astype(int)
            grid_radius = int(np.ceil(radius / self.grid_spacing))
            
            # 遍历邻近格点
            for i in range(max(0, grid_coord[0] - grid_radius),
                          min(grid['size'][0], grid_coord[0] + grid_radius + 1)):
                for j in range(max(0, grid_coord[1] - grid_radius),
                              min(grid['size'][1], grid_coord[1] + grid_radius + 1)):
                    for k in range(max(0, grid_coord[2] - grid_radius),
                                  min(grid['size'][2], grid_coord[2] + grid_radius + 1)):
                        # 计算实际距离
                        grid_pos = origin + np.array([i, j, k]) * self.grid_spacing
                        dist = np.linalg.norm(grid_pos - coord)
                        if dist <= radius + self.probe_radius:
                            occupied[i, j, k] = True
        
        return occupied
    
    def _mark_solvent(self, occupied):
        """标记外部溶剂可及区域（flood fill）"""
        solvent = np.zeros(occupied.shape, dtype=bool)
        
        # 从边界Start flood fill
        # using scipy 的连通性分析
        external = ~occupied
        labeled, num = ndimage.label(external)
        
        # 找出与边界相连的最大连通域
        border_labels = set()
        # 六个面
        border_labels.update(labeled[0, :, :].flat)
        border_labels.update(labeled[-1, :, :].flat)
        border_labels.update(labeled[:, 0, :].flat)
        border_labels.update(labeled[:, -1, :].flat)
        border_labels.update(labeled[:, :, 0].flat)
        border_labels.update(labeled[:, :, -1].flat)
        border_labels.discard(0)
        
        # 标记所有与边界相连的区域为溶剂
        for label in border_labels:
            solvent |= (labeled == label)
        
        return solvent
    
    def _analyze_pocket(self, pocket_id, pocket_mask, grid, origin, 
                       occupied, solvent, atoms):
        """分析单个口袋的性质"""
        pocket_coords = np.argwhere(pocket_mask)
        num_voxels = len(pocket_coords)
        
        # 几何性质
        volume = num_voxels * (self.grid_spacing ** 3)
        
        # 表面积：统计与蛋白接触的边界格点
        # using形态学膨胀
        struct = ndimage.generate_binary_structure(3, 1)  # 6-连通
        dilated = ndimage.binary_dilation(pocket_mask, structure=struct)
        surface_mask = dilated & occupied
        surface_area = np.sum(surface_mask) * (self.grid_spacing ** 2)
        
        # 中心点
        center_grid = pocket_coords.mean(axis=0)
        center = origin + center_grid * self.grid_spacing
        
        # 深degrees：口袋中心到最近溶剂格点的距离
        if np.any(solvent):
            solvent_coords = np.argwhere(solvent)
            distances = cdist([center_grid], solvent_coords)[0]
            depth = distances.min() * self.grid_spacing
        else:
            depth = 0.0
        
        # 开口Size：估计为与溶剂接触的表面积
        mouth_mask = dilated & solvent
        mouth_size = np.sqrt(np.sum(mouth_mask) * (self.grid_spacing ** 2) / np.pi) * 2
        
        # 球形degrees：实际体积 / 等效球体体积
        equivalent_radius = (3 * volume / (4 * np.pi)) ** (1/3)
        bounding_radius = np.max(cdist([center_grid], pocket_coords)) * self.grid_spacing
        sphericity = equivalent_radius / bounding_radius if bounding_radius > 0 else 0
        
        # 溶剂可及degrees
        boundary_mask = dilated & (~pocket_mask)
        boundary_coords = np.argwhere(boundary_mask)
        solvent_contact = np.sum(dilated & solvent)
        protein_contact = np.sum(surface_mask)
        total_contact = solvent_contact + protein_contact
        solvent_access = solvent_contact / total_contact if total_contact > 0 else 1.0
        
        # 邻近残基
        residues = self._find_pocket_residues(pocket_coords, origin, grid, atoms)
        
        return {
            'id': pocket_id,
            'volume': volume,
            'surface_area': surface_area,
            'depth': depth,
            'mouth_size': mouth_size,
            'sphericity': sphericity,
            'solvent_access': solvent_access,
            'center': tuple(center),
            'grid_coords': pocket_coords.tolist(),
            'residues': residues
        }
    
    def _find_pocket_residues(self, pocket_coords, origin, grid, atoms, cutoff=4.0):
        """找出口袋周围的残基"""
        # 口袋格点的实际坐标
        pocket_positions = origin + pocket_coords * self.grid_spacing
        
        # 原子坐标
        atom_coords = np.array([a['coord'] for a in atoms])
        
        # 计算距离矩阵
        distances = cdist(pocket_positions, atom_coords)
        
        # 找出距离 < cutoff 的原子
        close_atom_indices = np.unique(np.where(distances < cutoff)[1])
        
        # 提取残基Information（去重）
        residues = []
        seen = set()
        for idx in close_atom_indices:
            atom = atoms[idx]
            key = (atom['chain'], atom['resn'], atom['resi'])
            if key not in seen:
                seen.add(key)
                residues.append({
                    'chain': atom['chain'],
                    'resn': atom['resn'],
                    'resi': atom['resi']
                })
        
        return residues


class PocketAnalyzer:
    """口袋化学性质分析"""
    
    @staticmethod
    def analyze_properties(pocket, atoms):
        """
        分析口袋化学性质
        
        Return：{
            'hydrophobicity': float,
            'polarity': float,
            'net_charge': int,
            'aromaticity': float,
            'hbond_donors': int,
            'hbond_acceptors': int,
            'druggability_score': float
        }
        """
        residues = pocket['residues']
        resn_list = [r['resn'] for r in residues]
        
        if len(resn_list) == 0:
            return {
                'hydrophobicity': 0.0,
                'polarity': 0.0,
                'net_charge': 0,
                'aromaticity': 0.0,
                'hbond_donors': 0,
                'hbond_acceptors': 0,
                'druggability_score': 0.0
            }
        
        # 疏水性
        hydrophobic_count = sum(1 for r in resn_list if r in HYDROPHOBIC_RESIDUES)
        hydrophobicity = hydrophobic_count / len(resn_list)
        
        # 极性
        polar_count = sum(1 for r in resn_list if r in POLAR_RESIDUES)
        polarity = polar_count / len(resn_list)
        
        # 电荷
        positive_count = sum(1 for r in resn_list if r in POSITIVE_RESIDUES)
        negative_count = sum(1 for r in resn_list if r in NEGATIVE_RESIDUES)
        net_charge = positive_count - negative_count
        
        # 芳香性
        aromatic_count = sum(1 for r in resn_list if r in AROMATIC_RESIDUES)
        aromaticity = aromatic_count / len(resn_list)
        
        # 氢Key供受体（简化统计：从原子Type推断）
        hbond_donors = 0
        hbond_acceptors = 0
        for atom in atoms:
            # 检查是否属于口袋残基
            atom_key = (atom['chain'], atom['resn'], atom['resi'])
            if any(atom_key == (r['chain'], r['resn'], r['resi']) for r in residues):
                element = atom['element']
                if element in HBOND_DONOR_ATOMS:
                    hbond_donors += 1
                if element in HBOND_ACCEPTOR_ATOMS:
                    hbond_acceptors += 1
        
        # 可成药性评分（经验公式）
        # 参考 Fpocket 的评分：考虑体积、疏水性、深degrees、开口
        volume_score = min(pocket['volume'] / 500.0, 1.0)  # 理想体积 300-500 ų
        hydro_score = hydrophobicity * 0.7  # 适degrees疏水性
        depth_score = min(pocket['depth'] / 10.0, 1.0)  # 深degrees > 5Å
        mouth_score = 1.0 - min(pocket['mouth_size'] / 20.0, 1.0)  # 开口不宜过大
        
        druggability_score = (
            volume_score * 0.3 +
            hydro_score * 0.3 +
            depth_score * 0.2 +
            mouth_score * 0.2
        )
        
        return {
            'hydrophobicity': hydrophobicity,
            'polarity': polarity,
            'net_charge': net_charge,
            'aromaticity': aromaticity,
            'hbond_donors': hbond_donors,
            'hbond_acceptors': hbond_acceptors,
            'druggability_score': druggability_score
        }


def detect_pockets(obj_name=None, pdb_file=None, selection='all',
                   grid_spacing=0.5, probe_radius=1.4,
                   min_volume=30.0, min_depth=2.5, max_solvent_access=0.2,
                   output_csv=None):
    """
    检测蛋白口袋（严格标准Parameters）
    
    Parameters：
        obj_name: PyMOL 对象名
        pdb_file: PDB FilePath
        selection: PyMOL Select语法
        grid_spacing: 网格间距 (Å, 默认 0.5 高精degrees)
        probe_radius: 探针半径 (Å, 1.4 模拟水分子)
        min_volume: 最小口袋体积 (Å³, 默认 30)
        min_depth: 最小埋藏深degrees (Å, 默认 2.5)
        max_solvent_access: 最大溶剂可及degrees (0-1, 默认 0.2 高Select性)
        output_csv: 输出 CSV FilePath
    
    Return：口袋列表
    """
    
    # 检测
    detector = PocketDetector(
        grid_spacing=grid_spacing,
        probe_radius=probe_radius,
        min_volume=min_volume,
        min_depth=min_depth,
        max_solvent_access=max_solvent_access
    )
    
    pockets = detector.detect_pockets(obj_name, pdb_file, selection)
    
    if len(pockets) == 0:
        _info("未检测到符合条件的口袋", "No pockets found")
        return []
    
    # 分析化学性质
    atoms = detector._get_atoms(obj_name, pdb_file, selection)
    analyzer = PocketAnalyzer()
    
    for pocket in pockets:
        props = analyzer.analyze_properties(pocket, atoms)
        pocket.update(props)
    
    # 输出 CSV
    if output_csv:
        _export_csv(pockets, output_csv)
    
    # Print摘要
    _info(f"\n检测到 {len(pockets)} 个口袋：", f"\nDetected {len(pockets)} pockets:")
    for p in pockets:
        print(f"  Pocket {p['id']}: "
              f"V={p['volume']:.1f}ų, "
              f"SA={p['surface_area']:.1f}ų, "
              f"Depth={p['depth']:.1f}Å, "
              f"Hydro={p['hydrophobicity']:.2f}, "
              f"Drug={p['druggability_score']:.2f}")
    
    return pockets


def _export_csv(pockets, output_csv):
    """Export CSV"""
    fieldnames = [
        'Pocket_ID', 'Volume_A3', 'Surface_Area_A2', 'Depth_A',
        'Mouth_Size_A', 'Sphericity', 'Solvent_Access',
        'Center_X', 'Center_Y', 'Center_Z',
        'Hydrophobicity', 'Polarity', 'Net_Charge', 'Aromaticity',
        'HBond_Donors', 'HBond_Acceptors', 'Druggability_Score',
        'Residues'
    ]
    
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for p in pockets:
            residues_str = ';'.join([f"{r['chain']}:{r['resn']}:{r['resi']}" 
                                    for r in p['residues']])
            writer.writerow({
                'Pocket_ID': p['id'],
                'Volume_A3': f"{p['volume']:.2f}",
                'Surface_Area_A2': f"{p['surface_area']:.2f}",
                'Depth_A': f"{p['depth']:.2f}",
                'Mouth_Size_A': f"{p['mouth_size']:.2f}",
                'Sphericity': f"{p['sphericity']:.3f}",
                'Solvent_Access': f"{p['solvent_access']:.3f}",
                'Center_X': f"{p['center'][0]:.2f}",
                'Center_Y': f"{p['center'][1]:.2f}",
                'Center_Z': f"{p['center'][2]:.2f}",
                'Hydrophobicity': f"{p['hydrophobicity']:.3f}",
                'Polarity': f"{p['polarity']:.3f}",
                'Net_Charge': p['net_charge'],
                'Aromaticity': f"{p['aromaticity']:.3f}",
                'HBond_Donors': p['hbond_donors'],
                'HBond_Acceptors': p['hbond_acceptors'],
                'Druggability_Score': f"{p['druggability_score']:.3f}",
                'Residues': residues_str
            })
    
    _info(f"已Save到 {output_csv}", f"Saved to {output_csv}")


def compare_pockets(obj_a, obj_b, align=True, output_csv=None, **kwargs):
    """
    对比两个结构的口袋（如分子胶前后）
    
    Parameters：
        obj_a, obj_b: PyMOL 对象名或 PDB File
        align: 是否先对齐
        output_csv: 输出对比Results CSV
        **kwargs: 传递给 detect_pockets 的Parameters
    
    Return：(pockets_a, pockets_b, comparison)
    """
    # 对齐
    if align and cmd:
        try:
            cmd.align(obj_b, obj_a)
            _info(f"已对齐 {obj_b} 到 {obj_a}", f"Aligned {obj_b} to {obj_a}")
        except Exception:  # PyMOL 对齐可能Failed
            pass
    
    # 检测口袋
    _info(f"\n分析 {obj_a}...", f"\nAnalyzing {obj_a}...")
    pockets_a = detect_pockets(obj_name=obj_a, **kwargs)
    
    _info(f"\n分析 {obj_b}...", f"\nAnalyzing {obj_b}...")
    pockets_b = detect_pockets(obj_name=obj_b, **kwargs)
    
    # 匹配口袋（基于中心距离）
    comparison = _match_pockets(pockets_a, pockets_b)
    
    # 输出对比
    if output_csv:
        _export_comparison_csv(comparison, output_csv)
    
    # Print摘要
    _info(f"\n对比Results：", f"\nComparison results:")
    print(f"  {obj_a}: {len(pockets_a)} pockets")
    print(f"  {obj_b}: {len(pockets_b)} pockets")
    print(f"  Matched: {sum(1 for c in comparison if c['match_type'] == 'matched')}")
    print(f"  New in {obj_b}: {sum(1 for c in comparison if c['match_type'] == 'new')}")
    print(f"  Lost from {obj_a}: {sum(1 for c in comparison if c['match_type'] == 'lost')}")
    
    return pockets_a, pockets_b, comparison


def _match_pockets(pockets_a, pockets_b, distance_threshold=5.0):
    """匹配两组口袋"""
    comparison = []
    matched_b = set()
    
    # A 中的口袋找 B 中的匹配
    for pa in pockets_a:
        best_match = None
        best_dist = float('inf')
        
        for i, pb in enumerate(pockets_b):
            if i in matched_b:
                continue
            dist = np.linalg.norm(np.array(pa['center']) - np.array(pb['center']))
            if dist < best_dist and dist < distance_threshold:
                best_dist = dist
                best_match = i
        
        if best_match is not None:
            pb = pockets_b[best_match]
            matched_b.add(best_match)
            comparison.append({
                'match_type': 'matched',
                'pocket_a_id': pa['id'],
                'pocket_b_id': pb['id'],
                'distance': best_dist,
                'delta_volume': pb['volume'] - pa['volume'],
                'delta_surface_area': pb['surface_area'] - pa['surface_area'],
                'delta_druggability': pb['druggability_score'] - pa['druggability_score']
            })
        else:
            comparison.append({
                'match_type': 'lost',
                'pocket_a_id': pa['id'],
                'pocket_b_id': None,
                'distance': None,
                'delta_volume': -pa['volume'],
                'delta_surface_area': -pa['surface_area'],
                'delta_druggability': -pa['druggability_score']
            })
    
    # B 中新出现的口袋
    for i, pb in enumerate(pockets_b):
        if i not in matched_b:
            comparison.append({
                'match_type': 'new',
                'pocket_a_id': None,
                'pocket_b_id': pb['id'],
                'distance': None,
                'delta_volume': pb['volume'],
                'delta_surface_area': pb['surface_area'],
                'delta_druggability': pb['druggability_score']
            })
    
    return comparison


def _export_comparison_csv(comparison, output_csv):
    """Export对比Results"""
    fieldnames = [
        'Match_Type', 'Pocket_A_ID', 'Pocket_B_ID', 'Distance_A',
        'Delta_Volume_A3', 'Delta_Surface_Area_A2', 'Delta_Druggability'
    ]
    
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for c in comparison:
            writer.writerow({
                'Match_Type': c['match_type'],
                'Pocket_A_ID': c['pocket_a_id'] if c['pocket_a_id'] is not None else '',
                'Pocket_B_ID': c['pocket_b_id'] if c['pocket_b_id'] is not None else '',
                'Distance_A': f"{c['distance']:.2f}" if c['distance'] is not None else '',
                'Delta_Volume_A3': f"{c['delta_volume']:.2f}",
                'Delta_Surface_Area_A2': f"{c['delta_surface_area']:.2f}",
                'Delta_Druggability': f"{c['delta_druggability']:.3f}"
            })
    
    _info(f"对比Results已Save到 {output_csv}", f"Comparison saved to {output_csv}")


if __name__ == '__main__':
    # 测试（需要 PyMOL 环境）
    print("GLINT Pocket Detector - 请在 PyMOL 中using")
    print("示例：detect_pockets('protein', output_csv='pockets.csv')")
