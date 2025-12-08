# -*- coding: utf-8 -*-
"""
ligand_ligand_analyzer.py
小分子-小分子相互作用分析模块

功能：
1. 分析两个选定的小分子（或任意原子集合）之间的相互作用
2. 支持氢键、卤素键、金属配位（基于原子几何）
3. 如果安装了 RDKit，支持 π-π 堆积、π-阳离子、疏水相互作用（基于拓扑）
4. 生成可视化结果
"""

from __future__ import print_function
import os
import csv
import math
import tempfile
from collections import defaultdict
from pymol import cmd

# 复用 interaction_analyzer 的基础几何函数
from .interaction_analyzer import (
    parse_pdb_structure, distance, is_hbond_precise, is_halogen_bond,
    is_metal_coordination, calculate_confidence_score, centroid,
    get_element_from_atom_name, normal_vector, angle_between,
    calculate_angle_three_points
)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

def analyze_ligand_ligand_interactions(obj_name, sel1, sel2, cutoff=4.5, output_csv=None,
                                     visualize=True):
    """
    分析两个选择区域（通常是两个小分子）之间的相互作用
    
    参数:
        obj_name: PyMOL对象名
        sel1: 选择区域1 (PyMOL selection string, e.g. "resn LIG and chain A")
        sel2: 选择区域2 (PyMOL selection string, e.g. "resn DRG and chain A")
        cutoff: 距离截断 (Å)
        output_csv: 输出CSV路径
        visualize: 是否自动在PyMOL中显示
    """
    
    # 1. 获取原子信息
    # 我们使用 cmd.get_model 来获取选择区域的原子，这比 parse_pdb_structure 更灵活
    model1 = cmd.get_model(sel1)
    model2 = cmd.get_model(sel2)
    
    atoms1 = _model_to_list(model1)
    atoms2 = _model_to_list(model2)
    
    if not atoms1 or not atoms2:
        print(f"[GlueTK] ⚠️ Empty selection. Sel1: {len(atoms1)} atoms, Sel2: {len(atoms2)} atoms.")
        return []

    print(f"[GlueTK] Analyzing interactions between {len(atoms1)} atoms (Sel1) and {len(atoms2)} atoms (Sel2)...")

    interactions = []
    
    # 构建 lookup table 用于 is_hbond_precise
    # interaction_analyzer 期望的格式: {(chain, resn, resi): [atoms...]}
    all_atoms_dict = defaultdict(list)
    for a in atoms1 + atoms2:
        key = (a[0], a[1], a[2]) # chain, resn, resi
        all_atoms_dict[key].append(a)
        
    # 2. RDKit 高级分析 (如果可用)
    rdkit_features1 = None
    rdkit_features2 = None
    
    if RDKIT_AVAILABLE:
        try:
            rdkit_features1 = _analyze_with_rdkit(sel1)
            rdkit_features2 = _analyze_with_rdkit(sel2)
            print("[GlueTK] ✅ RDKit analysis successful (Rings, Aromaticity, Hydrophobicity detected)")
        except Exception as e:
            print(f"[GlueTK] ⚠️ RDKit analysis failed: {e}. Falling back to geometric only.")
            
    # 3. 遍历原子对进行分析
    # 优化：使用空间划分或简单距离预筛选
    
    # 提取 RDKit 特征 (如果可用)
    rings1 = rdkit_features1.get('rings', []) if rdkit_features1 else []
    rings2 = rdkit_features2.get('rings', []) if rdkit_features2 else []
    
    # --- A. 基于原子的相互作用 (HBond, Halogen, Metal, Basic VdW) ---
    for a1 in atoms1:
        for a2 in atoms2:
            d = distance(a1[4], a2[4])
            
            if d > cutoff:
                continue
                
            # 3.1 氢键 (使用精确模式)
            # is_hbond_precise 需要 (chain, resn, resi, name, coord)
            is_hb, hb_dist, hb_angle = is_hbond_precise(a1, a2, all_atoms_dict)
            if not is_hb:
                is_hb, hb_dist, hb_angle = is_hbond_precise(a2, a1, all_atoms_dict) # 反向
                
            if is_hb:
                conf = calculate_confidence_score("氢键", hb_dist, hb_angle)
                interactions.append({
                    "Type": "Hydrogen Bond",
                    "Atom1": a1, "Atom2": a2, "Distance": hb_dist, "Confidence": conf,
                    "Details": f"Angle: {hb_dist:.1f}°"
                })
                continue # 也就是优先判定为氢键
                
            # 3.2 卤素键
            is_hal, hal_dist, hal_angle = is_halogen_bond(a1, a2, d)
            if is_hal:
                conf = calculate_confidence_score("卤素键", hal_dist, hal_angle)
                interactions.append({
                    "Type": "Halogen Bond",
                    "Atom1": a1, "Atom2": a2, "Distance": hal_dist, "Confidence": conf,
                    "Details": ""
                })
                continue
                
            # 3.3 金属配位
            is_met, met_dist = is_metal_coordination(a1, a2)
            if is_met:
                interactions.append({
                    "Type": "Metal Coordination",
                    "Atom1": a1, "Atom2": a2, "Distance": met_dist, "Confidence": 1.0,
                    "Details": ""
                })
                continue
                
            # 3.4 弱极性相互作用 (Polar Contact) - Catch-all for non-ideal H-bonds
            # 当未通过严格氢键判定，但距离较近的极性原子对
            if d <= 3.6:
                e1 = get_element_from_atom_name(a1[3])
                e2 = get_element_from_atom_name(a2[3])
                polar_atoms = {'N', 'O', 'S', 'F', 'CL', 'BR', 'I'}
                if e1 in polar_atoms and e2 in polar_atoms:
                     interactions.append({
                        "Type": "Polar Contact",
                        "Atom1": a1, "Atom2": a2, "Distance": d, "Confidence": 0.4,
                        "Details": "Dipole-Dipole"
                    })
                     continue

            # 3.5 疏水相互作用 (如果有 RDKit 支持)
            if rdkit_features1 and rdkit_features2:
                if _is_hydrophobic_atom(a1, rdkit_features1) and _is_hydrophobic_atom(a2, rdkit_features2):
                    if d <= 4.0: # 疏水截断通常 3.8-4.0
                         interactions.append({
                            "Type": "Hydrophobic",
                            "Atom1": a1, "Atom2": a2, "Distance": d, "Confidence": 0.8,
                            "Details": "RDKit-verified"
                        })
            # 无 RDKit 时的疏水回退 (仅 C-C)
            elif not RDKIT_AVAILABLE:
                elem1 = get_element_from_atom_name(a1[3])
                elem2 = get_element_from_atom_name(a2[3])
                if elem1 == 'C' and elem2 == 'C' and d <= 3.6:  # 收紧阈值 3.8 -> 3.6
                    interactions.append({
                        "Type": "Hydrophobic (Generic)",
                        "Atom1": a1, "Atom2": a2, "Distance": d, "Confidence": 0.5,
                        "Details": "C-C contact"
                    })

    # --- B. 基于基团的相互作用 (Pi-Pi, Pi-Cation) ---
    if RDKIT_AVAILABLE and rings1 and rings2:
        # Pi-Pi
        for r1 in rings1:
            if not r1['is_aromatic']: continue
            for r2 in rings2:
                if not r2['is_aromatic']: continue
                
                c1 = r1['centroid']
                c2 = r2['centroid']
                d = distance(c1, c2)
                
                if 3.0 <= d <= 5.5:
                    n1 = r1['normal']
                    n2 = r2['normal']
                    angle = angle_between(n1, n2)
                    
                    # 面-面 (0-30 或 150-180) 或 边-面 (60-120)
                    itype = None
                    if angle <= 30 or angle >= 150:
                        itype = "Pi-Pi Stacking (Face-to-Face)"
                    elif 60 <= angle <= 120:
                        itype = "Pi-Pi Stacking (Edge-to-Face)"
                        
                    if itype:
                        interactions.append({
                            "Type": itype,
                            "Atom1": _dummy_atom_from_ring(r1, atoms1), # 用于可视化的代表原子
                            "Atom2": _dummy_atom_from_ring(r2, atoms2),
                            "Distance": d,
                            "Confidence": 0.9,
                            "Details": f"Angle: {angle:.1f}°",
                            "IsGroup": True,
                            "Group1": r1, "Group2": r2
                        })

    # 4. 输出与可视化
    if output_csv:
        _save_csv(interactions, output_csv)
        
    if visualize:
        visualize_ligand_interactions(obj_name, interactions)
        
    return interactions

def visualize_ligand_interactions(obj_name, interactions):
    """在 PyMOL 中可视化结果"""
    cmd.delete(f"{obj_name}_LL_inter_*")
    
    for i, inter in enumerate(interactions):
        name = f"{obj_name}_LL_inter_{i+1}"
        
        if inter.get("IsGroup"):
            # 绘制环中心连线
            p1 = inter['Group1']['centroid']
            p2 = inter['Group2']['centroid']
            # 创建伪原子
            cmd.pseudoatom("p1", pos=p1)
            cmd.pseudoatom("p2", pos=p2)
            cmd.distance(name, "p1", "p2")
            cmd.delete("p1")
            cmd.delete("p2")
            cmd.color("magenta", name)
        else:
            # 原子连线
            a1 = inter['Atom1']
            a2 = inter['Atom2']
            # selection syntax: chain and resi and name ...
            # 为了精准，最好使用 ID (但 atom tuple 没有 id)。这里使用 chain/resi/name
            sel1 = f"chain {a1[0]} and resn {a1[1]} and resi {a1[2]} and name {a1[3]}"
            sel2 = f"chain {a2[0]} and resn {a2[1]} and resi {a2[2]} and name {a2[3]}"
            
            cmd.distance(name, sel1, sel2)
            
            # 颜色
            color_map = {
                "Hydrogen Bond": "yellow",
                "Halogen Bond": "cyan",
                "Metal Coordination": "magenta",
                "Polar Contact": "orange",
                "Hydrophobic": "green",
                "Hydrophobic (Generic)": "lime",
            }
            col = color_map.get(inter['Type'], "white")
            cmd.color(col, name)
            
    print(f"[GlueTK] Visualized {len(interactions)} interactions.")

def _model_to_list(model):
    """将 PyMOL model 对象转换为 list 格式 [(chain, resn, resi, name, coord), ...]"""
    atoms = []
    for at in model.atom:
        chain = at.chain
        resn = at.resn
        resi = at.resi
        name = at.name
        coord = (at.coord[0], at.coord[1], at.coord[2])
        atoms.append((chain, resn, resi, name, coord))
    return atoms

def _save_csv(interactions, filepath):
    try:
        with open(filepath, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.writer(f)
            writer.writerow(["Type", "Distance", "Confidence", "Details", 
                           "Atom1_Chain", "Atom1_Resn", "Atom1_Resi", "Atom1_Name",
                           "Atom2_Chain", "Atom2_Resn", "Atom2_Resi", "Atom2_Name"])
            for i in interactions:
                a1 = i['Atom1']
                a2 = i['Atom2']
                writer.writerow([
                    i['Type'], f"{i['Distance']:.2f}", i['Confidence'], i['Details'],
                    a1[0], a1[1], a1[2], a1[3],
                    a2[0], a2[1], a2[2], a2[3]
                ])
        print(f"[GlueTK] Results saved to {filepath}")
    except Exception as e:
        print(f"[GlueTK] Failed to save CSV: {e}")

# --- RDKit Helpers ---

def _analyze_with_rdkit(selection):
    """使用 RDKit 分析选择区域的拓扑特征 (Ring, Aromaticity, Hydrophobicity)"""
    # 1. 保存为 SDF/PDB
    # 为了保留原子名称以便后续映射，PDB 比较好。但是 PDB 对小分子连接性支持有时不好。
    # 尝试使用 SDF， PyMOL save 支持
    
    tmp_pdb = tempfile.mktemp(suffix=".pdb")
    cmd.save(tmp_pdb, selection)
    
    mol = Chem.MolFromPDBFile(tmp_pdb, removeHs=False)
    os.remove(tmp_pdb)
    
    if not mol:
        return None
        
    # 2. 获取环信息
    rings = []
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    
    # 获取原子坐标 (Conformer)
    conf = mol.GetConformer()
    
    # 建立 PyMOL atom name 到 RDKit atom idx 的映射? 
    # PDB parser 会保留 PDBInfo。
    
    for ring_atom_indices in atom_rings:
        # 检查芳香性
        is_aromatic = all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring_atom_indices)
        
        # 计算质心和法向量
        coords = [conf.GetAtomPosition(i) for i in ring_atom_indices]
        coords_tuples = [(p.x, p.y, p.z) for p in coords]
        cen = centroid(coords_tuples)
        
        # 法向量 (取前三个点，假设平面)
        if len(coords) >= 3:
            norm = normal_vector(coords_tuples[0], coords_tuples[1], coords_tuples[2])
        else:
            norm = (0,0,1)
            
        # 记录原子名以便回溯
        atom_names = []
        for i in ring_atom_indices:
            info = mol.GetAtomWithIdx(i).GetPDBResidueInfo()
            if info:
                atom_names.append(info.GetName().strip())
        
        rings.append({
            "is_aromatic": is_aromatic,
            "centroid": cen,
            "normal": norm,
            "atom_indices": ring_atom_indices,
            "atom_names": atom_names # 用于后续匹配
        })
        
    # 3. 获取疏水原子
    # 简单定义: C atom not bonded to N, O, F, Cl, Br, I, S(polar?)
    hydrophobic_indices = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            is_polar = False
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() in ['N', 'O', 'F', 'Cl', 'Br', 'I']: # S?
                    is_polar = True
                    break
            if not is_polar:
                hydrophobic_indices.append(atom.GetIdx())

    # 建立 Atom Name -> Feature 映射
    # key: atom_name (assuming unique within the selection/residue context?)
    # 实际上 PyMOL selection 可能包含多个 residues。 RDKit PDB parser creates one molecule.
    # 我们需要一种方法将 RDKit atom index 映射回 atoms list 中的 entry。
    # 使用 (chain, resn, resi, name) 匹配。
    
    features = {
        "rings": rings,
        "hydrophobic_indices": set(hydrophobic_indices),
        "mol": mol
    }
    return features

def _is_hydrophobic_atom(atom_tuple, features):
    """检查原子是否为疏水 (利用 RDKit 结果)"""
    # atom_tuple: (chain, resn, resi, name, coord)
    mol = features['mol']
    hyd_indices = features['hydrophobic_indices']
    
    # 查找对应原子
    # 这是一个耗时操作，但在小分子场景下可接受
    tgt_name = atom_tuple[3]
    tgt_resi = atom_tuple[2]
    tgt_chain = atom_tuple[0]
    
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info:
            # PDB Residue Info parsing might vary
            # RDKit PDB parser usually strips spaces
            if (info.GetName().strip() == tgt_name and 
                str(info.GetResidueNumber()) == str(tgt_resi)): # chain check?
                # Check index
                return atom.GetIdx() in hyd_indices
    return False

def _dummy_atom_from_ring(ring, atom_list):
    """创建一个用于可视化的伪原子tuple，位于环中心"""
    # 我们需要找到环中的一个真实原子作为 'Atom1' 的引用信息 (用于 Chain/Res info)
    # 然后坐标改为 centroid
    
    # Find reference atom
    ref_name = ring['atom_names'][0]
    ref_atom = None
    for a in atom_list:
        if a[3] == ref_name:
            ref_atom = a
            break
            
    if ref_atom:
        # (chain, resn, resi, name, coord)
        return (ref_atom[0], ref_atom[1], ref_atom[2], "RING", ring['centroid'])
    return None
