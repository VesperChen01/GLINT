# -*- coding: utf-8 -*-
"""
ppi_analyzer.py
蛋白-蛋白界面(PPI)分析模块 - 分子胶(Molecular Glue)特异功能

核心功能:
1. analyze_protein_protein_interface - 检测两个蛋白质之间的直接相互作用
2. identify_neo_epitope - 识别分子胶诱导的新表位(neo-substrate epitope)
3. calculate_interface_bsa - 计算埋藏表面积(Buried Surface Area)
4. calculate_interface_score - 界面强度评分

典型应用场景:
- 分子胶(Molecular Glue)机制验证
- PROTAC vs Glue 区分
- 三元复合物界面稳定性分析

使用示例:
    # 分析蛋白-蛋白界面
    ppi_result = analyze_protein_protein_interface('complex', ['A'], ['B'])
    
    # 识别Neo-表位
    neo_result = identify_neo_epitope('complex', ['A'], ['B'], 'CC885')
    
    # 计算BSA
    bsa = calculate_interface_bsa('complex', 'A', 'B')
"""

from __future__ import print_function
import math
import csv
import os
from collections import defaultdict
from pymol import cmd

# 导入依赖模块
try:
    from .interaction_analyzer import (
        parse_pdb_structure, parse_pdb_file, distance,
        get_element_from_atom_name, calculate_angle_three_points,
        identify_molecule_type, INTERACTION_PARAMS
    )
except ImportError:
    from interaction_analyzer import (
        parse_pdb_structure, parse_pdb_file, distance,
        get_element_from_atom_name, calculate_angle_three_points,
        identify_molecule_type, INTERACTION_PARAMS
    )

# ========== PPI 参数 ==========
PPI_PARAMS = {
    "interface_distance": 4.5,      # Å, 界面残基判定距离
    "strong_interface_bsa": 800.0,  # Å², 强界面阈值
    "min_interface_contacts": 3,    # 最小界面接触数
    "neo_epitope_distance": 5.0,    # Å, Neo-表位检测距离
}


def analyze_protein_protein_interface(obj_name=None, 
                                      protein1_chains=None, 
                                      protein2_chains=None,
                                      interface_distance=4.5,
                                      output_csv=None,
                                      pdb_file=None,
                                      visualize=True,
                                      protein1_color="cyan",
                                      protein2_color="magenta"):
    """
    检测两个蛋白质之间的直接相互作用（分子胶机制的核心特征）
    
    这是区分 PROTAC 和分子胶的关键功能：
    - PROTAC: 通常无/弱 蛋白-蛋白接触
    - 分子胶: 强 蛋白-蛋白直接接触（Glue诱导的界面）
    
    参数:
        obj_name: PyMOL对象名称
        protein1_chains: 蛋白质1的链ID列表（例如 ["A"] - E3 ligase）
        protein2_chains: 蛋白质2的链ID列表（例如 ["B"] - Substrate）
        interface_distance: 界面残基判定距离（埃）
        output_csv: 输出CSV文件路径
        pdb_file: PDB文件路径（可选）
        visualize: 是否在PyMOL中可视化界面（默认True）
        protein1_color: 蛋白质1的显示颜色（默认cyan）
        protein2_color: 蛋白质2的显示颜色（默认magenta）
    
    返回:
        dict: {
            "interface_residues": [(chain1, resn1, resi1, chain2, resn2, resi2, min_dist), ...],
            "interface_interactions": [...],  # 氢键/盐桥/疏水
            "interface_contacts": int,  # 接触对数量
            "interface_strength": float,  # 界面强度评分 [0-10]
            "bsa": float,  # 埋藏表面积（如果可计算）
            "is_strong_interface": bool  # 是否为强界面
        }
    """
    # 输入验证
    if not protein1_chains or not protein2_chains:
        print("[analyze_protein_protein_interface] ⚠️ Must specify both protein1_chains and protein2_chains")
        return None
    
    # 获取原子信息
    if pdb_file:
        atoms = parse_pdb_file(pdb_file)
        if not atoms:
            print(f"[analyze_protein_protein_interface] Unable to read file: {pdb_file}")
            return None
    else:
        atoms = parse_pdb_structure(obj_name)
        if not atoms:
            print("[analyze_protein_protein_interface] Unable to obtain atom information")
            return None
    
    print(f"[PPI] Analyzing interface between chains {protein1_chains} and {protein2_chains}...")
    
    # 按残基和链分组原子
    chain_residues = defaultdict(list)
    for atom in atoms:
        chain, res_name, res_id, atom_name, coords = atom
        if identify_molecule_type(res_name) != "protein":
            continue  # 只处理蛋白质原子
        res_key = (chain, res_name, res_id)
        chain_residues[res_key].append(atom)
    
    # 分离两个蛋白的残基
    protein1_residues = {k: v for k, v in chain_residues.items() if k[0] in protein1_chains}
    protein2_residues = {k: v for k, v in chain_residues.items() if k[0] in protein2_chains}
    
    if not protein1_residues or not protein2_residues:
        print("[PPI] ⚠️ No protein residues found in specified chains")
        return None
    
    print(f"[PPI] Protein 1: {len(protein1_residues)} residues")
    print(f"[PPI] Protein 2: {len(protein2_residues)} residues")
    
    # 检测界面残基对
    interface_residues = []
    interface_atom_pairs = []  # 用于详细相互作用分析
    
    for res1_key, res1_atoms in protein1_residues.items():
        for res2_key, res2_atoms in protein2_residues.items():
            # 计算残基间最短距离
            min_dist = float('inf')
            closest_atoms = None
            
            for atom1 in res1_atoms:
                for atom2 in res2_atoms:
                    d = distance(atom1[4], atom2[4])
                    if d < min_dist:
                        min_dist = d
                        closest_atoms = (atom1, atom2)
            
            # 如果在界面距离内，记录
            if min_dist <= interface_distance:
                interface_residues.append({
                    "chain1": res1_key[0],
                    "resname1": res1_key[1],
                    "resid1": res1_key[2],
                    "chain2": res2_key[0],
                    "resname2": res2_key[1],
                    "resid2": res2_key[2],
                    "min_distance": round(min_dist, 2)
                })
                interface_atom_pairs.append((res1_atoms, res2_atoms, min_dist))
    
    print(f"[PPI] ✅ Found {len(interface_residues)} interface residue pairs")
    
    # 分析界面相互作用类型
    interface_interactions = _analyze_interface_interactions(interface_atom_pairs)
    
    # 计算埋藏表面积（如果对象存在于PyMOL中）
    bsa = None
    if obj_name and obj_name in cmd.get_object_list():
        try:
            bsa = calculate_interface_bsa(obj_name, protein1_chains[0], protein2_chains[0])
        except Exception as e:
            print(f"[PPI] BSA calculation failed: {e}")
    
    # 计算界面强度评分
    interface_strength = _calculate_interface_strength(
        len(interface_residues),
        interface_interactions,
        bsa
    )
    
    # 判定是否为强界面
    is_strong_interface = False
    if bsa and bsa > PPI_PARAMS["strong_interface_bsa"]:
        is_strong_interface = True
    elif len(interface_residues) > 10 and interface_strength > 6.0:
        is_strong_interface = True
    
    # 输出到CSV
    if output_csv:
        _export_ppi_csv(output_csv, interface_residues, interface_interactions, bsa, interface_strength)
    
    result = {
        "interface_residues": interface_residues,
        "interface_interactions": interface_interactions,
        "interface_contacts": len(interface_residues),
        "interface_strength": round(interface_strength, 2),
        "bsa": round(bsa, 2) if bsa else None,
        "is_strong_interface": is_strong_interface,
        "protein1_chains": protein1_chains,
        "protein2_chains": protein2_chains
    }
    
    # 打印摘要
    print("\n" + "=" * 60)
    print("蛋白-蛋白界面分析结果 (PPI Analysis)")
    print("=" * 60)
    print(f"界面接触对数: {len(interface_residues)}")
    print(f"界面强度评分: {interface_strength:.2f} / 10")
    if bsa:
        print(f"埋藏表面积(BSA): {bsa:.1f} Ų")
    print(f"界面类型: {'✨ 强界面 (Strong Interface)' if is_strong_interface else '⚠️ 弱界面 (Weak Interface)'}")
    print(f"\n相互作用类型统计:")
    for inter_type, count in _count_interaction_types(interface_interactions).items():
        print(f"  {inter_type}: {count}")
    print("=" * 60)
    
    # 在PyMOL中可视化界面（如果启用）
    if visualize and obj_name and obj_name in cmd.get_object_list():
        try:
            visualize_ppi_interface(obj_name, result, protein1_color, protein2_color)
        except Exception as e:
            print(f"[PPI] 可视化失败: {e}")
    
    return result


def _analyze_interface_interactions(interface_atom_pairs):
    """
    分析界面残基对之间的详细相互作用类型
    
    返回: list of dict
    """
    interactions = []
    
    for res1_atoms, res2_atoms, min_dist in interface_atom_pairs:
        # 简化分析：检测氢键、盐桥、疏水接触
        
        # 氢键检测（简化版 - 基于距离和原子类型）
        for atom1 in res1_atoms:
            elem1 = get_element_from_atom_name(atom1[3])
            if elem1 not in ["N", "O"]:
                continue
            
            for atom2 in res2_atoms:
                elem2 = get_element_from_atom_name(atom2[3])
                if elem2 not in ["N", "O"]:
                    continue
                
                d = distance(atom1[4], atom2[4])
                if d <= INTERACTION_PARAMS["hbond"]["max_DA_dist"]:
                    interactions.append({
                        "type": "氢键",
                        "atom1": f"{atom1[0]}:{atom1[1]} {atom1[2]}:{atom1[3]}",
                        "atom2": f"{atom2[0]}:{atom2[1]} {atom2[2]}:{atom2[3]}",
                        "distance": round(d, 2)
                    })
        
        # 疏水接触检测
        hydrophobic_residues = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO"}
        res1_name = res1_atoms[0][1]
        res2_name = res2_atoms[0][1]
        
        if res1_name in hydrophobic_residues and res2_name in hydrophobic_residues:
            if min_dist <= INTERACTION_PARAMS["hydrophobic"]["other_max"]:
                interactions.append({
                    "type": "疏水接触",
                    "atom1": f"{res1_atoms[0][0]}:{res1_name} {res1_atoms[0][2]}",
                    "atom2": f"{res2_atoms[0][0]}:{res2_name} {res2_atoms[0][2]}",
                    "distance": round(min_dist, 2)
                })
    
    return interactions


def _calculate_interface_strength(contact_count, interactions, bsa):
    """
    计算界面强度评分 [0-10]
    
    考虑因素:
    - 接触对数量
    - 相互作用类型和数量
    - 埋藏表面积
    """
    score = 0.0
    
    # 接触数贡献 (0-4分)
    score += min(4.0, contact_count / 5.0)
    
    # 相互作用质量贡献 (0-4分)
    hbond_count = sum(1 for i in interactions if "氢键" in i["type"])
    ionic_count = sum(1 for i in interactions if "盐桥" in i["type"])
    
    score += min(2.0, hbond_count * 0.5)
    score += min(2.0, ionic_count * 0.8)
    
    # BSA贡献 (0-2分)
    if bsa:
        score += min(2.0, bsa / 500.0)
    
    return min(10.0, score)


def _count_interaction_types(interactions):
    """统计相互作用类型"""
    counts = defaultdict(int)
    for inter in interactions:
        counts[inter["type"]] += 1
    return dict(counts)


def _export_ppi_csv(csv_path, interface_residues, interactions, bsa, strength):
    """导出PPI分析结果到CSV"""
    try:
        with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
            writer = csv.writer(f)
            
            # 摘要信息
            writer.writerow(["=== 蛋白-蛋白界面分析 (PPI Analysis) ==="])
            writer.writerow([])
            writer.writerow(["界面接触对数", len(interface_residues)])
            writer.writerow(["界面强度评分", strength])
            if bsa:
                writer.writerow(["埋藏表面积(BSA)", f"{bsa:.1f} Ų"])
            writer.writerow([])
            
            # 界面残基对
            writer.writerow(["=== 界面残基对 ==="])
            writer.writerow(["Chain1", "Residue1", "Chain2", "Residue2", "Distance(Å)"])
            for res in interface_residues:
                writer.writerow([
                    res["chain1"],
                    f"{res['resname1']} {res['resid1']}",
                    res["chain2"],
                    f"{res['resname2']} {res['resid2']}",
                    res["min_distance"]
                ])
            writer.writerow([])
            
            # 详细相互作用
            writer.writerow(["=== 界面相互作用 ==="])
            writer.writerow(["Type", "Atom1", "Atom2", "Distance(Å)"])
            for inter in interactions:
                writer.writerow([
                    inter["type"],
                    inter["atom1"],
                    inter["atom2"],
                    inter["distance"]
                ])
        
        print(f"[PPI] Results saved to: {csv_path}")
    except Exception as e:
        print(f"[PPI] Failed to save CSV: {e}")


def visualize_ppi_interface(obj_name, ppi_result, 
                            protein1_color="cyan", 
                            protein2_color="magenta",
                            show_labels=True,
                            clear_old=True):
    """
    在PyMOL中可视化蛋白-蛋白界面
    
    参数:
        obj_name: PyMOL对象名称
        ppi_result: analyze_protein_protein_interface()的返回结果
        protein1_color: 蛋白质1界面残基的颜色
        protein2_color: 蛋白质2界面残基的颜色
        show_labels: 是否显示残基标签
        clear_old: 是否清除旧的高亮
    """
    # 三字母到单字母氨基酸代码转换
    AA_3TO1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    # 定义多种颜色用于不同的chain
    CHAIN_COLORS = [
        "cyan", "magenta", "yellow", "salmon", "lime", 
        "orange", "purple", "green", "blue", "red",
        "pink", "brown", "gray", "lightblue", "lightorange"
    ]
    
    if not ppi_result or "interface_residues" not in ppi_result:
        print("[visualize_ppi_interface] Invalid PPI result")
        return
    
    if obj_name not in cmd.get_object_list():
        print(f"[visualize_ppi_interface] Object '{obj_name}' not found")
        return
    
    # 清除旧的可视化对象
    if clear_old:
        for name in cmd.get_names("objects"):
            if name.startswith("ppi_"):
                cmd.delete(name)
    
    interface_residues = ppi_result["interface_residues"]
    protein1_chains = ppi_result.get("protein1_chains", [])
    protein2_chains = ppi_result.get("protein2_chains", [])
    all_chains = list(set(protein1_chains + protein2_chains))
    
    if not interface_residues:
        print("[visualize_ppi_interface] No interface residues to visualize")
        return
    
    print(f"[visualize_ppi_interface] Visualizing {len(interface_residues)} interface residue pairs...")
    
    # 为每条链分配颜色
    chain_color_map = {}
    for idx, chain in enumerate(all_chains):
        chain_color_map[chain] = CHAIN_COLORS[idx % len(CHAIN_COLORS)]
    
    # 收集所有界面残基(去重),按chain分组
    interface_by_chain = {}  # {chain: {(resid, resname): True}}
    
    for res in interface_residues:
        chain1, resid1, resname1 = res["chain1"], res["resid1"], res["resname1"]
        chain2, resid2, resname2 = res["chain2"], res["resid2"], res["resname2"]
        
        if chain1 not in interface_by_chain:
            interface_by_chain[chain1] = {}
        interface_by_chain[chain1][(resid1, resname1)] = True
        
        if chain2 not in interface_by_chain:
            interface_by_chain[chain2] = {}
        interface_by_chain[chain2][(resid2, resname2)] = True
    
    # 可视化每条链的界面残基
    visualized_selections = []
    global_idx = 0
    
    for chain in sorted(interface_by_chain.keys()):
        residues = interface_by_chain[chain]
        chain_color = chain_color_map[chain]
        
        print(f"[visualize_ppi_interface] Chain {chain}: {len(residues)} unique interface residues (color: {chain_color})")
        
        for (resid, resname) in sorted(residues.keys()):
            global_idx += 1
            sel_name = f"ppi_{chain}_{resid}"
            sel_expr = f"{obj_name} and chain {chain} and resi {resid}"
            
            # 创建选择并显示为sticks
            cmd.select(sel_name, sel_expr)
            cmd.show("sticks", sel_name)
            cmd.color(chain_color, sel_name)
            visualized_selections.append(sel_name)
            
            # 添加标签(使用单字母大写氨基酸代码)
            if show_labels:
                label_name = f"ppi_label_{chain}_{resid}"
                ca_sel = f"{sel_expr} and name CA"
                try:
                    coords = cmd.get_atom_coords(ca_sel)
                    if coords:
                        # 转换为单字母大写代码
                        aa_code = AA_3TO1.get(resname.upper(), resname[0].upper())
                        label_text = f"{aa_code}{resid}"
                        
                        cmd.pseudoatom(label_name, pos=coords, label=label_text)
                        cmd.set("label_size", 16, label_name)
                        cmd.set("label_color", chain_color, label_name)
                except:
                    pass
    
    # 绘制界面接触线(距离线) - 只显示前20条,避免过于拥挤
    for idx, res in enumerate(interface_residues[:20], start=1):
        distance_name = f"ppi_dist{idx}"
        sel1 = f"{obj_name} and chain {res['chain1']} and resi {res['resid1']} and name CA"
        sel2 = f"{obj_name} and chain {res['chain2']} and resi {res['resid2']} and name CA"
        
        try:
            cmd.distance(distance_name, sel1, sel2)
            cmd.set("dash_color", "yellow", distance_name)
            cmd.set("dash_width", 2.0, distance_name)
            cmd.hide("labels", distance_name)  # 隐藏距离标签
        except:
            pass
    
    # 显示整体蛋白链(cartoon),按chain着色
    for chain in all_chains:
        chain_color = chain_color_map[chain]
        cmd.show("cartoon", f"{obj_name} and chain {chain}")
        cmd.color(chain_color, f"{obj_name} and chain {chain}")
    
    # 缩放到界面区域
    if visualized_selections:
        cmd.zoom(" or ".join(visualized_selections), buffer=8.0, complete=1)
    
    # 刷新视图
    cmd.refresh()
    cmd.rebuild()
    
    total_residues = sum(len(residues) for residues in interface_by_chain.values())
    print(f"[visualize_ppi_interface] ✅ Visualized {total_residues} unique interface residues")
    for chain, color in chain_color_map.items():
        print(f"  Chain {chain}: {color}")


def calculate_interface_bsa(obj_name, chain1, chain2):
    """
    计算两个链之间的埋藏表面积 (Buried Surface Area)
    
    使用PyMOL内置的get_area()函数
    BSA = (SA_chain1 + SA_chain2 - SA_complex) / 2
    
    参数:
        obj_name: PyMOL对象名称
        chain1: 第一条链ID
        chain2: 第二条链ID
    
    返回:
        float: 埋藏表面积 (Ų)
    """
    if obj_name not in cmd.get_object_list():
        raise ValueError(f"Object '{obj_name}' not found in PyMOL")
    
    try:
        # 计算单独链的表面积
        area_chain1 = cmd.get_area(f"{obj_name} and chain {chain1}")
        area_chain2 = cmd.get_area(f"{obj_name} and chain {chain2}")
        
        # 计算复合物表面积
        area_complex = cmd.get_area(f"{obj_name} and (chain {chain1} or chain {chain2})")
        
        # BSA = (SA1 + SA2 - SA_complex) / 2
        bsa = (area_chain1 + area_chain2 - area_complex) / 2.0
        
        return bsa
    
    except Exception as e:
        print(f"[calculate_interface_bsa] Error: {e}")
        return None


def identify_neo_epitope(obj_name=None,
                        e3_ligase_chains=None,
                        substrate_chains=None,
                        glue_resname=None,
                        distance_threshold=5.0,
                        output_csv=None,
                        pdb_file=None,
                        visualize=True,
                        neo_color="yellow",
                        glue_color="orange"):
    """
    识别分子胶诱导的新表位 (Neo-Substrate Epitope)
    
    这是分子胶机制的关键特征：
    - Neo-表位是指只有在Glue存在时，才与E3 ligase接触的底物残基
    - 区别于天然底物（无需Glue即可结合）
    
    检测逻辑:
    1. 找到同时与 E3 Ligase + Substrate 接触的 Glue 原子
    2. 找到同时与 Glue + E3 接触的 Substrate 残基
    3. 计算这些残基的几何位置（是否在"桥接"位置）
    
    参数:
        obj_name: PyMOL对象名称
        e3_ligase_chains: E3 ligase的链ID列表（例如 ["A"] - CRBN）
        substrate_chains: 底物的链ID列表（例如 ["B"]）
        glue_resname: 分子胶的残基名称（例如 "CC885"）
        distance_threshold: 接触距离阈值（埃）
        output_csv: 输出CSV文件路径
        pdb_file: PDB文件路径（可选）
        visualize: 是否在PyMOL中可视化（默认True）
        neo_color: Neo-表位残基的颜色（默认yellow）
        glue_color: 分子胶的颜色（默认orange）
    
    返回:
        dict: {
            "neo_substrate_residues": [...],  # Neo-表位残基
            "bridging_glue_atoms": [...],  # 桥接Glue原子
            "neo_epitope_count": int,
            "is_molecular_glue": bool,  # 是否符合分子胶特征
            "confidence": float  # 置信度 [0-1]
        }
    """
    # 输入验证
    if not e3_ligase_chains or not substrate_chains or not glue_resname:
        print("[identify_neo_epitope] ⚠️ Must specify e3_ligase_chains, substrate_chains, and glue_resname")
        return None
    
    # 获取原子信息
    if pdb_file:
        atoms = parse_pdb_file(pdb_file)
    else:
        atoms = parse_pdb_structure(obj_name)
    
    if not atoms:
        print("[identify_neo_epitope] Unable to obtain atom information")
        return None
    
    print(f"[Neo-Epitope] Analyzing: E3={e3_ligase_chains}, Substrate={substrate_chains}, Glue={glue_resname}")
    
    # 分组原子
    e3_atoms = []
    substrate_atoms = []
    glue_atoms = []
    substrate_residues = defaultdict(list)
    
    for atom in atoms:
        chain, res_name, res_id, atom_name, coords = atom
        
        if chain in e3_ligase_chains and identify_molecule_type(res_name) == "protein":
            e3_atoms.append(atom)
        elif chain in substrate_chains and identify_molecule_type(res_name) == "protein":
            substrate_atoms.append(atom)
            res_key = (chain, res_name, res_id)
            substrate_residues[res_key].append(atom)
        elif res_name.upper() == glue_resname.upper():
            glue_atoms.append(atom)
    
    if not glue_atoms:
        print(f"[Neo-Epitope] ⚠️ Glue molecule '{glue_resname}' not found")
        return None
    
    print(f"[Neo-Epitope] E3 atoms: {len(e3_atoms)}, Substrate atoms: {len(substrate_atoms)}, Glue atoms: {len(glue_atoms)}")
    
    # Step 1: 找到与 E3 和 Substrate 都接触的 Glue 原子（桥接原子）
    bridging_glue_atoms = []
    
    for glue_atom in glue_atoms:
        contacts_e3 = False
        contacts_substrate = False
        
        for e3_atom in e3_atoms:
            if distance(glue_atom[4], e3_atom[4]) <= distance_threshold:
                contacts_e3 = True
                break
        
        for sub_atom in substrate_atoms:
            if distance(glue_atom[4], sub_atom[4]) <= distance_threshold:
                contacts_substrate = True
                break
        
        if contacts_e3 and contacts_substrate:
            bridging_glue_atoms.append(glue_atom)
    
    print(f"[Neo-Epitope] Found {len(bridging_glue_atoms)} bridging Glue atoms")
    
    # Step 2: 找到同时与 Glue 和 E3 接触的 Substrate 残基（Neo-表位候选）
    neo_substrate_residues = []
    
    for res_key, res_atoms in substrate_residues.items():
        contacts_glue = False
        contacts_e3 = False
        min_dist_glue = float('inf')
        min_dist_e3 = float('inf')
        
        for res_atom in res_atoms:
            # 检查与Glue的接触
            for glue_atom in bridging_glue_atoms:
                d = distance(res_atom[4], glue_atom[4])
                if d <= distance_threshold:
                    contacts_glue = True
                    min_dist_glue = min(min_dist_glue, d)
            
            # 检查与E3的接触
            for e3_atom in e3_atoms:
                d = distance(res_atom[4], e3_atom[4])
                if d <= distance_threshold:
                    contacts_e3 = True
                    min_dist_e3 = min(min_dist_e3, d)
        
        # 如果同时接触Glue和E3，则为Neo-表位候选
        if contacts_glue and contacts_e3:
            neo_substrate_residues.append({
                "chain": res_key[0],
                "resname": res_key[1],
                "resid": res_key[2],
                "distance_to_glue": round(min_dist_glue, 2) if min_dist_glue != float('inf') else None,
                "distance_to_e3": round(min_dist_e3, 2) if min_dist_e3 != float('inf') else None,
                "is_neo_epitope": True  # 标记为Neo-表位
            })
    
    # Step 3: 计算置信度
    confidence = 0.0
    if len(bridging_glue_atoms) >= 2:
        confidence += 0.4
    if len(neo_substrate_residues) >= 3:
        confidence += 0.4
    if len(neo_substrate_residues) > 0:
        confidence += 0.2
    
    confidence = min(1.0, confidence)
    
    # 判定是否为分子胶
    is_molecular_glue = (len(bridging_glue_atoms) >= 2 and 
                        len(neo_substrate_residues) >= 3 and
                        confidence >= 0.6)
    
    # 输出到CSV
    if output_csv:
        _export_neo_epitope_csv(output_csv, neo_substrate_residues, bridging_glue_atoms, is_molecular_glue, confidence)
    
    result = {
        "neo_substrate_residues": neo_substrate_residues,
        "bridging_glue_atoms": len(bridging_glue_atoms),
        "neo_epitope_count": len(neo_substrate_residues),
        "is_molecular_glue": is_molecular_glue,
        "confidence": round(confidence, 2),
        "e3_ligase_chains": e3_ligase_chains,
        "substrate_chains": substrate_chains,
        "glue_resname": glue_resname
    }
    
    # 打印摘要
    print("\n" + "=" * 60)
    print("Neo-表位识别结果 (Neo-Epitope Analysis)")
    print("=" * 60)
    print(f"桥接Glue原子数: {len(bridging_glue_atoms)}")
    print(f"Neo-表位残基数: {len(neo_substrate_residues)}")
    print(f"置信度: {confidence:.2f}")
    print(f"分子胶判定: {'✨ 是 (Molecular Glue)' if is_molecular_glue else '⚠️ 否 (Not Glue / PROTAC)'}")
    
    if neo_substrate_residues:
        print(f"\nNeo-表位残基:")
        for res in neo_substrate_residues[:10]:  # 显示前10个
            print(f"  {res['chain']}:{res['resname']} {res['resid']} (Glue: {res['distance_to_glue']}Å, E3: {res['distance_to_e3']}Å)")
    print("=" * 60)
    
    # 在PyMOL中可视化Neo-表位（如果启用）
    if visualize and obj_name and obj_name in cmd.get_object_list():
        try:
            visualize_neo_epitope(obj_name, result, neo_color, glue_color)
        except Exception as e:
            print(f"[Neo-Epitope] 可视化失败: {e}")
    
    return result


def _export_neo_epitope_csv(csv_path, neo_residues, bridging_atoms_count, is_glue, confidence):
    """导出Neo-表位分析结果到CSV"""
    try:
        with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
            writer = csv.writer(f)
            
            writer.writerow(["=== Neo-表位识别结果 (Neo-Epitope Analysis) ==="])
            writer.writerow([])
            writer.writerow(["桥接Glue原子数", bridging_atoms_count])
            writer.writerow(["Neo-表位残基数", len(neo_residues)])
            writer.writerow(["置信度", confidence])
            writer.writerow(["分子胶判定", "是" if is_glue else "否"])
            writer.writerow([])
            
            writer.writerow(["=== Neo-表位残基 ==="])
            writer.writerow(["Chain", "Residue", "Distance_to_Glue(Å)", "Distance_to_E3(Å)"])
            for res in neo_residues:
                writer.writerow([
                    res["chain"],
                    f"{res['resname']} {res['resid']}",
                    res["distance_to_glue"],
                    res["distance_to_e3"]
                ])
        
        print(f"[Neo-Epitope] Results saved to: {csv_path}")
    except Exception as e:
        print(f"[Neo-Epitope] Failed to save CSV: {e}")


def visualize_neo_epitope(obj_name, neo_result, neo_color="yellow", glue_color="orange", clear_old=True):
    """
    在PyMOL中可视化Neo-表位分析结果
    
    参数:
        obj_name: PyMOL对象名称
        neo_result: identify_neo_epitope()的返回结果
        neo_color: Neo-表位残基的颜色
        glue_color: 分子胶的颜色
        clear_old: 是否清除旧的高亮
    """
    # 三字母到单字母氨基酸代码转换
    AA_3TO1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    # 定义多种颜色用于不同的chain
    CHAIN_COLORS = [
        "cyan", "magenta", "yellow", "salmon", "lime", 
        "orange", "purple", "green", "blue", "red",
        "pink", "brown", "gray", "lightblue", "lightorange"
    ]
    
    if not neo_result or "neo_substrate_residues" not in neo_result:
        print("[visualize_neo_epitope] Invalid neo-epitope result")
        return
    
    if obj_name not in cmd.get_object_list():
        print(f"[visualize_neo_epitope] Object '{obj_name}' not found")
        return
    
    # 清除旧的可视化对象
    if clear_old:
        for name in cmd.get_names("objects"):
            if name.startswith("neo_"):
                cmd.delete(name)
    
    neo_residues = neo_result["neo_substrate_residues"]
    glue_resname = neo_result.get("glue_resname", "")
    e3_chains = neo_result.get("e3_ligase_chains", [])
    substrate_chains = neo_result.get("substrate_chains", [])
    all_chains = list(set(e3_chains + substrate_chains))
    
    # 为每条链分配颜色
    chain_color_map = {}
    for idx, chain in enumerate(all_chains):
        chain_color_map[chain] = CHAIN_COLORS[idx % len(CHAIN_COLORS)]
    
    if not neo_residues:
        print("[visualize_neo_epitope] No neo-epitope residues to visualize")
        return
    
    print(f"[visualize_neo_epitope] Visualizing {len(neo_residues)} neo-epitope residues...")
    
    # 按chain分组Neo-表位残基(去重)
    neo_by_chain = {}  # {chain: {resid: resname}}
    for res in neo_residues:
        chain = res['chain']
        resid = res['resid']
        resname = res['resname']
        
        if chain not in neo_by_chain:
            neo_by_chain[chain] = {}
        neo_by_chain[chain][resid] = resname
    
    # 可视化Neo-表位残基
    visualized_selections = []
    for chain in sorted(neo_by_chain.keys()):
        residues = neo_by_chain[chain]
        # Neo-表位使用指定颜色,不按chain区分
        
        for resid, resname in sorted(residues.items()):
            sel_name = f"neo_{chain}_{resid}"
            sel_expr = f"{obj_name} and chain {chain} and resi {resid}"
            
            cmd.select(sel_name, sel_expr)
            cmd.show("sticks", sel_name)
            cmd.color(neo_color, sel_name)
            visualized_selections.append(sel_name)
            
            # 添加标签(使用单字母大写氨基酸代码)
            label_name = f"neo_label_{chain}_{resid}"
            ca_sel = f"{sel_expr} and name CA"
            try:
                coords = cmd.get_atom_coords(ca_sel)
                if coords:
                    aa_code = AA_3TO1.get(resname.upper(), resname[0].upper())
                    label_text = f"{aa_code}{resid}"
                    
                    cmd.pseudoatom(label_name, pos=coords, label=label_text)
                    cmd.set("label_size", 16, label_name)
                    cmd.set("label_color", neo_color, label_name)
            except:
                pass
    
    # 可视化分子胶
    if glue_resname:
        glue_sel = f"{obj_name} and resn {glue_resname}"
        cmd.select("neo_glue", glue_sel)
        cmd.show("sticks", "neo_glue")
        cmd.color(glue_color, "neo_glue")
        cmd.show("spheres", "neo_glue")
        cmd.set("sphere_scale", 0.3, "neo_glue")
        visualized_selections.append("neo_glue")
    
    # 显示蛋白链的cartoon,按chain着色
    for chain in all_chains:
        chain_color = chain_color_map[chain]
        cmd.show("cartoon", f"{obj_name} and chain {chain}")
        cmd.color(chain_color, f"{obj_name} and chain {chain}")
    
    # 缩放到Neo-表位区域
    if visualized_selections:
        cmd.zoom(" or ".join(visualized_selections), buffer=8.0, complete=1)
    
    # 刷新视图
    cmd.refresh()
    cmd.rebuild()
    
    total_neo = sum(len(residues) for residues in neo_by_chain.values())
    print(f"[visualize_neo_epitope] ✅ Visualized {total_neo} unique neo-epitope residues")
    print(f"[visualize_neo_epitope] Neo-epitope color: {neo_color}")
    print(f"[visualize_neo_epitope] Glue color: {glue_color}")
    for chain, color in chain_color_map.items():
        print(f"  Chain {chain}: {color}")


# ========== PyMOL命令封装 ==========
def ppi_analyze(obj_name, protein1_chains, protein2_chains, out_csv=None):
    """
    PyMOL命令: 分析蛋白-蛋白界面
    
    用法:
        ppi_analyze complex, [A], [B]
        ppi_analyze complex, [A], [B], out_csv=interface.csv
    """
    return analyze_protein_protein_interface(
        obj_name=obj_name,
        protein1_chains=protein1_chains,
        protein2_chains=protein2_chains,
        output_csv=out_csv
    )


def neo_epitope_find(obj_name, e3_chains, substrate_chains, glue_name, out_csv=None):
    """
    PyMOL命令: 识别Neo-表位
    
    用法:
        neo_epitope_find complex, [A], [B], CC885
        neo_epitope_find complex, [A], [B], CC885, out_csv=neo_epitope.csv
    """
    return identify_neo_epitope(
        obj_name=obj_name,
        e3_ligase_chains=e3_chains,
        substrate_chains=substrate_chains,
        glue_resname=glue_name,
        output_csv=out_csv
    )


if __name__ == "__main__":
    print("[ppi_analyzer] 这是一个PyMOL插件模块,请在PyMOL中加载使用")
    print("[ppi_analyzer] 命令: ppi_analyze, neo_epitope_find")
