# -*- coding: utf-8 -*-
"""
分子胶设计分析模块 (Molecular Glue Design Analyzer)

功能：
1. Ternary complex 建模：G-loop 替换与对齐
2. 碰撞检测与界面互补性分析
3. Exit vector 识别（PROTAC linker 连接点）
4. 静电环境与化学修饰建议

基于已有功能：
- G-motif 验证（validate_g_motif_geometry, validate_crbn_hbonds）
- 口袋检测（pocket_detector, pocket_visualizer）
- 静电分析（APBS 或 Quick 模式）
- PPI 界面分析（ppi_analyzer）

Author: Vesper
Date: 2025-11-11
"""

from __future__ import print_function
import os, csv, tempfile
from collections import defaultdict
from pymol import cmd
import numpy as np


# ====== 1) G-loop 对齐与替换 ======
def align_gloop_for_modeling(template_obj, template_chain, template_gloop_range,
                               target_obj, target_chain, target_gloop_range,
                               output_obj=None):
    """
    将目标蛋白的 G-loop 对齐到模板 ternary complex 的 G-loop
    
    用途：建模新的 CRBN-ligand-neosubstrate 复合物
    
    参数:
        template_obj: 模板结构（如 CRBN-IMiD-GSPT1 的 6H0G）
        template_chain: 模板 G-loop 链
        template_gloop_range: 模板 G-loop 残基范围 "60-67"
        target_obj: 目标蛋白（新底物）
        target_chain: 目标 G-loop 链
        target_gloop_range: 目标 G-loop 残基范围
        output_obj: 输出对象名（默认 target_obj + "_aligned"）
    
    返回:
        dict: {
            'rmsd': float,
            'aligned_obj': str,
            'template_gloop_sel': str,
            'target_gloop_sel': str
        }
    """
    if output_obj is None:
        output_obj = f"{target_obj}_aligned"
    
    # 构建选择
    template_sel = f"{template_obj} and chain {template_chain} and resi {template_gloop_range} and name CA"
    target_sel = f"{target_obj} and chain {target_chain} and resi {target_gloop_range} and name CA"
    
    # 检查原子数
    template_count = cmd.count_atoms(template_sel)
    target_count = cmd.count_atoms(target_sel)
    
    if template_count != target_count or template_count != 8:
        print(f"⚠️ 警告: G-loop 长度不匹配")
        print(f"   模板: {template_count} Cα, 目标: {target_count} Cα (预期 8)")
    
    # 复制目标结构
    cmd.create(output_obj, target_obj)
    
    # 对齐（仅使用 G-loop backbone）
    result = cmd.align(
        f"{output_obj} and chain {target_chain} and resi {target_gloop_range}",
        template_sel,
        cycles=0  # 不进行迭代优化，保持初始对齐
    )
    
    rmsd = result[0]
    
    print(f"\n{'='*70}")
    print(f"G-loop 对齐 (G-loop Alignment)")
    print(f"{'='*70}")
    print(f"模板: {template_obj} {template_chain}:{template_gloop_range}")
    print(f"目标: {target_obj} {target_chain}:{target_gloop_range}")
    print(f"RMSD: {rmsd:.3f} Å")
    print(f"输出: {output_obj}")
    print(f"{'='*70}\n")
    
    return {
        'rmsd': rmsd,
        'aligned_obj': output_obj,
        'template_gloop_sel': template_sel,
        'target_gloop_sel': f"{output_obj} and chain {target_chain} and resi {target_gloop_range}"
    }


# ====== 2) 碰撞检测 ======
def detect_clashes_at_interface(aligned_obj, aligned_chain, aligned_gloop_range,
                                  crbn_obj, crbn_chain,
                                  clash_threshold=2.0, contact_threshold=4.5):
    """
    检测新底物 G-loop 与 CRBN 的碰撞（clashes）和接触
    
    参数:
        aligned_obj: 对齐后的新底物
        aligned_chain: G-loop 链
        aligned_gloop_range: G-loop 范围
        crbn_obj: CRBN 对象（通常与模板同一对象）
        crbn_chain: CRBN 链
        clash_threshold: 碰撞阈值（Å）默认 2.0
        contact_threshold: 接触阈值（Å）默认 4.5
    
    返回:
        dict: {
            'clashes': [{'gloop_res': str, 'crbn_res': str, 'distance': float}, ...],
            'contacts': [...],
            'has_severe_clashes': bool,
            'contact_count': int
        }
    """
    # 获取原子
    gloop_sel = f"{aligned_obj} and chain {aligned_chain} and resi {aligned_gloop_range}"
    crbn_sel = f"{crbn_obj} and chain {crbn_chain} and polymer.protein"
    
    gloop_atoms = []
    for atom in cmd.get_model(gloop_sel).atom:
        gloop_atoms.append({
            'coord': np.array([atom.coord[0], atom.coord[1], atom.coord[2]]),
            'resn': atom.resn.strip(),
            'resi': atom.resi.strip(),
            'name': atom.name.strip(),
            'elem': atom.elem.strip()
        })
    
    crbn_atoms = []
    for atom in cmd.get_model(crbn_sel).atom:
        crbn_atoms.append({
            'coord': np.array([atom.coord[0], atom.coord[1], atom.coord[2]]),
            'resn': atom.resn.strip(),
            'resi': atom.resi.strip(),
            'name': atom.name.strip(),
            'elem': atom.elem.strip()
        })
    
    # 检测碰撞和接触
    clashes = []
    contacts = []
    
    for g_atom in gloop_atoms:
        for c_atom in crbn_atoms:
            dist = np.linalg.norm(g_atom['coord'] - c_atom['coord'])
            
            # 排除氢键 backbone atoms 的正常接触
            is_backbone_pair = (
                g_atom['name'] in ('N', 'CA', 'C', 'O') and
                c_atom['name'] in ('N', 'CA', 'C', 'O')
            )
            
            if dist < clash_threshold and not is_backbone_pair:
                clashes.append({
                    'gloop_res': f"{g_atom['resn']} {g_atom['resi']}",
                    'gloop_atom': g_atom['name'],
                    'crbn_res': f"{c_atom['resn']} {c_atom['resi']}",
                    'crbn_atom': c_atom['name'],
                    'distance': round(dist, 2)
                })
            elif clash_threshold <= dist <= contact_threshold:
                contacts.append({
                    'gloop_res': f"{g_atom['resn']} {g_atom['resi']}",
                    'gloop_atom': g_atom['name'],
                    'crbn_res': f"{c_atom['resn']} {c_atom['resi']}",
                    'crbn_atom': c_atom['name'],
                    'distance': round(dist, 2)
                })
    
    has_severe_clashes = any(c['distance'] < 1.5 for c in clashes)
    
    # 打印结果
    print(f"\n{'='*70}")
    print(f"碰撞检测 (Clash Detection)")
    print(f"{'='*70}")
    print(f"碰撞数: {len(clashes)} (< {clash_threshold} Å)")
    if clashes:
        print(f"\n⚠️ 碰撞列表:")
        for clash in clashes[:10]:  # 前10个
            print(f"  {clash['gloop_res']:15s} ({clash['gloop_atom']:4s}) ↔ "
                  f"{clash['crbn_res']:15s} ({clash['crbn_atom']:4s})  {clash['distance']} Å")
        if len(clashes) > 10:
            print(f"  ... 还有 {len(clashes)-10} 个碰撞")
    
    print(f"\n接触数: {len(contacts)} ({clash_threshold}-{contact_threshold} Å)")
    
    if has_severe_clashes:
        print(f"\n❌ 严重碰撞 (< 1.5 Å): 结构不可行")
    elif len(clashes) > 5:
        print(f"\n⚠️ 多处碰撞: 需要优化侧链")
    else:
        print(f"\n✅ 碰撞可接受")
    
    print(f"{'='*70}\n")
    
    return {
        'clashes': clashes,
        'contacts': contacts,
        'has_severe_clashes': has_severe_clashes,
        'contact_count': len(contacts),
        'clash_count': len(clashes)
    }


# ====== 3) Exit Vector 识别 ======
def identify_exit_vectors(ligand_obj, ligand_resname, 
                           crbn_obj, crbn_chain,
                           gloop_obj, gloop_chain, gloop_range,
                           solvent_radius=5.0):
    """
    识别分子胶上的 exit vectors（可用于 PROTAC linker 连接的位点）
    
    策略:
    1. 找到配体上远离 CRBN 和 G-loop 核心界面的原子
    2. 检查这些原子周围的溶剂暴露度
    3. 推荐最适合连接 linker 的位置
    
    参数:
        ligand_obj: 配体对象名
        ligand_resname: 配体残基名（如 CC9, LEN）
        crbn_obj: CRBN 对象
        crbn_chain: CRBN 链
        gloop_obj: G-loop 对象（对齐后的）
        gloop_chain: G-loop 链
        gloop_range: G-loop 范围
        solvent_radius: 溶剂暴露检测半径（Å）
    
    返回:
        list: [{'atom_name': str, 'score': float, 'reason': str}, ...]
    """
    # 获取配体原子
    lig_sel = f"{ligand_obj} and resn {ligand_resname}"
    lig_model = cmd.get_model(lig_sel)
    
    if len(lig_model.atom) == 0:
        print(f"⚠️ 未找到配体 {ligand_resname}")
        return []
    
    # 获取 CRBN 和 G-loop 原子（用于计算距离）
    crbn_coords = []
    for atom in cmd.get_model(f"{crbn_obj} and chain {crbn_chain} and polymer.protein").atom:
        crbn_coords.append(np.array([atom.coord[0], atom.coord[1], atom.coord[2]]))
    
    gloop_coords = []
    for atom in cmd.get_model(f"{gloop_obj} and chain {gloop_chain} and resi {gloop_range}").atom:
        gloop_coords.append(np.array([atom.coord[0], atom.coord[1], atom.coord[2]]))
    
    # 评估每个配体原子
    candidates = []
    
    for atom in lig_model.atom:
        coord = np.array([atom.coord[0], atom.coord[1], atom.coord[2]])
        aname = atom.name.strip()
        
        # 跳过核心结合原子（如 phthalimide/glutarimide 的羰基 O）
        if aname in ('O1', 'O2', 'O3', 'O4', 'N1'):  # 常见 IMiD 核心原子
            continue
        
        # 计算到 CRBN 和 G-loop 的最小距离
        min_dist_crbn = min([np.linalg.norm(coord - c) for c in crbn_coords])
        min_dist_gloop = min([np.linalg.norm(coord - g) for g in gloop_coords])
        
        # 评分标准：
        # 1. 远离核心界面 (+)
        # 2. 溶剂暴露 (+)
        # 3. 是否为末端/支链原子 (+)
        
        score = 0.0
        reasons = []
        
        # 距离评分
        if min_dist_crbn > 5.0:
            score += 3.0
            reasons.append("远离CRBN")
        elif min_dist_crbn > 4.0:
            score += 1.5
        
        if min_dist_gloop > 5.0:
            score += 3.0
            reasons.append("远离G-loop")
        elif min_dist_gloop > 4.0:
            score += 1.5
        
        # 溶剂暴露（简化检测：检查半径内蛋白原子数）
        nearby_protein = cmd.count_atoms(
            f"({crbn_obj} or {gloop_obj}) and polymer.protein within {solvent_radius} of "
            f"{lig_sel} and name {aname}"
        )
        if nearby_protein < 5:
            score += 2.0
            reasons.append("溶剂暴露")
        
        # 元素类型（C > N > O，避免破坏极性相互作用）
        elem = atom.elem.strip().upper()
        if elem == 'C':
            score += 1.0
            reasons.append("碳原子")
        
        if score > 3.0:
            candidates.append({
                'atom_name': aname,
                'element': elem,
                'score': round(score, 1),
                'min_dist_crbn': round(min_dist_crbn, 1),
                'min_dist_gloop': round(min_dist_gloop, 1),
                'reasons': ', '.join(reasons)
            })
    
    # 排序
    candidates.sort(key=lambda x: x['score'], reverse=True)
    
    # 打印结果
    print(f"\n{'='*70}")
    print(f"Exit Vector 识别 (Exit Vector Identification)")
    print(f"{'='*70}")
    print(f"配体: {ligand_resname}")
    print(f"候选位点: {len(candidates)}\n")
    
    if candidates:
        print(f"{'原子':<10} {'元素':<6} {'评分':<8} {'→CRBN':<10} {'→G-loop':<10} {'原因'}")
        print(f"{'-'*70}")
        for cand in candidates[:5]:
            print(f"{cand['atom_name']:<10} {cand['element']:<6} {cand['score']:<8} "
                  f"{cand['min_dist_crbn']:<10} {cand['min_dist_gloop']:<10} {cand['reasons']}")
        
        print(f"\n推荐: {candidates[0]['atom_name']} (评分 {candidates[0]['score']})")
    else:
        print("⚠️ 未找到合适的 exit vector")
    
    print(f"{'='*70}\n")
    
    return candidates


# ====== 4) 静电环境分析 ======
def analyze_electrostatic_environment(obj_name, gloop_chain, gloop_range,
                                       crbn_chain, sampling_radius=10.0):
    """
    分析 G-loop 周围的静电环境
    
    用途: 指导分子胶上的化学修饰（带电 vs 疏水）
    
    参数:
        obj_name: 对象名
        gloop_chain: G-loop 链
        gloop_range: G-loop 范围
        crbn_chain: CRBN 链
        sampling_radius: 采样半径（Å）
    
    返回:
        dict: {
            'net_charge': int,
            'charged_residues': {'positive': [...], 'negative': [...]},
            'hydrophobic_ratio': float,
            'recommendation': str
        }
    """
    # 获取 G-loop 中心
    gloop_sel = f"{obj_name} and chain {gloop_chain} and resi {gloop_range} and name CA"
    gloop_coords = []
    for atom in cmd.get_model(gloop_sel).atom:
        gloop_coords.append([atom.coord[0], atom.coord[1], atom.coord[2]])
    
    if not gloop_coords:
        print("⚠️ 无法获取 G-loop 坐标")
        return None
    
    center = np.mean(gloop_coords, axis=0)
    
    # 获取周围残基
    nearby_sel = (
        f"({obj_name} and chain {gloop_chain} and resi {gloop_range}) or "
        f"({obj_name} and chain {crbn_chain} and polymer.protein within {sampling_radius} of "
        f"{obj_name} and chain {gloop_chain} and resi {gloop_range})"
    )
    
    residues = defaultdict(int)
    charged_pos = []
    charged_neg = []
    hydrophobic_count = 0
    total_count = 0
    
    for atom in cmd.get_model(nearby_sel).atom:
        resn = atom.resn.strip().upper()
        resi = atom.resi.strip()
        chain = atom.chain.strip()
        res_key = f"{chain}:{resn}{resi}"
        
        if res_key in residues:
            continue
        
        residues[res_key] = 1
        total_count += 1
        
        # 分类
        if resn in ('ARG', 'LYS', 'HIS'):
            charged_pos.append(res_key)
        elif resn in ('ASP', 'GLU'):
            charged_neg.append(res_key)
        elif resn in ('ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'MET', 'PRO'):
            hydrophobic_count += 1
    
    net_charge = len(charged_pos) - len(charged_neg)
    hydrophobic_ratio = hydrophobic_count / total_count if total_count > 0 else 0
    
    # 推荐
    if net_charge > 2:
        recommendation = "环境偏正电：考虑引入负电基团（羧酸、磺酸）"
    elif net_charge < -2:
        recommendation = "环境偏负电：考虑引入正电基团（胺、胍）"
    elif hydrophobic_ratio > 0.5:
        recommendation = "环境疏水：保持芳香/脂肪族取代，避免极性基团"
    else:
        recommendation = "环境中性/混合：平衡疏水和亲水取代"
    
    # 打印结果
    print(f"\n{'='*70}")
    print(f"静电环境分析 (Electrostatic Environment Analysis)")
    print(f"{'='*70}")
    print(f"采样半径: {sampling_radius} Å")
    print(f"总残基数: {total_count}")
    print(f"\n带正电残基 (+): {len(charged_pos)}")
    if charged_pos:
        print(f"  {', '.join(charged_pos[:5])}" + (f" ... (+{len(charged_pos)-5})" if len(charged_pos) > 5 else ""))
    
    print(f"\n带负电残基 (-): {len(charged_neg)}")
    if charged_neg:
        print(f"  {', '.join(charged_neg[:5])}" + (f" ... (+{len(charged_neg)-5})" if len(charged_neg) > 5 else ""))
    
    print(f"\n净电荷: {net_charge:+d}")
    print(f"疏水比例: {hydrophobic_ratio:.1%}")
    print(f"\n💡 设计建议: {recommendation}")
    print(f"{'='*70}\n")
    
    return {
        'net_charge': net_charge,
        'charged_residues': {'positive': charged_pos, 'negative': charged_neg},
        'hydrophobic_ratio': hydrophobic_ratio,
        'recommendation': recommendation
    }


# ====== 5) 综合分析工作流 ======
def comprehensive_glue_design_analysis(
    template_obj, template_crbn_chain, template_gloop_chain, template_gloop_range,
    target_obj, target_chain, target_gloop_range,
    ligand_resname,
    output_report=None
):
    """
    综合分析工作流：从 G-loop 对齐到分子胶优化建议
    
    步骤:
    1. G-loop 对齐
    2. 碰撞检测
    3. G-motif 几何和 CRBN 氢键验证
    4. Exit vector 识别
    5. 静电环境分析
    6. 生成设计报告
    
    参数:
        template_obj: 模板 ternary complex (如 6H0G)
        template_crbn_chain: 模板 CRBN 链
        template_gloop_chain: 模板 G-loop 链
        template_gloop_range: 模板 G-loop 范围
        target_obj: 目标蛋白
        target_chain: 目标 G-loop 链
        target_gloop_range: 目标 G-loop 范围
        ligand_resname: 配体名称
        output_report: 报告输出路径（默认自动生成）
    
    返回:
        dict: 完整分析结果
    """
    print("\n" + "🔬"*35)
    print("分子胶设计分析 - 综合工作流")
    print("Molecular Glue Design Analysis - Comprehensive Workflow")
    print("🔬"*35 + "\n")
    
    results = {}
    
    # Step 1: G-loop 对齐
    print("📍 Step 1: G-loop 对齐")
    align_result = align_gloop_for_modeling(
        template_obj, template_gloop_chain, template_gloop_range,
        target_obj, target_chain, target_gloop_range
    )
    results['alignment'] = align_result
    aligned_obj = align_result['aligned_obj']
    
    # Step 2: 碰撞检测
    print("📍 Step 2: 碰撞检测")
    clash_result = detect_clashes_at_interface(
        aligned_obj, target_chain, target_gloop_range,
        template_obj, template_crbn_chain
    )
    results['clashes'] = clash_result
    
    # Step 3: G-motif 验证
    print("📍 Step 3: G-motif 几何验证")
    try:
        from g_motif_analyzer import validate_g_motif_geometry, validate_crbn_hbonds
        
        geom_result = validate_g_motif_geometry(
            aligned_obj, target_chain, target_gloop_range
        )
        results['geometry'] = geom_result
        
        hb_result = validate_crbn_hbonds(
            template_obj, template_gloop_chain, template_gloop_range,
            template_crbn_chain
        )
        results['crbn_hbonds'] = hb_result
    except Exception as e:
        print(f"⚠️ G-motif 验证失败: {e}")
        results['geometry'] = None
        results['crbn_hbonds'] = None
    
    # Step 4: Exit vector 识别
    print("📍 Step 4: Exit Vector 识别")
    exit_vectors = identify_exit_vectors(
        template_obj, ligand_resname,
        template_obj, template_crbn_chain,
        aligned_obj, target_chain, target_gloop_range
    )
    results['exit_vectors'] = exit_vectors
    
    # Step 5: 静电环境
    print("📍 Step 5: 静电环境分析")
    electro_result = analyze_electrostatic_environment(
        aligned_obj, target_chain, target_gloop_range,
        template_crbn_chain
    )
    results['electrostatics'] = electro_result
    
    # Step 6: 生成报告
    print("📍 Step 6: 生成设计报告")
    if output_report is None:
        output_report = f"glue_design_report_{target_obj}.txt"
    
    with open(output_report, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("分子胶设计分析报告\n")
        f.write("Molecular Glue Design Analysis Report\n")
        f.write("="*80 + "\n\n")
        
        f.write(f"模板: {template_obj} ({template_gloop_chain}:{template_gloop_range})\n")
        f.write(f"目标: {target_obj} ({target_chain}:{target_gloop_range})\n")
        f.write(f"配体: {ligand_resname}\n\n")
        
        f.write("-"*80 + "\n")
        f.write("1. G-loop 对齐\n")
        f.write("-"*80 + "\n")
        f.write(f"RMSD: {align_result['rmsd']:.3f} Å\n")
        f.write(f"对齐质量: {'✅ 优秀 (<1Å)' if align_result['rmsd'] < 1.0 else '⚠️ 一般'}\n\n")
        
        f.write("-"*80 + "\n")
        f.write("2. 碰撞分析\n")
        f.write("-"*80 + "\n")
        f.write(f"碰撞数: {clash_result['clash_count']}\n")
        f.write(f"接触数: {clash_result['contact_count']}\n")
        f.write(f"严重碰撞: {'❌ 是' if clash_result['has_severe_clashes'] else '✅ 否'}\n")
        if clash_result['clashes']:
            f.write("\n主要碰撞:\n")
            for clash in clash_result['clashes'][:5]:
                f.write(f"  {clash['gloop_res']} ↔ {clash['crbn_res']}: {clash['distance']} Å\n")
        f.write("\n")
        
        if results.get('geometry'):
            f.write("-"*80 + "\n")
            f.write("3. G-motif 验证\n")
            f.write("-"*80 + "\n")
            geom = results['geometry']
            f.write(f"中心 Gly: {'✅' if geom['central_gly_confirmed'] else '❌'}\n")
            f.write(f"α-turn 氢键: {'✅' if geom['has_alpha_turn_hbond'] else '❌'} ({geom['alpha_turn_distance']} Å)\n")
            f.write(f"RMSD: {geom['rmsd_to_gspt1']} Å\n")
            f.write(f"有效几何: {'✅' if geom['is_valid_geometry'] else '❌'}\n\n")
        
        if exit_vectors:
            f.write("-"*80 + "\n")
            f.write("4. Exit Vector 推荐\n")
            f.write("-"*80 + "\n")
            for i, ev in enumerate(exit_vectors[:3], 1):
                f.write(f"{i}. {ev['atom_name']} ({ev['element']}) - 评分 {ev['score']}\n")
                f.write(f"   原因: {ev['reasons']}\n")
            f.write("\n")
        
        if electro_result:
            f.write("-"*80 + "\n")
            f.write("5. 静电环境与设计建议\n")
            f.write("-"*80 + "\n")
            f.write(f"净电荷: {electro_result['net_charge']:+d}\n")
            f.write(f"疏水比例: {electro_result['hydrophobic_ratio']:.1%}\n")
            f.write(f"\n💡 建议: {electro_result['recommendation']}\n\n")
        
        f.write("="*80 + "\n")
        f.write("分析完成\n")
        f.write("="*80 + "\n")
    
    print(f"\n✅ 报告已生成: {output_report}\n")
    print("🔬"*35 + "\n")
    
    return results


# ====== 导出到命名空间 ======
if __name__ == "__main__":
    print("分子胶设计分析模块已加载")
    print("主要功能:")
    print("  - align_gloop_for_modeling()")
    print("  - detect_clashes_at_interface()")
    print("  - identify_exit_vectors()")
    print("  - analyze_electrostatic_environment()")
    print("  - comprehensive_glue_design_analysis()")
