# -*- coding: utf-8 -*-
"""
GlueTK Pocket-Glue Integration
================================
口袋分析与分子胶功能联动模块

功能：
- 口袋与 PPI 界面关联分析
- 口袋与相互作用联动
- 三元复合体口袋对比（有无分子胶）
- 与 APBS 静电势叠加分析
- 综合评分与可视化
"""

import os
import csv
import numpy as np

try:
    from pymol import cmd
except ImportError:
    cmd = None

from .pocket_detector import detect_pockets, compare_pockets, PocketDetector
from .pocket_visualizer import (
    visualize_pockets, 
    overlay_pocket_electrostatics,
    visualize_pocket_comparison
)


def _zh():
    """检测中文环境"""
    import locale
    try:
        lang = locale.getdefaultlocale()[0]
        return lang and lang.startswith('zh')
    except:
        return False


def _info(cn, en):
    """双语信息输出"""
    print(cn if _zh() else en)


def analyze_pockets_in_ppi_interface(obj_name, chain_a, chain_b, 
                                     output_csv=None, visualize=True,
                                     **pocket_kwargs):
    """
    分析 PPI 界面上的口袋
    
    参数：
        obj_name: PyMOL 对象名
        chain_a, chain_b: 两个蛋白链
        output_csv: 输出 CSV
        visualize: 是否可视化
        **pocket_kwargs: 传递给 detect_pockets 的参数
    
    返回：界面口袋列表
    """
    if not cmd:
        _info("错误：需要 PyMOL 环境", "Error: PyMOL required")
        return []
    
    _info(f"\n分析 {obj_name} 中链 {chain_a} 和 {chain_b} 的界面口袋...",
          f"\nAnalyzing interface pockets in {obj_name} between chains {chain_a} and {chain_b}...")
    
    # 检测整体口袋
    all_pockets = detect_pockets(obj_name=obj_name, **pocket_kwargs)
    
    # 获取界面残基（使用 PPI 分析模块）
    try:
        from .ppi_analyzer import analyze_protein_protein_interface
        ppi_result = analyze_protein_protein_interface(
            obj_name, chain_a, chain_b, cutoff=5.0
        )
        interface_residues = set()
        for res in ppi_result.get('interface_residues_A', []):
            interface_residues.add((chain_a, res))
        for res in ppi_result.get('interface_residues_B', []):
            interface_residues.add((chain_b, res))
    except Exception as e:
        _info(f"警告：无法获取界面残基信息：{e}", 
              f"Warning: Failed to get interface residues: {e}")
        interface_residues = set()
    
    # 筛选界面口袋（口袋残基与界面残基有重叠）
    interface_pockets = []
    for pocket in all_pockets:
        pocket_res_set = set((r['chain'], r['resi']) for r in pocket['residues'])
        overlap = pocket_res_set & interface_residues
        
        if len(overlap) > 0:
            # 计算界面重叠度
            interface_overlap = len(overlap) / len(pocket_res_set) if len(pocket_res_set) > 0 else 0
            pocket['interface_overlap'] = interface_overlap
            pocket['interface_residues'] = list(overlap)
            interface_pockets.append(pocket)
            
            _info(f"  口袋 {pocket['id']}: 界面重叠度 {interface_overlap:.2%}",
                  f"  Pocket {pocket['id']}: Interface overlap {interface_overlap:.2%}")
    
    # 输出 CSV
    if output_csv:
        _export_interface_pockets_csv(interface_pockets, output_csv)
    
    # 可视化
    if visualize:
        visualize_pockets(interface_pockets, obj_name='interface_pockets',
                         color_by='druggability', show_spheres=True)
    
    _info(f"\n在界面上发现 {len(interface_pockets)} 个口袋（共 {len(all_pockets)} 个）",
          f"\nFound {len(interface_pockets)} pockets at interface (out of {len(all_pockets)})")
    
    return interface_pockets


def _export_interface_pockets_csv(pockets, output_csv):
    """导出界面口袋 CSV"""
    fieldnames = [
        'Pocket_ID', 'Volume_A3', 'Druggability_Score', 'Interface_Overlap',
        'Hydrophobicity', 'Net_Charge', 'Center_X', 'Center_Y', 'Center_Z',
        'Interface_Residues'
    ]
    
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for p in pockets:
            interface_res_str = ';'.join([f"{chain}:{resi}" 
                                         for chain, resi in p['interface_residues']])
            writer.writerow({
                'Pocket_ID': p['id'],
                'Volume_A3': f"{p['volume']:.2f}",
                'Druggability_Score': f"{p['druggability_score']:.3f}",
                'Interface_Overlap': f"{p['interface_overlap']:.3f}",
                'Hydrophobicity': f"{p['hydrophobicity']:.3f}",
                'Net_Charge': p['net_charge'],
                'Center_X': f"{p['center'][0]:.2f}",
                'Center_Y': f"{p['center'][1]:.2f}",
                'Center_Z': f"{p['center'][2]:.2f}",
                'Interface_Residues': interface_res_str
            })
    
    _info(f"界面口袋已保存到 {output_csv}", f"Interface pockets saved to {output_csv}")


def analyze_pockets_with_glue(obj_with_glue, obj_without_glue=None,
                              glue_chain=None, protein_a_chain=None, protein_b_chain=None,
                              output_dir=None, visualize=True,
                              **pocket_kwargs):
    """
    分析分子胶对口袋的影响
    
    参数：
        obj_with_glue: 含分子胶的对象名
        obj_without_glue: 不含分子胶的对象名（可选，用于对比）
        glue_chain: 分子胶所在链
        protein_a_chain, protein_b_chain: 两个蛋白链
        output_dir: 输出目录
        visualize: 是否可视化
        **pocket_kwargs: 传递给 detect_pockets 的参数
    
    返回：(pockets_with, pockets_without, comparison)
    """
    if not cmd:
        return None, None, None
    
    _info("\n=== 分子胶口袋分析 ===", "\n=== Molecular Glue Pocket Analysis ===")
    
    # 分析含分子胶的结构
    _info(f"\n1. 分析 {obj_with_glue}（含分子胶）...", 
          f"\n1. Analyzing {obj_with_glue} (with glue)...")
    
    # 如果指定了界面链，优先分析界面口袋
    if protein_a_chain and protein_b_chain:
        pockets_with = analyze_pockets_in_ppi_interface(
            obj_with_glue, protein_a_chain, protein_b_chain,
            output_csv=os.path.join(output_dir, 'pockets_with_glue.csv') if output_dir else None,
            visualize=False,
            **pocket_kwargs
        )
    else:
        pockets_with = detect_pockets(
            obj_name=obj_with_glue,
            output_csv=os.path.join(output_dir, 'pockets_with_glue.csv') if output_dir else None,
            **pocket_kwargs
        )
    
    pockets_without = None
    comparison = None
    
    # 对比分析
    if obj_without_glue:
        _info(f"\n2. 分析 {obj_without_glue}（无分子胶）...", 
              f"\n2. Analyzing {obj_without_glue} (without glue)...")
        
        if protein_a_chain and protein_b_chain:
            pockets_without = analyze_pockets_in_ppi_interface(
                obj_without_glue, protein_a_chain, protein_b_chain,
                output_csv=os.path.join(output_dir, 'pockets_without_glue.csv') if output_dir else None,
                visualize=False,
                **pocket_kwargs
            )
        else:
            pockets_without = detect_pockets(
                obj_name=obj_without_glue,
                output_csv=os.path.join(output_dir, 'pockets_without_glue.csv') if output_dir else None,
                **pocket_kwargs
            )
        
        # 对比
        _info("\n3. 对比口袋变化...", "\n3. Comparing pocket changes...")
        pockets_with_full, pockets_without_full, comparison = compare_pockets(
            obj_without_glue, obj_with_glue,
            align=True,
            output_csv=os.path.join(output_dir, 'pocket_comparison.csv') if output_dir else None,
            **pocket_kwargs
        )
        
        # 分析新出现/扩大的口袋
        new_pockets = [c for c in comparison if c['match_type'] == 'new']
        expanded_pockets = [c for c in comparison if c['match_type'] == 'matched' and c['delta_volume'] > 10]
        
        if len(new_pockets) > 0 or len(expanded_pockets) > 0:
            _info(f"\n⭐ 分子胶效应：", f"\n⭐ Molecular Glue Effect:")
            _info(f"  - 新增口袋：{len(new_pockets)} 个", 
                  f"  - New pockets: {len(new_pockets)}")
            _info(f"  - 扩大口袋：{len(expanded_pockets)} 个", 
                  f"  - Expanded pockets: {len(expanded_pockets)}")
    
    # 可视化
    if visualize:
        if pockets_without:
            visualize_pocket_comparison(comparison, pockets_without, pockets_with)
        else:
            visualize_pockets(pockets_with, obj_name='glue_pockets',
                             color_by='druggability', show_spheres=True)
    
    return pockets_with, pockets_without, comparison


def correlate_pockets_with_interactions(pockets, interaction_csv, output_csv=None):
    """
    关联口袋与相互作用分析
    
    参数：
        pockets: detect_pockets 返回的口袋列表
        interaction_csv: 相互作用 CSV 文件
        output_csv: 输出 CSV
    
    返回：口袋-相互作用关联列表
    """
    if not os.path.exists(interaction_csv):
        _info(f"错误：未找到相互作用文件 {interaction_csv}",
              f"Error: Interaction file not found {interaction_csv}")
        return []
    
    _info(f"\n关联口袋与相互作用（{interaction_csv}）...",
          f"\nCorrelating pockets with interactions ({interaction_csv})...")
    
    # 读取相互作用
    interactions = []
    with open(interaction_csv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            interactions.append(row)
    
    # 为每个口袋统计相互作用
    correlations = []
    for pocket in pockets:
        pocket_residues = set((r['chain'], r['resi']) for r in pocket['residues'])
        
        # 统计在此口袋内的相互作用
        interactions_in_pocket = []
        interaction_types = {}
        
        for inter in interactions:
            # 提取残基信息（支持多种 CSV 格式）
            res_key = None
            inter_type = inter.get('Interaction', inter.get('InteractionType', 'Unknown'))
            
            if 'Protein_Chain' in inter and 'Protein_Residue' in inter:
                chain = inter['Protein_Chain']
                resi = inter['Protein_Residue'].split()[-1]
                res_key = (chain, resi)
            elif 'Chain1' in inter and 'Residue1' in inter:
                chain = inter['Chain1']
                resi = inter['Residue1'].split()[-1]
                res_key = (chain, resi)
            
            if res_key and res_key in pocket_residues:
                interactions_in_pocket.append(inter)
                interaction_types[inter_type] = interaction_types.get(inter_type, 0) + 1
        
        correlation = {
            'pocket_id': pocket['id'],
            'pocket_volume': pocket['volume'],
            'druggability_score': pocket['druggability_score'],
            'num_interactions': len(interactions_in_pocket),
            'interaction_types': interaction_types,
            'interaction_density': len(interactions_in_pocket) / pocket['volume'] if pocket['volume'] > 0 else 0
        }
        correlations.append(correlation)
        
        _info(f"  口袋 {pocket['id']}: {len(interactions_in_pocket)} 个相互作用",
              f"  Pocket {pocket['id']}: {len(interactions_in_pocket)} interactions")
        for itype, count in interaction_types.items():
            print(f"    - {itype}: {count}")
    
    # 输出 CSV
    if output_csv:
        _export_correlation_csv(correlations, output_csv)
    
    return correlations


def _export_correlation_csv(correlations, output_csv):
    """导出口袋-相互作用关联 CSV"""
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        fieldnames = [
            'Pocket_ID', 'Volume_A3', 'Druggability_Score',
            'Num_Interactions', 'Interaction_Density', 'Interaction_Types'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for corr in correlations:
            types_str = ';'.join([f"{k}:{v}" for k, v in corr['interaction_types'].items()])
            writer.writerow({
                'Pocket_ID': corr['pocket_id'],
                'Volume_A3': f"{corr['pocket_volume']:.2f}",
                'Druggability_Score': f"{corr['druggability_score']:.3f}",
                'Num_Interactions': corr['num_interactions'],
                'Interaction_Density': f"{corr['interaction_density']:.4f}",
                'Interaction_Types': types_str
            })
    
    _info(f"关联结果已保存到 {output_csv}", f"Correlation saved to {output_csv}")


def integrate_pockets_with_electrostatics(pockets, obj_name, apbs_dx_file=None,
                                          output_csv=None, visualize=True):
    """
    整合口袋与 APBS 静电势分析
    
    参数：
        pockets: detect_pockets 返回的口袋列表
        obj_name: PyMOL 对象名
        apbs_dx_file: APBS 输出的 .dx 文件
        output_csv: 输出 CSV
        visualize: 是否可视化
    
    返回：口袋静电势分析结果
    """
    if not apbs_dx_file or not os.path.exists(apbs_dx_file):
        _info("警告：未提供 APBS 地图文件，跳过静电势分析",
              "Warning: APBS map not provided, skipping electrostatics")
        return []
    
    _info(f"\n整合口袋与静电势（{apbs_dx_file}）...",
          f"\nIntegrating pockets with electrostatics ({apbs_dx_file})...")
    
    # TODO: 实现从 .dx 文件读取静电势值
    # 这里简化：仅根据 net_charge 估计
    electro_results = []
    for pocket in pockets:
        # 简化版：使用口袋的净电荷作为静电势指标
        electro_potential = pocket['net_charge'] * 10  # 粗略估计（kT/e）
        
        result = {
            'pocket_id': pocket['id'],
            'electro_potential': electro_potential,
            'net_charge': pocket['net_charge'],
            'polarity': pocket['polarity'],
            'center': pocket['center']
        }
        electro_results.append(result)
        
        _info(f"  口袋 {pocket['id']}: 静电势 ≈ {electro_potential:.1f} kT/e",
              f"  Pocket {pocket['id']}: Electrostatic potential ≈ {electro_potential:.1f} kT/e")
        
        # 可视化
        if visualize and cmd:
            overlay_pocket_electrostatics(pocket, apbs_dx_file, obj_name)
    
    # 输出 CSV
    if output_csv:
        with open(output_csv, 'w', newline='', encoding='utf-8') as f:
            fieldnames = ['Pocket_ID', 'Electro_Potential_kT_e', 'Net_Charge', 'Polarity']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for res in electro_results:
                writer.writerow({
                    'Pocket_ID': res['pocket_id'],
                    'Electro_Potential_kT_e': f"{res['electro_potential']:.2f}",
                    'Net_Charge': res['net_charge'],
                    'Polarity': f"{res['polarity']:.3f}"
                })
        
        _info(f"静电势结果已保存到 {output_csv}", f"Electrostatics saved to {output_csv}")
    
    return electro_results


def comprehensive_glue_pocket_analysis(obj_with_glue, obj_without_glue=None,
                                      protein_a_chain=None, protein_b_chain=None,
                                      interaction_csv=None, apbs_dx_file=None,
                                      output_dir='glue_pocket_analysis',
                                      **pocket_kwargs):
    """
    综合分子胶口袋分析（一键式）
    
    参数：
        obj_with_glue: 含分子胶的对象
        obj_without_glue: 不含分子胶的对象（对比用）
        protein_a_chain, protein_b_chain: 蛋白链
        interaction_csv: 相互作用 CSV
        apbs_dx_file: APBS 静电势文件
        output_dir: 输出目录
        **pocket_kwargs: detect_pockets 参数
    
    返回：综合分析结果字典
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    _info("\n" + "="*60, "\n" + "="*60)
    _info("  综合分子胶口袋分析", "  Comprehensive Glue Pocket Analysis")
    _info("="*60, "="*60)
    
    results = {}
    
    # 1. 口袋检测与对比
    pockets_with, pockets_without, comparison = analyze_pockets_with_glue(
        obj_with_glue, obj_without_glue,
        protein_a_chain=protein_a_chain,
        protein_b_chain=protein_b_chain,
        output_dir=output_dir,
        visualize=True,
        **pocket_kwargs
    )
    results['pockets_with_glue'] = pockets_with
    results['pockets_without_glue'] = pockets_without
    results['comparison'] = comparison
    
    # 2. 口袋-相互作用关联
    if interaction_csv and pockets_with:
        correlations = correlate_pockets_with_interactions(
            pockets_with, interaction_csv,
            output_csv=os.path.join(output_dir, 'pocket_interaction_correlation.csv')
        )
        results['interaction_correlations'] = correlations
    
    # 3. 口袋-静电势整合
    if apbs_dx_file and pockets_with:
        electro_results = integrate_pockets_with_electrostatics(
            pockets_with, obj_with_glue, apbs_dx_file,
            output_csv=os.path.join(output_dir, 'pocket_electrostatics.csv'),
            visualize=True
        )
        results['electrostatics'] = electro_results
    
    # 4. 生成总结报告
    _generate_summary_report(results, os.path.join(output_dir, 'SUMMARY.txt'))
    
    _info("\n✅ 综合分析完成，结果已保存到：" + output_dir,
          "\n✅ Comprehensive analysis complete, results saved to: " + output_dir)
    
    return results


def _generate_summary_report(results, output_file):
    """生成总结报告"""
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("="*60 + "\n")
        f.write("  GlueTK - Comprehensive Glue Pocket Analysis Report\n")
        f.write("="*60 + "\n\n")
        
        # 口袋统计
        f.write("1. Pocket Statistics\n")
        f.write("-" * 40 + "\n")
        if results.get('pockets_with_glue'):
            f.write(f"  Pockets with glue: {len(results['pockets_with_glue'])}\n")
            avg_volume = np.mean([p['volume'] for p in results['pockets_with_glue']])
            avg_drug = np.mean([p['druggability_score'] for p in results['pockets_with_glue']])
            f.write(f"  Average volume: {avg_volume:.2f} Ų\n")
            f.write(f"  Average druggability: {avg_drug:.3f}\n")
        
        if results.get('pockets_without_glue'):
            f.write(f"  Pockets without glue: {len(results['pockets_without_glue'])}\n")
        
        f.write("\n")
        
        # 对比结果
        if results.get('comparison'):
            f.write("2. Pocket Comparison\n")
            f.write("-" * 40 + "\n")
            comp = results['comparison']
            new_count = sum(1 for c in comp if c['match_type'] == 'new')
            lost_count = sum(1 for c in comp if c['match_type'] == 'lost')
            expanded_count = sum(1 for c in comp if c['match_type'] == 'matched' and c['delta_volume'] > 10)
            f.write(f"  New pockets: {new_count}\n")
            f.write(f"  Lost pockets: {lost_count}\n")
            f.write(f"  Expanded pockets: {expanded_count}\n")
            f.write("\n")
        
        # 相互作用统计
        if results.get('interaction_correlations'):
            f.write("3. Interaction Correlations\n")
            f.write("-" * 40 + "\n")
            for corr in results['interaction_correlations']:
                f.write(f"  Pocket {corr['pocket_id']}: "
                       f"{corr['num_interactions']} interactions, "
                       f"density {corr['interaction_density']:.4f}\n")
            f.write("\n")
        
        f.write("="*60 + "\n")
        f.write("Analysis completed successfully.\n")
    
    _info(f"总结报告已保存：{output_file}", f"Summary report saved: {output_file}")


def comprehensive_gmotif_pocket_analysis(obj_name, e3_chain, substrate_chain, 
                                         glue_chain=None, 
                                         output_dir='gmotif_pocket_analysis'):
    """
    G-motif 与口袋综合分析（一键式）
    
    参数：
        obj_name: PyMOL 对象名
        e3_chain: E3 ligase 链（如 CRBN 的链 A）
        substrate_chain: 底物链（如 GSPT1 的链 B）
        glue_chain: 分子胶链（可选）
        output_dir: 输出目录
    
    返回：综合分析结果字典
    
    功能：
    1. 识别 G-motif
    2. 检测 G-motif 周围口袋
    3. 分析分子胶与 G-motif 口袋的相互作用
    4. 关联口袋与相互作用
    5. 分析 PPI 界面口袋
    6. 可视化所有结果
    7. 生成综合报告
    """
    if not cmd:
        _info("错误：需要 PyMOL 环境", "Error: PyMOL required")
        return None
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    _info("\n" + "="*60, "\n" + "="*60)
    _info("  G-motif 口袋综合分析", "  Comprehensive G-motif Pocket Analysis")
    _info("="*60, "="*60)
    
    results = {}
    
    # 1. G-motif 识别
    _info("\n=== Step 1: G-motif 识别 ===", "\n=== Step 1: G-motif Detection ===")
    try:
        from .g_motif_analyzer import find_crbn_g_motif, analyze_g_motif_glue_binding
        
        gmotif_results = find_crbn_g_motif(
            obj_name, e3_chain, substrate_chain,
            template_mode='builtin',
            output_csv=os.path.join(output_dir, 'gmotif_detection.csv')
        )
        results['gmotif'] = gmotif_results
        
        if len(gmotif_results) == 0:
            _info("⚠️  未检测到 G-motif，分析终止", "⚠️  No G-motif detected, analysis terminated")
            return results
        
        best_gmotif = gmotif_results[0]
        _info(f"✅ 检测到 G-motif: {best_gmotif['start_resi']}-{best_gmotif['end_resi']}",
              f"✅ G-motif detected: {best_gmotif['start_resi']}-{best_gmotif['end_resi']}")
        _info(f"   RMSD: {best_gmotif['rmsd']:.2f} Å", f"   RMSD: {best_gmotif['rmsd']:.2f} Å")
    except Exception as e:
        _info(f"错误：G-motif 识别失败：{e}", f"Error: G-motif detection failed: {e}")
        return results
    
    # 2. G-motif 周围口袋检测
    _info("\n=== Step 2: G-motif 口袋检测 ===", "\n=== Step 2: G-motif Pocket Detection ===")
    gmotif_range = f"{best_gmotif['start_resi']}-{best_gmotif['end_resi']}"
    
    pockets_gmotif = detect_pockets(
        obj_name=obj_name,
        selection=f'byres ((chain {substrate_chain} and resi {gmotif_range}) around 15)',
        min_volume=15.0,  # G-motif 口袋通常较小
        output_csv=os.path.join(output_dir, 'gmotif_pockets.csv')
    )
    results['pockets'] = pockets_gmotif
    _info(f"✅ 检测到 {len(pockets_gmotif)} 个 G-motif 周围口袋",
          f"✅ Detected {len(pockets_gmotif)} pockets around G-motif")
    
    # 3. G-motif 与分子胶相互作用
    if glue_chain:
        _info("\n=== Step 3: G-motif-分子胶相互作用 ===", 
              "\n=== Step 3: G-motif-Glue Interactions ===")
        try:
            gmotif_glue_result = analyze_g_motif_glue_binding(
                obj_name, e3_chain, substrate_chain, glue_chain,
                output_csv=os.path.join(output_dir, 'gmotif_glue_interactions.csv')
            )
            results['gmotif_glue_interactions'] = gmotif_glue_result
            
            # 4. 口袋-相互作用关联
            _info("\n=== Step 4: 口袋-相互作用关联 ===", 
                  "\n=== Step 4: Pocket-Interaction Correlation ===")
            correlations = correlate_pockets_with_interactions(
                pockets_gmotif,
                os.path.join(output_dir, 'gmotif_glue_interactions.csv'),
                output_csv=os.path.join(output_dir, 'gmotif_pocket_correlations.csv')
            )
            results['correlations'] = correlations
            
            # 找出主要结合口袋
            if len(correlations) > 0:
                main_pocket = max(correlations, key=lambda c: c['num_interactions'])
                _info(f"\n✅ 主要结合口袋: #{main_pocket['pocket_id']}",
                      f"\n✅ Main binding pocket: #{main_pocket['pocket_id']}")
                print(f"   相互作用数 / Interactions: {main_pocket['num_interactions']}")
                print(f"   相互作用密度 / Density: {main_pocket['interaction_density']:.4f} /Ų")
                print(f"   可成药性 / Druggability: {main_pocket['druggability_score']:.3f}")
                results['main_binding_pocket'] = main_pocket
        except Exception as e:
            _info(f"警告：G-motif-分子胶分析失败：{e}", 
                  f"Warning: G-motif-glue analysis failed: {e}")
    
    # 5. 界面口袋分析
    _info("\n=== Step 5: PPI 界面口袋分析 ===", "\n=== Step 5: PPI Interface Pockets ===")
    interface_pockets = analyze_pockets_in_ppi_interface(
        obj_name, e3_chain, substrate_chain,
        output_csv=os.path.join(output_dir, 'interface_pockets.csv'),
        visualize=False
    )
    
    # 标记与 G-motif 重叠的口袋
    gmotif_resis = set(range(int(best_gmotif['start_resi']), 
                             int(best_gmotif['end_resi'])+1))
    
    gmotif_interface_pockets = []
    for pocket in interface_pockets:
        pocket_resis = set()
        for r in pocket['residues']:
            if r['chain'] == substrate_chain:
                try:
                    pocket_resis.add(int(r['resi']))
                except:
                    pass
        
        overlap = pocket_resis & gmotif_resis
        
        if len(overlap) > 0:
            pocket['gmotif_overlap_count'] = len(overlap)
            pocket['gmotif_overlap_ratio'] = len(overlap) / len(gmotif_resis)
            gmotif_interface_pockets.append(pocket)
    
    results['gmotif_interface_pockets'] = gmotif_interface_pockets
    _info(f"✅ 与 G-motif 相关的界面口袋: {len(gmotif_interface_pockets)} 个",
          f"✅ Interface pockets related to G-motif: {len(gmotif_interface_pockets)}")
    
    # 6. 可视化
    _info("\n=== Step 6: 可视化 ===", "\n=== Step 6: Visualization ===")
    
    # 可视化 G-motif 周围口袋
    if len(pockets_gmotif) > 0:
        visualize_pockets(pockets_gmotif, obj_name='gmotif_pockets',
                         color_by='druggability', show_spheres=True)
    
    # 高亮 G-motif 区域
    try:
        cmd.select('gmotif_region', 
                   f'chain {substrate_chain} and resi {gmotif_range}')
        cmd.show('cartoon', 'gmotif_region')
        cmd.color('red', 'gmotif_region')
        _info("✅ G-motif 区域已高亮（红色）", "✅ G-motif region highlighted (red)")
    except Exception as e:
        _info(f"警告：G-motif 高亮失败：{e}", f"Warning: G-motif highlighting failed: {e}")
    
    # 7. 生成报告
    _generate_gmotif_pocket_report(results, output_dir)
    
    _info(f"\n✅ 分析完成！结果保存在：{output_dir}",
          f"\n✅ Analysis complete! Results saved to: {output_dir}")
    
    return results


def _generate_gmotif_pocket_report(results, output_dir):
    """生成 G-motif 口袋分析报告"""
    report_path = os.path.join(output_dir, 'GMOTIF_POCKET_REPORT.txt')
    
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write("="*60 + "\n")
        f.write("  GlueTK - G-motif Pocket Analysis Report\n")
        f.write("="*60 + "\n\n")
        
        # G-motif 信息
        if results.get('gmotif'):
            gmotif = results['gmotif'][0]
            f.write("1. G-motif Detection\n")
            f.write("-" * 40 + "\n")
            f.write(f"  Position: {gmotif['start_resi']}-{gmotif['end_resi']}\n")
            f.write(f"  RMSD: {gmotif['rmsd']:.2f} Å\n")
            f.write(f"  Template: {gmotif.get('template', 'builtin')}\n")
            f.write(f"  Has Gly6: {gmotif.get('has_gly6', 'Unknown')}\n")
            f.write("\n")
        
        # 口袋统计
        if results.get('pockets'):
            pockets = results['pockets']
            f.write("2. Pocket Statistics\n")
            f.write("-" * 40 + "\n")
            f.write(f"  Total pockets around G-motif: {len(pockets)}\n")
            
            if len(pockets) > 0:
                avg_vol = np.mean([p['volume'] for p in pockets])
                avg_drug = np.mean([p['druggability_score'] for p in pockets])
                f.write(f"  Average volume: {avg_vol:.2f} Ų\n")
                f.write(f"  Average druggability: {avg_drug:.3f}\n")
                
                # 最大口袋
                max_pocket = max(pockets, key=lambda p: p['volume'])
                f.write(f"\n  Largest pocket:\n")
                f.write(f"    ID: {max_pocket['id']}\n")
                f.write(f"    Volume: {max_pocket['volume']:.2f} Ų\n")
                f.write(f"    Druggability: {max_pocket['druggability_score']:.3f}\n")
                f.write(f"    Hydrophobicity: {max_pocket['hydrophobicity']:.3f}\n")
                f.write(f"    Net charge: {max_pocket['net_charge']}\n")
            f.write("\n")
        
        # 主要结合口袋
        if results.get('main_binding_pocket'):
            main = results['main_binding_pocket']
            f.write("3. Main Binding Pocket\n")
            f.write("-" * 40 + "\n")
            f.write(f"  Pocket ID: {main['pocket_id']}\n")
            f.write(f"  Volume: {main['pocket_volume']:.2f} Ų\n")
            f.write(f"  Druggability: {main['druggability_score']:.3f}\n")
            f.write(f"  Interactions: {main['num_interactions']}\n")
            f.write(f"  Interaction density: {main['interaction_density']:.4f} /Ų\n")
            
            if 'interaction_types' in main:
                f.write(f"  Interaction types:\n")
                for itype, count in main['interaction_types'].items():
                    f.write(f"    - {itype}: {count}\n")
            f.write("\n")
        
        # G-motif 界面口袋
        if results.get('gmotif_interface_pockets'):
            gip = results['gmotif_interface_pockets']
            f.write("4. G-motif Interface Pockets\n")
            f.write("-" * 40 + "\n")
            f.write(f"  Pockets overlapping with G-motif: {len(gip)}\n")
            
            for pocket in gip:
                f.write(f"\n  Pocket {pocket['id']}:\n")
                f.write(f"    G-motif overlap: {pocket['gmotif_overlap_count']} residues ")
                f.write(f"({pocket['gmotif_overlap_ratio']:.1%})\n")
                f.write(f"    Interface overlap: {pocket.get('interface_overlap', 0):.1%}\n")
                f.write(f"    Volume: {pocket['volume']:.2f} Ų\n")
                f.write(f"    Druggability: {pocket['druggability_score']:.3f}\n")
                f.write(f"    Hydrophobicity: {pocket['hydrophobicity']:.3f}\n")
            f.write("\n")
        
        # 综合评估
        f.write("5. Comprehensive Assessment\n")
        f.write("-" * 40 + "\n")
        
        score = 0
        assessments = []
        
        # 评估 G-motif 质量
        if results.get('gmotif'):
            rmsd = results['gmotif'][0]['rmsd']
            if rmsd < 1.5:
                score += 3
                assessments.append(f"✅ Excellent G-motif match (RMSD {rmsd:.2f} Å)")
            elif rmsd < 2.5:
                score += 2
                assessments.append(f"✅ Good G-motif match (RMSD {rmsd:.2f} Å)")
            else:
                score += 1
                assessments.append(f"⚠️  Moderate G-motif match (RMSD {rmsd:.2f} Å)")
        
        # 评估口袋数量与质量
        if results.get('pockets'):
            num_pockets = len(results['pockets'])
            if num_pockets > 0:
                avg_drug = np.mean([p['druggability_score'] for p in results['pockets']])
                if num_pockets >= 2 and avg_drug > 0.5:
                    score += 3
                    assessments.append(f"✅ Multiple druggable pockets ({num_pockets}, avg drug {avg_drug:.2f})")
                elif num_pockets >= 1 and avg_drug > 0.3:
                    score += 2
                    assessments.append(f"✅ At least one pocket with moderate druggability")
                else:
                    score += 1
                    assessments.append(f"⚠️  Limited druggable pockets")
        
        # 评估相互作用
        if results.get('main_binding_pocket'):
            num_int = results['main_binding_pocket']['num_interactions']
            if num_int >= 5:
                score += 3
                assessments.append(f"✅ Strong glue-pocket interactions ({num_int})")
            elif num_int >= 3:
                score += 2
                assessments.append(f"✅ Moderate glue-pocket interactions ({num_int})")
            else:
                score += 1
                assessments.append(f"⚠️  Weak glue-pocket interactions ({num_int})")
        
        # 评估界面口袋
        if results.get('gmotif_interface_pockets'):
            num_gip = len(results['gmotif_interface_pockets'])
            if num_gip >= 2:
                score += 1
                assessments.append(f"✅ Multiple G-motif interface pockets ({num_gip})")
            elif num_gip >= 1:
                assessments.append(f"✅ G-motif interface pocket identified")
        
        f.write(f"  Overall Score: {score}/10\n\n")
        for assessment in assessments:
            f.write(f"  {assessment}\n")
        
        f.write("\n")
        f.write(f"  Conclusion: ")
        if score >= 8:
            f.write("Strong molecular glue characteristics\n")
        elif score >= 5:
            f.write("Moderate molecular glue characteristics\n")
        else:
            f.write("Weak molecular glue characteristics\n")
        
        f.write("\n")
        f.write("="*60 + "\n")
        f.write("Analysis completed successfully.\n")
    
    _info(f"报告已保存：{report_path}", f"Report saved: {report_path}")


if __name__ == '__main__':
    print("GlueTK Pocket-Glue Integration - 请在 PyMOL 中使用")
    print("示例：comprehensive_glue_pocket_analysis('complex_with_glue', 'complex_without_glue')")
    print("      comprehensive_gmotif_pocket_analysis('6bn7', 'A', 'B', 'L')")
