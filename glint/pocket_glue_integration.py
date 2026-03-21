# -*- coding: utf-8 -*-
"""
pocket_glue_integration.py
口袋分析与分子胶功能联动Module（精简版）

核心功能：
1. 口袋与 PPI interface关联分析
2. 分子胶前后口袋对比
3. 口袋-相互作用关联
4. G-motif 口袋综合分析
"""

import os
import csv
import numpy as np

try:
    from pymol import cmd
except ImportError:
    cmd = None

from .pocket_detector import detect_pockets, compare_pockets


# 从Package级别Import统一的中文检测Function，避免重复定义
from . import _zh


def _info(cn, en):
    """双语Information输出"""
    print(cn if _zh() else en)


def analyze_pockets_in_ppi_interface(obj_name, chain_a, chain_b, 
                                     output_csv=None, visualize=True,
                                     **pocket_kwargs):
    """
    分析 PPI interface上的口袋
    
    Parameters：
        obj_name: PyMOL 对象名
        chain_a, chain_b: 两个蛋白链
        output_csv: 输出 CSV
        visualize: 是否可视化
        **pocket_kwargs: 传递给 detect_pockets 的Parameters
    
    Return：interface口袋列表
    """
    if not cmd:
        _info("Error：需要 PyMOL 环境", "Error: PyMOL required")
        return []
    
    _info(f"\n分析 {obj_name} 中链 {chain_a} 和 {chain_b} 的interface口袋...",
          f"\nAnalyzing interface pockets in {obj_name} between chains {chain_a} and {chain_b}...")
    
    all_pockets = detect_pockets(obj_name=obj_name, **pocket_kwargs)
    
    # 获取interface残基
    try:
        from .ppi_analyzer import analyze_protein_protein_interface
        ppi_result = analyze_protein_protein_interface(obj_name, [chain_a], [chain_b], interface_distance=5.0, visualize=False)
        
        interface_residues = set()
        for res in ppi_result.get('interface_residues', []):
            interface_residues.add((res['chain1'], res['resid1']))
            interface_residues.add((res['chain2'], res['resid2']))
    except Exception as e:
        _info(f"Warning：无法获取interface残基：{e}", f"Warning: Failed to get interface residues: {e}")
        interface_residues = set()
    
    # 筛选interface口袋
    interface_pockets = []
    for pocket in all_pockets:
        pocket_res_set = set((r['chain'], r['resi']) for r in pocket['residues'])
        overlap = pocket_res_set & interface_residues
        
        if len(overlap) > 0:
            interface_overlap = len(overlap) / len(pocket_res_set) if len(pocket_res_set) > 0 else 0
            pocket['interface_overlap'] = interface_overlap
            pocket['interface_residues'] = list(overlap)
            interface_pockets.append(pocket)
            
            _info(f"  口袋 {pocket['id']}: interface重叠degrees {interface_overlap:.2%}",
                  f"  Pocket {pocket['id']}: Interface overlap {interface_overlap:.2%}")
    
    if output_csv:
        _export_interface_pockets_csv(interface_pockets, output_csv)
    
    if visualize:
        from .pocket_visualizer import visualize_pockets
        visualize_pockets(interface_pockets, obj_name='interface_pockets',
                         color_by='druggability', show_spheres=True)
    
    _info(f"\n在interface上发现 {len(interface_pockets)} 个口袋（共 {len(all_pockets)} 个）",
          f"\nFound {len(interface_pockets)} pockets at interface (out of {len(all_pockets)})")
    
    return interface_pockets


def _export_interface_pockets_csv(pockets, output_csv):
    """Exportinterface口袋 CSV"""
    fieldnames = [
        'Pocket_ID', 'Volume_A3', 'Druggability_Score', 'Interface_Overlap',
        'Hydrophobicity', 'Net_Charge', 'Center_X', 'Center_Y', 'Center_Z'
    ]
    
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for p in pockets:
            writer.writerow({
                'Pocket_ID': p['id'],
                'Volume_A3': f"{p['volume']:.2f}",
                'Druggability_Score': f"{p['druggability_score']:.3f}",
                'Interface_Overlap': f"{p['interface_overlap']:.3f}",
                'Hydrophobicity': f"{p['hydrophobicity']:.3f}",
                'Net_Charge': p['net_charge'],
                'Center_X': f"{p['center'][0]:.2f}",
                'Center_Y': f"{p['center'][1]:.2f}",
                'Center_Z': f"{p['center'][2]:.2f}"
            })


def analyze_pockets_with_glue(obj_with_glue, obj_without_glue=None,
                              protein_a_chain=None, protein_b_chain=None,
                              output_dir=None, visualize=True,
                              **pocket_kwargs):
    """
    分析分子胶对口袋的影响
    
    Parameters：
        obj_with_glue: 含分子胶的对象名
        obj_without_glue: 不含分子胶的对象名（可选）
        protein_a_chain, protein_b_chain: 两个蛋白链
        output_dir: 输出Directory
        visualize: 是否可视化
        **pocket_kwargs: 传递给 detect_pockets 的Parameters
    
    Return：(pockets_with, pockets_without, comparison)
    """
    if output_dir: output_dir = os.path.normpath(output_dir)

    if not cmd:
        return None, None, None
    
    _info("\n=== 分子胶口袋分析 ===", "\n=== Molecular Glue Pocket Analysis ===")
    
    # 分析含分子胶的结构
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
        
        pockets_with_full, pockets_without_full, comparison = compare_pockets(
            obj_without_glue, obj_with_glue, align=True,
            output_csv=os.path.join(output_dir, 'pocket_comparison.csv') if output_dir else None,
            **pocket_kwargs
        )
        
        new_pockets = [c for c in comparison if c['match_type'] == 'new']
        expanded_pockets = [c for c in comparison if c['match_type'] == 'matched' and c['delta_volume'] > 10]
        
        if len(new_pockets) > 0 or len(expanded_pockets) > 0:
            _info(f"\n⭐ 分子胶效应：", f"\n⭐ Molecular Glue Effect:")
            _info(f"  - 新增口袋：{len(new_pockets)} 个", f"  - New pockets: {len(new_pockets)}")
            _info(f"  - 扩大口袋：{len(expanded_pockets)} 个", f"  - Expanded pockets: {len(expanded_pockets)}")
    
    if visualize:
        from .pocket_visualizer import visualize_pocket_comparison, visualize_pockets
        if pockets_without:
            visualize_pocket_comparison(comparison, pockets_without, pockets_with)
        else:
            visualize_pockets(pockets_with, obj_name='glue_pockets', color_by='druggability')
    
    return pockets_with, pockets_without, comparison


def correlate_pockets_with_interactions(pockets, interaction_csv, output_csv=None):
    """
    关联口袋与相互作用分析
    
    Parameters：
        pockets: detect_pockets Return的口袋列表
        interaction_csv: 相互作用 CSV File
        output_csv: 输出 CSV
    
    Return：口袋-相互作用关联列表
    """
    if not os.path.exists(interaction_csv):
        _info(f"Error：未找到相互作用File {interaction_csv}", f"Error: Interaction file not found")
        return []
    
    interactions = []
    with open(interaction_csv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            interactions.append(row)
    
    correlations = []
    for pocket in pockets:
        pocket_residues = set((r['chain'], r['resi']) for r in pocket['residues'])
        
        interactions_in_pocket = []
        interaction_types = {}
        
        for inter in interactions:
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
    
    if output_csv:
        with open(output_csv, 'w', newline='', encoding='utf-8') as f:
            fieldnames = ['Pocket_ID', 'Volume_A3', 'Druggability_Score',
                         'Num_Interactions', 'Interaction_Density']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for corr in correlations:
                writer.writerow({
                    'Pocket_ID': corr['pocket_id'],
                    'Volume_A3': f"{corr['pocket_volume']:.2f}",
                    'Druggability_Score': f"{corr['druggability_score']:.3f}",
                    'Num_Interactions': corr['num_interactions'],
                    'Interaction_Density': f"{corr['interaction_density']:.4f}"
                })
    
    return correlations


def integrate_pockets_with_electrostatics(obj_name, pockets=None,
                                          output_csv=None, visualize=True,
                                          **pocket_kwargs):
    """
    整合口袋检测与静电势分析

    Parameters：
        obj_name: PyMOL 对象名
        pockets: 预先检测的口袋列表（可选，None 则自动检测）
        output_csv: 输出 CSV FilePath
        visualize: 是否可视化
        **pocket_kwargs: 传递给 detect_pockets 的Parameters

    Return：
        带有静电势Information的口袋列表
    """
    if not cmd:
        return None

    # 检测口袋
    if pockets is None:
        pockets = detect_pockets(obj_name, **pocket_kwargs)

    if not pockets:
        _info("⚠️  未检测到口袋", "⚠️  No pockets detected")
        return []

    _info(f"\n🔬 整合 {len(pockets)} 个口袋的静电势Information...",
          f"\n🔬 Integrating electrostatics for {len(pockets)} pockets...")

    # 为每个口袋计算静电势特征
    for pocket in pockets:
        center = pocket.get('center', [0, 0, 0])

        # 计算口袋内残基的净电荷
        pocket_residues = pocket.get('residues', [])
        positive_count = 0
        negative_count = 0

        for res in pocket_residues:
            resn = res.get('resn', '')
            if resn in ['ARG', 'LYS', 'HIS']:
                positive_count += 1
            elif resn in ['ASP', 'GLU']:
                negative_count += 1

        pocket['positive_residues'] = positive_count
        pocket['negative_residues'] = negative_count
        pocket['net_charge'] = positive_count - negative_count
        pocket['charge_ratio'] = positive_count / max(negative_count, 1)

        # 静电势极性评估
        if pocket['net_charge'] > 2:
            pocket['electrostatic_character'] = 'positive'
        elif pocket['net_charge'] < -2:
            pocket['electrostatic_character'] = 'negative'
        else:
            pocket['electrostatic_character'] = 'neutral'

    # 输出 CSV
    if output_csv:
        import csv
        fieldnames = [
            'Pocket_ID', 'Volume_A3', 'Druggability_Score',
            'Positive_Residues', 'Negative_Residues', 'Net_Charge',
            'Electrostatic_Character', 'Center_X', 'Center_Y', 'Center_Z'
        ]
        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for i, p in enumerate(pockets, 1):
                writer.writerow({
                    'Pocket_ID': i,
                    'Volume_A3': f"{p.get('volume', 0):.1f}",
                    'Druggability_Score': f"{p.get('druggability_score', 0):.2f}",
                    'Positive_Residues': p.get('positive_residues', 0),
                    'Negative_Residues': p.get('negative_residues', 0),
                    'Net_Charge': p.get('net_charge', 0),
                    'Electrostatic_Character': p.get('electrostatic_character', 'unknown'),
                    'Center_X': f"{p['center'][0]:.2f}",
                    'Center_Y': f"{p['center'][1]:.2f}",
                    'Center_Z': f"{p['center'][2]:.2f}"
                })
        _info(f"✅ Results已Save：{output_csv}", f"✅ Results saved: {output_csv}")

    # 可视化
    if visualize:
        from .pocket_visualizer import visualize_pockets
        visualize_pockets(pockets, obj_name='electrostatic_pockets', color_by='druggability')

    return pockets


def comprehensive_glue_pocket_analysis(obj_name, chain_a, chain_b,
                                       glue_selection=None,
                                       output_dir='glue_pocket_analysis'):
    """
    分子胶口袋综合分析（一Key式）

    Parameters：
        obj_name: PyMOL 对象名
        chain_a: 蛋白 A 链（如 E3）
        chain_b: 蛋白 B 链（如底物）
        glue_selection: 分子胶Select表达式（可选）
        output_dir: 输出Directory

    Return：综合分析Results字典
    """
    if output_dir: output_dir = os.path.normpath(output_dir)

    if not cmd:
        return None

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    _info("\n" + "="*60, "\n" + "="*60)
    _info("  分子胶口袋综合分析", "  Comprehensive Glue Pocket Analysis")
    _info("="*60, "="*60)

    results = {}

    # 1. 全局口袋检测
    _info("\n📍 步骤 1: 全局口袋检测...", "\n📍 Step 1: Global pocket detection...")
    all_pockets = detect_pockets(
        obj_name=obj_name,
        output_csv=os.path.join(output_dir, 'all_pockets.csv')
    )
    results['all_pockets'] = all_pockets
    _info(f"   检测到 {len(all_pockets)} 个口袋", f"   Detected {len(all_pockets)} pockets")

    # 2. interface口袋分析
    _info("\n📍 步骤 2: interface口袋分析...", "\n📍 Step 2: Interface pocket analysis...")
    interface_pockets = analyze_pockets_in_ppi_interface(
        obj_name, chain_a, chain_b,
        output_csv=os.path.join(output_dir, 'interface_pockets.csv'),
        visualize=False
    )
    results['interface_pockets'] = interface_pockets
    _info(f"   interface口袋: {len(interface_pockets)} 个", f"   Interface pockets: {len(interface_pockets)}")

    # 3. 静电势整合
    _info("\n📍 步骤 3: 静电势整合...", "\n📍 Step 3: Electrostatics integration...")
    electrostatic_pockets = integrate_pockets_with_electrostatics(
        obj_name, pockets=all_pockets,
        output_csv=os.path.join(output_dir, 'electrostatic_pockets.csv'),
        visualize=False
    )
    results['electrostatic_pockets'] = electrostatic_pockets

    # 4. 分子胶影响分析（如果提供）
    if glue_selection:
        _info("\n📍 步骤 4: 分子胶影响分析...", "\n📍 Step 4: Glue impact analysis...")
        # 这里可以Extension更多分析
        results['glue_selection'] = glue_selection

    # 5. 可视化
    from .pocket_visualizer import visualize_pockets
    if interface_pockets:
        visualize_pockets(interface_pockets, obj_name='interface_pockets', color_by='druggability')

    _info(f"\n✅ 分析Completed！ResultsSave在：{output_dir}", f"\n✅ Analysis complete! Results saved to: {output_dir}")

    return results


def comprehensive_gmotif_pocket_analysis(obj_name, e3_chain, substrate_chain,
                                         glue_chain=None, 
                                         output_dir='gmotif_pocket_analysis'):
    """
    G-motif 与口袋综合分析（一Key式）
    
    Parameters：
        obj_name: PyMOL 对象名
        e3_chain: E3 ligase 链
        substrate_chain: 底物链
        glue_chain: 分子胶链（可选）
        output_dir: 输出Directory
    
    Return：综合分析Results字典
    """
    if output_dir: output_dir = os.path.normpath(output_dir)

    if not cmd:
        return None
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    _info("\n" + "="*60, "\n" + "="*60)
    _info("  G-motif 口袋综合分析", "  Comprehensive G-motif Pocket Analysis")
    _info("="*60, "="*60)
    
    results = {}
    
    # 1. G-motif 识别
    try:
        from .g_motif_analyzer import find_crbn_g_motif, analyze_g_motif_glue_binding
        
        gmotif_results = find_crbn_g_motif(
            obj_name, e3_chain, substrate_chain,
            template_mode='builtin',
            output_csv=os.path.join(output_dir, 'gmotif_detection.csv')
        )
        results['gmotif'] = gmotif_results
        
        if len(gmotif_results) == 0:
            _info("⚠️  未检测到 G-motif", "⚠️  No G-motif detected")
            return results
        
        best_gmotif = gmotif_results[0]
        gmotif_range = f"{best_gmotif['start_resi']}-{best_gmotif['end_resi']}"
    except Exception as e:
        _info(f"Error：G-motif 识别Failed：{e}", f"Error: G-motif detection failed: {e}")
        return results
    
    # 2. G-motif 周围口袋检测
    pockets_gmotif = detect_pockets(
        obj_name=obj_name,
        selection=f'byres ((chain {substrate_chain} and resi {gmotif_range}) around 15)',
        output_csv=os.path.join(output_dir, 'gmotif_pockets.csv')
    )
    results['pockets'] = pockets_gmotif
    
    # 3. G-motif 与分子胶相互作用
    if glue_chain:
        try:
            gmotif_glue_result = analyze_g_motif_glue_binding(
                obj_name, e3_chain, substrate_chain, glue_chain,
                output_csv=os.path.join(output_dir, 'gmotif_glue_interactions.csv')
            )
            results['gmotif_glue_interactions'] = gmotif_glue_result
            
            correlations = correlate_pockets_with_interactions(
                pockets_gmotif,
                os.path.join(output_dir, 'gmotif_glue_interactions.csv'),
                output_csv=os.path.join(output_dir, 'gmotif_pocket_correlations.csv')
            )
            results['correlations'] = correlations
        except Exception as e:
            _info(f"Warning：G-motif-分子胶分析Failed：{e}", f"Warning: G-motif-glue analysis failed: {e}")
    
    # 4. interface口袋分析
    interface_pockets = analyze_pockets_in_ppi_interface(
        obj_name, e3_chain, substrate_chain,
        output_csv=os.path.join(output_dir, 'interface_pockets.csv'),
        visualize=False
    )
    results['interface_pockets'] = interface_pockets
    
    # 5. 可视化
    from .pocket_visualizer import visualize_pockets
    if len(pockets_gmotif) > 0:
        visualize_pockets(pockets_gmotif, obj_name='gmotif_pockets', color_by='druggability')
    
    try:
        cmd.select('gmotif_region', f'chain {substrate_chain} and resi {gmotif_range}')
        cmd.show('cartoon', 'gmotif_region')
        cmd.color('red', 'gmotif_region')
    except Exception:  # PyMOL 操作可能Failed
        pass
    
    _info(f"\n✅ 分析Completed！ResultsSave在：{output_dir}", f"\n✅ Analysis complete! Results saved to: {output_dir}")
    
    return results


if __name__ == '__main__':
    print("GLINT Pocket-Glue Integration - 请在 PyMOL 中using")
