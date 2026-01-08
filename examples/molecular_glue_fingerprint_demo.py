#!/usr/bin/env python3
"""
Molecular Glue Fingerprint (MGF) Demo
======================================

演示如何使用MGF系统分析分子胶三元复合物

示例结构：
- 5FQD: Lenalidomide-CRBN-CK1α
- 6H0G: CC-885-CRBN-GSPT1
"""

import sys
import os
from pymol import cmd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 添加GLINT路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

# 导入MGF模块（假设已实现）
# from glint.molecular_glue_fingerprint import (
#     MolecularGlueFingerprint,
#     plot_molecular_glue_fingerprint_heatmap,
#     plot_cooperativity_radar
# )


def demo_basic_usage():
    """演示1: 基础用法 - 分析单个分子胶"""
    print("\n" + "="*70)
    print("  Demo 1: Basic Usage - Analyzing Lenalidomide")
    print("="*70)
    
    # 加载结构
    cmd.fetch('5fqd', async_=0)
    
    # 配置参数
    obj_name = '5fqd'
    e3_chains = ['A', 'B']  # CRBN二聚体
    substrate_chains = ['C']  # CK1α
    glue_resname = '1NH'  # Lenalidomide
    
    print(f"\n📥 Loaded structure: {obj_name}")
    print(f"   E3 (CRBN): chains {e3_chains}")
    print(f"   Substrate (CK1α): chains {substrate_chains}")
    print(f"   Glue: {glue_resname}")
    
    # 生成指纹（模拟输出）
    print("\n🔬 Generating molecular glue fingerprint...")
    
    # 模拟结果
    mock_result = {
        'e3_glue_fp': np.array([0.6, 0.4, 0.2, 0.5, 0.3, 0.1, 0.0, 0.2]),
        'glue_substrate_fp': np.array([0.5, 0.3, 0.1, 0.4, 0.2, 0.0, 0.1, 0.1]),
        'neo_epitope_fp': np.array([0.75, 0.30, 0.65]),
        'bridging_atoms': [
            {'atom_name': 'C10', 'resi': '1501'},
            {'atom_name': 'C14', 'resi': '1501'},
            {'atom_name': 'N3', 'resi': '1501'},
            {'atom_name': 'O2', 'resi': '1501'}
        ],
        'cooperativity_score': 0.652,
        'metadata': {
            'e3_chains': e3_chains,
            'substrate_chains': substrate_chains,
            'glue_resname': glue_resname,
            'n_e3_glue_interactions': 12,
            'n_glue_substrate_interactions': 8,
            'n_neo_epitope_contacts': 15
        }
    }
    
    # 打印结果
    print("\n📊 Results:")
    print(f"   E3-Glue interactions: {mock_result['metadata']['n_e3_glue_interactions']}")
    print(f"   Glue-Substrate interactions: {mock_result['metadata']['n_glue_substrate_interactions']}")
    print(f"   Neo-epitope contacts: {mock_result['metadata']['n_neo_epitope_contacts']}")
    print(f"   Bridging atoms: {len(mock_result['bridging_atoms'])}")
    print(f"   Cooperativity score: {mock_result['cooperativity_score']:.3f}")
    
    print("\n✅ Bridging atoms (key for glue function):")
    for i, atom in enumerate(mock_result['bridging_atoms'], 1):
        print(f"   {i}. {atom['atom_name']} (resi {atom['resi']})")
    
    # 评估
    coop = mock_result['cooperativity_score']
    if coop > 0.7:
        rating = "Excellent"
    elif coop > 0.5:
        rating = "Good"
    elif coop > 0.3:
        rating = "Fair"
    else:
        rating = "Poor"
    
    print(f"\n🎯 Overall Rating: {rating}")
    print(f"   Lenalidomide shows {rating.lower()} molecular glue properties")
    
    cmd.delete('5fqd')
    return mock_result


def demo_comparison():
    """演示2: 比较多个分子胶"""
    print("\n" + "="*70)
    print("  Demo 2: Comparing Multiple Molecular Glues")
    print("="*70)
    
    # 定义要比较的分子胶
    compounds = [
        {'pdb': '5fqd', 'name': 'Lenalidomide', 'e3': ['A', 'B'], 'sub': ['C'], 'lig': '1NH'},
        {'pdb': '6h0g', 'name': 'CC-885', 'e3': ['A'], 'sub': ['B'], 'lig': 'CC9'},
    ]
    
    results = []
    
    for comp in compounds:
        print(f"\n📥 Analyzing {comp['name']} ({comp['pdb']})...")
        cmd.fetch(comp['pdb'], async_=0)
        
        # 模拟指纹生成
        if comp['name'] == 'Lenalidomide':
            coop = 0.652
            bridging = 4
            e3_glue = 12
            glue_sub = 8
            neo = 15
        else:  # CC-885
            coop = 0.724
            bridging = 5
            e3_glue = 15
            glue_sub = 12
            neo = 18
        
        results.append({
            'Compound': comp['name'],
            'PDB': comp['pdb'],
            'Cooperativity': coop,
            'Bridging Atoms': bridging,
            'E3-Glue': e3_glue,
            'Glue-Substrate': glue_sub,
            'Neo-epitope': neo,
            'Total Interactions': e3_glue + glue_sub + neo
        })
        
        cmd.delete(comp['pdb'])
    
    # 创建比较表
    df = pd.DataFrame(results)
    print("\n📊 Comparison Table:")
    print(df.to_string(index=False))
    
    # 排序
    df_sorted = df.sort_values('Cooperativity', ascending=False)
    print(f"\n🏆 Best molecular glue: {df_sorted.iloc[0]['Compound']}")
    print(f"   Cooperativity: {df_sorted.iloc[0]['Cooperativity']:.3f}")
    
    return df


def demo_visualization():
    """演示3: 可视化指纹"""
    print("\n" + "="*70)
    print("  Demo 3: Fingerprint Visualization")
    print("="*70)
    
    # 模拟数据
    interaction_types = ['HBond', 'Hydrophobic', 'SaltBridge', 'PiStack', 
                        'PiCation', 'Metal', 'Halogen', 'Water']
    
    data = np.array([
        [0.6, 0.4, 0.2, 0.5, 0.3, 0.1, 0.0, 0.2],  # E3-Glue
        [0.5, 0.3, 0.1, 0.4, 0.2, 0.0, 0.1, 0.1],  # Glue-Substrate
        [0.4, 0.5, 0.3, 0.2, 0.1, 0.0, 0.0, 0.3],  # Neo-epitope
    ])
    
    # 创建热图
    fig, ax = plt.subplots(figsize=(12, 4))
    
    sns.heatmap(data,
                xticklabels=interaction_types,
                yticklabels=['E3-Glue', 'Glue-Substrate', 'Neo-epitope'],
                cmap='YlOrRd',
                annot=True,
                fmt='.2f',
                cbar_kws={'label': 'Normalized Strength'},
                ax=ax)
    
    ax.set_title('Molecular Glue Three-Interface Fingerprint (Lenalidomide)', 
                 fontsize=14, fontweight='bold')
    ax.set_xlabel('Interaction Type', fontsize=12)
    ax.set_ylabel('Interface', fontsize=12)
    
    plt.tight_layout()
    
    output_file = 'lenalidomide_fingerprint_heatmap.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\n✅ Heatmap saved to: {output_file}")
    
    plt.close()
    
    return output_file


def demo_sar_analysis():
    """演示4: SAR分析 - 比较类似物"""
    print("\n" + "="*70)
    print("  Demo 4: SAR Analysis - Comparing Analogs")
    print("="*70)

    # 模拟一系列类似物数据
    analogs = [
        {'name': 'Lead Compound', 'coop': 0.650, 'bridging': 4, 'activity': 8.5},
        {'name': 'Analog 1 (R1=Me)', 'coop': 0.720, 'bridging': 5, 'activity': 9.2},
        {'name': 'Analog 2 (R1=Et)', 'coop': 0.580, 'bridging': 3, 'activity': 7.8},
        {'name': 'Analog 3 (R2=Cl)', 'coop': 0.690, 'bridging': 4, 'activity': 8.9},
        {'name': 'Analog 4 (R2=F)', 'coop': 0.710, 'bridging': 5, 'activity': 9.0},
    ]

    df = pd.DataFrame(analogs)
    df['Coop_Change'] = df['coop'] - df.iloc[0]['coop']
    df['Activity_Change'] = df['activity'] - df.iloc[0]['activity']

    print("\n📊 SAR Analysis Results:")
    print(df.to_string(index=False))

    # 相关性分析
    from scipy.stats import pearsonr
    corr, p_value = pearsonr(df['coop'], df['activity'])

    print(f"\n📈 Correlation Analysis:")
    print(f"   Cooperativity vs Activity: r = {corr:.3f}, p = {p_value:.4f}")

    if corr > 0.7:
        print(f"   ✅ Strong positive correlation!")
        print(f"   → Cooperativity score is a good predictor of activity")

    # SAR洞察
    print("\n💡 SAR Insights:")
    for i, row in df.iterrows():
        if i == 0:
            continue
        if row['Coop_Change'] > 0.05:
            print(f"   ✅ {row['name']}: Improved cooperativity (+{row['Coop_Change']:.3f})")
            print(f"      → Activity increased by {row['Activity_Change']:.1f} units")
        elif row['Coop_Change'] < -0.05:
            print(f"   ⚠️  {row['name']}: Decreased cooperativity ({row['Coop_Change']:.3f})")
            print(f"      → Activity decreased by {abs(row['Activity_Change']):.1f} units")

    return df


def demo_virtual_screening():
    """演示5: 虚拟筛选 - 候选物排序"""
    print("\n" + "="*70)
    print("  Demo 5: Virtual Screening - Ranking Candidates")
    print("="*70)

    # 模拟100个候选化合物
    np.random.seed(42)
    n_compounds = 100

    candidates = []
    for i in range(1, n_compounds + 1):
        coop = np.random.beta(2, 5)  # 偏向较低值，模拟真实筛选
        bridging = np.random.poisson(3)
        neo = np.random.poisson(12)

        # 综合评分
        score = coop * 0.4 + (bridging / 10) * 0.3 + (neo / 50) * 0.3

        candidates.append({
            'Compound_ID': f'CMPD_{i:03d}',
            'Cooperativity': coop,
            'Bridging_Atoms': bridging,
            'Neo_epitope': neo,
            'Score': score
        })

    df = pd.DataFrame(candidates)
    df = df.sort_values('Score', ascending=False)

    print(f"\n📊 Screened {n_compounds} compounds")
    print(f"\n🏆 Top 10 Candidates:")
    print(df.head(10).to_string(index=False))

    # 统计
    top_10_pct = (df.head(10)['Score'].mean() / df['Score'].mean() - 1) * 100
    print(f"\n📈 Statistics:")
    print(f"   Top 10 average score: {df.head(10)['Score'].mean():.3f}")
    print(f"   Overall average score: {df['Score'].mean():.3f}")
    print(f"   Top 10 enrichment: +{top_10_pct:.1f}%")

    # 推荐
    print(f"\n💊 Recommendation:")
    print(f"   Prioritize synthesis of: {', '.join(df.head(5)['Compound_ID'].tolist())}")
    print(f"   Expected success rate: ~50% (based on score distribution)")

    return df


def main():
    """运行所有演示"""
    print("\n" + "="*70)
    print("  Molecular Glue Fingerprint (MGF) System Demo")
    print("  GLINT Platform - Advanced Molecular Glue Analysis")
    print("="*70)

    try:
        # Demo 1: 基础用法
        result1 = demo_basic_usage()

        # Demo 2: 比较分析
        result2 = demo_comparison()

        # Demo 3: 可视化
        result3 = demo_visualization()

        # Demo 4: SAR分析
        result4 = demo_sar_analysis()

        # Demo 5: 虚拟筛选
        result5 = demo_virtual_screening()

        print("\n" + "="*70)
        print("  All Demos Completed Successfully! ✅")
        print("="*70)
        print("\n📁 Generated Files:")
        print(f"   - {result3}")
        print("\n📚 Next Steps:")
        print("   1. Review the analysis results")
        print("   2. Try with your own molecular glue structures")
        print("   3. Integrate MGF into your drug discovery pipeline")
        print("\n🚀 Happy Molecular Glue Discovery!")

    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == '__main__':
    main()

