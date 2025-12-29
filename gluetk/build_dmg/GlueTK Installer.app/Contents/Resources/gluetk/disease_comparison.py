# -*- coding: utf-8 -*-
"""
Disease Comparison and Multi-Disease Analysis
多疾病比较分析模块

支持同时分析多个疾病，找出共同靶点，生成可视化图表

Author: Vesper
"""

from typing import Dict, List, Optional, Set
import os

try:
    import pandas as pd
except ImportError:
    print("⚠️ Warning: pandas not installed")
    pd = None

try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')  # 使用非交互式后端
except ImportError:
    print("⚠️ Warning: matplotlib not installed")
    plt = None

try:
    import seaborn as sns
except ImportError:
    print("⚠️ Warning: seaborn not installed")
    sns = None

try:
    from matplotlib_venn import venn2, venn3
    _VENN_AVAILABLE = True
except ImportError:
    print("⚠️ Warning: matplotlib-venn not installed")
    _VENN_AVAILABLE = False

from .open_targets_api import search_disease, get_disease_targets, enrich_targets_with_e3_scores
from .disease_config import (
    DEFAULT_OUTPUT_DIR,
    MAX_DISEASES_FOR_VENN,
    HEATMAP_COLOR_SCHEME,
    COMPARISON_OUTPUT_FORMATS
)


# ============================================================================
# 多疾病靶点查询
# ============================================================================

def compare_diseases(
    disease_names: List[str],
    top_n: int = 30,
    include_e3_score: bool = True,
    e3_symbol: str = "CRBN",
    use_cache: bool = True
) -> Dict[str, pd.DataFrame]:
    """
    比较多个疾病的关联靶点
    
    Args:
        disease_names: 疾病名称列表（英文）
        top_n: 每个疾病返回的靶点数
        include_e3_score: 是否包含 E3 评分
        e3_symbol: E3 连接酶符号
        use_cache: 是否使用缓存
    
    Returns:
        {disease_name: targets_df} 字典
    
    Example:
        >>> results = compare_diseases(
        ...     ["multiple myeloma", "acute myeloid leukemia"],
        ...     top_n=30
        ... )
        >>> print(results.keys())
        dict_keys(['multiple myeloma', 'acute myeloid leukemia'])
    """
    if pd is None:
        raise RuntimeError("pandas is required for disease comparison")
    
    print("\n" + "=" * 100)
    print("🔬 MULTI-DISEASE COMPARISON ANALYSIS")
    print("=" * 100)
    print(f"Diseases to compare: {len(disease_names)}")
    for i, name in enumerate(disease_names, 1):
        print(f"  {i}. {name}")
    print("=" * 100)
    
    results = {}
    
    for disease_name in disease_names:
        print(f"\n📊 Processing: {disease_name}")
        print("-" * 80)
        
        try:
            # 搜索疾病
            candidates = search_disease(disease_name, max_results=5, use_cache=use_cache)
            
            if not candidates:
                print(f"⚠️ No disease found for '{disease_name}', skipping...")
                continue
            
            # 使用第一个候选
            disease_id = candidates[0]['id']
            disease_display_name = candidates[0]['name']
            
            print(f"✅ Selected: {disease_display_name} ({disease_id})")
            
            # 获取靶点
            df = get_disease_targets(disease_id, top_n=top_n, use_cache=use_cache)
            
            if df is None or df.empty:
                print(f"⚠️ No targets found for {disease_display_name}, skipping...")
                continue
            
            # 添加 E3 评分
            if include_e3_score:
                try:
                    df = enrich_targets_with_e3_scores(df, e3_symbol=e3_symbol)
                    print(f"✅ Added E3 scores ({e3_symbol})")
                except Exception as e:
                    print(f"⚠️ E3 scoring failed: {e}")
            
            # 添加疾病信息
            df['disease_name'] = disease_display_name
            df['disease_id'] = disease_id
            
            results[disease_display_name] = df
            print(f"✅ Retrieved {len(df)} targets")
            
        except Exception as e:
            print(f"❌ Error processing '{disease_name}': {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 100)
    print(f"✅ Comparison complete: {len(results)}/{len(disease_names)} diseases processed")
    print("=" * 100 + "\n")
    
    return results


# ============================================================================
# 共同靶点分析
# ============================================================================

def find_common_targets(
    disease_targets: Dict[str, pd.DataFrame],
    min_diseases: int = 2,
    score_threshold: float = 0.0
) -> pd.DataFrame:
    """
    找出多个疾病的共同靶点
    
    Args:
        disease_targets: compare_diseases() 的输出
        min_diseases: 最少出现在几个疾病中
        score_threshold: 最小评分阈值
    
    Returns:
        共同靶点 DataFrame，包含每个疾病的评分
    
    Example:
        >>> results = compare_diseases(["disease1", "disease2"])
        >>> common = find_common_targets(results, min_diseases=2)
        >>> print(common[['symbol', 'num_diseases', 'mean_score']])
    """
    if pd is None:
        raise RuntimeError("pandas is required")
    
    if not disease_targets:
        print("⚠️ No disease data provided")
        return pd.DataFrame()
    
    print("\n🔍 Finding common targets...")
    print(f"Criteria: Present in ≥{min_diseases} diseases, Score ≥{score_threshold}")
    
    # 收集所有靶点
    all_targets = {}  # {symbol: {disease: score}}
    
    for disease_name, df in disease_targets.items():
        # 过滤低分靶点
        df_filtered = df[df['score'] >= score_threshold]
        
        for _, row in df_filtered.iterrows():
            symbol = row['symbol']
            score = row['score']
            
            if symbol not in all_targets:
                all_targets[symbol] = {}
            
            all_targets[symbol][disease_name] = {
                'score': score,
                'e3_score': row.get('e3_score', None),
                'composite_score': row.get('composite_score', None)
            }
    
    # 过滤共同靶点
    common_targets = []
    
    for symbol, disease_scores in all_targets.items():
        num_diseases = len(disease_scores)
        
        if num_diseases >= min_diseases:
            # 计算平均分数
            scores = [d['score'] for d in disease_scores.values()]
            mean_score = sum(scores) / len(scores)
            
            # E3 评分
            e3_scores = [d['e3_score'] for d in disease_scores.values() if d['e3_score'] is not None]
            mean_e3_score = sum(e3_scores) / len(e3_scores) if e3_scores else None
            
            # 综合评分
            comp_scores = [d['composite_score'] for d in disease_scores.values() if d['composite_score'] is not None]
            mean_composite = sum(comp_scores) / len(comp_scores) if comp_scores else None
            
            common_targets.append({
                'symbol': symbol,
                'num_diseases': num_diseases,
                'diseases': list(disease_scores.keys()),
                'mean_score': mean_score,
                'mean_e3_score': mean_e3_score,
                'mean_composite_score': mean_composite,
                'disease_scores': disease_scores
            })
    
    # 创建 DataFrame
    if not common_targets:
        print("⚠️ No common targets found")
        return pd.DataFrame()
    
    df_common = pd.DataFrame(common_targets)
    
    # 按出现疾病数和平均分数排序
    df_common = df_common.sort_values(
        ['num_diseases', 'mean_score'],
        ascending=[False, False]
    ).reset_index(drop=True)
    
    print(f"✅ Found {len(df_common)} common targets")
    print(f"   Present in all {len(disease_targets)} diseases: {len(df_common[df_common['num_diseases'] == len(disease_targets)])}")
    
    return df_common


# ============================================================================
# Venn 图生成
# ============================================================================

def generate_venn_diagram(
    disease_targets: Dict[str, pd.DataFrame],
    output_path: str,
    score_threshold: float = 0.0,
    title: str = None
) -> bool:
    """
    生成 Venn 图（最多 3 个疾病）
    
    Args:
        disease_targets: 疾病靶点字典
        output_path: 输出图片路径
        score_threshold: 最小评分阈值
        title: 图表标题
    
    Returns:
        是否成功生成
    """
    if not _VENN_AVAILABLE or plt is None:
        print("❌ matplotlib-venn is required for Venn diagrams")
        return False
    
    num_diseases = len(disease_targets)
    
    if num_diseases < 2:
        print("❌ Need at least 2 diseases for Venn diagram")
        return False
    
    if num_diseases > MAX_DISEASES_FOR_VENN:
        print(f"⚠️ Venn diagram supports max {MAX_DISEASES_FOR_VENN} diseases, using first {MAX_DISEASES_FOR_VENN}")
        disease_targets = dict(list(disease_targets.items())[:MAX_DISEASES_FOR_VENN])
        num_diseases = MAX_DISEASES_FOR_VENN
    
    print(f"\n📊 Generating Venn diagram for {num_diseases} diseases...")
    
    # 提取靶点集合
    disease_sets = {}
    for disease_name, df in disease_targets.items():
        df_filtered = df[df['score'] >= score_threshold]
        disease_sets[disease_name] = set(df_filtered['symbol'].tolist())
    
    # 创建图表
    fig, ax = plt.subplots(figsize=(10, 8))
    
    disease_names = list(disease_sets.keys())
    sets = [disease_sets[name] for name in disease_names]
    
    try:
        if num_diseases == 2:
            venn = venn2(sets, set_labels=disease_names, ax=ax)
        elif num_diseases == 3:
            venn = venn3(sets, set_labels=disease_names, ax=ax)
        
        # 设置标题
        if title is None:
            title = f"Common Targets Across {num_diseases} Diseases"
        ax.set_title(title, fontsize=14, fontweight='bold')
        
        # 保存
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✅ Venn diagram saved to: {output_path}")
        return True
        
    except Exception as e:
        print(f"❌ Failed to generate Venn diagram: {e}")
        plt.close()
        return False


# ============================================================================
# 热图生成
# ============================================================================

def generate_heatmap(
    disease_targets: Dict[str, pd.DataFrame],
    output_path: str,
    top_n: int = 20,
    score_type: str = "composite",
    title: str = None
) -> bool:
    """
    生成疾病-靶点关联热图
    
    Args:
        disease_targets: 疾病靶点字典
        output_path: 输出图片路径
        top_n: 显示的靶点数量
        score_type: 评分类型 ("disease", "e3", "composite")
        title: 图表标题
    
    Returns:
        是否成功生成
    """
    if pd is None or plt is None or sns is None:
        print("❌ pandas, matplotlib, and seaborn are required for heatmaps")
        return False
    
    print(f"\n📊 Generating heatmap (top {top_n} targets, {score_type} scores)...")
    
    # 找出共同靶点
    common_df = find_common_targets(disease_targets, min_diseases=1)
    
    if common_df.empty:
        print("❌ No targets found for heatmap")
        return False
    
    # 选择 top N 靶点
    if score_type == "composite" and 'mean_composite_score' in common_df.columns:
        sort_col = 'mean_composite_score'
    elif score_type == "e3" and 'mean_e3_score' in common_df.columns:
        sort_col = 'mean_e3_score'
    else:
        sort_col = 'mean_score'
    
    top_targets = common_df.nlargest(top_n, sort_col)['symbol'].tolist()
    
    # 构建矩阵
    matrix_data = []
    disease_names = list(disease_targets.keys())
    
    for symbol in top_targets:
        row = []
        for disease_name in disease_names:
            df = disease_targets[disease_name]
            target_row = df[df['symbol'] == symbol]
            
            if not target_row.empty:
                if score_type == "composite" and 'composite_score' in target_row.columns:
                    score = target_row.iloc[0]['composite_score']
                elif score_type == "e3" and 'e3_score' in target_row.columns:
                    score = target_row.iloc[0]['e3_score']
                else:
                    score = target_row.iloc[0]['score']
            else:
                score = 0.0
            
            row.append(score)
        
        matrix_data.append(row)
    
    # 创建 DataFrame
    heatmap_df = pd.DataFrame(
        matrix_data,
        index=top_targets,
        columns=disease_names
    )
    
    # 绘制热图
    try:
        fig, ax = plt.subplots(figsize=(max(10, len(disease_names) * 2), max(8, top_n * 0.4)))
        
        sns.heatmap(
            heatmap_df,
            annot=True,
            fmt='.2f',
            cmap=HEATMAP_COLOR_SCHEME,
            cbar_kws={'label': f'{score_type.capitalize()} Score'},
            ax=ax,
            vmin=0,
            vmax=1
        )
        
        # 设置标题
        if title is None:
            title = f"Disease-Target Association Heatmap ({score_type.capitalize()} Scores)"
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        
        ax.set_xlabel("Disease", fontsize=12)
        ax.set_ylabel("Target Gene", fontsize=12)
        
        # 旋转标签
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        
        # 保存
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✅ Heatmap saved to: {output_path}")
        return True
        
    except Exception as e:
        print(f"❌ Failed to generate heatmap: {e}")
        import traceback
        traceback.print_exc()
        plt.close()
        return False


# ============================================================================
# 比较报告导出
# ============================================================================

def export_comparison_report(
    disease_targets: Dict[str, pd.DataFrame],
    output_dir: str = None,
    generate_plots: bool = True,
    top_n_heatmap: int = 20
) -> Dict[str, str]:
    """
    导出完整的比较报告
    
    Args:
        disease_targets: 疾病靶点字典
        output_dir: 输出目录
        generate_plots: 是否生成图表
        top_n_heatmap: 热图显示的靶点数
    
    Returns:
        生成的文件路径字典
    
    Output:
        - comparison_summary.csv: 汇总表
        - common_targets.csv: 共同靶点列表
        - venn_diagram.png: Venn 图（如果适用）
        - heatmap.png: 热图
    """
    if output_dir is None:
        output_dir = os.path.join(DEFAULT_OUTPUT_DIR, "comparison")
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n" + "=" * 100)
    print("📝 EXPORTING COMPARISON REPORT")
    print("=" * 100)
    print(f"Output directory: {output_dir}")
    
    generated_files = {}
    
    # 1. 保存共同靶点
    common_df = find_common_targets(disease_targets, min_diseases=1)
    
    if not common_df.empty:
        common_path = os.path.join(output_dir, "common_targets.csv")
        
        # 简化输出（移除复杂的嵌套字典）
        export_df = common_df[['symbol', 'num_diseases', 'mean_score']].copy()
        if 'mean_e3_score' in common_df.columns:
            export_df['mean_e3_score'] = common_df['mean_e3_score']
        if 'mean_composite_score' in common_df.columns:
            export_df['mean_composite_score'] = common_df['mean_composite_score']
        
        export_df.to_csv(common_path, index=False)
        generated_files['common_targets'] = common_path
        print(f"✅ Common targets saved: {common_path}")
    
    # 2. 保存各疾病的靶点列表
    for disease_name, df in disease_targets.items():
        safe_name = "".join(c if c.isalnum() or c in " _-" else "_" for c in disease_name)
        disease_path = os.path.join(output_dir, f"targets_{safe_name}.csv")
        df.to_csv(disease_path, index=False)
        generated_files[f'targets_{safe_name}'] = disease_path
    
    print(f"✅ Individual disease targets saved ({len(disease_targets)} files)")
    
    # 3. 生成图表
    if generate_plots:
        # Venn 图
        if len(disease_targets) >= 2 and len(disease_targets) <= MAX_DISEASES_FOR_VENN:
            venn_path = os.path.join(output_dir, "venn_diagram.png")
            if generate_venn_diagram(disease_targets, venn_path):
                generated_files['venn_diagram'] = venn_path
        
        # 热图
        heatmap_path = os.path.join(output_dir, "heatmap.png")
        if generate_heatmap(disease_targets, heatmap_path, top_n=top_n_heatmap):
            generated_files['heatmap'] = heatmap_path
    
    print("=" * 100)
    print(f"✅ Report export complete: {len(generated_files)} files generated")
    print("=" * 100 + "\n")
    
    return generated_files


# ============================================================================
# 打印比较摘要
# ============================================================================

def print_comparison_summary(disease_targets: Dict[str, pd.DataFrame]) -> None:
    """
    打印多疾病比较摘要
    """
    if not disease_targets:
        print("No data to summarize")
        return
    
    print("\n" + "=" * 100)
    print("📊 MULTI-DISEASE COMPARISON SUMMARY")
    print("=" * 100)
    
    # 疾病列表
    print(f"\nDiseases analyzed: {len(disease_targets)}")
    for i, (disease_name, df) in enumerate(disease_targets.items(), 1):
        print(f"  {i}. {disease_name}: {len(df)} targets")
    
    # 共同靶点统计
    common_df = find_common_targets(disease_targets, min_diseases=1)
    
    if not common_df.empty:
        print(f"\nCommon targets statistics:")
        for n in range(len(disease_targets), 1, -1):
            count = len(common_df[common_df['num_diseases'] == n])
            if count > 0:
                print(f"  Present in {n}/{len(disease_targets)} diseases: {count} targets")
        
        # Top 共同靶点
        print(f"\nTop 10 common targets (by mean score):")
        top_common = common_df.head(10)
        for i, row in top_common.iterrows():
            diseases_str = ", ".join(row['diseases'][:3])
            if len(row['diseases']) > 3:
                diseases_str += f", +{len(row['diseases']) - 3} more"
            print(f"  {i+1}. {row['symbol']}: {row['num_diseases']} diseases, mean score {row['mean_score']:.3f}")
            print(f"      ({diseases_str})")
    
    print("=" * 100 + "\n")
