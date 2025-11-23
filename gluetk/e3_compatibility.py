# -*- coding: utf-8 -*-
"""
E3 Ligase Compatibility Scoring
E3 连接酶兼容度评分模块

计算靶点蛋白与 E3 连接酶（CRBN、VHL 等）的兼容度评分
基于已知的分子胶数据和文献报道

Author: Vesper
"""

from typing import Dict, List, Optional
import os

try:
    import pandas as pd
except ImportError:
    print("⚠️ Warning: pandas not installed")
    pd = None

from .disease_config import (
    KNOWN_E3_COMPATIBILITY,
    DEFAULT_E3_LIGASE,
    E3_LIGASES,
    GENE_LEVEL_SCORES_CSV
)


# ============================================================================
# 已知 E3-靶点兼容性数据（基于文献）
# ============================================================================

# 扩展的已知兼容性数据
EXTENDED_E3_COMPATIBILITY = {
    "CRBN": {
        # 高兼容性靶点（已验证的分子胶靶点）
        "IKZF1": 0.95,    # Lenalidomide, Pomalidomide
        "IKZF3": 0.92,    # Lenalidomide, Pomalidomide
        "GSPT1": 0.88,    # CC-885, CC-90009
        "CSNK1A1": 0.85,  # Lenalidomide (CK1α)
        "ZFP91": 0.80,    # Lenalidomide
        "SALL4": 0.78,    # Thalidomide derivatives
        "IKZF2": 0.75,    # Helios
        "ZNF692": 0.72,   # Thalidomide
        "ZNF276": 0.70,   # Thalidomide
        
        # 中等兼容性（潜在靶点）
        "MEIS2": 0.65,
        "ZMYM2": 0.62,
        "PLZF": 0.60,
    },
    
    "VHL": {
        # 高兼容性靶点
        "HIF1A": 0.95,    # VHL-HIF1α 经典靶点
        "HIF2A": 0.93,    # EPAS1
        "TERT": 0.85,     # Telomerase
        
        # PROTAC 靶点（通过 VHL）
        "BRD4": 0.80,     # ARV-825, ARV-771
        "BCR-ABL": 0.78,  # DT2216
        "STAT3": 0.75,
    },
    
    "MDM2": {
        # 高兼容性靶点
        "TP53": 0.95,     # MDM2-p53 经典靶点
        "MDM4": 0.85,
        
        # PROTAC 靶点
        "BRD4": 0.75,
        "CDK9": 0.70,
    },
    
    "XIAP": {
        # 高兼容性靶点
        "CASP3": 0.90,    # Caspase-3
        "CASP7": 0.88,
        "CASP9": 0.85,
        "SMAC": 0.82,
    }
}


# ============================================================================
# E3 兼容度评分函数
# ============================================================================

def get_known_e3_score(target_symbol: str, e3_symbol: str = "CRBN") -> Optional[float]:
    """
    从已知兼容性表中获取评分
    
    Args:
        target_symbol: 靶点基因符号（如 "IKZF1"）
        e3_symbol: E3 连接酶符号（如 "CRBN"）
    
    Returns:
        兼容度评分 (0-1)，如果未知则返回 None
    """
    # 优先使用扩展表
    if e3_symbol in EXTENDED_E3_COMPATIBILITY:
        if target_symbol in EXTENDED_E3_COMPATIBILITY[e3_symbol]:
            return EXTENDED_E3_COMPATIBILITY[e3_symbol][target_symbol]
    
    # 回退到配置文件中的表
    if e3_symbol in KNOWN_E3_COMPATIBILITY:
        if target_symbol in KNOWN_E3_COMPATIBILITY[e3_symbol]:
            return KNOWN_E3_COMPATIBILITY[e3_symbol][target_symbol]
    
    return None


def calculate_e3_compatibility(
    target_symbol: str,
    e3_symbol: str = None,
    method: str = "known_targets",
    default_score: float = 0.5
) -> float:
    """
    计算靶点与 E3 连接酶的兼容度评分
    
    Args:
        target_symbol: 靶点基因符号
        e3_symbol: E3 连接酶符号（默认使用 DEFAULT_E3_LIGASE）
        method: 计算方法
            - "known_targets": 使用已知兼容性表（默认）
            - "coexpression": 共表达分析（未实现）
            - "literature": 文献挖掘（未实现）
        default_score: 未知靶点的默认评分（默认 0.5）
    
    Returns:
        兼容度评分 (0-1)
    
    Example:
        >>> score = calculate_e3_compatibility("IKZF1", "CRBN")
        >>> print(score)  # 0.95
    """
    if e3_symbol is None:
        e3_symbol = DEFAULT_E3_LIGASE
    
    # 验证 E3 连接酶
    if e3_symbol not in E3_LIGASES:
        print(f"⚠️ Warning: Unknown E3 ligase '{e3_symbol}', using default")
        e3_symbol = DEFAULT_E3_LIGASE
    
    # 方法 1: 已知靶点表
    if method == "known_targets":
        score = get_known_e3_score(target_symbol, e3_symbol)
        
        if score is not None:
            return score
        else:
            # 未知靶点，返回默认分数
            return default_score
    
    # 方法 2: 共表达分析（未实现）
    elif method == "coexpression":
        print(f"⚠️ Method 'coexpression' not implemented yet, using default score")
        return default_score
    
    # 方法 3: 文献挖掘（未实现）
    elif method == "literature":
        print(f"⚠️ Method 'literature' not implemented yet, using default score")
        return default_score
    
    else:
        raise ValueError(f"Unknown method: {method}")


def batch_calculate_e3_scores(
    targets_df: "pd.DataFrame",
    e3_symbol: str = None,
    method: str = "known_targets",
    default_score: float = 0.5
) -> "pd.DataFrame":
    """
    批量计算 E3 兼容度评分
    
    Args:
        targets_df: 靶点列表 DataFrame（来自 get_disease_targets）
        e3_symbol: E3 连接酶符号
        method: 计算方法
        default_score: 未知靶点的默认评分
    
    Returns:
        原 DataFrame + e3_score 列
    
    Example:
        >>> from gluetk.open_targets_api import get_disease_targets
        >>> df = get_disease_targets("EFO_0001378", top_n=10)
        >>> df_with_e3 = batch_calculate_e3_scores(df, "CRBN")
        >>> print(df_with_e3[['symbol', 'score', 'e3_score']])
    """
    if pd is None:
        raise RuntimeError("pandas is required for batch processing")
    
    if e3_symbol is None:
        e3_symbol = DEFAULT_E3_LIGASE
    
    # 复制 DataFrame
    result_df = targets_df.copy()
    
    # 计算每个靶点的 E3 评分
    e3_scores = []
    for symbol in result_df['symbol']:
        score = calculate_e3_compatibility(
            symbol,
            e3_symbol=e3_symbol,
            method=method,
            default_score=default_score
        )
        e3_scores.append(score)
    
    # 添加 E3 评分列
    result_df['e3_score'] = e3_scores
    
    # 添加 E3 连接酶信息
    result_df['e3_ligase'] = e3_symbol
    
    return result_df


# ============================================================================
# 基因级评分文件管理
# ============================================================================

def save_gene_level_scores(
    df: "pd.DataFrame",
    output_path: str = None,
    include_composite: bool = False
) -> None:
    """
    保存基因级评分到 CSV 文件
    
    Args:
        df: 包含评分的 DataFrame
        output_path: 输出路径（默认使用 GENE_LEVEL_SCORES_CSV）
        include_composite: 是否包含综合评分列
    
    Output format:
        symbol, disease_score, e3_score, [composite_score], e3_ligase
    """
    if pd is None:
        raise RuntimeError("pandas is required")
    
    if output_path is None:
        output_path = GENE_LEVEL_SCORES_CSV
    
    # 确保目录存在
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # 选择要保存的列
    columns = ['symbol', 'score', 'e3_score', 'e3_ligase']
    
    if include_composite and 'composite_score' in df.columns:
        columns.insert(3, 'composite_score')
    
    # 重命名列以便更清晰
    save_df = df[columns].copy()
    save_df.columns = ['symbol', 'disease_score', 'e3_score'] + \
                      (['composite_score'] if include_composite and 'composite_score' in df.columns else []) + \
                      ['e3_ligase']
    
    # 保存
    save_df.to_csv(output_path, index=False)
    print(f"💾 Saved gene-level scores to: {output_path}")


def load_gene_level_scores(input_path: str = None) -> Optional["pd.DataFrame"]:
    """
    从 CSV 文件加载基因级评分
    
    Args:
        input_path: 输入路径（默认使用 GENE_LEVEL_SCORES_CSV）
    
    Returns:
        DataFrame 或 None（如果文件不存在）
    """
    if pd is None:
        raise RuntimeError("pandas is required")
    
    if input_path is None:
        input_path = GENE_LEVEL_SCORES_CSV
    
    if not os.path.exists(input_path):
        return None
    
    try:
        df = pd.read_csv(input_path)
        print(f"📂 Loaded gene-level scores from: {input_path}")
        return df
    except Exception as e:
        print(f"❌ Failed to load gene-level scores: {e}")
        return None


# ============================================================================
# 统计与分析
# ============================================================================

def get_e3_compatibility_stats(e3_symbol: str = "CRBN") -> Dict:
    """
    获取 E3 兼容性数据的统计信息
    
    Args:
        e3_symbol: E3 连接酶符号
    
    Returns:
        统计信息字典
    """
    if e3_symbol not in EXTENDED_E3_COMPATIBILITY:
        return {"error": f"No data for E3 ligase: {e3_symbol}"}
    
    targets = EXTENDED_E3_COMPATIBILITY[e3_symbol]
    scores = list(targets.values())
    
    return {
        "e3_ligase": e3_symbol,
        "num_known_targets": len(targets),
        "known_targets": list(targets.keys()),
        "score_range": (min(scores), max(scores)),
        "mean_score": sum(scores) / len(scores),
        "high_compatibility": [k for k, v in targets.items() if v >= 0.8],
        "medium_compatibility": [k for k, v in targets.items() if 0.6 <= v < 0.8],
    }


def print_e3_stats(e3_symbol: str = "CRBN") -> None:
    """
    打印 E3 兼容性统计信息
    """
    stats = get_e3_compatibility_stats(e3_symbol)
    
    if "error" in stats:
        print(f"❌ {stats['error']}")
        return
    
    print("\n" + "=" * 80)
    print(f"E3 LIGASE COMPATIBILITY STATISTICS: {e3_symbol}")
    print("=" * 80)
    print(f"Known Targets: {stats['num_known_targets']}")
    print(f"Score Range: {stats['score_range'][0]:.2f} - {stats['score_range'][1]:.2f}")
    print(f"Mean Score: {stats['mean_score']:.2f}")
    print(f"\nHigh Compatibility (≥0.8): {len(stats['high_compatibility'])} targets")
    print(f"  {', '.join(stats['high_compatibility'][:10])}")
    if len(stats['high_compatibility']) > 10:
        print(f"  ... and {len(stats['high_compatibility']) - 10} more")
    print(f"\nMedium Compatibility (0.6-0.8): {len(stats['medium_compatibility'])} targets")
    print(f"  {', '.join(stats['medium_compatibility'][:10])}")
    print("=" * 80 + "\n")


# ============================================================================
# 测试函数
# ============================================================================

def test_e3_compatibility():
    """
    测试 E3 兼容度评分功能
    """
    print("\n🧪 Testing E3 Compatibility Module\n")
    
    # 测试 1: 已知靶点
    print("Test 1: Known Target (IKZF1 + CRBN)")
    print("-" * 50)
    score = calculate_e3_compatibility("IKZF1", "CRBN")
    print(f"Score: {score:.2f}")
    assert score > 0.8, "IKZF1 should have high CRBN compatibility"
    print("✅ Pass\n")
    
    # 测试 2: 未知靶点
    print("Test 2: Unknown Target (RANDOM_GENE + CRBN)")
    print("-" * 50)
    score = calculate_e3_compatibility("RANDOM_GENE", "CRBN")
    print(f"Score: {score:.2f} (default)")
    assert 0 <= score <= 1, "Score should be in valid range"
    print("✅ Pass\n")
    
    # 测试 3: 不同 E3 连接酶
    print("Test 3: Different E3 Ligases")
    print("-" * 50)
    for e3 in ["CRBN", "VHL", "MDM2"]:
        stats = get_e3_compatibility_stats(e3)
        print(f"{e3}: {stats['num_known_targets']} known targets")
    print("✅ Pass\n")
    
    # 测试 4: 统计信息
    print("Test 4: E3 Statistics")
    print("-" * 50)
    print_e3_stats("CRBN")
    print("✅ Pass\n")
    
    print("🎉 All tests passed!")


if __name__ == "__main__":
    # 直接运行此文件时执行测试
    test_e3_compatibility()
