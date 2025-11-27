# -*- coding: utf-8 -*-
"""
Disease Target Analysis - Configuration
配置文件：Open Targets API 和数据路径设置
"""

import os

# ============================================================================
# API 配置
# ============================================================================

# Open Targets Platform GraphQL API 端点
OPEN_TARGETS_GRAPHQL_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"

# HTTP 请求配置
REQUEST_TIMEOUT = 30  # 秒
MAX_RETRIES = 3       # 最大重试次数
RETRY_DELAY = 1       # 重试延迟（秒）

# ============================================================================
# 数据路径配置
# ============================================================================

# 默认输出目录（相对于当前工作目录）
DEFAULT_OUTPUT_DIR = "./disease_data"

# 残基评分文件目录
RESIDUE_SCORES_DIR = os.path.join(DEFAULT_OUTPUT_DIR, "residue_scores")

# 基因级评分文件（可选，用于 E3_score 等扩展功能）
GENE_LEVEL_SCORES_CSV = os.path.join(DEFAULT_OUTPUT_DIR, "gene_level_scores.csv")

# 查询结果缓存目录
CACHE_DIR = os.path.join(DEFAULT_OUTPUT_DIR, "cache")

# ============================================================================
# 可视化配置
# ============================================================================

# 热点残基阈值（SurfHot_score 默认阈值）
DEFAULT_HOTSPOT_THRESHOLD = 0.7

# 热点残基颜色
HOTSPOT_COLOR = "red"

# 热点残基表示方式
HOTSPOT_REPRESENTATION = "spheres"

# 非热点残基颜色（背景色）
BACKGROUND_COLOR = "gray80"

# ============================================================================
# V2 功能配置：E3 兼容度评分
# ============================================================================

# 支持的 E3 连接酶列表
E3_LIGASES = ["CRBN", "VHL", "MDM2", "XIAP"]

# 默认 E3 连接酶
DEFAULT_E3_LIGASE = "CRBN"

# 综合评分权重配置
COMPOSITE_SCORE_WEIGHTS = {
    "disease": 0.6,  # 疾病关联权重
    "e3": 0.4        # E3 兼容度权重
}

# 已知 E3-靶点兼容性（基于文献，简化版）
# 更完整的数据在 e3_compatibility.py 中
KNOWN_E3_COMPATIBILITY = {
    "CRBN": {
        "IKZF1": 0.95,
        "IKZF3": 0.92,
        "GSPT1": 0.88,
        "CSNK1A1": 0.85,
        "ZFP91": 0.80,
    },
    "VHL": {
        "HIF1A": 0.95,
        "HIF2A": 0.93,
        "BRD4": 0.80,
    },
    "MDM2": {
        "TP53": 0.95,
        "MDM4": 0.85,
    }
}

# E3 评分默认值（未知靶点）
DEFAULT_E3_SCORE = 0.5

# ============================================================================
# 多疾病比较配置
# ============================================================================

# Venn 图最多支持的疾病数量
MAX_DISEASES_FOR_VENN = 3

# 热图颜色方案
HEATMAP_COLOR_SCHEME = "YlOrRd"

# 比较报告输出格式
COMPARISON_OUTPUT_FORMATS = ["csv", "html", "png"]

# ============================================================================
# GraphQL 查询模板
# ============================================================================

# 疾病搜索查询
DISEASE_SEARCH_QUERY = """
query SearchDisease($queryString: String!, $size: Int!) {
  search(queryString: $queryString, entityNames: ["disease"], page: {index: 0, size: $size}) {
    hits {
      id
      name
      description
      entity
    }
  }
}
"""

# 疾病关联靶点查询
DISEASE_TARGETS_QUERY = """
query DiseaseTargets($efoId: String!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets(page: {index: 0, size: $size}) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
          proteinIds {
            id
            source
          }
        }
        score
      }
    }
  }
}
"""

# ============================================================================
# 文件命名约定
# ============================================================================

def get_targets_csv_path(disease_id: str, output_dir: str = None) -> str:
    """
    获取靶点列表 CSV 文件路径
    
    Args:
        disease_id: 疾病 EFO ID（如 EFO_0001378）
        output_dir: 输出目录（默认使用 DEFAULT_OUTPUT_DIR）
    
    Returns:
        CSV 文件完整路径
    """
    if output_dir is None:
        output_dir = DEFAULT_OUTPUT_DIR
    
    os.makedirs(output_dir, exist_ok=True)
    return os.path.join(output_dir, f"targets_{disease_id}.csv")


def get_residue_scores_path(gene_symbol: str, scores_dir: str = None) -> str:
    """
    获取残基评分文件路径
    
    Args:
        gene_symbol: 基因符号（如 IKZF1）
        scores_dir: 评分文件目录（默认使用 RESIDUE_SCORES_DIR）
    
    Returns:
        残基评分 CSV 文件完整路径
    """
    if scores_dir is None:
        scores_dir = RESIDUE_SCORES_DIR
    
    return os.path.join(scores_dir, f"residue_scores_{gene_symbol}.csv")


def get_cache_path(cache_key: str, cache_dir: str = None) -> str:
    """
    获取缓存文件路径
    
    Args:
        cache_key: 缓存键（如 disease_search_multiple_myeloma）
        cache_dir: 缓存目录（默认使用 CACHE_DIR）
    
    Returns:
        缓存 JSON 文件完整路径
    """
    if cache_dir is None:
        cache_dir = CACHE_DIR
    
    os.makedirs(cache_dir, exist_ok=True)
    # 清理缓存键中的特殊字符
    safe_key = "".join(c if c.isalnum() or c in "_-" else "_" for c in cache_key)
    return os.path.join(cache_dir, f"{safe_key}.json")
