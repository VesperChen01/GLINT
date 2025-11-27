# -*- coding: utf-8 -*-
"""
Open Targets Platform API Wrapper
封装 Open Targets GraphQL API，用于疾病-靶点关联查询

Author: Vesper
"""

import json
import time
from typing import List, Dict, Optional
import os

try:
    import requests
except ImportError:
    print("⚠️ Warning: 'requests' library not found. Please install: pip install requests")
    requests = None

try:
    import pandas as pd
except ImportError:
    print("⚠️ Warning: 'pandas' library not found. Please install: pip install pandas")
    pd = None

from .disease_config import (
    OPEN_TARGETS_GRAPHQL_ENDPOINT,
    REQUEST_TIMEOUT,
    MAX_RETRIES,
    RETRY_DELAY,
    DISEASE_SEARCH_QUERY,
    DISEASE_TARGETS_QUERY,
    get_cache_path
)


# ============================================================================
# 核心 API 调用函数
# ============================================================================

def _make_graphql_request(query: str, variables: dict, use_cache: bool = True) -> dict:
    """
    发送 GraphQL 请求到 Open Targets Platform API
    
    Args:
        query: GraphQL 查询字符串
        variables: 查询变量字典
        use_cache: 是否使用缓存（默认 True）
    
    Returns:
        API 响应的 JSON 数据
    
    Raises:
        RuntimeError: 当请求失败或达到最大重试次数时
    """
    if requests is None:
        raise RuntimeError("requests library is not installed. Please run: pip install requests")
    
    # 生成缓存键
    cache_key = f"{hash(query)}_{hash(json.dumps(variables, sort_keys=True))}"
    cache_path = get_cache_path(cache_key)
    
    # 尝试从缓存读取
    if use_cache and os.path.exists(cache_path):
        try:
            with open(cache_path, 'r', encoding='utf-8') as f:
                print(f"📦 Using cached data from {cache_path}")
                return json.load(f)
        except Exception as e:
            print(f"⚠️ Cache read failed: {e}, fetching from API...")
    
    # 准备请求
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    
    payload = {
        "query": query,
        "variables": variables
    }
    
    # 重试机制
    for attempt in range(MAX_RETRIES):
        try:
            print(f"🌐 Sending request to Open Targets API (attempt {attempt + 1}/{MAX_RETRIES})...")
            
            response = requests.post(
                OPEN_TARGETS_GRAPHQL_ENDPOINT,
                json=payload,
                headers=headers,
                timeout=REQUEST_TIMEOUT
            )
            
            # 检查 HTTP 状态码
            response.raise_for_status()
            
            # 解析 JSON
            data = response.json()
            
            # 检查 GraphQL 错误
            if "errors" in data:
                error_msg = "; ".join([e.get("message", str(e)) for e in data["errors"]])
                raise RuntimeError(f"GraphQL errors: {error_msg}")
            
            # 保存到缓存
            if use_cache:
                try:
                    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
                    with open(cache_path, 'w', encoding='utf-8') as f:
                        json.dump(data, f, indent=2)
                    print(f"💾 Cached response to {cache_path}")
                except Exception as e:
                    print(f"⚠️ Cache write failed: {e}")
            
            print("✅ Request successful")
            return data
            
        except requests.exceptions.Timeout:
            print(f"⏱️ Request timeout (attempt {attempt + 1}/{MAX_RETRIES})")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY * (attempt + 1))  # 指数退避
            else:
                raise RuntimeError(f"Request timeout after {MAX_RETRIES} attempts")
        
        except requests.exceptions.RequestException as e:
            print(f"❌ Request failed: {e} (attempt {attempt + 1}/{MAX_RETRIES})")
            if attempt < MAX_RETRIES - 1:
                time.sleep(RETRY_DELAY * (attempt + 1))
            else:
                raise RuntimeError(f"Request failed after {MAX_RETRIES} attempts: {e}")
        
        except Exception as e:
            print(f"❌ Unexpected error: {e}")
            raise


# ============================================================================
# 疾病搜索
# ============================================================================

def search_disease(disease_name: str, max_results: int = 10, use_cache: bool = True) -> List[Dict]:
    """
    根据疾病名称搜索候选疾病
    
    Args:
        disease_name: 疾病名称（英文，如 "multiple myeloma"）
        max_results: 返回的最大结果数（默认 10）
        use_cache: 是否使用缓存（默认 True）
    
    Returns:
        候选疾病列表，每个元素包含：
        - id: 疾病 EFO ID（如 "EFO_0001378"）
        - name: 疾病名称
        - description: 疾病描述
        - entity: 实体类型（应为 "disease"）
    
    Example:
        >>> results = search_disease("multiple myeloma")
        >>> print(results[0])
        {
            'id': 'EFO_0001378',
            'name': 'multiple myeloma',
            'description': 'A plasma cell neoplasm...',
            'entity': 'disease'
        }
    """
    print(f"\n🔍 Searching for disease: '{disease_name}'")
    
    variables = {
        "queryString": disease_name,
        "size": max_results
    }
    
    try:
        response = _make_graphql_request(DISEASE_SEARCH_QUERY, variables, use_cache)
        
        # 提取搜索结果
        hits = response.get("data", {}).get("search", {}).get("hits", [])
        
        if not hits:
            print(f"⚠️ No diseases found for '{disease_name}'")
            return []
        
        print(f"✅ Found {len(hits)} candidate disease(s)")
        
        # 格式化结果
        results = []
        for hit in hits:
            results.append({
                "id": hit.get("id", ""),
                "name": hit.get("name", ""),
                "description": hit.get("description", ""),
                "entity": hit.get("entity", "")
            })
        
        return results
    
    except Exception as e:
        print(f"❌ Disease search failed: {e}")
        raise


# ============================================================================
# 疾病-靶点关联查询
# ============================================================================

def get_disease_targets(disease_id: str, top_n: int = 30, use_cache: bool = True) -> "Optional[pd.DataFrame]":
    """
    获取指定疾病的关联靶点列表
    
    Args:
        disease_id: 疾病 EFO ID（如 "EFO_0001378"）
        top_n: 返回的靶点数量（默认 30）
        use_cache: 是否使用缓存（默认 True）
    
    Returns:
        pandas DataFrame，包含以下列：
        - symbol: 基因符号（如 "IKZF1"）
        - name: 基因全名
        - score: 疾病-靶点关联分数（0-1）
        - uniprot_id: UniProt ID（如果可用）
        - ensembl_id: Ensembl ID
        
        如果查询失败或无结果，返回 None
    
    Example:
        >>> df = get_disease_targets("EFO_0001378", top_n=10)
        >>> print(df.head())
           symbol                          name     score uniprot_id ensembl_id
        0  IKZF1   IKAROS family zinc finger 1  0.82      Q13422     ENSG00000185811
        1  GSPT1   G1 to S phase transition 1  0.78      P15170     ENSG00000103342
    """
    if pd is None:
        raise RuntimeError("pandas library is not installed. Please run: pip install pandas")
    
    print(f"\n🎯 Fetching targets for disease: {disease_id}")
    
    variables = {
        "efoId": disease_id,
        "size": top_n
    }
    
    try:
        response = _make_graphql_request(DISEASE_TARGETS_QUERY, variables, use_cache)
        
        # 提取疾病信息
        disease_data = response.get("data", {}).get("disease", {})
        
        if not disease_data:
            print(f"⚠️ No data found for disease ID: {disease_id}")
            return None
        
        disease_name = disease_data.get("name", "Unknown")
        print(f"📋 Disease: {disease_name}")
        
        # 提取靶点列表
        associated_targets = disease_data.get("associatedTargets", {})
        total_count = associated_targets.get("count", 0)
        rows = associated_targets.get("rows", [])
        
        if not rows:
            print(f"⚠️ No associated targets found for {disease_name}")
            return None
        
        print(f"✅ Found {total_count} total targets, returning top {len(rows)}")
        
        # 解析靶点数据
        targets_list = []
        for row in rows:
            target = row.get("target", {})
            score = row.get("score", 0.0)
            
            # 提取 UniProt ID
            uniprot_id = None
            protein_ids = target.get("proteinIds", [])
            for pid in protein_ids:
                if pid.get("source", "").lower() == "uniprot_swissprot":
                    uniprot_id = pid.get("id")
                    break
            
            targets_list.append({
                "symbol": target.get("approvedSymbol", ""),
                "name": target.get("approvedName", ""),
                "score": score,
                "uniprot_id": uniprot_id or "",
                "ensembl_id": target.get("id", "")
            })
        
        # 创建 DataFrame
        df = pd.DataFrame(targets_list)
        
        # 按分数降序排序
        df = df.sort_values("score", ascending=False).reset_index(drop=True)
        
        return df
    
    except Exception as e:
        print(f"❌ Target query failed: {e}")
        raise


# ============================================================================
# 辅助函数
# ============================================================================

def print_disease_candidates(candidates: List[Dict]) -> None:
    """
    打印疾病候选列表（格式化输出）
    
    Args:
        candidates: search_disease() 返回的候选疾病列表
    """
    if not candidates:
        print("No candidates to display.")
        return
    
    print("\n" + "=" * 80)
    print("DISEASE CANDIDATES")
    print("=" * 80)
    
    for i, disease in enumerate(candidates, 1):
        print(f"\n[{i}] {disease['name']}")
        print(f"    ID: {disease['id']}")
        if disease.get('description'):
            desc = disease['description'][:150] + "..." if len(disease['description']) > 150 else disease['description']
            print(f"    Description: {desc}")
    
    print("\n" + "=" * 80)


def print_targets_table(df: "pd.DataFrame", top_n: int = 20) -> None:
    """
    打印靶点列表表格（格式化输出）
    
    Args:
        df: get_disease_targets() 返回的 DataFrame
        top_n: 显示的行数（默认 20）
    """
    if df is None or df.empty:
        print("No targets to display.")
        return
    
    print("\n" + "=" * 100)
    print("DISEASE-ASSOCIATED TARGETS")
    print("=" * 100)
    
    # 显示前 N 行
    display_df = df.head(top_n).copy()
    
    # 格式化分数
    display_df['score'] = display_df['score'].apply(lambda x: f"{x:.4f}")
    
    # 截断长名称
    display_df['name'] = display_df['name'].apply(lambda x: x[:40] + "..." if len(x) > 40 else x)
    
    # 打印表格
    print(display_df.to_string(index=True))
    
    if len(df) > top_n:
        print(f"\n... and {len(df) - top_n} more targets")
    
    print("=" * 100 + "\n")


# ============================================================================
# 测试函数
# ============================================================================

def test_api():
    """
    测试 API 功能（用于开发调试）
    """
    print("🧪 Testing Open Targets API Wrapper\n")
    
    # 测试 1: 疾病搜索
    print("Test 1: Disease Search")
    print("-" * 50)
    candidates = search_disease("multiple myeloma", max_results=5)
    print_disease_candidates(candidates)
    
    if not candidates:
        print("❌ Test 1 failed: No candidates found")
        return
    
    # 测试 2: 靶点查询
    print("\nTest 2: Target Query")
    print("-" * 50)
    disease_id = candidates[0]['id']
    df = get_disease_targets(disease_id, top_n=10)
    print_targets_table(df, top_n=10)
    
    if df is not None and not df.empty:
        print("✅ All tests passed!")
    else:
        print("❌ Test 2 failed: No targets found")


if __name__ == "__main__":
    # 直接运行此文件时执行测试
    test_api()


# ============================================================================
# V2 功能：综合评分系统
# ============================================================================

def calculate_composite_score(
    disease_score: float,
    e3_score: float,
    weights: tuple = None
) -> float:
    """
    计算综合评分
    
    Args:
        disease_score: 疾病关联分数 (0-1)
        e3_score: E3 兼容度分数 (0-1)
        weights: (disease_weight, e3_weight)，默认从配置读取
    
    Returns:
        综合评分 (0-1)
    
    Example:
        >>> score = calculate_composite_score(0.8, 0.6)
        >>> print(f"{score:.2f}")  # 0.72 (0.8*0.6 + 0.6*0.4)
    """
    if weights is None:
        from .disease_config import COMPOSITE_SCORE_WEIGHTS
        weights = (
            COMPOSITE_SCORE_WEIGHTS["disease"],
            COMPOSITE_SCORE_WEIGHTS["e3"]
        )
    
    # 验证权重
    if len(weights) != 2:
        raise ValueError("weights must be a tuple of 2 values")
    
    if not (0 <= weights[0] <= 1 and 0 <= weights[1] <= 1):
        raise ValueError("weights must be in range [0, 1]")
    
    # 归一化权重（确保总和为 1）
    total_weight = weights[0] + weights[1]
    if total_weight == 0:
        raise ValueError("weights cannot both be 0")
    
    normalized_weights = (weights[0] / total_weight, weights[1] / total_weight)
    
    # 计算综合评分
    composite = (
        disease_score * normalized_weights[0] +
        e3_score * normalized_weights[1]
    )
    
    return composite


def enrich_targets_with_e3_scores(
    targets_df: "pd.DataFrame",
    e3_symbol: str = None,
    include_composite: bool = True,
    weights: tuple = None
) -> "pd.DataFrame":
    """
    为靶点列表添加 E3 评分和综合评分
    
    Args:
        targets_df: get_disease_targets() 的输出
        e3_symbol: E3 连接酶符号（默认使用 DEFAULT_E3_LIGASE）
        include_composite: 是否计算综合评分（默认 True）
        weights: 综合评分权重（默认从配置读取）
    
    Returns:
        扩展的 DataFrame，包含 e3_score 和 composite_score 列
    
    Example:
        >>> from gluetk.open_targets_api import get_disease_targets, enrich_targets_with_e3_scores
        >>> df = get_disease_targets("EFO_0001378", top_n=10)
        >>> df_enriched = enrich_targets_with_e3_scores(df, "CRBN")
        >>> print(df_enriched[['symbol', 'score', 'e3_score', 'composite_score']])
    """
    if pd is None:
        raise RuntimeError("pandas is required")
    
    # 导入 E3 兼容度模块
    try:
        from .e3_compatibility import batch_calculate_e3_scores
    except ImportError as e:
        print(f"❌ Failed to import e3_compatibility: {e}")
        raise
    
    # 批量计算 E3 评分
    enriched_df = batch_calculate_e3_scores(
        targets_df,
        e3_symbol=e3_symbol
    )
    
    # 计算综合评分
    if include_composite:
        composite_scores = []
        for _, row in enriched_df.iterrows():
            composite = calculate_composite_score(
                row['score'],
                row['e3_score'],
                weights=weights
            )
            composite_scores.append(composite)
        
        enriched_df['composite_score'] = composite_scores
        
        # 按综合评分重新排序
        enriched_df = enriched_df.sort_values(
            'composite_score',
            ascending=False
        ).reset_index(drop=True)
    
    return enriched_df


def print_enriched_targets_table(df: "pd.DataFrame", top_n: int = 20) -> None:
    """
    打印包含 E3 评分和综合评分的靶点列表
    
    Args:
        df: enrich_targets_with_e3_scores() 的输出
        top_n: 显示的行数（默认 20）
    """
    if df is None or df.empty:
        print("No targets to display.")
        return
    
    print("\n" + "=" * 120)
    print("DISEASE-ASSOCIATED TARGETS WITH E3 COMPATIBILITY SCORES")
    print("=" * 120)
    
    # 显示前 N 行
    display_df = df.head(top_n).copy()
    
    # 选择要显示的列
    display_columns = ['symbol', 'score', 'e3_score']
    if 'composite_score' in display_df.columns:
        display_columns.append('composite_score')
    display_columns.append('e3_ligase')
    
    # 格式化分数
    for col in ['score', 'e3_score', 'composite_score']:
        if col in display_df.columns:
            display_df[col] = display_df[col].apply(lambda x: f"{x:.4f}")
    
    # 截断长名称
    if 'name' in display_df.columns:
        display_df['name'] = display_df['name'].apply(
            lambda x: x[:35] + "..." if len(x) > 35 else x
        )
    
    # 打印表格
    print(display_df[display_columns].to_string(index=True))
    
    if len(df) > top_n:
        print(f"\n... and {len(df) - top_n} more targets")
    
    print("=" * 120 + "\n")
