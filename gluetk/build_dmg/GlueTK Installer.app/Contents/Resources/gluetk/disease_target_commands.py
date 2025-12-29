# -*- coding: utf-8 -*-
"""
Disease Target Commands for PyMOL
PyMOL 命令接口：疾病-靶点分析和分子胶热点可视化

Author: Vesper
"""

import os
from typing import Optional

try:
    from pymol import cmd
    _PYMOL_AVAILABLE = True
except ImportError:
    print("⚠️ Warning: PyMOL not available")
    _PYMOL_AVAILABLE = False
    cmd = None

try:
    import pandas as pd
except ImportError:
    print("⚠️ Warning: pandas not installed")
    pd = None

from .open_targets_api import (
    search_disease,
    get_disease_targets,
    print_disease_candidates,
    print_targets_table
)

from .disease_config import (
    DEFAULT_OUTPUT_DIR,
    RESIDUE_SCORES_DIR,
    DEFAULT_HOTSPOT_THRESHOLD,
    HOTSPOT_COLOR,
    HOTSPOT_REPRESENTATION,
    BACKGROUND_COLOR,
    get_targets_csv_path,
    get_residue_scores_path
)


# ============================================================================
# 命令 1: 疾病 → 靶点列表查询
# ============================================================================

def ot_disease_targets(
    disease_name: str,
    top_n: int = 30,
    output_dir: str = None,
    auto_select: bool = True,
    include_e3_score: bool = True,
    e3_symbol: str = "CRBN",
    mining: bool = False
):
    """
    查询指定疾病的关联靶点列表（V2: 支持 E3 评分和综合评分）
    
    Args:
        disease_name: 疾病名称（英文），如 "multiple myeloma"
        top_n: 返回的靶点数量（默认 30）
        output_dir: 输出目录（默认 ./disease_data）
        auto_select: 如果搜索结果唯一，自动选择（默认 True）
        include_e3_score: 是否计算 E3 兼容度评分（默认 True）
        e3_symbol: E3 连接酶符号（默认 CRBN）
        mining: 是否启用数据挖掘模式（显示证据类型分数）（默认 False）
    
    Example:
        ot_disease_targets disease_name="multiple myeloma", top_n=30
        ot_disease_targets disease_name="multiple myeloma", mining=True
        ot_disease_targets disease_name="acute myeloid leukemia", top_n=50, include_e3_score=True
    
    Output:
        - 在 PyMOL 控制台打印靶点列表
        - 保存 CSV 文件到 {output_dir}/targets_{disease_id}.csv
    """
    if output_dir is None:
        output_dir = DEFAULT_OUTPUT_DIR
    
    print("\n" + "=" * 100)
    print("🧬 OPEN TARGETS DISEASE-TARGET ANALYSIS")
    print("=" * 100)
    
    try:
        # 步骤 1: 搜索疾病
        candidates = search_disease(disease_name, max_results=10)
        
        if not candidates:
            print(f"\n❌ No diseases found for '{disease_name}'")
            print("💡 Tip: Try different keywords or check spelling")
            return
        
        # 步骤 2: 选择疾病
        selected_disease = None
        
        if len(candidates) == 1 and auto_select:
            # 唯一结果，自动选择
            selected_disease = candidates[0]
            print(f"\n✅ Auto-selected: {selected_disease['name']} ({selected_disease['id']})")
        else:
            # 多个结果，显示候选列表
            print_disease_candidates(candidates)
            print("\n💡 Multiple candidates found. Using the first one.")
            print("   To select a different disease, note its ID and use get_disease_targets() directly.")
            selected_disease = candidates[0]
        
        disease_id = selected_disease['id']
        disease_display_name = selected_disease['name']
        
        # 步骤 3: 查询靶点
        df = get_disease_targets(disease_id, top_n=top_n)
        
        if df is None or df.empty:
            print(f"\n❌ No targets found for {disease_display_name}")
            return
        
        if mining:
            print("\n⛏️  DATA MINING MODE: ENABLED")
            print("   Showing evidence datatype scores (Genetic, Drug, etc.)")
            # 打印详细表格
            # 获取所有可能的列（除去基本列）
            basic_cols = ['symbol', 'name', 'score', 'uniprot_id', 'ensembl_id']
            mining_cols = [col for col in df.columns if col not in basic_cols and col != 'e3_score' and col != 'composite_score']
            
            print("\n" + "=" * 120)
            print("EVIDENCE MINING REPORT")
            print("=" * 120)
            
            display_df = df.head(min(top_n, 30)).copy()
            cols_to_show = ['symbol', 'score'] + mining_cols
            
            # 格式化
            for col in cols_to_show:
                if col in display_df.columns and display_df[col].dtype == float:
                    display_df[col] = display_df[col].apply(lambda x: f"{x:.3f}")
            
            print(display_df[cols_to_show].to_string(index=True))
            print("=" * 120 + "\n")

        # 步骤 3.5: 添加 E3 评分（V2 功能）
        if include_e3_score and not mining: # Mining 模式下可能不需要混淆太多信息，或者由用户决定
            try:
                from .open_targets_api import enrich_targets_with_e3_scores, print_enriched_targets_table
                print(f"\n🧬 Calculating E3 compatibility scores (E3: {e3_symbol})...")
                df = enrich_targets_with_e3_scores(df, e3_symbol=e3_symbol)
                print(f"✅ E3 scores calculated, sorted by composite score")
                
                # 使用增强版表格显示
                print_enriched_targets_table(df, top_n=min(top_n, 30))
            except Exception as e:
                print(f"⚠️ E3 scoring failed: {e}")
                print("   Falling back to disease scores only")
                print_targets_table(df, top_n=min(top_n, 30))
        else:
            # 步骤 4: 显示结果（仅疾病评分）
            print_targets_table(df, top_n=min(top_n, 30))  # 控制台最多显示 30 行
        
        # 步骤 5: 保存 CSV
        csv_path = get_targets_csv_path(disease_id, output_dir)
        os.makedirs(output_dir, exist_ok=True)
        df.to_csv(csv_path, index=False)
        print(f"💾 Saved {len(df)} targets to: {csv_path}")
        
        # 步骤 6: 使用提示
        print("\n" + "=" * 100)
        print("📌 NEXT STEPS")
        print("=" * 100)
        print("1. Load a target structure in PyMOL:")
        print(f"   fetch <PDB_ID>, <object_name>")
        print("   # or load AlphaFold structure")
        print("")
        print("2. Visualize molecular glue hotspots:")
        print(f"   ot_glue_insight protein_obj=<object_name>, gene_symbol=<SYMBOL>, threshold=0.7")
        print("")
        print(f"3. Available targets saved in: {csv_path}")
        print("=" * 100 + "\n")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


# ============================================================================
# 命令 2: 单靶点分子胶热点可视化
# ============================================================================

def ot_glue_insight(
    protein_obj: str,
    gene_symbol: str,
    threshold: float = None,
    scores_dir: str = None,
    disease_data_dir: str = None
):
    """
    可视化单个靶点的分子胶热点
    
    Args:
        protein_obj: PyMOL 对象名称
        gene_symbol: 基因符号（如 IKZF1）
        threshold: 热点阈值（默认 0.7）
        scores_dir: 残基评分文件目录（默认 ./disease_data/residue_scores）
        disease_data_dir: 疾病数据目录（默认 ./disease_data）
    
    Example:
        # 先加载结构
        fetch 6h0f, IKZF1
        
        # 可视化热点
        ot_glue_insight protein_obj="IKZF1", gene_symbol="IKZF1", threshold=0.7
    
    Requirements:
        - 残基评分文件: {scores_dir}/residue_scores_{gene_symbol}.csv
        - 列: residue_index, SurfHot_score
    """
    if not _PYMOL_AVAILABLE:
        print("❌ PyMOL is not available")
        return
    
    if threshold is None:
        threshold = DEFAULT_HOTSPOT_THRESHOLD
    
    if scores_dir is None:
        scores_dir = RESIDUE_SCORES_DIR
    
    if disease_data_dir is None:
        disease_data_dir = DEFAULT_OUTPUT_DIR
    
    print("\n" + "=" * 100)
    print("🔬 MOLECULAR GLUE HOTSPOT VISUALIZATION")
    print("=" * 100)
    print(f"Target: {gene_symbol}")
    print(f"PyMOL Object: {protein_obj}")
    print(f"Threshold: {threshold}")
    print("=" * 100)
    
    try:
        # 步骤 1: 检查 PyMOL 对象是否存在
        if protein_obj not in cmd.get_object_list():
            print(f"\n❌ Object '{protein_obj}' not found in PyMOL")
            print("💡 Tip: Load the structure first using 'fetch' or 'load' command")
            return
        
        # 步骤 2: 读取疾病关联分数（可选）
        disease_score = None
        e3_score = None
        
        # 尝试从最近的 targets CSV 中读取（包括 E3 评分）
        csv_files = []
        if os.path.exists(disease_data_dir):
            csv_files = [f for f in os.listdir(disease_data_dir) if f.startswith("targets_") and f.endswith(".csv")]
        
        if csv_files and pd is not None:
            # 使用最新的文件
            latest_csv = sorted(csv_files)[-1]
            csv_path = os.path.join(disease_data_dir, latest_csv)
            
            try:
                targets_df = pd.read_csv(csv_path)
                target_row = targets_df[targets_df['symbol'] == gene_symbol]
                
                if not target_row.empty:
                    disease_score = target_row.iloc[0]['score']
                    print(f"📊 Disease Association Score: {disease_score:.4f}")
                    
                    # V2: 读取 E3 评分和综合评分
                    if 'e3_score' in target_row.columns:
                        e3_score = target_row.iloc[0]['e3_score']
                        print(f"🧬 E3 Compatibility Score: {e3_score:.4f}")
                    
                    if 'composite_score' in target_row.columns:
                        composite_score = target_row.iloc[0]['composite_score']
                        print(f"⭐ Composite Score: {composite_score:.4f}")
                    
                    print(f"   (from {latest_csv})")
            except Exception as e:
                print(f"⚠️ Could not read scores: {e}")
        
        # 步骤 3: 读取残基评分文件
        residue_scores_path = get_residue_scores_path(gene_symbol, scores_dir)
        
        if not os.path.exists(residue_scores_path):
            print(f"\n❌ Residue scores file not found: {residue_scores_path}")
            print("💡 Tip: Generate residue scores using G-loop analysis first")
            print(f"   Expected format: residue_index, SurfHot_score")
            return
        
        print(f"\n📂 Loading residue scores from: {residue_scores_path}")
        
        if pd is None:
            print("❌ pandas is not installed")
            return
        
        scores_df = pd.read_csv(residue_scores_path)
        
        # 验证列名
        required_cols = ['residue_index', 'SurfHot_score']
        missing_cols = [col for col in required_cols if col not in scores_df.columns]
        
        if missing_cols:
            print(f"❌ Missing required columns: {missing_cols}")
            print(f"   Available columns: {list(scores_df.columns)}")
            return
        
        # 步骤 4: 过滤热点残基
        hotspots = scores_df[scores_df['SurfHot_score'] >= threshold]
        
        if hotspots.empty:
            print(f"\n⚠️ No hotspots found with threshold >= {threshold}")
            print(f"   Max score in file: {scores_df['SurfHot_score'].max():.3f}")
            print("💡 Tip: Try lowering the threshold")
            return
        
        print(f"✅ Found {len(hotspots)} hotspot residues (score >= {threshold})")
        
        # 步骤 5: PyMOL 可视化
        selection_name = f"{protein_obj}_{gene_symbol}_hotspots"
        
        # 创建选择集
        residue_list = hotspots['residue_index'].tolist()
        residue_str = "+".join([str(int(r)) for r in residue_list])
        
        cmd.select(selection_name, f"{protein_obj} and resi {residue_str}")
        
        # 设置背景色（所有残基）
        cmd.color(BACKGROUND_COLOR, protein_obj)
        
        # 高亮热点
        cmd.color(HOTSPOT_COLOR, selection_name)
        cmd.show(HOTSPOT_REPRESENTATION, selection_name)
        
        # 设置视图
        cmd.zoom(protein_obj)
        cmd.center(selection_name)
        
        print(f"\n🎨 Visualization applied:")
        print(f"   Selection: {selection_name}")
        print(f"   Color: {HOTSPOT_COLOR}")
        print(f"   Representation: {HOTSPOT_REPRESENTATION}")
        
        # 步骤 6: 打印汇总信息
        print("\n" + "=" * 100)
        print("📈 SUMMARY")
        print("=" * 100)
        print(f"Gene: {gene_symbol}")
        
        if disease_score is not None:
            print(f"Disease Association Score: {disease_score:.4f}")
        
        if e3_score is not None:
            print(f"E3 Compatibility Score: {e3_score:.4f}")
        
        print(f"Hotspot Threshold: {threshold}")
        print(f"Number of Hotspot Residues: {len(hotspots)}")
        print(f"Hotspot Residues: {residue_str}")
        
        # 显示 Top 5 热点
        top_hotspots = hotspots.nlargest(5, 'SurfHot_score')
        print(f"\nTop 5 Hotspots:")
        for idx, row in top_hotspots.iterrows():
            print(f"   Residue {int(row['residue_index'])}: {row['SurfHot_score']:.3f}")
        
        print("=" * 100 + "\n")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


# ============================================================================
# 命令 3: 多疾病比较（V2 功能）
# ============================================================================

def ot_compare_diseases(
    disease_names: str,
    top_n: int = 30,
    output_dir: str = None,
    generate_plots: bool = True,
    include_e3_score: bool = True,
    e3_symbol: str = "CRBN"
):
    """
    比较多个疾病的关联靶点
    
    Args:
        disease_names: 疾病名称（逗号分隔），如 "multiple myeloma,acute myeloid leukemia"
        top_n: 每个疾病返回的靶点数（默认 30）
        output_dir: 输出目录（默认 ./disease_data/comparison）
        generate_plots: 是否生成图表（默认 True）
        include_e3_score: 是否包含 E3 评分（默认 True）
        e3_symbol: E3 连接酶符号（默认 CRBN）
    
    Example:
        ot_compare_diseases disease_names="multiple myeloma,acute myeloid leukemia", top_n=30
        ot_compare_diseases disease_names="breast cancer,ovarian cancer,lung cancer", generate_plots=True
    
    Output:
        - 共同靶点分析
        - Venn 图（2-3 个疾病）
        - 热图
        - CSV 报告
    """
    if output_dir is None:
        output_dir = os.path.join(DEFAULT_OUTPUT_DIR, "comparison")
    
    try:
        # 解析疾病名称
        disease_list = [name.strip() for name in disease_names.split(',')]
        
        if len(disease_list) < 2:
            print("❌ Need at least 2 diseases for comparison")
            print("💡 Tip: Separate disease names with commas")
            print("   Example: disease_names=\"multiple myeloma,acute myeloid leukemia\"")
            return
        
        # 导入比较模块
        try:
            from .disease_comparison import (
                compare_diseases,
                find_common_targets,
                export_comparison_report,
                print_comparison_summary
            )
        except ImportError as e:
            print(f"❌ Failed to import disease_comparison module: {e}")
            return
        
        # 执行比较
        results = compare_diseases(
            disease_list,
            top_n=top_n,
            include_e3_score=include_e3_score,
            e3_symbol=e3_symbol
        )
        
        if not results:
            print("❌ No results obtained")
            return
        
        # 打印摘要
        print_comparison_summary(results)
        
        # 导出报告
        generated_files = export_comparison_report(
            results,
            output_dir=output_dir,
            generate_plots=generate_plots
        )
        
        # 打印生成的文件
        print("\n" + "=" * 100)
        print("📁 GENERATED FILES")
        print("=" * 100)
        for file_type, file_path in generated_files.items():
            print(f"  • {file_type}: {file_path}")
        print("=" * 100 + "\n")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()


# ============================================================================
# 命令注册（在 __init__.py 中调用）
# ============================================================================

def register_disease_commands():
    """
    注册疾病分析命令到 PyMOL
    """
    if not _PYMOL_AVAILABLE:
        print("⚠️ PyMOL not available, skipping disease command registration")
        return
    
    cmd.extend("ot_disease_targets", ot_disease_targets)
    cmd.extend("ot_glue_insight", ot_glue_insight)
    cmd.extend("ot_compare_diseases", ot_compare_diseases)  # V2 新增
    
    print("✅ Registered disease analysis commands: ot_disease_targets, ot_glue_insight, ot_compare_diseases")


# ============================================================================
# 测试函数
# ============================================================================

def test_commands():
    """
    测试命令功能（用于开发调试）
    """
    print("🧪 Testing Disease Target Commands\n")
    
    # 测试命令 1
    print("Test 1: ot_disease_targets")
    print("-" * 50)
    ot_disease_targets("multiple myeloma", top_n=10)
    
    print("\n" + "=" * 100)
    print("✅ Test completed. Check output above.")
    print("=" * 100)


if __name__ == "__main__":
    # 直接运行此文件时执行测试
    test_commands()
