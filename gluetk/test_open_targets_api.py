# -*- coding: utf-8 -*-
"""
Test Open Targets API Integration
测试 Open Targets API 封装功能

运行方式:
    python test_open_targets_api.py
"""

import sys
import os

# 添加 gluetk 到路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_api_imports():
    """测试模块导入"""
    print("\n" + "=" * 80)
    print("TEST 1: Module Imports")
    print("=" * 80)
    
    try:
        from gluetk.disease_config import (
            OPEN_TARGETS_GRAPHQL_ENDPOINT,
            DEFAULT_OUTPUT_DIR
        )
        print("✅ disease_config imported successfully")
        print(f"   API Endpoint: {OPEN_TARGETS_GRAPHQL_ENDPOINT}")
        print(f"   Output Dir: {DEFAULT_OUTPUT_DIR}")
    except Exception as e:
        print(f"❌ Failed to import disease_config: {e}")
        return False
    
    try:
        from gluetk.open_targets_api import (
            search_disease,
            get_disease_targets
        )
        print("✅ open_targets_api imported successfully")
    except Exception as e:
        print(f"❌ Failed to import open_targets_api: {e}")
        return False
    
    try:
        from gluetk.disease_target_commands import (
            ot_disease_targets,
            ot_glue_insight
        )
        print("✅ disease_target_commands imported successfully")
    except Exception as e:
        print(f"❌ Failed to import disease_target_commands: {e}")
        return False
    
    return True


def test_disease_search():
    """测试疾病搜索功能"""
    print("\n" + "=" * 80)
    print("TEST 2: Disease Search")
    print("=" * 80)
    
    try:
        from gluetk.open_targets_api import search_disease, print_disease_candidates
        
        # 测试搜索
        results = search_disease("multiple myeloma", max_results=5)
        
        if not results:
            print("❌ No results returned")
            return False
        
        print(f"✅ Found {len(results)} disease candidates")
        print_disease_candidates(results)
        
        # 验证结果格式
        first = results[0]
        required_keys = ['id', 'name', 'description', 'entity']
        
        for key in required_keys:
            if key not in first:
                print(f"❌ Missing key in result: {key}")
                return False
        
        print("✅ Result format validated")
        return True
        
    except Exception as e:
        print(f"❌ Disease search failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_target_query():
    """测试靶点查询功能"""
    print("\n" + "=" * 80)
    print("TEST 3: Target Query")
    print("=" * 80)
    
    try:
        from gluetk.open_targets_api import get_disease_targets, print_targets_table
        
        # 使用多发性骨髓瘤的 EFO ID
        disease_id = "EFO_0001378"
        
        df = get_disease_targets(disease_id, top_n=10)
        
        if df is None or df.empty:
            print("❌ No targets returned")
            return False
        
        print(f"✅ Retrieved {len(df)} targets")
        print_targets_table(df, top_n=10)
        
        # 验证 DataFrame 格式
        required_cols = ['symbol', 'name', 'score', 'uniprot_id', 'ensembl_id']
        
        for col in required_cols:
            if col not in df.columns:
                print(f"❌ Missing column: {col}")
                return False
        
        print("✅ DataFrame format validated")
        
        # 验证分数范围
        if not all(0 <= score <= 1 for score in df['score']):
            print("❌ Invalid score values (should be 0-1)")
            return False
        
        print("✅ Score values validated")
        return True
        
    except Exception as e:
        print(f"❌ Target query failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_file_output():
    """测试文件输出功能"""
    print("\n" + "=" * 80)
    print("TEST 4: File Output")
    print("=" * 80)
    
    try:
        from gluetk.open_targets_api import get_disease_targets
        from gluetk.disease_config import get_targets_csv_path
        import pandas as pd
        
        disease_id = "EFO_0001378"
        output_dir = "./test_output"
        
        # 获取数据
        df = get_disease_targets(disease_id, top_n=10)
        
        if df is None or df.empty:
            print("❌ No data to save")
            return False
        
        # 保存 CSV
        csv_path = get_targets_csv_path(disease_id, output_dir)
        os.makedirs(output_dir, exist_ok=True)
        df.to_csv(csv_path, index=False)
        
        print(f"✅ Saved to: {csv_path}")
        
        # 验证文件存在
        if not os.path.exists(csv_path):
            print("❌ File not created")
            return False
        
        # 验证文件可读
        df_read = pd.read_csv(csv_path)
        
        if len(df_read) != len(df):
            print("❌ File content mismatch")
            return False
        
        print("✅ File validated")
        
        # 清理测试文件
        os.remove(csv_path)
        if os.path.exists(output_dir) and not os.listdir(output_dir):
            os.rmdir(output_dir)
        
        return True
        
    except Exception as e:
        print(f"❌ File output test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def create_sample_residue_scores():
    """创建示例残基评分文件"""
    print("\n" + "=" * 80)
    print("Creating Sample Residue Scores")
    print("=" * 80)
    
    try:
        import pandas as pd
        from gluetk.disease_config import RESIDUE_SCORES_DIR
        
        # 创建目录
        os.makedirs(RESIDUE_SCORES_DIR, exist_ok=True)
        
        # 创建示例数据（IKZF1）
        sample_data = {
            'residue_index': [23, 45, 67, 89, 102, 125, 148, 167, 189, 203],
            'SurfHot_score': [0.85, 0.92, 0.65, 0.78, 0.88, 0.71, 0.95, 0.62, 0.81, 0.73],
            'motif_flag': [1, 1, 0, 1, 1, 0, 1, 0, 1, 0]
        }
        
        df = pd.DataFrame(sample_data)
        
        # 保存文件
        csv_path = os.path.join(RESIDUE_SCORES_DIR, "residue_scores_IKZF1.csv")
        df.to_csv(csv_path, index=False)
        
        print(f"✅ Created sample file: {csv_path}")
        print(f"   Contains {len(df)} residues")
        print(f"   Score range: {df['SurfHot_score'].min():.2f} - {df['SurfHot_score'].max():.2f}")
        
        return True
        
    except Exception as e:
        print(f"❌ Failed to create sample data: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """运行所有测试"""
    print("\n" + "🧪" * 40)
    print("OPEN TARGETS API INTEGRATION TEST SUITE")
    print("🧪" * 40)
    
    results = []
    
    # 运行测试
    results.append(("Module Imports", test_api_imports()))
    results.append(("Disease Search", test_disease_search()))
    results.append(("Target Query", test_target_query()))
    results.append(("File Output", test_file_output()))
    results.append(("Sample Data Creation", create_sample_residue_scores()))
    
    # 汇总结果
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {test_name}")
    
    print("=" * 80)
    print(f"Total: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed!")
        print("\n📌 Next Steps:")
        print("1. Install requests if not already: pip install requests")
        print("2. Test in PyMOL:")
        print("   run /path/to/gluetk/disease_target_commands.py")
        print("   ot_disease_targets disease_name=\"multiple myeloma\", top_n=30")
        return 0
    else:
        print("⚠️ Some tests failed. Please check the errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
