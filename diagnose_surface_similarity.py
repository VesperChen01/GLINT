#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
诊断表面相似性分析问题的脚本

检查为什么不同蛋白的相似性分析结果都差不多一模一样
"""

import os
import sys
import numpy as np
from pathlib import Path

# 添加 gluetk 到路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'gluetk'))

print("\n" + "="*70)
print("表面相似性分析诊断工具")
print("="*70)

# ========== 检查依赖 ==========

print("\n[1/5] 检查依赖...")

try:
    import numpy as np
    print(f"  ✅ NumPy {np.__version__}")
except ImportError:
    print("  ❌ NumPy 未安装")
    sys.exit(1)

try:
    from scipy.spatial import cKDTree
    print(f"  ✅ SciPy (cKDTree)")
except ImportError:
    print("  ❌ SciPy 未安装")
    sys.exit(1)

try:
    from surface_similarity import SurfaceSimilarityAnalyzer, FeatureExtractor
    print(f"  ✅ surface_similarity 模块")
except ImportError:
    print("  ❌ surface_similarity 模块导入失败")
    sys.exit(1)

# ========== 测试特征提取 ==========

print("\n[2/5] 测试特征提取...")

# 创建两个不同的测试蛋白
def create_test_protein(name, size=50, center_offset=None):
    """创建测试蛋白"""
    atoms = []
    if center_offset is None:
        center_offset = np.array([0, 0, 0])
    
    # 创建球形蛋白
    for i in range(size):
        angle1 = 2 * np.pi * i / size
        angle2 = np.pi * (i % 10) / 10
        
        x = center_offset[0] + 10 * np.sin(angle2) * np.cos(angle1)
        y = center_offset[1] + 10 * np.sin(angle2) * np.sin(angle1)
        z = center_offset[2] + 10 * np.cos(angle2)
        
        # 交替使用不同的残基类型
        resn = ['ALA', 'GLY', 'VAL', 'LEU', 'ILE'][i % 5]
        
        atoms.append({
            'coord': np.array([x, y, z]),
            'element': 'C',
            'radius': 1.7,
            'chain': 'A',
            'resn': resn,
            'resi': str(i),
            'name': 'CA'
        })
    
    return atoms

# 创建两个不同的蛋白
protein1_atoms = create_test_protein("Protein1", size=50, center_offset=np.array([0, 0, 0]))
protein2_atoms = create_test_protein("Protein2", size=50, center_offset=np.array([20, 0, 0]))

print(f"  ✅ 创建测试蛋白 1: {len(protein1_atoms)} 个原子")
print(f"  ✅ 创建测试蛋白 2: {len(protein2_atoms)} 个原子")

# ========== 提取特征 ==========

print("\n[3/5] 提取表面特征...")

try:
    extractor1 = FeatureExtractor(protein1_atoms)
    extractor2 = FeatureExtractor(protein2_atoms)
    
    # 检查预计算的残基属性
    print(f"\n  蛋白 1 残基属性:")
    for key, hydro in list(extractor1.residue_hydrophobicity.items())[:3]:
        print(f"    {key}: hydrophobicity={hydro:.3f}")
    
    print(f"\n  蛋白 2 残基属性:")
    for key, hydro in list(extractor2.residue_hydrophobicity.items())[:3]:
        print(f"    {key}: hydrophobicity={hydro:.3f}")
    
    print(f"\n  ✅ 特征提取成功")
    
except Exception as e:
    print(f"  ❌ 特征提取失败: {e}")
    import traceback
    traceback.print_exc()

# ========== 检查相似性计算 ==========

print("\n[4/5] 检查相似性计算逻辑...")

# 测试直方图相关性
def test_histogram_correlation():
    """测试直方图相关性计算"""
    
    # 创建两个不同的分布
    dist1 = np.random.normal(0, 1, 100)  # 均值 0
    dist2 = np.random.normal(5, 1, 100)  # 均值 5
    dist3 = np.random.normal(0, 1, 100)  # 均值 0（与 dist1 相同）
    
    # 计算相关性
    from surface_similarity import SurfaceComparator
    comparator = SurfaceComparator()
    
    corr_1_2 = comparator._histogram_correlation(dist1, dist2)
    corr_1_3 = comparator._histogram_correlation(dist1, dist3)
    
    print(f"  分布 1 vs 分布 2 (不同): {corr_1_2:.4f}")
    print(f"  分布 1 vs 分布 3 (相同): {corr_1_3:.4f}")
    
    if corr_1_2 > 0.8:
        print(f"  ⚠️  警告: 不同分布的相关性过高 ({corr_1_2:.4f})")
        return False
    
    if corr_1_3 < 0.7:
        print(f"  ⚠️  警告: 相同分布的相关性过低 ({corr_1_3:.4f})")
        return False
    
    print(f"  ✅ 直方图相关性计算正常")
    return True

try:
    test_histogram_correlation()
except Exception as e:
    print(f"  ❌ 直方图相关性测试失败: {e}")

# ========== 建议 ==========

print("\n[5/5] 诊断建议...")

print("""
可能的原因和解决方案:

1. 【特征分布相似】
   - 问题: 所有蛋白的表面特征分布都很相似
   - 解决: 增加特征的多样性
     * 使用更多的几何特征 (例如: 高斯曲率、平均曲率)
     * 使用更精确的电势计算 (APBS 而不是 Coulomb)
     * 添加更多化学特征 (疏水性、H键密度等)

2. 【相似性度量不敏感】
   - 问题: 相似性计算对特征变化不敏感
   - 解决: 改进相似性度量
     * 使用欧几里得距离而不是余弦相似度
     * 使用 Wasserstein 距离而不是 Bhattacharyya 系数
     * 添加权重以强调重要特征

3. 【曲率计算失败】
   - 问题: 曲率计算可能都返回 0
   - 解决: 使用更好的曲率估计方法
     * 安装 Open3D 以获得更精确的曲率
     * 使用更大的邻域进行曲率估计
     * 验证曲率值是否真的为 0

4. 【表面生成问题】
   - 问题: 所有蛋白的表面生成方式相同
   - 解决: 改进表面生成
     * 使用 MSMS 而不是 EDTSurf
     * 调整表面密度参数
     * 验证生成的表面是否正确

建议的调试步骤:
1. 导出表面特征到 CSV 文件
2. 检查特征值的分布是否真的相似
3. 比较不同蛋白的曲率值
4. 检查相似性矩阵是否真的都接近 1.0
""")

print("\n" + "="*70 + "\n")
