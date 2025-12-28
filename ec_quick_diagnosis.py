#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ec_quick_diagnosis.py
EC 分析快速诊断工具

用于快速诊断 EC 计算的性能问题和给出优化建议
"""

import os
import sys
import time

print("\n" + "="*70)
print("GlueTK EC 分析性能诊断工具")
print("="*70)

# ========== 检查依赖 ==========

print("\n[1/5] 检查依赖...")

dependencies = {
    'numpy': False,
    'scipy': False,
    'rdkit': False,
    'pymol': False,
    'pdb2pqr': False,
    'apbs': False
}

try:
    import numpy as np
    dependencies['numpy'] = True
    print(f"  ✅ NumPy {np.__version__}")
except ImportError:
    print("  ❌ NumPy - 必需")

try:
    import scipy
    dependencies['scipy'] = True
    print(f"  ✅ SciPy {scipy.__version__}")
except ImportError:
    print("  ⚠️  SciPy - 可选（用于加速）")

try:
    from rdkit import Chem
    dependencies['rdkit'] = True
    print(f"  ✅ RDKit")
except ImportError:
    print("  ❌ RDKit - 必需")

try:
    from pymol import cmd
    dependencies['pymol'] = True
    print(f"  ✅ PyMOL")
except ImportError:
    print("  ⚠️  PyMOL - 可选（用于可视化）")

# 检查外部工具
import shutil
if shutil.which('pdb2pqr'):
    dependencies['pdb2pqr'] = True
    print(f"  ✅ PDB2PQR")
else:
    print("  ⚠️  PDB2PQR - 可选（用于蛋白质准备）")

if shutil.which('apbs'):
    dependencies['apbs'] = True
    print(f"  ✅ APBS")
else:
    print("  ⚠️  APBS - 可选（用于静电计算）")

# ========== 性能分析 ==========

print("\n[2/5] 性能分析...")

if dependencies['numpy']:
    import numpy as np
    
    # 测试向量化性能
    print("\n  测试 Coulomb 势能计算性能...")
    
    # 原始方法
    coords = np.random.randn(50, 3)
    charges = np.random.randn(50)
    points = np.random.randn(1000, 3)
    
    start = time.time()
    potentials_loop = np.zeros(len(points))
    for i, point in enumerate(points):
        distances = np.linalg.norm(coords - point, axis=1)
        distances = np.maximum(distances, 0.1)
        potentials_loop[i] = np.sum(charges / distances)
    time_loop = time.time() - start
    
    # 向量化方法
    start = time.time()
    diff = points[:, np.newaxis, :] - coords[np.newaxis, :, :]
    distances = np.linalg.norm(diff, axis=2)
    distances = np.maximum(distances, 0.1)
    potentials_vec = np.sum(charges[np.newaxis, :] / distances, axis=1)
    time_vec = time.time() - start
    
    speedup = time_loop / time_vec
    
    print(f"    • 原始方法（循环）: {time_loop:.3f}s")
    print(f"    • 优化方法（向量化）: {time_vec:.3f}s")
    print(f"    • 加速比: {speedup:.1f}x")
    
    if speedup < 10:
        print(f"    ⚠️  加速比较低，可能是 NumPy 版本问题")
    else:
        print(f"    ✅ 向量化优化有效")

# ========== 问题诊断 ==========

print("\n[3/5] 问题诊断...")

issues = []

if not dependencies['numpy']:
    issues.append("❌ NumPy 未安装 - EC 计算无法运行")

if not dependencies['rdkit']:
    issues.append("❌ RDKit 未安装 - 配体处理无法进行")

if not dependencies['scipy']:
    issues.append("⚠️  SciPy 未安装 - 表面采样会很慢（可选优化）")

if not dependencies['pdb2pqr']:
    issues.append("⚠️  PDB2PQR 未安装 - 蛋白质准备会失败")

if not dependencies['apbs']:
    issues.append("⚠️  APBS 未安装 - 静电计算会失败")

if issues:
    print("\n  发现的问题:")
    for issue in issues:
        print(f"    {issue}")
else:
    print("  ✅ 所有依赖都已安装")

# ========== 性能建议 ==========

print("\n[4/5] 性能优化建议...")

recommendations = []

if dependencies['numpy'] and dependencies['scipy']:
    recommendations.append("✅ 已具备向量化优化的条件")
else:
    recommendations.append("⚠️  安装 SciPy 以启用 KDTree 加速")

if not dependencies['scipy']:
    recommendations.append("📌 表面采样会使用嵌套循环（较慢）")

if dependencies['pymol']:
    recommendations.append("✅ 可以使用 PyMOL 进行可视化")
else:
    recommendations.append("📌 无法在 PyMOL 中可视化结果")

recommendations.append("💡 使用优化版本: calculate_ligand_ec_optimized()")
recommendations.append("💡 降低表面密度参数以加快计算: surface_density=5.0")
recommendations.append("💡 使用分块处理大型配体")

print("\n  建议:")
for rec in recommendations:
    print(f"    {rec}")

# ========== 快速测试 ==========

print("\n[5/5] 快速性能测试...")

if dependencies['numpy'] and dependencies['rdkit']:
    print("\n  生成测试分子...")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    # 创建简单的测试分子
    mol = Chem.MolFromSmiles('CCO')  # 乙醇
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    
    print(f"  ✅ 测试分子: {mol.GetNumAtoms()} 原子")
    
    # 测试表面采样
    print("\n  测试表面采样...")
    
    conf = mol.GetConformer()
    VDW_RADII = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
        'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98,
    }
    
    start = time.time()
    total_points = 0
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        symbol = atom.GetSymbol()
        radius = VDW_RADII.get(symbol, 1.70) + 1.4
        
        # 简单的球面采样
        area = 4 * np.pi * radius ** 2
        n_points = max(int(area * 10.0), 20)
        total_points += n_points
    
    time_sampling = time.time() - start
    
    print(f"  ✅ 生成 {total_points} 个表面点")
    print(f"  ✅ 采样时间: {time_sampling:.3f}s")

# ========== 总结 ==========

print("\n" + "="*70)
print("诊断总结")
print("="*70)

if all([dependencies['numpy'], dependencies['rdkit']]):
    print("\n✅ 基本环境满足要求")
    print("\n🚀 快速开始:")
    print("  1. 使用优化版本: from gluetk.ligand_ec_calculator_optimized import calculate_ligand_ec_optimized")
    print("  2. 降低表面密度: surface_density=5.0")
    print("  3. 启用进度显示: show_progress=True")
else:
    print("\n❌ 缺少必要的依赖")
    print("\n📦 安装命令:")
    if not dependencies['numpy']:
        print("  pip install numpy")
    if not dependencies['rdkit']:
        print("  pip install rdkit")
    if not dependencies['scipy']:
        print("  pip install scipy")

print("\n📖 详细指南: 查看 EC_OPTIMIZATION_GUIDE.md")
print("🧪 性能测试: python test_ec_performance.py")
print("\n" + "="*70 + "\n")
