#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test_ec_performance.py
EC 计算性能测试脚本

用于对比原始版本和优化版本的性能差异
"""

import os
import sys
import time
import tempfile
import numpy as np
from typing import Dict, Tuple

# 添加 gluetk 到路径
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'gluetk'))

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("⚠️ RDKit not available - skipping tests")

try:
    from scipy.spatial import cKDTree
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("⚠️ SciPy not available - skipping KDTree tests")


class PerformanceBenchmark:
    """性能基准测试"""
    
    def __init__(self):
        self.results = {}
    
    def benchmark(self, name: str, func, *args, **kwargs) -> float:
        """运行基准测试"""
        print(f"\n[Benchmark] {name}...", end='', flush=True)
        
        start_time = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start_time
        
        self.results[name] = elapsed
        print(f" ✓ {elapsed:.3f}s")
        
        return elapsed
    
    def print_summary(self):
        """打印摘要"""
        print("\n" + "="*60)
        print("Performance Benchmark Summary")
        print("="*60)
        
        for name, elapsed in sorted(self.results.items(), key=lambda x: x[1], reverse=True):
            print(f"{name:<40} {elapsed:>8.3f}s")
        
        print("="*60)


# ========== 测试用例 ==========

def test_coulomb_potential_original(coords: np.ndarray, charges: np.ndarray,
                                    points: np.ndarray, dielectric: float = 4.0) -> np.ndarray:
    """原始版本 - Python 循环"""
    potentials = np.zeros(len(points))
    
    for i, point in enumerate(points):
        distances = np.linalg.norm(coords - point, axis=1)
        distances = np.maximum(distances, 0.1)
        potentials[i] = np.sum(charges / (dielectric * distances))
    
    potentials *= 332.0637 / 0.593
    return potentials


def test_coulomb_potential_vectorized(coords: np.ndarray, charges: np.ndarray,
                                      points: np.ndarray, dielectric: float = 4.0) -> np.ndarray:
    """优化版本 - 向量化计算"""
    n_points = len(points)
    chunk_size = 1000
    potentials = np.zeros(n_points)
    
    for chunk_start in range(0, n_points, chunk_size):
        chunk_end = min(chunk_start + chunk_size, n_points)
        chunk_points = points[chunk_start:chunk_end]
        
        # 向量化计算距离
        diff = chunk_points[:, np.newaxis, :] - coords[np.newaxis, :, :]
        distances = np.linalg.norm(diff, axis=2)
        distances = np.maximum(distances, 0.1)
        
        chunk_potentials = np.sum(charges[np.newaxis, :] / (dielectric * distances), axis=1)
        potentials[chunk_start:chunk_end] = chunk_potentials
    
    potentials *= 332.0637 / 0.593
    return potentials


def test_surface_sampling_original(atoms: list) -> Tuple[np.ndarray, np.ndarray]:
    """原始版本 - 嵌套循环"""
    all_points = []
    all_normals = []
    
    for atom in atoms:
        points, normals = _sample_sphere(atom['pos'], atom['radius'])
        
        # 嵌套循环检查内部点
        valid_mask = np.ones(len(points), dtype=bool)
        for other_atom in atoms:
            if np.allclose(atom['pos'], other_atom['pos']):
                continue
            distances = np.linalg.norm(points - other_atom['pos'], axis=1)
            valid_mask &= (distances >= other_atom['radius'] - 0.1)
        
        all_points.append(points[valid_mask])
        all_normals.append(normals[valid_mask])
    
    return np.vstack(all_points), np.vstack(all_normals)


def test_surface_sampling_kdtree(atoms: list) -> Tuple[np.ndarray, np.ndarray]:
    """优化版本 - KDTree 加速"""
    if not SCIPY_AVAILABLE:
        print("⚠️ SciPy not available, skipping KDTree test")
        return None
    
    all_points = []
    all_normals = []
    
    atom_positions = np.array([a['pos'] for a in atoms])
    atom_radii = np.array([a['radius'] for a in atoms])
    tree = cKDTree(atom_positions)
    
    for i, atom in enumerate(atoms):
        points, normals = _sample_sphere(atom['pos'], atom['radius'])
        
        # 使用 KDTree 快速查询
        valid_mask = np.ones(len(points), dtype=bool)
        distances, indices = tree.query(points, k=len(atoms))
        
        for j, (dists, inds) in enumerate(zip(distances, indices)):
            for dist, atom_idx in zip(dists, inds):
                if atom_idx != i:
                    if dist < atom_radii[atom_idx] - 0.1:
                        valid_mask[j] = False
                        break
        
        all_points.append(points[valid_mask])
        all_normals.append(normals[valid_mask])
    
    return np.vstack(all_points), np.vstack(all_normals)


def _sample_sphere(center: np.ndarray, radius: float, density: float = 10.0) -> Tuple[np.ndarray, np.ndarray]:
    """在球面上生成均匀分布的点"""
    area = 4 * np.pi * radius ** 2
    n_points = max(int(area * density), 20)
    
    points = []
    normals = []
    
    golden_ratio = (1 + np.sqrt(5)) / 2
    
    for i in range(n_points):
        theta = 2 * np.pi * i / golden_ratio
        phi = np.arccos(1 - 2 * (i + 0.5) / n_points)
        
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        
        normal = np.array([x, y, z])
        point = center + radius * normal
        
        points.append(point)
        normals.append(normal)
    
    return np.array(points), np.array(normals)


# ========== 生成测试数据 ==========

def generate_test_data(n_atoms: int = 50, n_points: int = 1000) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """生成测试数据"""
    print(f"\n[Setup] Generating test data: {n_atoms} atoms, {n_points} surface points")
    
    # 随机原子坐标
    coords = np.random.randn(n_atoms, 3) * 5
    
    # 随机电荷
    charges = np.random.randn(n_atoms) * 0.5
    
    # 随机表面点
    points = np.random.randn(n_points, 3) * 10
    
    return coords, charges, points


def generate_test_atoms(n_atoms: int = 10) -> list:
    """生成测试原子"""
    print(f"\n[Setup] Generating test atoms: {n_atoms}")
    
    VDW_RADII = {
        'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'H': 1.20
    }
    
    atoms = []
    for i in range(n_atoms):
        element = list(VDW_RADII.keys())[i % len(VDW_RADII)]
        atoms.append({
            'pos': np.random.randn(3) * 5,
            'radius': VDW_RADII[element] + 1.4,  # 加上探针半径
            'symbol': element
        })
    
    return atoms


# ========== 主测试函数 ==========

def run_coulomb_benchmark():
    """运行 Coulomb 势能计算基准测试"""
    print("\n" + "="*60)
    print("Coulomb Potential Calculation Benchmark")
    print("="*60)
    
    benchmark = PerformanceBenchmark()
    
    # 测试不同规模
    test_cases = [
        (50, 500),    # 50 原子, 500 点
        (50, 1000),   # 50 原子, 1000 点
        (100, 1000),  # 100 原子, 1000 点
        (100, 2000),  # 100 原子, 2000 点
    ]
    
    for n_atoms, n_points in test_cases:
        coords, charges, points = generate_test_data(n_atoms, n_points)
        
        # 原始版本
        name_orig = f"Original ({n_atoms} atoms, {n_points} points)"
        t_orig = benchmark.benchmark(name_orig, test_coulomb_potential_original,
                                     coords, charges, points)
        
        # 优化版本
        name_opt = f"Vectorized ({n_atoms} atoms, {n_points} points)"
        t_opt = benchmark.benchmark(name_opt, test_coulomb_potential_vectorized,
                                    coords, charges, points)
        
        # 计算加速比
        speedup = t_orig / t_opt
        print(f"  → Speedup: {speedup:.1f}x")
    
    benchmark.print_summary()


def run_surface_sampling_benchmark():
    """运行表面采样基准测试"""
    print("\n" + "="*60)
    print("Surface Sampling Benchmark")
    print("="*60)
    
    if not SCIPY_AVAILABLE:
        print("⚠️ SciPy not available, skipping surface sampling benchmark")
        return
    
    benchmark = PerformanceBenchmark()
    
    # 测试不同规模
    test_cases = [5, 10, 15, 20]
    
    for n_atoms in test_cases:
        atoms = generate_test_atoms(n_atoms)
        
        # 原始版本
        name_orig = f"Original ({n_atoms} atoms)"
        t_orig = benchmark.benchmark(name_orig, test_surface_sampling_original, atoms)
        
        # 优化版本
        name_opt = f"KDTree ({n_atoms} atoms)"
        result = test_surface_sampling_kdtree(atoms)
        if result is not None:
            t_opt = benchmark.benchmark(name_opt, test_surface_sampling_kdtree, atoms)
            
            # 计算加速比
            speedup = t_orig / t_opt
            print(f"  → Speedup: {speedup:.1f}x")
    
    benchmark.print_summary()


def run_memory_benchmark():
    """运行内存使用基准测试"""
    print("\n" + "="*60)
    print("Memory Usage Benchmark")
    print("="*60)
    
    import sys
    
    # 测试不同规模的数据
    test_cases = [
        (50, 500),
        (50, 1000),
        (100, 1000),
        (100, 2000),
    ]
    
    for n_atoms, n_points in test_cases:
        coords, charges, points = generate_test_data(n_atoms, n_points)
        
        # 计算内存使用
        coords_size = sys.getsizeof(coords) / 1024 / 1024  # MB
        charges_size = sys.getsizeof(charges) / 1024 / 1024
        points_size = sys.getsizeof(points) / 1024 / 1024
        total_size = coords_size + charges_size + points_size
        
        print(f"\n{n_atoms} atoms, {n_points} points:")
        print(f"  Coordinates: {coords_size:.2f} MB")
        print(f"  Charges: {charges_size:.2f} MB")
        print(f"  Points: {points_size:.2f} MB")
        print(f"  Total: {total_size:.2f} MB")


# ========== 主函数 ==========

def main():
    """主函数"""
    print("\n" + "="*60)
    print("GlueTK EC Calculator Performance Benchmark")
    print("="*60)
    print(f"NumPy version: {np.__version__}")
    print(f"SciPy available: {SCIPY_AVAILABLE}")
    print(f"RDKit available: {RDKIT_AVAILABLE}")
    
    # 运行基准测试
    run_coulomb_benchmark()
    run_surface_sampling_benchmark()
    run_memory_benchmark()
    
    print("\n" + "="*60)
    print("Benchmark Complete")
    print("="*60)


if __name__ == '__main__':
    main()
