# -*- coding: utf-8 -*-
"""
ec_visualization.py
Enhanced EC Surface Visualization for PyMOL

This module provides publication-quality EC surface visualization,
creating smooth molecular surfaces colored by electrostatic complementarity
values with a red-white-green gradient.

Features:
- Smooth Gaussian/solvent-accessible surface generation
- Red-white-green color gradient (red=clash, white=neutral, green=complementary)
- Ligand sticks visible inside the surface
- Protein residues shown as lines around binding site
- High-quality rendering settings for publication figures

Usage:
    from glint.ec_visualization import visualize_ec_smooth_surface
    
    # After running calculate_ligand_ec:
    result = calculate_ligand_ec('complex', 'LIG', output_dir='./ec_output')
    visualize_ec_smooth_surface('complex', 'LIG', ec_result=result)

Author: GLINT Team
"""

from __future__ import print_function
import os
import tempfile
from typing import Dict, List, Optional, Any

# NumPy
try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

# SciPy
try:
    from scipy.spatial import cKDTree
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False

# PyMOL
try:
    from pymol import cmd
    from pymol.cgo import (
        BEGIN, END, VERTEX, COLOR, NORMAL, SPHERE, 
        TRIANGLES, ALPHA, LINEWIDTH, LINES
    )
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False


def _enhance_ec_contrast(ec_values: np.ndarray, power: float = 0.3) -> np.ndarray:
    """对 EC ValueApply非线性变换增强颜色对比degrees。

    通过 sign(x) * |x|^power 变换，使中间 EC Value（如 ±0.3）也能Display出
    明显的红/绿色，而非被线性映射压缩到接近白色。

    Args:
        ec_values: 原始 EC Value数组
        power: 变换指数，< 1 时增强中间Value饱和degrees（默认 0.3，
               进一步增强中间Value的颜色饱和degrees，让轻微互补/冲突更易辨识）

    Returns:
        增强后的 EC Value数组
    """
    return np.sign(ec_values) * np.abs(ec_values) ** power


def _auto_ec_range(ec_values: np.ndarray, default_range: tuple = (-1.0, 1.0),
                   percentile: float = 5.0) -> tuple:
    """根据实际 EC Value分布自适应计算颜色映射范围。

    如果数据的 5th/95th 百分位数远小于默认范围 (-1, 1)，则缩小范围
    以充分利用颜色梯degrees，避免大片白色。

    Args:
        ec_values: EC Value数组
        default_range: 默认映射范围
        percentile: 用于确定范围的百分位数

    Returns:
        (ec_min, ec_max) 元组
    """
    if len(ec_values) == 0:
        return default_range

    p_low = np.percentile(ec_values, percentile)
    p_high = np.percentile(ec_values, 100 - percentile)

    # 对称化范围，取绝对Value较大的一端
    abs_max = max(abs(p_low), abs(p_high), 0.1)  # 最小 0.1 防止范围过窄

    # 仅当实际范围远小于默认范围时才自适应（阈Value：默认范围的 70%）
    default_abs = max(abs(default_range[0]), abs(default_range[1]))
    if abs_max < default_abs * 0.7:
        print(f"[EC] 自适应颜色范围: ({-abs_max:.2f}, {abs_max:.2f})  "
              f"(原始数据 5th/95th 百分位: {p_low:.3f} / {p_high:.3f})")
        return (-abs_max, abs_max)

    return default_range



def _ec_5color_interpolate(t: float):
    """将归一化进degrees t ∈ [0,1] 映射到 5 色渐变的 RGB Value。
    
    5 色锚点：
      t=0.00 → 深红 [0.6, 0.0, 0.0]  (强 clash)
      t=0.25 → 浅红 [1.0, 0.4, 0.4]  (轻微 clash)
      t=0.50 → 白   [1.0, 1.0, 1.0]  (中性)
      t=0.75 → 浅绿 [0.4, 1.0, 0.4]  (轻微互补)
      t=1.00 → 深绿 [0.0, 0.6, 0.0]  (强互补)
    
    相邻锚点之间线性插Value。
    
    Args:
        t: 归一化进degrees [0, 1]
    
    Returns:
        (r, g, b) 元组，每个分量 ∈ [0, 1]
    """
    # 5 色锚点定义（增强饱和degrees，让颜色区分更明显）
    anchors = [
        (0.00, (0.5, 0.0, 0.0)),   # 深红（更暗，强 clash）
        (0.25, (0.95, 0.25, 0.25)), # 浅红（更饱和，轻微 clash）
        (0.50, (1.0, 1.0, 1.0)),    # 白（中性）
        (0.75, (0.25, 0.95, 0.25)), # 浅绿（更饱和，轻微互补）
        (1.00, (0.0, 0.5, 0.0)),    # 深绿（更暗，强互补）
    ]
    
    # 边界处理
    if t <= 0.0:
        return anchors[0][1]
    if t >= 1.0:
        return anchors[-1][1]
    
    # 找到 t 所在的区间并线性插Value
    for i in range(len(anchors) - 1):
        t0, c0 = anchors[i]
        t1, c1 = anchors[i + 1]
        if t0 <= t <= t1:
            s = (t - t0) / (t1 - t0)  # 区间内归一化进degrees
            r = c0[0] + s * (c1[0] - c0[0])
            g = c0[1] + s * (c1[1] - c0[1])
            b = c0[2] + s * (c1[2] - c0[2])
            return (r, g, b)
    
    # 不应该到达这里
    return (1.0, 1.0, 1.0)


def _apply_5color_ec_gradient(surface_obj: str, b_min: float, b_max: float):
    """在 PyMOL 中Apply 5 色 EC 渐变（深红-浅红-白-浅绿-深绿）。
    
    using cmd.set_color 定义 5 个自定义颜色锚点，然后用 cmd.spectrum
    的自定义色板Name进行渐变着色。比内置 red_white_green (3色) 更有区分degrees，
    能让用户直观分辨强clash、轻微clash、中性、轻微互补、强互补五个等级。
    
    颜色定义：
      深红 [0.5, 0.0, 0.0] = 强 clash    (EC < -0.5)
      浅红 [0.95, 0.25, 0.25] = 轻微 clash  (-0.5 < EC < -0.1)
      白   [1.0, 1.0, 1.0] = 中性        (-0.1 < EC < 0.1)
      浅绿 [0.25, 0.95, 0.25] = 轻微互补    (0.1 < EC < 0.5)
      深绿 [0.0, 0.5, 0.0] = 强互补      (EC > 0.5)
    
    Args:
        surface_obj: PyMOL surface 对象Name
        b_min: B-factor 最小Value（对应 EC 最小Value * 100）
        b_max: B-factor 最大Value（对应 EC 最大Value * 100）
    """
    # 定义 5 个自定义颜色锚点（增强饱和degrees，让颜色更鲜明）
    cmd.set_color('ec_dark_red',  [0.5, 0.0, 0.0])    # 强 clash（更深暗红）
    cmd.set_color('ec_light_red', [0.95, 0.25, 0.25])  # 轻微 clash（更饱和红）
    cmd.set_color('ec_white',     [1.0, 1.0, 1.0])     # 中性
    cmd.set_color('ec_light_green', [0.25, 0.95, 0.25]) # 轻微互补（更饱和绿）
    cmd.set_color('ec_dark_green',  [0.0, 0.5, 0.0])   # 强互补（更深暗绿）
    
    # using PyMOL 自定义色板Name列表进行 spectrum 着色
    # spectrum 会在这些颜色之间线性插Value
    palette = 'ec_dark_red ec_light_red ec_white ec_light_green ec_dark_green'
    cmd.spectrum('b', palette, surface_obj, minimum=b_min, maximum=b_max)
    
    print(f"[_apply_5color_ec_gradient] ✅ 已Apply 5 色 EC 渐变: "
          f"深红→浅红→白→浅绿→深绿 (B: {b_min:.0f}~{b_max:.0f})")



def visualize_ec_smooth_surface(obj_name: str, ligand_resname: str,
                                 ec_values: np.ndarray = None,
                                 surface_points: np.ndarray = None,
                                 ec_result: Dict = None,
                                 surface_type: str = 'gaussian',
                                 transparency: float = 0.0,
                                 show_ligand_sticks: bool = True,
                                 show_protein_lines: bool = True,
                                 protein_distance: float = 5.0,
                                 color_scheme: str = 'rwg',
                                 ec_range: tuple = (-1.0, 1.0),
                                 surface_quality: int = 2,
                                 ray_trace: bool = False,
                                 enhance_contrast: bool = True) -> bool:
    """
    Create a smooth EC-colored molecular surface visualization.
    
    This function generates a publication-quality surface visualization
    similar to electrostatic potential surfaces commonly shown in papers.
    
    Args:
        obj_name: PyMOL object name containing the complex
        ligand_resname: Ligand residue name (e.g., 'LIG')
        ec_values: EC values at surface points (from calculate_ligand_ec)
        surface_points: Surface point coordinates (from calculate_ligand_ec)
        ec_result: Full result dict from calculate_ligand_ec (alternative input)
        surface_type: Surface type - 'gaussian', 'solvent', or 'molecular'
        transparency: Surface transparency (0.0=opaque, 1.0=transparent)
        show_ligand_sticks: Show ligand as sticks inside surface
        show_protein_lines: Show nearby protein residues as lines
        protein_distance: Distance cutoff for showing protein residues (Å)
        color_scheme: Color scheme - 'rwg' (red-white-green) or 'bwr' (blue-white-red)
        ec_range: (min, max) EC values for color mapping
        surface_quality: PyMOL surface quality (0-4, higher=smoother)
        ray_trace: Whether to ray trace after visualization
        enhance_contrast: 是否Enable非线性对比degrees增强（默认 True）
        
    Returns:
        True if successful
    """
    if not PYMOL_AVAILABLE:
        print("[visualize_ec_smooth_surface] ❌ PyMOL required")
        return False
    
    if not NUMPY_AVAILABLE:
        print("[visualize_ec_smooth_surface] ❌ NumPy required")
        return False
    
    # Extract data from ec_result if provided
    if ec_result is not None:
        ec_values = ec_result.get('ec_values', ec_values)
        surface_points = ec_result.get('surface_points', surface_points)
    
    if ec_values is None or surface_points is None:
        print("[visualize_ec_smooth_surface] ❌ Need ec_values and surface_points")
        print("[visualize_ec_smooth_surface] 💡 Run calculate_ligand_ec first")
        return False
    
    print(f"[visualize_ec_smooth_surface] Creating smooth EC surface...")
    print(f"[visualize_ec_smooth_surface] Surface type: {surface_type}")
    print(f"[visualize_ec_smooth_surface] EC range: {ec_range}")
    
    # 自适应颜色范围：根据实际数据分布调整
    ec_range = _auto_ec_range(ec_values, default_range=ec_range)
    
    ligand_sel = f"{obj_name} and resn {ligand_resname}"
    surface_obj = f"ec_surface_{ligand_resname}"
    
    # Step 1: Create a copy of the ligand for surface generation
    cmd.create(surface_obj, ligand_sel)
    
    # Step 2: Map EC values to atom B-factors（含非线性增强）
    if SCIPY_AVAILABLE:
        _map_ec_to_bfactors_kdtree(surface_obj, surface_points, ec_values,
                                    enhance_contrast=enhance_contrast)
    else:
        _map_ec_to_bfactors_simple(surface_obj, surface_points, ec_values,
                                    enhance_contrast=enhance_contrast)
    
    # Step 3: Rebuild to apply B-factor changes
    cmd.rebuild(surface_obj)
    
    # Step 4: Generate and configure surface
    cmd.show('surface', surface_obj)
    
    # Set surface type
    surface_type_map = {
        'molecular': 0,
        'solvent': 1, 
        'gaussian': 2
    }
    cmd.set('surface_type', surface_type_map.get(surface_type, 2), surface_obj)
    
    # Set surface quality
    cmd.set('surface_quality', surface_quality, surface_obj)
    
    # Step 5: Apply color spectrum based on B-factors
    ec_min, ec_max = ec_range
    b_min = ec_min * 100  # Scale to B-factor range
    b_max = ec_max * 100
    
    if color_scheme == 'rwg':
        # === 5 色渐变方案：深红-浅红-白-浅绿-深绿 ===
        # 比原始 3 色 (red_white_green) 更有区分degrees，
        # 用户可直观分辨强clash/轻微clash/中性/轻微互补/强互补
        _apply_5color_ec_gradient(surface_obj, b_min, b_max)
    elif color_scheme == 'bwr':
        # Blue-White-Red gradient
        cmd.spectrum('b', 'blue_white_red', surface_obj, minimum=b_min, maximum=b_max)
    else:
        # 默认也using 5 色渐变
        _apply_5color_ec_gradient(surface_obj, b_min, b_max)
    
    # Step 6: Set surface properties for smooth appearance
    cmd.set('transparency', transparency, surface_obj)
    cmd.set('surface_color_smoothing', 1)
    cmd.set('surface_color_smoothing_threshold', 0.5)
    
    # Step 7: Show ligand sticks inside surface
    if show_ligand_sticks:
        cmd.show('sticks', ligand_sel)
        cmd.color('gray50', f"{ligand_sel} and elem C")
        cmd.color('blue', f"{ligand_sel} and elem N")
        cmd.color('red', f"{ligand_sel} and elem O")
        cmd.color('yellow', f"{ligand_sel} and elem S")
        cmd.color('green', f"{ligand_sel} and elem Cl")
        cmd.color('orange', f"{ligand_sel} and elem Br")
        cmd.set('stick_radius', 0.15, ligand_sel)
    
    # Step 8: Show nearby protein residues
    if show_protein_lines:
        protein_sel = f"{obj_name} and polymer within {protein_distance} of {ligand_sel}"
        cmd.show('lines', protein_sel)
        cmd.color('gray70', f"{protein_sel} and elem C")
        cmd.set('line_width', 1.5, protein_sel)
    
    # Step 9: Set protein cartoon transparency
    cmd.set('cartoon_transparency', 0.8, obj_name)
    
    # Step 10: Configure rendering settings
    _set_publication_rendering()
    
    # Step 11: Zoom to ligand
    cmd.zoom(ligand_sel, buffer=8)
    
    # Step 12: Optional ray tracing
    if ray_trace:
        print("[visualize_ec_smooth_surface] Ray tracing...")
        cmd.ray()
    
    print(f"[visualize_ec_smooth_surface] ✅ Created surface: {surface_obj}")
    print("[visualize_ec_smooth_surface] 🎨 5 色 EC 颜色图例:")
    print("  🟥 深红 = 强 clash        (EC < -0.5)   → 非常不利于结合，需优化")
    print("  🔴 浅红 = 轻微 clash      (-0.5 < EC < -0.1) → 不利，可考虑修饰")
    print("  ⚪ 白色 = 中性            (-0.1 < EC < 0.1)  → 无显著影响")
    print("  🟢 浅绿 = 轻微互补        (0.1 < EC < 0.5)   → 有利于结合")
    print("  🟩 深绿 = 强互补          (EC > 0.5)   → 非常有利，结合驱动力区域")
    print("")
    print("[visualize_ec_smooth_surface] 📊 Results解读指南:")
    print("  • 绿色越多 → 静电互补性越好 → 结合越有利")
    print("  • 红色越多 → 静电冲突越严重 → 结合越不利")
    print("  • 关注红色区域可指导药物化学优化方向")
    print("[visualize_ec_smooth_surface] 💡 Use 'ray' command for publication quality")
    
    # 自动在 PyMOL 3D 视图中Create色标条
    create_ec_legend(ec_range=ec_range, color_scheme=color_scheme)
    
    return True


def _map_ec_to_bfactors_kdtree(surface_obj: str, surface_points: np.ndarray, 
                                ec_values: np.ndarray,
                                enhance_contrast: bool = True):
    """Map EC values to atom B-factors using KDTree for efficiency.
    
    Args:
        enhance_contrast: 是否对 EC ValueApply非线性增强
    """
    # 非线性对比degrees增强
    if enhance_contrast:
        ec_values = _enhance_ec_contrast(ec_values)
    
    tree = cKDTree(surface_points)
    
    # Get atom coordinates and indices
    atom_coords = []
    atom_indices = []
    cmd.iterate_state(1, surface_obj,
                     "atom_coords.append([x,y,z]); atom_indices.append(index)",
                     space={'atom_coords': atom_coords, 'atom_indices': atom_indices})
    
    if not atom_coords:
        return
    
    atom_coords = np.array(atom_coords)
    
    # For each atom, find nearby surface points and average their EC values
    for coord, atom_idx in zip(atom_coords, atom_indices):
        # Find surface points within 2.5Å of this atom
        nearby_indices = tree.query_ball_point(coord, r=2.5)
        
        if nearby_indices:
            # Average EC of nearby surface points
            avg_ec = np.mean(ec_values[nearby_indices])
        else:
            # Use nearest point
            _, nearest_idx = tree.query(coord)
            avg_ec = ec_values[nearest_idx]
        
        # Scale EC to B-factor range (-99 to 99)
        bfactor = np.clip(avg_ec * 100, -99.99, 99.99)
        
        # Set B-factor for this atom
        cmd.alter(f"{surface_obj} and index {atom_idx}", f"b={bfactor}")


def _map_ec_to_bfactors_simple(surface_obj: str, surface_points: np.ndarray,
                                ec_values: np.ndarray,
                                enhance_contrast: bool = True):
    """Map EC values to atom B-factors using simple nearest neighbor.
    
    Args:
        enhance_contrast: 是否对 EC ValueApply非线性增强
    """
    # 非线性对比degrees增强
    if enhance_contrast:
        ec_values = _enhance_ec_contrast(ec_values)
    
    # Get atom coordinates and indices
    atom_coords = []
    atom_indices = []
    cmd.iterate_state(1, surface_obj,
                     "atom_coords.append([x,y,z]); atom_indices.append(index)",
                     space={'atom_coords': atom_coords, 'atom_indices': atom_indices})
    
    if not atom_coords:
        return
    
    atom_coords = np.array(atom_coords)
    
    for coord, atom_idx in zip(atom_coords, atom_indices):
        distances = np.linalg.norm(surface_points - coord, axis=1)
        nearest_idx = np.argmin(distances)
        bfactor = np.clip(ec_values[nearest_idx] * 100, -99.99, 99.99)
        cmd.alter(f"{surface_obj} and index {atom_idx}", f"b={bfactor}")


def _set_publication_rendering():
    """Set PyMOL rendering settings for publication-quality figures."""
    # Lighting
    cmd.set('ray_shadow', 0)
    cmd.set('ambient', 0.4)
    cmd.set('direct', 0.6)
    cmd.set('spec_reflect', 0.5)
    cmd.set('spec_power', 200)
    
    # Anti-aliasing
    cmd.set('antialias', 2)
    
    # Surface appearance
    cmd.set('surface_smooth_edges', 1)
    
    # Background
    cmd.bg_color('white')
    
    # Depth cueing
    cmd.set('depth_cue', 0)


def visualize_ec_cgo_surface(obj_name: str, ligand_resname: str,
                              ec_values: np.ndarray,
                              surface_points: np.ndarray,
                              surface_normals: np.ndarray = None,
                              sphere_scale: float = 0.25,
                              color_scheme: str = 'rwg',
                              ec_range: tuple = (-1.0, 1.0),
                              enhance_contrast: bool = True) -> bool:
    """
    Create EC visualization using CGO (Compiled Graphics Objects).
    
    This creates a point cloud representation where each surface point
    is rendered as a small sphere colored by its EC value.
    
    Args:
        obj_name: PyMOL object name
        ligand_resname: Ligand residue name
        ec_values: EC values at surface points
        surface_points: Surface point coordinates
        surface_normals: Surface normals (optional)
        sphere_scale: Size of each sphere
        color_scheme: 'rwg' (red-white-green) or 'bwr' (blue-white-red)
        ec_range: (min, max) EC values for color mapping
        enhance_contrast: 是否Enable非线性对比degrees增强（默认 True）
        
    Returns:
        True if successful
    """
    if not PYMOL_AVAILABLE:
        print("[visualize_ec_cgo_surface] ❌ PyMOL required")
        return False
    
    print(f"[visualize_ec_cgo_surface] Creating CGO surface with {len(surface_points)} points...")
    
    # 自适应颜色范围
    ec_range = _auto_ec_range(ec_values, default_range=ec_range)
    
    # 非线性对比degrees增强
    if enhance_contrast:
        ec_values = _enhance_ec_contrast(ec_values)
    
    ec_min, ec_max = ec_range
    
    # Color mapping function
    def ec_to_rgb(ec_val):
        """Map EC value to RGB color."""
        # Normalize to [0, 1]
        t = (ec_val - ec_min) / (ec_max - ec_min + 1e-10)
        t = np.clip(t, 0, 1)
        
        if color_scheme == 'rwg':
            # Red-White-Green
            if t < 0.5:
                # Red to White
                s = t * 2
                r, g, b = 1.0, s, s
            else:
                # White to Green
                s = (t - 0.5) * 2
                r, g, b = 1.0 - s, 1.0, 1.0 - s
        else:
            # Blue-White-Red
            if t < 0.5:
                s = t * 2
                r, g, b = s, s, 1.0
            else:
                s = (t - 0.5) * 2
                r, g, b = 1.0, 1.0 - s, 1.0 - s
        
        return (r, g, b)
    
    # Build CGO object
    cgo_obj = []
    
    for point, ec_val in zip(surface_points, ec_values):
        r, g, b = ec_to_rgb(ec_val)
        cgo_obj.extend([
            COLOR, r, g, b,
            SPHERE, point[0], point[1], point[2], sphere_scale
        ])
    
    # Load CGO
    cgo_name = f"ec_cgo_{ligand_resname}"
    cmd.load_cgo(cgo_obj, cgo_name)
    
    # Show ligand sticks
    ligand_sel = f"{obj_name} and resn {ligand_resname}"
    cmd.show('sticks', ligand_sel)
    cmd.color('gray50', f"{ligand_sel} and elem C")
    
    # Zoom
    cmd.zoom(ligand_sel, buffer=8)
    
    print(f"[visualize_ec_cgo_surface] ✅ Created CGO object: {cgo_name}")
    
    return True


def create_ec_legend(position: str = 'bottom_right',
                     ec_range: tuple = (-1.0, 1.0),
                     color_scheme: str = 'rwg',
                     title: str = 'EC Score') -> bool:
    """
    在 PyMOL 3D 视图中Create CGO 色标条。

    using一系列小方块组成竖直渐变色标，并用 pseudoatom 标注 EC Value和含义。
    如果 CGO CreateFailed则 fallback 到控制台文字输出。

    Args:
        position: 色标条位置（保留Parameters，暂只支持 'bottom_right'）
        ec_range: (min, max) EC Value范围
        color_scheme: using的颜色方案 ('rwg' 或 'bwr')
        title: 色标条标题

    Returns:
        True if successful
    """
    if not PYMOL_AVAILABLE:
        return False

    ec_min, ec_max = ec_range

    # ── 尝试Create CGO 色标条 ──
    try:
        n_steps = 50  # 色标条分段数
        bar_width = 1.5   # 色条宽degrees
        bar_height = 20.0  # 色条总高degrees
        step_h = bar_height / n_steps

        # 固定屏幕位置（右侧shift）
        x_base = 25.0
        y_base = -bar_height / 2

        cgo_obj = []

        for i in range(n_steps):
            # 归一化进degrees [0, 1]
            t = i / (n_steps - 1)

            # 颜色映射：5 色渐变（深红→浅红→白→浅绿→深绿）
            if color_scheme == 'rwg':
                r, g, b = _ec_5color_interpolate(t)
            else:  # bwr
                if t < 0.5:
                    s = t * 2
                    r, g, b = s, s, 1.0
                else:
                    s = (t - 0.5) * 2
                    r, g, b = 1.0, 1.0 - s, 1.0 - s

            y0 = y_base + i * step_h
            y1 = y0 + step_h

            # 两个三角形组成一个矩形色块
            cgo_obj.extend([
                COLOR, r, g, b,
                BEGIN, TRIANGLES,
                VERTEX, x_base, y0, 0.0,
                VERTEX, x_base + bar_width, y0, 0.0,
                VERTEX, x_base + bar_width, y1, 0.0,
                VERTEX, x_base, y0, 0.0,
                VERTEX, x_base + bar_width, y1, 0.0,
                VERTEX, x_base, y1, 0.0,
                END,
            ])

        # Delete旧的色标条对象（如果存在）
        try:
            cmd.delete("ec_legend_bar")
        except Exception:
            pass

        cmd.load_cgo(cgo_obj, "ec_legend_bar")
        cmd.set("cgo_line_width", 1.0, "ec_legend_bar")

        # ── Create 5 个标注 pseudoatom（对应 5 色等级，含好坏判断说明） ──
        if color_scheme == 'rwg':
            label_items = [
                ("ec_label_1", 0.00, f"Strong clash (<{ec_min * 0.5:.1f}) Bad"),
                ("ec_label_2", 0.25, f"Mild clash ({ec_min * 0.5:.1f}~{ec_min * 0.1:.1f})"),
                ("ec_label_3", 0.50, "Neutral (~0)"),
                ("ec_label_4", 0.75, f"Mild compl. ({ec_max * 0.1:.1f}~{ec_max * 0.5:.1f})"),
                ("ec_label_5", 1.00, f"Strong compl. (>{ec_max * 0.5:.1f}) Good"),
            ]
        else:
            label_items = [
                ("ec_label_1", 0.00, f"Clash ({ec_min:.2f})"),
                ("ec_label_3", 0.50, "Neutral (0)"),
                ("ec_label_5", 1.00, f"Compl. ({ec_max:.2f})"),
            ]

        label_x = x_base + bar_width + 1.0
        for label_name, frac, text in label_items:
            y_pos = y_base + frac * bar_height
            try:
                cmd.delete(label_name)
            except Exception:
                pass
            cmd.pseudoatom(label_name, pos=[label_x, y_pos, 0.0],
                           label=text)
            cmd.set("label_size", 12, label_name)
            cmd.set("label_color", "black", label_name)
            cmd.hide("everything", label_name)
            cmd.show("labels", label_name)

        print(f"[create_ec_legend] ✅ 色标条已Create (ec_legend_bar, 5色渐变)")
        return True

    except Exception as e:
        print(f"[create_ec_legend] ⚠️ CGO 色标条CreateFailed: {e}")
        print("[create_ec_legend] 回退到文字输出...")

    # ── Fallback：纯文字输出 ──
    print(f"\n{'='*50}")
    print(f"EC Color Legend: {title}")
    print(f"{'='*50}")

    if color_scheme == 'rwg':
        print(f"🟥 深红:   EC < {ec_min * 0.5:.1f}     (强 clash → 非常不利)")
        print(f"🔴 浅红:   {ec_min * 0.5:.1f} < EC < {ec_min * 0.1:.1f} (轻微 clash → 不利)")
        print(f"⚪ 白色:   {ec_min * 0.1:.1f} < EC < {ec_max * 0.1:.1f}   (中性)")
        print(f"🟢 浅绿:   {ec_max * 0.1:.1f} < EC < {ec_max * 0.5:.1f}  (轻微互补 → 有利)")
        print(f"🟩 深绿:   EC > {ec_max * 0.5:.1f}     (强互补 → 非常有利)")
    else:
        print(f"🔵 Blue:  EC = {ec_min:.1f} (Clash → Unfavorable)")
        print(f"⚪ White: EC = 0.0 (Neutral)")
        print(f"🔴 Red:   EC = {ec_max:.1f} (Complementary → Favorable)")

    print(f"\n📊 简要：绿色=有利于结合 | 红色=不利于结合 | 白色=中性")
    print(f"{'='*50}")

    return True


def save_ec_visualization(filename: str, width: int = 2400, height: int = 2400,
                          dpi: int = 300, ray: bool = True) -> bool:
    """
    Save current EC visualization as a high-resolution image.
    
    Args:
        filename: Output filename (PNG, TIFF, etc.)
        width: Image width in pixels
        height: Image height in pixels
        dpi: DPI for the output image
        ray: Whether to ray trace before saving
        
    Returns:
        True if successful
    """
    if not PYMOL_AVAILABLE:
        return False
    
    print(f"[save_ec_visualization] Saving to {filename}...")
    
    if ray:
        print("[save_ec_visualization] Ray tracing...")
        cmd.ray(width, height)
    
    cmd.png(filename, width=width, height=height, dpi=dpi)
    
    print(f"[save_ec_visualization] ✅ Saved: {filename}")
    
    return True


def visualize_ternary_ec_surfaces(obj_name: str, glue_resname: str,
                                   ternary_result: Dict,
                                   show_overlap: bool = True) -> bool:
    """
    Visualize EC surfaces for ternary complex (molecular glue) analysis.
    
    Creates separate surfaces for each protein interface, optionally
    highlighting the overlap region where the glue bridges both proteins.
    
    Args:
        obj_name: PyMOL object name
        glue_resname: Molecular glue residue name
        ternary_result: Result from analyze_ternary_ec
        show_overlap: Highlight overlap region
        
    Returns:
        True if successful
    """
    if not PYMOL_AVAILABLE:
        return False
    
    if not NUMPY_AVAILABLE:
        return False
    
    surface_points = ternary_result.get('surface_points')
    if surface_points is None:
        print("[visualize_ternary_ec_surfaces] ❌ No surface points in result")
        return False
    
    interfaces = ternary_result.get('interfaces', {})
    
    # Visualize each interface
    for interface_key, interface_data in interfaces.items():
        if interface_key in ['A_glue', 'B_glue']:
            ec_values = interface_data.get('ec_values')
            if ec_values is not None:
                # Create CGO surface for this interface
                cgo_name = f"ec_{interface_key}"
                
                # Use different color schemes for different interfaces
                if interface_key == 'A_glue':
                    _create_interface_cgo(surface_points, ec_values, cgo_name, 
                                         base_color='cyan')
                else:
                    _create_interface_cgo(surface_points, ec_values, cgo_name,
                                         base_color='magenta')
    
    # Highlight overlap region
    if show_overlap and 'A_glue' in interfaces and 'B_glue' in interfaces:
        ec_a = interfaces['A_glue'].get('ec_values')
        ec_b = interfaces['B_glue'].get('ec_values')
        
        if ec_a is not None and ec_b is not None:
            # Find overlap points (non-zero EC for both interfaces)
            overlap_mask = (ec_a != 0) & (ec_b != 0)
            
            if np.any(overlap_mask):
                overlap_points = surface_points[overlap_mask]
                # Average EC in overlap region
                overlap_ec = (ec_a[overlap_mask] + ec_b[overlap_mask]) / 2
                
                # Create overlap CGO
                _create_interface_cgo(overlap_points, overlap_ec, 'ec_overlap',
                                     base_color='yellow', sphere_scale=0.35)
                
                print(f"[visualize_ternary_ec_surfaces] Overlap region: {len(overlap_points)} points")
    
    # Show glue molecule
    glue_sel = f"{obj_name} and resn {glue_resname}"
    cmd.show('sticks', glue_sel)
    cmd.color('orange', f"{glue_sel} and elem C")
    
    # Color proteins
    protein_a_chains = ternary_result.get('protein_a_chains', [])
    protein_b_chains = ternary_result.get('protein_b_chains', [])
    
    for chain in protein_a_chains:
        cmd.color('cyan', f"{obj_name} and chain {chain} and elem C")
    
    for chain in protein_b_chains:
        cmd.color('magenta', f"{obj_name} and chain {chain} and elem C")
    
    cmd.set('cartoon_transparency', 0.7, obj_name)
    
    # Zoom
    cmd.zoom(glue_sel, buffer=12)
    
    print("[visualize_ternary_ec_surfaces] ✅ Created ternary EC visualization")
    print("  Cyan surface: Protein A - Glue interface")
    print("  Magenta surface: Protein B - Glue interface")
    if show_overlap:
        print("  Yellow surface: Overlap (bridging) region")
    
    return True


def _create_interface_cgo(points: np.ndarray, ec_values: np.ndarray,
                          cgo_name: str, base_color: str = 'white',
                          sphere_scale: float = 0.25):
    """Create CGO for an interface with EC-based coloring.
    
    using非线性增强后的 EC Value进行着色：
    - 正 EC（互补）：向饱和绿色过渡
    - 负 EC（冲突）：向明亮红色过渡
    """
    cgo_obj = []

    # 非线性增强 EC Value
    enhanced = _enhance_ec_contrast(ec_values)

    # 基础色 RGB
    base_colors = {
        'cyan': (0.0, 1.0, 1.0),
        'magenta': (1.0, 0.0, 1.0),
        'yellow': (1.0, 1.0, 0.0),
        'white': (1.0, 1.0, 1.0)
    }

    base_rgb = base_colors.get(base_color, (1.0, 1.0, 1.0))

    for point, ec_val in zip(points, enhanced):
        # Skip零Value点
        if ec_val == 0:
            continue

        if ec_val > 0:
            # 正 EC（互补）：从基础色向饱和绿色过渡
            t = min(1.0, abs(ec_val))
            # 混合基础色和饱和绿色（增强饱和degrees）
            r = base_rgb[0] * (1 - t) + 0.0 * t
            g = base_rgb[1] * (1 - t) + 0.9 * t
            b = base_rgb[2] * (1 - t) + 0.05 * t
        else:
            # 负 EC（冲突）：从基础色向明亮红色过渡
            t = min(1.0, abs(ec_val))
            r = base_rgb[0] * (1 - t) + 1.0 * t
            g = base_rgb[1] * (1 - t) + 0.1 * t
            b = base_rgb[2] * (1 - t) + 0.05 * t

        cgo_obj.extend([
            COLOR, r, g, b,
            SPHERE, point[0], point[1], point[2], sphere_scale
        ])

    if cgo_obj:
        cmd.load_cgo(cgo_obj, cgo_name)


# Register PyMOL commands
if PYMOL_AVAILABLE:
    cmd.extend('visualize_ec_smooth_surface', visualize_ec_smooth_surface)
    cmd.extend('visualize_ec_cgo_surface', visualize_ec_cgo_surface)
    cmd.extend('save_ec_visualization', save_ec_visualization)
    cmd.extend('visualize_ternary_ec_surfaces', visualize_ternary_ec_surfaces)
    cmd.extend('create_ec_legend', create_ec_legend)


# Module info
if __name__ == '__main__':
    print("="*60)
    print("GLINT EC Visualization Module")
    print("="*60)
    print("\nThis module provides publication-quality EC surface visualization.")
    print("\nUsage in PyMOL:")
    print("  from glint.ec_visualization import visualize_ec_smooth_surface")
    print("  ")
    print("  # After running calculate_ligand_ec:")
    print("  result = calculate_ligand_ec('complex', 'LIG')")
    print("  visualize_ec_smooth_surface('complex', 'LIG', ec_result=result)")
    print("\nDependencies:")
    print(f"  NumPy: {'✅' if NUMPY_AVAILABLE else '❌'}")
    print(f"  SciPy: {'✅' if SCIPY_AVAILABLE else '❌'}")
    print(f"  PyMOL: {'✅' if PYMOL_AVAILABLE else '❌'}")