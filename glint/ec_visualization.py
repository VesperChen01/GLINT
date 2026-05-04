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


# Surface mesh generation
try:
    from .surface_similarity import SurfaceGenerator
    SURFACE_MESH_AVAILABLE = True
except Exception:
    SURFACE_MESH_AVAILABLE = False


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


def _prepare_ec_color_values(ec_values: np.ndarray,
                             ec_range: tuple = (-1.0, 1.0),
                             enhance_contrast: bool = True) -> tuple:
    """Return EC values and range used only for visualization coloring."""
    color_values = np.asarray(ec_values, dtype=float)
    if enhance_contrast:
        color_values = _enhance_ec_contrast(color_values)
    color_range = _auto_ec_range(color_values, default_range=ec_range)
    return color_values, color_range


def _ec_5color_interpolate(t: float):
    """将归一化进degrees t ∈ [0,1] 映射到 5 色渐变的 RGB Value。

    5 色锚点：
      t=0.00 → 柔和深红 [0.72, 0.25, 0.25]  (强 clash)
      t=0.25 → 柔和浅红 [0.93, 0.63, 0.63]  (轻微 clash)
      t=0.50 → 柔白   [0.98, 0.98, 0.96]  (中性)
      t=0.75 → 柔和浅绿 [0.66, 0.86, 0.66]  (轻微互补)
      t=1.00 → 柔和深绿 [0.28, 0.55, 0.30]  (强互补)

    相邻锚点之间线性插Value。

    Args:
        t: 归一化进degrees [0, 1]

    Returns:
        (r, g, b) 元组，每个分量 ∈ [0, 1]
    """
    # 中文注释：使用更柔和、论文图风格更强的红白绿锚点，
    # 降低饱和度与纯色冲击，减少当前表面的"塑料感"和硬切感。
    anchors = [
        (0.00, (0.72, 0.25, 0.25)),
        (0.25, (0.93, 0.63, 0.63)),
        (0.50, (0.98, 0.98, 0.96)),
        (0.75, (0.66, 0.86, 0.66)),
        (1.00, (0.28, 0.55, 0.30)),
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

def _ensure_ec_colors_registered():
    """确保 EC 渐变与贡献残基着色所需的命名颜色已注册。"""
    cmd.set_color('ec_dark_red', [0.72, 0.25, 0.25])
    cmd.set_color('ec_light_red', [0.93, 0.63, 0.63])
    cmd.set_color('ec_white', [0.98, 0.98, 0.96])
    cmd.set_color('ec_light_green', [0.66, 0.86, 0.66])
    cmd.set_color('ec_dark_green', [0.28, 0.55, 0.30])




def _apply_5color_ec_gradient(surface_obj: str, b_min: float, b_max: float):
    """在 PyMOL 中Apply 5 色 EC 渐变（柔和深红-柔和浅红-柔白-柔和浅绿-柔和深绿）。

    using cmd.set_color 定义 5 个自定义颜色锚点，然后用 cmd.spectrum
    的自定义色板Name进行渐变着色。比内置 red_white_green (3色) 更有区分degrees，
    能让用户直观分辨强clash、轻微clash、中性、轻微互补、强互补五个等级。

    Args:
        surface_obj: PyMOL surface 对象Name
        b_min: B-factor 最小Value（对应 EC 最小Value * 100）
        b_max: B-factor 最大Value（对应 EC 最大Value * 100）
    """
    # 中文注释：采用更柔和的论文风配色，保留红=冲突、绿=互补的语义，
    # 同时降低纯色冲击，让 surface 看起来更高级、更自然。
    cmd.set_color('ec_dark_red', [0.72, 0.25, 0.25])
    cmd.set_color('ec_light_red', [0.93, 0.63, 0.63])
    cmd.set_color('ec_white', [0.98, 0.98, 0.96])
    cmd.set_color('ec_light_green', [0.66, 0.86, 0.66])
    cmd.set_color('ec_dark_green', [0.28, 0.55, 0.30])

    # using PyMOL 自定义色板Name列表进行 spectrum 着色
    # spectrum 会在这些颜色之间线性插Value
    palette = 'ec_dark_red ec_light_red ec_white ec_light_green ec_dark_green'
    cmd.spectrum('b', palette, surface_obj, minimum=b_min, maximum=b_max)

    print(f"[_apply_5color_ec_gradient] ✅ 已Apply 5 色 EC 渐变: "
          f"深红→浅红→白→浅绿→深绿 (B: {b_min:.0f}~{b_max:.0f})")



def _build_surface_atoms_from_selection(selection: str) -> List[Dict[str, Any]]:
    """从 PyMOL selection 提取真实表面生成所需原子信息。"""
    model = cmd.get_model(selection)
    atoms = []
    fallback_radii = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
        'P': 1.80, 'S': 1.80, 'CL': 1.75, 'BR': 1.85, 'I': 1.98
    }
    for atom in model.atom:
        elem = (atom.symbol or atom.name[:1] or 'C').upper()
        radius = getattr(atom, 'vdw', 0.0) or fallback_radii.get(elem, 1.70)
        atoms.append({
            'coord': np.array(atom.coord, dtype=float),
            'element': elem,
            'radius': float(radius),
        })
    return atoms


def _create_ec_real_surface_cgo(surface_name: str,
                                ligand_sel: str,
                                surface_points: np.ndarray,
                                ec_values: np.ndarray,
                                ec_range: tuple,
                                color_scheme: str = 'rwg',
                                alpha: float = 0.82) -> bool:
    """基于真实 mesh 重建 EC surface，并用 CGO 三角面着色。"""
    if not SURFACE_MESH_AVAILABLE or not SCIPY_AVAILABLE:
        return False

    atoms = _build_surface_atoms_from_selection(ligand_sel)
    if len(atoms) < 3:
        return False

    generator = SurfaceGenerator(method='auto')
    mesh = generator.generate(atoms, probe_radius=1.2)
    if mesh is None or len(mesh.vertices) == 0 or len(mesh.faces) == 0:
        return False

    tree = cKDTree(surface_points)
    k_neighbors = min(12, len(surface_points))
    distances, indices = tree.query(mesh.vertices, k=k_neighbors)

    # 中文注释：对每个 mesh 顶点使用多个最近邻 EC 采样点做反距离加权平均，
    # 减少单点投影造成的斑点/拼贴感，让 Solid Surface 更接近连续实体表面。
    if k_neighbors == 1:
        vertex_ec = ec_values[indices]
    else:
        distances = np.asarray(distances, dtype=float)
        indices = np.asarray(indices, dtype=int)
        weights = 1.0 / np.maximum(distances, 1e-6)
        weights = weights / np.sum(weights, axis=1, keepdims=True)
        vertex_ec = np.sum(ec_values[indices] * weights, axis=1)

        # 再做一轮轻度邻域平滑，进一步降低色块噪声感
        mesh_tree = cKDTree(mesh.vertices)
        smooth_k = min(8, len(mesh.vertices))
        _, mesh_neighbors = mesh_tree.query(mesh.vertices, k=smooth_k)
        if smooth_k > 1:
            vertex_ec = np.mean(vertex_ec[mesh_neighbors], axis=1)

    ec_min, ec_max = ec_range
    denom = (ec_max - ec_min) + 1e-10

    def ec_to_rgb(ec_val: float):
        t = np.clip((ec_val - ec_min) / denom, 0.0, 1.0)
        if color_scheme == 'rwg':
            return _ec_5color_interpolate(float(t))
        if t < 0.5:
            s = t * 2
            return (s, s, 1.0)
        s = (t - 0.5) * 2
        return (1.0, 1.0 - s, 1.0 - s)

    cgo_obj = []
    for face in mesh.faces:
        cgo_obj.extend([BEGIN, TRIANGLES, ALPHA, alpha])
        for vidx in face:
            nx, ny, nz = mesh.normals[vidx] if len(mesh.normals) > vidx else (0.0, 0.0, 1.0)
            x, y, z = mesh.vertices[vidx]
            r, g, b = ec_to_rgb(vertex_ec[vidx])
            cgo_obj.extend([COLOR, r, g, b, NORMAL, nx, ny, nz, VERTEX, x, y, z])
        cgo_obj.extend([END])

    try:
        cmd.delete(surface_name)
    except Exception:
        pass
    cmd.load_cgo(cgo_obj, surface_name)
    cmd.enable(surface_name)
    return True


def _create_ec_real_mesh_cgo(mesh_name: str,
                             ligand_sel: str,
                             surface_points: np.ndarray,
                             ec_values: np.ndarray,
                             ec_range: tuple,
                             color_scheme: str = 'rwg',
                             line_width: float = 0.7,
                             edge_stride: int = 2,
                             color_softness: float = 0.35) -> bool:
    """基于表面采样点生成点状/泡状 EC 外壳。"""
    if surface_points is None or ec_values is None or len(surface_points) == 0:
        return False

    ec_min, ec_max = ec_range
    denom = (ec_max - ec_min) + 1e-10

    def ec_to_rgb(ec_val: float):
        t = np.clip((ec_val - ec_min) / denom, 0.0, 1.0)
        if color_scheme == 'rwg':
            base = _ec_5color_interpolate(float(t))
        elif t < 0.5:
            s = t * 2
            base = (s, s, 1.0)
        else:
            s = (t - 0.5) * 2
            base = (1.0, 1.0 - s, 1.0 - s)

        # 中文注释：向白色轻微混合，让点状外壳更接近参考图中的柔和电势云层。
        soft = float(np.clip(color_softness, 0.0, 0.8))
        return tuple((1.0 - soft) * c + soft * 1.0 for c in base)

    # 中文注释：规则抽稀表面点，避免外壳过厚、过密，形成更轻盈的 dotted shell。
    stride = max(int(edge_stride), 1)
    visible_points = surface_points[::stride]
    visible_ec = ec_values[::stride]

    # 中文注释：小球半径使用较小固定值，使视觉更接近"点状 mesh"，而非连续 surface。
    sphere_radius = 0.18
    cgo_obj = []
    for point, ec_val in zip(visible_points, visible_ec):
        x, y, z = [float(v) for v in point]
        r, g, b = ec_to_rgb(float(ec_val))
        cgo_obj.extend([ALPHA, 0.55, COLOR, r, g, b, SPHERE, x, y, z, sphere_radius])

    try:
        cmd.delete(mesh_name)
    except Exception:
        pass
    cmd.load_cgo(cgo_obj, mesh_name)
    cmd.enable(mesh_name)
    return True



def _ec_residue_color_name(mean_ec: float) -> str:
    """根据残基平均 EC Value 返回对应颜色名。"""
    if mean_ec <= -0.03:
        return 'ec_dark_red' if mean_ec <= -0.30 else 'ec_light_red'
    if mean_ec >= 0.03:
        return 'ec_dark_green' if mean_ec >= 0.30 else 'ec_light_green'
    return 'ec_white'


def _make_residue_selection(obj_name: str, chain: str, resi: str, resn: str) -> str:
    """Build a PyMOL residue selection that also works for blank chain IDs."""
    chain = str(chain or '').strip()
    parts = [obj_name, f"resi {resi}", f"resn {resn}"]
    if chain:
        parts.insert(1, f"chain {chain}")
    return "(" + " and ".join(parts) + ")"


def _make_atom_selection(obj_name: str, chain: str, resi: str, resn: str, atom_name: str) -> str:
    """Build a PyMOL atom selection that also works for blank chain IDs."""
    residue_sel = _make_residue_selection(obj_name, chain, resi, resn)
    return f"({residue_sel} and name {atom_name})"


def _interaction_distance_limit(interaction_type: str) -> float:
    """Return a strict visualization distance cutoff for a detected interaction."""
    limits = {
        'Hydrogen Bond': 3.8,
        'Salt Bridge': 5.0,
        'Hydrophobic': 4.5,
        'Halogen Bond': 4.2,
        'Metal Coordination': 3.5,
        'Water Bridge': 4.0,
    }
    return limits.get(interaction_type, 5.0)


def _interaction_coords(inter: Dict[str, Any]):
    """Return ligand/protein coordinates stored in an interaction result."""
    try:
        lig = (
            float(inter['Ligand_Atom_X']),
            float(inter['Ligand_Atom_Y']),
            float(inter['Ligand_Atom_Z']),
        )
        prot = (
            float(inter['Protein_Atom_X']),
            float(inter['Protein_Atom_Y']),
            float(inter['Protein_Atom_Z']),
        )
        return lig, prot
    except (KeyError, TypeError, ValueError):
        return None, None


def _get_protein_ligand_interaction_map(obj_name: str,
                                        ligand_resname: str) -> Dict[tuple, List[str]]:
    """调用现有蛋白-配体相互作用模块，返回按残基索引的相互作用类型映射。"""
    try:
        from .interaction_analyzer import analyze_protein_ligand_interactions
    except Exception as e:
        print(f"[EC Residues] ⚠️ 无法导入 protein-ligand interaction 模块: {e}")
        return {}

    try:
        result = analyze_protein_ligand_interactions(
            obj_name=obj_name,
            ligand_resname=ligand_resname,
            distance_cutoff=4.5,
            key_interactions_only=False,
        )
    except Exception as e:
        print(f"[EC Residues] ⚠️ protein-ligand interaction 分析失败: {e}")
        return {}

    interaction_map: Dict[tuple, List[str]] = {}
    for inter in (result or {}).get('interactions', []):
        chain = str(inter.get('Protein_Chain', '')).strip()
        residue_text = str(inter.get('Protein_Residue', '')).strip()
        interaction_type = str(inter.get('Interaction', '')).strip()
        if not chain or not residue_text or not interaction_type:
            continue
        parts = residue_text.split()
        if len(parts) < 2:
            continue
        resn, resi = parts[0], parts[1]
        key = (chain, resi, resn)
        interaction_map.setdefault(key, [])
        if interaction_type not in interaction_map[key]:
            interaction_map[key].append(interaction_type)
    return interaction_map

def _draw_ec_residue_interaction_dashes(obj_name: str,
                                       ligand_resname: str,
                                       red_residue_keys: List[tuple]):
    """为红色冲突残基绘制与配体之间的相互作用虚线，不在 label 中显示 interaction 文本。"""
    if not ligand_resname or not red_residue_keys:
        return

    try:
        from .interaction_analyzer import analyze_protein_ligand_interactions, apply_interaction_dash_style
    except Exception as e:
        print(f"[EC Residues] ⚠️ 无法导入 interaction 可视化模块: {e}")
        return

    try:
        result = analyze_protein_ligand_interactions(
            obj_name=obj_name,
            ligand_resname=ligand_resname,
            distance_cutoff=4.5,
            key_interactions_only=False,
        )
    except Exception as e:
        print(f"[EC Residues] ⚠️ 相互作用分析失败，无法绘制 dash: {e}")
        return

    color_map = {
        'Hydrogen Bond': 'glue_hbond',
        'Salt Bridge': 'orange',
        'Hydrophobic': 'gray50',
        'Pi-Pi': 'violet',
        'Pi-Pi Stacking': 'violet',
        'Pi-Cation': 'magenta',
        'Halogen Bond': 'cyan',
        'Metal Coordination': 'yellow',
        'Water Bridge': 'teal',
    }

    red_residue_keys = set(red_residue_keys)
    drawn = 0
    skipped_long = 0
    skipped_ambiguous = 0
    for i, inter in enumerate((result or {}).get('interactions', []), 1):
        chain = str(inter.get('Protein_Chain', '')).strip()
        residue_text = str(inter.get('Protein_Residue', '')).strip()
        prot_atom = str(inter.get('Protein_Atom', '')).strip()
        lig_chain = str(inter.get('Ligand_Chain', '')).strip()
        lig_residue = str(inter.get('Ligand_Residue', '')).strip()
        lig_atom = str(inter.get('Ligand_Atom', '')).strip()
        interaction_type = str(inter.get('Interaction', '')).strip()
        parts = residue_text.split()
        lig_parts = lig_residue.split()
        if len(parts) < 2 or len(lig_parts) < 2:
            continue
        resn, resi = parts[0], parts[1]
        if (chain, resi, resn) not in red_residue_keys:
            continue

        if not prot_atom or not lig_atom or prot_atom.startswith('Ring(') or lig_atom.startswith('Ring('):
            skipped_ambiguous += 1
            continue

        # 清理配体 resn/resi 中的非法字符（如尖括号）
        lig_resn_clean = lig_parts[0].replace('<', '').replace('>', '')
        lig_resi_clean = lig_parts[1].replace('<', '').replace('>', '')
        lig_coords, prot_coords = _interaction_coords(inter)
        if lig_coords is None or prot_coords is None:
            skipped_ambiguous += 1
            continue

        dash_name = f"ec_interact_{i}"
        prot_pseudo = f"{dash_name}_prot"
        lig_pseudo = f"{dash_name}_lig"
        try:
            actual_distance = float(np.linalg.norm(np.asarray(prot_coords) - np.asarray(lig_coords)))
            max_distance = _interaction_distance_limit(interaction_type)
            reported_distance = inter.get('Distance', None)
            try:
                reported_distance = float(reported_distance)
            except (TypeError, ValueError):
                reported_distance = actual_distance

            if actual_distance > max_distance or reported_distance > max_distance:
                skipped_long += 1
                print(
                    f"[EC Residues] ⚠️ Skip long/invalid dash {interaction_type}: "
                    f"{resn}{resi}:{prot_atom} - {lig_resn_clean}{lig_resi_clean}:{lig_atom}, "
                    f"actual={actual_distance:.2f} Å, reported={reported_distance:.2f} Å"
                )
                continue

            cmd.pseudoatom(prot_pseudo, pos=prot_coords, name='ECI')
            cmd.pseudoatom(lig_pseudo, pos=lig_coords, name='ECI')
            cmd.hide('everything', prot_pseudo)
            cmd.hide('everything', lig_pseudo)
            cmd.distance(dash_name, prot_pseudo, lig_pseudo, mode=0)
            cmd.enable(dash_name)
            apply_interaction_dash_style(
                cmd,
                dash_name,
                color_name=color_map.get(interaction_type, 'gray60'),
                dash_width=2.2,
                dash_gap=0.28,
                dash_length=0.22,
                dash_radius=0.08,
                hide_labels=True,
            )
            cmd.show('dashes', dash_name)
            drawn += 1
        except Exception:
            try:
                cmd.delete(prot_pseudo)
                cmd.delete(lig_pseudo)
            except Exception:
                pass
            continue

    print(
        f"[EC Residues] ✅ 已绘制 {drawn} 条红色残基相关相互作用 dash"
        f" (跳过 ambiguous={skipped_ambiguous}, long={skipped_long})"
    )
def _cleanup_ec_interaction_objects():
    """清理旧的 EC interaction dash/label 对象，避免重复叠加。"""
    try:
        for name in cmd.get_names('all'):
            if name.startswith('ec_interact_') or name.startswith('ec_residue_label_'):
                try:
                    cmd.delete(name)
                except Exception:
                    pass
    except Exception:
        pass








def _show_ec_contributing_residues(obj_name: str,
                                   ligand_sel: str,
                                   surface_points: np.ndarray,
                                   ec_values: np.ndarray,
                                   protein_distance: float = 5.0,
                                   assignment_cutoff: float = 6.5,
                                   label_top_n: int = 12):
    """显示对 EC 有局部贡献的附近蛋白残基，并按冲突/互补着色与标注。"""
    if surface_points is None or ec_values is None or len(surface_points) == 0:
        return

    # 清理 ligand_sel 中解析出的 resn 尖括号
    ligand_resname_raw = ligand_sel.split('resn', 1)[1].strip() if 'resn' in ligand_sel else ''
    ligand_resname = ligand_resname_raw.replace('<', '').replace('>', '') if ligand_resname_raw else ''
    interaction_map = _get_protein_ligand_interaction_map(obj_name, ligand_resname) if ligand_resname else {}

    try:
        cmd.delete('ec_contrib_residues')
    except Exception:
        pass

    # 中文注释：EC residue 视图只允许红色冲突残基以 sticks 显示。
    # 先隐藏所有蛋白 sticks/licorice，避免全蛋白背景紫色或历史显示状态漏出来。
    try:
        cmd.hide('sticks', f"({obj_name} and polymer)")
        cmd.hide('licorice', f"({obj_name} and polymer)")
    except Exception:
        pass

    protein_sel = f"({obj_name} and polymer within {protein_distance} of ({ligand_sel}))"
    model = cmd.get_model(protein_sel)
    if not model.atom:
        print('[EC Residues] ⚠️ 未找到配体附近蛋白原子，跳过贡献残基显示')
        return

    atom_coords = []
    atom_meta = []
    for atom in model.atom:
        coord = getattr(atom, 'coord', None)
        if coord is None or len(coord) != 3:
            continue
        atom_coords.append(coord)
        atom_meta.append((atom.chain, atom.resi, atom.resn))


    # 中文注释：贡献残基着色前确保 EC 命名颜色可用，避免 Unknown color 异常。
    _ensure_ec_colors_registered()

    if not atom_coords:
        return

    atom_coords = np.asarray(atom_coords, dtype=float)
    tree = cKDTree(atom_coords)
    distances, nearest_idx = tree.query(surface_points, k=1)

    residue_scores = {}
    for i, dist in enumerate(np.asarray(distances, dtype=float)):
        if dist > assignment_cutoff:
            continue
        key = atom_meta[int(nearest_idx[i])]
        residue_scores.setdefault(key, []).append(float(ec_values[i]))

    if not residue_scores:
        print('[EC Residues] ⚠️ 表面点未直接分配到残基，改用附近蛋白原子反向匹配')
        surface_tree = cKDTree(surface_points)
        atom_to_surface_dist, atom_to_surface_idx = surface_tree.query(atom_coords, k=1)
        reverse_cutoff = max(float(assignment_cutoff), float(protein_distance) + 2.0)
        for atom_i, dist in enumerate(np.asarray(atom_to_surface_dist, dtype=float)):
            if dist > reverse_cutoff:
                continue
            key = atom_meta[atom_i]
            residue_scores.setdefault(key, []).append(float(ec_values[int(atom_to_surface_idx[atom_i])]))

    if not residue_scores:
        print('[EC Residues] ⚠️ 未分配到有效贡献残基，跳过显示')
        return

    ranked = []
    for key, values in residue_scores.items():
        mean_ec = float(np.mean(values))
        strength = float(np.mean(np.abs(values)))
        ranked.append((key, mean_ec, strength, len(values)))
    ranked.sort(key=lambda x: x[2], reverse=True)

    label_count = 0
    shown_count = 0
    red_residue_keys = []
    for (chain, resi, resn), mean_ec, strength, n_points in ranked:
        # 中文注释：EC 残基提示只突出红色冲突残基；绿色互补已由 EC surface 表达。
        if mean_ec > -0.03:
            continue
        residue_sel = _make_residue_selection(obj_name, chain, resi, resn)
        color_name = _ec_residue_color_name(mean_ec)
        cmd.show('sticks', residue_sel)
        cmd.set('stick_radius', 0.22, residue_sel)
        cmd.color(color_name, f"{residue_sel} and elem C")
        cmd.color('blue', f"{residue_sel} and elem N")
        cmd.color('red', f"{residue_sel} and elem O")
        cmd.color('yellow', f"{residue_sel} and elem S")
        cmd.color('green', f"{residue_sel} and elem Cl")
        cmd.color('orange', f"{residue_sel} and elem Br")
        shown_count += 1
        red_residue_keys.append((chain, resi, resn))
        if label_count < label_top_n and n_points >= 3:
            label_sel = f"{residue_sel} and name CA"
            if cmd.count_atoms(label_sel) == 0:
                label_sel = residue_sel
            cmd.label(label_sel, f'"{resn}{resi} ({mean_ec:+.2f})"')
            cmd.set('label_color', color_name, label_sel)
            label_count += 1

    if ligand_resname:
        _draw_ec_residue_interaction_dashes(obj_name, ligand_resname, red_residue_keys)

    cmd.set('label_size', 14)
    print(f"[EC Residues] ✅ 已显示 {shown_count} 个红色 EC 冲突残基，标注前 {label_count} 个")


def visualize_ec_smooth_surface(obj_name: str, ligand_resname: str,
                                 ec_values: np.ndarray = None,
                                 surface_points: np.ndarray = None,
                                 ec_result: Dict = None,
                                 surface_type: str = 'molecular',
                                 transparency: float = 0.30,
                                 show_ligand_sticks: bool = True,
                                 show_protein_lines: bool = False,
                                 protein_distance: float = 5.0,
                                 color_scheme: str = 'rwg',
                                 ec_range: tuple = (-1.0, 1.0),
                                 surface_quality: int = 4,
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

    # 自适应颜色范围：根据实际数据分布调整。这里的增强只用于显示，不改变 EC 分数。
    color_ec_values, color_ec_range = _prepare_ec_color_values(
        ec_values,
        ec_range=ec_range,
        enhance_contrast=enhance_contrast,
    )

    ligand_sel = f"{obj_name} and resn {ligand_resname}"

    # 中文注释：清理旧 legend，避免 3D 图例对象干扰主分子视图。
    for legend_obj in [
        "ec_legend_bar", "ec_legend_title",
        "ec_label_1", "ec_label_2", "ec_label_3", "ec_label_4", "ec_label_5"
    ]:
        try:
            cmd.delete(legend_obj)
        except Exception:
            pass

    surface_obj = f"ec_surface_{ligand_resname}"
    mesh_obj = f"ec_mesh_{ligand_resname}"

    # 中文注释：切换到 solid 时顺便清理旧 mesh 对象，避免样式切换后彼此覆盖。
    try:
        cmd.delete(mesh_obj)
    except Exception:
        pass


    # 中文注释：先删除旧对象，避免重复创建后的显示状态残留。
    try:
        cmd.delete(surface_obj)
    except Exception:
        pass

    # Step 1: Prefer a real per-vertex CGO molecular surface. This preserves local
    # EC color variation instead of averaging many surface samples onto a few atoms.
    if surface_type == 'molecular' and _create_ec_real_surface_cgo(
        surface_name=surface_obj,
        ligand_sel=ligand_sel,
        surface_points=surface_points,
        ec_values=color_ec_values,
        ec_range=color_ec_range,
        color_scheme=color_scheme,
        alpha=max(0.0, min(1.0, 1.0 - float(transparency))),
    ):
        if show_ligand_sticks:
            cmd.show('sticks', ligand_sel)
            cmd.set('stick_radius', 0.18, ligand_sel)
            cmd.color('gray70', f"{ligand_sel} and elem C")
            cmd.color('blue', f"{ligand_sel} and elem N")
            cmd.color('red', f"{ligand_sel} and elem O")
            cmd.color('yellow', f"{ligand_sel} and elem S")
            cmd.color('green', f"{ligand_sel} and elem Cl")
            cmd.color('orange', f"{ligand_sel} and elem Br")
            cmd.set('stick_transparency', 0.0, ligand_sel)

        cmd.set_color('glint_protein_purple', [0.67, 0.55, 0.86])
        if show_protein_lines:
            protein_sel = f"{obj_name} and polymer within {protein_distance} of {ligand_sel}"
            cmd.show('lines', protein_sel)
            cmd.color('glint_protein_purple', f"{protein_sel} and elem C")
            cmd.set('line_width', 1.5, protein_sel)

        cmd.hide('cartoon', obj_name)
        cmd.color('glint_protein_purple', f"{obj_name} and polymer")
        _show_ec_contributing_residues(
            obj_name=obj_name,
            ligand_sel=ligand_sel,
            surface_points=surface_points,
            ec_values=ec_values,
            protein_distance=protein_distance,
        )
        _set_publication_rendering()
        cmd.zoom(ligand_sel, buffer=8)

        if ray_trace:
            print("[visualize_ec_smooth_surface] Ray tracing...")
            cmd.ray()

        print(f"[visualize_ec_smooth_surface] ✅ Created surface: {surface_obj}")
        return True

    # Step 1b: Create a copy of the ligand for PyMOL built-in surface generation
    cmd.create(surface_obj, ligand_sel)
    cmd.enable(surface_obj)

    # Step 2: Map EC values to atom B-factors（含非线性增强）
    if SCIPY_AVAILABLE:
        _map_ec_to_bfactors_kdtree(surface_obj, surface_points, color_ec_values,
                                    enhance_contrast=False)
    else:
        _map_ec_to_bfactors_simple(surface_obj, surface_points, color_ec_values,
                                    enhance_contrast=False)

    # Step 3: Rebuild to apply B-factor changes
    cmd.rebuild(surface_obj)

    # Step 4: Generate and configure surface
    cmd.hide('everything', surface_obj)
    cmd.show('surface', surface_obj)
    cmd.set('surface_mode', 0, surface_obj)
    cmd.set('two_sided_lighting', 1, surface_obj)

    # 中文注释：solid 默认改为更接近参考图的 molecular surface，
    # 并适当降低透明度，减少 PyMOL 视口中的点阵/抖动感。
    surface_type_map = {
        'molecular': 0,
        'solvent': 1,
        'gaussian': 2
    }
    cmd.set('surface_type', surface_type_map.get(surface_type, 2), surface_obj)

    # 中文注释：solid 路线只保留平滑半透明 surface，不使用 mesh/三角壳逻辑。
    cmd.set('surface_quality', max(surface_quality, 4), surface_obj)
    cmd.set('surface_smooth_edges', 1, surface_obj)
    cmd.set('surface_proximity', 1, surface_obj)

    # Step 5: Apply color spectrum based on B-factors
    ec_min, ec_max = color_ec_range
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
    cmd.show('surface', surface_obj)
    cmd.set('transparency', transparency, surface_obj)
    cmd.set('surface_color_smoothing', 1)
    cmd.set('surface_color_smoothing_threshold', 0.22)
    cmd.set('surface_quality', max(surface_quality, 4), surface_obj)

    # 中文注释：solid 参考图里 sticks 是清楚但不过粗，因此适当减细并弱化碳颜色。
    if show_ligand_sticks:
        cmd.show('sticks', ligand_sel)
        cmd.set('stick_radius', 0.18, ligand_sel)
        cmd.color('gray70', f"{ligand_sel} and elem C")
        cmd.color('blue', f"{ligand_sel} and elem N")
        cmd.color('red', f"{ligand_sel} and elem O")
        cmd.color('yellow', f"{ligand_sel} and elem S")
        cmd.color('green', f"{ligand_sel} and elem Cl")
        cmd.color('orange', f"{ligand_sel} and elem Br")

    # Step 7: Show ligand sticks inside surface
    if show_ligand_sticks:
        cmd.show('sticks', ligand_sel)
        cmd.set('stick_transparency', 0.0, ligand_sel)


    # 中文注释：在实际着色前先确保自定义紫色已注册，避免 Unknown color 异常。
    cmd.set_color('glint_protein_purple', [0.67, 0.55, 0.86])

    # Step 8: Show nearby protein residues
    if show_protein_lines:
        protein_sel = f"{obj_name} and polymer within {protein_distance} of {ligand_sel}"
        cmd.show('lines', protein_sel)
        cmd.color('glint_protein_purple', f"{protein_sel} and elem C")
        cmd.set('line_width', 1.5, protein_sel)

    # Step 9: EC 视图默认隐藏蛋白 cartoon，只保留局部线框/残基高亮与 EC 表面。
    cmd.hide('cartoon', obj_name)
    cmd.color('glint_protein_purple', f"{obj_name} and polymer")
    _show_ec_contributing_residues(
        obj_name=obj_name,
        ligand_sel=ligand_sel,
        surface_points=surface_points,
        ec_values=ec_values,
        protein_distance=protein_distance,
    )

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

    # 中文注释：使用更适合半透明实体表面的显示参数，避免视口里出现明显点阵/抖动感。
    cmd.set('surface_smooth_edges', 1)
    cmd.set('two_sided_lighting', 1)
    cmd.set('transparency_mode', 2)

    # Background
    cmd.bg_color('white')

    # 中文注释：将蛋白主体统一设为柔和紫色，提升和红绿 EC surface 的区分度。
    cmd.set_color('glint_protein_purple', [0.67, 0.55, 0.86])

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

    # 自适应颜色范围。增强只用于显示，不改变 EC 分数。
    color_ec_values, color_ec_range = _prepare_ec_color_values(
        ec_values,
        ec_range=ec_range,
        enhance_contrast=enhance_contrast,
    )

    ec_values = color_ec_values
    ec_min, ec_max = color_ec_range

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

    using一系列小方块组成水平渐变色标，并用 pseudoatom 标注 EC Value和含义。
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
        n_steps = 60  # 色标条分段数
        bar_width = 18.0   # 水平色条总宽度
        bar_height = 1.2   # 水平色条高度
        step_w = bar_width / n_steps

        # 中文注释：改为水平图例并向右下角偏移，降低透视倾斜和文字裁切问题。
        x_base = 18.0
        y_base = -12.0
        z_base = 0.0

        cgo_obj = []

        for i in range(n_steps):
            t = i / (n_steps - 1)

            if color_scheme == 'rwg':
                r, g, b = _ec_5color_interpolate(t)
            else:
                if t < 0.5:
                    s = t * 2
                    r, g, b = s, s, 1.0
                else:
                    s = (t - 0.5) * 2
                    r, g, b = 1.0, 1.0 - s, 1.0 - s

            x0 = x_base + i * step_w
            x1 = x0 + step_w
            y0 = y_base
            y1 = y_base + bar_height

            cgo_obj.extend([
                COLOR, r, g, b,
                BEGIN, TRIANGLES,
                VERTEX, x0, y0, z_base,
                VERTEX, x1, y0, z_base,
                VERTEX, x1, y1, z_base,
                VERTEX, x0, y0, z_base,
                VERTEX, x1, y1, z_base,
                VERTEX, x0, y1, z_base,
                END,
            ])

        try:
            cmd.delete("ec_legend_bar")
            cmd.delete("ec_label_1")
            cmd.delete("ec_label_2")
            cmd.delete("ec_label_3")
            cmd.delete("ec_label_4")
            cmd.delete("ec_label_5")
            cmd.delete("ec_legend_title")
        except Exception:
            pass

        cmd.load_cgo(cgo_obj, "ec_legend_bar")
        cmd.set("cgo_line_width", 1.0, "ec_legend_bar")

        if color_scheme == 'rwg':
            label_items = [
                ("ec_label_1", 0.00, f"Strong clash (<{ec_min * 0.5:.1f})"),
                ("ec_label_2", 0.25, f"Mild clash ({ec_min * 0.5:.1f}~{ec_min * 0.1:.1f})"),
                ("ec_label_3", 0.50, "Neutral (~0)"),
                ("ec_label_4", 0.75, f"Mild compl. ({ec_max * 0.1:.1f}~{ec_max * 0.5:.1f})"),
                ("ec_label_5", 1.00, f"Strong compl. (>{ec_max * 0.5:.1f})"),
            ]
        else:
            label_items = [
                ("ec_label_1", 0.00, f"Clash ({ec_min:.2f})"),
                ("ec_label_3", 0.50, "Neutral (0)"),
                ("ec_label_5", 1.00, f"Compl. ({ec_max:.2f})"),
            ]

        cmd.pseudoatom("ec_legend_title", pos=[x_base + bar_width / 2, y_base + 2.6, z_base], label=title)
        cmd.set("label_size", 14, "ec_legend_title")
        cmd.set("label_color", "black", "ec_legend_title")
        cmd.hide("everything", "ec_legend_title")
        cmd.show("labels", "ec_legend_title")

        for label_name, frac, text in label_items:
            x_pos = x_base + frac * bar_width
            y_pos = y_base - 1.8
            cmd.pseudoatom(label_name, pos=[x_pos, y_pos, z_base], label=text)
            cmd.set("label_size", 11, label_name)
            cmd.set("label_color", "black", label_name)
            cmd.set("label_font_id", 7, label_name)
            cmd.hide("everything", label_name)
            cmd.show("labels", label_name)

        print(f"[create_ec_legend] ✅ 色标条已Create (ec_legend_bar, 水平5色渐变)")
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

    cmd.hide('cartoon', obj_name)

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




# ============================================================================
# 双面板 FMO 风格可视化函数（实体表面 + 半透明 mesh 对比展示）
# ============================================================================

def _calculate_mol_extent(selection: str) -> float:
    """计算分子选区的包围盒最大跨度（X 方向），用于双面板平移距离计算。

    Args:
        selection: PyMOL 选区字符串

    Returns:
        X 方向跨度（Å）
    """
    coords = []
    cmd.iterate_state(1, selection,
                     "coords.append(x)",
                     space={'coords': coords})
    if not coords:
        return 20.0  # 默认值
    return max(coords) - min(coords)


def _apply_sticks_coloring(selection: str, carbon_color: str = 'green'):
    """为分子骨架设置元素着色（碳用指定颜色，其他元素用标准着色）。

    模仿参考图中绿色碳 + 元素着色的风格。

    Args:
        selection: PyMOL 选区字符串
        carbon_color: 碳原子颜色（默认 green，模仿参考图风格）
    """
    cmd.show('sticks', selection)
    cmd.color(carbon_color, f"{selection} and elem C")
    cmd.color('blue', f"{selection} and elem N")
    cmd.color('red', f"{selection} and elem O")
    cmd.color('yellow', f"{selection} and elem S")
    cmd.color('green', f"{selection} and elem Cl")
    cmd.color('orange', f"{selection} and elem Br")
    cmd.set('stick_radius', 0.15, selection)


def _add_ec_stat_labels(ec_values: np.ndarray, center_pos: list,
                         label_prefix: str = '',
                         y_offset: float = 0.0) -> List[str]:
    """在 PyMOL 视图中用 pseudoatom 标注 EC 统计信息。

    类似参考图中 E_HOMO / E_LUMO 的文字标注风格，显示关键 EC 统计值。

    Args:
        ec_values: EC 值数组
        center_pos: 标注中心位置 [x, y, z]
        label_prefix: 标注名称前缀（避免多次调用冲突）
        y_offset: Y 方向偏移量

    Returns:
        创建的 pseudoatom 名称列表
    """
    if not NUMPY_AVAILABLE:
        return []

    # 计算统计信息
    ec_mean = float(np.mean(ec_values))
    ec_median = float(np.median(ec_values))
    ec_positive_ratio = float(np.sum(ec_values > 0) / len(ec_values) * 100)

    # 标注项：(名称后缀, Y 偏移, 标注文本)
    label_items = [
        ('mean',   2.0, f"EC_mean = {ec_mean:.3f}"),
        ('median', 0.0, f"EC_median = {ec_median:.3f}"),
        ('ratio', -2.0, f"Positive = {ec_positive_ratio:.1f}%"),
    ]

    created_labels = []
    x, y, z = center_pos[0], center_pos[1] + y_offset, center_pos[2]

    for suffix, dy, text in label_items:
        label_name = f"ec_stat_{label_prefix}{suffix}"
        # 清理已存在的同名对象
        try:
            cmd.delete(label_name)
        except Exception:
            pass

        cmd.pseudoatom(label_name, pos=[x, y + dy, z], label=text)
        cmd.set('label_size', 14, label_name)
        cmd.set('label_color', 'black', label_name)
        cmd.set('label_font_id', 7, label_name)  # 等宽字体
        cmd.hide('everything', label_name)
        cmd.show('labels', label_name)
        created_labels.append(label_name)

    return created_labels


def visualize_ec_mesh_surface(obj_name: str, ligand_resname: str,
                               ec_values: np.ndarray = None,
                               surface_points: np.ndarray = None,
                               ec_result: Dict = None,
                               mesh_width: float = 0.55,
                               surface_type: str = 'gaussian',
                               show_ligand_sticks: bool = True,
                               show_protein_lines: bool = False,
                               protein_distance: float = 5.0,
                               color_scheme: str = 'rwg',
                               ec_range: tuple = (-1.0, 1.0),
                               surface_quality: int = 1,
                               ray_trace: bool = False,
                               enhance_contrast: bool = True) -> bool:
    """创建 EC 着色的半透明 mesh（网格）表面可视化。

    类似参考图右侧 Electrostatic Potential 风格：蓝/红/白半透明网格
    包裹分子，分子骨架清晰可见。

    Args:
        obj_name: PyMOL 对象名
        ligand_resname: 配体残基名（如 'LIG'）
        ec_values: EC 值数组（来自 calculate_ligand_ec）
        surface_points: 表面点坐标数组
        ec_result: calculate_ligand_ec 的完整结果字典（替代输入）
        mesh_width: mesh 线宽（默认 1.0）
        surface_type: 表面类型 - 'gaussian', 'solvent', 'molecular'
        show_ligand_sticks: 是否显示配体 sticks 模型
        show_protein_lines: 是否显示蛋白 lines
        protein_distance: 蛋白显示距离截断（Å）
        color_scheme: 颜色方案 - 'rwg'（红白绿）或 'bwr'（蓝白红）
        ec_range: (min, max) EC 值映射范围
        surface_quality: 表面质量（0-4）
        ray_trace: 是否进行光线追踪
        enhance_contrast: 是否启用非线性对比度增强

    Returns:
        True if successful
    """
    if not PYMOL_AVAILABLE:
        print("[visualize_ec_mesh_surface] ❌ PyMOL required")
        return False
    if not NUMPY_AVAILABLE:
        print("[visualize_ec_mesh_surface] ❌ NumPy required")
        return False

    # 从 ec_result 提取数据
    if ec_result is not None:
        ec_values = ec_result.get('ec_values', ec_values)
        surface_points = ec_result.get('surface_points', surface_points)

    # 中文注释：清理旧 legend，避免 3D 图例对象干扰主分子视图。
    for legend_obj in [
        "ec_legend_bar", "ec_legend_title",
        "ec_label_1", "ec_label_2", "ec_label_3", "ec_label_4", "ec_label_5"
    ]:
        try:
            cmd.delete(legend_obj)
        except Exception:
            pass


    if ec_values is None or surface_points is None:
        print("[visualize_ec_mesh_surface] ❌ Need ec_values and surface_points")
        return False
    surface_obj = f"ec_surface_{ligand_resname}"

    # 中文注释：切换到 mesh 时顺便清理旧 solid surface 对象，避免旧壳层残留覆盖 mesh。
    try:
        cmd.delete(surface_obj)
    except Exception:
        pass


    print(f"[visualize_ec_mesh_surface] 创建 EC mesh 表面...")

    # 自适应颜色范围。增强只用于显示，不改变 EC 分数。
    color_ec_values, color_ec_range = _prepare_ec_color_values(
        ec_values,
        ec_range=ec_range,
        enhance_contrast=enhance_contrast,
    )

    ligand_sel = f"{obj_name} and resn {ligand_resname}"
    mesh_obj = f"ec_mesh_{ligand_resname}"

    # 中文注释：删除旧 mesh 对象，避免历史显示状态干扰本次绘制。
    try:
        cmd.delete(mesh_obj)
    except Exception:
        pass

    # 创建配体副本用于 mesh 表面
    cmd.create(mesh_obj, ligand_sel)
    cmd.enable(mesh_obj)

    # 映射 EC 值到 B-factor
    if SCIPY_AVAILABLE:
        _map_ec_to_bfactors_kdtree(mesh_obj, surface_points, color_ec_values,
                                    enhance_contrast=False)
    else:
        _map_ec_to_bfactors_simple(mesh_obj, surface_points, color_ec_values,
                                    enhance_contrast=False)

    cmd.rebuild(mesh_obj)

    # 中文注释：优先使用真实 CGO 线框 mesh，确保在当前图形环境下也能稳定看到外壳。
    if _create_ec_real_mesh_cgo(
        mesh_name=mesh_obj,
        ligand_sel=ligand_sel,
        surface_points=surface_points,
        ec_values=color_ec_values,
        ec_range=color_ec_range,
        color_scheme=color_scheme,
        line_width=max(mesh_width, 0.7),
        edge_stride=3,
        color_softness=0.40,
    ):
        _set_publication_rendering()
        cmd.set_color('glint_protein_purple', [0.67, 0.55, 0.86])
        if show_ligand_sticks:
            _apply_sticks_coloring(ligand_sel, carbon_color='green')
            cmd.set('stick_radius', 0.16, ligand_sel)
        if show_protein_lines:
            protein_sel = f"{obj_name} and polymer within {protein_distance} of {ligand_sel}"
            cmd.show('lines', protein_sel)
            cmd.color('glint_protein_purple', f"{protein_sel} and elem C")
            cmd.set('line_width', 1.5, protein_sel)
        cmd.hide('cartoon', obj_name)
        cmd.color('glint_protein_purple', f"{obj_name} and polymer")
        _show_ec_contributing_residues(
            obj_name=obj_name,
            ligand_sel=ligand_sel,
            surface_points=surface_points,
            ec_values=ec_values,
            protein_distance=protein_distance,
        )
        cmd.zoom(ligand_sel, buffer=8)

        print(f"[visualize_ec_mesh_surface] ✅ Created mesh surface: {mesh_obj}")
        return True

    # Fallback: use PyMOL's built-in mesh representation
    cmd.hide('everything', mesh_obj)
    cmd.show('mesh', mesh_obj)
    cmd.set('mesh_width', mesh_width, mesh_obj)

    # Apply coloring
    ec_min, ec_max = color_ec_range
    b_min = ec_min * 100
    b_max = ec_max * 100

    if color_scheme == 'rwg':
        _apply_5color_ec_gradient(mesh_obj, b_min, b_max)
    elif color_scheme == 'bwr':
        cmd.spectrum('b', 'blue_white_red', mesh_obj, minimum=b_min, maximum=b_max)
    else:
        _apply_5color_ec_gradient(mesh_obj, b_min, b_max)

    cmd.set('transparency', 0.5, mesh_obj)
    cmd.set('surface_color_smoothing', 1)
    cmd.set('surface_color_smoothing_threshold', 0.22)

    if show_ligand_sticks:
        _apply_sticks_coloring(ligand_sel, carbon_color='green')
        cmd.set('stick_transparency', 0.0, ligand_sel)

    cmd.set_color('glint_protein_purple', [0.67, 0.55, 0.86])

    if show_protein_lines:
        protein_sel = f"{obj_name} and polymer within {protein_distance} of {ligand_sel}"
        cmd.show('lines', protein_sel)
        cmd.color('glint_protein_purple', f"{protein_sel} and elem C")
        cmd.set('line_width', 1.5, protein_sel)

    cmd.hide('cartoon', obj_name)
    cmd.color('glint_protein_purple', f"{obj_name} and polymer")
    _show_ec_contributing_residues(
        obj_name=obj_name,
        ligand_sel=ligand_sel,
        surface_points=surface_points,
        ec_values=ec_values,
        protein_distance=protein_distance,
    )

    _set_publication_rendering()
    cmd.zoom(ligand_sel, buffer=8)

    if ray_trace:
        print("[visualize_ec_mesh_surface] Ray tracing...")
        cmd.ray()

    print(f"[visualize_ec_mesh_surface] ✅ Created mesh surface: {mesh_obj}")
    return True
