# -*- coding: utf-8 -*-
"""
color_scheme.py
Unified color scheme for GLINT - Schrödinger-inspired professional palette

This module provides a consistent color scheme across all GLINT modules including:
- 2D interaction diagrams
- 3D PyMOL visualizations
- PPI analysis
- GUI components

All colors are defined in multiple formats for compatibility:
- Hex codes (for matplotlib, 2D plotting)
- RGB tuples 0-1 range (for PyMOL)
- RGB tuples 0-255 range (for Qt/GUI)
"""

from dataclasses import dataclass
from typing import Dict, Optional


@dataclass
class VisualizationSettings:
    """统一的3D可视化配置对象（用于 PPI / 蛋白-配体 渲染）。"""
    # 通用渲染参数
    background: str = 'white'
    label_size: int = 16
    label_font_id: int = 5

    # 交互显示参数
    display_mode: str = 'cartoon_surface_interaction'
    show_labels: bool = True
    show_hydrophobic: bool = False
    min_confidence: float = 0.8
    max_interactions_per_type: Optional[Dict[str, int]] = None

    # PPI 相关参数
    protein1_color: str = 'cyan'
    protein2_color: str = 'magenta'
    clear_old: bool = True


# ============================================================================
# Core Interaction Colors (Schrödinger-inspired)
# ============================================================================

INTERACTION_COLORS_HEX = {
    'hbond':       '#2979FF',  # Blue - Hydrogen bonds
    'salt':        '#D32F2F',  # Red - Salt bridges
    'pipi':        '#7B1FA2',  # Purple - π-π stacking
    'pication':    '#F57C00',  # Orange - π-cation
    'hydrophobic': '#388E3C',  # Green - Hydrophobic interactions
    'halogen':     '#FFB300',  # Amber/Gold - Halogen bonds
    'metal':       '#795548',  # Brown - Metal coordination
    'water':       '#00ACC1',  # Cyan - Water bridges
    'other':       '#9E9E9E',  # Gray - Other/unknown
}

# PyMOL-compatible RGB (0-1 range)
INTERACTION_COLORS_PYMOL = {
    'hbond':       [0.161, 0.475, 1.000],  # #2979FF
    'salt':        [0.827, 0.184, 0.184],  # #D32F2F
    'pipi':        [0.482, 0.122, 0.635],  # #7B1FA2
    'pication':    [0.961, 0.486, 0.000],  # #F57C00
    'hydrophobic': [0.220, 0.557, 0.235],  # #388E3C
    'halogen':     [1.000, 0.702, 0.000],  # #FFB300
    'metal':       [0.475, 0.333, 0.282],  # #795548
    'water':       [0.000, 0.675, 0.757],  # #00ACC1
    'other':       [0.620, 0.620, 0.620],  # #9E9E9E
}

# Qt/GUI-compatible RGB (0-255 range)
INTERACTION_COLORS_QT = {
    'hbond':       (41, 121, 255),   # #2979FF
    'salt':        (211, 47, 47),    # #D32F2F
    'pipi':        (123, 31, 162),   # #7B1FA2
    'pication':    (245, 124, 0),    # #F57C00
    'hydrophobic': (56, 142, 60),    # #388E3C
    'halogen':     (255, 179, 0),    # #FFB300
    'metal':       (121, 85, 72),   # #795548
    'water':       (0, 172, 193),    # #00ACC1
    'other':       (158, 158, 158),  # #9E9E9E
}

# ============================================================================
# Residue Type Colors (Discovery Studio style)
# ============================================================================

RESIDUE_COLORS_HEX = {
    'hydrophobic': {'face': '#C5E1A5', 'edge': '#388E3C', 'text': '#2E7D32'},  # Light green
    'nonpolar':    {'face': '#FFE0B2', 'edge': '#FFB300', 'text': '#E65100'},  # Light orange
    'polar':       {'face': '#BBDEFB', 'edge': '#2979FF', 'text': '#1565C0'},  # Light blue
    'negative':    {'face': '#FFCDD2', 'edge': '#F44336', 'text': '#C62828'},  # Light red (ASP/GLU)
    'positive':    {'face': '#E1BEE7', 'edge': '#7B1FA2', 'text': '#6A1B9A'},  # Light purple (LYS/ARG/HIS)
}

# ============================================================================
# PyMOL Custom Color Names
# ============================================================================

PYMOL_COLOR_NAMES = {
    'hbond':       'glue_hbond',
    'salt':        'glue_salt',
    'pipi':        'glue_pipi',
    'pication':    'glue_pication',
    'hydrophobic': 'glue_hydrophobic',
    'halogen':     'glue_halogen',
    'metal':       'glue_metal',
    'water':       'glue_water',
    'other':       'glue_other',
}

# ============================================================================
# Interaction Line Styles (for 2D and 3D visualization)
# ============================================================================

INTERACTION_LINE_STYLES = {
    'hbond': {
        'color': INTERACTION_COLORS_HEX['hbond'],
        'linewidth': 2.0,
        'linestyle': '--',
        'label': 'Hydrogen Bond',
        'dash_width': 2.0,
        'dash_gap': 0.3,
    },
    'salt': {
        'color': INTERACTION_COLORS_HEX['salt'],
        'linewidth': 2.5,
        'linestyle': '-',
        'label': 'Salt Bridge',
        'dash_width': 2.5,
        'dash_gap': 0.25,
    },
    'pipi': {
        'color': INTERACTION_COLORS_HEX['pipi'],
        'linewidth': 2.0,
        'linestyle': '--',
        'label': 'Pi-Pi Stacking',
        'dash_width': 2.0,
        'dash_gap': 0.3,
    },
    'pication': {
        'color': INTERACTION_COLORS_HEX['pication'],
        'linewidth': 2.0,
        'linestyle': '--',
        'label': 'Pi-Cation',
        'dash_width': 2.0,
        'dash_gap': 0.3,
    },
    'hydrophobic': {
        'color': INTERACTION_COLORS_HEX['hydrophobic'],
        'linewidth': 1.5,
        'linestyle': ':',
        'label': 'Hydrophobic',
        'dash_width': 1.5,
        'dash_gap': 0.35,
    },
    'halogen': {
        'color': INTERACTION_COLORS_HEX['halogen'],
        'linewidth': 2.0,
        'linestyle': '--',
        'label': 'Halogen Bond',
        'dash_width': 2.0,
        'dash_gap': 0.3,
    },
    'metal': {
        'color': INTERACTION_COLORS_HEX['metal'],
        'linewidth': 2.5,
        'linestyle': '-',
        'label': 'Metal Coordination',
        'dash_width': 2.5,
        'dash_gap': 0.25,
    },
    'water': {
        'color': INTERACTION_COLORS_HEX['water'],
        'linewidth': 2.0,
        'linestyle': '--',
        'label': 'Water Bridge',
        'dash_width': 2.0,
        'dash_gap': 0.3,
    },
    'other': {
        'color': INTERACTION_COLORS_HEX['other'],
        'linewidth': 1.5,
        'linestyle': ':',
        'label': 'Other',
        'dash_width': 1.5,
        'dash_gap': 0.35,
    },
}

# ============================================================================
# Helper Functions
# ============================================================================

def hex_to_rgb_0_1(hex_color):
    """
    Convert hex color to RGB tuple (0-1 range) for PyMOL
    
    Args:
        hex_color: Hex color string (e.g., '#2979FF')
    
    Returns:
        tuple: (r, g, b) with values in 0-1 range
    """
    hex_color = hex_color.lstrip('#')
    r, g, b = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    return (r/255.0, g/255.0, b/255.0)

def hex_to_rgb_0_255(hex_color):
    """
    Convert hex color to RGB tuple (0-255 range) for Qt
    
    Args:
        hex_color: Hex color string (e.g., '#2979FF')
    
    Returns:
        tuple: (r, g, b) with values in 0-255 range
    """
    hex_color = hex_color.lstrip('#')
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))

def register_pymol_colors(cmd):
    """
    Register custom colors in PyMOL
    
    Args:
        cmd: PyMOL cmd module
    """
    for key, rgb in INTERACTION_COLORS_PYMOL.items():
        color_name = PYMOL_COLOR_NAMES.get(key, f'glue_{key}')
        try:
            cmd.set_color(color_name, rgb)
        except Exception as e:
            print(f"[color_scheme] Warning: Failed to register color {color_name}: {e}")

def get_interaction_color(interaction_type, format='hex'):
    """
    Get color for an interaction type in the specified format
    
    Args:
        interaction_type: Interaction type key (e.g., 'hbond', 'salt')
        format: Color format - 'hex', 'pymol', 'qt', or 'name'
    
    Returns:
        Color in the requested format, or gray if type not found
    """
    itype = interaction_type.lower()
    
    if format == 'hex':
        return INTERACTION_COLORS_HEX.get(itype, INTERACTION_COLORS_HEX['other'])
    elif format == 'pymol':
        return INTERACTION_COLORS_PYMOL.get(itype, INTERACTION_COLORS_PYMOL['other'])
    elif format == 'qt':
        return INTERACTION_COLORS_QT.get(itype, INTERACTION_COLORS_QT['other'])
    elif format == 'name':
        return PYMOL_COLOR_NAMES.get(itype, 'gray')
    else:
        raise ValueError(f"Unknown format: {format}")

def get_interaction_style(interaction_type):
    """
    Get complete line style for an interaction type

    Args:
        interaction_type: Interaction type key

    Returns:
        dict: Style dictionary with color, linewidth, linestyle, etc.
    """
    return INTERACTION_LINE_STYLES.get(
        interaction_type.lower(),
        INTERACTION_LINE_STYLES['other']
    )


def apply_publication_pymol_style(cmd, background='white', label_size=16, label_font_id=5):
    """统一应用出版级 PyMOL 渲染参数（3D 图通用）。"""
    # 背景与光照：白底、较柔和光照，减少阴影干扰，适合论文插图
    cmd.bg_color(background)
    cmd.set('ray_opaque_background', 1)
    cmd.set('ray_shadow', 0)
    # 方案B：保持纯白底，同时降低整体光照强度，避免画面过亮
    cmd.set('ambient', 0.22)
    cmd.set('direct', 0.48)
    cmd.set('spec_reflect', 0.20)
    cmd.set('spec_power', 80)

    # 抗锯齿与景深：保证边缘清晰、避免景深导致细节发灰
    cmd.set('antialias', 2)
    cmd.set('depth_cue', 0)
    cmd.set('ray_trace_mode', 1)
    cmd.set('orthoscopic', 1)

    # 标签字体统一：Times-like（font_id=5）+ 黑色，便于出版阅读
    cmd.set('label_size', label_size)
    cmd.set('label_font_id', label_font_id)
    cmd.set('label_color', 'black')


def apply_interaction_dash_style(cmd, obj_name, color_name=None,
                                 dash_width=2.0, dash_gap=0.3,
                                 dash_length=0.25, dash_radius=0.08,
                                 hide_labels=True):
    """统一应用相互作用虚线样式，避免不同模块风格不一致。"""
    cmd.show('dashes', obj_name)
    cmd.set('dash_width', dash_width, obj_name)
    cmd.set('dash_gap', dash_gap, obj_name)
    cmd.set('dash_length', dash_length, obj_name)
    cmd.set('dash_radius', dash_radius, obj_name)

    if color_name:
        cmd.color(color_name, obj_name)
        cmd.set('dash_color', color_name, obj_name)

    if hide_labels:
        cmd.hide('labels', obj_name)