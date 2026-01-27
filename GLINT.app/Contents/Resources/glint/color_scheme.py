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

# ============================================================================
# Core Interaction Colors (Schrödinger-inspired)
# ============================================================================

INTERACTION_COLORS_HEX = {
    'hbond':       '#2196F3',  # Blue - Hydrogen bonds
    'salt':        '#FF5722',  # Orange-red - Salt bridges
    'pipi':        '#9C27B0',  # Purple - π-π stacking
    'pication':    '#E91E63',  # Pink - π-cation
    'hydrophobic': '#4CAF50',  # Green - Hydrophobic interactions
    'halogen':     '#FF9800',  # Orange - Halogen bonds
    'metal':       '#673AB7',  # Deep purple - Metal coordination
    'water':       '#00BCD4',  # Cyan - Water bridges
    'other':       '#9E9E9E',  # Gray - Other/unknown
}

# PyMOL-compatible RGB (0-1 range)
INTERACTION_COLORS_PYMOL = {
    'hbond':       [0.129, 0.588, 0.953],  # #2196F3
    'salt':        [1.000, 0.341, 0.133],  # #FF5722
    'pipi':        [0.612, 0.153, 0.690],  # #9C27B0
    'pication':    [0.914, 0.118, 0.388],  # #E91E63
    'hydrophobic': [0.298, 0.686, 0.314],  # #4CAF50
    'halogen':     [1.000, 0.596, 0.000],  # #FF9800
    'metal':       [0.404, 0.227, 0.718],  # #673AB7
    'water':       [0.000, 0.737, 0.831],  # #00BCD4
    'other':       [0.620, 0.620, 0.620],  # #9E9E9E
}

# Qt/GUI-compatible RGB (0-255 range)
INTERACTION_COLORS_QT = {
    'hbond':       (33, 150, 243),   # #2196F3
    'salt':        (255, 87, 34),    # #FF5722
    'pipi':        (156, 39, 176),   # #9C27B0
    'pication':    (233, 30, 99),    # #E91E63
    'hydrophobic': (76, 175, 80),    # #4CAF50
    'halogen':     (255, 152, 0),    # #FF9800
    'metal':       (103, 58, 183),   # #673AB7
    'water':       (0, 188, 212),    # #00BCD4
    'other':       (158, 158, 158),  # #9E9E9E
}

# ============================================================================
# Residue Type Colors (Discovery Studio style)
# ============================================================================

RESIDUE_COLORS_HEX = {
    'hydrophobic': {'face': '#C5E1A5', 'edge': '#4CAF50', 'text': '#2E7D32'},  # Light green
    'nonpolar':    {'face': '#FFE0B2', 'edge': '#FF9800', 'text': '#E65100'},  # Light orange
    'polar':       {'face': '#BBDEFB', 'edge': '#2196F3', 'text': '#1565C0'},  # Light blue
    'negative':    {'face': '#FFCDD2', 'edge': '#F44336', 'text': '#C62828'},  # Light red (ASP/GLU)
    'positive':    {'face': '#E1BEE7', 'edge': '#9C27B0', 'text': '#6A1B9A'},  # Light purple (LYS/ARG/HIS)
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
        hex_color: Hex color string (e.g., '#2196F3')
    
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
        hex_color: Hex color string (e.g., '#2196F3')
    
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