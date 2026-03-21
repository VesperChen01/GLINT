# -*- coding: utf-8 -*-
"""
pymol_styles.py
Professional PyMOL visualization styles for GLINT

This module provides publication-quality PyMOL rendering styles inspired by
professional molecular visualization standards. It includes:
- Protein representation styles (cartoon, surface, sticks)
- Ligand representation styles
- Color schemes for different molecular components
- Ray-tracing and rendering settings
"""

from typing import Optional, List, Tuple
try:
    from pymol import cmd
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False
    cmd = None


# ============================================================================
# Color Definitions (Professional Palette)
# ============================================================================

# Protein colors (soft, publication-friendly)
PROTEIN_COLORS = {
    'chain_a': [0.53, 0.81, 0.92],      # Light sky blue
    'chain_b': [1.00, 0.71, 0.76],      # Light pink
    'chain_c': [0.68, 0.85, 0.68],      # Light green
    'chain_d': [1.00, 0.85, 0.56],      # Light yellow
    'helix': [0.94, 0.20, 0.20],        # Red for helices
    'sheet': [1.00, 0.78, 0.05],        # Yellow for sheets
    'loop': [0.13, 0.55, 0.13],         # Forest green for loops
}

# Ligand colors (vibrant but professional)
LIGAND_COLORS = {
    'carbon': [0.20, 0.60, 0.20],       # Forest green
    'nitrogen': [0.20, 0.40, 0.80],     # Blue
    'oxygen': [0.80, 0.20, 0.20],       # Red
    'sulfur': [0.90, 0.78, 0.00],       # Yellow
    'phosphorus': [1.00, 0.50, 0.00],   # Orange
    'halogen': [0.12, 0.94, 0.12],      # Bright green
    'hydrogen': [0.90, 0.90, 0.90],     # Light gray
}

# Additional colors for Science-style visualization
SCIENCE_COLORS = {
    'deeporange': [1.00, 0.27, 0.00],   # Deep orange for carbon
    'wheat': [0.96, 0.87, 0.70],        # Wheat for hydrogen
    'tan': [0.82, 0.71, 0.55],          # Tan for nitrogen
    'limon': [0.75, 1.00, 0.00],        # Lime green for halogens
}

# Surface colors
SURFACE_COLORS = {
    'hydrophobic': [0.95, 0.78, 0.00],  # Gold
    'hydrophilic': [0.53, 0.81, 0.92],  # Light blue
    'positive': [0.20, 0.40, 0.80],     # Blue
    'negative': [0.80, 0.20, 0.20],     # Red
    'neutral': [0.90, 0.90, 0.90],      # Light gray
}


# ============================================================================
# Core Style Application Functions
# ============================================================================

def apply_professional_style(background='white', ray_trace_mode=1):
    """
    Apply professional publication-quality rendering settings.
    
    This sets up optimal lighting, shadows, and anti-aliasing for
    high-quality molecular visualization.
    
    Args:
        background: Background color ('white', 'black', or custom)
        ray_trace_mode: Ray tracing mode (0=no shadows, 1=normal, 3=black outline)
    """
    if not PYMOL_AVAILABLE:
        return
    
    # Background
    cmd.bg_color(background)
    cmd.set('ray_opaque_background', 1)
    
    # Lighting and shadows
    cmd.set('ambient', 0.2)
    cmd.set('direct', 0.6)
    cmd.set('reflect', 0.4)
    cmd.set('specular', 0.5)
    cmd.set('shininess', 10)
    cmd.set('spec_power', 200)
    cmd.set('spec_reflect', 0.5)
    
    # Ray tracing
    cmd.set('ray_trace_mode', ray_trace_mode)
    cmd.set('ray_shadow', 1 if ray_trace_mode > 0 else 0)
    cmd.set('ray_trace_fog', 0)
    
    # Anti-aliasing
    cmd.set('antialias', 2)
    cmd.set('hash_max', 300)
    
    # Depth and perspective
    cmd.set('depth_cue', 0)
    try:
        cmd.set('ray_depth_cue', 0)
    except:
        # ray_depth_cue not available in all PyMOL versions
        pass
    cmd.set('orthoscopic', 1)
    
    # Surface quality
    cmd.set('surface_quality', 1)
    cmd.set('surface_type', 0)  # 0=solid, 1=mesh, 2=dots
    cmd.set('transparency_mode', 1)
    
    print("[pymol_styles] ✓ Professional rendering style applied")


def setup_protein_cartoon(selection, color='auto', transparency=0.0, 
                         fancy_helices=True, smooth_loops=True, color_scheme='science'):
    """
    Set up professional cartoon representation for proteins.
    
    Args:
        selection: PyMOL selection string
        color: Color scheme ('auto', 'ss' for secondary structure, or color name)
        transparency: Cartoon transparency (0.0-1.0)
        fancy_helices: Use fancy helix representation
        smooth_loops: Smooth loop regions
        color_scheme: 'science' (deep blue/green) or 'classic' (default PyMOL colors)
    """
    if not PYMOL_AVAILABLE:
        return
    
    # Hide everything first
    cmd.hide('everything', selection)
    
    # Show cartoon
    cmd.show('cartoon', selection)
    
    # Cartoon settings
    cmd.set('cartoon_fancy_helices', 1 if fancy_helices else 0, selection)
    cmd.set('cartoon_smooth_loops', 1 if smooth_loops else 0, selection)
    cmd.set('cartoon_loop_radius', 0.2, selection)
    cmd.set('cartoon_tube_radius', 0.3, selection)
    cmd.set('cartoon_oval_length', 1.2, selection)
    cmd.set('cartoon_oval_width', 0.25, selection)
    cmd.set('cartoon_rect_length', 1.4, selection)
    cmd.set('cartoon_rect_width', 0.4, selection)
    
    # Transparency
    if transparency > 0:
        cmd.set('cartoon_transparency', transparency, selection)
    
    # Coloring
    if color == 'auto':
        if color_scheme == 'science':
            # Science cover style: deep blue + deep lake-green
            cmd.set_color('FOG_BLUE_DEEP', [0.45, 0.62, 0.78])
            cmd.set_color('LAKE_GREEN_DEEP', [0.60, 0.76, 0.68])
            
            try:
                chains = cmd.get_chains(selection)
                if len(chains) >= 1:
                    cmd.color('FOG_BLUE_DEEP', f'{selection} and chain {chains[0]}')
                if len(chains) >= 2:
                    cmd.color('LAKE_GREEN_DEEP', f'{selection} and chain {chains[1]}')
                # Additional chains cycle through these colors
                for i in range(2, len(chains)):
                    color_name = 'FOG_BLUE_DEEP' if i % 2 == 0 else 'LAKE_GREEN_DEEP'
                    cmd.color(color_name, f'{selection} and chain {chains[i]}')
            except:
                cmd.color('FOG_BLUE_DEEP', selection)
        else:
            cmd.util.cbc(selection)  # Color by chain (classic)
    elif color == 'ss':
        cmd.color('red', f'{selection} and ss h')      # Helices
        cmd.color('yellow', f'{selection} and ss s')   # Sheets
        cmd.color('green', f'{selection} and ss l+')   # Loops
    else:
        cmd.color(color, selection)
    
    print(f"[pymol_styles] ✓ Cartoon representation applied to {selection}")


def setup_protein_surface(selection, color='auto', transparency=0.5,
                         surface_type='solid', quality=1, color_scheme='science'):
    """
    Set up professional surface representation for proteins.
    
    Args:
        selection: PyMOL selection string
        color: Color scheme ('auto', 'hydrophobicity', 'electrostatics', or color name)
        transparency: Surface transparency (0.0-1.0)
        surface_type: 'solid', 'mesh', or 'dots'
        quality: Surface quality (0=low, 1=medium, 2=high)
        color_scheme: 'science' (deep blue/green) or 'classic' (default PyMOL colors)
    """
    if not PYMOL_AVAILABLE:
        return
    
    # Show surface
    cmd.show('surface', selection)
    
    # Surface settings
    cmd.set('surface_quality', quality, selection)
    
    surface_type_map = {'solid': 0, 'mesh': 1, 'dots': 2}
    cmd.set('surface_type', surface_type_map.get(surface_type, 0), selection)
    
    # Transparency
    cmd.set('transparency', transparency, selection)
    
    # Coloring
    if color == 'auto':
        if color_scheme == 'science':
            # Science cover style: deep blue + deep lake-green
            cmd.set_color('FOG_BLUE_DEEP', [0.45, 0.62, 0.78])
            cmd.set_color('LAKE_GREEN_DEEP', [0.60, 0.76, 0.68])
            
            try:
                chains = cmd.get_chains(selection)
                if len(chains) >= 1:
                    cmd.color('FOG_BLUE_DEEP', f'{selection} and chain {chains[0]}')
                if len(chains) >= 2:
                    cmd.color('LAKE_GREEN_DEEP', f'{selection} and chain {chains[1]}')
                # Additional chains cycle through these colors
                for i in range(2, len(chains)):
                    color_name = 'FOG_BLUE_DEEP' if i % 2 == 0 else 'LAKE_GREEN_DEEP'
                    cmd.color(color_name, f'{selection} and chain {chains[i]}')
            except:
                cmd.color('FOG_BLUE_DEEP', selection)
        else:
            cmd.util.cbc(selection)
    elif color == 'hydrophobicity':
        # Color by hydrophobicity (requires pre-calculated B-factors)
        cmd.spectrum('b', 'blue_white_yellow', selection)
    elif color == 'electrostatics':
        # Color by electrostatics (requires APBS or similar)
        cmd.spectrum('b', 'red_white_blue', selection)
    else:
        cmd.color(color, selection)
    
    print(f"[pymol_styles] ✓ Surface representation applied to {selection}")


def setup_ligand_sticks(selection, carbon_color='deeporange', show_polar_h=True,
                       stick_radius=0.25, ball_and_stick=False, color_scheme='orange'):
    """
    Set up professional stick representation for ligands.
    
    Args:
        selection: PyMOL selection string
        carbon_color: Color for carbon atoms ('green', 'cyan', 'yellow', 'deeporange', etc.)
        show_polar_h: Show only polar hydrogens (NH, OH, SH)
        stick_radius: Stick radius (default 0.25)
        ball_and_stick: Use ball-and-stick representation
        color_scheme: 'orange' (Science cover style) or 'classic' (traditional green)
    """
    if not PYMOL_AVAILABLE:
        return
    
    # Hide everything first
    cmd.hide('everything', selection)
    
    # Show sticks
    cmd.show('sticks', selection)
    cmd.set('stick_radius', stick_radius, selection)
    
    # Color by element - Science cover style (orange family)
    if color_scheme == 'orange':
        # Define custom deep orange color
        cmd.set_color('deep_orange', [0.90, 0.45, 0.20])
        cmd.color('deep_orange', f'{selection} and elem C')
        cmd.color('wheat', f'{selection} and elem H')
        cmd.color('tan', f'{selection} and elem N')
        cmd.color('firebrick', f'{selection} and elem O')
        cmd.color('yellow', f'{selection} and elem S')
        cmd.color('orange', f'{selection} and elem P')
        cmd.color('limon', f'{selection} and elem F,Cl,Br,I')
    else:
        # Classic color scheme
        cmd.color(carbon_color, f'{selection} and elem C')
        cmd.color('blue', f'{selection} and elem N')
        cmd.color('red', f'{selection} and elem O')
        cmd.color('yellow', f'{selection} and elem S')
        cmd.color('orange', f'{selection} and elem P')
        cmd.color('limon', f'{selection} and elem F,Cl,Br,I')
    
    # Hydrogen display
    if show_polar_h:
        cmd.hide('sticks', f'{selection} and elem H')
        cmd.show('sticks', f'{selection} and elem H and (neighbor elem N+O+S)')
    
    # Ball-and-stick mode
    if ball_and_stick:
        cmd.show('spheres', selection)
        cmd.set('sphere_scale', 0.25, selection)
    
    print(f"[pymol_styles] ✓ Stick representation applied to {selection}")


def setup_binding_site(protein_sel, ligand_sel, distance=5.0,
                      show_protein_surface=True, show_residue_sticks=True,
                      color_scheme='science'):
    """
    Set up professional binding site visualization.
    
    Args:
        protein_sel: Protein selection string
        ligand_sel: Ligand selection string
        distance: Distance cutoff for binding site residues (Å)
        show_protein_surface: Show protein surface around binding site
        show_residue_sticks: Show binding site residues as sticks
        color_scheme: 'science' (deep blue/green + orange) or 'classic' (green)
    """
    if not PYMOL_AVAILABLE:
        return
    
    # Create binding site selection
    binding_site = f'binding_site_{ligand_sel}'
    cmd.select(binding_site, f'{protein_sel} within {distance} of {ligand_sel}')
    
    # Protein cartoon (full protein)
    setup_protein_cartoon(protein_sel, color='auto', transparency=0.3, color_scheme=color_scheme)
    
    # Binding site surface
    if show_protein_surface:
        setup_protein_surface(binding_site, color='auto', transparency=0.5, color_scheme=color_scheme)
    
    # Binding site residues as sticks
    if show_residue_sticks:
        cmd.show('sticks', f'{binding_site} and (sidechain or name CA)')
        cmd.set('stick_radius', 0.15, binding_site)
        cmd.util.cbag(binding_site)  # Color by atom (with green carbons)
    
    # Ligand as sticks
    setup_ligand_sticks(ligand_sel, ball_and_stick=False, color_scheme=color_scheme)
    
    # Zoom to binding site
    cmd.zoom(ligand_sel, buffer=8)
    
    print(f"[pymol_styles] ✓ Binding site visualization complete")


# ============================================================================
# Specialized Visualization Functions
# ============================================================================

def apply_publication_figure_style(obj_name, ligand_resname=None,
                                   protein_color='lightblue', 
                                   ligand_carbon_color='deeporange',
                                   show_surface=False,
                                   background='white',
                                   color_scheme='science'):
    """
    Apply complete publication-quality figure style.
    
    This is a high-level function that sets up everything for a
    publication-ready molecular visualization.
    
    Args:
        obj_name: PyMOL object name
        ligand_resname: Ligand residue name (if present)
        protein_color: Color for protein cartoon
        ligand_carbon_color: Color for ligand carbons
        show_surface: Show protein surface
        background: Background color
        color_scheme: 'science' (deep blue/green + orange) or 'classic' (green)
    """
    if not PYMOL_AVAILABLE:
        return
    
    # Apply professional rendering settings
    apply_professional_style(background=background, ray_trace_mode=1)
    
    # Protein selection
    protein_sel = f'{obj_name} and polymer.protein'
    
    # Set up protein
    if color_scheme == 'science':
        setup_protein_cartoon(protein_sel, color='auto', transparency=0.0, color_scheme='science')
    else:
        setup_protein_cartoon(protein_sel, color=protein_color, transparency=0.0, color_scheme='classic')
    
    if show_surface:
        if color_scheme == 'science':
            setup_protein_surface(protein_sel, color='auto', transparency=0.36, color_scheme='science')
        else:
            setup_protein_surface(protein_sel, color=protein_color, transparency=0.5, color_scheme='classic')
    
    # Set up ligand if present
    if ligand_resname:
        ligand_sel = f'{obj_name} and resn {ligand_resname}'
        setup_ligand_sticks(ligand_sel, color_scheme=color_scheme)
        
        # Highlight binding site
        binding_site = f'binding_site_{ligand_resname}'
        cmd.select(binding_site, f'{protein_sel} within 5.0 of {ligand_sel}')
        cmd.show('sticks', f'{binding_site} and (sidechain or name CA)')
        cmd.set('stick_radius', 0.15, binding_site)
        cmd.util.cbag(binding_site)
        
        # Zoom to ligand
        cmd.zoom(ligand_sel, buffer=8)
    else:
        cmd.zoom(obj_name)
    
    print(f"[pymol_styles] ✓ Publication figure style applied to {obj_name}")


def register_custom_colors():
    """Register all custom colors in PyMOL."""
    if not PYMOL_AVAILABLE:
        return
    
    # Register protein colors
    for name, rgb in PROTEIN_COLORS.items():
        try:
            cmd.set_color(f'glint_{name}', rgb)
        except:
            pass
    
    # Register ligand colors
    for name, rgb in LIGAND_COLORS.items():
        try:
            cmd.set_color(f'glint_lig_{name}', rgb)
        except:
            pass
    
    # Register Science-style colors
    for name, rgb in SCIENCE_COLORS.items():
        try:
            cmd.set_color(name, rgb)
        except:
            pass
    
    # Register surface colors
    for name, rgb in SURFACE_COLORS.items():
        try:
            cmd.set_color(f'glint_surf_{name}', rgb)
        except:
            pass
    
    print("[pymol_styles] ✓ Custom colors registered")


# ============================================================================
# Convenience Functions
# ============================================================================

def quick_protein_view(obj_name, style='cartoon', color='auto', color_scheme='science'):
    """Quick protein visualization with sensible defaults."""
    if not PYMOL_AVAILABLE:
        return
    
    apply_professional_style()
    
    if style == 'cartoon':
        setup_protein_cartoon(obj_name, color=color, color_scheme=color_scheme)
    elif style == 'surface':
        setup_protein_surface(obj_name, color=color, color_scheme=color_scheme)
    elif style == 'both':
        setup_protein_cartoon(obj_name, color=color, transparency=0.3, color_scheme=color_scheme)
        setup_protein_surface(obj_name, color=color, transparency=0.5, color_scheme=color_scheme)
    
    cmd.zoom(obj_name)


def quick_ligand_view(obj_name, ligand_resname, show_binding_site=True, color_scheme='science'):
    """Quick ligand visualization with binding site."""
    if not PYMOL_AVAILABLE:
        return
    
    apply_professional_style()
    
    protein_sel = f'{obj_name} and polymer.protein'
    ligand_sel = f'{obj_name} and resn {ligand_resname}'
    
    if show_binding_site:
        setup_binding_site(protein_sel, ligand_sel, color_scheme=color_scheme)
    else:
        setup_protein_cartoon(protein_sel, color_scheme=color_scheme)
        setup_ligand_sticks(ligand_sel, color_scheme=color_scheme)
        cmd.zoom(ligand_sel, buffer=8)


# ============================================================================
# PyMOL Command Extensions
# ============================================================================

if PYMOL_AVAILABLE:
    # Register commands
    cmd.extend('glint_style', apply_professional_style)
    cmd.extend('glint_protein', setup_protein_cartoon)
    cmd.extend('glint_surface', setup_protein_surface)
    cmd.extend('glint_ligand', setup_ligand_sticks)
    cmd.extend('glint_binding_site', setup_binding_site)
    cmd.extend('glint_publication', apply_publication_figure_style)
    cmd.extend('glint_quick_protein', quick_protein_view)
    cmd.extend('glint_quick_ligand', quick_ligand_view)
    
    # Register colors on import
    register_custom_colors()
