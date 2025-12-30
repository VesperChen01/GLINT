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
    from gluetk.ec_visualization import visualize_ec_smooth_surface
    
    # After running calculate_ligand_ec:
    result = calculate_ligand_ec('complex', 'LIG', output_dir='./ec_output')
    visualize_ec_smooth_surface('complex', 'LIG', ec_result=result)

Author: GlueTK Team
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
                                 ray_trace: bool = False) -> bool:
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
    
    ligand_sel = f"{obj_name} and resn {ligand_resname}"
    surface_obj = f"ec_surface_{ligand_resname}"
    
    # Step 1: Create a copy of the ligand for surface generation
    cmd.create(surface_obj, ligand_sel)
    
    # Step 2: Map EC values to atom B-factors
    if SCIPY_AVAILABLE:
        _map_ec_to_bfactors_kdtree(surface_obj, surface_points, ec_values)
    else:
        _map_ec_to_bfactors_simple(surface_obj, surface_points, ec_values)
    
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
        # Red-White-Green gradient
        cmd.spectrum('b', 'red_white_green', surface_obj, minimum=b_min, maximum=b_max)
    elif color_scheme == 'bwr':
        # Blue-White-Red gradient
        cmd.spectrum('b', 'blue_white_red', surface_obj, minimum=b_min, maximum=b_max)
    else:
        # Default to red-white-green
        cmd.spectrum('b', 'red_white_green', surface_obj, minimum=b_min, maximum=b_max)
    
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
    print("[visualize_ec_smooth_surface] Color scheme:")
    print("  🔴 Red = Electrostatic clash (EC < 0)")
    print("  ⚪ White = Neutral (EC ≈ 0)")
    print("  🟢 Green = Electrostatic complementarity (EC > 0)")
    print("[visualize_ec_smooth_surface] 💡 Use 'ray' command for publication quality")
    
    return True


def _map_ec_to_bfactors_kdtree(surface_obj: str, surface_points: np.ndarray, 
                                ec_values: np.ndarray):
    """Map EC values to atom B-factors using KDTree for efficiency."""
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
                                ec_values: np.ndarray):
    """Map EC values to atom B-factors using simple nearest neighbor."""
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
                              ec_range: tuple = (-1.0, 1.0)) -> bool:
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
        
    Returns:
        True if successful
    """
    if not PYMOL_AVAILABLE:
        print("[visualize_ec_cgo_surface] ❌ PyMOL required")
        return False
    
    print(f"[visualize_ec_cgo_surface] Creating CGO surface with {len(surface_points)} points...")
    
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
    Create a color legend for EC visualization.
    
    Args:
        position: Legend position ('bottom_right', 'bottom_left', etc.)
        ec_range: (min, max) EC values
        color_scheme: Color scheme used
        title: Legend title
        
    Returns:
        True if successful
    """
    if not PYMOL_AVAILABLE:
        return False
    
    # Create a simple text-based legend using PyMOL's label feature
    # This is a simplified version - for publication, use external tools
    
    ec_min, ec_max = ec_range
    
    print(f"\n{'='*40}")
    print(f"EC Color Legend: {title}")
    print(f"{'='*40}")
    
    if color_scheme == 'rwg':
        print(f"🔴 Red:   EC = {ec_min:.1f} (Clash)")
        print(f"⚪ White: EC = 0.0 (Neutral)")
        print(f"🟢 Green: EC = {ec_max:.1f} (Complementary)")
    else:
        print(f"🔵 Blue:  EC = {ec_min:.1f} (Clash)")
        print(f"⚪ White: EC = 0.0 (Neutral)")
        print(f"🔴 Red:   EC = {ec_max:.1f} (Complementary)")
    
    print(f"{'='*40}")
    
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
    """Create CGO for an interface with EC-based coloring."""
    cgo_obj = []
    
    # Base color RGB values
    base_colors = {
        'cyan': (0.0, 1.0, 1.0),
        'magenta': (1.0, 0.0, 1.0),
        'yellow': (1.0, 1.0, 0.0),
        'white': (1.0, 1.0, 1.0)
    }
    
    base_rgb = base_colors.get(base_color, (1.0, 1.0, 1.0))
    
    for point, ec_val in zip(points, ec_values):
        # Skip masked points (EC = 0)
        if ec_val == 0:
            continue
        
        # Modulate color by EC value
        # Positive EC: more saturated base color
        # Negative EC: shift toward red
        if ec_val > 0:
            # Positive: base color with intensity based on EC
            intensity = min(1.0, 0.5 + ec_val * 0.5)
            r = base_rgb[0] * intensity
            g = base_rgb[1] * intensity
            b = base_rgb[2] * intensity
        else:
            # Negative: shift toward red
            t = min(1.0, abs(ec_val))
            r = 1.0
            g = base_rgb[1] * (1 - t)
            b = base_rgb[2] * (1 - t)
        
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


# Module info
if __name__ == '__main__':
    print("="*60)
    print("GlueTK EC Visualization Module")
    print("="*60)
    print("\nThis module provides publication-quality EC surface visualization.")
    print("\nUsage in PyMOL:")
    print("  from gluetk.ec_visualization import visualize_ec_smooth_surface")
    print("  ")
    print("  # After running calculate_ligand_ec:")
    print("  result = calculate_ligand_ec('complex', 'LIG')")
    print("  visualize_ec_smooth_surface('complex', 'LIG', ec_result=result)")
    print("\nDependencies:")
    print(f"  NumPy: {'✅' if NUMPY_AVAILABLE else '❌'}")
    print(f"  SciPy: {'✅' if SCIPY_AVAILABLE else '❌'}")
    print(f"  PyMOL: {'✅' if PYMOL_AVAILABLE else '❌'}")