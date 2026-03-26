# -*- coding: utf-8 -*-
"""
ligand_ec_calculator.py
Electrostatic Complementarity (EC) Calculator for Molecular Glue Ligands

This module implements the EC analysis workflow described in apbs_ec_open_source_workflow.md:
1. PDB2PQR integration for protein structure preparation
2. APBS integration for protein electrostatic potential calculation
3. Ligand surface sampling
4. Ligand electrostatic potential calculation (Gasteiger charges)
5. EC score calculation and EC map generation
6. Ternary complex (molecular glue) EC analysis

Core Concepts:
- PIP (Protein Interaction Potential) ≈ Protein PB electrostatic field
- EC (Electrostatic Complementarity): Compare protein potential with ligand potential
  at ligand surface points. Opposite signs = complementary; same signs = clash.

Advanced Features (v2.0):
- σ-hole virtual points for halogen bonds (Cl/Br/I)
- Lone pair virtual points for carbonyl O
- Bridging water filter for structural waters

Usage:
    # Basic EC analysis
    result = calculate_ligand_ec('complex', 'LIG', output_dir='./ec_output')
    
    # Ternary complex (molecular glue) analysis
    result = analyze_ternary_ec('complex', 'GLUE', ['A'], ['B'], output_dir='./ec_output')
    
    # Enhanced EC with σ-hole (for halogen-containing ligands)
    result = calculate_ligand_ec('complex', 'LIG', use_sigma_holes=True)

Author: GLINT Team
"""

from __future__ import print_function
import os
import sys
import math
import tempfile
import subprocess
import shutil
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Any, Union
import csv

# NumPy is required for this module
try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False
    print("[ligand_ec_calculator] ⚠️ NumPy not installed - EC calculations unavailable")

# RDKit for ligand handling
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("[ligand_ec_calculator] ⚠️ RDKit not installed - ligand EC calculations unavailable")

# SciPy for interpolation
try:
    from scipy.interpolate import RegularGridInterpolator
    from scipy.spatial import cKDTree
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("[ligand_ec_calculator] ⚠️ SciPy not installed - interpolation unavailable")

# PyMOL integration
try:
    from pymol import cmd
    PYMOL_AVAILABLE = True
except ImportError:
    PYMOL_AVAILABLE = False

# Advanced patches (σ-hole, lone pairs, bridging waters)
try:
    from glint.ec_advanced_patches import (
        SigmaHoleGenerator,
        LonePairGenerator,
        BridgingWaterFilter,
        EnhancedChargeCalculator,
    )
    ADVANCED_PATCHES_AVAILABLE = True
except ImportError:
    ADVANCED_PATCHES_AVAILABLE = False
    print("[ligand_ec_calculator] ⚠️ Advanced patches not available (σ-hole, lone pairs)")

# Enhanced EC visualization module
try:
    from glint.ec_visualization import (
        visualize_ec_smooth_surface,
        visualize_ec_cgo_surface,
        save_ec_visualization,
        visualize_ternary_ec_surfaces,
    )
    EC_VISUALIZATION_AVAILABLE = True
except ImportError:
    EC_VISUALIZATION_AVAILABLE = False

# ========== Constants ==========
# Van der Waals radii (Å) for surface sampling
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47,
    'P': 1.80, 'S': 1.80, 'Cl': 1.75, 'Br': 1.85, 'I': 1.98,
    'Si': 2.10, 'B': 1.92, 'Se': 1.90, 'As': 1.85,
    'default': 1.70
}

# Probe radius for solvent accessible surface
PROBE_RADIUS = 1.4  # Å (water molecule)

# EC calculation parameters
EC_PARAMS = {
    'epsilon': 1e-6,           # Small value to avoid division by zero
    'phi_clip_min': -10.0,     # Minimum potential clip value (kT/e)
    'phi_clip_max': 10.0,      # Maximum potential clip value (kT/e)
    'surface_density': 10.0,   # Points per Å² for surface sampling
    'grid_spacing': 0.5,       # Å for APBS grid
}

# APBS default parameters
APBS_PARAMS = {
    'pdie': 2.0,      # Protein dielectric constant
    'sdie': 78.54,    # Solvent dielectric constant
    'temp': 298.15,   # Temperature (K)
    'ion_conc': 0.15, # Ionic concentration (M)
    'ion_charge': 1,  # Ion charge
    'ion_radius': 2.0 # Ion radius (Å)
}


class DXGrid:
    """
    Class to read and interpolate APBS OpenDX format potential grids.
    
    The OpenDX format contains:
    - Grid dimensions (nx, ny, nz)
    - Origin coordinates
    - Grid spacing (delta)
    - Potential values at each grid point
    """
    
    def __init__(self, filename: str = None):
        self.origin = np.zeros(3)
        self.delta = np.zeros((3, 3))  # Can be non-orthogonal
        self.dims = np.zeros(3, dtype=int)
        self.data = None
        self.interpolator = None
        
        if filename:
            self.read(filename)
    
    def read(self, filename: str):
        """Read an OpenDX format file from APBS."""
        if not os.path.exists(filename):
            raise FileNotFoundError(f"DX file not found: {filename}")
        
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        data_values = []
        reading_data = False
        
        for line in lines:
            line = line.strip()
            
            if line.startswith('#') or not line:
                continue
            
            # Parse grid dimensions
            if 'object 1 class gridpositions counts' in line:
                parts = line.split()
                self.dims = np.array([int(parts[-3]), int(parts[-2]), int(parts[-1])])
            
            # Parse origin
            elif line.startswith('origin'):
                parts = line.split()
                self.origin = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
            
            # Parse delta (grid spacing)
            elif line.startswith('delta'):
                parts = line.split()
                delta_row = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
                # Find which row this is based on non-zero element
                for i in range(3):
                    if abs(delta_row[i]) > 1e-10:
                        self.delta[i] = delta_row
                        break
            
            # Start reading data after this line
            elif 'object 3 class array' in line:
                reading_data = True
                continue
            
            # Stop reading at attribute line
            elif line.startswith('attribute') or line.startswith('object'):
                reading_data = False
            
            # Read data values
            elif reading_data:
                values = line.split()
                for v in values:
                    try:
                        data_values.append(float(v))
                    except ValueError:
                        pass
        
        # Reshape data to 3D grid
        expected_size = self.dims[0] * self.dims[1] * self.dims[2]
        if len(data_values) != expected_size:
            raise ValueError(f"Data size mismatch: got {len(data_values)}, expected {expected_size}")
        
        # APBS uses Fortran ordering (column-major)
        self.data = np.array(data_values).reshape(self.dims, order='F')
        
        # Create interpolator
        self._create_interpolator()
        
        print(f"[DXGrid] Loaded grid: dims={self.dims}, origin={self.origin}")
        print(f"[DXGrid] Grid spacing: {np.diag(self.delta)}")
        print(f"[DXGrid] Potential range: [{self.data.min():.3f}, {self.data.max():.3f}] kT/e")
    
    def _create_interpolator(self):
        """Create a scipy interpolator for the grid."""
        if not SCIPY_AVAILABLE:
            return
        
        # Create coordinate arrays for each axis
        x = self.origin[0] + np.arange(self.dims[0]) * self.delta[0, 0]
        y = self.origin[1] + np.arange(self.dims[1]) * self.delta[1, 1]
        z = self.origin[2] + np.arange(self.dims[2]) * self.delta[2, 2]
        
        self.interpolator = RegularGridInterpolator(
            (x, y, z), self.data,
            method='linear',
            bounds_error=False,
            fill_value=0.0
        )
    
    def interpolate(self, points: np.ndarray) -> np.ndarray:
        """
        Interpolate potential values at given points.
        
        Args:
            points: Nx3 array of coordinates
            
        Returns:
            N array of interpolated potential values
        """
        if self.interpolator is None:
            return self._interpolate_manual(points)
        
        return self.interpolator(points)
    
    def _interpolate_manual(self, points: np.ndarray) -> np.ndarray:
        """Manual trilinear interpolation (fallback if scipy unavailable)."""
        results = np.zeros(len(points))
        
        for i, point in enumerate(points):
            # Convert to grid coordinates
            grid_coords = (point - self.origin) / np.diag(self.delta)
            
            # Get integer indices
            i0 = int(np.floor(grid_coords[0]))
            j0 = int(np.floor(grid_coords[1]))
            k0 = int(np.floor(grid_coords[2]))
            
            # Check bounds
            if (i0 < 0 or i0 >= self.dims[0] - 1 or
                j0 < 0 or j0 >= self.dims[1] - 1 or
                k0 < 0 or k0 >= self.dims[2] - 1):
                results[i] = 0.0
                continue
            
            # Fractional parts
            xd = grid_coords[0] - i0
            yd = grid_coords[1] - j0
            zd = grid_coords[2] - k0
            
            # Trilinear interpolation
            c000 = self.data[i0, j0, k0]
            c001 = self.data[i0, j0, k0 + 1]
            c010 = self.data[i0, j0 + 1, k0]
            c011 = self.data[i0, j0 + 1, k0 + 1]
            c100 = self.data[i0 + 1, j0, k0]
            c101 = self.data[i0 + 1, j0, k0 + 1]
            c110 = self.data[i0 + 1, j0 + 1, k0]
            c111 = self.data[i0 + 1, j0 + 1, k0 + 1]
            
            c00 = c000 * (1 - xd) + c100 * xd
            c01 = c001 * (1 - xd) + c101 * xd
            c10 = c010 * (1 - xd) + c110 * xd
            c11 = c011 * (1 - xd) + c111 * xd
            
            c0 = c00 * (1 - yd) + c10 * yd
            c1 = c01 * (1 - yd) + c11 * yd
            
            results[i] = c0 * (1 - zd) + c1 * zd
        
        return results



def _get_element_charge(element: str) -> float:
    """
    Get approximate partial charge for an element (fallback for simple calculations).
    
    Args:
        element: Element symbol
        
    Returns:
        Approximate partial charge
    """
    charges = {
        'H': 0.1, 'C': 0.0, 'N': -0.3, 'O': -0.3, 'S': -0.1,
        'F': -0.2, 'Cl': -0.2, 'Br': -0.1, 'I': -0.1,
        'P': 0.1, 'FE': 2.0, 'ZN': 2.0, 'MG': 2.0, 'CA': 2.0
    }
    return charges.get(element.upper(), 0.0)

class LigandSurfaceSampler:
    """
    Generate surface sampling points for a ligand molecule.
    
    Uses a simple sphere-based approach:
    1. Place points on sphere around each atom (vdW radius + probe)
    2. Remove points that are inside other atoms
    3. Optionally cluster/downsample for efficiency
    """
    
    def __init__(self, density: float = 10.0, probe_radius: float = 1.4):
        """
        Args:
            density: Target points per Å² on the surface
            probe_radius: Solvent probe radius (Å)
        """
        self.density = density
        self.probe_radius = probe_radius
    
    def sample_molecule(self, mol: 'Chem.Mol', conformer_id: int = 0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate surface points for a molecule.
        
        Args:
            mol: RDKit molecule with 3D coordinates
            conformer_id: Which conformer to use
            
        Returns:
            Tuple of (surface_points, surface_normals)
        """
        if not RDKIT_AVAILABLE:
            raise RuntimeError("RDKit required for surface sampling")
        
        conf = mol.GetConformer(conformer_id)
        
        # Get atom coordinates and radii
        atoms = []
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            radius = VDW_RADII.get(symbol, VDW_RADII['default'])
            atoms.append({
                'pos': np.array([pos.x, pos.y, pos.z]),
                'radius': radius + self.probe_radius,
                'symbol': symbol
            })
        
        # Generate surface points
        all_points = []
        all_normals = []
        
        for atom in atoms:
            points, normals = self._sample_sphere(atom['pos'], atom['radius'])
            
            # Filter points inside other atoms
            valid_mask = np.ones(len(points), dtype=bool)
            for other_atom in atoms:
                if np.allclose(atom['pos'], other_atom['pos']):
                    continue
                
                distances = np.linalg.norm(points - other_atom['pos'], axis=1)
                valid_mask &= (distances >= other_atom['radius'] - 0.1)  # Small tolerance
            
            all_points.append(points[valid_mask])
            all_normals.append(normals[valid_mask])
        
        surface_points = np.vstack(all_points)
        surface_normals = np.vstack(all_normals)
        
        print(f"[LigandSurfaceSampler] Generated {len(surface_points)} surface points")
        
        return surface_points, surface_normals
    
    def _sample_sphere(self, center: np.ndarray, radius: float) -> Tuple[np.ndarray, np.ndarray]:
        """Generate evenly distributed points on a sphere."""
        # Calculate number of points based on surface area and density
        area = 4 * np.pi * radius ** 2
        n_points = max(int(area * self.density), 20)
        
        # Use Fibonacci sphere for even distribution
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
    
    def sample_from_coords(self, coords: np.ndarray, elements: List[str]) -> Tuple[np.ndarray, np.ndarray]:
        """
        Generate surface points from raw coordinates and elements.
        
        Args:
            coords: Nx3 array of atom coordinates
            elements: List of element symbols
            
        Returns:
            Tuple of (surface_points, surface_normals)
        """
        atoms = []
        for i, (coord, elem) in enumerate(zip(coords, elements)):
            radius = VDW_RADII.get(elem.upper(), VDW_RADII['default'])
            atoms.append({
                'pos': np.array(coord),
                'radius': radius + self.probe_radius,
                'symbol': elem
            })
        
        all_points = []
        all_normals = []
        
        for atom in atoms:
            points, normals = self._sample_sphere(atom['pos'], atom['radius'])
            
            valid_mask = np.ones(len(points), dtype=bool)
            for other_atom in atoms:
                if np.allclose(atom['pos'], other_atom['pos']):
                    continue
                
                distances = np.linalg.norm(points - other_atom['pos'], axis=1)
                valid_mask &= (distances >= other_atom['radius'] - 0.1)
            
            all_points.append(points[valid_mask])
            all_normals.append(normals[valid_mask])
        
        surface_points = np.vstack(all_points)
        surface_normals = np.vstack(all_normals)
        
        return surface_points, surface_normals


class GasteigerChargeCalculator:
    """
    Calculate Gasteiger partial charges and electrostatic potential for ligands.
    
    Uses RDKit's Gasteiger charge implementation and computes Coulomb potential
    at surface points.
    
    Enhanced Features (v2.0):
    - Optional σ-hole virtual points for halogens (Cl/Br/I)
    - Optional lone pair virtual points for carbonyl O
    """
    
    def __init__(self, dielectric: float = 4.0,
                 use_sigma_holes: bool = False,
                 use_lone_pairs: bool = False):
        """
        Args:
            dielectric: Effective dielectric constant for Coulomb calculation
            use_sigma_holes: Add σ-hole virtual points for Cl/Br/I
            use_lone_pairs: Add lone pair virtual points for carbonyl O
        """
        self.dielectric = dielectric
        self.use_sigma_holes = use_sigma_holes
        self.use_lone_pairs = use_lone_pairs
        
        # Initialize enhanced calculator if patches available
        if (use_sigma_holes or use_lone_pairs) and ADVANCED_PATCHES_AVAILABLE:
            self.enhanced_calc = EnhancedChargeCalculator(
                use_sigma_holes=use_sigma_holes,
                use_lone_pairs=use_lone_pairs,
                dielectric=dielectric
            )
        else:
            self.enhanced_calc = None
            if use_sigma_holes or use_lone_pairs:
                print("[GasteigerChargeCalculator] ⚠️ Advanced patches not available, using standard charges")
    
    def calculate_charges(self, mol: 'Chem.Mol') -> np.ndarray:
        """
        Calculate Gasteiger charges for a molecule.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Array of partial charges
        """
        if not RDKIT_AVAILABLE:
            raise RuntimeError("RDKit required for charge calculation")
        
        # Compute Gasteiger charges
        AllChem.ComputeGasteigerCharges(mol)
        
        charges = []
        for atom in mol.GetAtoms():
            charge = atom.GetDoubleProp('_GasteigerCharge')
            # Handle NaN values (can occur for some atoms)
            if np.isnan(charge):
                charge = 0.0
            charges.append(charge)
        
        return np.array(charges)
    
    def calculate_potential(self, mol: 'Chem.Mol', points: np.ndarray,
                           conformer_id: int = 0) -> np.ndarray:
        """
        Calculate electrostatic potential at given points using Coulomb's law.
        
        φ(r) = Σ q_i / (ε * |r - r_i|)
        
        If σ-hole or lone pair patches are enabled, virtual points are included
        in the potential calculation.
        
        Args:
            mol: RDKit molecule with Gasteiger charges
            points: Nx3 array of points to evaluate
            conformer_id: Which conformer to use
            
        Returns:
            Array of potential values at each point
        """
        if not RDKIT_AVAILABLE:
            raise RuntimeError("RDKit required for potential calculation")
        
        # Use enhanced calculator if available
        if self.enhanced_calc is not None:
            return self.enhanced_calc.calculate_potential(mol, points, conformer_id)
        
        # Standard calculation (no virtual points)
        # Get charges
        charges = self.calculate_charges(mol)
        
        # Get atom coordinates
        conf = mol.GetConformer(conformer_id)
        coords = []
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            coords.append([pos.x, pos.y, pos.z])
        coords = np.array(coords)
        
        # Calculate Coulomb potential at each point
        # Using kT/e units (multiply by 332.0637 for kcal/mol·e)
        # Here we use a simplified form
        potentials = np.zeros(len(points))
        
        for i, point in enumerate(points):
            distances = np.linalg.norm(coords - point, axis=1)
            # Avoid division by zero
            distances = np.maximum(distances, 0.1)
            potentials[i] = np.sum(charges / (self.dielectric * distances))
        
        # Convert to kT/e units (approximate)
        # 1 kT/e ≈ 0.0257 V at 298K
        # Coulomb constant k = 332.0637 kcal·Å/(mol·e²)
        # kT at 298K ≈ 0.593 kcal/mol
        # So 1 kT/e ≈ 0.593/332.0637 ≈ 0.00179 in our units
        # We scale to match APBS output range
        potentials *= 332.0637 / 0.593  # Convert to kT/e
        
        return potentials
    
    def calculate_potential_from_coords(self, coords: np.ndarray, charges: np.ndarray,
                                        points: np.ndarray) -> np.ndarray:
        """
        Calculate potential from raw coordinates and charges.

        Args:
            coords: Nx3 array of atom coordinates
            charges: N array of partial charges
            points: Mx3 array of evaluation points

        Returns:
            M array of potential values
        """
        potentials = np.zeros(len(points))

        for i, point in enumerate(points):
            distances = np.linalg.norm(coords - point, axis=1)
            distances = np.maximum(distances, 0.1)
            potentials[i] = np.sum(charges / (self.dielectric * distances))

        potentials *= 332.0637 / 0.593

        return potentials

    def calculate_potential_from_pqr(self, pqr_file: str, points: np.ndarray) -> np.ndarray:
        """
        Calculate potential from PQR file (fallback when APBS unavailable).

        Args:
            pqr_file: Path to PQR file
            points: Mx3 array of evaluation points

        Returns:
            M array of potential values
        """
        coords = []
        charges = []

        try:
            with open(pqr_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        # PQR format: x, y, z, charge, radius
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        charge = float(line[54:62])

                        coords.append([x, y, z])
                        charges.append(charge)
        except Exception as e:
            print(f"[calculate_potential_from_pqr] Error reading PQR: {e}")
            return np.zeros(len(points))

        if not coords:
            print(f"[calculate_potential_from_pqr] No atoms found in {pqr_file}")
            return np.zeros(len(points))

        coords = np.array(coords)
        charges = np.array(charges)

        return self.calculate_potential_from_coords(coords, charges, points)


# Try to import pdb2pqr Python module
try:
    import pdb2pqr
    from pdb2pqr import run as pdb2pqr_run
    from pdb2pqr.main import main_driver as pdb2pqr_main
    PDB2PQR_PYTHON_AVAILABLE = True
except ImportError:
    PDB2PQR_PYTHON_AVAILABLE = False


class PDB2PQRRunner:
    """
    Interface to run PDB2PQR for protein structure preparation.
    
    PDB2PQR adds hydrogens, assigns charges and radii based on force field.
    Supports both Python API and command-line interface.
    """
    
    def __init__(self, pdb2pqr_path: str = None, force_field: str = 'AMBER'):
        """
        Args:
            pdb2pqr_path: Path to pdb2pqr executable (auto-detect if None)
            force_field: Force field to use (AMBER, CHARMM, PARSE, etc.)
        """
        self.pdb2pqr_path = pdb2pqr_path or self._find_pdb2pqr()
        self.force_field = force_field
        self.use_python_api = PDB2PQR_PYTHON_AVAILABLE
    
    def _find_pdb2pqr(self) -> str:
        """Try to find pdb2pqr in PATH or common locations."""
        # Try common names
        for name in ['pdb2pqr', 'pdb2pqr30', 'pdb2pqr.py']:
            path = shutil.which(name)
            if path:
                return path

        # Try Python module
        if PDB2PQR_PYTHON_AVAILABLE:
            return 'python_api'

        # Try conda environment
        conda_env = os.environ.get('CONDA_PREFIX')
        if conda_env:
            conda_path = os.path.join(conda_env, 'bin', 'pdb2pqr')
            if os.path.exists(conda_path):
                return conda_path

        # Try common conda locations
        for conda_dir in [os.path.expanduser('~/miniconda3'), os.path.expanduser('~/anaconda3')]:
            if os.path.exists(conda_dir):
                conda_path = os.path.join(conda_dir, 'bin', 'pdb2pqr')
                if os.path.exists(conda_path):
                    return conda_path

        # Try homebrew on macOS
        homebrew_path = '/usr/local/opt/pdb2pqr/bin/pdb2pqr'
        if os.path.exists(homebrew_path):
            return homebrew_path

        return 'pdb2pqr'  # Hope it's in PATH
    
    def run(self, input_pdb: str, output_pqr: str, ph: float = 7.4,
            keep_chain: bool = True, remove_water: bool = True) -> bool:
        """
        Run PDB2PQR on a PDB file.
        
        Args:
            input_pdb: Input PDB file path
            output_pqr: Output PQR file path
            ph: pH for protonation state assignment
            keep_chain: Keep chain IDs in output
            remove_water: Remove water molecules
            
        Returns:
            True if successful
        """
        # Try Python API first (faster and more reliable)
        if self.use_python_api and PDB2PQR_PYTHON_AVAILABLE:
            return self._run_python_api(input_pdb, output_pqr, ph, keep_chain, remove_water)
        
        # Fall back to command line
        return self._run_command_line(input_pdb, output_pqr, ph, keep_chain, remove_water)
    
    def _run_python_api(self, input_pdb: str, output_pqr: str, ph: float = 7.4,
                        keep_chain: bool = True, remove_water: bool = True) -> bool:
        """Run PDB2PQR using Python API."""
        print(f"[PDB2PQR] Running via Python API...")
        print(f"[PDB2PQR] Input: {input_pdb}")
        print(f"[PDB2PQR] Output: {output_pqr}")
        print(f"[PDB2PQR] Force field: {self.force_field}, pH: {ph}")
        
        try:
            # Build arguments for pdb2pqr
            args = [
                input_pdb,
                output_pqr,
                f'--ff={self.force_field}',
                f'--with-ph={ph}',
                '--titration-state-method=propka',
            ]
            
            if keep_chain:
                args.append('--keep-chain')
            
            if remove_water:
                args.append('--drop-water')
            
            # Run pdb2pqr
            from pdb2pqr.main import main_driver, build_main_parser
            
            parser = build_main_parser()
            parsed_args = parser.parse_args(args)
            main_driver(parsed_args)
            
            if os.path.exists(output_pqr):
                # 修复 PQR 格式：将固定列格式转为空格分隔，防止 APBS 解析Failed
                _fix_pqr_format(output_pqr)
                print(f"[PDB2PQR] ✅ Generated: {output_pqr}")
                return True
            else:
                print(f"[PDB2PQR] ❌ Output file not created")
                return False
                
        except Exception as e:
            print(f"[PDB2PQR] Python API error: {e}")
            print("[PDB2PQR] Falling back to command line...")
            return self._run_command_line(input_pdb, output_pqr, ph, keep_chain, remove_water)
    
    def _run_command_line(self, input_pdb: str, output_pqr: str, ph: float = 7.4,
                          keep_chain: bool = True, remove_water: bool = True) -> bool:
        """Run PDB2PQR using command line."""
        cmd_parts = [self.pdb2pqr_path]
        
        # Add options
        cmd_parts.extend(['--ff=' + self.force_field])
        cmd_parts.extend(['--with-ph=' + str(ph)])
        
        if keep_chain:
            cmd_parts.append('--keep-chain')
        
        if remove_water:
            cmd_parts.append('--drop-water')
        
        # Input and output
        cmd_parts.append(input_pdb)
        cmd_parts.append(output_pqr)
        
        # Run
        cmd_str = ' '.join(cmd_parts)
        print(f"[PDB2PQR] Running: {cmd_str}")
        
        try:
            result = subprocess.run(
                cmd_parts,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode != 0:
                # Check for "no such option: --keep-chain" and retry without it
                if "--keep-chain" in result.stderr and "no such option" in result.stderr:
                    print("[PDB2PQR] ⚠️ '--keep-chain' not supported by this version, retrying without it...")
                    new_cmd_parts = [p for p in cmd_parts if p != '--keep-chain']
                    result = subprocess.run(
                        new_cmd_parts,
                        capture_output=True,
                        text=True,
                        timeout=300
                    )
                    if result.returncode == 0 and os.path.exists(output_pqr):
                        # 修复 PQR 格式：将固定列格式转为空格分隔，防止 APBS 解析Failed
                        _fix_pqr_format(output_pqr)
                        print(f"[PDB2PQR] ✅ Generated (without keep-chain): {output_pqr}")
                        return True
                
                print(f"[PDB2PQR] Error: {result.stderr}")
                return False
            
            if os.path.exists(output_pqr):
                # 修复 PQR 格式：将固定列格式转为空格分隔，防止 APBS 解析Failed
                _fix_pqr_format(output_pqr)
                print(f"[PDB2PQR] ✅ Generated: {output_pqr}")
                return True
            else:
                print(f"[PDB2PQR] ❌ Output file not created")
                return False
                
        except subprocess.TimeoutExpired:
            print("[PDB2PQR] ❌ Timeout")
            return False
        except FileNotFoundError:
            print(f"[PDB2PQR] ❌ pdb2pqr not found at: {self.pdb2pqr_path}")
            print("[PDB2PQR] 💡 Install with: pip install pdb2pqr")
            return False
        except Exception as e:
            print(f"[PDB2PQR] ❌ Error: {e}")
            return False


# 尝试Import APBS 原生 Python Module（如 conda install apbs）
APBS_PYTHON_AVAILABLE = False
_apbs_module = None
try:
    import apbs
    _apbs_module = apbs
    APBS_PYTHON_AVAILABLE = True
except ImportError:
    pass

# 尝试Import apbs-binary pip Package（安装器通过 pip install apbs-binary 安装）
# 该Package提供 run_apbs / popen_apbs Function，内部调用打Package的 apbs 可执行File
APBS_BINARY_AVAILABLE = False
_apbs_binary_run = None
try:
    from apbs_binary import run_apbs as _apbs_binary_run_fn
    _apbs_binary_run = _apbs_binary_run_fn
    APBS_BINARY_AVAILABLE = True
except ImportError:
    pass


class APBSRunner:
    """
    Interface to run APBS for electrostatic potential calculation.
    
    APBS solves the Poisson-Boltzmann equation to compute electrostatic
    potential around a molecule.
    
    Supports both Python API and command-line interface.
    """
    
    def __init__(self, apbs_path: str = None):
        """
        Args:
            apbs_path: Path to APBS executable (auto-detect if None)
        """
        self.apbs_path = apbs_path or self._find_apbs()
        self.use_python_api = APBS_PYTHON_AVAILABLE
        # 标记 apbs-binary pip Package是否可用，作为第二优先级运行方式
        self.use_apbs_binary = APBS_BINARY_AVAILABLE
    
    def _find_apbs(self) -> str:
        """尝试在 PATH、conda 环境、apbs-binary Package等多个位置Find APBS。"""
        # 1. 系统 PATH 中Find
        path = shutil.which('apbs')
        if path:
            return path

        # 2. 原生 Python API 可用时直接using
        if APBS_PYTHON_AVAILABLE:
            return 'python_api'

        # 3. 通过 apbs-binary pip PackageLocate其内置的 apbs 可执行File
        if APBS_BINARY_AVAILABLE:
            try:
                import apbs_binary
                pkg_dir = os.path.dirname(apbs_binary.__file__)
                # apbs-binary Package通常在PackageDirectory或 bin 子Directory下放置可执行File
                for candidate in [
                    os.path.join(pkg_dir, 'apbs'),
                    os.path.join(pkg_dir, 'apbs.exe'),
                    os.path.join(pkg_dir, 'bin', 'apbs'),
                    os.path.join(pkg_dir, 'bin', 'apbs.exe'),
                ]:
                    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
                        return candidate
            except Exception:
                pass

        # 4. 通过当前 Python 解释器Path推断 conda 环境的 bin Directory
        #    从 .app 束启动 PyMOL 时 CONDA_PREFIX 可能未Settings，
        #    但 sys.executable 仍指向 conda 环境内的 python
        try:
            python_dir = os.path.dirname(os.path.realpath(sys.executable))
            # python 通常在 envs/glint/bin/ 下
            inferred_apbs = os.path.join(python_dir, 'apbs')
            if os.path.isfile(inferred_apbs) and os.access(inferred_apbs, os.X_OK):
                return inferred_apbs
        except Exception:
            pass

        # 5. CONDA_PREFIX 环境变量（终端启动时可用）
        conda_env = os.environ.get('CONDA_PREFIX')
        if conda_env:
            conda_path = os.path.join(conda_env, 'bin', 'apbs')
            if os.path.exists(conda_path):
                return conda_path

        # 6. GLINT 安装器Create的 glint conda 环境（硬编码Path兜底）
        for conda_root in [os.path.expanduser('~/miniconda3'), os.path.expanduser('~/anaconda3')]:
            if os.path.exists(conda_root):
                # 优先检查 glint 虚拟环境
                glint_env_apbs = os.path.join(conda_root, 'envs', 'glint', 'bin', 'apbs')
                if os.path.exists(glint_env_apbs):
                    return glint_env_apbs
                # 再检查 base 环境
                base_apbs = os.path.join(conda_root, 'bin', 'apbs')
                if os.path.exists(base_apbs):
                    return base_apbs

        # 7. macOS Homebrew Path
        for homebrew_path in ['/usr/local/opt/apbs/bin/apbs', '/opt/homebrew/bin/apbs']:
            if os.path.exists(homebrew_path):
                return homebrew_path

        # 8. 其他常见安装Path
        common_paths = [
            '/usr/local/bin/apbs',
            '/usr/bin/apbs',
            os.path.expanduser('~/apbs/bin/apbs'),
            'C:\\Program Files\\APBS\\apbs.exe'
        ]

        for p in common_paths:
            if os.path.exists(p):
                return p

        # 如果 apbs-binary 可用，可以通过其 run_apbs Function运行，无需可执行FilePath
        if APBS_BINARY_AVAILABLE:
            return 'apbs_binary_api'

        return 'apbs'  # 兜底：期望它在 PATH 中
    
    def generate_input(self, pqr_file: str, output_prefix: str,
                      grid_center: np.ndarray = None,
                      grid_dims: Tuple[int, int, int] = (161, 161, 161),
                      coarse_len: Tuple[float, float, float] = (60, 60, 60),
                      fine_len: Tuple[float, float, float] = (40, 40, 40)) -> str:
        """
        Generate APBS input file.
        
        Args:
            pqr_file: Input PQR file
            output_prefix: Prefix for output files
            grid_center: Center of grid (use molecule center if None)
            grid_dims: Grid dimensions
            coarse_len: Coarse grid lengths (Å)
            fine_len: Fine grid lengths (Å)
            
        Returns:
            Path to generated input file
        """
        # 网格中心：如提供 grid_center 则以坐标Locate，否则以分子质心Locate
        if grid_center is not None:
            cgcent_line = f"  cgcent {grid_center[0]:.4f} {grid_center[1]:.4f} {grid_center[2]:.4f}"
            fgcent_line = f"  fgcent {grid_center[0]:.4f} {grid_center[1]:.4f} {grid_center[2]:.4f}"
            print(f"[APBS] using自定义网格中心: ({grid_center[0]:.2f}, {grid_center[1]:.2f}, {grid_center[2]:.2f})")
        else:
            cgcent_line = "  cgcent mol 1"
            fgcent_line = "  fgcent mol 1"
        
        input_content = f"""# APBS input file generated by GLINT
read
  mol pqr {pqr_file}
end

elec
  mg-auto
  mol 1
  lpbe
  bcfl sdh
  pdie {APBS_PARAMS['pdie']}
  sdie {APBS_PARAMS['sdie']}
  chgm spl2
  srfm smol
  srad {PROBE_RADIUS}
  swin 0.3
  sdens 10.0
  temp {APBS_PARAMS['temp']}
  calcenergy no
  calcforce no

  dime {grid_dims[0]} {grid_dims[1]} {grid_dims[2]}
  cglen {coarse_len[0]} {coarse_len[1]} {coarse_len[2]}
  fglen {fine_len[0]} {fine_len[1]} {fine_len[2]}
{cgcent_line}
{fgcent_line}

  ion charge {APBS_PARAMS['ion_charge']} conc {APBS_PARAMS['ion_conc']} radius {APBS_PARAMS['ion_radius']}
  ion charge -{APBS_PARAMS['ion_charge']} conc {APBS_PARAMS['ion_conc']} radius {APBS_PARAMS['ion_radius']}

  write pot dx {output_prefix}
end
quit
"""
        
        input_file = output_prefix + '.in'
        with open(input_file, 'w') as f:
            f.write(input_content)
        
        print(f"[APBS] Generated input file: {input_file}")
        return input_file
    
    def run(self, input_file: str, working_dir: str = None) -> Optional[str]:
        """
        运行 APBS 计算。
        
        优先级：原生 Python API > apbs-binary Package > CLI 命令行
        
        Args:
            input_file: APBS 输入File
            working_dir: 工作Directory（默认using输入File所在Directory）
            
        Returns:
            输出 DX FilePath，Failed则Return None
        """
        if working_dir is None:
            working_dir = os.path.dirname(input_file) or '.'
        
        # 方式1: 原生 Python API（如 conda install apbs 提供的Module）
        if self.use_python_api and APBS_PYTHON_AVAILABLE:
            result = self._run_python_api(input_file, working_dir)
            if result:
                return result
            print("[APBS] Python API 执行Failed，尝试其他方式...")
        
        # 方式2: apbs-binary pip Package（安装器通过 pip install apbs-binary 安装）
        if self.use_apbs_binary and APBS_BINARY_AVAILABLE:
            result = self._run_apbs_binary(input_file, working_dir)
            if result:
                return result
            print("[APBS] apbs-binary Package执行Failed，尝试命令行...")
        
        # 方式3: 命令行 CLI（直接调用 apbs 可执行File）
        return self._run_command_line(input_file, working_dir)
    
    def _run_python_api(self, input_file: str, working_dir: str) -> Optional[str]:
        """Run APBS using Python API."""
        print(f"[APBS] Running via Python API...")
        print(f"[APBS] Input file: {input_file}")
        
        try:
            # Read input file content
            with open(input_file, 'r') as f:
                input_content = f.read()
            
            # Try different APBS Python API approaches
            # Method 1: Using apbs.run() if available
            if hasattr(apbs, 'run'):
                old_cwd = os.getcwd()
                try:
                    os.chdir(working_dir)
                    apbs.run(input_file)
                finally:
                    os.chdir(old_cwd)
            
            # Method 2: Using apbs.Apbs class if available
            elif hasattr(apbs, 'Apbs'):
                solver = apbs.Apbs(input_file)
                solver.run()
            
            # Method 3: Using apbs.input_file if available
            elif hasattr(apbs, 'input_file'):
                old_cwd = os.getcwd()
                try:
                    os.chdir(working_dir)
                    apbs.input_file(input_file)
                finally:
                    os.chdir(old_cwd)
            
            else:
                print("[APBS] Python API available but no known interface found")
                return None
            
            # Find output DX file
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            dx_file = os.path.join(working_dir, base_name + '.dx')
            
            if os.path.exists(dx_file):
                print(f"[APBS] ✅ Generated: {dx_file}")
                return dx_file
            
            # Try alternative naming
            for f in os.listdir(working_dir):
                if f.endswith('.dx'):
                    dx_file = os.path.join(working_dir, f)
                    print(f"[APBS] ✅ Found output: {dx_file}")
                    return dx_file
            
            print("[APBS] ❌ No output DX file found")
            return None
            
        except Exception as e:
            print(f"[APBS] Python API error: {e}")
            return None
    
    def _run_apbs_binary(self, input_file: str, working_dir: str) -> Optional[str]:
        """通过 apbs-binary pip Package运行 APBS。
        
        apbs-binary Package提供 run_apbs() Function，内部调用打Package的 apbs 可执行File。
        Return subprocess.CompletedProcess 对象。
        """
        print("[APBS] 通过 apbs-binary Package运行...")
        print(f"[APBS] 输入File: {input_file}")
        
        try:
            from apbs_binary import run_apbs
            
            # run_apbs 在当前工作Directory下运行，需要切换到工作Directory
            old_cwd = os.getcwd()
            try:
                os.chdir(working_dir)
                # run_apbs 接受输入FilePath，Return subprocess.CompletedProcess
                completed = run_apbs(input_file)
                
                if completed.returncode != 0:
                    stderr_text = completed.stderr if hasattr(completed, 'stderr') and completed.stderr else ''
                    print(f"[APBS] apbs-binary Return非零Exit码: {completed.returncode}")
                    if stderr_text:
                        print(f"[APBS] Error输出:\n{stderr_text}")
                    # APBS 有时Return非零但仍生成输出，ContinueFind
            finally:
                os.chdir(old_cwd)
            
            # Find输出 DX File（复用与其他方式相同的Find逻辑）
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            dx_file = os.path.join(working_dir, base_name + '.dx')
            
            if os.path.exists(dx_file):
                print(f"[APBS] ✅ 生成: {dx_file}")
                return dx_file
            
            # 尝试其他命名方式（APBS 可能using不同的输出File名）
            for f in os.listdir(working_dir):
                if f.endswith('.dx'):
                    dx_file = os.path.join(working_dir, f)
                    print(f"[APBS] ✅ 找到输出: {dx_file}")
                    return dx_file
            
            print("[APBS] ❌ 未找到输出 DX File")
            return None
            
        except Exception as e:
            print(f"[APBS] apbs-binary 执行Error: {e}")
            return None
    
    def _run_command_line(self, input_file: str, working_dir: str) -> Optional[str]:
        """Run APBS using command line."""
        cmd_parts = [self.apbs_path, input_file]

        print(f"[APBS] Running: {' '.join(cmd_parts)}")
        print(f"[APBS] Working directory: {working_dir}")

        # 为 APBS 进程构造运行时库路径，修复 macOS 下 libmaloc 等动态库找不到的问题
        run_env = os.environ.copy()
        lib_dirs = []
        apbs_real = os.path.realpath(self.apbs_path)
        apbs_bin_dir = os.path.dirname(apbs_real)

        # 优先：apbs 所在环境的 lib / Frameworks
        for d in [
            os.path.join(apbs_bin_dir, '..', 'lib'),
            os.path.join(apbs_bin_dir, '..', 'Frameworks'),
        ]:
            d = os.path.realpath(d)
            if os.path.isdir(d):
                lib_dirs.append(d)

        # 次优先：当前 conda 环境（若存在）
        conda_prefix = run_env.get('CONDA_PREFIX')
        if conda_prefix:
            for d in [
                os.path.join(conda_prefix, 'lib'),
                os.path.join(conda_prefix, 'Frameworks'),
            ]:
                d = os.path.realpath(d)
                if os.path.isdir(d):
                    lib_dirs.append(d)

        # 去重并合并到 DYLD_LIBRARY_PATH / LD_LIBRARY_PATH
        dedup_lib_dirs = []
        for d in lib_dirs:
            if d not in dedup_lib_dirs:
                dedup_lib_dirs.append(d)

        if dedup_lib_dirs:
            # macOS 动态库搜索路径
            old_dyld = run_env.get('DYLD_LIBRARY_PATH', '')
            run_env['DYLD_LIBRARY_PATH'] = ':'.join(dedup_lib_dirs + ([old_dyld] if old_dyld else []))
            # 兼容 Linux 路径变量（无副作用）
            old_ld = run_env.get('LD_LIBRARY_PATH', '')
            run_env['LD_LIBRARY_PATH'] = ':'.join(dedup_lib_dirs + ([old_ld] if old_ld else []))

        try:
            result = subprocess.run(
                cmd_parts,
                capture_output=True,
                text=True,
                cwd=working_dir,
                timeout=600,  # 10 minute timeout
                env=run_env,
            )

            if result.returncode != 0:
                print(f"[APBS] Error output:\n{result.stderr}")
                # APBS sometimes returns non-zero but still produces output

            # Find output DX file
            # APBS adds .dx extension to the prefix
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            dx_file = os.path.join(working_dir, base_name + '.dx')

            if os.path.exists(dx_file):
                print(f"[APBS] ✅ Generated: {dx_file}")
                return dx_file

            # Try alternative naming
            for f in os.listdir(working_dir):
                if f.endswith('.dx'):
                    dx_file = os.path.join(working_dir, f)
                    print(f"[APBS] ✅ Found output: {dx_file}")
                    return dx_file

            print("[APBS] ❌ No output DX file found")
            return None

        except subprocess.TimeoutExpired:
            print("[APBS] ❌ Timeout")
            return None
        except FileNotFoundError:
            print(f"[APBS] ❌ APBS not found at: {self.apbs_path}")
            print("[APBS] 💡 Install APBS from: https://github.com/Electrostatics/apbs")
            return None
        except Exception as e:
            print(f"[APBS] ❌ Error: {e}")
            return None


class ECCalculator:
    """
    Main class for Electrostatic Complementarity (EC) calculation.
    
    EC measures how well the electrostatic potential of a ligand matches
    (complements) the potential from the protein at the ligand surface.
    
    EC_i = -2 * φ_protein * φ_ligand / (φ_protein² + φ_ligand² + ε)
    
    - Opposite signs (complementary): EC > 0 (protein provides what ligand surface needs)
    - Same signs (clash): EC < 0 (electrostatic repulsion)
    - Perfect match: EC → +1
    """
    
    def __init__(self, epsilon: float = 1e-6, clip_range: Tuple[float, float] = (-10.0, 10.0)):
        """
        Args:
            epsilon: Small value to avoid division by zero
            clip_range: (min, max) for potential clipping
        """
        self.epsilon = epsilon
        self.clip_min, self.clip_max = clip_range
    
    def calculate_ec_local(self, phi_protein: np.ndarray, phi_ligand: np.ndarray) -> np.ndarray:
        """
        Calculate local EC values at each surface point.
        
        Args:
            phi_protein: Protein potential at surface points
            phi_ligand: Ligand potential at surface points
            
        Returns:
            Array of EC values at each point
        """
        # Clip potentials to avoid extreme values
        phi_p = np.clip(phi_protein, self.clip_min, self.clip_max)
        phi_l = np.clip(phi_ligand, self.clip_min, self.clip_max)
        
        # Calculate EC
        
        # EC formula: EC = -2 * φ_protein * φ_ligand / (φ_protein² + φ_ligand² + ε)
        # When protein and ligand potentials have OPPOSITE signs, EC > 0 (complementary)
        # Opposite signs mean protein provides what ligand surface needs
        numerator = -2.0 * phi_p * phi_l  # Positive when opposite signs (complementary)
        denominator = phi_p ** 2 + phi_l ** 2 + self.epsilon
        
        ec = numerator / denominator
        
        return ec
    
    def calculate_ec_score(self, ec_local: np.ndarray, weights: np.ndarray = None) -> float:
        """
        Calculate overall EC score from local values.
        
        Args:
            ec_local: Array of local EC values
            weights: Optional weights for each point (e.g., surface area)
            
        Returns:
            Scalar EC score
        """
        # Exclude points where EC == 0 (masked by distance cutoff)
        nonzero_mask = ec_local != 0
        if np.sum(nonzero_mask) == 0:
            return 0.0  # All points masked
        
        ec_nonzero = ec_local[nonzero_mask]
        print(f"[EC Score] Using {len(ec_nonzero)}/{len(ec_local)} non-zero points")
        
        if weights is None:
            return float(np.mean(ec_nonzero))
        else:
            weights_nonzero = weights[nonzero_mask] if weights is not None else None
            return float(np.average(ec_nonzero, weights=weights_nonzero))
    
    def calculate_ec_statistics(self, ec_local: np.ndarray, exclude_zeros: bool = False) -> Dict[str, float]:
        """
        Calculate various statistics for EC distribution.
        
        Args:
            ec_local: Array of local EC values
            exclude_zeros: 是否排除零Value点（三元复合物 mask 后的点）。默认 False 保持向后兼容。
            
        Returns:
            Dictionary of statistics
        """
        # 修复: 当 exclude_zeros=True 时，先Filter掉被 mask 为 0 的点，避免统计被稀释
        if exclude_zeros:
            nonzero_mask = ec_local != 0
            if np.sum(nonzero_mask) == 0:
                return {
                    'ec_mean': 0.0, 'ec_median': 0.0, 'ec_std': 0.0,
                    'ec_min': 0.0, 'ec_max': 0.0,
                    'ec_positive_fraction': 0.0, 'ec_negative_fraction': 0.0,
                    'ec_q25': 0.0, 'ec_q75': 0.0
                }
            data = ec_local[nonzero_mask]
        else:
            data = ec_local
        
        return {
            'ec_mean': float(np.mean(data)),
            'ec_median': float(np.median(data)),
            'ec_std': float(np.std(data)),
            'ec_min': float(np.min(data)),
            'ec_max': float(np.max(data)),
            'ec_positive_fraction': float(np.mean(data > 0)),
            'ec_negative_fraction': float(np.mean(data < 0)),
            'ec_q25': float(np.percentile(data, 25)),
            'ec_q75': float(np.percentile(data, 75))
        }


class ECMapWriter:
    """
    Write EC values to various formats for visualization.
    """
    
    @staticmethod
    def write_pseudo_pdb(filename: str, points: np.ndarray, values: np.ndarray,
                        scale: float = 100.0, atom_name: str = 'EC'):
        """
        Write surface points as pseudo-PDB with EC values in B-factor.
        
        Args:
            filename: Output PDB file
            points: Nx3 array of coordinates
            values: N array of EC values
            scale: Scale factor for B-factor (EC * scale)
            atom_name: Atom name to use
        """
        with open(filename, 'w') as f:
            f.write("REMARK EC map generated by GLINT\n")
            f.write(f"REMARK EC values scaled by {scale} in B-factor column\n")
            
            for i, (point, value) in enumerate(zip(points, values)):
                # Scale EC to B-factor range
                bfactor = value * scale
                bfactor = max(-99.99, min(99.99, bfactor))  # PDB B-factor limits
                
                f.write(
                    f"HETATM{i+1:5d}  {atom_name:3s} ECM A   1    "
                    f"{point[0]:8.3f}{point[1]:8.3f}{point[2]:8.3f}"
                    f"  1.00{bfactor:6.2f}           C\n"
                )
            
            f.write("END\n")
        
        print(f"[ECMapWriter] Wrote {len(points)} points to {filename}")
    
    @staticmethod
    def write_csv(filename: str, points: np.ndarray, ec_values: np.ndarray,
                 phi_protein: np.ndarray = None, phi_ligand: np.ndarray = None):
        """
        Write EC data to CSV file.
        
        Args:
            filename: Output CSV file
            points: Nx3 array of coordinates
            ec_values: N array of EC values
            phi_protein: Optional protein potential values
            phi_ligand: Optional ligand potential values
        """
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Header
            header = ['x', 'y', 'z', 'EC']
            if phi_protein is not None:
                header.append('phi_protein')
            if phi_ligand is not None:
                header.append('phi_ligand')
            writer.writerow(header)
            
            # Data
            for i in range(len(points)):
                row = [points[i, 0], points[i, 1], points[i, 2], ec_values[i]]
                if phi_protein is not None:
                    row.append(phi_protein[i])
                if phi_ligand is not None:
                    row.append(phi_ligand[i])
                writer.writerow(row)
        
        print(f"[ECMapWriter] Wrote EC data to {filename}")


# ========== Main Analysis Functions ==========

def calculate_ligand_ec(obj_name: str = None, ligand_resname: str = None,
                       protein_chains: List[str] = None,
                       output_dir: str = None,
                       ph: float = 7.4,
                       surface_density: float = 10.0,
                       visualize: bool = True,
                       pdb_file: str = None,
                       use_sigma_holes: bool = False,
                       use_lone_pairs: bool = False,
                       keep_bridging_waters: bool = False,
                       contact_cutoff: float = 10.0) -> Optional[Dict[str, Any]]:
    """
    Calculate Electrostatic Complementarity for a protein-ligand complex.
    
    This is the main entry point for EC analysis. It:
    1. Extracts protein and ligand from the complex
    2. Runs PDB2PQR on the protein
    3. Runs APBS to calculate protein electrostatic potential
    4. Generates ligand surface points
    5. Calculates ligand potential using Gasteiger charges
    6. Computes EC at each surface point
    7. Generates EC map for visualization
    
    Args:
        obj_name: PyMOL object name
        ligand_resname: Ligand residue name (e.g., 'LIG')
        protein_chains: List of protein chain IDs (auto-detect if None)
        output_dir: Directory for output files (temp dir if None)
        ph: pH for protonation state
        surface_density: Surface sampling density (points/Å²)
        visualize: Whether to visualize results in PyMOL
        pdb_file: Alternative: use PDB file instead of PyMOL object
        use_sigma_holes: Add σ-hole virtual points for Cl/Br/I (improves halogen bond EC)
        use_lone_pairs: Add lone pair virtual points for carbonyl O (improves H-bond EC)
        keep_bridging_waters: Keep structurally important bridging waters
        contact_cutoff: 接触区截断距离 (Å)，仅提取配体周围此距离内的蛋白残基进行计算（默认 10.0）
        
    Returns:
        Dictionary containing:
        - ec_score: Overall EC score
        - ec_statistics: Detailed statistics
        - surface_points: Surface point coordinates
        - ec_values: EC value at each point
        - output_files: Paths to generated files
        - advanced_features: Info about σ-hole/lone pair usage
    """
    # Check dependencies
    if not NUMPY_AVAILABLE:
        print("[calculate_ligand_ec] ❌ NumPy required")
        return None
    
    if not RDKIT_AVAILABLE:
        print("[calculate_ligand_ec] ❌ RDKit required")
        return None
    

    if output_dir: output_dir = os.path.normpath(output_dir)

    # Create output directory
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='glint_ec_')
    else:
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"[calculate_ligand_ec] Output directory: {output_dir}")
    
    # Step 1: Extract protein and ligand
    print("\n[Step 1] Extracting protein and ligand...")
    
    if pdb_file:
        # Use provided PDB file
        complex_pdb = pdb_file
    elif PYMOL_AVAILABLE and obj_name:
        # Export from PyMOL
        complex_pdb = os.path.join(output_dir, 'complex.pdb')
        cmd.save(complex_pdb, obj_name)
    else:
        print("[calculate_ligand_ec] ❌ Need either PyMOL object or PDB file")
        return None
    
    # Extract protein (remove ligand, water, ions)
    protein_pdb = os.path.join(output_dir, 'protein.pdb')
    ligand_sdf = os.path.join(output_dir, 'ligand.sdf')
    
    # 追踪保留的水分子Count（用于后续 advanced_features 记录）
    n_waters_kept = 0
    
    if PYMOL_AVAILABLE and obj_name:
        # 接触区残基提取：仅提取配体周围 contact_cutoff 距离内的蛋白残基
        ligand_within_sel = f"{obj_name} and resn {ligand_resname}" if ligand_resname else f"{obj_name} and organic and not polymer"
        if protein_chains:
            chain_sel = ' or '.join([f"chain {c}" for c in protein_chains])
            protein_sel = f"byres (({obj_name} and polymer and ({chain_sel})) within {contact_cutoff} of ({ligand_within_sel}))"
        else:
            protein_sel = f"byres (({obj_name} and polymer) within {contact_cutoff} of ({ligand_within_sel}))"
        
        # Print接触区Information
        try:
            n_contact_atoms = cmd.count_atoms(protein_sel)
            # 统计残基数
            contact_residues = set()
            cmd.iterate(protein_sel, "contact_residues.add((chain, resi))",
                       space={'contact_residues': contact_residues})
            n_contact_res = len(contact_residues)
            print(f"[calculate_ligand_ec] 📍 接触区提取 (cutoff={contact_cutoff} Å):")
            print(f"  蛋白残基: {n_contact_res} 个, 原子: {n_contact_atoms} 个")
        except Exception:
            pass
        
        # 如果需要保留桥联水，using BridgingWaterFilter 筛选结构性水分子
        if keep_bridging_waters and ADVANCED_PATCHES_AVAILABLE:
            print("[calculate_ligand_ec] ✨ 筛选结构性桥联水...")
            try:
                water_filter = BridgingWaterFilter()
                # 先SavePackage含水的复合物临时File
                complex_with_water_pdb = os.path.join(output_dir, 'complex_with_water.pdb')
                cmd.save(complex_with_water_pdb, f"{obj_name} and (polymer or resn HOH)")
                
                # 筛选桥联水（ReturnValue为 (chain, resnum) 元组列表）
                water_result = water_filter.filter_pdb(complex_with_water_pdb, ligand_resname)
                bridging_waters = water_result.get('bridging_waters', [])
                protein_bound_waters = water_result.get('protein_bound_waters', [])
                
                # 合并桥联水和蛋白结合水
                waters_to_keep = set()
                for chain, resnum in bridging_waters:
                    waters_to_keep.add((chain, resnum))
                for chain, resnum in protein_bound_waters:
                    waters_to_keep.add((chain, resnum))
                
                n_waters_kept = len(waters_to_keep)
                
                if waters_to_keep:
                    # 构建Package含关Key水分子的 PyMOL Select器
                    water_resi_sel = ' or '.join(
                        [f"(chain {c} and resi {r})" for c, r in waters_to_keep]
                    )
                    water_sel = f"({obj_name} and resn HOH and ({water_resi_sel}))"
                    protein_sel = f"({protein_sel}) or ({water_sel})"
                    print(f"[calculate_ligand_ec] 保留 {n_waters_kept} 个结构性/桥联水分子")
                else:
                    print("[calculate_ligand_ec] 未发现符合条件的桥联水")
            except Exception as e:
                print(f"[calculate_ligand_ec] ⚠️ 桥联水筛选Failed: {e}")
                print("[calculate_ligand_ec] Continue分析（不含水）")
        
        cmd.save(protein_pdb, protein_sel)
        
        # Build ligand selection and validate it exists
        if ligand_resname:
            ligand_sel = f"{obj_name} and resn {ligand_resname}"
        else:
            # Auto-detect ligand (organic molecules that are not polymer)
            ligand_sel = f"{obj_name} and organic and not polymer"
        
        # Check if ligand selection has atoms
        ligand_atom_count = cmd.count_atoms(ligand_sel)
        if ligand_atom_count == 0:
            if ligand_resname:
                print(f"[calculate_ligand_ec] ❌ No atoms found for ligand '{ligand_resname}' in object '{obj_name}'")
                print(f"[calculate_ligand_ec] 💡 Available residue names in structure:")
                # List available organic residues
                try:
                    organic_residues = set()
                    cmd.iterate(f"{obj_name} and organic and not polymer",
                               "organic_residues.add(resn)",
                               space={'organic_residues': organic_residues})
                    if organic_residues:
                        print(f"[calculate_ligand_ec]    Organic residues: {', '.join(sorted(organic_residues))}")
                    else:
                        print(f"[calculate_ligand_ec]    No organic residues found in structure")
                except Exception as e:
                    print(f"[calculate_ligand_ec]    Could not list residues: {e}")
            else:
                print(f"[calculate_ligand_ec] ❌ No organic molecules found in object '{obj_name}'")
            return None
        
        print(f"[calculate_ligand_ec] Found {ligand_atom_count} atoms in ligand selection")
        
        # Save ligand to SDF format
        try:
            cmd.save(ligand_sdf, ligand_sel, format='sdf')
        except Exception as e:
            print(f"[calculate_ligand_ec] ❌ Failed to save ligand to SDF: {e}")
            # Try alternative: save as MOL2 then convert
            try:
                ligand_mol2 = os.path.join(output_dir, 'ligand.mol2')
                cmd.save(ligand_mol2, ligand_sel, format='mol2')
                print(f"[calculate_ligand_ec] Saved ligand as MOL2, converting to SDF...")
                mol = Chem.MolFromMol2File(ligand_mol2, removeHs=False)
                if mol:
                    writer = Chem.SDWriter(ligand_sdf)
                    writer.write(mol)
                    writer.close()
                    print(f"[calculate_ligand_ec] ✅ Converted MOL2 to SDF successfully")
            except Exception as e2:
                print(f"[calculate_ligand_ec] ❌ Alternative conversion also failed: {e2}")
                return None
    else:
        # Parse PDB file manually
        _extract_protein_ligand_from_pdb(complex_pdb, protein_pdb, ligand_sdf, ligand_resname)
    
    # Verify SDF file exists and is not empty
    if not os.path.exists(ligand_sdf):
        print(f"[calculate_ligand_ec] ❌ Ligand SDF file was not created: {ligand_sdf}")
        return None
    
    sdf_size = os.path.getsize(ligand_sdf)
    if sdf_size == 0:
        print(f"[calculate_ligand_ec] ❌ Ligand SDF file is empty: {ligand_sdf}")
        print(f"[calculate_ligand_ec] 💡 This usually means the ligand selection matched no atoms")
        return None
    
    print(f"[calculate_ligand_ec] Ligand SDF file size: {sdf_size} bytes")
    
    # using RDKit Load配体 SDF（含 sanitize=False 降级策略）
    ligand_mol = None
    try:
        # 优先尝试正常Load（sanitize=True，默认Value）
        supplier = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
        ligand_mol = supplier[0] if len(supplier) > 0 else None
    except Exception as e:
        print(f"[calculate_ligand_ec] ⚠️ RDKit standard read failed: {e}")
    
    # 降级策略：Close sanitize 重新Load（常见于 MD 快照或 PyMOL Export的 SDF）
    if ligand_mol is None:
        try:
            print(f"[calculate_ligand_ec] 🔄 Retrying with sanitize=False (common for MD snapshots)...")
            supplier = Chem.SDMolSupplier(ligand_sdf, removeHs=False, sanitize=False)
            ligand_mol = supplier[0] if len(supplier) > 0 else None
            if ligand_mol is not None:
                # 尝试部分清理：Skip严格的价态检查，保留其他校验
                try:
                    Chem.SanitizeMol(ligand_mol, 
                        sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
                    print(f"[calculate_ligand_ec] ✅ Loaded with partial sanitization (valence check skipped)")
                except Exception:
                    # 部分清理也Failed，using未清理的分子（坐标和原子Type仍然可用）
                    print(f"[calculate_ligand_ec] ⚠️ Partial sanitization failed, using unsanitized molecule")
                    print(f"[calculate_ligand_ec] 💡 Coordinates and atom types are valid; charges may be approximate")
        except Exception as e2:
            print(f"[calculate_ligand_ec] ❌ Fallback read also failed: {e2}")
    
    # 最终检查：如果仍然无法Load，输出调试Information并Return
    if ligand_mol is None:
        print(f"[calculate_ligand_ec] ❌ Failed to load ligand from SDF file")
        print(f"[calculate_ligand_ec] 💡 The SDF file may have invalid format")
        try:
            with open(ligand_sdf, 'r') as f:
                content = f.read(500)
                print(f"[calculate_ligand_ec] SDF file content preview:\n{content}")
        except OSError:
            pass
        return None
    
    print(f"[calculate_ligand_ec] Ligand: {ligand_mol.GetNumAtoms()} atoms")
    
    # Step 2: Run PDB2PQR
    # 如果保留了桥联水，不要让 PDB2PQR Delete它们
    print("\n[Step 2] Running PDB2PQR...")
    protein_pqr = os.path.join(output_dir, 'protein.pqr')
    
    should_remove_water = not (keep_bridging_waters and ADVANCED_PATCHES_AVAILABLE and n_waters_kept > 0)
    pdb2pqr = PDB2PQRRunner()
    if not pdb2pqr.run(protein_pdb, protein_pqr, ph=ph, remove_water=should_remove_water):
        print("[calculate_ligand_ec] ⚠️ PDB2PQR failed, trying alternative...")
        # Try to continue without PDB2PQR (use simple charge assignment)
        protein_pqr = _simple_pdb_to_pqr(protein_pdb, protein_pqr)
        if protein_pqr is None:
            return None
    
    # Step 3: Run APBS
    print("\n[Step 3] Running APBS...")
    apbs = APBSRunner()
    
    # Determine grid size based on ligand position
    ligand_center = _get_molecule_center(ligand_mol)
    
    apbs_prefix = os.path.join(output_dir, 'protein_pot')
    apbs_input = apbs.generate_input(
        protein_pqr, apbs_prefix,
        grid_center=ligand_center,
        grid_dims=(161, 161, 161),
        coarse_len=(60, 60, 60),
        fine_len=(40, 40, 40)
    )
    
    dx_file = apbs.run(apbs_input, output_dir)

    # 追踪是否using了 Coulomb 近似回退
    _coulomb_fallback = False
    _coulomb_fallback_reason = ""
    
    if dx_file is None:
        print("[calculate_ligand_ec] ⚠️ APBS failed, using Coulomb approximation...")
        # Fallback: Use Coulomb potential from protein charges
        protein_grid = None
        _coulomb_fallback = True
        _coulomb_fallback_reason = "APBS 执行Failed或未安装，已using Coulomb 近似计算蛋白静电势。Results精degrees可能较低。"
    else:
        # Step 4: Load protein potential grid
        print("\n[Step 4] Loading protein potential...")
        protein_grid = DXGrid(dx_file)
    
    # Step 5: Generate ligand surface points
    print("\n[Step 5] Generating ligand surface...")
    sampler = LigandSurfaceSampler(density=surface_density)
    surface_points, surface_normals = sampler.sample_molecule(ligand_mol)

    # Step 6: Calculate potentials at surface points
    print("\n[Step 6] Calculating potentials...")

    # Protein potential (from APBS grid or Coulomb approximation)
    if protein_grid is not None:
        phi_protein = protein_grid.interpolate(surface_points)
    else:
        # Fallback: Use Coulomb potential from protein charges
        print("[calculate_ligand_ec] Using Coulomb approximation for protein potential...")
        charge_calc_protein = GasteigerChargeCalculator()
        phi_protein = charge_calc_protein.calculate_potential_from_pqr(protein_pqr, surface_points)

    # Ligand potential (Gasteiger charges + Coulomb)
    # Use enhanced calculator if σ-hole or lone pairs requested
    charge_calc = GasteigerChargeCalculator(
        use_sigma_holes=use_sigma_holes,
        use_lone_pairs=use_lone_pairs
    )
    phi_ligand = charge_calc.calculate_potential(ligand_mol, surface_points)
    
    # 记录高级功能using情况
    advanced_features = {
        'sigma_holes_enabled': use_sigma_holes,
        'lone_pairs_enabled': use_lone_pairs,
        'bridging_waters_enabled': keep_bridging_waters,
        'n_waters_kept': n_waters_kept,
        'advanced_patches_available': ADVANCED_PATCHES_AVAILABLE,
    }
    
    if use_sigma_holes and ADVANCED_PATCHES_AVAILABLE:
        print("[calculate_ligand_ec] ✨ Using σ-hole virtual points for halogens")
    if use_lone_pairs and ADVANCED_PATCHES_AVAILABLE:
        print("[calculate_ligand_ec] ✨ Using lone pair virtual points")
    
    print(f"[calculate_ligand_ec] φ_protein range: [{phi_protein.min():.3f}, {phi_protein.max():.3f}]")
    print(f"[calculate_ligand_ec] φ_ligand range: [{phi_ligand.min():.3f}, {phi_ligand.max():.3f}]")
    
    # Step 7: Calculate EC
    print("\n[Step 7] Calculating EC...")
    ec_calc = ECCalculator()
    ec_values = ec_calc.calculate_ec_local(phi_protein, phi_ligand)
    ec_score = ec_calc.calculate_ec_score(ec_values)
    ec_stats = ec_calc.calculate_ec_statistics(ec_values)
    
    print(f"\n{'='*60}")
    print("Electrostatic Complementarity Analysis Results")
    print('='*60)
    print(f"EC Score: {ec_score:.4f}")
    print(f"EC Mean: {ec_stats['ec_mean']:.4f}")
    print(f"EC Median: {ec_stats['ec_median']:.4f}")
    print(f"EC Std: {ec_stats['ec_std']:.4f}")
    print(f"Positive EC fraction: {ec_stats['ec_positive_fraction']*100:.1f}%")
    print(f"Negative EC fraction: {ec_stats['ec_negative_fraction']*100:.1f}%")
    print('='*60)
    
    # Step 8: Generate output files
    print("\n[Step 8] Generating output files...")
    
    ec_pdb = os.path.join(output_dir, 'ec_map.pdb')
    ec_csv = os.path.join(output_dir, 'ec_data.csv')
    
    ECMapWriter.write_pseudo_pdb(ec_pdb, surface_points, ec_values)
    ECMapWriter.write_csv(ec_csv, surface_points, ec_values, phi_protein, phi_ligand)
    
    # Step 9: Visualize in PyMOL
    if visualize and PYMOL_AVAILABLE:
        print("\n[Step 9] Visualizing in PyMOL...")
        # Use enhanced smooth surface visualization if available
        if EC_VISUALIZATION_AVAILABLE:
            visualize_ec_smooth_surface(
                obj_name, ligand_resname,
                ec_values=ec_values,
                surface_points=surface_points,
                surface_type='gaussian',
                transparency=0.0
            )
        else:
            # Fallback to basic visualization
            _visualize_ec_map(obj_name, ec_pdb, ligand_resname)
    
    result = {
        'ec_score': ec_score,
        'ec_statistics': ec_stats,
        'surface_points': surface_points,
        'ec_values': ec_values,
        'phi_protein': phi_protein,
        'phi_ligand': phi_ligand,
        'output_files': {
            'ec_map_pdb': ec_pdb,
            'ec_data_csv': ec_csv,
            'protein_pqr': protein_pqr,
            'apbs_dx': dx_file
        },
        'output_dir': output_dir,
        'advanced_features': advanced_features,
        # Coulomb 回退标记：当 APBS 不可用时通知 GUI 层
        'coulomb_fallback': _coulomb_fallback,
        'coulomb_fallback_reason': _coulomb_fallback_reason,
    }
    
    return result


def analyze_ternary_ec(obj_name: str = None, glue_resname: str = None,
                      protein_a_chains: List[str] = None,
                      protein_b_chains: List[str] = None,
                      output_dir: str = None,
                      ph: float = 7.4,
                      surface_density: float = 10.0,
                      visualize: bool = True,
                      pdb_file: str = None,
                      keep_bridging_waters: bool = False,
                      contact_cutoff: float = 10.0) -> Optional[Dict[str, Any]]:
    """
    Analyze Electrostatic Complementarity for a ternary complex (molecular glue).
    
    For molecular glue systems, we calculate EC at multiple interfaces:
    1. EC(A-glue): Protein A potential vs glue surface
    2. EC(B-glue): Protein B potential vs glue surface
    3. EC(A-B): Optional PPI interface analysis
    
    Args:
        obj_name: PyMOL object name
        glue_resname: Molecular glue residue name
        protein_a_chains: Chains for protein A (e.g., E3 ligase)
        protein_b_chains: Chains for protein B (e.g., substrate)
        output_dir: Directory for output files
        ph: pH for protonation state
        surface_density: Surface sampling density
        visualize: Whether to visualize results
        pdb_file: Alternative PDB file input
        keep_bridging_waters: 是否保留结构性桥联水分子（默认 False，向后兼容）
        contact_cutoff: 接触区截断距离 (Å)，仅提取 Glue 周围此距离内的蛋白残基（默认 10.0）
        
    Returns:
        Dictionary containing EC results for each interface
    """
    if not NUMPY_AVAILABLE or not RDKIT_AVAILABLE:
        print("[analyze_ternary_ec] ❌ NumPy and RDKit required")
        return None
    
    if not protein_a_chains or not protein_b_chains:
        print("[analyze_ternary_ec] ❌ Must specify both protein_a_chains and protein_b_chains")
        return None
    

    if output_dir: output_dir = os.path.normpath(output_dir)

    # Create output directory
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='glint_ternary_ec_')
    else:
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"[analyze_ternary_ec] Output directory: {output_dir}")
    print(f"[analyze_ternary_ec] Protein A chains: {protein_a_chains}")
    print(f"[analyze_ternary_ec] Protein B chains: {protein_b_chains}")
    print(f"[analyze_ternary_ec] Glue: {glue_resname}")
    print(f"[analyze_ternary_ec] Contact cutoff: {contact_cutoff} Å")
    
    results = {
        'glue_resname': glue_resname,
        'protein_a_chains': protein_a_chains,
        'protein_b_chains': protein_b_chains,
        'interfaces': {},
        # Coulomb 回退追踪
        'coulomb_fallback': False,
        'coulomb_fallback_reason': '',
    }
    
    # Extract complex
    if pdb_file:
        complex_pdb = pdb_file
    elif PYMOL_AVAILABLE and obj_name:
        complex_pdb = os.path.join(output_dir, 'complex.pdb')
        cmd.save(complex_pdb, obj_name)
    else:
        print("[analyze_ternary_ec] ❌ Need either PyMOL object or PDB file")
        return None
    
    # Extract glue molecule
    glue_sdf = os.path.join(output_dir, 'glue.sdf')
    if PYMOL_AVAILABLE and obj_name:
        cmd.save(glue_sdf, f"{obj_name} and resn {glue_resname}", format='sdf')
    
    # Load glue with RDKit
    supplier = Chem.SDMolSupplier(glue_sdf, removeHs=False)
    glue_mol = supplier[0] if len(supplier) > 0 else None
    
    if glue_mol is None:
        print("[analyze_ternary_ec] ❌ Failed to load glue molecule")
        return None
    
    # 生成 Glue 分子表面点（两个interface共用基础点集）
    sampler = LigandSurfaceSampler(density=surface_density)
    surface_points, surface_normals = sampler.sample_molecule(glue_mol)
    
    # 计算 Glue 电势（共用）
    charge_calc = GasteigerChargeCalculator()
    phi_glue = charge_calc.calculate_potential(glue_mol, surface_points)
    
    # 链Name字符串
    chains_a_str = "_".join(protein_a_chains)
    chains_b_str = "_".join(protein_b_chains)
    
    # ====== 三元复合物关Key改进：Glue 表面朝向分区 ======
    # 提前获取 Protein A / B 坐标，用于面朝向判断
    INTERFACE_CUTOFF = 5.0  # Å — 与蛋白表面的距离截断
    face_a_mask = None  # Glue 表面朝向 Protein A 的点
    face_b_mask = None  # Glue 表面朝向 Protein B 的点
    overlap_face_mask = None  # 桥接区（同时靠近两个蛋白）
    prot_a_coords = None
    prot_b_coords = None
    tree_a = None
    tree_b = None
    
    if PYMOL_AVAILABLE and obj_name and SCIPY_AVAILABLE:
        try:
            # 提取 Protein A 坐标
            coords_a_list = []
            chain_sel_a = ' or '.join([f"chain {c}" for c in protein_a_chains])
            cmd.iterate_state(1, f"{obj_name} and polymer and ({chain_sel_a})",
                             "coords_a_list.append([x,y,z])",
                             space={'coords_a_list': coords_a_list})
            
            # 提取 Protein B 坐标
            coords_b_list = []
            chain_sel_b = ' or '.join([f"chain {c}" for c in protein_b_chains])
            cmd.iterate_state(1, f"{obj_name} and polymer and ({chain_sel_b})",
                             "coords_b_list.append([x,y,z])",
                             space={'coords_b_list': coords_b_list})
            
            if coords_a_list and coords_b_list:
                prot_a_coords = np.array(coords_a_list)
                prot_b_coords = np.array(coords_b_list)
                tree_a = cKDTree(prot_a_coords)
                tree_b = cKDTree(prot_b_coords)
                
                # 计算每个 Glue 表面点到两个蛋白的最近距离
                dists_to_a, _ = tree_a.query(surface_points)
                dists_to_b, _ = tree_b.query(surface_points)
                
                # 面朝向分区
                face_a_mask = (dists_to_a <= INTERFACE_CUTOFF)
                face_b_mask = (dists_to_b <= INTERFACE_CUTOFF)
                overlap_face_mask = face_a_mask & face_b_mask  # 桥接区
                a_only_mask = face_a_mask & ~face_b_mask
                b_only_mask = face_b_mask & ~face_a_mask
                
                n_face_a = int(np.sum(face_a_mask))
                n_face_b = int(np.sum(face_b_mask))
                n_overlap = int(np.sum(overlap_face_mask))
                n_neither = int(np.sum(~face_a_mask & ~face_b_mask))
                
                print(f"[analyze_ternary_ec] 🔬 Glue 表面朝向分区:")
                print(f"  朝 Protein A: {n_face_a} 点 (含桥接区)")
                print(f"  朝 Protein B: {n_face_b} 点 (含桥接区)")
                print(f"  桥接区 (重叠): {n_overlap} 点")
                print(f"  非interface区: {n_neither} 点")
            else:
                print("[analyze_ternary_ec] ⚠️ 无法获取蛋白坐标，降级为全表面分析")
        except Exception as e:
            print(f"[analyze_ternary_ec] ⚠️ 面朝向分区Failed: {e}")
            print("[analyze_ternary_ec] 降级为全表面分析")
    
    # 桥联水筛选（三元复合物共用，如Enable）
    ternary_waters_to_keep = set()
    if keep_bridging_waters and ADVANCED_PATCHES_AVAILABLE:
        print("[analyze_ternary_ec] ✨ 筛选结构性桥联水...")
        try:
            water_filter = BridgingWaterFilter()
            # SavePackage含水的复合物临时File
            complex_with_water_pdb = os.path.join(output_dir, 'complex_with_water.pdb')
            if PYMOL_AVAILABLE and obj_name:
                cmd.save(complex_with_water_pdb, f"{obj_name} and (polymer or resn HOH)")
            
            # 筛选桥联水（ReturnValue为 (chain, resnum) 元组列表）
            water_result = water_filter.filter_pdb(complex_with_water_pdb, glue_resname)
            for chain, resnum in water_result.get('bridging_waters', []):
                ternary_waters_to_keep.add((chain, resnum))
            for chain, resnum in water_result.get('protein_bound_waters', []):
                ternary_waters_to_keep.add((chain, resnum))
            
            if ternary_waters_to_keep:
                print(f"[analyze_ternary_ec] 保留 {len(ternary_waters_to_keep)} 个结构性/桥联水分子")
            else:
                print("[analyze_ternary_ec] 未发现符合条件的桥联水")
        except Exception as e:
            print(f"[analyze_ternary_ec] ⚠️ 桥联水筛选Failed: {e}")
            print("[analyze_ternary_ec] Continue分析（不含水）")
    
    # 如果保留了桥联水，PDB2PQR 不Delete水
    should_remove_water = not (keep_bridging_waters and ADVANCED_PATCHES_AVAILABLE and len(ternary_waters_to_keep) > 0)
    
    # Initialize共享变量（避免变量作用域问题 —— 当某个 Interface 的 APBS Failed时仍可用）
    grid_a = None
    grid_b = None
    ec_calc = ECCalculator()  # 修复: 提前Initialize，避免 Interface A Failed时 Interface B 找不到 ec_calc
    dx_file_a = None  # APBS 输出FilePath
    dx_file_b = None

    # Interface 1: Protein A - Glue
    print("\n" + "="*60)
    print(f"Interface 1: Protein A ({chains_a_str}) - Glue")
    print("="*60)
    
    protein_a_dir = os.path.join(output_dir, f'protein_{chains_a_str}')
    os.makedirs(protein_a_dir, exist_ok=True)
    
    # Extract protein A（接触区残基 + 桥联水支持）
    protein_a_pdb = os.path.join(protein_a_dir, 'protein_a.pdb')
    if PYMOL_AVAILABLE and obj_name:
        chain_sel = ' or '.join([f"chain {c}" for c in protein_a_chains])
        glue_sel = f"{obj_name} and resn {glue_resname}"
        # 接触区残基提取：仅提取 Glue 周围 contact_cutoff 距离内的蛋白残基
        protein_a_sel = f"byres (({obj_name} and polymer and ({chain_sel})) within {contact_cutoff} of ({glue_sel}))"
        
        # Print接触区Information
        try:
            n_contact_atoms_a = cmd.count_atoms(protein_a_sel)
            contact_residues_a = set()
            cmd.iterate(protein_a_sel, "contact_residues_a.add((chain, resi))",
                       space={'contact_residues_a': contact_residues_a})
            print(f"[Interface A] 📍 接触区提取 (cutoff={contact_cutoff} Å): {len(contact_residues_a)} 残基, {n_contact_atoms_a} 原子")
        except Exception:
            pass
        
        # 如果有桥联水，ExtensionSelect器以Package含关Key水分子
        if ternary_waters_to_keep:
            water_resi_sel = ' or '.join(
                [f"(chain {c} and resi {r})" for c, r in ternary_waters_to_keep]
            )
            water_sel = f"({obj_name} and resn HOH and ({water_resi_sel}))"
            protein_a_sel = f"({protein_a_sel}) or ({water_sel})"
        
        cmd.save(protein_a_pdb, protein_a_sel)
    
    # Run PDB2PQR and APBS for protein A
    protein_a_pqr = os.path.join(protein_a_dir, 'protein_a.pqr')
    pdb2pqr = PDB2PQRRunner()
    pdb2pqr.run(protein_a_pdb, protein_a_pqr, ph=ph, remove_water=should_remove_water)
    
    apbs = APBSRunner()
    apbs_prefix_a = os.path.join(protein_a_dir, 'protein_a_pot')
    # using Glue 几何中心作为 APBS 网格中心，确保网格覆盖 Glue 表面区域
    glue_center = _get_molecule_center(glue_mol)
    # 修复: 根据 contact_cutoff 动态调整网格尺寸，确保覆盖整个接触区
    fine_size = max(30, int(contact_cutoff * 2 + 10))
    ternary_fine_len = (fine_size, fine_size, fine_size)
    apbs_input_a = apbs.generate_input(
        protein_a_pqr, apbs_prefix_a,
        grid_center=glue_center,
        fine_len=ternary_fine_len
    )
    dx_file_a = apbs.run(apbs_input_a, protein_a_dir)
    
    # APBS Success时using DX 网格插Value，Failed时回退到 Coulomb 近似
    if dx_file_a:
        grid_a = DXGrid(dx_file_a)
        phi_protein_a = grid_a.interpolate(surface_points)
    else:
        # Coulomb 近似回退：当 APBS Failed（如 PQR 解析Error）时，using Gasteiger 电荷直接计算静电势
        print("[analyze_ternary_ec] ⚠️ APBS failed for Protein A, using Coulomb approximation...")
        grid_a = None
        charge_calc_a = GasteigerChargeCalculator()
        phi_protein_a = charge_calc_a.calculate_potential_from_pqr(protein_a_pqr, surface_points)
        # 标记 Coulomb 回退
        results['coulomb_fallback'] = True
        results['coulomb_fallback_reason'] = "Protein A APBS 计算Failed，已using Coulomb 近似。"
    
    # using面朝向 mask（优先）或距离截断来筛选有效interface点（无论 APBS 是否Success都执行）
    phi_protein_a_masked = phi_protein_a.copy()
    if face_a_mask is not None:
        # 三元复合物改进：只保留朝向 Protein A 的 Glue 表面点
        mask_a_invalid = ~face_a_mask
        phi_protein_a_masked[mask_a_invalid] = 0.0
        print(f"[Interface A-Glue] using面朝向 mask: {int(np.sum(face_a_mask))}/{len(surface_points)} 有效点")
    elif SCIPY_AVAILABLE:
        # 降级：using PDB 坐标 + 距离截断
        prot_a_mol = Chem.MolFromPDBFile(protein_a_pdb, removeHs=False)
        if prot_a_mol:
            prot_coords = prot_a_mol.GetConformer().GetPositions()
            tree = cKDTree(prot_coords)
            dists, _ = tree.query(surface_points)
            mask = dists > INTERFACE_CUTOFF
            phi_protein_a_masked[mask] = 0.0

    ec_calc = ECCalculator()
    ec_a_glue = ec_calc.calculate_ec_local(phi_protein_a_masked, phi_glue)
    ec_score_a = ec_calc.calculate_ec_score(ec_a_glue)
    ec_stats_a = ec_calc.calculate_ec_statistics(ec_a_glue, exclude_zeros=True)
    
    results['interfaces']['A_glue'] = {
        'ec_score': ec_score_a,
        'ec_statistics': ec_stats_a,
        'ec_values': ec_a_glue,
        'phi_protein': phi_protein_a,
        'chains': protein_a_chains,
        'n_face_points': int(np.sum(face_a_mask)) if face_a_mask is not None else len(surface_points)
    }
    
    print(f"EC(A-Glue) Score: {ec_score_a:.4f}")
    
    # Write EC map with chain-based name
    ec_pdb_a = os.path.join(protein_a_dir, f'ec_map_{chains_a_str}_glue.pdb')
    ECMapWriter.write_pseudo_pdb(ec_pdb_a, surface_points, ec_a_glue)
    
    # Interface 2: Protein B - Glue
    print("\n" + "="*60)
    print(f"Interface 2: Protein B ({chains_b_str}) - Glue")
    print("="*60)
    
    protein_b_dir = os.path.join(output_dir, f'protein_{chains_b_str}')
    os.makedirs(protein_b_dir, exist_ok=True)
    
    # Extract protein B（接触区残基 + 桥联水支持）
    protein_b_pdb = os.path.join(protein_b_dir, 'protein_b.pdb')
    if PYMOL_AVAILABLE and obj_name:
        chain_sel = ' or '.join([f"chain {c}" for c in protein_b_chains])
        glue_sel = f"{obj_name} and resn {glue_resname}"
        # 接触区残基提取：仅提取 Glue 周围 contact_cutoff 距离内的蛋白残基
        protein_b_sel = f"byres (({obj_name} and polymer and ({chain_sel})) within {contact_cutoff} of ({glue_sel}))"
        
        # Print接触区Information
        try:
            n_contact_atoms_b = cmd.count_atoms(protein_b_sel)
            contact_residues_b = set()
            cmd.iterate(protein_b_sel, "contact_residues_b.add((chain, resi))",
                       space={'contact_residues_b': contact_residues_b})
            print(f"[Interface B] 📍 接触区提取 (cutoff={contact_cutoff} Å): {len(contact_residues_b)} 残基, {n_contact_atoms_b} 原子")
        except Exception:
            pass
        
        # 如果有桥联水，ExtensionSelect器以Package含关Key水分子
        if ternary_waters_to_keep:
            water_resi_sel = ' or '.join(
                [f"(chain {c} and resi {r})" for c, r in ternary_waters_to_keep]
            )
            water_sel = f"({obj_name} and resn HOH and ({water_resi_sel}))"
            protein_b_sel = f"({protein_b_sel}) or ({water_sel})"
        
        cmd.save(protein_b_pdb, protein_b_sel)
    
    # Run PDB2PQR and APBS for protein B
    protein_b_pqr = os.path.join(protein_b_dir, 'protein_b.pqr')
    pdb2pqr.run(protein_b_pdb, protein_b_pqr, ph=ph, remove_water=should_remove_water)
    
    apbs_prefix_b = os.path.join(protein_b_dir, 'protein_b_pot')
    # 同样using Glue 中心和动态网格
    apbs_input_b = apbs.generate_input(
        protein_b_pqr, apbs_prefix_b,
        grid_center=glue_center,
        fine_len=ternary_fine_len
    )
    dx_file_b = apbs.run(apbs_input_b, protein_b_dir)
    
    # APBS Success时using DX 网格插Value，Failed时回退到 Coulomb 近似
    if dx_file_b:
        grid_b = DXGrid(dx_file_b)
        phi_protein_b = grid_b.interpolate(surface_points)
    else:
        # Coulomb 近似回退：当 APBS Failed（如 PQR 解析Error）时，using Gasteiger 电荷直接计算静电势
        print("[analyze_ternary_ec] ⚠️ APBS failed for Protein B, using Coulomb approximation...")
        grid_b = None
        charge_calc_b = GasteigerChargeCalculator()
        phi_protein_b = charge_calc_b.calculate_potential_from_pqr(protein_b_pqr, surface_points)
        # 标记 Coulomb 回退
        results['coulomb_fallback'] = True
        reason_b = "Protein B APBS 计算Failed，已using Coulomb 近似。"
        existing_reason = results.get('coulomb_fallback_reason', '')
        results['coulomb_fallback_reason'] = (existing_reason + " " + reason_b).strip()
    
    # using面朝向 mask（优先）或距离截断来筛选有效interface点（无论 APBS 是否Success都执行）
    phi_protein_b_masked = phi_protein_b.copy()
    if face_b_mask is not None:
        # 三元复合物改进：只保留朝向 Protein B 的 Glue 表面点
        mask_b_invalid = ~face_b_mask
        phi_protein_b_masked[mask_b_invalid] = 0.0
        print(f"[Interface B-Glue] using面朝向 mask: {int(np.sum(face_b_mask))}/{len(surface_points)} 有效点")
    elif SCIPY_AVAILABLE:
        # 降级：using PDB 坐标 + 距离截断
        prot_b_mol = Chem.MolFromPDBFile(protein_b_pdb, removeHs=False)
        if prot_b_mol:
            prot_coords = prot_b_mol.GetConformer().GetPositions()
            tree = cKDTree(prot_coords)
            dists, _ = tree.query(surface_points)
            mask = dists > INTERFACE_CUTOFF
            phi_protein_b_masked[mask] = 0.0

    ec_b_glue = ec_calc.calculate_ec_local(phi_protein_b_masked, phi_glue)
    ec_score_b = ec_calc.calculate_ec_score(ec_b_glue)
    ec_stats_b = ec_calc.calculate_ec_statistics(ec_b_glue, exclude_zeros=True)
    
    results['interfaces']['B_glue'] = {
        'ec_score': ec_score_b,
        'ec_statistics': ec_stats_b,
        'ec_values': ec_b_glue,
        'phi_protein': phi_protein_b,
        'chains': protein_b_chains,
        'n_face_points': int(np.sum(face_b_mask)) if face_b_mask is not None else len(surface_points)
    }
    
    print(f"EC(B-Glue) Score: {ec_score_b:.4f}")
    
    # Write EC map with chain-based name
    ec_pdb_b = os.path.join(protein_b_dir, f'ec_map_{chains_b_str}_glue.pdb')
    ECMapWriter.write_pseudo_pdb(ec_pdb_b, surface_points, ec_b_glue)
    
    # Combined analysis
    print("\n" + "="*60)
    print("Combined Ternary EC Analysis")
    print("="*60)
    
    if 'A_glue' in results['interfaces'] and 'B_glue' in results['interfaces']:
        ec_a_values = results['interfaces']['A_glue']['ec_values']
        ec_b_values = results['interfaces']['B_glue']['ec_values']
        
        # Identify interface regions
        interface_a = ec_a_values != 0  # Points near Protein A
        interface_b = ec_b_values != 0  # Points near Protein B
        overlap = interface_a & interface_b  # Points near both proteins
        only_a = interface_a & ~interface_b
        only_b = ~interface_a & interface_b
        
        # Calculate EC for different regions
        ec_calc = ECCalculator()
        
        # Overall scores (using non-zero points)
        ec_a = ec_calc.calculate_ec_score(ec_a_values)
        ec_b = ec_calc.calculate_ec_score(ec_b_values)
        
        # Overlap region scores (most important for molecular glue!)
        if np.sum(overlap) > 0:
            ec_a_overlap = float(np.mean(ec_a_values[overlap]))
            ec_b_overlap = float(np.mean(ec_b_values[overlap]))
            ec_overlap_combined = (ec_a_overlap + ec_b_overlap) / 2
            # Positive EC fraction (more intuitive metric)
            pos_frac_a = float(np.mean(ec_a_values[overlap] > 0))
            pos_frac_b = float(np.mean(ec_b_values[overlap] > 0))
            pos_frac_combined = (pos_frac_a + pos_frac_b) / 2
        else:
            ec_a_overlap = 0.0
            ec_b_overlap = 0.0
            ec_overlap_combined = 0.0
            pos_frac_a = 0.0
            pos_frac_b = 0.0
            pos_frac_combined = 0.0
        
        # Non-overlap region scores
        ec_a_only = float(np.mean(ec_a_values[only_a])) if np.sum(only_a) > 0 else 0.0
        ec_b_only = float(np.mean(ec_b_values[only_b])) if np.sum(only_b) > 0 else 0.0
        
        # Combined score (average of overall)
        ec_combined = (ec_a + ec_b) / 2
        
        # Asymmetry (difference between interfaces)
        ec_asymmetry = abs(ec_a - ec_b)
        
        results['combined'] = {
            'ec_combined_score': ec_combined,
            'ec_asymmetry': ec_asymmetry,
            'ec_a_glue': ec_a,
            'ec_b_glue': ec_b,
            'ec_overlap': {
                'ec_a_overlap': ec_a_overlap,
                'ec_b_overlap': ec_b_overlap,
                'ec_combined': ec_overlap_combined,
                'pos_frac_a': pos_frac_a,
                'pos_frac_b': pos_frac_b,
                'pos_frac_combined': pos_frac_combined,
                'n_points': int(np.sum(overlap))
            },
            'ec_non_overlap': {
                'ec_a_only': ec_a_only,
                'ec_b_only': ec_b_only,
                'n_points_a': int(np.sum(only_a)),
                'n_points_b': int(np.sum(only_b))
            }
        }
        
        print(f"\nOverall EC (all interface points):")
        print(f"  EC(A-Glue): {ec_a:.4f}")
        print(f"  EC(B-Glue): {ec_b:.4f}")
        print(f"  Combined EC: {ec_combined:.4f}")
        
        print(f"\n⭐ Overlap Region (glue bridging zone - {np.sum(overlap)} points):")
        print(f"  Positive EC fraction (A-Glue): {pos_frac_a*100:.1f}%")
        print(f"  Positive EC fraction (B-Glue): {pos_frac_b*100:.1f}%")
        print(f"  Combined positive fraction: {pos_frac_combined*100:.1f}%")
        print(f"  (Mean EC values: A={ec_a_overlap:.4f}, B={ec_b_overlap:.4f})")
        
        print(f"\nNon-overlap Region EC:")
        print(f"  EC(A-Glue) A-only region ({np.sum(only_a)} points): {ec_a_only:.4f}")
        print(f"  EC(B-Glue) B-only region ({np.sum(only_b)} points): {ec_b_only:.4f}")
        
        print(f"\nAsymmetry: {ec_asymmetry:.4f}")
        
        # Interpretation based on positive EC fraction (more intuitive!)
        if pos_frac_combined > 0.6:
            print("\n✨ Good electrostatic complementarity in bridging zone (>60% positive)")
        elif pos_frac_combined > 0.5:
            print("\n✓ Moderate electrostatic complementarity in bridging zone (50-60% positive)")
        else:
            print("\n⚠️ Poor electrostatic complementarity in bridging zone (<50% positive)")

    # Interface 3: Protein A - Protein B (PPI interface with glue)
    # 注意：PPI 分析using APBS 网格或 Coulomb 近似来计算 Protein B 在interface的电势
    print("\n" + "="*60)
    print("Interface 3: Protein A - Protein B (PPI with Glue)")
    print("="*60)
    
    ppi_dir = os.path.join(output_dir, 'ppi_interface')
    os.makedirs(ppi_dir, exist_ok=True)
    
    # Get PPI interface residues
    if PYMOL_AVAILABLE and obj_name:
        # Find interface residues between A and B
        chain_a_sel = ' or '.join([f"chain {c}" for c in protein_a_chains])
        chain_b_sel = ' or '.join([f"chain {c}" for c in protein_b_chains])
        
        # Get interface atoms (within 5Å of each other) — 修复: Add polymer Filter，排除配体/水/离子
        interface_a_sel = f"({obj_name} and polymer and ({chain_a_sel})) within 5.0 of ({obj_name} and polymer and ({chain_b_sel}))"
        interface_b_sel = f"({obj_name} and polymer and ({chain_b_sel})) within 5.0 of ({obj_name} and polymer and ({chain_a_sel}))"
        
        # Sample PPI interface surface
        try:
            # Get interface coordinates
            interface_a_coords = []
            interface_a_elements = []
            cmd.iterate_state(1, interface_a_sel,
                             "interface_a_coords.append([x,y,z]); interface_a_elements.append(elem)",
                             space={'interface_a_coords': interface_a_coords,
                                   'interface_a_elements': interface_a_elements})
            
            interface_b_coords = []
            interface_b_elements = []
            cmd.iterate_state(1, interface_b_sel,
                             "interface_b_coords.append([x,y,z]); interface_b_elements.append(elem)",
                             space={'interface_b_coords': interface_b_coords,
                                   'interface_b_elements': interface_b_elements})
            
            if interface_a_coords and interface_b_coords:
                # Sample interface surface
                ppi_sampler = LigandSurfaceSampler(density=surface_density / 2)  # Lower density for PPI
                
                # Sample from interface A (facing B)
                ppi_points_a, _ = ppi_sampler.sample_from_coords(
                    np.array(interface_a_coords), interface_a_elements)
                
                # 计算 Protein B 在 PPI interface表面点的静电势
                # 优先using APBS 网格插Value，APBS Failed时回退到 Coulomb 近似
                if grid_b is not None:
                    phi_b_at_interface = grid_b.interpolate(ppi_points_a)
                else:
                    # Coulomb 近似回退：using Protein B 的 PQR 电荷计算 PPI interface电势
                    print("[analyze_ternary_ec] ⚠️ PPI: Using Coulomb approximation for Protein B potential...")
                    charge_calc_ppi_b = GasteigerChargeCalculator()
                    phi_b_at_interface = charge_calc_ppi_b.calculate_potential_from_pqr(protein_b_pqr, ppi_points_a)
                
                # Calculate "potential" from interface A atoms (simplified)
                # 修复: 预计算 charges 和坐标数组，避免循环内重复Create
                charges_a = np.array([_get_element_charge(e) for e in interface_a_elements])
                interface_a_np = np.array(interface_a_coords)
                phi_a_at_interface = np.zeros(len(ppi_points_a))
                for i, point in enumerate(ppi_points_a):
                    distances = np.linalg.norm(interface_a_np - point, axis=1)
                    distances = np.maximum(distances, 0.5)
                    phi_a_at_interface[i] = np.sum(charges_a / distances) * 332.0637 / 0.593
                
                ec_ppi = ec_calc.calculate_ec_local(phi_b_at_interface, phi_a_at_interface)
                ec_score_ppi = ec_calc.calculate_ec_score(ec_ppi)
                ec_stats_ppi = ec_calc.calculate_ec_statistics(ec_ppi)
                
                results['interfaces']['A_B_ppi'] = {
                    'ec_score': ec_score_ppi,
                    'ec_statistics': ec_stats_ppi,
                    'ec_values': ec_ppi,
                    'surface_points': ppi_points_a
                }
                
                print(f"EC(A-B PPI) Score: {ec_score_ppi:.4f}")
                
                # Write PPI EC map
                ec_pdb_ppi = os.path.join(ppi_dir, 'ec_map_ppi.pdb')
                ECMapWriter.write_pseudo_pdb(ec_pdb_ppi, ppi_points_a, ec_ppi)
        except Exception as e:
            # PPI interface分析异常处理（修复：补全缺失的 except 子句）
            print(f"[analyze_ternary_ec] PPI interface analysis failed: {e}")
    
    results['surface_points'] = surface_points
    results['phi_glue'] = phi_glue
    results['output_dir'] = output_dir
    
    # Visualize
    if visualize and PYMOL_AVAILABLE and obj_name:
        _visualize_ternary_ec(obj_name, results)
    
    return results


# ========== Helper Functions ==========

def _extract_protein_ligand_from_pdb(complex_pdb: str, protein_pdb: str, 
                                     ligand_sdf: str, ligand_resname: str = None):
    """Extract protein and ligand from a PDB file without PyMOL."""
    # Standard amino acids
    standard_aa = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    
    # Solvent/ions to exclude
    exclude = {'HOH', 'WAT', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE'}
    
    protein_lines = []
    ligand_lines = []
    
    with open(complex_pdb, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                resname = line[17:20].strip()
                
                if resname in standard_aa:
                    protein_lines.append(line)
                elif resname not in exclude:
                    if ligand_resname is None or resname == ligand_resname:
                        ligand_lines.append(line)
    
    # Write protein PDB
    with open(protein_pdb, 'w') as f:
        f.writelines(protein_lines)
        f.write('END\n')
    
    # Write ligand (as PDB, then convert to SDF with RDKit)
    ligand_pdb = ligand_sdf.replace('.sdf', '.pdb')
    with open(ligand_pdb, 'w') as f:
        f.writelines(ligand_lines)
        f.write('END\n')
    
    # Convert to SDF
    if RDKIT_AVAILABLE:
        mol = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
        if mol:
            writer = Chem.SDWriter(ligand_sdf)
            writer.write(mol)
            writer.close()


def _simple_pdb_to_pqr(pdb_file: str, pqr_file: str) -> Optional[str]:
    """Simple PDB to PQR conversion with basic charges (fallback).
    
    using空格分隔的 PQR 格式输出，避免固定列拼接导致 APBS 解析Failed。
    """
    # 基于残基Type的简单电荷分配
    charges = {
        'ARG': {'NH1': 0.5, 'NH2': 0.5, 'NE': 0.0},
        'LYS': {'NZ': 1.0},
        'ASP': {'OD1': -0.5, 'OD2': -0.5},
        'GLU': {'OE1': -0.5, 'OE2': -0.5},
        'HIS': {'ND1': 0.25, 'NE2': 0.25}
    }
    
    # 默认原子半径
    radii = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2}
    
    try:
        with open(pdb_file, 'r') as f_in, open(pqr_file, 'w') as f_out:
            for line in f_in:
                if line.startswith(('ATOM', 'HETATM')):
                    # 从 PDB 固定列位置提取各字段
                    record_type = line[0:6].strip()
                    serial = line[6:11].strip()
                    atom_name = line[12:16].strip()
                    resname = line[17:20].strip()
                    chain = line[21:22].strip()
                    resseq = line[22:26].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                    
                    # 获取电荷
                    charge = 0.0
                    if resname in charges and atom_name in charges[resname]:
                        charge = charges[resname][atom_name]
                    
                    # 获取半径
                    radius = radii.get(element.upper(), 1.7)
                    
                    # 格式化原子Name（保持 PDB 对齐规则：1字符元素右移一位）
                    if len(atom_name) < 4:
                        atom_name_fmt = f" {atom_name:<3s}"
                    else:
                        atom_name_fmt = f"{atom_name:<4s}"
                    
                    # using空格分隔的 PQR 格式写入，确保 APBS 可正确解析
                    pqr_line = (
                        f"{record_type:<6s}{int(serial):>5d} {atom_name_fmt} "
                        f"{resname:<3s} {chain:1s}{resseq:>4s}    "
                        f"{x:8.3f} {y:8.3f} {z:8.3f} "
                        f"{charge:7.4f} {radius:6.4f}\n"
                    )
                    f_out.write(pqr_line)
                elif line.startswith(('END', 'TER', 'REMARK')):
                    # 保留结构标记行
                    f_out.write(line)
        
        return pqr_file
    except Exception as e:
        print(f"[_simple_pdb_to_pqr] Error: {e}")
        return None


def _fix_pqr_format(pqr_file: str) -> bool:
    """修复 PQR File格式，将固定列格式转为空格分隔格式。
    
    PDB2PQR（Python API 或命令行）生成的 PQR File可能using PDB 式固定列格式，
    当坐标Value较大（如 y=-104.24）时，相邻字段会连在一起（字段拼接），
    导致 APBS 无法解析并报错：
      "Valist_readPQR: Error parsing atom! ...no concatenated fields."
    
    本Function依次尝试三种解析策略：
      1. PDB 固定列解析（标准 PDB 列位置）
      2. 空格分隔解析（适用于已经是空格分隔的行）
      3. 正则表达式兜底解析（从行中提取所有浮点数）
    
    修复Completed后会验证第一条 ATOM 行的字段数是否符合 APBS 要求。
    
    Args:
        pqr_file: PQR FilePath
        
    Returns:
        True 表示修复Success，False 表示修复Failed（原File不受影响）
    """
    import re
    
    try:
        with open(pqr_file, 'r') as f:
            lines = f.readlines()
        
        # 调试日志：Print修复前第一条 ATOM/HETATM 行
        first_atom_line = next((l for l in lines if l.startswith(('ATOM', 'HETATM'))), None)
        if first_atom_line:
            print(f"[_fix_pqr_format] 修复前第一行: {first_atom_line.rstrip()}")
        
        fixed_lines = []
        first_fixed_logged = False
        parse_fail_count = 0
        total_atom_count = 0
        
        for line in lines:
            # 仅处理 ATOM/HETATM 记录行
            if not line.startswith(('ATOM', 'HETATM')):
                fixed_lines.append(line)
                continue
            
            total_atom_count += 1
            
            # 依次尝试三种解析Method
            # Method一：按 PDB 固定列位置解析
            parsed = _parse_pqr_fixed_columns(line)
            
            # Method二：按空格分割解析
            if parsed is None:
                parsed = _parse_pqr_whitespace(line)
            
            # Method三：正则表达式兜底（从拼接的行中提取浮点数）
            if parsed is None:
                parsed = _parse_pqr_regex(line)
            
            if parsed is None:
                # 三种方式全部Failed，记录并保留原行
                parse_fail_count += 1
                print(f"[_fix_pqr_format] ⚠️ 无法解析行 (三种Method均Failed)，保留原样: {line.rstrip()}")
                fixed_lines.append(line)
                continue
            
            record_type, serial, atom_name, resname, chain, resseq, x, y, z, charge, radius = parsed
            
            # APBS PQR 格式：纯空格分隔，每个字段间至少一个空格
            # 格式: RECORD serial atom_name resname chain resseq x y z charge radius
            # 注意：APBS 对 chain 字段不敏感，但空格分隔必须正确
            chain_str = chain.strip() if chain.strip() else "A"
            fixed_line = (
                f"{record_type:<6s} {serial:>5d} {atom_name:<4s} {resname:<3s} "
                f"{chain_str:>1s} {resseq:>4s} "
                f"{x:>8.3f} {y:>8.3f} {z:>8.3f} "
                f"{charge:>7.4f} {radius:>6.4f}\n"
            )
            fixed_lines.append(fixed_line)
            
            # 调试日志：Print修复后第一条行，方便对比
            if not first_fixed_logged:
                print(f"[_fix_pqr_format] 修复后第一行: {fixed_line.rstrip()}")
                first_fixed_logged = True
        
        # === 验证步骤：检查修复后第一条 ATOM 行是否可被正确解析 ===
        first_fixed_atom = next(
            (l for l in fixed_lines if l.startswith(('ATOM', 'HETATM'))), None
        )
        if first_fixed_atom:
            parts = first_fixed_atom.split()
            # APBS 期望空格分隔后至少有 10-11 个字段
            if len(parts) < 10:
                print(f"[_fix_pqr_format] ❌ 验证Failed: 修复后第一行仅有 {len(parts)} 个字段 "
                      f"(期望 ≥10): {first_fixed_atom.rstrip()}")
                # 不写回File，ReturnFailed
                return False
            
            # 验证坐标和电荷/半径字段是否为有效浮点数
            try:
                # 尝试从末尾解析 5 个浮点数 (x, y, z, charge, radius)
                for idx in range(-5, 0):
                    float(parts[idx])
            except (ValueError, IndexError):
                print(f"[_fix_pqr_format] ❌ 验证Failed: 修复后行的数Value字段无法解析: "
                      f"{first_fixed_atom.rstrip()}")
                return False
            
            print(f"[_fix_pqr_format] ✅ 验证通过: 修复后第一行有 {len(parts)} 个字段，数Value字段正常")
        
        if parse_fail_count > 0:
            print(f"[_fix_pqr_format] ⚠️ 有 {parse_fail_count}/{total_atom_count} 行解析Failed")
        
        # 写回File
        with open(pqr_file, 'w') as f:
            f.writelines(fixed_lines)
        
        print(f"[_fix_pqr_format] ✅ PQR 格式已修复: {pqr_file} ({total_atom_count} 个原子)")
        return True
        
    except Exception as e:
        print(f"[_fix_pqr_format] ❌ 修复Failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def _parse_pqr_fixed_columns(line: str):
    """按 PDB 固定列位置解析 PQR 行。
    
    PDB/PQR 固定列定义（标准 PDB 格式）：
      record_type: [0:6], serial: [6:11], atom_name: [12:16],
      resname: [17:20], chain: [21:22], resseq: [22:26],
      x: [30:38], y: [38:46], z: [46:54],
      charge: [54:62], radius: [62:70] (PQR Extension字段)
    
    增强处理：
      - 当坐标Value较大（如 -104.240）导致字段拼接时，尝试智能拆分
      - 支持负坐标导致的列宽溢出（如 -1004.240 溢出 8 字符列宽）
      - 支持 PDB2PQR Python API 输出的微小列shift
    
    Returns:
        解析SuccessReturn元组 (record_type, serial, atom_name, resname, chain, resseq, x, y, z, charge, radius)
        解析FailedReturn None
    """
    import re
    
    try:
        # 行长degrees不够则无法按固定列解析
        if len(line) < 54:
            return None
        
        record_type = line[0:6].strip()
        if record_type not in ('ATOM', 'HETATM'):
            return None
        
        serial_str = line[6:11].strip()
        if not serial_str:
            return None
        serial = int(serial_str)
        
        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        # 有时 resname 会延伸到第 21 列（如 4 字符残基名），Package容处理
        if not resname:
            resname = line[16:21].strip()
        
        chain = line[21:22].strip() if len(line) > 21 and line[21:22].strip() else ' '
        resseq = line[22:26].strip()
        
        # ── 坐标解析（增强版：处理负坐标导致的列溢出和字段拼接） ──
        # 先尝试标准 PDB 列位置 [30:38], [38:46], [46:54]
        coords_parsed = False
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords_parsed = True
        except ValueError:
            pass
        
        if not coords_parsed:
            # 标准列位置解析Failed，using增强正则从坐标区域提取
            # 正则匹配：限制小数位数为 1-4 位，防止贪婪匹配吃掉相邻字段
            # 例如 "27.340-104.240" 正确拆分为 ["27.340", "-104.240"]
            enhanced_float = re.compile(r'[+-]?(?:\d+\.\d{1,4}|\.\d{1,4}|\d+(?:[eE][+-]?\d+))')
            
            # 从第 26 列StartSearch（覆盖可能的列shift），范围扩大到第 60 列
            coord_region = line[26:60] if len(line) >= 60 else line[26:]
            coord_matches = enhanced_float.findall(coord_region)
            
            # Filter掉可能混入的残基Number（通常是纯整数且位于区域开头）
            # 坐标Value通常含小数点，残基Number通常不含
            if len(coord_matches) > 3:
                # 优先Select含小数点的匹配项
                decimal_matches = [m for m in coord_matches if '.' in m]
                if len(decimal_matches) >= 3:
                    coord_matches = decimal_matches[:3]
                else:
                    # 取最后 3 个（Skip前面可能的残基Number）
                    # 但需确保取出的是坐标而非 charge/radius
                    coord_matches = coord_matches[:3]
            
            if len(coord_matches) < 3:
                return None
            
            x = float(coord_matches[0])
            y = float(coord_matches[1])
            z = float(coord_matches[2])
        
        # ── 电荷和半径解析（增强版：处理列shift） ──
        charge = 0.0
        radius = 0.0
        
        # 尝试多个起始位置提取 charge/radius（应对坐标列溢出导致的shift）
        # 优先从标准位置 [54:] Start，如果Failed则向后shift 1-2 列
        tail_candidates = []
        for start_col in [54, 55, 56, 52, 53]:
            if len(line) > start_col:
                tail_candidates.append(line[start_col:].strip())
        
        # 匹配含小数点的浮点数（charge/radius 一定有小数点，限制 1-4 位小数防止贪婪溢出）
        tail_float = re.compile(r'[+-]?(?:\d+\.\d{1,4}|\.\d{1,4})')
        
        for tail in tail_candidates:
            if not tail:
                continue
            tail_matches = tail_float.findall(tail)
            if len(tail_matches) >= 2:
                c_val = float(tail_matches[0])
                r_val = float(tail_matches[1])
                # 合理性检查：charge 通常在 [-3, 3]，radius 通常在 [0, 5]
                if abs(c_val) <= 10 and 0 <= r_val <= 10:
                    charge = c_val
                    radius = r_val
                    break
            elif len(tail_matches) == 1:
                charge = float(tail_matches[0])
                break
        
        # 基本合理性检查：坐标不应超过 ±9999，半径应为正数或零
        if abs(x) > 9999 or abs(y) > 9999 or abs(z) > 9999:
            return None
        if radius < 0:
            return None
        
        return (record_type, serial, atom_name, resname, chain, resseq, x, y, z, charge, radius)
        
    except (ValueError, IndexError):
        return None


def _parse_pqr_whitespace(line: str):
    """按空格分割解析 PQR 行。
    
    PQR 空格分隔格式通常为：
      ATOM serial atom_name resname [chain] resseq x y z charge radius
    
    注意：chain 字段可能缺失，需要根据字段Count判断。
    
    Returns:
        解析SuccessReturn元组 (record_type, serial, atom_name, resname, chain, resseq, x, y, z, charge, radius)
        解析FailedReturn None
    """
    try:
        parts = line.split()
        
        # 最少需要 10 个字段（无 chain 时）：record serial name resname resseq x y z charge radius
        # 有 chain 时为 11 个字段：record serial name resname chain resseq x y z charge radius
        if len(parts) < 10:
            return None
        
        record_type = parts[0]
        serial = int(parts[1])
        atom_name = parts[2]
        resname = parts[3]
        
        if len(parts) == 11:
            # Package含 chain 字段
            chain = parts[4]
            resseq = parts[5]
            x = float(parts[6])
            y = float(parts[7])
            z = float(parts[8])
            charge = float(parts[9])
            radius = float(parts[10])
        elif len(parts) == 10:
            # 无 chain 字段
            chain = ' '
            resseq = parts[4]
            x = float(parts[5])
            y = float(parts[6])
            z = float(parts[7])
            charge = float(parts[8])
            radius = float(parts[9])
        else:
            # 字段过多，尝试取最后 5 个为坐标+电荷+半径
            chain = parts[4] if not parts[4].lstrip('-').replace('.', '').isdigit() else ' '
            idx = 5 if chain != ' ' else 4
            resseq = parts[idx]
            # 从末尾往前取 5 个浮点数
            radius = float(parts[-1])
            charge = float(parts[-2])
            z = float(parts[-3])
            y = float(parts[-4])
            x = float(parts[-5])
        
        return (record_type, serial, atom_name, resname, chain, resseq, x, y, z, charge, radius)
        
    except (ValueError, IndexError):
        return None


def _parse_pqr_regex(line: str):
    """正则表达式兜底方案：从 PQR 行中提取所有浮点数。
    
    当固定列和空格分隔Method都Failed时（通常因为大坐标Value导致字段拼接），
    using正则表达式从行中提取所有浮点数Value，然后根据Count推断各字段含义。
    
    PQR 行应Package含 5 个浮点数：x, y, z, charge, radius
    
    Returns:
        解析SuccessReturn元组 (record_type, serial, atom_name, resname, chain, resseq, x, y, z, charge, radius)
        解析FailedReturn None
    """
    import re
    
    try:
        # Confirm是 ATOM/HETATM 行
        if not line.startswith(('ATOM', 'HETATM')):
            return None
        
        record_type = line[0:6].strip()
        
        # 从行中提取所有浮点数（含可选的正负号和小数点）
        # 限制小数位数为 1-4 位，防止贪婪匹配吃掉相邻拼接字段的数字
        # 例如 "0.14501.8240" → ["0.1450", "1.8240"] 而非 ["0.14501"]
        float_pattern = re.compile(r'[+-]?\d+\.\d{1,4}')
        all_floats = float_pattern.findall(line)
        
        if len(all_floats) < 5:
            # 浮点数不足 5 个，无法确定 x,y,z,charge,radius
            return None
        
        # 最后 5 个浮点数应该是 x, y, z, charge, radius
        # （前面可能有混入的数Value如残基Number等，但残基Number通常是整数）
        x = float(all_floats[-5])
        y = float(all_floats[-4])
        z = float(all_floats[-3])
        charge = float(all_floats[-2])
        radius = float(all_floats[-1])
        
        # 从行首提取Index、原子名、残基名等文本字段
        # using另一个正则提取 record_type 后面的整数和文本
        # 典型格式: "ATOM     1  N   MET A   1   ..."
        header_match = re.match(
            r'(ATOM|HETATM)\s+(\d+)\s+(\S+)\s+(\S+)\s*(\S?)\s*(\d+)',
            line
        )
        
        if header_match:
            serial = int(header_match.group(2))
            atom_name = header_match.group(3)
            resname = header_match.group(4)
            chain = header_match.group(5) if header_match.group(5) else ' '
            resseq = header_match.group(6)
        else:
            # 如果连 header 都无法解析，尝试从固定列提取（宽容模式）
            try:
                serial = int(re.search(r'(\d+)', line[6:12]).group(1))
            except (AttributeError, ValueError):
                serial = 1
            atom_name = line[12:16].strip() if len(line) > 16 else 'X'
            resname = line[17:21].strip() if len(line) > 20 else 'UNK'
            chain = line[21:22].strip() if len(line) > 21 and line[21:22].strip() else ' '
            try:
                resseq = re.search(r'(\d+)', line[22:27]).group(1) if len(line) > 26 else '1'
            except (AttributeError, ValueError):
                resseq = '1'
        
        # 合理性校验
        if abs(x) > 9999 or abs(y) > 9999 or abs(z) > 9999:
            return None
        if radius < 0:
            return None
        
        print(f"[_parse_pqr_regex] 🔧 正则兜底解析Success: {record_type} {serial} {atom_name} "
              f"x={x:.3f} y={y:.3f} z={z:.3f} q={charge:.4f} r={radius:.4f}")
        
        return (record_type, serial, atom_name, resname, chain, resseq, x, y, z, charge, radius)
        
    except (ValueError, IndexError):
        return None



def _get_molecule_center(mol: 'Chem.Mol', conformer_id: int = 0) -> np.ndarray:
    """Get the geometric center of a molecule."""
    conf = mol.GetConformer(conformer_id)
    coords = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        coords.append([pos.x, pos.y, pos.z])
    return np.mean(coords, axis=0)


def _visualize_ec_map(obj_name: str, ec_pdb: str, ligand_resname: str = None):
    """Visualize EC map in PyMOL."""
    if not PYMOL_AVAILABLE:
        return
    
    # Load EC map
    ec_obj = 'ec_map'
    cmd.load(ec_pdb, ec_obj)
    
    # Color by B-factor (EC values)
    # minimum/maximum are scaled to B-factor range (usually -1 to 1 or -100 to 100)
    cmd.spectrum('b', 'red_white_green', ec_obj, minimum=-100, maximum=100)
    
    # Show as spheres with transparency
    cmd.show('spheres', ec_obj)
    cmd.set('sphere_scale', 0.15, ec_obj)
    cmd.set('sphere_transparency', 0.4, ec_obj)
    
    # Highlight ligand
    if ligand_resname:
        cmd.show('sticks', f"{obj_name} and resn {ligand_resname}")
        cmd.color('yellow', f"{obj_name} and resn {ligand_resname} and elem C")
    
    # Set transparency for protein
    cmd.set('cartoon_transparency', 0.5, obj_name)
    
    # Zoom to ligand region
    if ligand_resname:
        cmd.zoom(f"{obj_name} and resn {ligand_resname}", buffer=10)
    
    print("[Visualization] EC map loaded as 'ec_map'")
    print("[Visualization] Green = complementary (EC > 0), Red = clash (EC < 0)")


def visualize_ec_surface(obj_name: str, ligand_resname: str,
                         ec_values: np.ndarray = None,
                         surface_points: np.ndarray = None,
                         output_dir: str = None,
                         surface_type: str = 'gaussian',
                         transparency: float = 0.0,
                         ec_result: Dict = None) -> bool:
    """
    Visualize EC as a smooth colored molecular surface (like the reference image).
    
    This creates a smooth surface representation with EC values mapped as colors,
    similar to electrostatic potential surfaces commonly shown in publications.
    
    Args:
        obj_name: PyMOL object name
        ligand_resname: Ligand residue name
        ec_values: EC values at surface points (from calculate_ligand_ec)
        surface_points: Surface point coordinates (from calculate_ligand_ec)
        output_dir: Output directory for temporary files
        surface_type: Surface type ('gaussian', 'solvent', 'molecular')
        transparency: Surface transparency (0.0 = opaque, 1.0 = fully transparent)
        ec_result: Full result dict from calculate_ligand_ec (alternative to ec_values/surface_points)
        
    Returns:
        True if successful
    """
    if not PYMOL_AVAILABLE:
        print("[visualize_ec_surface] ❌ PyMOL required")
        return False
    
    if not NUMPY_AVAILABLE:
        print("[visualize_ec_surface] ❌ NumPy required")
        return False
    
    # Get EC data from result dict if provided
    if ec_result is not None:
        ec_values = ec_result.get('ec_values')
        surface_points = ec_result.get('surface_points')
        output_dir = ec_result.get('output_dir', output_dir)
    
    if ec_values is None or surface_points is None:
        print("[visualize_ec_surface] ❌ Need ec_values and surface_points")
        print("[visualize_ec_surface] 💡 Run calculate_ligand_ec first and pass the result")
        return False
    

    if output_dir: output_dir = os.path.normpath(output_dir)

    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='glint_ec_vis_')
    
    print(f"[visualize_ec_surface] Creating smooth EC surface visualization...")
    
    # Method 1: Create a CGO (Compiled Graphics Object) surface
    # This creates a smooth interpolated surface colored by EC values
    
    ligand_sel = f"{obj_name} and resn {ligand_resname}"
    surface_obj = f"ec_surface_{ligand_resname}"
    
    # First, copy the ligand to a new object for surface generation
    cmd.create(surface_obj, ligand_sel)
    
    # Assign EC values to B-factors of the ligand atoms
    # We need to map surface point EC values back to atoms
    if SCIPY_AVAILABLE:
        # Use KDTree to find nearest surface point for each atom
        tree = cKDTree(surface_points)
        
        # Get ligand atom coordinates
        atom_coords = []
        atom_indices = []
        cmd.iterate_state(1, surface_obj,
                         "atom_coords.append([x,y,z]); atom_indices.append(index)",
                         space={'atom_coords': atom_coords, 'atom_indices': atom_indices})
        
        if atom_coords:
            atom_coords = np.array(atom_coords)
            
            # For each atom, find nearby surface points and average their EC values
            for i, (coord, atom_idx) in enumerate(zip(atom_coords, atom_indices)):
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
    else:
        # Fallback: simple nearest neighbor assignment
        for i, (coord, atom_idx) in enumerate(zip(atom_coords, atom_indices)):
            distances = np.linalg.norm(surface_points - coord, axis=1)
            nearest_idx = np.argmin(distances)
            bfactor = np.clip(ec_values[nearest_idx] * 100, -99.99, 99.99)
            cmd.alter(f"{surface_obj} and index {atom_idx}", f"b={bfactor}")
    
    # Rebuild to apply B-factor changes
    cmd.rebuild(surface_obj)
    
    # Generate surface
    cmd.show('surface', surface_obj)
    
    # Set surface quality
    cmd.set('surface_quality', 1, surface_obj)  # Higher quality
    
    # Color by B-factor (EC values) with red-white-green gradient
    # Red = negative EC (clash), White = neutral, Green = positive EC (complementary)
    cmd.spectrum('b', 'red_white_green', surface_obj, minimum=-100, maximum=100)
    
    # Set surface properties
    cmd.set('surface_type', 0 if surface_type == 'molecular' else
                           1 if surface_type == 'solvent' else 2, surface_obj)  # 2 = gaussian
    cmd.set('transparency', transparency, surface_obj)
    cmd.set('surface_color_smoothing', 1)  # Smooth color transitions
    cmd.set('surface_color_smoothing_threshold', 0.5)
    
    # Also show ligand sticks inside the surface
    cmd.show('sticks', ligand_sel)
    cmd.color('gray50', f"{ligand_sel} and elem C")
    cmd.set('stick_transparency', 0.3, ligand_sel)
    
    # Show protein as lines/cartoon with transparency
    cmd.set('cartoon_transparency', 0.7, obj_name)
    cmd.show('lines', f"{obj_name} and polymer within 5 of {ligand_sel}")
    
    # Zoom to ligand
    cmd.zoom(ligand_sel, buffer=8)
    
    # Set nice rendering options
    cmd.set('ray_shadow', 0)
    cmd.set('antialias', 2)
    cmd.set('ambient', 0.4)
    cmd.set('spec_reflect', 0.5)
    
    print(f"[visualize_ec_surface] ✅ Created surface object: {surface_obj}")
    print("[visualize_ec_surface] Color scheme: Red = clash (EC < 0), White = neutral, Green = complementary (EC > 0)")
    print("[visualize_ec_surface] 💡 Use 'ray' command for high-quality rendering")
    
    return True


def visualize_ec_surface_cgo(obj_name: str, ligand_resname: str,
                             ec_values: np.ndarray,
                             surface_points: np.ndarray,
                             surface_normals: np.ndarray = None,
                             point_size: float = 0.3) -> bool:
    """
    Visualize EC as a CGO (Compiled Graphics Object) point cloud surface.
    
    This creates a dense point cloud that approximates a smooth surface,
    with each point colored by its EC value.
    
    Args:
        obj_name: PyMOL object name
        ligand_resname: Ligand residue name
        ec_values: EC values at surface points
        surface_points: Surface point coordinates
        surface_normals: Surface normals (optional, for lighting)
        point_size: Size of each point
        
    Returns:
        True if successful
    """
    if not PYMOL_AVAILABLE:
        print("[visualize_ec_surface_cgo] ❌ PyMOL required")
        return False
    
    from pymol.cgo import BEGIN, END, VERTEX, COLOR, NORMAL, SPHERE, TRIANGLES
    from pymol.cgo import ALPHA
    
    print(f"[visualize_ec_surface_cgo] Creating CGO surface with {len(surface_points)} points...")
    
    # Create CGO object
    cgo_obj = []
    
    # Color mapping function: EC -> RGB
    def ec_to_rgb(ec_val):
        """Map EC value to RGB color (red-white-green gradient)."""
        # Clip EC to [-1, 1] range
        ec_clipped = np.clip(ec_val, -1.0, 1.0)
        
        if ec_clipped < 0:
            # Negative EC: red to white
            t = 1.0 + ec_clipped  # 0 to 1 as EC goes from -1 to 0
            r = 1.0
            g = t
            b = t
        else:
            # Positive EC: white to green
            t = ec_clipped  # 0 to 1 as EC goes from 0 to 1
            r = 1.0 - t
            g = 1.0
            b = 1.0 - t
        
        return (r, g, b)
    
    # Add spheres for each surface point
    for i, (point, ec_val) in enumerate(zip(surface_points, ec_values)):
        r, g, b = ec_to_rgb(ec_val)
        
        # Add colored sphere
        cgo_obj.extend([
            COLOR, r, g, b,
            SPHERE, point[0], point[1], point[2], point_size
        ])
    
    # Load CGO object
    cgo_name = f"ec_cgo_{ligand_resname}"
    cmd.load_cgo(cgo_obj, cgo_name)
    
    # Show ligand sticks
    ligand_sel = f"{obj_name} and resn {ligand_resname}"
    cmd.show('sticks', ligand_sel)
    cmd.color('gray50', f"{ligand_sel} and elem C")
    
    # Zoom
    cmd.zoom(ligand_sel, buffer=8)
    
    print(f"[visualize_ec_surface_cgo] ✅ Created CGO object: {cgo_name}")
    print("[visualize_ec_surface_cgo] Color scheme: Red = clash, White = neutral, Green = complementary")
    
    return True


def create_ec_surface_mesh(surface_points: np.ndarray, ec_values: np.ndarray,
                           output_file: str, mesh_resolution: float = 0.5) -> bool:
    """
    Create a triangulated mesh surface colored by EC values.
    
    This uses Delaunay triangulation to create a proper mesh surface
    that can be rendered smoothly.
    
    Args:
        surface_points: Surface point coordinates
        ec_values: EC values at surface points
        output_file: Output file path (.obj or .ply format)
        mesh_resolution: Resolution for mesh generation
        
    Returns:
        True if successful
    """
    if not NUMPY_AVAILABLE:
        return False
    
    try:
        from scipy.spatial import Delaunay, ConvexHull
        
        print(f"[create_ec_surface_mesh] Creating mesh from {len(surface_points)} points...")
        
        # Create Delaunay triangulation
        tri = Delaunay(surface_points)
        
        # Filter triangles to keep only surface triangles
        # (remove internal triangles based on edge length)
        max_edge_length = mesh_resolution * 3
        
        valid_simplices = []
        for simplex in tri.simplices:
            pts = surface_points[simplex]
            # Calculate edge lengths
            edges = [
                np.linalg.norm(pts[0] - pts[1]),
                np.linalg.norm(pts[1] - pts[2]),
                np.linalg.norm(pts[2] - pts[0]),
                np.linalg.norm(pts[0] - pts[3]),
                np.linalg.norm(pts[1] - pts[3]),
                np.linalg.norm(pts[2] - pts[3]),
            ]
            if max(edges) < max_edge_length:
                valid_simplices.append(simplex)
        
        # Write to OBJ file with vertex colors
        if output_file.endswith('.obj'):
            with open(output_file, 'w') as f:
                f.write("# EC surface mesh generated by GLINT\n")
                
                # Write vertices with colors
                for i, (point, ec_val) in enumerate(zip(surface_points, ec_values)):
                    # Map EC to color
                    ec_clipped = np.clip(ec_val, -1.0, 1.0)
                    if ec_clipped < 0:
                        r, g, b = 1.0, 1.0 + ec_clipped, 1.0 + ec_clipped
                    else:
                        r, g, b = 1.0 - ec_clipped, 1.0, 1.0 - ec_clipped
                    
                    f.write(f"v {point[0]:.4f} {point[1]:.4f} {point[2]:.4f} {r:.4f} {g:.4f} {b:.4f}\n")
                
                # Write faces (triangles from valid simplices)
                for simplex in valid_simplices:
                    # OBJ uses 1-indexed vertices
                    f.write(f"f {simplex[0]+1} {simplex[1]+1} {simplex[2]+1}\n")
                    f.write(f"f {simplex[0]+1} {simplex[2]+1} {simplex[3]+1}\n")
                    f.write(f"f {simplex[0]+1} {simplex[1]+1} {simplex[3]+1}\n")
                    f.write(f"f {simplex[1]+1} {simplex[2]+1} {simplex[3]+1}\n")
        
        print(f"[create_ec_surface_mesh] ✅ Wrote mesh to {output_file}")
        return True
        
    except Exception as e:
        print(f"[create_ec_surface_mesh] ❌ Error: {e}")
        return False


def _visualize_ternary_ec(obj_name: str, results: Dict):
    """Visualize ternary EC analysis in PyMOL."""
    if not PYMOL_AVAILABLE:
        return
    
    output_dir = results.get('output_dir', '.')
    
    # Load EC maps for both interfaces
    for interface_key in ['A_glue', 'B_glue']:
        if interface_key in results['interfaces']:
            interface_data = results['interfaces'][interface_key]
            chains = interface_data.get('chains', [])
            chains_str = "_".join(chains)
            
            subdir = f'protein_{chains_str}'
            filename = f'ec_map_{chains_str}_glue.pdb'
            ec_pdb = os.path.join(output_dir, subdir, filename)
            
            if os.path.exists(ec_pdb):
                ec_obj = f'ec_{interface_key.lower()}'
                cmd.load(ec_pdb, ec_obj)
                cmd.spectrum('b', 'red_white_green', ec_obj, minimum=-100, maximum=100)
                cmd.show('spheres', ec_obj)
                cmd.set('sphere_scale', 0.15, ec_obj)
                cmd.set('sphere_transparency', 0.4, ec_obj)
    
    # Color proteins
    protein_a_chains = results.get('protein_a_chains', [])
    protein_b_chains = results.get('protein_b_chains', [])
    
    for chain in protein_a_chains:
        cmd.color('cyan', f"{obj_name} and chain {chain}")
    
    for chain in protein_b_chains:
        cmd.color('magenta', f"{obj_name} and chain {chain}")
    
    # Highlight glue
    glue_resname = results.get('glue_resname')
    if glue_resname:
        cmd.show('sticks', f"{obj_name} and resn {glue_resname}")
        cmd.color('orange', f"{obj_name} and resn {glue_resname}")
    
    print("[Visualization] Ternary EC maps loaded")
    print("[Visualization] Cyan = Protein A, Magenta = Protein B, Orange = Glue")
    print("[Visualization] Green = complementary (EC > 0), Red = clash (EC < 0)")


def _get_element_charge(element: str) -> float:
    """Get approximate partial charge for an element (simplified)."""
    # Simplified charge assignment for PPI interface analysis
    charges = {
        'N': -0.3,   # Typically negative (backbone N, Lys NH3+)
        'O': -0.5,   # Typically negative (carbonyl, carboxyl)
        'C': 0.0,    # Neutral
        'S': -0.1,   # Slightly negative
        'H': 0.1,    # Slightly positive
    }
    return charges.get(element.upper(), 0.0)


def compare_ligand_ec(obj_name: str, ligand_resnames: List[str],
                      protein_chains: List[str] = None,
                      output_dir: str = None,
                      ph: float = 7.4,
                      surface_density: float = 10.0) -> Optional[Dict[str, Any]]:
    """
    Compare EC scores for multiple ligands (SAR analysis).
    
    This function calculates EC for multiple ligands in the same binding site,
    useful for understanding structure-activity relationships (SAR) driven by
    electrostatic effects.
    
    Args:
        obj_name: PyMOL object name containing all ligands
        ligand_resnames: List of ligand residue names to compare
        protein_chains: Protein chain IDs
        output_dir: Output directory
        ph: pH for protonation
        surface_density: Surface sampling density
        
    Returns:
        Dictionary with EC comparison results
    """
    if not NUMPY_AVAILABLE or not RDKIT_AVAILABLE:
        print("[compare_ligand_ec] ❌ NumPy and RDKit required")
        return None
    

    if output_dir: output_dir = os.path.normpath(output_dir)

    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='glint_ec_compare_')
    else:
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"[compare_ligand_ec] Comparing {len(ligand_resnames)} ligands")
    print(f"[compare_ligand_ec] Output directory: {output_dir}")
    
    results = {
        'ligands': {},
        'comparison': {},
        'output_dir': output_dir
    }
    
    # Calculate EC for each ligand
    for i, ligand_resname in enumerate(ligand_resnames):
        print(f"\n{'='*60}")
        print(f"Ligand {i+1}/{len(ligand_resnames)}: {ligand_resname}")
        print('='*60)
        
        ligand_dir = os.path.join(output_dir, ligand_resname)
        
        ec_result = calculate_ligand_ec(
            obj_name=obj_name,
            ligand_resname=ligand_resname,
            protein_chains=protein_chains,
            output_dir=ligand_dir,
            ph=ph,
            surface_density=surface_density,
            visualize=False  # Don't visualize each one
        )
        
        if ec_result:
            results['ligands'][ligand_resname] = {
                'ec_score': ec_result['ec_score'],
                'ec_statistics': ec_result['ec_statistics'],
                'output_dir': ligand_dir
            }
    
    # Generate comparison
    if len(results['ligands']) > 1:
        scores = [(name, data['ec_score']) for name, data in results['ligands'].items()]
        scores.sort(key=lambda x: x[1], reverse=True)  # Sort by EC score (higher = better)
        
        results['comparison'] = {
            'ranking': scores,
            'best_ligand': scores[0][0],
            'best_ec_score': scores[0][1],
            'worst_ligand': scores[-1][0],
            'worst_ec_score': scores[-1][1],
            'ec_range': scores[0][1] - scores[-1][1]
        }
        
        print("\n" + "="*60)
        print("EC Comparison Summary (SAR Analysis)")
        print("="*60)
        print(f"{'Rank':<6}{'Ligand':<15}{'EC Score':<12}{'Interpretation'}")
        print("-"*60)
        
        for rank, (name, score) in enumerate(scores, 1):
            if score > 0.3:
                interp = "Strong complementarity"
            elif score > 0:
                interp = "Moderate"
            else:
                interp = "Poor/Clash"
            print(f"{rank:<6}{name:<15}{score:<12.4f}{interp}")
        
        print("-"*60)
        print(f"EC Range: {results['comparison']['ec_range']:.4f}")
        print(f"Best: {results['comparison']['best_ligand']} (EC={results['comparison']['best_ec_score']:.4f})")
    
    # Write comparison CSV
    csv_path = os.path.join(output_dir, 'ec_comparison.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Ligand', 'EC_Score', 'EC_Mean', 'EC_Median', 'EC_Std',
                        'Positive_Fraction', 'Negative_Fraction'])
        for name, data in results['ligands'].items():
            stats = data['ec_statistics']
            writer.writerow([
                name, data['ec_score'], stats['ec_mean'], stats['ec_median'],
                stats['ec_std'], stats['ec_positive_fraction'], stats['ec_negative_fraction']
            ])
    
    print(f"\n[compare_ligand_ec] Comparison saved to: {csv_path}")
    
    return results


def analyze_multiconformer_ec(obj_name: str, ligand_resname: str,
                              conformer_states: List[int] = None,
                              protein_chains: List[str] = None,
                              output_dir: str = None,
                              ph: float = 7.4,
                              surface_density: float = 10.0) -> Optional[Dict[str, Any]]:
    """
    Analyze EC across multiple conformations (e.g., docking poses or MD frames).
    
    This function calculates EC statistics across multiple conformations,
    providing mean, median, and distribution of EC scores.
    
    Args:
        obj_name: PyMOL object name
        ligand_resname: Ligand residue name
        conformer_states: List of PyMOL states to analyze (1-indexed)
        protein_chains: Protein chain IDs
        output_dir: Output directory
        ph: pH for protonation
        surface_density: Surface sampling density
        
    Returns:
        Dictionary with multi-conformer EC analysis
    """
    if not NUMPY_AVAILABLE or not RDKIT_AVAILABLE:
        print("[analyze_multiconformer_ec] ❌ NumPy and RDKit required")
        return None
    
    if not PYMOL_AVAILABLE:
        print("[analyze_multiconformer_ec] ❌ PyMOL required for multi-state analysis")
        return None
    
    # Get number of states if not specified
    if conformer_states is None:
        n_states = cmd.count_states(obj_name)
        if n_states <= 1:
            print("[analyze_multiconformer_ec] Only 1 state found, running single EC analysis")
            return calculate_ligand_ec(obj_name, ligand_resname, protein_chains,
                                       output_dir, ph, surface_density)
        conformer_states = list(range(1, n_states + 1))
    

    if output_dir: output_dir = os.path.normpath(output_dir)

    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='glint_ec_multiconf_')
    else:
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"[analyze_multiconformer_ec] Analyzing {len(conformer_states)} conformations")
    print(f"[analyze_multiconformer_ec] Output directory: {output_dir}")
    
    results = {
        'conformers': {},
        'statistics': {},
        'output_dir': output_dir
    }
    
    ec_scores = []
    
    for state in conformer_states:
        print(f"\n{'='*60}")
        print(f"Conformer State {state}/{len(conformer_states)}")
        print('='*60)
        
        # Set state
        cmd.frame(state)
        
        state_dir = os.path.join(output_dir, f'state_{state}')
        
        # Export current state
        state_pdb = os.path.join(state_dir, 'complex.pdb')
        os.makedirs(state_dir, exist_ok=True)
        cmd.save(state_pdb, obj_name, state=state)
        
        # Run EC analysis
        ec_result = calculate_ligand_ec(
            obj_name=None,  # Use PDB file instead
            ligand_resname=ligand_resname,
            protein_chains=protein_chains,
            output_dir=state_dir,
            ph=ph,
            surface_density=surface_density,
            visualize=False,
            pdb_file=state_pdb
        )
        
        if ec_result:
            results['conformers'][state] = {
                'ec_score': ec_result['ec_score'],
                'ec_statistics': ec_result['ec_statistics']
            }
            ec_scores.append(ec_result['ec_score'])
    
    # Calculate overall statistics
    if ec_scores:
        ec_scores = np.array(ec_scores)
        results['statistics'] = {
            'n_conformers': len(ec_scores),
            'ec_mean': float(np.mean(ec_scores)),
            'ec_median': float(np.median(ec_scores)),
            'ec_std': float(np.std(ec_scores)),
            'ec_min': float(np.min(ec_scores)),
            'ec_max': float(np.max(ec_scores)),
            'ec_q25': float(np.percentile(ec_scores, 25)),
            'ec_q75': float(np.percentile(ec_scores, 75)),
            'best_conformer': int(conformer_states[np.argmax(ec_scores)]),
            'worst_conformer': int(conformer_states[np.argmin(ec_scores)])
        }
        
        print("\n" + "="*60)
        print("Multi-Conformer EC Analysis Summary")
        print("="*60)
        print(f"Number of conformers: {results['statistics']['n_conformers']}")
        print(f"EC Mean: {results['statistics']['ec_mean']:.4f}")
        print(f"EC Median: {results['statistics']['ec_median']:.4f}")
        print(f"EC Std: {results['statistics']['ec_std']:.4f}")
        print(f"EC Range: [{results['statistics']['ec_min']:.4f}, {results['statistics']['ec_max']:.4f}]")
        print(f"EC IQR: [{results['statistics']['ec_q25']:.4f}, {results['statistics']['ec_q75']:.4f}]")
        print(f"Best conformer: State {results['statistics']['best_conformer']} (EC={results['statistics']['ec_max']:.4f})")
        print(f"Worst conformer: State {results['statistics']['worst_conformer']} (EC={results['statistics']['ec_min']:.4f})")
        
        # Interpretation
        if results['statistics']['ec_mean'] > 0.3:
            print("\n✨ Strong average electrostatic complementarity")
        elif results['statistics']['ec_mean'] > 0:
            print("\n✓ Moderate average electrostatic complementarity")
        else:
            print("\n⚠️ Poor average electrostatic complementarity")
        
        if results['statistics']['ec_std'] > 0.2:
            print("⚠️ High variability across conformers - binding may be pose-dependent")
    
    # Write summary CSV
    csv_path = os.path.join(output_dir, 'multiconformer_ec.csv')
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['State', 'EC_Score', 'EC_Mean', 'EC_Positive_Fraction'])
        for state, data in results['conformers'].items():
            writer.writerow([
                state, data['ec_score'],
                data['ec_statistics']['ec_mean'],
                data['ec_statistics']['ec_positive_fraction']
            ])
    
    print(f"\n[analyze_multiconformer_ec] Results saved to: {csv_path}")
    
    return results


def calculate_ec_hotspots(obj_name: str, ligand_resname: str,
                         protein_chains: List[str] = None,
                         output_dir: str = None,
                         ph: float = 7.4,
                         surface_density: float = 15.0,
                         hotspot_threshold: float = 0.5) -> Optional[Dict[str, Any]]:
    """
    Identify EC hotspots - regions of strong electrostatic complementarity.
    
    This function identifies specific regions on the ligand surface that
    contribute most to electrostatic complementarity, useful for guiding
    SAR optimization.
    
    Args:
        obj_name: PyMOL object name
        ligand_resname: Ligand residue name
        protein_chains: Protein chain IDs
        output_dir: Output directory
        ph: pH for protonation
        surface_density: Higher density for better hotspot resolution
        hotspot_threshold: EC threshold for hotspot identification
        
    Returns:
        Dictionary with hotspot analysis
    """
    if not NUMPY_AVAILABLE or not RDKIT_AVAILABLE:
        print("[calculate_ec_hotspots] ❌ NumPy and RDKit required")
        return None
    

    if output_dir: output_dir = os.path.normpath(output_dir)

    # First run standard EC analysis with higher density
    ec_result = calculate_ligand_ec(
        obj_name=obj_name,
        ligand_resname=ligand_resname,
        protein_chains=protein_chains,
        output_dir=output_dir,
        ph=ph,
        surface_density=surface_density,
        visualize=False
    )
    
    if ec_result is None:
        return None
    
    surface_points = ec_result['surface_points']
    ec_values = ec_result['ec_values']
    phi_protein = ec_result['phi_protein']
    phi_ligand = ec_result['phi_ligand']
    output_dir = ec_result['output_dir']
    
    # Identify hotspots (high positive EC)
    positive_hotspots = ec_values > hotspot_threshold
    negative_hotspots = ec_values < -hotspot_threshold
    
    results = {
        'ec_result': ec_result,
        'hotspots': {
            'positive': {
                'count': int(np.sum(positive_hotspots)),
                'fraction': float(np.mean(positive_hotspots)),
                'points': surface_points[positive_hotspots],
                'ec_values': ec_values[positive_hotspots],
                'mean_ec': float(np.mean(ec_values[positive_hotspots])) if np.any(positive_hotspots) else 0.0
            },
            'negative': {
                'count': int(np.sum(negative_hotspots)),
                'fraction': float(np.mean(negative_hotspots)),
                'points': surface_points[negative_hotspots],
                'ec_values': ec_values[negative_hotspots],
                'mean_ec': float(np.mean(ec_values[negative_hotspots])) if np.any(negative_hotspots) else 0.0
            }
        },
        'output_dir': output_dir
    }
    
    # Cluster hotspots to find distinct regions
    if SCIPY_AVAILABLE and np.sum(positive_hotspots) > 5:
        from scipy.cluster.hierarchy import fcluster, linkage
        
        pos_points = surface_points[positive_hotspots]
        if len(pos_points) > 1:
            # Hierarchical clustering
            Z = linkage(pos_points, method='average')
            clusters = fcluster(Z, t=3.0, criterion='distance')  # 3Å cluster distance
            
            n_clusters = len(np.unique(clusters))
            results['hotspots']['positive']['n_clusters'] = n_clusters
            
            # Get cluster centers
            cluster_centers = []
            for c in np.unique(clusters):
                mask = clusters == c
                center = np.mean(pos_points[mask], axis=0)
                cluster_centers.append(center)
            results['hotspots']['positive']['cluster_centers'] = np.array(cluster_centers)
    
    # Write hotspot PDBs
    if np.sum(positive_hotspots) > 0:
        pos_pdb = os.path.join(output_dir, 'ec_hotspots_positive.pdb')
        ECMapWriter.write_pseudo_pdb(pos_pdb,
                                     surface_points[positive_hotspots],
                                     ec_values[positive_hotspots],
                                     atom_name='POS')
    
    if np.sum(negative_hotspots) > 0:
        neg_pdb = os.path.join(output_dir, 'ec_hotspots_negative.pdb')
        ECMapWriter.write_pseudo_pdb(neg_pdb,
                                     surface_points[negative_hotspots],
                                     ec_values[negative_hotspots],
                                     atom_name='NEG')
    
    print("\n" + "="*60)
    print("EC Hotspot Analysis")
    print("="*60)
    print(f"Hotspot threshold: |EC| > {hotspot_threshold}")
    print(f"\nPositive hotspots (complementary):")
    print(f"  Count: {results['hotspots']['positive']['count']} points")
    print(f"  Fraction: {results['hotspots']['positive']['fraction']*100:.1f}%")
    print(f"  Mean EC: {results['hotspots']['positive']['mean_ec']:.4f}")
    if 'n_clusters' in results['hotspots']['positive']:
        print(f"  Distinct regions: {results['hotspots']['positive']['n_clusters']}")
    
    print(f"\nNegative hotspots (clash):")
    print(f"  Count: {results['hotspots']['negative']['count']} points")
    print(f"  Fraction: {results['hotspots']['negative']['fraction']*100:.1f}%")
    print(f"  Mean EC: {results['hotspots']['negative']['mean_ec']:.4f}")
    
    # Visualize in PyMOL
    if PYMOL_AVAILABLE:
        _visualize_ec_hotspots(obj_name, results, ligand_resname)
    
    return results


def _visualize_ec_hotspots(obj_name: str, results: Dict, ligand_resname: str = None):
    """Visualize EC hotspots in PyMOL."""
    if not PYMOL_AVAILABLE:
        return
    
    output_dir = results.get('output_dir', '.')
    
    # Load positive hotspots
    pos_pdb = os.path.join(output_dir, 'ec_hotspots_positive.pdb')
    if os.path.exists(pos_pdb):
        cmd.load(pos_pdb, 'ec_hotspots_pos')
        cmd.color('green', 'ec_hotspots_pos')
        cmd.show('spheres', 'ec_hotspots_pos')
        cmd.set('sphere_scale', 0.2, 'ec_hotspots_pos')
        cmd.set('sphere_transparency', 0.3, 'ec_hotspots_pos')
    
    # Load negative hotspots
    neg_pdb = os.path.join(output_dir, 'ec_hotspots_negative.pdb')
    if os.path.exists(neg_pdb):
        cmd.load(neg_pdb, 'ec_hotspots_neg')
        cmd.color('red', 'ec_hotspots_neg')
        cmd.show('spheres', 'ec_hotspots_neg')
        cmd.set('sphere_scale', 0.2, 'ec_hotspots_neg')
    
    # Highlight ligand
    if ligand_resname:
        cmd.show('sticks', f"{obj_name} and resn {ligand_resname}")
        cmd.color('yellow', f"{obj_name} and resn {ligand_resname} and elem C")
    
    print("[Visualization] EC hotspots loaded")
    print("[Visualization] Blue = complementary hotspots, Red = clash hotspots")


def analyze_substituent_ec_effect(obj_name: str, ligand_resname: str,
                                  substituent_atoms: List[str],
                                  protein_chains: List[str] = None,
                                  output_dir: str = None,
                                  ph: float = 7.4) -> Optional[Dict[str, Any]]:
    """
    Analyze EC contribution of specific substituent atoms.
    
    This function calculates the EC contribution of specific atoms/groups
    on the ligand, useful for understanding SAR (e.g., CF3 vs CH3 effect).
    
    Args:
        obj_name: PyMOL object name
        ligand_resname: Ligand residue name
        substituent_atoms: List of atom names to analyze (e.g., ['F1', 'F2', 'F3'])
        protein_chains: Protein chain IDs
        output_dir: Output directory
        ph: pH for protonation
        
    Returns:
        Dictionary with substituent EC analysis
    """
    if not NUMPY_AVAILABLE or not RDKIT_AVAILABLE:
        print("[analyze_substituent_ec_effect] ❌ NumPy and RDKit required")
        return None
    

    if output_dir: output_dir = os.path.normpath(output_dir)

    # First run full EC analysis
    ec_result = calculate_ligand_ec(
        obj_name=obj_name,
        ligand_resname=ligand_resname,
        protein_chains=protein_chains,
        output_dir=output_dir,
        ph=ph,
        surface_density=15.0,  # Higher density for better resolution
        visualize=False
    )
    
    if ec_result is None:
        return None
    
    surface_points = ec_result['surface_points']
    ec_values = ec_result['ec_values']
    output_dir = ec_result['output_dir']
    
    # Get substituent atom coordinates from PyMOL
    if not PYMOL_AVAILABLE:
        print("[analyze_substituent_ec_effect] ❌ PyMOL required")
        return None
    
    substituent_coords = []
    for atom_name in substituent_atoms:
        sel = f"{obj_name} and resn {ligand_resname} and name {atom_name}"
        coords = []
        cmd.iterate_state(1, sel, "coords.append([x,y,z])", space={'coords': coords})
        if coords:
            substituent_coords.append(np.array(coords[0]))
    
    if not substituent_coords:
        print("[analyze_substituent_ec_effect] ❌ No substituent atoms found")
        return None
    
    substituent_coords = np.array(substituent_coords)
    substituent_center = np.mean(substituent_coords, axis=0)
    
    # Find surface points near substituent (within 4Å of any substituent atom)
    near_substituent = np.zeros(len(surface_points), dtype=bool)
    for coord in substituent_coords:
        distances = np.linalg.norm(surface_points - coord, axis=1)
        near_substituent |= (distances < 4.0)
    
    results = {
        'full_ec_result': ec_result,
        'substituent_atoms': substituent_atoms,
        'substituent_center': substituent_center,
        'substituent_analysis': {
            'n_surface_points': int(np.sum(near_substituent)),
            'ec_mean': float(np.mean(ec_values[near_substituent])) if np.any(near_substituent) else 0.0,
            'ec_std': float(np.std(ec_values[near_substituent])) if np.any(near_substituent) else 0.0,
            'ec_positive_fraction': float(np.mean(ec_values[near_substituent] > 0)) if np.any(near_substituent) else 0.0,
            'contribution_to_total': float(np.sum(ec_values[near_substituent]) / np.sum(np.abs(ec_values))) if np.sum(np.abs(ec_values)) > 0 else 0.0
        },
        'output_dir': output_dir
    }
    
    # Compare to rest of ligand
    rest_of_ligand = ~near_substituent
    if np.any(rest_of_ligand):
        results['rest_of_ligand'] = {
            'n_surface_points': int(np.sum(rest_of_ligand)),
            'ec_mean': float(np.mean(ec_values[rest_of_ligand])),
            'ec_std': float(np.std(ec_values[rest_of_ligand])),
            'ec_positive_fraction': float(np.mean(ec_values[rest_of_ligand] > 0))
        }
    
    print("\n" + "="*60)
    print("Substituent EC Effect Analysis")
    print("="*60)
    print(f"Substituent atoms: {substituent_atoms}")
    print(f"\nSubstituent region:")
    print(f"  Surface points: {results['substituent_analysis']['n_surface_points']}")
    print(f"  EC Mean: {results['substituent_analysis']['ec_mean']:.4f}")
    print(f"  EC Std: {results['substituent_analysis']['ec_std']:.4f}")
    print(f"  Positive EC fraction: {results['substituent_analysis']['ec_positive_fraction']*100:.1f}%")
    print(f"  Contribution to total: {results['substituent_analysis']['contribution_to_total']*100:.1f}%")
    
    if 'rest_of_ligand' in results:
        print(f"\nRest of ligand:")
        print(f"  EC Mean: {results['rest_of_ligand']['ec_mean']:.4f}")
        print(f"  Positive EC fraction: {results['rest_of_ligand']['ec_positive_fraction']*100:.1f}%")
        
        # Interpretation
        sub_ec = results['substituent_analysis']['ec_mean']
        rest_ec = results['rest_of_ligand']['ec_mean']
        if sub_ec > rest_ec + 0.1:
            print(f"\n✨ Substituent has BETTER EC than rest of ligand (+{sub_ec - rest_ec:.3f})")
        elif sub_ec < rest_ec - 0.1:
            print(f"\n⚠️ Substituent has WORSE EC than rest of ligand ({sub_ec - rest_ec:.3f})")
        else:
            print(f"\n→ Substituent has similar EC to rest of ligand")
    
    # Write substituent surface points
    if np.any(near_substituent):
        sub_pdb = os.path.join(output_dir, 'ec_substituent.pdb')
        ECMapWriter.write_pseudo_pdb(sub_pdb,
                                     surface_points[near_substituent],
                                     ec_values[near_substituent],
                                     atom_name='SUB')
        
        # Visualize
        if PYMOL_AVAILABLE:
            cmd.load(sub_pdb, 'ec_substituent')
            cmd.spectrum('b', 'red_white_green', 'ec_substituent', minimum=-100, maximum=100)
            cmd.show('spheres', 'ec_substituent')
            cmd.set('sphere_scale', 0.2, 'ec_substituent')
            cmd.set('sphere_transparency', 0.4, 'ec_substituent')
    
    return results


# ========== PyMOL Command Registration ==========

if PYMOL_AVAILABLE:
    cmd.extend('calculate_ligand_ec', calculate_ligand_ec)
    cmd.extend('analyze_ternary_ec', analyze_ternary_ec)
    cmd.extend('compare_ligand_ec', compare_ligand_ec)
    cmd.extend('analyze_multiconformer_ec', analyze_multiconformer_ec)
    cmd.extend('calculate_ec_hotspots', calculate_ec_hotspots)
    cmd.extend('analyze_substituent_ec_effect', analyze_substituent_ec_effect)
    cmd.extend('visualize_ec_surface', visualize_ec_surface)
    cmd.extend('visualize_ec_surface_cgo', visualize_ec_surface_cgo)
    
    # Register enhanced visualization commands if available
    if EC_VISUALIZATION_AVAILABLE:
        cmd.extend('visualize_ec_smooth_surface', visualize_ec_smooth_surface)
        cmd.extend('save_ec_visualization', save_ec_visualization)
        cmd.extend('visualize_ternary_ec_surfaces', visualize_ternary_ec_surfaces)


# ========== Module Info ==========

if __name__ == '__main__':
    print("="*60)
    print("GLINT Ligand EC Calculator")
    print("="*60)
    print("\nThis module calculates Electrostatic Complementarity (EC)")
    print("for protein-ligand and ternary (molecular glue) complexes.")
    print("\nUsage in PyMOL:")
    print("  calculate_ligand_ec 'complex', 'LIG', output_dir='./ec_output'")
    print("  analyze_ternary_ec 'complex', 'GLUE', ['A'], ['B']")
    print("  compare_ligand_ec 'complex', ['LIG1', 'LIG2', 'LIG3']")
    print("  analyze_multiconformer_ec 'complex', 'LIG'")
    print("  calculate_ec_hotspots 'complex', 'LIG'")
    print("  analyze_substituent_ec_effect 'complex', 'LIG', ['F1', 'F2', 'F3']")
    print("\nCore Dependencies:")
    print(f"  NumPy: {'✅' if NUMPY_AVAILABLE else '❌'}")
    print(f"  RDKit: {'✅' if RDKIT_AVAILABLE else '❌'}")
    print(f"  SciPy: {'✅' if SCIPY_AVAILABLE else '❌'}")
    print(f"  PyMOL: {'✅' if PYMOL_AVAILABLE else '❌'}")
    print("\nElectrostatics Tools:")
    print(f"  PDB2PQR Python API: {'✅' if PDB2PQR_PYTHON_AVAILABLE else '❌ (pip install pdb2pqr)'}")
    print(f"  APBS Python API: {'✅' if APBS_PYTHON_AVAILABLE else '❌ (pip install apbs or use CLI)'}")
    print(f"  APBS Binary Pkg: {'✅' if APBS_BINARY_AVAILABLE else '❌ (pip install apbs-binary)'}")
    print("\nNote: If Python APIs are not available, command-line tools will be used as fallback.")