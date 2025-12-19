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

Usage:
    # Basic EC analysis
    result = calculate_ligand_ec('complex', 'LIG', output_dir='./ec_output')
    
    # Ternary complex (molecular glue) analysis
    result = analyze_ternary_ec('complex', 'GLUE', ['A'], ['B'], output_dir='./ec_output')

Author: GlueTK Team
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
    """
    
    def __init__(self, dielectric: float = 4.0):
        """
        Args:
            dielectric: Effective dielectric constant for Coulomb calculation
        """
        self.dielectric = dielectric
    
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
        
        Args:
            mol: RDKit molecule with Gasteiger charges
            points: Nx3 array of points to evaluate
            conformer_id: Which conformer to use
            
        Returns:
            Array of potential values at each point
        """
        if not RDKIT_AVAILABLE:
            raise RuntimeError("RDKit required for potential calculation")
        
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
                print(f"[PDB2PQR] Error: {result.stderr}")
                return False
            
            if os.path.exists(output_pqr):
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


# Try to import APBS Python module
try:
    import apbs
    APBS_PYTHON_AVAILABLE = True
except ImportError:
    APBS_PYTHON_AVAILABLE = False


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
    
    def _find_apbs(self) -> str:
        """Try to find APBS in PATH or common locations."""
        path = shutil.which('apbs')
        if path:
            return path
        
        # Try Python module
        if APBS_PYTHON_AVAILABLE:
            return 'python_api'
        
        # Common installation locations
        common_paths = [
            '/usr/local/bin/apbs',
            '/usr/bin/apbs',
            os.path.expanduser('~/apbs/bin/apbs'),
            'C:\\Program Files\\APBS\\apbs.exe'
        ]
        
        for p in common_paths:
            if os.path.exists(p):
                return p
        
        return 'apbs'  # Hope it's in PATH
    
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
        input_content = f"""# APBS input file generated by GlueTK
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
  swin 0.3
  sdens 10.0
  temp {APBS_PARAMS['temp']}
  calcenergy no
  calcforce no

  dime {grid_dims[0]} {grid_dims[1]} {grid_dims[2]}
  cglen {coarse_len[0]} {coarse_len[1]} {coarse_len[2]}
  fglen {fine_len[0]} {fine_len[1]} {fine_len[2]}
  cgcent mol 1
  fgcent mol 1

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
        Run APBS calculation.
        
        Args:
            input_file: APBS input file
            working_dir: Working directory (use input file dir if None)
            
        Returns:
            Path to output DX file, or None if failed
        """
        if working_dir is None:
            working_dir = os.path.dirname(input_file) or '.'
        
        # Try Python API first
        if self.use_python_api and APBS_PYTHON_AVAILABLE:
            result = self._run_python_api(input_file, working_dir)
            if result:
                return result
            print("[APBS] Python API failed, falling back to command line...")
        
        # Fall back to command line
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
    
    def _run_command_line(self, input_file: str, working_dir: str) -> Optional[str]:
        """Run APBS using command line."""
        cmd_parts = [self.apbs_path, input_file]
        
        print(f"[APBS] Running: {' '.join(cmd_parts)}")
        print(f"[APBS] Working directory: {working_dir}")
        
        try:
            result = subprocess.run(
                cmd_parts,
                capture_output=True,
                text=True,
                cwd=working_dir,
                timeout=600  # 10 minute timeout
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
    
    - Opposite signs (complementary): EC > 0
    - Same signs (clash): EC < 0
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
        numerator = -2.0 * phi_p * phi_l
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
        if weights is None:
            return float(np.mean(ec_local))
        else:
            return float(np.average(ec_local, weights=weights))
    
    def calculate_ec_statistics(self, ec_local: np.ndarray) -> Dict[str, float]:
        """
        Calculate various statistics for EC distribution.
        
        Args:
            ec_local: Array of local EC values
            
        Returns:
            Dictionary of statistics
        """
        return {
            'ec_mean': float(np.mean(ec_local)),
            'ec_median': float(np.median(ec_local)),
            'ec_std': float(np.std(ec_local)),
            'ec_min': float(np.min(ec_local)),
            'ec_max': float(np.max(ec_local)),
            'ec_positive_fraction': float(np.mean(ec_local > 0)),
            'ec_negative_fraction': float(np.mean(ec_local < 0)),
            'ec_q25': float(np.percentile(ec_local, 25)),
            'ec_q75': float(np.percentile(ec_local, 75))
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
            f.write("REMARK EC map generated by GlueTK\n")
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
                       pdb_file: str = None) -> Optional[Dict[str, Any]]:
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
        
    Returns:
        Dictionary containing:
        - ec_score: Overall EC score
        - ec_statistics: Detailed statistics
        - surface_points: Surface point coordinates
        - ec_values: EC value at each point
        - output_files: Paths to generated files
    """
    # Check dependencies
    if not NUMPY_AVAILABLE:
        print("[calculate_ligand_ec] ❌ NumPy required")
        return None
    
    if not RDKIT_AVAILABLE:
        print("[calculate_ligand_ec] ❌ RDKit required")
        return None
    
    # Create output directory
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='gluetk_ec_')
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
    
    if PYMOL_AVAILABLE and obj_name:
        # Use PyMOL selections
        protein_sel = f"{obj_name} and polymer"
        if protein_chains:
            chain_sel = ' or '.join([f"chain {c}" for c in protein_chains])
            protein_sel = f"({protein_sel}) and ({chain_sel})"
        
        cmd.save(protein_pdb, protein_sel)
        
        if ligand_resname:
            cmd.save(ligand_sdf, f"{obj_name} and resn {ligand_resname}", format='sdf')
        else:
            # Auto-detect ligand
            cmd.save(ligand_sdf, f"{obj_name} and organic and not polymer", format='sdf')
    else:
        # Parse PDB file manually
        _extract_protein_ligand_from_pdb(complex_pdb, protein_pdb, ligand_sdf, ligand_resname)
    
    # Load ligand with RDKit
    supplier = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
    ligand_mol = supplier[0] if len(supplier) > 0 else None
    
    if ligand_mol is None:
        print("[calculate_ligand_ec] ❌ Failed to load ligand")
        return None
    
    print(f"[calculate_ligand_ec] Ligand: {ligand_mol.GetNumAtoms()} atoms")
    
    # Step 2: Run PDB2PQR
    print("\n[Step 2] Running PDB2PQR...")
    protein_pqr = os.path.join(output_dir, 'protein.pqr')
    
    pdb2pqr = PDB2PQRRunner()
    if not pdb2pqr.run(protein_pdb, protein_pqr, ph=ph):
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
    
    if dx_file is None:
        print("[calculate_ligand_ec] ❌ APBS failed")
        return None
    
    # Step 4: Load protein potential grid
    print("\n[Step 4] Loading protein potential...")
    protein_grid = DXGrid(dx_file)
    
    # Step 5: Generate ligand surface points
    print("\n[Step 5] Generating ligand surface...")
    sampler = LigandSurfaceSampler(density=surface_density)
    surface_points, surface_normals = sampler.sample_molecule(ligand_mol)
    
    # Step 6: Calculate potentials at surface points
    print("\n[Step 6] Calculating potentials...")
    
    # Protein potential (from APBS grid)
    phi_protein = protein_grid.interpolate(surface_points)
    
    # Ligand potential (Gasteiger charges + Coulomb)
    charge_calc = GasteigerChargeCalculator()
    phi_ligand = charge_calc.calculate_potential(ligand_mol, surface_points)
    
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
        'output_dir': output_dir
    }
    
    return result


def analyze_ternary_ec(obj_name: str = None, glue_resname: str = None,
                      protein_a_chains: List[str] = None,
                      protein_b_chains: List[str] = None,
                      output_dir: str = None,
                      ph: float = 7.4,
                      surface_density: float = 10.0,
                      visualize: bool = True,
                      pdb_file: str = None) -> Optional[Dict[str, Any]]:
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
        
    Returns:
        Dictionary containing EC results for each interface
    """
    if not NUMPY_AVAILABLE or not RDKIT_AVAILABLE:
        print("[analyze_ternary_ec] ❌ NumPy and RDKit required")
        return None
    
    if not protein_a_chains or not protein_b_chains:
        print("[analyze_ternary_ec] ❌ Must specify both protein_a_chains and protein_b_chains")
        return None
    
    # Create output directory
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='gluetk_ternary_ec_')
    else:
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"[analyze_ternary_ec] Output directory: {output_dir}")
    print(f"[analyze_ternary_ec] Protein A chains: {protein_a_chains}")
    print(f"[analyze_ternary_ec] Protein B chains: {protein_b_chains}")
    print(f"[analyze_ternary_ec] Glue: {glue_resname}")
    
    results = {
        'glue_resname': glue_resname,
        'protein_a_chains': protein_a_chains,
        'protein_b_chains': protein_b_chains,
        'interfaces': {}
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
    
    # Generate glue surface points (shared for both interfaces)
    sampler = LigandSurfaceSampler(density=surface_density)
    surface_points, surface_normals = sampler.sample_molecule(glue_mol)
    
    # Calculate glue potential (shared)
    charge_calc = GasteigerChargeCalculator()
    phi_glue = charge_calc.calculate_potential(glue_mol, surface_points)
    
    # Interface 1: Protein A - Glue
    print("\n" + "="*60)
    print("Interface 1: Protein A - Glue")
    print("="*60)
    
    protein_a_dir = os.path.join(output_dir, 'protein_a')
    os.makedirs(protein_a_dir, exist_ok=True)
    
    # Extract protein A
    protein_a_pdb = os.path.join(protein_a_dir, 'protein_a.pdb')
    if PYMOL_AVAILABLE and obj_name:
        chain_sel = ' or '.join([f"chain {c}" for c in protein_a_chains])
        cmd.save(protein_a_pdb, f"{obj_name} and polymer and ({chain_sel})")
    
    # Run PDB2PQR and APBS for protein A
    protein_a_pqr = os.path.join(protein_a_dir, 'protein_a.pqr')
    pdb2pqr = PDB2PQRRunner()
    pdb2pqr.run(protein_a_pdb, protein_a_pqr, ph=ph)
    
    apbs = APBSRunner()
    apbs_prefix_a = os.path.join(protein_a_dir, 'protein_a_pot')
    apbs_input_a = apbs.generate_input(protein_a_pqr, apbs_prefix_a)
    dx_file_a = apbs.run(apbs_input_a, protein_a_dir)
    
    if dx_file_a:
        grid_a = DXGrid(dx_file_a)
        phi_protein_a = grid_a.interpolate(surface_points)
        
        ec_calc = ECCalculator()
        ec_a_glue = ec_calc.calculate_ec_local(phi_protein_a, phi_glue)
        ec_score_a = ec_calc.calculate_ec_score(ec_a_glue)
        ec_stats_a = ec_calc.calculate_ec_statistics(ec_a_glue)
        
        results['interfaces']['A_glue'] = {
            'ec_score': ec_score_a,
            'ec_statistics': ec_stats_a,
            'ec_values': ec_a_glue,
            'phi_protein': phi_protein_a
        }
        
        print(f"EC(A-Glue) Score: {ec_score_a:.4f}")
        
        # Write EC map
        ec_pdb_a = os.path.join(protein_a_dir, 'ec_map_a_glue.pdb')
        ECMapWriter.write_pseudo_pdb(ec_pdb_a, surface_points, ec_a_glue)
    
    # Interface 2: Protein B - Glue
    print("\n" + "="*60)
    print("Interface 2: Protein B - Glue")
    print("="*60)
    
    protein_b_dir = os.path.join(output_dir, 'protein_b')
    os.makedirs(protein_b_dir, exist_ok=True)
    
    # Extract protein B
    protein_b_pdb = os.path.join(protein_b_dir, 'protein_b.pdb')
    if PYMOL_AVAILABLE and obj_name:
        chain_sel = ' or '.join([f"chain {c}" for c in protein_b_chains])
        cmd.save(protein_b_pdb, f"{obj_name} and polymer and ({chain_sel})")
    
    # Run PDB2PQR and APBS for protein B
    protein_b_pqr = os.path.join(protein_b_dir, 'protein_b.pqr')
    pdb2pqr.run(protein_b_pdb, protein_b_pqr, ph=ph)
    
    apbs_prefix_b = os.path.join(protein_b_dir, 'protein_b_pot')
    apbs_input_b = apbs.generate_input(protein_b_pqr, apbs_prefix_b)
    dx_file_b = apbs.run(apbs_input_b, protein_b_dir)
    
    if dx_file_b:
        grid_b = DXGrid(dx_file_b)
        phi_protein_b = grid_b.interpolate(surface_points)
        
        ec_b_glue = ec_calc.calculate_ec_local(phi_protein_b, phi_glue)
        ec_score_b = ec_calc.calculate_ec_score(ec_b_glue)
        ec_stats_b = ec_calc.calculate_ec_statistics(ec_b_glue)
        
        results['interfaces']['B_glue'] = {
            'ec_score': ec_score_b,
            'ec_statistics': ec_stats_b,
            'ec_values': ec_b_glue,
            'phi_protein': phi_protein_b
        }
        
        print(f"EC(B-Glue) Score: {ec_score_b:.4f}")
        
        # Write EC map
        ec_pdb_b = os.path.join(protein_b_dir, 'ec_map_b_glue.pdb')
        ECMapWriter.write_pseudo_pdb(ec_pdb_b, surface_points, ec_b_glue)
    
    # Combined analysis
    print("\n" + "="*60)
    print("Combined Ternary EC Analysis")
    print("="*60)
    
    if 'A_glue' in results['interfaces'] and 'B_glue' in results['interfaces']:
        ec_a = results['interfaces']['A_glue']['ec_score']
        ec_b = results['interfaces']['B_glue']['ec_score']
        
        # Combined score (average)
        ec_combined = (ec_a + ec_b) / 2
        
        # Asymmetry (difference between interfaces)
        ec_asymmetry = abs(ec_a - ec_b)
        
        results['combined'] = {
            'ec_combined_score': ec_combined,
            'ec_asymmetry': ec_asymmetry,
            'ec_a_glue': ec_a,
            'ec_b_glue': ec_b
        }
        
        print(f"EC(A-Glue): {ec_a:.4f}")
        print(f"EC(B-Glue): {ec_b:.4f}")
        print(f"Combined EC: {ec_combined:.4f}")
        print(f"Asymmetry: {ec_asymmetry:.4f}")
        
        # Interpretation
        if ec_combined > 0.3:
            print("\n✨ Strong electrostatic complementarity - favorable glue binding")
        elif ec_combined > 0:
            print("\n✓ Moderate electrostatic complementarity")
        else:
            print("\n⚠️ Poor electrostatic complementarity - potential clash")
    
    # Interface 3: Protein A - Protein B (PPI interface with glue)
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
        
        # Get interface atoms (within 5Å of each other)
        interface_a_sel = f"({obj_name} and ({chain_a_sel})) within 5.0 of ({obj_name} and ({chain_b_sel}))"
        interface_b_sel = f"({obj_name} and ({chain_b_sel})) within 5.0 of ({obj_name} and ({chain_a_sel}))"
        
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
                
                # Calculate EC at PPI interface
                # Use protein B potential at interface A surface
                if dx_file_b:
                    phi_b_at_interface = grid_b.interpolate(ppi_points_a)
                    
                    # Calculate "potential" from interface A atoms (simplified)
                    # Use distance-weighted charge approximation
                    phi_a_at_interface = np.zeros(len(ppi_points_a))
                    for i, point in enumerate(ppi_points_a):
                        distances = np.linalg.norm(np.array(interface_a_coords) - point, axis=1)
                        distances = np.maximum(distances, 0.5)
                        # Simplified: assume partial charges based on element
                        charges = np.array([_get_element_charge(e) for e in interface_a_elements])
                        phi_a_at_interface[i] = np.sum(charges / distances) * 332.0637 / 0.593
                    
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
    """Simple PDB to PQR conversion with basic charges (fallback)."""
    # Simple charge assignment based on residue type
    charges = {
        'ARG': {'NH1': 0.5, 'NH2': 0.5, 'NE': 0.0},
        'LYS': {'NZ': 1.0},
        'ASP': {'OD1': -0.5, 'OD2': -0.5},
        'GLU': {'OE1': -0.5, 'OE2': -0.5},
        'HIS': {'ND1': 0.25, 'NE2': 0.25}
    }
    
    # Default radii
    radii = {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8, 'H': 1.2}
    
    try:
        with open(pdb_file, 'r') as f_in, open(pqr_file, 'w') as f_out:
            for line in f_in:
                if line.startswith(('ATOM', 'HETATM')):
                    resname = line[17:20].strip()
                    atomname = line[12:16].strip()
                    element = line[76:78].strip() if len(line) > 76 else atomname[0]
                    
                    # Get charge
                    charge = 0.0
                    if resname in charges and atomname in charges[resname]:
                        charge = charges[resname][atomname]
                    
                    # Get radius
                    radius = radii.get(element.upper(), 1.7)
                    
                    # Write PQR line
                    pqr_line = (
                        f"{line[:54]}"
                        f"{charge:8.4f}"
                        f"{radius:7.4f}\n"
                    )
                    f_out.write(pqr_line)
        
        return pqr_file
    except Exception as e:
        print(f"[_simple_pdb_to_pqr] Error: {e}")
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
    cmd.spectrum('b', 'blue_white_red', ec_obj, minimum=-100, maximum=100)
    
    # Show as spheres
    cmd.show('spheres', ec_obj)
    cmd.set('sphere_scale', 0.15, ec_obj)
    
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
    print("[Visualization] Blue = complementary (EC > 0), Red = clash (EC < 0)")


def _visualize_ternary_ec(obj_name: str, results: Dict):
    """Visualize ternary EC analysis in PyMOL."""
    if not PYMOL_AVAILABLE:
        return
    
    output_dir = results.get('output_dir', '.')
    
    # Load EC maps for both interfaces
    for interface in ['A_glue', 'B_glue']:
        if interface in results['interfaces']:
            subdir = 'protein_a' if interface == 'A_glue' else 'protein_b'
            ec_pdb = os.path.join(output_dir, subdir, f'ec_map_{interface.lower()}.pdb')
            
            if os.path.exists(ec_pdb):
                ec_obj = f'ec_{interface}'
                cmd.load(ec_pdb, ec_obj)
                cmd.spectrum('b', 'blue_white_red', ec_obj, minimum=-100, maximum=100)
                cmd.show('spheres', ec_obj)
                cmd.set('sphere_scale', 0.15, ec_obj)
    
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
    
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='gluetk_ec_compare_')
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
    
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix='gluetk_ec_multiconf_')
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
        cmd.color('blue', 'ec_hotspots_pos')
        cmd.show('spheres', 'ec_hotspots_pos')
        cmd.set('sphere_scale', 0.2, 'ec_hotspots_pos')
    
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
            cmd.spectrum('b', 'blue_white_red', 'ec_substituent', minimum=-100, maximum=100)
            cmd.show('spheres', 'ec_substituent')
            cmd.set('sphere_scale', 0.2, 'ec_substituent')
    
    return results


# ========== PyMOL Command Registration ==========

if PYMOL_AVAILABLE:
    cmd.extend('calculate_ligand_ec', calculate_ligand_ec)
    cmd.extend('analyze_ternary_ec', analyze_ternary_ec)
    cmd.extend('compare_ligand_ec', compare_ligand_ec)
    cmd.extend('analyze_multiconformer_ec', analyze_multiconformer_ec)
    cmd.extend('calculate_ec_hotspots', calculate_ec_hotspots)
    cmd.extend('analyze_substituent_ec_effect', analyze_substituent_ec_effect)


# ========== Module Info ==========

if __name__ == '__main__':
    print("="*60)
    print("GlueTK Ligand EC Calculator")
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
    print("\nNote: If Python APIs are not available, command-line tools will be used as fallback.")