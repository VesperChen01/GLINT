# -*- coding: utf-8 -*-
"""
Protein Surface Analyzer
========================
Analyzes protein surface properties including electrostatic potential,
hydrophobicity, and geometric features to identify potential binding sites.

Features:
- Electrostatic patch detection (APBS/PDB2PQR or Coulomb fallback)
- Hydrophobic patch detection
- Geometric feature analysis (curvature, shape index)
- Comprehensive druggability scoring

Updated: Now supports APBS/PDB2PQR for accurate electrostatic potential calculation.
"""

import numpy as np
from scipy.spatial.distance import cdist
from scipy import ndimage
import math
import os
import tempfile
import shutil
from typing import List, Dict, Optional, Tuple, Any

try:
    from pymol import cmd
except ImportError:
    cmd = None

from .feature_extractor import MolecularFeatureExtractor
from .pocket_detector import PocketDetector

# Residue-level chemistry used for lightweight surface summaries.
HYDROPHOBIC_RESIDUES = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO', 'GLY'}
POLAR_RESIDUES = {'SER', 'THR', 'CYS', 'TYR', 'ASN', 'GLN'}
POSITIVE_RESIDUES = {'LYS', 'ARG', 'HIS'}
NEGATIVE_RESIDUES = {'ASP', 'GLU'}
AROMATIC_RESIDUES = {'PHE', 'TRP', 'TYR', 'HIS'}

HBOND_DONOR_ATOMS = {
    ('LYS', 'NZ'), ('ARG', 'NE'), ('ARG', 'NH1'), ('ARG', 'NH2'),
    ('HIS', 'ND1'), ('HIS', 'NE2'), ('SER', 'OG'), ('THR', 'OG1'),
    ('TYR', 'OH'), ('ASN', 'ND2'), ('GLN', 'NE2'), ('TRP', 'NE1'),
    ('CYS', 'SG'),
}
HBOND_ACCEPTOR_ATOMS = {
    ('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), ('GLU', 'OE2'),
    ('ASN', 'OD1'), ('GLN', 'OE1'), ('SER', 'OG'), ('THR', 'OG1'),
    ('TYR', 'OH'), ('HIS', 'ND1'), ('HIS', 'NE2'), ('CYS', 'SG'),
}

# Try to import APBS/PDB2PQR runners
try:
    from .ligand_ec_calculator import PDB2PQRRunner, APBSRunner, DXGrid
    APBS_AVAILABLE = True
except ImportError:
    APBS_AVAILABLE = False
    PDB2PQRRunner = None
    APBSRunner = None
    DXGrid = None

class SurfacePatch:
    """Represents a continuous patch on the protein surface"""
    
    def __init__(self, patch_id: int, patch_type: str):
        self.id = patch_id
        self.type = patch_type  # 'electrostatic_pos', 'electrostatic_neg', 'hydrophobic', 'concave'
        self.points: List[Tuple[float, float, float]] = []  # Surface points (x,y,z)
        self.residues: List[Dict] = []  # Associated residues
        self.area: float = 0.0
        self.center: Tuple[float, float, float] = (0.0, 0.0, 0.0)
        
        # Properties
        self.avg_potential: float = 0.0  # For electrostatic patches
        self.avg_hydrophobicity: float = 0.0  # For hydrophobic patches
        self.curvature: float = 0.0
        self.avg_curvature: float = 0.0
        self.avg_signed_curvature: float = 0.0
        self.avg_concavity: float = 0.0
        self.shape_index: float = 0.0  # -1 (cup) to +1 (cap)
        self.score: float = 0.0  # Druggability/Importance score
        self.grid_indices: List[Tuple[int, int, int]] = []

        # Residue chemistry summary for this patch.
        self.hydrophobic_fraction: float = 0.0
        self.polar_fraction: float = 0.0
        self.positive_fraction: float = 0.0
        self.negative_fraction: float = 0.0
        self.aromatic_fraction: float = 0.0
        self.hbond_donors: int = 0
        self.hbond_acceptors: int = 0
        self.avg_hbond_role: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            'id': self.id,
            'type': self.type,
            'area': float(self.area),
            'center': [float(x) for x in self.center],
            'avg_potential': float(self.avg_potential),
            'avg_hydrophobicity': float(self.avg_hydrophobicity),
            'curvature': float(self.curvature),
            'avg_curvature': float(self.avg_curvature),
            'avg_signed_curvature': float(self.avg_signed_curvature),
            'avg_concavity': float(self.avg_concavity),
            'shape_index': float(self.shape_index),
            'score': float(self.score),
            'hydrophobic_fraction': float(self.hydrophobic_fraction),
            'polar_fraction': float(self.polar_fraction),
            'positive_fraction': float(self.positive_fraction),
            'negative_fraction': float(self.negative_fraction),
            'aromatic_fraction': float(self.aromatic_fraction),
            'hbond_donors': int(self.hbond_donors),
            'hbond_acceptors': int(self.hbond_acceptors),
            'avg_hbond_role': float(self.avg_hbond_role),
            'num_residues': len(self.residues),
            'residues': ";".join([f"{r['chain']}:{r['resn']}:{r['resi']}" for r in self.residues])
        }

class SurfaceAnalyzer:
    """Analyzes protein surface for features"""
    
    def __init__(self, obj_name: str, grid_spacing: float = 0.5, 
                 use_apbs: bool = True, ph: float = 7.4,
                 output_dir: str = None, selection: str = "all"):
        """
        Initialize surface analyzer.
        
        Args:
            obj_name: PyMOL object name
            grid_spacing: Grid spacing for surface generation (Å)
            use_apbs: Use APBS/PDB2PQR for accurate electrostatics (default: True)
            ph: pH for PDB2PQR protonation state (default: 7.4)
            output_dir: Directory for temporary files (auto-created if None)
            selection: PyMOL selection inside obj_name to analyze/render
        """
        self.obj_name = obj_name
        self.selection = selection or "all"
        self.grid_spacing = grid_spacing
        self.use_apbs = use_apbs and APBS_AVAILABLE
        self.ph = ph
        self.output_dir = output_dir
        self.atoms = []
        self.grid = None
        self.origin = None
        self.dx_data = None  # APBS electrostatic potential grid
        self._temp_dir = None
        self._surface_cache = None
        
    def analyze(self) -> List[SurfacePatch]:
        """Main analysis pipeline"""
        if not cmd:
            print("Error: PyMOL not available")
            return []
            
        print(f"Starting surface analysis for {self.obj_name}...")
        print(f"  Use APBS: {self.use_apbs} (available: {APBS_AVAILABLE})")
        
        # 1. Get Atoms
        self.atoms = self._get_atoms()
        if not self.atoms:
            return []
        
        # 2. Run APBS/PDB2PQR if enabled
        if self.use_apbs:
            print("  Running APBS/PDB2PQR for accurate electrostatics...")
            self.dx_data = self._run_apbs_pipeline()
            if self.dx_data:
                print(f"  ✅ APBS completed: grid shape {self.dx_data['data'].shape}")
            else:
                print("  ⚠️ APBS failed, falling back to residue-based approximation")
            
        # 3. Build Grid & Surface
        # Use PocketDetector's robust grid generation
        pd = PocketDetector(grid_spacing=self.grid_spacing)
        self.grid_info, self.origin, atom_coords = pd._build_grid(self.atoms)
        occupied = pd._mark_occupied(self.grid_info, atom_coords, self.origin)
        
        # Generate surface points (solvent accessible surface)
        # Dilate occupied grid to find surface
        struct = ndimage.generate_binary_structure(3, 1)
        dilated = ndimage.binary_dilation(occupied, structure=struct)
        surface_mask = dilated & (~occupied)
        surface_points_indices = np.argwhere(surface_mask)
        
        if len(surface_points_indices) == 0:
            return []
            
        surface_coords = self.origin + surface_points_indices * self.grid_spacing
        
        print(f"Generated {len(surface_coords)} surface points")

        # 4. Calculate Properties for Surface Points
        # We need to map surface points to nearest atoms to get properties
        atom_kdtree = None
        try:
            from scipy.spatial import cKDTree
            atom_kdtree = cKDTree(atom_coords)
        except ImportError:
            pass
            
        properties = self._map_properties_to_surface(surface_coords, surface_points_indices, atom_kdtree, self.atoms)
        
        # 5. Detect Patches
        patches = []
        patches.extend(self._detect_electrostatic_patches(surface_points_indices, properties, surface_coords))
        patches.extend(self._detect_hydrophobic_patches(surface_points_indices, properties, surface_coords))
        patches.extend(self._detect_curvature_patches(surface_points_indices, properties, surface_coords))
        patches.extend(self._detect_hbond_patches(surface_points_indices, properties, surface_coords))
        
        # 6. Calculate Geometric Features & Score
        self._calculate_geometry_and_score(patches, surface_coords)
        
        # 7. Cleanup temporary files
        self._cleanup()
        
        print(f"Found {len(patches)} surface patches")
        return patches
    
    def _run_apbs_pipeline(self) -> Optional[Dict]:
        """
        Run PDB2PQR and APBS to calculate electrostatic potential.
        
        This method automatically calculates optimal grid parameters based on
        protein size, similar to PyMOL's APBS Electrostatics plugin.
        
        Returns:
            Dict with 'data', 'origin', 'spacing' keys, or None if failed
        """
        if not APBS_AVAILABLE:
            return None
        
        try:
            # Create temporary directory
            if self.output_dir:
                work_dir = self.output_dir
                os.makedirs(work_dir, exist_ok=True)
            else:
                self._temp_dir = tempfile.mkdtemp(prefix='glint_surface_')
                work_dir = self._temp_dir
            
            # Step 1: Export protein to PDB
            protein_pdb = os.path.join(work_dir, 'protein.pdb')
            cmd.save(protein_pdb, self._selection_expr())
            
            if not os.path.exists(protein_pdb) or os.path.getsize(protein_pdb) == 0:
                print("  ❌ Failed to export protein PDB")
                return None
            
            # Step 2: Run PDB2PQR
            protein_pqr = os.path.join(work_dir, 'protein.pqr')
            pdb2pqr = PDB2PQRRunner()
            
            print(f"  Running PDB2PQR (pH={self.ph})...")
            if not pdb2pqr.run(protein_pdb, protein_pqr, ph=self.ph):
                print("  ❌ PDB2PQR failed")
                return None
            
            if not os.path.exists(protein_pqr) or os.path.getsize(protein_pqr) == 0:
                print("  ❌ PDB2PQR output file not found")
                return None
            
            # Step 3: Calculate grid parameters based on protein size
            # This mimics PyMOL APBS plugin behavior for optimal results
            coords = np.array([a['coord'] for a in self.atoms])
            center = coords.mean(axis=0)
            
            # Calculate protein extent (bounding box)
            min_coords = coords.min(axis=0)
            max_coords = coords.max(axis=0)
            extent = max_coords - min_coords
            max_extent = max(extent)
            
            print(f"  Protein extent: {extent[0]:.1f} x {extent[1]:.1f} x {extent[2]:.1f} Å")
            print(f"  Max dimension: {max_extent:.1f} Å")
            
            # Calculate grid parameters like PyMOL APBS plugin
            # Fine grid should cover the protein with some padding
            # Coarse grid should be ~2x the fine grid
            # Grid spacing should be ~0.5 Å for good resolution
            
            # Add padding (20 Å on each side for solvent)
            padding = 20.0
            fine_len = extent + 2 * padding
            
            # Coarse grid is typically 2x fine grid
            coarse_len = fine_len * 2.0
            
            # Calculate grid dimensions for ~0.5 Å spacing
            # APBS requires odd dimensions
            grid_spacing = 0.5
            grid_dims = np.ceil(fine_len / grid_spacing).astype(int)
            # Make dimensions odd (APBS requirement)
            grid_dims = grid_dims + (1 - grid_dims % 2)
            # Ensure minimum size of 65 and maximum of 225
            grid_dims = np.clip(grid_dims, 65, 225)
            
            print(f"  Grid dimensions: {grid_dims[0]} x {grid_dims[1]} x {grid_dims[2]}")
            print(f"  Fine grid: {fine_len[0]:.1f} x {fine_len[1]:.1f} x {fine_len[2]:.1f} Å")
            print(f"  Coarse grid: {coarse_len[0]:.1f} x {coarse_len[1]:.1f} x {coarse_len[2]:.1f} Å")
            
            # Step 4: Run APBS with calculated parameters
            apbs = APBSRunner()
            apbs_prefix = os.path.join(work_dir, 'protein_pot')
            
            print("  Generating APBS input...")
            apbs_input = apbs.generate_input(
                protein_pqr, apbs_prefix,
                grid_center=center.tolist(),
                grid_dims=tuple(grid_dims.tolist()),
                coarse_len=tuple(coarse_len.tolist()),
                fine_len=tuple(fine_len.tolist())
            )
            
            print("  Running APBS...")
            dx_file = apbs.run(apbs_input, work_dir)
            
            if not dx_file or not os.path.exists(dx_file):
                print("  ❌ APBS failed to generate DX file")
                return None
            
            # Step 5: Read DX file using DXGrid class
            print(f"  Reading DX file: {os.path.basename(dx_file)}")
            dx_grid = DXGrid(dx_file)
            
            # Convert DXGrid to dict format expected by the rest of the code
            dx_data = {
                'data': dx_grid.data,
                'origin': dx_grid.origin,
                'spacing': np.diag(dx_grid.delta),  # Extract diagonal (spacing in each dimension)
                'dims': dx_grid.dims,
                'interpolator': dx_grid  # Keep reference to DXGrid for interpolation
            }
            
            print(f"  ✅ DX file loaded: shape={dx_grid.dims}, range=[{dx_grid.data.min():.2f}, {dx_grid.data.max():.2f}] kT/e")
            
            return dx_data
            
        except Exception as e:
            print(f"  ❌ APBS pipeline error: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def _cleanup(self):
        """Clean up temporary files"""
        if self._temp_dir and os.path.exists(self._temp_dir):
            try:
                shutil.rmtree(self._temp_dir)
            except Exception as e:
                print(f"Warning: Failed to cleanup temp dir: {e}")

    def _get_atoms(self):
        """Extract atoms from PyMOL"""
        atoms = []
        try:
            model = cmd.get_model(self._selection_expr())
            for atom in model.atom:
                atoms.append({
                    'coord': np.array(atom.coord),
                    'element': atom.symbol.upper(),
                    'resn': atom.resn,
                    'resi': atom.resi,
                    'chain': atom.chain,
                    'partial_charge': getattr(atom, 'partial_charge', 0.0), # Requires PDB2PQR or similar usually
                    'name': atom.name
                })
        except Exception as e:
            print(f"Error getting atoms: {e}")
        return atoms

    def _selection_expr(self) -> str:
        """Return the PyMOL expression for the active analysis scope."""
        selection = (self.selection or "all").strip()
        if not selection or selection == "all":
            return self.obj_name
        return f"({self.obj_name}) and ({selection})"

    def _map_properties_to_surface(self, surface_coords, surface_indices, kdtree, atoms):
        """
        Map electrostatic potential and hydrophobicity to surface points.
        
        If APBS data is available, uses accurate electrostatic potential from DX file.
        Otherwise, falls back to residue-based approximation.
        
        Returns dict with 'potential' and 'hydrophobicity' arrays.
        """
        n_points = len(surface_coords)
        potentials = np.zeros(n_points)
        hydrophobicities = np.zeros(n_points)  # 1.0 = hydrophobic, 0.0 = hydrophilic
        hbond_roles = np.zeros(n_points)  # +1 donor, -1 acceptor, 0 neutral/unknown
        nearest_atom_indices = np.zeros(n_points, dtype=int)
        
        # Hydrophobicity scale (Kyte-Doolittle normalized to 0-1)
        hydro_scale = {
            'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
            'MET': 1.9, 'ALA': 1.8, 'GLY': -0.4, 'THR': -0.7, 'SER': -0.8,
            'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6, 'HIS': -3.2, 'GLU': -3.5,
            'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5
        }
        
        # Normalize to 0-1
        for k in hydro_scale:
            hydro_scale[k] = (hydro_scale[k] + 4.5) / 9.0

        # === Method 1: Use APBS electrostatic potential (accurate) ===
        if self.dx_data is not None:
            print("  Using APBS electrostatic potential (accurate)")
            potentials = self._sample_dx_at_points(surface_coords)
            
            # Map hydrophobicity from nearest atoms
            if kdtree:
                dists, indices = kdtree.query(surface_coords, k=1)
                nearest_atom_indices = np.asarray(indices, dtype=int)
                for i, idx in enumerate(indices):
                    atom = atoms[idx]
                    resn = atom['resn']
                    hydrophobicities[i] = hydro_scale.get(resn, 0.5)
                    hbond_roles[i] = self._get_hbond_role(atom)
        
        # === Method 2: Fallback to residue-based approximation ===
        else:
            print("  Using residue-based electrostatic approximation (fallback)")
            if kdtree:
                # Find nearest atom for each surface point
                dists, indices = kdtree.query(surface_coords, k=1)
                nearest_atom_indices = np.asarray(indices, dtype=int)
                
                for i, idx in enumerate(indices):
                    atom = atoms[idx]
                    resn = atom['resn']
                    
                    # Hydrophobicity
                    hydrophobicities[i] = hydro_scale.get(resn, 0.5)
                    hbond_roles[i] = self._get_hbond_role(atom)
                    
                    # Electrostatics (residue-based approximation)
                    charge = atom.get('partial_charge', 0.0)
                    if charge == 0.0:
                        # Use residue type to infer charge
                        if resn in ['ASP', 'GLU']:
                            charge = -1.0
                        elif resn in ['LYS', 'ARG']:
                            charge = 1.0
                        elif resn == 'HIS':
                            charge = 0.5  # pH dependent
                    
                    potentials[i] = charge

        curvature = self._compute_surface_curvature(
            surface_coords,
            np.array([a['coord'] for a in atoms]),
            nearest_atom_indices,
        )

        self._surface_cache = {
            'surface_coords': surface_coords,
            'surface_indices': surface_indices,
            'nearest_atom_indices': nearest_atom_indices,
            'properties': {
                'potential': potentials,
                'hydrophobicity': hydrophobicities,
                'hbond_role': hbond_roles,
                **curvature,
            },
        }

        return self._surface_cache['properties']

    def _get_hbond_role(self, atom: Dict[str, Any]) -> float:
        """Approximate whether an atom behaves as a donor or acceptor on the surface."""
        resn = str(atom.get('resn', '')).upper()
        name = str(atom.get('name', '')).upper()
        element = str(atom.get('element', '')).upper()

        donor = (resn, name) in HBOND_DONOR_ATOMS or (name == 'N' and element == 'N')
        acceptor = (resn, name) in HBOND_ACCEPTOR_ATOMS or (name == 'O' and element == 'O')

        if donor and not acceptor:
            return 1.0
        if acceptor and not donor:
            return -1.0
        return 0.0

    def _compute_surface_curvature(self, surface_coords: np.ndarray,
                                   atom_coords: np.ndarray,
                                   nearest_atom_indices: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Estimate local surface curvature from neighboring surface points.

        The signed value is a local geometric heuristic:
        negative values indicate concave/cup-like regions, positive values indicate
        convex/cap-like regions, and values near zero indicate flatter surface.
        """
        n_points = len(surface_coords)
        curvature = np.zeros(n_points)
        signed_curvature = np.zeros(n_points)
        shape_index = np.zeros(n_points)

        if n_points == 0:
            return {
                'curvature': curvature,
                'signed_curvature': signed_curvature,
                'concavity': curvature,
                'shape_index': shape_index,
            }

        try:
            from scipy.spatial import cKDTree
            tree = cKDTree(surface_coords)
            radius = max(2.0, self.grid_spacing * 4.0)

            for i, point in enumerate(surface_coords):
                neighbor_ids = tree.query_ball_point(point, r=radius)
                if len(neighbor_ids) < 6:
                    _, neighbor_ids = tree.query(point, k=min(12, n_points))
                    neighbor_ids = np.atleast_1d(neighbor_ids).tolist()

                neighbors = surface_coords[neighbor_ids]
                centered = neighbors - neighbors.mean(axis=0)
                cov = np.cov(centered, rowvar=False)
                eigvals = np.linalg.eigvalsh(cov)
                eigvals = np.clip(eigvals, 0.0, None)
                denom = float(eigvals.sum())
                if denom <= 1e-12:
                    continue

                # Surface variation: small on flat patches, larger in highly curved regions.
                curv_mag = float(eigvals[0] / denom)
                curvature[i] = curv_mag

                nearest_atom = atom_coords[nearest_atom_indices[i]]
                normal = point - nearest_atom
                norm = np.linalg.norm(normal)
                if norm > 1e-8:
                    normal = normal / norm
                    local_offset = neighbors.mean(axis=0) - point
                    # Convex surfaces have neighbor centroids slightly inward.
                    sign = -1.0 if np.dot(local_offset, normal) > 0 else 1.0
                else:
                    sign = 1.0

                signed_curvature[i] = sign * curv_mag
                shape_index[i] = float(np.tanh(signed_curvature[i] * 24.0))

        except Exception as e:
            print(f"  Warning: curvature estimation failed: {e}")

        concavity = np.maximum(-signed_curvature, 0.0)
        return {
            'curvature': curvature,
            'signed_curvature': signed_curvature,
            'concavity': concavity,
            'shape_index': shape_index,
        }
    
    def _sample_dx_at_points(self, points: np.ndarray) -> np.ndarray:
        """
        Sample electrostatic potential from DX grid at given points.
        
        Uses DXGrid's interpolate method which leverages scipy's RegularGridInterpolator
        for efficient and accurate trilinear interpolation.
        
        Args:
            points: Nx3 array of coordinates
            
        Returns:
            N-length array of potential values (in kT/e)
        """
        if self.dx_data is None:
            return np.zeros(len(points))
        
        # Use DXGrid's interpolate method if available (more efficient)
        if 'interpolator' in self.dx_data and self.dx_data['interpolator'] is not None:
            dx_grid = self.dx_data['interpolator']
            return dx_grid.interpolate(points)
        
        # Fallback to manual trilinear interpolation
        data = self.dx_data['data']
        origin = np.array(self.dx_data['origin'])
        spacing = np.array(self.dx_data['spacing'])
        
        # Convert world coordinates to grid indices
        grid_coords = (points - origin) / spacing
        
        # Get grid dimensions
        nx, ny, nz = data.shape
        
        potentials = np.zeros(len(points))
        
        for i, gc in enumerate(grid_coords):
            # Check bounds
            if (gc[0] < 0 or gc[0] >= nx - 1 or
                gc[1] < 0 or gc[1] >= ny - 1 or
                gc[2] < 0 or gc[2] >= nz - 1):
                potentials[i] = 0.0
                continue
            
            # Trilinear interpolation
            x0, y0, z0 = int(gc[0]), int(gc[1]), int(gc[2])
            xd, yd, zd = gc[0] - x0, gc[1] - y0, gc[2] - z0
            
            # Get 8 corner values
            c000 = data[x0, y0, z0]
            c001 = data[x0, y0, z0 + 1]
            c010 = data[x0, y0 + 1, z0]
            c011 = data[x0, y0 + 1, z0 + 1]
            c100 = data[x0 + 1, y0, z0]
            c101 = data[x0 + 1, y0, z0 + 1]
            c110 = data[x0 + 1, y0 + 1, z0]
            c111 = data[x0 + 1, y0 + 1, z0 + 1]
            
            # Interpolate
            c00 = c000 * (1 - xd) + c100 * xd
            c01 = c001 * (1 - xd) + c101 * xd
            c10 = c010 * (1 - xd) + c110 * xd
            c11 = c011 * (1 - xd) + c111 * xd
            
            c0 = c00 * (1 - yd) + c10 * yd
            c1 = c01 * (1 - yd) + c11 * yd
            
            potentials[i] = c0 * (1 - zd) + c1 * zd
        
        return potentials

    def _detect_electrostatic_patches(self, surface_indices, properties, surface_coords):
        """
        Cluster positive and negative electrostatic regions.
        
        Thresholds are adaptive based on whether APBS or fallback method was used:
        - APBS: Uses kT/e units, typical range -5 to +5
        - Fallback: Uses discrete charges -1, 0, +1
        """
        patches = []
        pot = properties['potential']
        
        # Adaptive thresholds based on data source
        if self.dx_data is not None:
            # APBS data: values in kT/e, typically -5 to +5
            # Use ±1 kT/e as threshold for significant electrostatic regions
            pos_threshold = 1.0   # kT/e
            neg_threshold = -1.0  # kT/e
            print(f"  Electrostatic thresholds (APBS): pos > {pos_threshold} kT/e, neg < {neg_threshold} kT/e")
        else:
            # Fallback: discrete charges
            pos_threshold = 0.5
            neg_threshold = -0.5
            print(f"  Electrostatic thresholds (fallback): pos > {pos_threshold}, neg < {neg_threshold}")
        
        pos_mask = pot > pos_threshold
        neg_mask = pot < neg_threshold
        
        # Calculate average potential for patches
        pos_patches = self._cluster_surface(surface_indices[pos_mask], "electrostatic_pos", surface_coords[pos_mask])
        neg_patches = self._cluster_surface(surface_indices[neg_mask], "electrostatic_neg", surface_coords[neg_mask])
        
        # Store average potential in patches
        for p in pos_patches:
            self._assign_patch_property_averages(p, surface_indices, properties)
        for p in neg_patches:
            self._assign_patch_property_averages(p, surface_indices, properties)
        
        patches.extend(pos_patches)
        patches.extend(neg_patches)
        
        return patches

    def _detect_hydrophobic_patches(self, surface_indices, properties, surface_coords):
        """Cluster hydrophobic regions"""
        hy = properties['hydrophobicity']
        mask = hy > 0.6 # High hydrophobicity

        patches = self._cluster_surface(surface_indices[mask], "hydrophobic", surface_coords[mask])
        for p in patches:
            self._assign_patch_property_averages(p, surface_indices, properties)
        return patches

    def _detect_curvature_patches(self, surface_indices, properties, surface_coords):
        """Cluster concave high-curvature regions."""
        concavity = properties.get('concavity')
        if concavity is None or len(concavity) == 0:
            return []

        nonzero = concavity[concavity > 0]
        if len(nonzero) == 0:
            return []

        threshold = max(float(np.percentile(nonzero, 75)), 0.015)
        mask = concavity >= threshold
        print(f"  Curvature threshold: concavity >= {threshold:.4f}")

        patches = self._cluster_surface(surface_indices[mask], "concave", surface_coords[mask])
        for p in patches:
            self._assign_patch_property_averages(p, surface_indices, properties)
        return patches

    def _detect_hbond_patches(self, surface_indices, properties, surface_coords):
        """Cluster approximate donor and acceptor surface regions."""
        role = properties.get('hbond_role')
        if role is None or len(role) == 0:
            return []

        donor_mask = role > 0.5
        acceptor_mask = role < -0.5
        donor_patches = self._cluster_surface(surface_indices[donor_mask], "hbond_donor", surface_coords[donor_mask])
        acceptor_patches = self._cluster_surface(surface_indices[acceptor_mask], "hbond_acceptor", surface_coords[acceptor_mask])

        for p in donor_patches + acceptor_patches:
            self._assign_patch_property_averages(p, surface_indices, properties)

        return donor_patches + acceptor_patches

    def _cluster_surface(self, indices, ptype, coords):
        """Generic clustering for surface points using distance"""
        if len(coords) < 10: return []
        
        patches = []
        # Use DBSCAN or simple connectivity if we have grid indices
        # Since we have grid indices, connectivity is fast
        
        # Convert list of indices back to dense grid for labeling
        # Find bounds to minimize grid size
        if len(indices) == 0: return []
        
        mins = indices.min(axis=0)
        maxs = indices.max(axis=0)
        shape = maxs - mins + 1
        
        local_grid = np.zeros(shape, dtype=bool)
        for idx in indices:
            local_pos = idx - mins
            local_grid[tuple(local_pos)] = True
            
        # Label connected components
        labeled, num_features = ndimage.label(local_grid, structure=ndimage.generate_binary_structure(3, 2)) # 18-connectivity
        
        for lbl in range(1, num_features + 1):
            mask = (labeled == lbl)
            # Count points
            points_count = np.sum(mask)
            
            # Min patch size (approx 10 points ~ small patch)
            if points_count < 10: continue
            
            # Create patch object
            patch = SurfacePatch(patch_id=len(patches), patch_type=ptype)
            
            # Get real coords
            loc_indices = np.argwhere(mask)
            glob_indices = loc_indices + mins
            real_coords = self.origin + glob_indices * self.grid_spacing
            patch.grid_indices = [tuple(int(v) for v in idx) for idx in glob_indices]
            patch.points = [tuple(float(v) for v in coord) for coord in real_coords]
            
            # Map back to subset coords to filter them
            # This is tricky because we need original coords match
            # Simplified: calculate center and area
            patch.area = points_count * (self.grid_spacing ** 2) # Approx area
            
            # Calculate center
            patch_center_idx = glob_indices.mean(axis=0)
            patch.center = self.origin + patch_center_idx * self.grid_spacing
            
            # Find nearest residues (simplified)
            # We will populate residues later
            
            patches.append(patch)
            
        return patches

    def _assign_patch_property_averages(self, patch: SurfacePatch,
                                        surface_indices: np.ndarray,
                                        properties: Dict[str, np.ndarray]) -> None:
        """Populate patch-level average electrostatic, hydrophobic, and curvature values."""
        if not patch.grid_indices:
            return

        index_to_pos = {tuple(int(v) for v in idx): i for i, idx in enumerate(surface_indices)}
        positions = [index_to_pos[idx] for idx in patch.grid_indices if idx in index_to_pos]
        if not positions:
            return

        patch.avg_potential = float(np.mean(properties['potential'][positions]))
        patch.avg_hydrophobicity = float(np.mean(properties['hydrophobicity'][positions]))
        patch.avg_curvature = float(np.mean(properties['curvature'][positions]))
        patch.avg_signed_curvature = float(np.mean(properties['signed_curvature'][positions]))
        patch.avg_concavity = float(np.mean(properties['concavity'][positions]))
        patch.avg_hbond_role = float(np.mean(properties['hbond_role'][positions]))
        patch.curvature = patch.avg_curvature
        patch.shape_index = float(np.mean(properties['shape_index'][positions]))

    def _calculate_geometry_and_score(self, patches, all_surface_coords):
        """
        Calculate curvature and other geometric features.
        Score patches based on area and geometry.
        """
        for p in patches:
            # Scoring
            # Area is good, concavity is good for binding
            # Hydrophobic: larger area = better
            # Electrostatic: extreme charge = better
            
            size_score = min(p.area / 100.0, 1.0) # 100 A^2 max score
            concavity_score = min(p.avg_concavity / 0.05, 1.0)
            if p.type == "concave":
                p.score = 0.6 * size_score + 0.4 * concavity_score
            else:
                p.score = size_score * (1.0 + concavity_score)
            
            # Identify residues
            # Need to search atoms near p.center
            # Using simple distance check
            p.residues = self._find_residues_near(p.center, 6.0) # 6A radius
            self._assign_patch_residue_summary(p)

    def _assign_patch_residue_summary(self, patch: SurfacePatch) -> None:
        """Summarize residue chemistry classes around a patch."""
        residues = patch.residues or []
        n_res = len(residues)
        if n_res == 0:
            return

        resn_list = [str(r.get('resn', '')).upper() for r in residues]
        patch.hydrophobic_fraction = sum(r in HYDROPHOBIC_RESIDUES for r in resn_list) / n_res
        patch.polar_fraction = sum(r in POLAR_RESIDUES for r in resn_list) / n_res
        patch.positive_fraction = sum(r in POSITIVE_RESIDUES for r in resn_list) / n_res
        patch.negative_fraction = sum(r in NEGATIVE_RESIDUES for r in resn_list) / n_res
        patch.aromatic_fraction = sum(r in AROMATIC_RESIDUES for r in resn_list) / n_res

        residue_keys = {(r['chain'], r['resn'], r['resi']) for r in residues}
        donors = 0
        acceptors = 0
        for atom in self.atoms:
            key = (atom['chain'], atom['resn'], atom['resi'])
            if key not in residue_keys:
                continue
            role = self._get_hbond_role(atom)
            if role > 0:
                donors += 1
            elif role < 0:
                acceptors += 1
        patch.hbond_donors = donors
        patch.hbond_acceptors = acceptors

    def _find_residues_near(self, center, radius):
        """Find residues within radius of center"""
        found = set()
        residues = []
        
        center_arr = np.array(center)
        
        for atom in self.atoms:
            if np.linalg.norm(atom['coord'] - center_arr) < radius:
                key = (atom['chain'], atom['resn'], atom['resi'])
                if key not in found:
                    found.add(key)
                    residues.append({
                        'chain': atom['chain'],
                        'resn': atom['resn'],
                        'resi': atom['resi']
                    })
        return residues

    def _get_residue_charge(self, resn: str) -> float:
        """Get approximate charge for a residue type"""
        charge_map = {
            'ASP': -1.0, 'GLU': -1.0,  # Negative
            'LYS': 1.0, 'ARG': 1.0,    # Positive
            'HIS': 0.5,                 # pH dependent
        }
        return charge_map.get(resn, 0.0)
    
    def _get_hydrophobicity(self, resn: str) -> float:
        """Get normalized hydrophobicity for a residue type (0-1 scale)"""
        # Kyte-Doolittle scale normalized to 0-1
        hydro_scale = {
            'ILE': 1.0, 'VAL': 0.97, 'LEU': 0.92, 'PHE': 0.81, 'CYS': 0.78,
            'MET': 0.71, 'ALA': 0.70, 'GLY': 0.46, 'THR': 0.42, 'SER': 0.41,
            'TRP': 0.40, 'TYR': 0.36, 'PRO': 0.32, 'HIS': 0.14, 'GLU': 0.11,
            'GLN': 0.11, 'ASP': 0.11, 'ASN': 0.11, 'LYS': 0.07, 'ARG': 0.0
        }
        return hydro_scale.get(resn, 0.5)

    def render_full_surface(self, color_by: str = 'potential',
                           surface_name: str = None,
                           transparency: float = 0.0,
                           color_range: tuple = None,
                           save_files: bool = False,
                           output_dir: str = None) -> bool:
        """
        Render the complete protein surface with continuous color mapping.
        
        For electrostatic potential: Uses PyMOL's native ramp_new method with DX grid
        for smooth, high-quality visualization (same as PyMOL APBS plugin).
        
        For hydrophobicity: Uses B-factor coloring with Kyte-Doolittle scale.
        
        Args:
            color_by: 'potential', 'hydrophobicity', 'curvature', 'hbond', or 'combined'
            surface_name: Name for the surface object (default: obj_name + '_ESP')
            transparency: Surface transparency (0.0 = opaque)
            color_range: (min, max) for color scale, or None for auto
            save_files: If True, save intermediate files (PDB, PQR, DX) to output_dir
            output_dir: Directory to save files (required if save_files=True)
            
        Returns:
            True if successful
        """
        # If save_files is requested, set output_dir
        if save_files and output_dir:
            self.output_dir = output_dir
            os.makedirs(output_dir, exist_ok=True)
            print(f"  Output files will be saved to: {output_dir}")
        if not cmd:
            print("Error: PyMOL not available")
            return False

        cmd.bg_color('white')
        
        # Use different default names for different surface types
        # ESP = Electrostatic Potential
        if surface_name is None:
            if color_by == 'potential':
                surface_name = f"{self.obj_name}_ESP"
            elif color_by == 'hydrophobicity':
                surface_name = f"{self.obj_name}_hydro"
            elif color_by == 'curvature':
                surface_name = f"{self.obj_name}_curvature"
            elif color_by == 'hbond':
                surface_name = f"{self.obj_name}_hbond"
            else:
                surface_name = f"{self.obj_name}_surface"
        
        print(f"Rendering full surface for {self.obj_name}...")
        print(f"  Color by: {color_by}")
        
        # Get atoms if not already loaded
        if not self.atoms:
            self.atoms = self._get_atoms()
        
        if not self.atoms:
            print("  Error: No atoms found")
            return False
        
        if color_by == 'potential':
            # Use PyMOL native method with DX grid for electrostatic potential
            return self._render_potential_surface_native(surface_name, transparency, color_range)
            
        elif color_by == 'hydrophobicity':
            # Use B-factor method for hydrophobicity
            return self._render_hydrophobicity_surface(surface_name, transparency, color_range)

        elif color_by == 'curvature':
            # Use B-factor method for signed curvature / concavity
            return self._render_curvature_surface(surface_name, transparency, color_range)

        elif color_by == 'hbond':
            # Use B-factor method for approximate donor/acceptor surface roles
            return self._render_hbond_surface(surface_name, transparency, color_range)
            
        elif color_by == 'combined':
            # Combined mode
            return self._render_combined_surface_v2(surface_name, transparency, color_range)
        else:
            print(f"  Unknown color mode: {color_by}")
            return False

    def _interpolate_rgb(self, colors: List[Tuple[float, float, float]], t: float) -> Tuple[float, float, float]:
        """Interpolate across a list of RGB color stops."""
        if not colors:
            return (1.0, 1.0, 1.0)
        if len(colors) == 1:
            return colors[0]
        t = max(0.0, min(1.0, float(t)))
        pos = t * (len(colors) - 1)
        idx = min(int(math.floor(pos)), len(colors) - 2)
        frac = pos - idx
        c0 = np.array(colors[idx], dtype=float)
        c1 = np.array(colors[idx + 1], dtype=float)
        c = c0 * (1.0 - frac) + c1 * frac
        return tuple(float(x) for x in c)

    def _create_surface_colorbar(self, group_name: str,
                                 colors: List[Tuple[float, float, float]],
                                 value_min: float,
                                 value_max: float,
                                 title: str,
                                 value_format: str = "{:.2f}") -> Optional[str]:
        """Create a compact 3D colorbar next to the rendered molecular surface."""
        if not cmd:
            return None

        try:
            from pymol.cgo import COLOR, CYLINDER
        except Exception:
            return None

        try:
            cmd.delete(f"{group_name}_colorbar*")
        except Exception:
            pass

        try:
            min_corner, max_corner = cmd.get_extent(self._selection_expr())
        except Exception:
            min_corner, max_corner = ([0.0, 0.0, 0.0], [20.0, 20.0, 20.0])

        min_corner = np.asarray(min_corner, dtype=float)
        max_corner = np.asarray(max_corner, dtype=float)
        span = np.maximum(max_corner - min_corner, 1.0)

        height = max(float(span[2]) * 0.45, 10.0)
        radius = max(height * 0.018, 0.12)
        x = float(max_corner[0] + max(span[0] * 0.12, 4.0))
        y = float(max_corner[1] + max(span[1] * 0.05, 1.5))
        z0 = float(min_corner[2] + span[2] * 0.25)
        z1 = z0 + height

        n_segments = 28
        cgo_obj = []
        for i in range(n_segments):
            t0 = i / n_segments
            t1 = (i + 1) / n_segments
            c0 = self._interpolate_rgb(colors, t0)
            c1 = self._interpolate_rgb(colors, t1)
            za = z0 + height * t0
            zb = z0 + height * t1
            cgo_obj.extend([
                CYLINDER,
                x, y, za,
                x, y, zb,
                radius,
                c0[0], c0[1], c0[2],
                c1[0], c1[1], c1[2],
            ])

        bar_obj = f"{group_name}_colorbar_bar"
        cmd.load_cgo(cgo_obj, bar_obj)

        label_specs = [
            (f"{group_name}_colorbar_max", value_format.format(value_max), z1),
            (f"{group_name}_colorbar_mid", value_format.format((value_min + value_max) / 2.0), (z0 + z1) / 2.0),
            (f"{group_name}_colorbar_min", value_format.format(value_min), z0),
            (f"{group_name}_colorbar_title", title, z1 + height * 0.10),
        ]
        label_objs = []
        for obj_name, text, z in label_specs:
            cmd.pseudoatom(obj_name, pos=[x + radius * 3.5, y, z])
            cmd.label(obj_name, f'"{text}"')
            cmd.set('label_color', 'black', obj_name)
            cmd.set('label_size', 14, obj_name)
            cmd.hide('nonbonded', obj_name)
            label_objs.append(obj_name)

        colorbar_group = f"{group_name}_colorbar"
        cmd.group(colorbar_group, " ".join([bar_obj] + label_objs))
        return colorbar_group
    
    def _render_potential_surface_native(self, surface_name: str,
                                          transparency: float = 0.0,
                                          color_range: tuple = None) -> bool:
        """
        Render electrostatic potential surface using PyMOL native method.
        
        This uses ramp_new to directly map DX grid values to surface colors,
        exactly matching PyMOL's APBS Electrostatics plugin behavior.
        
        Reference: https://pymolwiki.org/index.php/Ramp_new
        Example from PyMOL Wiki:
            load 1ubq_apbs.dx, e_pot_map
            fetch 1ubq, async=0
            as surface
            ramp_new e_pot_color, e_pot_map, [-5, 0, 5], [red, white, blue]
            set surface_color, e_pot_color
            set surface_ramp_above_mode
        """
        # Run APBS if enabled and not already done
        if self.use_apbs and self.dx_data is None:
            print("  Running APBS for accurate electrostatics...")
            self.dx_data = self._run_apbs_pipeline()
        
        if self.dx_data is None:
            print("  ⚠️ APBS not available, falling back to residue-based method")
            return self._render_potential_surface_fallback(surface_name, transparency, color_range)
        
        # Get the DX file path from temp directory or output directory
        dx_file = None
        import glob
        
        # Check output_dir first (if save_files was requested)
        if self.output_dir and os.path.exists(self.output_dir):
            dx_files = glob.glob(os.path.join(self.output_dir, '*.dx'))
            if dx_files:
                dx_file = dx_files[0]
        
        # Fall back to temp directory
        if not dx_file and self._temp_dir and os.path.exists(self._temp_dir):
            dx_files = glob.glob(os.path.join(self._temp_dir, '*.dx'))
            if dx_files:
                dx_file = dx_files[0]
        
        if not dx_file or not os.path.exists(dx_file):
            print("  ⚠️ DX file not found, falling back to residue-based method")
            return self._render_potential_surface_fallback(surface_name, transparency, color_range)
        
        # Set color range - PyMOL APBS plugin default is ±5 kT/e
        # This is the standard range for electrostatic visualization
        if color_range:
            vmin = color_range[0]
            vmax = color_range[1]
        else:
            # Use standard ±5 kT/e range (same as PyMOL APBS plugin default)
            vmin, vmax = -5.0, 5.0
        
        print(f"  Loading DX file: {os.path.basename(dx_file)}")
        print(f"  Color scale: [{vmin:.1f}, {vmax:.1f}] kT/e")
        
        # === PyMOL Wiki official method ===
        # Use surface_name as the group name (e.g., "8D7Z_ESP")
        group_name = surface_name
        
        # Step 1: Load DX file as a map object
        map_name = f"{group_name}_map"
        cmd.load(dx_file, map_name)
        
        # Step 2: Create a copy of the object for the surface
        # This avoids modifying the original object's surface display
        surface_obj = f"{group_name}_surf"
        cmd.create(surface_obj, self._selection_expr())
        
        # Step 3: Hide everything on the surface object, then show only surface
        cmd.hide('everything', surface_obj)
        cmd.show('surface', surface_obj)
        
        # Step 4: Create color ramp using PyMOL's ramp_new
        # PyMOL Wiki example: ramp_new e_pot_color, e_pot_map, [-5, 0, 5], [red, white, blue]
        ramp_name = f"{group_name}_ramp"
        cmd.ramp_new(ramp_name, map_name, [vmin, 0.0, vmax], ['red', 'white', 'blue'])
        
        # Step 5: Apply the ramp to color the surface object (not the original)
        # PyMOL Wiki: set surface_color, e_pot_color
        cmd.set('surface_color', ramp_name, surface_obj)
        
        # Step 6: Enable surface_ramp_above_mode (critical for proper coloring)
        # PyMOL Wiki: set surface_ramp_above_mode
        cmd.set('surface_ramp_above_mode', 1)
        colorbar_group = self._create_surface_colorbar(
            group_name,
            colors=[(1.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0.0, 0.0, 1.0)],
            value_min=vmin,
            value_max=vmax,
            title="Electrostatic kT/e",
            value_format="{:.1f}",
        )
        
        # === Advanced rendering settings ===
        # Depth and lighting
        cmd.set('depth_cue', 0)
        cmd.set('ambient', 0.3)
        cmd.set('direct', 0.7)
        cmd.set('light_count', 3)
        
        # Cartoon style for protein backbone
        cmd.set_color('color1', [220, 223, 228])
        cmd.set_color('color2', [144, 148, 159])
        cmd.set('cartoon_color', 'color1', self.obj_name)
        cmd.set('cartoon_highlight_color', 'color2', self.obj_name)
        cmd.set('cartoon_oval_length', 1, self.obj_name)
        cmd.set('cartoon_oval_width', 0.2, self.obj_name)
        cmd.set('cartoon_transparency', 0.3, self.obj_name)
        
        # Shadow and specular highlights
        cmd.set('ray_shadow', 'off')
        cmd.set('specular', 'on')
        cmd.set('specular_intensity', 0.6)
        cmd.set('shininess', 80)
        cmd.set('reflect', 0.25)

        
        # Set surface quality for smooth rendering
        cmd.set('surface_quality', 1, surface_obj)
        
        # Set transparency if requested
        if transparency > 0:
            cmd.set('transparency', transparency, surface_obj)
        
        # Hide surface on original object (keep cartoon visible with normal colors)
        cmd.hide('surface', self.obj_name)
        
        # Step 7: Group all ESP-related objects together
        # This keeps the PyMOL object list clean and organized
        group_members = f"{map_name} {ramp_name} {surface_obj}"
        if colorbar_group:
            group_members += f" {colorbar_group}"
        cmd.group(group_name, group_members)
        
        # Don't cleanup temp dir yet - we need the DX file for the ramp
        # The map will be deleted when PyMOL session ends
        
        print(f"✅ Electrostatic surface rendered as '{group_name}'")
        print("   Method: PyMOL native ramp_new (same as APBS plugin)")
        print("   Color scheme: Red (negative) → White (neutral) → Blue (positive)")
        print(f"   Group contains: {surface_obj}, {map_name}, {ramp_name}")
        if colorbar_group:
            print(f"   Colorbar: {colorbar_group}")
        
        # Export per-residue electrostatic potential to CSV if output_dir is set
        if self.output_dir:
            self._export_electrostatic_csv()
        
        # Report saved files if output_dir was specified
        work_dir = self.output_dir or self._temp_dir
        if work_dir and os.path.exists(work_dir):
            saved_files = []
            for f in os.listdir(work_dir):
                fpath = os.path.join(work_dir, f)
                if os.path.isfile(fpath):
                    saved_files.append(f)
            if saved_files and self.output_dir:
                print(f"   Output files saved to: {self.output_dir}")
                for f in saved_files:
                    print(f"     • {f}")
        
        return True
    
    def _export_electrostatic_csv(self) -> str:
        """
        Export per-residue electrostatic potential values to CSV file.
        
        For each residue, samples the electrostatic potential at multiple surface-exposed
        atom positions and takes the average to get a more representative value.
        
        Returns:
            Path to the exported CSV file
        """
        import csv
        
        if not self.output_dir:
            return None
        
        if self.dx_data is None:
            print("   ⚠️ No APBS data available for electrostatic CSV export")
            return None
        
        csv_path = os.path.join(self.output_dir, 'electrostatic_potential.csv')
        
        # Collect unique residues with their atoms
        residue_atoms = {}
        for atom in self.atoms:
            key = (atom['chain'], atom['resi'], atom['resn'])
            if key not in residue_atoms:
                residue_atoms[key] = []
            residue_atoms[key].append(atom)
        
        # For each residue, sample at multiple surface-exposed positions
        residue_data = {}
        
        # Surface-exposed atoms to sample (in order of preference)
        surface_atoms = ['CB', 'CG', 'CD', 'CE', 'NZ', 'OD1', 'OD2', 'OE1', 'OE2',
                         'ND1', 'ND2', 'NE', 'NH1', 'NH2', 'OH', 'SG', 'SD']
        
        for key, atoms in residue_atoms.items():
            chain, resi, resn = key
            
            # Collect sample coordinates from surface-exposed atoms
            sample_coords = []
            
            # First try surface-exposed side chain atoms
            for atom in atoms:
                if atom['name'] in surface_atoms:
                    sample_coords.append(atom['coord'])
            
            # If no surface atoms found, use Cα
            if not sample_coords:
                for atom in atoms:
                    if atom['name'] == 'CA':
                        sample_coords.append(atom['coord'])
                        break
            
            # If still no coords, use centroid
            if not sample_coords:
                coords = np.array([a['coord'] for a in atoms])
                sample_coords.append(coords.mean(axis=0))
            
            # Sample electrostatic potential at all positions
            sample_coords = np.array(sample_coords)
            potentials = self._sample_dx_at_points(sample_coords)
            
            # Filter out extreme values (likely grid boundary artifacts)
            # Valid APBS values are typically in range [-10, +10] kT/e
            valid_potentials = potentials[(potentials > -10) & (potentials < 10)]
            
            if len(valid_potentials) > 0:
                # Use median to be robust against outliers
                avg_potential = np.median(valid_potentials)
            else:
                # If all values are extreme, use the one closest to zero
                avg_potential = potentials[np.argmin(np.abs(potentials))]
                # Clamp to reasonable range
                avg_potential = np.clip(avg_potential, -10.0, 10.0)
            
            residue_data[key] = {
                'chain': chain,
                'resi': resi,
                'resn': resn,
                'potential': avg_potential,
                'n_samples': len(sample_coords)
            }
        
        # Sort by chain and residue number
        sorted_residues = sorted(residue_data.values(),
                                  key=lambda x: (x['chain'], int(x['resi']) if x['resi'].isdigit() else 0))
        
        # Write CSV
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            fieldnames = ['Chain', 'Residue_Number', 'Residue_Name', 'Electrostatic_Potential_kTe']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for res in sorted_residues:
                writer.writerow({
                    'Chain': res['chain'],
                    'Residue_Number': res['resi'],
                    'Residue_Name': res['resn'],
                    'Electrostatic_Potential_kTe': f"{res['potential']:.4f}"
                })
        
        print(f"   Electrostatic CSV saved: {os.path.basename(csv_path)}")
        print(f"   Total residues: {len(sorted_residues)}")
        
        # Report statistics
        all_potentials = [r['potential'] for r in sorted_residues]
        print(f"   Potential range: [{min(all_potentials):.2f}, {max(all_potentials):.2f}] kT/e")
        
        return csv_path

    def _ensure_surface_cache(self) -> bool:
        """Build surface points and mapped properties if they are not already cached."""
        if self._surface_cache is not None:
            return True

        if not self.atoms:
            self.atoms = self._get_atoms()
        if not self.atoms:
            return False

        pd = PocketDetector(grid_spacing=self.grid_spacing)
        self.grid_info, self.origin, atom_coords = pd._build_grid(self.atoms)
        occupied = pd._mark_occupied(self.grid_info, atom_coords, self.origin)
        struct = ndimage.generate_binary_structure(3, 1)
        dilated = ndimage.binary_dilation(occupied, structure=struct)
        surface_mask = dilated & (~occupied)
        surface_indices = np.argwhere(surface_mask)
        if len(surface_indices) == 0:
            return False

        surface_coords = self.origin + surface_indices * self.grid_spacing
        try:
            from scipy.spatial import cKDTree
            atom_kdtree = cKDTree(atom_coords)
        except ImportError:
            atom_kdtree = None

        self._map_properties_to_surface(surface_coords, surface_indices, atom_kdtree, self.atoms)
        return self._surface_cache is not None
    
    def _render_potential_surface_fallback(self, surface_name: str,
                                            transparency: float = 0.0,
                                            color_range: tuple = None) -> bool:
        """
        Fallback method for electrostatic potential using residue-based approximation.
        Used when APBS is not available.
        """
        atom_coords = np.array([a['coord'] for a in self.atoms])
        values = np.array([self._get_residue_charge(a['resn']) for a in self.atoms])
        
        vmin = color_range[0] if color_range else -1.0
        vmax = color_range[1] if color_range else 1.0
        
        print(f"  Using residue-based approximation")
        print(f"  Charge range: [{values.min():.2f}, {values.max():.2f}]")
        
        # Create surface and color by B-factor
        cmd.create(surface_name, self._selection_expr())
        
        model = cmd.get_model(surface_name)
        for i, atom in enumerate(model.atom):
            if i < len(values):
                atom.b = values[i]
        cmd.load_model(model, surface_name, state=1)
        
        cmd.hide('everything', surface_name)
        cmd.show('surface', surface_name)
        cmd.spectrum('b', 'red_white_blue', surface_name, minimum=vmin, maximum=vmax)
        self._create_surface_colorbar(
            surface_name,
            colors=[(1.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0.0, 0.0, 1.0)],
            value_min=vmin,
            value_max=vmax,
            title="Charge",
            value_format="{:.1f}",
        )
        
        if transparency > 0:
            cmd.set('transparency', transparency, surface_name)
        
        cmd.set('surface_quality', 1, surface_name)
        
        # Hide surface on original object (keep cartoon visible with normal colors)
        cmd.hide('surface', self.obj_name)
        
        self._cleanup()
        
        print(f"✅ Electrostatic surface rendered as '{surface_name}'")
        print("   Method: Residue-based approximation (APBS not available)")
        print("   Color scheme: Red (negative) → White (neutral) → Blue (positive)")
        
        return True

    def _render_curvature_surface(self, surface_name: str,
                                  transparency: float = 0.0,
                                  color_range: tuple = None) -> bool:
        """
        Render signed surface curvature using B-factor coloring.

        Negative values are concave/cup-like regions, positive values are convex/cap-like
        regions. This is intended to expose grooves, shallow pockets, and curved surface
        patches alongside electrostatic and hydrophobic surface analysis.
        """
        if not self._ensure_surface_cache():
            print("  Error: failed to calculate surface curvature")
            return False

        cache = self._surface_cache
        props = cache['properties']
        nearest_atom_indices = cache['nearest_atom_indices']
        signed = props['signed_curvature']

        atom_values = np.zeros(len(self.atoms), dtype=float)
        atom_counts = np.zeros(len(self.atoms), dtype=int)
        for value, atom_idx in zip(signed, nearest_atom_indices):
            atom_values[atom_idx] += value
            atom_counts[atom_idx] += 1
        has_values = atom_counts > 0
        atom_values[has_values] /= atom_counts[has_values]

        if color_range:
            vmin, vmax = color_range
        else:
            populated = atom_values[has_values]
            if len(populated):
                max_abs = float(np.percentile(np.abs(populated), 95))
                max_abs = max(max_abs, 0.02)
            else:
                max_abs = 0.05
            vmin, vmax = -max_abs, max_abs

        print(f"  Signed curvature range: [{atom_values.min():.4f}, {atom_values.max():.4f}]")
        print(f"  Color scale: [{vmin:.4f}, {vmax:.4f}]")

        cmd.create(surface_name, self._selection_expr())

        model = cmd.get_model(surface_name)
        for i, atom in enumerate(model.atom):
            if i < len(atom_values):
                atom.b = float(atom_values[i])
        cmd.load_model(model, surface_name, state=1)

        cmd.hide('everything', surface_name)
        cmd.show('surface', surface_name)
        cmd.set_color('glint_curv_concave', [0.50, 0.18, 0.78])
        cmd.set_color('glint_curv_flat', [0.88, 0.88, 0.84])
        cmd.set_color('glint_curv_convex', [0.95, 0.66, 0.14])
        try:
            cmd.spectrum(
                'b',
                'glint_curv_concave glint_curv_flat glint_curv_convex',
                surface_name,
                minimum=vmin,
                maximum=vmax
            )
        except Exception:
            cmd.spectrum('b', 'purple_white_orange', surface_name, minimum=vmin, maximum=vmax)
        self._create_surface_colorbar(
            surface_name,
            colors=[(0.50, 0.18, 0.78), (0.88, 0.88, 0.84), (0.95, 0.66, 0.14)],
            value_min=vmin,
            value_max=vmax,
            title="Signed Curvature",
            value_format="{:.3f}",
        )

        if transparency > 0:
            cmd.set('transparency', transparency, surface_name)

        cmd.set('surface_quality', 1, surface_name)
        cmd.hide('surface', self.obj_name)

        if self.output_dir:
            self._export_curvature_csv(atom_values)

        self._cleanup()

        print(f"✅ Curvature surface rendered as '{surface_name}'")
        print("   Method: local surface-neighborhood curvature")
        print("   Color scheme: Purple (concave) → light gray (flat) → gold (convex)")

        return True

    def _export_curvature_csv(self, atom_values: np.ndarray) -> str:
        """Export per-residue signed curvature and concavity values to CSV."""
        import csv

        if not self.output_dir:
            return None

        csv_path = os.path.join(self.output_dir, 'surface_curvature.csv')
        residue_values = {}
        for atom, value in zip(self.atoms, atom_values):
            key = (atom['chain'], atom['resi'], atom['resn'])
            residue_values.setdefault(key, []).append(float(value))

        rows = []
        for (chain, resi, resn), values in residue_values.items():
            vals = np.array(values, dtype=float)
            mean_signed = float(np.mean(vals))
            rows.append({
                'chain': chain,
                'resi': resi,
                'resn': resn,
                'signed_curvature': mean_signed,
                'curvature': float(np.mean(np.abs(vals))),
                'concavity': max(-mean_signed, 0.0),
            })

        rows = sorted(rows, key=lambda x: (x['chain'], int(x['resi']) if str(x['resi']).isdigit() else 0))

        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            fieldnames = [
                'Chain', 'Residue_Number', 'Residue_Name',
                'Signed_Curvature', 'Curvature_Magnitude', 'Concavity'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow({
                    'Chain': row['chain'],
                    'Residue_Number': row['resi'],
                    'Residue_Name': row['resn'],
                    'Signed_Curvature': f"{row['signed_curvature']:.5f}",
                    'Curvature_Magnitude': f"{row['curvature']:.5f}",
                    'Concavity': f"{row['concavity']:.5f}",
                })

        print(f"   Curvature CSV saved: {os.path.basename(csv_path)}")
        print(f"   Total residues: {len(rows)}")

        return csv_path

    def _render_hbond_surface(self, surface_name: str,
                              transparency: float = 0.0,
                              color_range: tuple = None) -> bool:
        """Render approximate hydrogen-bond donor/acceptor distribution."""
        if not self._ensure_surface_cache():
            print("  Error: failed to calculate donor/acceptor surface roles")
            return False

        cache = self._surface_cache
        role_values = cache['properties']['hbond_role']
        nearest_atom_indices = cache['nearest_atom_indices']

        atom_values = np.zeros(len(self.atoms), dtype=float)
        atom_counts = np.zeros(len(self.atoms), dtype=int)
        for value, atom_idx in zip(role_values, nearest_atom_indices):
            atom_values[atom_idx] += value
            atom_counts[atom_idx] += 1
        has_values = atom_counts > 0
        atom_values[has_values] /= atom_counts[has_values]

        vmin, vmax = color_range if color_range else (-1.0, 1.0)
        print(f"  Donor/acceptor role range: [{atom_values.min():.2f}, {atom_values.max():.2f}]")

        cmd.create(surface_name, self._selection_expr())
        model = cmd.get_model(surface_name)
        for i, atom in enumerate(model.atom):
            if i < len(atom_values):
                atom.b = float(atom_values[i])
        cmd.load_model(model, surface_name, state=1)

        cmd.hide('everything', surface_name)
        cmd.show('surface', surface_name)
        cmd.set_color('glint_hbond_acceptor', [0.78, 0.16, 0.54])
        cmd.set_color('glint_hbond_neutral', [0.88, 0.88, 0.84])
        cmd.set_color('glint_hbond_donor', [0.00, 0.62, 0.72])
        try:
            cmd.spectrum(
                'b',
                'glint_hbond_acceptor glint_hbond_neutral glint_hbond_donor',
                surface_name,
                minimum=vmin,
                maximum=vmax
            )
        except Exception:
            cmd.spectrum('b', 'magenta_white_cyan', surface_name, minimum=vmin, maximum=vmax)
        self._create_surface_colorbar(
            surface_name,
            colors=[(0.78, 0.16, 0.54), (0.88, 0.88, 0.84), (0.00, 0.62, 0.72)],
            value_min=vmin,
            value_max=vmax,
            title="Acceptor / Donor",
            value_format="{:.1f}",
        )

        if transparency > 0:
            cmd.set('transparency', transparency, surface_name)
        cmd.set('surface_quality', 1, surface_name)
        cmd.hide('surface', self.obj_name)

        if self.output_dir:
            self._export_hbond_csv(atom_values)

        self._cleanup()
        print(f"✅ Donor/acceptor surface rendered as '{surface_name}'")
        print("   Color scheme: magenta (acceptor) → light gray (neutral) → cyan (donor)")
        return True

    def _export_hbond_csv(self, atom_values: np.ndarray) -> str:
        """Export per-residue approximate donor/acceptor role values."""
        import csv

        if not self.output_dir:
            return None

        csv_path = os.path.join(self.output_dir, 'hbond_donor_acceptor.csv')
        residue_data = {}
        for atom, value in zip(self.atoms, atom_values):
            key = (atom['chain'], atom['resi'], atom['resn'])
            if key not in residue_data:
                residue_data[key] = {'values': [], 'donors': 0, 'acceptors': 0}
            residue_data[key]['values'].append(float(value))
            role = self._get_hbond_role(atom)
            if role > 0:
                residue_data[key]['donors'] += 1
            elif role < 0:
                residue_data[key]['acceptors'] += 1

        rows = []
        for (chain, resi, resn), data in residue_data.items():
            rows.append({
                'chain': chain,
                'resi': resi,
                'resn': resn,
                'role': float(np.mean(data['values'])),
                'donors': data['donors'],
                'acceptors': data['acceptors'],
            })
        rows = sorted(rows, key=lambda x: (x['chain'], int(x['resi']) if str(x['resi']).isdigit() else 0))

        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            fieldnames = [
                'Chain', 'Residue_Number', 'Residue_Name',
                'Donor_Acceptor_Role', 'Approx_Donor_Atoms', 'Approx_Acceptor_Atoms'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow({
                    'Chain': row['chain'],
                    'Residue_Number': row['resi'],
                    'Residue_Name': row['resn'],
                    'Donor_Acceptor_Role': f"{row['role']:.4f}",
                    'Approx_Donor_Atoms': row['donors'],
                    'Approx_Acceptor_Atoms': row['acceptors'],
                })

        print(f"   Donor/acceptor CSV saved: {os.path.basename(csv_path)}")
        print(f"   Total residues: {len(rows)}")
        return csv_path
    
    def _render_hydrophobicity_surface(self, surface_name: str,
                                        transparency: float = 0.0,
                                        color_range: tuple = None) -> bool:
        """
        Render hydrophobicity surface using B-factor coloring.
        Also exports per-residue hydrophobicity values to CSV if output_dir is set.
        """
        values = np.array([self._get_hydrophobicity(a['resn']) for a in self.atoms])
        
        vmin = color_range[0] if color_range else 0.0
        vmax = color_range[1] if color_range else 1.0
        
        print(f"  Hydrophobicity range: [{values.min():.2f}, {values.max():.2f}]")
        
        # Create surface and color by B-factor
        cmd.create(surface_name, self._selection_expr())
        
        model = cmd.get_model(surface_name)
        for i, atom in enumerate(model.atom):
            if i < len(values):
                atom.b = values[i]
        cmd.load_model(model, surface_name, state=1)
        
        cmd.hide('everything', surface_name)
        cmd.show('surface', surface_name)
        cmd.spectrum('b', 'white_green', surface_name, minimum=vmin, maximum=vmax)
        self._create_surface_colorbar(
            surface_name,
            colors=[(1.0, 1.0, 1.0), (0.0, 0.55, 0.15)],
            value_min=vmin,
            value_max=vmax,
            title="Hydrophobicity",
            value_format="{:.2f}",
        )
        
        if transparency > 0:
            cmd.set('transparency', transparency, surface_name)
        
        cmd.set('surface_quality', 1, surface_name)
        
        # Hide surface on original object (keep cartoon visible with normal colors)
        cmd.hide('surface', self.obj_name)
        
        # Export per-residue hydrophobicity to CSV if output_dir is set
        if self.output_dir:
            self._export_hydrophobicity_csv()
        
        self._cleanup()
        
        print(f"✅ Hydrophobicity surface rendered as '{surface_name}'")
        print("   Method: Kyte-Doolittle scale")
        print("   Color scheme: White (hydrophilic) → Green (hydrophobic)")
        
        return True
    
    def _export_hydrophobicity_csv(self) -> str:
        """
        Export per-residue hydrophobicity values to CSV file.
        
        Returns:
            Path to the exported CSV file
        """
        import csv
        
        if not self.output_dir:
            return None
        
        csv_path = os.path.join(self.output_dir, 'hydrophobicity.csv')
        
        # Collect unique residues with their hydrophobicity values
        residue_data = {}
        for atom in self.atoms:
            key = (atom['chain'], atom['resi'], atom['resn'])
            if key not in residue_data:
                hydro = self._get_hydrophobicity(atom['resn'])
                residue_data[key] = {
                    'chain': atom['chain'],
                    'resi': atom['resi'],
                    'resn': atom['resn'],
                    'hydrophobicity': hydro
                }
        
        # Sort by chain and residue number
        sorted_residues = sorted(residue_data.values(),
                                  key=lambda x: (x['chain'], int(x['resi']) if x['resi'].isdigit() else 0))
        
        # Write CSV
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            fieldnames = ['Chain', 'Residue_Number', 'Residue_Name', 'Hydrophobicity_Normalized']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for res in sorted_residues:
                writer.writerow({
                    'Chain': res['chain'],
                    'Residue_Number': res['resi'],
                    'Residue_Name': res['resn'],
                    'Hydrophobicity_Normalized': f"{res['hydrophobicity']:.4f}"
                })
        
        print(f"   Hydrophobicity CSV saved: {os.path.basename(csv_path)}")
        print(f"   Total residues: {len(sorted_residues)}")
        
        return csv_path
    
    def _render_combined_surface_v2(self, surface_name: str,
                                     transparency: float = 0.0,
                                     color_range: tuple = None) -> bool:
        """
        Render combined surface showing both electrostatic and hydrophobicity.
        Creates two semi-transparent surfaces.
        """
        print("  Creating combined visualization (two surfaces)...")
        
        # Create electrostatic surface
        esp_name = f"{surface_name}_esp"
        self._render_potential_surface_native(esp_name, transparency=0.5, color_range=color_range)
        
        # Create hydrophobicity surface
        hydro_name = f"{surface_name}_hydro"
        self._render_hydrophobicity_surface(hydro_name, transparency=0.5, color_range=None)
        
        # Group them
        cmd.group(surface_name, f"{esp_name} {hydro_name}")
        
        print(f"✅ Combined surface rendered as '{surface_name}'")
        print("   Contains: {esp_name} (electrostatic) + {hydro_name} (hydrophobicity)")
        
        return True
    
    def _render_combined_surface(self, surface_name: str, atom_coords: np.ndarray,
                                  potentials: np.ndarray, hydrophobicities: np.ndarray,
                                  transparency: float = 0.0) -> bool:
        """
        Render surface with combined electrostatic + hydrophobicity coloring.
        
        Color scheme:
        - Hydrophobic regions (hydrophobicity > 0.6): Green shades
        - Positive charged (potential > threshold): Blue shades
        - Negative charged (potential < -threshold): Red shades
        - Neutral hydrophilic: White
        """
        # Determine thresholds
        if self.dx_data is not None:
            pot_threshold = 1.0  # kT/e for APBS
        else:
            pot_threshold = 0.5  # For residue-based
        
        hydro_threshold = 0.6  # Kyte-Doolittle normalized
        
        # Create a copy of the object
        cmd.create(surface_name, self._selection_expr())
        
        # Calculate combined color values
        # Priority: charged > hydrophobic > neutral
        # Encode as: negative=-10 to -1, neutral=0, positive=1 to 10, hydrophobic=20 to 30
        combined_values = np.zeros(len(potentials))
        
        for i in range(len(potentials)):
            pot = potentials[i]
            hydro = hydrophobicities[i]
            
            if pot > pot_threshold:
                # Positive charged - map to 1-10 (will be blue)
                combined_values[i] = min(pot / pot_threshold * 5, 10)
            elif pot < -pot_threshold:
                # Negative charged - map to -10 to -1 (will be red)
                combined_values[i] = max(pot / pot_threshold * 5, -10)
            elif hydro > hydro_threshold:
                # Hydrophobic - map to 20-30 (will be green)
                combined_values[i] = 20 + (hydro - hydro_threshold) / (1 - hydro_threshold) * 10
            else:
                # Neutral hydrophilic - 0 (will be white)
                combined_values[i] = 0
        
        # Set B-factors
        model = cmd.get_model(surface_name)
        for i, atom in enumerate(model.atom):
            if i < len(combined_values):
                atom.b = combined_values[i]
        cmd.load_model(model, surface_name, state=1)
        
        # Show surface
        cmd.hide('everything', surface_name)
        cmd.show('surface', surface_name)
        
        # Apply custom coloring using selections
        # Negative (red)
        cmd.select('_neg_charged', f'{surface_name} and b < -0.5')
        cmd.color('red', '_neg_charged')
        
        # Positive (blue)
        cmd.select('_pos_charged', f'{surface_name} and b > 0.5 and b < 15')
        cmd.color('blue', '_pos_charged')
        
        # Hydrophobic (green)
        cmd.select('_hydrophobic', f'{surface_name} and b > 15')
        cmd.color('green', '_hydrophobic')
        
        # Neutral (white)
        cmd.select('_neutral', f'{surface_name} and b > -0.5 and b < 0.5')
        cmd.color('white', '_neutral')
        
        # Clean up selections
        cmd.delete('_neg_charged')
        cmd.delete('_pos_charged')
        cmd.delete('_hydrophobic')
        cmd.delete('_neutral')
        
        # Set transparency
        if transparency > 0:
            cmd.set('transparency', transparency, surface_name)
        
        # Set surface quality
        cmd.set('surface_quality', 1, surface_name)
        
        # Hide surface on original object (keep cartoon visible with normal colors)
        cmd.hide('surface', self.obj_name)
        
        # Cleanup
        self._cleanup()
        
        print(f"✅ Combined surface rendered as '{surface_name}'")
        print("   Color scheme: Green (疏水) + Red (负电) + Blue (正电) + White (中性亲水)")
        
        return True
    
def render_protein_surface(obj_name: str, use_apbs: bool = True, ph: float = 7.4,
                          color_by: str = 'potential', transparency: float = 0.0,
                          color_min: float = None, color_max: float = None,
                          save_files: bool = False, output_dir: str = None,
                          selection: str = "all") -> bool:
    """
    Render complete protein surface with APBS electrostatics using PyMOL native surface.
    
    This creates a smooth, high-quality molecular surface colored by electrostatic
    potential or hydrophobicity.
    
    Args:
        obj_name: PyMOL object name
        use_apbs: Use APBS for accurate electrostatics (default: True)
        ph: pH for PDB2PQR protonation (default: 7.4)
        color_by: 'potential', 'hydrophobicity', 'curvature', or 'hbond'
        transparency: Surface transparency (0.0 = opaque)
        color_min: Minimum value for color scale (default: auto)
        color_max: Maximum value for color scale (default: auto)
        save_files: If True, save intermediate files (PDB, PQR, DX) to output_dir
        output_dir: Directory to save files (required if save_files=True)
        selection: PyMOL selection inside obj_name to analyze/render
        
    Returns:
        True if successful
        
    Usage in PyMOL:
        render_protein_surface '8D7Z', use_apbs=True, ph=7.4
        render_protein_surface '8D7Z', color_by='hydrophobicity'
        render_protein_surface '8D7Z', color_by='curvature'
        render_protein_surface '8D7Z', color_by='hbond'
        render_protein_surface '8D7Z', save_files=True, output_dir='/path/to/output'
    """
    color_range = None
    if color_min is not None and color_max is not None:
        color_range = (color_min, color_max)
    
    analyzer = SurfaceAnalyzer(obj_name, grid_spacing=0.5,
                               use_apbs=use_apbs, ph=ph,
                               selection=selection)
    return analyzer.render_full_surface(color_by=color_by,
                                        transparency=transparency,
                                        color_range=color_range,
                                        save_files=save_files,
                                        output_dir=output_dir)


# Register PyMOL command
if cmd:
    cmd.extend('render_protein_surface', render_protein_surface)
