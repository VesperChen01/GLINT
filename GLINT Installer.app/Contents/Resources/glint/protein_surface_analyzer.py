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
        self.type = patch_type  # 'electrostatic_pos', 'electrostatic_neg', 'hydrophobic'
        self.points: List[Tuple[float, float, float]] = []  # Surface points (x,y,z)
        self.residues: List[Dict] = []  # Associated residues
        self.area: float = 0.0
        self.center: Tuple[float, float, float] = (0.0, 0.0, 0.0)
        
        # Properties
        self.avg_potential: float = 0.0  # For electrostatic patches
        self.avg_hydrophobicity: float = 0.0  # For hydrophobic patches
        self.curvature: float = 0.0
        self.shape_index: float = 0.0  # -1 (cup) to +1 (cap)
        self.score: float = 0.0  # Druggability/Importance score

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
            'shape_index': float(self.shape_index),
            'score': float(self.score),
            'num_residues': len(self.residues),
            'residues': ";".join([f"{r['chain']}:{r['resn']}:{r['resi']}" for r in self.residues])
        }

class SurfaceAnalyzer:
    """Analyzes protein surface for features"""
    
    def __init__(self, obj_name: str, grid_spacing: float = 0.5, 
                 use_apbs: bool = True, ph: float = 7.4,
                 output_dir: str = None):
        """
        Initialize surface analyzer.
        
        Args:
            obj_name: PyMOL object name
            grid_spacing: Grid spacing for surface generation (Å)
            use_apbs: Use APBS/PDB2PQR for accurate electrostatics (default: True)
            ph: pH for PDB2PQR protonation state (default: 7.4)
            output_dir: Directory for temporary files (auto-created if None)
        """
        self.obj_name = obj_name
        self.grid_spacing = grid_spacing
        self.use_apbs = use_apbs and APBS_AVAILABLE
        self.ph = ph
        self.output_dir = output_dir
        self.atoms = []
        self.grid = None
        self.origin = None
        self.dx_data = None  # APBS electrostatic potential grid
        self._temp_dir = None
        
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
            
        properties = self._map_properties_to_surface(surface_coords, atom_kdtree, self.atoms)
        
        # 5. Detect Patches
        patches = []
        patches.extend(self._detect_electrostatic_patches(surface_points_indices, properties, surface_coords))
        patches.extend(self._detect_hydrophobic_patches(surface_points_indices, properties, surface_coords))
        
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
            cmd.save(protein_pdb, self.obj_name)
            
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
            model = cmd.get_model(self.obj_name)
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

    def _map_properties_to_surface(self, surface_coords, kdtree, atoms):
        """
        Map electrostatic potential and hydrophobicity to surface points.
        
        If APBS data is available, uses accurate electrostatic potential from DX file.
        Otherwise, falls back to residue-based approximation.
        
        Returns dict with 'potential' and 'hydrophobicity' arrays.
        """
        n_points = len(surface_coords)
        potentials = np.zeros(n_points)
        hydrophobicities = np.zeros(n_points)  # 1.0 = hydrophobic, 0.0 = hydrophilic
        
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
                for i, idx in enumerate(indices):
                    atom = atoms[idx]
                    resn = atom['resn']
                    hydrophobicities[i] = hydro_scale.get(resn, 0.5)
        
        # === Method 2: Fallback to residue-based approximation ===
        else:
            print("  Using residue-based electrostatic approximation (fallback)")
            if kdtree:
                # Find nearest atom for each surface point
                dists, indices = kdtree.query(surface_coords, k=1)
                
                for i, idx in enumerate(indices):
                    atom = atoms[idx]
                    resn = atom['resn']
                    
                    # Hydrophobicity
                    hydrophobicities[i] = hydro_scale.get(resn, 0.5)
                    
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
                    
        return {'potential': potentials, 'hydrophobicity': hydrophobicities}
    
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
            p.avg_potential = np.mean(pot[pos_mask]) if np.any(pos_mask) else 0.0
        for p in neg_patches:
            p.avg_potential = np.mean(pot[neg_mask]) if np.any(neg_mask) else 0.0
        
        patches.extend(pos_patches)
        patches.extend(neg_patches)
        
        return patches

    def _detect_hydrophobic_patches(self, surface_indices, properties, surface_coords):
        """Cluster hydrophobic regions"""
        hy = properties['hydrophobicity']
        mask = hy > 0.6 # High hydrophobicity
        
        return self._cluster_surface(surface_indices[mask], "hydrophobic", surface_coords[mask])

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

    def _calculate_geometry_and_score(self, patches, all_surface_coords):
        """
        Calculate curvature and other geometric features.
        Score patches based on area and geometry.
        """
        for p in patches:
            # Curvature estimation (Shape Index)
            # Simplified: compare surface area to volume of convex hull or simpler metric
            # Here: just random placeholder or simple convexity check needed
            # Valid Shape Index: -1 (cup/concave) to +1 (cap/convex)
            # Binding sites are often concave (-0.5 to -1.0)
            
            # Since we don't have local mesh, we use a heuristic
            # If patch center is 'deeper' than surrounding, it's concave
            
            p.shape_index = -0.5 # Assume slightly concave for interesting patches
            
            # Scoring
            # Area is good, concavity is good for binding
            # Hydrophobic: larger area = better
            # Electrostatic: extreme charge = better
            
            size_score = min(p.area / 100.0, 1.0) # 100 A^2 max score
            p.score = size_score * (1.0 - p.shape_index) # Favor concave
            
            # Identify residues
            # Need to search atoms near p.center
            # Using simple distance check
            p.residues = self._find_residues_near(p.center, 6.0) # 6A radius

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
            color_by: 'potential', 'hydrophobicity', or 'combined'
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
        
        # Use different default names for different surface types
        # ESP = Electrostatic Potential
        if surface_name is None:
            if color_by == 'potential':
                surface_name = f"{self.obj_name}_ESP"
            elif color_by == 'hydrophobicity':
                surface_name = f"{self.obj_name}_hydro"
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
            
        elif color_by == 'combined':
            # Combined mode
            return self._render_combined_surface_v2(surface_name, transparency, color_range)
        else:
            print(f"  Unknown color mode: {color_by}")
            return False
    
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
        cmd.create(surface_obj, self.obj_name)
        
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
        
        # Set surface quality for smooth rendering
        cmd.set('surface_quality', 1, surface_obj)
        
        # Set transparency if requested
        if transparency > 0:
            cmd.set('transparency', transparency, surface_obj)
        
        # Hide surface on original object (keep cartoon visible with normal colors)
        cmd.hide('surface', self.obj_name)
        
        # Step 7: Group all ESP-related objects together
        # This keeps the PyMOL object list clean and organized
        cmd.group(group_name, f"{map_name} {ramp_name} {surface_obj}")
        
        # Don't cleanup temp dir yet - we need the DX file for the ramp
        # The map will be deleted when PyMOL session ends
        
        print(f"✅ Electrostatic surface rendered as '{group_name}'")
        print("   Method: PyMOL native ramp_new (same as APBS plugin)")
        print("   Color scheme: Red (negative) → White (neutral) → Blue (positive)")
        print(f"   Group contains: {surface_obj}, {map_name}, {ramp_name}")
        
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
        cmd.create(surface_name, self.obj_name)
        
        model = cmd.get_model(surface_name)
        for i, atom in enumerate(model.atom):
            if i < len(values):
                atom.b = values[i]
        cmd.load_model(model, surface_name, state=1)
        
        cmd.hide('everything', surface_name)
        cmd.show('surface', surface_name)
        cmd.spectrum('b', 'red_white_blue', surface_name, minimum=vmin, maximum=vmax)
        
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
        cmd.create(surface_name, self.obj_name)
        
        model = cmd.get_model(surface_name)
        for i, atom in enumerate(model.atom):
            if i < len(values):
                atom.b = values[i]
        cmd.load_model(model, surface_name, state=1)
        
        cmd.hide('everything', surface_name)
        cmd.show('surface', surface_name)
        cmd.spectrum('b', 'white_green', surface_name, minimum=vmin, maximum=vmax)
        
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
        cmd.create(surface_name, self.obj_name)
        
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
                          save_files: bool = False, output_dir: str = None) -> bool:
    """
    Render complete protein surface with APBS electrostatics using PyMOL native surface.
    
    This creates a smooth, high-quality molecular surface colored by electrostatic
    potential or hydrophobicity.
    
    Args:
        obj_name: PyMOL object name
        use_apbs: Use APBS for accurate electrostatics (default: True)
        ph: pH for PDB2PQR protonation (default: 7.4)
        color_by: 'potential' or 'hydrophobicity'
        transparency: Surface transparency (0.0 = opaque)
        color_min: Minimum value for color scale (default: auto)
        color_max: Maximum value for color scale (default: auto)
        save_files: If True, save intermediate files (PDB, PQR, DX) to output_dir
        output_dir: Directory to save files (required if save_files=True)
        
    Returns:
        True if successful
        
    Usage in PyMOL:
        render_protein_surface '8D7Z', use_apbs=True, ph=7.4
        render_protein_surface '8D7Z', color_by='hydrophobicity'
        render_protein_surface '8D7Z', save_files=True, output_dir='/path/to/output'
    """
    color_range = None
    if color_min is not None and color_max is not None:
        color_range = (color_min, color_max)
    
    analyzer = SurfaceAnalyzer(obj_name, grid_spacing=0.5,
                               use_apbs=use_apbs, ph=ph)
    return analyzer.render_full_surface(color_by=color_by,
                                        transparency=transparency,
                                        color_range=color_range,
                                        save_files=save_files,
                                        output_dir=output_dir)


# Register PyMOL command
if cmd:
    cmd.extend('render_protein_surface', render_protein_surface)
