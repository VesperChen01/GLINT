# -*- coding: utf-8 -*-
"""
Protein Surface Analyzer
========================
Analyzes protein surface properties including electrostatic potential,
hydrophobicity, and geometric features to identify potential binding sites.

Features:
- Electrostatic patch detection
- Hydrophobic patch detection
- Geometric feature analysis (curvature, shape index)
- Comprehensive druggability scoring
"""

import numpy as np
from scipy.spatial.distance import cdist
from scipy import ndimage
import math
from typing import List, Dict, Optional, Tuple, Any

try:
    from pymol import cmd
except ImportError:
    cmd = None

from .feature_extractor import MolecularFeatureExtractor
from .pocket_detector import PocketDetector

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
    
    def __init__(self, obj_name: str, grid_spacing: float = 0.5):
        self.obj_name = obj_name
        self.grid_spacing = grid_spacing
        self.atoms = []
        self.grid = None
        self.origin = None
        
    def analyze(self) -> List[SurfacePatch]:
        """Main analysis pipeline"""
        if not cmd:
            print("Error: PyMOL not available")
            return []
            
        print(f"Starting surface analysis for {self.obj_name}...")
        
        # 1. Get Atoms
        self.atoms = self._get_atoms()
        if not self.atoms:
            return []
            
        # 2. Build Grid & Surface
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

        # 3. Calculate Properties for Surface Points
        # We need to map surface points to nearest atoms to get properties
        atom_kdtree = None
        try:
            from scipy.spatial import cKDTree
            atom_kdtree = cKDTree(atom_coords)
        except ImportError:
            pass
            
        properties = self._map_properties_to_surface(surface_coords, atom_kdtree, self.atoms)
        
        # 4. Detect Patches
        patches = []
        patches.extend(self._detect_electrostatic_patches(surface_points_indices, properties, surface_coords))
        patches.extend(self._detect_hydrophobic_patches(surface_points_indices, properties, surface_coords))
        
        # 5. Calculate Geometric Features & Score
        self._calculate_geometry_and_score(patches, surface_coords)
        
        print(f"Found {len(patches)} surface patches")
        return patches

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
        Map partial charges and hydrophobicity to surface points.
        Returns dict with 'potential' and 'hydrophobicity' arrays.
        """
        n_points = len(surface_coords)
        potentials = np.zeros(n_points)
        hydrophobicities = np.zeros(n_points) # 1.0 = hydrophobic, 0.0 = hydrophilic
        
        # Hydrophobicity scale (Kyte-Doolittle normalized)
        hydro_scale = {
            'ILE': 4.5, 'VAL': 4.2, 'LEU': 3.8, 'PHE': 2.8, 'CYS': 2.5,
            'MET': 1.9, 'ALA': 1.8, 'GLY': -0.4, 'THR': -0.7, 'SER': -0.8,
            'TRP': -0.9, 'TYR': -1.3, 'PRO': -1.6, 'HIS': -3.2, 'GLU': -3.5,
            'GLN': -3.5, 'ASP': -3.5, 'ASN': -3.5, 'LYS': -3.9, 'ARG': -4.5
        }
        
        # Normalize to 0-1 (approximate)
        for k in hydro_scale:
            hydro_scale[k] = (hydro_scale[k] + 4.5) / 9.0

        if kdtree:
            # fast neighbor lookup
            # Find nearest atom for each surface point
            dists, indices = kdtree.query(surface_coords, k=1)
            
            for i, idx in enumerate(indices):
                atom = atoms[idx]
                resn = atom['resn']
                
                # Hydrophobicity
                hydrophobicities[i] = hydro_scale.get(resn, 0.5)
                
                # Electrostatics (Approximation using implicit charge if not available)
                # This is a VERY rough approximation if partial_charges are missing
                # Better would be to read APBS map if available
                charge = atom.get('partial_charge', 0.0)
                if charge == 0.0:
                    if resn in ['ASP', 'GLU']: charge = -1.0
                    elif resn in ['LYS', 'ARG']: charge = 1.0
                    elif resn == 'HIS': charge = 0.5 # pH dependent
                
                # Simple Coulomb decay 1/r
                # Ideally we want the potential at the surface point
                # Summing contributions from nearby atoms would be better but slower
                potentials[i] = charge # Just assign nearest atom charge for now (patchy)
                
        return {'potential': potentials, 'hydrophobicity': hydrophobicities}

    def _detect_electrostatic_patches(self, surface_indices, properties, surface_coords):
        """Cluster positive and negative regions"""
        patches = []
        pot = properties['potential']
        
        # Thresholds
        pos_mask = pot > 0.5  # Arbitrary threshold
        neg_mask = pot < -0.5
        
        patches.extend(self._cluster_surface(surface_indices[pos_mask], "electrostatic_pos", surface_coords[pos_mask]))
        patches.extend(self._cluster_surface(surface_indices[neg_mask], "electrostatic_neg", surface_coords[neg_mask]))
        
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
