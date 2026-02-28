# -*- coding: utf-8 -*-
"""
Surface Similarity and Complementarity Analysis Module
======================================================
MaSIF-inspired protein surface analysis for GLINT.

Features:
- Surface mesh generation (MSMS, EDTSurf, or built-in)
- Geometric features: curvature, shape index, curvedness
- Chemical features: electrostatics, hydrophobicity, H-bond capacity
- Surface similarity comparison
- Surface complementarity analysis for PPI
- Batch analysis for multiple structures

Dependencies:
- Required: numpy, scipy
- Optional: open3d (for advanced mesh processing)
- Optional: MSMS binary (for accurate surface generation)
- Optional: APBS (for accurate electrostatics)

Cross-platform support: macOS, Windows, Linux
"""

from __future__ import annotations
import os
import sys
import math
import tempfile
import subprocess
import platform
from typing import List, Dict, Tuple, Optional, Any, Union
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree

# Import ligand features module for small molecule support
try:
    from .ligand_features import is_ligand, get_ligand_hydrophobicity, get_ligand_charge
    HAS_LIGAND_FEATURES = True
except ImportError:
    HAS_LIGAND_FEATURES = False
from scipy.spatial.distance import cdist
from scipy import ndimage

try:
    from pymol import cmd
    HAS_PYMOL = True
except ImportError:
    cmd = None
    HAS_PYMOL = False

# Optional: Open3D for advanced mesh processing
try:
    import open3d as o3d
    HAS_OPEN3D = True
except ImportError:
    o3d = None
    HAS_OPEN3D = False


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class SurfacePoint:
    """A single point on the protein surface with features."""
    coord: np.ndarray  # (x, y, z)
    normal: np.ndarray  # surface normal vector
    
    # Geometric features
    gaussian_curvature: float = 0.0  # K = κ1 * κ2
    mean_curvature: float = 0.0      # H = (κ1 + κ2) / 2
    shape_index: float = 0.0         # SI ∈ [-1, 1]
    curvedness: float = 0.0          # C = sqrt((κ1² + κ2²)/2)
    
    # Chemical features
    electrostatic_potential: float = 0.0  # in kT/e
    hydrophobicity: float = 0.0           # normalized [0, 1]
    hbond_donor_density: float = 0.0      # donor atoms nearby
    hbond_acceptor_density: float = 0.0   # acceptor atoms nearby
    
    # Metadata
    nearest_residue: Optional[str] = None  # "A:ALA:123"
    
    def to_feature_vector(self, include_geometry: bool = True, 
                          include_chemistry: bool = True) -> np.ndarray:
        """Convert to feature vector for comparison."""
        features = []
        if include_geometry:
            features.extend([
                self.gaussian_curvature,
                self.mean_curvature,
                self.shape_index,
                self.curvedness
            ])
        if include_chemistry:
            features.extend([
                self.electrostatic_potential,
                self.hydrophobicity,
                self.hbond_donor_density,
                self.hbond_acceptor_density
            ])
        return np.array(features)


@dataclass
class SurfacePatch:
    """A patch of surface points centered at a specific location."""
    center: np.ndarray
    radius: float
    points: List[SurfacePoint] = field(default_factory=list)
    
    # Aggregated features
    mean_features: Optional[np.ndarray] = None
    feature_histogram: Optional[np.ndarray] = None
    
    def compute_mean_features(self):
        """Compute mean feature vector for the patch."""
        if not self.points:
            return
        vectors = [p.to_feature_vector() for p in self.points]
        self.mean_features = np.mean(vectors, axis=0)
    
    def compute_feature_histogram(self, n_bins: int = 10):
        """Compute feature histogram for more detailed comparison."""
        if not self.points:
            return
        vectors = np.array([p.to_feature_vector() for p in self.points])
        histograms = []
        for i in range(vectors.shape[1]):
            hist, _ = np.histogram(vectors[:, i], bins=n_bins, density=True)
            histograms.append(hist)
        self.feature_histogram = np.concatenate(histograms)


@dataclass
class SurfaceMesh:
    """Complete surface mesh with vertices, faces, and features."""
    vertices: np.ndarray      # (N, 3) vertex coordinates
    faces: np.ndarray         # (M, 3) face indices
    normals: np.ndarray       # (N, 3) vertex normals
    
    # Per-vertex features
    points: List[SurfacePoint] = field(default_factory=list)
    
    # Metadata
    source_object: str = ""
    generation_method: str = ""
    
    @property
    def n_vertices(self) -> int:
        return len(self.vertices)
    
    @property
    def n_faces(self) -> int:
        return len(self.faces)
    
    @property
    def surface_area(self) -> float:
        """Calculate total surface area."""
        if len(self.faces) == 0:
            return 0.0
        total = 0.0
        for face in self.faces:
            v0, v1, v2 = self.vertices[face]
            # Triangle area = 0.5 * |cross(v1-v0, v2-v0)|
            edge1 = v1 - v0
            edge2 = v2 - v0
            cross = np.cross(edge1, edge2)
            total += 0.5 * np.linalg.norm(cross)
        return total


@dataclass
class SimilarityResult:
    """Result of surface similarity comparison."""
    score: float  # Overall similarity score [0, 1]
    
    # Component scores
    geometric_similarity: float = 0.0
    chemical_similarity: float = 0.0
    
    # Detailed metrics
    shape_index_correlation: float = 0.0
    curvature_correlation: float = 0.0
    electrostatic_correlation: float = 0.0
    hydrophobicity_correlation: float = 0.0
    
    # Matched patches
    matched_patches: List[Tuple[int, int, float]] = field(default_factory=list)


@dataclass
class PatchSearchResult:
    """Result of searching similar patches in a target protein using a template patch."""
    template_patch_idx: int  # Index of the template patch
    template_center: np.ndarray  # Center of template patch
    
    # List of matching patches in target, sorted by similarity (descending)
    # Each tuple: (target_patch_idx, target_center, similarity_score, geometric_sim, chemical_sim)
    matches: List[Tuple[int, np.ndarray, float, float, float]] = field(default_factory=list)
    
    # Best match info
    best_match_idx: int = -1
    best_match_score: float = 0.0
    best_match_center: Optional[np.ndarray] = None


@dataclass
class SimilaritySearchResult:
    """Result of similarity search: using one protein as template to search another."""
    template_object: str  # Name of template protein
    target_object: str    # Name of target protein
    
    # Overall search statistics
    n_template_patches: int = 0
    n_target_patches: int = 0
    n_matches_found: int = 0
    
    # Per-template-patch search results
    patch_results: List[PatchSearchResult] = field(default_factory=list)
    
    # Aggregated scores
    mean_best_similarity: float = 0.0
    max_similarity: float = 0.0
    
    # Best overall match
    best_template_patch_idx: int = -1
    best_target_patch_idx: int = -1
    best_overall_score: float = 0.0
    
    # Residue mapping for best matches
    # List of (template_residue, target_residue, similarity) tuples
    residue_matches: List[Tuple[str, str, float]] = field(default_factory=list)


@dataclass
class ComplementarityResult:
    """Result of surface complementarity analysis."""
    score: float  # Overall complementarity score [0, 1]
    
    # Component scores
    geometric_complementarity: float = 0.0  # Shape matching
    electrostatic_complementarity: float = 0.0  # Charge matching
    hydrophobic_complementarity: float = 0.0  # Hydrophobic matching
    
    # Interface analysis
    interface_area: float = 0.0
    n_contacts: int = 0
    gap_volume: float = 0.0
    
    # Residue pairs at interface
    interface_residues: List[Tuple[str, str]] = field(default_factory=list)


# =============================================================================
# Surface Generation
# =============================================================================

class SurfaceGenerator:
    """
    Generate molecular surface using various methods.
    
    Supports:
    - MSMS (most accurate, requires binary)
    - EDTSurf (built-in approximation)
    - PyMOL surface (if available)
    - Open3D ball pivoting (if available)
    """
    
    def __init__(self, method: str = "auto", msms_path: Optional[str] = None):
        """
        Initialize surface generator.
        
        Args:
            method: "msms", "edtsurf", "pymol", "open3d", or "auto"
            msms_path: Path to MSMS binary (optional)
        """
        self.method = method
        self.msms_path = msms_path or self._find_msms()
        
    def _find_msms(self) -> Optional[str]:
        """Find MSMS binary in common locations."""
        system = platform.system().lower()
        
        # Common installation paths
        search_paths = []
        
        if system == "darwin":  # macOS
            search_paths = [
                "/usr/local/bin/msms",
                os.path.expanduser("~/bin/msms"),
                "/opt/homebrew/bin/msms",
                os.path.expanduser("~/.local/bin/msms"),
            ]
        elif system == "linux":
            search_paths = [
                "/usr/bin/msms",
                "/usr/local/bin/msms",
                os.path.expanduser("~/bin/msms"),
                os.path.expanduser("~/.local/bin/msms"),
            ]
        elif system == "windows":
            search_paths = [
                r"C:\Program Files\MSMS\msms.exe",
                r"C:\msms\msms.exe",
                os.path.expanduser(r"~\msms\msms.exe"),
            ]
        
        for path in search_paths:
            if os.path.isfile(path) and os.access(path, os.X_OK):
                return path
        
        return None
    
    def generate(self, atoms: List[Dict], probe_radius: float = 1.4) -> SurfaceMesh:
        """
        Generate surface mesh from atom coordinates.
        
        Args:
            atoms: List of atom dicts with 'coord', 'element', 'radius' keys
            probe_radius: Solvent probe radius (default 1.4 Å for water)
        
        Returns:
            SurfaceMesh object
        """
        method = self.method
        
        if method == "auto":
            if self.msms_path:
                method = "msms"
            elif HAS_OPEN3D:
                method = "open3d"
            else:
                method = "edtsurf"
        
        if method == "msms":
            return self._generate_msms(atoms, probe_radius)
        elif method == "open3d":
            return self._generate_open3d(atoms, probe_radius)
        else:
            return self._generate_edtsurf(atoms, probe_radius)
    
    def _generate_msms(self, atoms: List[Dict], probe_radius: float) -> SurfaceMesh:
        """Generate surface using MSMS."""
        if not self.msms_path:
            raise RuntimeError("MSMS binary not found")
        
        # Create temporary files
        with tempfile.TemporaryDirectory() as tmpdir:
            xyzr_file = os.path.join(tmpdir, "input.xyzr")
            vert_file = os.path.join(tmpdir, "output.vert")
            face_file = os.path.join(tmpdir, "output.face")
            
            # Write XYZR file
            with open(xyzr_file, 'w') as f:
                for atom in atoms:
                    x, y, z = atom['coord']
                    r = atom.get('radius', 1.7)  # Default VDW radius
                    f.write(f"{x:.3f} {y:.3f} {z:.3f} {r:.3f}\n")
            
            # Run MSMS
            cmd_args = [
                self.msms_path,
                "-if", xyzr_file,
                "-of", os.path.join(tmpdir, "output"),
                "-probe_radius", str(probe_radius),
                "-density", "3.0",  # Vertex density
                "-no_header"
            ]
            
            try:
                subprocess.run(cmd_args, check=True, capture_output=True)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"MSMS failed: {e.stderr.decode()}")
            
            # Read vertices
            vertices = []
            normals = []
            with open(vert_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 6:
                        vertices.append([float(parts[0]), float(parts[1]), float(parts[2])])
                        normals.append([float(parts[3]), float(parts[4]), float(parts[5])])
            
            # Read faces
            faces = []
            with open(face_file, 'r') as f:
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 3:
                        # MSMS uses 1-based indexing
                        faces.append([int(parts[0])-1, int(parts[1])-1, int(parts[2])-1])
        
        return SurfaceMesh(
            vertices=np.array(vertices),
            faces=np.array(faces),
            normals=np.array(normals),
            generation_method="msms"
        )
    
    def _generate_open3d(self, atoms: List[Dict], probe_radius: float) -> SurfaceMesh:
        """Generate surface using Open3D."""
        if not HAS_OPEN3D:
            raise RuntimeError("Open3D not available")
        
        # Create point cloud from atom positions
        coords = np.array([a['coord'] for a in atoms])
        radii = np.array([a.get('radius', 1.7) for a in atoms])
        
        # Expand points to surface using spherical sampling
        surface_points = []
        for coord, radius in zip(coords, radii):
            # Sample points on sphere surface
            n_samples = max(10, int(4 * np.pi * (radius + probe_radius)**2))
            phi = np.random.uniform(0, 2*np.pi, n_samples)
            costheta = np.random.uniform(-1, 1, n_samples)
            theta = np.arccos(costheta)
            
            r = radius + probe_radius
            x = coord[0] + r * np.sin(theta) * np.cos(phi)
            y = coord[1] + r * np.sin(theta) * np.sin(phi)
            z = coord[2] + r * np.cos(theta)
            
            surface_points.extend(zip(x, y, z))
        
        surface_points = np.array(surface_points)
        
        # Create point cloud
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(surface_points)
        
        # Estimate normals
        pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(
            radius=2.0, max_nn=30))
        
        # Ball pivoting reconstruction
        radii_list = [0.5, 1.0, 1.5, 2.0]
        mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(
            pcd, o3d.utility.DoubleVector(radii_list))
        
        # Extract vertices, faces, normals
        vertices = np.asarray(mesh.vertices)
        faces = np.asarray(mesh.triangles)
        
        # Compute vertex normals
        mesh.compute_vertex_normals()
        normals = np.asarray(mesh.vertex_normals)
        
        return SurfaceMesh(
            vertices=vertices,
            faces=faces,
            normals=normals,
            generation_method="open3d"
        )
    
    def _generate_edtsurf(self, atoms: List[Dict], probe_radius: float) -> SurfaceMesh:
        """
        Generate surface using EDT (Euclidean Distance Transform) method.
        Pure NumPy/SciPy implementation.
        """
        coords = np.array([a['coord'] for a in atoms])
        radii = np.array([a.get('radius', 1.7) for a in atoms])
        
        # Grid parameters
        grid_spacing = 0.5
        padding = probe_radius + 3.0
        
        min_coords = coords.min(axis=0) - padding - radii.max()
        max_coords = coords.max(axis=0) + padding + radii.max()
        
        grid_size = np.ceil((max_coords - min_coords) / grid_spacing).astype(int)
        
        # Create distance field
        distance_field = np.full(grid_size, np.inf)
        
        # Mark atom interiors
        for coord, radius in zip(coords, radii):
            grid_coord = ((coord - min_coords) / grid_spacing).astype(int)
            r_grid = int(np.ceil((radius + probe_radius) / grid_spacing)) + 1
            
            for i in range(max(0, grid_coord[0]-r_grid), min(grid_size[0], grid_coord[0]+r_grid+1)):
                for j in range(max(0, grid_coord[1]-r_grid), min(grid_size[1], grid_coord[1]+r_grid+1)):
                    for k in range(max(0, grid_coord[2]-r_grid), min(grid_size[2], grid_coord[2]+r_grid+1)):
                        point = min_coords + np.array([i, j, k]) * grid_spacing
                        dist = np.linalg.norm(point - coord) - radius - probe_radius
                        distance_field[i, j, k] = min(distance_field[i, j, k], dist)
        
        # Extract isosurface at distance = 0 using marching cubes
        vertices, faces, normals = self._marching_cubes(
            distance_field, 0.0, min_coords, grid_spacing)
        
        return SurfaceMesh(
            vertices=vertices,
            faces=faces,
            normals=normals,
            generation_method="edtsurf"
        )
    
    def _marching_cubes(self, volume: np.ndarray, level: float,
                        origin: np.ndarray, spacing: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Simple marching cubes implementation.
        For production, consider using skimage.measure.marching_cubes
        """
        try:
            from skimage.measure import marching_cubes
            verts, faces, normals, _ = marching_cubes(volume, level, spacing=(spacing, spacing, spacing))
            verts = verts + origin
            return verts, faces, normals
        except ImportError:
            # Fallback: simple surface extraction
            return self._simple_surface_extraction(volume, level, origin, spacing)
    
    def _simple_surface_extraction(self, volume: np.ndarray, level: float,
                                   origin: np.ndarray, spacing: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Simple surface point extraction (fallback)."""
        # Find surface voxels (sign change)
        surface_mask = np.zeros(volume.shape, dtype=bool)

        # Check neighbors for sign change
        for axis in range(3):
            shifted_pos = np.roll(volume, 1, axis=axis)
            shifted_neg = np.roll(volume, -1, axis=axis)
            surface_mask |= ((volume * shifted_pos) < 0) | ((volume * shifted_neg) < 0)

        # Extract surface points
        indices = np.argwhere(surface_mask & (np.abs(volume) < spacing))
        vertices = origin + indices * spacing

        # Estimate normals from gradient
        grad = np.gradient(volume)
        normals = np.zeros((len(vertices), 3))
        for i, idx in enumerate(indices):
            n = np.array([grad[0][tuple(idx)], grad[1][tuple(idx)], grad[2][tuple(idx)]])
            norm = np.linalg.norm(n)
            normals[i] = n / norm if norm > 0 else [0, 0, 1]

        # Create triangulation using nearest neighbors
        faces = self._triangulate_surface_points(vertices, normals, spacing)

        return vertices, faces, normals

    def _triangulate_surface_points(self, vertices: np.ndarray, normals: np.ndarray,
                                    spacing: float) -> np.ndarray:
        """
        Create triangular faces from surface points using local Delaunay-like triangulation.
        This is a simplified approach that works for molecular surfaces.
        """
        if len(vertices) < 3:
            return np.array([]).reshape(0, 3).astype(int)

        # Use KDTree for efficient neighbor finding
        tree = cKDTree(vertices)

        faces = []
        used_triangles = set()

        # For each vertex, try to form triangles with nearby vertices
        # Use a distance threshold based on grid spacing
        max_edge_length = spacing * 2.5

        for i in range(len(vertices)):
            # Find nearby vertices
            nearby_indices = tree.query_ball_point(vertices[i], max_edge_length)
            nearby_indices = [j for j in nearby_indices if j != i]

            if len(nearby_indices) < 2:
                continue

            # Sort by distance
            distances = [np.linalg.norm(vertices[j] - vertices[i]) for j in nearby_indices]
            sorted_neighbors = [x for _, x in sorted(zip(distances, nearby_indices))]

            # Try to form triangles with pairs of nearby vertices
            for j_idx, j in enumerate(sorted_neighbors[:8]):  # Limit to 8 nearest
                for k in sorted_neighbors[j_idx+1:8]:
                    # Check if j and k are also close to each other
                    jk_dist = np.linalg.norm(vertices[j] - vertices[k])
                    if jk_dist > max_edge_length:
                        continue

                    # Create canonical triangle representation (sorted indices)
                    tri = tuple(sorted([i, j, k]))
                    if tri in used_triangles:
                        continue

                    # Check normal consistency - triangle normal should roughly align with vertex normals
                    v0, v1, v2 = vertices[i], vertices[j], vertices[k]
                    edge1 = v1 - v0
                    edge2 = v2 - v0
                    tri_normal = np.cross(edge1, edge2)
                    tri_norm_len = np.linalg.norm(tri_normal)

                    if tri_norm_len < 1e-10:
                        continue  # Degenerate triangle

                    tri_normal = tri_normal / tri_norm_len

                    # Check if triangle normal aligns with average vertex normal
                    avg_normal = (normals[i] + normals[j] + normals[k]) / 3
                    avg_norm_len = np.linalg.norm(avg_normal)
                    if avg_norm_len > 0:
                        avg_normal = avg_normal / avg_norm_len
                        dot = np.dot(tri_normal, avg_normal)

                        # Accept if normals are roughly aligned (allow some tolerance)
                        if abs(dot) > 0.3:
                            # Ensure consistent winding order
                            if dot < 0:
                                faces.append([i, k, j])  # Flip winding
                            else:
                                faces.append([i, j, k])
                            used_triangles.add(tri)

        if len(faces) == 0:
            return np.array([]).reshape(0, 3).astype(int)

        return np.array(faces, dtype=int)


# =============================================================================
# Feature Extraction
# =============================================================================

class FeatureExtractor:
    """Extract geometric and chemical features for surface points."""
    
    # Kyte-Doolittle hydrophobicity scale (normalized to [0, 1])
    HYDROPHOBICITY_SCALE = {
        'ILE': 1.000, 'VAL': 0.967, 'LEU': 0.922, 'PHE': 0.811, 'CYS': 0.722,
        'MET': 0.656, 'ALA': 0.622, 'GLY': 0.456, 'THR': 0.411, 'SER': 0.400,
        'TRP': 0.389, 'TYR': 0.356, 'PRO': 0.322, 'HIS': 0.144, 'GLU': 0.111,
        'GLN': 0.111, 'ASP': 0.111, 'ASN': 0.111, 'LYS': 0.067, 'ARG': 0.000
    }
    
    # Partial charges for common atoms (simplified)
    PARTIAL_CHARGES = {
        'N': -0.4, 'O': -0.5, 'C': 0.1, 'S': -0.2, 'H': 0.2,
        'NZ': 1.0,   # Lys
        'NH1': 0.5, 'NH2': 0.5,  # Arg
        'OD1': -0.5, 'OD2': -0.5,  # Asp
        'OE1': -0.5, 'OE2': -0.5,  # Glu
    }
    
    def __init__(self, atoms: List[Dict], apbs_path: Optional[str] = None):
        """
        Initialize feature extractor.
        
        Args:
            atoms: List of atom dicts with 'coord', 'element', 'resn', 'resi', 'chain', 'name'
            apbs_path: Path to APBS binary for accurate electrostatics (optional)
        """
        self.atoms = atoms
        self.apbs_path = apbs_path or self._find_apbs()
        
        # Build spatial index
        self.coords = np.array([a['coord'] for a in atoms])
        self.kdtree = cKDTree(self.coords)
        
        # Precompute residue properties
        self._precompute_residue_properties()
    
    def _find_apbs(self) -> Optional[str]:
        """Find APBS binary."""
        system = platform.system().lower()
        
        search_paths = []
        if system == "darwin":
            search_paths = [
                "/usr/local/bin/apbs",
                os.path.expanduser("~/bin/apbs"),
                "/opt/homebrew/bin/apbs",
            ]
        elif system == "linux":
            search_paths = [
                "/usr/bin/apbs",
                "/usr/local/bin/apbs",
                os.path.expanduser("~/bin/apbs"),
            ]
        elif system == "windows":
            search_paths = [
                r"C:\Program Files\APBS\apbs.exe",
                r"C:\APBS\apbs.exe",
            ]
        
        for path in search_paths:
            if os.path.isfile(path) and os.access(path, os.X_OK):
                return path
        
        return None
    
    def _precompute_residue_properties(self):
        """Precompute properties for each residue, including ligands."""
        self.residue_hydrophobicity = {}
        self.residue_charge = {}
        self.atom_hydrophobicity = {}  # Per-atom for ligands
        self.atom_charge = {}  # Per-atom for ligands
        
        for i, atom in enumerate(self.atoms):
            key = (atom.get('chain', ''), atom.get('resn', ''), atom.get('resi', ''))
            resn = atom.get('resn', 'ALA')
            name = atom.get('name', '')
            elem = atom.get('element', 'C')
            
            # Check if this is a ligand atom
            is_lig = HAS_LIGAND_FEATURES and is_ligand(resn)
            
            if is_lig:
                # Use atom-type based features for ligands
                self.atom_hydrophobicity[i] = get_ligand_hydrophobicity(name, elem, resn)
                self.atom_charge[i] = get_ligand_charge(name, elem, resn)
                # Also set residue-level defaults
                if key not in self.residue_hydrophobicity:
                    self.residue_hydrophobicity[key] = 0.5  # Neutral default
                if key not in self.residue_charge:
                    self.residue_charge[key] = 0.0
            else:
                # Standard amino acid handling
                if key not in self.residue_hydrophobicity:
                    self.residue_hydrophobicity[key] = self.HYDROPHOBICITY_SCALE.get(resn, 0.5)
                
                if key not in self.residue_charge:
                    if resn in ['LYS', 'ARG', 'HIS']:
                        self.residue_charge[key] = 1.0
                    elif resn in ['ASP', 'GLU']:
                        self.residue_charge[key] = -1.0
                    else:
                        self.residue_charge[key] = 0.0
    
    def extract_features(self, mesh: SurfaceMesh, 
                        use_apbs: bool = False) -> List[SurfacePoint]:
        """
        Extract features for all surface vertices.
        
        Args:
            mesh: Surface mesh
            use_apbs: Use APBS for electrostatics (slower but more accurate)
        
        Returns:
            List of SurfacePoint objects with features
        """
        points = []
        
        # Compute curvatures
        curvatures = self._compute_curvatures(mesh)
        
        # Compute electrostatics
        if use_apbs and self.apbs_path:
            electrostatics = self._compute_apbs_electrostatics(mesh)
        else:
            electrostatics = self._compute_coulomb_electrostatics(mesh)
        
        # Extract features for each vertex
        for i, (vertex, normal) in enumerate(zip(mesh.vertices, mesh.normals)):
            point = SurfacePoint(
                coord=vertex,
                normal=normal,
                gaussian_curvature=curvatures['gaussian'][i],
                mean_curvature=curvatures['mean'][i],
                shape_index=curvatures['shape_index'][i],
                curvedness=curvatures['curvedness'][i],
                electrostatic_potential=electrostatics[i],
            )
            
            # Find nearest residue and compute local properties
            self._compute_local_properties(point)
            
            points.append(point)
        
        return points
    
    def _compute_curvatures(self, mesh: SurfaceMesh) -> Dict[str, np.ndarray]:
        """Compute curvature features for all vertices."""
        n = len(mesh.vertices)
        
        gaussian = np.zeros(n)
        mean = np.zeros(n)
        shape_index = np.zeros(n)
        curvedness = np.zeros(n)
        
        if HAS_OPEN3D and len(mesh.faces) > 0:
            # Use Open3D for accurate curvature computation
            o3d_mesh = o3d.geometry.TriangleMesh()
            o3d_mesh.vertices = o3d.utility.Vector3dVector(mesh.vertices)
            o3d_mesh.triangles = o3d.utility.Vector3iVector(mesh.faces)
            o3d_mesh.compute_vertex_normals()
            
            # Compute curvatures using quadric fitting
            for i in range(n):
                k1, k2 = self._estimate_principal_curvatures_open3d(
                    mesh.vertices, mesh.normals, i, mesh.faces)
                
                gaussian[i] = k1 * k2
                mean[i] = (k1 + k2) / 2
                
                # Shape index: SI = 2/π * arctan((k1+k2)/(k1-k2))
                if abs(k1 - k2) > 1e-10:
                    shape_index[i] = (2/np.pi) * np.arctan((k1 + k2) / (k1 - k2))
                else:
                    shape_index[i] = 0.0
                
                # Curvedness: C = sqrt((k1² + k2²)/2)
                curvedness[i] = np.sqrt((k1**2 + k2**2) / 2)
        else:
            # Fallback: estimate from local neighborhood
            vertex_tree = cKDTree(mesh.vertices)
            
            for i in range(n):
                k1, k2 = self._estimate_principal_curvatures_local(
                    mesh.vertices, mesh.normals, i, vertex_tree)
                
                gaussian[i] = k1 * k2
                mean[i] = (k1 + k2) / 2
                
                if abs(k1 - k2) > 1e-10:
                    shape_index[i] = (2/np.pi) * np.arctan((k1 + k2) / (k1 - k2))
                
                curvedness[i] = np.sqrt((k1**2 + k2**2) / 2)
        
        return {
            'gaussian': gaussian,
            'mean': mean,
            'shape_index': shape_index,
            'curvedness': curvedness
        }
    
    def _estimate_principal_curvatures_open3d(self, vertices: np.ndarray, 
                                               normals: np.ndarray, 
                                               idx: int,
                                               faces: np.ndarray) -> Tuple[float, float]:
        """Estimate principal curvatures using Open3D."""
        # Find neighboring vertices
        vertex = vertices[idx]
        normal = normals[idx]
        
        # Get 1-ring neighbors from faces
        neighbors = set()
        for face in faces:
            if idx in face:
                neighbors.update(face)
        neighbors.discard(idx)
        
        if len(neighbors) < 3:
            return 0.0, 0.0
        
        neighbor_verts = vertices[list(neighbors)]
        
        # Fit quadric surface
        return self._fit_quadric_curvature(vertex, normal, neighbor_verts)
    
    def _estimate_principal_curvatures_local(self, vertices: np.ndarray,
                                              normals: np.ndarray,
                                              idx: int,
                                              kdtree: cKDTree) -> Tuple[float, float]:
        """Estimate principal curvatures from local neighborhood."""
        vertex = vertices[idx]
        normal = normals[idx]
        
        # Find k nearest neighbors
        k = min(20, len(vertices) - 1)
        _, neighbor_indices = kdtree.query(vertex, k=k+1)
        neighbor_indices = neighbor_indices[1:]  # Exclude self
        
        neighbor_verts = vertices[neighbor_indices]
        
        return self._fit_quadric_curvature(vertex, normal, neighbor_verts)
    
    def _fit_quadric_curvature(self, vertex: np.ndarray, normal: np.ndarray,
                                neighbors: np.ndarray) -> Tuple[float, float]:
        """Fit quadric surface to estimate principal curvatures."""
        # Transform to local coordinate system
        # z-axis aligned with normal
        z = normal / np.linalg.norm(normal)
        
        # Find perpendicular vectors
        if abs(z[0]) < 0.9:
            x = np.cross(z, [1, 0, 0])
        else:
            x = np.cross(z, [0, 1, 0])
        x = x / np.linalg.norm(x)
        y = np.cross(z, x)
        
        # Transform neighbors to local coordinates
        local_coords = []
        for n in neighbors:
            d = n - vertex
            local_coords.append([np.dot(d, x), np.dot(d, y), np.dot(d, z)])
        local_coords = np.array(local_coords)
        
        if len(local_coords) < 5:
            return 0.0, 0.0
        
        # Fit z = ax² + bxy + cy² (paraboloid)
        # Using least squares
        X = local_coords[:, :2]
        Z = local_coords[:, 2]
        
        # Design matrix: [x², xy, y²]
        A = np.column_stack([X[:, 0]**2, X[:, 0]*X[:, 1], X[:, 1]**2])
        
        try:
            coeffs, _, _, _ = np.linalg.lstsq(A, Z, rcond=None)
            a, b, c = coeffs
            
            # Principal curvatures from Hessian
            # H = [[2a, b], [b, 2c]]
            # Eigenvalues give principal curvatures
            H = np.array([[2*a, b], [b, 2*c]])
            eigenvalues = np.linalg.eigvalsh(H)
            
            return float(eigenvalues[1]), float(eigenvalues[0])  # k1 >= k2
        except (np.linalg.LinAlgError, ValueError):  # 线性代数计算可能Failed
            return 0.0, 0.0
    
    def _compute_coulomb_electrostatics(self, mesh: SurfaceMesh) -> np.ndarray:
        """Compute electrostatic potential using Coulomb's law."""
        n = len(mesh.vertices)
        potentials = np.zeros(n)
        
        # Dielectric constant
        epsilon = 80.0  # Water
        k_coulomb = 332.0  # kcal/mol·Å·e²
        
        # Assign charges to atoms
        charges = []
        for atom in self.atoms:
            resn = atom.get('resn', '')
            name = atom.get('name', '')
            
            # Get charge from lookup or estimate
            if name in self.PARTIAL_CHARGES:
                q = self.PARTIAL_CHARGES[name]
            elif resn in ['LYS', 'ARG']:
                q = 0.2 if name.startswith('N') else 0.0
            elif resn in ['ASP', 'GLU']:
                q = -0.2 if name.startswith('O') else 0.0
            else:
                q = 0.0
            charges.append(q)
        
        charges = np.array(charges)
        
        # Compute potential at each surface point
        for i, vertex in enumerate(mesh.vertices):
            distances = np.linalg.norm(self.coords - vertex, axis=1)
            distances = np.maximum(distances, 1.0)  # Avoid division by zero
            
            # Coulomb potential: V = k * q / (ε * r)
            potentials[i] = np.sum(k_coulomb * charges / (epsilon * distances))
        
        # Convert to kT/e units (at 300K)
        kT = 0.592  # kcal/mol at 300K
        potentials = potentials / kT
        
        return potentials
    
    def _compute_apbs_electrostatics(self, mesh: SurfaceMesh) -> np.ndarray:
        """Compute electrostatic potential using APBS."""
        if not self.apbs_path:
            return self._compute_coulomb_electrostatics(mesh)
        
        # This would require PDB2PQR and APBS setup
        # For now, fall back to Coulomb
        print("APBS integration not yet implemented, using Coulomb approximation")
        return self._compute_coulomb_electrostatics(mesh)
    
    def _compute_local_properties(self, point: SurfacePoint):
        """Compute local chemical properties for a surface point."""
        # Find nearest atoms
        distances, indices = self.kdtree.query(point.coord, k=min(10, len(self.atoms)))
        
        if isinstance(indices, np.integer):
            indices = [indices]
            distances = [distances]
        
        # Nearest residue
        nearest_atom = self.atoms[indices[0]]
        point.nearest_residue = f"{nearest_atom.get('chain', '')}:{nearest_atom.get('resn', '')}:{nearest_atom.get('resi', '')}"
        
        # Hydrophobicity (weighted average of nearby atoms/residues)
        hydro_sum = 0.0
        weight_sum = 0.0
        
        for idx, dist in zip(indices, distances):
            if dist > 6.0:
                continue
            
            # Check if we have atom-level hydrophobicity (for ligands)
            if idx in self.atom_hydrophobicity:
                hydro = self.atom_hydrophobicity[idx]
            else:
                atom = self.atoms[idx]
                key = (atom.get('chain', ''), atom.get('resn', ''), atom.get('resi', ''))
                hydro = self.residue_hydrophobicity.get(key, 0.5)
            
            weight = 1.0 / (dist + 0.1)
            hydro_sum += hydro * weight
            weight_sum += weight
        
        point.hydrophobicity = hydro_sum / weight_sum if weight_sum > 0 else 0.5
        
        # H-bond donor/acceptor density
        donor_count = 0
        acceptor_count = 0
        
        for idx, dist in zip(indices, distances):
            if dist > 4.0:
                continue
            atom = self.atoms[idx]
            elem = atom.get('element', 'C').upper()
            name = atom.get('name', '').upper()
            
            # Donors: N-H, O-H
            if elem == 'N' or (elem == 'O' and 'H' in name):
                donor_count += 1
            
            # Acceptors: O, N (sp2)
            if elem == 'O' or (elem == 'N' and atom.get('resn', '') not in ['LYS', 'ARG']):
                acceptor_count += 1
        
        point.hbond_donor_density = donor_count / 10.0  # Normalize
        point.hbond_acceptor_density = acceptor_count / 10.0


# =============================================================================
# Similarity and Complementarity Analysis
# =============================================================================

class SurfaceComparator:
    """Compare surfaces for similarity and complementarity."""
    
    def __init__(self, patch_radius: float = 12.0, n_patches: int = 100):
        """
        Initialize comparator.
        
        Args:
            patch_radius: Radius of surface patches for comparison (Å)
            n_patches: Number of patches to sample
        """
        self.patch_radius = patch_radius
        self.n_patches = n_patches
    
    def compute_similarity(self, mesh1: SurfaceMesh, mesh2: SurfaceMesh,
                          points1: List[SurfacePoint], 
                          points2: List[SurfacePoint]) -> SimilarityResult:
        """
        Compute similarity between two surfaces.
        
        Args:
            mesh1, mesh2: Surface meshes
            points1, points2: Surface points with features
        
        Returns:
            SimilarityResult
        """
        # Extract patches from both surfaces
        patches1 = self._extract_patches(mesh1, points1)
        patches2 = self._extract_patches(mesh2, points2)
        
        # Compute patch features
        for patch in patches1 + patches2:
            patch.compute_mean_features()
            patch.compute_feature_histogram()
        
        # Match patches and compute similarity
        matched_patches = []
        similarities = []
        
        for i, p1 in enumerate(patches1):
            best_match = -1
            best_sim = -1
            
            for j, p2 in enumerate(patches2):
                sim = self._patch_similarity(p1, p2)
                if sim > best_sim:
                    best_sim = sim
                    best_match = j
            
            if best_match >= 0:
                matched_patches.append((i, best_match, best_sim))
                similarities.append(best_sim)
        
        # Compute component scores
        geo_sims = []
        chem_sims = []
        
        for i, j, _ in matched_patches:
            p1, p2 = patches1[i], patches2[j]
            geo_sims.append(self._geometric_similarity(p1, p2))
            chem_sims.append(self._chemical_similarity(p1, p2))
        
        # Compute correlations
        si1 = np.array([p.shape_index for p in points1])
        si2 = np.array([p.shape_index for p in points2])
        
        curv1 = np.array([p.curvedness for p in points1])
        curv2 = np.array([p.curvedness for p in points2])
        
        esp1 = np.array([p.electrostatic_potential for p in points1])
        esp2 = np.array([p.electrostatic_potential for p in points2])
        
        hydro1 = np.array([p.hydrophobicity for p in points1])
        hydro2 = np.array([p.hydrophobicity for p in points2])
        
        return SimilarityResult(
            score=np.mean(similarities) if similarities else 0.0,
            geometric_similarity=np.mean(geo_sims) if geo_sims else 0.0,
            chemical_similarity=np.mean(chem_sims) if chem_sims else 0.0,
            shape_index_correlation=self._histogram_correlation(si1, si2),
            curvature_correlation=self._histogram_correlation(curv1, curv2),
            electrostatic_correlation=self._histogram_correlation(esp1, esp2),
            hydrophobicity_correlation=self._histogram_correlation(hydro1, hydro2),
            matched_patches=matched_patches
        )
    
    def compute_complementarity(self, mesh1: SurfaceMesh, mesh2: SurfaceMesh,
                                points1: List[SurfacePoint],
                                points2: List[SurfacePoint],
                                interface_distance: float = 4.0) -> ComplementarityResult:
        """
        Compute complementarity between two surfaces (for PPI analysis).

        Args:
            mesh1, mesh2: Surface meshes
            points1, points2: Surface points with features
            interface_distance: Distance threshold for interface contacts (Å)

        Returns:
            ComplementarityResult

        Note:
            This analysis is designed for surfaces that are in contact (e.g.,
            receptor-ligand interface within the same complex). If the two
            surfaces are from different structures that are not spatially
            aligned, no interface will be found.
        """
        # Find interface points
        coords1 = np.array([p.coord for p in points1])
        coords2 = np.array([p.coord for p in points2])

        tree2 = cKDTree(coords2)

        interface_pairs = []
        interface_points1 = []
        interface_points2 = []

        # Also track minimum distance for diagnostic purposes
        min_distance = float('inf')

        for i, coord in enumerate(coords1):
            distances, indices = tree2.query(coord, k=1)
            if distances < min_distance:
                min_distance = distances
            if distances < interface_distance:
                interface_pairs.append((i, indices))
                interface_points1.append(points1[i])
                interface_points2.append(points2[indices])

        if not interface_pairs:
            # Provide diagnostic information
            print(f"\n⚠️  No interface found between the two surfaces!")
            print(f"   Minimum distance between surfaces: {min_distance:.1f} Å")
            print(f"   Interface distance threshold: {interface_distance:.1f} Å")
            if min_distance > 50:
                print(f"\n💡 The surfaces appear to be from different structures that are not spatially aligned.")
                print(f"   For complementarity analysis, use chains from the SAME complex (e.g., 'chain A' vs 'chain B').")
                print(f"   For comparing different structures, use 'Similarity Search' instead.")
            elif min_distance > interface_distance:
                print(f"\n💡 Try increasing the 'Interface Distance' parameter to {min_distance + 2:.1f} Å or more.")

            return ComplementarityResult(
                score=0.0,
                geometric_complementarity=0.0,
                electrostatic_complementarity=0.0,
                hydrophobic_complementarity=0.0,
                interface_area=0.0,
                n_contacts=0
            )
        
        # Compute complementarity scores
        geo_comp = self._geometric_complementarity(interface_points1, interface_points2)
        esp_comp = self._electrostatic_complementarity(interface_points1, interface_points2)
        hydro_comp = self._hydrophobic_complementarity(interface_points1, interface_points2)
        
        # Estimate interface area
        interface_area = len(interface_pairs) * (interface_distance ** 2) * 0.5
        
        # Get interface residues
        interface_residues = []
        seen = set()
        for p1, p2 in zip(interface_points1, interface_points2):
            pair = (p1.nearest_residue, p2.nearest_residue)
            if pair not in seen:
                seen.add(pair)
                interface_residues.append(pair)
        
        # Overall score (weighted average)
        overall = 0.4 * geo_comp + 0.35 * esp_comp + 0.25 * hydro_comp
        
        return ComplementarityResult(
            score=overall,
            geometric_complementarity=geo_comp,
            electrostatic_complementarity=esp_comp,
            hydrophobic_complementarity=hydro_comp,
            interface_area=interface_area,
            n_contacts=len(interface_pairs),
            interface_residues=interface_residues
        )
    
    def _extract_patches(self, mesh: SurfaceMesh, 
                        points: List[SurfacePoint]) -> List[SurfacePatch]:
        """Extract surface patches."""
        patches = []
        coords = np.array([p.coord for p in points])
        
        # Sample patch centers
        if len(coords) <= self.n_patches:
            centers = coords
        else:
            # Farthest point sampling for uniform coverage
            centers = self._farthest_point_sampling(coords, self.n_patches)
        
        tree = cKDTree(coords)
        
        for center in centers:
            indices = tree.query_ball_point(center, self.patch_radius)
            patch_points = [points[i] for i in indices]
            
            if len(patch_points) >= 5:  # Minimum points for meaningful patch
                patches.append(SurfacePatch(
                    center=center,
                    radius=self.patch_radius,
                    points=patch_points
                ))
        
        return patches
    
    def _farthest_point_sampling(self, points: np.ndarray, n: int) -> np.ndarray:
        """Sample n points using farthest point sampling."""
        if len(points) <= n:
            return points
        
        selected = [0]
        distances = np.full(len(points), np.inf)
        
        for _ in range(n - 1):
            last = points[selected[-1]]
            new_distances = np.linalg.norm(points - last, axis=1)
            distances = np.minimum(distances, new_distances)
            selected.append(np.argmax(distances))
        
        return points[selected]
    
    def _patch_similarity(self, p1: SurfacePatch, p2: SurfacePatch) -> float:
        """Compute similarity between two patches."""
        if p1.mean_features is None or p2.mean_features is None:
            return 0.0
        
        # Cosine similarity of mean features
        dot = np.dot(p1.mean_features, p2.mean_features)
        norm1 = np.linalg.norm(p1.mean_features)
        norm2 = np.linalg.norm(p2.mean_features)
        
        if norm1 == 0 or norm2 == 0:
            return 0.0
        
        return (dot / (norm1 * norm2) + 1) / 2  # Normalize to [0, 1]
    
    def _geometric_similarity(self, p1: SurfacePatch, p2: SurfacePatch) -> float:
        """Compute geometric similarity."""
        si1 = np.array([pt.shape_index for pt in p1.points])
        si2 = np.array([pt.shape_index for pt in p2.points])
        
        curv1 = np.array([pt.curvedness for pt in p1.points])
        curv2 = np.array([pt.curvedness for pt in p2.points])
        
        si_corr = self._histogram_correlation(si1, si2)
        curv_corr = self._histogram_correlation(curv1, curv2)
        
        return (si_corr + curv_corr) / 2
    
    def _chemical_similarity(self, p1: SurfacePatch, p2: SurfacePatch) -> float:
        """Compute chemical similarity."""
        esp1 = np.array([pt.electrostatic_potential for pt in p1.points])
        esp2 = np.array([pt.electrostatic_potential for pt in p2.points])
        
        hydro1 = np.array([pt.hydrophobicity for pt in p1.points])
        hydro2 = np.array([pt.hydrophobicity for pt in p2.points])
        
        esp_corr = self._histogram_correlation(esp1, esp2)
        hydro_corr = self._histogram_correlation(hydro1, hydro2)
        
        return (esp_corr + hydro_corr) / 2
    
    def _geometric_complementarity(self, points1: List[SurfacePoint],
                                   points2: List[SurfacePoint]) -> float:
        """
        Compute geometric complementarity.
        Complementary surfaces should have opposite shape indices (convex-concave).
        """
        if not points1 or not points2:
            return 0.0
        
        complementarity_scores = []
        
        for p1, p2 in zip(points1, points2):
            # Shape index complementarity: SI1 ≈ -SI2
            si_comp = 1.0 - abs(p1.shape_index + p2.shape_index) / 2
            
            # Normal alignment: normals should be roughly opposite
            dot = np.dot(p1.normal, p2.normal)
            normal_comp = (1 - dot) / 2  # -1 (opposite) -> 1, +1 (same) -> 0
            
            complementarity_scores.append(0.6 * si_comp + 0.4 * normal_comp)
        
        return np.mean(complementarity_scores)
    
    def _electrostatic_complementarity(self, points1: List[SurfacePoint],
                                       points2: List[SurfacePoint]) -> float:
        """
        Compute electrostatic complementarity.
        Complementary surfaces should have opposite charges.
        """
        if not points1 or not points2:
            return 0.0
        
        complementarity_scores = []
        
        for p1, p2 in zip(points1, points2):
            esp1 = p1.electrostatic_potential
            esp2 = p2.electrostatic_potential
            
            # Opposite charges are complementary
            # Normalize to [-1, 1] range first
            esp1_norm = np.tanh(esp1 / 5.0)
            esp2_norm = np.tanh(esp2 / 5.0)
            
            # Complementarity: opposite signs
            comp = (1 - esp1_norm * esp2_norm) / 2
            complementarity_scores.append(comp)
        
        return np.mean(complementarity_scores)
    
    def _hydrophobic_complementarity(self, points1: List[SurfacePoint],
                                     points2: List[SurfacePoint]) -> float:
        """
        Compute hydrophobic complementarity.
        Hydrophobic patches should match hydrophobic patches.
        """
        if not points1 or not points2:
            return 0.0
        
        complementarity_scores = []
        
        for p1, p2 in zip(points1, points2):
            h1 = p1.hydrophobicity
            h2 = p2.hydrophobicity
            
            # Hydrophobic-hydrophobic or hydrophilic-hydrophilic matching
            comp = 1.0 - abs(h1 - h2)
            complementarity_scores.append(comp)
        
        return np.mean(complementarity_scores)
    
    def _histogram_correlation(self, values1: np.ndarray,
                               values2: np.ndarray,
                               n_bins: int = 20) -> float:
        """Compute correlation between two distributions using histograms."""
        if len(values1) == 0 or len(values2) == 0:
            return 0.0
        
        # Determine common range
        all_values = np.concatenate([values1, values2])
        vmin, vmax = np.percentile(all_values, [1, 99])
        
        if vmax - vmin < 1e-10:
            return 1.0  # Both distributions are essentially constant
        
        bins = np.linspace(vmin, vmax, n_bins + 1)
        
        hist1, _ = np.histogram(values1, bins=bins, density=True)
        hist2, _ = np.histogram(values2, bins=bins, density=True)
        
        # Normalize
        hist1 = hist1 / (np.sum(hist1) + 1e-10)
        hist2 = hist2 / (np.sum(hist2) + 1e-10)
        
        # Bhattacharyya coefficient
        bc = np.sum(np.sqrt(hist1 * hist2))
        
        return bc
    
    def search_similar_patches(self,
                               template_mesh: SurfaceMesh,
                               target_mesh: SurfaceMesh,
                               template_points: List[SurfacePoint],
                               target_points: List[SurfacePoint],
                               template_name: str = "template",
                               target_name: str = "target",
                               similarity_threshold: float = 0.5,
                               top_k: int = 5) -> SimilaritySearchResult:
        """
        Search for similar patches in target protein using template protein patches.
        
        This is the correct similarity search logic:
        1. Generate patches from the template protein
        2. For each template patch, search ALL patches in the target protein
        3. Return ranked matches for each template patch
        
        Args:
            template_mesh: Surface mesh of template protein
            target_mesh: Surface mesh of target protein
            template_points: Surface points of template protein
            target_points: Surface points of target protein
            template_name: Name of template protein
            target_name: Name of target protein
            similarity_threshold: Minimum similarity score to consider a match
            top_k: Number of top matches to return per template patch
        
        Returns:
            SimilaritySearchResult with detailed match information
        """
        # Extract patches from template (these are our query patches)
        template_patches = self._extract_patches(template_mesh, template_points)
        
        # Extract ALL patches from target (these are what we search through)
        # Use more patches for target to ensure good coverage
        original_n_patches = self.n_patches
        self.n_patches = max(self.n_patches * 2, 200)  # More patches for thorough search
        target_patches = self._extract_patches(target_mesh, target_points)
        self.n_patches = original_n_patches  # Restore
        
        # Compute features for all patches
        for patch in template_patches:
            patch.compute_mean_features()
            patch.compute_feature_histogram()
        
        for patch in target_patches:
            patch.compute_mean_features()
            patch.compute_feature_histogram()
        
        # Search: for each template patch, find similar patches in target
        patch_results = []
        all_best_scores = []
        best_overall_score = 0.0
        best_template_idx = -1
        best_target_idx = -1
        
        for t_idx, template_patch in enumerate(template_patches):
            # Compare this template patch against ALL target patches
            matches = []
            
            for tgt_idx, target_patch in enumerate(target_patches):
                # Compute similarity
                overall_sim = self._patch_similarity(template_patch, target_patch)
                geo_sim = self._geometric_similarity(template_patch, target_patch)
                chem_sim = self._chemical_similarity(template_patch, target_patch)
                
                if overall_sim >= similarity_threshold:
                    matches.append((
                        tgt_idx,
                        target_patch.center.copy(),
                        overall_sim,
                        geo_sim,
                        chem_sim
                    ))
            
            # Sort matches by similarity (descending)
            matches.sort(key=lambda x: x[2], reverse=True)
            
            # Keep top_k matches
            top_matches = matches[:top_k]
            
            # Create result for this template patch
            best_match_idx = top_matches[0][0] if top_matches else -1
            best_match_score = top_matches[0][2] if top_matches else 0.0
            best_match_center = top_matches[0][1] if top_matches else None
            
            patch_result = PatchSearchResult(
                template_patch_idx=t_idx,
                template_center=template_patch.center.copy(),
                matches=top_matches,
                best_match_idx=best_match_idx,
                best_match_score=best_match_score,
                best_match_center=best_match_center
            )
            patch_results.append(patch_result)
            
            if best_match_score > 0:
                all_best_scores.append(best_match_score)
            
            # Track overall best
            if best_match_score > best_overall_score:
                best_overall_score = best_match_score
                best_template_idx = t_idx
                best_target_idx = best_match_idx
        
        # Compute residue matches for top results
        residue_matches = []
        for pr in patch_results:
            if pr.best_match_idx >= 0 and pr.best_match_score >= similarity_threshold:
                template_patch = template_patches[pr.template_patch_idx]
                target_patch = target_patches[pr.best_match_idx]
                
                # Get representative residues from each patch
                template_residues = set()
                for pt in template_patch.points:
                    if pt.nearest_residue:
                        template_residues.add(pt.nearest_residue)
                
                target_residues = set()
                for pt in target_patch.points:
                    if pt.nearest_residue:
                        target_residues.add(pt.nearest_residue)
                
                # Add residue pairs
                for t_res in list(template_residues)[:3]:  # Limit to 3 per patch
                    for tgt_res in list(target_residues)[:3]:
                        residue_matches.append((t_res, tgt_res, pr.best_match_score))
        
        # Create final result
        result = SimilaritySearchResult(
            template_object=template_name,
            target_object=target_name,
            n_template_patches=len(template_patches),
            n_target_patches=len(target_patches),
            n_matches_found=sum(1 for pr in patch_results if pr.best_match_score >= similarity_threshold),
            patch_results=patch_results,
            mean_best_similarity=np.mean(all_best_scores) if all_best_scores else 0.0,
            max_similarity=best_overall_score,
            best_template_patch_idx=best_template_idx,
            best_target_patch_idx=best_target_idx,
            best_overall_score=best_overall_score,
            residue_matches=residue_matches
        )
        
        return result


# =============================================================================
# Main Analysis Class
# =============================================================================

class SurfaceSimilarityAnalyzer:
    """
    Main class for surface similarity and complementarity analysis.
    
    Usage:
        analyzer = SurfaceSimilarityAnalyzer()
        
        # Analyze single surface
        result = analyzer.analyze_surface("protein_obj")
        
        # Compare two surfaces
        similarity = analyzer.compare_surfaces("protein1", "protein2")
        
        # Analyze PPI complementarity
        complementarity = analyzer.analyze_complementarity("receptor", "ligand")
    """
    
    def __init__(self, 
                 surface_method: str = "auto",
                 msms_path: Optional[str] = None,
                 apbs_path: Optional[str] = None,
                 patch_radius: float = 12.0):
        """
        Initialize analyzer.
        
        Args:
            surface_method: "msms", "open3d", "edtsurf", or "auto"
            msms_path: Path to MSMS binary
            apbs_path: Path to APBS binary
            patch_radius: Radius for surface patches (Å)
        """
        self.surface_generator = SurfaceGenerator(surface_method, msms_path)
        self.apbs_path = apbs_path
        self.comparator = SurfaceComparator(patch_radius=patch_radius)
        
        # Cache for analyzed surfaces
        self._cache: Dict[str, Tuple[SurfaceMesh, List[SurfacePoint]]] = {}
    
    def analyze_surface(self, obj_name: str = None, 
                       pdb_file: str = None,
                       selection: str = "all",
                       use_apbs: bool = False) -> Tuple[SurfaceMesh, List[SurfacePoint]]:
        """
        Analyze a protein surface.
        
        Args:
            obj_name: PyMOL object name
            pdb_file: PDB file path (alternative to obj_name)
            selection: PyMOL selection (if using obj_name)
            use_apbs: Use APBS for electrostatics
        
        Returns:
            (SurfaceMesh, List[SurfacePoint])
        """
        # Get atoms
        atoms = self._get_atoms(obj_name, pdb_file, selection)
        if not atoms:
            raise ValueError("No atoms found")
        
        print(f"Analyzing surface: {len(atoms)} atoms")
        
        # Generate surface
        mesh = self.surface_generator.generate(atoms)
        mesh.source_object = obj_name or pdb_file or "unknown"
        
        print(f"Generated surface: {mesh.n_vertices} vertices, {mesh.n_faces} faces")
        
        # Extract features
        extractor = FeatureExtractor(atoms, self.apbs_path)
        points = extractor.extract_features(mesh, use_apbs=use_apbs)
        
        mesh.points = points
        
        # Cache result
        cache_key = obj_name or pdb_file or "default"
        self._cache[cache_key] = (mesh, points)
        
        return mesh, points
    
    def compare_surfaces(self, obj1: str, obj2: str,
                        selection1: str = "all",
                        selection2: str = "all") -> SimilarityResult:
        """
        Compare two protein surfaces for similarity.
        
        Args:
            obj1, obj2: PyMOL object names or PDB files
            selection1, selection2: PyMOL selections
        
        Returns:
            SimilarityResult
        """
        # Analyze both surfaces
        mesh1, points1 = self._get_or_analyze(obj1, selection1)
        mesh2, points2 = self._get_or_analyze(obj2, selection2)
        
        # Compare
        result = self.comparator.compute_similarity(mesh1, mesh2, points1, points2)
        
        return result
    
    def analyze_complementarity(self, receptor: str, ligand: str,
                               receptor_sel: str = "all",
                               ligand_sel: str = "all",
                               interface_distance: float = 4.0) -> ComplementarityResult:
        """
        Analyze surface complementarity between receptor and ligand.
        
        Args:
            receptor, ligand: PyMOL object names or PDB files
            receptor_sel, ligand_sel: PyMOL selections
            interface_distance: Distance threshold for interface (Å)
        
        Returns:
            ComplementarityResult
        """
        # Analyze both surfaces
        mesh1, points1 = self._get_or_analyze(receptor, receptor_sel)
        mesh2, points2 = self._get_or_analyze(ligand, ligand_sel)
        
        # Compute complementarity
        result = self.comparator.compute_complementarity(
            mesh1, mesh2, points1, points2, interface_distance)
        
        return result
    
    def _get_or_analyze(self, obj: str, selection: str) -> Tuple[SurfaceMesh, List[SurfacePoint]]:
        """Get cached analysis or perform new analysis."""
        cache_key = f"{obj}_{selection}"
        
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Check if it's a file path
        if os.path.isfile(obj):
            return self.analyze_surface(pdb_file=obj)
        else:
            return self.analyze_surface(obj_name=obj, selection=selection)
    
    def search_similar_surfaces(self,
                                template_obj: str,
                                target_obj: str,
                                template_sel: str = "all",
                                target_sel: str = "all",
                                similarity_threshold: float = 0.5,
                                top_k: int = 5) -> SimilaritySearchResult:
        """
        Search for similar surface patches in target protein using template protein.
        
        This is the correct similarity search logic:
        - Template protein: the protein with known binding site (e.g., CRBN)
        - Target protein: the protein to search for similar sites
        
        The method generates patches from the template and searches for
        similar patches in the target protein.
        
        Args:
            template_obj: PyMOL object name or PDB file for template protein
            target_obj: PyMOL object name or PDB file for target protein
            template_sel: PyMOL selection for template
            target_sel: PyMOL selection for target
            similarity_threshold: Minimum similarity score (0-1) to consider a match
            top_k: Number of top matches to return per template patch
        
        Returns:
            SimilaritySearchResult with detailed match information
        """
        # Analyze both surfaces
        template_mesh, template_points = self._get_or_analyze(template_obj, template_sel)
        target_mesh, target_points = self._get_or_analyze(target_obj, target_sel)
        
        # Perform search
        result = self.comparator.search_similar_patches(
            template_mesh=template_mesh,
            target_mesh=target_mesh,
            template_points=template_points,
            target_points=target_points,
            template_name=template_obj,
            target_name=target_obj,
            similarity_threshold=similarity_threshold,
            top_k=top_k
        )

        return result

    def _normalize_selection(self, obj_name: str, selection: str) -> str:
        """
        Normalize selection string - convert short chain names to full PyMOL selection.

        Examples:
            'A' -> 'chain A'
            'A B' -> 'chain A or chain B'
            'A+B' -> 'chain A or chain B'
            'all' -> 'all'
            'chain A' -> 'chain A' (unchanged)
            'resi 1-100' -> 'resi 1-100' (unchanged)
        """
        if not selection:
            return "all"

        selection = selection.strip()

        # If it's already a valid PyMOL selection keyword, return as-is
        pymol_keywords = ['all', 'chain', 'resi', 'resn', 'name', 'elem', 'hetatm',
                         'polymer', 'solvent', 'organic', 'inorganic', 'and', 'or', 'not']
        first_word = selection.split()[0].lower()
        if first_word in pymol_keywords:
            return selection

        # Check if it looks like chain identifiers (single letters/numbers, possibly separated)
        # Pattern: single chars separated by space, comma, or plus
        import re
        chain_pattern = re.compile(r'^[A-Za-z0-9](\s*[,+\s]\s*[A-Za-z0-9])*$')

        if chain_pattern.match(selection):
            # Split by common separators
            chains = re.split(r'[\s,+]+', selection)
            chains = [c.strip().upper() for c in chains if c.strip()]

            if len(chains) == 1:
                # Single chain - verify it exists
                chain = chains[0]
                if HAS_PYMOL:
                    try:
                        count = cmd.count_atoms(f"{obj_name} and chain {chain}")
                        if count > 0:
                            print(f"💡 Selection '{selection}' interpreted as 'chain {chain}' ({count} atoms)")
                            return f"chain {chain}"
                    except Exception:  # PyMOL count_atoms 可能Failed
                        pass
            else:
                # Multiple chains
                chain_sels = [f"chain {c}" for c in chains]
                combined = " or ".join(chain_sels)
                if HAS_PYMOL:
                    try:
                        count = cmd.count_atoms(f"{obj_name} and ({combined})")
                        if count > 0:
                            print(f"💡 Selection '{selection}' interpreted as '({combined})' ({count} atoms)")
                            return f"({combined})"
                    except Exception:  # PyMOL count_atoms 可能Failed
                        pass

        # Return original selection if no conversion needed
        return selection

    def _get_atoms(self, obj_name: Optional[str],
                   pdb_file: Optional[str],
                   selection: str) -> List[Dict]:
        """Get atom information from PyMOL or PDB file."""
        atoms = []

        # VDW radii
        vdw_radii = {
            'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
            'P': 1.80, 'H': 1.20, 'F': 1.47, 'CL': 1.75,
            'BR': 1.85, 'I': 1.98, 'MG': 1.73, 'CA': 2.31,
            'FE': 2.00, 'ZN': 1.39, 'CU': 1.40, 'MN': 1.61
        }

        if obj_name and HAS_PYMOL:
            # Normalize selection - convert short chain names to full selection
            selection = self._normalize_selection(obj_name, selection)

            # Validate selection first
            full_selection = f"{obj_name} and {selection}"
            try:
                atom_count = cmd.count_atoms(full_selection)
                if atom_count == 0:
                    raise ValueError(
                        f"Selection '{selection}' in object '{obj_name}' contains no atoms.\n"
                        f"Please check:\n"
                        f"  - Is the selection syntax correct? (e.g., 'all', 'chain A', 'A', 'resi 1-100')\n"
                        f"  - Does the object '{obj_name}' exist in PyMOL?"
                    )
            except Exception as e:
                if "Invalid selection" in str(e) or "Selector-Error" in str(e):
                    raise ValueError(
                        f"Invalid selection: '{selection}'\n"
                        f"Common selections: 'all', 'A' (chain A), 'chain A', 'resi 1-100'\n"
                        f"Did you mean 'all' instead of '{selection}'?"
                    )
                raise

            model = cmd.get_model(full_selection)
            for atom in model.atom:
                elem = atom.symbol.upper()
                atoms.append({
                    'coord': np.array(atom.coord),
                    'element': elem,
                    'radius': vdw_radii.get(elem, 1.70),
                    'chain': atom.chain,
                    'resn': atom.resn,
                    'resi': str(atom.resi),
                    'name': atom.name
                })
        elif pdb_file:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        try:
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            elem = line[76:78].strip().upper()
                            if not elem:
                                elem = line[12:14].strip()[0].upper()
                            
                            atoms.append({
                                'coord': np.array([x, y, z]),
                                'element': elem,
                                'radius': vdw_radii.get(elem, 1.70),
                                'chain': line[21:22].strip(),
                                'resn': line[17:20].strip(),
                                'resi': line[22:27].strip(),
                                'name': line[12:16].strip()
                            })
                        except (ValueError, IndexError):  # PDB 行解析可能Failed
                            continue
        
        return atoms
    
    def export_features(self, obj_name: str, output_csv: str,
                       selection: str = "all"):
        """
        Export surface features to CSV.
        
        Args:
            obj_name: PyMOL object name
            output_csv: Output CSV file path
            selection: PyMOL selection
        """
        mesh, points = self._get_or_analyze(obj_name, selection)
        
        import csv
        
        with open(output_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                'X', 'Y', 'Z',
                'Normal_X', 'Normal_Y', 'Normal_Z',
                'Gaussian_Curvature', 'Mean_Curvature',
                'Shape_Index', 'Curvedness',
                'Electrostatic_Potential', 'Hydrophobicity',
                'HBond_Donor_Density', 'HBond_Acceptor_Density',
                'Nearest_Residue'
            ])
            
            for p in points:
                writer.writerow([
                    f"{p.coord[0]:.3f}", f"{p.coord[1]:.3f}", f"{p.coord[2]:.3f}",
                    f"{p.normal[0]:.3f}", f"{p.normal[1]:.3f}", f"{p.normal[2]:.3f}",
                    f"{p.gaussian_curvature:.4f}", f"{p.mean_curvature:.4f}",
                    f"{p.shape_index:.4f}", f"{p.curvedness:.4f}",
                    f"{p.electrostatic_potential:.4f}", f"{p.hydrophobicity:.4f}",
                    f"{p.hbond_donor_density:.4f}", f"{p.hbond_acceptor_density:.4f}",
                    p.nearest_residue or ""
                ])
        
        print(f"Exported {len(points)} surface points to {output_csv}")


# =============================================================================
# PyMOL Command Interface
# =============================================================================

def analyze_surface_similarity(obj1: str, obj2: str = None,
                               selection1: str = "all",
                               selection2: str = "all",
                               output_csv: str = None,
                               patch_radius: float = 12.0) -> Union[SimilarityResult, Tuple[SurfaceMesh, List[SurfacePoint]]]:
    """
    PyMOL command: Analyze surface similarity.
    
    Usage:
        analyze_surface_similarity protein1, protein2
        analyze_surface_similarity protein1  # Single surface analysis
    """
    analyzer = SurfaceSimilarityAnalyzer(patch_radius=patch_radius)
    
    if obj2:
        # Compare two surfaces
        result = analyzer.compare_surfaces(obj1, obj2, selection1, selection2)
        
        print("\n" + "="*50)
        print("Surface Similarity Analysis")
        print("="*50)
        print(f"Object 1: {obj1}")
        print(f"Object 2: {obj2}")
        print(f"\nOverall Similarity: {result.score:.3f}")
        print(f"  Geometric: {result.geometric_similarity:.3f}")
        print(f"  Chemical:  {result.chemical_similarity:.3f}")
        print(f"\nFeature Correlations:")
        print(f"  Shape Index:    {result.shape_index_correlation:.3f}")
        print(f"  Curvature:      {result.curvature_correlation:.3f}")
        print(f"  Electrostatic:  {result.electrostatic_correlation:.3f}")
        print(f"  Hydrophobicity: {result.hydrophobicity_correlation:.3f}")
        print("="*50)
        
        return result
    else:
        # Single surface analysis
        mesh, points = analyzer.analyze_surface(obj_name=obj1, selection=selection1)
        
        if output_csv:
            analyzer.export_features(obj1, output_csv, selection1)
        
        print("\n" + "="*50)
        print("Surface Analysis")
        print("="*50)
        print(f"Object: {obj1}")
        print(f"Vertices: {mesh.n_vertices}")
        print(f"Faces: {mesh.n_faces}")
        print(f"Surface Area: {mesh.surface_area:.1f} Ų")
        print("="*50)
        
        return mesh, points


def analyze_surface_complementarity(receptor: str, ligand: str,
                                    receptor_sel: str = "all",
                                    ligand_sel: str = "all",
                                    interface_distance: float = 4.0,
                                    output_csv: str = None) -> ComplementarityResult:
    """
    PyMOL command: Analyze surface complementarity for PPI.
    
    Usage:
        analyze_surface_complementarity receptor, ligand
    """
    analyzer = SurfaceSimilarityAnalyzer()
    
    result = analyzer.analyze_complementarity(
        receptor, ligand, receptor_sel, ligand_sel, interface_distance)
    
    print("\n" + "="*50)
    print("Surface Complementarity Analysis")
    print("="*50)
    print(f"Receptor: {receptor}")
    print(f"Ligand: {ligand}")
    print(f"\nOverall Complementarity: {result.score:.3f}")
    print(f"  Geometric:     {result.geometric_complementarity:.3f}")
    print(f"  Electrostatic: {result.electrostatic_complementarity:.3f}")
    print(f"  Hydrophobic:   {result.hydrophobic_complementarity:.3f}")
    print(f"\nInterface:")
    print(f"  Area: {result.interface_area:.1f} Ų")
    print(f"  Contacts: {result.n_contacts}")
    print(f"  Residue pairs: {len(result.interface_residues)}")
    print("="*50)
    
    if output_csv:
        import csv
        with open(output_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['Receptor_Residue', 'Ligand_Residue'])
            for r1, r2 in result.interface_residues:
                writer.writerow([r1, r2])
        print(f"Interface residues saved to {output_csv}")
    
    return result


# =============================================================================
# Batch Analysis
# =============================================================================

class BatchSurfaceAnalyzer:
    """
    Batch analysis for multiple protein structures.
    
    Supports:
    - Analyzing multiple PDB files
    - Pairwise similarity matrix
    - Clustering based on surface features
    """
    
    def __init__(self,
                 surface_method: str = "auto",
                 patch_radius: float = 12.0,
                 n_jobs: int = 1):
        """
        Initialize batch analyzer.
        
        Args:
            surface_method: Surface generation method
            patch_radius: Patch radius for comparison
            n_jobs: Number of parallel jobs (1 = sequential)
        """
        self.surface_method = surface_method
        self.patch_radius = patch_radius
        self.n_jobs = n_jobs
        self.analyzer = SurfaceSimilarityAnalyzer(
            surface_method=surface_method,
            patch_radius=patch_radius
        )
        
        # Results cache
        self._surfaces: Dict[str, Tuple[SurfaceMesh, List[SurfacePoint]]] = {}
        self._similarity_matrix: Optional[np.ndarray] = None
        self._structure_names: List[str] = []
    
    def add_structure(self, name: str, obj_name: str = None,
                     pdb_file: str = None, selection: str = "all"):
        """
        Add a structure for batch analysis.
        
        Args:
            name: Unique identifier for this structure
            obj_name: PyMOL object name
            pdb_file: PDB file path
            selection: PyMOL selection
        """
        print(f"Adding structure: {name}")
        
        mesh, points = self.analyzer.analyze_surface(
            obj_name=obj_name,
            pdb_file=pdb_file,
            selection=selection
        )
        
        self._surfaces[name] = (mesh, points)
        if name not in self._structure_names:
            self._structure_names.append(name)
        
        print(f"  Added: {mesh.n_vertices} vertices")
    
    def add_structures_from_directory(self, directory: str,
                                      pattern: str = "*.pdb",
                                      recursive: bool = False):
        """
        Add all PDB files from a directory.
        
        Args:
            directory: Directory path
            pattern: File pattern (glob)
            recursive: Search recursively
        """
        dir_path = Path(directory)
        
        if recursive:
            files = list(dir_path.rglob(pattern))
        else:
            files = list(dir_path.glob(pattern))
        
        print(f"Found {len(files)} files matching '{pattern}'")
        
        for pdb_file in files:
            name = pdb_file.stem
            try:
                self.add_structure(name, pdb_file=str(pdb_file))
            except Exception as e:
                print(f"  Error processing {name}: {e}")
    
    def compute_similarity_matrix(self,
                                  metric: str = "overall") -> np.ndarray:
        """
        Compute pairwise similarity matrix.
        
        Args:
            metric: "overall", "geometric", or "chemical"
        
        Returns:
            np.ndarray: N x N similarity matrix
        """
        n = len(self._structure_names)
        if n < 2:
            raise ValueError("Need at least 2 structures for comparison")
        
        print(f"Computing {n}x{n} similarity matrix...")
        
        matrix = np.zeros((n, n))
        
        for i in range(n):
            matrix[i, i] = 1.0  # Self-similarity
            
            for j in range(i + 1, n):
                name_i = self._structure_names[i]
                name_j = self._structure_names[j]
                
                mesh_i, points_i = self._surfaces[name_i]
                mesh_j, points_j = self._surfaces[name_j]
                
                result = self.analyzer.comparator.compute_similarity(
                    mesh_i, mesh_j, points_i, points_j)
                
                if metric == "geometric":
                    sim = result.geometric_similarity
                elif metric == "chemical":
                    sim = result.chemical_similarity
                else:
                    sim = result.score
                
                matrix[i, j] = sim
                matrix[j, i] = sim
                
                print(f"  {name_i} vs {name_j}: {sim:.3f}")
        
        self._similarity_matrix = matrix
        return matrix
    
    def cluster_structures(self, n_clusters: int = None,
                          method: str = "hierarchical") -> Dict[str, int]:
        """
        Cluster structures based on surface similarity.
        
        Args:
            n_clusters: Number of clusters (auto if None)
            method: "hierarchical" or "kmeans"
        
        Returns:
            Dict[str, int]: {structure_name: cluster_id}
        """
        if self._similarity_matrix is None:
            self.compute_similarity_matrix()
        
        # Convert similarity to distance
        distance_matrix = 1.0 - self._similarity_matrix
        
        try:
            from scipy.cluster.hierarchy import linkage, fcluster
            from scipy.spatial.distance import squareform
            
            # Hierarchical clustering
            condensed = squareform(distance_matrix)
            Z = linkage(condensed, method='average')
            
            if n_clusters is None:
                # Auto-determine number of clusters
                n_clusters = max(2, len(self._structure_names) // 3)
            
            labels = fcluster(Z, n_clusters, criterion='maxclust')
            
            clusters = {}
            for name, label in zip(self._structure_names, labels):
                clusters[name] = int(label)
            
            return clusters
            
        except ImportError:
            print("scipy.cluster not available, returning single cluster")
            return {name: 1 for name in self._structure_names}
    
    def find_most_similar(self, query_name: str, top_k: int = 5) -> List[Tuple[str, float]]:
        """
        Find most similar structures to a query.
        
        Args:
            query_name: Name of query structure
            top_k: Number of results to return
        
        Returns:
            List of (structure_name, similarity_score) tuples
        """
        if self._similarity_matrix is None:
            self.compute_similarity_matrix()
        
        if query_name not in self._structure_names:
            raise ValueError(f"Unknown structure: {query_name}")
        
        query_idx = self._structure_names.index(query_name)
        similarities = self._similarity_matrix[query_idx]
        
        # Sort by similarity (descending), exclude self
        results = []
        for i, sim in enumerate(similarities):
            if i != query_idx:
                results.append((self._structure_names[i], sim))
        
        results.sort(key=lambda x: x[1], reverse=True)
        return results[:top_k]
    
    def export_results(self, output_dir: str):
        """
        Export all results to files.
        
        Args:
            output_dir: Output directory
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Export similarity matrix
        if self._similarity_matrix is not None:
            matrix_file = os.path.join(output_dir, "similarity_matrix.csv")
            with open(matrix_file, 'w') as f:
                # Header
                f.write("," + ",".join(self._structure_names) + "\n")
                # Data
                for i, name in enumerate(self._structure_names):
                    row = [name] + [f"{self._similarity_matrix[i, j]:.4f}"
                                   for j in range(len(self._structure_names))]
                    f.write(",".join(row) + "\n")
            print(f"Saved similarity matrix to {matrix_file}")
        
        # Export individual surface features
        features_dir = os.path.join(output_dir, "features")
        os.makedirs(features_dir, exist_ok=True)
        
        for name, (mesh, points) in self._surfaces.items():
            feature_file = os.path.join(features_dir, f"{name}_features.csv")
            
            with open(feature_file, 'w') as f:
                f.write("X,Y,Z,Shape_Index,Curvedness,ESP,Hydrophobicity\n")
                for p in points:
                    f.write(f"{p.coord[0]:.3f},{p.coord[1]:.3f},{p.coord[2]:.3f},"
                           f"{p.shape_index:.4f},{p.curvedness:.4f},"
                           f"{p.electrostatic_potential:.4f},{p.hydrophobicity:.4f}\n")
        
        print(f"Saved {len(self._surfaces)} feature files to {features_dir}")
    
    def get_summary(self) -> Dict:
        """Get analysis summary."""
        return {
            "n_structures": len(self._structure_names),
            "structure_names": self._structure_names,
            "has_similarity_matrix": self._similarity_matrix is not None,
            "surface_method": self.surface_method,
            "patch_radius": self.patch_radius
        }


def batch_analyze_surfaces(pdb_files: List[str] = None,
                          pdb_directory: str = None,
                          output_dir: str = None,
                          surface_method: str = "auto",
                          patch_radius: float = 12.0) -> BatchSurfaceAnalyzer:
    """
    PyMOL command: Batch analyze multiple protein surfaces.
    
    Usage:
        batch_analyze_surfaces pdb_directory="/path/to/pdbs", output_dir="/path/to/output"
    """
    analyzer = BatchSurfaceAnalyzer(
        surface_method=surface_method,
        patch_radius=patch_radius
    )
    
    if pdb_files:
        for pdb_file in pdb_files:
            name = os.path.splitext(os.path.basename(pdb_file))[0]
            analyzer.add_structure(name, pdb_file=pdb_file)
    
    if pdb_directory:
        analyzer.add_structures_from_directory(pdb_directory)
    
    if len(analyzer._structure_names) >= 2:
        analyzer.compute_similarity_matrix()
    
    if output_dir:
        analyzer.export_results(output_dir)
    
    # Print summary
    summary = analyzer.get_summary()
    print("\n" + "="*50)
    print("Batch Surface Analysis Summary")
    print("="*50)
    print(f"Structures analyzed: {summary['n_structures']}")
    print(f"Surface method: {summary['surface_method']}")
    print(f"Patch radius: {summary['patch_radius']} Å")
    print("="*50)
    
    return analyzer


# Register PyMOL commands
if HAS_PYMOL:
    cmd.extend("analyze_surface_similarity", analyze_surface_similarity)
    cmd.extend("analyze_surface_complementarity", analyze_surface_complementarity)
    cmd.extend("batch_analyze_surfaces", batch_analyze_surfaces)