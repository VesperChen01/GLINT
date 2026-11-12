# -*- coding: utf-8 -*-
"""
Surface Similarity and Complementarity Analysis Module
======================================================
MaSIF-inspired protein surface analysis for GLINT.

Features:
- Surface mesh generation with MSMS
- Geometric features: curvature, shape index, curvedness
- Chemical features: electrostatics, hydrophobicity, H-bond capacity
- Surface similarity comparison
- Surface complementarity analysis for PPI
- Batch analysis for multiple structures

Dependencies:
- Required: numpy, scipy
- Required for surface analysis: MSMS executable
- Optional: APBS (for accurate electrostatics)

Cross-platform support: macOS, Windows, Linux
"""

from __future__ import annotations
import os
import math
import platform
import subprocess
import tempfile
import shutil
import gzip
from typing import List, Dict, Tuple, Optional, Any, Union
from dataclasses import dataclass, field
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from glint.path_utils import find_executable


# Import ligand features module for small molecule support
try:
    from .ligand_features import is_ligand, get_ligand_hydrophobicity, get_ligand_charge
    HAS_LIGAND_FEATURES = True
except ImportError:
    HAS_LIGAND_FEATURES = False

try:
    from .ligand_ec_calculator import PDB2PQRRunner, APBSRunner, DXGrid
    HAS_APBS_RUNNERS = True
except ImportError:
    PDB2PQRRunner = None
    APBSRunner = None
    DXGrid = None
    HAS_APBS_RUNNERS = False
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
    surface_residues: List[str] = field(default_factory=list)
    
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
    directional_electrostatic_complementarity: float = 0.0  # Spatially organized charge matching
    effective_electrostatic_complementarity: float = 0.0  # Charge score used by gate
    hydrophobic_complementarity: float = 0.0  # Hydrophobic matching
    hbond_complementarity: float = 0.0  # Donor/acceptor matching
    
    # Interface analysis
    interface_area: float = 0.0
    n_contacts: int = 0
    gap_volume: float = 0.0
    min_distance: float = 0.0
    analysis_mode: str = "patch"
    n_patch_matches: int = 0
    patch_radius: float = 0.0
    patch_score_std: float = 0.0
    electrostatic_gate_threshold: float = 0.0
    electrostatic_gate_multiplier: float = 1.0
    core_overall_score: float = 0.0
    main_model_name: str = "HybridPatchScore"
    template_class: str = "imid"
    ablation_rankings: Dict[str, List["ComplementarityPatchMatch"]] = field(default_factory=dict)

    # Representative residue pairs from best patch matches
    interface_residues: List[Tuple[str, str]] = field(default_factory=list)
    patch_matches: List["ComplementarityPatchMatch"] = field(default_factory=list)


@dataclass
class ComplementarityPatchMatch:
    """Detailed patch-level complementarity match for visualization and reporting."""
    center1: np.ndarray
    center2: np.ndarray
    overall_score: float
    geometric_score: float
    electrostatic_score: float
    hydrophobic_score: float
    core_overall_score: float = 0.0
    directional_electrostatic_score: float = 0.0
    effective_electrostatic_score: float = 0.0
    hbond_score: float = 0.0
    interface_score: float = 0.0
    focused_interface_score: float = 0.0
    background_score: float = 0.0
    interface_delta: float = 0.0
    fingerprint_score: float = 0.0
    fingerprint_shape_score: float = 0.0
    fingerprint_curvature_score: float = 0.0
    fingerprint_esp_score: float = 0.0
    fingerprint_hydro_score: float = 0.0
    fingerprint_hbond_score: float = 0.0
    fingerprint_hbond_total_score: float = 0.0
    specificity_score: float = 0.0
    target_pocketness_score: float = 0.0
    template_match_score: float = 0.0
    template_only_score: float = 0.0
    pocket_only_score: float = 0.0
    template_specificity_score: float = 0.0
    template_pocket_score: float = 0.0
    template_anchor_distance: float = float("nan")
    electrostatic_gate_multiplier: float = 1.0
    residues1: List[str] = field(default_factory=list)
    residues2: List[str] = field(default_factory=list)


# =============================================================================
# Surface Generation
# =============================================================================

class SurfaceGenerator:
    """
    Generate molecular surface using MSMS.
    
    Supports:
    - MSMS molecular surface
    """
    
    def __init__(self, method: str = "auto"):
        """
        Initialize surface generator.
        
        Args:
            method: "auto" or "msms". Both use MSMS.
        """
        self.method = (method or "auto").lower()
        self.last_requested_method = self.method
        self.last_resolved_method = self.method
        self.msms_path = self._find_msms()

    def _find_msms(self) -> Optional[str]:
        """Find an MSMS executable."""
        resolved = find_executable("msms")
        if resolved and self._is_usable_msms(resolved):
            return resolved

        candidates = [
            "/opt/homebrew/Caskroom/miniconda/base/envs/glint/bin/msms",
            "/opt/homebrew/Caskroom/miniconda/base/bin/msms",
            os.path.expanduser("~/miniconda3/envs/glint/bin/msms"),
            os.path.expanduser("~/anaconda3/envs/glint/bin/msms"),
            os.path.expanduser("~/ADFRsuite/bin/msms"),
            "/usr/local/bin/msms",
            "/opt/homebrew/bin/msms",
        ]
        for candidate in candidates:
            if os.path.exists(candidate) and os.access(candidate, os.X_OK) and self._is_usable_msms(candidate):
                return candidate
        return None

    def _is_usable_msms(self, path: str) -> bool:
        """Return True only if MSMS can be launched on this machine."""
        try:
            result = subprocess.run([path], capture_output=True, text=True, timeout=5)
            combined = f"{result.stdout}\n{result.stderr}".lower()
            if "bad cpu type" in combined or "exec format" in combined:
                return False
            return True
        except OSError:
            return False
        except subprocess.TimeoutExpired:
            return True

    def _resolve_method(self, requested: str) -> str:
        """Resolve a requested method. GLINT surface complementarity uses MSMS only."""
        method = (requested or "auto").lower()

        if method in {"auto", "msms"}:
            if self.msms_path:
                return "msms"
            raise RuntimeError(
                "MSMS is required for surface analysis, but no compatible executable is available. "
                "Install it in the GLINT conda environment with: "
                "conda install -n glint -c bioconda -c conda-forge scripps-msms"
            )

        raise RuntimeError(
            f"Unsupported surface method '{requested}'. GLINT surface complementarity uses MSMS only."
        )
    
    def generate(self, atoms: List[Dict], probe_radius: float = 1.4) -> SurfaceMesh:
        """
        Generate surface mesh from atom coordinates.
        
        Args:
            atoms: List of atom dicts with 'coord', 'element', 'radius' keys
            probe_radius: Solvent probe radius (default 1.4 Å for water)
        
        Returns:
            SurfaceMesh object
        """
        method = self._resolve_method(self.method)
        self.last_requested_method = self.method
        self.last_resolved_method = method

        generators = {
            "msms": self._generate_msms,
        }

        try:
            mesh = generators[method](atoms, probe_radius)
            if mesh.n_vertices == 0:
                raise RuntimeError(f"{method} generated an empty surface")
            mesh.generation_method = method
            return mesh
        except Exception as exc:
            raise RuntimeError(f"MSMS surface generation failed: {exc}")

    def _generate_msms(self, atoms: List[Dict], probe_radius: float) -> SurfaceMesh:
        """Generate a molecular surface with external MSMS."""
        if not self.msms_path:
            raise RuntimeError("MSMS executable not found")

        with tempfile.TemporaryDirectory(prefix="glint_msms_") as work_dir:
            xyzr_path = os.path.join(work_dir, "surface.xyzr")
            out_prefix = os.path.join(work_dir, "surface")

            with open(xyzr_path, "w", encoding="utf-8") as fh:
                for atom in atoms:
                    coord = atom["coord"]
                    radius = atom.get("radius", 1.7)
                    fh.write(f"{coord[0]:.4f} {coord[1]:.4f} {coord[2]:.4f} {radius:.4f}\n")

            cmd_parts = [
                self.msms_path,
                "-if", xyzr_path,
                "-of", out_prefix,
                "-probe_radius", str(probe_radius),
                "-density", "3.0",
            ]
            result = subprocess.run(cmd_parts, capture_output=True, text=True, timeout=180)
            if result.returncode != 0:
                raise RuntimeError(result.stderr.strip() or result.stdout.strip() or "MSMS failed")

            vert_path = out_prefix + ".vert"
            face_path = out_prefix + ".face"
            if not os.path.exists(vert_path) or not os.path.exists(face_path):
                raise RuntimeError("MSMS did not generate .vert/.face files")

            vertices = []
            normals = []
            with open(vert_path, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) < 6:
                        continue
                    try:
                        vals = [float(x) for x in parts[:6]]
                    except ValueError:
                        continue
                    vertices.append(vals[:3])
                    normal = np.array(vals[3:6], dtype=float)
                    norm = np.linalg.norm(normal)
                    normals.append((normal / norm).tolist() if norm > 0 else [0.0, 0.0, 1.0])

            faces = []
            with open(face_path, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    try:
                        face = [int(parts[0]) - 1, int(parts[1]) - 1, int(parts[2]) - 1]
                    except ValueError:
                        continue
                    if min(face) >= 0:
                        faces.append(face)

        if not vertices:
            raise RuntimeError("MSMS generated an empty surface")

        return SurfaceMesh(
            vertices=np.asarray(vertices, dtype=float),
            faces=np.asarray(faces, dtype=int) if faces else np.empty((0, 3), dtype=int),
            normals=np.asarray(normals, dtype=float),
            generation_method="msms"
        )
    
# =============================================================================
# Feature Extraction
# =============================================================================

class FeatureExtractor:
    """Extract geometric and chemical features for surface points."""

    STANDARD_AMINO_ACIDS = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'MSE',
    }
    
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
    
    def __init__(self, atoms: List[Dict], apbs_path: Optional[str] = None,
                 pdb_source: Optional[str] = None):
        """
        Initialize feature extractor.
        
        Args:
            atoms: List of atom dicts with 'coord', 'element', 'resn', 'resi', 'chain', 'name'
            apbs_path: Path to APBS binary for accurate electrostatics (optional)
            pdb_source: PDB file for PDB2PQR/APBS preparation.
        """
        self.atoms = atoms
        self.apbs_path = apbs_path or self._find_apbs()
        self.pdb_source = pdb_source
        
        # Build spatial index
        self.coords = np.array([a['coord'] for a in atoms])
        self.kdtree = cKDTree(self.coords)
        
        # Precompute residue properties
        self._precompute_residue_properties()
    
    def _find_apbs(self) -> Optional[str]:
        """Find APBS binary."""
        resolved = find_executable("apbs")
        if resolved:
            return resolved

        system = platform.system().lower()
        
        search_paths = []
        if system == "darwin":
            search_paths = [
                "/opt/homebrew/Caskroom/miniconda/base/envs/glint/bin/apbs",
                "/opt/homebrew/Caskroom/miniconda/base/bin/apbs",
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
        if use_apbs:
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
        if not self.pdb_source:
            return self._compute_coulomb_electrostatics(mesh)

        try:
            try:
                from .ligand_ec_calculator import PDB2PQRRunner as _PDB2PQRRunner
                from .ligand_ec_calculator import APBSRunner as _APBSRunner
                from .ligand_ec_calculator import DXGrid as _DXGrid
            except Exception as exc:
                print(f"APBS runner import failed ({exc}), using Coulomb approximation")
                return self._compute_coulomb_electrostatics(mesh)

            with tempfile.TemporaryDirectory(prefix="glint_apbs_") as work_dir:
                input_pdb = os.path.join(work_dir, "input.pdb")
                protein_pqr = os.path.join(work_dir, "protein.pqr")
                apbs_prefix = os.path.join(work_dir, "protein_pot")

                if str(self.pdb_source).endswith(".gz"):
                    with gzip.open(self.pdb_source, "rt", encoding="utf-8", errors="ignore") as src:
                        with open(input_pdb, "w", encoding="utf-8") as dst:
                            shutil.copyfileobj(src, dst)
                else:
                    shutil.copyfile(self.pdb_source, input_pdb)

                pdb2pqr = _PDB2PQRRunner()
                print("Running PDB2PQR for surface electrostatics...")
                if not pdb2pqr.run(input_pdb, protein_pqr, ph=7.4):
                    print("PDB2PQR failed, using Coulomb approximation")
                    return self._compute_coulomb_electrostatics(mesh)

                coords = self.coords
                center = coords.mean(axis=0)
                min_coords = coords.min(axis=0)
                max_coords = coords.max(axis=0)
                extent = max_coords - min_coords
                padding = 20.0
                fine_len = extent + 2.0 * padding
                coarse_len = fine_len * 2.0
                grid_spacing = 0.5
                grid_dims = np.ceil(fine_len / grid_spacing).astype(int)
                grid_dims = grid_dims + (1 - grid_dims % 2)
                grid_dims = np.clip(grid_dims, 65, 225)

                apbs = _APBSRunner(apbs_path=self.apbs_path)
                print("Running APBS for surface electrostatics...")
                apbs_input = apbs.generate_input(
                    protein_pqr,
                    apbs_prefix,
                    grid_center=center.tolist(),
                    grid_dims=tuple(grid_dims.tolist()),
                    coarse_len=tuple(coarse_len.tolist()),
                    fine_len=tuple(fine_len.tolist()),
                )
                dx_file = apbs.run(apbs_input, work_dir)
                if not dx_file or not os.path.exists(dx_file):
                    print("APBS failed to generate DX grid, using Coulomb approximation")
                    return self._compute_coulomb_electrostatics(mesh)

                dx_grid = _DXGrid(dx_file)
                potentials = dx_grid.interpolate(mesh.vertices)
                print(
                    "Using APBS electrostatics for surface points "
                    f"(range {np.nanmin(potentials):.2f} to {np.nanmax(potentials):.2f} kT/e)"
                )
                return np.nan_to_num(potentials, nan=0.0, posinf=0.0, neginf=0.0)

        except Exception as exc:
            print(f"APBS electrostatics failed ({exc}), using Coulomb approximation")
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
        point.surface_residues = self._surface_residues_near_point(point.coord, radius=4.5)
        if not point.surface_residues and nearest_atom.get('resn', '') in self.STANDARD_AMINO_ACIDS:
            point.surface_residues = [point.nearest_residue]
        
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

    def _surface_residues_near_point(self, coord: np.ndarray, radius: float = 4.5) -> List[str]:
        """Collect protein residues with atoms close to an MSMS surface point."""
        indices = self.kdtree.query_ball_point(coord, radius)
        residues = set()
        for idx in indices:
            atom = self.atoms[idx]
            resn = atom.get('resn', '')
            if resn not in self.STANDARD_AMINO_ACIDS:
                continue
            residues.add(f"{atom.get('chain', '')}:{resn}:{atom.get('resi', '')}")
        return sorted(residues)


# =============================================================================
# Similarity and Complementarity Analysis
# =============================================================================

class SurfaceComparator:
    """Compare surfaces for similarity and complementarity."""

    COMPLEMENTARITY_WEIGHTS = {
        "geometric": 0.40,
        "hydrophobic": 0.30,
        "hbond": 0.30,
    }
    ELECTROSTATIC_DIRECTIONAL_WEIGHT = 0.40
    ELECTROSTATIC_GATE_THRESHOLD = 0.60
    ELECTROSTATIC_GATE_MIN_MULTIPLIER = 0.35
    
    def __init__(self, patch_radius: float = 9.0, n_patches: int = 100):
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
                                interface_distance: float = 4.0,
                                reference_anchor_coords: Optional[np.ndarray] = None,
                                reference_anchor_radius: float = 8.0,
                                reference_anchor_mode: str = "include",
                                reference_template_center: Optional[np.ndarray] = None,
                                template_class: str = "imid") -> ComplementarityResult:
        """
        Compute patch-search complementarity between two independent surfaces.

        Args:
            mesh1, mesh2: Surface meshes
            points1, points2: Surface points with features
            interface_distance: Kept for backward compatibility; independent
                patch-search mode does not use coordinate-space contacts.

        Returns:
            ComplementarityResult

        Note:
            Object coordinates are not interpreted as a physical complex here.
            Two independently loaded proteins can overlap in PyMOL coordinates
            by accident, so spatial distance is deliberately ignored.
        """
        print("\nPatch-search surface complementarity")
        print("   Mode: independent template patches vs target patches")
        print("   Coordinate-space distances are ignored.")
        if reference_anchor_coords is not None and len(reference_anchor_coords) > 0:
            anchor_mode = "excluding" if reference_anchor_mode == "exclude" else "keeping"
            print(
                f"   Reference anchor: {anchor_mode} template patch centers within "
                f"{reference_anchor_radius:.1f} Å of {len(reference_anchor_coords)} anchor atoms."
            )
        if reference_template_center is not None:
            center = np.asarray(reference_template_center, dtype=float).reshape(-1)
            if center.size == 3:
                print(
                    "   Fixed template patch center requested: "
                    f"[{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}]"
                )
            print(f"   Template class: {self._normalize_template_class(template_class)}")
        return self._compute_patch_complementarity(
            mesh1,
            mesh2,
            points1,
            points2,
            reference_anchor_coords=reference_anchor_coords,
            reference_anchor_radius=reference_anchor_radius,
            reference_anchor_mode=reference_anchor_mode,
            reference_template_center=reference_template_center,
            template_class=template_class,
        )

    def _compute_interface_complementarity(self, points1: List[SurfacePoint],
                                           points2: List[SurfacePoint],
                                           min_distance: float) -> ComplementarityResult:
        """Compute complementarity from nearest surface point pairs at an interface."""
        geo = self._geometric_complementarity(points1, points2)
        esp = self._electrostatic_complementarity(points1, points2)
        directional_esp = esp
        effective_esp = self._effective_electrostatic_complementarity(esp, directional_esp)
        hydro = self._hydrophobic_complementarity(points1, points2)
        hbond = self._hbond_complementarity(points1, points2)
        overall = self._combined_complementarity_score(geo, effective_esp, hydro, hbond)

        residue_pairs = []
        seen = set()
        for p1, p2 in zip(points1, points2):
            pair = (p1.nearest_residue, p2.nearest_residue)
            if pair[0] and pair[1] and pair not in seen:
                seen.add(pair)
                residue_pairs.append(pair)

        return ComplementarityResult(
            score=float(overall),
            geometric_complementarity=float(geo),
            electrostatic_complementarity=float(esp),
            directional_electrostatic_complementarity=float(directional_esp),
            effective_electrostatic_complementarity=float(effective_esp),
            hydrophobic_complementarity=float(hydro),
            hbond_complementarity=float(hbond),
            interface_area=float(len(points1)) * 0.5,
            n_contacts=len(points1),
            min_distance=min_distance,
            analysis_mode="interface",
            n_patch_matches=0,
            patch_radius=self.patch_radius,
            electrostatic_gate_threshold=self.ELECTROSTATIC_GATE_THRESHOLD,
            electrostatic_gate_multiplier=self._electrostatic_gate_multiplier(effective_esp),
            interface_residues=residue_pairs[:50],
            patch_matches=[],
        )

    def _compute_patch_complementarity(self, mesh1: SurfaceMesh, mesh2: SurfaceMesh,
                                       points1: List[SurfacePoint],
                                       points2: List[SurfacePoint],
                                       min_distance: float = float('nan'),
                                       reference_anchor_coords: Optional[np.ndarray] = None,
                                       reference_anchor_radius: float = 8.0,
                                       reference_anchor_mode: str = "include",
                                       reference_template_center: Optional[np.ndarray] = None,
                                       template_class: str = "imid") -> ComplementarityResult:
        """Compute cross-structure complementarity using best matching surface patches."""
        template_class = self._normalize_template_class(template_class)
        patches1 = self._extract_patches(mesh1, points1)
        original_n_patches = self.n_patches
        if reference_template_center is not None:
            fixed_patch = self._extract_patch_at_center(points1, reference_template_center)
            patches1 = [fixed_patch] if fixed_patch is not None else []
            self.n_patches = max(self.n_patches * 4, 400)
            if patches1:
                center = patches1[0].center
                print(
                    "   Fixed template patch: nearest protein surface center "
                    f"[{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}]"
                )
                print(f"   Target global patch sampling: {self.n_patches} centers")
            else:
                print("   WARNING: fixed template patch could not be created.")
        elif reference_anchor_coords is not None and len(reference_anchor_coords) > 0:
            if reference_anchor_mode == "exclude":
                patches1 = self._filter_patches_near_anchor(
                    patches1,
                    reference_anchor_coords,
                    reference_anchor_radius,
                    exclude=True,
                )
            else:
                patches1 = self._extract_patches_near_anchor(
                    points1,
                    reference_anchor_coords,
                    reference_anchor_radius,
                )
            if not patches1:
                print("   WARNING: no template patches remained after reference anchor selection.")
        patches2 = self._extract_patches(mesh2, points2)
        self.n_patches = original_n_patches

        patch_scores = []
        core_overall_scores = []
        geo_scores = []
        esp_scores = []
        directional_esp_scores = []
        effective_esp_scores = []
        hydro_scores = []
        hbond_scores = []
        residue_pairs = []
        best_patch_matches = []
        fixed_template_mode = reference_template_center is not None

        def _build_patch_match(p1: SurfacePatch, p2: SurfacePatch, components) -> ComplementarityPatchMatch:
            (
                main_score,
                core_overall,
                geo,
                esp,
                directional_esp,
                effective_esp,
                hydro,
                hbond,
                interface_score,
                focused_interface_score,
                background_score,
                specificity_score,
                fingerprint_scores,
                target_pocketness,
                template_match_score,
                template_only_score,
                pocket_only_score,
                template_specificity_score,
                template_pocket_score,
            ) = components
            residues1 = self._surface_residues_for_patch(p1)
            residues2 = self._surface_residues_for_patch(p2) if p2 else []
            template_anchor_distance = self._distance_to_anchor(p1.center, reference_anchor_coords)
            return ComplementarityPatchMatch(
                center1=p1.center.copy(),
                center2=p2.center.copy() if p2 is not None else np.zeros(3, dtype=float),
                overall_score=main_score,
                core_overall_score=core_overall,
                geometric_score=geo,
                electrostatic_score=esp,
                hydrophobic_score=hydro,
                directional_electrostatic_score=directional_esp,
                effective_electrostatic_score=effective_esp,
                hbond_score=hbond,
                interface_score=interface_score,
                focused_interface_score=focused_interface_score,
                background_score=background_score,
                interface_delta=interface_score - background_score,
                fingerprint_score=fingerprint_scores.get("combined", 0.0),
                fingerprint_shape_score=fingerprint_scores.get("shape", 0.0),
                fingerprint_curvature_score=fingerprint_scores.get("curvature", 0.0),
                fingerprint_esp_score=fingerprint_scores.get("esp", 0.0),
                fingerprint_hydro_score=fingerprint_scores.get("hydro", 0.0),
                fingerprint_hbond_score=fingerprint_scores.get("hbond", 0.0),
                fingerprint_hbond_total_score=fingerprint_scores.get("hbond_total", 0.0),
                specificity_score=specificity_score,
                target_pocketness_score=target_pocketness,
                template_match_score=template_match_score,
                template_only_score=template_only_score,
                pocket_only_score=pocket_only_score,
                template_specificity_score=template_specificity_score,
                template_pocket_score=template_pocket_score,
                template_anchor_distance=template_anchor_distance,
                electrostatic_gate_multiplier=self._electrostatic_gate_multiplier(effective_esp),
                residues1=residues1,
                residues2=residues2,
            )

        for p1 in patches1:
            best = None
            best_components = None
            best_pair = (None, None)
            best_patch = None
            target_candidate_records = []
            for p2 in patches2:
                geo = self._geometric_patch_complementarity(p1, p2)
                esp = self._electrostatic_patch_complementarity(p1, p2)
                directional_esp = self._directional_electrostatic_patch_complementarity(p1, p2)
                effective_esp = self._effective_electrostatic_complementarity(esp, directional_esp)
                hydro = self._hydrophobic_patch_complementarity(p1, p2)
                hbond = self._hbond_patch_complementarity(p1, p2)
                core_overall = self._combined_complementarity_score(geo, effective_esp, hydro, hbond)
                interface_score = self._interface_diagnostic_score(hydro, hbond, directional_esp)
                focused_interface_score = self._focused_interface_diagnostic_score(hydro, hbond, esp, effective_esp)
                fingerprint_scores = self._patch_fingerprint_complementarity(p1, p2)
                if fixed_template_mode:
                    fixed_interface_score = self._fixed_template_interface_score(
                        geo,
                        esp,
                        fingerprint_scores,
                        template_class=template_class,
                    )
                    background_score = self._fixed_template_background_score(geo)
                    target_pocketness = self._target_pocketness_score(p2)
                    specificity_score = self._specificity_residual_score(fixed_interface_score - background_score)
                    soft_pocket_gate = self._fixed_template_soft_pocket_gate(
                        target_pocketness,
                        template_class=template_class,
                    )
                    template_only_score = fixed_interface_score
                    pocket_only_score = target_pocketness
                    template_specificity_score = float(np.clip(
                        fixed_interface_score * specificity_score,
                        0.0,
                        1.0,
                    ))
                    template_pocket_score = float(np.clip(
                        fixed_interface_score * soft_pocket_gate,
                        0.0,
                        1.0,
                    ))
                    main_score = self._fixed_template_target_patch_score(
                        fixed_interface_score,
                        target_pocketness,
                        specificity_score,
                        template_class=template_class,
                    )
                else:
                    background_score = self._background_diagnostic_score(geo, hydro)
                    target_pocketness = 1.0
                    fixed_interface_score = fingerprint_scores.get("combined", 0.0)
                    specificity_score = self._specificity_residual_score(interface_score - background_score)
                    template_only_score = fixed_interface_score
                    pocket_only_score = target_pocketness
                    template_specificity_score = fixed_interface_score * specificity_score
                    template_pocket_score = fixed_interface_score
                    main_score = self._hybrid_main_patch_score(
                        core_overall,
                        focused_interface_score,
                        fingerprint_scores.get("combined", 0.0),
                        specificity_score,
                    )
                components = (
                    main_score,
                    core_overall,
                    geo,
                    esp,
                    directional_esp,
                    effective_esp,
                    hydro,
                    hbond,
                    interface_score,
                    focused_interface_score,
                    background_score,
                    specificity_score,
                    fingerprint_scores,
                    target_pocketness,
                    fixed_interface_score,
                    template_only_score,
                    pocket_only_score,
                    template_specificity_score,
                    template_pocket_score,
                )
                res1 = p1.points[0].nearest_residue if p1.points else None
                res2 = p2.points[0].nearest_residue if p2.points else None
                target_candidate_records.append((main_score, p2, components, (res1, res2)))
                if best is None or main_score > best:
                    best = main_score
                    best_components = components
                    best_patch = p2

                    best_pair = (res1, res2)

            if best is not None:
                (
                    main_score,
                    core_overall,
                    geo,
                    esp,
                    directional_esp,
                    effective_esp,
                    hydro,
                    hbond,
                    interface_score,
                    focused_interface_score,
                    background_score,
                    specificity_score,
                    fingerprint_scores,
                    target_pocketness,
                    template_match_score,
                    template_only_score,
                    pocket_only_score,
                    template_specificity_score,
                    template_pocket_score,
                ) = best_components
                patch_scores.append(main_score)
                core_overall_scores.append(core_overall)
                geo_scores.append(geo)
                esp_scores.append(esp)
                directional_esp_scores.append(directional_esp)
                effective_esp_scores.append(effective_esp)
                hydro_scores.append(hydro)
                hbond_scores.append(hbond)
                if fixed_template_mode:
                    target_candidate_records.sort(key=lambda item: item[0], reverse=True)
                    target_matches = [
                        _build_patch_match(p1, candidate_patch, components)
                        for _, candidate_patch, components, _ in target_candidate_records
                    ]
                    best_patch_matches.extend(self._diversify_patch_matches(target_matches, max_matches=20))
                else:
                    best_patch_matches.append(_build_patch_match(p1, best_patch, best_components))
                if best_pair[0] and best_pair[1]:
                    residue_pairs.append(best_pair)

        unique_pairs = []
        seen = set()
        for pair in residue_pairs:
            if pair not in seen:
                seen.add(pair)
                unique_pairs.append(pair)

        best_patch_matches.sort(key=lambda x: x.overall_score, reverse=True)
        diverse_patch_matches = self._diversify_patch_matches(best_patch_matches, max_matches=20)
        ablation_rankings = {}
        if fixed_template_mode and best_patch_matches:
            ablation_rankings = {
                "Full_GLINT": self._diversify_patch_matches(
                    sorted(best_patch_matches, key=lambda x: x.overall_score, reverse=True),
                    max_matches=20,
                ),
                "TemplateMatch_only": self._diversify_patch_matches(
                    sorted(best_patch_matches, key=lambda x: x.template_only_score, reverse=True),
                    max_matches=20,
                ),
                "PocketGate_only": self._diversify_patch_matches(
                    sorted(best_patch_matches, key=lambda x: x.pocket_only_score, reverse=True),
                    max_matches=20,
                ),
                "TemplateMatch_x_SpecificityGate": self._diversify_patch_matches(
                    sorted(best_patch_matches, key=lambda x: x.template_specificity_score, reverse=True),
                    max_matches=20,
                ),
                "TemplateMatch_x_SoftPocketGate": self._diversify_patch_matches(
                    sorted(best_patch_matches, key=lambda x: x.template_pocket_score, reverse=True),
                    max_matches=20,
                ),
            }

        return ComplementarityResult(
            score=float(np.mean(patch_scores)) if patch_scores else 0.0,
            geometric_complementarity=float(np.mean(geo_scores)) if geo_scores else 0.0,
            electrostatic_complementarity=float(np.mean(esp_scores)) if esp_scores else 0.0,
            directional_electrostatic_complementarity=float(np.mean(directional_esp_scores)) if directional_esp_scores else 0.0,
            effective_electrostatic_complementarity=float(np.mean(effective_esp_scores)) if effective_esp_scores else 0.0,
            hydrophobic_complementarity=float(np.mean(hydro_scores)) if hydro_scores else 0.0,
            hbond_complementarity=float(np.mean(hbond_scores)) if hbond_scores else 0.0,
            interface_area=0.0,
            n_contacts=0,
            min_distance=min_distance,
            analysis_mode="patch_search",
            n_patch_matches=len(diverse_patch_matches),
            patch_radius=self.patch_radius,
            patch_score_std=float(np.std(patch_scores, ddof=1)) if len(patch_scores) > 1 else 0.0,
            electrostatic_gate_threshold=self.ELECTROSTATIC_GATE_THRESHOLD,
            electrostatic_gate_multiplier=float(np.mean([
                self._electrostatic_gate_multiplier(score) for score in effective_esp_scores
            ])) if effective_esp_scores else 1.0,
            core_overall_score=float(np.mean(core_overall_scores)) if core_overall_scores else 0.0,
            main_model_name="FixedTemplateTargetScore" if fixed_template_mode else "HybridPatchScore",
            template_class=template_class,
            ablation_rankings=ablation_rankings,
            interface_residues=unique_pairs[:20],
            patch_matches=diverse_patch_matches,
        )

    def _filter_patches_near_anchor(self, patches: List[SurfacePatch],
                                    anchor_coords: np.ndarray,
                                    radius: float,
                                    exclude: bool = False) -> List[SurfacePatch]:
        """Keep or exclude patches whose centers are close to at least one anchor atom."""
        anchors = np.asarray(anchor_coords, dtype=float)
        if anchors.ndim != 2 or anchors.shape[1] != 3:
            return patches
        tree = cKDTree(anchors)
        kept = []
        for patch in patches:
            is_near = float(tree.query(patch.center, k=1)[0]) <= radius
            if is_near != exclude:
                kept.append(patch)
        mode = "exclusion" if exclude else "inclusion"
        print(f"   Template patches after anchor {mode} filter: {len(kept)}/{len(patches)}")
        return kept

    def _extract_patches_near_anchor(self, points: List[SurfacePoint],
                                     anchor_coords: np.ndarray,
                                     radius: float) -> List[SurfacePatch]:
        """Extract template patches from surface centers close to anchor atoms."""
        anchors = np.asarray(anchor_coords, dtype=float)
        if anchors.ndim != 2 or anchors.shape[1] != 3 or not points:
            return []

        coords = np.asarray([p.coord for p in points], dtype=float)
        anchor_tree = cKDTree(anchors)
        distances, _ = anchor_tree.query(coords, k=1)
        local_centers = coords[distances <= radius]
        if len(local_centers) == 0:
            return []
        if len(local_centers) > self.n_patches:
            centers = self._farthest_point_sampling(local_centers, self.n_patches)
        else:
            centers = local_centers

        surface_tree = cKDTree(coords)
        patches = []
        for center in centers:
            indices = surface_tree.query_ball_point(center, self.patch_radius)
            patch_points = [points[i] for i in indices]
            if len(patch_points) >= 5:
                patches.append(SurfacePatch(
                    center=center,
                    radius=self.patch_radius,
                    points=patch_points,
                ))
        print(f"   Template anchor-local patches: {len(patches)} centers within {radius:.1f} Å")
        return patches

    def _extract_patch_at_center(self, points: List[SurfacePoint],
                                 requested_center: np.ndarray) -> Optional[SurfacePatch]:
        """Extract one patch centered at the surface point nearest to requested_center."""
        if not points:
            return None
        coords = np.asarray([p.coord for p in points], dtype=float)
        center = np.asarray(requested_center, dtype=float).reshape(-1)
        if center.size != 3:
            return None
        surface_tree = cKDTree(coords)
        _, center_idx = surface_tree.query(center, k=1)
        surface_center = coords[int(center_idx)]
        indices = surface_tree.query_ball_point(surface_center, self.patch_radius)
        patch_points = [points[i] for i in indices]
        if len(patch_points) < 5:
            return None
        return SurfacePatch(
            center=surface_center,
            radius=self.patch_radius,
            points=patch_points,
        )

    def _distance_to_anchor(self, center: np.ndarray,
                            anchor_coords: Optional[np.ndarray]) -> float:
        """Nearest distance from a patch center to anchor atoms."""
        if anchor_coords is None:
            return float("nan")
        anchors = np.asarray(anchor_coords, dtype=float)
        if anchors.ndim != 2 or anchors.shape[1] != 3 or len(anchors) == 0:
            return float("nan")
        return float(np.min(np.linalg.norm(anchors - np.asarray(center, dtype=float), axis=1)))

    def _diversify_patch_matches(self, matches: List[ComplementarityPatchMatch],
                                 max_matches: int = 20) -> List[ComplementarityPatchMatch]:
        """
        Keep high-scoring patch pairs while suppressing near-duplicate locations.

        Patch extraction uses overlapping neighborhoods. Without this step, the top
        ranked complementarity hits can all describe the same physical region.
        """
        if not matches:
            return []

        min_center_distance = max(2.0, self.patch_radius * 0.90)
        selected: List[ComplementarityPatchMatch] = []
        overflow: List[ComplementarityPatchMatch] = []

        for match in matches:
            is_duplicate = False
            for kept in selected:
                dist1 = float(np.linalg.norm(match.center1 - kept.center1))
                dist2 = float(np.linalg.norm(match.center2 - kept.center2))
                # In fixed-template searches all template centers are identical;
                # target-side spatial diversity is what prevents one region from
                # occupying the whole ranked list.
                if dist2 < min_center_distance or (dist1 < min_center_distance and dist2 < min_center_distance):
                    is_duplicate = True
                    break

            if is_duplicate:
                overflow.append(match)
                continue

            selected.append(match)
            if len(selected) >= max_matches:
                return selected

        if len(selected) < max_matches:
            selected.extend(overflow[: max_matches - len(selected)])

        return selected

    def _combined_complementarity_score(self, geo: float, esp: float,
                                        hydro: float, hbond: float) -> float:
        """
        Combine patch component scores for ranking.

        Combine geometry, hydrophobicity, and donor/acceptor matching, then
        use electrostatic complementarity as a gate. This keeps electrostatics
        from being counted twice while still penalizing weak charge matching.
        """
        weights = self.COMPLEMENTARITY_WEIGHTS
        weighted_score = (
            weights["geometric"] * geo
            + weights["hydrophobic"] * hydro
            + weights["hbond"] * hbond
        )
        if esp >= self.ELECTROSTATIC_GATE_THRESHOLD:
            return weighted_score

        return weighted_score * self._electrostatic_gate_multiplier(esp)

    def _effective_electrostatic_complementarity(self, esp: float, directional_esp: float) -> float:
        """Blend distributional and radial electrostatic complementarity for gating."""
        directional_weight = self.ELECTROSTATIC_DIRECTIONAL_WEIGHT
        effective = (1.0 - directional_weight) * esp + directional_weight * directional_esp
        return float(np.clip(effective, 0.0, 1.0))

    def _electrostatic_gate_multiplier(self, esp: float) -> float:
        """Return the score multiplier applied for weak electrostatic complementarity."""
        if esp >= self.ELECTROSTATIC_GATE_THRESHOLD:
            return 1.0

        multiplier = max(
            self.ELECTROSTATIC_GATE_MIN_MULTIPLIER,
            esp / self.ELECTROSTATIC_GATE_THRESHOLD,
        )
        return float(multiplier)

    def _interface_diagnostic_score(self, hydro: float, hbond: float, directional_esp: float) -> float:
        """Chemistry-first patch score used for interface diagnostics, not GUI ranking."""
        return float(np.clip(0.70 * hbond + 0.15 * hydro + 0.15 * directional_esp, 0.0, 1.0))

    def _focused_interface_diagnostic_score(self, hydro: float, hbond: float, esp: float, effective_esp: float) -> float:
        """ESP-focused interface diagnostic for separating real interface patches from background."""
        return float(np.clip(0.45 * hbond + 0.15 * hydro + 0.25 * effective_esp + 0.15 * esp, 0.0, 1.0))

    def _background_diagnostic_score(self, geo: float, hydro: float) -> float:
        """Geometry/hydrophobic background score that can dominate broad surface matches."""
        return float(np.clip(0.60 * geo + 0.40 * hydro, 0.0, 1.0))

    def _fixed_template_background_score(self, geo: float) -> float:
        """
        Background penalty for fixed ligand-pocket template searches.

        Hydrophobic agreement is meaningful for molecular-glue recruitment
        patches, so it should not be treated as non-specific background here.
        Penalize only broad geometry matches, which are common across unrelated
        protein surfaces.
        """
        return float(np.clip(0.70 * geo, 0.0, 1.0))

    def _specificity_residual_score(self, interface_delta: float) -> float:
        """Convert interface-vs-background residual into a bounded weak specificity term."""
        return float(np.clip(0.5 + interface_delta, 0.0, 1.0))

    def _hybrid_main_patch_score(
        self,
        core_overall: float,
        focused_interface: float,
        fingerprint: float,
        specificity: float,
    ) -> float:
        """Primary GUI patch score: conservative hybrid of legacy and spatial evidence."""
        score = (
            0.50 * core_overall
            + 0.25 * focused_interface
            + 0.20 * fingerprint
            + 0.05 * specificity
        )
        return float(np.clip(score, 0.0, 1.0))

    def _normalize_template_class(self, template_class: Optional[str]) -> str:
        """Normalize ligand template class names used by fixed-template scoring."""
        value = str(template_class or "imid").strip().lower()
        if value in {"macro", "macrocycle", "large", "large_glue", "macrocycle_contact"}:
            return "macrocycle"
        return "imid"

    def _fixed_template_interface_score(
        self,
        geo: float,
        esp: float,
        fingerprint_scores: Dict[str, float],
        template_class: str = "imid",
    ) -> float:
        """
        Multi-channel template-matching score for fixed ligand-pocket search.

        This score compares the fixed template pocket against each target patch
        after channel-specific complementarity transforms: shape and ESP are
        inverted, hydrophobicity is matched, and donor/acceptor channels are
        swapped in _patch_fingerprint_complementarity().
        """
        esp_balance = 1.0 - abs(float(esp) - 0.75) / 0.75
        esp_balance = float(np.clip(esp_balance, 0.0, 1.0))
        if self._normalize_template_class(template_class) == "macrocycle":
            score = (
                0.40 * fingerprint_scores.get("hydro", 0.0)
                + 0.20 * fingerprint_scores.get("curvature", 0.0)
                + 0.15 * fingerprint_scores.get("shape", geo)
                + 0.15 * fingerprint_scores.get("esp", 0.0)
                + 0.05 * fingerprint_scores.get("hbond", 0.0)
                + 0.05 * esp_balance
            )
        else:
            score = (
                0.35 * fingerprint_scores.get("hydro", 0.0)
                + 0.25 * fingerprint_scores.get("esp", 0.0)
                + 0.15 * fingerprint_scores.get("hbond", 0.0)
                + 0.10 * fingerprint_scores.get("shape", geo)
                + 0.10 * esp_balance
                + 0.05 * fingerprint_scores.get("curvature", 0.0)
            )
        return float(np.clip(score, 0.0, 1.0))

    def _fixed_template_target_patch_score(
        self,
        interface_score: float,
        pocket_gate: float,
        specificity: float,
        template_class: str = "imid",
    ) -> float:
        """
        Target-side score after the ligand-pocket template has been fixed.

        Fixed-template search is treated as multi-channel surface template
        matching with two gates: target pocketness and interface specificity.
        The multiplicative gates keep flat/background patches from outranking
        true pocket-like regions just because ESP or H-bond histograms match.
        """
        pocket_factor = self._fixed_template_soft_pocket_gate(pocket_gate, template_class)
        score = interface_score * pocket_factor * specificity
        return float(np.clip(score, 0.0, 1.0))

    def _fixed_template_soft_pocket_gate(
        self,
        pocket_gate: float,
        template_class: str = "imid",
    ) -> float:
        """Class-aware softened pocket gate used by fixed-template ranking."""
        if self._normalize_template_class(template_class) == "macrocycle":
            return 0.72 + 0.28 * float(np.clip(pocket_gate, 0.0, 1.0))
        return 0.55 + 0.45 * float(np.clip(pocket_gate, 0.0, 1.0))

    def _target_pocketness_score(self, patch: SurfacePatch) -> float:
        """
        Estimate whether a target patch is pocket-like enough for ranking.

        Shape-index sign can vary with surface normal conventions, so this uses
        both the dominant concavity signal and the stronger of the two signed
        alternatives. It also rewards radial depth and curvature contrast, while
        penalizing nearly flat patches.
        """
        if not patch.points:
            return 0.0

        shape = np.asarray([pt.shape_index for pt in patch.points], dtype=float)
        curv = np.asarray([pt.curvedness for pt in patch.points], dtype=float)
        finite = np.isfinite(shape) & np.isfinite(curv)
        if not np.any(finite):
            return 0.0
        shape = shape[finite]
        curv = curv[finite]

        concave_negative = float(np.mean(np.clip(-shape, 0.0, 1.0)))
        concave_positive = float(np.mean(np.clip(shape, 0.0, 1.0)))
        concavity_score = max(concave_negative, concave_positive)

        mean_curv = float(np.mean(np.abs(curv))) + 1e-8
        curvature_contrast = float(np.clip(np.std(curv) / (mean_curv + np.std(curv) + 1e-8), 0.0, 1.0))

        coords = np.asarray([pt.coord for pt in patch.points], dtype=float)[finite]
        distances = np.linalg.norm(coords - patch.center, axis=1)
        scaled = np.clip(distances / max(float(patch.radius), 1e-8), 0.0, 0.999999)
        center_mask = scaled <= 0.35
        rim_mask = scaled >= 0.70
        radial_depth = 0.0
        if np.any(center_mask) and np.any(rim_mask):
            signed_shape = shape if concave_positive >= concave_negative else -shape
            center_concavity = float(np.mean(np.clip(signed_shape[center_mask], 0.0, 1.0)))
            rim_concavity = float(np.mean(np.clip(signed_shape[rim_mask], 0.0, 1.0)))
            radial_depth = float(np.clip(center_concavity - rim_concavity + 0.5, 0.0, 1.0))

        flat_shape = np.abs(shape) < 0.12
        low_curv = np.abs(curv) < np.percentile(np.abs(curv), 40)
        flatness_penalty = float(np.mean(flat_shape & low_curv))

        raw = (
            2.5 * concavity_score
            + 1.5 * curvature_contrast
            + 1.0 * radial_depth
            - 1.5 * flatness_penalty
            - 1.2
        )
        return float(np.clip(1.0 / (1.0 + np.exp(-raw)), 0.0, 1.0))

    def _patch_fingerprint(self, patch: SurfacePatch, n_bins: int = 4) -> Dict[str, np.ndarray]:
        """Return radial multi-channel fingerprint profiles for one patch."""
        if not patch.points:
            zero = np.zeros(n_bins, dtype=float)
            return {
                "shape": zero.copy(),
                "curvature": zero.copy(),
                "esp": zero.copy(),
                "hydro": zero.copy(),
                "donor": zero.copy(),
                "acceptor": zero.copy(),
                "hbond_total": zero.copy(),
            }

        shape = np.asarray([pt.shape_index for pt in patch.points], dtype=float)
        curvature = np.asarray([pt.curvedness for pt in patch.points], dtype=float)
        esp = np.tanh(np.asarray([pt.electrostatic_potential for pt in patch.points], dtype=float) / 5.0)
        hydro = np.asarray([pt.hydrophobicity for pt in patch.points], dtype=float)
        donor = np.asarray([pt.hbond_donor_density for pt in patch.points], dtype=float)
        acceptor = np.asarray([pt.hbond_acceptor_density for pt in patch.points], dtype=float)

        return {
            "shape": self._patch_radial_profile(patch, shape, n_bins),
            "curvature": self._patch_radial_profile(patch, curvature, n_bins),
            "esp": self._patch_radial_profile(patch, esp, n_bins),
            "hydro": self._patch_radial_profile(patch, hydro, n_bins),
            "donor": self._patch_radial_profile(patch, donor, n_bins),
            "acceptor": self._patch_radial_profile(patch, acceptor, n_bins),
            "hbond_total": self._patch_radial_profile(patch, donor + acceptor, n_bins),
        }

    def _patch_radial_profile(self, patch: SurfacePatch, values: np.ndarray, n_bins: int) -> np.ndarray:
        """Average patch values in radial bins from center to patch edge."""
        if not patch.points:
            return np.zeros(n_bins, dtype=float)
        coords = np.asarray([pt.coord for pt in patch.points], dtype=float)
        distances = np.linalg.norm(coords - patch.center, axis=1)
        if float(np.max(distances, initial=0.0)) <= 1e-8:
            return np.full(n_bins, self._finite_mean_or_zero(values), dtype=float)
        scaled = np.clip(distances / max(float(patch.radius), 1e-8), 0.0, 0.999999)
        bin_ids = np.floor(scaled * n_bins).astype(int)
        profile = np.zeros(n_bins, dtype=float)
        for bin_id in range(n_bins):
            profile[bin_id] = self._finite_mean_or_zero(values[bin_ids == bin_id])
        return profile

    def _patch_fingerprint_complementarity(self, patch1: SurfacePatch, patch2: SurfacePatch) -> Dict[str, float]:
        """Compare patch fingerprints using channel-specific complementarity rules."""
        fp1 = self._patch_fingerprint(patch1)
        fp2 = self._patch_fingerprint(patch2)
        shape = self._profile_match(fp1["shape"], fp2["shape"], opposite=True)
        curvature = self._profile_match(fp1["curvature"], fp2["curvature"], opposite=False)
        esp = self._profile_match(fp1["esp"], fp2["esp"], opposite=True)
        hydro = self._profile_match(fp1["hydro"], fp2["hydro"], opposite=False)
        donor_acceptor = self._profile_match(fp1["donor"], fp2["acceptor"], opposite=False)
        acceptor_donor = self._profile_match(fp1["acceptor"], fp2["donor"], opposite=False)
        hbond = 0.5 * donor_acceptor + 0.5 * acceptor_donor
        hbond_total = self._profile_similarity(fp1["hbond_total"], fp2["hbond_total"])
        combined = (
            0.20 * shape
            + 0.10 * curvature
            + 0.25 * esp
            + 0.15 * hydro
            + 0.25 * hbond
            + 0.05 * hbond_total
        )
        return {
            "shape": float(shape),
            "curvature": float(curvature),
            "esp": float(esp),
            "hydro": float(hydro),
            "hbond": float(hbond),
            "hbond_total": float(hbond_total),
            "combined": float(np.clip(combined, 0.0, 1.0)),
        }

    def _profile_match(self, a: np.ndarray, b: np.ndarray, opposite: bool = False) -> float:
        target_b = -b if opposite else b
        diff = float(np.mean(np.abs(a - target_b))) if a.size and b.size else 2.0
        return float(np.clip(1.0 - diff / 2.0, 0.0, 1.0))

    def _profile_similarity(self, a: np.ndarray, b: np.ndarray) -> float:
        denom = float(np.linalg.norm(a) * np.linalg.norm(b))
        if denom <= 1e-8:
            return 0.5
        return float(np.clip((float(np.dot(a, b)) / denom + 1.0) / 2.0, 0.0, 1.0))

    def _finite_mean_or_zero(self, values: np.ndarray) -> float:
        finite = values[np.isfinite(values)]
        return float(np.mean(finite)) if finite.size else 0.0

    def _surface_residues_for_patch(self, patch: SurfacePatch) -> List[str]:
        """Return protein residues represented by surface points in a patch."""
        residues = set()
        for point in patch.points:
            for residue in getattr(point, "surface_residues", []) or []:
                if residue:
                    residues.add(residue)
            if not getattr(point, "surface_residues", None) and point.nearest_residue:
                residues.add(point.nearest_residue)
        return sorted(residues)
    
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

    def _geometric_patch_complementarity(self, p1: SurfacePatch, p2: SurfacePatch) -> float:
        """Compute patch-level geometric complementarity from feature distributions."""
        si1 = np.array([pt.shape_index for pt in p1.points])
        si2 = -np.array([pt.shape_index for pt in p2.points])
        curv1 = np.array([pt.curvedness for pt in p1.points])
        curv2 = np.array([pt.curvedness for pt in p2.points])
        si_comp = self._histogram_correlation(si1, si2)
        curv_comp = self._histogram_correlation(curv1, curv2)
        return (si_comp + curv_comp) / 2
    
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

    def _electrostatic_patch_complementarity(self, p1: SurfacePatch, p2: SurfacePatch) -> float:
        """Compute patch-level electrostatic complementarity from opposite-charge distributions."""
        esp1 = np.tanh(np.array([pt.electrostatic_potential for pt in p1.points]) / 5.0)
        esp2 = -np.tanh(np.array([pt.electrostatic_potential for pt in p2.points]) / 5.0)
        return self._histogram_correlation(esp1, esp2)

    def _directional_electrostatic_patch_complementarity(self, p1: SurfacePatch, p2: SurfacePatch) -> float:
        """
        Compare spatially organized electrostatic complementarity.

        The regular electrostatic score compares only the distribution of charge
        values. This keeps a lightweight notion of direction by comparing
        center-to-rim radial electrostatic profiles across the two patches.
        """
        profile1 = self._radial_electrostatic_profile(p1)
        profile2 = self._radial_electrostatic_profile(p2)
        valid = np.isfinite(profile1) & np.isfinite(profile2)
        if not np.any(valid):
            return self._electrostatic_patch_complementarity(p1, p2)
        # Opposite potentials in the same radial shell are complementary.
        shell_scores = 1.0 - np.abs(profile1[valid] + profile2[valid]) / 2.0
        return float(np.clip(np.mean(shell_scores), 0.0, 1.0))

    def _radial_electrostatic_profile(self, patch: SurfacePatch, n_bins: int = 4) -> np.ndarray:
        """Mean normalized electrostatic potential from patch center to rim."""
        if not patch.points:
            return np.full(n_bins, np.nan, dtype=float)
        coords = np.asarray([pt.coord for pt in patch.points], dtype=float)
        esp = np.tanh(np.asarray([pt.electrostatic_potential for pt in patch.points], dtype=float) / 5.0)
        distances = np.linalg.norm(coords - patch.center, axis=1)
        if np.max(distances, initial=0.0) <= 1e-8:
            return np.full(n_bins, np.nan, dtype=float)
        scaled = np.clip(distances / self.patch_radius, 0.0, 0.999999)
        bin_ids = np.floor(scaled * n_bins).astype(int)
        profile = np.full(n_bins, np.nan, dtype=float)
        for bin_id in range(n_bins):
            values = esp[bin_ids == bin_id]
            if len(values) > 0:
                profile[bin_id] = float(np.mean(values))
        return profile
    
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

    def _hydrophobic_patch_complementarity(self, p1: SurfacePatch, p2: SurfacePatch) -> float:
        """Compute patch-level hydrophobic complementarity as hydrophobicity matching."""
        hydro1 = np.array([pt.hydrophobicity for pt in p1.points])
        hydro2 = np.array([pt.hydrophobicity for pt in p2.points])
        return self._histogram_correlation(hydro1, hydro2)

    def _hbond_patch_complementarity(self, p1: SurfacePatch, p2: SurfacePatch) -> float:
        """Compute patch-level donor/acceptor complementarity."""
        donor1 = np.array([pt.hbond_donor_density for pt in p1.points])
        acceptor1 = np.array([pt.hbond_acceptor_density for pt in p1.points])
        donor2 = np.array([pt.hbond_donor_density for pt in p2.points])
        acceptor2 = np.array([pt.hbond_acceptor_density for pt in p2.points])

        if (
            np.max(donor1 + acceptor1, initial=0.0) < 1e-8
            or np.max(donor2 + acceptor2, initial=0.0) < 1e-8
        ):
            return 0.5

        directed_scores = [
            self._directed_hbond_complementarity(donor1, acceptor2),
            self._directed_hbond_complementarity(acceptor1, donor2),
        ]
        weighted_sum = sum(score * weight for score, weight in directed_scores)
        total_weight = sum(weight for _, weight in directed_scores)
        return weighted_sum / total_weight if total_weight > 1e-8 else 0.5

    def _directed_hbond_complementarity(self, values1: np.ndarray, values2: np.ndarray) -> Tuple[float, float]:
        """Compare one donor/acceptor direction without rewarding absent signal as perfect."""
        has_signal1 = np.max(values1, initial=0.0) >= 1e-8
        has_signal2 = np.max(values2, initial=0.0) >= 1e-8
        weight = float(np.mean(values1) + np.mean(values2))
        if not has_signal1 and not has_signal2:
            return 0.5, 0.0
        if has_signal1 != has_signal2:
            return 0.0, weight
        return self._histogram_correlation(values1, values2), weight
    
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
                 apbs_path: Optional[str] = None,
                 use_apbs: bool = True,
                 patch_radius: float = 9.0,
                 max_feature_vertices: int = 12000):
        """
        Initialize analyzer.
        
        Args:
            surface_method: "auto" or "msms". Both require MSMS.
            apbs_path: Path to APBS binary
            use_apbs: Prefer APBS electrostatics when available
            patch_radius: Radius for surface patches (Å)
        """
        self.surface_generator = SurfaceGenerator(surface_method)
        self.apbs_path = apbs_path
        self.use_apbs = use_apbs
        self.comparator = SurfaceComparator(patch_radius=patch_radius)
        self.max_feature_vertices = max(1000, int(max_feature_vertices))
        
        # Cache for analyzed surfaces
        self._cache: Dict[str, Tuple[SurfaceMesh, List[SurfacePoint]]] = {}

    def _reduce_mesh_for_feature_extraction(self, mesh: SurfaceMesh) -> Tuple[SurfaceMesh, int]:
        """Downsample very dense meshes before per-vertex feature extraction."""
        n_vertices = len(mesh.vertices)
        if n_vertices <= self.max_feature_vertices:
            return mesh, n_vertices

        sample_indices = np.linspace(0, n_vertices - 1, self.max_feature_vertices, dtype=int)
        reduced_mesh = SurfaceMesh(
            vertices=mesh.vertices[sample_indices].copy(),
            faces=np.empty((0, 3), dtype=int),
            normals=mesh.normals[sample_indices].copy(),
            source_object=mesh.source_object,
            generation_method=mesh.generation_method,
        )
        return reduced_mesh, len(sample_indices)
    
    def analyze_surface(self, obj_name: str = None, 
                       pdb_file: str = None,
                       selection: str = "all",
                       use_apbs: Optional[bool] = None) -> Tuple[SurfaceMesh, List[SurfacePoint]]:
        """
        Analyze a protein surface.
        
        Args:
            obj_name: PyMOL object name
            pdb_file: PDB file path (alternative to obj_name)
            selection: PyMOL selection (if using obj_name)
            use_apbs: Use APBS for electrostatics. If None, use analyzer default.
        
        Returns:
            (SurfaceMesh, List[SurfacePoint])
        """
        # Get atoms
        atoms = self._get_atoms(obj_name, pdb_file, selection)
        if not atoms:
            raise ValueError("No atoms found")
        
        print(f"Analyzing surface: {len(atoms)} atoms")
        if use_apbs is None:
            use_apbs = self.use_apbs
        pdb_source = pdb_file
        temp_pdb_source = None
        if use_apbs and obj_name and HAS_PYMOL:
            normalized_selection = self._normalize_selection(obj_name, selection)
            full_selection = obj_name if normalized_selection == "all" else f"({obj_name}) and ({normalized_selection})"
            tmp = tempfile.NamedTemporaryFile(prefix="glint_surface_", suffix=".pdb", delete=False)
            temp_pdb_source = tmp.name
            tmp.close()
            cmd.save(temp_pdb_source, full_selection)
            pdb_source = temp_pdb_source
        
        # Generate surface
        mesh = self.surface_generator.generate(atoms)
        mesh.source_object = obj_name or pdb_file or "unknown"
        
        print(
            f"Generated surface: {mesh.n_vertices} vertices, {mesh.n_faces} faces "
            f"(method={mesh.generation_method})"
        )

        feature_mesh, feature_vertex_count = self._reduce_mesh_for_feature_extraction(mesh)
        if feature_vertex_count != mesh.n_vertices:
            print(
                f"Downsampling feature extraction mesh: {mesh.n_vertices} -> "
                f"{feature_vertex_count} vertices"
            )
        
        # Extract features
        try:
            extractor = FeatureExtractor(atoms, self.apbs_path, pdb_source=pdb_source)
            points = extractor.extract_features(feature_mesh, use_apbs=use_apbs)
        finally:
            if temp_pdb_source and os.path.exists(temp_pdb_source):
                try:
                    os.remove(temp_pdb_source)
                except OSError:
                    pass
        
        mesh.points = points
        mesh.feature_vertex_count = feature_vertex_count
        
        # Cache result
        cache_key = f"{obj_name or pdb_file or 'default'}|{selection}|apbs={int(bool(use_apbs))}"
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
                               interface_distance: float = 4.0,
                               receptor_anchor_resn: Optional[str] = None,
                               receptor_anchor_coords: Optional[np.ndarray] = None,
                               receptor_anchor_radius: float = 8.0,
                               receptor_anchor_mode: str = "include",
                               receptor_template_center: Optional[np.ndarray] = None,
                               template_class: str = "imid") -> ComplementarityResult:
        """
        Analyze surface complementarity between receptor and ligand.
        
        Args:
            receptor, ligand: PyMOL object names or PDB files
            receptor_sel, ligand_sel: PyMOL selections
            interface_distance: Backward-compatible argument; patch mode does not
                use a direct-contact threshold.
        
        Returns:
            ComplementarityResult
        """
        # Analyze both surfaces
        mesh1, points1 = self._get_or_analyze(receptor, receptor_sel)
        mesh2, points2 = self._get_or_analyze(ligand, ligand_sel)
        anchor_coords = None
        if receptor_anchor_resn:
            anchor_coords = self._anchor_atom_coords(receptor, receptor_sel, receptor_anchor_resn)
            if anchor_coords.size == 0:
                print(f"WARNING: no anchor atoms found for receptor resn {receptor_anchor_resn!r}.")
        if receptor_anchor_coords is not None:
            supplied_coords = np.asarray(receptor_anchor_coords, dtype=float)
            if supplied_coords.ndim == 2 and supplied_coords.shape[1] == 3 and len(supplied_coords) > 0:
                if anchor_coords is not None and len(anchor_coords) > 0:
                    anchor_coords = np.vstack([anchor_coords, supplied_coords])
                else:
                    anchor_coords = supplied_coords
            else:
                print("WARNING: supplied receptor_anchor_coords had invalid shape; ignoring.")
        
        # Compute complementarity
        result = self.comparator.compute_complementarity(
            mesh1,
            mesh2,
            points1,
            points2,
            interface_distance,
            reference_anchor_coords=anchor_coords,
            reference_anchor_radius=receptor_anchor_radius,
            reference_anchor_mode=receptor_anchor_mode,
            reference_template_center=receptor_template_center,
            template_class=template_class,
        )
        
        return result

    def _anchor_atom_coords(self, obj: str, selection: str, resn: str) -> np.ndarray:
        """Return atom coordinates for a residue name used as an anchor."""
        if os.path.isfile(obj):
            atoms = self._get_atoms(None, obj, selection)
        else:
            atoms = self._get_atoms(obj, None, selection)
        wanted = resn.strip().upper()
        coords = [
            atom["coord"] for atom in atoms
            if str(atom.get("resn", "")).strip().upper() == wanted
        ]
        return np.asarray(coords, dtype=float)
    
    def _get_or_analyze(self, obj: str, selection: str) -> Tuple[SurfaceMesh, List[SurfacePoint]]:
        """Get cached analysis or perform new analysis."""
        cache_key = f"{obj}|{selection}|apbs={int(bool(self.use_apbs))}"
        
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Check if it's a file path
        if os.path.isfile(obj):
            return self.analyze_surface(pdb_file=obj, use_apbs=self.use_apbs)
        else:
            return self.analyze_surface(obj_name=obj, selection=selection, use_apbs=self.use_apbs)
    
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
            opener = gzip.open if str(pdb_file).endswith(".gz") else open
            with opener(pdb_file, 'rt') as f:
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
                               patch_radius: float = 9.0) -> Union[SimilarityResult, Tuple[SurfaceMesh, List[SurfacePoint]]]:
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
    print(f"  H-bond D/A:    {result.hbond_complementarity:.3f}")
    print(f"\nPatch-based matching:")
    print(f"  Patch matches: {result.n_patch_matches}")
    print(f"  Patch radius: {result.patch_radius:.1f} Å")
    print(f"  Minimum surface distance: {result.min_distance:.1f} Å")
    print("="*50)
    
    if output_csv:
        import csv
        with open(output_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Rank', 'Overall_Score', 'Geometric_Score',
                'Electrostatic_Score', 'Hydrophobic_Score', 'Hbond_DA_Score',
                'Interface_Score', 'Focused_Interface_Score',
                'Background_Score', 'Interface_Delta',
                'Template_Match_Score', 'Target_Pocketness_Gate',
                'Specificity_Score',
                'Receptor_Center_X', 'Receptor_Center_Y', 'Receptor_Center_Z',
                'Ligand_Center_X', 'Ligand_Center_Y', 'Ligand_Center_Z',
                'Receptor_Patch_Residues', 'Ligand_Patch_Residues'
            ])
            for idx, match in enumerate(result.patch_matches, start=1):
                writer.writerow([
                    idx,
                    f"{match.overall_score:.4f}",
                    f"{match.geometric_score:.4f}",
                    f"{match.electrostatic_score:.4f}",
                    f"{match.hydrophobic_score:.4f}",
                    f"{match.hbond_score:.4f}",
                    f"{getattr(match, 'interface_score', 0.0):.4f}",
                    f"{getattr(match, 'focused_interface_score', 0.0):.4f}",
                    f"{getattr(match, 'background_score', 0.0):.4f}",
                    f"{getattr(match, 'interface_delta', 0.0):.4f}",
                    f"{getattr(match, 'template_match_score', 0.0):.4f}",
                    f"{getattr(match, 'target_pocketness_score', 0.0):.4f}",
                    f"{getattr(match, 'specificity_score', 0.0):.4f}",
                    f"{match.center1[0]:.3f}",
                    f"{match.center1[1]:.3f}",
                    f"{match.center1[2]:.3f}",
                    f"{match.center2[0]:.3f}",
                    f"{match.center2[1]:.3f}",
                    f"{match.center2[2]:.3f}",
                    ";".join(match.residues1),
                    ";".join(match.residues2),
                ])
        print(f"Patch complementarity results saved to {output_csv}")
    
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
                 patch_radius: float = 9.0,
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
                          patch_radius: float = 9.0) -> BatchSurfaceAnalyzer:
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
