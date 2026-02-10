# -*- coding: utf-8 -*-
"""
Worker threads for GLINT GUI.

This module provides QThread-based worker classes for running
long-running analysis tasks in the background without blocking
the GUI.

Each worker emits signals for:
- progress: Status updates during execution
- finished: Results when complete
- error: Error messages if something goes wrong
"""
import os
import traceback
from typing import Optional, List, Tuple, Dict, Any, Union

from .qt_adapter import QThread, Signal as pyqtSignal

from .utils import (
    get_lang, t,
    analyze_pdb_interactions,
    analyze_protein_nucleic_interactions,
    find_crbn_g_motif
)


class AnalysisWorker(QThread):
    """
    Worker thread for general interaction analysis.
    
    Analyzes protein-protein interactions in a structure and optionally
    exports results to CSV.
    
    Signals:
        progress(str): Emitted with status messages during analysis
        finished(list): Emitted with list of interactions when complete
        error(str): Emitted with error message if analysis fails
    """
    progress = pyqtSignal(str)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)
    
    def __init__(self, obj_name: str, pdb_file: Optional[str], output_csv: Optional[str]):
        """
        Initialize the AnalysisWorker.
        
        Args:
            obj_name: Name of the PyMOL object to analyze
            pdb_file: Optional path to PDB file (uses PyMOL object if None)
            output_csv: Optional path for CSV output
        """
        super().__init__()
        self.obj_name = obj_name
        self.pdb_file = pdb_file
        self.output_csv = output_csv
        
    def run(self) -> None:
        """Execute the analysis in a background thread."""
        try:
            self.progress.emit(t("log_start"))
            
            if analyze_pdb_interactions is None:
                raise RuntimeError("analyze_pdb_interactions not available. Check module imports.")
            
            interactions = analyze_pdb_interactions(
                obj_name=self.obj_name,
                pdb_file=self.pdb_file,
                only_between_chains=True,
                output_csv=self.output_csv,
                auto_highlight=True,
            )
            self.progress.emit(t("log_done").format(n=len(interactions)))
            self.finished.emit(interactions)
        except Exception as e:
            self.error.emit(f"{str(e)}\n{traceback.format_exc()}")

class PNAnalysisWorker(QThread):
    """
    Worker thread for Protein-Nucleic Acid interaction analysis.
    
    Analyzes interactions between protein and nucleic acid chains
    in a structure.
    
    Signals:
        progress(str): Emitted with status messages during analysis
        finished(list): Emitted with list of interactions when complete
        error(str): Emitted with error message if analysis fails
    """
    progress = pyqtSignal(str)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)
    
    def __init__(
        self,
        obj_name: str,
        nucleic_chains: List[str],
        protein_chains: List[str],
        output_csv: Optional[str],
        distance_cutoff: float,
        pdb_file: Optional[str] = None
    ):
        """
        Initialize the PNAnalysisWorker.
        
        Args:
            obj_name: Name of the PyMOL object to analyze
            nucleic_chains: List of chain IDs for nucleic acid
            protein_chains: List of chain IDs for protein
            output_csv: Optional path for CSV output
            distance_cutoff: Distance threshold for interactions (Å)
            pdb_file: Optional path to PDB file
        """
        super().__init__()
        self.obj_name = obj_name
        self.nucleic_chains = nucleic_chains
        self.protein_chains = protein_chains
        self.output_csv = output_csv
        self.distance_cutoff = distance_cutoff
        self.pdb_file = pdb_file
        
    def run(self) -> None:
        """Execute the analysis in a background thread."""
        try:
            self.progress.emit("Starting Protein-Nucleic Acid Analysis...")
            
            if analyze_protein_nucleic_interactions is None:
                raise RuntimeError("analyze_protein_nucleic_interactions not available.")
            
            result = analyze_protein_nucleic_interactions(
                obj_name=self.obj_name,
                nucleic_chains=self.nucleic_chains,
                protein_chains=self.protein_chains,
                output_csv=self.output_csv,
                distance_cutoff=self.distance_cutoff,
                pdb_file=self.pdb_file
            )
            interactions = result.get("interactions", [])
            self.progress.emit(f"Analysis complete. Found {len(interactions)} interactions.")
            self.finished.emit(interactions)
        except Exception as e:
            self.error.emit(f"{str(e)}\n{traceback.format_exc()}")

class GMotifWorker(QThread):
    """
    Worker thread for G-Motif (CRBN G-loop) detection.
    
    Detects G-motif patterns in protein structures that may be
    suitable for molecular glue-mediated degradation.
    
    Signals:
        progress(str): Emitted with status messages during detection
        finished(list, str, object, object): Emitted with (hits, csv_path, coords_data, surface_info)
        error(str): Emitted with error message if detection fails
    """
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, str, object, object)  # (hits, csv_path, coords_data, surface_info)
    error = pyqtSignal(str)
    
    def __init__(
        self,
        obj_name: str,
        pdb_file: Optional[str],
        rmsd: float,
        require_gly_pos: str,
        out_csv: Optional[str],
        template_mode: str,
        template_sel: Optional[str],
        template_builtin: Optional[str],
        highlight_surface: bool = False,
        export_coords: bool = False,
        exclude_proline: bool = False,
        check_surface_exposure: bool = True,
        min_sasa: float = 15.0
    ):
        """
        Initialize the GMotifWorker.
        
        Args:
            obj_name: Name of the PyMOL object to analyze
            pdb_file: Optional path to PDB file
            rmsd: RMSD cutoff for template matching (Å)
            require_gly_pos: Glycine position requirement ("6,3", "6", "3", or None)
            out_csv: Optional path for CSV output
            template_mode: Template selection mode ('builtin' or 'selection')
            template_sel: PyMOL selection for custom template
            template_builtin: Name of built-in template (GSPT1, CK1α, VAV1)
            highlight_surface: Whether to highlight G-loop surface
            export_coords: Whether to export coordinates
            exclude_proline: Whether to exclude proline-containing loops (default: False)
            check_surface_exposure: Whether to check surface accessibility
            min_sasa: Minimum SASA per residue (Å²)
        """
        super().__init__()
        self.obj_name = obj_name
        self.pdb_file = pdb_file
        self.rmsd = rmsd
        self.require_gly_pos = require_gly_pos
        self.out_csv = out_csv
        self.template_mode = template_mode
        self.template_sel = template_sel
        self.template_builtin = template_builtin
        self.highlight_surface = highlight_surface
        self.export_coords = export_coords
        self.exclude_proline = exclude_proline
        self.check_surface_exposure = check_surface_exposure
        self.min_sasa = min_sasa
    def run(self):
        try:
            if find_crbn_g_motif is None:
                raise RuntimeError("find_crbn_g_motif not found; ensure g_motif_analyzer.py exists in the plugin directory.")
            self.progress.emit("[G-Motif] Detecting...")
            out_csv_path = self.out_csv
            if not out_csv_path:
                import tempfile, os
                fd, out_csv_path = tempfile.mkstemp(suffix="_gmotif.csv"); os.close(fd)
            
            # Call find_crbn_g_motif with surface highlight and coordinate export support
            result = find_crbn_g_motif(
                obj_name=self.obj_name,
                pdb_file=self.pdb_file,
                template_mode=self.template_mode,
                template_sel=self.template_sel,
                template_builtin=self.template_builtin,
                rmsd_cutoff=float(self.rmsd),
                out_csv=out_csv_path,
                auto_highlight=1,  # Enable auto-highlight
                require_gly_pos=self.require_gly_pos,
                exclude_proline=bool(self.exclude_proline),
                check_surface_exposure=bool(self.check_surface_exposure),
                min_sasa_per_residue=float(self.min_sasa),
                highlight_surface=bool(self.highlight_surface),
                export_coords=bool(self.export_coords),
            )
            
            # Process return result (GMotifResult object)
            if hasattr(result, 'hits'):
                # New version returns GMotifResult object
                hits = list(result.hits) if result.hits else []
                coords_data = result.coordinates
                surface_info = result.surface_info
            else:
                # Old version returns list (backward compatible)
                hits = list(result) if result else []
                coords_data = None
                surface_info = None
            
            self.progress.emit(f"[G-Motif] Done, {len(hits)} hits")
            self.finished.emit(hits, out_csv_path, coords_data, surface_info)
        except Exception as e:
            self.error.emit(str(e))



class SurfaceAnalysisWorker(QThread):
    """Worker for protein surface analysis"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, str) # (patches, out_csv)
    error = pyqtSignal(str)

    def __init__(self, obj_name: str, out_csv: str, use_apbs: bool = True, ph: float = 7.4):
        """
        Initialize the SurfaceAnalysisWorker.
        
        Args:
            obj_name: Name of the PyMOL object to analyze
            out_csv: Optional path for CSV output
            use_apbs: Use APBS/PDB2PQR for accurate electrostatics (default: True)
            ph: pH for PDB2PQR protonation state (default: 7.4)
        """
        super().__init__()
        self.obj_name = obj_name
        self.out_csv = out_csv
        self.use_apbs = use_apbs
        self.ph = ph
        
    def run(self):
        try:
            self.progress.emit(f"Starting surface analysis for {self.obj_name}...")
            if self.use_apbs:
                self.progress.emit(f"  Using APBS/PDB2PQR for accurate electrostatics (pH={self.ph})")
            else:
                self.progress.emit("  Using residue-based electrostatic approximation")
            
            # Late import to avoid circular dependencies if any
            from ..protein_surface_analyzer import SurfaceAnalyzer
            
            analyzer = SurfaceAnalyzer(
                self.obj_name, 
                use_apbs=self.use_apbs, 
                ph=self.ph
            )
            patches = analyzer.analyze()
            
            self.progress.emit(f"Analysis complete. Found {len(patches)} patches.")
            
            # Export to CSV if requested
            if self.out_csv:
                import csv
                with open(self.out_csv, 'w', newline='', encoding='utf-8') as f:
                    fieldnames = ['id', 'type', 'area', 'score', 'center', 'avg_potential', 'residues']
                    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
                    writer.writeheader()
                    for p in patches:
                        writer.writerow(p.to_dict())
                self.progress.emit(f"Report saved to {self.out_csv}")
                
            self.finished.emit(patches, self.out_csv)
            
        except Exception as e:
            self.error.emit(str(e))


class SurfaceSimilarityWorker(QThread):
    """Worker for surface similarity and complementarity analysis"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(object, str)  # (result, out_csv)
    error = pyqtSignal(str)
    
    def __init__(self, obj1: str, obj2: str = None,
                 selection1: str = "all", selection2: str = "all",
                 analysis_type: str = "similarity",  # "similarity", "complementarity", or "search"
                 patch_radius: float = 12.0,
                 interface_distance: float = 4.0,
                 out_csv: str = None,
                 surface_method: str = "auto",
                 use_apbs: bool = False,
                 similarity_threshold: float = 0.5,
                 top_k: int = 5):
        super().__init__()
        self.obj1 = obj1
        self.obj2 = obj2
        self.selection1 = selection1
        self.selection2 = selection2
        self.analysis_type = analysis_type
        self.patch_radius = patch_radius
        self.interface_distance = interface_distance
        self.out_csv = out_csv
        self.surface_method = surface_method
        self.use_apbs = use_apbs
        self.similarity_threshold = similarity_threshold
        self.top_k = top_k
        
    def run(self):
        try:
            from ..surface_similarity import SurfaceSimilarityAnalyzer
            
            self.progress.emit(f"Initializing surface analyzer (method={self.surface_method})...")
            
            analyzer = SurfaceSimilarityAnalyzer(
                surface_method=self.surface_method,
                patch_radius=self.patch_radius
            )
            
            if self.analysis_type == "complementarity" and self.obj2:
                # Complementarity analysis
                self.progress.emit(f"Analyzing surface complementarity: {self.obj1} vs {self.obj2}...")
                
                result = analyzer.analyze_complementarity(
                    self.obj1, self.obj2,
                    self.selection1, self.selection2,
                    self.interface_distance
                )
                
                self.progress.emit(f"Complementarity analysis complete. Score: {result.score:.3f}")
                
                # Export interface residues if CSV requested
                if self.out_csv:
                    import csv
                    with open(self.out_csv, 'w', newline='', encoding='utf-8') as f:
                        writer = csv.writer(f)
                        writer.writerow(['Metric', 'Value'])
                        writer.writerow(['Overall_Score', f"{result.score:.4f}"])
                        writer.writerow(['Geometric_Complementarity', f"{result.geometric_complementarity:.4f}"])
                        writer.writerow(['Electrostatic_Complementarity', f"{result.electrostatic_complementarity:.4f}"])
                        writer.writerow(['Hydrophobic_Complementarity', f"{result.hydrophobic_complementarity:.4f}"])
                        writer.writerow(['Interface_Area', f"{result.interface_area:.2f}"])
                        writer.writerow(['N_Contacts', result.n_contacts])
                        writer.writerow(['', ''])
                        writer.writerow(['Receptor_Residue', 'Ligand_Residue'])
                        for r1, r2 in result.interface_residues:
                            writer.writerow([r1, r2])
                    self.progress.emit(f"Results saved to {self.out_csv}")
            
            elif self.analysis_type == "search" and self.obj2:
                # Similarity search: use obj1 as template, search in obj2
                self.progress.emit(f"Searching similar patches: template={self.obj1}, target={self.obj2}...")
                self.progress.emit(f"  Template selection: {self.selection1}")
                self.progress.emit(f"  Target selection: {self.selection2}")
                self.progress.emit(f"  Similarity threshold: {self.similarity_threshold}")
                
                result = analyzer.search_similar_surfaces(
                    template_obj=self.obj1,
                    target_obj=self.obj2,
                    template_sel=self.selection1,
                    target_sel=self.selection2,
                    similarity_threshold=self.similarity_threshold,
                    top_k=self.top_k
                )
                
                self.progress.emit(f"Search complete:")
                self.progress.emit(f"  Template patches: {result.n_template_patches}")
                self.progress.emit(f"  Target patches: {result.n_target_patches}")
                self.progress.emit(f"  Matches found: {result.n_matches_found}")
                self.progress.emit(f"  Best similarity: {result.max_similarity:.3f}")
                self.progress.emit(f"  Mean best similarity: {result.mean_best_similarity:.3f}")
                
                # Export results if CSV requested
                if self.out_csv:
                    import csv
                    with open(self.out_csv, 'w', newline='', encoding='utf-8') as f:
                        writer = csv.writer(f)
                        # Summary section
                        writer.writerow(['=== Similarity Search Results ===', ''])
                        writer.writerow(['Template', result.template_object])
                        writer.writerow(['Target', result.target_object])
                        writer.writerow(['Template_Patches', result.n_template_patches])
                        writer.writerow(['Target_Patches', result.n_target_patches])
                        writer.writerow(['Matches_Found', result.n_matches_found])
                        writer.writerow(['Max_Similarity', f"{result.max_similarity:.4f}"])
                        writer.writerow(['Mean_Best_Similarity', f"{result.mean_best_similarity:.4f}"])
                        writer.writerow(['', ''])
                        
                        # Per-patch results
                        writer.writerow(['=== Per-Patch Results ===', ''])
                        writer.writerow(['Template_Patch_Idx', 'Template_Center_X', 'Template_Center_Y', 'Template_Center_Z',
                                        'Best_Match_Idx', 'Best_Match_Score', 'Match_Center_X', 'Match_Center_Y', 'Match_Center_Z'])
                        for pr in result.patch_results:
                            if pr.best_match_score >= self.similarity_threshold:
                                tc = pr.template_center
                                mc = pr.best_match_center if pr.best_match_center is not None else [0, 0, 0]
                                writer.writerow([
                                    pr.template_patch_idx,
                                    f"{tc[0]:.2f}", f"{tc[1]:.2f}", f"{tc[2]:.2f}",
                                    pr.best_match_idx,
                                    f"{pr.best_match_score:.4f}",
                                    f"{mc[0]:.2f}", f"{mc[1]:.2f}", f"{mc[2]:.2f}"
                                ])
                        
                        writer.writerow(['', ''])
                        writer.writerow(['=== Residue Matches ===', ''])
                        writer.writerow(['Template_Residue', 'Target_Residue', 'Similarity'])
                        for t_res, tgt_res, sim in result.residue_matches:
                            writer.writerow([t_res, tgt_res, f"{sim:.4f}"])
                    
                    self.progress.emit(f"Results saved to {self.out_csv}")
                    
            elif self.obj2:
                # Similarity comparison (legacy symmetric comparison)
                self.progress.emit(f"Comparing surfaces: {self.obj1} vs {self.obj2}...")
                
                result = analyzer.compare_surfaces(
                    self.obj1, self.obj2,
                    self.selection1, self.selection2
                )
                
                self.progress.emit(f"Similarity analysis complete. Score: {result.score:.3f}")
                
                # Export results if CSV requested
                if self.out_csv:
                    import csv
                    with open(self.out_csv, 'w', newline='', encoding='utf-8') as f:
                        writer = csv.writer(f)
                        writer.writerow(['Metric', 'Value'])
                        writer.writerow(['Overall_Similarity', f"{result.score:.4f}"])
                        writer.writerow(['Geometric_Similarity', f"{result.geometric_similarity:.4f}"])
                        writer.writerow(['Chemical_Similarity', f"{result.chemical_similarity:.4f}"])
                        writer.writerow(['Shape_Index_Correlation', f"{result.shape_index_correlation:.4f}"])
                        writer.writerow(['Curvature_Correlation', f"{result.curvature_correlation:.4f}"])
                        writer.writerow(['Electrostatic_Correlation', f"{result.electrostatic_correlation:.4f}"])
                        writer.writerow(['Hydrophobicity_Correlation', f"{result.hydrophobicity_correlation:.4f}"])
                    self.progress.emit(f"Results saved to {self.out_csv}")
                    
            else:
                # Single surface analysis
                self.progress.emit(f"Analyzing surface: {self.obj1}...")
                
                mesh, points = analyzer.analyze_surface(
                    obj_name=self.obj1,
                    selection=self.selection1,
                    use_apbs=self.use_apbs
                )
                
                self.progress.emit(f"Surface analysis complete. {mesh.n_vertices} vertices, {mesh.n_faces} faces")
                
                # Export features if CSV requested
                if self.out_csv:
                    analyzer.export_features(self.obj1, self.out_csv, self.selection1)
                    self.progress.emit(f"Features exported to {self.out_csv}")
                
                # Create a simple result object for single surface
                result = {
                    'mesh': mesh,
                    'points': points,
                    'n_vertices': mesh.n_vertices,
                    'n_faces': mesh.n_faces,
                    'surface_area': mesh.surface_area
                }
            
            self.finished.emit(result, self.out_csv or "")
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.error.emit(str(e))
