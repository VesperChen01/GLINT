# -*- coding: utf-8 -*-
"""
Worker threads for GlueTK GUI.
"""
import os
from typing import Optional, List, Tuple, Dict, Any

try:
    from PyQt5.QtCore import QThread, pyqtSignal
except ImportError:
    try:
        from PyQt6.QtCore import QThread, pyqtSignal
    except ImportError:
        raise RuntimeError("PyQt5 or PyQt6 must be installed.")

from .utils import (
    get_lang, t,
    analyze_pdb_interactions,
    analyze_protein_nucleic_interactions,
    find_crbn_g_motif
)

class AnalysisWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)
    def __init__(self, obj_name: str, pdb_file: Optional[str], output_csv: Optional[str]):
        super().__init__()
        self.obj_name = obj_name
        self.pdb_file = pdb_file
        self.output_csv = output_csv
    def run(self):
        try:
            self.progress.emit(t("log_start"))
            interactions = analyze_pdb_interactions(
                obj_name=self.obj_name,
                pdb_file=self.pdb_file,
                only_between_chains=True,  # 默认启用链间分析
                output_csv=self.output_csv,
                auto_highlight=True,  # 启用自动高亮
            )
            self.progress.emit(t("log_done").format(n=len(interactions)))
            self.finished.emit(interactions)
        except Exception as e:
            self.error.emit(str(e))

class PNAnalysisWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)
    def __init__(self, obj_name, nucleic_chains, protein_chains, output_csv, distance_cutoff, pdb_file=None):
        super().__init__()
        self.obj_name = obj_name
        self.nucleic_chains = nucleic_chains
        self.protein_chains = protein_chains
        self.output_csv = output_csv
        self.distance_cutoff = distance_cutoff
        self.pdb_file = pdb_file
        
    def run(self):
        try:
            self.progress.emit("Starting Protein-Nucleic Acid Analysis...")
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
            self.error.emit(str(e))

class GMotifWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, str)
    error = pyqtSignal(str)
    def __init__(self, obj_name: str, pdb_file: Optional[str], rmsd: float, require_gly: bool,
                 out_csv: Optional[str], template_mode: str, template_sel: Optional[str], template_builtin: Optional[str],
                 exclude_proline: bool = True, check_surface_exposure: bool = True, min_sasa: float = 15.0):
        super().__init__()
        self.obj_name = obj_name; self.pdb_file = pdb_file
        self.rmsd = rmsd; self.require_gly = require_gly; self.out_csv = out_csv
        self.template_mode = template_mode; self.template_sel = template_sel; self.template_builtin = template_builtin
        self.exclude_proline = exclude_proline
        self.check_surface_exposure = check_surface_exposure
        self.min_sasa = min_sasa
    def run(self):
        try:
            if find_crbn_g_motif is None:
                raise RuntimeError("find_crbn_g_motif not found; ensure g_motif_analyzer.py exists in the plugin directory.")
            self.progress.emit("[G-Motif] " + ("开始识别…" if get_lang()=="zh" else "Detecting…"))
            out_csv_path = self.out_csv
            if not out_csv_path:
                import tempfile, os
                fd, out_csv_path = tempfile.mkstemp(suffix="_gmotif.csv"); os.close(fd)
            hits = find_crbn_g_motif(
                obj_name=self.obj_name,
                pdb_file=self.pdb_file,
                template_mode=self.template_mode,
                template_sel=self.template_sel,
                template_builtin=self.template_builtin,
                rmsd_cutoff=float(self.rmsd),
                out_csv=out_csv_path,
                auto_highlight=1,  # 启用自动高亮
                require_gly_pos6=bool(self.require_gly),
                exclude_proline=bool(self.exclude_proline),
                check_surface_exposure=bool(self.check_surface_exposure),
                min_sasa_per_residue=float(self.min_sasa),
            ) or []
            self.progress.emit("[G-Motif] " + (f"完成，命中 {len(hits)} 条" if get_lang()=="zh" else f"Done, {len(hits)} hits"))
            self.finished.emit(hits, out_csv_path)
        except Exception as e:
            self.error.emit(str(e))
