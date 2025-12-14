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
    find_crbn_g_motif,
    find_c2h2_domains
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


class C2H2Worker(QThread):
    """C2H2 锌指蛋白检测 Worker"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, str)  # (domains, csv_path)
    error = pyqtSignal(str)
    
    def __init__(self, obj_name: str, pdb_file: Optional[str], 
                 turn_rmsd: float, global_rmsd: float,
                 require_turn_gly: bool, skip_low_complexity: bool,
                 out_csv: Optional[str], hmm_profile: Optional[str] = None):
        super().__init__()
        self.obj_name = obj_name
        self.pdb_file = pdb_file
        self.turn_rmsd = turn_rmsd
        self.global_rmsd = global_rmsd
        self.require_turn_gly = require_turn_gly
        self.skip_low_complexity = skip_low_complexity
        self.out_csv = out_csv
        self.hmm_profile = hmm_profile
        
    def run(self):
        try:
            if find_c2h2_domains is None:
                raise RuntimeError("find_c2h2_domains not found; ensure c2h2_finder.py exists.")
            
            self.progress.emit("[C2H2] " + ("开始检测锌指域…" if get_lang()=="zh" else "Detecting zinc fingers…"))
            
            out_csv_path = self.out_csv
            if not out_csv_path:
                import tempfile
                fd, out_csv_path = tempfile.mkstemp(suffix="_c2h2.csv")
                os.close(fd)
            
            domains = find_c2h2_domains(
                obj_name=self.obj_name,
                pdb_file=self.pdb_file,
                hmm_profile=self.hmm_profile,
                turn_rmsd_cutoff=float(self.turn_rmsd),
                global_rmsd_cutoff=float(self.global_rmsd),
                require_turn_gly=bool(self.require_turn_gly),
                skip_low_complexity=bool(self.skip_low_complexity),
                out_csv=out_csv_path,
                auto_highlight=True,
            ) or []
            
            # 转换为简单 tuple 列表以便 GUI 处理
            hits = [(d.chain, d.sequence, d.domain_start, d.domain_end, 
                     d.turn_rmsd if d.turn_rmsd else 0.0, d.priority_score, d.status) 
                    for d in domains]
            
            self.progress.emit("[C2H2] " + (f"完成，发现 {len(hits)} 个锌指域" if get_lang()=="zh" else f"Done, {len(hits)} zinc fingers found"))
            self.finished.emit(hits, out_csv_path)
        except Exception as e:
            self.error.emit(str(e))

class SurfaceAnalysisWorker(QThread):
    """Worker for protein surface analysis"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, str) # (patches, out_csv)
    error = pyqtSignal(str)

    def __init__(self, obj_name: str, out_csv: str):
        super().__init__()
        self.obj_name = obj_name
        self.out_csv = out_csv
        
    def run(self):
        try:
            self.progress.emit(f"Starting surface analysis for {self.obj_name}...")
            
            # Late import to avoid circular dependencies if any
            from ..protein_surface_analyzer import SurfaceAnalyzer
            
            analyzer = SurfaceAnalyzer(self.obj_name)
            patches = analyzer.analyze()
            
            self.progress.emit(f"Analysis complete. Found {len(patches)} patches.")
            
            # Export to CSV if requested
            if self.out_csv:
                import csv
                with open(self.out_csv, 'w', newline='', encoding='utf-8') as f:
                    fieldnames = ['id', 'type', 'area', 'score', 'center', 'residues']
                    writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
                    writer.writeheader()
                    for p in patches:
                        writer.writerow(p.to_dict())
                self.progress.emit(f"Report saved to {self.out_csv}")
                
            self.finished.emit(patches, self.out_csv)
            
        except Exception as e:
            self.error.emit(str(e))
