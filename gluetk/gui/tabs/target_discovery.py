# -*- coding: utf-8 -*-
"""
Target Discovery Tab: G-Motif, Disease, Pocket
"""
import os
from typing import Optional, List, Tuple, Dict, Any

try:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
        QCheckBox, QComboBox, QGroupBox, QGridLayout, QScrollArea, QFrame,
        QFileDialog, QMessageBox
    )
except ImportError:
    try:
        from PyQt6.QtCore import Qt
        from PyQt6.QtWidgets import (
            QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
            QCheckBox, QComboBox, QGroupBox, QGridLayout, QScrollArea, QFrame,
            QFileDialog, QMessageBox
        )
    except ImportError:
        raise RuntimeError("PyQt5 or PyQt6 must be installed.")

from ..utils import t
from .common import CommonTab
from ..workers import GMotifWorker

try:
    from ...open_targets_api import search_disease
except (ImportError, ValueError):
    try:
        from ...open_targets_api import search_disease
    except Exception:
        search_disease = None

class TargetDiscoveryTab(CommonTab):
    def __init__(self, parent):
        super().__init__(parent)
        self._gmotif_hits = []
        self._last_gmotif_csv = None
        
        # Remove proxies for methods implemented here to avoid shadowing
        for attr in ['start_gmotif', 'render_gmotif_with_esp', 'browse_gm_pdb', 'browse_gm_out_csv']:
            if attr in self.__dict__:
                del self.__dict__[attr]
        
        self.init_ui()
        
    def init_ui(self):
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        
        content_widget = QWidget()
        content_widget.setObjectName("scroll_content")
        self.parent_window._target_scroll_content = content_widget
        bg_color = "#0d1117" if getattr(self.parent_window, "_dark_mode", True) else "#f8fafc"
        content_widget.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")
        
        layout = QVBoxLayout(content_widget)
        layout.setSpacing(14)
        layout.setContentsMargins(12, 12, 12, 12)
        
        # 1. Disease Analysis
        if search_disease is not None:
            try:
                # Try to import DiseaseAnalysisTab
                from ...disease_analysis_gui import DiseaseAnalysisTab
                disease_card = DiseaseAnalysisTab(self.parent_window)
                grp_disease = QGroupBox("Disease Target Analysis")
                d_layout = QVBoxLayout(grp_disease)
                d_layout.addWidget(disease_card)
                layout.addWidget(grp_disease)
            except Exception as e:
                self.log(f"Failed to load Disease Analysis: {e}")
        
        # 2. G-Motif Detection
        grp_gm = QGroupBox("G-Motif (CRBN G-loop) Detection")
        gm_grid = QGridLayout(grp_gm)
        gm_grid.setColumnStretch(1, 1); gm_grid.setColumnStretch(3, 1)
        gm_grid.setHorizontalSpacing(8); gm_grid.setVerticalSpacing(10)
        
        # Row 0
        gm_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.obj_combo_gm = QComboBox(); self.parent_window.obj_combo_gm.setMinimumHeight(32)
        self.parent_window.refresh_obj_gm = QPushButton(t("refresh")); self.parent_window.refresh_obj_gm.clicked.connect(self.refresh_objects)
        r0 = QHBoxLayout(); r0.addWidget(self.parent_window.obj_combo_gm, 1); r0.addWidget(self.parent_window.refresh_obj_gm)
        gm_grid.addLayout(r0, 0, 1)
        
        gm_grid.addWidget(QLabel("PDB File (opt):"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.gm_pdb = QLineEdit(); self.parent_window.gm_pdb_browse = QPushButton(t("browse"))
        self.parent_window.gm_pdb_browse.clicked.connect(self.browse_gm_pdb)
        r0b = QHBoxLayout(); r0b.addWidget(self.parent_window.gm_pdb, 1); r0b.addWidget(self.parent_window.gm_pdb_browse)
        gm_grid.addLayout(r0b, 0, 3)
        
        # Row 1
        gm_grid.addWidget(QLabel("RMSD cutoff (Å):"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.gm_rmsd = QLineEdit("3.5")
        gm_grid.addWidget(self.parent_window.gm_rmsd, 1, 1)
        
        self.parent_window.gm_require_gly = QCheckBox(t("require_gly")); self.parent_window.gm_require_gly.setChecked(True)
        gm_grid.addWidget(self.parent_window.gm_require_gly, 1, 3)
        
        # Row 2
        gm_grid.addWidget(QLabel("Template:"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.gm_template_mode = QComboBox()
        self.parent_window.gm_template_mode.addItems(["Idealized (8×Cα)", "Built-in: GSPT1", "Built-in: CK1α", "Built-in: VAV1", "From Selection"])
        gm_grid.addWidget(self.parent_window.gm_template_mode, 2, 1)
        
        gm_grid.addWidget(QLabel("Selection:"), 2, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.gm_template_sel = QLineEdit()
        self.parent_window.gm_template_pick = QPushButton("Pick (sele)"); self.parent_window.gm_template_pick.clicked.connect(lambda: self.parent_window.gm_template_sel.setText("sele"))
        r2b = QHBoxLayout(); r2b.addWidget(self.parent_window.gm_template_sel, 1); r2b.addWidget(self.parent_window.gm_template_pick)
        gm_grid.addLayout(r2b, 2, 3)
        
        def _toggle_template_inputs(idx):
            use_sel = (idx == 4)
            self.parent_window.gm_template_sel.setEnabled(use_sel); self.parent_window.gm_template_pick.setEnabled(use_sel)
        self.parent_window.gm_template_mode.currentIndexChanged.connect(_toggle_template_inputs)
        _toggle_template_inputs(0)
        
        # Row 3
        gm_grid.addWidget(QLabel("Output CSV:"), 3, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.gm_out_csv = QLineEdit()
        self.parent_window.gm_out_browse = QPushButton(t("browse")); self.parent_window.gm_out_browse.clicked.connect(self.browse_gm_out_csv)
        r3 = QHBoxLayout(); r3.addWidget(self.parent_window.gm_out_csv, 1); r3.addWidget(self.parent_window.gm_out_browse)
        gm_grid.addLayout(r3, 3, 1, 1, 3)
        
        # Buttons
        gm_btn_row = QHBoxLayout()
        self.parent_window.gm_btn = QPushButton("Detect POI"); self.parent_window.gm_btn.setObjectName("highlight_btn")
        self.parent_window.gm_btn.clicked.connect(self.start_gmotif)
        self.parent_window.gm_btn_render = QPushButton("Render All (POI + ESP + PNG)"); self.parent_window.gm_btn_render.setObjectName("highlight_btn")
        self.parent_window.gm_btn_render.clicked.connect(self.render_gmotif_with_esp)
        gm_btn_row.addWidget(self.parent_window.gm_btn); gm_btn_row.addWidget(self.parent_window.gm_btn_render); gm_btn_row.addStretch(1)
        
        layout.addWidget(grp_gm)
        layout.addLayout(gm_btn_row)
        
        layout.addStretch(1)
        scroll_area.setWidget(content_widget)
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(scroll_area)

    # --- G-Motif Logic ---
    def browse_gm_pdb(self):
        fn, _ = QFileDialog.getOpenFileName(self, t("select_pdb"), "", "PDB (*.pdb *.cif);;All Files (*)")
        if fn: self.parent_window.gm_pdb.setText(fn); self.parent_window.update_enablement()

    def browse_gm_out_csv(self):
        fn, _ = QFileDialog.getSaveFileName(self, t("select_outcsv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.parent_window.gm_out_csv.setText(fn); self.parent_window.update_enablement()

    def start_gmotif(self):
        obj = self.parent_window.obj_combo_gm.currentText().strip()
        pdb = self.parent_window.gm_pdb.text().strip() or None
        outcsv = self.parent_window.gm_out_csv.text().strip() or None
        if not obj or obj == t("no_object"):
            QMessageBox.warning(self, t("title"), t("no_object")); return
        try:
            rmsd = float(self.parent_window.gm_rmsd.text().strip() or "3.5")
        except Exception:
            rmsd = 3.5
        require_gly = self.parent_window.gm_require_gly.isChecked()

        idx = self.parent_window.gm_template_mode.currentIndex()
        if idx == 0:
            template_mode, template_sel, template_builtin = "ideal", None, None
        elif idx in (1, 2, 3):
            template_mode = "builtin"
            template_sel = None
            template_builtin = [
                "GSPT1 (6H0G A:60-67)",
                "CK1α (3M51 A:36-43)",
                "VAV1 (2MC1 A:95-102)",
            ][idx - 1]
        else:
            template_mode, template_sel, template_builtin = "selection", (self.parent_window.gm_template_sel.text().strip() or None), None

        self.parent_window.gm_btn.setEnabled(False)
        self.parent_window.progress_bar.setVisible(True); self.parent_window.progress_bar.setRange(0, 0)
        self.parent_window.gmotif_thread = GMotifWorker(obj, pdb, rmsd, require_gly, outcsv,
                                        template_mode, template_sel, template_builtin)
        self.parent_window.gmotif_thread.progress.connect(self.log)
        self.parent_window.gmotif_thread.error.connect(self.on_error)
        self.parent_window.gmotif_thread.finished.connect(self.on_finished_gmotif)
        self.parent_window.gmotif_thread.start()

    def on_finished_gmotif(self, hits: List[Tuple], out_csv_path: str):
        self._gmotif_hits = hits or []
        self._last_gmotif_csv = out_csv_path
        self.log(f"✅ G-Motif detection complete: {len(hits)} hits found")
        if os.path.exists(out_csv_path):
            self.parent_window.gm_out_csv.setText(out_csv_path)
            self.log(f"   Saved to: {os.path.basename(out_csv_path)}")
        self.parent_window.progress_bar.setVisible(False); self.parent_window.progress_bar.setRange(0, 1)
        self.parent_window.gm_btn.setEnabled(True)

    def render_gmotif_with_esp(self):
        # Similar to VisualizationTab methods but specific to G-Motif results
        try:
            obj = self.parent_window.obj_combo_gm.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object")); return

            from pymol import cmd
            
            try:
                from ...highlight_residues import highlight_gmotif_loops
            except ImportError:
                try: from highlight_residues import highlight_gmotif_loops
                except ImportError: highlight_gmotif_loops = None

            if self._last_gmotif_csv and os.path.exists(self._last_gmotif_csv) and highlight_gmotif_loops:
                try:
                    highlight_gmotif_loops(self._last_gmotif_csv, obj, color="yellow", show_labels=True)
                    self.log("Applied G-Motif highlight (latest CSV)")
                except Exception as e:
                    self.log(f"G-Motif highlight skipped: {e}")

            # Generate ESP (using APBS tab settings if available, else default)
            grid = 1.0
            vmin, v0, vmax = -5.0, 0.0, 5.0
            if hasattr(self.parent_window, "apbs_grid"):
                 try: grid = float(self.parent_window.apbs_grid.text().strip() or "1.0")
                 except: pass
            if hasattr(self.parent_window, "apbs_range"):
                 try: vmin, v0, vmax = [float(x) for x in self.parent_window.apbs_range.text().split(",")]
                 except: pass

            map_name = f"{obj}_esp_map"
            ramp_name = f"{obj}_esp_ramp"
            cmd.map_new(map_name, "coulomb", grid, obj)
            cmd.ramp_new(ramp_name, map_name, [vmin, v0, vmax], ["blue", "white", "red"])
            
            cmd.hide("everything", obj)
            cmd.show("cartoon", obj)
            cmd.set("cartoon_transparency", 0.3, obj)

            # 4) Show surface only around G-loop
            if self._last_gmotif_csv and os.path.exists(self._last_gmotif_csv):
                try:
                    import csv
                    gloop_regions = []
                    with open(self._last_gmotif_csv, "r", encoding="utf-8") as f:
                        reader = csv.DictReader(f)
                        for r in reader:
                            chain = r.get("Chain", "")
                            start = r.get("Start", "")
                            end = r.get("End", "")
                            if chain and start and end:
                                gloop_regions.append((chain, start, end))

                    if gloop_regions:
                        gloop_selections = [f"(chain {ch} and resi {st}-{ed})" for ch, st, ed in gloop_regions]
                        all_gloop_sel = " or ".join(gloop_selections)
                        surface_sel_name = "gloop_surface_area"
                        cmd.select(surface_sel_name, f"byres ({obj} within 10 of ({all_gloop_sel}))")
                        cmd.show("surface", surface_sel_name)
                        cmd.set("surface_quality", 1, surface_sel_name)
                        cmd.set("surface_color_smoothing", 1, surface_sel_name)
                        cmd.set("transparency", 0.2, surface_sel_name)
                        cmd.color(ramp_name, surface_sel_name)

                        for ch, st, ed in gloop_regions:
                            gloop_sel = f"{obj} and chain {ch} and resi {st}-{ed}"
                            cmd.hide("surface", gloop_sel)
                        self.log("Electrostatic surface shown only around G-loop (10 Å region)")
                    else:
                         self._show_full_surface_esp(obj, ramp_name)
                except Exception as e:
                    self.log(f"ESP around G-loop failed: {e}")
                    self._show_full_surface_esp(obj, ramp_name)
            else:
                self._show_full_surface_esp(obj, ramp_name)

            cmd.set("ambient", 0.2)
            cmd.set("spec_power", 80)
            cmd.set("spec_reflect", 0.3)
            cmd.set("depth_cue", 1)
            cmd.set("fog_start", 0.45)
            cmd.orient(obj)
            
            self.log(f"Rendered: G-Motif + ESP (grid={grid} Å)")
            
        except Exception as e:
            self.on_error(str(e))

    def _show_full_surface_esp(self, obj: str, ramp_name: str):
        from pymol import cmd
        cmd.show("surface", obj)
        cmd.set("surface_quality", 1, obj)
        cmd.set("surface_color_smoothing", 1, obj)
        cmd.set("transparency", 0.2, obj)
        cmd.color(ramp_name, obj)
        self.log("Showing full-protein electrostatic surface")


