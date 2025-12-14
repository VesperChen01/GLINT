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
from ..workers import GMotifWorker, C2H2Worker, SurfaceAnalysisWorker

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
        self._c2h2_hits = []
        self._last_c2h2_csv = None
        
        # Remove proxies for methods implemented here to avoid shadowing
        for attr in ['start_gmotif', 'render_gmotif_with_esp', 'browse_gm_pdb', 'browse_gm_out_csv',
                     'start_c2h2', 'browse_c2h2_pdb', 'browse_c2h2_out_csv']:
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
        
        # 3. C2H2 Zinc Finger Detection
        grp_c2h2 = QGroupBox("C2H2 Zinc Finger Detection")
        c2h2_grid = QGridLayout(grp_c2h2)
        c2h2_grid.setColumnStretch(1, 1); c2h2_grid.setColumnStretch(3, 1)
        c2h2_grid.setHorizontalSpacing(8); c2h2_grid.setVerticalSpacing(10)
        
        # Row 0: Target Object / PDB File
        c2h2_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.obj_combo_c2h2 = QComboBox(); self.parent_window.obj_combo_c2h2.setMinimumHeight(32)
        self.parent_window.refresh_obj_c2h2 = QPushButton(t("refresh")); self.parent_window.refresh_obj_c2h2.clicked.connect(self.refresh_objects)
        r0_c2h2 = QHBoxLayout(); r0_c2h2.addWidget(self.parent_window.obj_combo_c2h2, 1); r0_c2h2.addWidget(self.parent_window.refresh_obj_c2h2)
        c2h2_grid.addLayout(r0_c2h2, 0, 1)
        
        c2h2_grid.addWidget(QLabel("PDB File (opt):"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.c2h2_pdb = QLineEdit(); self.parent_window.c2h2_pdb_browse = QPushButton(t("browse"))
        self.parent_window.c2h2_pdb_browse.clicked.connect(self.browse_c2h2_pdb)
        r0b_c2h2 = QHBoxLayout(); r0b_c2h2.addWidget(self.parent_window.c2h2_pdb, 1); r0b_c2h2.addWidget(self.parent_window.c2h2_pdb_browse)
        c2h2_grid.addLayout(r0b_c2h2, 0, 3)
        
        # Row 1: Turn RMSD / Global RMSD
        c2h2_grid.addWidget(QLabel("Turn RMSD (\u00c5):"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.c2h2_turn_rmsd = QLineEdit("2.0")
        self.parent_window.c2h2_turn_rmsd.setToolTip("局部 turn 对齐 RMSD 阈值（更敏感）")
        c2h2_grid.addWidget(self.parent_window.c2h2_turn_rmsd, 1, 1)
        
        c2h2_grid.addWidget(QLabel("Global RMSD (\u00c5):"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.c2h2_global_rmsd = QLineEdit("3.5")
        self.parent_window.c2h2_global_rmsd.setToolTip("全局 fold check RMSD 阈值")
        c2h2_grid.addWidget(self.parent_window.c2h2_global_rmsd, 1, 3)
        
        # Row 2: Checkboxes
        self.parent_window.c2h2_require_turn_gly = QCheckBox("Require Turn Gly")
        self.parent_window.c2h2_require_turn_gly.setChecked(False)  # 默认不做 hard filter
        self.parent_window.c2h2_require_turn_gly.setToolTip("要求 turn 区域有关键 Gly（启用会降低召回）")
        c2h2_grid.addWidget(self.parent_window.c2h2_require_turn_gly, 2, 1)
        
        self.parent_window.c2h2_skip_low_complexity = QCheckBox("Skip Low-Complexity")
        self.parent_window.c2h2_skip_low_complexity.setChecked(False)  # 默认保留但降权
        self.parent_window.c2h2_skip_low_complexity.setToolTip("跳过 polyQ 等低复杂度区域（启用会降低召回）")
        c2h2_grid.addWidget(self.parent_window.c2h2_skip_low_complexity, 2, 3)
        
        # Row 3: Output CSV
        c2h2_grid.addWidget(QLabel("Output CSV:"), 3, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.c2h2_out_csv = QLineEdit()
        self.parent_window.c2h2_out_browse = QPushButton(t("browse")); self.parent_window.c2h2_out_browse.clicked.connect(self.browse_c2h2_out_csv)
        r3_c2h2 = QHBoxLayout(); r3_c2h2.addWidget(self.parent_window.c2h2_out_csv, 1); r3_c2h2.addWidget(self.parent_window.c2h2_out_browse)
        c2h2_grid.addLayout(r3_c2h2, 3, 1, 1, 3)
        
        # C2H2 Buttons
        c2h2_btn_row = QHBoxLayout()
        self.parent_window.c2h2_btn = QPushButton("Find Zinc Fingers"); self.parent_window.c2h2_btn.setObjectName("highlight_btn")
        self.parent_window.c2h2_btn.clicked.connect(self.start_c2h2)
        self.parent_window.c2h2_btn_render = QPushButton("Render C2H2 + ESP"); self.parent_window.c2h2_btn_render.setObjectName("highlight_btn")
        self.parent_window.c2h2_btn_render.clicked.connect(self.render_c2h2_with_esp)
        c2h2_btn_row.addWidget(self.parent_window.c2h2_btn); c2h2_btn_row.addWidget(self.parent_window.c2h2_btn_render); c2h2_btn_row.addStretch(1)
        
        layout.addWidget(grp_c2h2)
        layout.addLayout(c2h2_btn_row)
        
        layout.addWidget(grp_c2h2)
        layout.addLayout(c2h2_btn_row)
        
        # 4. Surface Analysis
        grp_surf = QGroupBox("Protein Surface Analysis")
        surf_grid = QGridLayout(grp_surf)
        surf_grid.setColumnStretch(1, 1); surf_grid.setColumnStretch(3, 1)
        surf_grid.setHorizontalSpacing(8); surf_grid.setVerticalSpacing(10)
        
        # Row 0: Target Object
        surf_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.obj_combo_surf = QComboBox(); self.parent_window.obj_combo_surf.setMinimumHeight(32)
        self.parent_window.refresh_obj_surf = QPushButton(t("refresh")); self.parent_window.refresh_obj_surf.clicked.connect(self.refresh_objects)
        r0_surf = QHBoxLayout(); r0_surf.addWidget(self.parent_window.obj_combo_surf, 1); r0_surf.addWidget(self.parent_window.refresh_obj_surf)
        surf_grid.addLayout(r0_surf, 0, 1)

        # Row 0: Output CSV
        surf_grid.addWidget(QLabel("Output CSV:"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.surf_out_csv = QLineEdit()
        self.parent_window.surf_out_browse = QPushButton(t("browse"))
        self.parent_window.surf_out_browse.clicked.connect(self.browse_surf_out_csv)
        r0b_surf = QHBoxLayout(); r0b_surf.addWidget(self.parent_window.surf_out_csv, 1); r0b_surf.addWidget(self.parent_window.surf_out_browse)
        surf_grid.addLayout(r0b_surf, 0, 3)
        
        # Row 1: Description
        desc_label = QLabel("Detects electrostatic and hydrophobic patches on the surface.")
        desc_label.setStyleSheet("color: gray; font-style: italic;")
        surf_grid.addWidget(desc_label, 1, 1, 1, 3)
        
        # Buttons
        surf_btn_row = QHBoxLayout()
        self.parent_window.surf_btn = QPushButton("Analyze Surface")
        self.parent_window.surf_btn.setObjectName("highlight_btn")
        self.parent_window.surf_btn.clicked.connect(self.start_surface_analysis)
        
        self.parent_window.surf_vis_btn = QPushButton("Visualize Patches")
        self.parent_window.surf_vis_btn.setObjectName("highlight_btn")
        self.parent_window.surf_vis_btn.clicked.connect(self.render_surface_patches)
        
        surf_btn_row.addWidget(self.parent_window.surf_btn)
        surf_btn_row.addWidget(self.parent_window.surf_vis_btn)
        surf_btn_row.addStretch(1)
        
        layout.addWidget(grp_surf)
        layout.addLayout(surf_btn_row)

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

    # --- C2H2 Zinc Finger Logic ---
    def browse_c2h2_pdb(self):
        fn, _ = QFileDialog.getOpenFileName(self, t("select_pdb"), "", "PDB (*.pdb *.cif);;All Files (*)")
        if fn: self.parent_window.c2h2_pdb.setText(fn); self.parent_window.update_enablement()

    def browse_c2h2_out_csv(self):
        fn, _ = QFileDialog.getSaveFileName(self, t("select_outcsv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.parent_window.c2h2_out_csv.setText(fn); self.parent_window.update_enablement()

    def start_c2h2(self):
        """Start C2H2 zinc finger detection"""
        obj = self.parent_window.obj_combo_c2h2.currentText().strip()
        pdb = self.parent_window.c2h2_pdb.text().strip() or None
        outcsv = self.parent_window.c2h2_out_csv.text().strip() or None
        if not obj or obj == t("no_object"):
            QMessageBox.warning(self, t("title"), t("no_object")); return
        
        try:
            turn_rmsd = float(self.parent_window.c2h2_turn_rmsd.text().strip() or "2.0")
        except Exception:
            turn_rmsd = 2.0
        try:
            global_rmsd = float(self.parent_window.c2h2_global_rmsd.text().strip() or "3.5")
        except Exception:
            global_rmsd = 3.5
        
        require_turn_gly = self.parent_window.c2h2_require_turn_gly.isChecked()
        skip_low_complexity = self.parent_window.c2h2_skip_low_complexity.isChecked()
        
        self.parent_window.c2h2_btn.setEnabled(False)
        self.parent_window.progress_bar.setVisible(True); self.parent_window.progress_bar.setRange(0, 0)
        self.parent_window.c2h2_thread = C2H2Worker(
            obj, pdb, turn_rmsd, global_rmsd, require_turn_gly, skip_low_complexity, outcsv
        )
        self.parent_window.c2h2_thread.progress.connect(self.log)
        self.parent_window.c2h2_thread.error.connect(self.on_error)
        self.parent_window.c2h2_thread.finished.connect(self.on_finished_c2h2)
        self.parent_window.c2h2_thread.start()

    def on_finished_c2h2(self, hits: List[Tuple], out_csv_path: str):
        """Handle C2H2 detection completion"""
        self._c2h2_hits = hits or []
        self._last_c2h2_csv = out_csv_path
        
        # Count by status
        n_pass = sum(1 for h in hits if h[6] == "pass")
        n_candidate = sum(1 for h in hits if h[6] == "candidate")
        
        self.log(f"✅ C2H2 detection complete: {len(hits)} domains found ({n_pass} pass, {n_candidate} candidate)")
        if os.path.exists(out_csv_path):
            self.parent_window.c2h2_out_csv.setText(out_csv_path)
            self.log(f"   Saved to: {os.path.basename(out_csv_path)}")
        self.parent_window.progress_bar.setVisible(False); self.parent_window.progress_bar.setRange(0, 1)
        self.parent_window.c2h2_btn.setEnabled(True)

    def render_c2h2_with_esp(self):
        """Render C2H2 domains with electrostatic surface"""
        try:
            obj = self.parent_window.obj_combo_c2h2.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object")); return

            from pymol import cmd
            
            # Generate ESP
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

            # Show surface around C2H2 domains
            if self._last_c2h2_csv and os.path.exists(self._last_c2h2_csv):
                try:
                    import csv
                    c2h2_regions = []
                    with open(self._last_c2h2_csv, "r", encoding="utf-8") as f:
                        reader = csv.DictReader(f)
                        for r in reader:
                            chain = r.get("Chain", "")
                            start = r.get("Domain_Start", "")
                            end = r.get("Domain_End", "")
                            if chain and start and end:
                                c2h2_regions.append((chain, start, end))

                    if c2h2_regions:
                        c2h2_selections = [f"(chain {ch} and resi {st}-{ed})" for ch, st, ed in c2h2_regions]
                        all_c2h2_sel = " or ".join(c2h2_selections)
                        surface_sel_name = "c2h2_surface_area"
                        cmd.select(surface_sel_name, f"byres ({obj} within 8 of ({all_c2h2_sel}))")
                        cmd.show("surface", surface_sel_name)
                        cmd.set("surface_quality", 1, surface_sel_name)
                        cmd.set("transparency", 0.2, surface_sel_name)
                        cmd.color(ramp_name, surface_sel_name)
                        self.log(f"Electrostatic surface shown around {len(c2h2_regions)} C2H2 domains (8 \u00c5 region)")
                    else:
                        self._show_full_surface_esp(obj, ramp_name)
                except Exception as e:
                    self.log(f"ESP around C2H2 failed: {e}")
                    self._show_full_surface_esp(obj, ramp_name)
            else:
                self._show_full_surface_esp(obj, ramp_name)

            cmd.set("ambient", 0.2)
            cmd.set("spec_power", 80)
            cmd.orient(obj)
            
            self.log(f"Rendered: C2H2 + ESP (grid={grid} \u00c5)")
            
        except Exception as e:
            self.on_error(str(e))

    # --- Surface Analysis Logic ---
    def browse_surf_out_csv(self):
        fn, _ = QFileDialog.getSaveFileName(self, t("select_outcsv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.parent_window.surf_out_csv.setText(fn)

    def start_surface_analysis(self):
        """Start Surface Analysis"""
        obj = self.parent_window.obj_combo_surf.currentText().strip()
        outcsv = self.parent_window.surf_out_csv.text().strip() or None
        
        if not obj or obj == t("no_object"):
            QMessageBox.warning(self, t("title"), t("no_object")); return
            
        self.parent_window.surf_btn.setEnabled(False)
        self.parent_window.progress_bar.setVisible(True); self.parent_window.progress_bar.setRange(0, 0)
        
        self.parent_window.surf_thread = SurfaceAnalysisWorker(obj, outcsv)
        self.parent_window.surf_thread.progress.connect(self.log)
        self.parent_window.surf_thread.error.connect(self.on_error)
        self.parent_window.surf_thread.finished.connect(self.on_finished_surface)
        self.parent_window.surf_thread.start()

    def on_finished_surface(self, patches: List, out_csv_path: str):
        self._last_surf_csv = out_csv_path
        self._last_surf_patches = patches
        
        self.log(f"✅ Surface analysis complete: {len(patches)} patches found")
        if out_csv_path and os.path.exists(out_csv_path):
             self.parent_window.surf_out_csv.setText(out_csv_path)
             
        self.parent_window.progress_bar.setVisible(False); self.parent_window.progress_bar.setRange(0, 1)
        self.parent_window.surf_btn.setEnabled(True)

    def render_surface_patches(self):
        """Visualize detected patches"""
        try:
            from pymol import cmd
            obj = self.parent_window.obj_combo_surf.currentText().strip()
            
            if not getattr(self, '_last_surf_patches', None):
                QMessageBox.warning(self, "Warning", "No analysis results to visualize. Please run analysis first.")
                return

            cmd.hide("everything", obj)
            cmd.show("surface", obj)
            cmd.color("white", obj)
            cmd.set("transparency", 0.3, obj)
            
            # Create selections for each patch and color them
            for p in self._last_surf_patches:
                # Use residues to select patch area
                # p.residues is list of dicts {chain, resn, resi}
                if not p.residues: continue
                
                sel_str = " or ".join([f"(chain {r['chain']} and resi {r['resi']})" for r in p.residues])
                patch_name = f"patch_{p.id}_{p.type}"
                
                # Expand selection slightly to cover surface
                cmd.select(patch_name, f"byres ({obj} and ({sel_str}))")
                
                # Color based on type
                color = "blue" # neg
                if p.type == "electrostatic_pos": color = "red"
                elif p.type == "hydrophobic": color = "green"
                elif p.type == "electrostatic_neg": color = "blue"
                
                cmd.color(color, patch_name)
                
            self.log(f"Visualized {len(self._last_surf_patches)} patches on {obj}")
            
        except Exception as e:
            self.on_error(str(e))

