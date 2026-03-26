# -*- coding: utf-8 -*-
"""
Target Discovery Tab: G-Motif, Disease, Pocket
"""
import os
from typing import Optional, List, Tuple, Dict, Any

from ..qt_adapter import (
    Qt, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QCheckBox, QComboBox, QGroupBox, QGridLayout, QScrollArea, QFrame,
    QFileDialog, QMessageBox
)

from ..utils import t, show_message_box
from .common import CommonTab
from ..workers import GMotifWorker, SurfaceAnalysisWorker, SurfaceSimilarityWorker

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
        self._last_similarity_result = None
        self._last_complementarity_result = None
        
        # Remove proxies for methods implemented here to avoid shadowing
        for attr in ['start_gmotif', 'browse_gm_out_csv']:
            if attr in self.__dict__:
                del self.__dict__[attr]
        
        self.init_ui()
        
    def init_ui(self):
        """InitializeUI - 现代卡片式布局"""
        self.setObjectName("scroll_content")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.parent_window._target_scroll_content = self

        is_dark = getattr(self.parent_window, "_dark_mode", False)
        bg_color = "#161b22" if is_dark else "#f8fafc"
        self.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")

        layout = QVBoxLayout(self)
        layout.setSpacing(16)
        layout.setContentsMargins(20, 20, 20, 20)

        # === 页面标题 ===
        header = QHBoxLayout()
        title = QLabel("Target Discovery")
        title.setStyleSheet("""
            font-size: 20px; font-weight: 600;
            color: #3b82f6; padding: 4px 0;
        """)
        header.addWidget(title)
        header.addStretch(1)
        layout.addLayout(header)

        # === G-Motif Detection 卡片 ===
        grp_gm = QFrame()
        grp_gm.setStyleSheet(self._get_card_style(is_dark))
        gm_layout = QVBoxLayout(grp_gm)
        gm_layout.setSpacing(12)
        gm_layout.setContentsMargins(16, 14, 16, 14)

        gm_title = QLabel("G-Motif (CRBN G-loop) Detection")
        gm_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        gm_layout.addWidget(gm_title)

        gm_grid = QGridLayout()
        gm_grid.setContentsMargins(12, 8, 12, 8)
        gm_grid.setColumnStretch(1, 1); gm_grid.setColumnStretch(3, 1)
        gm_grid.setHorizontalSpacing(12); gm_grid.setVerticalSpacing(8)

        # Row 0: Target Object | Template
        gm_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.obj_combo_gm = QComboBox(); self.parent_window.obj_combo_gm.setMinimumHeight(32)
        self.parent_window.refresh_obj_gm = QPushButton(t("refresh")); self.parent_window.refresh_obj_gm.setMinimumHeight(32); self.parent_window.refresh_obj_gm.clicked.connect(self.refresh_objects)
        r0 = QHBoxLayout(); r0.addWidget(self.parent_window.obj_combo_gm, 1); r0.addWidget(self.parent_window.refresh_obj_gm)
        gm_grid.addLayout(r0, 0, 1)

        gm_grid.addWidget(QLabel("Template:"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.gm_template_mode = QComboBox(); self.parent_window.gm_template_mode.setMinimumHeight(32)
        self.parent_window.gm_template_mode.addItems(["GSPT1 (5HXB)", "CK1α (5FQD)", "From Selection"])
        gm_grid.addWidget(self.parent_window.gm_template_mode, 0, 3)

        # Row 1: RMSD cutoff | Require Gly
        gm_grid.addWidget(QLabel("RMSD cutoff (Å):"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.gm_rmsd = QLineEdit("3.5"); self.parent_window.gm_rmsd.setMinimumHeight(32)
        gm_grid.addWidget(self.parent_window.gm_rmsd, 1, 1)

        gm_grid.addWidget(QLabel("Require Gly:"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.gm_require_gly_pos = QComboBox()
        self.parent_window.gm_require_gly_pos.addItems([
            "Pos 6 or 3 (default)",  # "6,3"
            "Pos 6 only",             # "6"
            "Pos 3 only",             # "3"
            "No requirement"          # None
        ])
        self.parent_window.gm_require_gly_pos.setToolTip("Require glycine at specific position(s) in the 8-residue window")
        self.parent_window.gm_require_gly_pos.setMinimumHeight(32)
        gm_grid.addWidget(self.parent_window.gm_require_gly_pos, 1, 3)

        # Row 2: Output CSV (full width)
        gm_grid.addWidget(QLabel("Output CSV:"), 2, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.gm_out_csv = QLineEdit(); self.parent_window.gm_out_csv.setMinimumHeight(32)
        self.parent_window.gm_out_browse = QPushButton(t("browse")); self.parent_window.gm_out_browse.setMinimumHeight(32); self.parent_window.gm_out_browse.clicked.connect(self.browse_gm_out_csv)
        r2 = QHBoxLayout(); r2.addWidget(self.parent_window.gm_out_csv, 1); r2.addWidget(self.parent_window.gm_out_browse)
        gm_grid.addLayout(r2, 2, 1, 1, 3)

        # Row 3: Checkboxes | Max Patches
        self.parent_window.gm_highlight_surface = QCheckBox("Highlight Surface")
        self.parent_window.gm_highlight_surface.setChecked(True)
        self.parent_window.gm_highlight_surface.setToolTip("Highlight molecular surface of G-loop region")
        gm_grid.addWidget(self.parent_window.gm_highlight_surface, 3, 1)

        self.parent_window.gm_export_coords = QCheckBox("Export Coordinates")
        self.parent_window.gm_export_coords.setChecked(False)
        self.parent_window.gm_export_coords.setToolTip("Export G-loop atom coordinates for downstream analysis")
        gm_grid.addWidget(self.parent_window.gm_export_coords, 3, 3)

        # Row 4: Max Patches
        gm_grid.addWidget(QLabel("Max Patches:"), 4, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.gm_max_patches = QComboBox()
        self.parent_window.gm_max_patches.addItems(["5", "10", "15", "20", "All"])
        self.parent_window.gm_max_patches.setCurrentIndex(1)  # Default: 10
        self.parent_window.gm_max_patches.setMinimumHeight(32)
        self.parent_window.gm_max_patches.setToolTip("Maximum number of surface patches to display in analysis")
        gm_grid.addWidget(self.parent_window.gm_max_patches, 4, 1)

        gm_layout.addLayout(gm_grid)

        # Buttons inside card
        gm_btn_row = QHBoxLayout()
        gm_btn_row.setSpacing(10)

        self.parent_window.gm_btn = QPushButton("Detect POI")
        self.parent_window.gm_btn.setMinimumHeight(36)
        self.parent_window.gm_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.gm_btn.clicked.connect(self.start_gmotif)
        gm_btn_row.addWidget(self.parent_window.gm_btn)
        gm_btn_row.addStretch(1)

        gm_layout.addLayout(gm_btn_row)

        layout.addWidget(grp_gm)

        # 4. Surface Analysis
        grp_surf = QFrame()
        grp_surf.setStyleSheet(self._get_card_style(is_dark))
        surf_layout = QVBoxLayout(grp_surf)
        surf_layout.setSpacing(12)
        surf_layout.setContentsMargins(16, 14, 16, 14)

        surf_title = QLabel("Protein Surface Analysis")
        surf_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        surf_layout.addWidget(surf_title)

        surf_grid = QGridLayout()
        surf_grid.setContentsMargins(12, 16, 12, 8)
        surf_grid.setColumnStretch(1, 1); surf_grid.setColumnStretch(3, 1)
        surf_grid.setHorizontalSpacing(12); surf_grid.setVerticalSpacing(4)
        
        
        # Row 0: Target Object
        surf_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.obj_combo_surf = QComboBox(); self.parent_window.obj_combo_surf.setMinimumHeight(32)
        self.parent_window.refresh_obj_surf = QPushButton(t("refresh")); self.parent_window.refresh_obj_surf.setMinimumHeight(32); self.parent_window.refresh_obj_surf.clicked.connect(self.refresh_objects)
        r0_surf = QHBoxLayout(); r0_surf.addWidget(self.parent_window.obj_combo_surf, 1); r0_surf.addWidget(self.parent_window.refresh_obj_surf)
        surf_grid.addLayout(r0_surf, 0, 1)

        # Row 0: Output CSV
        surf_grid.addWidget(QLabel("Output CSV:"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.surf_out_csv = QLineEdit(); self.parent_window.surf_out_csv.setMinimumHeight(32)
        self.parent_window.surf_out_browse = QPushButton(t("browse")); self.parent_window.surf_out_browse.setMinimumHeight(32)
        self.parent_window.surf_out_browse.clicked.connect(self.browse_surf_out_csv)
        r0b_surf = QHBoxLayout(); r0b_surf.addWidget(self.parent_window.surf_out_csv, 1); r0b_surf.addWidget(self.parent_window.surf_out_browse)
        surf_grid.addLayout(r0b_surf, 0, 3)

        # Row 1: pH option (APBS is always used)
        surf_grid.addWidget(QLabel("pH:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.surf_ph = QLineEdit("7.4")
        self.parent_window.surf_ph.setMinimumHeight(32)
        self.parent_window.surf_ph.setToolTip("pH value for PDB2PQR protonation state (default: 7.4)")
        surf_grid.addWidget(self.parent_window.surf_ph, 1, 1)
        
        # Row 1 continued: Surface color mode
        surf_grid.addWidget(QLabel("Surface Property:"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.surf_color_mode = QComboBox()
        self.parent_window.surf_color_mode.setMinimumHeight(32)
        self.parent_window.surf_color_mode.addItems([
            "Electrostatic Potential",
            "Hydrophobicity"
        ])
        self.parent_window.surf_color_mode.setToolTip(
            "Select surface coloring mode:\n"
            "• Electrostatic Potential: Red (negative) → White (neutral) → Blue (positive)\n"
            "• Hydrophobicity: White (hydrophilic) → Green (hydrophobic)"
        )
        surf_grid.addWidget(self.parent_window.surf_color_mode, 1, 3)

        surf_layout.addLayout(surf_grid)

        # Buttons - using Ternary Evaluation 风格
        surf_btn_row = QHBoxLayout()
        surf_btn_row.setSpacing(10)

        self.parent_window.surf_btn = QPushButton("Analyze Surface")
        self.parent_window.surf_btn.setMinimumHeight(36)
        self.parent_window.surf_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.surf_btn.clicked.connect(self.start_surface_analysis)

        self.parent_window.surf_vis_btn = QPushButton("Visualize Patches")
        self.parent_window.surf_vis_btn.setMinimumHeight(36)
        self.parent_window.surf_vis_btn.setStyleSheet(self._get_green_btn_style())
        self.parent_window.surf_vis_btn.clicked.connect(self.render_surface_patches)

        self.parent_window.surf_full_btn = QPushButton("Full Surface")
        self.parent_window.surf_full_btn.setMinimumHeight(36)
        self.parent_window.surf_full_btn.setStyleSheet(self._get_purple_btn_style())
        self.parent_window.surf_full_btn.setToolTip(
            "Render complete protein surface with continuous color mapping.\n"
            "Red = negative potential, White = neutral, Blue = positive potential"
        )
        self.parent_window.surf_full_btn.clicked.connect(self.render_full_surface)

        surf_btn_row.addWidget(self.parent_window.surf_btn)
        surf_btn_row.addWidget(self.parent_window.surf_vis_btn)
        surf_btn_row.addWidget(self.parent_window.surf_full_btn)
        surf_btn_row.addStretch(1)

        layout.addWidget(grp_surf)
        layout.addLayout(surf_btn_row)

        # 5. Surface Similarity & Complementarity Analysis
        grp_sim = QFrame()
        grp_sim.setStyleSheet(self._get_card_style(is_dark))
        sim_layout = QVBoxLayout(grp_sim)
        sim_layout.setSpacing(12)
        sim_layout.setContentsMargins(16, 14, 16, 14)

        sim_title = QLabel("Surface Similarity & Complementarity")
        sim_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        sim_layout.addWidget(sim_title)

        sim_grid = QGridLayout()
        sim_grid.setContentsMargins(12, 16, 12, 8)
        sim_grid.setColumnStretch(1, 1); sim_grid.setColumnStretch(3, 1)
        sim_grid.setHorizontalSpacing(12); sim_grid.setVerticalSpacing(4)
        
        
        # Row 0: Object 1 / Object 2
        sim_grid.addWidget(QLabel("Object 1:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.obj_combo_sim1 = QComboBox(); self.parent_window.obj_combo_sim1.setMinimumHeight(32)
        self.parent_window.refresh_obj_sim1 = QPushButton(t("refresh")); self.parent_window.refresh_obj_sim1.setMinimumHeight(32); self.parent_window.refresh_obj_sim1.clicked.connect(self.refresh_objects)
        r0_sim = QHBoxLayout(); r0_sim.addWidget(self.parent_window.obj_combo_sim1, 1); r0_sim.addWidget(self.parent_window.refresh_obj_sim1)
        sim_grid.addLayout(r0_sim, 0, 1)
        
        sim_grid.addWidget(QLabel("Object 2:"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.obj_combo_sim2 = QComboBox(); self.parent_window.obj_combo_sim2.setMinimumHeight(32)
        self.parent_window.obj_combo_sim2.addItem("(None - Single Surface)")
        sim_grid.addWidget(self.parent_window.obj_combo_sim2, 0, 3)
        
        # Row 1: Selection 1 / Selection 2
        sim_grid.addWidget(QLabel("Selection 1:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.sim_sel1 = QLineEdit("all"); self.parent_window.sim_sel1.setMinimumHeight(32)
        self.parent_window.sim_sel1.setToolTip(
            "Selection for Object 1. Examples:\n"
            "  • 'all' - entire structure\n"
            "  • 'A' or 'chain A' - chain A\n"
            "  • 'A B' or 'A+B' - chains A and B\n"
            "  • 'resi 1-100' - residues 1-100"
        )
        sim_grid.addWidget(self.parent_window.sim_sel1, 1, 1)

        sim_grid.addWidget(QLabel("Selection 2:"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.sim_sel2 = QLineEdit("all"); self.parent_window.sim_sel2.setMinimumHeight(32)
        self.parent_window.sim_sel2.setToolTip(
            "Selection for Object 2. Examples:\n"
            "  • 'all' - entire structure\n"
            "  • 'B' or 'chain B' - chain B\n"
            "  • 'C D' or 'C+D' - chains C and D\n"
            "  • 'resi 1-100' - residues 1-100"
        )
        sim_grid.addWidget(self.parent_window.sim_sel2, 1, 3)
        
        # Row 2: Analysis Type / Surface Method
        sim_grid.addWidget(QLabel("Analysis Type:"), 2, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.sim_analysis_type = QComboBox(); self.parent_window.sim_analysis_type.setMinimumHeight(32)
        self.parent_window.sim_analysis_type.addItems(["Similarity Search", "Complementarity (PPI)"])
        self.parent_window.sim_analysis_type.setToolTip(
            "Similarity Search: Compare surface features between two different structures.\n"
            "  - Use for finding similar binding sites across proteins\n"
            "  - Objects can be from different PDB files\n\n"
            "Complementarity (PPI): Analyze interface fit between receptor and ligand.\n"
            "  - Use for chains within the SAME complex (e.g., 'chain A' vs 'chain B')\n"
            "  - Objects must be spatially close (in contact)"
        )
        sim_grid.addWidget(self.parent_window.sim_analysis_type, 2, 1)
        
        sim_grid.addWidget(QLabel("Surface Method:"), 2, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.sim_surface_method = QComboBox(); self.parent_window.sim_surface_method.setMinimumHeight(32)
        self.parent_window.sim_surface_method.addItems(["auto", "msms", "open3d", "edtsurf"])
        self.parent_window.sim_surface_method.setToolTip("auto: use best available\nmsms: most accurate (requires MSMS)\nopen3d: good quality (requires open3d)\nedtsurf: built-in fallback")
        sim_grid.addWidget(self.parent_window.sim_surface_method, 2, 3)
        
        # Row 3: Patch Radius / Interface Distance
        sim_grid.addWidget(QLabel("Patch Radius (Å):"), 3, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.sim_patch_radius = QLineEdit("12.0"); self.parent_window.sim_patch_radius.setMinimumHeight(32)
        self.parent_window.sim_patch_radius.setToolTip("Radius of surface patches for comparison (default: 12 Å)")
        sim_grid.addWidget(self.parent_window.sim_patch_radius, 3, 1)
        
        sim_grid.addWidget(QLabel("Interface Dist (Å):"), 3, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.sim_interface_dist = QLineEdit("4.0"); self.parent_window.sim_interface_dist.setMinimumHeight(32)
        self.parent_window.sim_interface_dist.setToolTip("Distance threshold for interface contacts (default: 4 Å)")
        sim_grid.addWidget(self.parent_window.sim_interface_dist, 3, 3)
        
        # Row 4: Output CSV / Use APBS
        sim_grid.addWidget(QLabel("Output CSV:"), 4, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.sim_out_csv = QLineEdit(); self.parent_window.sim_out_csv.setMinimumHeight(32)
        self.parent_window.sim_out_browse = QPushButton(t("browse")); self.parent_window.sim_out_browse.setMinimumHeight(32); self.parent_window.sim_out_browse.clicked.connect(self.browse_sim_out_csv)
        r4_sim = QHBoxLayout(); r4_sim.addWidget(self.parent_window.sim_out_csv, 1); r4_sim.addWidget(self.parent_window.sim_out_browse)
        sim_grid.addLayout(r4_sim, 4, 1)
        
        self.parent_window.sim_use_apbs = QCheckBox("Use APBS (accurate ESP)")
        self.parent_window.sim_use_apbs.setChecked(False)
        self.parent_window.sim_use_apbs.setToolTip("Use APBS for accurate electrostatics (slower, requires APBS)")
        sim_grid.addWidget(self.parent_window.sim_use_apbs, 4, 3)
        
        # Row 5: Include Ligands option
        self.parent_window.sim_include_ligands = QCheckBox("Include Ligands")
        self.parent_window.sim_include_ligands.setChecked(True)
        self.parent_window.sim_include_ligands.setToolTip("Include small molecules (HETATM) in surface analysis using atom-type based features")
        sim_grid.addWidget(self.parent_window.sim_include_ligands, 5, 1)

        sim_layout.addLayout(sim_grid)

        # Buttons - using Ternary Evaluation 风格
        sim_btn_row = QHBoxLayout()
        sim_btn_row.setSpacing(10)

        self.parent_window.sim_analyze_btn = QPushButton("Analyze Surface")
        self.parent_window.sim_analyze_btn.setMinimumHeight(36)
        self.parent_window.sim_analyze_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.sim_analyze_btn.clicked.connect(self.start_surface_similarity)

        self.parent_window.sim_compare_btn = QPushButton("Compare / Complementarity")
        self.parent_window.sim_compare_btn.setMinimumHeight(36)
        self.parent_window.sim_compare_btn.setStyleSheet(self._get_purple_btn_style())
        self.parent_window.sim_compare_btn.clicked.connect(self.start_surface_comparison)

        self.parent_window.sim_visualize_btn = QPushButton("Visualize Features")
        self.parent_window.sim_visualize_btn.setMinimumHeight(36)
        self.parent_window.sim_visualize_btn.setStyleSheet(self._get_green_btn_style())
        self.parent_window.sim_visualize_btn.clicked.connect(self.visualize_surface_features)
        
        sim_btn_row.addWidget(self.parent_window.sim_analyze_btn)
        sim_btn_row.addWidget(self.parent_window.sim_compare_btn)
        sim_btn_row.addWidget(self.parent_window.sim_visualize_btn)
        sim_btn_row.addStretch(1)
        
        layout.addWidget(grp_sim)
        layout.addLayout(sim_btn_row)
        
        # Add stretch at the end to prevent compression
        layout.addStretch(1)

    # --- G-Motif Logic ---
    def browse_gm_out_csv(self):
        fn, _ = QFileDialog.getSaveFileName(self, t("select_outcsv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.parent_window.gm_out_csv.setText(fn); self.parent_window.update_enablement()

    def start_gmotif(self):
        obj = self.parent_window.obj_combo_gm.currentText().strip()
        pdb = None
        outcsv = self.parent_window.gm_out_csv.text().strip() or None
        if not obj or obj == t("no_object"):
            show_message_box(self, t("title"), t("no_object"), "warning"); return
        try:
            rmsd = float(self.parent_window.gm_rmsd.text().strip() or "3.5")
        except Exception:
            rmsd = 3.5
        
        # Parse glycine position requirement from dropdown
        gly_pos_idx = self.parent_window.gm_require_gly_pos.currentIndex()
        require_gly_pos_map = {
            0: "6,3",   # Pos 6 or 3 (default)
            1: "6",     # Pos 6 only
            2: "3",     # Pos 3 only
            3: None     # No requirement
        }
        require_gly_pos = require_gly_pos_map.get(gly_pos_idx, "6,3")
        
        # New options
        highlight_surface = self.parent_window.gm_highlight_surface.isChecked()
        export_coords = self.parent_window.gm_export_coords.isChecked()

        idx = self.parent_window.gm_template_mode.currentIndex()
        # Template options: 0=GSPT1, 1=CK1α, 2=From Selection
        if idx in (0, 1):
            template_mode = "builtin"
            template_sel = None
            # Use simplified template names (consistent with BUILTIN_TEMPLATES in g_motif_analyzer.py)
            template_builtin = ["GSPT1", "CK1α"][idx]
        else:
            # idx == 2: From Selection
            template_mode = "selection"
            template_sel = None
            template_builtin = None

        self.parent_window.gm_btn.setEnabled(False)
        self.parent_window.progress_bar.setVisible(True); self.parent_window.progress_bar.setRange(0, 0)
        self.parent_window.gmotif_thread = GMotifWorker(obj, pdb, rmsd, require_gly_pos, outcsv,
                                        template_mode, template_sel, template_builtin,
                                        highlight_surface, export_coords)
        self.parent_window.gmotif_thread.progress.connect(self.log)
        self.parent_window.gmotif_thread.error.connect(self.on_error)
        self.parent_window.gmotif_thread.finished.connect(self.on_finished_gmotif)
        self.parent_window.gmotif_thread.start()

    def on_finished_gmotif(self, hits: List[Tuple], out_csv_path: str, coords_data: dict = None, surface_info: dict = None):
        self._gmotif_hits = hits or []
        self._last_gmotif_csv = out_csv_path
        self._last_gmotif_coords = coords_data
        self._last_gmotif_surface = surface_info
        
        self.log(f"✅ G-Motif detection complete: {len(hits)} hits found")
        
        # Print sequence for each hit
        if hits:
            for i, hit in enumerate(hits):
                ch, resi_s, resi_e, seq8, rmsd = hit
                self.log(f"   Hit {i+1}: Chain {ch}, Resi {resi_s}-{resi_e}, Seq: {seq8}, RMSD: {rmsd:.2f}Å")
        
        if os.path.exists(out_csv_path):
            self.parent_window.gm_out_csv.setText(out_csv_path)
            self.log(f"   Saved to: {os.path.basename(out_csv_path)}")
        
        # Show surface info
        if surface_info:
            self.log(f"   Surface: {surface_info.get('atom_count', 0)} atoms, "
                    f"area={surface_info.get('surface_area', 0):.1f} Å²")
        
        # Show coordinate info
        if coords_data:
            centroid = coords_data.get('centroid', (0, 0, 0))
            n_ca = len(coords_data.get('ca_coords', []))
            self.log(f"   Coordinates: {n_ca} Cα atoms, centroid=({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f})")
            if coords_data.get('csv_path'):
                self.log(f"   Coords CSV: {os.path.basename(coords_data['csv_path'])}")
        
        self.parent_window.progress_bar.setVisible(False); self.parent_window.progress_bar.setRange(0, 1)
        self.parent_window.gm_btn.setEnabled(True)

    def show_gloop_surface(self):
        """Manually show G-loop surface"""
        if not self._gmotif_hits:
            show_message_box(self, "Warning", "No G-loop hits found. Please run detection first.", "warning")
            return

        try:
            from pymol import cmd
            try:
                from ...highlight_residues import highlight_gloop_surface
            except ImportError:
                try:
                    from highlight_residues import highlight_gloop_surface
                except ImportError:
                    show_message_box(self, "Error", "highlight_gloop_surface not available", "warning")
                    return
            
            obj = self.parent_window.obj_combo_gm.currentText().strip()
            
            # Highlight surface for all hits
            for i, hit in enumerate(self._gmotif_hits):
                ch, resi_s, resi_e, seq8, rmsd = hit
                result = highlight_gloop_surface(
                    obj=obj, chain=ch, start_resi=resi_s, end_resi=resi_e,
                    surface_color="yellow", surface_transparency=0.3,
                    selection_name=f"gloop_surf_{i+1}",
                    clear_old=(i == 0)  # Only clear old on first one
                )
                if result:
                    self.log(f"✅ Surface highlighted: {ch}:{resi_s}-{resi_e} ({seq8})")
            
            # Zoom to all G-loops
            if self._gmotif_hits:
                all_sel = " or ".join([f"gloop_surf_{i+1}" for i in range(len(self._gmotif_hits))])
                cmd.zoom(all_sel, buffer=8.0)
                
        except Exception as e:
            self.on_error(f"Surface highlight failed: {e}")

    def analyze_gloop_surface(self):
        """Analyze protein surface properties around G-loop patches"""
        if not self._gmotif_hits:
            show_message_box(self, "Warning", "No G-loop hits found. Please run detection first.", "warning")
            return
        
        try:
            from pymol import cmd

            obj = self.parent_window.obj_combo_gm.currentText().strip()
            if not obj or obj == t("no_object"):
                show_message_box(self, t("title"), t("no_object"), "warning")
                return

            self.log("🔬 Starting G-loop surface analysis...")

            # Import surface analyzer
            try:
                from ..protein_surface_analyzer import SurfaceAnalyzer, SurfacePatch
            except ImportError:
                try:
                    from protein_surface_analyzer import SurfaceAnalyzer, SurfacePatch
                except ImportError:
                    show_message_box(self, "Error", "SurfaceAnalyzer not available", "warning")
                    return

            # Build selection for G-loop regions + surrounding area
            gloop_regions = []
            for hit in self._gmotif_hits:
                ch, resi_s, resi_e, seq8, rmsd = hit
                gloop_regions.append((ch, resi_s, resi_e))

            if not gloop_regions:
                show_message_box(self, "Warning", "No valid G-loop regions found.", "warning")
                return
            
            # Create selection for G-loop + 10Å surrounding
            gloop_selections = [f"(chain {ch} and resi {st}-{ed})" for ch, st, ed in gloop_regions]
            all_gloop_sel = " or ".join(gloop_selections)
            
            # Create a temporary object with just the G-loop region + surrounding
            gloop_env_name = f"{obj}_gloop_env"
            cmd.select("_gloop_core", f"{obj} and ({all_gloop_sel})")
            cmd.select("_gloop_env", f"byres ({obj} within 12 of _gloop_core)")
            cmd.create(gloop_env_name, "_gloop_env")
            cmd.delete("_gloop_core")
            cmd.delete("_gloop_env")
            
            self.log(f"   Created analysis region: {gloop_env_name}")
            
            # Run surface analysis on the G-loop environment
            analyzer = SurfaceAnalyzer(gloop_env_name, grid_spacing=0.5)
            patches = analyzer.analyze()
            
            # Store results
            self._gloop_surface_patches = patches
            
            # Summarize results
            n_electro_pos = sum(1 for p in patches if p.type == "electrostatic_pos")
            n_electro_neg = sum(1 for p in patches if p.type == "electrostatic_neg")
            n_hydrophobic = sum(1 for p in patches if p.type == "hydrophobic")
            
            self.log(f"✅ G-loop surface analysis complete:")
            self.log(f"   Total patches: {len(patches)}")
            self.log(f"   • Positive electrostatic: {n_electro_pos}")
            self.log(f"   • Negative electrostatic: {n_electro_neg}")
            self.log(f"   • Hydrophobic: {n_hydrophobic}")
            
            # Get max patches setting
            max_patches_text = self.parent_window.gm_max_patches.currentText()
            if max_patches_text == "All":
                max_patches = len(patches)
            else:
                max_patches = int(max_patches_text)
            
            # Report top patches by score
            if patches:
                sorted_patches = sorted(patches, key=lambda p: p.score, reverse=True)
                display_patches = sorted_patches[:max_patches]
                self.log(f"   Top {len(display_patches)} patches by druggability score:")
                for i, p in enumerate(display_patches):
                    res_str = ", ".join([f"{r['chain']}:{r['resn']}{r['resi']}" for r in p.residues[:3]])
                    if len(p.residues) > 3:
                        res_str += f"... (+{len(p.residues)-3} more)"
                    self.log(f"     {i+1}. {p.type}: area={p.area:.1f}Å², score={p.score:.2f}")
                    self.log(f"        Residues: {res_str}")
                
                # Store only the top patches for visualization
                self._gloop_surface_patches = display_patches
            
            # Visualize patches on the G-loop environment (use filtered patches)
            self._visualize_gloop_patches(gloop_env_name, self._gloop_surface_patches)
            
            # Zoom to the analysis region
            cmd.zoom(gloop_env_name, buffer=5.0)
            
            # Optionally export to CSV
            if self.parent_window.gm_out_csv.text().strip():
                base_csv = self.parent_window.gm_out_csv.text().strip()
                surf_csv = base_csv.replace(".csv", "_surface_patches.csv")
                self._export_gloop_patches_csv(patches, surf_csv)
                self.log(f"   Surface patches saved to: {os.path.basename(surf_csv)}")
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.on_error(f"G-loop surface analysis failed: {e}")

    def _visualize_gloop_patches(self, obj_name: str, patches: list):
        """Visualize surface patches on the G-loop environment"""
        try:
            from pymol import cmd
            
            # Prepare main object visualization
            cmd.show("cartoon", obj_name)
            cmd.hide("surface", obj_name)
            
            # Delete old patch objects
            for name in cmd.get_names("objects"):
                if name.startswith("patch_"):
                    cmd.delete(name)
            
            # Create separate objects for each patch
            patch_names = []
            for p in patches:
                if not p.residues:
                    continue
                
                # Build selection for patch residues
                sel_parts = [f"(chain {r['chain']} and resi {r['resi']})" for r in p.residues]
                sel_str = " or ".join(sel_parts)
                patch_name = f"patch_{p.id}_{p.type[:4]}"
                patch_names.append(patch_name)
                
                # Create separate object
                cmd.create(patch_name, f"{obj_name} and ({sel_str})")
                
                # Set up visualization
                cmd.hide("everything", patch_name)
                cmd.show("surface", patch_name)
                cmd.set("transparency", 0.3, patch_name)
                
                # Color based on type
                if p.type == "electrostatic_pos":
                    cmd.color("blue", patch_name)
                elif p.type == "electrostatic_neg":
                    cmd.color("red", patch_name)
                elif p.type == "hydrophobic":
                    cmd.color("green", patch_name)
                
                # Reduce transparency for patches
                cmd.set("transparency", 0.2, patch_name)
            
            # Group patches
            if patch_names:
                cmd.group("gloop_patches", "patch_*")
            
            # Also show the G-loop as sticks
            if self._gmotif_hits:
                for i, hit in enumerate(self._gmotif_hits):
                    ch, resi_s, resi_e, seq8, rmsd = hit
                    gloop_sel = f"{obj_name} and chain {ch} and resi {resi_s}-{resi_e}"
                    cmd.show("sticks", gloop_sel)
                    cmd.color("yellow", gloop_sel)
                    cmd.set("stick_radius", 0.2, gloop_sel)
            
            self.log(f"   Visualized {len(patches)} top patches on {obj_name}")
            self.log("   Colors: Blue=positive, Red=negative, Green=hydrophobic")
            for p in patches:
                self.log(f"     • patch_{p.id}_{p.type[:4]}: score={p.score:.2f}")
            
        except Exception as e:
            self.log(f"   Visualization warning: {e}")

    def _export_gloop_patches_csv(self, patches: list, csv_path: str):
        """Export G-loop surface patches to CSV"""
        import csv
        
        with open(csv_path, 'w', newline='', encoding='utf-8') as f:
            fieldnames = ['Patch_ID', 'Type', 'Area_A2', 'Score', 'Center_X', 'Center_Y', 'Center_Z',
                         'Avg_Potential', 'Avg_Hydrophobicity', 'Shape_Index', 'N_Residues', 'Residues']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for p in patches:
                res_str = ";".join([f"{r['chain']}:{r['resn']}{r['resi']}" for r in p.residues])
                writer.writerow({
                    'Patch_ID': p.id,
                    'Type': p.type,
                    'Area_A2': f"{p.area:.2f}",
                    'Score': f"{p.score:.3f}",
                    'Center_X': f"{p.center[0]:.2f}",
                    'Center_Y': f"{p.center[1]:.2f}",
                    'Center_Z': f"{p.center[2]:.2f}",
                    'Avg_Potential': f"{p.avg_potential:.3f}",
                    'Avg_Hydrophobicity': f"{p.avg_hydrophobicity:.3f}",
                    'Shape_Index': f"{p.shape_index:.3f}",
                    'N_Residues': len(p.residues),
                    'Residues': res_str
                })

    # --- Surface Analysis Logic ---
    def browse_surf_out_csv(self):
        fn, _ = QFileDialog.getSaveFileName(self, t("select_outcsv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.parent_window.surf_out_csv.setText(fn)

    def start_surface_analysis(self):
        """Start Surface Analysis"""
        obj = self.parent_window.obj_combo_surf.currentText().strip()
        outcsv = self.parent_window.surf_out_csv.text().strip() or None

        if not obj or obj == t("no_object"):
            show_message_box(self, t("title"), t("no_object"), "warning"); return
        
        # Get pH value (APBS is always used)
        try:
            ph = float(self.parent_window.surf_ph.text().strip() or "7.4")
        except (ValueError, TypeError):  # float() 转换可能Failed
            ph = 7.4
            
        # Get surface property mode and map to internal key
        color_mode_idx = self.parent_window.surf_color_mode.currentIndex()
        # 中文注释：Analyze Surface 需要根据用户选择筛选 patch 类型，避免两个选项返回相同结果
        surface_property = 'electrostatic' if color_mode_idx == 0 else 'hydrophobicity'

        self.parent_window.surf_btn.setEnabled(False)
        self.parent_window.progress_bar.setVisible(True); self.parent_window.progress_bar.setRange(0, 0)

        self.log(f"🔬 Starting surface analysis with APBS (pH={ph})...")
        self.log(f"   Surface property: {surface_property}")

        self.parent_window.surf_thread = SurfaceAnalysisWorker(
            obj,
            outcsv,
            use_apbs=True,
            ph=ph,
            surface_property=surface_property,
        )
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
                show_message_box(self, "Warning", "No analysis results to visualize. Please run analysis first.", "warning")
                return

            # Prepare main object visualization (show cartoon, hide old patches)
            cmd.show("cartoon", obj)
            cmd.hide("surface", obj)
            cmd.hide("lines", obj)
            cmd.hide("sticks", obj)
            
            # Delete old patch objects
            for name in cmd.get_names("objects"):
                if name.startswith("patch_"):
                    cmd.delete(name)
            
            # Create separate objects for each patch
            patch_names = []
            for p in self._last_surf_patches:
                # Use residues to select patch area
                if not p.residues: continue
                
                sel_str = " or ".join([f"(chain {r['chain']} and resi {r['resi']})" for r in p.residues])
                patch_name = f"patch_{p.id}_{p.type}"
                patch_names.append(patch_name)
                
                # Create a new object from the selection
                # We use 'byres' to ensure we get full residues for the surface
                cmd.create(patch_name, f"byres ({obj} and ({sel_str}))")
                
                # Set up visualization for this patch object
                cmd.hide("everything", patch_name)
                cmd.show("surface", patch_name)
                cmd.set("transparency", 0.3, patch_name)
                
                # Color based on type
                color = "red" # neg
                if p.type == "electrostatic_pos": color = "blue"
                elif p.type == "hydrophobic": color = "green"
                elif p.type == "electrostatic_neg": color = "red"
                
                cmd.color(color, patch_name)
                
            # Group patches for cleaner object list
            if patch_names:
                cmd.group("surface_patches", "patch_*")
                
            self.log(f"Visualized {len(self._last_surf_patches)} patches on {obj} (created separate objects)")
        
        except Exception as e:
            self.on_error(str(e))

    def render_full_surface(self):
        """Render complete protein surface with continuous color mapping"""
        try:
            from pymol import cmd
            obj = self.parent_window.obj_combo_surf.currentText().strip()

            if not obj or obj == t("no_object"):
                show_message_box(self, "Warning", "Please select an object first.", "warning")
                return

            # Get pH value (APBS is always used)
            try:
                ph = float(self.parent_window.surf_ph.text().strip() or "7.4")
            except (ValueError, TypeError):  # float() 转换可能Failed
                ph = 7.4
            
            # Get color mode: 0=potential, 1=hydrophobicity
            color_mode_idx = self.parent_window.surf_color_mode.currentIndex()
            color_modes = ['potential', 'hydrophobicity']
            color_by = color_modes[color_mode_idx]

            self.log(f"🎨 Rendering full surface for {obj}...")
            if color_by == 'potential':
                self.log(f"   Using APBS electrostatics (pH={ph})")
            elif color_by == 'hydrophobicity':
                self.log("   Coloring by hydrophobicity (Kyte-Doolittle scale)")

            # Import and run
            try:
                from ...protein_surface_analyzer import render_protein_surface
            except ImportError:
                try:
                    from protein_surface_analyzer import render_protein_surface
                except ImportError:
                    show_message_box(self, "Error", "render_protein_surface not available", "warning")
                    return

            # Check if output directory is specified (use Output CSV directory)
            output_csv = self.parent_window.surf_out_csv.text().strip()
            save_files = False
            output_dir = None
            if output_csv:
                output_dir = os.path.dirname(output_csv)
                if output_dir:
                    save_files = True
                    self.log(f"   Output files will be saved to: {output_dir}")

            success = render_protein_surface(
                obj_name=obj,
                use_apbs=True,  # Always use APBS
                ph=ph,
                color_by=color_by,
                transparency=0.0,
                save_files=save_files,
                output_dir=output_dir
            )

            if success:
                self.log(f"✅ Full surface rendered for {obj}")
                if color_by == 'potential':
                    self.log("   Color scheme: Red (negative) → White (neutral) → Blue (positive)")
                    # ESP group contains: {obj}_ESP_surf, {obj}_ESP_map, {obj}_ESP_ramp
                    surface_name = f"{obj}_ESP_surf"
                elif color_by == 'hydrophobicity':
                    self.log("   Color scheme: White (hydrophilic) → Green (hydrophobic)")
                    surface_name = f"{obj}_hydro"
                else:
                    surface_name = f"{obj}_surface"
                # Zoom to the surface
                cmd.zoom(surface_name, buffer=5.0)
            else:
                self.log("❌ Surface rendering failed")

        except Exception as e:
            import traceback
            traceback.print_exc()
            self.on_error(f"Full surface rendering failed: {e}")

    # --- Surface Similarity & Complementarity Logic ---
    def browse_sim_out_csv(self):
        fn, _ = QFileDialog.getSaveFileName(self, t("select_outcsv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.parent_window.sim_out_csv.setText(fn)

    def start_surface_similarity(self):
        """Start single surface analysis"""
        obj1 = self.parent_window.obj_combo_sim1.currentText().strip()
        if not obj1 or obj1 == t("no_object"):
            show_message_box(self, t("title"), t("no_object"), "warning"); return
        
        sel1 = self.parent_window.sim_sel1.text().strip() or "all"
        outcsv = self.parent_window.sim_out_csv.text().strip() or None
        surface_method = self.parent_window.sim_surface_method.currentText()
        use_apbs = self.parent_window.sim_use_apbs.isChecked()
        
        try:
            patch_radius = float(self.parent_window.sim_patch_radius.text().strip() or "12.0")
        except (ValueError, TypeError):  # float() 转换可能Failed
            patch_radius = 12.0
        
        self.parent_window.sim_analyze_btn.setEnabled(False)
        self.parent_window.progress_bar.setVisible(True); self.parent_window.progress_bar.setRange(0, 0)
        
        self.parent_window.sim_thread = SurfaceSimilarityWorker(
            obj1=obj1, obj2=None,
            selection1=sel1, selection2="all",
            analysis_type="single",
            patch_radius=patch_radius,
            out_csv=outcsv,
            surface_method=surface_method,
            use_apbs=use_apbs
        )
        self.parent_window.sim_thread.progress.connect(self.log)
        self.parent_window.sim_thread.error.connect(self.on_error)
        self.parent_window.sim_thread.finished.connect(self.on_finished_similarity)
        self.parent_window.sim_thread.start()

    def start_surface_comparison(self):
        """Start surface comparison or complementarity analysis
        
        For Similarity Search:
        - Object 1 is the TEMPLATE (e.g., known binding site like CRBN)
        - Object 2 is the TARGET (protein to search for similar sites)
        - The method generates patches from template and searches in target
        """
        obj1 = self.parent_window.obj_combo_sim1.currentText().strip()
        obj2_text = self.parent_window.obj_combo_sim2.currentText().strip()
        
        if not obj1 or obj1 == t("no_object"):
            show_message_box(self, t("title"), t("no_object"), "warning"); return
        
        # Check if obj2 is selected
        if obj2_text == "(None - Single Surface)" or not obj2_text:
            # Single surface analysis
            self.start_surface_similarity()
            return
        
        obj2 = obj2_text
        sel1 = self.parent_window.sim_sel1.text().strip() or "all"
        sel2 = self.parent_window.sim_sel2.text().strip() or "all"
        outcsv = self.parent_window.sim_out_csv.text().strip() or None
        surface_method = self.parent_window.sim_surface_method.currentText()
        
        # Determine analysis type
        analysis_idx = self.parent_window.sim_analysis_type.currentIndex()
        # 0 = Similarity Search (use obj1 as template, search in obj2)
        # 1 = Complementarity (PPI)
        analysis_type = "complementarity" if analysis_idx == 1 else "search"
        
        try:
            patch_radius = float(self.parent_window.sim_patch_radius.text().strip() or "12.0")
        except (ValueError, TypeError):  # float() 转换可能Failed
            patch_radius = 12.0
        
        try:
            interface_dist = float(self.parent_window.sim_interface_dist.text().strip() or "4.0")
        except (ValueError, TypeError):  # float() 转换可能Failed
            interface_dist = 4.0
        
        self.parent_window.sim_compare_btn.setEnabled(False)
        self.parent_window.progress_bar.setVisible(True); self.parent_window.progress_bar.setRange(0, 0)
        
        # Log the search direction for clarity
        if analysis_type == "search":
            self.log(f"🔍 Starting similarity search:")
            self.log(f"   Template (Object 1): {obj1} [{sel1}]")
            self.log(f"   Target (Object 2): {obj2} [{sel2}]")
            self.log(f"   → Searching for patches in {obj2} similar to {obj1}")
        
        self.parent_window.sim_thread = SurfaceSimilarityWorker(
            obj1=obj1, obj2=obj2,
            selection1=sel1, selection2=sel2,
            analysis_type=analysis_type,
            patch_radius=patch_radius,
            interface_distance=interface_dist,
            out_csv=outcsv,
            surface_method=surface_method,
            similarity_threshold=0.5,  # Default threshold
            top_k=5  # Return top 5 matches per template patch
        )
        self.parent_window.sim_thread.progress.connect(self.log)
        self.parent_window.sim_thread.error.connect(self.on_error)
        self.parent_window.sim_thread.finished.connect(self.on_finished_comparison)
        self.parent_window.sim_thread.start()

    def on_finished_similarity(self, result, out_csv_path: str):
        """Handle single surface analysis completion"""
        self._last_similarity_result = result
        
        if isinstance(result, dict):
            # Single surface result
            self.log(f"✅ Surface analysis complete:")
            self.log(f"   Vertices: {result.get('n_vertices', 0)}")
            self.log(f"   Faces: {result.get('n_faces', 0)}")
            self.log(f"   Surface Area: {result.get('surface_area', 0):.1f} Å²")
        
        if out_csv_path and os.path.exists(out_csv_path):
            self.parent_window.sim_out_csv.setText(out_csv_path)
            self.log(f"   Saved to: {os.path.basename(out_csv_path)}")
        
        self.parent_window.progress_bar.setVisible(False); self.parent_window.progress_bar.setRange(0, 1)
        self.parent_window.sim_analyze_btn.setEnabled(True)

    def on_finished_comparison(self, result, out_csv_path: str):
        """Handle comparison/complementarity/search analysis completion"""
        analysis_type = self.parent_window.sim_analysis_type.currentText()
        
        if "Complementarity" in analysis_type:
            self._last_complementarity_result = result
            self.log(f"✅ Complementarity analysis complete:")
            self.log(f"   Overall Score: {result.score:.3f}")
            self.log(f"   Geometric: {result.geometric_complementarity:.3f}")
            self.log(f"   Electrostatic: {result.electrostatic_complementarity:.3f}")
            self.log(f"   Hydrophobic: {result.hydrophobic_complementarity:.3f}")
            self.log(f"   Interface Area: {result.interface_area:.1f} Å²")
            self.log(f"   Contacts: {result.n_contacts}")
        elif hasattr(result, 'n_template_patches'):
            # SimilaritySearchResult - new search logic
            self._last_search_result = result
            self.log(f"✅ Similarity search complete:")
            self.log(f"   Template: {result.template_object}")
            self.log(f"   Target: {result.target_object}")
            self.log(f"   Template patches: {result.n_template_patches}")
            self.log(f"   Target patches: {result.n_target_patches}")
            self.log(f"   Matches found: {result.n_matches_found}")
            self.log(f"   Best similarity: {result.max_similarity:.3f}")
            self.log(f"   Mean best similarity: {result.mean_best_similarity:.3f}")
            
            # Show top matches
            if result.patch_results:
                self.log(f"   Top matching patches:")
                sorted_results = sorted(result.patch_results,
                                        key=lambda x: x.best_match_score, reverse=True)
                for i, pr in enumerate(sorted_results[:5]):
                    if pr.best_match_score > 0:
                        self.log(f"     {i+1}. Template patch {pr.template_patch_idx} → "
                                f"Target patch {pr.best_match_idx} (score: {pr.best_match_score:.3f})")
        else:
            # Legacy SimilarityResult
            self._last_similarity_result = result
            self.log(f"✅ Similarity analysis complete:")
            self.log(f"   Overall Score: {result.score:.3f}")
            self.log(f"   Geometric: {result.geometric_similarity:.3f}")
            self.log(f"   Chemical: {result.chemical_similarity:.3f}")
            self.log(f"   Shape Index Corr: {result.shape_index_correlation:.3f}")
            self.log(f"   ESP Corr: {result.electrostatic_correlation:.3f}")
        
        if out_csv_path and os.path.exists(out_csv_path):
            self.parent_window.sim_out_csv.setText(out_csv_path)
            self.log(f"   Saved to: {os.path.basename(out_csv_path)}")
        
        self.parent_window.progress_bar.setVisible(False); self.parent_window.progress_bar.setRange(0, 1)
        self.parent_window.sim_compare_btn.setEnabled(True)

    def visualize_surface_features(self):
        """Visualize surface features in PyMOL"""
        try:
            from pymol import cmd

            obj1 = self.parent_window.obj_combo_sim1.currentText().strip()
            if not obj1 or obj1 == t("no_object"):
                show_message_box(self, "Warning", "Please select an object first.", "warning")
                return

            # Check if we have analysis results
            if not self._last_similarity_result:
                show_message_box(self, "Warning", "Please run surface analysis first.", "warning")
                return
            
            result = self._last_similarity_result
            
            # If it's a dict with mesh/points, visualize features
            if isinstance(result, dict) and 'points' in result:
                points = result['points']
                
                # Create CGO objects for visualization
                # Color by shape index
                self._visualize_shape_index(obj1, points)
                
                self.log(f"✅ Visualized surface features for {obj1}")
                self.log("   Shape Index: blue (concave) → white (saddle) → red (convex)")
            else:
                # For comparison results, show matched regions
                self.log("Visualization for comparison results not yet implemented")
                
        except Exception as e:
            self.on_error(str(e))

    def _visualize_shape_index(self, obj_name: str, points):
        """Visualize shape index on surface"""
        try:
            from pymol import cmd
            from pymol import cgo

            # Create a color ramp based on shape index
            # SI: -1 (concave/blue) to +1 (convex/red)

            cgo_obj = []
            sphere_radius = 0.3

            for p in points:
                si = p.shape_index

                # Color mapping: -1 -> blue, 0 -> white, +1 -> red
                if si < 0:
                    r = 1.0 + si  # 0 to 1
                    g = 1.0 + si
                    b = 1.0
                else:
                    r = 1.0
                    g = 1.0 - si
                    b = 1.0 - si

                cgo_obj.extend([
                    cgo.COLOR, r, g, b,
                    cgo.SPHERE, p.coord[0], p.coord[1], p.coord[2], sphere_radius
                ])

            cgo_name = f"{obj_name}_shape_index"
            cmd.load_cgo(cgo_obj, cgo_name)

            # Also show the original structure
            cmd.show("cartoon", obj_name)
            cmd.set("cartoon_transparency", 0.7, obj_name)

        except Exception as e:
            self.log(f"Shape index visualization failed: {e}")

    def _get_card_style(self, is_dark: bool) -> str:
        """获取卡片样式"""
        if is_dark:
            return """
                QFrame {
                    background: #1e2530;
                    border: 1px solid #30363d;
                    border-radius: 10px;
                    padding: 12px;
                }
                QLabel {
                    border: none;
                    background: transparent;
                }
            """
        return """
            QFrame {
                background: white;
                border: 1px solid #e2e8f0;
                border-radius: 10px;
                padding: 12px;
            }
            QLabel {
                border: none;
                background: transparent;
            }
        """

    def _create_card(self, title: str, is_dark: bool) -> QFrame:
        """Create带标题的卡片"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(12)
        layout.setContentsMargins(16, 14, 16, 14)

        # 卡片标题
        title_lbl = QLabel(title)
        title_lbl.setStyleSheet("""
            font-size: 15px; font-weight: 600;
            color: #1e293b; padding-bottom: 4px;
        """ if not is_dark else """
            font-size: 15px; font-weight: 600;
            color: #e2e8f0; padding-bottom: 4px;
        """)
        layout.addWidget(title_lbl)

        return card
