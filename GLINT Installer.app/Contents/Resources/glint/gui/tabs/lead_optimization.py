# -*- coding: utf-8 -*-
"""
Lead Optimization: PPI, Glue, Ternary, Mutation, Electrostatic Complementarity
"""
import os
import traceback
from typing import Optional

from ..qt_adapter import (
    Qt, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QCheckBox, QComboBox, QGroupBox, QGridLayout, QScrollArea, QFrame,
    QFileDialog, QMessageBox
)

from ..utils import t, show_message_box, show_question_box
from .common import CommonTab

class LeadOptimizationTab(CommonTab):
    def __init__(self, parent):
        super().__init__(parent)

        # Remove proxies for methods implemented here
        for attr in ['run_mutation', 'run_minimize', 'run_mutation_analysis']:
            if attr in self.__dict__:
                del self.__dict__[attr]

        self.init_ui()

    def init_ui(self):
        """初始化UI - 现代卡片式布局"""
        self.setObjectName("scroll_content")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.parent_window._lead_scroll_content = self

        is_dark = getattr(self.parent_window, "_dark_mode", False)
        bg_color = "#161b22" if is_dark else "#f8fafc"
        self.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")

        layout = QVBoxLayout(self)
        layout.setSpacing(16)
        layout.setContentsMargins(20, 20, 20, 20)

        # === 页面标题 ===
        header = QHBoxLayout()
        title = QLabel("Lead Optimization")
        title.setStyleSheet("""
            font-size: 20px; font-weight: 600;
            color: #3b82f6; padding: 4px 0;
        """)
        header.addWidget(title)
        header.addStretch(1)
        layout.addLayout(header)

        # 1. Protein-Protein Interface (PPI) Analysis
        grp_ppi = QFrame()
        grp_ppi.setStyleSheet(self._get_card_style(is_dark))
        ppi_layout = QVBoxLayout(grp_ppi)
        ppi_layout.setSpacing(12)
        ppi_layout.setContentsMargins(16, 14, 16, 14)

        ppi_title = QLabel("Protein-Protein Interface (PPI) Analysis")
        ppi_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        ppi_layout.addWidget(ppi_title)

        ppi_grid = QGridLayout()
        ppi_grid.setContentsMargins(12, 16, 12, 8)
        ppi_grid.setColumnStretch(1, 1); ppi_grid.setColumnStretch(3, 1)
        ppi_grid.setHorizontalSpacing(12); ppi_grid.setVerticalSpacing(4)
        
        
        ppi_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_obj_combo = QComboBox(); self.parent_window.ppi_obj_combo.setMinimumHeight(32)
        self.parent_window.ppi_refresh_btn = QPushButton(t("refresh")); self.parent_window.ppi_refresh_btn.setMinimumHeight(32); self.parent_window.ppi_refresh_btn.clicked.connect(self.refresh_objects)
        r0 = QHBoxLayout(); r0.addWidget(self.parent_window.ppi_obj_combo, 1); r0.addWidget(self.parent_window.ppi_refresh_btn)
        ppi_grid.addLayout(r0, 0, 1)
        
        ppi_grid.addWidget(QLabel("Protein1 Chains:"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_protein1_chains = QLineEdit(); self.parent_window.ppi_protein1_chains.setPlaceholderText("e.g. A"); self.parent_window.ppi_protein1_chains.setMinimumHeight(32)
        ppi_grid.addWidget(self.parent_window.ppi_protein1_chains, 0, 3)
        
        ppi_grid.addWidget(QLabel("Protein2 Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_protein2_chains = QLineEdit(); self.parent_window.ppi_protein2_chains.setPlaceholderText("e.g. B"); self.parent_window.ppi_protein2_chains.setMinimumHeight(32)
        ppi_grid.addWidget(self.parent_window.ppi_protein2_chains, 1, 1)
        
        ppi_grid.addWidget(QLabel("Interface Dist (Å):"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_interface_dist = QLineEdit("4.5"); self.parent_window.ppi_interface_dist.setMinimumHeight(32)
        ppi_grid.addWidget(self.parent_window.ppi_interface_dist, 1, 3)
        
        ppi_grid.addWidget(QLabel("Output CSV:"), 2, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_csv = QLineEdit(); self.parent_window.ppi_csv.setPlaceholderText("Optional"); self.parent_window.ppi_csv.setMinimumHeight(32)
        self.parent_window.ppi_csv_btn = QPushButton(t("browse")); self.parent_window.ppi_csv_btn.setMinimumHeight(32)
        self.parent_window.ppi_csv_btn.clicked.connect(lambda: self._browse_save_file(self.parent_window.ppi_csv, "CSV (*.csv)"))
        r2 = QHBoxLayout(); r2.addWidget(self.parent_window.ppi_csv, 1); r2.addWidget(self.parent_window.ppi_csv_btn)
        ppi_grid.addLayout(r2, 2, 1, 1, 3)

        # Row 3: Visualization Options
        ppi_grid.addWidget(QLabel("Display Mode:"), 3, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_display_mode = QComboBox()
        self.parent_window.ppi_display_mode.setMinimumHeight(32)
        self.parent_window.ppi_display_mode.addItems(["Cartoon + Surface + Interaction", "Cartoon + Interaction"])
        ppi_grid.addWidget(self.parent_window.ppi_display_mode, 3, 1)
        
        self.parent_window.ppi_show_labels = QCheckBox("Show Distance Labels")
        self.parent_window.ppi_show_labels.setChecked(True)
        ppi_grid.addWidget(self.parent_window.ppi_show_labels, 3, 3)

        ppi_layout.addLayout(ppi_grid)

        ppi_btn_row = QHBoxLayout()
        ppi_btn_row.setSpacing(10)

        self.parent_window.ppi_analyze_btn = QPushButton("Analyze PPI Interface")
        self.parent_window.ppi_analyze_btn.setMinimumHeight(36)
        self.parent_window.ppi_analyze_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.ppi_analyze_btn.clicked.connect(self.run_ppi_analysis)
        ppi_btn_row.addWidget(self.parent_window.ppi_analyze_btn)
        ppi_btn_row.addStretch(1)

        layout.addWidget(grp_ppi)
        layout.addLayout(ppi_btn_row)

        # 2. Protein-Ligand Interactions
        grp_pl = QFrame()
        grp_pl.setStyleSheet(self._get_card_style(is_dark))
        pl_layout = QVBoxLayout(grp_pl)
        pl_layout.setSpacing(12)
        pl_layout.setContentsMargins(16, 14, 16, 14)

        pl_title = QLabel("Protein-Ligand Interactions")
        pl_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        pl_layout.addWidget(pl_title)

        pl_grid = QGridLayout()
        pl_grid.setContentsMargins(12, 16, 12, 8)
        pl_grid.setColumnStretch(1, 1); pl_grid.setColumnStretch(3, 1)
        pl_grid.setHorizontalSpacing(12); pl_grid.setVerticalSpacing(4)
        
        
        pl_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_obj_combo = QComboBox(); self.parent_window.pl_obj_combo.setMinimumWidth(150); self.parent_window.pl_obj_combo.setMinimumHeight(32)
        self.parent_window.pl_refresh_btn = QPushButton(t("refresh")); self.parent_window.pl_refresh_btn.setMinimumHeight(32); self.parent_window.pl_refresh_btn.clicked.connect(self.refresh_objects)
        r0_pl = QHBoxLayout(); r0_pl.addWidget(self.parent_window.pl_obj_combo, 1); r0_pl.addWidget(self.parent_window.pl_refresh_btn)
        pl_grid.addLayout(r0_pl, 0, 1)
        
        pl_grid.addWidget(QLabel("Ligand Name:"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_ligand_name = QLineEdit(); self.parent_window.pl_ligand_name.setPlaceholderText("e.g. CC885, Auto-detect if blank"); self.parent_window.pl_ligand_name.setMinimumHeight(32)
        pl_grid.addWidget(self.parent_window.pl_ligand_name, 0, 3)
        
        pl_grid.addWidget(QLabel("Protein Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_protein_chains = QLineEdit(); self.parent_window.pl_protein_chains.setPlaceholderText("e.g. A,B (optional)"); self.parent_window.pl_protein_chains.setMinimumHeight(32)
        pl_grid.addWidget(self.parent_window.pl_protein_chains, 1, 1)
        
        pl_grid.addWidget(QLabel("Distance (Å):"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_distance = QLineEdit("4.5"); self.parent_window.pl_distance.setMinimumHeight(32)
        pl_grid.addWidget(self.parent_window.pl_distance, 1, 3)
        
        # Output CSV moved to row 3 to make space for 3D options
        # See below for new layout positioning
        
        pl_btn_row = QHBoxLayout()
        pl_btn_row.setSpacing(10)

        self.parent_window.pl_analyze_btn = QPushButton("Analyze Protein-Ligand")
        self.parent_window.pl_analyze_btn.setMinimumHeight(36)
        self.parent_window.pl_analyze_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.pl_analyze_btn.clicked.connect(self.run_pl_analysis)

        self.parent_window.pl_2d_btn = QPushButton("Generate 2D Diagram")
        self.parent_window.pl_2d_btn.setMinimumHeight(36)
        self.parent_window.pl_2d_btn.setStyleSheet(self._get_secondary_btn_style())
        self.parent_window.pl_2d_btn.clicked.connect(self.run_pl_2d_diagram)
        
        pl_btn_row.addWidget(self.parent_window.pl_analyze_btn)
        pl_btn_row.addWidget(self.parent_window.pl_2d_btn)
        pl_btn_row.addStretch(1)
        
        layout.addWidget(grp_pl)
        layout.addLayout(pl_btn_row)
        
        # 3D Visualization Options
        pl_grid.addWidget(QLabel("Min Confidence:"), 2, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_min_confidence = QComboBox(); self.parent_window.pl_min_confidence.setMinimumHeight(32)
        self.parent_window.pl_min_confidence.addItems(["0.0 (Show All)", "0.5", "0.6", "0.7", "0.8 (High)", "0.9"])
        self.parent_window.pl_min_confidence.setCurrentText("0.8 (High)")
        pl_grid.addWidget(self.parent_window.pl_min_confidence, 2, 1)

        self.parent_window.pl_show_hydrophobic = QCheckBox("Show Hydrophobic Interactions")
        self.parent_window.pl_show_hydrophobic.setChecked(False)
        pl_grid.addWidget(self.parent_window.pl_show_hydrophobic, 2, 3)
        
        # Output CSV row
        pl_grid.addWidget(QLabel("Output CSV:"), 3, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        
        self.parent_window.pl_csv = QLineEdit(); self.parent_window.pl_csv.setPlaceholderText("Optional"); self.parent_window.pl_csv.setMinimumHeight(32)
        self.parent_window.pl_csv_btn = QPushButton(t("browse")); self.parent_window.pl_csv_btn.setMinimumHeight(32)
        self.parent_window.pl_csv_btn.clicked.connect(lambda: self._browse_save_file(self.parent_window.pl_csv, "CSV (*.csv)"))
        
        r2_pl = QHBoxLayout(); r2_pl.addWidget(self.parent_window.pl_csv, 1); r2_pl.addWidget(self.parent_window.pl_csv_btn)
        pl_grid.addLayout(r2_pl, 3, 1, 1, 3)

        pl_layout.addLayout(pl_grid)

        # 3. Ligand-Ligand Interactions
        grp_ll = QFrame()
        grp_ll.setStyleSheet(self._get_card_style(is_dark))
        ll_layout = QVBoxLayout(grp_ll)
        ll_layout.setSpacing(12)
        ll_layout.setContentsMargins(16, 14, 16, 14)

        ll_title = QLabel("Ligand-Ligand Interactions (Small Molecule - Small Molecule)")
        ll_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        ll_layout.addWidget(ll_title)

        ll_grid = QGridLayout()
        ll_grid.setContentsMargins(12, 16, 12, 8)
        ll_grid.setColumnStretch(1, 1); ll_grid.setColumnStretch(3, 1)
        ll_grid.setHorizontalSpacing(12); ll_grid.setVerticalSpacing(4)
        
        
        ll_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ll_obj_combo = QComboBox(); self.parent_window.ll_obj_combo.setMinimumWidth(150); self.parent_window.ll_obj_combo.setMinimumHeight(32)
        self.parent_window.ll_refresh_btn = QPushButton(t("refresh")); self.parent_window.ll_refresh_btn.setMinimumHeight(32); self.parent_window.ll_refresh_btn.clicked.connect(self.refresh_objects)
        r0_ll = QHBoxLayout(); r0_ll.addWidget(self.parent_window.ll_obj_combo, 1); r0_ll.addWidget(self.parent_window.ll_refresh_btn)
        ll_grid.addLayout(r0_ll, 0, 1, 1, 3)
        
        ll_grid.addWidget(QLabel("Selection 1:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ll_sel1 = QLineEdit(); self.parent_window.ll_sel1.setPlaceholderText("e.g. resn LIG1 or resi 100"); self.parent_window.ll_sel1.setMinimumHeight(32)
        ll_grid.addWidget(self.parent_window.ll_sel1, 1, 1)
        
        ll_grid.addWidget(QLabel("Selection 2:"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ll_sel2 = QLineEdit(); self.parent_window.ll_sel2.setPlaceholderText("e.g. resn LIG2 or resi 200"); self.parent_window.ll_sel2.setMinimumHeight(32)
        ll_grid.addWidget(self.parent_window.ll_sel2, 1, 3)
        
        ll_grid.addWidget(QLabel("Distance (Å):"), 2, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ll_dist = QLineEdit("4.5"); self.parent_window.ll_dist.setMinimumHeight(32)
        ll_grid.addWidget(self.parent_window.ll_dist, 2, 1)
        
        ll_grid.addWidget(QLabel("Output CSV:"), 2, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ll_csv = QLineEdit(); self.parent_window.ll_csv.setPlaceholderText("Optional"); self.parent_window.ll_csv.setMinimumHeight(32)
        self.parent_window.ll_csv_btn = QPushButton(t("browse")); self.parent_window.ll_csv_btn.setMinimumHeight(32)
        self.parent_window.ll_csv_btn.clicked.connect(lambda: self._browse_save_file(self.parent_window.ll_csv, "CSV (*.csv)"))
        r2_ll = QHBoxLayout(); r2_ll.addWidget(self.parent_window.ll_csv, 1); r2_ll.addWidget(self.parent_window.ll_csv_btn)
        ll_grid.addLayout(r2_ll, 2, 3)

        ll_layout.addLayout(ll_grid)

        ll_btn_row = QHBoxLayout()
        ll_btn_row.setSpacing(10)

        self.parent_window.ll_analyze_btn = QPushButton("Analyze Ligand-Ligand")
        self.parent_window.ll_analyze_btn.setMinimumHeight(36)
        self.parent_window.ll_analyze_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.ll_analyze_btn.clicked.connect(self.run_ll_analysis)
        ll_btn_row.addWidget(self.parent_window.ll_analyze_btn)
        ll_btn_row.addStretch(1)

        layout.addWidget(grp_ll)
        layout.addLayout(ll_btn_row)

        # 4. Electrostatic Complementarity (EC) Analysis
        grp_ec = QFrame()
        grp_ec.setStyleSheet(self._get_card_style(is_dark))
        ec_layout = QVBoxLayout(grp_ec)
        ec_layout.setSpacing(12)
        ec_layout.setContentsMargins(16, 14, 16, 14)

        ec_title = QLabel("Electrostatic Complementarity (EC) Analysis")
        ec_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        ec_layout.addWidget(ec_title)

        ec_grid = QGridLayout()
        ec_grid.setContentsMargins(12, 16, 12, 8)
        ec_grid.setColumnStretch(1, 1); ec_grid.setColumnStretch(3, 1)
        ec_grid.setHorizontalSpacing(12); ec_grid.setVerticalSpacing(4)
        
        
        ec_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_obj_combo = QComboBox(); self.parent_window.ec_obj_combo.setMinimumWidth(150); self.parent_window.ec_obj_combo.setMinimumHeight(32)
        self.parent_window.ec_refresh_btn = QPushButton(t("refresh")); self.parent_window.ec_refresh_btn.setMinimumHeight(32); self.parent_window.ec_refresh_btn.clicked.connect(self.refresh_objects)
        r0_ec = QHBoxLayout(); r0_ec.addWidget(self.parent_window.ec_obj_combo, 1); r0_ec.addWidget(self.parent_window.ec_refresh_btn)
        ec_grid.addLayout(r0_ec, 0, 1)
        
        ec_grid.addWidget(QLabel("Ligand/Glue Name:"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_ligand_name = QLineEdit(); self.parent_window.ec_ligand_name.setPlaceholderText("e.g. LIG, CC885"); self.parent_window.ec_ligand_name.setMinimumHeight(32)
        ec_grid.addWidget(self.parent_window.ec_ligand_name, 0, 3)
        
        ec_grid.addWidget(QLabel("Analysis Mode:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_mode_combo = QComboBox(); self.parent_window.ec_mode_combo.setMinimumHeight(32)
        self.parent_window.ec_mode_combo.addItems(["Protein-Ligand EC", "Ternary Complex EC (Molecular Glue)"])
        self.parent_window.ec_mode_combo.currentIndexChanged.connect(self._on_ec_mode_changed)
        ec_grid.addWidget(self.parent_window.ec_mode_combo, 1, 1)
        
        ec_grid.addWidget(QLabel("pH:"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_ph = QLineEdit("7.4"); self.parent_window.ec_ph.setMinimumHeight(32)
        ec_grid.addWidget(self.parent_window.ec_ph, 1, 3)
        
        # Ternary-specific options (initially hidden)
        ec_grid.addWidget(QLabel("Protein A Chains:"), 2, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_protein_a_chains = QLineEdit(); self.parent_window.ec_protein_a_chains.setPlaceholderText("e.g. A (E3 ligase)"); self.parent_window.ec_protein_a_chains.setMinimumHeight(32)
        ec_grid.addWidget(self.parent_window.ec_protein_a_chains, 2, 1)
        
        ec_grid.addWidget(QLabel("Protein B Chains:"), 2, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_protein_b_chains = QLineEdit(); self.parent_window.ec_protein_b_chains.setPlaceholderText("e.g. B (Substrate)"); self.parent_window.ec_protein_b_chains.setMinimumHeight(32)
        ec_grid.addWidget(self.parent_window.ec_protein_b_chains, 2, 3)
        
        ec_grid.addWidget(QLabel("Surface Density:"), 3, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_surface_density = QLineEdit("10.0"); self.parent_window.ec_surface_density.setPlaceholderText("Points/Ų"); self.parent_window.ec_surface_density.setMinimumHeight(32)
        ec_grid.addWidget(self.parent_window.ec_surface_density, 3, 1)
        
        ec_grid.addWidget(QLabel("Output Directory:"), 3, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_output_dir = QLineEdit(); self.parent_window.ec_output_dir.setPlaceholderText("Optional (temp dir if blank)"); self.parent_window.ec_output_dir.setMinimumHeight(32)
        self.parent_window.ec_output_btn = QPushButton(t("browse")); self.parent_window.ec_output_btn.setMinimumHeight(32)
        self.parent_window.ec_output_btn.clicked.connect(self._browse_ec_output_dir)
        r3_ec = QHBoxLayout(); r3_ec.addWidget(self.parent_window.ec_output_dir, 1); r3_ec.addWidget(self.parent_window.ec_output_btn)
        ec_grid.addLayout(r3_ec, 3, 3)
        
        # Advanced EC Options (σ-hole, Lone Pairs, Bridging Waters)
        ec_grid.addWidget(QLabel("Advanced Options:"), 4, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        
        self.parent_window.ec_use_sigma_holes = QCheckBox("σ-hole (Cl/Br/I)")
        self.parent_window.ec_use_sigma_holes.setChecked(True)
        self.parent_window.ec_use_sigma_holes.setToolTip(
            "Add σ-hole virtual points for halogen atoms (Cl, Br, I).\n"
            "Improves EC accuracy for halogen bond interactions.\n"
            "Recommended for ligands containing halogens."
        )
        ec_grid.addWidget(self.parent_window.ec_use_sigma_holes, 4, 1)
        
        self.parent_window.ec_use_lone_pairs = QCheckBox("Lone Pairs (C=O)")
        self.parent_window.ec_use_lone_pairs.setChecked(True)
        self.parent_window.ec_use_lone_pairs.setToolTip(
            "Add lone pair virtual points for carbonyl oxygen.\n"
            "Improves H-bond directionality prediction.\n"
            "May increase computation time."
        )
        ec_grid.addWidget(self.parent_window.ec_use_lone_pairs, 4, 3)
        
        self.parent_window.ec_visualize = QCheckBox("Visualize EC Map in PyMOL")
        self.parent_window.ec_visualize.setChecked(True)
        ec_grid.addWidget(self.parent_window.ec_visualize, 5, 1)
        
        self.parent_window.ec_keep_waters = QCheckBox("Keep Bridging Waters")
        self.parent_window.ec_keep_waters.setToolTip(
            "Retain structurally important bridging water molecules.\n"
            "Waters that form ≥2 H-bonds with protein or bridge protein-ligand.\n"
            "Useful for high-resolution crystal structures."
        )
        ec_grid.addWidget(self.parent_window.ec_keep_waters, 5, 3)

        ec_layout.addLayout(ec_grid)

        ec_btn_row = QHBoxLayout()
        ec_btn_row.setSpacing(10)

        self.parent_window.ec_analyze_btn = QPushButton("Calculate EC")
        self.parent_window.ec_analyze_btn.setMinimumHeight(36)
        self.parent_window.ec_analyze_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.ec_analyze_btn.clicked.connect(self.run_ec_analysis)

        self.parent_window.ec_hotspots_btn = QPushButton("Find EC Hotspots")
        self.parent_window.ec_hotspots_btn.setMinimumHeight(36)
        self.parent_window.ec_hotspots_btn.setStyleSheet(self._get_secondary_btn_style())
        self.parent_window.ec_hotspots_btn.setToolTip("Identify regions of strong electrostatic complementarity")
        self.parent_window.ec_hotspots_btn.clicked.connect(self.run_ec_hotspots)
        
        ec_btn_row.addWidget(self.parent_window.ec_analyze_btn)
        ec_btn_row.addWidget(self.parent_window.ec_hotspots_btn)
        ec_btn_row.addStretch(1)
        
        layout.addWidget(grp_ec)
        layout.addLayout(ec_btn_row)
        
        # Initially hide ternary-specific fields
        self._on_ec_mode_changed(0)
        
        
        layout.addStretch(1)

    # --- PPI Analysis Logic ---
    def run_ppi_analysis(self):
        """Analyze protein-protein interface (PPI Interface)"""
        try:
            obj = self.parent_window.ppi_obj_combo.currentText().strip()
            if not obj or obj == t("no_object"):
                show_message_box(self, "Warning", "Please select a structure object", "warning")
                return

            p1_chains = self.parent_window.ppi_protein1_chains.text().strip()
            p2_chains = self.parent_window.ppi_protein2_chains.text().strip()

            if not p1_chains or not p2_chains:
                show_message_box(self, "Warning", "Please specify both protein chain groups", "warning")
                return
            
            self.log(f"Starting PPI analysis for {obj}...")
            
            try: from ...ppi_analyzer import analyze_protein_protein_interface
            except ImportError: from ppi_analyzer import analyze_protein_protein_interface
            
            p1_list = [c.strip() for c in p1_chains.split(",")]
            p2_list = [c.strip() for c in p2_chains.split(",")]
            interface_dist = float(self.parent_window.ppi_interface_dist.text())
            output_csv = self.parent_window.ppi_csv.text().strip() or None
            
            # Get Visualization Options
            display_mode_idx = self.parent_window.ppi_display_mode.currentIndex()
            display_mode = "cartoon_surface_interaction" if display_mode_idx == 0 else "cartoon_interaction"
            show_labels = self.parent_window.ppi_show_labels.isChecked()

            result = analyze_protein_protein_interface(
                obj_name=obj,
                protein1_chains=p1_list,
                protein2_chains=p2_list,
                interface_distance=interface_dist,
                output_csv=output_csv,
                visualize=True,
                display_mode=display_mode,
                show_labels=show_labels
            )
            
            if result:
                contacts = result.get('interface_contacts', 0)
                bsa = result.get('bsa')
                strength = result.get('interface_strength', 0)
                is_strong = result.get('is_strong_interface', False)
                
                self.log(f"PPI Analysis Complete:")
                self.log(f"  Interface Contacts: {contacts}")
                if bsa:
                    self.log(f"  BSA: {bsa:.1f} Ų")
                self.log(f"  Interface Strength: {strength:.1f}/10")
                self.log(f"  Classification: {'Strong Interface' if is_strong else 'Weak Interface'}")

                show_message_box(self, "PPI Analysis Complete",
                    f"Interface Contacts: {contacts}\n"
                    f"{'BSA: ' + str(round(bsa, 1)) + ' Ų' if bsa else 'BSA: N/A'}\n"
                    f"Interface Strength: {strength:.1f}/10\n\n"
                    f"{'Strong Interface' if is_strong else 'Weak Interface'}"
                )
            else:
                self.log("PPI analysis failed")
        except Exception as e:
            self.on_error(str(e))
            import traceback; traceback.print_exc()
    


    # --- PL Logic ---
    def run_pl_analysis(self):
        """Analyze protein-ligand interactions"""
        try:
            obj_name = self.parent_window.pl_obj_combo.currentText()
            if obj_name == t("no_object"):
                show_message_box(self, "Warning", "Please select a structure object", "warning")
                return

            ligand_name = self.parent_window.pl_ligand_name.text().strip() or None
            protein_chains_str = self.parent_window.pl_protein_chains.text().strip()
            protein_chains = [c.strip() for c in protein_chains_str.split(",")] if protein_chains_str else None
            distance = float(self.parent_window.pl_distance.text())
            output_csv = self.parent_window.pl_csv.text().strip() or None

            # Get 3D Visualization Options
            show_hydrophobic = self.parent_window.pl_show_hydrophobic.isChecked()
            min_conf_str = self.parent_window.pl_min_confidence.currentText().split()[0]
            try:
                min_confidence = float(min_conf_str)
            except:
                min_confidence = 0.0

            self.log(f"Starting Protein-Ligand analysis for {obj_name}...")
            if ligand_name:
                self.log(f"  Ligand: {ligand_name}")
            else:
                self.log(f"  Ligand: Auto-detect")

            try: from ...interaction_analyzer import analyze_protein_ligand_interactions
            except ImportError: from interaction_analyzer import analyze_protein_ligand_interactions

            result = analyze_protein_ligand_interactions(
                obj_name=obj_name,
                ligand_resname=ligand_name,
                protein_chains=protein_chains,
                distance_cutoff=distance,
                output_csv=output_csv
            )

            if result:
                n = len(result.get("interactions", []))
                ligand_used = result.get("ligand_resname", ligand_name or "Unknown")
                self.log(f"Protein-Ligand Analysis Complete:")
                self.log(f"  Ligand: {ligand_used}")
                self.log(f"  Interactions found: {n}")
                self.parent_window.current_pl_result = result

                # Add 3D visualization
                if n > 0:
                    try:
                        try: from ...interaction_analyzer import visualize_protein_ligand_3d
                        except ImportError: from interaction_analyzer import visualize_protein_ligand_3d

                        # Get actual ligand name
                        actual_ligand = ligand_name
                        if result.get("ligand_residues"):
                            actual_ligand = result["ligand_residues"][0].get("resname", ligand_name)

                        visualize_protein_ligand_3d(
                            obj_name,
                            result,
                            actual_ligand,
                            show_hydrophobic=show_hydrophobic,
                            min_confidence=min_confidence
                        )
                        self.log(f"  3D visualization generated (Conf>={min_confidence}, Hydrophobic={show_hydrophobic})")
                    except Exception as viz_e:
                        self.log(f"  3D visualization failed: {viz_e}")
                        import traceback; traceback.print_exc()

                show_message_box(self, "Analysis Complete",
                    f"Ligand: {ligand_used}\nInteractions found: {n}")
            else:
                self.log("No interactions found")
        except Exception as e:
            self.on_error(str(e))
            import traceback; traceback.print_exc()


    def run_pl_2d_diagram(self):
        """Generate 2D Interaction Diagram"""
        try:
            # Check for RDKit
            try:
                import rdkit
            except ImportError:
                show_message_box(self, "Error", "RDKit is required for 2D diagrams.\nPlease install it: pip install rdkit", "critical")
                return

            # Get parameters
            obj_name = self.parent_window.pl_obj_combo.currentText()
            ligand_name = self.parent_window.pl_ligand_name.text().strip()
            
            # Try to get data from current result or CSV
            csv_path = self.parent_window.pl_csv.text().strip()
            
            # If no CSV path provided, try to use a temp file from current result
            if not csv_path and hasattr(self.parent_window, "current_pl_result") and self.parent_window.current_pl_result:
                # We need to save the current result to a temp CSV if not saved yet
                # For simplicity, let's ask user to analyze first/provide CSV if they haven't
                pass

            if not csv_path or not os.path.exists(csv_path):
                # Try to use the last analysis result if available
                # But generate_2d_interaction_diagram requires a CSV path currently

                # If we have a result object, maybe we can save it to temp
                if hasattr(self.parent_window, "current_pl_result") and self.parent_window.current_pl_result:
                    import tempfile, csv
                    result = self.parent_window.current_pl_result

                    # Check if result has interactions
                    interactions = result.get("interactions", [])
                    if not interactions:
                        show_message_box(self, "Warning", "No interactions to plot.", "warning")
                        return

                    # Infer ligand name from result if not provided
                    if not ligand_name and result.get("ligand_residues"):
                         ligand_name = result["ligand_residues"][0].get("resname", "LIG")

                    # Create temp CSV
                    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".csv", mode='w', encoding='utf-8-sig')
                    # Use keys from first item, but ignore extras in others
                    if interactions:
                        fieldnames = list(interactions[0].keys())
                        writer = csv.DictWriter(tmp, fieldnames=fieldnames, extrasaction='ignore')
                        writer.writeheader()
                        writer.writerows(interactions)
                    tmp.close()
                    csv_path = tmp.name
                    self.log(f"Using temporary CSV: {csv_path}")
                else:
                    show_message_box(self, "Missing Data", "Please run analysis first (with output CSV) or select an existing CSV file.", "warning")
                    return

            if not ligand_name:
                show_message_box(self, "Missing Input", "Please specify the Ligand Name (Residue Name).", "warning")
                return

            # Use the same visualization options as 3D view
            show_hydrophobic = self.parent_window.pl_show_hydrophobic.isChecked()
            min_conf_str = self.parent_window.pl_min_confidence.currentText().split()[0]
            try:
                min_confidence = float(min_conf_str)
            except Exception:
                min_confidence = 0.8

            # Output path
            default_name = f"{ligand_name}_2d.png"
            out_path, _ = QFileDialog.getSaveFileName(self, "Save 2D Diagram", default_name, "PNG Image (*.png)")
            if not out_path: return

            self.log(f"Generating 2D diagram for {ligand_name}...")

            try: from ...interaction_2d_plot import generate_2d_interaction_diagram
            except ImportError: from interaction_2d_plot import generate_2d_interaction_diagram

            final_path = generate_2d_interaction_diagram(
                csv_path=csv_path,
                ligand_resname=ligand_name,
                obj_name=obj_name,
                output_path=out_path,
                min_confidence=min_confidence
            )

            if final_path and os.path.exists(final_path):
                self.log(f"2D Diagram saved: {final_path}")

                show_message_box(self, "Success",
                    f"2D Diagram saved successfully!\n\nFile location:\n{final_path}")

                 # Try to open the file (Mac/Linux/Windows)
                try:
                    import subprocess, platform
                    if platform.system() == 'Darwin':       # macOS
                        subprocess.call(('open', final_path))
                    elif platform.system() == 'Windows':    # Windows
                        os.startfile(final_path)
                    else:                                   # linux variants
                        subprocess.call(('xdg-open', final_path))
                except:
                    pass
            else:
                self.log("Failed to generate 2D diagram.")
                show_message_box(self, "Error", "Failed to generate diagram. See log for details.", "warning")
    #
        except Exception as e:
            self.on_error(f"2D Diagram Error: {str(e)}")
            import traceback; traceback.print_exc()

    # --- Browse file helpers ---
    def _browse_save_file(self, line_edit, file_filter):
        """Browse and select save file path"""
        fn, _ = QFileDialog.getSaveFileName(self, "Save File", "", file_filter)
        if fn:
            line_edit.setText(fn)
    
    def _browse_file(self, line_edit, file_filter):
        """Browse and select file to open"""
        fn, _ = QFileDialog.getOpenFileName(self, "Open File", "", file_filter)
        if fn:
            line_edit.setText(fn)
    
    # --- LL (Ligand-Ligand) Logic ---
    def run_ll_analysis(self):
        """Run ligand-ligand interaction analysis"""
        try:
            # Try to import analysis function
            try: from ...ligand_ligand_analyzer import analyze_ligand_ligand_interactions
            except ImportError:
                try: from ligand_ligand_analyzer import analyze_ligand_ligand_interactions
                except ImportError:
                    show_message_box(self, "Error", "Ligand-Ligand analysis module not found.", "critical")
                    return
            
            obj = self.parent_window.ll_obj_combo.currentText()
            sel1 = self.parent_window.ll_sel1.text().strip()
            sel2 = self.parent_window.ll_sel2.text().strip()
            dist_str = self.parent_window.ll_dist.text().strip()
            csv_path = self.parent_window.ll_csv.text().strip() or None
            
            if not obj or obj == t("no_object"):
                show_message_box(self, "Missing Input", "Please select a target object.", "warning")
                return
            if not sel1 or not sel2:
                show_message_box(self, "Missing Input", "Please define both Selection 1 and Selection 2.", "warning")
                return

            try:
                dist = float(dist_str) if dist_str else 4.5
            except ValueError:
                show_message_box(self, "Invalid Input", "Distance must be a number.", "warning")
                return

            self.log(f"Starting Ligand-Ligand analysis for {obj}...")
            self.log(f"  Selection 1: {sel1}")
            self.log(f"  Selection 2: {sel2}")
            self.log(f"  Distance cutoff: {dist} Å")

            try:
                import pymol

                # Build complete PyMOL selection
                full_sel1 = f"({obj}) and ({sel1})"
                full_sel2 = f"({obj}) and ({sel2})"

                interactions = analyze_ligand_ligand_interactions(
                    obj_name=obj,
                    sel1=full_sel1,
                    sel2=full_sel2,
                    cutoff=dist,
                    output_csv=csv_path,
                    visualize=True
                )

                self.log(f"Ligand-Ligand Analysis Complete:")
                self.log(f"  Interactions found: {len(interactions)}")

                msg = f"Analysis complete.\nFound {len(interactions)} interactions."
                if csv_path:
                    msg += f"\n\nSaved to: {csv_path}"

                show_message_box(self, "Success", msg)

            except pymol.CmdException as e:
                msg = str(e)
                if "Invalid selection" in msg:
                    show_message_box(self, "Selection Error",
                        f"PyMOL could not understand your selection.\n\n"
                        f"Error: {msg}\n\n"
                        f"Tip: Please use valid PyMOL selection syntax.\n"
                        f"Examples:\n"
                        f"\u2022 resn LIG (by residue name)\n"
                        f"\u2022 resi 900 (by residue index)\n"
                        f"\u2022 chain A (by chain)\n\n"
                        f"You entered: '{sel1}' and '{sel2}'", "critical")
                else:
                    show_message_box(self, "PyMOL Error", str(e), "critical")
                return
                
        except Exception as e:
            self.on_error(str(e))
            import traceback; traceback.print_exc()
    
    # --- EC Analysis Logic ---
    def _on_ec_mode_changed(self, index):
        """Show/hide ternary-specific fields based on mode selection"""
        is_ternary = (index == 1)
        self.parent_window.ec_protein_a_chains.setEnabled(is_ternary)
        self.parent_window.ec_protein_b_chains.setEnabled(is_ternary)
        
        if is_ternary:
            self.parent_window.ec_protein_a_chains.setPlaceholderText("e.g. A (E3 ligase) - Required")
            self.parent_window.ec_protein_b_chains.setPlaceholderText("e.g. B (Substrate) - Required")
        else:
            self.parent_window.ec_protein_a_chains.setPlaceholderText("Not used in this mode")
            self.parent_window.ec_protein_b_chains.setPlaceholderText("Not used in this mode")
    
    def _browse_ec_output_dir(self):
        """Browse for EC output directory"""
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.parent_window.ec_output_dir.setText(dir_path)
    
    def run_ec_analysis(self):
        """Run Electrostatic Complementarity analysis"""
        try:
            obj_name = self.parent_window.ec_obj_combo.currentText()
            if obj_name == t("no_object") or not obj_name:
                show_message_box(self, "Warning", "Please select a structure object", "warning")
                return

            ligand_name = self.parent_window.ec_ligand_name.text().strip()
            if not ligand_name:
                show_message_box(self, "Warning", "Please specify the ligand/glue residue name", "warning")
                return

            mode = self.parent_window.ec_mode_combo.currentIndex()
            ph = float(self.parent_window.ec_ph.text().strip() or "7.4")
            surface_density = float(self.parent_window.ec_surface_density.text().strip() or "10.0")
            output_dir = self.parent_window.ec_output_dir.text().strip() or None
            visualize = self.parent_window.ec_visualize.isChecked()

            # Get advanced options
            use_sigma_holes = self.parent_window.ec_use_sigma_holes.isChecked()
            use_lone_pairs = self.parent_window.ec_use_lone_pairs.isChecked()
            keep_bridging_waters = self.parent_window.ec_keep_waters.isChecked()

            self.log(f"Starting EC analysis for {obj_name}...")
            self.log(f"  Ligand/Glue: {ligand_name}")
            self.log(f"  Mode: {'Ternary Complex' if mode == 1 else 'Protein-Ligand'}")
            self.log(f"  pH: {ph}")

            # Log advanced options
            if use_sigma_holes:
                self.log(f"  ✨ σ-hole virtual points: ENABLED")
            if use_lone_pairs:
                self.log(f"  ✨ Lone pair virtual points: ENABLED")
            if keep_bridging_waters:
                self.log(f"  ✨ Bridging waters: ENABLED")

            # Import EC calculator
            try:
                from ...ligand_ec_calculator import calculate_ligand_ec, analyze_ternary_ec
            except ImportError:
                try:
                    from ligand_ec_calculator import calculate_ligand_ec, analyze_ternary_ec
                except ImportError:
                    show_message_box(self, "Error",
                        "EC Calculator module not found.\n\n"
                        "Please ensure ligand_ec_calculator.py is installed.", "critical")
                    return

            if mode == 0:
                # Protein-Ligand EC
                result = calculate_ligand_ec(
                    obj_name=obj_name,
                    ligand_resname=ligand_name,
                    output_dir=output_dir,
                    ph=ph,
                    surface_density=surface_density,
                    visualize=visualize,
                    use_sigma_holes=use_sigma_holes,
                    use_lone_pairs=use_lone_pairs,
                    keep_bridging_waters=keep_bridging_waters
                )

                if result:
                    ec_score = result.get('ec_score', 0)
                    ec_stats = result.get('ec_statistics', {})

                    self.log(f"EC Analysis Complete:")
                    self.log(f"  EC Score: {ec_score:.4f}")
                    self.log(f"  EC Mean: {ec_stats.get('ec_mean', 0):.4f}")
                    self.log(f"  Positive EC fraction: {ec_stats.get('ec_positive_fraction', 0)*100:.1f}%")

                    # Interpretation
                    if ec_score > 0.3:
                        interpretation = "Strong electrostatic complementarity - favorable binding"
                    elif ec_score > 0:
                        interpretation = "Moderate electrostatic complementarity"
                    else:
                        interpretation = "Poor electrostatic complementarity - potential clash"

                    self.log(f"  Interpretation: {interpretation}")

                    show_message_box(self, "EC Analysis Complete",
                        f"EC Score: {ec_score:.4f}\n"
                        f"EC Mean: {ec_stats.get('ec_mean', 0):.4f}\n"
                        f"Positive EC: {ec_stats.get('ec_positive_fraction', 0)*100:.1f}%\n\n"
                        f"{interpretation}\n\n"
                        f"Output: {result.get('output_dir', 'N/A')}")
                else:
                    self.log("EC analysis failed")
                    show_message_box(self, "Error", "EC analysis failed. Check the log for details.", "warning")

            else:
                # Ternary Complex EC (Molecular Glue)
                protein_a_chains = self.parent_window.ec_protein_a_chains.text().strip()
                protein_b_chains = self.parent_window.ec_protein_b_chains.text().strip()

                if not protein_a_chains or not protein_b_chains:
                    show_message_box(self, "Warning",
                        "For ternary complex analysis, please specify both Protein A and Protein B chains", "warning")
                    return

                protein_a_list = [c.strip() for c in protein_a_chains.split(",")]
                protein_b_list = [c.strip() for c in protein_b_chains.split(",")]

                self.log(f"  Protein A chains: {protein_a_list}")
                self.log(f"  Protein B chains: {protein_b_list}")
                
                result = analyze_ternary_ec(
                    obj_name=obj_name,
                    glue_resname=ligand_name,
                    protein_a_chains=protein_a_list,
                    protein_b_chains=protein_b_list,
                    output_dir=output_dir,
                    ph=ph,
                    surface_density=surface_density,
                    visualize=visualize
                )
                
                if result and 'combined' in result:
                    combined = result['combined']
                    ec_a = combined.get('ec_a_glue', 0)
                    ec_b = combined.get('ec_b_glue', 0)
                    ec_combined = combined.get('ec_combined_score', 0)
                    
                    # Get overlap region data (most important for molecular glue!)
                    overlap_data = combined.get('ec_overlap', {})
                    pos_frac_a = overlap_data.get('pos_frac_a', 0)
                    pos_frac_b = overlap_data.get('pos_frac_b', 0)
                    pos_frac_combined = overlap_data.get('pos_frac_combined', 0)
                    n_overlap = overlap_data.get('n_points', 0)
                    
                    self.log(f"Ternary EC Analysis Complete:")
                    self.log(f"  Bridging Zone ({n_overlap} points):")
                    self.log(f"    Positive EC (A-Glue): {pos_frac_a*100:.1f}%")
                    self.log(f"    Positive EC (B-Glue): {pos_frac_b*100:.1f}%")
                    self.log(f"    Combined: {pos_frac_combined*100:.1f}%")
                    
                    # Interpretation based on positive EC fraction
                    if pos_frac_combined > 0.6:
                        interpretation = "Good complementarity (>60% positive)"
                    elif pos_frac_combined > 0.5:
                        interpretation = "Moderate complementarity (50-60% positive)"
                    else:
                        interpretation = "Poor complementarity (<50% positive)"
                    
                    self.log(f"  Interpretation: {interpretation}")

                    show_message_box(self, "Ternary EC Analysis Complete",
                        f"Bridging Zone Analysis ({n_overlap} points):\n\n"
                        f"Positive EC (A-Glue): {pos_frac_a*100:.1f}%\n"
                        f"Positive EC (B-Glue): {pos_frac_b*100:.1f}%\n"
                        f"Combined: {pos_frac_combined*100:.1f}%\n\n"
                        f"{interpretation}\n\n"
                        f"Output: {result.get('output_dir', 'N/A')}")
                else:
                    self.log("Ternary EC analysis failed or incomplete")
                    show_message_box(self, "Error", "Ternary EC analysis failed. Check the log for details.", "warning")
                    
        except Exception as e:
            self.on_error(str(e))
            import traceback; traceback.print_exc()
    
    def run_ec_hotspots(self):
        """Run EC hotspot analysis to identify regions of strong complementarity"""
        try:
            obj_name = self.parent_window.ec_obj_combo.currentText()
            if obj_name == t("no_object") or not obj_name:
                show_message_box(self, "Warning", "Please select a structure object", "warning")
                return

            ligand_name = self.parent_window.ec_ligand_name.text().strip()
            if not ligand_name:
                show_message_box(self, "Warning", "Please specify the ligand residue name", "warning")
                return

            ph = float(self.parent_window.ec_ph.text().strip() or "7.4")
            output_dir = self.parent_window.ec_output_dir.text().strip() or None

            # Get advanced options
            use_sigma_holes = self.parent_window.ec_use_sigma_holes.isChecked()
            use_lone_pairs = self.parent_window.ec_use_lone_pairs.isChecked()

            self.log(f"Starting EC hotspot analysis for {obj_name}...")
            self.log(f"  Ligand: {ligand_name}")

            # Import EC calculator
            try:
                from ...ligand_ec_calculator import calculate_ec_hotspots
            except ImportError:
                try:
                    from ligand_ec_calculator import calculate_ec_hotspots
                except ImportError:
                    show_message_box(self, "Error",
                        "EC Calculator module not found.\n\n"
                        "Please ensure ligand_ec_calculator.py is installed.", "critical")
                    return

            result = calculate_ec_hotspots(
                obj_name=obj_name,
                ligand_resname=ligand_name,
                output_dir=output_dir,
                ph=ph,
                surface_density=15.0,  # Higher density for hotspot detection
                hotspot_threshold=0.5
            )

            if result and 'hotspots' in result:
                pos_hotspots = result['hotspots']['positive']
                neg_hotspots = result['hotspots']['negative']

                self.log(f"EC Hotspot Analysis Complete:")
                self.log(f"  Positive hotspots (complementary): {pos_hotspots['count']} points ({pos_hotspots['fraction']*100:.1f}%)")
                self.log(f"  Negative hotspots (clash): {neg_hotspots['count']} points ({neg_hotspots['fraction']*100:.1f}%)")

                if 'n_clusters' in pos_hotspots:
                    self.log(f"  Distinct complementary regions: {pos_hotspots['n_clusters']}")

                show_message_box(self, "EC Hotspot Analysis Complete",
                    f"Positive Hotspots (Complementary):\n"
                    f"  Count: {pos_hotspots['count']} points\n"
                    f"  Fraction: {pos_hotspots['fraction']*100:.1f}%\n"
                    f"  Mean EC: {pos_hotspots['mean_ec']:.4f}\n\n"
                    f"Negative Hotspots (Clash):\n"
                    f"  Count: {neg_hotspots['count']} points\n"
                    f"  Fraction: {neg_hotspots['fraction']*100:.1f}%\n"
                    f"  Mean EC: {neg_hotspots['mean_ec']:.4f}\n\n"
                    f"Output: {result.get('output_dir', 'N/A')}")
            else:
                self.log("EC hotspot analysis failed")
                show_message_box(self, "Error", "EC hotspot analysis failed. Check the log for details.", "warning")
                
        except Exception as e:
            self.on_error(str(e))
            import traceback; traceback.print_exc()
    
    # --- Mutation Logic ---
    def run_mutation(self):
        obj_name = self.parent_window.mut_obj_combo.currentText()
        mutations_str = self.parent_window.mut_input.text().strip()
        if not obj_name or not mutations_str: return
        
        try:
            mutations = []
            for mut in mutations_str.split(','):
                parts = mut.strip().split(':')
                if len(parts) == 3: mutations.append((parts[0], parts[1], parts[2]))
            
            if not mutations:
                self.log("Invalid mutation format (e.g. A:23:ALA)")
                return
                
            from pymol import cmd
            # Assuming perform_mutation is monkey-patched or available
            # If not, we need to implement it or check if it's a custom method
            if hasattr(cmd, 'perform_mutation'):
                cmd.perform_mutation(obj_name, mutations, method='pymol')
                self.log("Mutation performed")
            else:
                # Fallback basic mutation using PyMOL wizard or simple command
                for chain, resi, resn in mutations:
                    sel = f"/{obj_name}//{chain}/{resi}"
                    cmd.wizard("mutagenesis")
                    cmd.get_wizard().do_select(sel)
                    cmd.get_wizard().set_mode(resn)
                    cmd.get_wizard().apply()
                    cmd.set_wizard() # close
                self.log("Mutation applied via PyMOL Wizard")
        except Exception as e:
            self.on_error(str(e))

    def run_minimize(self):
        obj = self.parent_window.mut_obj_combo.currentText()
        if not obj: return
        try:
            from pymol import cmd
            # minimize_energy is likely a custom function or part of a plugin
            # Standard PyMOL has 'minimize' or 'clean'
            # Checking unified_gui original code, it calls cmd.minimize_energy
            if hasattr(cmd, 'minimize_energy'):
                cmd.minimize_energy(obj)
            else:
                # Fallback
                cmd.protect(f"not {obj}")
                cmd.sculpt_activate(obj)
                cmd.sculpt_iterate(obj, cycles=100)
                self.log(" minimization done (sculpt)")
        except Exception as e:
            self.on_error(str(e))

    def run_mutation_analysis(self):
        """Run full mutation ΔΔG analysis using FoldX"""
        obj_name = self.parent_window.mut_obj_combo.currentText()
        if not obj_name or obj_name == t("no_object"):
            show_message_box(self, "Warning", "Please select a structure object", "warning")
            return

        # Check if FoldX is available
        try:
            from ...mutation_analyzer import _detect_foldx
        except ImportError:
            try:
                from mutation_analyzer import _detect_foldx
            except ImportError:
                self.log("❌ mutation_analyzer module not found")
                return

        foldx_path = _detect_foldx()
        if not foldx_path:
            show_message_box(self, "FoldX Not Found",
                "FoldX is required for ΔΔG analysis.\n\n"
                "Installation:\n"
                "1. Download FoldX: https://foldxsuite.crg.eu/\n"
                "2. Set environment variable: export FOLDX=/path/to/foldx\n"
                "   Or add FoldX to your PATH", "warning")
            self.log("❌ FoldX not found. Please install FoldX for ΔΔG analysis.")
            return
        
        self.log(f"✅ FoldX detected: {foldx_path}")
        self.log("💡 For ΔΔG heatmap analysis, use the PyMOL command:")
        self.log("   ddg_heatmap('CRBN_selection', 'POI_selection')")
        self.log("")
        self.log("FoldX Input:")
        self.log("  - PDB file: CRBN + POI complex structure")
        self.log("  - Mutation file: <chain><resi><icode><WT><Mut> format")
        self.log("")
        self.log("FoldX Output:")
        self.log("  - DifferencesBetweenMutantAndWildType_fxout.csv")
        self.log("  - Contains ΔΔG values (kcal/mol)")

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

    def _get_primary_btn_style(self) -> str:
        """主按钮样式 - 蓝色渐变"""
        return """
            QPushButton {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #3b82f6, stop:1 #2563eb);
                color: white; border: none; border-radius: 8px;
                font-weight: 600; font-size: 13px; padding: 8px 20px;
            }
            QPushButton:hover { background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #2563eb, stop:1 #1d4ed8); }
            QPushButton:pressed { background: #1d4ed8; }
            QPushButton:disabled { background: #94a3b8; }
        """

    def _get_secondary_btn_style(self) -> str:
        """次要按钮样式 - 浅灰色"""
        return """
            QPushButton {
                background: #f1f5f9; color: #475569; border: 1px solid #e2e8f0;
                border-radius: 8px; padding: 8px 16px; font-weight: 500;
            }
            QPushButton:hover { background: #e2e8f0; }
        """
