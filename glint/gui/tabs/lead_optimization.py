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
    QFileDialog, QMessageBox, QThread, Signal as pyqtSignal
)

from ..utils import t, show_message_box, show_question_box
from .common import CommonTab


class _Diagram2DWorker(QThread):
    """Background worker that generates a 2D interaction diagram safely."""
    finished = pyqtSignal(str)
    error = pyqtSignal(str)

    TOTAL_TIMEOUT = 90

    def __init__(self, csv_path: str, ligand_resname: str,
                 pdb_file: str, output_path: str, min_confidence: float):
        super().__init__()
        self.csv_path = csv_path
        self.ligand_resname = ligand_resname
        self.pdb_file = pdb_file
        self.output_path = output_path
        self.min_confidence = min_confidence
        self._timed_out = False

    def run(self) -> None:
        """Run diagram generation in a worker thread with timeout protection."""
        import threading

        result_holder = {'path': None, 'error': None}

        def _generate():
            try:
                import matplotlib
                matplotlib.use('Agg')

                try:
                    from ...interaction_2d_plot import generate_2d_interaction_diagram
                except ImportError:
                    from interaction_2d_plot import generate_2d_interaction_diagram

                final_path = generate_2d_interaction_diagram(
                    csv_path=self.csv_path,
                    ligand_resname=self.ligand_resname,
                    pdb_file=self.pdb_file,
                    output_path=self.output_path,
                    min_confidence=self.min_confidence,
                )
                result_holder['path'] = final_path
            except Exception as e:
                result_holder['error'] = str(e)

        gen_thread = threading.Thread(target=_generate, daemon=True)
        gen_thread.start()
        gen_thread.join(timeout=self.TOTAL_TIMEOUT)

        if gen_thread.is_alive():
            self._timed_out = True
            self.error.emit(
                f"Diagram generation timed out after {self.TOTAL_TIMEOUT}s.\n"
                "Try a simpler ligand structure or verify that Open Babel is available."
            )
        elif result_holder['error']:
            self.error.emit(result_holder['error'])
        elif result_holder['path'] and os.path.exists(result_holder['path']):
            self.finished.emit(result_holder['path'])
        else:
            self.error.emit("Failed to generate diagram. See log for details.")

        try:
            if self.pdb_file and os.path.exists(self.pdb_file):
                os.remove(self.pdb_file)
        except OSError:
            pass


class LeadOptimizationTab(CommonTab):
    def __init__(self, parent):
        super().__init__(parent)

        # Remove proxies for methods implemented here
        for attr in ['run_mutation', 'run_minimize', 'run_mutation_analysis']:
            if attr in self.__dict__:
                del self.__dict__[attr]

        self.init_ui()

    def init_ui(self):
        """Initialize the tab with a modern card layout."""
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

        # 1. Electrostatic Complementarity (EC) Analysis
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
        self.parent_window.ec_ligand_name = QLineEdit()
        self.parent_window.ec_ligand_name.setPlaceholderText("e.g. LIG, CC885")
        self.parent_window.ec_ligand_name.setMinimumHeight(32)
        ec_grid.addWidget(self.parent_window.ec_ligand_name, 0, 3)
# EC surface style dropdown: Solid Surface calls visualize_ec_smooth_surface, Mesh Surface calls visualize_ec_mesh_surface
# Row 1: Surface Style and Execute Button
        ec_grid.addWidget(QLabel("Surface Style:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_surface_style = QComboBox()
        self.parent_window.ec_surface_style.setMinimumHeight(32)
        self.parent_window.ec_surface_style.addItems(["Solid Surface", "Mesh Surface"])
        self.parent_window.ec_surface_style.setToolTip("Solid: Render EC as solid surface | Mesh: Render EC as semi-transparent mesh surface")
        ec_grid.addWidget(self.parent_window.ec_surface_style, 1, 1)
        
        self.parent_window.ec_analyze_btn = QPushButton("Analyze EC")
        self.parent_window.ec_analyze_btn.setMinimumHeight(36)
        self.parent_window.ec_analyze_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.ec_analyze_btn.clicked.connect(self.run_ec_analysis)
        ec_grid.addWidget(self.parent_window.ec_analyze_btn, 1, 2, 1, 2)

        ec_layout.addLayout(ec_grid)
        layout.addWidget(grp_ec)

        # 2. Protein-Protein Interface (PPI) Analysis
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
        self.parent_window.ppi_display_mode.addItems(["Surface + Interaction", "Cartoon + Interaction"])
        ppi_grid.addWidget(self.parent_window.ppi_display_mode, 3, 1)
        
        self.parent_window.ppi_show_residue_labels = QCheckBox("Show Residue Labels")
        self.parent_window.ppi_show_residue_labels.setChecked(True)
        ppi_grid.addWidget(self.parent_window.ppi_show_residue_labels, 3, 2)
        
        self.parent_window.ppi_show_distance_labels = QCheckBox("Show Distance Labels")
        self.parent_window.ppi_show_distance_labels.setChecked(False)
        ppi_grid.addWidget(self.parent_window.ppi_show_distance_labels, 3, 3)


        self.parent_window.ppi_show_hydrophobic = QCheckBox("Show Hydrophobic Interactions")
        self.parent_window.ppi_show_hydrophobic.setChecked(False)
        ppi_grid.addWidget(self.parent_window.ppi_show_hydrophobic, 4, 1)
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

        # 3. Protein-Ligand Interactions
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
        
        # Distance labels option
        pl_grid.addWidget(QLabel("3D Display:"), 3, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_show_distance_labels = QCheckBox("Show Distance Labels")
        self.parent_window.pl_show_distance_labels.setChecked(False)
        self.parent_window.pl_show_distance_labels.setToolTip("Display distance values on interaction lines in 3D view")
        pl_grid.addWidget(self.parent_window.pl_show_distance_labels, 3, 1)
        
        # Output CSV row
        pl_grid.addWidget(QLabel("Output CSV:"), 4, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        
        self.parent_window.pl_csv = QLineEdit(); self.parent_window.pl_csv.setPlaceholderText("Optional"); self.parent_window.pl_csv.setMinimumHeight(32)
        self.parent_window.pl_csv_btn = QPushButton(t("browse")); self.parent_window.pl_csv_btn.setMinimumHeight(32)
        self.parent_window.pl_csv_btn.clicked.connect(lambda: self._browse_save_file(self.parent_window.pl_csv, "CSV (*.csv)"))
        
        r2_pl = QHBoxLayout(); r2_pl.addWidget(self.parent_window.pl_csv, 1); r2_pl.addWidget(self.parent_window.pl_csv_btn)
        pl_grid.addLayout(r2_pl, 4, 1, 1, 3)

        pl_layout.addLayout(pl_grid)

        # 4. Ligand-Ligand Interactions
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

        layout.addStretch(1)

    def refresh_objects(self):
        """刷新 PyMOL 对象列表到本 Tab 的所有 combo box"""
        try:
            from pymol import cmd
            objects = cmd.get_object_list()
            for combo_attr in ("ppi_obj_combo", "pl_obj_combo", "ll_obj_combo", "ec_obj_combo"):
                combo = getattr(self.parent_window, combo_attr, None)
                if combo is None:
                    continue
                current = combo.currentText()
                combo.clear()
                combo.addItems(objects)
                if current in objects:
                    combo.setCurrentText(current)
            print(f"[GLINT] ✓ Lead Optimization: refreshed objects ({len(objects)} found)")
        except Exception as e:
            print(f"[GLINT] ❌ Failed to refresh objects: {e}")


    # ------------------------------------------------------------------ #
    #  业务方法：PPI / PL / LL / 2D Diagram
    # ------------------------------------------------------------------ #

    def run_ppi_analysis(self):
        """运行蛋白-蛋白界面（PPI）分析"""
        try:
            from pymol import cmd
        except ImportError:
            show_message_box(self, "Error", "PyMOL not available", "critical")
            return

        obj_name = self.parent_window.ppi_obj_combo.currentText().strip()
        if not obj_name:
            show_message_box(self, "Warning", "Please select a target object.", "warning")
            return

        p1 = self.parent_window.ppi_protein1_chains.text().strip()
        p2 = self.parent_window.ppi_protein2_chains.text().strip()
        if not p1 or not p2:
            show_message_box(self, "Warning", "Please specify both Protein1 and Protein2 chain IDs.", "warning")
            return

        protein1_chains = [c.strip() for c in p1.split(",") if c.strip()]
        protein2_chains = [c.strip() for c in p2.split(",") if c.strip()]

        try:
            interface_dist = float(self.parent_window.ppi_interface_dist.text())
        except ValueError:
            show_message_box(self, "Warning", "Invalid interface distance value.", "warning")
            return

        output_csv = self.parent_window.ppi_csv.text().strip() or None

        # 显示模式映射
        mode_map = {"Surface + Interaction": "surface_interaction", "Cartoon + Interaction": "cartoon_interaction"}
        display_mode = mode_map.get(self.parent_window.ppi_display_mode.currentText(), "surface_interaction")

        show_labels = self.parent_window.ppi_show_residue_labels.isChecked()
        show_dist_labels = self.parent_window.ppi_show_distance_labels.isChecked()
        show_hydrophobic = self.parent_window.ppi_show_hydrophobic.isChecked()

        try:
            from ...ppi_analyzer import analyze_protein_protein_interface
            print(f"[GLINT] 🔬 Running PPI analysis: {obj_name} ({p1} vs {p2})")
            result = analyze_protein_protein_interface(
                obj_name=obj_name,
                protein1_chains=protein1_chains,
                protein2_chains=protein2_chains,
                interface_distance=interface_dist,
                output_csv=output_csv,
                visualize=True,
                display_mode=display_mode,
                show_labels=show_labels,
                show_distance_labels=show_dist_labels,
                show_hydrophobic=show_hydrophobic,
            )
            if result:
                contacts = result.get("interface_contacts", 0)
                strength = result.get("interface_strength", 0)
                strong = result.get("is_strong_interface", False)
                msg = (f"PPI Analysis Complete\n\n"
                       f"Interface contacts: {contacts}\n"
                       f"Interface strength: {strength:.1f}/10\n"
                       f"Strong interface: {'Yes' if strong else 'No'}")
                if output_csv:
                    msg += f"\nResults saved to: {output_csv}"
                show_message_box(self, "PPI Analysis", msg, "info")
            else:
                show_message_box(self, "PPI Analysis", "No interface detected between the specified chains.", "info")
        except Exception as e:
            traceback.print_exc()
            show_message_box(self, "Error", f"PPI analysis failed:\n{str(e)}", "critical")

    def run_pl_analysis(self):
        """运行蛋白-配体相互作用分析"""
        try:
            from pymol import cmd
        except ImportError:
            show_message_box(self, "Error", "PyMOL not available", "critical")
            return

        obj_name = self.parent_window.pl_obj_combo.currentText().strip()
        if not obj_name:
            show_message_box(self, "Warning", "Please select a target object.", "warning")
            return

        ligand_name = self.parent_window.pl_ligand_name.text().strip() or None
        chains_text = self.parent_window.pl_protein_chains.text().strip()
        protein_chains = [c.strip() for c in chains_text.split(",") if c.strip()] if chains_text else None

        try:
            distance = float(self.parent_window.pl_distance.text())
        except ValueError:
            show_message_box(self, "Warning", "Invalid distance value.", "warning")
            return

        output_csv = self.parent_window.pl_csv.text().strip() or None

        # 解析 min_confidence
        conf_text = self.parent_window.pl_min_confidence.currentText()
        try:
            min_confidence = float(conf_text.split("(")[0].strip())
        except ValueError:
            min_confidence = 0.8

        show_hydrophobic = self.parent_window.pl_show_hydrophobic.isChecked()
        show_dist_labels = self.parent_window.pl_show_distance_labels.isChecked()

        try:
            from ...interaction_analyzer import (
                analyze_protein_ligand_interactions,
                visualize_protein_ligand_3d,
            )
            print(f"[GLINT] 🔬 Running Protein-Ligand analysis: {obj_name} (ligand={ligand_name})")
            result = analyze_protein_ligand_interactions(
                obj_name=obj_name,
                ligand_resname=ligand_name,
                protein_chains=protein_chains,
                output_csv=output_csv,
                distance_cutoff=distance,
            )
            if result and result.get("interactions"):
                n = len(result["interactions"])
                print(f"[GLINT] ✅ Found {n} interactions")

                # 3D 可视化
                visualize_protein_ligand_3d(
                    obj_name=obj_name,
                    interactions_result=result,
                    ligand_resname=ligand_name,
                    show_hydrophobic=show_hydrophobic,
                    min_confidence=min_confidence,
                    show_distance_labels=show_dist_labels,
                )

                msg = f"Found {n} protein-ligand interactions."
                if output_csv:
                    msg += f"\nResults saved to: {output_csv}"
                show_message_box(self, "Protein-Ligand Analysis", msg, "info")
            else:
                show_message_box(self, "Protein-Ligand Analysis", "No interactions detected.", "info")
        except Exception as e:
            traceback.print_exc()
            show_message_box(self, "Error", f"Protein-Ligand analysis failed:\n{str(e)}", "critical")

    def run_pl_2d_diagram(self):
        """生成蛋白-配体 2D 相互作用图"""
        try:
            from pymol import cmd
        except ImportError:
            show_message_box(self, "Error", "PyMOL not available", "critical")
            return

        obj_name = self.parent_window.pl_obj_combo.currentText().strip()
        if not obj_name:
            show_message_box(self, "Warning", "Please select a target object.", "warning")
            return

        ligand_name = self.parent_window.pl_ligand_name.text().strip()
        if not ligand_name:
            show_message_box(self, "Warning", "Please specify a ligand name for 2D diagram.", "warning")
            return

        csv_path = self.parent_window.pl_csv.text().strip()
        if not csv_path or not os.path.exists(csv_path):
            show_message_box(self, "Warning",
                             "Please run Protein-Ligand analysis first and specify a valid output CSV path.",
                             "warning")
            return

        conf_text = self.parent_window.pl_min_confidence.currentText()
        try:
            min_confidence = float(conf_text.split("(")[0].strip())
        except ValueError:
            min_confidence = 0.8

        # 在主线程导出临时 PDB 供后台线程使用（避免后台调用 PyMOL）
        import tempfile
        tmp_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb").name
        try:
            cmd.save(tmp_pdb, obj_name)
        except Exception as e:
            show_message_box(self, "Error", f"Failed to export PDB: {e}", "critical")
            return

        # 输出路径：与 CSV 同目录
        output_path = os.path.splitext(csv_path)[0] + "_2d_diagram.png"

        self.parent_window.pl_2d_btn.setEnabled(False)
        self.parent_window.pl_2d_btn.setText("Generating...")

        self._diagram_worker = _Diagram2DWorker(
            csv_path=csv_path,
            ligand_resname=ligand_name,
            pdb_file=tmp_pdb,
            output_path=output_path,
            min_confidence=min_confidence,
        )
        self._diagram_worker.finished.connect(self._on_2d_diagram_finished)
        self._diagram_worker.error.connect(self._on_2d_diagram_error)
        self._diagram_worker.start()

    def _on_2d_diagram_finished(self, path: str):
        """2D 图表生成成功回调"""
        self.parent_window.pl_2d_btn.setEnabled(True)
        self.parent_window.pl_2d_btn.setText("Generate 2D Diagram")
        print(f"[GLINT] ✅ 2D diagram saved: {path}")
        show_message_box(self, "2D Diagram", f"2D interaction diagram saved:\n{path}", "info")

    def _on_2d_diagram_error(self, err: str):
        """2D 图表生成失败回调"""
        self.parent_window.pl_2d_btn.setEnabled(True)
        self.parent_window.pl_2d_btn.setText("Generate 2D Diagram")
        print(f"[GLINT] ❌ 2D diagram failed: {err}")
        show_message_box(self, "Error", f"2D diagram generation failed:\n{err}", "critical")

    def run_ll_analysis(self):
        """运行配体-配体相互作用分析"""
        try:
            from pymol import cmd
        except ImportError:
            show_message_box(self, "Error", "PyMOL not available", "critical")
            return

        obj_name = self.parent_window.ll_obj_combo.currentText().strip()
        if not obj_name:
            show_message_box(self, "Warning", "Please select a target object.", "warning")
            return

        sel1 = self.parent_window.ll_sel1.text().strip()
        sel2 = self.parent_window.ll_sel2.text().strip()
        if not sel1 or not sel2:
            show_message_box(self, "Warning", "Please specify both Selection 1 and Selection 2.", "warning")
            return

        try:
            cutoff = float(self.parent_window.ll_dist.text())
        except ValueError:
            show_message_box(self, "Warning", "Invalid distance value.", "warning")
            return

        output_csv = self.parent_window.ll_csv.text().strip() or None

        try:
            from ..ligand_ligand_analyzer import analyze_ligand_ligand_interactions
            print(f"[GLINT] 🔬 Running Ligand-Ligand analysis: {obj_name} ({sel1} vs {sel2})")
            interactions = analyze_ligand_ligand_interactions(
                obj_name=obj_name,
                sel1=sel1,
                sel2=sel2,
                cutoff=cutoff,
                output_csv=output_csv,
                visualize=True,
            )
            n = len(interactions) if interactions else 0
            msg = f"Found {n} ligand-ligand interactions."
            if output_csv:
                msg += f"\nResults saved to: {output_csv}"
            show_message_box(self, "Ligand-Ligand Analysis", msg, "info")
        except Exception as e:
            traceback.print_exc()
            show_message_box(self, "Error", f"Ligand-Ligand analysis failed:\n{str(e)}", "critical")

    def run_ec_analysis(self):
        """运行 EC (Electrostatic Complementarity) 分析"""
        try:
            from pymol import cmd
        except ImportError:
            show_message_box(self, "Error", "PyMOL not available", "critical")
            return

        obj_name = self.parent_window.ec_obj_combo.currentText().strip()
        if not obj_name:
            show_message_box(self, "Warning", "Please select a target object.", "warning")
            return

        ligand_name = self.parent_window.ec_ligand_name.text().strip()
        if not ligand_name:
            show_message_box(self, "Warning", "Please specify a ligand/glue name.", "warning")
            return

        # 中文注释：清理用户输入中的异常符号（如 Y˙70 里的中点），
        # 仅保留残基名常见的字母数字字符，避免 PyMOL resn 匹配失败。
        ligand_name_clean = ''.join(ch for ch in ligand_name if ch.isalnum())
        if ligand_name_clean and ligand_name_clean != ligand_name:
            print(f"[GLINT] ⚠️ Normalized ligand name: '{ligand_name}' -> '{ligand_name_clean}'")
            ligand_name = ligand_name_clean
            self.parent_window.ec_ligand_name.setText(ligand_name)

        surface_style = self.parent_window.ec_surface_style.currentText()

        try:
            from ...ligand_ec_calculator import calculate_ligand_ec
            from ...ec_visualization import visualize_ec_smooth_surface, visualize_ec_mesh_surface
            
            print(f"[GLINT] 🔬 Running EC analysis: {obj_name} (ligand={ligand_name}, style={surface_style})")
            
            # 计算 EC
            result = calculate_ligand_ec(
                obj_name=obj_name,
                ligand_resname=ligand_name,
                visualize=False  # 先计算，再根据样式可视化
            )
            
            if not result:
                show_message_box(self, "EC Analysis", "EC calculation failed. Check console for details.", "warning")
                return
            
            ec_score = result.get('ec_score', 0)
            ec_stats = result.get('ec_statistics', {})
            
            # 根据用户选择的样式进行可视化
            if surface_style == "Solid Surface":
                # 中文注释：Solid Surface 应使用光滑的 molecular surface，
                # 避免 GUI 层误传 gaussian 导致表面重新变成颗粒/泡状外观。
                visualize_ec_smooth_surface(
                    obj_name, ligand_name,
                    ec_result=result,
                    surface_type='molecular',
                    transparency=0.30
                )
            else:  # Mesh Surface
                visualize_ec_mesh_surface(
                    obj_name, ligand_name,
                    ec_result=result
                )
            
            msg = f"EC Analysis completed.\n\n"
            msg += f"EC Score: {ec_score:.4f}\n"
            if ec_stats:
                msg += f"EC Mean: {ec_stats.get('ec_mean', 0):.4f}\n"
                msg += f"EC Median: {ec_stats.get('ec_median', 0):.4f}\n"
                msg += f"Positive EC: {ec_stats.get('ec_positive_fraction', 0)*100:.1f}%\n"
                msg += f"Negative EC: {ec_stats.get('ec_negative_fraction', 0)*100:.1f}%"
            
            output_dir = result.get('output_dir')
            if output_dir:
                msg += f"\n\nOutput saved to: {output_dir}"
            
            show_message_box(self, "EC Analysis", msg, "info")
            
        except Exception as e:
            traceback.print_exc()
            show_message_box(self, "Error", f"EC analysis failed:\n{str(e)}", "critical")