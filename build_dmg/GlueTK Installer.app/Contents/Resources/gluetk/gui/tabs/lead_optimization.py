# -*- coding: utf-8 -*-
"""
Lead Optimization: PPI, Glue, Ternary, Mutation
"""
import os
import traceback
from typing import Optional

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

class LeadOptimizationTab(CommonTab):
    def __init__(self, parent):
        super().__init__(parent)
        
        # Remove proxies for methods implemented here
        for attr in ['run_mutation', 'run_minimize', 'run_mutation_analysis']:
            if attr in self.__dict__:
                del self.__dict__[attr]
                
        self.init_ui()
        
    def init_ui(self):
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        
        content_widget = QWidget()
        content_widget.setObjectName("scroll_content")
        self.parent_window._lead_scroll_content = content_widget
        bg_color = "#0d1117" if getattr(self.parent_window, "_dark_mode", True) else "#f8fafc"
        content_widget.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")
        
        layout = QVBoxLayout(content_widget)
        layout.setSpacing(14)
        layout.setContentsMargins(12, 12, 12, 12)
        
        # 1. Molecular Glue Specifics
        grp_glue = QGroupBox("Molecular Glue Analysis (PPI + Neo-Epitope)")
        glue_grid = QGridLayout(grp_glue)
        glue_grid.setColumnStretch(1, 1); glue_grid.setColumnStretch(3, 1)
        glue_grid.setHorizontalSpacing(8); glue_grid.setVerticalSpacing(10)
        
        glue_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.glue_obj_combo = QComboBox(); self.parent_window.glue_obj_combo.setMinimumHeight(36)
        self.parent_window.glue_refresh_btn = QPushButton(t("refresh")); self.parent_window.glue_refresh_btn.clicked.connect(self.refresh_objects)
        r0 = QHBoxLayout(); r0.addWidget(self.parent_window.glue_obj_combo, 1); r0.addWidget(self.parent_window.glue_refresh_btn)
        glue_grid.addLayout(r0, 0, 1)
        
        glue_grid.addWidget(QLabel("Glue Resname:"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.glue_resname = QLineEdit(); self.parent_window.glue_resname.setPlaceholderText("e.g. CC885")
        glue_grid.addWidget(self.parent_window.glue_resname, 0, 3)
        
        glue_grid.addWidget(QLabel("E3 Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.glue_e3_chains = QLineEdit(); self.parent_window.glue_e3_chains.setPlaceholderText("e.g. A")
        glue_grid.addWidget(self.parent_window.glue_e3_chains, 1, 1)
        
        glue_grid.addWidget(QLabel("Substrate Chains:"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.glue_sub_chains = QLineEdit(); self.parent_window.glue_sub_chains.setPlaceholderText("e.g. B")
        glue_grid.addWidget(self.parent_window.glue_sub_chains, 1, 3)
        
        glue_grid.addWidget(QLabel("Interface Dist (Å):"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.glue_interface_dist = QLineEdit("4.5")
        glue_grid.addWidget(self.parent_window.glue_interface_dist, 2, 1)
        
        glue_grid.addWidget(QLabel("Neo-Epitope Dist (Å):"), 2, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.glue_neo_dist = QLineEdit("5.0")
        glue_grid.addWidget(self.parent_window.glue_neo_dist, 2, 3)
        
        glue_btn_row = QHBoxLayout()
        self.parent_window.glue_full_btn = QPushButton("Full Glue Analysis"); self.parent_window.glue_full_btn.setObjectName("highlight_btn")
        self.parent_window.glue_full_btn.clicked.connect(self.run_glue_full_analysis)
        glue_btn_row.addWidget(self.parent_window.glue_full_btn)
        glue_btn_row.addStretch(1)
        
        layout.addWidget(grp_glue)
        layout.addLayout(glue_btn_row)
        
        # 2. Ternary Complex
        grp_ternary = QGroupBox("Ternary Complex Analysis")
        t_grid = QGridLayout(grp_ternary)
        t_grid.setColumnStretch(1, 1); t_grid.setColumnStretch(3, 1)
        
        t_grid.addWidget(QLabel("Target Object:"), 0, 0)
        self.parent_window.tc_obj_combo = QComboBox()
        t_grid.addWidget(self.parent_window.tc_obj_combo, 0, 1)
        
        t_grid.addWidget(QLabel("Ligand:"), 0, 2)
        self.parent_window.tc_ligand_name = QLineEdit(); self.parent_window.tc_ligand_name.setPlaceholderText("Auto")
        t_grid.addWidget(self.parent_window.tc_ligand_name, 0, 3)
        
        t_grid.addWidget(QLabel("E3 Chains:"), 1, 0)
        self.parent_window.tc_protein1_chains = QLineEdit(); self.parent_window.tc_protein1_chains.setPlaceholderText("e.g. A")
        t_grid.addWidget(self.parent_window.tc_protein1_chains, 1, 1)
        
        t_grid.addWidget(QLabel("POI Chains:"), 1, 2)
        self.parent_window.tc_protein2_chains = QLineEdit(); self.parent_window.tc_protein2_chains.setPlaceholderText("e.g. B")
        t_grid.addWidget(self.parent_window.tc_protein2_chains, 1, 3)
        
        t_btn_row = QHBoxLayout()
        self.parent_window.tc_analyze_btn = QPushButton("Analyze Complex"); self.parent_window.tc_analyze_btn.setObjectName("highlight_btn")
        self.parent_window.tc_analyze_btn.clicked.connect(self.run_tc_analysis)
        self.parent_window.tc_render_btn = QPushButton("Render All"); self.parent_window.tc_render_btn.setObjectName("highlight_btn")
        self.parent_window.tc_render_btn.clicked.connect(self.run_tc_render)
        t_btn_row.addWidget(self.parent_window.tc_analyze_btn); t_btn_row.addWidget(self.parent_window.tc_render_btn); t_btn_row.addStretch(1)
        
        layout.addWidget(grp_ternary)
        layout.addLayout(t_btn_row)

        # 3. Protein-Protein Docking (HDOCK)
        grp_dock = QGroupBox("Protein-Protein Docking (Ternary Complex Prediction)")
        dock_grid = QGridLayout(grp_dock)
        dock_grid.setColumnStretch(1, 1)
        
        dock_grid.addWidget(QLabel("Receptor (E3+Glue):"), 0, 0)
        self.parent_window.hdock_receptor = QLineEdit(); self.parent_window.hdock_receptor.setPlaceholderText("PDB File")
        self.parent_window.hdock_rec_browse = QPushButton("Browse"); self.parent_window.hdock_rec_browse.clicked.connect(self.browse_hdock_receptor)
        r_box = QHBoxLayout(); r_box.addWidget(self.parent_window.hdock_receptor); r_box.addWidget(self.parent_window.hdock_rec_browse)
        dock_grid.addLayout(r_box, 0, 1)
        
        dock_grid.addWidget(QLabel("Ligand (POI):"), 1, 0)
        self.parent_window.hdock_ligand = QLineEdit(); self.parent_window.hdock_ligand.setPlaceholderText("PDB File")
        self.parent_window.hdock_lig_browse = QPushButton("Browse"); self.parent_window.hdock_lig_browse.clicked.connect(self.browse_hdock_ligand)
        l_box = QHBoxLayout(); l_box.addWidget(self.parent_window.hdock_ligand); l_box.addWidget(self.parent_window.hdock_lig_browse)
        dock_grid.addLayout(l_box, 1, 1)
        
        self.parent_window.hdock_run_btn = QPushButton("Run HDOCK Prediction")
        self.parent_window.hdock_run_btn.setObjectName("highlight_btn")
        self.parent_window.hdock_run_btn.clicked.connect(self.run_hdock)
        dock_grid.addWidget(self.parent_window.hdock_run_btn, 2, 1)
        
        layout.addWidget(grp_dock)
        
        # 4. General Interactions (PP, PL)
        grp_pl = QGroupBox("Protein-Ligand Interactions")
        pl_layout = QHBoxLayout(grp_pl)
        self.parent_window.pl_obj_combo = QComboBox(); self.parent_window.pl_obj_combo.setMinimumWidth(150)
        self.parent_window.pl_analyze_btn = QPushButton("Analyze PL"); self.parent_window.pl_analyze_btn.clicked.connect(self.run_pl_analysis)
        self.parent_window.pl_refresh_btn = QPushButton("R"); self.parent_window.pl_refresh_btn.clicked.connect(self.refresh_objects)
        
        pl_layout.addWidget(QLabel("Object:")); pl_layout.addWidget(self.parent_window.pl_obj_combo)
        pl_layout.addWidget(self.parent_window.pl_refresh_btn)
        pl_layout.addWidget(self.parent_window.pl_analyze_btn)
        layout.addWidget(grp_pl)
        
        # Hidden fields required for PL logic but not fully exposed in this simplified view
        self.parent_window.pl_ligand_name = QLineEdit(); self.parent_window.pl_ligand_name.setVisible(False)
        self.parent_window.pl_protein_chains = QLineEdit(); self.parent_window.pl_protein_chains.setVisible(False)
        self.parent_window.pl_distance = QLineEdit("4.5"); self.parent_window.pl_distance.setVisible(False)
        self.parent_window.pl_csv = QLineEdit(); self.parent_window.pl_csv.setVisible(False)
        self.parent_window.pl_show_hydrophobic = QCheckBox(); self.parent_window.pl_show_hydrophobic.setVisible(False)
        self.parent_window.pl_min_confidence = QComboBox(); self.parent_window.pl_min_confidence.setVisible(False)
        
        # 4. Mutation Analysis
        layout.addWidget(self._create_mutation_analysis_card())
        
        layout.addStretch(1)
        scroll_area.setWidget(content_widget)
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(scroll_area)

    # --- Glue Logic ---
    def run_glue_full_analysis(self):
        # Placeholder logic based on original unified_gui.py
        # Actual implementation needs imports from pocket_glue_integration etc.
        # For now, we log and import what we can
        try:
            obj = self.parent_window.glue_obj_combo.currentText().strip()
            glue_res = self.parent_window.glue_resname.text().strip()
            e3_chain = self.parent_window.glue_e3_chains.text().strip()
            sub_chain = self.parent_window.glue_sub_chains.text().strip()
            
            if not obj or obj == t("no_object") or not glue_res or not e3_chain or not sub_chain:
                QMessageBox.warning(self, "Warning", "Please fill all fields")
                return
                
            self.log(f"Starting Glue Analysis for {obj}...")
            try: from ...pocket_glue_integration import analyze_molecular_glue
            except ImportError: from pocket_glue_integration import analyze_molecular_glue
            
            results = analyze_molecular_glue(obj, glue_res, e3_chain, sub_chain,
                                           interface_dist_cutoff=float(self.parent_window.glue_interface_dist.text()),
                                           neo_epitope_dist_cutoff=float(self.parent_window.glue_neo_dist.text()))
            self.log("Glue Analysis Complete. See results in PyMOL and logs.")
        except Exception as e:
            self.on_error(str(e))

    # --- Ternary Complex Logic ---
    def run_tc_analysis(self):
        try:
            obj = self.parent_window.tc_obj_combo.currentText().strip()
            if not obj or obj == t("no_object"): return
            
            try: from ...interaction_analyzer import analyze_ternary_complex
            except ImportError: from interaction_analyzer import analyze_ternary_complex
            
            analyze_ternary_complex(
                obj, 
                ligand_resname=self.parent_window.tc_ligand_name.text().strip() or None,
                protein1_chains=self.parent_window.tc_protein1_chains.text().strip() or None,
                protein2_chains=self.parent_window.tc_protein2_chains.text().strip() or None
            )
            self.log(f"Ternary analysis for {obj} complete.")
        except Exception as e:
            self.on_error(str(e))

    def run_tc_render(self):
        # Just call analysis then maybe some visualization
        self.run_tc_analysis()

    # --- HDOCK Logic ---
    def browse_hdock_receptor(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Receptor PDB (E3)", "", "PDB (*.pdb)")
        if fn: self.parent_window.hdock_receptor.setText(fn)

    def browse_hdock_ligand(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Ligand PDB (POI)", "", "PDB (*.pdb)")
        if fn: self.parent_window.hdock_ligand.setText(fn)

    def run_hdock(self):
        rec = self.parent_window.hdock_receptor.text().strip()
        lig = self.parent_window.hdock_ligand.text().strip()
        
        if not rec or not lig:
            QMessageBox.warning(self, "Warning", "Please select both Receptor and Ligand PDB files.")
            return
            
        if not os.path.exists(rec) or not os.path.exists(lig):
            QMessageBox.warning(self, "Warning", "Selected files do not exist.")
            return
            
        self.log(f"Starting HDOCK...")
        self.log(f"   Receptor: {os.path.basename(rec)}")
        self.log(f"   Ligand: {os.path.basename(lig)}")
        
        try:
            from ...hdock_integration import check_hdock_available, HDockRunner
            hdock_path = check_hdock_available()
            
            if not hdock_path:
                QMessageBox.warning(self, "Error", "HDOCKlite not found. Please ensure 'HDOCKlite-v1.1' is in the plugin directory.")
                return
                
            runner = HDockRunner(hdock_path)
            
            # Run in background thread ideally, but for now synchronous with UI update
            self.parent_window.hdock_run_btn.setEnabled(False)
            self.parent_window.hdock_run_btn.setText("Running HDOCK...")
            self.parent_window.repaint() # Force UI update
            
            result = runner.run_docking(rec, lig)
            
            self.parent_window.hdock_run_btn.setEnabled(True)
            self.parent_window.hdock_run_btn.setText("Run HDOCK Prediction")
            
            if result['success']:
                self.log("HDOCK complete!")
                self.log(f"   Models: {result['models_pdb']}")
                
                # Load into PyMOL
                from pymol import cmd
                cmd.load(result['models_pdb'], "hdock_models")
                self.log("   Loaded 'hdock_models' into PyMOL")
                QMessageBox.information(self, "Success", "Docking complete. Top 10 models loaded.")
            else:
                self.log(f"HDOCK failed: {result['error']}")
                QMessageBox.critical(self, "Error", f"HDOCK failed: {result['error']}")
                
        except Exception as e:
            self.log(f"Error: {e}")
            self.on_error(str(e))
            self.parent_window.hdock_run_btn.setEnabled(True)
            self.parent_window.hdock_run_btn.setText("Run HDOCK Prediction")

    # --- PL Logic ---
    def run_pl_analysis(self):
        try:
            obj_name = self.parent_window.pl_obj_combo.currentText()
            if obj_name == t("no_object"): return
            
            try: from ...interaction_analyzer import analyze_protein_ligand_interactions
            except ImportError: from interaction_analyzer import analyze_protein_ligand_interactions
            
            result = analyze_protein_ligand_interactions(
                obj_name=obj_name,
                ligand_resname=self.parent_window.pl_ligand_name.text().strip() or None,
                protein_chains=None,
                distance_cutoff=4.5,
                output_csv=None
            )
            
            if result:
                n = len(result.get("interactions", []))
                self.log(f"Found {n} interactions for {obj_name}")
                self.parent_window.current_pl_result = result
                QMessageBox.information(self, "Done", f"Found {n} interactions")
            else:
                self.log("No interactions found")
        except Exception as e:
            self.on_error(str(e))

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
        # Placeholder for full analysis
        self.log("Full mutation analysis requires external tools (FoldX/PyRosetta).")
