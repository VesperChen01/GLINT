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
        
        # 1. Protein-Protein Interface (PPI) Analysis
        grp_ppi = QGroupBox("Protein-Protein Interface (PPI) Analysis")
        ppi_grid = QGridLayout(grp_ppi)
        ppi_grid.setColumnStretch(1, 1); ppi_grid.setColumnStretch(3, 1)
        ppi_grid.setHorizontalSpacing(8); ppi_grid.setVerticalSpacing(10)
        
        ppi_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ppi_obj_combo = QComboBox(); self.parent_window.ppi_obj_combo.setMinimumHeight(36)
        self.parent_window.ppi_refresh_btn = QPushButton(t("refresh")); self.parent_window.ppi_refresh_btn.clicked.connect(self.refresh_objects)
        r0 = QHBoxLayout(); r0.addWidget(self.parent_window.ppi_obj_combo, 1); r0.addWidget(self.parent_window.ppi_refresh_btn)
        ppi_grid.addLayout(r0, 0, 1)
        
        ppi_grid.addWidget(QLabel("Protein1 Chains:"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ppi_protein1_chains = QLineEdit(); self.parent_window.ppi_protein1_chains.setPlaceholderText("e.g. A")
        ppi_grid.addWidget(self.parent_window.ppi_protein1_chains, 0, 3)
        
        ppi_grid.addWidget(QLabel("Protein2 Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ppi_protein2_chains = QLineEdit(); self.parent_window.ppi_protein2_chains.setPlaceholderText("e.g. B")
        ppi_grid.addWidget(self.parent_window.ppi_protein2_chains, 1, 1)
        
        ppi_grid.addWidget(QLabel("Interface Dist (Å):"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ppi_interface_dist = QLineEdit("4.5")
        ppi_grid.addWidget(self.parent_window.ppi_interface_dist, 1, 3)
        
        ppi_grid.addWidget(QLabel("Output CSV:"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ppi_csv = QLineEdit(); self.parent_window.ppi_csv.setPlaceholderText("Optional")
        self.parent_window.ppi_csv_btn = QPushButton(t("browse"))
        self.parent_window.ppi_csv_btn.clicked.connect(lambda: self._browse_save_file(self.parent_window.ppi_csv, "CSV (*.csv)"))
        r2 = QHBoxLayout(); r2.addWidget(self.parent_window.ppi_csv, 1); r2.addWidget(self.parent_window.ppi_csv_btn)
        ppi_grid.addLayout(r2, 2, 1, 1, 3)

        # Row 3: Visualization Options
        ppi_grid.addWidget(QLabel("Display Mode:"), 3, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ppi_display_mode = QComboBox()
        self.parent_window.ppi_display_mode.addItems(["Cartoon + Surface + Interactions", "Surface + Interactions"])
        ppi_grid.addWidget(self.parent_window.ppi_display_mode, 3, 1)
        
        self.parent_window.ppi_show_labels = QCheckBox("Show Distance Labels")
        self.parent_window.ppi_show_labels.setChecked(True)
        ppi_grid.addWidget(self.parent_window.ppi_show_labels, 3, 3)
        
        ppi_btn_row = QHBoxLayout()
        self.parent_window.ppi_analyze_btn = QPushButton("Analyze PPI Interface"); self.parent_window.ppi_analyze_btn.setObjectName("highlight_btn")
        self.parent_window.ppi_analyze_btn.clicked.connect(self.run_ppi_analysis)
        ppi_btn_row.addWidget(self.parent_window.ppi_analyze_btn)
        ppi_btn_row.addStretch(1)
        
        layout.addWidget(grp_ppi)
        layout.addLayout(ppi_btn_row)
        
        # 2. Protein-Ligand Interactions
        grp_pl = QGroupBox("Protein-Ligand Interactions")
        pl_grid = QGridLayout(grp_pl)
        pl_grid.setColumnStretch(1, 1); pl_grid.setColumnStretch(3, 1)
        pl_grid.setHorizontalSpacing(8); pl_grid.setVerticalSpacing(10)
        
        pl_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.pl_obj_combo = QComboBox(); self.parent_window.pl_obj_combo.setMinimumWidth(150); self.parent_window.pl_obj_combo.setMinimumHeight(36)
        self.parent_window.pl_refresh_btn = QPushButton(t("refresh")); self.parent_window.pl_refresh_btn.clicked.connect(self.refresh_objects)
        r0_pl = QHBoxLayout(); r0_pl.addWidget(self.parent_window.pl_obj_combo, 1); r0_pl.addWidget(self.parent_window.pl_refresh_btn)
        pl_grid.addLayout(r0_pl, 0, 1)
        
        pl_grid.addWidget(QLabel("Ligand Name:"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.pl_ligand_name = QLineEdit(); self.parent_window.pl_ligand_name.setPlaceholderText("e.g. CC885, Auto-detect if blank")
        pl_grid.addWidget(self.parent_window.pl_ligand_name, 0, 3)
        
        pl_grid.addWidget(QLabel("Protein Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.pl_protein_chains = QLineEdit(); self.parent_window.pl_protein_chains.setPlaceholderText("e.g. A,B (optional)")
        pl_grid.addWidget(self.parent_window.pl_protein_chains, 1, 1)
        
        pl_grid.addWidget(QLabel("Distance (Å):"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.pl_distance = QLineEdit("4.5")
        pl_grid.addWidget(self.parent_window.pl_distance, 1, 3)
        
        pl_grid.addWidget(QLabel("Output CSV:"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.pl_csv = QLineEdit(); self.parent_window.pl_csv.setPlaceholderText("Optional")
        self.parent_window.pl_csv_btn = QPushButton(t("browse"))
        self.parent_window.pl_csv_btn.clicked.connect(lambda: self._browse_save_file(self.parent_window.pl_csv, "CSV (*.csv)"))
        r2_pl = QHBoxLayout(); r2_pl.addWidget(self.parent_window.pl_csv, 1); r2_pl.addWidget(self.parent_window.pl_csv_btn)
        pl_grid.addLayout(r2_pl, 2, 1, 1, 3)
        
        pl_btn_row = QHBoxLayout()
        self.parent_window.pl_analyze_btn = QPushButton("Analyze Protein-Ligand"); self.parent_window.pl_analyze_btn.setObjectName("highlight_btn")
        self.parent_window.pl_analyze_btn.clicked.connect(self.run_pl_analysis)
        
        self.parent_window.pl_2d_btn = QPushButton("Generate 2D Diagram")
        self.parent_window.pl_2d_btn.clicked.connect(self.run_pl_2d_diagram)
        
        pl_btn_row.addWidget(self.parent_window.pl_analyze_btn)
        pl_btn_row.addWidget(self.parent_window.pl_2d_btn)
        pl_btn_row.addStretch(1)
        
        layout.addWidget(grp_pl)
        layout.addLayout(pl_btn_row)
        
        # Hidden fields for backward compatibility
        self.parent_window.pl_show_hydrophobic = QCheckBox(); self.parent_window.pl_show_hydrophobic.setVisible(False)
        self.parent_window.pl_min_confidence = QComboBox(); self.parent_window.pl_min_confidence.setVisible(False)
        
        # 3. Ligand-Ligand Interactions
        grp_ll = QGroupBox("Ligand-Ligand Interactions (Small Molecule - Small Molecule)")
        ll_grid = QGridLayout(grp_ll)
        ll_grid.setColumnStretch(1, 1); ll_grid.setColumnStretch(3, 1)
        ll_grid.setHorizontalSpacing(8); ll_grid.setVerticalSpacing(10)
        
        ll_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ll_obj_combo = QComboBox(); self.parent_window.ll_obj_combo.setMinimumWidth(150); self.parent_window.ll_obj_combo.setMinimumHeight(36)
        self.parent_window.ll_refresh_btn = QPushButton(t("refresh")); self.parent_window.ll_refresh_btn.clicked.connect(self.refresh_objects)
        r0_ll = QHBoxLayout(); r0_ll.addWidget(self.parent_window.ll_obj_combo, 1); r0_ll.addWidget(self.parent_window.ll_refresh_btn)
        ll_grid.addLayout(r0_ll, 0, 1)
        
        ll_grid.addWidget(QLabel("PDB File (opt):"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ll_pdb = QLineEdit(); self.parent_window.ll_pdb.setPlaceholderText("Load from file...")
        self.parent_window.ll_pdb_browse = QPushButton(t("browse"))
        self.parent_window.ll_pdb_browse.clicked.connect(lambda: self._browse_file(self.parent_window.ll_pdb, "PDB Files (*.pdb *.cif *.sdf)"))
        r0b_ll = QHBoxLayout(); r0b_ll.addWidget(self.parent_window.ll_pdb, 1); r0b_ll.addWidget(self.parent_window.ll_pdb_browse)
        ll_grid.addLayout(r0b_ll, 0, 3)
        
        ll_grid.addWidget(QLabel("Selection 1:"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ll_sel1 = QLineEdit(); self.parent_window.ll_sel1.setPlaceholderText("e.g. resn LIG1 or resi 100")
        ll_grid.addWidget(self.parent_window.ll_sel1, 1, 1)
        
        ll_grid.addWidget(QLabel("Selection 2:"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ll_sel2 = QLineEdit(); self.parent_window.ll_sel2.setPlaceholderText("e.g. resn LIG2 or resi 200")
        ll_grid.addWidget(self.parent_window.ll_sel2, 1, 3)
        
        ll_grid.addWidget(QLabel("Distance (Å):"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ll_dist = QLineEdit("4.5")
        ll_grid.addWidget(self.parent_window.ll_dist, 2, 1)
        
        ll_grid.addWidget(QLabel("Output CSV:"), 2, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ll_csv = QLineEdit(); self.parent_window.ll_csv.setPlaceholderText("Optional")
        self.parent_window.ll_csv_btn = QPushButton(t("browse"))
        self.parent_window.ll_csv_btn.clicked.connect(lambda: self._browse_save_file(self.parent_window.ll_csv, "CSV (*.csv)"))
        r2_ll = QHBoxLayout(); r2_ll.addWidget(self.parent_window.ll_csv, 1); r2_ll.addWidget(self.parent_window.ll_csv_btn)
        ll_grid.addLayout(r2_ll, 2, 3)
        
        ll_btn_row = QHBoxLayout()
        self.parent_window.ll_analyze_btn = QPushButton("Analyze Ligand-Ligand"); self.parent_window.ll_analyze_btn.setObjectName("highlight_btn")
        self.parent_window.ll_analyze_btn.clicked.connect(self.run_ll_analysis)
        ll_btn_row.addWidget(self.parent_window.ll_analyze_btn)
        ll_btn_row.addStretch(1)
        
        layout.addWidget(grp_ll)
        layout.addLayout(ll_btn_row)
        
        # 4. Mutation Analysis
        layout.addWidget(self._create_mutation_analysis_card())
        
        layout.addStretch(1)
        scroll_area.setWidget(content_widget)
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(scroll_area)

    # --- PPI Analysis Logic ---
    def run_ppi_analysis(self):
        """分析蛋白-蛋白界面 (PPI Interface)"""
        try:
            obj = self.parent_window.ppi_obj_combo.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select a structure object")
                return
            
            p1_chains = self.parent_window.ppi_protein1_chains.text().strip()
            p2_chains = self.parent_window.ppi_protein2_chains.text().strip()
            
            if not p1_chains or not p2_chains:
                QMessageBox.warning(self, "Warning", "Please specify both protein chain groups")
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
            display_mode = "backbone_surface" if display_mode_idx == 0 else "surface_only"
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
                
                QMessageBox.information(self, "PPI Analysis Complete",
                    f"Interface Contacts: {contacts}\n"
                    f"{'BSA: ' + str(round(bsa, 1)) + ' Ų' if bsa else 'BSA: N/A'}\n"
                    f"Interface Strength: {strength:.1f}/10\n\n"
                    f"{'Strong Interface' if is_strong else 'Weak Interface'}")
            else:
                self.log("PPI analysis failed")
        except Exception as e:
            self.on_error(str(e))
            import traceback; traceback.print_exc()
    


    # --- PL Logic ---
    def run_pl_analysis(self):
        """分析蛋白-配体相互作用"""
        try:
            obj_name = self.parent_window.pl_obj_combo.currentText()
            if obj_name == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select a structure object")
                return
            
            ligand_name = self.parent_window.pl_ligand_name.text().strip() or None
            protein_chains_str = self.parent_window.pl_protein_chains.text().strip()
            protein_chains = [c.strip() for c in protein_chains_str.split(",")] if protein_chains_str else None
            distance = float(self.parent_window.pl_distance.text())
            output_csv = self.parent_window.pl_csv.text().strip() or None
            
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
                
                # ✅ 添加 3D 可视化
                if n > 0:
                    try:
                        try: from ...interaction_analyzer import visualize_protein_ligand_3d
                        except ImportError: from interaction_analyzer import visualize_protein_ligand_3d
                        
                        # 获取实际配体名称
                        actual_ligand = ligand_name
                        if result.get("ligand_residues"):
                            actual_ligand = result["ligand_residues"][0].get("resname", ligand_name)
                        
                        visualize_protein_ligand_3d(obj_name, result, actual_ligand)
                        self.log("  3D visualization generated")
                    except Exception as viz_e:
                        self.log(f"  3D visualization failed: {viz_e}")
                        import traceback; traceback.print_exc()
                
                QMessageBox.information(self, "Analysis Complete", 
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
                QMessageBox.critical(self, "Error", "RDKit is required for 2D diagrams.\nPlease install it: pip install rdkit")
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
                        QMessageBox.warning(self, "Warning", "No interactions to plot.")
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
                    QMessageBox.warning(self, "Missing Data", "Please run analysis first (with output CSV) or select an existing CSV file.")
                    return

            if not ligand_name:
                QMessageBox.warning(self, "Missing Input", "Please specify the Ligand Name (Residue Name).")
                return

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
                output_path=out_path
            )
            
            if final_path and os.path.exists(final_path):
                self.log(f"2D Diagram saved: {final_path}")
                QMessageBox.information(self, "Success", f"2D Diagram saved to:\n{final_path}")
                
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
                QMessageBox.warning(self, "Error", "Failed to generate diagram. See log for details.")

        except Exception as e:
            self.on_error(f"2D Diagram Error: {str(e)}")
            import traceback; traceback.print_exc()

    # --- Browse file helpers ---
    def _browse_save_file(self, line_edit, file_filter):
        """浏览并选择保存文件路径"""
        fn, _ = QFileDialog.getSaveFileName(self, "Save File", "", file_filter)
        if fn:
            line_edit.setText(fn)
    
    def _browse_file(self, line_edit, file_filter):
        """浏览并选择打开文件路径"""
        fn, _ = QFileDialog.getOpenFileName(self, "Open File", "", file_filter)
        if fn:
            line_edit.setText(fn)
    
    # --- LL (Ligand-Ligand) Logic ---
    def run_ll_analysis(self):
        """运行配体-配体相互作用分析"""
        try:
            # 尝试导入分析函数
            try: from ...ligand_ligand_analyzer import analyze_ligand_ligand_interactions
            except ImportError:
                try: from ligand_ligand_analyzer import analyze_ligand_ligand_interactions
                except ImportError:
                    QMessageBox.critical(self, "Error", "Ligand-Ligand analysis module not found.")
                    return
            
            # 检查PDB文件
            pdb_file = self.parent_window.ll_pdb.text().strip()
            obj = self.parent_window.ll_obj_combo.currentText()
            
            # 如果提供了PDB文件，先加载它
            if pdb_file and os.path.exists(pdb_file):
                try:
                    from pymol import cmd
                    loaded_obj = os.path.basename(pdb_file).split('.')[0]
                    # 确保唯一名称
                    if hasattr(cmd, 'get_unused_name'):
                        loaded_obj = cmd.get_unused_name(loaded_obj)
                    cmd.load(pdb_file, loaded_obj)
                    obj = loaded_obj
                    self.log(f"Loaded {pdb_file} as {obj}")
                    self.refresh_objects()
                except Exception as e:
                    QMessageBox.critical(self, "Error", f"Failed to load PDB file: {e}")
                    return
            
            sel1 = self.parent_window.ll_sel1.text().strip()
            sel2 = self.parent_window.ll_sel2.text().strip()
            dist_str = self.parent_window.ll_dist.text().strip()
            csv_path = self.parent_window.ll_csv.text().strip() or None
            
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, "Missing Input", "Please select a target object or load a PDB file.")
                return
            if not sel1 or not sel2:
                QMessageBox.warning(self, "Missing Input", "Please define both Selection 1 and Selection 2.")
                return
            
            try:
                dist = float(dist_str) if dist_str else 4.5
            except ValueError:
                QMessageBox.warning(self, "Invalid Input", "Distance must be a number.")
                return
            
            self.log(f"Starting Ligand-Ligand analysis for {obj}...")
            self.log(f"  Selection 1: {sel1}")
            self.log(f"  Selection 2: {sel2}")
            self.log(f"  Distance cutoff: {dist} Å")
            
            try:
                import pymol
                
                # 构建完整的PyMOL selection
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
                
                QMessageBox.information(self, "Success", msg)
                
            except pymol.CmdException as e:
                msg = str(e)
                if "Invalid selection" in msg:
                    QMessageBox.critical(self, "Selection Error",
                        f"PyMOL could not understand your selection.\n\n"
                        f"Error: {msg}\n\n"
                        f"Tip: Please use valid PyMOL selection syntax.\n"
                        f"Examples:\n"
                        f"\u2022 resn LIG (by residue name)\n"
                        f"\u2022 resi 900 (by residue index)\n"
                        f"\u2022 chain A (by chain)\n\n"
                        f"You entered: '{sel1}' and '{sel2}'")
                else:
                    QMessageBox.critical(self, "PyMOL Error", str(e))
                return
                
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
        # Placeholder for full analysis
        self.log("Full mutation analysis requires external tools (FoldX/PyRosetta).")
