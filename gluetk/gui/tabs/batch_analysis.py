# -*- coding: utf-8 -*-
"""
Batch Analysis Tab: Multi-structure analysis
"""
import os
from typing import Optional, List

from ..qt_adapter import (
    Qt, QThread, Signal as pyqtSignal,
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QGroupBox, QScrollArea, QFrame, QFileDialog, QMessageBox,
    QComboBox, QCheckBox, QTextEdit, QProgressBar, QListWidget,
    QListWidgetItem, QAbstractItemView, QSpinBox
)

from ..utils import t
from .common import CommonTab


class BatchWorker(QThread):
    """Batch analysis worker thread"""
    progress = pyqtSignal(str)
    error = pyqtSignal(str)
    finished = pyqtSignal(dict)
    
    def __init__(self, analysis_type: str, pdb_sources: List[str], params: dict):
        super().__init__()
        self.analysis_type = analysis_type
        self.pdb_sources = pdb_sources
        self.params = params
        
    def run(self):
        try:
            from ...batch_analyzer import BatchAnalyzer
            
            analyzer = BatchAnalyzer(
                output_dir=self.params.get("output_dir"),
                max_workers=self.params.get("max_workers", 4)
            )
            
            if self.analysis_type == "gmotif":
                result = analyzer.batch_gmotif_detection(
                    self.pdb_sources,
                    rmsd_cutoff=self.params.get("rmsd_cutoff", 3.5),
                    require_gly=self.params.get("require_gly", True),
                    template=self.params.get("template", "GSPT1"),
                    output_csv=self.params.get("output_csv")
                )
            elif self.analysis_type == "ppi":
                result = analyzer.batch_ppi_analysis(
                    self.pdb_sources,
                    chain_pairs=self.params.get("chain_pairs", []),
                    interface_distance=self.params.get("interface_distance", 4.5),
                    output_csv=self.params.get("output_csv")
                )
            elif self.analysis_type == "pockets":
                result = analyzer.batch_pocket_detection(
                    self.pdb_sources,
                    min_volume=self.params.get("min_volume", 30.0),
                    min_depth=self.params.get("min_depth", 2.5),
                    output_csv=self.params.get("output_csv")
                )
            elif self.analysis_type == "interactions":
                result = analyzer.batch_interaction_analysis(
                    self.pdb_sources,
                    ligand_resnames=self.params.get("ligand_resnames"),
                    distance_cutoff=self.params.get("distance_cutoff", 4.5),
                    output_csv=self.params.get("output_csv")
                )
            elif self.analysis_type == "comprehensive":
                result = analyzer.comprehensive_batch_analysis(
                    self.pdb_sources,
                    chain_pairs=self.params.get("chain_pairs"),
                    ligand_resnames=self.params.get("ligand_resnames"),
                    output_dir=self.params.get("output_dir")
                )
            else:
                result = {"success": False, "error": f"Unknown analysis type: {self.analysis_type}"}
            
            self.finished.emit(result)
            
        except Exception as e:
            import traceback
            self.error.emit(f"{str(e)}\n{traceback.format_exc()}")


class BatchAnalysisTab(CommonTab):
    """Batch analysis tab"""
    
    def __init__(self, parent):
        super().__init__(parent)
        self._pdb_sources = []
        self._batch_thread = None
        self.init_ui()
        
    def init_ui(self):
        self.setObjectName("scroll_content")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)

        bg_color = "#0d1117" if getattr(self.parent_window, "_dark_mode", True) else "#f8fafc"
        self.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")

        layout = QVBoxLayout(self)
        layout.setSpacing(16)
        layout.setContentsMargins(12, 12, 12, 12)
        
        # 1. Input file selection
        grp_input = QGroupBox("Input Structures")
        input_layout = QVBoxLayout(grp_input)
        input_layout.setSpacing(12)
        
        # PDB ID input
        pdb_id_row = QHBoxLayout()
        pdb_id_row.addWidget(QLabel("PDB IDs:"))
        self.pdb_id_input = QLineEdit()
        self.pdb_id_input.setPlaceholderText("Enter PDB IDs separated by comma (e.g., 6H0G,6H0F,5FQD)")
        self.pdb_id_input.setMinimumHeight(32)
        pdb_id_row.addWidget(self.pdb_id_input, 1)
        
        self.add_pdb_ids_btn = QPushButton("Add")
        self.add_pdb_ids_btn.setMinimumHeight(32)
        self.add_pdb_ids_btn.clicked.connect(self._add_pdb_ids)
        pdb_id_row.addWidget(self.add_pdb_ids_btn)
        input_layout.addLayout(pdb_id_row)
        
        # File selection
        file_row = QHBoxLayout()
        file_row.addWidget(QLabel("PDB Files:"))
        self.browse_files_btn = QPushButton("Browse Files...")
        self.browse_files_btn.setMinimumHeight(32)
        self.browse_files_btn.clicked.connect(self._browse_files)
        file_row.addWidget(self.browse_files_btn)
        
        self.browse_folder_btn = QPushButton("Browse Folder...")
        self.browse_folder_btn.setMinimumHeight(32)
        self.browse_folder_btn.clicked.connect(self._browse_folder)
        file_row.addWidget(self.browse_folder_btn)
        file_row.addStretch()
        input_layout.addLayout(file_row)
        
        # Selected structures list
        input_layout.addWidget(QLabel("Selected Structures:"))
        self.structure_list = QListWidget()
        self.structure_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.structure_list.setMinimumHeight(120)
        self.structure_list.setMaximumHeight(200)
        input_layout.addWidget(self.structure_list)
        
        # List operation buttons
        list_btn_row = QHBoxLayout()
        self.clear_list_btn = QPushButton("Clear All")
        self.clear_list_btn.clicked.connect(self._clear_list)
        list_btn_row.addWidget(self.clear_list_btn)
        
        self.remove_selected_btn = QPushButton("Remove Selected")
        self.remove_selected_btn.clicked.connect(self._remove_selected)
        list_btn_row.addWidget(self.remove_selected_btn)
        
        self.structure_count_label = QLabel("0 structures")
        list_btn_row.addStretch()
        list_btn_row.addWidget(self.structure_count_label)
        input_layout.addLayout(list_btn_row)
        
        layout.addWidget(grp_input)
        
        # 2. Analysis type selection
        grp_analysis = QGroupBox("Analysis Type")
        analysis_layout = QVBoxLayout(grp_analysis)
        
        type_row = QHBoxLayout()
        type_row.addWidget(QLabel("Analysis:"))
        self.analysis_type_combo = QComboBox()
        self.analysis_type_combo.addItems([
            "G-motif Detection",
            "PPI Interface Analysis",
            "Pocket Detection",
            "Protein-Ligand Interactions",
            "Comprehensive (All)"
        ])
        self.analysis_type_combo.setMinimumHeight(32)
        self.analysis_type_combo.currentIndexChanged.connect(self._on_analysis_type_changed)
        type_row.addWidget(self.analysis_type_combo, 1)
        analysis_layout.addLayout(type_row)
        
        layout.addWidget(grp_analysis)
        
        # 3. Parameter settings (dynamically displayed)
        self.params_group = QGroupBox("Parameters")
        self.params_layout = QVBoxLayout(self.params_group)
        self._create_gmotif_params()  # Default: show G-motif parameters
        layout.addWidget(self.params_group)

        # 4. Output settings
        grp_output = QGroupBox("Output")
        output_layout = QVBoxLayout(grp_output)
        
        output_row = QHBoxLayout()
        output_row.addWidget(QLabel("Output Directory:"))
        self.output_dir_input = QLineEdit()
        self.output_dir_input.setPlaceholderText("Default: ./batch_results")
        self.output_dir_input.setMinimumHeight(32)
        output_row.addWidget(self.output_dir_input, 1)
        
        self.browse_output_btn = QPushButton("Browse...")
        self.browse_output_btn.setMinimumHeight(32)
        self.browse_output_btn.clicked.connect(self._browse_output_dir)
        output_row.addWidget(self.browse_output_btn)
        output_layout.addLayout(output_row)
        
        layout.addWidget(grp_output)
        
        # 5. Run buttons
        run_row = QHBoxLayout()
        self.run_btn = QPushButton("🚀 Run Batch Analysis")
        self.run_btn.setObjectName("primary_btn")
        self.run_btn.setMinimumHeight(48)
        self.run_btn.clicked.connect(self._run_analysis)
        run_row.addWidget(self.run_btn)
        
        self.stop_btn = QPushButton("⏹ Stop")
        self.stop_btn.setMinimumHeight(48)
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self._stop_analysis)
        run_row.addWidget(self.stop_btn)
        layout.addLayout(run_row)
        
        # 6. Progress and logs
        grp_progress = QGroupBox("Progress")
        progress_layout = QVBoxLayout(grp_progress)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimumHeight(24)
        progress_layout.addWidget(self.progress_bar)
        
        self.log_output = QTextEdit()
        self.log_output.setReadOnly(True)
        self.log_output.setMinimumHeight(150)
        self.log_output.setMaximumHeight(250)
        progress_layout.addWidget(self.log_output)
        
        layout.addWidget(grp_progress)
        
        layout.addStretch(1)
    
    def _create_gmotif_params(self):
        """Create G-motif parameter controls"""
        self._clear_params_layout()
        
        row1 = QHBoxLayout()
        row1.addWidget(QLabel("RMSD Cutoff (Å):"))
        self.gmotif_rmsd = QLineEdit("3.5")
        self.gmotif_rmsd.setMaximumWidth(100)
        self.gmotif_rmsd.setMinimumHeight(32)
        row1.addWidget(self.gmotif_rmsd)
        
        row1.addWidget(QLabel("Template:"))
        self.gmotif_template = QComboBox()
        self.gmotif_template.addItems(["GSPT1", "CK1α", "VAV1"])
        self.gmotif_template.setMinimumHeight(32)
        row1.addWidget(self.gmotif_template)
        
        self.gmotif_require_gly = QCheckBox("Require Central Gly")
        self.gmotif_require_gly.setChecked(True)
        row1.addWidget(self.gmotif_require_gly)
        row1.addStretch()
        
        self.params_layout.addLayout(row1)
    
    def _create_ppi_params(self):
        """Create PPI parameter controls"""
        self._clear_params_layout()
        
        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Protein 1 Chains:"))
        self.ppi_chains1 = QLineEdit("A")
        self.ppi_chains1.setPlaceholderText("e.g., A or A,C")
        self.ppi_chains1.setMaximumWidth(150)
        self.ppi_chains1.setMinimumHeight(32)
        row1.addWidget(self.ppi_chains1)
        
        row1.addWidget(QLabel("Protein 2 Chains:"))
        self.ppi_chains2 = QLineEdit("B")
        self.ppi_chains2.setPlaceholderText("e.g., B or B,D")
        self.ppi_chains2.setMaximumWidth(150)
        self.ppi_chains2.setMinimumHeight(32)
        row1.addWidget(self.ppi_chains2)
        
        row1.addWidget(QLabel("Interface Distance (Å):"))
        self.ppi_distance = QLineEdit("4.5")
        self.ppi_distance.setMaximumWidth(80)
        self.ppi_distance.setMinimumHeight(32)
        row1.addWidget(self.ppi_distance)
        row1.addStretch()
        
        self.params_layout.addLayout(row1)
    
    def _create_pocket_params(self):
        """Create pocket detection parameter controls"""
        self._clear_params_layout()
        
        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Min Volume (ų):"))
        self.pocket_min_volume = QLineEdit("30.0")
        self.pocket_min_volume.setMaximumWidth(100)
        self.pocket_min_volume.setMinimumHeight(32)
        row1.addWidget(self.pocket_min_volume)
        
        row1.addWidget(QLabel("Min Depth (Å):"))
        self.pocket_min_depth = QLineEdit("2.5")
        self.pocket_min_depth.setMaximumWidth(100)
        self.pocket_min_depth.setMinimumHeight(32)
        row1.addWidget(self.pocket_min_depth)
        row1.addStretch()
        
        self.params_layout.addLayout(row1)
    
    def _create_interaction_params(self):
        """Create interaction analysis parameter controls"""
        self._clear_params_layout()
        
        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Ligand Names:"))
        self.interaction_ligands = QLineEdit()
        self.interaction_ligands.setPlaceholderText("Optional: e.g., CC9,LEN (auto-detect if empty)")
        self.interaction_ligands.setMinimumHeight(32)
        row1.addWidget(self.interaction_ligands, 1)
        
        row1.addWidget(QLabel("Distance (Å):"))
        self.interaction_distance = QLineEdit("4.5")
        self.interaction_distance.setMaximumWidth(80)
        self.interaction_distance.setMinimumHeight(32)
        row1.addWidget(self.interaction_distance)
        
        self.params_layout.addLayout(row1)
    
    def _create_comprehensive_params(self):
        """Create comprehensive analysis parameter controls"""
        self._clear_params_layout()

        # PPI parameters
        row1 = QHBoxLayout()
        row1.addWidget(QLabel("PPI Chains 1:"))
        self.comp_chains1 = QLineEdit("A")
        self.comp_chains1.setMaximumWidth(100)
        self.comp_chains1.setMinimumHeight(32)
        row1.addWidget(self.comp_chains1)
        
        row1.addWidget(QLabel("PPI Chains 2:"))
        self.comp_chains2 = QLineEdit("B")
        self.comp_chains2.setMaximumWidth(100)
        self.comp_chains2.setMinimumHeight(32)
        row1.addWidget(self.comp_chains2)
        
        row1.addWidget(QLabel("Ligands:"))
        self.comp_ligands = QLineEdit()
        self.comp_ligands.setPlaceholderText("Optional")
        self.comp_ligands.setMaximumWidth(150)
        self.comp_ligands.setMinimumHeight(32)
        row1.addWidget(self.comp_ligands)
        row1.addStretch()
        
        self.params_layout.addLayout(row1)
        
        note = QLabel("Note: Comprehensive analysis runs all analysis types and generates a summary report.")
        note.setStyleSheet("color: gray; font-style: italic;")
        self.params_layout.addWidget(note)
    
    def _clear_params_layout(self):
        """Clear all widgets from parameter layout"""
        while self.params_layout.count():
            item = self.params_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
            elif item.layout():
                self._clear_layout(item.layout())
    
    def _clear_layout(self, layout):
        """Recursively clear layout"""
        while layout.count():
            item = layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
            elif item.layout():
                self._clear_layout(item.layout())
    
    def _on_analysis_type_changed(self, index):
        """Update parameter panel when analysis type changes"""
        if index == 0:  # G-motif
            self._create_gmotif_params()
        elif index == 1:  # PPI
            self._create_ppi_params()
        elif index == 2:  # Pocket
            self._create_pocket_params()
        elif index == 3:  # Interactions
            self._create_interaction_params()
        elif index == 4:  # Comprehensive
            self._create_comprehensive_params()
    
    def _add_pdb_ids(self):
        """Add PDB IDs to list"""
        text = self.pdb_id_input.text().strip()
        if not text:
            return
        
        pdb_ids = [p.strip().upper() for p in text.split(",") if p.strip()]
        for pdb_id in pdb_ids:
            if len(pdb_id) == 4 and pdb_id.isalnum():
                if pdb_id not in self._pdb_sources:
                    self._pdb_sources.append(pdb_id)
                    item = QListWidgetItem(f"📦 {pdb_id} (PDB)")
                    self.structure_list.addItem(item)
        
        self.pdb_id_input.clear()
        self._update_count()
    
    def _browse_files(self):
        """Browse and select PDB files"""
        files, _ = QFileDialog.getOpenFileNames(
            self, "Select PDB Files", "",
            "PDB Files (*.pdb *.cif *.ent);;All Files (*)"
        )
        
        for f in files:
            if f not in self._pdb_sources:
                self._pdb_sources.append(f)
                item = QListWidgetItem(f"📄 {os.path.basename(f)}")
                item.setToolTip(f)
                self.structure_list.addItem(item)
        
        self._update_count()
    
    def _browse_folder(self):
        """Browse folder and add all PDB files"""
        folder = QFileDialog.getExistingDirectory(self, "Select Folder with PDB Files")
        if not folder:
            return
        
        count = 0
        for f in os.listdir(folder):
            if f.lower().endswith(('.pdb', '.cif', '.ent')):
                full_path = os.path.join(folder, f)
                if full_path not in self._pdb_sources:
                    self._pdb_sources.append(full_path)
                    item = QListWidgetItem(f"📄 {f}")
                    item.setToolTip(full_path)
                    self.structure_list.addItem(item)
                    count += 1
        
        if count > 0:
            self.log(f"Added {count} PDB files from {folder}")
        self._update_count()
    
    def _clear_list(self):
        """Clear list"""
        self._pdb_sources.clear()
        self.structure_list.clear()
        self._update_count()
    
    def _remove_selected(self):
        """Remove selected items"""
        for item in self.structure_list.selectedItems():
            row = self.structure_list.row(item)
            self.structure_list.takeItem(row)
            if row < len(self._pdb_sources):
                self._pdb_sources.pop(row)
        self._update_count()
    
    def _update_count(self):
        """Update structure count"""
        count = len(self._pdb_sources)
        self.structure_count_label.setText(f"{count} structure{'s' if count != 1 else ''}")
    
    def _browse_output_dir(self):
        """Browse output directory"""
        folder = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if folder:
            self.output_dir_input.setText(folder)
    
    def _run_analysis(self):
        """Run batch analysis"""
        if not self._pdb_sources:
            QMessageBox.warning(self, "Warning", "Please add at least one structure to analyze.")
            return
        
        # Get analysis type
        analysis_types = ["gmotif", "ppi", "pockets", "interactions", "comprehensive"]
        analysis_type = analysis_types[self.analysis_type_combo.currentIndex()]
        
        # Collect parameters
        params = {
            "output_dir": self.output_dir_input.text().strip() or None,
            "output_csv": None
        }
        
        if analysis_type == "gmotif":
            try:
                params["rmsd_cutoff"] = float(self.gmotif_rmsd.text())
            except:
                params["rmsd_cutoff"] = 3.5
            params["require_gly"] = self.gmotif_require_gly.isChecked()
            params["template"] = self.gmotif_template.currentText()
            
        elif analysis_type == "ppi":
            chains1 = [c.strip() for c in self.ppi_chains1.text().split(",") if c.strip()]
            chains2 = [c.strip() for c in self.ppi_chains2.text().split(",") if c.strip()]
            if not chains1 or not chains2:
                QMessageBox.warning(self, "Warning", "Please specify chain IDs for both proteins.")
                return
            params["chain_pairs"] = [(chains1, chains2)]
            try:
                params["interface_distance"] = float(self.ppi_distance.text())
            except:
                params["interface_distance"] = 4.5
                
        elif analysis_type == "pockets":
            try:
                params["min_volume"] = float(self.pocket_min_volume.text())
            except:
                params["min_volume"] = 30.0
            try:
                params["min_depth"] = float(self.pocket_min_depth.text())
            except:
                params["min_depth"] = 2.5
                
        elif analysis_type == "interactions":
            ligands = self.interaction_ligands.text().strip()
            params["ligand_resnames"] = [l.strip() for l in ligands.split(",") if l.strip()] if ligands else None
            try:
                params["distance_cutoff"] = float(self.interaction_distance.text())
            except:
                params["distance_cutoff"] = 4.5
                
        elif analysis_type == "comprehensive":
            chains1 = [c.strip() for c in self.comp_chains1.text().split(",") if c.strip()]
            chains2 = [c.strip() for c in self.comp_chains2.text().split(",") if c.strip()]
            params["chain_pairs"] = [(chains1, chains2)] if chains1 and chains2 else None
            ligands = self.comp_ligands.text().strip()
            params["ligand_resnames"] = [l.strip() for l in ligands.split(",") if l.strip()] if ligands else None
        
        # Set output CSV
        if params["output_dir"]:
            params["output_csv"] = os.path.join(params["output_dir"], f"{analysis_type}_results.csv")
        
        # Update UI
        self.run_btn.setEnabled(False)
        self.stop_btn.setEnabled(True)
        self.progress_bar.setRange(0, 0)  # Indeterminate
        self.log_output.clear()
        self.log(f"Starting {analysis_type} analysis on {len(self._pdb_sources)} structures...")
        
        # Start worker thread
        self._batch_thread = BatchWorker(analysis_type, self._pdb_sources.copy(), params)
        self._batch_thread.progress.connect(self.log)
        self._batch_thread.error.connect(self._on_error)
        self._batch_thread.finished.connect(self._on_finished)
        self._batch_thread.start()
    
    def _stop_analysis(self):
        """Stop analysis"""
        if self._batch_thread and self._batch_thread.isRunning():
            self._batch_thread.terminate()
            self._batch_thread.wait()
            self.log("Analysis stopped by user.")
            self._reset_ui()
    
    def _on_error(self, error_msg: str):
        """Handle error"""
        self.log(f"❌ Error: {error_msg}")
        QMessageBox.critical(self, "Error", error_msg)
        self._reset_ui()
    
    def _on_finished(self, result: dict):
        """Analysis complete"""
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(1)
        
        if result.get("success", True):
            total = result.get("total_structures", len(self._pdb_sources))
            successful = result.get("successful", 0)
            self.log(f"\n✅ Analysis complete!")
            self.log(f"   Total structures: {total}")
            self.log(f"   Successful: {successful}")
            
            # Show specific results
            if "total_hits" in result:
                self.log(f"   Total G-motifs: {result['total_hits']}")
            if "total_pockets" in result:
                self.log(f"   Total pockets: {result['total_pockets']}")
            if "avg_contacts" in result:
                self.log(f"   Avg contacts: {result['avg_contacts']:.1f}")
            
            QMessageBox.information(self, "Complete", 
                f"Batch analysis complete!\n\n"
                f"Analyzed: {total} structures\n"
                f"Successful: {successful}")
        else:
            self.log(f"❌ Analysis failed: {result.get('error', 'Unknown error')}")
        
        self._reset_ui()
    
    def _reset_ui(self):
        """Reset UI state"""
        self.run_btn.setEnabled(True)
        self.stop_btn.setEnabled(False)
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)