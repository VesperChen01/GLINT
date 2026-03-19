# -*- coding: utf-8 -*-
"""
Hit Identification Tab: Vina Docking, HADDOCK3, Mutation Analysis
"""
import os
from typing import Optional

from ..qt_adapter import (
    Qt, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QGroupBox, QScrollArea, QFrame, QFileDialog, QMessageBox,
    QComboBox, QCheckBox, QSpinBox, QGridLayout
)

from ..utils import t, show_message_box
from .common import CommonTab

class HitIdentificationTab(CommonTab):
    def __init__(self, parent):
        super().__init__(parent)
        self.init_ui()

    def init_ui(self):
        """InitializeUI - 现代卡片式布局"""
        self.setObjectName("scroll_content")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.parent_window._hit_scroll_content = self

        is_dark = getattr(self.parent_window, "_dark_mode", False)
        bg_color = "#161b22" if is_dark else "#f8fafc"
        self.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")

        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(16)
        main_layout.setContentsMargins(20, 20, 20, 20)

        # === 页面标题 ===
        header = QHBoxLayout()
        title = QLabel("Hit Identification")
        title.setStyleSheet("""
            font-size: 20px; font-weight: 600;
            color: #3b82f6; padding: 4px 0;
        """)
        header.addWidget(title)
        header.addStretch(1)
        main_layout.addLayout(header)

        # 1. Vina Docking (Small Molecule)
        main_layout.addWidget(self._create_vina_card(is_dark))
        
        # Vina buttons row
        vina_btn_row = QHBoxLayout()
        vina_btn_row.setSpacing(10)

        self.parent_window.vina_dock_btn = QPushButton("Run Docking")
        self.parent_window.vina_dock_btn.setMinimumHeight(36)
        self.parent_window.vina_dock_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.vina_dock_btn.clicked.connect(self.run_vina_docking)

        self.parent_window.vina_load_result_btn = QPushButton("Load Result")
        self.parent_window.vina_load_result_btn.setMinimumHeight(36)
        self.parent_window.vina_load_result_btn.setStyleSheet(self._get_secondary_btn_style())
        self.parent_window.vina_load_result_btn.clicked.connect(self.load_vina_result)
        
        vina_btn_row.addWidget(self.parent_window.vina_dock_btn)
        vina_btn_row.addWidget(self.parent_window.vina_load_result_btn)
        vina_btn_row.addStretch(1)
        main_layout.addLayout(vina_btn_row)

        # 2. HADDOCK3 (Protein-Protein Docking)
        main_layout.addWidget(self._create_hdock_card(is_dark))

        # HADDOCK3 buttons row
        hdock_btn_row = QHBoxLayout()
        hdock_btn_row.setSpacing(10)

        self.parent_window.hdock_run_btn = QPushButton("Run HADDOCK3")
        self.parent_window.hdock_run_btn.setMinimumHeight(36)
        self.parent_window.hdock_run_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.hdock_run_btn.clicked.connect(self.run_hdock)

        hdock_btn_row.addWidget(self.parent_window.hdock_run_btn)
        hdock_btn_row.addStretch(1)
        main_layout.addLayout(hdock_btn_row)

        # 3. Protein Mutation & ΔΔG Analysis (independent card under HADDOCK3)
        main_layout.addWidget(self._create_mutation_analysis_card(is_dark))

        main_layout.addStretch(1)

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

    def _create_vina_card(self, is_dark: bool = False) -> QWidget:
        """Create compact Vina docking card"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(12)
        layout.setContentsMargins(16, 14, 16, 14)

        title = QLabel("AutoDock Vina - Small Molecule Docking")
        title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        layout.addWidget(title)
        
        # Row 1: Receptor + Ligand
        row1 = QHBoxLayout()
        row1.setSpacing(6)
        
        self.parent_window.vina_receptor_combo = QComboBox()
        self.parent_window.vina_receptor_combo.setMinimumHeight(28)
        self.parent_window.vina_receptor_combo.setMinimumWidth(120)
        refresh_btn = QPushButton("Refresh")
        refresh_btn.setMinimumHeight(28)
        refresh_btn.clicked.connect(self.refresh_objects)
        
        self.parent_window.vina_ligand = QLineEdit()
        self.parent_window.vina_ligand.setPlaceholderText("Ligand file or folder")
        self.parent_window.vina_ligand.setMinimumHeight(28)
        self.parent_window.vina_ligand.setMinimumWidth(150)
        
        lig_browse_btn = QPushButton("Browse")
        lig_browse_btn.setMinimumHeight(28)
        lig_browse_btn.setToolTip("Select ligand file or folder for batch docking")
        lig_browse_btn.clicked.connect(self.browse_vina_ligand)
        
        row1.addWidget(QLabel("Receptor:"))
        row1.addWidget(self.parent_window.vina_receptor_combo)
        row1.addWidget(refresh_btn)
        row1.addSpacing(10)
        row1.addWidget(QLabel("Ligand:"))
        row1.addWidget(self.parent_window.vina_ligand)
        row1.addWidget(lig_browse_btn)
        row1.addStretch()
        layout.addLayout(row1)
        
        # Row 2: Selection + Get Center + Center/Size/Exhaust/Modes
        row2 = QHBoxLayout()
        row2.setSpacing(6)
        
        self.parent_window.vina_selection = QLineEdit()
        self.parent_window.vina_selection.setPlaceholderText("resn LIG")
        self.parent_window.vina_selection.setMinimumHeight(28)
        self.parent_window.vina_selection.setMinimumWidth(100)
        
        get_center_btn = QPushButton("Get Center")
        get_center_btn.setMinimumHeight(28)
        get_center_btn.clicked.connect(self.get_vina_center)
        
        self.parent_window.vina_center_x = QLineEdit("0")
        self.parent_window.vina_center_x.setMinimumHeight(28)
        self.parent_window.vina_center_x.setMaximumWidth(50)
        self.parent_window.vina_center_y = QLineEdit("0")
        self.parent_window.vina_center_y.setMinimumHeight(28)
        self.parent_window.vina_center_y.setMaximumWidth(50)
        self.parent_window.vina_center_z = QLineEdit("0")
        self.parent_window.vina_center_z.setMinimumHeight(28)
        self.parent_window.vina_center_z.setMaximumWidth(50)
        
        self.parent_window.vina_size_x = QLineEdit("20")
        self.parent_window.vina_size_x.setMinimumHeight(28)
        self.parent_window.vina_size_x.setMaximumWidth(50)
        self.parent_window.vina_size_y = QLineEdit("20")
        self.parent_window.vina_size_y.setMinimumHeight(28)
        self.parent_window.vina_size_y.setMaximumWidth(50)
        self.parent_window.vina_size_z = QLineEdit("20")
        self.parent_window.vina_size_z.setMinimumHeight(28)
        self.parent_window.vina_size_z.setMaximumWidth(50)
        
        self.parent_window.vina_exhaustiveness = QSpinBox()
        self.parent_window.vina_exhaustiveness.setRange(1, 32)
        self.parent_window.vina_exhaustiveness.setValue(8)
        self.parent_window.vina_exhaustiveness.setMinimumHeight(28)
        self.parent_window.vina_exhaustiveness.setMaximumWidth(60)
        
        self.parent_window.vina_num_modes = QSpinBox()
        self.parent_window.vina_num_modes.setRange(1, 20)
        self.parent_window.vina_num_modes.setValue(9)
        self.parent_window.vina_num_modes.setMinimumHeight(28)
        self.parent_window.vina_num_modes.setMaximumWidth(60)
        
        row2.addWidget(QLabel("Selection:"))
        row2.addWidget(self.parent_window.vina_selection)
        row2.addWidget(get_center_btn)
        row2.addSpacing(6)
        row2.addWidget(QLabel("Center:"))
        row2.addWidget(self.parent_window.vina_center_x)
        row2.addWidget(self.parent_window.vina_center_y)
        row2.addWidget(self.parent_window.vina_center_z)
        row2.addSpacing(6)
        row2.addWidget(QLabel("Size:"))
        row2.addWidget(self.parent_window.vina_size_x)
        row2.addWidget(self.parent_window.vina_size_y)
        row2.addWidget(self.parent_window.vina_size_z)
        row2.addSpacing(6)
        row2.addWidget(QLabel("Exhaust:"))
        row2.addWidget(self.parent_window.vina_exhaustiveness)
        row2.addWidget(QLabel("Modes:"))
        row2.addWidget(self.parent_window.vina_num_modes)
        row2.addStretch()
        layout.addLayout(row2)
        
        return card

    def _create_hdock_card(self, is_dark: bool = False) -> QWidget:
        """Create compact HADDOCK3 card"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(12)
        layout.setContentsMargins(16, 14, 16, 14)

        title = QLabel("HADDOCK3 - Protein-Protein Docking")
        title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        layout.addWidget(title)
        
        # Row 1: Receptor + Ligand
        row1 = QHBoxLayout()
        row1.setSpacing(6)
        
        self.parent_window.hdock_receptor_combo = QComboBox()
        self.parent_window.hdock_receptor_combo.setMinimumHeight(28)
        self.parent_window.hdock_receptor_combo.setMinimumWidth(120)
        
        refresh_btn = QPushButton("Refresh")
        refresh_btn.setMinimumHeight(28)
        refresh_btn.clicked.connect(self.refresh_objects)
        
        self.parent_window.hdock_ligand_combo = QComboBox()
        self.parent_window.hdock_ligand_combo.setMinimumHeight(28)
        self.parent_window.hdock_ligand_combo.setMinimumWidth(120)
        
        row1.addWidget(QLabel("Receptor:"))
        row1.addWidget(self.parent_window.hdock_receptor_combo)
        row1.addWidget(refresh_btn)
        row1.addSpacing(10)
        row1.addWidget(QLabel("Ligand:"))
        row1.addWidget(self.parent_window.hdock_ligand_combo)
        row1.addStretch()
        layout.addLayout(row1)
        
        # Row 2: Active/Passive residues
        row2 = QHBoxLayout()
        row2.setSpacing(6)
        
        self.parent_window.hdock_active_receptor = QLineEdit()
        self.parent_window.hdock_active_receptor.setPlaceholderText("e.g., 10,15,20")
        self.parent_window.hdock_active_receptor.setMinimumHeight(28)
        self.parent_window.hdock_active_receptor.setMinimumWidth(100)
        
        self.parent_window.hdock_passive_receptor = QLineEdit()
        self.parent_window.hdock_passive_receptor.setPlaceholderText("e.g., 11,16,21")
        self.parent_window.hdock_passive_receptor.setMinimumHeight(28)
        self.parent_window.hdock_passive_receptor.setMinimumWidth(100)
        
        self.parent_window.hdock_active_ligand = QLineEdit()
        self.parent_window.hdock_active_ligand.setPlaceholderText("e.g., 5,8,12")
        self.parent_window.hdock_active_ligand.setMinimumHeight(28)
        self.parent_window.hdock_active_ligand.setMinimumWidth(100)
        
        self.parent_window.hdock_passive_ligand = QLineEdit()
        self.parent_window.hdock_passive_ligand.setPlaceholderText("e.g., 6,9,13")
        self.parent_window.hdock_passive_ligand.setMinimumHeight(28)
        self.parent_window.hdock_passive_ligand.setMinimumWidth(100)
        
        row2.addWidget(QLabel("Receptor Active:"))
        row2.addWidget(self.parent_window.hdock_active_receptor)
        row2.addWidget(QLabel("Passive:"))
        row2.addWidget(self.parent_window.hdock_passive_receptor)
        row2.addSpacing(10)
        row2.addWidget(QLabel("Ligand Active:"))
        row2.addWidget(self.parent_window.hdock_active_ligand)
        row2.addWidget(QLabel("Passive:"))
        row2.addWidget(self.parent_window.hdock_passive_ligand)
        row2.addStretch()
        layout.addLayout(row2)
        
        return card

    def _create_mutation_analysis_card(self, is_dark: bool = False) -> QWidget:
        """Create compact mutation analysis card"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(12)
        layout.setContentsMargins(16, 14, 16, 14)

        title = QLabel("Protein Mutation & ΔΔG Analysis")
        title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        layout.addWidget(title)
        
        # Row 1: Structure + Mutation List
        row1 = QHBoxLayout()
        row1.setSpacing(6)
        
        self.parent_window.mutation_structure_combo = QComboBox()
        self.parent_window.mutation_structure_combo.setMinimumHeight(28)
        self.parent_window.mutation_structure_combo.setMinimumWidth(120)
        
        refresh_btn = QPushButton("Refresh")
        refresh_btn.setMinimumHeight(28)
        refresh_btn.clicked.connect(self.refresh_objects)
        
        self.parent_window.mutation_list = QLineEdit()
        self.parent_window.mutation_list.setPlaceholderText("e.g., A50G,A51V,A52L")
        self.parent_window.mutation_list.setMinimumHeight(28)
        self.parent_window.mutation_list.setMinimumWidth(200)
        
        row1.addWidget(QLabel("Structure:"))
        row1.addWidget(self.parent_window.mutation_structure_combo)
        row1.addWidget(refresh_btn)
        row1.addSpacing(10)
        row1.addWidget(QLabel("Mutations:"))
        row1.addWidget(self.parent_window.mutation_list)
        row1.addStretch()
        layout.addLayout(row1)
        
        # Row 2: Run button
        row2 = QHBoxLayout()
        row2.setSpacing(10)
        
        self.parent_window.mutation_run_btn = QPushButton("Run ΔΔG Analysis")
        self.parent_window.mutation_run_btn.setMinimumHeight(36)
        self.parent_window.mutation_run_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.mutation_run_btn.clicked.connect(self.run_mutation_analysis)
        
        row2.addWidget(self.parent_window.mutation_run_btn)
        row2.addStretch()
        layout.addLayout(row2)
        
        return card

    def _get_primary_btn_style(self) -> str:
        """获取主按钮样式"""
        return """
            QPushButton {
                background: #3b82f6;
                color: white;
                border: none;
                border-radius: 6px;
                padding: 8px 16px;
                font-weight: 500;
            }
            QPushButton:hover {
                background: #2563eb;
            }
            QPushButton:pressed {
                background: #1d4ed8;
            }
            QPushButton:disabled {
                background: #94a3b8;
            }
        """

    def _get_secondary_btn_style(self) -> str:
        """获取次要按钮样式"""
        return """
            QPushButton {
                background: #64748b;
                color: white;
                border: none;
                border-radius: 6px;
                padding: 8px 16px;
                font-weight: 500;
            }
            QPushButton:hover {
                background: #475569;
            }
            QPushButton:pressed {
                background: #334155;
            }
            QPushButton:disabled {
                background: #cbd5e1;
            }
        """

    def refresh_objects(self):
        """刷新PyMOL对象列表"""
        try:
            from pymol import cmd
            objects = cmd.get_object_list()
            
            # Update Vina receptor
            current_vina = self.parent_window.vina_receptor_combo.currentText()
            self.parent_window.vina_receptor_combo.clear()
            self.parent_window.vina_receptor_combo.addItems(objects)
            if current_vina in objects:
                self.parent_window.vina_receptor_combo.setCurrentText(current_vina)
            
            # Update HADDOCK3 receptor and ligand
            current_hdock_rec = self.parent_window.hdock_receptor_combo.currentText()
            current_hdock_lig = self.parent_window.hdock_ligand_combo.currentText()
            self.parent_window.hdock_receptor_combo.clear()
            self.parent_window.hdock_receptor_combo.addItems(objects)
            self.parent_window.hdock_ligand_combo.clear()
            self.parent_window.hdock_ligand_combo.addItems(objects)
            if current_hdock_rec in objects:
                self.parent_window.hdock_receptor_combo.setCurrentText(current_hdock_rec)
            if current_hdock_lig in objects:
                self.parent_window.hdock_ligand_combo.setCurrentText(current_hdock_lig)
            
            # Update mutation structure
            current_mutation = self.parent_window.mutation_structure_combo.currentText()
            self.parent_window.mutation_structure_combo.clear()
            self.parent_window.mutation_structure_combo.addItems(objects)
            if current_mutation in objects:
                self.parent_window.mutation_structure_combo.setCurrentText(current_mutation)
            
            self.log(f"✓ Refreshed objects: {len(objects)} found")
        except Exception as e:
            self.log(f"❌ Failed to refresh objects: {str(e)}")

    def browse_vina_ligand(self):
        """Browse for ligand file or folder"""
        from ..qt_adapter import QFileDialog
        
        # Ask user to choose file or folder
        choice = show_message_box(
            self,
            "Select Ligand Input",
            "Choose input type:\n\n"
            "• Single File: Select one ligand file (.pdbqt, .sdf, .mol2)\n"
            "• Folder: Select folder containing multiple ligands for batch docking",
            "question",
            buttons=["Single File", "Folder", "Cancel"]
        )
        
        if choice == "Single File":
            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Select Ligand File",
                "",
                "Ligand Files (*.pdbqt *.sdf *.mol2);;All Files (*)"
            )
            if file_path:
                self.parent_window.vina_ligand.setText(file_path)
                self.log(f"✓ Selected ligand file: {file_path}")
        
        elif choice == "Folder":
            folder_path = QFileDialog.getExistingDirectory(
                self,
                "Select Ligand Folder"
            )
            if folder_path:
                self.parent_window.vina_ligand.setText(folder_path)
                self.log(f"✓ Selected ligand folder: {folder_path}")

    def get_vina_center(self):
        """Get center coordinates from PyMOL selection"""
        try:
            from pymol import cmd
            selection = self.parent_window.vina_selection.text().strip()
            if not selection:
                show_message_box(self, "Error", "Please enter a selection", "warning")
                return
            
            # Get center of mass
            model = cmd.get_model(selection)
            if not model.atom:
                show_message_box(self, "Error", f"No atoms found in selection: {selection}", "warning")
                return
            
            x = sum(atom.coord[0] for atom in model.atom) / len(model.atom)
            y = sum(atom.coord[1] for atom in model.atom) / len(model.atom)
            z = sum(atom.coord[2] for atom in model.atom) / len(model.atom)
            
            self.parent_window.vina_center_x.setText(f"{x:.2f}")
            self.parent_window.vina_center_y.setText(f"{y:.2f}")
            self.parent_window.vina_center_z.setText(f"{z:.2f}")
            
            self.log(f"✓ Center calculated: ({x:.2f}, {y:.2f}, {z:.2f})")
        except Exception as e:
            show_message_box(self, "Error", f"Failed to get center: {str(e)}", "critical")
            self.log(f"❌ Failed to get center: {str(e)}")

    def run_vina_docking(self):
        """Run AutoDock Vina docking"""
        try:
            from pymol import cmd
            import subprocess
            import tempfile
            import shutil
            
            # Get parameters
            receptor_obj = self.parent_window.vina_receptor_combo.currentText()
            ligand_path = self.parent_window.vina_ligand.text().strip()
            
            if not receptor_obj or not ligand_path:
                show_message_box(self, "Error", "Please select receptor and ligand", "warning")
                return
            
            if not os.path.exists(ligand_path):
                show_message_box(self, "Error", f"Ligand path not found: {ligand_path}", "warning")
                return
            
            # Get box parameters
            try:
                center_x = float(self.parent_window.vina_center_x.text())
                center_y = float(self.parent_window.vina_center_y.text())
                center_z = float(self.parent_window.vina_center_z.text())
                size_x = float(self.parent_window.vina_size_x.text())
                size_y = float(self.parent_window.vina_size_y.text())
                size_z = float(self.parent_window.vina_size_z.text())
                exhaustiveness = int(self.parent_window.vina_exhaustiveness.value())
                num_modes = int(self.parent_window.vina_num_modes.value())
            except ValueError as e:
                show_message_box(self, "Error", f"Invalid numeric parameter: {str(e)}", "warning")
                return
            
            self.log("🚀 Starting Vina docking...")
            self.log(f"   Receptor: {receptor_obj}")
            self.log(f"   Ligand: {ligand_path}")
            self.log(f"   Box center: ({center_x}, {center_y}, {center_z})")
            self.log(f"   Box size: ({size_x}, {size_y}, {size_z})")
            
            # Create temp directory
            temp_dir = tempfile.mkdtemp(prefix="vina_")
            receptor_pdbqt = os.path.join(temp_dir, "receptor.pdbqt")
            
            # Save receptor as PDBQT
            cmd.save(receptor_pdbqt, receptor_obj)
            
            # Check if ligand is file or folder
            is_batch = os.path.isdir(ligand_path)
            
            if is_batch:
                # Batch docking
                ligand_files = [f for f in os.listdir(ligand_path) 
                              if f.endswith(('.pdbqt', '.sdf', '.mol2'))]
                if not ligand_files:
                    show_message_box(self, "Error", "No ligand files found in folder", "warning")
                    shutil.rmtree(temp_dir)
                    return
                
                self.log(f"📦 Batch docking: {len(ligand_files)} ligands")
                results = []
                
                for i, lig_file in enumerate(ligand_files, 1):
                    lig_path = os.path.join(ligand_path, lig_file)
                    output_pdbqt = os.path.join(temp_dir, f"output_{i}.pdbqt")
                    
                    self.log(f"   [{i}/{len(ligand_files)}] Docking {lig_file}...")
                    
                    # Run Vina
                    cmd_list = [
                        "vina",
                        "--receptor", receptor_pdbqt,
                        "--ligand", lig_path,
                        "--out", output_pdbqt,
                        "--center_x", str(center_x),
                        "--center_y", str(center_y),
                        "--center_z", str(center_z),
                        "--size_x", str(size_x),
                        "--size_y", str(size_y),
                        "--size_z", str(size_z),
                        "--exhaustiveness", str(exhaustiveness),
                        "--num_modes", str(num_modes)
                    ]
                    
                    result = subprocess.run(cmd_list, capture_output=True, text=True)
                    
                    if result.returncode == 0 and os.path.exists(output_pdbqt):
                        # Parse affinity from output
                        affinity = self._parse_vina_affinity(result.stdout)
                        results.append((lig_file, affinity, output_pdbqt))
                        self.log(f"      ✓ Affinity: {affinity} kcal/mol")
                    else:
                        self.log(f"      ❌ Failed: {result.stderr}")
                
                # Sort by affinity and load top results
                results.sort(key=lambda x: x[1])
                top_n = min(5, len(results))
                
                self.log(f"\n📊 Top {top_n} results:")
                for i, (lig_file, affinity, output_path) in enumerate(results[:top_n], 1):
                    obj_name = f"vina_result_{i}"
                    cmd.load(output_path, obj_name)
                    self.log(f"   {i}. {lig_file}: {affinity} kcal/mol → {obj_name}")
                
                self.log(f"✅ Batch docking completed: {len(results)} ligands processed")
                
            else:
                # Single ligand docking
                output_pdbqt = os.path.join(temp_dir, "output.pdbqt")

                self.log("📦 Single ligand docking...")

                cmd_list = [
                    "vina",
                    "--receptor", receptor_pdbqt,
                    "--ligand", ligand_path,
                    "--out", output_pdbqt,
                    "--center_x", str(center_x),
                    "--center_y", str(center_y),
                    "--center_z", str(center_z),
                    "--size_x", str(size_x),
                    "--size_y", str(size_y),
                    "--size_z", str(size_z),
                    "--exhaustiveness", str(exhaustiveness),
                    "--num_modes", str(num_modes)
                ]

                result = subprocess.run(cmd_list, capture_output=True, text=True)

                if result.returncode == 0 and os.path.exists(output_pdbqt):
                    affinity = self._parse_vina_affinity(result.stdout)
                    obj_name = "vina_result"
                    cmd.load(output_pdbqt, obj_name)
                    self.log(f"✅ Docking completed: Affinity {affinity} kcal/mol → {obj_name}")
                else:
                    error_msg = result.stderr.strip() or result.stdout.strip() or "Unknown error"
                    self.log(f"❌ Docking failed: {error_msg}")

        except Exception as e:
            show_message_box(self, "Error", f"Vina docking failed: {str(e)}", "critical")
            self.log(f"❌ Vina docking failed: {str(e)}")
        finally:
            try:
                import shutil as _shutil
                if 'temp_dir' in locals() and temp_dir and os.path.exists(temp_dir):
                    _shutil.rmtree(temp_dir)
            except Exception as e:
                self.log(f"⚠️ Failed to clean temp dir: {str(e)}")

    def load_vina_result(self):
        """Load Vina docking result"""
        try:
            from pymol import cmd

            file_path, _ = QFileDialog.getOpenFileName(
                self,
                "Select Vina Result",
                "",
                "Vina Results (*.pdbqt *.pdb);;All Files (*)"
            )

            if not file_path:
                return

            if not os.path.exists(file_path):
                show_message_box(self, "Error", f"File not found: {file_path}", "warning")
                return

            obj_name = os.path.splitext(os.path.basename(file_path))[0] or "vina_result"
            cmd.load(file_path, obj_name)
            self.log(f"✓ Loaded Vina result: {file_path} → {obj_name}")
        except Exception as e:
            show_message_box(self, "Error", f"Failed to load Vina result: {str(e)}", "critical")
            self.log(f"❌ Failed to load Vina result: {str(e)}")


    def _parse_vina_affinity(self, output: str) -> Optional[float]:
        """Parse affinity from Vina output"""
        if not output:
            return None

        for line in output.splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit() and parts[0] == "1":
                try:
                    return float(parts[1])
                except ValueError:
                    continue

        self.log("未能从 Vina 输出解析出亲和力数值")

    def _parse_mutation_list(self, mutation_text: str):
        """Parse mutation list like A50G,A51V into tuples"""
        import re

        raw = mutation_text.replace(";", ",")
        items = [item.strip() for item in raw.split(",") if item.strip()]
        if not items:
            raise ValueError("突变列表为空，请输入如 A50G,A51V")

        mutations = []
        pattern = re.compile(r"^([A-Za-z])(\d+)([A-Za-z]{1,3})$")
        for item in items:
            match = pattern.match(item)
            if not match:
                raise ValueError(f"突变格式错误：{item}（示例：A50G）")
            chain, resi, target = match.groups()
            mutations.append((chain.upper(), resi, target.upper()))
        return mutations

    def run_mutation_analysis(self):
        """Run mutation analysis using mutation_analyzer"""
        try:
            from pymol import cmd
            from ...mutation_analyzer import analyze_mutation_effects
        except Exception as e:
            show_message_box(self, "Error", f"Mutation analyzer not available: {str(e)}", "critical")
            self.log(f"❌ Mutation analyzer not available: {str(e)}")
            return

        obj_name = self.parent_window.mutation_structure_combo.currentText().strip()
        mutation_text = self.parent_window.mutation_list.text().strip()

        if not obj_name:
            show_message_box(self, "Warning", "请选择要分析的结构对象", "warning")
            return

        if not mutation_text:
            show_message_box(self, "Warning", "请输入突变列表（如 A50G,A51V）", "warning")
            return

        try:
            if obj_name not in cmd.get_names("objects"):
                show_message_box(self, "Warning", f"结构对象不存在：{obj_name}", "warning")
                return
        except Exception as e:
            show_message_box(self, "Error", f"无法获取对象列表：{str(e)}", "critical")
            return

        try:
            mutations = self._parse_mutation_list(mutation_text)
        except ValueError as e:
            show_message_box(self, "Warning", str(e), "warning")
            return

        self.log(f"🧬 Starting mutation analysis: {obj_name}")
        self.log(f"   Mutations: {', '.join([f'{c}{r}{a}' for c, r, a in mutations])}")

        results = analyze_mutation_effects(obj_name, mutations, partner_sel=None, output_csv=None, method="auto")
        if not results:
            show_message_box(self, "Mutation Analysis", "突变分析失败，请检查日志输出", "warning")
            self.log("❌ Mutation analysis failed")
            return

        ddg_result = results.get("ddg_result") if isinstance(results, dict) else None
        if isinstance(ddg_result, dict) and ddg_result.get("ddg") is not None:
            self.log(f"✅ ΔΔG: {ddg_result.get('ddg'):.3f} kcal/mol")
        else:
            self.log("⚠️ ΔΔG 结果不可用或未计算")

        show_message_box(self, "Mutation Analysis", "突变分析完成，详情请查看日志", "info")

        return None