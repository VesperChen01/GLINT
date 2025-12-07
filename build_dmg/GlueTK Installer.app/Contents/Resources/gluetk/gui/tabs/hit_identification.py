# -*- coding: utf-8 -*-
"""
Hit Identification Tab: Binding Site Detection, Vina Docking, HDOCK
"""
import os
from typing import Optional

try:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
        QGroupBox, QScrollArea, QFrame, QFileDialog, QMessageBox
    )
except ImportError:
    try:
        from PyQt6.QtCore import Qt
        from PyQt6.QtWidgets import (
            QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
            QGroupBox, QScrollArea, QFrame, QFileDialog, QMessageBox
        )
    except ImportError:
        raise RuntimeError("PyQt5 or PyQt6 must be installed.")

from ..utils import t
from .common import CommonTab

class HitIdentificationTab(CommonTab):
    def __init__(self, parent):
        super().__init__(parent)
        
        self._detected_pockets = []
        
        # Remove proxies
        for attr in ['browse_vina_ligand', 'run_vina_docking', 'load_vina_result',
                     'run_pocket_detection', 'run_pocket_visualization',
                     'browse_hdock_receptor', 'browse_hdock_ligand', 'run_hdock']:
            if attr in self.__dict__:
                del self.__dict__[attr]

        self.init_ui()
        
    def init_ui(self):
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        
        content_widget = QWidget()
        content_widget.setObjectName("scroll_content")
        self.parent_window._hit_scroll_content = content_widget
        bg_color = "#0d1117" if getattr(self.parent_window, "_dark_mode", True) else "#f8fafc"
        content_widget.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")
        
        layout = QVBoxLayout(content_widget)
        layout.setSpacing(14)
        layout.setContentsMargins(12, 12, 12, 12)
        
        # 1. Binding Site Detection
        layout.addWidget(self._create_pocket_detection_card())
        
        # 2. Vina Docking
        layout.addWidget(self._create_vina_docking_card())
        
        # 3. HDOCK (Protein-Protein Docking)
        layout.addWidget(self._create_hdock_card())
        
        layout.addStretch(1)
        scroll_area.setWidget(content_widget)
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(scroll_area)

    def _create_hdock_card(self) -> QWidget:
        """创建 HDOCK 蛋白-蛋白对接卡片"""
        card = QGroupBox("HDOCK - Protein-Protein Docking")
        layout = QVBoxLayout(card)
        layout.setSpacing(10)
        layout.setContentsMargins(16, 20, 16, 16)
        
        # Receptor (E3+Glue)
        rec_layout = QHBoxLayout()
        rec_layout.setSpacing(8)
        
        self.parent_window.hdock_receptor = QLineEdit()
        self.parent_window.hdock_receptor.setPlaceholderText("Receptor PDB (e.g., E3+Glue)")
        self.parent_window.hdock_receptor.setMinimumHeight(36)
        
        self.parent_window.hdock_rec_browse = QPushButton("Browse")
        self.parent_window.hdock_rec_browse.setObjectName("browse_btn")
        self.parent_window.hdock_rec_browse.setMinimumHeight(36)
        self.parent_window.hdock_rec_browse.clicked.connect(self.browse_hdock_receptor)
        
        rec_layout.addWidget(QLabel("Receptor:"))
        rec_layout.addWidget(self.parent_window.hdock_receptor, 1)
        rec_layout.addWidget(self.parent_window.hdock_rec_browse)
        layout.addLayout(rec_layout)
        
        # Ligand (POI)
        lig_layout = QHBoxLayout()
        lig_layout.setSpacing(8)
        
        self.parent_window.hdock_ligand = QLineEdit()
        self.parent_window.hdock_ligand.setPlaceholderText("Ligand PDB (e.g., POI/Substrate)")
        self.parent_window.hdock_ligand.setMinimumHeight(36)
        
        self.parent_window.hdock_lig_browse = QPushButton("Browse")
        self.parent_window.hdock_lig_browse.setObjectName("browse_btn")
        self.parent_window.hdock_lig_browse.setMinimumHeight(36)
        self.parent_window.hdock_lig_browse.clicked.connect(self.browse_hdock_ligand)
        
        lig_layout.addWidget(QLabel("Ligand:"))
        lig_layout.addWidget(self.parent_window.hdock_ligand, 1)
        lig_layout.addWidget(self.parent_window.hdock_lig_browse)
        layout.addLayout(lig_layout)
        
        # Run button
        self.parent_window.hdock_run_btn = QPushButton("Run HDOCK")
        self.parent_window.hdock_run_btn.setObjectName("primary_btn")
        self.parent_window.hdock_run_btn.setMinimumHeight(36)
        self.parent_window.hdock_run_btn.clicked.connect(self.run_hdock)
        layout.addWidget(self.parent_window.hdock_run_btn)
        
        return card

    # --- Pocket Detection Logic ---
    def run_pocket_detection(self):
        try:
            obj = self.parent_window.pocket_obj_combo.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select an object")
                return
            try:
                grid_spacing = float(self.parent_window.pocket_grid_spacing.text().strip() or "0.6")
                min_volume = float(self.parent_window.pocket_min_volume.text().strip() or "20")
            except ValueError:
                QMessageBox.warning(self, "Warning", "Invalid parameter values")
                return
            
            self.log(f"\n🔍 Detecting pockets in {obj}...")
            try: from ...pocket_detector import detect_pockets
            except ImportError: from pocket_detector import detect_pockets
            
            pockets = detect_pockets(obj, grid_spacing, min_volume)
            if not pockets:
                self.log("⚠️  No pockets detected")
                QMessageBox.information(self, "Result", "No pockets detected.")
                return
            
            self._detected_pockets = pockets
            self.log(f"\n✅ Detected {len(pockets)} pockets")
            self.run_pocket_visualization()
            
        except Exception as e:
            self.on_error(str(e))

    def run_pocket_visualization(self):
        try:
            if not self._detected_pockets:
                QMessageBox.warning(self, "Warning", "Please run pocket detection first")
                return
            
            color_by_map = {"Volume": "volume", "Druggability": "druggability", "Hydrophobicity": "hydrophobicity", "Depth": "depth"}
            color_by = color_by_map.get(self.parent_window.pocket_color_by.currentText(), "volume")
            
            try: from ...pocket_visualizer import visualize_pockets
            except ImportError: from pocket_visualizer import visualize_pockets
            
            visualize_pockets(self._detected_pockets, 'gluetk_pockets', color_by=color_by, show_spheres=True, sphere_radius=1.5)
            self.log("Pockets visualized")
        except Exception as e:
            self.on_error(str(e))

    # --- Vina Logic ---
    def browse_vina_ligand(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Ligand", "", "MOL2/SDF/PDBQT (*.mol2 *.sdf *.pdbqt);;All Files (*)")
        if fn:
            self.parent_window.vina_ligand.setText(fn)

    def run_vina_docking(self):
        try:
            from pymol import cmd
            
            if not hasattr(self.parent_window, 'pocket_obj_combo'):
                 QMessageBox.warning(self, "Warning", "Receptor selector not found (init Target Tab first)")
                 return

            receptor_obj = self.parent_window.pocket_obj_combo.currentText().strip()
            if not receptor_obj or receptor_obj == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select a receptor object from Pocket Detection")
                return
            
            try:
                names = cmd.get_names("objects")
            except:
                try:
                    names = cmd.get_object_list()
                except:
                    names = []
                
            if receptor_obj not in names:
                QMessageBox.warning(self, "Warning", f"Object '{receptor_obj}' not found in PyMOL")
                return
            
            ligand = self.parent_window.vina_ligand.text().strip()
            if not ligand or not os.path.exists(ligand):
                QMessageBox.warning(self, "Warning", "Please select a valid ligand file")
                return
            
            max_pockets = self.parent_window.vina_max_pockets.value()
            exhaustiveness = self.parent_window.vina_exhaustiveness.value()
            
            self.log(f"\nStarting Vina docking...")
            self.log(f"   Receptor: {receptor_obj}")
            self.log(f"   Ligand: {os.path.basename(ligand)}")
            
            try: from ...vina_integration import pocket_based_docking, manual_box_docking, check_vina_available
            except ImportError: from vina_integration import pocket_based_docking, manual_box_docking, check_vina_available
            
            if not check_vina_available():
                QMessageBox.warning(self, "Vina Not Found", "AutoDock Vina is not installed.")
                return
            
            if self.parent_window.vina_use_custom_box.isChecked():
                try:
                    box = {
                        'center_x': float(self.parent_window.vina_cx.text().strip()),
                        'center_y': float(self.parent_window.vina_cy.text().strip()),
                        'center_z': float(self.parent_window.vina_cz.text().strip()),
                        'size_x': float(self.parent_window.vina_sx.text().strip()),
                        'size_y': float(self.parent_window.vina_sy.text().strip()),
                        'size_z': float(self.parent_window.vina_sz.text().strip())
                    }
                except ValueError:
                    QMessageBox.warning(self, "Warning", "Invalid custom box parameters")
                    return
                
                self.log("   Using custom docking box")
                result = manual_box_docking(receptor_obj, ligand, box, exhaustiveness)
            else:
                self.log(f"   Auto-detecting max {max_pockets} pockets...")
                result = pocket_based_docking(receptor_obj, ligand, max_pockets, exhaustiveness)
            
            if result.get('success'):
                self.log("Docking complete!")
                if 'results' in result:
                    sorted_results = sorted([r for r in result['results'] if r.get('success')], key=lambda r: r.get('affinity', 0))
                    for i, r in enumerate(sorted_results[:3], 1):
                        self.log(f"   {i}. Pocket {r.get('pocket_id','Custom')}: {r.get('affinity'):.2f} kcal/mol")
                QMessageBox.information(self, "Success", f"Docking completed.\nSaved to: {result.get('output_dir')}")
            else:
                self.log(f"Docking failed: {result.get('error')}")
                QMessageBox.warning(self, "Error", f"Docking failed: {result.get('error')}")
                
        except Exception as e:
            self.on_error(str(e))

    def load_vina_result(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Docking Result", "", "PDBQT (*.pdbqt);;All Files (*)")
        if fn:
            try:
                from pymol import cmd
                obj_name = os.path.splitext(os.path.basename(fn))[0]
                cmd.load(fn, obj_name)
                self.log(f"Loaded docking result: {obj_name}")
                self.refresh_objects()
            except Exception as e:
                self.on_error(str(e))

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
            
            self.parent_window.hdock_run_btn.setEnabled(False)
            self.parent_window.hdock_run_btn.setText("Running HDOCK...")
            self.parent_window.repaint()
            
            result = runner.run_docking(rec, lig)
            
            self.parent_window.hdock_run_btn.setEnabled(True)
            self.parent_window.hdock_run_btn.setText("Run HDOCK")
            
            if result['success']:
                self.log("HDOCK complete!")
                self.log(f"   Models: {result['models_pdb']}")
                
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
            self.parent_window.hdock_run_btn.setText("Run HDOCK")
