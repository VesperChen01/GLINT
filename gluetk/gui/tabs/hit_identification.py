# -*- coding: utf-8 -*-
"""
Hit Identification Tab: Binding Site Detection, Vina Docking, HADDOCK3
"""
import os
from typing import Optional

from ..qt_adapter import (
    Qt, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QGroupBox, QScrollArea, QFrame, QFileDialog, QMessageBox,
    QComboBox, QCheckBox
)

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
        self.setObjectName("scroll_content")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.parent_window._hit_scroll_content = self

        bg_color = "#0d1117" if getattr(self.parent_window, "_dark_mode", True) else "#f8fafc"
        self.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")

        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(14)
        main_layout.setContentsMargins(12, 12, 12, 12)

        # 1. Binding Site Detection
        main_layout.addWidget(self._create_pocket_detection_card())

        # 2. Vina Docking
        main_layout.addWidget(self._create_vina_docking_card())

        # 3. HADDOCK3 (Protein-Protein Docking)
        main_layout.addWidget(self._create_hdock_card())

        main_layout.addStretch(1)

    def _create_hdock_card(self) -> QWidget:
        """Create HADDOCK3 protein-protein docking card"""
        card = QGroupBox("HADDOCK3 - Protein-Protein Docking")
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
        
        # Docking mode
        mode_layout = QHBoxLayout()
        mode_layout.setSpacing(8)
        self.parent_window.haddock_mode = QComboBox()
        self.parent_window.haddock_mode.addItems([
            "Blind docking (Random AIR)",
            "Blind docking (Centroid)",
            "Blind docking (Surface)",
            "Pocket-constrained (Residue AIR)"
        ])
        self.parent_window.haddock_mode.setMinimumHeight(32)
        mode_layout.addWidget(QLabel("Mode:"))
        mode_layout.addWidget(self.parent_window.haddock_mode, 1)
        layout.addLayout(mode_layout)
        
        # Auto passive expansion (only for pocket-constrained)
        self.parent_window.haddock_auto_passive = QCheckBox("Auto expand passive residues (6.5 Å)")
        self.parent_window.haddock_auto_passive.setChecked(True)
        layout.addWidget(self.parent_window.haddock_auto_passive)
        
        # Receptor/Ligand site residues
        rsite_layout = QHBoxLayout()
        rsite_layout.setSpacing(8)
        self.parent_window.hdock_rsite = QLineEdit()
        self.parent_window.hdock_rsite.setPlaceholderText("Rec site residues, e.g., 195:A,203-206:A")
        self.parent_window.hdock_rsite.setMinimumHeight(32)
        rsite_layout.addWidget(QLabel("Rec site:"))
        rsite_layout.addWidget(self.parent_window.hdock_rsite, 1)
        layout.addLayout(rsite_layout)

        lsite_layout = QHBoxLayout()
        lsite_layout.setSpacing(8)
        self.parent_window.hdock_lsite = QLineEdit()
        self.parent_window.hdock_lsite.setPlaceholderText("Lig site residues, e.g., 108:B,120-123:B")
        self.parent_window.hdock_lsite.setMinimumHeight(32)
        lsite_layout.addWidget(QLabel("Lig site:"))
        lsite_layout.addWidget(self.parent_window.hdock_lsite, 1)
        layout.addLayout(lsite_layout)
        
        # Output directory (optional)
        out_layout = QHBoxLayout()
        out_layout.setSpacing(8)
        self.parent_window.hdock_output = QLineEdit()
        self.parent_window.hdock_output.setPlaceholderText("Output dir (default: same as receptor)")
        self.parent_window.hdock_output.setMinimumHeight(32)
        self.parent_window.hdock_out_browse = QPushButton("Browse")
        self.parent_window.hdock_out_browse.setObjectName("browse_btn")
        self.parent_window.hdock_out_browse.setMinimumHeight(32)
        self.parent_window.hdock_out_browse.clicked.connect(self.browse_hdock_output)
        out_layout.addWidget(QLabel("Output:"))
        out_layout.addWidget(self.parent_window.hdock_output, 1)
        out_layout.addWidget(self.parent_window.hdock_out_browse)
        layout.addLayout(out_layout)
        
        # Toggle visibility by mode
        def _toggle_mode(idx: int):
            text = self.parent_window.haddock_mode.currentText()
            pocket_mode = text.startswith("Pocket-constrained")
            for w in (self.parent_window.hdock_rsite, self.parent_window.hdock_lsite, self.parent_window.haddock_auto_passive):
                w.setEnabled(pocket_mode)
        self.parent_window.haddock_mode.currentIndexChanged.connect(_toggle_mode)
        _toggle_mode(self.parent_window.haddock_mode.currentIndex())
        
        # Run button
        self.parent_window.hdock_run_btn = QPushButton("Run HADDOCK3")
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

    def browse_hdock_output(self):
        fn = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if fn: self.parent_window.hdock_output.setText(fn)

    def run_hdock(self):
        rec = self.parent_window.hdock_receptor.text().strip()
        lig = self.parent_window.hdock_ligand.text().strip()
        
        if not rec or not lig:
            QMessageBox.warning(self, "Warning", "Please select both Receptor and Ligand PDB files.")
            return
            
        if not os.path.exists(rec) or not os.path.exists(lig):
            QMessageBox.warning(self, "Warning", "Selected files do not exist.")
            return
            
        self.log(f"Starting HADDOCK3...")
        self.log(f"   Receptor: {os.path.basename(rec)}")
        self.log(f"   Ligand: {os.path.basename(lig)}")
        
        try:
            from ...haddock3_integration import check_haddock3_available, Haddock3Runner
            info = check_haddock3_available()
            if not info.get('available'):
                QMessageBox.warning(
                    self,
                    "Error",
                    "HADDOCK3 not detected. Please install in the same Python/environment as PyMOL and restart:\n\n"
                    "pip install -U haddock3"
                )
                return
                
            runner = Haddock3Runner()
            
            self.parent_window.hdock_run_btn.setEnabled(False)
            self.parent_window.hdock_run_btn.setText("Running HADDOCK3...")
            self.parent_window.repaint()
            
            rsite = getattr(self.parent_window, 'hdock_rsite', None)
            rsite_txt = rsite.text().strip() if rsite else ""
            lsite = getattr(self.parent_window, 'hdock_lsite', None)
            lsite_txt = lsite.text().strip() if lsite else ""

            # Map UI mode to runner mode
            mode_map = {
                "Blind docking (Random AIR)": "blind_ranair",
                "Blind docking (Centroid)": "blind_cm",
                "Blind docking (Surface)": "blind_surf",
                "Pocket-constrained (Residue AIR)": "air_from_residues",
            }
            mode_txt = self.parent_window.haddock_mode.currentText()
            mode = mode_map.get(mode_txt, 'blind_ranair')
            expand_passive = self.parent_window.haddock_auto_passive.isChecked()
            
            # Get output directory (optional)
            output_dir = self.parent_window.hdock_output.text().strip() or None
            
            result = runner.run_docking(
                rec, lig,
                output_dir=output_dir,
                rsite=(rsite_txt or None), lsite=(lsite_txt or None),
                mode=mode, expand_passive=expand_passive
            )
            
            self.parent_window.hdock_run_btn.setEnabled(True)
            self.parent_window.hdock_run_btn.setText("Run HADDOCK3")
            
            if result.get('success'):
                self.log("HADDOCK3 complete!")
                # Show path fallback warning (if any)
                if result.get('warning'):
                    self.log(f"   ⚠️  {result['warning']}")
                self.log(f"   Models: {result['models_pdb']}")
                
                from pymol import cmd
                cmd.load(result['models_pdb'], "haddock3_models")
                self.log("   Loaded 'haddock3_models' into PyMOL")
                
                msg = "Docking complete. Top models loaded."
                if result.get('warning'):
                    msg += f"\n\n⚠️  {result['warning']}"
                QMessageBox.information(self, "Success", msg)
            else:
                self.log(f"HADDOCK3 failed: {result.get('error')}")
                QMessageBox.critical(self, "Error", f"HADDOCK3 failed: {result.get('error')}")
                
        except Exception as e:
            self.log(f"Error: {e}")
            self.on_error(str(e))
            self.parent_window.hdock_run_btn.setEnabled(True)
            self.parent_window.hdock_run_btn.setText("Run HADDOCK3")
