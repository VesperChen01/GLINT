# -*- coding: utf-8 -*-
"""
Hit Identification Tab: Vina Docking
"""
import os
from typing import Optional

try:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox, QScrollArea, QFrame,
        QFileDialog, QMessageBox
    )
except ImportError:
    try:
        from PyQt6.QtCore import Qt
        from PyQt6.QtWidgets import (
            QWidget, QVBoxLayout, QHBoxLayout, QLabel, QGroupBox, QScrollArea, QFrame,
            QFileDialog, QMessageBox
        )
    except ImportError:
        raise RuntimeError("PyQt5 or PyQt6 must be installed.")

from ..utils import t
from .common import CommonTab

class HitIdentificationTab(CommonTab):
    def __init__(self, parent):
        super().__init__(parent)
        
        # Remove proxies
        for attr in ['browse_vina_ligand', 'run_vina_docking', 'load_vina_result']:
            if attr in self.__dict__:
                del self.__dict__[attr]

        self.init_ui()
        
    def init_ui(self):
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        
        content_widget = QWidget()
        content_widget.setObjectName("scroll_content")
        bg_color = "#0d1117" if getattr(self.parent_window, "_dark_mode", True) else "#f8fafc"
        content_widget.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")
        
        layout = QVBoxLayout(content_widget)
        layout.setSpacing(14)
        layout.setContentsMargins(12, 12, 12, 12)
        
        # Vina Docking Card
        layout.addWidget(self._create_vina_docking_card())
        
        # Placeholder for future Virtual Screening
        grp_vs = QGroupBox("Virtual Screening (Coming Soon)")
        vs_layout = QVBoxLayout(grp_vs)
        vs_layout.addWidget(QLabel("Batch docking and scoring functionality will be available in future updates."))
        layout.addWidget(grp_vs)
        
        layout.addStretch(1)
        scroll_area.setWidget(content_widget)
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(scroll_area)

    # --- Vina Logic ---
    def browse_vina_ligand(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Ligand", "", "MOL2/SDF/PDBQT (*.mol2 *.sdf *.pdbqt);;All Files (*)")
        if fn:
            self.parent_window.vina_ligand.setText(fn)

    def run_vina_docking(self):
        try:
            from pymol import cmd
            
            # Try to get receptor object from parent's pocket combo
            # This assumes TargetDiscoveryTab has initialized this widget on parent
            if not hasattr(self.parent_window, 'pocket_obj_combo'):
                 QMessageBox.warning(self, "Warning", "Receptor selector not found (init Target Tab first)")
                 return

            receptor_obj = self.parent_window.pocket_obj_combo.currentText().strip()
            if not receptor_obj or receptor_obj == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select a receptor object from Pocket tab (Target Discovery)")
                return
            
            # Check if object exists in PyMOL
            # We can use get_names if available or fallback
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
            
            self.log(f"\\nStarting Vina docking...")
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
                QMessageBox.information(self, "Success", f"Docking completed.\\nSaved to: {result.get('output_dir')}")
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
