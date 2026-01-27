
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
        """初始化UI - 现代卡片式布局"""
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
        self.parent_window.vina_selection.setMaximumWidth(80)
        get_center_btn = QPushButton("Get Center")
        get_center_btn.setMinimumHeight(28)
        get_center_btn.clicked.connect(self.get_center_from_selection)
        
        self.parent_window.vina_cx = QLineEdit()
        self.parent_window.vina_cx.setPlaceholderText("X")
        self.parent_window.vina_cx.setFixedWidth(60)
        self.parent_window.vina_cx.setMinimumHeight(28)
        self.parent_window.vina_cy = QLineEdit()
        self.parent_window.vina_cy.setPlaceholderText("Y")
        self.parent_window.vina_cy.setFixedWidth(60)
        self.parent_window.vina_cy.setMinimumHeight(28)
        self.parent_window.vina_cz = QLineEdit()
        self.parent_window.vina_cz.setPlaceholderText("Z")
        self.parent_window.vina_cz.setFixedWidth(60)
        self.parent_window.vina_cz.setMinimumHeight(28)
        
        self.parent_window.vina_sx = QLineEdit("20")
        self.parent_window.vina_sx.setFixedWidth(35)
        self.parent_window.vina_sx.setMinimumHeight(28)
        self.parent_window.vina_sy = QLineEdit("20")
        self.parent_window.vina_sy.setFixedWidth(35)
        self.parent_window.vina_sy.setMinimumHeight(28)
        self.parent_window.vina_sz = QLineEdit("20")
        self.parent_window.vina_sz.setFixedWidth(35)
        self.parent_window.vina_sz.setMinimumHeight(28)
        
        self.parent_window.vina_exhaustiveness = QSpinBox()
        self.parent_window.vina_exhaustiveness.setRange(1, 32)
        self.parent_window.vina_exhaustiveness.setValue(8)
        self.parent_window.vina_exhaustiveness.setFixedWidth(50)
        self.parent_window.vina_exhaustiveness.setMinimumHeight(28)
        
        self.parent_window.vina_num_modes = QSpinBox()
        self.parent_window.vina_num_modes.setRange(1, 20)
        self.parent_window.vina_num_modes.setValue(9)
        self.parent_window.vina_num_modes.setFixedWidth(50)
        self.parent_window.vina_num_modes.setMinimumHeight(28)
        
        row2.addWidget(self.parent_window.vina_selection)
        row2.addWidget(get_center_btn)
        row2.addWidget(QLabel("Center:"))
        row2.addWidget(self.parent_window.vina_cx)
        row2.addWidget(self.parent_window.vina_cy)
        row2.addWidget(self.parent_window.vina_cz)
        row2.addWidget(QLabel("Size:"))
        row2.addWidget(self.parent_window.vina_sx)
        row2.addWidget(self.parent_window.vina_sy)
        row2.addWidget(self.parent_window.vina_sz)
        row2.addWidget(QLabel("Exhaust:"))
        row2.addWidget(self.parent_window.vina_exhaustiveness)
        row2.addWidget(QLabel("Modes:"))
        row2.addWidget(self.parent_window.vina_num_modes)
        row2.addStretch()
        layout.addLayout(row2)
        
        # Row 3: Output directory
        row3 = QHBoxLayout()
        row3.setSpacing(6)
        
        self.parent_window.vina_output_dir = QLineEdit()
        self.parent_window.vina_output_dir.setPlaceholderText("Output directory (optional)")
        self.parent_window.vina_output_dir.setMinimumHeight(28)
        self.parent_window.vina_output_dir.setMinimumWidth(200)
        out_browse = QPushButton("Browse")
        out_browse.setMinimumHeight(28)
        out_browse.clicked.connect(lambda: self._browse_directory(self.parent_window.vina_output_dir, "Select Output Directory"))
        
        row3.addWidget(QLabel("Output:"))
        row3.addWidget(self.parent_window.vina_output_dir)
        row3.addWidget(out_browse)
        row3.addStretch()
        layout.addLayout(row3)
        
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

        # Row 1: Receptor + Ligand files
        row1 = QHBoxLayout()
        row1.setSpacing(6)

        self.parent_window.hdock_receptor = QLineEdit()
        self.parent_window.hdock_receptor.setPlaceholderText("Receptor PDB")
        self.parent_window.hdock_receptor.setMinimumHeight(28)
        self.parent_window.hdock_receptor.setMinimumWidth(150)
        rec_browse = QPushButton("Browse")
        rec_browse.setMinimumHeight(28)
        rec_browse.clicked.connect(self.browse_hdock_receptor)

        self.parent_window.hdock_ligand = QLineEdit()
        self.parent_window.hdock_ligand.setPlaceholderText("Ligand PDB")
        self.parent_window.hdock_ligand.setMinimumHeight(28)
        self.parent_window.hdock_ligand.setMinimumWidth(150)
        lig_browse = QPushButton("Browse")
        lig_browse.setMinimumHeight(28)
        lig_browse.clicked.connect(self.browse_hdock_ligand)

        row1.addWidget(QLabel("Receptor:"))
        row1.addWidget(self.parent_window.hdock_receptor)
        row1.addWidget(rec_browse)
        row1.addSpacing(10)
        row1.addWidget(QLabel("Ligand:"))
        row1.addWidget(self.parent_window.hdock_ligand)
        row1.addWidget(lig_browse)
        row1.addStretch()
        layout.addLayout(row1)

        # Row 2: Mode + site residues
        row2 = QHBoxLayout()
        row2.setSpacing(6)

        self.parent_window.haddock_mode = QComboBox()
        self.parent_window.haddock_mode.addItems([
            "Blind (Random AIR)", "Blind (Centroid)", "Blind (Surface)", "Pocket-constrained"
        ])
        self.parent_window.haddock_mode.setMinimumHeight(28)
        self.parent_window.haddock_mode.setMinimumWidth(130)

        self.parent_window.hdock_rsite = QLineEdit()
        self.parent_window.hdock_rsite.setPlaceholderText("Rec site: 195:A,203-206:A")
        self.parent_window.hdock_rsite.setMinimumHeight(28)
        self.parent_window.hdock_rsite.setMinimumWidth(140)
        self.parent_window.hdock_rsite.setEnabled(False)

        self.parent_window.hdock_lsite = QLineEdit()
        self.parent_window.hdock_lsite.setPlaceholderText("Lig site: 108:B,120-123:B")
        self.parent_window.hdock_lsite.setMinimumHeight(28)
        self.parent_window.hdock_lsite.setMinimumWidth(140)
        self.parent_window.hdock_lsite.setEnabled(False)

        row2.addWidget(QLabel("Mode:"))
        row2.addWidget(self.parent_window.haddock_mode)
        row2.addWidget(self.parent_window.hdock_rsite)
        row2.addWidget(self.parent_window.hdock_lsite)
        row2.addStretch()
        layout.addLayout(row2)

        # Row 3: Options + output
        row3 = QHBoxLayout()
        row3.setSpacing(6)

        self.parent_window.haddock_auto_passive = QCheckBox("Auto passive (6.5Å)")
        self.parent_window.haddock_auto_passive.setChecked(True)
        self.parent_window.haddock_auto_passive.setEnabled(False)

        self.parent_window.hdock_output = QLineEdit()
        self.parent_window.hdock_output.setPlaceholderText("Output dir (optional)")
        self.parent_window.hdock_output.setMinimumHeight(28)
        self.parent_window.hdock_output.setMinimumWidth(150)
        out_browse = QPushButton("Browse")
        out_browse.setMinimumHeight(28)
        out_browse.clicked.connect(self.browse_hdock_output)

        row3.addWidget(self.parent_window.haddock_auto_passive)
        row3.addWidget(QLabel("Output:"))
        row3.addWidget(self.parent_window.hdock_output)
        row3.addWidget(out_browse)
        row3.addStretch()
        layout.addLayout(row3)

        # Toggle site fields by mode
        def _toggle_mode(idx: int):
            pocket_mode = self.parent_window.haddock_mode.currentText().startswith("Pocket")
            self.parent_window.hdock_rsite.setEnabled(pocket_mode)
            self.parent_window.hdock_lsite.setEnabled(pocket_mode)
            self.parent_window.haddock_auto_passive.setEnabled(pocket_mode)
        self.parent_window.haddock_mode.currentIndexChanged.connect(_toggle_mode)

        return card

    # --- Vina Logic ---
    def browse_vina_ligand(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Ligand (Cancel for Folder)", "", 
            "Ligand Files (*.mol2 *.sdf *.pdbqt *.mol *.pdb);;All (*)")
        if fn:
            self.parent_window.vina_ligand.setText(fn)
        else:
            folder = QFileDialog.getExistingDirectory(self, "Select Ligand Folder (Batch)")
            if folder:
                self.parent_window.vina_ligand.setText(folder)

    def get_center_from_selection(self):
        try:
            from pymol import cmd
            selection = self.parent_window.vina_selection.text().strip()
            if not selection:
                show_message_box(self, "Warning", "Enter a PyMOL selection (e.g. resn LIG)", "warning")
                return
            try:
                com = cmd.centerofmass(selection)
                if com:
                    self.parent_window.vina_cx.setText(f"{com[0]:.2f}")
                    self.parent_window.vina_cy.setText(f"{com[1]:.2f}")
                    self.parent_window.vina_cz.setText(f"{com[2]:.2f}")
                    self.log(f"Center: ({com[0]:.2f}, {com[1]:.2f}, {com[2]:.2f})")
                else:
                    show_message_box(self, "Warning", f"Could not get center for '{selection}'", "warning")
            except Exception as e:
                show_message_box(self, "Warning", f"Error: {e}", "warning")
        except Exception as e:
            self.on_error(str(e))

    def run_vina_docking(self):
        try:
            from pymol import cmd

            receptor_obj = self.parent_window.vina_receptor_combo.currentText().strip()
            if not receptor_obj or receptor_obj == t("no_object"):
                show_message_box(self, "Warning", "Select a receptor object", "warning")
                return

            ligand_path = self.parent_window.vina_ligand.text().strip()
            if not ligand_path or not os.path.exists(ligand_path):
                show_message_box(self, "Warning", "Select a valid ligand file or folder", "warning")
                return

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
                show_message_box(self, "Warning", "Invalid box parameters. Use 'Get Center' first.", "warning")
                return

            exhaustiveness = self.parent_window.vina_exhaustiveness.value()
            num_modes = self.parent_window.vina_num_modes.value()
            output_dir = self.parent_window.vina_output_dir.text().strip() or None
            remove_selection = self.parent_window.vina_selection.text().strip() or None

            try:
                from ...vina_integration import batch_docking, get_ligand_files, check_vina_available
            except ImportError:
                from vina_integration import batch_docking, get_ligand_files, check_vina_available

            if not check_vina_available():
                show_message_box(self, "Error", "AutoDock Vina not found", "warning")
                return

            ligand_files = get_ligand_files(ligand_path)
            if not ligand_files:
                show_message_box(self, "Warning", "No valid ligand files found", "warning")
                return

            is_batch = len(ligand_files) > 1
            
            if is_batch:
                self.log(f"Batch docking: {len(ligand_files)} ligands")
                if remove_selection:
                    self.log(f"Removing co-crystal: {remove_selection}")
            else:
                self.log(f"Vina: {receptor_obj} + {os.path.basename(ligand_files[0])}")
            
            self.parent_window.vina_dock_btn.setEnabled(False)
            self.parent_window.vina_dock_btn.setText("Running...")
            self.parent_window.repaint()
            
            def progress_callback(current, total, name):
                self.log(f"  [{current}/{total}] {name}")
                self.parent_window.repaint()
            
            result = batch_docking(
                receptor_obj, ligand_path, box,
                output_dir=output_dir,
                exhaustiveness=exhaustiveness,
                num_modes=num_modes,
                remove_selection=remove_selection,
                progress_callback=progress_callback if is_batch else None
            )
            
            self.parent_window.vina_dock_btn.setEnabled(True)
            self.parent_window.vina_dock_btn.setText("Run Docking")
            
            if result.get('success'):
                if is_batch:
                    self.log(f"✅ Batch complete: {result['successful']}/{result['total']} successful")
                    self.log(f"   Results: {result['csv_path']}")
                    show_message_box(self, "Success",
                        f"Batch docking complete.\n"
                        f"Successful: {result['successful']}/{result['total']}\n"
                        f"Results saved to: {result['csv_path']}")
                else:
                    affinity = result['results'][0].get('affinity', 0) if result.get('results') else 0
                    self.log(f"✅ Done! Best: {affinity:.2f} kcal/mol")
                    show_message_box(self, "Success", f"Docking complete.\nBest: {affinity:.2f} kcal/mol")
            else:
                self.log(f"❌ Failed: {result.get('error')}")
                show_message_box(self, "Error", f"Docking failed: {result.get('error')}", "warning")
        except Exception as e:
            self.on_error(str(e))
            self.parent_window.vina_dock_btn.setEnabled(True)
            self.parent_window.vina_dock_btn.setText("Run Docking")

    def load_vina_result(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Result", "", "PDBQT (*.pdbqt);;All (*)")
        if fn:
            try:
                from pymol import cmd
                obj_name = os.path.splitext(os.path.basename(fn))[0]
                cmd.load(fn, obj_name)
                self.log(f"Loaded: {obj_name}")
                self.refresh_objects()
            except Exception as e:
                self.on_error(str(e))

    # --- HADDOCK3 Logic ---
    def browse_hdock_receptor(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Receptor PDB", "", "PDB (*.pdb)")
        if fn: self.parent_window.hdock_receptor.setText(fn)

    def browse_hdock_ligand(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Ligand PDB", "", "PDB (*.pdb)")
        if fn: self.parent_window.hdock_ligand.setText(fn)

    def browse_hdock_output(self):
        fn = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if fn: self.parent_window.hdock_output.setText(fn)

    def run_hdock(self):
        rec = self.parent_window.hdock_receptor.text().strip()
        lig = self.parent_window.hdock_ligand.text().strip()

        if not rec or not lig:
            show_message_box(self, "Warning", "Select both Receptor and Ligand PDB files", "warning")
            return
        if not os.path.exists(rec) or not os.path.exists(lig):
            show_message_box(self, "Warning", "Selected files do not exist", "warning")
            return

        self.log(f"HADDOCK3: {os.path.basename(rec)} + {os.path.basename(lig)}")

        try:
            from ...haddock3_integration import check_haddock3_available, Haddock3Runner
            info = check_haddock3_available()
            if not info.get('available'):
                show_message_box(self, "Error", "HADDOCK3 not found. Install: pip install -U haddock3", "warning")
                return

            runner = Haddock3Runner()
            self.parent_window.hdock_run_btn.setEnabled(False)
            self.parent_window.hdock_run_btn.setText("Running...")
            self.parent_window.repaint()

            mode_map = {
                "Blind (Random AIR)": "blind_ranair",
                "Blind (Centroid)": "blind_cm",
                "Blind (Surface)": "blind_surf",
                "Pocket-constrained": "air_from_residues",
            }
            mode = mode_map.get(self.parent_window.haddock_mode.currentText(), 'blind_ranair')

            result = runner.run_docking(
                rec, lig,
                output_dir=self.parent_window.hdock_output.text().strip() or None,
                rsite=self.parent_window.hdock_rsite.text().strip() or None,
                lsite=self.parent_window.hdock_lsite.text().strip() or None,
                mode=mode,
                expand_passive=self.parent_window.haddock_auto_passive.isChecked()
            )

            self.parent_window.hdock_run_btn.setEnabled(True)
            self.parent_window.hdock_run_btn.setText("Run HADDOCK3")

            if result.get('success'):
                self.log("✅ HADDOCK3 complete!")
                from pymol import cmd
                cmd.load(result['models_pdb'], "haddock3_models")
                show_message_box(self, "Success", "Docking complete. Models loaded.")
            else:
                self.log(f"❌ Failed: {result.get('error')}")
                show_message_box(self, "Error", f"HADDOCK3 failed: {result.get('error')}", "critical")

        except Exception as e:
            self.on_error(str(e))
            self.parent_window.hdock_run_btn.setEnabled(True)
            self.parent_window.hdock_run_btn.setText("Run HADDOCK3")

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
            if hasattr(cmd, 'perform_mutation'):
                cmd.perform_mutation(obj_name, mutations, method='pymol')
                self.log("Mutation performed")
            else:
                for chain, resi, resn in mutations:
                    sel = f"/{obj_name}//{chain}/{resi}"
                    cmd.wizard("mutagenesis")
                    cmd.get_wizard().do_select(sel)
                    cmd.get_wizard().set_mode(resn)
                    cmd.get_wizard().apply()
                    cmd.set_wizard()
                self.log("Mutation applied via PyMOL Wizard")
        except Exception as e:
            self.on_error(str(e))

    def run_minimize(self):
        obj = self.parent_window.mut_obj_combo.currentText()
        if not obj: return
        try:
            from pymol import cmd
            if hasattr(cmd, 'minimize_energy'):
                cmd.minimize_energy(obj)
            else:
                cmd.protect(f"not {obj}")
                cmd.sculpt_activate(obj)
                cmd.sculpt_iterate(obj, cycles=100)
                self.log("minimization done (sculpt)")
        except Exception as e:
            self.on_error(str(e))

    def run_mutation_analysis(self):
        """Run full mutation ΔΔG analysis using FoldX"""
        obj_name = self.parent_window.mut_obj_combo.currentText()
        if not obj_name or obj_name == t("no_object"):
            show_message_box(self, "Warning", "Please select a structure object", "warning")
            return

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
