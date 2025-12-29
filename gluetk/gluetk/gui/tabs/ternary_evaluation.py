# -*- coding: utf-8 -*-
"""
Ternary Complex Evaluation Tab
三元复合物评估标签页 - 包含界面模块、配体模块、三元几何模块
"""
import os
from typing import Optional, Dict, Any

from ..qt_adapter import (
    Qt, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QComboBox, QGroupBox, QGridLayout, QScrollArea, QFrame, QTabWidget,
    QFileDialog, QMessageBox, QTextEdit, QDoubleSpinBox, QCheckBox
)

from ..utils import t
from .common import CommonTab


class TernaryEvaluationTab(CommonTab):
    """三元复合物评估标签页 - 与Target Discovery并列"""
    
    def __init__(self, parent):
        super().__init__(parent)
        self._last_result = None
        self.init_ui()
        
    def init_ui(self):
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        
        content_widget = QWidget()
        content_widget.setObjectName("scroll_content")
        bg_color = "#161b22" if getattr(self.parent_window, "_dark_mode", False) else "#ffffff"
        content_widget.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")
        
        layout = QVBoxLayout(content_widget)
        layout.setSpacing(16)
        layout.setContentsMargins(16, 16, 16, 16)
        
        # 标题
        title = QLabel("Ternary Complex Evaluation")
        title.setStyleSheet("font-size: 20px; font-weight: bold; color: #3b82f6;")
        layout.addWidget(title)
        
        desc = QLabel("Evaluate molecular glue ternary complexes: Interface, Ligand, and Geometry features")
        desc.setStyleSheet("color: #64748b; margin-bottom: 10px;")
        layout.addWidget(desc)
        
        # === 输入区域 ===
        grp_input = QGroupBox("Structure Input")
        input_layout = QGridLayout(grp_input)
        input_layout.setColumnStretch(1, 1)
        input_layout.setColumnStretch(3, 1)
        
        # Row 0: PyMOL Object / PDB File
        input_layout.addWidget(QLabel("PyMOL Object:"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.ternary_obj_combo = QComboBox()
        self.ternary_obj_combo.setMinimumHeight(28)
        # 同时保存到parent_window以便refresh_objects能找到
        self.parent_window.ternary_obj_combo = self.ternary_obj_combo
        
        refresh_btn = QPushButton(t("refresh"))
        refresh_btn.clicked.connect(self._do_refresh)
        r0 = QHBoxLayout()
        r0.addWidget(self.ternary_obj_combo, 1)
        r0.addWidget(refresh_btn)
        input_layout.addLayout(r0, 0, 1)
        
        input_layout.addWidget(QLabel("PDB File:"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ternary_pdb_path = QLineEdit()
        browse_btn = QPushButton(t("browse"))
        browse_btn.clicked.connect(self.browse_pdb)
        r0b = QHBoxLayout()
        r0b.addWidget(self.parent_window.ternary_pdb_path, 1)
        r0b.addWidget(browse_btn)
        input_layout.addLayout(r0b, 0, 3)
        
        # Row 1: Chain IDs
        input_layout.addWidget(QLabel("E3 Chain:"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ternary_e3_chain = QLineEdit("A")
        self.parent_window.ternary_e3_chain.setToolTip("E3 ligase chain ID (e.g., CRBN)")
        input_layout.addWidget(self.parent_window.ternary_e3_chain, 1, 1)
        
        input_layout.addWidget(QLabel("POI Chain:"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ternary_poi_chain = QLineEdit("B")
        self.parent_window.ternary_poi_chain.setToolTip("Protein of Interest chain ID")
        input_layout.addWidget(self.parent_window.ternary_poi_chain, 1, 3)
        
        # Row 2: Ligand
        input_layout.addWidget(QLabel("Ligand Chain:"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ternary_lig_chain = QLineEdit("C")
        self.parent_window.ternary_lig_chain.setToolTip("Molecular glue chain/residue ID")
        input_layout.addWidget(self.parent_window.ternary_lig_chain, 2, 1)
        
        input_layout.addWidget(QLabel("SMILES (optional):"), 2, 2, Qt.AlignmentFlag.AlignRight)
        self.parent_window.ternary_smiles = QLineEdit()
        self.parent_window.ternary_smiles.setToolTip("SMILES for ligand property calculation")
        input_layout.addWidget(self.parent_window.ternary_smiles, 2, 3)
        
        layout.addWidget(grp_input)
        
        # === 三个子模块标签页 ===
        sub_tabs = QTabWidget()
        sub_tabs.setStyleSheet("""
            QTabWidget::pane { border: 1px solid #e5e7eb; background: white; }
            QTabBar::tab { padding: 8px 16px; }
            QTabBar::tab:selected { background: #3b82f6; color: white; }
        """)
        
        # 1. 界面模块
        interface_tab = self._create_interface_tab()
        sub_tabs.addTab(interface_tab, "Interface Module")
        
        # 2. 配体模块
        ligand_tab = self._create_ligand_tab()
        sub_tabs.addTab(ligand_tab, "Ligand Module")
        
        # 3. 三元几何模块
        geometry_tab = self._create_geometry_tab()
        sub_tabs.addTab(geometry_tab, "Geometry Module")
        
        layout.addWidget(sub_tabs)
        
        # === 运行按钮 ===
        btn_row = QHBoxLayout()
        self.parent_window.ternary_run_all_btn = QPushButton("Run Full Evaluation")
        self.parent_window.ternary_run_all_btn.setObjectName("primary_btn")
        self.parent_window.ternary_run_all_btn.clicked.connect(self.run_full_evaluation)
        btn_row.addWidget(self.parent_window.ternary_run_all_btn)
        
        self.parent_window.ternary_visualize_btn = QPushButton("Visualize")
        self.parent_window.ternary_visualize_btn.setObjectName("highlight_btn")
        self.parent_window.ternary_visualize_btn.clicked.connect(self.visualize_geometry)
        btn_row.addWidget(self.parent_window.ternary_visualize_btn)
        
        export_btn = QPushButton("Export CSV")
        export_btn.clicked.connect(self.export_results)
        btn_row.addWidget(export_btn)
        
        btn_row.addStretch(1)
        layout.addLayout(btn_row)
        
        # === 结果显示 ===
        grp_result = QGroupBox("Evaluation Results")
        result_layout = QVBoxLayout(grp_result)
        self.parent_window.ternary_result_text = QTextEdit()
        self.parent_window.ternary_result_text.setReadOnly(True)
        self.parent_window.ternary_result_text.setMinimumHeight(250)
        self.parent_window.ternary_result_text.setStyleSheet(
            "font-family: 'Consolas', 'Monaco', monospace; font-size: 12px;"
        )
        result_layout.addWidget(self.parent_window.ternary_result_text)
        layout.addWidget(grp_result)
        
        layout.addStretch(1)
        scroll_area.setWidget(content_widget)
        
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(scroll_area)

    def _do_refresh(self):
        """刷新PyMOL对象列表"""
        # 直接获取PyMOL对象，不依赖parent_window.refresh_objects
        names = []
        try:
            from pymol import cmd
            names = cmd.get_names("objects") if hasattr(cmd, "get_names") else cmd.get_object_list()
            self.log(f"PyMOL objects found: {names}")
        except Exception as e:
            self.log(f"PyMOL connection error: {e}")
        
        if not names:
            names = [t("no_object")]
        
        # 更新本地combo
        self.ternary_obj_combo.blockSignals(True)
        self.ternary_obj_combo.clear()
        self.ternary_obj_combo.addItems(names)
        self.ternary_obj_combo.blockSignals(False)
        
        self.log(f"Refreshed: {len(names)} objects")

    def _create_interface_tab(self) -> QWidget:
        """创建界面模块子标签页"""
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(12, 12, 12, 12)
        
        desc = QLabel("Calculate Buried Surface Area (BSA) and contact counts between components")
        desc.setWordWrap(True)
        desc.setStyleSheet("color: #64748b;")
        layout.addWidget(desc)
        
        grid = QGridLayout()
        grid.addWidget(QLabel("Probe Radius (Å):"), 0, 0)
        self.parent_window.ternary_probe_radius = QDoubleSpinBox()
        self.parent_window.ternary_probe_radius.setRange(0.5, 3.0)
        self.parent_window.ternary_probe_radius.setValue(1.4)
        self.parent_window.ternary_probe_radius.setSingleStep(0.1)
        grid.addWidget(self.parent_window.ternary_probe_radius, 0, 1)
        
        grid.addWidget(QLabel("Contact Distance (Å):"), 1, 0)
        self.parent_window.ternary_contact_dist = QDoubleSpinBox()
        self.parent_window.ternary_contact_dist.setRange(3.0, 8.0)
        self.parent_window.ternary_contact_dist.setValue(4.5)
        self.parent_window.ternary_contact_dist.setSingleStep(0.5)
        grid.addWidget(self.parent_window.ternary_contact_dist, 1, 1)
        
        layout.addLayout(grid)
        
        run_btn = QPushButton("Calculate Interface")
        run_btn.clicked.connect(self.run_interface_only)
        layout.addWidget(run_btn)
        
        layout.addStretch(1)
        return w

    def _create_ligand_tab(self) -> QWidget:
        """创建配体模块子标签页"""
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(12, 12, 12, 12)
        
        desc = QLabel("Calculate molecular glue properties: MW, LogP, TPSA, HBD/HBA, Fsp3, etc.")
        desc.setWordWrap(True)
        desc.setStyleSheet("color: #64748b;")
        layout.addWidget(desc)
        
        grid = QGridLayout()
        grid.addWidget(QLabel("SMILES:"), 0, 0)
        self.parent_window.ternary_lig_smiles = QLineEdit()
        self.parent_window.ternary_lig_smiles.setPlaceholderText("Enter SMILES or extract from structure")
        grid.addWidget(self.parent_window.ternary_lig_smiles, 0, 1)
        
        extract_btn = QPushButton("Extract from PDB")
        extract_btn.clicked.connect(self.extract_smiles)
        grid.addWidget(extract_btn, 0, 2)
        
        layout.addLayout(grid)
        
        # 属性显示
        props_grp = QGroupBox("Calculated Properties")
        props_layout = QGridLayout(props_grp)
        self.lig_prop_labels = {}
        props = ["MW", "LogP", "TPSA", "HBD", "HBA", "RotBonds", "Fsp3", "Rings"]
        for i, p in enumerate(props):
            props_layout.addWidget(QLabel(f"{p}:"), i // 4, (i % 4) * 2)
            lbl = QLabel("-")
            lbl.setStyleSheet("font-weight: bold;")
            self.lig_prop_labels[p] = lbl
            props_layout.addWidget(lbl, i // 4, (i % 4) * 2 + 1)
        layout.addWidget(props_grp)
        
        run_btn = QPushButton("Calculate Ligand Properties")
        run_btn.clicked.connect(self.run_ligand_only)
        layout.addWidget(run_btn)
        
        layout.addStretch(1)
        return w

    def _create_geometry_tab(self) -> QWidget:
        """创建三元几何模块子标签页"""
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(12, 12, 12, 12)
        
        desc = QLabel("Calculate ternary complex geometry: center of gravity, angles, distances, cooperativity")
        desc.setWordWrap(True)
        desc.setStyleSheet("color: #64748b;")
        layout.addWidget(desc)
        
        # 几何参数显示
        geom_grp = QGroupBox("Geometry Features")
        geom_layout = QGridLayout(geom_grp)
        self.geom_labels = {}
        features = [
            ("COG Shift", "Å"), ("Angle", "°"), 
            ("E3-POI Dist", "Å"), ("E3-MG Dist", "Å"), ("POI-MG Dist", "Å"),
            ("Cooperativity", "kcal/mol"), ("Hook Risk", ""), ("Duality", "")
        ]
        for i, (name, unit) in enumerate(features):
            geom_layout.addWidget(QLabel(f"{name}:"), i // 2, (i % 2) * 3)
            lbl = QLabel("-")
            lbl.setStyleSheet("font-weight: bold;")
            self.geom_labels[name] = lbl
            geom_layout.addWidget(lbl, i // 2, (i % 2) * 3 + 1)
            geom_layout.addWidget(QLabel(unit), i // 2, (i % 2) * 3 + 2)
        layout.addWidget(geom_grp)
        
        run_btn = QPushButton("Calculate Geometry")
        run_btn.clicked.connect(self.run_geometry_only)
        layout.addWidget(run_btn)
        
        layout.addStretch(1)
        return w

    def browse_pdb(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select PDB File", "", "PDB (*.pdb *.cif);;All (*)")
        if fn:
            self.parent_window.ternary_pdb_path.setText(fn)

    def _get_pdb_path(self) -> Optional[str]:
        """获取PDB路径，如果没有则从PyMOL导出"""
        pdb_path = self.parent_window.ternary_pdb_path.text().strip()
        if pdb_path and os.path.exists(pdb_path):
            return pdb_path
        
        obj = self.ternary_obj_combo.currentText().strip()
        if obj and obj != t("no_object"):
            try:
                from pymol import cmd
                import tempfile
                fd, pdb_path = tempfile.mkstemp(suffix=".pdb")
                os.close(fd)
                cmd.save(pdb_path, obj)
                return pdb_path
            except Exception as e:
                self.on_error(f"Failed to export PyMOL object: {e}")
        return None

    def run_full_evaluation(self):
        """运行完整评估"""
        pdb_path = self._get_pdb_path()
        if not pdb_path:
            QMessageBox.warning(self, "Warning", "Please provide a PDB file or select a PyMOL object.")
            return
        
        e3 = self.parent_window.ternary_e3_chain.text().strip()
        poi = self.parent_window.ternary_poi_chain.text().strip()
        lig = self.parent_window.ternary_lig_chain.text().strip()
        smiles = self.parent_window.ternary_smiles.text().strip() or self.parent_window.ternary_lig_smiles.text().strip()
        
        if not all([e3, poi, lig]):
            QMessageBox.warning(self, "Warning", "Please specify all chain IDs.")
            return
        
        self.parent_window.ternary_run_all_btn.setEnabled(False)
        self.parent_window.progress_bar.setVisible(True)
        self.parent_window.progress_bar.setRange(0, 0)
        
        from ..workers_ternary import TernaryEvaluationWorker
        self.parent_window.ternary_thread = TernaryEvaluationWorker(
            pdb_path, e3, poi, lig, smiles or None, None
        )
        self.parent_window.ternary_thread.progress.connect(self.log)
        self.parent_window.ternary_thread.error.connect(self.on_error)
        self.parent_window.ternary_thread.finished.connect(self._on_finished)
        self.parent_window.ternary_thread.start()

    def _on_finished(self, result: Dict[str, Any], out_csv: str):
        """评估完成"""
        self._last_result = result
        self.parent_window.progress_bar.setVisible(False)
        self.parent_window.progress_bar.setRange(0, 1)
        self.parent_window.ternary_run_all_btn.setEnabled(True)
        
        # 更新结果显示
        self.parent_window.ternary_result_text.setText(self._format_results(result))
        
        # 更新配体属性标签
        for key, lbl in self.lig_prop_labels.items():
            map_key = f"ligand_{key.lower()}"
            if key == "RotBonds":
                map_key = "ligand_rotatable_bonds"
            val = result.get(map_key, 0)
            lbl.setText(f"{val:.2f}" if isinstance(val, float) else str(val))
        
        # 更新几何标签
        geom_map = {
            "COG Shift": "geom_cog_shift",
            "Angle": "geom_angle_deg",
            "E3-POI Dist": "dist_e3_poi",
            "E3-MG Dist": "dist_e3_mg",
            "POI-MG Dist": "dist_poi_mg",
            "Cooperativity": "cooperativity_energy",
            "Hook Risk": "hook_risk_score",
            "Duality": "duality_index"
        }
        for name, key in geom_map.items():
            val = result.get(key, 0)
            self.geom_labels[name].setText(f"{val:.2f}" if isinstance(val, float) else str(val))
        
        self.log("✅ Full evaluation complete")

    def _format_results(self, result: Dict[str, Any]) -> str:
        """格式化结果"""
        lines = [
            "=" * 50,
            "  TERNARY COMPLEX EVALUATION",
            "=" * 50,
            "",
            "【Interface】",
            f"  Total BSA:     {result.get('bsa_total', 0):.1f} Å²",
            f"  MG-E3 BSA:     {result.get('bsa_mg_e3', 0):.1f} Å²",
            f"  MG-POI BSA:    {result.get('bsa_mg_poi', 0):.1f} Å²",
            f"  E3-POI BSA:    {result.get('bsa_e3_poi', 0):.1f} Å²",
            f"  Contacts:      {result.get('contact_count_45', 0)}",
            "",
            "【Geometry】",
            f"  COG Shift:     {result.get('geom_cog_shift', 0):.2f} Å",
            f"  Angle:         {result.get('geom_angle_deg', 0):.1f}°",
            f"  Cooperativity: {result.get('cooperativity_energy', 0):.2f} kcal/mol",
            f"  Hook Risk:     {result.get('hook_risk_score', 0):.2f}",
        ]
        
        if result.get('ligand_mw', 0) > 0:
            lines.extend([
                "",
                "【Ligand】",
                f"  MW:            {result.get('ligand_mw', 0):.1f} Da",
                f"  LogP:          {result.get('ligand_logp', 0):.2f}",
                f"  TPSA:          {result.get('ligand_tpsa', 0):.1f} Å²",
                f"  Fsp3:          {result.get('ligand_fsp3', 0):.2f}",
            ])
        
        return "\n".join(lines)

    def run_interface_only(self):
        """仅运行界面计算"""
        self.log("Running interface calculation...")
        # 简化版：直接调用完整评估
        self.run_full_evaluation()

    def run_ligand_only(self):
        """仅运行配体计算"""
        smiles = self.parent_window.ternary_lig_smiles.text().strip()
        if not smiles:
            QMessageBox.warning(self, "Warning", "Please enter SMILES.")
            return
        
        try:
            from ..ternary_complex_evaluator import LigandCalculator
            calc = LigandCalculator()
            props = calc.calculate(smiles)
            
            for key, lbl in self.lig_prop_labels.items():
                map_key = key.lower()
                if key == "RotBonds":
                    map_key = "rotatable_bonds"
                val = props.get(map_key, 0)
                lbl.setText(f"{val:.2f}" if isinstance(val, float) else str(val))
            
            self.log("✅ Ligand properties calculated")
        except Exception as e:
            self.on_error(str(e))

    def run_geometry_only(self):
        """仅运行几何计算"""
        self.log("Running geometry calculation...")
        self.run_full_evaluation()

    def extract_smiles(self):
        """从PDB提取SMILES"""
        QMessageBox.information(self, "Info", "SMILES extraction from PDB requires OpenBabel. Please enter SMILES manually.")

    def visualize_geometry(self):
        """可视化几何"""
        if not self._last_result:
            QMessageBox.warning(self, "Warning", "Please run evaluation first.")
            return
        
        try:
            from pymol import cmd, cgo
            
            obj = self.ternary_obj_combo.currentText().strip()
            e3 = self.parent_window.ternary_e3_chain.text().strip()
            poi = self.parent_window.ternary_poi_chain.text().strip()
            lig = self.parent_window.ternary_lig_chain.text().strip()
            
            cmd.show("cartoon", obj)
            cmd.color("lightblue", f"{obj} and chain {e3}")
            cmd.color("lightgreen", f"{obj} and chain {poi}")
            cmd.show("sticks", f"{obj} and chain {lig}")
            cmd.color("yellow", f"{obj} and chain {lig}")
            
            # 质心球
            cog_e3 = (self._last_result.get('cog_e3_x', 0),
                      self._last_result.get('cog_e3_y', 0),
                      self._last_result.get('cog_e3_z', 0))
            cog_poi = (self._last_result.get('cog_poi_x', 0),
                       self._last_result.get('cog_poi_y', 0),
                       self._last_result.get('cog_poi_z', 0))
            cog_mg = (self._last_result.get('cog_mg_x', 0),
                      self._last_result.get('cog_mg_y', 0),
                      self._last_result.get('cog_mg_z', 0))
            
            cgo_obj = [
                cgo.COLOR, 0.3, 0.5, 1.0,
                cgo.SPHERE, *cog_e3, 2.0,
                cgo.COLOR, 0.3, 0.8, 0.3,
                cgo.SPHERE, *cog_poi, 2.0,
                cgo.COLOR, 1.0, 0.8, 0.0,
                cgo.SPHERE, *cog_mg, 1.5,
            ]
            
            cmd.delete("ternary_cog")
            cmd.load_cgo(cgo_obj, "ternary_cog")
            cmd.zoom(obj, buffer=5.0)
            
            self.log("✅ Geometry visualized")
        except Exception as e:
            self.on_error(str(e))

    def export_results(self):
        """导出结果到CSV"""
        if not self._last_result:
            QMessageBox.warning(self, "Warning", "No results to export.")
            return
        
        fn, _ = QFileDialog.getSaveFileName(self, "Save Results", "", "CSV (*.csv)")
        if fn:
            import csv
            with open(fn, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(['Feature', 'Value'])
                for k, v in self._last_result.items():
                    writer.writerow([k, v])
            self.log(f"✅ Results exported to {fn}")