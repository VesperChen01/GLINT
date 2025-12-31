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
        self.setObjectName("scroll_content")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)

        bg_color = "#161b22" if getattr(self.parent_window, "_dark_mode", False) else "#ffffff"
        self.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")

        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        layout.setContentsMargins(12, 12, 12, 12)
        
        # 标题行
        title_row = QHBoxLayout()
        title = QLabel("Ternary Complex Evaluation")
        title.setStyleSheet("font-size: 16px; font-weight: bold; color: #3b82f6;")
        title_row.addWidget(title)
        title_row.addStretch(1)
        layout.addLayout(title_row)
        
        # === 输入区域 - 紧凑的单行布局 ===
        grp_input = QGroupBox("Structure Input")
        input_layout = QHBoxLayout(grp_input)
        input_layout.setContentsMargins(8, 10, 8, 8)
        input_layout.setSpacing(6)
        
        # PyMOL Object
        input_layout.addWidget(QLabel("Object:"))
        self.ternary_obj_combo = QComboBox()
        self.ternary_obj_combo.setMinimumHeight(26)
        self.ternary_obj_combo.setMinimumWidth(120)
        self.ternary_obj_combo.setMaximumWidth(180)
        self.ternary_obj_combo.setToolTip("Select PyMOL object")
        self.parent_window.ternary_obj_combo = self.ternary_obj_combo
        input_layout.addWidget(self.ternary_obj_combo)
        
        refresh_btn = QPushButton(t("refresh"))
        refresh_btn.setFixedWidth(70)
        refresh_btn.setMinimumHeight(26)
        refresh_btn.clicked.connect(self._do_refresh)
        input_layout.addWidget(refresh_btn)
        
        input_layout.addSpacing(10)
        
        # E3 Chain
        input_layout.addWidget(QLabel("E3:"))
        self.parent_window.ternary_e3_chain = QLineEdit("A")
        self.parent_window.ternary_e3_chain.setFixedWidth(40)
        self.parent_window.ternary_e3_chain.setMinimumHeight(26)
        self.parent_window.ternary_e3_chain.setToolTip("E3 ligase chain ID")
        input_layout.addWidget(self.parent_window.ternary_e3_chain)
        
        # POI Chain
        input_layout.addWidget(QLabel("POI:"))
        self.parent_window.ternary_poi_chain = QLineEdit("B")
        self.parent_window.ternary_poi_chain.setFixedWidth(40)
        self.parent_window.ternary_poi_chain.setMinimumHeight(26)
        self.parent_window.ternary_poi_chain.setToolTip("POI chain ID")
        input_layout.addWidget(self.parent_window.ternary_poi_chain)
        
        # Ligand Residue Name
        input_layout.addWidget(QLabel("Ligand:"))
        self.parent_window.ternary_lig_resn = QLineEdit("UNL")
        self.parent_window.ternary_lig_resn.setFixedWidth(50)
        self.parent_window.ternary_lig_resn.setMinimumHeight(26)
        self.parent_window.ternary_lig_resn.setToolTip("Ligand residue name (e.g., UNL, LIG)")
        input_layout.addWidget(self.parent_window.ternary_lig_resn)
        
        input_layout.addStretch(1)
        
        layout.addWidget(grp_input)
        
        # === 三个子模块标签页 ===
        sub_tabs = QTabWidget()
        sub_tabs.setStyleSheet("""
            QTabWidget::pane { border: 1px solid #e5e7eb; background: white; padding: 4px; }
            QTabBar::tab { padding: 4px 12px; }
            QTabBar::tab:selected { background: #3b82f6; color: white; }
        """)
        sub_tabs.setMaximumHeight(180)
        
        # 1. 界面模块
        interface_tab = self._create_interface_tab()
        sub_tabs.addTab(interface_tab, "Interface")
        
        # 2. 配体模块
        ligand_tab = self._create_ligand_tab()
        sub_tabs.addTab(ligand_tab, "Ligand")
        
        # 3. 三元几何模块
        geometry_tab = self._create_geometry_tab()
        sub_tabs.addTab(geometry_tab, "Geometry")
        
        layout.addWidget(sub_tabs)
        
        # === 运行按钮 ===
        btn_row = QHBoxLayout()
        btn_row.setSpacing(6)
        self.parent_window.ternary_run_all_btn = QPushButton("Run Evaluation")
        self.parent_window.ternary_run_all_btn.setObjectName("primary_btn")
        self.parent_window.ternary_run_all_btn.setMinimumHeight(28)
        self.parent_window.ternary_run_all_btn.clicked.connect(self.run_full_evaluation)
        btn_row.addWidget(self.parent_window.ternary_run_all_btn)
        
        self.parent_window.ternary_visualize_btn = QPushButton("Visualize")
        self.parent_window.ternary_visualize_btn.setObjectName("highlight_btn")
        self.parent_window.ternary_visualize_btn.setMinimumHeight(28)
        self.parent_window.ternary_visualize_btn.clicked.connect(self.visualize_geometry)
        btn_row.addWidget(self.parent_window.ternary_visualize_btn)
        
        export_btn = QPushButton("Export")
        export_btn.setMinimumHeight(28)
        export_btn.clicked.connect(self.export_results)
        btn_row.addWidget(export_btn)
        
        btn_row.addStretch(1)
        layout.addLayout(btn_row)
        
        # === 结果显示 ===
        grp_result = QGroupBox("Results")
        result_layout = QVBoxLayout(grp_result)
        result_layout.setContentsMargins(6, 10, 6, 6)
        self.parent_window.ternary_result_text = QTextEdit()
        self.parent_window.ternary_result_text.setReadOnly(True)
        self.parent_window.ternary_result_text.setMinimumHeight(150)
        self.parent_window.ternary_result_text.setStyleSheet(
            "font-family: 'Consolas', 'Monaco', monospace; font-size: 11px;"
        )
        result_layout.addWidget(self.parent_window.ternary_result_text)
        layout.addWidget(grp_result, 1)

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
        layout = QHBoxLayout(w)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(12)
        
        # 左侧：参数
        left = QVBoxLayout()
        left.setSpacing(4)
        
        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Probe (Å):"))
        self.parent_window.ternary_probe_radius = QDoubleSpinBox()
        self.parent_window.ternary_probe_radius.setRange(0.5, 3.0)
        self.parent_window.ternary_probe_radius.setValue(1.4)
        self.parent_window.ternary_probe_radius.setSingleStep(0.1)
        self.parent_window.ternary_probe_radius.setFixedWidth(70)
        row1.addWidget(self.parent_window.ternary_probe_radius)
        row1.addStretch(1)
        left.addLayout(row1)
        
        row2 = QHBoxLayout()
        row2.addWidget(QLabel("Contact (Å):"))
        self.parent_window.ternary_contact_dist = QDoubleSpinBox()
        self.parent_window.ternary_contact_dist.setRange(3.0, 8.0)
        self.parent_window.ternary_contact_dist.setValue(4.5)
        self.parent_window.ternary_contact_dist.setSingleStep(0.5)
        self.parent_window.ternary_contact_dist.setFixedWidth(70)
        row2.addWidget(self.parent_window.ternary_contact_dist)
        row2.addStretch(1)
        left.addLayout(row2)
        
        run_btn = QPushButton("Calculate")
        run_btn.setMinimumHeight(26)
        run_btn.clicked.connect(self.run_interface_only)
        left.addWidget(run_btn)
        
        layout.addLayout(left)
        
        # 右侧：说明
        desc = QLabel("Calculate BSA and contacts\nbetween E3, POI, and ligand")
        desc.setStyleSheet("color: #64748b; font-size: 11px;")
        layout.addWidget(desc)
        layout.addStretch(1)
        
        return w

    def _create_ligand_tab(self) -> QWidget:
        """创建配体模块子标签页"""
        w = QWidget()
        layout = QHBoxLayout(w)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)
        
        # 属性显示 - 紧凑的网格
        self.lig_prop_labels = {}
        props = ["MW", "LogP", "TPSA", "HBD", "HBA", "RotB", "Fsp3", "Ring"]
        
        grid = QGridLayout()
        grid.setSpacing(4)
        for i, p in enumerate(props):
            grid.addWidget(QLabel(f"{p}:"), i // 4, (i % 4) * 2)
            lbl = QLabel("-")
            lbl.setStyleSheet("font-weight: bold; min-width: 40px;")
            self.lig_prop_labels[p if p != "RotB" else "RotBonds"] = lbl
            self.lig_prop_labels[p if p != "Ring" else "Rings"] = lbl
            grid.addWidget(lbl, i // 4, (i % 4) * 2 + 1)
        
        layout.addLayout(grid)
        
        run_btn = QPushButton("Calculate")
        run_btn.setMinimumHeight(26)
        run_btn.setFixedWidth(80)
        run_btn.clicked.connect(self.run_ligand_only)
        layout.addWidget(run_btn)
        layout.addStretch(1)
        
        return w

    def _create_geometry_tab(self) -> QWidget:
        """创建三元几何模块子标签页"""
        w = QWidget()
        layout = QHBoxLayout(w)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(8)
        
        # 几何参数显示 - 紧凑网格
        self.geom_labels = {}
        features = [
            ("COG", "Å"), ("Angle", "°"), ("E3-POI", "Å"), ("E3-MG", "Å"),
            ("POI-MG", "Å"), ("Coop", "kcal"), ("Hook", ""), ("Dual", "")
        ]
        full_names = ["COG Shift", "Angle", "E3-POI Dist", "E3-MG Dist",
                      "POI-MG Dist", "Cooperativity", "Hook Risk", "Duality"]
        
        grid = QGridLayout()
        grid.setSpacing(4)
        for i, ((name, unit), full) in enumerate(zip(features, full_names)):
            grid.addWidget(QLabel(f"{name}:"), i // 4, (i % 4) * 2)
            lbl = QLabel("-")
            lbl.setStyleSheet("font-weight: bold; min-width: 45px;")
            lbl.setToolTip(f"{full} ({unit})" if unit else full)
            self.geom_labels[full] = lbl
            grid.addWidget(lbl, i // 4, (i % 4) * 2 + 1)
        
        layout.addLayout(grid)
        
        run_btn = QPushButton("Calculate")
        run_btn.setMinimumHeight(26)
        run_btn.setFixedWidth(80)
        run_btn.clicked.connect(self.run_geometry_only)
        layout.addWidget(run_btn)
        layout.addStretch(1)
        
        return w

    def _get_pdb_path(self) -> Optional[str]:
        """从PyMOL对象导出PDB文件"""
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
            QMessageBox.warning(self, "Warning", "Please select a PyMOL object.")
            return
        
        e3 = self.parent_window.ternary_e3_chain.text().strip()
        poi = self.parent_window.ternary_poi_chain.text().strip()
        lig_resn = self.parent_window.ternary_lig_resn.text().strip()
        
        if not all([e3, poi, lig_resn]):
            QMessageBox.warning(self, "Warning", "Please specify E3 chain, POI chain, and ligand residue name.")
            return
        
        self.parent_window.ternary_run_all_btn.setEnabled(False)
        self.parent_window.progress_bar.setVisible(True)
        self.parent_window.progress_bar.setRange(0, 0)
        
        from ..workers_ternary import TernaryEvaluationWorker
        # 使用 ligand residue name 作为配体标识
        obj_name = self.ternary_obj_combo.currentText().strip()
        self.parent_window.ternary_thread = TernaryEvaluationWorker(
            pdb_path, e3, poi, lig_resn, None, None, obj_name
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

    def _safe_get(self, result: Dict[str, Any], key: str, default: float = 0.0) -> float:
        """安全获取字典值，处理 None 的情况"""
        val = result.get(key, default)
        return default if val is None else val
    
    def _format_results(self, result: Dict[str, Any]) -> str:
        """格式化结果 - 包含公式说明和中间计算值"""
        
        # 获取各项数值
        bsa_total = self._safe_get(result, 'bsa_total')
        bsa_mg_e3 = self._safe_get(result, 'bsa_mg_e3')
        bsa_mg_poi = self._safe_get(result, 'bsa_mg_poi')
        bsa_e3_poi = self._safe_get(result, 'bsa_e3_poi')
        contacts = result.get('contact_count_45', 0) or 0
        
        # 获取中间计算值
        sa_e3_mg = self._safe_get(result, 'sa_e3_mg')
        sa_poi_mg = self._safe_get(result, 'sa_poi_mg')
        sa_ternary = self._safe_get(result, 'sa_ternary')
        sa_e3 = self._safe_get(result, 'sa_e3')
        sa_poi = self._safe_get(result, 'sa_poi')
        sa_e3_poi_complex = self._safe_get(result, 'sa_e3_poi_complex')
        contact_mg_e3 = result.get('contact_mg_e3', 0) or 0
        contact_mg_poi = result.get('contact_mg_poi', 0) or 0
        
        cog_shift = self._safe_get(result, 'geom_cog_shift')
        angle = self._safe_get(result, 'geom_angle_deg')
        coop = self._safe_get(result, 'cooperativity_energy')
        hook = self._safe_get(result, 'hook_risk_score')
        duality = self._safe_get(result, 'duality_index')
        
        # 计算 BSA_MG_total
        bsa_mg_total = bsa_mg_e3 + bsa_mg_poi
        total_contact = contact_mg_e3 + contact_mg_poi
        ratio_e3 = contact_mg_e3 / total_contact if total_contact > 0 else 0.5
        ratio_poi = contact_mg_poi / total_contact if total_contact > 0 else 0.5
        
        lines = [
            "=" * 70,
            "  TERNARY COMPLEX EVALUATION RESULTS",
            "=" * 70,
            "",
            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
            "【1. Interface Analysis - BSA (Buried Surface Area)】",
            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
            "",
            "  ▶ MG (Molecular Glue) BSA Calculation:",
            "",
            "    Formula: BSA_MG = SA(E3-MG) + SA(POI-MG) - SA(E3-POI-MG)",
            "",
            "    Surface Areas (measured):",
            f"      SA(E3-MG)      = {sa_e3_mg:>10.1f} Å²  (E3 + Ligand binary complex)",
            f"      SA(POI-MG)     = {sa_poi_mg:>10.1f} Å²  (POI + Ligand binary complex)",
            f"      SA(E3-POI-MG)  = {sa_ternary:>10.1f} Å²  (Ternary complex)",
            "",
            "    Calculation:",
            f"      BSA_MG_total   = {sa_e3_mg:.1f} + {sa_poi_mg:.1f} - {sa_ternary:.1f}",
            f"                     = {bsa_mg_total:>10.1f} Å²",
            "",
            "    Distribution (by contact atoms):",
            f"      MG-E3 contacts = {contact_mg_e3:>5d} atoms  (ratio: {ratio_e3:.2f})",
            f"      MG-POI contacts= {contact_mg_poi:>5d} atoms  (ratio: {ratio_poi:.2f})",
            f"      BSA_MG-E3      = {bsa_mg_total:.1f} × {ratio_e3:.2f} = {bsa_mg_e3:>8.1f} Å²",
            f"      BSA_MG-POI     = {bsa_mg_total:.1f} × {ratio_poi:.2f} = {bsa_mg_poi:>8.1f} Å²",
            "",
            "  ▶ E3-POI BSA Calculation (standard formula):",
            "",
            "    Formula: BSA = (SA_E3 + SA_POI - SA_E3-POI) / 2",
            "",
            "    Surface Areas (measured):",
            f"      SA(E3)         = {sa_e3:>10.1f} Å²  (E3 alone)",
            f"      SA(POI)        = {sa_poi:>10.1f} Å²  (POI alone)",
            f"      SA(E3-POI)     = {sa_e3_poi_complex:>10.1f} Å²  (E3-POI binary complex)",
            "",
            "    Calculation:",
            f"      BSA_E3-POI     = ({sa_e3:.1f} + {sa_poi:.1f} - {sa_e3_poi_complex:.1f}) / 2",
            f"                     = {bsa_e3_poi:>10.1f} Å²",
            "",
            "  ▶ Summary:",
            f"    ┌────────────────────────────────────────────────────────────┐",
            f"    │  Total BSA:      {bsa_total:>8.1f} Å²  (= MG-E3 + MG-POI + E3-POI) │",
            f"    │  MG-E3 BSA:      {bsa_mg_e3:>8.1f} Å²  (Ligand ↔ E3 interface)    │",
            f"    │  MG-POI BSA:     {bsa_mg_poi:>8.1f} Å²  (Ligand ↔ POI interface)   │",
            f"    │  E3-POI BSA:     {bsa_e3_poi:>8.1f} Å²  (E3 ↔ POI direct contact)  │",
            f"    │  Contacts (4.5Å): {contacts:>7d}     (atom pairs within 4.5Å)   │",
            f"    └────────────────────────────────────────────────────────────┘",
            "",
            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
            "【2. Geometry Analysis】",
            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
            "",
            f"    COG Shift:      {cog_shift:>8.2f} Å   (MG deviation from E3-POI axis)",
            f"    Angle:          {angle:>8.1f} °   (E3-MG-POI angle)",
            "",
            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
            "【3. Cooperativity Metrics】",
            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
            "",
            "  Cooperativity Energy (estimated):",
            f"    ΔG_coop:        {coop:>8.2f} kcal/mol",
            "    (Negative = favorable ternary complex formation)",
            "",
            "  Hook Effect Risk:",
            "  ┌─────────────────────────────────────────────────────────┐",
            "  │ Hook_Risk = max(BSA_E3, BSA_POI) / BSA_total           │",
            "  └─────────────────────────────────────────────────────────┘",
            f"    Hook Risk:      {hook:>8.2f}     (0-1, higher = more risk)",
            "",
            "  Duality Index:",
            "  ┌─────────────────────────────────────────────────────────┐",
            "  │ Duality = min(BSA_E3, BSA_POI) / max(BSA_E3, BSA_POI)  │",
            "  └─────────────────────────────────────────────────────────┘",
            f"    Duality:        {duality:>8.2f}     (0-1, higher = more balanced)",
        ]
        
        ligand_mw = self._safe_get(result, 'ligand_mw')
        if ligand_mw > 0:
            lines.extend([
                "",
                "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
                "【4. Ligand Properties】",
                "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
                "",
                f"    MW:             {ligand_mw:>8.1f} Da",
                f"    LogP:           {self._safe_get(result, 'ligand_logp'):>8.2f}",
                f"    TPSA:           {self._safe_get(result, 'ligand_tpsa'):>8.1f} Å²",
                f"    Fsp3:           {self._safe_get(result, 'ligand_fsp3'):>8.2f}",
            ])
        
        lines.extend([
            "",
            "=" * 60,
        ])
        
        return "\n".join(lines)

    def run_interface_only(self):
        """仅运行界面计算"""
        self.log("Running interface calculation...")
        # 简化版：直接调用完整评估
        self.run_full_evaluation()

    def run_ligand_only(self):
        """仅运行配体计算 - 从结构中提取配体信息"""
        lig_resn = self.parent_window.ternary_lig_resn.text().strip()
        if not lig_resn:
            QMessageBox.warning(self, "Warning", "Please specify ligand residue name.")
            return
        
        obj = self.ternary_obj_combo.currentText().strip()
        if not obj or obj == t("no_object"):
            QMessageBox.warning(self, "Warning", "Please select a PyMOL object.")
            return
        
        try:
            from pymol import cmd
            import tempfile
            
            # 只导出配体部分到临时 PDB 文件
            fd, lig_pdb_path = tempfile.mkstemp(suffix=".pdb")
            os.close(fd)
            
            # 选择并保存配体
            lig_sel = f"{obj} and resn {lig_resn}"
            lig_count = cmd.count_atoms(lig_sel)
            
            if lig_count == 0:
                QMessageBox.warning(self, "Warning", f"No atoms found for ligand residue '{lig_resn}'.")
                return
            
            cmd.save(lig_pdb_path, lig_sel)
            self.log(f"Exported ligand ({lig_count} atoms) to temp file")
            
            # 尝试从配体 PDB 提取 SMILES
            smiles = None
            try:
                from openbabel import openbabel as ob
                obConversion = ob.OBConversion()
                obConversion.SetInAndOutFormats("pdb", "smi")
                mol = ob.OBMol()
                obConversion.ReadFile(mol, lig_pdb_path)
                smiles = obConversion.WriteString(mol).strip().split()[0]  # 只取第一个 SMILES
                self.log(f"Extracted SMILES: {smiles[:50]}..." if len(smiles) > 50 else f"Extracted SMILES: {smiles}")
            except ImportError:
                self.log("⚠️ OpenBabel not available, trying RDKit...")
                # 尝试使用 RDKit 从 PDB 读取
                try:
                    from rdkit import Chem
                    mol = Chem.MolFromPDBFile(lig_pdb_path, removeHs=False)
                    if mol:
                        smiles = Chem.MolToSmiles(mol)
                        self.log(f"Extracted SMILES (RDKit): {smiles[:50]}..." if len(smiles) > 50 else f"Extracted SMILES: {smiles}")
                except ImportError:
                    pass
            except Exception as e:
                self.log(f"⚠️ OpenBabel error: {e}")
            
            # 清理临时文件
            try:
                os.remove(lig_pdb_path)
            except:
                pass
            
            if smiles:
                # 计算配体属性
                from ...ternary_complex_evaluator import LigandCalculator
                calc = LigandCalculator()
                props = calc.calculate(smiles)
                
                # 更新 GUI 标签
                attr_map = {
                    "MW": "molecular_weight",
                    "LogP": "logp",
                    "TPSA": "tpsa",
                    "HBD": "hbd_count",
                    "HBA": "hba_count",
                    "RotB": "rotatable_bonds",
                    "RotBonds": "rotatable_bonds",
                    "Fsp3": "fsp3",
                    "Ring": "num_rings",
                    "Rings": "num_rings",
                }
                for key, lbl in self.lig_prop_labels.items():
                    attr_name = attr_map.get(key, key.lower())
                    val = getattr(props, attr_name, 0)
                    lbl.setText(f"{val:.2f}" if isinstance(val, float) else str(val))
                
                self.log(f"✅ Ligand properties: MW={props.molecular_weight:.1f}, LogP={props.logp:.2f}")
            else:
                QMessageBox.information(self, "Info",
                    "Could not extract SMILES from ligand structure.\n"
                    "Please ensure OpenBabel or RDKit is installed.")
                self.log("⚠️ Could not extract SMILES from ligand")
                
        except Exception as e:
            self.on_error(str(e))

    def run_geometry_only(self):
        """仅运行几何计算"""
        self.log("Running geometry calculation...")
        self.run_full_evaluation()

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
            lig_resn = self.parent_window.ternary_lig_resn.text().strip()
            
            cmd.show("cartoon", obj)
            # 使用 PyMOL 标准颜色名称
            cmd.color("slate", f"{obj} and chain {e3}")      # 浅蓝色
            cmd.color("palegreen", f"{obj} and chain {poi}") # 浅绿色
            cmd.show("sticks", f"{obj} and resn {lig_resn}")
            cmd.color("tv_yellow", f"{obj} and resn {lig_resn}")  # 亮黄色
            
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