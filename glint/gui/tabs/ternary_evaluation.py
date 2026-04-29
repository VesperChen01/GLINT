# -*- coding: utf-8 -*-
"""
Ternary Complex Evaluation Tab
三元复合物评估Label页 - Package含interfaceModule、配体Module、三元几何Module
"""
import os
from typing import Optional, Dict, Any

from ..qt_adapter import (
    Qt, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QComboBox, QGroupBox, QGridLayout, QScrollArea, QFrame, QTabWidget,
    QFileDialog, QMessageBox, QTextEdit, QDoubleSpinBox, QCheckBox
)

from ..utils import t, show_message_box
from .common import CommonTab


class TernaryEvaluationTab(CommonTab):
    """三元复合物评估Label页 - 与Target Discovery并列"""
    
    def __init__(self, parent):
        super().__init__(parent)
        self._last_result = None
        self.init_ui()
        
    def init_ui(self):
        """InitializeUI - 现代卡片式布局"""
        self.setObjectName("scroll_content")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)

        is_dark = getattr(self.parent_window, "_dark_mode", False)
        bg_color = "#161b22" if is_dark else "#ffffff"
        self.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")

        layout = QVBoxLayout(self)
        layout.setSpacing(16)
        layout.setContentsMargins(20, 20, 20, 20)

        # === 页面标题 ===
        header = QHBoxLayout()
        title = QLabel("Ternary Complex Evaluation")
        title.setStyleSheet("""
            font-size: 20px; font-weight: 600;
            color: #3b82f6; padding: 4px 0;
        """)
        header.addWidget(title)
        header.addStretch(1)
        layout.addLayout(header)

        # === 输入卡片 - 两行清晰布局 ===
        input_card = self._create_input_card(is_dark)
        layout.addWidget(input_card)

        # === 三列分析卡片 ===
        analysis_row = QHBoxLayout()
        analysis_row.setSpacing(12)

        # Interface 卡片
        interface_card = self._create_interface_card(is_dark)
        analysis_row.addWidget(interface_card, 1)

        # Ligand 卡片
        ligand_card = self._create_ligand_card(is_dark)
        analysis_row.addWidget(ligand_card, 1)

        # Geometry 卡片
        geometry_card = self._create_geometry_card(is_dark)
        analysis_row.addWidget(geometry_card, 1)

        layout.addLayout(analysis_row)

        # === 操作按钮行 ===
        btn_row = QHBoxLayout()
        btn_row.setSpacing(10)

        self.parent_window.ternary_run_all_btn = QPushButton("Run Full Evaluation")
        self.parent_window.ternary_run_all_btn.setObjectName("primary_btn")
        self.parent_window.ternary_run_all_btn.setMinimumHeight(36)
        self.parent_window.ternary_run_all_btn.setMinimumWidth(180)
        self.parent_window.ternary_run_all_btn.setStyleSheet("""
            QPushButton {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #3b82f6, stop:1 #2563eb);
                color: white; border: none; border-radius: 8px;
                font-weight: 600; font-size: 13px; padding: 8px 20px;
            }
            QPushButton:hover { background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #2563eb, stop:1 #1d4ed8); }
            QPushButton:pressed { background: #1d4ed8; }
            QPushButton:disabled { background: #94a3b8; }
        """)
        self.parent_window.ternary_run_all_btn.clicked.connect(self.run_full_evaluation)
        btn_row.addWidget(self.parent_window.ternary_run_all_btn)

        self.parent_window.ternary_visualize_btn = QPushButton("Visualize")
        self.parent_window.ternary_visualize_btn.setMinimumHeight(36)
        self.parent_window.ternary_visualize_btn.setStyleSheet("""
            QPushButton {
                background: #10b981; color: white; border: none;
                border-radius: 8px; font-weight: 500; padding: 8px 16px;
            }
            QPushButton:hover { background: #059669; }
        """)
        self.parent_window.ternary_visualize_btn.clicked.connect(self.visualize_geometry)
        btn_row.addWidget(self.parent_window.ternary_visualize_btn)

        export_btn = QPushButton("Export")
        export_btn.setMinimumHeight(36)
        export_btn.setStyleSheet("""
            QPushButton {
                background: #6366f1; color: white; border: none;
                border-radius: 8px; font-weight: 500; padding: 8px 16px;
            }
            QPushButton:hover { background: #4f46e5; }
        """)
        export_btn.clicked.connect(self.export_results)
        btn_row.addWidget(export_btn)

        btn_row.addStretch(1)
        layout.addLayout(btn_row)

        # === Results摘要卡片 ===
        summary_card = self._create_summary_card(is_dark)
        layout.addWidget(summary_card)

        # === 详细Results区域 ===
        result_card = self._create_result_card(is_dark)
        layout.addWidget(result_card, 1)

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
            """
        return """
            QFrame {
                background: white;
                border: 1px solid #e2e8f0;
                border-radius: 10px;
                padding: 12px;
                box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            }
        """

    def _create_input_card(self, is_dark: bool) -> QFrame:
        """Create输入卡片 - 两行清晰布局"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(12)
        layout.setContentsMargins(16, 14, 16, 14)

        # 第一行：PyMOL对象Select
        row1 = QHBoxLayout()
        row1.setSpacing(10)

        obj_label = QLabel("PyMOL Object:")
        obj_label.setStyleSheet("font-weight: 500; color: #64748b;" if not is_dark else "font-weight: 500; color: #94a3b8;")
        row1.addWidget(obj_label)

        self.ternary_obj_combo = QComboBox()
        self.ternary_obj_combo.setMinimumHeight(32)
        self.ternary_obj_combo.setMinimumWidth(180)
        self.ternary_obj_combo.setToolTip("Select PyMOL object containing ternary complex")
        self.parent_window.ternary_obj_combo = self.ternary_obj_combo
        row1.addWidget(self.ternary_obj_combo)

        refresh_btn = QPushButton("Refresh")
        refresh_btn.setMinimumHeight(32)
        refresh_btn.setStyleSheet("""
            QPushButton {
                background: #f1f5f9; color: #475569; border: 1px solid #e2e8f0;
                border-radius: 6px; padding: 6px 12px; font-weight: 500;
            }
            QPushButton:hover { background: #e2e8f0; }
        """)
        refresh_btn.clicked.connect(self._do_refresh)
        row1.addWidget(refresh_btn)

        row1.addStretch(1)
        layout.addLayout(row1)

        # 分隔线
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("background: #e2e8f0;" if not is_dark else "background: #30363d;")
        sep.setFixedHeight(1)
        layout.addWidget(sep)

        # 第二行：链和配体Settings
        row2 = QHBoxLayout()
        row2.setSpacing(20)

        # E3 Chain
        e3_group = QHBoxLayout()
        e3_group.setSpacing(6)
        e3_label = QLabel("E3 Chain:")
        e3_label.setStyleSheet("font-weight: 500;")
        e3_group.addWidget(e3_label)
        self.parent_window.ternary_e3_chain = QLineEdit("A")
        self.parent_window.ternary_e3_chain.setFixedWidth(50)
        self.parent_window.ternary_e3_chain.setMinimumHeight(30)
        self.parent_window.ternary_e3_chain.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.parent_window.ternary_e3_chain.setToolTip("E3 ligase chain ID")
        e3_group.addWidget(self.parent_window.ternary_e3_chain)
        row2.addLayout(e3_group)

        # POI Chain
        poi_group = QHBoxLayout()
        poi_group.setSpacing(6)
        poi_label = QLabel("POI Chain:")
        poi_label.setStyleSheet("font-weight: 500;")
        poi_group.addWidget(poi_label)
        self.parent_window.ternary_poi_chain = QLineEdit("B")
        self.parent_window.ternary_poi_chain.setFixedWidth(50)
        self.parent_window.ternary_poi_chain.setMinimumHeight(30)
        self.parent_window.ternary_poi_chain.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.parent_window.ternary_poi_chain.setToolTip("POI chain ID")
        poi_group.addWidget(self.parent_window.ternary_poi_chain)
        row2.addLayout(poi_group)

        # Ligand Residue Name
        lig_group = QHBoxLayout()
        lig_group.setSpacing(6)
        lig_label = QLabel("Ligand Resn:")
        lig_label.setStyleSheet("font-weight: 500;")
        lig_group.addWidget(lig_label)
        self.parent_window.ternary_lig_resn = QLineEdit("UNL")
        self.parent_window.ternary_lig_resn.setFixedWidth(60)
        self.parent_window.ternary_lig_resn.setMinimumHeight(30)
        self.parent_window.ternary_lig_resn.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.parent_window.ternary_lig_resn.setToolTip("Ligand residue name (e.g., UNL, LIG)")
        lig_group.addWidget(self.parent_window.ternary_lig_resn)
        row2.addLayout(lig_group)

        row2.addStretch(1)
        layout.addLayout(row2)

        return card

    def _do_refresh(self):
        """RefreshPyMOL对象列表"""
        names = []
        try:
            from pymol import cmd
            names = cmd.get_names("objects") if hasattr(cmd, "get_names") else cmd.get_object_list()
            self.log(f"PyMOL objects found: {names}")
        except Exception as e:
            self.log(f"PyMOL connection error: {e}")

        if not names:
            names = [t("no_object")]

        self.ternary_obj_combo.blockSignals(True)
        self.ternary_obj_combo.clear()
        self.ternary_obj_combo.addItems(names)
        self.ternary_obj_combo.blockSignals(False)

        self.log(f"Refreshed: {len(names)} objects")

    def _create_interface_card(self, is_dark: bool) -> QFrame:
        """Createinterface分析卡片"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(10)
        layout.setContentsMargins(14, 12, 14, 12)

        # 卡片标题
        title = QLabel("Interface")
        title.setStyleSheet("font-size: 14px; font-weight: 600; color: #3b82f6;")
        layout.addWidget(title)

        desc = QLabel("BSA & Contact Analysis")
        desc.setStyleSheet("font-size: 11px; color: #94a3b8; margin-bottom: 6px;")
        layout.addWidget(desc)

        # Parameters区域
        params = QGridLayout()
        params.setSpacing(8)

        params.addWidget(QLabel("Probe (Å):"), 0, 0)
        self.parent_window.ternary_probe_radius = QDoubleSpinBox()
        self.parent_window.ternary_probe_radius.setRange(0.5, 3.0)
        self.parent_window.ternary_probe_radius.setValue(1.4)
        self.parent_window.ternary_probe_radius.setSingleStep(0.1)
        self.parent_window.ternary_probe_radius.setMinimumHeight(28)
        params.addWidget(self.parent_window.ternary_probe_radius, 0, 1)

        params.addWidget(QLabel("Contact (Å):"), 1, 0)
        self.parent_window.ternary_contact_dist = QDoubleSpinBox()
        self.parent_window.ternary_contact_dist.setRange(3.0, 8.0)
        self.parent_window.ternary_contact_dist.setValue(4.5)
        self.parent_window.ternary_contact_dist.setSingleStep(0.5)
        self.parent_window.ternary_contact_dist.setMinimumHeight(28)
        params.addWidget(self.parent_window.ternary_contact_dist, 1, 1)

        layout.addLayout(params)
        layout.addStretch(1)

        # 计算按钮
        run_btn = QPushButton("Calculate")
        run_btn.setMinimumHeight(30)
        run_btn.setStyleSheet("""
            QPushButton { background: #e0f2fe; color: #0369a1; border: none; border-radius: 6px; font-weight: 500; }
            QPushButton:hover { background: #bae6fd; }
        """)
        run_btn.clicked.connect(self.run_interface_only)
        layout.addWidget(run_btn)

        return card

    def _create_ligand_card(self, is_dark: bool) -> QFrame:
        """Create配体Property卡片"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(10)
        layout.setContentsMargins(14, 12, 14, 12)

        # 卡片标题
        title = QLabel("Ligand")
        title.setStyleSheet("font-size: 14px; font-weight: 600; color: #10b981;")
        layout.addWidget(title)

        desc = QLabel("Molecular Properties")
        desc.setStyleSheet("font-size: 11px; color: #94a3b8; margin-bottom: 6px;")
        layout.addWidget(desc)

        # Property网格
        self.lig_prop_labels = {}
        props = [("MW", "Da"), ("LogP", ""), ("TPSA", "Å²"), ("HBD", ""),
                 ("HBA", ""), ("RotB", ""), ("Fsp3", ""), ("Ring", "")]

        grid = QGridLayout()
        grid.setSpacing(6)
        for i, (p, unit) in enumerate(props):
            lbl_name = QLabel(f"{p}:")
            lbl_name.setStyleSheet("color: #64748b; font-size: 11px;")
            grid.addWidget(lbl_name, i // 2, (i % 2) * 2)

            lbl_val = QLabel("-")
            lbl_val.setStyleSheet("font-weight: 600; font-size: 12px; min-width: 45px;")
            lbl_val.setToolTip(f"{p} ({unit})" if unit else p)
            key = "RotBonds" if p == "RotB" else ("Rings" if p == "Ring" else p)
            self.lig_prop_labels[key] = lbl_val
            grid.addWidget(lbl_val, i // 2, (i % 2) * 2 + 1)

        layout.addLayout(grid)
        layout.addStretch(1)

        # 计算按钮
        run_btn = QPushButton("Calculate")
        run_btn.setMinimumHeight(30)
        run_btn.setStyleSheet("""
            QPushButton { background: #d1fae5; color: #047857; border: none; border-radius: 6px; font-weight: 500; }
            QPushButton:hover { background: #a7f3d0; }
        """)
        run_btn.clicked.connect(self.run_ligand_only)
        layout.addWidget(run_btn)

        return card

    def _create_geometry_card(self, is_dark: bool) -> QFrame:
        """Create几何分析卡片"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(10)
        layout.setContentsMargins(14, 12, 14, 12)

        # 卡片标题
        title = QLabel("Geometry")
        title.setStyleSheet("font-size: 14px; font-weight: 600; color: #8b5cf6;")
        layout.addWidget(title)

        desc = QLabel("Ternary Complex Metrics")
        desc.setStyleSheet("font-size: 11px; color: #94a3b8; margin-bottom: 6px;")
        layout.addWidget(desc)

        # 几何Parameters网格
        self.geom_labels = {}
        features = [
            ("COG", "Å", "COG Shift"), ("Angle", "°", "Angle"),
            ("E3-POI", "Å", "E3-POI Dist"), ("E3-MG", "Å", "E3-MG Dist"),
            ("POI-MG", "Å", "POI-MG Dist"), ("Dual", "", "Duality")
        ]

        grid = QGridLayout()
        grid.setSpacing(6)
        for i, (name, unit, full) in enumerate(features):
            lbl_name = QLabel(f"{name}:")
            lbl_name.setStyleSheet("color: #64748b; font-size: 11px;")
            grid.addWidget(lbl_name, i // 2, (i % 2) * 2)

            lbl_val = QLabel("-")
            lbl_val.setStyleSheet("font-weight: 600; font-size: 12px; min-width: 45px;")
            lbl_val.setToolTip(f"{full} ({unit})" if unit else full)
            self.geom_labels[full] = lbl_val
            grid.addWidget(lbl_val, i // 2, (i % 2) * 2 + 1)

        layout.addLayout(grid)
        layout.addStretch(1)

        # 计算按钮
        run_btn = QPushButton("Calculate")
        run_btn.setMinimumHeight(30)
        run_btn.setStyleSheet("""
            QPushButton { background: #ede9fe; color: #6d28d9; border: none; border-radius: 6px; font-weight: 500; }
            QPushButton:hover { background: #ddd6fe; }
        """)
        run_btn.clicked.connect(self.run_geometry_only)
        layout.addWidget(run_btn)

        return card

    def _create_summary_card(self, is_dark: bool) -> QFrame:
        """CreateResults摘要卡片"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(10)
        layout.setContentsMargins(16, 14, 16, 14)

        # 标题
        title = QLabel("Quick Summary")
        title.setStyleSheet("font-size: 14px; font-weight: 600; color: #f59e0b;")
        layout.addWidget(title)

        # 摘要指标网格
        self.summary_labels = {}
        metrics = [
            ("Total BSA", "Å²", "#3b82f6"),
            ("MG-E3 BSA", "Å²", "#3b82f6"),
            ("MG-POI BSA", "Å²", "#10b981"),
            ("E3-POI BSA", "Å²", "#8b5cf6"),
            ("Contacts", "", "#f59e0b"),
        ]

        grid = QHBoxLayout()
        grid.setSpacing(16)

        for name, unit, color in metrics:
            metric_box = QVBoxLayout()
            metric_box.setSpacing(2)

            val_lbl = QLabel("-")
            val_lbl.setStyleSheet(f"font-size: 18px; font-weight: 700; color: {color};")
            val_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.summary_labels[name] = val_lbl
            metric_box.addWidget(val_lbl)

            name_lbl = QLabel(f"{name}" + (f" ({unit})" if unit else ""))
            name_lbl.setStyleSheet("font-size: 10px; color: #94a3b8;")
            name_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            metric_box.addWidget(name_lbl)

            grid.addLayout(metric_box)

        layout.addLayout(grid)

        return card

    def _create_result_card(self, is_dark: bool) -> QFrame:
        """Create详细Results卡片"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(8)
        layout.setContentsMargins(14, 12, 14, 12)

        # 标题行
        header = QHBoxLayout()
        title = QLabel("Detailed Results")
        title.setStyleSheet("font-size: 14px; font-weight: 600; color: #64748b;")
        header.addWidget(title)
        header.addStretch(1)
        layout.addLayout(header)

        # Results文本区域
        self.parent_window.ternary_result_text = QTextEdit()
        self.parent_window.ternary_result_text.setReadOnly(True)
        self.parent_window.ternary_result_text.setMinimumHeight(180)

        text_style = """
            QTextEdit {
                font-family: 'Consolas', 'Courier New', 'SF Mono', 'Monaco', monospace;
                font-size: 12px;
                background: #0d1117;
                color: #c9d1d9;
                border: 1px solid #30363d;
                border-radius: 6px;
                padding: 10px;
                selection-background-color: #1f6feb;
                selection-color: #ffffff;
            }
        """ if is_dark else """
            QTextEdit {
                font-family: 'Consolas', 'Courier New', 'SF Mono', 'Monaco', monospace;
                font-size: 12px;
                background: #ffffff;
                color: #1e293b;
                border: 1px solid #e2e8f0;
                border-radius: 6px;
                padding: 10px;
                selection-background-color: #3b82f6;
                selection-color: #ffffff;
            }
        """
        self.parent_window.ternary_result_text.setStyleSheet(text_style)
        self.parent_window.ternary_result_text.setPlaceholderText(
            "Results will appear here after running evaluation...\n\n"
            "Click 'Run Full Evaluation' to analyze the ternary complex."
        )
        layout.addWidget(self.parent_window.ternary_result_text)

        return card

    def _get_pdb_path(self) -> Optional[str]:
        """从PyMOL对象ExportPDBFile"""
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
            show_message_box(self, "Warning", "Please select a PyMOL object.", "warning")
            return

        e3 = self.parent_window.ternary_e3_chain.text().strip()
        poi = self.parent_window.ternary_poi_chain.text().strip()
        lig_resn = self.parent_window.ternary_lig_resn.text().strip()

        if not all([e3, poi, lig_resn]):
            show_message_box(self, "Warning", "Please specify E3 chain, POI chain, and ligand residue name.", "warning")
            return
        
        self.parent_window.ternary_run_all_btn.setEnabled(False)
        self.parent_window.progress_bar.setVisible(True)
        self.parent_window.progress_bar.setRange(0, 0)
        
        from ..workers_ternary import TernaryEvaluationWorker
        # using ligand residue name 作为配体ID
        obj_name = self.ternary_obj_combo.currentText().strip()
        self.parent_window.ternary_thread = TernaryEvaluationWorker(
            pdb_path, e3, poi, lig_resn, None, None, obj_name
        )
        self.parent_window.ternary_thread.progress.connect(self.log)
        self.parent_window.ternary_thread.error.connect(self.on_error)
        self.parent_window.ternary_thread.finished.connect(self._on_finished)
        self.parent_window.ternary_thread.start()

    def _on_finished(self, result: Dict[str, Any], out_csv: str):
        """评估Completed"""
        self._last_result = result
        self.parent_window.progress_bar.setVisible(False)
        self.parent_window.progress_bar.setRange(0, 1)
        self.parent_window.ternary_run_all_btn.setEnabled(True)

        # UpdateResultsDisplay
        self.parent_window.ternary_result_text.setText(self._format_results(result))

        # Update摘要卡片
        summary_map = {
            "Total BSA": ("bsa_total", "Å²"),
            "MG-E3 BSA": ("bsa_mg_e3", "Å²"),
            "MG-POI BSA": ("bsa_mg_poi", "Å²"),
            "E3-POI BSA": ("bsa_e3_poi", "Å²"),
            "Contacts": ("contact_count_45", ""),
        }
        for name, (key, unit) in summary_map.items():
            if name in self.summary_labels:
                val = self._safe_get(result, key) if key != "contact_count_45" else (result.get(key, 0) or 0)
                if isinstance(val, float):
                    self.summary_labels[name].setText(f"{val:.1f}")
                else:
                    self.summary_labels[name].setText(str(val))

        # Update配体PropertyLabel
        for key, lbl in self.lig_prop_labels.items():
            map_key = f"ligand_{key.lower()}"
            if key == "RotBonds":
                map_key = "ligand_rotatable_bonds"
            val = result.get(map_key, 0)
            lbl.setText(f"{val:.2f}" if isinstance(val, float) else str(val))

        # Update几何Label
        geom_map = {
            "COG Shift": "geom_cog_shift",
            "Angle": "geom_angle_deg",
            "E3-POI Dist": "dist_e3_poi",
            "E3-MG Dist": "dist_e3_mg",
            "POI-MG Dist": "dist_poi_mg",
            "Duality": "duality_index"
        }
        for name, key in geom_map.items():
            val = result.get(key, 0)
            self.geom_labels[name].setText(f"{val:.2f}" if isinstance(val, float) else str(val))

        self.log("✅ Full evaluation complete")

    def _safe_get(self, result: Dict[str, Any], key: str, default: float = 0.0) -> float:
        """安全获取字典Value，处理 None 的情况"""
        val = result.get(key, default)
        return default if val is None else val

    def _parse_numeric_text(self, text: str) -> Optional[float]:
        """Parse a numeric label value from the GUI."""
        value = (text or "").strip()
        if not value or value == "-":
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def _snapshot_current_result(self) -> Dict[str, Any]:
        """Merge the latest GUI-visible values into a result dict for export/display."""
        result = dict(self._last_result or {})

        summary_map = {
            "Total BSA": "bsa_total",
            "MG-E3 BSA": "bsa_mg_e3",
            "MG-POI BSA": "bsa_mg_poi",
            "E3-POI BSA": "bsa_e3_poi",
            "Contacts": "contact_count_45",
        }
        for label_name, key in summary_map.items():
            label = self.summary_labels.get(label_name)
            numeric = self._parse_numeric_text(label.text()) if label is not None else None
            if numeric is None:
                continue
            result[key] = int(round(numeric)) if key == "contact_count_45" else numeric

        ligand_map = {
            "MW": "ligand_mw",
            "LogP": "ligand_logp",
            "TPSA": "ligand_tpsa",
            "HBD": "ligand_hbd",
            "HBA": "ligand_hba",
            "RotBonds": "ligand_rotatable_bonds",
            "Fsp3": "ligand_fsp3",
            "Rings": "ligand_rings",
        }
        for label_name, key in ligand_map.items():
            label = self.lig_prop_labels.get(label_name)
            numeric = self._parse_numeric_text(label.text()) if label is not None else None
            if numeric is None:
                continue
            if key in {"ligand_hbd", "ligand_hba", "ligand_rotatable_bonds", "ligand_rings"}:
                result[key] = int(round(numeric))
            else:
                result[key] = numeric

        geom_map = {
            "COG Shift": "geom_cog_shift",
            "Angle": "geom_angle_deg",
            "E3-POI Dist": "dist_e3_poi",
            "E3-MG Dist": "dist_e3_mg",
            "POI-MG Dist": "dist_poi_mg",
            "Duality": "duality_index",
        }
        for label_name, key in geom_map.items():
            label = self.geom_labels.get(label_name)
            numeric = self._parse_numeric_text(label.text()) if label is not None else None
            if numeric is None:
                continue
            result[key] = numeric
            if key == "duality_index":
                result["balance_index"] = numeric

        return result
    
    def _format_results(self, result: Dict[str, Any]) -> str:
        """格式化Results - Package含公式说明和中间计算Value"""
        
        # 获取各项数Value
        bsa_total = self._safe_get(result, 'bsa_total')
        bsa_mg_e3 = self._safe_get(result, 'bsa_mg_e3')
        bsa_mg_poi = self._safe_get(result, 'bsa_mg_poi')
        bsa_e3_poi = self._safe_get(result, 'bsa_e3_poi')
        contacts = result.get('contact_count_45', 0) or 0
        
        # 获取中间计算Value
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
            "  * MG (Molecular Glue) BSA Calculation:",
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
            "  * E3-POI BSA Calculation (standard formula):",
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
            "  * Summary:",
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
            "【3. Quality Metrics】",
            "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━",
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
        """仅运行interface计算"""
        self.log("Running interface calculation...")
        # 简化版：直接调用完整评估
        self.run_full_evaluation()

    def run_ligand_only(self):
        """仅运行配体计算 - 从结构中提取配体Information"""
        lig_resn = self.parent_window.ternary_lig_resn.text().strip()
        if not lig_resn:
            show_message_box(self, "Warning", "Please specify ligand residue name.", "warning")
            return

        obj = self.ternary_obj_combo.currentText().strip()
        if not obj or obj == t("no_object"):
            show_message_box(self, "Warning", "Please select a PyMOL object.", "warning")
            return

        try:
            from pymol import cmd
            import tempfile

            # 只Export配体部分到临时 PDB File
            fd, lig_pdb_path = tempfile.mkstemp(suffix=".pdb")
            os.close(fd)

            # Select并Save配体
            lig_sel = f"{obj} and resn {lig_resn}"
            lig_count = cmd.count_atoms(lig_sel)

            if lig_count == 0:
                show_message_box(self, "Warning", f"No atoms found for ligand residue '{lig_resn}'.", "warning")
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
                smiles = obConversion.WriteString(mol).strip().split()[0]  # 只取First SMILES
                self.log(f"Extracted SMILES: {smiles[:50]}..." if len(smiles) > 50 else f"Extracted SMILES: {smiles}")
            except ImportError:
                self.log("⚠️ OpenBabel not available, trying RDKit...")
                # 尝试using RDKit 从 PDB 读取
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
            
            # 清理临时File
            try:
                os.remove(lig_pdb_path)
            except OSError:  # 临时FileDelete可能Failed
                pass
            
            if smiles:
                # 计算配体Property
                from ...ternary_complex_evaluator import LigandCalculator
                calc = LigandCalculator()
                props = calc.calculate(smiles)
                
                # Update GUI Label
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

                if self._last_result is None:
                    self._last_result = {}
                self._last_result.update({
                    "ligand_mw": props.molecular_weight,
                    "ligand_logp": props.logp,
                    "ligand_tpsa": props.tpsa,
                    "ligand_hbd": props.hbd_count,
                    "ligand_hba": props.hba_count,
                    "ligand_rotatable_bonds": props.rotatable_bonds,
                    "ligand_fsp3": props.fsp3,
                    "ligand_rings": props.num_rings,
                })
                if self.parent_window.ternary_result_text.toPlainText().strip():
                    self.parent_window.ternary_result_text.setText(self._format_results(self._snapshot_current_result()))
                
                self.log(f"✅ Ligand properties: MW={props.molecular_weight:.1f}, LogP={props.logp:.2f}")
            else:
                show_message_box(self, "Info",
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
            show_message_box(self, "Warning", "Please run evaluation first.", "warning")
            return
        
        try:
            from pymol import cmd, cgo
            
            obj = self.ternary_obj_combo.currentText().strip()
            e3 = self.parent_window.ternary_e3_chain.text().strip()
            poi = self.parent_window.ternary_poi_chain.text().strip()
            lig_resn = self.parent_window.ternary_lig_resn.text().strip()
            
            cmd.show("cartoon", obj)
            # using PyMOL 标准颜色Name
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
        """ExportResults到CSV - 与GUI摘要卡片内容保持一致"""
        # 使用GUI摘要卡片中的显示值
        export_data = []

        # 导出摘要卡片数据（与GUI显示顺序一致）
        summary_items = [
            ("Total BSA", "Å²"),
            ("MG-E3 BSA", "Å²"),
            ("MG-POI BSA", "Å²"),
            ("E3-POI BSA", "Å²"),
            ("Contacts", ""),
        ]

        for name, unit in summary_items:
            label = self.summary_labels.get(name)
            if label:
                value = label.text()
                if value and value != "-":
                    export_data.append((name + (f" ({unit})" if unit else ""), value))

        # 导出配体性质
        ligand_items = [
            ("MW", "Da"),
            ("LogP", ""),
            ("TPSA", "Å²"),
            ("HBD", ""),
            ("HBA", ""),
            ("RotBonds", ""),
            ("Fsp3", ""),
            ("Rings", ""),
        ]

        for name, unit in ligand_items:
            label = self.lig_prop_labels.get(name)
            if label:
                value = label.text()
                if value and value != "-":
                    export_data.append((f"Ligand {name}" + (f" ({unit})" if unit else ""), value))

        # 导出几何参数
        geom_items = [
            ("COG Shift", "Å"),
            ("Angle", "°"),
            ("E3-POI Dist", "Å"),
            ("E3-MG Dist", "Å"),
            ("POI-MG Dist", "Å"),
            ("Duality", ""),
        ]

        for name, unit in geom_items:
            label = self.geom_labels.get(name)
            if label:
                value = label.text()
                if value and value != "-":
                    export_data.append((name + (f" ({unit})" if unit else ""), value))

        if not export_data:
            show_message_box(self, "Warning", "No results to export. Please run evaluation first.", "warning")
            return

        fn, _ = QFileDialog.getSaveFileName(self, "Save Results", "", "CSV (*.csv)")
        if fn:
            import csv
            from datetime import datetime

            with open(fn, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                # 写入标题
                writer.writerow(['GLINT Ternary Complex Evaluation Results'])
                writer.writerow(['Export Time:', datetime.now().strftime('%Y-%m-%d %H:%M:%S')])
                writer.writerow([])

                # 写入摘要信息
                writer.writerow(['=== Quick Summary ==='])
                for name, value in export_data:
                    writer.writerow([name, value])

                writer.writerow([])

                # 添加完整结果文本（如果存在）
                full_text = self.parent_window.ternary_result_text.toPlainText()
                if full_text and full_text.strip():
                    writer.writerow(['=== Detailed Results ==='])
                    # 将多行文本按行写入
                    for line in full_text.split('\n'):
                        writer.writerow([line])

            self.log(f"✅ Results exported to {fn}")
