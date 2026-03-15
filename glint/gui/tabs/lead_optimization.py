# -*- coding: utf-8 -*-
"""
Lead Optimization: PPI, Glue, Ternary, Mutation, Electrostatic Complementarity
"""
import os
import traceback
from typing import Optional

from ..qt_adapter import (
    Qt, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QCheckBox, QComboBox, QGroupBox, QGridLayout, QScrollArea, QFrame,
    QFileDialog, QMessageBox, QThread, Signal as pyqtSignal
)

from ..utils import t, show_message_box, show_question_box
from .common import CommonTab


class _Diagram2DWorker(QThread):
    """后台线程：生成 2D 相互作用图，避免阻塞 Qt 主线程。

    注意：PyMOL cmd 不是线程安全的，因此所有 PyMOL 操作必须在主线程中Completed。
    本 Worker 通过 pdb_file（主线程预先Export的临时 PDB File）获取配体数据，
    避免在后台线程中调用 PyMOL，防止 Qt 事件循环死锁。

    超时保护：总超时 90 秒，防止 RDKit/matplotlib 在复杂分子上无限期运行。

    Signals:
        finished(str): 生成Success时发射，携带输出FilePath
        error(str): 生成Failed时发射，携带ErrorInformation
    """
    finished = pyqtSignal(str)
    error = pyqtSignal(str)

    # 总超时Time（秒）
    TOTAL_TIMEOUT = 90

    def __init__(self, csv_path: str, ligand_resname: str,
                 pdb_file: str, output_path: str, min_confidence: float):
        super().__init__()
        self.csv_path = csv_path
        self.ligand_resname = ligand_resname
        self.pdb_file = pdb_file          # 主线程预Export的配体 PDB FilePath
        self.output_path = output_path
        self.min_confidence = min_confidence
        self._timed_out = False

    def run(self) -> None:
        """在后台线程中执行 generate_2d_interaction_diagram，带总超时保护。"""
        import threading

        result_holder = {'path': None, 'error': None}

        def _generate():
            """实际的图表生成逻辑，在子线程中运行。"""
            try:
                # 确保 matplotlib using非交互式后端
                import matplotlib
                matplotlib.use('Agg')

                try:
                    from ...interaction_2d_plot import generate_2d_interaction_diagram
                except ImportError:
                    from interaction_2d_plot import generate_2d_interaction_diagram

                # using pdb_file 而非 obj_name，确保后台线程不调用 PyMOL
                final_path = generate_2d_interaction_diagram(
                    csv_path=self.csv_path,
                    ligand_resname=self.ligand_resname,
                    pdb_file=self.pdb_file,
                    output_path=self.output_path,
                    min_confidence=self.min_confidence
                )
                result_holder['path'] = final_path

            except Exception as e:
                import traceback
                traceback.print_exc()
                result_holder['error'] = str(e)

        # 在独立线程中运行图表生成，带超时保护
        gen_thread = threading.Thread(target=_generate, daemon=True)
        gen_thread.start()
        gen_thread.join(timeout=self.TOTAL_TIMEOUT)

        if gen_thread.is_alive():
            # 超时：线程仍在运行
            self._timed_out = True
            print(f"[2D Diagram] ⚠️ 图表生成超时 ({self.TOTAL_TIMEOUT}s)，强制中止")
            self.error.emit(
                f"图表生成超时 ({self.TOTAL_TIMEOUT}s)。\n"
                "可能原因：配体结构过于复杂或 Open Babel 未响应。\n"
                "请尝试简化配体或Confirm Open Babel 已正确安装。"
            )
        elif result_holder['error']:
            self.error.emit(result_holder['error'])
        elif result_holder['path'] and os.path.exists(result_holder['path']):
            self.finished.emit(result_holder['path'])
        else:
            self.error.emit("Failed to generate diagram. See log for details.")

        # 清理主线程预Export的临时 PDB File
        try:
            if self.pdb_file and os.path.exists(self.pdb_file):
                os.remove(self.pdb_file)
        except OSError:
            pass


class LeadOptimizationTab(CommonTab):
    def __init__(self, parent):
        super().__init__(parent)

        # Remove proxies for methods implemented here
        for attr in ['run_mutation', 'run_minimize', 'run_mutation_analysis']:
            if attr in self.__dict__:
                del self.__dict__[attr]

        self.init_ui()

    def init_ui(self):
        """InitializeUI - 现代卡片式布局"""
        self.setObjectName("scroll_content")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.parent_window._lead_scroll_content = self

        is_dark = getattr(self.parent_window, "_dark_mode", False)
        bg_color = "#161b22" if is_dark else "#f8fafc"
        self.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")

        layout = QVBoxLayout(self)
        layout.setSpacing(16)
        layout.setContentsMargins(20, 20, 20, 20)

        # === 页面标题 ===
        header = QHBoxLayout()
        title = QLabel("Lead Optimization")
        title.setStyleSheet("""
            font-size: 20px; font-weight: 600;
            color: #3b82f6; padding: 4px 0;
        """)
        header.addWidget(title)
        header.addStretch(1)
        layout.addLayout(header)

        # 1. Protein-Protein Interface (PPI) Analysis
        grp_ppi = QFrame()
        grp_ppi.setStyleSheet(self._get_card_style(is_dark))
        ppi_layout = QVBoxLayout(grp_ppi)
        ppi_layout.setSpacing(12)
        ppi_layout.setContentsMargins(16, 14, 16, 14)

        ppi_title = QLabel("Protein-Protein Interface (PPI) Analysis")
        ppi_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        ppi_layout.addWidget(ppi_title)

        ppi_grid = QGridLayout()
        ppi_grid.setContentsMargins(12, 16, 12, 8)
        ppi_grid.setColumnStretch(1, 1); ppi_grid.setColumnStretch(3, 1)
        ppi_grid.setHorizontalSpacing(12); ppi_grid.setVerticalSpacing(4)
        
        
        ppi_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_obj_combo = QComboBox(); self.parent_window.ppi_obj_combo.setMinimumHeight(32)
        self.parent_window.ppi_refresh_btn = QPushButton(t("refresh")); self.parent_window.ppi_refresh_btn.setMinimumHeight(32); self.parent_window.ppi_refresh_btn.clicked.connect(self.refresh_objects)
        r0 = QHBoxLayout(); r0.addWidget(self.parent_window.ppi_obj_combo, 1); r0.addWidget(self.parent_window.ppi_refresh_btn)
        ppi_grid.addLayout(r0, 0, 1)
        
        ppi_grid.addWidget(QLabel("Protein1 Chains:"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_protein1_chains = QLineEdit(); self.parent_window.ppi_protein1_chains.setPlaceholderText("e.g. A"); self.parent_window.ppi_protein1_chains.setMinimumHeight(32)
        ppi_grid.addWidget(self.parent_window.ppi_protein1_chains, 0, 3)
        
        ppi_grid.addWidget(QLabel("Protein2 Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_protein2_chains = QLineEdit(); self.parent_window.ppi_protein2_chains.setPlaceholderText("e.g. B"); self.parent_window.ppi_protein2_chains.setMinimumHeight(32)
        ppi_grid.addWidget(self.parent_window.ppi_protein2_chains, 1, 1)
        
        ppi_grid.addWidget(QLabel("Interface Dist (Å):"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_interface_dist = QLineEdit("4.5"); self.parent_window.ppi_interface_dist.setMinimumHeight(32)
        ppi_grid.addWidget(self.parent_window.ppi_interface_dist, 1, 3)
        
        ppi_grid.addWidget(QLabel("Output CSV:"), 2, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_csv = QLineEdit(); self.parent_window.ppi_csv.setPlaceholderText("Optional"); self.parent_window.ppi_csv.setMinimumHeight(32)
        self.parent_window.ppi_csv_btn = QPushButton(t("browse")); self.parent_window.ppi_csv_btn.setMinimumHeight(32)
        self.parent_window.ppi_csv_btn.clicked.connect(lambda: self._browse_save_file(self.parent_window.ppi_csv, "CSV (*.csv)"))
        r2 = QHBoxLayout(); r2.addWidget(self.parent_window.ppi_csv, 1); r2.addWidget(self.parent_window.ppi_csv_btn)
        ppi_grid.addLayout(r2, 2, 1, 1, 3)

        # Row 3: Visualization Options
        ppi_grid.addWidget(QLabel("Display Mode:"), 3, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ppi_display_mode = QComboBox()
        self.parent_window.ppi_display_mode.setMinimumHeight(32)
        self.parent_window.ppi_display_mode.addItems(["Surface + Interaction", "Cartoon + Interaction"])
        ppi_grid.addWidget(self.parent_window.ppi_display_mode, 3, 1)
        
        self.parent_window.ppi_show_residue_labels = QCheckBox("Show Residue Labels")
        self.parent_window.ppi_show_residue_labels.setChecked(True)
        ppi_grid.addWidget(self.parent_window.ppi_show_residue_labels, 3, 2)
        
        self.parent_window.ppi_show_distance_labels = QCheckBox("Show Distance Labels")
        self.parent_window.ppi_show_distance_labels.setChecked(False)
        ppi_grid.addWidget(self.parent_window.ppi_show_distance_labels, 3, 3)


        self.parent_window.ppi_show_hydrophobic = QCheckBox("Show Hydrophobic Interactions")
        self.parent_window.ppi_show_hydrophobic.setChecked(False)
        ppi_grid.addWidget(self.parent_window.ppi_show_hydrophobic, 4, 1)
        ppi_layout.addLayout(ppi_grid)

        ppi_btn_row = QHBoxLayout()
        ppi_btn_row.setSpacing(10)

        self.parent_window.ppi_analyze_btn = QPushButton("Analyze PPI Interface")
        self.parent_window.ppi_analyze_btn.setMinimumHeight(36)
        self.parent_window.ppi_analyze_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.ppi_analyze_btn.clicked.connect(self.run_ppi_analysis)
        ppi_btn_row.addWidget(self.parent_window.ppi_analyze_btn)
        ppi_btn_row.addStretch(1)

        layout.addWidget(grp_ppi)
        layout.addLayout(ppi_btn_row)

        # 2. Protein-Ligand Interactions
        grp_pl = QFrame()
        grp_pl.setStyleSheet(self._get_card_style(is_dark))
        pl_layout = QVBoxLayout(grp_pl)
        pl_layout.setSpacing(12)
        pl_layout.setContentsMargins(16, 14, 16, 14)

        pl_title = QLabel("Protein-Ligand Interactions")
        pl_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        pl_layout.addWidget(pl_title)

        pl_grid = QGridLayout()
        pl_grid.setContentsMargins(12, 16, 12, 8)
        pl_grid.setColumnStretch(1, 1); pl_grid.setColumnStretch(3, 1)
        pl_grid.setHorizontalSpacing(12); pl_grid.setVerticalSpacing(4)
        
        
        pl_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_obj_combo = QComboBox(); self.parent_window.pl_obj_combo.setMinimumWidth(150); self.parent_window.pl_obj_combo.setMinimumHeight(32)
        self.parent_window.pl_refresh_btn = QPushButton(t("refresh")); self.parent_window.pl_refresh_btn.setMinimumHeight(32); self.parent_window.pl_refresh_btn.clicked.connect(self.refresh_objects)
        r0_pl = QHBoxLayout(); r0_pl.addWidget(self.parent_window.pl_obj_combo, 1); r0_pl.addWidget(self.parent_window.pl_refresh_btn)
        pl_grid.addLayout(r0_pl, 0, 1)
        
        pl_grid.addWidget(QLabel("Ligand Name:"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_ligand_name = QLineEdit(); self.parent_window.pl_ligand_name.setPlaceholderText("e.g. CC885, Auto-detect if blank"); self.parent_window.pl_ligand_name.setMinimumHeight(32)
        pl_grid.addWidget(self.parent_window.pl_ligand_name, 0, 3)
        
        pl_grid.addWidget(QLabel("Protein Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_protein_chains = QLineEdit(); self.parent_window.pl_protein_chains.setPlaceholderText("e.g. A,B (optional)"); self.parent_window.pl_protein_chains.setMinimumHeight(32)
        pl_grid.addWidget(self.parent_window.pl_protein_chains, 1, 1)
        
        pl_grid.addWidget(QLabel("Distance (Å):"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_distance = QLineEdit("4.5"); self.parent_window.pl_distance.setMinimumHeight(32)
        pl_grid.addWidget(self.parent_window.pl_distance, 1, 3)
        
        # Output CSV moved to row 3 to make space for 3D options
        # See below for new layout positioning
        
        pl_btn_row = QHBoxLayout()
        pl_btn_row.setSpacing(10)

        self.parent_window.pl_analyze_btn = QPushButton("Analyze Protein-Ligand")
        self.parent_window.pl_analyze_btn.setMinimumHeight(36)
        self.parent_window.pl_analyze_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.pl_analyze_btn.clicked.connect(self.run_pl_analysis)

        self.parent_window.pl_2d_btn = QPushButton("Generate 2D Diagram")
        self.parent_window.pl_2d_btn.setMinimumHeight(36)
        self.parent_window.pl_2d_btn.setStyleSheet(self._get_secondary_btn_style())
        self.parent_window.pl_2d_btn.clicked.connect(self.run_pl_2d_diagram)
        
        pl_btn_row.addWidget(self.parent_window.pl_analyze_btn)
        pl_btn_row.addWidget(self.parent_window.pl_2d_btn)
        pl_btn_row.addStretch(1)
        
        layout.addWidget(grp_pl)
        layout.addLayout(pl_btn_row)
        
        # 3D Visualization Options
        pl_grid.addWidget(QLabel("Min Confidence:"), 2, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_min_confidence = QComboBox(); self.parent_window.pl_min_confidence.setMinimumHeight(32)
        self.parent_window.pl_min_confidence.addItems(["0.0 (Show All)", "0.5", "0.6", "0.7", "0.8 (High)", "0.9"])
        self.parent_window.pl_min_confidence.setCurrentText("0.8 (High)")
        pl_grid.addWidget(self.parent_window.pl_min_confidence, 2, 1)

        self.parent_window.pl_show_hydrophobic = QCheckBox("Show Hydrophobic Interactions")
        self.parent_window.pl_show_hydrophobic.setChecked(False)
        pl_grid.addWidget(self.parent_window.pl_show_hydrophobic, 2, 3)
        
        # Distance labels option
        pl_grid.addWidget(QLabel("3D Display:"), 3, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pl_show_distance_labels = QCheckBox("Show Distance Labels")
        self.parent_window.pl_show_distance_labels.setChecked(False)
        self.parent_window.pl_show_distance_labels.setToolTip("Display distance values on interaction lines in 3D view")
        pl_grid.addWidget(self.parent_window.pl_show_distance_labels, 3, 1)
        
        # Output CSV row
        pl_grid.addWidget(QLabel("Output CSV:"), 4, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        
        self.parent_window.pl_csv = QLineEdit(); self.parent_window.pl_csv.setPlaceholderText("Optional"); self.parent_window.pl_csv.setMinimumHeight(32)
        self.parent_window.pl_csv_btn = QPushButton(t("browse")); self.parent_window.pl_csv_btn.setMinimumHeight(32)
        self.parent_window.pl_csv_btn.clicked.connect(lambda: self._browse_save_file(self.parent_window.pl_csv, "CSV (*.csv)"))
        
        r2_pl = QHBoxLayout(); r2_pl.addWidget(self.parent_window.pl_csv, 1); r2_pl.addWidget(self.parent_window.pl_csv_btn)
        pl_grid.addLayout(r2_pl, 4, 1, 1, 3)

        pl_layout.addLayout(pl_grid)

        # 3. Ligand-Ligand Interactions
        grp_ll = QFrame()
        grp_ll.setStyleSheet(self._get_card_style(is_dark))
        ll_layout = QVBoxLayout(grp_ll)
        ll_layout.setSpacing(12)
        ll_layout.setContentsMargins(16, 14, 16, 14)

        ll_title = QLabel("Ligand-Ligand Interactions (Small Molecule - Small Molecule)")
        ll_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        ll_layout.addWidget(ll_title)

        ll_grid = QGridLayout()
        ll_grid.setContentsMargins(12, 16, 12, 8)
        ll_grid.setColumnStretch(1, 1); ll_grid.setColumnStretch(3, 1)
        ll_grid.setHorizontalSpacing(12); ll_grid.setVerticalSpacing(4)
        
        
        ll_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ll_obj_combo = QComboBox(); self.parent_window.ll_obj_combo.setMinimumWidth(150); self.parent_window.ll_obj_combo.setMinimumHeight(32)
        self.parent_window.ll_refresh_btn = QPushButton(t("refresh")); self.parent_window.ll_refresh_btn.setMinimumHeight(32); self.parent_window.ll_refresh_btn.clicked.connect(self.refresh_objects)
        r0_ll = QHBoxLayout(); r0_ll.addWidget(self.parent_window.ll_obj_combo, 1); r0_ll.addWidget(self.parent_window.ll_refresh_btn)
        ll_grid.addLayout(r0_ll, 0, 1, 1, 3)
        
        ll_grid.addWidget(QLabel("Selection 1:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ll_sel1 = QLineEdit(); self.parent_window.ll_sel1.setPlaceholderText("e.g. resn LIG1 or resi 100"); self.parent_window.ll_sel1.setMinimumHeight(32)
        ll_grid.addWidget(self.parent_window.ll_sel1, 1, 1)
        
        ll_grid.addWidget(QLabel("Selection 2:"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ll_sel2 = QLineEdit(); self.parent_window.ll_sel2.setPlaceholderText("e.g. resn LIG2 or resi 200"); self.parent_window.ll_sel2.setMinimumHeight(32)
        ll_grid.addWidget(self.parent_window.ll_sel2, 1, 3)
        
        ll_grid.addWidget(QLabel("Distance (Å):"), 2, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ll_dist = QLineEdit("4.5"); self.parent_window.ll_dist.setMinimumHeight(32)
        ll_grid.addWidget(self.parent_window.ll_dist, 2, 1)
        
        ll_grid.addWidget(QLabel("Output CSV:"), 2, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ll_csv = QLineEdit(); self.parent_window.ll_csv.setPlaceholderText("Optional"); self.parent_window.ll_csv.setMinimumHeight(32)
        self.parent_window.ll_csv_btn = QPushButton(t("browse")); self.parent_window.ll_csv_btn.setMinimumHeight(32)
        self.parent_window.ll_csv_btn.clicked.connect(lambda: self._browse_save_file(self.parent_window.ll_csv, "CSV (*.csv)"))
        r2_ll = QHBoxLayout(); r2_ll.addWidget(self.parent_window.ll_csv, 1); r2_ll.addWidget(self.parent_window.ll_csv_btn)
        ll_grid.addLayout(r2_ll, 2, 3)

        ll_layout.addLayout(ll_grid)

        ll_btn_row = QHBoxLayout()
        ll_btn_row.setSpacing(10)

        self.parent_window.ll_analyze_btn = QPushButton("Analyze Ligand-Ligand")
        self.parent_window.ll_analyze_btn.setMinimumHeight(36)
        self.parent_window.ll_analyze_btn.setStyleSheet(self._get_primary_btn_style())
        self.parent_window.ll_analyze_btn.clicked.connect(self.run_ll_analysis)
        ll_btn_row.addWidget(self.parent_window.ll_analyze_btn)
        ll_btn_row.addStretch(1)

        layout.addWidget(grp_ll)
        layout.addLayout(ll_btn_row)

        # 4. Electrostatic Complementarity (EC) Analysis
        grp_ec = QFrame()
        grp_ec.setStyleSheet(self._get_card_style(is_dark))
        ec_layout = QVBoxLayout(grp_ec)
        ec_layout.setSpacing(12)
        ec_layout.setContentsMargins(16, 14, 16, 14)

        ec_title = QLabel("Electrostatic Complementarity (EC) Analysis")
        ec_title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        ec_layout.addWidget(ec_title)

        ec_grid = QGridLayout()
        ec_grid.setContentsMargins(12, 16, 12, 8)
        ec_grid.setColumnStretch(1, 1); ec_grid.setColumnStretch(3, 1)
        ec_grid.setHorizontalSpacing(12); ec_grid.setVerticalSpacing(4)
        
        
        ec_grid.addWidget(QLabel("Target Object:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_obj_combo = QComboBox(); self.parent_window.ec_obj_combo.setMinimumWidth(150); self.parent_window.ec_obj_combo.setMinimumHeight(32)
        self.parent_window.ec_refresh_btn = QPushButton(t("refresh")); self.parent_window.ec_refresh_btn.setMinimumHeight(32); self.parent_window.ec_refresh_btn.clicked.connect(self.refresh_objects)
        r0_ec = QHBoxLayout(); r0_ec.addWidget(self.parent_window.ec_obj_combo, 1); r0_ec.addWidget(self.parent_window.ec_refresh_btn)
        ec_grid.addLayout(r0_ec, 0, 1)
        
        ec_grid.addWidget(QLabel("Ligand/Glue Name:"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.ec_ligand_name = QLineEdit(); self.parent_window.ec_ligand_name.setPlaceholderText("e.g. LIG, CC885"); self.parent_window