# -*- coding: utf-8 -*-
"""
Common UI components for GlueTK.
"""
import os
from typing import Optional

try:
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
        QCheckBox, QComboBox, QGroupBox, QGridLayout, QTabWidget, QSpinBox,
        QFileDialog
    )
except ImportError:
    try:
        from PyQt6.QtCore import Qt
        from PyQt6.QtWidgets import (
            QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
            QCheckBox, QComboBox, QGroupBox, QGridLayout, QTabWidget, QSpinBox,
            QFileDialog
        )
    except ImportError:
        raise RuntimeError("PyQt5 or PyQt6 must be installed.")

from ..utils import t

class CommonTab(QWidget):
    """Base class for all workflow tabs"""
    def __init__(self, parent):
        super().__init__()
        self.parent_window = parent
        # Proxy logging to parent
        self.log = parent.log
        self.on_error = parent.on_error
        self.refresh_objects = parent.refresh_objects

    def _create_pocket_detection_card(self) -> QWidget:
        """创建口袋检测与可视化卡片"""
        card = QGroupBox("Pocket Detection")
        layout = QVBoxLayout(card)
        layout.setSpacing(10)
        layout.setContentsMargins(16, 20, 16, 16)
        
        # Target object
        obj_row = QHBoxLayout()
        obj_row.setSpacing(8)
        # We bind to parent's attributes to keep state
        self.parent_window.pocket_obj_combo = QComboBox()
        self.parent_window.pocket_obj_combo.setMinimumHeight(36)
        refresh_btn = QPushButton("Refresh")
        refresh_btn.setObjectName("refresh_btn")
        refresh_btn.setMinimumHeight(36)
        refresh_btn.setToolTip("Refresh objects")
        refresh_btn.clicked.connect(self.refresh_objects)
        obj_row.addWidget(QLabel("Target:"))
        obj_row.addWidget(self.parent_window.pocket_obj_combo, 1)
        obj_row.addWidget(refresh_btn)
        layout.addLayout(obj_row)
        
        # Parameters row
        param_row = QGridLayout()
        param_row.setSpacing(10)
        
        param_row.addWidget(QLabel("Grid Spacing:"), 0, 0)
        self.parent_window.pocket_grid_spacing = QLineEdit("0.5")
        self.parent_window.pocket_grid_spacing.setMinimumHeight(36)
        self.parent_window.pocket_grid_spacing.setMaximumWidth(80)
        param_row.addWidget(self.parent_window.pocket_grid_spacing, 0, 1)
        param_row.addWidget(QLabel("Å"), 0, 2)
        
        param_row.addWidget(QLabel("Min Volume:"), 0, 3)
        self.parent_window.pocket_min_volume = QLineEdit("30")
        self.parent_window.pocket_min_volume.setMinimumHeight(36)
        self.parent_window.pocket_min_volume.setMaximumWidth(80)
        param_row.addWidget(self.parent_window.pocket_min_volume, 0, 4)
        param_row.addWidget(QLabel("Ų"), 0, 5)
        
        param_row.addWidget(QLabel("Color by:"), 1, 0)
        self.parent_window.pocket_color_by = QComboBox()
        self.parent_window.pocket_color_by.addItems(["Volume", "Druggability", "Hydrophobicity", "Depth"])
        self.parent_window.pocket_color_by.setMinimumHeight(36)
        param_row.addWidget(self.parent_window.pocket_color_by, 1, 1, 1, 5)
        
        layout.addLayout(param_row)
        
        # Buttons
        btn_row = QHBoxLayout()
        btn_row.setSpacing(8)
        
        detect_btn = QPushButton("Detect Pockets")
        detect_btn.setObjectName("primary_btn")
        detect_btn.setMinimumHeight(36)
        detect_btn.clicked.connect(self.run_pocket_detection)
        
        viz_btn = QPushButton("Visualize")
        viz_btn.setObjectName("secondary_btn")
        viz_btn.setMinimumHeight(36)
        viz_btn.clicked.connect(self.run_pocket_visualization)
        
        btn_row.addWidget(detect_btn)
        btn_row.addWidget(viz_btn)
        layout.addLayout(btn_row)
        
        return card
    
    def _create_vina_docking_card(self) -> QWidget:
        """创建 Vina 对接卡片（支持自定义盒子）"""
        card = QGroupBox("AutoDock Vina")
        layout = QVBoxLayout(card)
        layout.setSpacing(10)
        layout.setContentsMargins(16, 20, 16, 16)
        
        # Ligand file
        ligand_layout = QHBoxLayout()
        ligand_layout.setSpacing(8)
        
        self.parent_window.vina_ligand = QLineEdit()
        self.parent_window.vina_ligand.setPlaceholderText("Select ligand file (MOL2/SDF/PDBQT)")
        self.parent_window.vina_ligand.setMinimumHeight(36)
        
        ligand_browse = QPushButton("Browse")
        ligand_browse.setObjectName("browse_btn")
        ligand_browse.setMinimumHeight(36)
        ligand_browse.setToolTip("Browse ligand file")
        ligand_browse.clicked.connect(self.browse_vina_ligand)
        
        ligand_layout.addWidget(QLabel("Ligand File:"))
        ligand_layout.addWidget(self.parent_window.vina_ligand, 1)
        ligand_layout.addWidget(ligand_browse)
        layout.addLayout(ligand_layout)
        
        # 参数设置
        param_layout = QGridLayout()
        param_layout.setSpacing(10)
        
        param_layout.addWidget(QLabel("Max Pockets:"), 0, 0)
        self.parent_window.vina_max_pockets = QSpinBox()
        self.parent_window.vina_max_pockets.setRange(1, 10)
        self.parent_window.vina_max_pockets.setValue(3)
        self.parent_window.vina_max_pockets.setMinimumHeight(36)
        self.parent_window.vina_max_pockets.setMaximumWidth(80)
        param_layout.addWidget(self.parent_window.vina_max_pockets, 0, 1)
        
        param_layout.addWidget(QLabel("Exhaustiveness:"), 0, 2)
        self.parent_window.vina_exhaustiveness = QSpinBox()
        self.parent_window.vina_exhaustiveness.setRange(1, 32)
        self.parent_window.vina_exhaustiveness.setValue(8)
        self.parent_window.vina_exhaustiveness.setMinimumHeight(36)
        self.parent_window.vina_exhaustiveness.setMaximumWidth(80)
        param_layout.addWidget(self.parent_window.vina_exhaustiveness, 0, 3)
        
        layout.addLayout(param_layout)
        
        # 自定义盒子选项
        self.parent_window.vina_use_custom_box = QCheckBox("Use custom box (skip pocket detection)")
        layout.addWidget(self.parent_window.vina_use_custom_box)
        
        box_grid = QGridLayout()
        box_grid.setSpacing(8)
        
        # Center
        box_grid.addWidget(QLabel("center_x"), 0, 0)
        self.parent_window.vina_cx = QLineEdit(); self.parent_window.vina_cx.setPlaceholderText("e.g. 10.0"); self.parent_window.vina_cx.setEnabled(False)
        box_grid.addWidget(self.parent_window.vina_cx, 0, 1)
        box_grid.addWidget(QLabel("center_y"), 0, 2)
        self.parent_window.vina_cy = QLineEdit(); self.parent_window.vina_cy.setPlaceholderText("e.g. 20.0"); self.parent_window.vina_cy.setEnabled(False)
        box_grid.addWidget(self.parent_window.vina_cy, 0, 3)
        box_grid.addWidget(QLabel("center_z"), 0, 4)
        self.parent_window.vina_cz = QLineEdit(); self.parent_window.vina_cz.setPlaceholderText("e.g. 30.0"); self.parent_window.vina_cz.setEnabled(False)
        box_grid.addWidget(self.parent_window.vina_cz, 0, 5)
        
        # Size
        box_grid.addWidget(QLabel("size_x"), 1, 0)
        self.parent_window.vina_sx = QLineEdit(); self.parent_window.vina_sx.setPlaceholderText("e.g. 20.0"); self.parent_window.vina_sx.setEnabled(False)
        box_grid.addWidget(self.parent_window.vina_sx, 1, 1)
        box_grid.addWidget(QLabel("size_y"), 1, 2)
        self.parent_window.vina_sy = QLineEdit(); self.parent_window.vina_sy.setPlaceholderText("e.g. 20.0"); self.parent_window.vina_sy.setEnabled(False)
        box_grid.addWidget(self.parent_window.vina_sy, 1, 3)
        box_grid.addWidget(QLabel("size_z"), 1, 4)
        self.parent_window.vina_sz = QLineEdit(); self.parent_window.vina_sz.setPlaceholderText("e.g. 20.0"); self.parent_window.vina_sz.setEnabled(False)
        box_grid.addWidget(self.parent_window.vina_sz, 1, 5)
        
        layout.addLayout(box_grid)
        
        def _toggle_box_fields(checked: bool):
            for w in (self.parent_window.vina_cx, self.parent_window.vina_cy, self.parent_window.vina_cz, 
                      self.parent_window.vina_sx, self.parent_window.vina_sy, self.parent_window.vina_sz):
                w.setEnabled(checked)
        self.parent_window.vina_use_custom_box.toggled.connect(_toggle_box_fields)
        
        # Buttons
        btn_row = QHBoxLayout()
        btn_row.setSpacing(8)
        
        self.parent_window.vina_dock_btn = QPushButton("Run Docking")
        self.parent_window.vina_dock_btn.setObjectName("primary_btn")
        self.parent_window.vina_dock_btn.setMinimumHeight(36)
        self.parent_window.vina_dock_btn.clicked.connect(self.run_vina_docking)
        
        self.parent_window.vina_load_result_btn = QPushButton("Load Result")
        self.parent_window.vina_load_result_btn.setObjectName("secondary_btn")
        self.parent_window.vina_load_result_btn.setMinimumHeight(36)
        self.parent_window.vina_load_result_btn.clicked.connect(self.load_vina_result)
        
        btn_row.addWidget(self.parent_window.vina_dock_btn)
        btn_row.addWidget(self.parent_window.vina_load_result_btn)
        layout.addLayout(btn_row)
        
        return card
    
    def _create_advanced_pocket_card(self) -> QWidget:
        """创建高级口袋分析卡片"""
        card = QGroupBox("Advanced Pocket Analysis")
        layout = QVBoxLayout(card)
        layout.setSpacing(8)
        
        # 创建标签页
        self.parent_window.pocket_advanced_tabs = QTabWidget()
        self.parent_window.pocket_advanced_tabs.setFixedHeight(220)
        
        # Tab 1: 口袋对比
        comparison_tab = self._create_pocket_comparison_tab()
        self.parent_window.pocket_advanced_tabs.addTab(comparison_tab, "Comparison")
        
        # Tab 2: 界面口袋
        interface_tab = self._create_pocket_interface_tab()
        self.parent_window.pocket_advanced_tabs.addTab(interface_tab, "PPI Interface")
        
        # Tab 3: 口袋-相互作用关联
        correlation_tab = self._create_pocket_correlation_tab()
        self.parent_window.pocket_advanced_tabs.addTab(correlation_tab, "Interactions")
        
        # Tab 4: G-motif 口袋
        gmotif_pocket_tab = self._create_gmotif_pocket_tab()
        self.parent_window.pocket_advanced_tabs.addTab(gmotif_pocket_tab, "G-motif")
        
        # 强制设置所有 tab 的背景色（macOS Qt 兼容性）
        bg_color = "#161b22" if getattr(self.parent_window, "_dark_mode", True) else "white"
        for i in range(self.parent_window.pocket_advanced_tabs.count()):
            tab_widget = self.parent_window.pocket_advanced_tabs.widget(i)
            if tab_widget:
                tab_widget.setObjectName("tab_content")
                tab_widget.setStyleSheet(f"#tab_content {{ background-color: {bg_color}; }}")
        
        layout.addWidget(self.parent_window.pocket_advanced_tabs)
        
        return card
    
    def _create_pocket_comparison_tab(self) -> QWidget:
        """创建口袋对比标签页"""
        w = QWidget()
        w.setObjectName("tab_content")
        w.setAutoFillBackground(True)
        layout = QVBoxLayout(w)
        layout.setSpacing(6)
        layout.setContentsMargins(6, 6, 6, 6)
        
        # Object A & B
        obj_grid = QGridLayout()
        obj_grid.setSpacing(6)
        
        self.parent_window.pocket_comp_obj_a = QComboBox()
        self.parent_window.pocket_comp_obj_a.setFixedHeight(28)
        self.parent_window.pocket_comp_obj_b = QComboBox()
        self.parent_window.pocket_comp_obj_b.setFixedHeight(28)
        
        refresh_a = QPushButton("Refresh")
        refresh_a.setObjectName("refresh_btn")
        refresh_a.setMinimumHeight(28)
        refresh_a.clicked.connect(self.refresh_objects)
        
        refresh_b = QPushButton("Refresh")
        refresh_b.setObjectName("refresh_btn")
        refresh_b.setMinimumHeight(28)
        refresh_b.clicked.connect(self.refresh_objects)
        
        obj_grid.addWidget(QLabel("Object A:"), 0, 0)
        obj_grid.addWidget(self.parent_window.pocket_comp_obj_a, 0, 1)
        obj_grid.addWidget(refresh_a, 0, 2)
        obj_grid.addWidget(QLabel("Object B:"), 1, 0)
        obj_grid.addWidget(self.parent_window.pocket_comp_obj_b, 1, 1)
        obj_grid.addWidget(refresh_b, 1, 2)
        layout.addLayout(obj_grid)
        
        # Align checkbox
        self.parent_window.pocket_comp_align = QCheckBox("Align structures before comparison")
        self.parent_window.pocket_comp_align.setChecked(True)
        layout.addWidget(self.parent_window.pocket_comp_align)
        
        # Buttons
        btn_row = QHBoxLayout()
        btn_row.setSpacing(6)
        
        compare_btn = QPushButton("Compare Pockets")
        compare_btn.setObjectName("primary_btn")
        compare_btn.setMinimumHeight(32)
        compare_btn.clicked.connect(self.run_pocket_comparison)
        
        viz_btn = QPushButton("Visualize")
        viz_btn.setObjectName("secondary_btn")
        viz_btn.setMinimumHeight(32)
        viz_btn.clicked.connect(self.visualize_pocket_comparison)
        
        btn_row.addWidget(compare_btn)
        btn_row.addWidget(viz_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)
        
        layout.addStretch()
        return w
    
    def _create_pocket_interface_tab(self) -> QWidget:
        """创建 PPI 界面口袋标签页"""
        w = QWidget()
        w.setObjectName("tab_content")
        w.setAutoFillBackground(True)
        layout = QVBoxLayout(w)
        layout.setSpacing(6)
        layout.setContentsMargins(6, 6, 6, 6)
        
        # Object
        obj_row = QHBoxLayout()
        obj_row.setSpacing(6)
        self.parent_window.pocket_interface_obj = QComboBox()
        self.parent_window.pocket_interface_obj.setFixedHeight(28)
        refresh_btn = QPushButton("Refresh")
        refresh_btn.setObjectName("refresh_btn")
        refresh_btn.setMinimumHeight(28)
        refresh_btn.clicked.connect(self.refresh_objects)
        obj_row.addWidget(QLabel("Object:"))
        obj_row.addWidget(self.parent_window.pocket_interface_obj, 1)
        obj_row.addWidget(refresh_btn)
        layout.addLayout(obj_row)
        
        # Chains
        chain_row = QHBoxLayout()
        chain_row.setSpacing(6)
        self.parent_window.pocket_interface_chain_a = QLineEdit()
        self.parent_window.pocket_interface_chain_a.setPlaceholderText("Chain A")
        self.parent_window.pocket_interface_chain_a.setFixedHeight(28)
        self.parent_window.pocket_interface_chain_a.setFixedWidth(60)
        
        self.parent_window.pocket_interface_chain_b = QLineEdit()
        self.parent_window.pocket_interface_chain_b.setPlaceholderText("Chain B")
        self.parent_window.pocket_interface_chain_b.setFixedHeight(28)
        self.parent_window.pocket_interface_chain_b.setFixedWidth(60)
        
        chain_row.addWidget(QLabel("Chains:"))
        chain_row.addWidget(self.parent_window.pocket_interface_chain_a)
        chain_row.addWidget(QLabel("+"))
        chain_row.addWidget(self.parent_window.pocket_interface_chain_b)
        chain_row.addStretch()
        layout.addLayout(chain_row)
        
        # Buttons
        btn_row = QHBoxLayout()
        btn_row.setSpacing(6)
        
        analyze_btn = QPushButton("Analyze Interface Pockets")
        analyze_btn.setObjectName("primary_btn")
        analyze_btn.setMinimumHeight(32)
        analyze_btn.clicked.connect(self.run_interface_pockets)
        
        btn_row.addWidget(analyze_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)
        
        layout.addStretch()
        return w
    
    def _create_pocket_correlation_tab(self) -> QWidget:
        """创建口袋-相互作用关联标签页"""
        w = QWidget()
        w.setObjectName("tab_content")
        w.setAutoFillBackground(True)
        layout = QVBoxLayout(w)
        layout.setSpacing(6)
        layout.setContentsMargins(6, 6, 6, 6)
        
        # Info label
        info_label = QLabel("First detect pockets, then select interaction CSV")
        info_label.setStyleSheet("color: #64748b; font-size: 11px;")
        layout.addWidget(info_label)
        
        # Interaction CSV
        csv_row = QHBoxLayout()
        csv_row.setSpacing(6)
        self.parent_window.pocket_corr_csv = QLineEdit()
        self.parent_window.pocket_corr_csv.setPlaceholderText("Interaction CSV file")
        self.parent_window.pocket_corr_csv.setFixedHeight(28)
        
        csv_browse = QPushButton("Browse")
        csv_browse.setObjectName("browse_btn")
        csv_browse.setMinimumHeight(28)
        csv_browse.clicked.connect(self.browse_pocket_corr_csv)
        
        csv_row.addWidget(self.parent_window.pocket_corr_csv, 1)
        csv_row.addWidget(csv_browse)
        layout.addLayout(csv_row)
        
        # Buttons
        btn_row = QHBoxLayout()
        btn_row.setSpacing(6)
        
        correlate_btn = QPushButton("Correlate with Pockets")
        correlate_btn.setObjectName("primary_btn")
        correlate_btn.setMinimumHeight(32)
        correlate_btn.clicked.connect(self.run_pocket_correlation)
        
        btn_row.addWidget(correlate_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)
        
        layout.addStretch()
        return w
    
    def _create_gmotif_pocket_tab(self) -> QWidget:
        """创建 G-motif 口袋分析标签页"""
        w = QWidget()
        w.setObjectName("tab_content")
        w.setAutoFillBackground(True)
        layout = QVBoxLayout(w)
        layout.setSpacing(6)
        layout.setContentsMargins(6, 6, 6, 6)
        
        # Object
        obj_row = QHBoxLayout()
        obj_row.setSpacing(6)
        self.parent_window.gmotif_pocket_obj = QComboBox()
        self.parent_window.gmotif_pocket_obj.setFixedHeight(28)
        refresh_btn = QPushButton("R")
        refresh_btn.setObjectName("refresh_btn")
        refresh_btn.setFixedSize(28, 28)
        refresh_btn.clicked.connect(self.refresh_objects)
        obj_row.addWidget(QLabel("Object:"))
        obj_row.addWidget(self.parent_window.gmotif_pocket_obj, 1)
        obj_row.addWidget(refresh_btn)
        layout.addLayout(obj_row)
        
        # Chains
        chain_grid = QGridLayout()
        chain_grid.setSpacing(6)
        
        self.parent_window.gmotif_pocket_e3 = QLineEdit()
        self.parent_window.gmotif_pocket_e3.setPlaceholderText("E3 (e.g., A)")
        self.parent_window.gmotif_pocket_e3.setFixedHeight(28)
        self.parent_window.gmotif_pocket_e3.setFixedWidth(70)
        
        self.parent_window.gmotif_pocket_sub = QLineEdit()
        self.parent_window.gmotif_pocket_sub.setPlaceholderText("Substrate (e.g., B)")
        self.parent_window.gmotif_pocket_sub.setFixedHeight(28)
        self.parent_window.gmotif_pocket_sub.setFixedWidth(70)
        
        self.parent_window.gmotif_pocket_glue = QLineEdit()
        self.parent_window.gmotif_pocket_glue.setPlaceholderText("Glue (optional)")
        self.parent_window.gmotif_pocket_glue.setFixedHeight(28)
        self.parent_window.gmotif_pocket_glue.setFixedWidth(70)
        
        chain_grid.addWidget(QLabel("E3:"), 0, 0)
        chain_grid.addWidget(self.parent_window.gmotif_pocket_e3, 0, 1)
        chain_grid.addWidget(QLabel("Sub:"), 0, 2)
        chain_grid.addWidget(self.parent_window.gmotif_pocket_sub, 0, 3)
        chain_grid.addWidget(QLabel("Glue:"), 1, 0)
        chain_grid.addWidget(self.parent_window.gmotif_pocket_glue, 1, 1)
        layout.addLayout(chain_grid)
        
        # Buttons
        btn_row = QHBoxLayout()
        btn_row.setSpacing(6)
        
        analyze_btn = QPushButton("Comprehensive Analysis")
        analyze_btn.setObjectName("highlight_btn")
        analyze_btn.setFixedHeight(28)
        analyze_btn.clicked.connect(self.run_gmotif_pocket_analysis)
        
        btn_row.addWidget(analyze_btn)
        btn_row.addStretch()
        layout.addLayout(btn_row)
        
        layout.addStretch()
        return w

    def _create_mutation_analysis_card(self) -> QWidget:
        """创建突变分析卡片"""
        card = QGroupBox("Protein Mutation & ΔΔG Analysis")
        layout = QVBoxLayout(card)
        layout.setSpacing(8)
        
        # 输入区域
        input_grid = QGridLayout()
        input_grid.setSpacing(6)
        
        # 对象选择
        input_grid.addWidget(QLabel("Object:"), 0, 0)
        self.parent_window.mut_obj_combo = QComboBox()
        self.parent_window.mut_obj_combo.setFixedHeight(28)
        self.parent_window.mut_refresh_btn = QPushButton("Refresh")
        self.parent_window.mut_refresh_btn.setObjectName("refresh_btn")
        self.parent_window.mut_refresh_btn.setFixedHeight(28)
        self.parent_window.mut_refresh_btn.clicked.connect(self.refresh_objects)
        obj_row = QHBoxLayout()
        obj_row.addWidget(self.parent_window.mut_obj_combo, 1)
        obj_row.addWidget(self.parent_window.mut_refresh_btn)
        input_grid.addLayout(obj_row, 0, 1)
        
        # 突变输入
        input_grid.addWidget(QLabel("Mutations:"), 1, 0)
        self.parent_window.mut_input = QLineEdit()
        self.parent_window.mut_input.setPlaceholderText("e.g., A:23:ALA, B:45:GLY")
        self.parent_window.mut_input.setFixedHeight(28)
        input_grid.addWidget(self.parent_window.mut_input, 1, 1)
        
        # 方法选择
        input_grid.addWidget(QLabel("Method:"), 2, 0)
        method_row = QHBoxLayout()
        self.parent_window.mut_method_combo = QComboBox()
        self.parent_window.mut_method_combo.addItems(["Auto", "FoldX", "PyRosetta"])
        self.parent_window.mut_method_combo.setFixedHeight(28)
        self.parent_window.mut_method_combo.setFixedWidth(120)
        method_row.addWidget(self.parent_window.mut_method_combo)
        method_row.addStretch()
        input_grid.addLayout(method_row, 2, 1)
        
        layout.addLayout(input_grid)
        
        # 按钮行
        btn_row = QHBoxLayout()
        btn_row.setSpacing(6)
        
        self.parent_window.mut_perform_btn = QPushButton("Perform Mutation")
        self.parent_window.mut_perform_btn.setObjectName("primary_btn")
        self.parent_window.mut_perform_btn.setFixedHeight(32)
        self.parent_window.mut_perform_btn.clicked.connect(self.run_mutation)
        
        self.parent_window.mut_minimize_btn = QPushButton("Minimize Energy")
        self.parent_window.mut_minimize_btn.setObjectName("highlight_btn")
        self.parent_window.mut_minimize_btn.setFixedHeight(32)
        self.parent_window.mut_minimize_btn.clicked.connect(self.run_minimize)
        
        self.parent_window.mut_analyze_btn = QPushButton("Full Analysis")
        self.parent_window.mut_analyze_btn.setObjectName("highlight_btn")
        self.parent_window.mut_analyze_btn.setFixedHeight(32)
        self.parent_window.mut_analyze_btn.clicked.connect(self.run_mutation_analysis)
        
        btn_row.addWidget(self.parent_window.mut_perform_btn)
        btn_row.addWidget(self.parent_window.mut_minimize_btn)
        btn_row.addWidget(self.parent_window.mut_analyze_btn)
        btn_row.addStretch()
        
        layout.addLayout(btn_row)
        
        return card
