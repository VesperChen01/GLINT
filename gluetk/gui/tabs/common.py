# -*- coding: utf-8 -*-
"""
Common UI components for GlueTK.

This module provides base classes and shared UI components used across
all workflow tabs in the GlueTK GUI.
"""
import os
from typing import Optional, Callable, Any, TYPE_CHECKING

if TYPE_CHECKING:
    from ..main_window import GlueTKDialog

from ..qt_adapter import (
    Qt, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton,
    QCheckBox, QComboBox, QGroupBox, QGridLayout, QTabWidget, QSpinBox,
    QFileDialog, QMessageBox
)

from ..utils import t


class CommonTab(QWidget):
    """
    Base class for all workflow tabs in GlueTK.
    
    Provides common functionality including:
    - Access to parent window for logging and error handling
    - Shared UI component creation methods
    - Common file browsing utilities
    
    Attributes:
        parent_window: Reference to the main GlueTKDialog window
        log: Logging function from parent
        on_error: Error handling function from parent
        refresh_objects: Function to refresh PyMOL objects
    """
    
    def __init__(self, parent: 'GlueTKDialog') -> None:
        """
        Initialize the CommonTab.
        
        Args:
            parent: The parent GlueTKDialog window
        """
        super().__init__()
        self.parent_window = parent
        # Proxy logging to parent
        self.log: Callable[[str], None] = parent.log
        self.on_error: Callable[[str], None] = parent.on_error
        self.refresh_objects: Callable[[], None] = parent.refresh_objects

    def _create_pocket_detection_card(self) -> QWidget:
        """Create pocket detection and visualization card"""
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
        """Create Vina docking card (supports custom box)"""
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
        
        # Parameter settings
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
        
        # Custom box options
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
    
    def _create_mutation_analysis_card(self) -> QWidget:
        """Create mutation analysis card"""
        card = QGroupBox("Protein Mutation & ΔΔG Analysis")
        layout = QVBoxLayout(card)
        layout.setSpacing(8)
        
        # Input area
        input_grid = QGridLayout()
        input_grid.setSpacing(6)
        
        # Object selection
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
        
        # Mutation input
        input_grid.addWidget(QLabel("Mutations:"), 1, 0)
        self.parent_window.mut_input = QLineEdit()
        self.parent_window.mut_input.setPlaceholderText("e.g., A:23:ALA, B:45:GLY")
        self.parent_window.mut_input.setFixedHeight(28)
        input_grid.addWidget(self.parent_window.mut_input, 1, 1)
        
        # Method selection
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
        
        # Button row
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

    # --- Common File Browsing Methods ---
    
    def _browse_file(self, line_edit: QLineEdit, file_filter: str, title: str = "Open File") -> Optional[str]:
        """
        Browse and select a file to open.
        
        Args:
            line_edit: The QLineEdit to populate with the selected path
            file_filter: File filter string (e.g., "PDB Files (*.pdb)")
            title: Dialog title
            
        Returns:
            The selected file path, or None if cancelled
        """
        fn, _ = QFileDialog.getOpenFileName(self, title, "", file_filter)
        if fn:
            line_edit.setText(fn)
        return fn if fn else None
    
    def _browse_save_file(self, line_edit: QLineEdit, file_filter: str, title: str = "Save File") -> Optional[str]:
        """
        Browse and select a file path for saving.
        
        Args:
            line_edit: The QLineEdit to populate with the selected path
            file_filter: File filter string (e.g., "CSV Files (*.csv)")
            title: Dialog title
            
        Returns:
            The selected file path, or None if cancelled
        """
        fn, _ = QFileDialog.getSaveFileName(self, title, "", file_filter)
        if fn:
            line_edit.setText(fn)
        return fn if fn else None
    
    def _browse_directory(self, line_edit: QLineEdit, title: str = "Select Directory") -> Optional[str]:
        """
        Browse and select a directory.
        
        Args:
            line_edit: The QLineEdit to populate with the selected path
            title: Dialog title
            
        Returns:
            The selected directory path, or None if cancelled
        """
        folder = QFileDialog.getExistingDirectory(self, title)
        if folder:
            line_edit.setText(folder)
        return folder if folder else None

    # --- Placeholder Methods for Subclasses ---
    # These should be overridden in subclasses that use the corresponding cards
    
    def run_pocket_detection(self) -> None:
        """Run pocket detection. Override in subclass."""
        self.log("Pocket detection not implemented in this tab.")
    
    def run_pocket_visualization(self) -> None:
        """Visualize detected pockets. Override in subclass."""
        self.log("Pocket visualization not implemented in this tab.")
    
    def browse_vina_ligand(self) -> None:
        """Browse for Vina ligand file. Override in subclass."""
        pass
    
    def run_vina_docking(self) -> None:
        """Run Vina docking. Override in subclass."""
        self.log("Vina docking not implemented in this tab.")
    
    def load_vina_result(self) -> None:
        """Load Vina docking result. Override in subclass."""
        pass
    
    def run_mutation(self) -> None:
        """Perform mutation. Override in subclass."""
        self.log("Mutation not implemented in this tab.")
    
    def run_minimize(self) -> None:
        """Run energy minimization. Override in subclass."""
        self.log("Energy minimization not implemented in this tab.")
    
    def run_mutation_analysis(self) -> None:
        """Run full mutation analysis. Override in subclass."""
        self.log("Mutation analysis not implemented in this tab.")
