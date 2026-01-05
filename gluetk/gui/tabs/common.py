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
    QFileDialog, QMessageBox, QFrame
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

    def _get_card_style_common(self, is_dark: bool) -> str:
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

    def _create_pocket_detection_card(self) -> QWidget:
        """Create pocket detection and visualization card"""
        card = QGroupBox("Pocket Detection")
        layout = QVBoxLayout(card)
        layout.setSpacing(14)
        layout.setContentsMargins(12, 16, 12, 8)
        
        # Target object
        obj_row = QHBoxLayout()
        obj_row.setSpacing(8)
        # We bind to parent's attributes to keep state
        self.parent_window.pocket_obj_combo = QComboBox()
        self.parent_window.pocket_obj_combo.setMinimumHeight(32)
        refresh_btn = QPushButton("Refresh")
        refresh_btn.setObjectName("refresh_btn")
        refresh_btn.setMinimumHeight(32)
        refresh_btn.setToolTip("Refresh objects")
        refresh_btn.clicked.connect(self.refresh_objects)
        obj_row.addWidget(QLabel("Target:"))
        obj_row.addWidget(self.parent_window.pocket_obj_combo, 1)
        obj_row.addWidget(refresh_btn)
        layout.addLayout(obj_row)
        
        # Parameters row
        param_row = QGridLayout()
        param_row.setSpacing(12)
        param_row.setVerticalSpacing(4)
        
        
        param_row.addWidget(QLabel("Grid Spacing:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pocket_grid_spacing = QLineEdit("0.5")
        self.parent_window.pocket_grid_spacing.setMinimumHeight(32)
        self.parent_window.pocket_grid_spacing.setMaximumWidth(80)
        param_row.addWidget(self.parent_window.pocket_grid_spacing, 0, 1)
        param_row.addWidget(QLabel("Å"), 0, 2)
        
        param_row.addWidget(QLabel("Min Volume:"), 0, 3, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pocket_min_volume = QLineEdit("30")
        self.parent_window.pocket_min_volume.setMinimumHeight(32)
        self.parent_window.pocket_min_volume.setMaximumWidth(80)
        param_row.addWidget(self.parent_window.pocket_min_volume, 0, 4)
        param_row.addWidget(QLabel("Ų"), 0, 5)
        
        param_row.addWidget(QLabel("Color by:"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.pocket_color_by = QComboBox()
        self.parent_window.pocket_color_by.addItems(["Volume", "Druggability", "Hydrophobicity", "Depth"])
        self.parent_window.pocket_color_by.setMinimumHeight(32)
        param_row.addWidget(self.parent_window.pocket_color_by, 1, 1, 1, 5)
        
        layout.addLayout(param_row)
        
        # Buttons
        btn_row = QHBoxLayout()
        btn_row.setSpacing(8)
        
        detect_btn = QPushButton("Detect Pockets")
        detect_btn.setObjectName("primary_btn")
        detect_btn.setMinimumHeight(32)
        detect_btn.clicked.connect(self.run_pocket_detection)
        
        viz_btn = QPushButton("Visualize")
        viz_btn.setObjectName("secondary_btn")
        viz_btn.setMinimumHeight(32)
        viz_btn.clicked.connect(self.run_pocket_visualization)
        
        btn_row.addWidget(detect_btn)
        btn_row.addWidget(viz_btn)
        layout.addLayout(btn_row)
        
        return card
    
    def _create_vina_docking_card(self) -> QWidget:
        """Create Vina docking card (custom box only, no pocket detection)"""
        card = QGroupBox("AutoDock Vina")
        layout = QVBoxLayout(card)
        layout.setSpacing(14)
        layout.setContentsMargins(12, 16, 12, 8)
        
        # Receptor object selection
        rec_layout = QHBoxLayout()
        rec_layout.setSpacing(8)
        self.parent_window.vina_receptor_combo = QComboBox()
        self.parent_window.vina_receptor_combo.setMinimumHeight(32)
        rec_refresh = QPushButton("Refresh")
        rec_refresh.setObjectName("refresh_btn")
        rec_refresh.setMinimumHeight(32)
        rec_refresh.clicked.connect(self.refresh_objects)
        rec_layout.addWidget(QLabel("Receptor:"))
        rec_layout.addWidget(self.parent_window.vina_receptor_combo, 1)
        rec_layout.addWidget(rec_refresh)
        layout.addLayout(rec_layout)
        
        # Ligand file
        ligand_layout = QHBoxLayout()
        ligand_layout.setSpacing(8)
        
        self.parent_window.vina_ligand = QLineEdit()
        self.parent_window.vina_ligand.setPlaceholderText("Select ligand file (MOL2/SDF/PDBQT)")
        self.parent_window.vina_ligand.setMinimumHeight(32)
        
        ligand_browse = QPushButton("Browse")
        ligand_browse.setObjectName("browse_btn")
        ligand_browse.setMinimumHeight(32)
        ligand_browse.setToolTip("Browse ligand file")
        ligand_browse.clicked.connect(self.browse_vina_ligand)
        
        ligand_layout.addWidget(QLabel("Ligand File:"))
        ligand_layout.addWidget(self.parent_window.vina_ligand, 1)
        ligand_layout.addWidget(ligand_browse)
        layout.addLayout(ligand_layout)
        
        # Output directory
        out_layout = QHBoxLayout()
        out_layout.setSpacing(8)
        self.parent_window.vina_output_dir = QLineEdit()
        self.parent_window.vina_output_dir.setPlaceholderText("Output directory (optional)")
        self.parent_window.vina_output_dir.setMinimumHeight(32)
        out_browse = QPushButton("Browse")
        out_browse.setObjectName("browse_btn")
        out_browse.setMinimumHeight(32)
        out_browse.clicked.connect(lambda: self._browse_directory(self.parent_window.vina_output_dir, "Select Output Directory"))
        out_layout.addWidget(QLabel("Output Dir:"))
        out_layout.addWidget(self.parent_window.vina_output_dir, 1)
        out_layout.addWidget(out_browse)
        layout.addLayout(out_layout)

        # Parameter settings
        param_layout = QGridLayout()
        param_layout.setSpacing(12)
        param_layout.setVerticalSpacing(4)
        
        param_layout.addWidget(QLabel("Exhaustiveness:"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.vina_exhaustiveness = QSpinBox()
        self.parent_window.vina_exhaustiveness.setRange(1, 32)
        self.parent_window.vina_exhaustiveness.setValue(8)
        self.parent_window.vina_exhaustiveness.setMinimumHeight(32)
        self.parent_window.vina_exhaustiveness.setMaximumWidth(80)
        param_layout.addWidget(self.parent_window.vina_exhaustiveness, 0, 1)
        
        param_layout.addWidget(QLabel("Num Modes:"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.vina_num_modes = QSpinBox()
        self.parent_window.vina_num_modes.setRange(1, 20)
        self.parent_window.vina_num_modes.setValue(9)
        self.parent_window.vina_num_modes.setMinimumHeight(32)
        self.parent_window.vina_num_modes.setMaximumWidth(80)
        param_layout.addWidget(self.parent_window.vina_num_modes, 0, 3)
        
        layout.addLayout(param_layout)
        
        # Box parameters (always enabled)
        box_label = QLabel("Docking Box (use 'Get from Selection' to auto-fill center)")
        layout.addWidget(box_label)
        
        box_grid = QGridLayout()
        box_grid.setSpacing(8)
        box_grid.setVerticalSpacing(4)
        
        # Center
        box_grid.addWidget(QLabel("center_x"), 0, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.vina_cx = QLineEdit()
        self.parent_window.vina_cx.setPlaceholderText("e.g. 10.0")
        self.parent_window.vina_cx.setMinimumHeight(32)
        box_grid.addWidget(self.parent_window.vina_cx, 0, 1)
        
        box_grid.addWidget(QLabel("center_y"), 0, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.vina_cy = QLineEdit()
        self.parent_window.vina_cy.setPlaceholderText("e.g. 20.0")
        self.parent_window.vina_cy.setMinimumHeight(32)
        box_grid.addWidget(self.parent_window.vina_cy, 0, 3)
        
        box_grid.addWidget(QLabel("center_z"), 0, 4, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.vina_cz = QLineEdit()
        self.parent_window.vina_cz.setPlaceholderText("e.g. 30.0")
        self.parent_window.vina_cz.setMinimumHeight(32)
        box_grid.addWidget(self.parent_window.vina_cz, 0, 5)
        
        # Size
        box_grid.addWidget(QLabel("size_x"), 1, 0, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.vina_sx = QLineEdit("20.0")
        self.parent_window.vina_sx.setMinimumHeight(32)
        box_grid.addWidget(self.parent_window.vina_sx, 1, 1)
        
        box_grid.addWidget(QLabel("size_y"), 1, 2, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.vina_sy = QLineEdit("20.0")
        self.parent_window.vina_sy.setMinimumHeight(32)
        box_grid.addWidget(self.parent_window.vina_sy, 1, 3)
        
        box_grid.addWidget(QLabel("size_z"), 1, 4, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.parent_window.vina_sz = QLineEdit("20.0")
        self.parent_window.vina_sz.setMinimumHeight(32)
        box_grid.addWidget(self.parent_window.vina_sz, 1, 5)
        
        layout.addLayout(box_grid)
        
        # Get center from selection
        sel_layout = QHBoxLayout()
        sel_layout.setSpacing(8)
        self.parent_window.vina_selection = QLineEdit()
        self.parent_window.vina_selection.setPlaceholderText("PyMOL selection (e.g. resn LIG)")
        self.parent_window.vina_selection.setMinimumHeight(32)
        get_center_btn = QPushButton("Get Center from Selection")
        get_center_btn.setObjectName("secondary_btn")
        get_center_btn.setMinimumHeight(32)
        get_center_btn.clicked.connect(self.get_center_from_selection)
        sel_layout.addWidget(self.parent_window.vina_selection, 1)
        sel_layout.addWidget(get_center_btn)
        layout.addLayout(sel_layout)
        
        # Buttons
        btn_row = QHBoxLayout()
        btn_row.setSpacing(8)
        
        self.parent_window.vina_dock_btn = QPushButton("Run Docking")
        self.parent_window.vina_dock_btn.setObjectName("primary_btn")
        self.parent_window.vina_dock_btn.setMinimumHeight(32)
        self.parent_window.vina_dock_btn.clicked.connect(self.run_vina_docking)
        
        self.parent_window.vina_load_result_btn = QPushButton("Load Result")
        self.parent_window.vina_load_result_btn.setObjectName("secondary_btn")
        self.parent_window.vina_load_result_btn.setMinimumHeight(32)
        self.parent_window.vina_load_result_btn.clicked.connect(self.load_vina_result)
        
        btn_row.addWidget(self.parent_window.vina_dock_btn)
        btn_row.addWidget(self.parent_window.vina_load_result_btn)
        layout.addLayout(btn_row)
        
        return card
    
    def _create_mutation_analysis_card(self, is_dark: bool = False) -> QWidget:
        """Create mutation analysis card"""
        card = QFrame()
        card.setStyleSheet(self._get_card_style_common(is_dark))
        layout = QVBoxLayout(card)
        layout.setSpacing(12)
        layout.setContentsMargins(16, 14, 16, 14)

        title = QLabel("Protein Mutation & ΔΔG Analysis")
        title.setStyleSheet("font-size: 15px; font-weight: 600; color: #1e293b; padding-bottom: 4px;" if not is_dark else "font-size: 15px; font-weight: 600; color: #e2e8f0; padding-bottom: 4px;")
        layout.addWidget(title)
        
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
        
        # Method selection (FoldX only)
        input_grid.addWidget(QLabel("Method:"), 2, 0)
        method_row = QHBoxLayout()
        self.parent_window.mut_method_combo = QComboBox()
        self.parent_window.mut_method_combo.addItems(["FoldX"])
        self.parent_window.mut_method_combo.setFixedHeight(28)
        self.parent_window.mut_method_combo.setFixedWidth(120)
        self.parent_window.mut_method_combo.setToolTip("ΔΔG calculation requires FoldX\nDownload: https://foldxsuite.crg.eu/")
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
    
    def get_center_from_selection(self) -> None:
        """Get box center from PyMOL selection using centerofmass. Override in subclass."""
        self.log("Get center from selection not implemented in this tab.")
    
    def run_mutation(self) -> None:
        """Perform mutation. Override in subclass."""
        self.log("Mutation not implemented in this tab.")
    
    def run_minimize(self) -> None:
        """Run energy minimization. Override in subclass."""
        self.log("Energy minimization not implemented in this tab.")
    
    def run_mutation_analysis(self) -> None:
        """Run full mutation analysis. Override in subclass."""
        self.log("Mutation analysis not implemented in this tab.")
