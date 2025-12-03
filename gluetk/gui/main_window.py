# -*- coding: utf-8 -*-
"""
GlueTK Main Window
Modularized version of the Unified GUI.
"""
import os
import sys
from typing import Optional, List, Dict, Tuple, Any

try:
    from PyQt5.QtCore import Qt, QTimer, QSize
    from PyQt5.QtWidgets import (
        QDialog, QVBoxLayout, QHBoxLayout, QListWidget, QStackedWidget,
        QWidget, QPushButton, QLabel, QFrame, QTextEdit, QProgressBar, QMessageBox
    )
    from PyQt5.QtGui import QIcon, QPixmap
except ImportError:
    try:
        from PyQt6.QtCore import Qt, QTimer, QSize
        from PyQt6.QtWidgets import (
            QDialog, QVBoxLayout, QHBoxLayout, QListWidget, QStackedWidget,
            QWidget, QPushButton, QLabel, QFrame, QTextEdit, QProgressBar, QMessageBox
        )
        from PyQt6.QtGui import QIcon, QPixmap
    except ImportError:
        raise RuntimeError("PyQt5 or PyQt6 must be installed.")

from .utils import t, get_lang
from .tabs.target_discovery import TargetDiscoveryTab
from .tabs.hit_identification import HitIdentificationTab
from .tabs.lead_optimization import LeadOptimizationTab
from .tabs.visualization import VisualizationTab

# Import worker classes if needed for type hinting or global usage
from .workers import AnalysisWorker, GMotifWorker

def _get_logo_path() -> Optional[str]:
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    logo_path = os.path.join(here, "assets", "logo.png")
    if os.path.exists(logo_path):
        return logo_path
    return None

class GlueTKDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(t("title"))
        
        # Window setup
        logo_path = _get_logo_path()
        if logo_path:
            self.setWindowIcon(QIcon(logo_path))
            
        min_w, min_h = 1280, 720
        self.setMinimumSize(min_w, min_h)
        self.resize(min_w, min_h)
        
        # State
        self._dark_mode = True
        self._interactions = []
        self._gmotif_hits = []
        self._esp_maps = {}
        
        # Workers
        self.analysis_thread = None
        self.gmotif_thread = None
        
        # Initialize UI
        self.build_ui()
        self.setup_style()
        
        # Init Timer
        QTimer.singleShot(100, self.refresh_objects)
        
        self.log(t("log_ready"))

    def build_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)
        
        # Top split: Nav | Content
        content_row = QHBoxLayout()
        content_row.setSpacing(0) # Spacing handled by widgets
        
        # Navigation
        nav_widget = self._create_nav_widget()
        content_row.addWidget(nav_widget, 0)
        
        # Content Stack
        self.content_stack = QStackedWidget()
        
        # 0. Welcome
        self.content_stack.addWidget(self.create_welcome_tab())
        
        # 1. Target Discovery
        self.target_tab = TargetDiscoveryTab(self)
        self.content_stack.addWidget(self.target_tab)
        
        # 2. Hit Identification
        self.hit_tab = HitIdentificationTab(self)
        self.content_stack.addWidget(self.hit_tab)
        
        # 3. Lead Optimization
        self.lead_tab = LeadOptimizationTab(self)
        self.content_stack.addWidget(self.lead_tab)
        
        # 4. Visualization
        self.viz_tab = VisualizationTab(self)
        self.content_stack.addWidget(self.viz_tab)
        
        # 5. README
        self.content_stack.addWidget(self.create_readme_tab())
        
        # 6. Contact
        self.content_stack.addWidget(self.create_contact_tab())
        
        content_row.addWidget(self.content_stack, 1)
        main_layout.addLayout(content_row)
        
        # Bottom: Status/Log (Hidden log_edit for compatibility)
        self.log_edit = QTextEdit()
        self.log_edit.setVisible(False)
        main_layout.addWidget(self.log_edit)
        
        # Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        main_layout.addWidget(self.progress_bar)

    def _create_nav_widget(self) -> QWidget:
        nav_widget = QWidget()
        nav_widget.setObjectName("nav_widget")
        self.nav_widget = nav_widget  # 保存引用以便切换主题时更新
        nav_widget.setFixedWidth(240)
        nav_widget.setStyleSheet("background-color: #f1f5f9; border-right: 1px solid #e2e8f0;" if not self._dark_mode else "background-color: #0d1117; border-right: 1px solid #30363d;")
        
        layout = QVBoxLayout(nav_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        
        # Logo Area
        logo_frame = QFrame()
        logo_frame.setFixedHeight(80)
        logo_layout = QVBoxLayout(logo_frame)
        logo_path = _get_logo_path()
        if logo_path:
            logo_lbl = QLabel()
            pix = QPixmap(logo_path)
            logo_lbl.setPixmap(pix.scaled(180, 60, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))
            logo_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            logo_layout.addWidget(logo_lbl)
        layout.addWidget(logo_frame)
        
        # Header (Theme)
        header = QWidget()
        h_layout = QHBoxLayout(header)
        self.theme_btn = QPushButton("☾" if self._dark_mode else "☀")
        self.theme_btn.setFixedSize(30, 30)
        self.theme_btn.setFlat(True)
        self.theme_btn.clicked.connect(self.toggle_theme)
        h_layout.addStretch()
        h_layout.addWidget(self.theme_btn)
        h_layout.setContentsMargins(10, 5, 10, 5)
        layout.addWidget(header)
        
        # Nav List
        self.nav_list = QListWidget()
        self.nav_list.setFrameShape(QFrame.Shape.NoFrame)
        self.nav_list.setStyleSheet("QListWidget { background: transparent; font-size: 14px; } QListWidget::item { padding: 12px; } QListWidget::item:selected { background: #3b82f6; color: white; border-radius: 4px; }")
        
        items = [
            "Welcome",
            "Target Discovery",
            "Hit Identification",
            "Lead Optimization",
            "Visualization"
        ]
        self.nav_list.addItems(items)
        self.nav_list.setCurrentRow(0)
        self.nav_list.currentRowChanged.connect(self.on_nav_changed)
        layout.addWidget(self.nav_list)
        
        layout.addStretch()
        
        # Bottom links
        bottom_bar = QWidget()
        b_layout = QVBoxLayout(bottom_bar)
        b_layout.setSpacing(2)
        
        btn_readme = QPushButton("README")
        btn_readme.setFlat(True)
        btn_readme.clicked.connect(lambda: self.content_stack.setCurrentIndex(5))
        b_layout.addWidget(btn_readme)
        
        btn_contact = QPushButton("Contact Us")
        btn_contact.setFlat(True)
        btn_contact.clicked.connect(lambda: self.content_stack.setCurrentIndex(6))
        b_layout.addWidget(btn_contact)
        
        btn_check = QPushButton("Check Env")
        btn_check.setFlat(True)
        btn_check.clicked.connect(self.check_environment)
        b_layout.addWidget(btn_check)
        
        layout.addWidget(bottom_bar)
        
        return nav_widget

    def on_nav_changed(self, index):
        self.content_stack.setCurrentIndex(index)

    def log(self, msg: str):
        print(f"[GlueTK] {msg}")
        # self.log_edit.append(msg) # Log edit is hidden

    def on_error(self, msg: str):
        self.log(f"Error: {msg}")
        QMessageBox.critical(self, "Error", msg)
        self.progress_bar.setVisible(False)
        # Re-enable buttons if needed
        if hasattr(self, 'gm_btn'): self.gm_btn.setEnabled(True)

    def refresh_objects(self):
        names = []
        try:
            from pymol import cmd
            if hasattr(cmd, "get_names"):
                names = cmd.get_names("objects")
            else:
                names = cmd.get_object_list()
        except:
            names = []
            
        if not names: names = [t("no_object")]
        
        # List of combos to update
        combos = [
            getattr(self, "obj_combo_gm", None),
            getattr(self, "obj_combo_apbs", None),
            getattr(self, "pocket_obj_combo", None),
            getattr(self, "glue_obj_combo", None),
            getattr(self, "tc_obj_combo", None),
            getattr(self, "pl_obj_combo", None),
            getattr(self, "mut_obj_combo", None)
        ]
        
        for cb in combos:
            if cb is not None:
                cb.blockSignals(True)
                cb.clear()
                cb.addItems(names)
                cb.blockSignals(False)
        
        self.log(f"Refreshed {len(names)} objects")

    def update_enablement(self):
        pass # Simplified for now

    def toggle_theme(self):
        self._dark_mode = not self._dark_mode
        self.setup_style()
        self.theme_btn.setText("☾" if self._dark_mode else "☀")
        self.log(f"Switched to {'Dark' if self._dark_mode else 'Light'} mode")

    def setup_style(self):
        if self._dark_mode:
            self.setStyleSheet("""
                QDialog, QWidget { background-color: #0d1117; color: #c9d1d9; }
                QGroupBox { border: 1px solid #30363d; border-radius: 6px; margin-top: 12px; padding-top: 10px; }
                QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left; padding: 0 5px; color: #58a6ff; font-weight: bold; }
                QLineEdit, QComboBox, QSpinBox { background-color: #010409; border: 1px solid #30363d; border-radius: 4px; padding: 4px; color: #c9d1d9; selection-background-color: #1f6feb; }
                QPushButton { background-color: #21262d; border: 1px solid #30363d; border-radius: 4px; padding: 6px 12px; color: #c9d1d9; }
                QPushButton:hover { background-color: #30363d; border-color: #8b949e; }
                QPushButton:pressed { background-color: #282e33; }
                QPushButton#highlight_btn { background-color: #238636; color: #ffffff; border: 1px solid #2ea043; }
                QPushButton#highlight_btn:hover { background-color: #2ea043; }
                QScrollArea { border: none; background-color: transparent; }
                QScrollBar:vertical { border: none; background: #0d1117; width: 10px; margin: 0px 0 0px 0; }
                QScrollBar::handle:vertical { background: #30363d; min-height: 20px; border-radius: 5px; }
                QListWidget { background: transparent; font-size: 14px; color: #c9d1d9; }
                QListWidget::item { padding: 12px; }
                QListWidget::item:selected { background: #3b82f6; color: white; border-radius: 4px; }
                #nav_widget { background-color: #0d1117; border-right: 1px solid #30363d; }
                #tab_content { background-color: #161b22; }
            """)
            # 更新导航栏样式
            if hasattr(self, 'nav_widget'):
                self.nav_widget.setStyleSheet("background-color: #0d1117; border-right: 1px solid #30363d;")
        else:
            # 浅色主题样式
            self.setStyleSheet("""
                QDialog, QWidget { background-color: #ffffff; color: #1f2937; }
                QGroupBox { border: 1px solid #e5e7eb; border-radius: 6px; margin-top: 12px; padding-top: 10px; background-color: #f9fafb; }
                QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left; padding: 0 5px; color: #3b82f6; font-weight: bold; }
                QLineEdit, QComboBox, QSpinBox { background-color: #ffffff; border: 1px solid #d1d5db; border-radius: 4px; padding: 4px; color: #1f2937; selection-background-color: #3b82f6; }
                QPushButton { background-color: #f3f4f6; border: 1px solid #d1d5db; border-radius: 4px; padding: 6px 12px; color: #374151; }
                QPushButton:hover { background-color: #e5e7eb; border-color: #9ca3af; }
                QPushButton:pressed { background-color: #d1d5db; }
                QPushButton#highlight_btn { background-color: #22c55e; color: #ffffff; border: 1px solid #16a34a; }
                QPushButton#highlight_btn:hover { background-color: #16a34a; }
                QPushButton#primary_btn { background-color: #3b82f6; color: #ffffff; border: 1px solid #2563eb; }
                QPushButton#primary_btn:hover { background-color: #2563eb; }
                QScrollArea { border: none; background-color: transparent; }
                QScrollBar:vertical { border: none; background: #f3f4f6; width: 10px; margin: 0px 0 0px 0; }
                QScrollBar::handle:vertical { background: #d1d5db; min-height: 20px; border-radius: 5px; }
                QListWidget { background: transparent; font-size: 14px; color: #1f2937; }
                QListWidget::item { padding: 12px; }
                QListWidget::item:selected { background: #3b82f6; color: white; border-radius: 4px; }
                #nav_widget { background-color: #f1f5f9; border-right: 1px solid #e2e8f0; }
                #tab_content { background-color: #ffffff; }
                QTabWidget::pane { border: 1px solid #e5e7eb; background-color: #ffffff; }
                QTabBar::tab { background-color: #f3f4f6; border: 1px solid #e5e7eb; padding: 6px 12px; }
                QTabBar::tab:selected { background-color: #ffffff; border-bottom-color: #ffffff; }
            """)
            # 更新导航栏样式
            if hasattr(self, 'nav_widget'):
                self.nav_widget.setStyleSheet("background-color: #f1f5f9; border-right: 1px solid #e2e8f0;")

    def check_environment(self):
        try:
            from ..env_checker import ensure_dependencies
            ensure_dependencies(silent=False)
        except ImportError:
            QMessageBox.warning(self, "Error", "Environment checker not found.")

    # --- Simple Pages ---
    def create_welcome_tab(self) -> QWidget:
        w = QWidget()
        layout = QVBoxLayout(w)
        title = QLabel("Welcome to GlueTK")
        title.setStyleSheet("font-size: 24px; font-weight: bold; color: #3b82f6;")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addStretch()
        layout.addWidget(title)
        layout.addWidget(QLabel("Molecular Glue Analysis Plugin for PyMOL", alignment=Qt.AlignmentFlag.AlignCenter))
        layout.addStretch()
        return w

    def create_readme_tab(self) -> QWidget:
        w = QWidget()
        layout = QVBoxLayout(w)
        text = QTextEdit()
        text.setReadOnly(True)
        text.setHtml("<h2>GlueTK Documentation</h2><p>Please refer to QUICK_START.md in the plugin directory.</p>")
        layout.addWidget(text)
        return w

    def create_contact_tab(self) -> QWidget:
        w = QWidget()
        layout = QVBoxLayout(w)
        text = QLabel("Contact Us\n\nEmail: support@gluetk.org\nGitHub: github.com/gluetk")
        text.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(text)
        return w
