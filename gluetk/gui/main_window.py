# -*- coding: utf-8 -*-
"""
GlueTK Main Window
Modularized version of the Unified GUI.
"""
import os
import sys
from typing import Optional, List, Dict, Tuple, Any

try:
    from PyQt5.QtCore import Qt, QTimer, QSize, QSettings
    from PyQt5.QtWidgets import (
        QDialog, QVBoxLayout, QHBoxLayout, QListWidget, QStackedWidget,
        QWidget, QPushButton, QLabel, QFrame, QTextEdit, QProgressBar, QMessageBox
    )
    from PyQt5.QtGui import QIcon, QPixmap, QColor, QBrush, QRadialGradient, QLinearGradient
except ImportError:
    try:
        from PyQt6.QtCore import Qt, QTimer, QSize, QSettings
        from PyQt6.QtWidgets import (
            QDialog, QVBoxLayout, QHBoxLayout, QListWidget, QStackedWidget,
            QWidget, QPushButton, QLabel, QFrame, QTextEdit, QProgressBar, QMessageBox
        )
        from PyQt6.QtGui import QIcon, QPixmap, QColor, QBrush, QRadialGradient, QLinearGradient
    except ImportError:
        raise RuntimeError("PyQt5 or PyQt6 must be installed.")

from .utils import t, get_lang
from .tabs.target_discovery import TargetDiscoveryTab
from .tabs.hit_identification import HitIdentificationTab
from .tabs.lead_optimization import LeadOptimizationTab
from .tabs.batch_analysis import BatchAnalysisTab

# Import worker classes if needed for type hinting or global usage
from .workers import AnalysisWorker, GMotifWorker

def _get_logo_path() -> Optional[str]:
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    logo_path = os.path.join(here, "assets", "logo.png")
    if os.path.exists(logo_path):
        return logo_path
    return None

def _get_welcome_bg_path() -> Optional[str]:
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bg_path = os.path.join(here, "assets", "welcome_bg.png")
    if os.path.exists(bg_path):
        return bg_path
    return None

class HeroHeader(QWidget):
    """Custom widget for the hero section with adaptive background and overlay."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedHeight(300)
        self.bg_path = _get_welcome_bg_path()
        
    def paintEvent(self, event):
        from PyQt5.QtGui import QPainter, QBrush, QColor, QRadialGradient, QLinearGradient
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        
        # 1. Draw Background Image
        if self.bg_path:
            pix = QPixmap(self.bg_path)
            if not pix.isNull():
                # Scale pixmap to cover the entire widget
                scaled_pix = pix.scaled(self.size(), Qt.AspectRatioMode.KeepAspectRatioByExpanding, Qt.TransformationMode.SmoothTransformation)
                # Center the pixmap
                x = (self.width() - scaled_pix.width()) // 2
                y = (self.height() - scaled_pix.height()) // 2
                painter.drawPixmap(x, y, scaled_pix)
        else:
            # Fallback gradient
            grad = QLinearGradient(0, 0, self.width(), self.height())
            grad.setColorAt(0, QColor("#1e293b"))
            grad.setColorAt(1, QColor("#0f172a"))
            painter.fillRect(self.rect(), grad)
            
        # 2. Draw Dark Overlay (Mask)
        painter.fillRect(self.rect(), QColor(15, 23, 42, 100)) # rgba(15, 23, 42, 0.4)
        
        # 3. Bottom border
        painter.setPen(QColor("#3b82f6"))
        painter.drawLine(0, self.height()-1, self.width(), self.height()-1)
        painter.end()

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
        self.settings = QSettings("GlueTK", "GlueTK_App")
        
        # State
        self._dark_mode = False
        self._interactions = []
        self._gmotif_hits = []
        self._esp_maps = {}
        
        # Workers
        self.analysis_thread = None
        self.gmotif_thread = None
        
        # Initialize UI
        self.build_ui()
        self.setup_style()
        self._update_content_bg()
        
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
        
        # 4. Batch Analysis
        self.batch_tab = BatchAnalysisTab(self)
        self.content_stack.addWidget(self.batch_tab)
        
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
        self.nav_widget = nav_widget  # Save reference for theme switching
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
        
        # Header (Theme) - REMOVED for specific white-only request
        header = QWidget()
        h_layout = QHBoxLayout(header)
        # self.theme_btn = QPushButton("☾" if self._dark_mode else "☀")
        # self.theme_btn.setFixedSize(30, 30)
        # self.theme_btn.setFlat(True)
        # self.theme_btn.clicked.connect(self.toggle_theme)
        h_layout.addStretch()
        # h_layout.addWidget(self.theme_btn)
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
            "Batch Analysis"
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
        
        btn_resources = QPushButton("Resources")
        btn_resources.setFlat(True)
        btn_resources.clicked.connect(self._show_resources)
        b_layout.addWidget(btn_resources)
        
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
            getattr(self, "obj_combo_c2h2", None),  # C2H2 Zinc Finger
            getattr(self, "obj_combo_apbs", None),
            getattr(self, "pocket_obj_combo", None),
            getattr(self, "ppi_obj_combo", None),
            getattr(self, "pl_obj_combo", None),
            getattr(self, "ll_obj_combo", None),
            getattr(self, "mut_obj_combo", None),
            getattr(self, "obj_combo_surf", None),  # Surface Analysis
            getattr(self, "obj_combo_sim1", None),  # Surface Similarity Object 1
            getattr(self, "obj_combo_sim2", None),  # Surface Similarity Object 2
            getattr(self, "ec_obj_combo", None),  # EC Analysis
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
        # self.theme_btn.setText("☾" if self._dark_mode else "☀")
        self._update_content_bg()
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
            # Update navigation bar style
            if hasattr(self, 'nav_widget'):
                self.nav_widget.setStyleSheet("background-color: #0d1117; border-right: 1px solid #30363d;")
        else:
            # Light theme style
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
            # Update navigation bar style
            if hasattr(self, 'nav_widget'):
                self.nav_widget.setStyleSheet("background-color: #f1f5f9; border-right: 1px solid #e2e8f0;")

    def _update_content_bg(self):
        """Update scroll content background for all tabs when theme changes."""
        bg = "#161b22" if self._dark_mode else "#ffffff"
        for name in ("_target_scroll_content", "_lead_scroll_content", "_hit_scroll_content"):
            w = getattr(self, name, None)
            if w is not None:
                w.setStyleSheet(f"#scroll_content {{ background-color: {bg}; }}")

    def check_environment(self):
        try:
            from ..env_checker import ensure_dependencies
            ensure_dependencies(silent=False)
        except ImportError:
            QMessageBox.warning(self, "Error", "Environment checker not found.")

    # --- Simple Pages ---
    def create_welcome_tab(self) -> QWidget:
        """Create a modern, visually stunning welcome page following the scientific showcase concept."""
        from PyQt5.QtWidgets import QGraphicsDropShadowEffect, QGridLayout, QScrollArea
        
        main_page = QWidget()
        main_page.setObjectName("welcome_page")
        main_layout = QVBoxLayout(main_page)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Scroll Area for the whole page
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        content_container = QWidget()
        layout = QVBoxLayout(content_container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        scroll.setWidget(content_container)
        main_layout.addWidget(scroll)

        # Hero Section (Top)
        hero_widget = HeroHeader()
        hero_layout = QVBoxLayout(hero_widget)
        hero_layout.setContentsMargins(40, 40, 40, 40)
        
        title_lbl = QLabel("GlueTK")
        title_lbl.setStyleSheet("font-size: 80px; font-weight: 900; color: #ffffff; background: transparent; letter-spacing: 4px;")
        
        # Add a subtle glow/shadow to Title
        title_shadow = QGraphicsDropShadowEffect()
        title_shadow.setBlurRadius(20)
        title_shadow.setColor(Qt.GlobalColor.black)
        title_shadow.setOffset(0, 0)
        title_lbl.setGraphicsEffect(title_shadow)
        
        slogan_lbl = QLabel("Professional Toolkit for Molecular Glue Engineering")
        slogan_lbl.setStyleSheet("font-size: 24px; color: #f8fafc; background: transparent; font-weight: 500; letter-spacing: 1px;")
        
        # Add shadow to Slogan
        slogan_shadow = QGraphicsDropShadowEffect()
        slogan_shadow.setBlurRadius(10)
        slogan_shadow.setColor(Qt.GlobalColor.black)
        slogan_shadow.setOffset(2, 2)
        slogan_lbl.setGraphicsEffect(slogan_shadow)
        
        hero_layout.addWidget(title_lbl)
        hero_layout.addWidget(slogan_lbl)
        hero_layout.addStretch()
        layout.addWidget(hero_widget)

        # Content Section
        content_area = QWidget()
        content_layout = QVBoxLayout(content_area)
        content_layout.setContentsMargins(40, 40, 40, 40)
        content_layout.setSpacing(40)

        # SECTION 1: CORE TOOLS
        tools_lbl = QLabel("CORE MODULES")
        tools_lbl.setStyleSheet("font-size: 14px; font-weight: bold; color: #3b82f6; letter-spacing: 1px;")
        content_layout.addWidget(tools_lbl)

        grid_layout = QGridLayout()
        grid_layout.setSpacing(20)

        def create_card(title: str, desc: str, icon_text: str, index: int):
            card = QFrame()
            card.setObjectName("feature_card")
            card.setCursor(Qt.CursorShape.PointingHandCursor)
            card.setStyleSheet("""
                #feature_card {
                    background-color: rgba(255, 255, 255, 0.05);
                    border: 1px solid rgba(255, 255, 255, 0.1);
                    border-radius: 12px; padding: 20px;
                }
                #feature_card:hover { background-color: rgba(59, 130, 246, 0.1); border: 1px solid #3b82f6; }
            """ if self._dark_mode else """
                #feature_card {
                    background-color: #ffffff; border: 1px solid #e2e8f0;
                    border-radius: 12px; padding: 20px;
                }
                #feature_card:hover { border-color: #3b82f6; background-color: #f8fafc; }
            """)
            
            shadow = QGraphicsDropShadowEffect()
            shadow.setBlurRadius(15)
            shadow.setXOffset(0)
            shadow.setYOffset(4)
            shadow.setColor(Qt.GlobalColor.black if self._dark_mode else Qt.GlobalColor.lightGray)
            card.setGraphicsEffect(shadow)

            c_layout = QVBoxLayout(card)
            ico = QLabel(icon_text); ico.setStyleSheet("font-size: 32px; background: transparent;")
            t_lbl = QLabel(title); t_lbl.setStyleSheet("font-size: 18px; font-weight: bold; color: #3b82f6; background: transparent;")
            d_lbl = QLabel(desc); d_lbl.setWordWrap(True); d_lbl.setStyleSheet("color: #64748b; background: transparent;")
            
            c_layout.addWidget(ico); c_layout.addWidget(t_lbl); c_layout.addWidget(d_lbl); c_layout.addStretch()
            card.mousePressEvent = lambda e: self.nav_list.setCurrentRow(index)
            return card

        grid_layout.addWidget(create_card("Target Discovery", "Identify promising residues for ternary complex formation.", "🎯", 1), 0, 0)
        grid_layout.addWidget(create_card("Hit Identification", "Analyze interactions and interface complementarity.", "🔍", 2), 0, 1)
        grid_layout.addWidget(create_card("Lead Optimization", "Fine-tune molecules for better binding affinity.", "⚡", 3), 1, 0)
        grid_layout.addWidget(create_card("Visualization", "High-quality rendering of interaction networks.", "🎨", 4), 1, 1)
        content_layout.addLayout(grid_layout)

        # SECTION 2: CLASSIC CASES (User Request)
        case_lbl = QLabel("FEATURED CASES")
        case_lbl.setStyleSheet("font-size: 14px; font-weight: bold; color: #3b82f6; margin-top: 20px; letter-spacing: 1px;")
        content_layout.addWidget(case_lbl)

        case_row = QHBoxLayout()
        case_row.setSpacing(20)

        def create_case_box(name: str, pdb_id: str, desc: str):
            box = QFrame()
            box.setCursor(Qt.CursorShape.PointingHandCursor)
            box.setStyleSheet("""
                QFrame { background: #f1f5f9; border-radius: 8px; padding: 15px; border: 1px solid #e2e8f0; }
                QFrame:hover { background: #e2e8f0; border-color: #3b82f6; }
            """ if not self._dark_mode else """
                QFrame { background: #1e293b; border-radius: 8px; padding: 15px; border: 1px solid #334155; }
                QFrame:hover { background: #334155; border-color: #3b82f6; }
            """)
            b_layout = QVBoxLayout(box)
            h_lbl = QLabel(f"{name} ({pdb_id})")
            h_lbl.setStyleSheet("font-weight: bold; color: #0f172a;" if not self._dark_mode else "font-weight: bold; color: #f1f5f9;")
            d_lbl = QLabel(desc); d_lbl.setStyleSheet("font-size: 12px; color: #64748b;"); d_lbl.setWordWrap(True)
            b_layout.addWidget(h_lbl); b_layout.addWidget(d_lbl)
            
            box.mousePressEvent = lambda e: self._load_case(pdb_id)
            return box

        case_row.addWidget(create_case_box("Thalidomide", "4CIW", "CRBN degrader, the classic prototype."))
        case_row.addWidget(create_case_box("Indisulam", "5S9M", "DCAF15-RBM39 splicing factor degrader."))
        case_row.addWidget(create_case_box("QS1", "8P1D", "GSPT1 degrader, next-gen target."))
        content_layout.addLayout(case_row)

        layout.addWidget(content_area)
        layout.addStretch()
        return main_page

    def _load_case(self, pdb_id: str):
        """Helper to load a PDB structure in PyMOL."""
        try:
            from pymol import cmd
            self.log(f"Fetching PDB: {pdb_id}")
            cmd.reinitialize()
            cmd.fetch(pdb_id)
            cmd.show_as("cartoon")
            cmd.show("sticks", f"resn {pdb_id}") # Simple attempt to show ligand if named like PDB
            cmd.zoom()
            QMessageBox.information(self, "Case Loaded", f"Successfully loaded {pdb_id} from PDB.")
            self.refresh_objects()
        except Exception as e:
            self.on_error(f"Failed to load PDB: {str(e)}")

    def create_readme_tab(self) -> QWidget:
        """Create a comprehensive README tab with GUI documentation."""
        from PyQt5.QtWidgets import QScrollArea, QTextBrowser
        
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(0, 0, 0, 0)
        
        browser = QTextBrowser()
        browser.setReadOnly(True)
        browser.setOpenExternalLinks(True)
        
        html_content = """
        <html>
        <head>
        <style>
            body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; padding: 30px; line-height: 1.6; color: #1f2937; }
            h1 { color: #1e40af; font-size: 28px; border-bottom: 3px solid #3b82f6; padding-bottom: 10px; }
            h2 { color: #3b82f6; font-size: 22px; margin-top: 30px; border-bottom: 1px solid #e5e7eb; padding-bottom: 8px; }
            h3 { color: #6366f1; font-size: 18px; margin-top: 20px; }
            h4 { color: #8b5cf6; font-size: 16px; margin-top: 15px; }
            p { margin: 10px 0; }
            ul { margin: 10px 0 10px 20px; }
            li { margin: 5px 0; }
            code { background: #f3f4f6; padding: 2px 6px; border-radius: 4px; font-family: 'Monaco', 'Consolas', monospace; }
            .feature-box { background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 8px; padding: 15px; margin: 15px 0; }
            .feature-title { color: #3b82f6; font-weight: bold; font-size: 16px; margin-bottom: 10px; }
            table { border-collapse: collapse; width: 100%; margin: 15px 0; }
            th { background: #3b82f6; color: white; padding: 10px; text-align: left; }
            td { border: 1px solid #e5e7eb; padding: 8px; }
            tr:nth-child(even) { background: #f9fafb; }
            .emoji { font-size: 20px; margin-right: 8px; }
            .tip { background: #ecfdf5; border-left: 4px solid #10b981; padding: 10px 15px; margin: 15px 0; }
            .warning { background: #fef3c7; border-left: 4px solid #f59e0b; padding: 10px 15px; margin: 15px 0; }
        </style>
        </head>
        <body>
        
        <h1>🧬 GlueTK User Guide</h1>
        <p>Welcome to GlueTK - Professional Toolkit for Molecular Glue Engineering</p>
        
        <h2><span class="emoji">🎯</span>Target Discovery Tab</h2>
        <p>Identify potential molecular glue targets and binding sites.</p>
        
        <div class="feature-box">
            <div class="feature-title">G-Motif (CRBN G-loop) Detection</div>
            <p>Detect CRBN G-loop binding motifs using RMSD-based template matching.</p>
            <table>
                <tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
                <tr><td>Target Object</td><td>PyMOL object to analyze</td><td>-</td></tr>
                <tr><td>RMSD Cutoff</td><td>Maximum RMSD for match (Å)</td><td>3.5</td></tr>
                <tr><td>Template</td><td>Reference structure (GSPT1/CK1α/VAV1)</td><td>GSPT1</td></tr>
                <tr><td>Require Central Gly</td><td>Filter for glycine at position 4</td><td>✓</td></tr>
                <tr><td>Highlight Surface</td><td>Show molecular surface of G-loop</td><td>✓</td></tr>
            </table>
            <p><b>Buttons:</b></p>
            <ul>
                <li><b>Detect POI</b> - Run G-motif detection</li>
                <li><b>Show Surface</b> - Highlight G-loop surface</li>
                <li><b>Export Coords</b> - Save coordinates to CSV</li>
                <li><b>Surface Analysis</b> - Analyze electrostatic/hydrophobic patches</li>
            </ul>
        </div>
        
        <div class="feature-box">
            <div class="feature-title">C2H2 Zinc Finger Detection</div>
            <p>Find C2H2 zinc finger domains (e.g., IKZF1/3 for lenalidomide targets).</p>
            <table>
                <tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
                <tr><td>Turn RMSD</td><td>Local turn alignment threshold (Å)</td><td>2.0</td></tr>
                <tr><td>Global RMSD</td><td>Global fold check threshold (Å)</td><td>3.5</td></tr>
                <tr><td>Require Turn Gly</td><td>Require key Gly in turn region</td><td>✗</td></tr>
            </table>
            <p><b>Buttons:</b> Find Zinc Fingers, Render C2H2 + ESP, Highlight Domains</p>
        </div>
        
        <div class="feature-box">
            <div class="feature-title">Surface Similarity & Complementarity (MaSIF-style)</div>
            <p>Compare binding site features using geometric and chemical descriptors.</p>
            <table>
                <tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
                <tr><td>Analysis Type</td><td>Similarity or Complementarity (PPI)</td><td>Similarity</td></tr>
                <tr><td>Surface Method</td><td>auto/msms/open3d/edtsurf</td><td>auto</td></tr>
                <tr><td>Patch Radius</td><td>Radius for comparison (Å)</td><td>12.0</td></tr>
                <tr><td>Interface Distance</td><td>Contact threshold (Å)</td><td>4.0</td></tr>
            </table>
        </div>
        
        <h2><span class="emoji">🔍</span>Hit Identification Tab</h2>
        <p>Identify binding sites and perform molecular docking.</p>
        
        <div class="feature-box">
            <div class="feature-title">Binding Site Detection</div>
            <p>Detect druggable binding pockets with volume, depth, and druggability scoring.</p>
            <table>
                <tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
                <tr><td>Grid Spacing</td><td>Detection resolution (Å)</td><td>0.6</td></tr>
                <tr><td>Min Volume</td><td>Minimum pocket volume (Å³)</td><td>20</td></tr>
                <tr><td>Color By</td><td>Volume/Druggability/Hydrophobicity/Depth</td><td>Volume</td></tr>
            </table>
        </div>
        
        <div class="feature-box">
            <div class="feature-title">AutoDock Vina Docking</div>
            <p>Integrated molecular docking to detected pockets.</p>
            <table>
                <tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
                <tr><td>Ligand File</td><td>MOL2/SDF/PDBQT file</td><td>-</td></tr>
                <tr><td>Max Pockets</td><td>Number of pockets to dock</td><td>3</td></tr>
                <tr><td>Exhaustiveness</td><td>Search thoroughness</td><td>8</td></tr>
            </table>
        </div>
        
        <div class="feature-box">
            <div class="feature-title">HADDOCK3 Protein-Protein Docking</div>
            <p>Model ternary complexes with protein-protein docking.</p>
            <table>
                <tr><th>Mode</th><th>Description</th></tr>
                <tr><td>Blind (Random AIR)</td><td>No restraints, random sampling</td></tr>
                <tr><td>Blind (Centroid)</td><td>Center-of-mass guided</td></tr>
                <tr><td>Blind (Surface)</td><td>Surface-based sampling</td></tr>
                <tr><td>Pocket-constrained</td><td>Residue-based AIR restraints</td></tr>
            </table>
        </div>
        
        <h2><span class="emoji">⚡</span>Lead Optimization Tab</h2>
        <p>Analyze and optimize molecular interactions.</p>
        
        <div class="feature-box">
            <div class="feature-title">Protein-Protein Interface (PPI) Analysis</div>
            <p>Analyze protein-protein interfaces with interaction classification.</p>
            <table>
                <tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
                <tr><td>Protein1 Chains</td><td>First protein chain IDs</td><td>-</td></tr>
                <tr><td>Protein2 Chains</td><td>Second protein chain IDs</td><td>-</td></tr>
                <tr><td>Interface Distance</td><td>Contact threshold (Å)</td><td>4.5</td></tr>
                <tr><td>Display Mode</td><td>Cartoon+Surface or Surface only</td><td>Cartoon+Surface</td></tr>
            </table>
            <p><b>Output:</b> Interface contacts, Buried Surface Area (BSA), Interface strength score (0-10)</p>
        </div>
        
        <div class="feature-box">
            <div class="feature-title">Protein-Ligand Interactions</div>
            <p>Comprehensive interaction analysis with visualization.</p>
            <table>
                <tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
                <tr><td>Ligand Name</td><td>Residue name (auto-detect if blank)</td><td>-</td></tr>
                <tr><td>Protein Chains</td><td>Specific chains (optional)</td><td>-</td></tr>
                <tr><td>Distance</td><td>Interaction cutoff (Å)</td><td>4.5</td></tr>
            </table>
            <p><b>Interaction Types Detected:</b></p>
            <ul>
                <li>Hydrogen bonds</li>
                <li>Salt bridges</li>
                <li>Hydrophobic contacts</li>
                <li>π-π stacking</li>
                <li>Cation-π interactions</li>
                <li>Halogen bonds</li>
            </ul>
            <p><b>Buttons:</b> Analyze Protein-Ligand (3D), Generate 2D Diagram (RDKit)</p>
        </div>
        
        <div class="feature-box">
            <div class="feature-title">Ligand-Ligand Interactions</div>
            <p>Analyze interactions between two small molecules.</p>
            <p>Use PyMOL selection syntax (e.g., <code>resn LIG1</code>, <code>resi 100</code>)</p>
        </div>
        
        <div class="feature-box">
            <div class="feature-title">Mutation Analysis</div>
            <p>Predict mutation effects on binding.</p>
            <ul>
                <li><b>Perform Mutation</b> - Apply mutations via PyMOL wizard</li>
                <li><b>Minimize</b> - Energy minimization with sculpting</li>
                <li><b>Full Analysis</b> - ΔΔG prediction (requires FoldX/PyRosetta)</li>
            </ul>
        </div>
        
        <h2><span class="emoji">📊</span>Batch Analysis Tab</h2>
        <p>Analyze multiple structures simultaneously.</p>
        
        <div class="feature-box">
            <div class="feature-title">Input Methods</div>
            <table>
                <tr><th>Method</th><th>Description</th></tr>
                <tr><td>PDB IDs</td><td>Enter comma-separated IDs (e.g., 6H0G,6H0F,5FQD)</td></tr>
                <tr><td>Browse Files</td><td>Select individual PDB/CIF files</td></tr>
                <tr><td>Browse Folder</td><td>Add all PDB files from a directory</td></tr>
            </table>
        </div>
        
        <div class="feature-box">
            <div class="feature-title">Analysis Types</div>
            <table>
                <tr><th>Type</th><th>Description</th></tr>
                <tr><td>G-motif Detection</td><td>Batch G-loop finding</td></tr>
                <tr><td>PPI Interface</td><td>Compare interfaces across structures</td></tr>
                <tr><td>Pocket Detection</td><td>Screen for druggable sites</td></tr>
                <tr><td>Protein-Ligand</td><td>Batch interaction analysis</td></tr>
                <tr><td>Comprehensive</td><td>Run all analyses</td></tr>
            </table>
            <p><b>Output:</b> CSV files with detailed results, summary statistics, JSON reports</p>
        </div>
        
        <h2><span class="emoji">🎨</span>Visualization (APBS)</h2>
        <p>Generate publication-quality images with electrostatic surfaces.</p>
        
        <div class="feature-box">
            <div class="feature-title">Electrostatic Surface</div>
            <table>
                <tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
                <tr><td>Grid Spacing</td><td>Resolution (Å)</td><td>1.0</td></tr>
                <tr><td>Range</td><td>Color scale (min,mid,max)</td><td>-5,0,5</td></tr>
            </table>
            <p><b>Buttons:</b></p>
            <ul>
                <li><b>Quick ESP</b> - Fast Coulomb-based calculation</li>
                <li><b>True APBS</b> - Full APBS calculation (requires APBS tools)</li>
            </ul>
        </div>
        
        <div class="feature-box">
            <div class="feature-title">Image Export</div>
            <table>
                <tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
                <tr><td>Width</td><td>Image width (px)</td><td>3000</td></tr>
                <tr><td>Height</td><td>Image height (px)</td><td>2000</td></tr>
                <tr><td>DPI</td><td>Resolution</td><td>300</td></tr>
                <tr><td>Background</td><td>White or Transparent</td><td>White</td></tr>
                <tr><td>Ray Trace</td><td>High-quality rendering</td><td>✓</td></tr>
            </table>
            <p><b>Buttons:</b> Export PNG, Export DX, Fill Viewport</p>
        </div>
        
        <h2><span class="emoji">⌨️</span>Quick Commands</h2>
        <div class="tip">
            <p>You can also use these commands in PyMOL command line:</p>
        </div>
        <table>
            <tr><th>Command</th><th>Description</th></tr>
            <tr><td><code>gluetk_gui</code></td><td>Launch GlueTK GUI</td></tr>
            <tr><td><code>find_crbn_g_motif obj, rmsd_cutoff=3.5</code></td><td>G-motif detection</td></tr>
            <tr><td><code>ppi_analyze obj, ['A'], ['B']</code></td><td>PPI interface analysis</td></tr>
            <tr><td><code>neo_epitope_find obj, ['A'], ['B'], LIG</code></td><td>Neo-epitope detection</td></tr>
            <tr><td><code>generate_2d_diagram obj, LIG</code></td><td>2D interaction diagram</td></tr>
            <tr><td><code>detect_pockets obj</code></td><td>Pocket detection</td></tr>
            <tr><td><code>batch_gmotif "6H0G,6H0F"</code></td><td>Batch G-motif</td></tr>
        </table>
        
        <h2><span class="emoji">📚</span>Resources</h2>
        <p>Click <b>Resources</b> button in the sidebar for databases and literature.</p>
        <p>GitHub: <a href="https://github.com/VesperChen01/GlueTK">https://github.com/VesperChen01/GlueTK</a></p>
        
        <div class="warning">
            <b>Note:</b> Some features require additional dependencies (RDKit for 2D diagrams, AutoDock Vina for docking, APBS for accurate electrostatics).
        </div>
        
        </body>
        </html>
        """
        
        browser.setHtml(html_content)
        layout.addWidget(browser)
        return w

    def create_contact_tab(self) -> QWidget:
        """Create a clean contact page with structured information."""
        from PyQt5.QtWidgets import QScrollArea, QGridLayout
        
        main_page = QWidget()
        main_layout = QVBoxLayout(main_page)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        
        # Scroll Area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        
        content = QWidget()
        content.setStyleSheet("background-color: #ffffff;")
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(60, 50, 60, 50)
        content_layout.setSpacing(30)
        
        # Header Section
        header_section = QWidget()
        header_layout = QVBoxLayout(header_section)
        header_layout.setSpacing(10)
        
        title_lbl = QLabel("Contact Us")
        title_lbl.setStyleSheet("font-size: 36px; font-weight: 800; color: #1e293b;")
        
        subtitle_lbl = QLabel("We'd love to hear from you.")
        subtitle_lbl.setStyleSheet("font-size: 20px; color: #3b82f6; font-weight: 500;")
        
        desc_lbl = QLabel("If you have any questions, suggestions, or collaboration ideas while using GlueTK, feel free to reach out anytime.")
        desc_lbl.setStyleSheet("font-size: 14px; color: #64748b;")
        desc_lbl.setWordWrap(True)
        
        header_layout.addWidget(title_lbl)
        header_layout.addWidget(subtitle_lbl)
        header_layout.addWidget(desc_lbl)
        content_layout.addWidget(header_section)
        
        # Divider
        divider1 = QFrame()
        divider1.setFrameShape(QFrame.Shape.HLine)
        divider1.setStyleSheet("background-color: #e2e8f0; max-height: 1px;")
        content_layout.addWidget(divider1)
        
        # About the Maintainer Section
        maintainer_section = QWidget()
        maintainer_layout = QVBoxLayout(maintainer_section)
        maintainer_layout.setSpacing(20)
        
        maintainer_title = QLabel("👤 About the Maintainer")
        maintainer_title.setStyleSheet("font-size: 18px; font-weight: bold; color: #1e293b;")
        maintainer_layout.addWidget(maintainer_title)
        
        # Info Grid
        info_grid = QGridLayout()
        info_grid.setSpacing(20)
        info_grid.setColumnStretch(1, 1)
        
        # Name & Affiliation
        name_card = self._create_info_card("Roufen Chen", "Zhejiang University", "#3b82f6")
        info_grid.addWidget(name_card, 0, 0)
        
        # Email
        email_card = self._create_info_card("📧 Email", "12319021@zju.edu.cn", "#10b981")
        info_grid.addWidget(email_card, 0, 1)
        
        # GitHub
        github_card = self._create_info_card("💻 GitHub", "VesperChen01 / GlueTK", "#8b5cf6")
        info_grid.addWidget(github_card, 1, 0, 1, 2)
        
        maintainer_layout.addLayout(info_grid)
        content_layout.addWidget(maintainer_section)
        
        # Divider
        divider2 = QFrame()
        divider2.setFrameShape(QFrame.Shape.HLine)
        divider2.setStyleSheet("background-color: #e2e8f0; max-height: 1px;")
        content_layout.addWidget(divider2)
        
        # Send Us a Message Section
        message_section = QWidget()
        message_layout = QVBoxLayout(message_section)
        message_layout.setSpacing(15)
        
        message_title = QLabel("✉️ Send Us a Message")
        message_title.setStyleSheet("font-size: 18px; font-weight: bold; color: #1e293b;")
        message_layout.addWidget(message_title)
        
        message_desc = QLabel("Please contact us through one of the following ways:")
        message_desc.setStyleSheet("font-size: 14px; color: #64748b;")
        message_layout.addWidget(message_desc)
        
        # Contact options
        options = [
            ("🐛 Bug Report", "Found a bug or unexpected behavior? Please report it with reproduction steps and environment info."),
            ("💡 Feature Request", "Have a new feature idea or improvement suggestion? We'd love to hear from you."),
            ("🤝 Collaboration", "Interested in GlueTK and want to collaborate or discuss? Feel free to reach out."),
            ("📝 Other Inquiries", "Any other questions or ideas."),
        ]
        
        for title, desc in options:
            option_widget = self._create_contact_option(title, desc)
            message_layout.addWidget(option_widget)
        
        content_layout.addWidget(message_section)
        content_layout.addStretch()
        
        scroll.setWidget(content)
        main_layout.addWidget(scroll)
        
        return main_page
    
    def _create_info_card(self, title: str, value: str, accent_color: str) -> QWidget:
        """Create an info card for contact page."""
        card = QFrame()
        card.setStyleSheet(f"""
            QFrame {{
                background-color: #f8fafc;
                border: 1px solid #e2e8f0;
                border-radius: 12px;
                padding: 20px;
            }}
        """)
        
        layout = QVBoxLayout(card)
        layout.setSpacing(5)
        
        title_lbl = QLabel(title)
        title_lbl.setStyleSheet(f"font-size: 16px; font-weight: bold; color: {accent_color};")
        
        value_lbl = QLabel(value)
        value_lbl.setStyleSheet("font-size: 14px; color: #64748b;")
        
        layout.addWidget(title_lbl)
        layout.addWidget(value_lbl)
        
        return card
    
    def _create_contact_option(self, title: str, desc: str) -> QWidget:
        """Create a contact option item."""
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(15, 10, 15, 10)
        layout.setSpacing(5)
        
        title_lbl = QLabel(title)
        title_lbl.setStyleSheet("font-size: 15px; font-weight: bold; color: #1e293b;")
        
        desc_lbl = QLabel(desc)
        desc_lbl.setStyleSheet("font-size: 13px; color: #64748b;")
        desc_lbl.setWordWrap(True)
        
        layout.addWidget(title_lbl)
        layout.addWidget(desc_lbl)
        
        widget.setStyleSheet("""
            QWidget {
                background-color: #f8fafc;
                border-radius: 8px;
            }
            QWidget:hover {
                background-color: #f1f5f9;
            }
        """)
        
        return widget
    
    def _show_resources(self):
        """Show molecular glue resources dialog"""
        try:
            # Try to find and open the resources markdown file
            here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            resources_path = os.path.join(here, "resources", "molecular_glue_resources.md")
            
            if os.path.exists(resources_path):
                # Create a dialog to show resources
                from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTextBrowser, QPushButton
                
                dialog = QDialog(self)
                dialog.setWindowTitle("Molecular Glue Resources")
                dialog.setMinimumSize(800, 600)
                
                layout = QVBoxLayout(dialog)
                
                browser = QTextBrowser()
                browser.setOpenExternalLinks(True)
                
                # Read and convert markdown to HTML (simple conversion)
                with open(resources_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                # Simple markdown to HTML conversion
                html_content = self._markdown_to_html(content)
                browser.setHtml(html_content)
                
                layout.addWidget(browser)
                
                close_btn = QPushButton("Close")
                close_btn.clicked.connect(dialog.close)
                layout.addWidget(close_btn)
                
                dialog.exec_()
            else:
                QMessageBox.information(self, "Resources",
                    "Resources file not found.\n\n"
                    "Key databases:\n"
                    "• PROTAC-DB: https://protacdb.weizmann.ac.il/\n"
                    "• Open Targets: https://www.opentargets.org/\n"
                    "• PDB: https://www.rcsb.org/\n\n"
                    "Key PDB structures:\n"
                    "• 6H0G - CRBN-CC885-GSPT1\n"
                    "• 5FQD - CRBN-Lenalidomide-CK1α\n"
                    "• 5S9M - DCAF15-Indisulam-RBM39")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to load resources: {e}")
    
    def _markdown_to_html(self, md_content: str) -> str:
        """Simple markdown to HTML conversion"""
        import re
        
        html = md_content
        
        # Headers
        html = re.sub(r'^### (.+)$', r'<h3>\1</h3>', html, flags=re.MULTILINE)
        html = re.sub(r'^## (.+)$', r'<h2>\1</h2>', html, flags=re.MULTILINE)
        html = re.sub(r'^# (.+)$', r'<h1>\1</h1>', html, flags=re.MULTILINE)
        
        # Bold
        html = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', html)
        
        # Links
        html = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<a href="\2">\1</a>', html)
        
        # Code blocks
        html = re.sub(r'```(\w+)?\n(.*?)```', r'<pre><code>\2</code></pre>', html, flags=re.DOTALL)
        html = re.sub(r'`([^`]+)`', r'<code>\1</code>', html)
        
        # Tables (simple)
        lines = html.split('\n')
        in_table = False
        new_lines = []
        for line in lines:
            if '|' in line and not line.strip().startswith('```'):
                if not in_table:
                    new_lines.append('<table border="1" cellpadding="5" cellspacing="0">')
                    in_table = True
                if line.strip().startswith('|---') or line.strip().startswith('| ---'):
                    continue  # Skip separator line
                cells = [c.strip() for c in line.split('|')[1:-1]]
                row = '<tr>' + ''.join(f'<td>{c}</td>' for c in cells) + '</tr>'
                new_lines.append(row)
            else:
                if in_table:
                    new_lines.append('</table>')
                    in_table = False
                new_lines.append(line)
        if in_table:
            new_lines.append('</table>')
        html = '\n'.join(new_lines)
        
        # Horizontal rules
        html = re.sub(r'^---+$', r'<hr>', html, flags=re.MULTILINE)
        
        # Line breaks
        html = html.replace('\n\n', '</p><p>')
        html = f'<p>{html}</p>'
        
        # Style
        html = f'''
        <html>
        <head>
        <style>
            body {{ font-family: Arial, sans-serif; padding: 20px; }}
            h1 {{ color: #1e40af; }}
            h2 {{ color: #3b82f6; border-bottom: 1px solid #e5e7eb; padding-bottom: 5px; }}
            h3 {{ color: #6366f1; }}
            a {{ color: #2563eb; }}
            code {{ background: #f3f4f6; padding: 2px 6px; border-radius: 4px; }}
            pre {{ background: #1e293b; color: #e2e8f0; padding: 15px; border-radius: 8px; overflow-x: auto; }}
            table {{ border-collapse: collapse; margin: 10px 0; }}
            td {{ border: 1px solid #e5e7eb; padding: 8px; }}
            tr:nth-child(even) {{ background: #f9fafb; }}
        </style>
        </head>
        <body>
        {html}
        </body>
        </html>
        '''
        
        return html
