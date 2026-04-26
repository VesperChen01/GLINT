# -*- coding: utf-8 -*-
"""
GLINT Main Window
Modularized version of the Unified GUI.
"""
import os
import sys
import platform
from typing import Optional, List, Dict, Tuple, Any

from .qt_adapter import QtCore, QtWidgets, QtGui, Qt, Signal, Slot, Property

if QtWidgets is None:
    raise RuntimeError("No Qt binding (PyQt5, PyQt6, PySide2, or PySide6) found.")

# Re-map commonly used classes from submodules for compatibility
from .qt_adapter import (
    QDialog, QVBoxLayout, QHBoxLayout, QListWidget, QStackedWidget,
    QWidget, QPushButton, QLabel, QFrame, QTextEdit, QProgressBar,
    QIcon, QPixmap, QColor, QLinearGradient,
    QTimer, QSize, QSettings, QScrollArea
)

from .utils import t, get_lang, show_message_box, show_question_box, _get_logo_path, _get_asset_path
from .tabs.target_discovery import TargetDiscoveryTab
from .tabs.hit_identification import HitIdentificationTab
from .tabs.lead_optimization import LeadOptimizationTab
from .tabs.ternary_evaluation import TernaryEvaluationTab

# Import worker classes if needed for type hinting or global usage
from .workers import AnalysisWorker, GMotifWorker
from .._version import __version__

# Navigation bar page index constants（Eliminate magic numbers）
NAV_INDEX_README = 5
NAV_INDEX_CONTACT = 6

def _get_platform_font() -> str:
    """Get platform-specific font family for best rendering."""
    system = platform.system()
    if system == "Darwin":  # macOS
        return "SF Pro Text, -apple-system, Helvetica Neue, sans-serif"
    elif system == "Windows":
        return "Segoe UI, Microsoft YaHei UI, sans-serif"
    else:  # Linux
        return "Ubuntu, Noto Sans, DejaVu Sans, sans-serif"

def _get_platform_adjustments() -> Dict[str, int]:
    """Get platform-specific size adjustments."""
    system = platform.system()
    if system == "Darwin":  # macOS - Retina displays
        return {"min_height": 24, "padding": 4, "font_size": 13}
    elif system == "Windows":  # Windows - needs slightly larger
        return {"min_height": 26, "padding": 5, "font_size": 12}
    else:  # Linux
        return {"min_height": 26, "padding": 5, "font_size": 12}

class HeroHeader(QWidget):
    """Custom widget for the hero section with adaptive background and overlay."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedHeight(220)
        # using通用的 _get_asset_path 替代重复的 _get_welcome_bg_path
        self.bg_path = _get_asset_path("welcome_bg.png")
        
    def paintEvent(self, event):
        from .qt_adapter import QPainter, QColor, QLinearGradient
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

class GLINTDialog(QDialog):
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
        self.settings = QSettings("GLINT", "GLINT_App")
        
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

        # Helper to wrap tab in scroll area
        def wrap_in_scroll(tab_widget):
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QFrame.Shape.NoFrame)
            scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            scroll.setWidget(tab_widget)
            return scroll

        # 1. Target Discovery
        self.target_tab = TargetDiscoveryTab(self)
        self.content_stack.addWidget(wrap_in_scroll(self.target_tab))

        # 2. Hit Identification
        self.hit_tab = HitIdentificationTab(self)
        self.content_stack.addWidget(wrap_in_scroll(self.hit_tab))

        # 3. Ternary Evaluation (moved before Lead Optimization)
        self.ternary_tab = TernaryEvaluationTab(self)
        self.content_stack.addWidget(wrap_in_scroll(self.ternary_tab))

        # 4. Lead Optimization
        self.lead_tab = LeadOptimizationTab(self)
        self.content_stack.addWidget(wrap_in_scroll(self.lead_tab))

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
        logo_frame.setFixedHeight(130)
        logo_layout = QVBoxLayout(logo_frame)
        logo_layout.setContentsMargins(0, 0, 0, 0)
        logo_path = _get_logo_path()
        if logo_path:
            logo_lbl = QLabel()
            pix = QPixmap(logo_path)
            logo_lbl.setPixmap(pix.scaledToHeight(120, Qt.TransformationMode.SmoothTransformation))
            logo_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            logo_layout.addWidget(logo_lbl)
        layout.addWidget(logo_frame)
        
        # 顶部间距区域（主题切换已Remove，仅保留布局）
        header = QWidget()
        h_layout = QHBoxLayout(header)
        h_layout.addStretch()
        h_layout.setContentsMargins(10, 5, 10, 5)
        layout.addWidget(header)
        
        # Nav List
        self.nav_list = QListWidget()
        self.nav_list.setFrameShape(QFrame.Shape.NoFrame)
        # 减少 padding，使菜单更紧凑
        self.nav_list.setStyleSheet("""
            QListWidget { background: transparent; font-size: 14px; }
            QListWidget::item { padding: 10px 16px; margin: 2px 8px; border-radius: 6px; }
            QListWidget::item:selected { background: #3b82f6; color: white; }
            QListWidget::item:hover:!selected { background: rgba(59, 130, 246, 0.1); }
        """)
        
        items = [
            "Welcome",
            "Target Discovery",
            "Hit Identification",
            "Ternary Evaluation",
            "Lead Optimization",
        ]
        self.nav_list.addItems(items)
        self.nav_list.setCurrentRow(0)
        self.nav_list.currentRowChanged.connect(self.on_nav_changed)
        layout.addWidget(self.nav_list, 1)  # stretch factor 1，让列表填充可用空间

        # Bottom links
        bottom_bar = QWidget()
        b_layout = QVBoxLayout(bottom_bar)
        b_layout.setSpacing(2)

        btn_readme = QPushButton("README")
        btn_readme.setFlat(True)
        btn_readme.clicked.connect(lambda: self.content_stack.setCurrentIndex(NAV_INDEX_README))
        b_layout.addWidget(btn_readme)

        btn_contact = QPushButton("Contact Us")
        btn_contact.setFlat(True)
        btn_contact.clicked.connect(lambda: self.content_stack.setCurrentIndex(NAV_INDEX_CONTACT))
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
        print(f"[GLINT] {msg}")
        # self.log_edit.append(msg) # Log edit is hidden

    def on_error(self, msg: str):
        self.log(f"Error: {msg}")
        show_message_box(self, "Error", msg, icon_type="critical")
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 1)
        # Re-enable buttons if needed
        if hasattr(self, 'gm_btn'): self.gm_btn.setEnabled(True)
        if hasattr(self, 'surf_btn'): self.surf_btn.setEnabled(True)
        if hasattr(self, 'sim_similarity_btn'): self.sim_similarity_btn.setEnabled(True)
        if hasattr(self, 'sim_complement_btn'): self.sim_complement_btn.setEnabled(True)

    def refresh_objects(self):
        names = []
        try:
            from pymol import cmd
            if hasattr(cmd, "get_names"):
                names = cmd.get_names("objects")
            else:
                names = cmd.get_object_list()
        except Exception:  # PyMOL API 调用可能Failed
            names = []
            
        if not names: names = [t("no_object")]
        
        current_by_combo = {}

        def remember_selection(cb):
            if cb is not None:
                current_by_combo[cb] = cb.currentText().strip()

        def restore_selection(cb, fallback_index=0):
            if cb is None:
                return
            previous = current_by_combo.get(cb, "")
            idx = cb.findText(previous) if previous else -1
            if idx < 0 and cb.count():
                idx = min(fallback_index, cb.count() - 1)
            if idx >= 0:
                cb.setCurrentIndex(idx)

        # List of combos to update
        combos = [
            getattr(self, "obj_combo_gm", None),
            getattr(self, "vina_receptor_combo", None),

            getattr(self, "obj_combo_apbs", None),
            getattr(self, "pocket_obj_combo", None),
            getattr(self, "ppi_obj_combo", None),
            getattr(self, "pl_obj_combo", None),
            getattr(self, "ll_obj_combo", None),
            getattr(self, "mut_obj_combo", None),
            getattr(self, "obj_combo_surf", None),  # Surface Analysis
            getattr(self, "obj_combo_sim1", None),  # Surface Similarity Object 1
            getattr(self, "ec_obj_combo", None),  # EC Analysis
            getattr(self, "ternary_obj_combo", None),  # Ternary Evaluation
        ]
        
        for cb in combos:
            if cb is not None:
                remember_selection(cb)
                cb.blockSignals(True)
                cb.clear()
                cb.addItems(names)
                restore_selection(cb)
                cb.blockSignals(False)

        sim2 = getattr(self, "obj_combo_sim2", None)
        if sim2 is not None:
            remember_selection(sim2)
            sim2.blockSignals(True)
            sim2.clear()
            sim2.addItem("(None - Single Surface)")
            sim2.addItems(names)
            restore_selection(sim2)
            sim1 = getattr(self, "obj_combo_sim1", None)
            sel1 = getattr(self, "sim_sel1", None)
            sel2 = getattr(self, "sim_sel2", None)
            if (
                sim1 is not None
                and sim1.currentText().strip() == sim2.currentText().strip()
                and (sel1 is None or sel2 is None or sel1.text().strip() == sel2.text().strip())
            ):
                sim2.setCurrentIndex(0)
            sim2.blockSignals(False)
        
        self.log(f"Refreshed {len(names)} objects")

    def update_enablement(self):
        """占位Method：按需Enable/Disable UI 控件。当前无需实现，子Class可覆写。"""
        pass

    def toggle_theme(self):
        """切换主题模式（保留兼容性，当前仅支持浅色主题）"""
        self._dark_mode = not self._dark_mode
        self.setup_style()
        self._update_content_bg()
        self.log(f"Switched to {'Dark' if self._dark_mode else 'Light'} mode")

    def setup_style(self):
        font = _get_platform_font()
        adj = _get_platform_adjustments()
        h = adj["min_height"]
        p = adj["padding"]
        fs = adj["font_size"]

        base_style = f"""
            * {{ font-family: {font}; font-size: {fs}px; }}
        """

        if self._dark_mode:
            self.setStyleSheet(base_style + f"""
                QDialog, QWidget {{ background-color: #0d1117; color: #c9d1d9; }}
                QGroupBox {{ border: 1px solid #30363d; border-radius: 6px; margin-top: 12px; padding: 16px 12px 12px 12px; }}
                QGroupBox::title {{ subcontrol-origin: margin; subcontrol-position: top left; padding: 0 8px; color: #58a6ff; font-weight: bold; }}
                QLabel {{ min-height: {h-2}px; padding: {p-2}px 0px; }}
                QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{ background-color: #010409; border: 1px solid #30363d; border-radius: 4px; padding: {p}px 6px; color: #c9d1d9; selection-background-color: #1f6feb; min-height: {h}px; }}
                QCheckBox {{ min-height: {h-2}px; padding: {p-2}px 0px; }}
                QPushButton {{ background-color: #21262d; border: 1px solid #30363d; border-radius: 4px; padding: {p}px 10px; color: #c9d1d9; min-height: {h}px; }}
                QPushButton:hover {{ background-color: #30363d; border-color: #8b949e; }}
                QPushButton:pressed {{ background-color: #282e33; }}
                QPushButton#highlight_btn {{ background-color: #238636; color: #ffffff; border: 1px solid #2ea043; }}
                QPushButton#highlight_btn:hover {{ background-color: #2ea043; }}
                QPushButton#primary_btn {{ background-color: #1f6feb; color: #ffffff; border: 1px solid #388bfd; }}
                QPushButton#primary_btn:hover {{ background-color: #388bfd; }}
                QScrollArea {{ border: none; background-color: transparent; }}
                QScrollBar:vertical {{ border: none; background: #0d1117; width: 10px; margin: 0px 0 0px 0; }}
                QScrollBar::handle:vertical {{ background: #30363d; min-height: 20px; border-radius: 5px; }}
                QListWidget {{ background: transparent; font-size: 14px; color: #c9d1d9; }}
                QListWidget::item {{ padding: 10px 16px; margin: 2px 8px; border-radius: 6px; }}
                QListWidget::item:selected {{ background: #3b82f6; color: white; }}
                QListWidget::item:hover:!selected {{ background: rgba(59, 130, 246, 0.1); }}
                #nav_widget {{ background-color: #0d1117; border-right: 1px solid #30363d; }}
                #tab_content {{ background-color: #161b22; }}
            """)
            if hasattr(self, 'nav_widget'):
                self.nav_widget.setStyleSheet("background-color: #0d1117; border-right: 1px solid #30363d;")
        else:
            self.setStyleSheet(base_style + f"""
                QDialog, QWidget {{ background-color: #ffffff; color: #1f2937; }}
                QGroupBox {{ border: 1px solid #e5e7eb; border-radius: 6px; margin-top: 12px; padding: 16px 12px 12px 12px; background-color: #f9fafb; }}
                QGroupBox::title {{ subcontrol-origin: margin; subcontrol-position: top left; padding: 0 8px; color: #3b82f6; font-weight: bold; }}
                QLabel {{ min-height: {h-2}px; padding: {p-2}px 0px; }}
                QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{ background-color: #ffffff; border: 1px solid #d1d5db; border-radius: 4px; padding: {p}px 6px; color: #1f2937; selection-background-color: #3b82f6; min-height: {h}px; }}
                QCheckBox {{ min-height: {h-2}px; padding: {p-2}px 0px; }}
                QPushButton {{ background-color: #f3f4f6; border: 1px solid #d1d5db; border-radius: 4px; padding: {p}px 10px; color: #374151; min-height: {h}px; }}
                QPushButton:hover {{ background-color: #e5e7eb; border-color: #9ca3af; }}
                QPushButton:pressed {{ background-color: #d1d5db; }}
                QPushButton#highlight_btn {{ background-color: #22c55e; color: #ffffff; border: 1px solid #16a34a; }}
                QPushButton#highlight_btn:hover {{ background-color: #16a34a; }}
                QPushButton#primary_btn {{ background-color: #3b82f6; color: #ffffff; border: 1px solid #2563eb; }}
                QPushButton#primary_btn:hover {{ background-color: #2563eb; }}
                QScrollArea {{ border: none; background-color: transparent; }}
                QScrollBar:vertical {{ border: none; background: #f3f4f6; width: 10px; margin: 0px 0 0px 0; }}
                QScrollBar::handle:vertical {{ background: #d1d5db; min-height: 20px; border-radius: 5px; }}
                QListWidget {{ background: transparent; font-size: 14px; color: #1f2937; }}
                QListWidget::item {{ padding: 10px 16px; margin: 2px 8px; border-radius: 6px; }}
                QListWidget::item:selected {{ background: #3b82f6; color: white; }}
                QListWidget::item:hover:!selected {{ background: rgba(59, 130, 246, 0.1); }}
                #nav_widget {{ background-color: #f1f5f9; border-right: 1px solid #e2e8f0; }}
                #tab_content {{ background-color: #ffffff; }}
                QTabWidget::pane {{ border: 1px solid #e5e7eb; background-color: #ffffff; }}
                QTabBar::tab {{ background-color: #f3f4f6; border: 1px solid #e5e7eb; padding: 6px 12px; }}
                QTabBar::tab:selected {{ background-color: #ffffff; border-bottom-color: #ffffff; }}
            """)
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
            show_message_box(self, "Error", "Environment checker not found.", icon_type="warning")

    # --- Simple Pages ---
    def create_welcome_tab(self) -> QWidget:
        """Create a modern, visually stunning welcome page following the scientific showcase concept."""
        from .qt_adapter import QGraphicsDropShadowEffect, QGridLayout, QScrollArea
        
        main_page = QWidget()
        main_page.setObjectName("welcome_page")
        main_layout = QVBoxLayout(main_page)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Direct layout without scroll
        content_container = QWidget()
        layout = QVBoxLayout(content_container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        main_layout.addWidget(content_container)

        # Hero Section (Top) - Compact
        hero_widget = HeroHeader()
        hero_layout = QVBoxLayout(hero_widget)
        hero_layout.setContentsMargins(40, 30, 40, 20)

        title_lbl = QLabel("GLINT")
        title_lbl.setStyleSheet("font-size: 64px; font-weight: 900; color: #ffffff; background: transparent; letter-spacing: 4px;")

        # Add a subtle glow/shadow to Title
        title_shadow = QGraphicsDropShadowEffect()
        title_shadow.setBlurRadius(20)
        title_shadow.setColor(Qt.GlobalColor.black)
        title_shadow.setOffset(0, 0)
        title_lbl.setGraphicsEffect(title_shadow)

        slogan_lbl = QLabel("A Toolkit for Molecular Glue Interface Analysis and Rational Design")
        slogan_lbl.setStyleSheet("font-size: 16px; color: #f8fafc; background: transparent; font-weight: 500; letter-spacing: 1px;")

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

        # Content Section - Compact
        content_area = QWidget()
        content_layout = QVBoxLayout(content_area)
        content_layout.setContentsMargins(30, 20, 30, 20)
        content_layout.setSpacing(12)

        # SECTION 1: WORKFLOW
        tools_lbl = QLabel("WORKFLOW")
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
                    border-radius: 12px; padding: 18px;
                }
                #feature_card:hover { background-color: rgba(59, 130, 246, 0.1); border: 1px solid #3b82f6; }
            """ if self._dark_mode else """
                #feature_card {
                    background-color: #ffffff; border: 1px solid #e2e8f0;
                    border-radius: 12px; padding: 18px;
                }
                #feature_card:hover { border-color: #3b82f6; background-color: #f8fafc; }
            """)

            shadow = QGraphicsDropShadowEffect()
            shadow.setBlurRadius(10)
            shadow.setXOffset(0)
            shadow.setYOffset(3)
            shadow.setColor(Qt.GlobalColor.black if self._dark_mode else Qt.GlobalColor.lightGray)
            card.setGraphicsEffect(shadow)

            c_layout = QHBoxLayout(card)
            c_layout.setSpacing(16)
            c_layout.setContentsMargins(16, 18, 16, 18)
            ico = QLabel(icon_text); ico.setStyleSheet("font-size: 28px; background: transparent;")
            text_container = QWidget()
            text_layout = QVBoxLayout(text_container)
            text_layout.setContentsMargins(0, 0, 0, 0)
            text_layout.setSpacing(6)
            t_lbl = QLabel(title); t_lbl.setStyleSheet("font-size: 15px; font-weight: bold; color: #3b82f6; background: transparent;")
            d_lbl = QLabel(desc); d_lbl.setWordWrap(True); d_lbl.setStyleSheet("font-size: 12px; color: #64748b; background: transparent;")
            text_layout.addWidget(t_lbl)
            text_layout.addWidget(d_lbl)

            c_layout.addWidget(ico)
            c_layout.addWidget(text_container, 1)
            card.mousePressEvent = lambda e: self.nav_list.setCurrentRow(index)
            return card

        grid_layout.addWidget(create_card("Target Discovery", "G-motif detection, surface analysis", "🎯", 1), 0, 0)
        grid_layout.addWidget(create_card("Hit Identification", "Pocket detection, docking integration", "🔍", 2), 0, 1)
        grid_layout.addWidget(create_card("Ternary Evaluation", "Neo-epitope mapping, interface scoring", "🔬", 3), 1, 0)
        grid_layout.addWidget(create_card("Lead Optimization", "PPI analysis, interaction profiling, mutation analysis", "⚡", 4), 1, 1)
        content_layout.addLayout(grid_layout)

        # Community Section
        community_section = QWidget()
        community_section.setStyleSheet("background: transparent;")
        community_layout = QVBoxLayout(community_section)
        community_layout.setContentsMargins(0, 20, 0, 0)
        community_layout.setSpacing(12)

        community_title = QLabel("COMMUNITY")
        community_title.setStyleSheet("font-size: 11px; font-weight: bold; color: #94a3b8; letter-spacing: 2px; background: transparent;")
        community_layout.addWidget(community_title)

        # Community buttons row
        btn_row = QHBoxLayout()
        btn_row.setSpacing(12)

        def create_community_btn(text, icon, url):
            btn = QPushButton(f"{icon}  {text}")
            btn.setCursor(Qt.CursorShape.PointingHandCursor)
            btn.setStyleSheet("""
                QPushButton {
                    background: #f1f5f9;
                    border: 1px solid #e2e8f0;
                    border-radius: 8px;
                    padding: 12px 20px;
                    font-size: 13px;
                    font-weight: 500;
                    color: #475569;
                }
                QPushButton:hover {
                    background: #e2e8f0;
                    border-color: #3b82f6;
                    color: #3b82f6;
                }
            """)
            btn.clicked.connect(lambda: QDesktopServices.openUrl(QUrl(url)))
            return btn

        btn_row.addWidget(create_community_btn("GitHub", "⭐", "https://github.com/Augus1999/GlueTK"))
        btn_row.addWidget(create_community_btn("Documentation", "📖", "https://github.com/Augus1999/GlueTK#readme"))
        btn_row.addWidget(create_community_btn("Report Issue", "🐛", "https://github.com/Augus1999/GlueTK/issues"))
        btn_row.addWidget(create_community_btn("Discussions", "💬", "https://github.com/Augus1999/GlueTK/discussions"))
        btn_row.addWidget(create_community_btn("Cite", "📝", "https://github.com/Augus1999/GlueTK#citation"))
        btn_row.addStretch()

        community_layout.addLayout(btn_row)
        content_layout.addWidget(community_section)

        layout.addWidget(content_area)
        layout.addStretch()
        return main_page

    def create_readme_tab(self) -> QWidget:
        """Create a comprehensive README tab with GUI documentation."""
        from .qt_adapter import QScrollArea, QTextBrowser
        
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
        
        <h1>🧬 GLINT User Guide</h1>
        <p>Welcome to GLINT - Professional Toolkit for Molecular Glue Interface Analysis and Rational Design</p>
        
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
                <li><b>Surface Analysis</b> - Analyze electrostatic/hydrophobic patches</li>
            </ul>
            <p><b>Options:</b></p>
            <ul>
                <li><b>Export Coordinates</b> - Check to save G-loop coordinates to CSV during detection</li>
            </ul>
        </div>
        
        <div class="feature-box">
            <div class="feature-title">Surface Similarity & Complementarity</div>
            <p>Compare binding site features using geometric and chemical descriptors.</p>
            <table>
                <tr><th>Parameter</th><th>Description</th><th>Default</th></tr>
                <tr><td>Analysis Type</td><td>Similarity or Complementarity (PPI)</td><td>Similarity</td></tr>
                <tr><td>Surface Method</td><td>auto/open3d/edtsurf</td><td>auto</td></tr>
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
                <tr><td>Pocket-constrained</td><td>Residue-based AIR restraints (requires active residues)</td></tr>
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
            <tr><td><code>glint_gui</code></td><td>Launch GLINT GUI</td></tr>
            <tr><td><code>find_crbn_g_motif obj, rmsd_cutoff=3.5</code></td><td>G-motif detection</td></tr>
            <tr><td><code>ppi_analyze obj, ['A'], ['B']</code></td><td>PPI interface analysis</td></tr>
            <tr><td><code>neo_epitope_find obj, ['A'], ['B'], LIG</code></td><td>Neo-epitope detection</td></tr>
            <tr><td><code>generate_2d_diagram obj, LIG</code></td><td>2D interaction diagram</td></tr>
            <tr><td><code>detect_pockets obj</code></td><td>Pocket detection</td></tr>
            <tr><td><code>batch_gmotif "6H0G,6H0F"</code></td><td>Batch G-motif</td></tr>
        </table>
        
        <h2><span class="emoji">📚</span>Resources</h2>
        <p>Click <b>Resources</b> button in the sidebar for databases and literature.</p>
        <p>GitHub: <a href="https://github.com/VesperChen01/GLINT">https://github.com/VesperChen01/GLINT</a></p>
        
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
        """Create a clean contact page with app info and maintainer cards."""
        from .qt_adapter import QScrollArea, QGridLayout

        main_page = QWidget()
        main_layout = QVBoxLayout(main_page)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Scroll Area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        content = QWidget()
        content.setStyleSheet("background-color: #ffffff;")
        content_layout = QVBoxLayout(content)
        content_layout.setContentsMargins(60, 50, 60, 50)
        content_layout.setSpacing(30)
        content_layout.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignHCenter)

        # === Header Section (居中) ===
        header_section = QWidget()
        header_layout = QVBoxLayout(header_section)
        header_layout.setSpacing(12)
        header_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)

        # App Icon (居中)
        logo_path = _get_logo_path()
        if logo_path:
            icon_lbl = QLabel()
            pix = QPixmap(logo_path)
            icon_lbl.setPixmap(pix.scaled(100, 110, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))
            icon_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            header_layout.addWidget(icon_lbl)

        # Title: Glint (H1, Bold, 高亮黑色)
        title_lbl = QLabel("GLue INTerface analyzer")
        title_lbl.setStyleSheet("font-size: 32px; font-weight: bold; color: #0f172a;")
        title_lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
        header_layout.addWidget(title_lbl)

        # Meta Info Row: Version号 Badge + Description (同行Display)
        meta_row = QWidget()
        meta_layout = QHBoxLayout(meta_row)
        meta_layout.setSpacing(12)
        meta_layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        meta_layout.setContentsMargins(0, 0, 0, 0)

        # Version号 Badge (蓝色描边/背景)
        version_badge = QLabel(f"v{__version__.lstrip('v')}")
        version_badge.setStyleSheet("""
            QLabel {
                background-color: #eff6ff;
                color: #2563eb;
                border: 1px solid #3b82f6;
                border-radius: 12px;
                padding: 4px 12px;
                font-size: 13px;
                font-weight: 600;
            }
        """)
        meta_layout.addWidget(version_badge)

        # Description文案
        desc_lbl = QLabel("Professional Toolkit for Molecular Glue Interface Analysis and Rational Design")
        desc_lbl.setStyleSheet("font-size: 14px; color: #64748b;")
        meta_layout.addWidget(desc_lbl)

        header_layout.addWidget(meta_row)
        content_layout.addWidget(header_section)

        # === Info Cards Section (两张等宽卡片横向排列) ===
        cards_container = QWidget()
        cards_layout = QHBoxLayout(cards_container)
        cards_layout.setSpacing(20)
        cards_layout.setContentsMargins(0, 20, 0, 0)

        # 卡片1: Maintainers
        card1 = self._create_contact_card(
            "Maintainers",
            [
                ("Roufen Chen", "Creator & Lead Maintainer", "12319021@zju.edu.cn"),
                ("Xinle Yang", "Co-Maintainer", "221124070197@zjut.edu.cn"),
            ]
        )
        cards_layout.addWidget(card1, 1)

        # 卡片2: Repository
        card2 = self._create_repo_card()
        cards_layout.addWidget(card2, 1)

        content_layout.addWidget(cards_container)
        content_layout.addStretch()

        scroll.setWidget(content)
        main_layout.addWidget(scroll)

        return main_page

    def _create_contact_card(self, title: str, contacts: list) -> QWidget:
        """Create a contact info card with multiple contacts."""
        card = QFrame()
        card.setStyleSheet("""
            QFrame {
                background-color: #f8fafc;
                border-radius: 12px;
            }
        """)

        layout = QVBoxLayout(card)
        layout.setContentsMargins(24, 20, 24, 20)
        layout.setSpacing(16)

        # Card Title
        title_lbl = QLabel(title)
        title_lbl.setStyleSheet("font-size: 16px; font-weight: bold; color: #1e293b;")
        layout.addWidget(title_lbl)

        # Contact Items
        for name, role, email in contacts:
            item = QWidget()
            item_layout = QVBoxLayout(item)
            item_layout.setContentsMargins(0, 0, 0, 0)
            item_layout.setSpacing(2)

            name_lbl = QLabel(name)
            name_lbl.setStyleSheet("font-size: 14px; font-weight: 600; color: #334155;")

            role_lbl = QLabel(role)
            role_lbl.setStyleSheet("font-size: 12px; color: #64748b;")

            email_lbl = QLabel(f"📧 {email}")
            email_lbl.setStyleSheet("font-size: 12px; color: #3b82f6;")

            item_layout.addWidget(name_lbl)
            item_layout.addWidget(role_lbl)
            item_layout.addWidget(email_lbl)
            layout.addWidget(item)

        layout.addStretch()
        return card

    def _create_repo_card(self) -> QWidget:
        """Create a repository info card."""
        card = QFrame()
        card.setStyleSheet("""
            QFrame {
                background-color: #f8fafc;
                border-radius: 12px;
            }
        """)

        layout = QVBoxLayout(card)
        layout.setContentsMargins(24, 20, 24, 20)
        layout.setSpacing(16)

        # Card Title
        title_lbl = QLabel("Repository")
        title_lbl.setStyleSheet("font-size: 16px; font-weight: bold; color: #1e293b;")
        layout.addWidget(title_lbl)

        # GitHub Info
        github_item = QWidget()
        github_layout = QVBoxLayout(github_item)
        github_layout.setContentsMargins(0, 0, 0, 0)
        github_layout.setSpacing(4)

        platform_lbl = QLabel("GitHub")
        platform_lbl.setStyleSheet("font-size: 14px; font-weight: 600; color: #334155;")

        url_lbl = QLabel("https://github.com/VesperChen01/GLINT")
        url_lbl.setStyleSheet("font-size: 12px; color: #3b82f6;")
        url_lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        github_layout.addWidget(platform_lbl)
        github_layout.addWidget(url_lbl)
        layout.addWidget(github_item)

        layout.addStretch()
        return card
