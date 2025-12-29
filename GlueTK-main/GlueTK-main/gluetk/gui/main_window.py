# -*- coding: utf-8 -*-
"""
GlueTK Main Window
Modularized version of the Unified GUI.
"""
import os
import sys
from typing import Optional, List, Dict, Tuple, Any

from .qt_adapter import QtCore, QtWidgets, QtGui, Qt, Signal, Slot, Property

if QtWidgets is None:
    raise RuntimeError("No Qt binding (PyQt5, PyQt6, PySide2, or PySide6) found.")

from .qt_adapter import (
    QDialog, QVBoxLayout, QHBoxLayout, QListWidget, QStackedWidget,
    QWidget, QPushButton, QLabel, QFrame, QTextEdit, QProgressBar, QMessageBox,
    QIcon, QPixmap, QColor, QBrush, QRadialGradient, QLinearGradient,
    QTimer, QSize, QSettings
)

from .utils import t, get_lang
from .tabs.target_discovery import TargetDiscoveryTab
from .tabs.hit_identification import HitIdentificationTab
from .tabs.lead_optimization import LeadOptimizationTab
from .tabs.batch_analysis import BatchAnalysisTab
from .tabs.ternary_evaluation import TernaryEvaluationTab

from .workers import AnalysisWorker, GMotifWorker

def _get_logo_path() -> Optional[str]:
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    logo_path = os.path.join(here, "assets", "logo.png")
    return logo_path if os.path.exists(logo_path) else None

def _get_welcome_bg_path() -> Optional[str]:
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    bg_path = os.path.join(here, "assets", "welcome_bg.png")
    return bg_path if os.path.exists(bg_path) else None

class HeroHeader(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedHeight(300)
        self.bg_path = _get_welcome_bg_path()
        
    def paintEvent(self, event):
        from .qt_adapter import QPainter, QLinearGradient
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        if self.bg_path:
            pix = QPixmap(self.bg_path)
            if not pix.isNull():
                scaled_pix = pix.scaled(self.size(), Qt.AspectRatioMode.KeepAspectRatioByExpanding, Qt.TransformationMode.SmoothTransformation)
                x = (self.width() - scaled_pix.width()) // 2
                y = (self.height() - scaled_pix.height()) // 2
                painter.drawPixmap(x, y, scaled_pix)
        else:
            grad = QLinearGradient(0, 0, self.width(), self.height())
            grad.setColorAt(0, QColor("#1e293b"))
            grad.setColorAt(1, QColor("#0f172a"))
            painter.fillRect(self.rect(), grad)
        painter.fillRect(self.rect(), QColor(15, 23, 42, 100))
        painter.setPen(QColor("#3b82f6"))
        painter.drawLine(0, self.height()-1, self.width(), self.height()-1)
        painter.end()

class GlueTKDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(t("title"))
        logo_path = _get_logo_path()
        if logo_path:
            self.setWindowIcon(QIcon(logo_path))
        self.setMinimumSize(1280, 720)
        self.resize(1280, 720)
        self.settings = QSettings("GlueTK", "GlueTK_App")
        self._dark_mode = False
        self._interactions = []
        self._gmotif_hits = []
        self._esp_maps = {}
        self.analysis_thread = None
        self.gmotif_thread = None
        self.ternary_thread = None
        self.build_ui()
        self.setup_style()
        self._update_content_bg()
        QTimer.singleShot(100, self.refresh_objects)
        self.log(t("log_ready"))

    def build_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)
        content_row = QHBoxLayout()
        content_row.setSpacing(0)
        nav_widget = self._create_nav_widget()
        content_row.addWidget(nav_widget, 0)
        self.content_stack = QStackedWidget()
        self.content_stack.addWidget(self.create_welcome_tab())
        self.target_tab = TargetDiscoveryTab(self)
        self.content_stack.addWidget(self.target_tab)
        self.ternary_tab = TernaryEvaluationTab(self)
        self.content_stack.addWidget(self.ternary_tab)
        self.hit_tab = HitIdentificationTab(self)
        self.content_stack.addWidget(self.hit_tab)
        self.lead_tab = LeadOptimizationTab(self)
        self.content_stack.addWidget(self.lead_tab)
        self.batch_tab = BatchAnalysisTab(self)
        self.content_stack.addWidget(self.batch_tab)
        self.content_stack.addWidget(self.create_readme_tab())
        self.content_stack.addWidget(self.create_contact_tab())
        content_row.addWidget(self.content_stack, 1)
        main_layout.addLayout(content_row)
        self.log_edit = QTextEdit()
        self.log_edit.setVisible(False)
        main_layout.addWidget(self.log_edit)
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        main_layout.addWidget(self.progress_bar)

    def _create_nav_widget(self) -> QWidget:
        nav_widget = QWidget()
        nav_widget.setObjectName("nav_widget")
        self.nav_widget = nav_widget
        nav_widget.setFixedWidth(240)
        nav_widget.setStyleSheet("background-color: #f1f5f9; border-right: 1px solid #e2e8f0;")
        layout = QVBoxLayout(nav_widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
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
        header = QWidget()
        h_layout = QHBoxLayout(header)
        h_layout.addStretch()
        h_layout.setContentsMargins(10, 5, 10, 5)
        layout.addWidget(header)
        self.nav_list = QListWidget()
        self.nav_list.setFrameShape(QFrame.Shape.NoFrame)
        self.nav_list.setStyleSheet("QListWidget { background: transparent; font-size: 14px; } QListWidget::item { padding: 12px; } QListWidget::item:selected { background: #3b82f6; color: white; border-radius: 4px; }")
        items = ["Welcome", "Target Discovery", "Ternary Evaluation", "Hit Identification", "Lead Optimization", "Batch Analysis"]
        self.nav_list.addItems(items)
        self.nav_list.setCurrentRow(0)
        self.nav_list.currentRowChanged.connect(self.on_nav_changed)
        layout.addWidget(self.nav_list)
        layout.addStretch()
        bottom_bar = QWidget()
        b_layout = QVBoxLayout(bottom_bar)
        b_layout.setSpacing(2)
        btn_readme = QPushButton("README")
        btn_readme.setFlat(True)
        btn_readme.clicked.connect(lambda: self.content_stack.setCurrentIndex(6))
        b_layout.addWidget(btn_readme)
        btn_contact = QPushButton("Contact Us")
        btn_contact.setFlat(True)
        btn_contact.clicked.connect(lambda: self.content_stack.setCurrentIndex(7))
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

    def on_error(self, msg: str):
        self.log(f"Error: {msg}")
        QMessageBox.critical(self, "Error", msg)
        self.progress_bar.setVisible(False)
        if hasattr(self, 'gm_btn'): self.gm_btn.setEnabled(True)
        if hasattr(self, 'ternary_run_all_btn'): self.ternary_run_all_btn.setEnabled(True)

    def refresh_objects(self):
        names = []
        try:
            from pymol import cmd
            names = cmd.get_names("objects") if hasattr(cmd, "get_names") else cmd.get_object_list()
        except:
            pass
        if not names: names = [t("no_object")]
        combos = [
            getattr(self, "obj_combo_gm", None), getattr(self, "obj_combo_apbs", None),
            getattr(self, "pocket_obj_combo", None), getattr(self, "ppi_obj_combo", None),
            getattr(self, "pl_obj_combo", None), getattr(self, "ll_obj_combo", None),
            getattr(self, "mut_obj_combo", None), getattr(self, "obj_combo_surf", None),
            getattr(self, "obj_combo_sim1", None), getattr(self, "obj_combo_sim2", None),
            getattr(self, "ec_obj_combo", None), getattr(self, "ternary_obj_combo", None),
        ]
        for cb in combos:
            if cb:
                cb.blockSignals(True)
                cb.clear()
                cb.addItems(names)
                cb.blockSignals(False)
        self.log(f"Refreshed {len(names)} objects")

    def update_enablement(self):
        pass

    def toggle_theme(self):
        self._dark_mode = not self._dark_mode
        self.setup_style()
        self._update_content_bg()

    def setup_style(self):
        if self._dark_mode:
            self.setStyleSheet("""
                QDialog, QWidget { background-color: #0d1117; color: #c9d1d9; }
                QGroupBox { border: 1px solid #30363d; border-radius: 6px; margin-top: 12px; padding-top: 10px; }
                QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left; padding: 0 5px; color: #58a6ff; font-weight: bold; }
                QLineEdit, QComboBox, QSpinBox { background-color: #010409; border: 1px solid #30363d; border-radius: 4px; padding: 4px; color: #c9d1d9; }
                QPushButton { background-color: #21262d; border: 1px solid #30363d; border-radius: 4px; padding: 6px 12px; color: #c9d1d9; }
                QPushButton:hover { background-color: #30363d; }
                QPushButton#highlight_btn { background-color: #238636; color: #ffffff; }
                QPushButton#primary_btn { background-color: #3b82f6; color: #ffffff; }
                QListWidget { background: transparent; font-size: 14px; color: #c9d1d9; }
                QListWidget::item { padding: 12px; }
                QListWidget::item:selected { background: #3b82f6; color: white; }
            """)
        else:
            self.setStyleSheet("""
                QDialog, QWidget { background-color: #ffffff; color: #1f2937; }
                QGroupBox { border: 1px solid #e5e7eb; border-radius: 6px; margin-top: 12px; padding-top: 10px; background-color: #f9fafb; }
                QGroupBox::title { subcontrol-origin: margin; subcontrol-position: top left; padding: 0 5px; color: #3b82f6; font-weight: bold; }
                QLineEdit, QComboBox, QSpinBox { background-color: #ffffff; border: 1px solid #d1d5db; border-radius: 4px; padding: 4px; color: #1f2937; }
                QPushButton { background-color: #f3f4f6; border: 1px solid #d1d5db; border-radius: 4px; padding: 6px 12px; color: #374151; }
                QPushButton:hover { background-color: #e5e7eb; }
                QPushButton#highlight_btn { background-color: #22c55e; color: #ffffff; }
                QPushButton#primary_btn { background-color: #3b82f6; color: #ffffff; }
                QListWidget { background: transparent; font-size: 14px; color: #1f2937; }
                QListWidget::item { padding: 12px; }
                QListWidget::item:selected { background: #3b82f6; color: white; }
                QTabWidget::pane { border: 1px solid #e5e7eb; background-color: #ffffff; }
                QTabBar::tab { background-color: #f3f4f6; border: 1px solid #e5e7eb; padding: 6px 12px; }
                QTabBar::tab:selected { background-color: #ffffff; }
            """)

    def _update_content_bg(self):
        bg = "#161b22" if self._dark_mode else "#ffffff"
        for name in ("_target_scroll_content", "_lead_scroll_content", "_hit_scroll_content"):
            w = getattr(self, name, None)
            if w: w.setStyleSheet(f"#scroll_content {{ background-color: {bg}; }}")

    def check_environment(self):
        try:
            from ..env_checker import ensure_dependencies
            ensure_dependencies(silent=False)
        except ImportError:
            QMessageBox.warning(self, "Error", "Environment checker not found.")

    def create_welcome_tab(self) -> QWidget:
        from .qt_adapter import QGraphicsDropShadowEffect, QGridLayout, QScrollArea
        main_page = QWidget()
        main_layout = QVBoxLayout(main_page)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        content_container = QWidget()
        layout = QVBoxLayout(content_container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        scroll.setWidget(content_container)
        main_layout.addWidget(scroll)
        hero_widget = HeroHeader()
        hero_layout = QVBoxLayout(hero_widget)
        hero_layout.setContentsMargins(40, 40, 40, 40)
        title_lbl = QLabel("GlueTK")
        title_lbl.setStyleSheet("font-size: 80px; font-weight: 900; color: #ffffff; background: transparent;")
        slogan_lbl = QLabel("Professional Toolkit for Molecular Glue Engineering")
        slogan_lbl.setStyleSheet("font-size: 24px; color: #f8fafc; background: transparent;")
        hero_layout.addWidget(title_lbl)
        hero_layout.addWidget(slogan_lbl)
        hero_layout.addStretch()
        layout.addWidget(hero_widget)
        content_area = QWidget()
        content_layout = QVBoxLayout(content_area)
        content_layout.setContentsMargins(40, 40, 40, 40)
        content_layout.setSpacing(40)
        tools_lbl = QLabel("CORE MODULES")
        tools_lbl.setStyleSheet("font-size: 14px; font-weight: bold; color: #3b82f6;")
        content_layout.addWidget(tools_lbl)
        grid_layout = QGridLayout()
        grid_layout.setSpacing(20)
        def create_card(title, desc, icon, idx):
            card = QFrame()
            card.setStyleSheet("QFrame { background: #ffffff; border: 1px solid #e2e8f0; border-radius: 12px; padding: 20px; } QFrame:hover { border-color: #3b82f6; }")
            card.setCursor(Qt.CursorShape.PointingHandCursor)
            c_layout = QVBoxLayout(card)
            c_layout.addWidget(QLabel(icon))
            t_lbl = QLabel(title); t_lbl.setStyleSheet("font-size: 18px; font-weight: bold; color: #3b82f6;")
            c_layout.addWidget(t_lbl)
            d_lbl = QLabel(desc); d_lbl.setWordWrap(True); d_lbl.setStyleSheet("color: #64748b;")
            c_layout.addWidget(d_lbl)
            c_layout.addStretch()
            card.mousePressEvent = lambda e: self.nav_list.setCurrentRow(idx)
            return card
        grid_layout.addWidget(create_card("Target Discovery", "Identify promising residues.", "🎯", 1), 0, 0)
        grid_layout.addWidget(create_card("Ternary Evaluation", "Evaluate interface, ligand, geometry.", "🔬", 2), 0, 1)
        grid_layout.addWidget(create_card("Hit Identification", "Analyze interactions.", "🔍", 3), 1, 0)
        grid_layout.addWidget(create_card("Lead Optimization", "Fine-tune molecules.", "⚡", 4), 1, 1)
        content_layout.addLayout(grid_layout)
        layout.addWidget(content_area)
        layout.addStretch()
        return main_page

    def _load_case(self, pdb_id: str):
        try:
            from pymol import cmd
            cmd.reinitialize()
            cmd.fetch(pdb_id)
            cmd.show_as("cartoon")
            cmd.zoom()
            QMessageBox.information(self, "Case Loaded", f"Loaded {pdb_id}")
            self.refresh_objects()
        except Exception as e:
            self.on_error(str(e))

    def create_readme_tab(self) -> QWidget:
        from .qt_adapter import QTextBrowser
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(0, 0, 0, 0)
        browser = QTextBrowser()
        browser.setReadOnly(True)
        browser.setOpenExternalLinks(True)
        browser.setHtml("""
        <html><body style="font-family: Arial; padding: 30px;">
        <h1 style="color: #1e40af;">GlueTK User Guide</h1>
        <h2 style="color: #3b82f6;">Target Discovery</h2>
        <p>Identify G-motif patterns for molecular glue targets.</p>
        <h2 style="color: #3b82f6;">Ternary Evaluation</h2>
        <p>Evaluate ternary complexes with Interface, Ligand, and Geometry modules.</p>
        <h2 style="color: #3b82f6;">Hit Identification</h2>
        <p>Detect binding pockets and perform docking.</p>
        <h2 style="color: #3b82f6;">Lead Optimization</h2>
        <p>Analyze PPI, protein-ligand interactions, and mutations.</p>
        </body></html>
        """)
        layout.addWidget(browser)
        return w

    def create_contact_tab(self) -> QWidget:
        from .qt_adapter import QScrollArea
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        content = QWidget()
        content.setStyleSheet("background-color: #ffffff;")
        c_layout = QVBoxLayout(content)
        c_layout.setContentsMargins(60, 50, 60, 50)
        c_layout.setSpacing(30)
        title = QLabel("Contact Us")
        title.setStyleSheet("font-size: 36px; font-weight: 800; color: #1e293b;")
        c_layout.addWidget(title)
        info = QLabel("Email: 12319021@zju.edu.cn\nGitHub: VesperChen01/GlueTK")
        info.setStyleSheet("font-size: 16px; color: #64748b;")
        c_layout.addWidget(info)
        c_layout.addStretch()
        scroll.setWidget(content)
        layout.addWidget(scroll)
        return w

    def _show_resources(self):
        QMessageBox.information(self, "Resources",
            "Key databases:\n• PROTAC-DB\n• Open Targets\n• PDB\n\n"
            "Key structures:\n• 6H0G - CRBN-CC885-GSPT1\n• 5FQD - CRBN-Lenalidomide-CK1α")
