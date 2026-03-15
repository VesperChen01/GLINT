# -*- coding: utf-8 -*-
"""
GLINT Unified GUI - PyMOL Plugin for Molecular Glue Analysis
- Molecular Glue vs PROTAC Classification
- PPI Interface & Neo-Epitope Detection
- G-Motif Recognition & CRBN Analysis
- Interaction Analysis & Visualization
- Electrostatics (APBS/Quick) + Export

Features:
1) Non-modal GUI integration with PyMOL
2) Advanced molecular glue analysis with PPI quantification
3) Neo-substrate epitope detection
4) CRBN G-motif recognition with glue binding validation
5) Quick electrostatics (no APBS required) + full APBS support
"""

from __future__ import annotations
import os, sys, csv, importlib.util
from typing import List, Dict, Any, Tuple, Optional

# Open Targets API module (integrated)
try:
    from .open_targets_api import (
        search_disease,
        get_disease_targets,
        enrich_targets_with_e3_scores
    )
    from .disease_config import DEFAULT_OUTPUT_DIR, E3_LIGASES
except ImportError:
    try:
        from open_targets_api import (
            search_disease,
            get_disease_targets,
            enrich_targets_with_e3_scores
        )
        from disease_config import DEFAULT_OUTPUT_DIR, E3_LIGASES
    except ImportError:
        search_disease = None
        get_disease_targets = None
        enrich_targets_with_e3_scores = None
        DEFAULT_OUTPUT_DIR = "./disease_data"
        E3_LIGASES = ["CRBN", "VHL", "MDM2", "XIAP"]

# -------- Qt compatibility (prefer PyQt5) --------
QT_LIB = None
try:
    from PyQt5.QtCore import Qt, QThread, pyqtSignal, QLocale, QTimer
    from PyQt5.QtWidgets import (
        QApplication, QDialog, QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QLabel,
        QLineEdit, QPushButton, QCheckBox, QComboBox, QFileDialog, QGroupBox,
        QFormLayout, QMessageBox, QTextEdit, QProgressBar, QFrame, QTabWidget,
        QTableWidget, QTableWidgetItem, QSizePolicy, QGridLayout, QListWidget, QStackedWidget,
        QSpinBox, QScrollArea
    )
    from PyQt5.QtGui import QIcon, QPixmap
    QT_LIB = "PyQt5"
except Exception:
    try:
        from PyQt6.QtCore import Qt, QThread, pyqtSignal, QLocale, QTimer
        from PyQt6.QtWidgets import (
            QApplication, QDialog, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
            QLineEdit, QPushButton, QCheckBox, QComboBox, QFileDialog, QGroupBox,
            QFormLayout, QMessageBox, QTextEdit, QProgressBar, QFrame, QTabWidget,
            QTableWidget, QTableWidgetItem, QSizePolicy, QGridLayout, QListWidget, QStackedWidget,
            QSpinBox, QScrollArea
        )
        from PyQt6.QtGui import QIcon, QPixmap
        QT_LIB = "PyQt6"
    except Exception as e:
        raise RuntimeError("PyQt5 or PyQt6 is required") from e

# -------- Logo path --------
def _get_logo_path() -> Optional[str]:
    """Get the path to the logo file"""
    here = os.path.dirname(os.path.abspath(__file__))
    logo_path = os.path.join(here, "assets", "logo.png")
    if os.path.exists(logo_path):
        return logo_path
    return None

# -------- Module import helper (relative → absolute → dynamic) --------
def _dynamic_load_by_filenames(names: List[str], symbol: str):
    here = os.path.dirname(os.path.abspath(__file__))
    for nm in names:
        path = os.path.join(here, nm)
        if os.path.exists(path):
            spec = importlib.util.spec_from_file_location(f"_dyn_{os.path.splitext(nm)[0]}", path)
            if spec and spec.loader:
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)  # type: ignore
                if hasattr(mod, symbol):
                    return getattr(mod, symbol)
    return None

def _import_helpers():
    """Returns: highlight_csv_residues, highlight_gmotif_loops, analyze_pdb_interactions, find_crbn_g_motif, render_interactions_beautifully, generate_2d_interaction_diagram, analyze_protein_ligand_interactions, visualize_protein_ligand_3d, generate_interaction_network_plot, analyze_ternary_complex, analyze_atom_pair_interactions, visualize_atom_pairs, analyze_ligand_ligand_interactions"""
    # Try package-relative import
    try:
        from .highlight_residues import highlight_csv_residues, highlight_gmotif_loops  # type: ignore
        from .interaction_analyzer import (  # type: ignore
            analyze_pdb_interactions, render_interactions_beautifully,
            analyze_protein_ligand_interactions, visualize_protein_ligand_3d,
            generate_interaction_network_plot, analyze_ternary_complex,
            analyze_atom_pair_interactions, visualize_atom_pairs,
            analyze_protein_nucleic_interactions
        )
        from .interaction_2d_plot import generate_2d_interaction_diagram  # type: ignore
        try:
            from .ligand_ligand_analyzer import analyze_ligand_ligand_interactions  # type: ignore
        except ImportError:
            analyze_ligand_ligand_interactions = None

        find_crbn_g_motif = None
        try:
            try:
                from .g_motif_analyzer import find_crbn_g_motif  # type: ignore
            except Exception:
                from .g_motif import find_crbn_g_motif  # type: ignore
        except Exception:
            find_crbn_g_motif = _dynamic_load_by_filenames(
                ["g_motif_analyzer.py", "g_motif.py", "g-motif.py"], "find_crbn_g_motif"
            )
        return highlight_csv_residues, highlight_gmotif_loops, analyze_pdb_interactions, find_crbn_g_motif, render_interactions_beautifully, generate_2d_interaction_diagram, analyze_protein_ligand_interactions, visualize_protein_ligand_3d, generate_interaction_network_plot, analyze_ternary_complex, analyze_atom_pair_interactions, visualize_atom_pairs, analyze_ligand_ligand_interactions, analyze_protein_nucleic_interactions
    except Exception:
        pass
    # Absolute import from same directory
    here = os.path.dirname(os.path.abspath(__file__))
    if here and here not in sys.path:
        sys.path.insert(0, here)
    try:
        from highlight_residues import highlight_csv_residues, highlight_gmotif_loops  # type: ignore
        from interaction_analyzer import (  # type: ignore
            analyze_pdb_interactions, render_interactions_beautifully,
            analyze_protein_ligand_interactions, visualize_protein_ligand_3d,
            generate_interaction_network_plot, analyze_ternary_complex,
            analyze_atom_pair_interactions, visualize_atom_pairs,
            analyze_protein_nucleic_interactions
        )
        from interaction_2d_plot import generate_2d_interaction_diagram  # type: ignore
        try:
            from ligand_ligand_analyzer import analyze_ligand_ligand_interactions  # type: ignore
        except ImportError:
            analyze_ligand_ligand_interactions = None
    except Exception as e:
        raise ModuleNotFoundError(
            "highlight_residues / interaction_analyzer not found; ensure they are in the same directory as unified_gui.py or inside the package."
        ) from e
    find_crbn_g_motif = None
    try:
        from g_motif_analyzer import find_crbn_g_motif  # type: ignore
    except Exception:
        try:
            from g_motif import find_crbn_g_motif  # type: ignore
        except Exception:
            find_crbn_g_motif = _dynamic_load_by_filenames(
                ["g_motif_analyzer.py", "g_motif.py", "g-motif.py"], "find_crbn_g_motif"
            )
    return highlight_csv_residues, highlight_gmotif_loops, analyze_pdb_interactions, find_crbn_g_motif, render_interactions_beautifully, generate_2d_interaction_diagram, analyze_protein_ligand_interactions, visualize_protein_ligand_3d, generate_interaction_network_plot, analyze_ternary_complex, analyze_atom_pair_interactions, visualize_atom_pairs, analyze_ligand_ligand_interactions, analyze_protein_nucleic_interactions

highlight_csv_residues, highlight_gmotif_loops, analyze_pdb_interactions, find_crbn_g_motif, render_interactions_beautifully, generate_2d_interaction_diagram, analyze_protein_ligand_interactions, visualize_protein_ligand_3d, generate_interaction_network_plot, analyze_ternary_complex, analyze_atom_pair_interactions, visualize_atom_pairs, analyze_ligand_ligand_interactions, analyze_protein_nucleic_interactions = _import_helpers()

# -------- Dependency check --------
def _check_and_install_deps():
    """Check and install dependencies, called at GUI startup"""
    try:
        # Try to import env_checker module
        here = os.path.dirname(os.path.abspath(__file__))
        sys.path.insert(0, here)
        try:
            from .env_checker import ensure_dependencies
        except ImportError:
            from env_checker import ensure_dependencies
        
        return ensure_dependencies(silent=False)
    except Exception as e:
        print(f"[GLINT] Dependency check failed: {e}")
        return False

# -------- Language & Text --------
LANG_FORCE = "en"  # Force English for all GUI
def get_lang() -> str:
    return "en"  # Always return English

T = {
    "title": {"zh": "GLINT 统一 GUI", "en": "GLINT - Molecular Glue Analyzer"},
    "tab_gmotif": {"zh": "G-Motif 识别", "en": "G-Motif Detection"},
    "tab_analysis": {"zh": "相互作用分析与高亮", "en": "Interaction Analysis & Highlight"},
    "tab_apbs": {"zh": "静电势（APBS/Quick）", "en": "Electrostatics (APBS/Quick)"},
    "grp_csv": {"zh": "CSV 高亮", "en": "CSV Highlight"},
    "csv_path": {"zh": "CSV FilePath", "en": "CSV File Path"},
    "browse": {"zh": "浏览…", "en": "Browse…"},
    "target_obj": {"zh": "目标对象", "en": "Target Object"},
    "refresh": {"zh": "Refresh", "en": "Refresh"},
    "btn_highlight": {"zh": "高亮Display", "en": "Highlight"},
    "btn_clear": {"zh": "Clear高亮", "en": "Clear"},
    "grp_analysis": {"zh": "相互作用分析", "en": "Interaction Analysis"},
    "pdb_file": {"zh": "PDB File（可选）", "en": "PDB File (optional)"},
    "output_csv": {"zh": "输出 CSV（可选）", "en": "Output CSV (optional)"},
    "btn_analyze": {"zh": "Start分析", "en": "Start"},
    "btn_render_interactions": {"zh": "一Key渲染（分析+美化+PNG）", "en": "Render (Analyze + Beautify + PNG)"},
    "right_results": {"zh": "Results与日志", "en": "Results & Logs"},
    "table_header": {
        "zh": ["链1", "残基1", "链2", "残基2", "距离", "相互作用"],
        "en": ["Chain1", "Residue1", "Chain2", "Residue2", "Distance", "Interaction"],
    },
    "table_header_gmotif": {
        "zh": ["链", "序列", "起始", "结束", "RMSD (Å)", "Type"],
        "en": ["Chain", "Sequence", "Start", "End", "RMSD (Å)", "Type"],
    },
    "no_object": {"zh": "(无对象)", "en": "(No object)"},
    "select_csv": {"zh": "Select CSV File", "en": "Select CSV"},
    "select_pdb": {"zh": "Select PDB File", "en": "Select PDB"},
    "select_outcsv": {"zh": "Select输出 CSV", "en": "Select output CSV"},
    "log_ready": {"zh": "就绪", "en": "Ready"},
    "log_highlight_ok": {"zh": "高亮Completed", "en": "Highlight done"},
    "log_clear_ok": {"zh": "已Clear高亮", "en": "Cleared"},
    "log_start": {"zh": "Start相互作用分析…", "en": "Starting interaction analysis…"},
    "log_done": {"zh": "Completed，发现 {n} 条记录", "en": "Done: {n} records"},
    "log_error": {"zh": "Error：{msg}", "en": "Error: {msg}"},
    "btn_close": {"zh": "Close", "en": "Close"},
    "btn_load_csv_to_table": {"zh": "载入到表格", "en": "Load to Table"},
    # G-Motif
    "grp_gmotif": {"zh": "G-Motif（CRBN G-loop）识别", "en": "G-Motif (CRBN G-loop) Detection"},
    "rmsd": {"zh": "RMSD 阈Value (Å)", "en": "RMSD cutoff (Å)"},
    "require_gly": {"zh": "第6位必须为 Gly", "en": "Require Gly at pos6"},
    "btn_gmotif": {"zh": "Start识别", "en": "Detect"},
    "btn_gmotif_render": {"zh": "一Key渲染（G-Motif + 电势 + PNG）", "en": "Render (G-Motif + ESP + PNG)"},
    # APBS/Quick + Export
    "grp_apbs": {"zh": "静电势Display", "en": "Electrostatics Display"},
    "apbs_target": {"zh": "目标对象", "en": "Target Object"},
    "apbs_grid": {"zh": "网格间距 (Å)", "en": "Grid spacing (Å)"},
    "apbs_range": {"zh": "色阶范围 (kT/e)", "en": "Color range (kT/e)"},
    "btn_quick": {"zh": "快速静电图（Coulomb）", "en": "Quick (Coulomb)"},
    "btn_apbs": {"zh": "运行 APBS（若可用）", "en": "Run APBS (if available)"},
    "grp_export": {"zh": "Export", "en": "Export"},
    "img_w": {"zh": "宽degrees (px)", "en": "Width (px)"},
    "img_h": {"zh": "高degrees (px)", "en": "Height (px)"},
    "img_dpi": {"zh": "DPI", "en": "DPI"},
    "img_bg": {"zh": "背景", "en": "Background"},
    "bg_white": {"zh": "白色", "en": "White"},
    "bg_trans": {"zh": "透明", "en": "Transparent"},
    "raytrace": {"zh": "光线追踪（高质量）", "en": "Ray trace (high quality)"},
    "btn_viewport": {"zh": "using当前视口尺寸", "en": "Use current viewport"},
    "btn_export_png": {"zh": "Export PNG", "en": "Export PNG"},
    "btn_export_dx": {"zh": "Export DX 网格", "en": "Export DX"},
}

def t(key: str) -> str:
    lang = get_lang()
    d = T.get(key)
    if isinstance(d, dict):
        return d.get(lang, list(d.values())[0])
    return str(d)

# -------- Worker threads --------
class AnalysisWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)
    def __init__(self, obj_name: str, pdb_file: Optional[str], output_csv: Optional[str]):
        super().__init__()
        self.obj_name = obj_name
        self.pdb_file = pdb_file
        self.output_csv = output_csv
    def run(self):
        try:
            self.progress.emit(t("log_start"))
            interactions = analyze_pdb_interactions(
                obj_name=self.obj_name,
                pdb_file=self.pdb_file,
                only_between_chains=True,  # Enable inter-chain analysis by default
                output_csv=self.output_csv,
                auto_highlight=True,  # Enable auto-highlight
            )
            self.progress.emit(t("log_done").format(n=len(interactions)))
            self.finished.emit(interactions)
        except Exception as e:
            self.error.emit(str(e))

class PNAnalysisWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)
    def __init__(self, obj_name, nucleic_chains, protein_chains, output_csv, distance_cutoff, pdb_file=None):
        super().__init__()
        self.obj_name = obj_name
        self.nucleic_chains = nucleic_chains
        self.protein_chains = protein_chains
        self.output_csv = output_csv
        self.distance_cutoff = distance_cutoff
        self.pdb_file = pdb_file
        
    def run(self):
        try:
            self.progress.emit("Starting Protein-Nucleic Acid Analysis...")
            result = analyze_protein_nucleic_interactions(
                obj_name=self.obj_name,
                nucleic_chains=self.nucleic_chains,
                protein_chains=self.protein_chains,
                output_csv=self.output_csv,
                distance_cutoff=self.distance_cutoff,
                pdb_file=self.pdb_file
            )
            interactions = result.get("interactions", [])
            self.progress.emit(f"Analysis complete. Found {len(interactions)} interactions.")
            self.finished.emit(interactions)
        except Exception as e:
            self.error.emit(str(e))

class GMotifWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, str)
    error = pyqtSignal(str)
    def __init__(self, obj_name: str, pdb_file: Optional[str], rmsd: float, require_gly: bool,
                 out_csv: Optional[str], template_mode: str, template_sel: Optional[str], template_builtin: Optional[str],
                 exclude_proline: bool = True, check_surface_exposure: bool = True, min_sasa: float = 15.0):
        super().__init__()
        self.obj_name = obj_name; self.pdb_file = pdb_file
        self.rmsd = rmsd; self.require_gly = require_gly; self.out_csv = out_csv
        self.template_mode = template_mode; self.template_sel = template_sel; self.template_builtin = template_builtin
        self.exclude_proline = exclude_proline
        self.check_surface_exposure = check_surface_exposure
        self.min_sasa = min_sasa
    def run(self):
        try:
            if find_crbn_g_motif is None:
                raise RuntimeError("find_crbn_g_motif not found; ensure g_motif_analyzer.py exists in the plugin directory.")
            self.progress.emit("[G-Motif] " + ("Start识别…" if get_lang()=="zh" else "Detecting…"))
            out_csv_path = self.out_csv
            if not out_csv_path:
                import tempfile, os
                fd, out_csv_path = tempfile.mkstemp(suffix="_gmotif.csv"); os.close(fd)
            hits = find_crbn_g_motif(
                obj_name=self.obj_name,
                pdb_file=self.pdb_file,
                template_mode=self.template_mode,
                template_sel=self.template_sel,
                template_builtin=self.template_builtin,
                rmsd_cutoff=float(self.rmsd),
                out_csv=out_csv_path,
                auto_highlight=1,  # Enable auto-highlight
                require_gly_pos6=bool(self.require_gly),
                exclude_proline=bool(self.exclude_proline),
                check_surface_exposure=bool(self.check_surface_exposure),
                min_sasa_per_residue=float(self.min_sasa),
            ) or []
            self.progress.emit("[G-Motif] " + (f"Completed，命中 {len(hits)} 条" if get_lang()=="zh" else f"Done, {len(hits)} hits"))
            self.finished.emit(hits, out_csv_path)
        except Exception as e:
            self.error.emit(str(e))

# -------- Main dialog --------
class GLINTDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(t("title"))
        
        # Set window icon
        logo_path = _get_logo_path()
        if logo_path:
            self.setWindowIcon(QIcon(logo_path))
        
        # Set minimum window size (allow user resizing)
        min_w, min_h = 1280, 720
        self.setMinimumSize(min_w, min_h)  # Allow resizing
        self.resize(min_w, min_h)  # Initial size
        
        # Center on screen
        try:
            from PyQt5.QtGui import QGuiApplication as _QGA  # type: ignore
        except Exception:
            try:
                from PyQt6.QtGui import QGuiApplication as _QGA  # type: ignore
            except Exception:
                _QGA = None
        if _QGA:
            scr = _QGA.primaryScreen()
            if scr:
                geom = scr.availableGeometry()
                x = (geom.width() - min_w) // 2
                y = (geom.height() - min_h) // 2
                self.move(x, y)

        self.analysis_thread: Optional[AnalysisWorker] = None
        self.gmotif_thread: Optional[GMotifWorker] = None

        self._interactions: List[Dict[str, Any]] = []
        self._gmotif_hits: List[Tuple] = []
        self._last_gmotif_csv: Optional[str] = None
        self.settings = QSettings("GLINT", "GLINT_App")
        self._dark_mode: bool = False  # Default light theme
        self._ui_scale: float = 1.0   # Auto-scaling ratio

        # Pre-create log_edit and progress_bar (before build_ui)
        # Log panel removed, but keep invisible controls for log API compatibility
        # Create invisible log_edit to prevent errors
        self.log_edit = QTextEdit()
        self.log_edit.setVisible(False)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)

        # Check and install dependencies
        _check_and_install_deps()
        
        self.build_ui()
        self.setup_style()
        # Initialize theme icon - using Unicode symbols
        # self.theme_toggle_btn.setText("☾" if self._dark_mode else "☀")
        # Disable auto-scaling to maintain fixed height
        # self.apply_auto_scaling()

        # ❗Key fix: defer first PyMOL call to avoid blocking during construction
        QTimer.singleShot(0, self.refresh_objects)

        self.update_enablement()
        self.update_modules_button_style()  # Initialize Modules button color
        self.log(t("log_ready"))

    # Uniformly adjust spacing and margins for all layouts, reducing clutter and improving consistency
    def _tune_layouts(self, widget: QWidget):
        def _tune_layout_obj(lay):
            if lay is None:
                return
            try:
                if isinstance(lay, (QVBoxLayout, QHBoxLayout)):
                    lay.setSpacing(12)
                    lay.setContentsMargins(12, 12, 12, 12)
                elif isinstance(lay, QGridLayout):
                    lay.setHorizontalSpacing(12)
                    lay.setVerticalSpacing(10)
                    lay.setContentsMargins(12, 12, 12, 12)
                elif isinstance(lay, QFormLayout):
                    try:
                        lay.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
                        lay.setFormAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
                    except Exception:
                        # PyQt5 enum names
                        lay.setLabelAlignment(Qt.AlignRight)
                        lay.setFormAlignment(Qt.AlignLeft | Qt.AlignTop)
                    try:
                        lay.setRowWrapPolicy(QFormLayout.RowWrapPolicy.DontWrapRows)
                        lay.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
                    except Exception:
                        pass
                    lay.setHorizontalSpacing(12)
                    lay.setVerticalSpacing(10)
                    lay.setContentsMargins(12, 12, 12, 12)
            except Exception:
                pass

        layout = widget.layout()
        if layout:
            _tune_layout_obj(layout)

            # Iterate child items, recursively tune (set sub-layouts directly, dive into child widget layouts)
            for i in range(layout.count()):
                item = layout.itemAt(i)
                if not item:
                    continue
                try:
                    sub_layout = item.layout()
                except Exception:
                    sub_layout = None
                if sub_layout is not None:
                    _tune_layout_obj(sub_layout)
                    # Continue深入子布局的子项
                    try:
                        for j in range(sub_layout.count()):
                            sub_item = sub_layout.itemAt(j)
                            if sub_item and sub_item.widget():
                                self._tune_layouts(sub_item.widget())
                    except Exception:
                        pass
                # Child widgets (may have internal layouts)
                w = item.widget()
                if w is not None:
                    self._tune_layouts(w)

    # --- UI 结构 ---
    def build_ui(self):
        # ========== Main layout ==========
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(15, 15, 15, 15)  # Comfortable outer margins
        main_layout.setSpacing(12)  # Comfortable spacing

        # Horizontal split: navigation | content (fullscreen)
        content_row = QHBoxLayout()
        content_row.setSpacing(10)

        # ========== Left: Navigation List ==========
        nav_widget = QWidget()
        nav_widget.setObjectName("nav_widget")
        nav_layout = QVBoxLayout(nav_widget)
        nav_layout.setSpacing(10)
        nav_layout.setContentsMargins(8, 8, 8, 8)

        # Logo at top of navigation
        logo_path = _get_logo_path()
        if logo_path:
            logo_container = QHBoxLayout()
            logo_container.setContentsMargins(0, 4, 0, 8)
            logo_label = QLabel()
            logo_label.setObjectName("nav_logo")
            logo_pixmap = QPixmap(logo_path)
            # Scale logo to fit navigation width (max 80px height)
            try:
                # PyQt6 style
                scaled_pixmap = logo_pixmap.scaledToHeight(80, Qt.TransformationMode.SmoothTransformation)
                logo_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            except AttributeError:
                # PyQt5 style
                scaled_pixmap = logo_pixmap.scaledToHeight(80, Qt.SmoothTransformation)
                logo_label.setAlignment(Qt.AlignCenter)
            logo_label.setPixmap(scaled_pixmap)
            logo_container.addStretch()
            logo_container.addWidget(logo_label)
            logo_container.addStretch()
            nav_layout.addLayout(logo_container)

        # Header: Theme icon + Modules label (centered)
        nav_header = QHBoxLayout()
        nav_header.setSpacing(6)
        
        # Theme toggle button - REMOVED for specific white-only request
        # self.theme_toggle_btn = QPushButton("")  # Icon set after setup_style
        # self.theme_toggle_btn.setObjectName("theme_toggle_btn")
        # self.theme_toggle_btn.setFlat(True)
        # self.theme_toggle_btn.setFixedSize(32, 32)
        # self.theme_toggle_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        # self.theme_toggle_btn.setToolTip("Toggle theme (Light/Dark)")
        # self.theme_toggle_btn.clicked.connect(self.toggle_theme)
        
        # Modules label (centered, borderless)
        self.modules_btn = QPushButton("Modules")
        self.modules_btn.setObjectName("modules_btn")
        self.modules_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.modules_btn.setFlat(True)
        self.modules_btn.setToolTip("Modules")
        self.modules_btn.setMinimumWidth(120)  # Settings最小宽degrees以Display完整文本
        
        nav_header.addStretch()
        # nav_header.addWidget(self.theme_toggle_btn)
        nav_header.addWidget(self.modules_btn)
        nav_header.addStretch()
        
        nav_header_widget = QWidget()
        nav_header_widget.setLayout(nav_header)
        # 去除外框，保持简洁
        nav_header_widget.setStyleSheet("""
            background-color: transparent;
            border: none;
            padding: 0px;
        """)
        nav_layout.addWidget(nav_header_widget)

        # 导航按钮列表
        self.nav_list = QListWidget()

        # Update为工作流导向的导航结构
        nav_items =