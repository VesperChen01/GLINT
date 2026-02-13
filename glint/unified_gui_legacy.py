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
    "csv_path": {"zh": "CSV 文件路径", "en": "CSV File Path"},
    "browse": {"zh": "浏览…", "en": "Browse…"},
    "target_obj": {"zh": "目标对象", "en": "Target Object"},
    "refresh": {"zh": "刷新", "en": "Refresh"},
    "btn_highlight": {"zh": "高亮显示", "en": "Highlight"},
    "btn_clear": {"zh": "清空高亮", "en": "Clear"},
    "grp_analysis": {"zh": "相互作用分析", "en": "Interaction Analysis"},
    "pdb_file": {"zh": "PDB 文件（可选）", "en": "PDB File (optional)"},
    "output_csv": {"zh": "输出 CSV（可选）", "en": "Output CSV (optional)"},
    "btn_analyze": {"zh": "开始分析", "en": "Start"},
    "btn_render_interactions": {"zh": "一键渲染（分析+美化+PNG）", "en": "Render (Analyze + Beautify + PNG)"},
    "right_results": {"zh": "结果与日志", "en": "Results & Logs"},
    "table_header": {
        "zh": ["链1", "残基1", "链2", "残基2", "距离", "相互作用"],
        "en": ["Chain1", "Residue1", "Chain2", "Residue2", "Distance", "Interaction"],
    },
    "table_header_gmotif": {
        "zh": ["链", "序列", "起始", "结束", "RMSD (Å)", "类型"],
        "en": ["Chain", "Sequence", "Start", "End", "RMSD (Å)", "Type"],
    },
    "no_object": {"zh": "(无对象)", "en": "(No object)"},
    "select_csv": {"zh": "选择 CSV 文件", "en": "Select CSV"},
    "select_pdb": {"zh": "选择 PDB 文件", "en": "Select PDB"},
    "select_outcsv": {"zh": "选择输出 CSV", "en": "Select output CSV"},
    "log_ready": {"zh": "就绪", "en": "Ready"},
    "log_highlight_ok": {"zh": "高亮完成", "en": "Highlight done"},
    "log_clear_ok": {"zh": "已清空高亮", "en": "Cleared"},
    "log_start": {"zh": "开始相互作用分析…", "en": "Starting interaction analysis…"},
    "log_done": {"zh": "完成，发现 {n} 条记录", "en": "Done: {n} records"},
    "log_error": {"zh": "错误：{msg}", "en": "Error: {msg}"},
    "btn_close": {"zh": "关闭", "en": "Close"},
    "btn_load_csv_to_table": {"zh": "载入到表��", "en": "Load to Table"},
    # G-Motif
    "grp_gmotif": {"zh": "G-Motif（CRBN G-loop）识别", "en": "G-Motif (CRBN G-loop) Detection"},
    "rmsd": {"zh": "RMSD 阈值 (Å)", "en": "RMSD cutoff (Å)"},
    "require_gly": {"zh": "第6位必须为 Gly", "en": "Require Gly at pos6"},
    "btn_gmotif": {"zh": "开始识别", "en": "Detect"},
    "btn_gmotif_render": {"zh": "一键渲染（G-Motif + 电势 + PNG）", "en": "Render (G-Motif + ESP + PNG)"},
    # APBS/Quick + 导出
    "grp_apbs": {"zh": "静电势显示", "en": "Electrostatics Display"},
    "apbs_target": {"zh": "目标对象", "en": "Target Object"},
    "apbs_grid": {"zh": "网格间距 (Å)", "en": "Grid spacing (Å)"},
    "apbs_range": {"zh": "色阶范围 (kT/e)", "en": "Color range (kT/e)"},
    "btn_quick": {"zh": "快速静电图（Coulomb）", "en": "Quick (Coulomb)"},
    "btn_apbs": {"zh": "运行 APBS（若可用）", "en": "Run APBS (if available)"},
    "grp_export": {"zh": "导出", "en": "Export"},
    "img_w": {"zh": "宽度 (px)", "en": "Width (px)"},
    "img_h": {"zh": "高度 (px)", "en": "Height (px)"},
    "img_dpi": {"zh": "DPI", "en": "DPI"},
    "img_bg": {"zh": "背景", "en": "Background"},
    "bg_white": {"zh": "白色", "en": "White"},
    "bg_trans": {"zh": "透明", "en": "Transparent"},
    "raytrace": {"zh": "光线追踪（高质量）", "en": "Ray trace (high quality)"},
    "btn_viewport": {"zh": "使用当前视口尺寸", "en": "Use current viewport"},
    "btn_export_png": {"zh": "导出 PNG", "en": "Export PNG"},
    "btn_export_dx": {"zh": "导出 DX 网格", "en": "Export DX"},
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
            self.progress.emit("[G-Motif] " + ("开始识别…" if get_lang()=="zh" else "Detecting…"))
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
            self.progress.emit("[G-Motif] " + (f"完成，命中 {len(hits)} 条" if get_lang()=="zh" else f"Done, {len(hits)} hits"))
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
                    # 继续深入子布局的子项
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
        self.modules_btn.setMinimumWidth(120)  # 设置最小宽度以显示完整文本
        
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

        # 更新为工作流导向的导航结构
        nav_items = [
            ("Welcome", "Welcome"),
            ("Target Discovery", "Target Discovery"),      # G-Motif, Disease, Pocket
            ("Hit Identification", "Hit Identification"),  # Docking, Virtual Screening
            ("Lead Optimization", "Lead Optimization"),    # Interaction, Ternary, Glue, Mutation
            ("Visualization", "Visualization"),            # Electrostatics, Ray Tracing
        ]

        for zh_text, en_text in nav_items:
            text = zh_text if get_lang() == "zh" else en_text
            self.nav_list.addItem(text)

        self.nav_list.setCurrentRow(0)
        self.nav_list.currentRowChanged.connect(self.on_nav_changed)
        nav_layout.addWidget(self.nav_list, 1)
        
        # Bottom items (Check Environment, README and Contact)
        from PyQt5.QtWidgets import QFrame
        self.nav_separator = QFrame()
        self.nav_separator.setFrameShape(QFrame.Shape.HLine)
        self.nav_separator.setFrameShadow(QFrame.Shadow.Sunken)
        self.nav_separator.setStyleSheet("background-color: #e2e8f0; margin: 8px 0;")
        nav_layout.addWidget(self.nav_separator)
        
        # Check Environment button
        check_env_btn = QPushButton("Check Env")
        check_env_btn.setObjectName("bottom_nav_btn")
        check_env_btn.setFlat(True)
        check_env_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        check_env_btn.clicked.connect(self.check_environment)
        check_env_btn.setToolTip("Check dependencies and environment status")
        nav_layout.addWidget(check_env_btn)
        
        # README button
        readme_btn = QPushButton("README")
        readme_btn.setObjectName("bottom_nav_btn")
        readme_btn.setFlat(True)
        readme_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        readme_btn.clicked.connect(lambda: self.content_stack.setCurrentIndex(5))
        nav_layout.addWidget(readme_btn)
        
        # Contact button
        contact_btn = QPushButton("Contact Us")
        contact_btn.setObjectName("bottom_nav_btn")
        contact_btn.setFlat(True)
        contact_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        contact_btn.clicked.connect(lambda: self.content_stack.setCurrentIndex(6))
        nav_layout.addWidget(contact_btn)

        # log_edit已在__init__中创建，这里不需要再设置
        
        # ========== 内容区域 ==========
        self.content_stack = QStackedWidget()

        # 创建各个页面（log_edit 必须先创建）
        # 0: Welcome
        self.content_stack.addWidget(self.create_welcome_tab())         
        
        # 1: Target Discovery (G-Motif + Disease + Pocket)
        self.content_stack.addWidget(self.create_target_discovery_tab())
        
        # 2: Hit Identification (Docking)
        self.content_stack.addWidget(self.create_hit_identification_tab())
        
        # 3: Lead Optimization (Interaction + Ternary + Glue + Mutation)
        self.content_stack.addWidget(self.create_lead_optimization_tab())
        
        # 4: Visualization (Electrostatics + Export)
        self.content_stack.addWidget(self.create_visualization_tab())
        
        # 5: Check Env (placeholder for nav compatibility, not used in main flow)
        # Note: Button click handlers below use direct index or method call, 
        # but we keep the stack clean.
        # Check Env button calls self.check_environment() directly.
        
        # 5: README
        self.content_stack.addWidget(self.create_readme_tab())          
        
        # 6: Contact
        self.content_stack.addWidget(self.create_contact_tab())

        # ========== Assemble horizontal layout (nav + full content) ==========
        content_row.addWidget(nav_widget, 0)  # Left nav: auto-size to content
        nav_widget.setMaximumWidth(240)
        nav_widget.setMinimumWidth(220)
        content_row.addWidget(self.content_stack, 1)  # Content: takes all remaining space

        main_layout.addLayout(content_row, 1)  # 内容区占全部空间

        # 全局布局调优，统一控件间距与外边距，避免拥挤
        try:
            self._tune_layouts(self)
        except Exception:
            pass

    def on_nav_changed(self, index):
        """导航切换"""
        self.content_stack.setCurrentIndex(index)

    def create_interaction_tab(self) -> QWidget:
        """创建整合的相互作用分析标签页（Protein-Protein + Protein-Ligand + Atom Pairs）"""
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        
        # 创建内容容器
        content_widget = QWidget()
        content_widget.setObjectName("scroll_content")
        self._interaction_scroll_content = content_widget  # 保存引用以便主题切换
        # 强制设置背景色（macOS 兼容性）
        bg_color = "#0d1117" if self._dark_mode else "#f8fafc"
        content_widget.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")
        main_layout = QVBoxLayout(content_widget)
        main_layout.setSpacing(14)
        main_layout.setContentsMargins(12, 12, 12, 12)
        
        # ========== 1. Protein-Protein Interaction ==========
        grp_pp = QGroupBox("Protein-Protein Interaction Analysis")
        pp_layout = QVBoxLayout(grp_pp)
        pp_layout.setSpacing(12)
        pp_layout.setContentsMargins(12, 12, 12, 12)
        
        # Row 1: Target Object
        row1 = QHBoxLayout()
        row1.addWidget(QLabel("Target Object:"), 0)
        self.obj_combo_analysis = QComboBox()
        self.obj_combo_analysis.setMinimumHeight(36)
        self.refresh_obj_analysis = QPushButton(t("refresh"))
        self.refresh_obj_analysis.setObjectName("refresh_btn")
        self.refresh_obj_analysis.setMinimumHeight(36)
        self.refresh_obj_analysis.clicked.connect(self.refresh_objects)
        row1.addWidget(self.obj_combo_analysis, 1)
        row1.addWidget(self.refresh_obj_analysis, 0)
        pp_layout.addLayout(row1)
        
        # Row 2: Output CSV
        row4 = QHBoxLayout()
        row4.addWidget(QLabel("Output CSV (optional):"), 0)
        self.out_csv = QLineEdit()
        self.out_csv.setMinimumHeight(36)
        self.out_browse = QPushButton(t("browse"))
        self.out_browse.setMinimumHeight(36)
        self.out_browse.setObjectName("save_btn")
        self.out_browse.clicked.connect(self.browse_out_csv)
        row4.addWidget(self.out_csv, 1)
        row4.addWidget(self.out_browse, 0)
        pp_layout.addLayout(row4)
        
        btn_pp_row = QHBoxLayout()
        self.analyze_btn = QPushButton("Analyze")
        self.analyze_btn.setObjectName("highlight_btn")
        self.analyze_btn.setMinimumHeight(36)
        self.analyze_btn.clicked.connect(self.start_analysis)
        
        self.pp_heatmap_btn = QPushButton("Heatmap")
        self.pp_heatmap_btn.setObjectName("highlight_btn")
        self.pp_heatmap_btn.setMinimumHeight(36)
        self.pp_heatmap_btn.clicked.connect(self.run_pp_heatmap)
        
        btn_pp_row.addWidget(self.analyze_btn)
        btn_pp_row.addWidget(self.pp_heatmap_btn)
        btn_pp_row.addStretch(1)
        pp_layout.addLayout(btn_pp_row)
        
        main_layout.addWidget(grp_pp)
        
        # ========== 2. Protein-Ligand Interaction ==========
        grp_pl = QGroupBox("Protein-Ligand Interaction Analysis")
        pl_layout = QVBoxLayout(grp_pl)
        pl_layout.setSpacing(12)
        pl_layout.setContentsMargins(12, 12, 12, 12)
        
        # Row 1: Target Object
        pl_row1 = QHBoxLayout()
        pl_row1.addWidget(QLabel("Target Object:"), 0)
        self.pl_obj_combo = QComboBox()
        self.pl_obj_combo.setMinimumHeight(36)
        self.pl_refresh_btn = QPushButton(t("refresh"))
        self.pl_refresh_btn.setObjectName("refresh_btn")
        self.pl_refresh_btn.setMinimumHeight(36)
        self.pl_refresh_btn.clicked.connect(self.refresh_objects)
        pl_row1.addWidget(self.pl_obj_combo, 1)
        pl_row1.addWidget(self.pl_refresh_btn, 0)
        pl_layout.addLayout(pl_row1)
        
        # Row 2: Ligand Resname
        pl_row2 = QHBoxLayout()
        pl_row2.addWidget(QLabel("Ligand Resname:"), 0)
        self.pl_ligand_name = QLineEdit()
        self.pl_ligand_name.setMinimumHeight(36)
        self.pl_ligand_name.setPlaceholderText("Auto-detect if blank")
        pl_row2.addWidget(self.pl_ligand_name, 1)
        pl_layout.addLayout(pl_row2)
        
        # Row 3: Protein Chains
        pl_row3 = QHBoxLayout()
        pl_row3.addWidget(QLabel("Protein Chains:"), 0)
        self.pl_protein_chains = QLineEdit()
        self.pl_protein_chains.setMinimumHeight(36)
        self.pl_protein_chains.setPlaceholderText("Auto-detect if blank")
        pl_row3.addWidget(self.pl_protein_chains, 1)
        pl_layout.addLayout(pl_row3)
        
        # Row 4: Distance cutoff
        pl_row4 = QHBoxLayout()
        pl_row4.addWidget(QLabel("Distance cutoff (Å):"), 0)
        self.pl_distance = QLineEdit("4.5")
        self.pl_distance.setMinimumHeight(36)
        pl_row4.addWidget(self.pl_distance, 1)
        pl_layout.addLayout(pl_row4)
        
        # Row 5: Output CSV
        pl_row5 = QHBoxLayout()
        pl_row5.addWidget(QLabel("Output CSV (optional):"), 0)
        self.pl_csv = QLineEdit()
        self.pl_csv.setMinimumHeight(36)
        self.pl_csv.setPlaceholderText("Optional")
        self.pl_csv_btn = QPushButton(t("browse"))
        self.pl_csv_btn.setMinimumHeight(36)
        self.pl_csv_btn.setObjectName("browse_btn")
        self.pl_csv_btn.clicked.connect(lambda: self._browse_save_file(self.pl_csv, "CSV (*.csv)"))
        pl_row5.addWidget(self.pl_csv, 1)
        pl_row5.addWidget(self.pl_csv_btn, 0)
        pl_layout.addLayout(pl_row5)
        
        # Row 6: Plot Style
        pl_row6 = QHBoxLayout()
        pl_row6.addWidget(QLabel("🎨 Plot Style:"), 0)
        self.pl_plot_style = QComboBox()
        self.pl_plot_style.setMinimumHeight(36)
        self.pl_plot_style.addItems([
            "Professional (Default)",
            "Hand-drawn (xkcd)",
            "Minimalist",
            "Publication",
            "Colorful"
        ])
        self.pl_plot_style.setToolTip(
            "Choose the visual style for network plots:\n"
            "- Professional: Clean and modern\n"
            "- Hand-drawn: Casual xkcd style\n"
            "- Minimalist: Simple and elegant\n"
            "- Publication: High-quality for papers\n"
            "- Colorful: Vibrant and eye-catching"
        )
        pl_row6.addWidget(self.pl_plot_style, 1)
        pl_layout.addLayout(pl_row6)
        
        # Row 7: Visualization Options
        pl_row7 = QHBoxLayout()
        pl_row7.addWidget(QLabel("📐 Visualization:"), 0)
        
        # Show hydrophobic checkbox
        self.pl_show_hydrophobic = QCheckBox("Show hydrophobic")
        self.pl_show_hydrophobic.setToolTip(
            "Include hydrophobic interactions in visualization.\n"
            "Note: Even with high confidence, hydrophobic interactions\n"
            "are hidden by default (professional mode)."
        )
        pl_row7.addWidget(self.pl_show_hydrophobic, 0)
        
        # Min confidence filter
        pl_row7.addWidget(QLabel("Min Confidence:"), 0)
        self.pl_min_confidence = QComboBox()
        self.pl_min_confidence.setMinimumHeight(36)
        self.pl_min_confidence.addItems(["0.0", "0.3", "0.5", "0.7", "0.8 (default)", "0.9", "1.0"])
        self.pl_min_confidence.setCurrentIndex(4)  # 0.8 by default
        self.pl_min_confidence.setToolTip(
            "Minimum confidence score for interactions to display.\n"
            "Lower values show more interactions; higher values\n"
            "show only high-confidence interactions."
        )
        self.pl_min_confidence.setFixedWidth(120)
        pl_row7.addWidget(self.pl_min_confidence, 0)
        
        pl_row7.addStretch(1)
        pl_layout.addLayout(pl_row7)
        
        btn_pl_row = QHBoxLayout()
        self.pl_analyze_btn = QPushButton("Analyze")
        self.pl_analyze_btn.setObjectName("highlight_btn")
        self.pl_analyze_btn.setMinimumHeight(36)
        self.pl_analyze_btn.clicked.connect(self.run_pl_analysis)
        
        self.pl_visualize_btn = QPushButton("3D Visualize")
        self.pl_visualize_btn.setObjectName("highlight_btn")
        self.pl_visualize_btn.setMinimumHeight(36)
        self.pl_visualize_btn.clicked.connect(self.run_pl_visualize)
        
        self.pl_network_btn = QPushButton("Network Plot")
        self.pl_network_btn.setObjectName("highlight_btn")
        self.pl_network_btn.setMinimumHeight(36)
        self.pl_network_btn.clicked.connect(self.run_pl_network)
        
        btn_pl_row.addWidget(self.pl_analyze_btn)
        btn_pl_row.addWidget(self.pl_visualize_btn)
        btn_pl_row.addWidget(self.pl_network_btn)
        btn_pl_row.addStretch(1)
        pl_layout.addLayout(btn_pl_row)
        
        main_layout.addWidget(grp_pl)
        
        # ========== 3. Ligand-Ligand Interaction (New Feature) ==========
        grp_ll = QGroupBox("Ligand-Ligand Interaction Analysis (Small Molecule - Small Molecule)")
        ll_layout = QVBoxLayout(grp_ll)
        ll_layout.setSpacing(12)
        ll_layout.setContentsMargins(12, 12, 12, 12)
        
        # Row 1: Target Object & PDB File
        ll_row1 = QHBoxLayout()
        
        # Left: Object Combo
        ll_row1.addWidget(QLabel("Target Object:"), 0)
        self.ll_obj_combo = QComboBox()
        self.ll_obj_combo.setMinimumHeight(36)
        self.ll_refresh_btn = QPushButton(t("refresh"))
        self.ll_refresh_btn.setObjectName("refresh_btn")
        self.ll_refresh_btn.setMinimumHeight(36)
        self.ll_refresh_btn.clicked.connect(self.refresh_objects)
        
        obj_container = QWidget()
        obj_layout = QHBoxLayout(obj_container)
        obj_layout.setContentsMargins(0,0,0,0)
        obj_layout.addWidget(self.ll_obj_combo, 1)
        obj_layout.addWidget(self.ll_refresh_btn, 0)
        
        ll_row1.addWidget(obj_container, 1)
        
        # Right: PDB File (Optional)
        ll_row1.addWidget(QLabel("PDB File (optional):"), 0)
        self.ll_pdb = QLineEdit()
        self.ll_pdb.setPlaceholderText("Load from file...")
        self.ll_pdb.setMinimumHeight(36)
        self.ll_pdb_browse = QPushButton(t("browse"))
        self.ll_pdb_browse.setObjectName("browse_btn")
        self.ll_pdb_browse.setMinimumHeight(36)
        self.ll_pdb_browse.clicked.connect(lambda: self._browse_file(self.ll_pdb, "PDB Files (*.pdb *.cif *.sdf)"))
        
        pdb_container = QWidget()
        pdb_layout = QHBoxLayout(pdb_container)
        pdb_layout.setContentsMargins(0,0,0,0)
        pdb_layout.addWidget(self.ll_pdb, 1)
        pdb_layout.addWidget(self.ll_pdb_browse, 0)
        
        ll_row1.addWidget(pdb_container, 1)
        
        ll_layout.addLayout(ll_row1)
        
        # Row 2: Selection 1 & 2
        ll_row2 = QHBoxLayout()
        
        ll_row2.addWidget(QLabel("Selection 1:"), 0)
        self.ll_sel1 = QLineEdit()
        self.ll_sel1.setPlaceholderText("e.g. resn LIG1 or resi 100")
        self.ll_sel1.setMinimumHeight(36)
        ll_row2.addWidget(self.ll_sel1, 1)
        
        ll_row2.addWidget(QLabel("Selection 2:"), 0)
        self.ll_sel2 = QLineEdit()
        self.ll_sel2.setPlaceholderText("e.g. resn LIG2 or resi 200")
        self.ll_sel2.setMinimumHeight(36)
        ll_row2.addWidget(self.ll_sel2, 1)
        
        ll_layout.addLayout(ll_row2)
        
        # Row 3: Distance & CSV
        ll_row3 = QHBoxLayout()
        
        ll_row3.addWidget(QLabel("Distance (Å):"), 0)
        self.ll_dist = QLineEdit("4.5")
        self.ll_dist.setFixedWidth(60)
        self.ll_dist.setMinimumHeight(36)
        ll_row3.addWidget(self.ll_dist, 0)
        
        ll_row3.addWidget(QLabel("Output CSV:"), 0)
        self.ll_csv = QLineEdit()
        self.ll_csv.setPlaceholderText("Optional")
        self.ll_csv.setMinimumHeight(36)
        ll_row3.addWidget(self.ll_csv, 1)
        
        self.ll_csv_btn = QPushButton(t("browse"))
        self.ll_csv_btn.setMinimumHeight(36)
        self.ll_csv_btn.setObjectName("browse_btn")
        self.ll_csv_btn.clicked.connect(lambda: self._browse_save_file(self.ll_csv, "CSV (*.csv)"))
        ll_row3.addWidget(self.ll_csv_btn, 0)
        
        ll_layout.addLayout(ll_row3)
        
        # Row 4: Analyze Button
        ll_row4 = QHBoxLayout()
        self.ll_analyze_btn = QPushButton("Analyze Ligand-Ligand Interactions")
        self.ll_analyze_btn.setObjectName("highlight_btn")
        self.ll_analyze_btn.setMinimumHeight(36)
        self.ll_analyze_btn.clicked.connect(self.run_ll_analysis)
        ll_row4.addWidget(self.ll_analyze_btn)
        ll_row4.addStretch(1)
        ll_layout.addLayout(ll_row4)
        
        main_layout.addWidget(grp_ll)
        
        # ========== 4. Protein-Nucleic Acid Interaction ==========
        grp_pn = QGroupBox("Protein-Nucleic Acid Interaction Analysis")
        pn_layout = QVBoxLayout(grp_pn)
        pn_layout.setSpacing(12)
        pn_layout.setContentsMargins(12, 12, 12, 12)
        
        # Row 1: Target Object
        pn_row1 = QHBoxLayout()
        pn_row1.addWidget(QLabel("Target Object:"), 0)
        self.pn_obj_combo = QComboBox()
        self.pn_obj_combo.setMinimumHeight(36)
        self.pn_refresh_btn = QPushButton(t("refresh"))
        self.pn_refresh_btn.setObjectName("refresh_btn")
        self.pn_refresh_btn.setMinimumHeight(36)
        self.pn_refresh_btn.clicked.connect(self.refresh_objects)
        pn_row1.addWidget(self.pn_obj_combo, 1)
        pn_row1.addWidget(self.pn_refresh_btn, 0)
        pn_layout.addLayout(pn_row1)
        
        # Row 1.5: PDB File (Optional)
        pn_row_pdb = QHBoxLayout()
        pn_row_pdb.addWidget(QLabel("PDB File (optional):"), 0)
        self.pn_pdb = QLineEdit()
        self.pn_pdb.setPlaceholderText("Load from file...")
        self.pn_pdb.setMinimumHeight(36)
        self.pn_pdb_browse = QPushButton(t("browse"))
        self.pn_pdb_browse.setObjectName("browse_btn")
        self.pn_pdb_browse.setMinimumHeight(36)
        self.pn_pdb_browse.clicked.connect(lambda: self._browse_file(self.pn_pdb, "PDB Files (*.pdb *.cif *.sdf)"))
        pn_row_pdb.addWidget(self.pn_pdb, 1)
        pn_row_pdb.addWidget(self.pn_pdb_browse, 0)
        pn_layout.addLayout(pn_row_pdb)
        
        # Row 2: Nucleic Chains
        pn_row2 = QHBoxLayout()
        pn_row2.addWidget(QLabel("Nucleic Chains:"), 0)
        self.pn_nucleic_chains = QLineEdit()
        self.pn_nucleic_chains.setMinimumHeight(36)
        self.pn_nucleic_chains.setPlaceholderText("Auto-detect if blank")
        pn_row2.addWidget(self.pn_nucleic_chains, 1)
        pn_layout.addLayout(pn_row2)
        
        # Row 3: Protein Chains
        pn_row3 = QHBoxLayout()
        pn_row3.addWidget(QLabel("Protein Chains:"), 0)
        self.pn_protein_chains = QLineEdit()
        self.pn_protein_chains.setMinimumHeight(36)
        self.pn_protein_chains.setPlaceholderText("Auto-detect if blank")
        pn_row3.addWidget(self.pn_protein_chains, 1)
        pn_layout.addLayout(pn_row3)
        
        # Row 4: Distance & CSV
        pn_row4 = QHBoxLayout()
        pn_row4.addWidget(QLabel("Distance (Å):"), 0)
        self.pn_distance = QLineEdit("4.5")
        self.pn_distance.setFixedWidth(60)
        self.pn_distance.setMinimumHeight(36)
        pn_row4.addWidget(self.pn_distance, 0)
        
        pn_row4.addWidget(QLabel("Output CSV:"), 0)
        self.pn_csv = QLineEdit()
        self.pn_csv.setPlaceholderText("Optional")
        self.pn_csv.setMinimumHeight(36)
        pn_row4.addWidget(self.pn_csv, 1)
        
        self.pn_csv_btn = QPushButton(t("browse"))
        self.pn_csv_btn.setMinimumHeight(36)
        self.pn_csv_btn.setObjectName("browse_btn")
        self.pn_csv_btn.clicked.connect(lambda: self._browse_save_file(self.pn_csv, "CSV (*.csv)"))
        pn_row4.addWidget(self.pn_csv_btn, 0)
        pn_layout.addLayout(pn_row4)
        
        # Row 5: Analyze Button
        pn_row5 = QHBoxLayout()
        self.pn_analyze_btn = QPushButton("Analyze Protein-Nucleic Acid")
        self.pn_analyze_btn.setObjectName("highlight_btn")
        self.pn_analyze_btn.setMinimumHeight(36)
        self.pn_analyze_btn.clicked.connect(self.run_pn_analysis)
        
        self.pn_heatmap_btn = QPushButton("Heatmap")
        self.pn_heatmap_btn.setObjectName("highlight_btn")
        self.pn_heatmap_btn.setMinimumHeight(36)
        self.pn_heatmap_btn.clicked.connect(self.run_pn_heatmap)
        
        pn_row5.addWidget(self.pn_analyze_btn)
        pn_row5.addWidget(self.pn_heatmap_btn)
        pn_row5.addStretch(1)
        pn_layout.addLayout(pn_row5)
        
        main_layout.addWidget(grp_pn)
        
        # ========== 3. Atom Pair Analysis ========== (DISABLED - 高级功能,一般用户不需要)
        # grp_ap = QGroupBox("Atom Pair Analysis (Atomic-Level Precision)")
        # ap_layout = QVBoxLayout(grp_ap)
        # ap_layout.setSpacing(12)
        # ap_layout.setContentsMargins(12, 12, 12, 12)
        
        # # Row 1: Target Object
        # ap_row1 = QHBoxLayout()
        # ap_row1.addWidget(QLabel("Target Object:"), 0)
        # self.ap_obj_combo = QComboBox()
        # self.ap_obj_combo.setMinimumHeight(36)
        # self.ap_refresh_btn = QPushButton(t("refresh"))
        # self.ap_refresh_btn.setObjectName("refresh_btn")
        # self.ap_refresh_btn.setMinimumHeight(36)
        # self.ap_refresh_btn.clicked.connect(self.refresh_objects)
        # ap_row1.addWidget(self.ap_obj_combo, 1)
        # ap_row1.addWidget(self.ap_refresh_btn, 0)
        # ap_layout.addLayout(ap_row1)
        # 
        # # Row 2: Atom1 Selection
        # ap_row2 = QHBoxLayout()
        # ap_row2.addWidget(QLabel("Atom1 Selection:"), 0)
        # self.ap_atom1 = QLineEdit()
        # self.ap_atom1.setMinimumHeight(36)
        # self.ap_atom1.setPlaceholderText('e.g.: "resn LIG and name N1"')
        # ap_row2.addWidget(self.ap_atom1, 1)
        # ap_layout.addLayout(ap_row2)
        # 
        # # Row 3: Atom2 Selection
        # ap_row3 = QHBoxLayout()
        # ap_row3.addWidget(QLabel("Atom2 Selection:"), 0)
        # self.ap_atom2 = QLineEdit()
        # self.ap_atom2.setMinimumHeight(36)
        # self.ap_atom2.setPlaceholderText('e.g.: "elem O"')
        # ap_row3.addWidget(self.ap_atom2, 1)
        # ap_layout.addLayout(ap_row3)
        # 
        # # Row 4: Distance cutoff
        # ap_row4 = QHBoxLayout()
        # ap_row4.addWidget(QLabel("Distance cutoff (Å):"), 0)
        # self.ap_distance = QLineEdit("5.0")
        # self.ap_distance.setMinimumHeight(36)
        # ap_row4.addWidget(self.ap_distance, 1)
        # ap_layout.addLayout(ap_row4)
        # 
        # # Row 5: Output CSV
        # ap_row5 = QHBoxLayout()
        # ap_row5.addWidget(QLabel("Output CSV (optional):"), 0)
        # self.ap_csv = QLineEdit()
        # self.ap_csv.setMinimumHeight(36)
        # self.ap_csv.setPlaceholderText("Optional")
        # self.ap_csv_btn = QPushButton(t("browse"))
        # self.ap_csv_btn.setMinimumHeight(36)
        # self.ap_csv_btn.setObjectName("browse_btn")
        # self.ap_csv_btn.clicked.connect(lambda: self._browse_save_file(self.ap_csv, "CSV (*.csv)"))
        # ap_row5.addWidget(self.ap_csv, 1)
        # ap_row5.addWidget(self.ap_csv_btn, 0)
        # ap_layout.addLayout(ap_row5)
        # 
        # main_layout.addWidget(grp_ap)
        # 
        # # Quick templates - 确保所有按钮在一行显示
        # main_layout.addWidget(QLabel("Quick Templates:"))
        # template_row = QHBoxLayout()
        # templates = [
        #     ("N-O H-bonds", '"elem N"', '"elem O"', "3.5"),
        #     ("Lig-SER", '"resn LIG"', '"resn SER"', "4.5"),
        #     ("S-S", '"name SG"', '"name SG"', "2.5"),
        #     ("π-Stacking", '"aromatic"', '"aromatic"', "4.5"),
        #     ("π-Cation", '"aromatic"', '"basic"', "4.0"),
        # ]
        # for label, atom1, atom2, dist in templates:
        #     btn = QPushButton(label)
        #     btn.setObjectName("browse_btn")
        #     btn.setMinimumHeight(36)
        #     btn.clicked.connect(lambda checked, a1=atom1, a2=atom2, d=dist: self.apply_ap_template(a1, a2, d))
        #     template_row.addWidget(btn)
        # template_row.addStretch(1)
        # main_layout.addLayout(template_row)
        # 
        # btn_ap_row = QHBoxLayout()
        # self.ap_analyze_btn = QPushButton("Analyze Pairs")
        # self.ap_analyze_btn.setObjectName("highlight_btn")
        # self.ap_analyze_btn.setMinimumHeight(36)
        # self.ap_analyze_btn.clicked.connect(self.run_ap_analysis)
        # 
        # self.ap_visualize_btn = QPushButton("Visualize")
        # self.ap_visualize_btn.setObjectName("highlight_btn")
        # self.ap_visualize_btn.setMinimumHeight(36)
        # self.ap_visualize_btn.clicked.connect(self.run_ap_visualize)
        # 
        # btn_ap_row.addWidget(self.ap_analyze_btn)
        # btn_ap_row.addWidget(self.ap_visualize_btn)
        # btn_ap_row.addStretch(1)
        # main_layout.addLayout(btn_ap_row)
        
        main_layout.addStretch(1)
        
        # 将内容设置到滚动区域
        scroll_area.setWidget(content_widget)
        
        # 返回包裹的滚动区域
        wrapper = QWidget()
        wrapper_layout = QVBoxLayout(wrapper)
        wrapper_layout.setContentsMargins(0, 0, 0, 0)
        wrapper_layout.addWidget(scroll_area)
        return wrapper
    
    def create_analysis_tab(self) -> QWidget:
        w = QWidget(); lay = QVBoxLayout(w); lay.setSpacing(10)

        # ===== 相互作用分析组 =====
        grp = QGroupBox(t("grp_analysis")); form = QFormLayout(grp); form.setLabelAlignment(Qt.AlignmentFlag.AlignRight); form.setSpacing(12)
        row0 = QHBoxLayout()
        self.obj_combo_analysis = QComboBox()
        self.refresh_obj_analysis = QPushButton(t("refresh")); self.refresh_obj_analysis.setObjectName("refresh_btn"); self.refresh_obj_analysis.clicked.connect(self.refresh_objects)
        row0.addWidget(self.obj_combo_analysis, 1); row0.addWidget(self.refresh_obj_analysis)
        form.addRow(QLabel(t("target_obj")), row0)
        row1 = QHBoxLayout()
        self.pdb_path = QLineEdit()
        self.pdb_browse = QPushButton(t("browse")); self.pdb_browse.setObjectName("browse_btn"); self.pdb_browse.clicked.connect(self.browse_pdb)
        row1.addWidget(self.pdb_path, 1); row1.addWidget(self.pdb_browse)
        form.addRow(QLabel(t("pdb_file")), row1)
        # Removed checkbox - default to only_between_chains=True
        row2 = QHBoxLayout()
        self.out_csv = QLineEdit()
        self.out_browse = QPushButton(t("browse")); self.out_browse.setObjectName("save_btn"); self.out_browse.clicked.connect(self.browse_out_csv)
        row2.addWidget(self.out_csv, 1); row2.addWidget(self.out_browse)
        form.addRow(QLabel(t("output_csv")), row2)

        # 按钮行
        btn_row = QHBoxLayout()
        self.analyze_btn = QPushButton(t("btn_analyze"))
        self.analyze_btn.setObjectName("highlight_btn")
        self.analyze_btn.clicked.connect(self.start_analysis)
        self.render_interact_btn = QPushButton(t("btn_render_interactions"))
        self.render_interact_btn.setObjectName("highlight_btn")
        self.render_interact_btn.clicked.connect(self.render_interactions_beautifully_clicked)
        btn_row.addWidget(self.analyze_btn)
        btn_row.addWidget(self.render_interact_btn)
        btn_row.addStretch(1)

        # ===== CSV 高亮组 =====
        grp_csv = QGroupBox(t("grp_csv")); form_csv = QFormLayout(grp_csv); form_csv.setLabelAlignment(Qt.AlignmentFlag.AlignRight); form_csv.setSpacing(12)
        row_csv = QHBoxLayout()
        self.csv_path = QLineEdit()
        self.csv_browse = QPushButton(t("browse")); self.csv_browse.setObjectName("browse_btn"); self.csv_browse.clicked.connect(self.browse_csv)
        self.csv_load_to_table = QPushButton(t("btn_load_csv_to_table")); self.csv_load_to_table.setObjectName("save_btn")
        self.csv_load_to_table.clicked.connect(self.load_csv_to_table)
        row_csv.addWidget(self.csv_path, 1); row_csv.addWidget(self.csv_browse); row_csv.addWidget(self.csv_load_to_table)
        form_csv.addRow(QLabel(t("csv_path")), row_csv)
        row2_csv = QHBoxLayout()
        self.obj_combo_csv = QComboBox()
        self.refresh_obj_csv = QPushButton(t("refresh")); self.refresh_obj_csv.setObjectName("refresh_btn"); self.refresh_obj_csv.clicked.connect(self.refresh_objects)
        row2_csv.addWidget(self.obj_combo_csv, 1); row2_csv.addWidget(self.refresh_obj_csv)
        form_csv.addRow(QLabel(t("target_obj")), row2_csv)

        lay.addWidget(grp)
        lay.addLayout(btn_row)
        lay.addWidget(grp_csv)
        lay.addStretch(1)
        return w
    
    def create_welcome_tab(self) -> QWidget:
        """Create welcome page with molecular glue introduction"""
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setSpacing(14)
        layout.setContentsMargins(20, 16, 20, 16)
        
        # Welcome title
        title_label = QLabel("GLINT - Molecular Glue Analyzer")
        title_label.setStyleSheet("""
            QLabel {
                font-size: 36px;
                font-weight: bold;
                color: #1e40af;
                padding: 10px;
            }
        """)
        title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(title_label)
        
        # Subtitle
        subtitle = QLabel("POI Discovery & Targeted Protein Degradation Analysis")
        subtitle.setStyleSheet("""
            QLabel {
                font-size: 18px;
                color: #475569;
                padding: 5px;
            }
        """)
        subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.addWidget(subtitle)
        
        
        layout.addStretch()
        return w

    def create_molecular_glue_tab(self) -> QWidget:
        """创建分子胶设计模块 - POI Discovery"""
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        
        # 创建内容容器
        content_widget = QWidget()
        content_widget.setObjectName("scroll_content")
        self._molecular_glue_scroll_content = content_widget  # 保存引用以便主题切换
        # 强制设置背景色（macOS 兼容性）
        bg_color = "#0d1117" if self._dark_mode else "#f8fafc"
        content_widget.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")
        main_layout = QVBoxLayout(content_widget)
        main_layout.setSpacing(12) # 恢复间距
        main_layout.setContentsMargins(10, 10, 10, 10) # 恢复边距
        
        # ========== POI Discovery - G-Motif Detection ==========
        grp_gm = QGroupBox("POI Discovery - G-Motif (CRBN G-loop) Detection")
        # 使用网格布局代替FormLayout，实现两栏式
        gm_grid = QGridLayout(grp_gm)
        gm_grid.setColumnStretch(0, 0)
        gm_grid.setColumnStretch(1, 1)
        gm_grid.setColumnStretch(2, 0)
        gm_grid.setColumnStretch(3, 1)
        gm_grid.setHorizontalSpacing(8)
        gm_grid.setVerticalSpacing(10)  # 减少垂直间距以节省空间
        
        # 第一行： Target Object | [combo + refresh] | PDB File | [input + browse]
        gm_grid.addWidget(QLabel("Target Object"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.obj_combo_gm = QComboBox()
        self.obj_combo_gm.setMinimumHeight(32)
        self.refresh_obj_gm = QPushButton(t("refresh"))
        self.refresh_obj_gm.setObjectName("refresh_btn")
        self.refresh_obj_gm.setMinimumHeight(32)
        self.refresh_obj_gm.clicked.connect(self.refresh_objects)
        obj_row_container = QWidget()
        obj_row = QHBoxLayout(obj_row_container)
        obj_row.setContentsMargins(0, 0, 0, 0)
        obj_row.addWidget(self.obj_combo_gm, 1)
        obj_row.addWidget(self.refresh_obj_gm)
        gm_grid.addWidget(obj_row_container, 0, 1)
        
        gm_grid.addWidget(QLabel("PDB File (optional)"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.gm_pdb = QLineEdit()
        self.gm_pdb.setMinimumHeight(32)
        self.gm_pdb_browse = QPushButton(t("browse"))
        self.gm_pdb_browse.setMinimumHeight(32)
        self.gm_pdb_browse.setObjectName("browse_btn")
        self.gm_pdb_browse.clicked.connect(self.browse_gm_pdb)
        pdb_row_container = QWidget()
        pdb_row = QHBoxLayout(pdb_row_container)
        pdb_row.setContentsMargins(0, 0, 0, 0)
        pdb_row.addWidget(self.gm_pdb, 1)
        pdb_row.addWidget(self.gm_pdb_browse)
        gm_grid.addWidget(pdb_row_container, 0, 3)
        
        # 第二行： RMSD cutoff | [input] | Require Gly | [checkbox]
        gm_grid.addWidget(QLabel("RMSD cutoff (Å)"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.gm_rmsd = QLineEdit("3.5")
        self.gm_rmsd.setMinimumHeight(32)
        gm_grid.addWidget(self.gm_rmsd, 1, 1)
        
        gm_grid.addWidget(QLabel(""), 1, 2)  # 空位
        self.gm_require_gly = QCheckBox(t("require_gly"))
        self.gm_require_gly.setChecked(True)
        gm_grid.addWidget(self.gm_require_gly, 1, 3)
        
        # 第三行： Template Source | [combo] | Template Selection | [input + button]
        gm_grid.addWidget(QLabel("Template Source"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.gm_template_mode = QComboBox()
        self.gm_template_mode.setMinimumHeight(32)
        self.gm_template_mode.addItems([
            "Idealized (8×Cα)",
            "Built-in: GSPT1 (6H0G A:60-67)",
            "Built-in: CK1α (3M51 A:36-43)",
            "Built-in: VAV1 (2MC1 A:95-102)",
            "From Selection (8×Cα)",
        ])
        gm_grid.addWidget(self.gm_template_mode, 2, 1)
        
        gm_grid.addWidget(QLabel("Template Selection"), 2, 2, Qt.AlignmentFlag.AlignRight)
        self.gm_template_sel = QLineEdit()
        self.gm_template_sel.setMinimumHeight(32)
        self.gm_template_sel.setPlaceholderText("Enter selection, e.g., sele")
        self.gm_template_pick = QPushButton("Get Current (sele)")
        self.gm_template_pick.setMinimumHeight(32)
        self.gm_template_pick.setObjectName("refresh_btn")
        self.gm_template_pick.clicked.connect(lambda: self.gm_template_sel.setText("sele"))
        temp_row_container = QWidget()
        temp_row = QHBoxLayout(temp_row_container)
        temp_row.setContentsMargins(0, 0, 0, 0)
        temp_row.addWidget(self.gm_template_sel, 1)
        temp_row.addWidget(self.gm_template_pick)
        gm_grid.addWidget(temp_row_container, 2, 3)
        
        def _toggle_template_inputs(idx):
            use_sel = (idx == 4)
            self.gm_template_sel.setEnabled(use_sel)
            self.gm_template_pick.setEnabled(use_sel)
        self.gm_template_mode.currentIndexChanged.connect(_toggle_template_inputs)
        _toggle_template_inputs(self.gm_template_mode.currentIndex())
        
        # 第四行： Output CSV | [input + browse] (跨两列)
        gm_grid.addWidget(QLabel("Output CSV (optional)"), 3, 0, Qt.AlignmentFlag.AlignRight)
        self.gm_out_csv = QLineEdit()
        self.gm_out_csv.setMinimumHeight(36)
        self.gm_out_csv.setPlaceholderText("Optional - leave blank for temp CSV")
        self.gm_out_browse = QPushButton(t("browse"))
        self.gm_out_browse.setMinimumHeight(36)
        self.gm_out_browse.setObjectName("save_btn")
        self.gm_out_browse.clicked.connect(self.browse_gm_out_csv)
        out_row_container = QWidget()
        out_row = QHBoxLayout(out_row_container)
        out_row.setContentsMargins(0, 0, 0, 0)
        out_row.addWidget(self.gm_out_csv, 1)
        out_row.addWidget(self.gm_out_browse)
        gm_grid.addWidget(out_row_container, 3, 1, 1, 3)  # 跨三列
        
        btn_gm_row = QHBoxLayout()
        self.gm_btn = QPushButton("Detect POI")
        self.gm_btn.setObjectName("highlight_btn")
        self.gm_btn.clicked.connect(self.start_gmotif)
        self.gm_btn_render = QPushButton("Render All (POI + ESP + PNG)")
        self.gm_btn_render.setObjectName("highlight_btn")
        self.gm_btn_render.clicked.connect(self.render_gmotif_with_esp)
        btn_gm_row.addWidget(self.gm_btn)
        btn_gm_row.addWidget(self.gm_btn_render)
        btn_gm_row.addStretch(1)
        
        main_layout.addWidget(grp_gm)
        main_layout.addLayout(btn_gm_row)
        
        # ========== Ternary Complex Analysis (with integrated interface tools) ==========
        grp_ternary = QGroupBox("Ternary Complex Analysis (E3-PROTAC-POI with Interface Tools)")
        # 使用网格布局实现两栏式
        ternary_grid = QGridLayout(grp_ternary)
        ternary_grid.setColumnStretch(0, 0)
        ternary_grid.setColumnStretch(1, 1)
        ternary_grid.setColumnStretch(2, 0)
        ternary_grid.setColumnStretch(3, 1)
        ternary_grid.setHorizontalSpacing(8)
        ternary_grid.setVerticalSpacing(10)  # 减少垂直间距以节省空间
        
        # 第一行: Target Object | [combo + refresh] | Ligand Resname | [input]
        ternary_grid.addWidget(QLabel("Target Object"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.tc_obj_combo = QComboBox()
        self.tc_obj_combo.setMinimumHeight(36)
        self.tc_refresh_btn = QPushButton(t("refresh"))
        self.tc_refresh_btn.setObjectName("refresh_btn")
        self.tc_refresh_btn.setMinimumHeight(36)
        self.tc_refresh_btn.clicked.connect(self.refresh_objects)
        tc_obj_row_container = QWidget()
        tc_obj_row = QHBoxLayout(tc_obj_row_container)
        tc_obj_row.setContentsMargins(0, 0, 0, 0)
        tc_obj_row.addWidget(self.tc_obj_combo, 1)
        tc_obj_row.addWidget(self.tc_refresh_btn)
        ternary_grid.addWidget(tc_obj_row_container, 0, 1)
        
        ternary_grid.addWidget(QLabel("Ligand Resname:"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.tc_ligand_name = QLineEdit()
        self.tc_ligand_name.setMinimumHeight(36)
        self.tc_ligand_name.setPlaceholderText("PROTAC/Glue name, auto-detect if blank")
        ternary_grid.addWidget(self.tc_ligand_name, 0, 3)
        
        # 第二行: E3 Chains (CRBN/VHL) | [input] | POI Chains | [input]
        ternary_grid.addWidget(QLabel("E3 Ligase Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.tc_protein1_chains = QLineEdit()
        self.tc_protein1_chains.setMinimumHeight(36)
        self.tc_protein1_chains.setPlaceholderText("e.g.: A (CRBN/VHL/IAP)")
        ternary_grid.addWidget(self.tc_protein1_chains, 1, 1)
        
        ternary_grid.addWidget(QLabel("POI Chains:"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.tc_protein2_chains = QLineEdit()
        self.tc_protein2_chains.setMinimumHeight(36)
        self.tc_protein2_chains.setPlaceholderText("e.g.: B (Target Protein)")
        ternary_grid.addWidget(self.tc_protein2_chains, 1, 3)
        
        # 第三行: Interface cutoff | [input] | Output CSV | [input + browse]
        ternary_grid.addWidget(QLabel("Interface cutoff (Å):"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.tc_distance = QLineEdit("4.5")
        self.tc_distance.setMinimumHeight(36)
        ternary_grid.addWidget(self.tc_distance, 2, 1)
        
        ternary_grid.addWidget(QLabel("Output CSV (optional)"), 2, 2, Qt.AlignmentFlag.AlignRight)
        self.tc_csv = QLineEdit()
        self.tc_csv.setMinimumHeight(36)
        self.tc_csv.setPlaceholderText("Optional")
        self.tc_csv_btn = QPushButton(t("browse"))
        self.tc_csv_btn.setMinimumHeight(36)
        self.tc_csv_btn.setObjectName("browse_btn")
        self.tc_csv_btn.clicked.connect(lambda: self._browse_save_file(self.tc_csv, "CSV (*.csv)"))
        tc_csv_row = QHBoxLayout()
        tc_csv_row.addWidget(self.tc_csv, 1)
        tc_csv_row.addWidget(self.tc_csv_btn)
        ternary_grid.addLayout(tc_csv_row, 2, 3)
        
        # 第四行: Analyze P-P Interface | [checkbox] | Include ΔΔG | [checkbox]
        ternary_grid.addWidget(QLabel(""), 3, 0)  # 空位
        self.tc_analyze_interface = QCheckBox("Analyze Protein-Protein Interface")
        self.tc_analyze_interface.setChecked(True)
        ternary_grid.addWidget(self.tc_analyze_interface, 3, 1, 1, 2)
        
        self.tc_include_ddg = QCheckBox("Calculate ΔΔG")
        ternary_grid.addWidget(self.tc_include_ddg, 3, 3)
        
        # Ternary buttons (expanded functionality)
        tc_btn_row = QHBoxLayout()
        self.tc_analyze_btn = QPushButton("Analyze Complex")
        self.tc_analyze_btn.setObjectName("highlight_btn")
        self.tc_analyze_btn.clicked.connect(self.run_tc_analysis)
        
        self.tc_interface_btn = QPushButton("Interface Map")
        self.tc_interface_btn.setObjectName("highlight_btn")
        self.tc_interface_btn.clicked.connect(self.run_tc_interface)
        
        self.tc_network_btn = QPushButton("Network Plot")
        self.tc_network_btn.setObjectName("highlight_btn")
        self.tc_network_btn.clicked.connect(self.run_tc_network)
        
        self.tc_render_btn = QPushButton("Render All")
        self.tc_render_btn.setObjectName("highlight_btn") 
        self.tc_render_btn.clicked.connect(self.run_tc_render)
        
        tc_btn_row.addWidget(self.tc_analyze_btn)
        tc_btn_row.addWidget(self.tc_interface_btn)
        tc_btn_row.addWidget(self.tc_network_btn)
        tc_btn_row.addWidget(self.tc_render_btn)
        tc_btn_row.addStretch(1)
        
        main_layout.addWidget(grp_ternary)
        main_layout.addLayout(tc_btn_row)
        
        # ========== 新增: Molecular Glue-Specific Analysis ==========
        grp_glue = QGroupBox("Molecular Glue Analysis (PPI + Neo-Epitope Detection)")
        glue_grid = QGridLayout(grp_glue)
        glue_grid.setColumnStretch(0, 0)
        glue_grid.setColumnStretch(1, 1)
        glue_grid.setColumnStretch(2, 0)
        glue_grid.setColumnStretch(3, 1)
        glue_grid.setHorizontalSpacing(8)
        glue_grid.setVerticalSpacing(10)  # 减少垂直间距以节省空间
        
        # 第一行: Target Object | [combo + refresh] | Glue Residue Name | [input]
        glue_grid.addWidget(QLabel("Target Object"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.glue_obj_combo = QComboBox()
        self.glue_obj_combo.setMinimumHeight(36)
        self.glue_refresh_btn = QPushButton(t("refresh"))
        self.glue_refresh_btn.setObjectName("refresh_btn")
        self.glue_refresh_btn.setMinimumHeight(36)
        self.glue_refresh_btn.clicked.connect(self.refresh_objects)
        glue_obj_row_container = QWidget()
        glue_obj_row = QHBoxLayout(glue_obj_row_container)
        glue_obj_row.setContentsMargins(0, 0, 0, 0)
        glue_obj_row.addWidget(self.glue_obj_combo, 1)
        glue_obj_row.addWidget(self.glue_refresh_btn)
        glue_grid.addWidget(glue_obj_row_container, 0, 1)
        
        glue_grid.addWidget(QLabel("Glue Residue Name:"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.glue_resname = QLineEdit()
        self.glue_resname.setMinimumHeight(36)
        self.glue_resname.setPlaceholderText("e.g., CC885, 1N6 (Lenalidomide)")
        glue_grid.addWidget(self.glue_resname, 0, 3)
        
        # 第二行: E3 Chains | [input] | Substrate Chains | [input]
        glue_grid.addWidget(QLabel("E3 Ligase Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.glue_e3_chains = QLineEdit()
        self.glue_e3_chains.setMinimumHeight(36)
        self.glue_e3_chains.setPlaceholderText("e.g., A (CRBN)")
        glue_grid.addWidget(self.glue_e3_chains, 1, 1)
        
        glue_grid.addWidget(QLabel("Substrate Chains:"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.glue_sub_chains = QLineEdit()
        self.glue_sub_chains.setMinimumHeight(36)
        self.glue_sub_chains.setPlaceholderText("e.g., B (Substrate)")
        glue_grid.addWidget(self.glue_sub_chains, 1, 3)
        
        # 第三行: Interface Distance | [input] | Neo-Epitope Distance | [input]
        glue_grid.addWidget(QLabel("Interface Distance (Å):"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.glue_interface_dist = QLineEdit("4.5")
        self.glue_interface_dist.setMinimumHeight(36)
        glue_grid.addWidget(self.glue_interface_dist, 2, 1)
        
        glue_grid.addWidget(QLabel("Neo-Epitope Dist (Å):"), 2, 2, Qt.AlignmentFlag.AlignRight)
        self.glue_neo_dist = QLineEdit("5.0")
        self.glue_neo_dist.setMinimumHeight(36)
        glue_grid.addWidget(self.glue_neo_dist, 2, 3)
        
        # 第四行: Output PPI CSV | [input + browse]
        glue_grid.addWidget(QLabel("Output PPI CSV:"), 3, 0, Qt.AlignmentFlag.AlignRight)
        self.glue_ppi_csv = QLineEdit()
        self.glue_ppi_csv.setMinimumHeight(36)
        self.glue_ppi_csv.setPlaceholderText("Optional")
        self.glue_ppi_browse = QPushButton(t("browse"))
        self.glue_ppi_browse.setMinimumHeight(36)
        self.glue_ppi_browse.setObjectName("save_btn")
        self.glue_ppi_browse.clicked.connect(lambda: self._browse_save_file(self.glue_ppi_csv, "CSV (*.csv)"))
        ppi_csv_row_container = QWidget()
        ppi_csv_row = QHBoxLayout(ppi_csv_row_container)
        ppi_csv_row.setContentsMargins(0, 0, 0, 0)
        ppi_csv_row.addWidget(self.glue_ppi_csv, 1)
        ppi_csv_row.addWidget(self.glue_ppi_browse)
        glue_grid.addWidget(ppi_csv_row_container, 3, 1)
        
        # 第四行右侧: Output Neo-Epitope CSV | [input + browse]
        glue_grid.addWidget(QLabel("Output Neo-Epitope CSV:"), 3, 2, Qt.AlignmentFlag.AlignRight)
        self.glue_neo_csv = QLineEdit()
        self.glue_neo_csv.setMinimumHeight(36)
        self.glue_neo_csv.setPlaceholderText("Optional")
        self.glue_neo_browse = QPushButton(t("browse"))
        self.glue_neo_browse.setMinimumHeight(36)
        self.glue_neo_browse.setObjectName("save_btn")
        self.glue_neo_browse.clicked.connect(lambda: self._browse_save_file(self.glue_neo_csv, "CSV (*.csv)"))
        neo_csv_row_container = QWidget()
        neo_csv_row = QHBoxLayout(neo_csv_row_container)
        neo_csv_row.setContentsMargins(0, 0, 0, 0)
        neo_csv_row.addWidget(self.glue_neo_csv, 1)
        neo_csv_row.addWidget(self.glue_neo_browse)
        glue_grid.addWidget(neo_csv_row_container, 3, 3)
        
        # Molecular Glue 按钮行
        glue_btn_row = QHBoxLayout()
        self.glue_ppi_btn = QPushButton("Analyze PPI Interface")
        self.glue_ppi_btn.setObjectName("highlight_btn")
        self.glue_ppi_btn.setToolTip("Detect protein-protein interface (key for glue vs PROTAC)")
        self.glue_ppi_btn.clicked.connect(self.run_glue_ppi_analysis)
        
        self.glue_neo_btn = QPushButton("Detect Neo-Epitope")
        self.glue_neo_btn.setObjectName("highlight_btn")
        self.glue_neo_btn.setToolTip("Identify neo-substrate epitope residues")
        self.glue_neo_btn.clicked.connect(self.run_glue_neo_epitope)
        
        self.glue_full_btn = QPushButton("Full Glue Analysis")
        self.glue_full_btn.setObjectName("highlight_btn")
        self.glue_full_btn.setToolTip("PPI + Neo-Epitope + Scoring + Classification")
        self.glue_full_btn.clicked.connect(self.run_glue_full_analysis)
        
        glue_btn_row.addWidget(self.glue_ppi_btn)
        glue_btn_row.addWidget(self.glue_neo_btn)
        glue_btn_row.addWidget(self.glue_full_btn)
        glue_btn_row.addStretch(1)
        
        main_layout.addWidget(grp_glue)
        main_layout.addLayout(glue_btn_row)
        
        main_layout.addStretch(1)
        
        # 将内容设置到滚动区域
        scroll_area.setWidget(content_widget)
        
        # 返回包裹的滚动区域
        wrapper = QWidget()
        wrapper_layout = QVBoxLayout(wrapper)
        wrapper_layout.setContentsMargins(0, 0, 0, 0)
        wrapper_layout.addWidget(scroll_area)
        return wrapper
    
    def create_gmotif_tab(self) -> QWidget:
        w = QWidget(); lay = QVBoxLayout(w); lay.setSpacing(10)
        grp = QGroupBox(t("grp_gmotif")); form = QFormLayout(grp); form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)

        r0 = QHBoxLayout()
        self.obj_combo_gm = QComboBox()
        self.refresh_obj_gm = QPushButton(t("refresh")); self.refresh_obj_gm.setObjectName("refresh_btn"); self.refresh_obj_gm.clicked.connect(self.refresh_objects)
        r0.addWidget(self.obj_combo_gm, 1); r0.addWidget(self.refresh_obj_gm)
        form.addRow(QLabel(t("target_obj")), r0)

        r1 = QHBoxLayout()
        self.gm_pdb = QLineEdit()
        self.gm_pdb_browse = QPushButton(t("browse")); self.gm_pdb_browse.setObjectName("browse_btn"); self.gm_pdb_browse.clicked.connect(self.browse_gm_pdb)
        r1.addWidget(self.gm_pdb, 1); r1.addWidget(self.gm_pdb_browse)
        form.addRow(QLabel(t("pdb_file")), r1)

        self.gm_rmsd = QLineEdit("3.5")              # Default cutoff 3.5 Å
        form.addRow(QLabel(t("rmsd")), self.gm_rmsd)
        self.gm_require_gly = QCheckBox(t("require_gly")); self.gm_require_gly.setChecked(True)
        form.addRow(QLabel(""), self.gm_require_gly)

        # --- Template Source --- (Ideal / Built-in / Custom Selection)
        self.gm_template_mode = QComboBox()
        self.gm_template_mode.addItems([
            "Idealized (8×Cα)",
            "Built-in: GSPT1 (6H0G A:60-67)",
            "Built-in: CK1α (3M51 A:36-43)",
            "Built-in: VAV1 (2MC1 A:95-102)",
            "From Selection (8×Cα)",
        ])
        form.addRow(QLabel("Template Source"), self.gm_template_mode)

        r_temp = QHBoxLayout()
        self.gm_template_sel = QLineEdit()
        self.gm_template_sel.setPlaceholderText("For 'From Selection': a selection name")
        self.gm_template_pick = QPushButton("Pick Current (sele)")
        self.gm_template_pick.setObjectName("refresh_btn")
        self.gm_template_pick.clicked.connect(lambda: self.gm_template_sel.setText("sele"))
        r_temp.addWidget(self.gm_template_sel, 1); r_temp.addWidget(self.gm_template_pick)
        form.addRow(QLabel("Template Selection"), r_temp)

        # Button row
        row_btns = QHBoxLayout()
        self.gm_btn = QPushButton(t("btn_gmotif"))
        self.gm_btn.setObjectName("highlight_btn")
        self.gm_btn.clicked.connect(self.start_gmotif)
        self.gm_btn_render = QPushButton(t("btn_gmotif_render"))
        self.gm_btn_render.setObjectName("highlight_btn")
        self.gm_btn_render.clicked.connect(self.render_gmotif_with_esp)
        row_btns.addWidget(self.gm_btn)
        row_btns.addWidget(self.gm_btn_render)
        row_btns.addStretch(1)

        # Output CSV
        r2 = QHBoxLayout()
        self.gm_out_csv = QLineEdit()
        self.gm_out_csv.setPlaceholderText("(Optional) A temp file is used if blank")
        self.gm_out_browse = QPushButton(t("browse")); self.gm_out_browse.setObjectName("save_btn"); self.gm_out_browse.clicked.connect(self.browse_gm_out_csv)
        r2.addWidget(self.gm_out_csv, 1); r2.addWidget(self.gm_out_browse)
        form.addRow(QLabel(t("output_csv")), r2)

        # Enable/disable logic for template inputs
        def _toggle_template_inputs(idx):
            use_sel = (idx == 4)
            self.gm_template_sel.setEnabled(use_sel)
            self.gm_template_pick.setEnabled(use_sel)
        self.gm_template_mode.currentIndexChanged.connect(_toggle_template_inputs)
        _toggle_template_inputs(self.gm_template_mode.currentIndex())


        lay.addWidget(grp)
        lay.addLayout(row_btns)
        lay.addStretch(1)
        return w

    def create_docking_scoring_tab(self) -> QWidget:
        """创建 Docking & Scoring 标签页 - 优化布局"""
        # 创建滚动区域
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setFrameShape(QFrame.Shape.NoFrame)
        
        # 创建内容容器
        content_widget = QWidget()
        content_widget.setObjectName("scroll_content")
        self._docking_scroll_content = content_widget  # 保存引用以便主题切换
        # 强制设置背景色（macOS 兼容性）
        bg_color = "#0d1117" if self._dark_mode else "#f8fafc"
        content_widget.setStyleSheet(f"#scroll_content {{ background-color: {bg_color}; }}")
        main_layout = QVBoxLayout(content_widget)
        main_layout.setSpacing(12)
        main_layout.setContentsMargins(12, 12, 12, 12)
        
        # ========== 口袋检测与对接（并排） ==========
        docking_row = QHBoxLayout()
        docking_row.setSpacing(12)
        
        # 口袋检测
        pocket_card = self._create_pocket_detection_card()
        docking_row.addWidget(pocket_card, 1)
        
        # Vina 对接
        docking_card = self._create_vina_docking_card()
        docking_row.addWidget(docking_card, 1)
        
        main_layout.addLayout(docking_row)
        
        # ========== 高级口袋分析（全宽） ==========
        advanced_pocket_card = self._create_advanced_pocket_card()
        main_layout.addWidget(advanced_pocket_card)
        
        # ========== 突变分析（全宽）========== ⭐ 新增
        mutation_card = self._create_mutation_analysis_card()
        main_layout.addWidget(mutation_card)
        
        main_layout.addStretch(1)
        
        # 结果显示在底部日志
        self.score_result_text = self.log_edit
        
        # 将内容设置到滚动区域
        scroll_area.setWidget(content_widget)
        
        # 返回包裹的滚动区域
        wrapper = QWidget()
        wrapper_layout = QVBoxLayout(wrapper)
        wrapper_layout.setContentsMargins(0, 0, 0, 0)
        wrapper_layout.addWidget(scroll_area)
        return wrapper
    def _create_pocket_detection_card(self) -> QWidget:
        """创建口袋检测与可视化卡片"""
        card = QGroupBox("Pocket Detection")
        layout = QVBoxLayout(card)
        layout.setSpacing(10)
        layout.setContentsMargins(16, 20, 16, 16)
        
        # Target object
        obj_row = QHBoxLayout()
        obj_row.setSpacing(8)
        self.pocket_obj_combo = QComboBox()
        self.pocket_obj_combo.setMinimumHeight(36)
        refresh_btn = QPushButton("Refresh")
        refresh_btn.setObjectName("refresh_btn")
        refresh_btn.setMinimumHeight(36)
        refresh_btn.setToolTip("Refresh objects")
        refresh_btn.clicked.connect(self.refresh_objects)
        obj_row.addWidget(QLabel("Target:"))
        obj_row.addWidget(self.pocket_obj_combo, 1)
        obj_row.addWidget(refresh_btn)
        layout.addLayout(obj_row)
        
        # Parameters row
        param_row = QGridLayout()
        param_row.setSpacing(10)
        
        param_row.addWidget(QLabel("Grid Spacing:"), 0, 0)
        self.pocket_grid_spacing = QLineEdit("0.5")
        self.pocket_grid_spacing.setMinimumHeight(36)
        self.pocket_grid_spacing.setMaximumWidth(80)
        param_row.addWidget(self.pocket_grid_spacing, 0, 1)
        param_row.addWidget(QLabel("Å"), 0, 2)
        
        param_row.addWidget(QLabel("Min Volume:"), 0, 3)
        self.pocket_min_volume = QLineEdit("30")
        self.pocket_min_volume.setMinimumHeight(36)
        self.pocket_min_volume.setMaximumWidth(80)
        param_row.addWidget(self.pocket_min_volume, 0, 4)
        param_row.addWidget(QLabel("Ų"), 0, 5)
        
        param_row.addWidget(QLabel("Color by:"), 1, 0)
        self.pocket_color_by = QComboBox()
        self.pocket_color_by.addItems(["Volume", "Druggability", "Hydrophobicity", "Depth"])
        self.pocket_color_by.setMinimumHeight(36)
        param_row.addWidget(self.pocket_color_by, 1, 1, 1, 5)
        
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
        
        self.vina_ligand = QLineEdit()
        self.vina_ligand.setPlaceholderText("Select ligand file (MOL2/SDF/PDBQT)")
        self.vina_ligand.setMinimumHeight(36)
        
        ligand_browse = QPushButton("Browse")
        ligand_browse.setObjectName("browse_btn")
        ligand_browse.setMinimumHeight(36)
        ligand_browse.setToolTip("Browse ligand file")
        ligand_browse.clicked.connect(self.browse_vina_ligand)
        
        ligand_layout.addWidget(QLabel("Ligand File:"))
        ligand_layout.addWidget(self.vina_ligand, 1)
        ligand_layout.addWidget(ligand_browse)
        layout.addLayout(ligand_layout)
        
        # 参数设置
        param_layout = QGridLayout()
        param_layout.setSpacing(10)
        
        param_layout.addWidget(QLabel("Max Pockets:"), 0, 0)
        self.vina_max_pockets = QSpinBox()
        self.vina_max_pockets.setRange(1, 10)
        self.vina_max_pockets.setValue(3)
        self.vina_max_pockets.setMinimumHeight(36)
        self.vina_max_pockets.setMaximumWidth(80)
        param_layout.addWidget(self.vina_max_pockets, 0, 1)
        
        param_layout.addWidget(QLabel("Exhaustiveness:"), 0, 2)
        self.vina_exhaustiveness = QSpinBox()
        self.vina_exhaustiveness.setRange(1, 32)
        self.vina_exhaustiveness.setValue(8)
        self.vina_exhaustiveness.setMinimumHeight(36)
        self.vina_exhaustiveness.setMaximumWidth(80)
        param_layout.addWidget(self.vina_exhaustiveness, 0, 3)
        
        layout.addLayout(param_layout)
        
        # 自定义盒子选项
        self.vina_use_custom_box = QCheckBox("Use custom box (skip pocket detection)")
        layout.addWidget(self.vina_use_custom_box)
        
        box_grid = QGridLayout()
        box_grid.setSpacing(8)
        
        # Center
        box_grid.addWidget(QLabel("center_x"), 0, 0)
        self.vina_cx = QLineEdit(); self.vina_cx.setPlaceholderText("e.g. 10.0"); self.vina_cx.setEnabled(False)
        box_grid.addWidget(self.vina_cx, 0, 1)
        box_grid.addWidget(QLabel("center_y"), 0, 2)
        self.vina_cy = QLineEdit(); self.vina_cy.setPlaceholderText("e.g. 20.0"); self.vina_cy.setEnabled(False)
        box_grid.addWidget(self.vina_cy, 0, 3)
        box_grid.addWidget(QLabel("center_z"), 0, 4)
        self.vina_cz = QLineEdit(); self.vina_cz.setPlaceholderText("e.g. 30.0"); self.vina_cz.setEnabled(False)
        box_grid.addWidget(self.vina_cz, 0, 5)
        
        # Size
        box_grid.addWidget(QLabel("size_x"), 1, 0)
        self.vina_sx = QLineEdit(); self.vina_sx.setPlaceholderText("e.g. 20.0"); self.vina_sx.setEnabled(False)
        box_grid.addWidget(self.vina_sx, 1, 1)
        box_grid.addWidget(QLabel("size_y"), 1, 2)
        self.vina_sy = QLineEdit(); self.vina_sy.setPlaceholderText("e.g. 20.0"); self.vina_sy.setEnabled(False)
        box_grid.addWidget(self.vina_sy, 1, 3)
        box_grid.addWidget(QLabel("size_z"), 1, 4)
        self.vina_sz = QLineEdit(); self.vina_sz.setPlaceholderText("e.g. 20.0"); self.vina_sz.setEnabled(False)
        box_grid.addWidget(self.vina_sz, 1, 5)
        
        layout.addLayout(box_grid)
        
        def _toggle_box_fields(checked: bool):
            for w in (self.vina_cx, self.vina_cy, self.vina_cz, self.vina_sx, self.vina_sy, self.vina_sz):
                w.setEnabled(checked)
        self.vina_use_custom_box.toggled.connect(_toggle_box_fields)
        
        # Buttons
        btn_row = QHBoxLayout()
        btn_row.setSpacing(8)
        
        self.vina_dock_btn = QPushButton("Run Docking")
        self.vina_dock_btn.setObjectName("primary_btn")
        self.vina_dock_btn.setMinimumHeight(36)
        self.vina_dock_btn.clicked.connect(self.run_vina_docking)
        
        self.vina_load_result_btn = QPushButton("Load Result")
        self.vina_load_result_btn.setObjectName("secondary_btn")
        self.vina_load_result_btn.setMinimumHeight(36)
        self.vina_load_result_btn.clicked.connect(self.load_vina_result)
        
        btn_row.addWidget(self.vina_dock_btn)
        btn_row.addWidget(self.vina_load_result_btn)
        layout.addLayout(btn_row)
        
        return card
    
    
    def _create_advanced_pocket_card(self) -> QWidget:
        """创建高级口袋分析卡片"""
        card = QGroupBox("Advanced Pocket Analysis")
        layout = QVBoxLayout(card)
        layout.setSpacing(8)
        
        # 创建标签页
        self.pocket_advanced_tabs = QTabWidget()
        self.pocket_advanced_tabs.setFixedHeight(220)
        
        # Tab 1: 口袋对比
        comparison_tab = self._create_pocket_comparison_tab()
        self.pocket_advanced_tabs.addTab(comparison_tab, "Comparison")
        
        # Tab 2: 界面口袋
        interface_tab = self._create_pocket_interface_tab()
        self.pocket_advanced_tabs.addTab(interface_tab, "PPI Interface")
        
        # Tab 3: 口袋-相互作用关联
        correlation_tab = self._create_pocket_correlation_tab()
        self.pocket_advanced_tabs.addTab(correlation_tab, "Interactions")
        
        # Tab 4: G-motif 口袋
        gmotif_pocket_tab = self._create_gmotif_pocket_tab()
        self.pocket_advanced_tabs.addTab(gmotif_pocket_tab, "G-motif")
        
        # 强制设置所有 tab 的背景色（macOS Qt 兼容性）
        bg_color = "#161b22" if self._dark_mode else "white"
        for i in range(self.pocket_advanced_tabs.count()):
            tab_widget = self.pocket_advanced_tabs.widget(i)
            if tab_widget:
                tab_widget.setObjectName("tab_content")
                tab_widget.setStyleSheet(f"#tab_content {{ background-color: {bg_color}; }}")
        
        layout.addWidget(self.pocket_advanced_tabs)
        
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
        
        self.pocket_comp_obj_a = QComboBox()
        self.pocket_comp_obj_a.setFixedHeight(28)
        self.pocket_comp_obj_b = QComboBox()
        self.pocket_comp_obj_b.setFixedHeight(28)
        
        refresh_a = QPushButton("Refresh")
        refresh_a.setObjectName("refresh_btn")
        refresh_a.setMinimumHeight(28)
        refresh_a.clicked.connect(self.refresh_objects)
        
        refresh_b = QPushButton("Refresh")
        refresh_b.setObjectName("refresh_btn")
        refresh_b.setMinimumHeight(28)
        refresh_b.clicked.connect(self.refresh_objects)
        
        obj_grid.addWidget(QLabel("Object A:"), 0, 0)
        obj_grid.addWidget(self.pocket_comp_obj_a, 0, 1)
        obj_grid.addWidget(refresh_a, 0, 2)
        obj_grid.addWidget(QLabel("Object B:"), 1, 0)
        obj_grid.addWidget(self.pocket_comp_obj_b, 1, 1)
        obj_grid.addWidget(refresh_b, 1, 2)
        layout.addLayout(obj_grid)
        
        # Align checkbox
        self.pocket_comp_align = QCheckBox("Align structures before comparison")
        self.pocket_comp_align.setChecked(True)
        layout.addWidget(self.pocket_comp_align)
        
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
        self.pocket_interface_obj = QComboBox()
        self.pocket_interface_obj.setFixedHeight(28)
        refresh_btn = QPushButton("Refresh")
        refresh_btn.setObjectName("refresh_btn")
        refresh_btn.setMinimumHeight(28)
        refresh_btn.clicked.connect(self.refresh_objects)
        obj_row.addWidget(QLabel("Object:"))
        obj_row.addWidget(self.pocket_interface_obj, 1)
        obj_row.addWidget(refresh_btn)
        layout.addLayout(obj_row)
        
        # Chains
        chain_row = QHBoxLayout()
        chain_row.setSpacing(6)
        self.pocket_interface_chain_a = QLineEdit()
        self.pocket_interface_chain_a.setPlaceholderText("Chain A")
        self.pocket_interface_chain_a.setFixedHeight(28)
        self.pocket_interface_chain_a.setFixedWidth(60)
        
        self.pocket_interface_chain_b = QLineEdit()
        self.pocket_interface_chain_b.setPlaceholderText("Chain B")
        self.pocket_interface_chain_b.setFixedHeight(28)
        self.pocket_interface_chain_b.setFixedWidth(60)
        
        chain_row.addWidget(QLabel("Chains:"))
        chain_row.addWidget(self.pocket_interface_chain_a)
        chain_row.addWidget(QLabel("+"))
        chain_row.addWidget(self.pocket_interface_chain_b)
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
        self.mut_obj_combo = QComboBox()
        self.mut_obj_combo.setFixedHeight(28)
        self.mut_refresh_btn = QPushButton("Refresh")
        self.mut_refresh_btn.setObjectName("refresh_btn")
        self.mut_refresh_btn.setFixedHeight(28)
        self.mut_refresh_btn.clicked.connect(self.refresh_objects)
        obj_row = QHBoxLayout()
        obj_row.addWidget(self.mut_obj_combo, 1)
        obj_row.addWidget(self.mut_refresh_btn)
        input_grid.addLayout(obj_row, 0, 1)
        
        # 突变输入
        input_grid.addWidget(QLabel("Mutations:"), 1, 0)
        self.mut_input = QLineEdit()
        self.mut_input.setPlaceholderText("e.g., A:23:ALA, B:45:GLY")
        self.mut_input.setFixedHeight(28)
        input_grid.addWidget(self.mut_input, 1, 1)
        
        # 方法选择 (FoldX only)
        input_grid.addWidget(QLabel("Method:"), 2, 0)
        method_row = QHBoxLayout()
        self.mut_method_combo = QComboBox()
        self.mut_method_combo.addItems(["FoldX"])
        self.mut_method_combo.setFixedHeight(28)
        self.mut_method_combo.setFixedWidth(120)
        self.mut_method_combo.setToolTip("ΔΔG calculation requires FoldX\nDownload: https://foldxsuite.crg.eu/")
        method_row.addWidget(self.mut_method_combo)
        method_row.addStretch()
        input_grid.addLayout(method_row, 2, 1)
        
        layout.addLayout(input_grid)
        
        # 按钮行
        btn_row = QHBoxLayout()
        btn_row.setSpacing(6)
        
        self.mut_perform_btn = QPushButton("Perform Mutation")
        self.mut_perform_btn.setObjectName("primary_btn")
        self.mut_perform_btn.setFixedHeight(32)
        self.mut_perform_btn.clicked.connect(self.run_mutation)
        
        self.mut_minimize_btn = QPushButton("Minimize Energy")
        self.mut_minimize_btn.setObjectName("highlight_btn")
        self.mut_minimize_btn.setFixedHeight(32)
        self.mut_minimize_btn.clicked.connect(self.run_minimize)
        
        self.mut_analyze_btn = QPushButton("Full Analysis")
        self.mut_analyze_btn.setObjectName("highlight_btn")
        self.mut_analyze_btn.setFixedHeight(32)
        self.mut_analyze_btn.clicked.connect(self.run_mutation_analysis)
        
        btn_row.addWidget(self.mut_perform_btn)
        btn_row.addWidget(self.mut_minimize_btn)
        btn_row.addWidget(self.mut_analyze_btn)
        btn_row.addStretch()
        
        layout.addLayout(btn_row)
        
        return card
    
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
        self.pocket_corr_csv = QLineEdit()
        self.pocket_corr_csv.setPlaceholderText("Interaction CSV file")
        self.pocket_corr_csv.setFixedHeight(28)
        
        csv_browse = QPushButton("Browse")
        csv_browse.setObjectName("browse_btn")
        csv_browse.setMinimumHeight(28)
        csv_browse.clicked.connect(self.browse_pocket_corr_csv)
        
        csv_row.addWidget(self.pocket_corr_csv, 1)
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
        self.gmotif_pocket_obj = QComboBox()
        self.gmotif_pocket_obj.setFixedHeight(28)
        refresh_btn = QPushButton("R")
        refresh_btn.setObjectName("refresh_btn")
        refresh_btn.setFixedSize(28, 28)
        refresh_btn.clicked.connect(self.refresh_objects)
        obj_row.addWidget(QLabel("Object:"))
        obj_row.addWidget(self.gmotif_pocket_obj, 1)
        obj_row.addWidget(refresh_btn)
        layout.addLayout(obj_row)
        
        # Chains
        chain_grid = QGridLayout()
        chain_grid.setSpacing(6)
        
        self.gmotif_pocket_e3 = QLineEdit()
        self.gmotif_pocket_e3.setPlaceholderText("E3 (e.g., A)")
        self.gmotif_pocket_e3.setFixedHeight(28)
        self.gmotif_pocket_e3.setFixedWidth(70)
        
        self.gmotif_pocket_sub = QLineEdit()
        self.gmotif_pocket_sub.setPlaceholderText("Substrate (e.g., B)")
        self.gmotif_pocket_sub.setFixedHeight(28)
        self.gmotif_pocket_sub.setFixedWidth(70)
        
        self.gmotif_pocket_glue = QLineEdit()
        self.gmotif_pocket_glue.setPlaceholderText("Glue (optional)")
        self.gmotif_pocket_glue.setFixedHeight(28)
        self.gmotif_pocket_glue.setFixedWidth(70)
        
        chain_grid.addWidget(QLabel("E3:"), 0, 0)
        chain_grid.addWidget(self.gmotif_pocket_e3, 0, 1)
        chain_grid.addWidget(QLabel("Sub:"), 0, 2)
        chain_grid.addWidget(self.gmotif_pocket_sub, 0, 3)
        chain_grid.addWidget(QLabel("Glue:"), 1, 0)
        chain_grid.addWidget(self.gmotif_pocket_glue, 1, 1)
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

    def create_apbs_tab(self) -> QWidget:
        w = QWidget(); lay = QVBoxLayout(w); lay.setSpacing(10)
        grp = QGroupBox(t("grp_apbs")); form = QFormLayout(grp); form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        r0 = QHBoxLayout()
        self.obj_combo_apbs = QComboBox()
        self.refresh_obj_apbs = QPushButton(t("refresh")); self.refresh_obj_apbs.setObjectName("refresh_btn"); self.refresh_obj_apbs.clicked.connect(self.refresh_objects)
        r0.addWidget(self.obj_combo_apbs, 1); r0.addWidget(self.refresh_obj_apbs)
        form.addRow(QLabel(t("apbs_target")), r0)
        self.apbs_grid = QLineEdit("1.0"); form.addRow(QLabel(t("apbs_grid")), self.apbs_grid)
        self.apbs_range = QLineEdit("-5,0,5"); form.addRow(QLabel(t("apbs_range")), self.apbs_range)
        
        btn_row = QHBoxLayout()
        self.btn_apbs_quick = QPushButton(t("btn_quick"))
        self.btn_apbs_quick.setObjectName("highlight_btn")
        self.btn_apbs_quick.clicked.connect(self.apbs_run_quick)
        self.btn_apbs_true = QPushButton(t("btn_apbs"))
        self.btn_apbs_true.setObjectName("highlight_btn")
        self.btn_apbs_true.clicked.connect(self.apbs_run_true)
        btn_row.addWidget(self.btn_apbs_quick)
        btn_row.addWidget(self.btn_apbs_true)
        btn_row.addStretch(1)
        
        grp_exp = QGroupBox(t("grp_export"))
        exp_form = QFormLayout(grp_exp)
        
        exp_row1 = QHBoxLayout()
        self.png_w = QLineEdit("3000"); self.png_h = QLineEdit("2000"); self.png_dpi = QLineEdit("300")
        exp_row1.addWidget(QLabel(t("img_w"))); exp_row1.addWidget(self.png_w)
        exp_row1.addWidget(QLabel(t("img_h"))); exp_row1.addWidget(self.png_h)
        exp_row1.addWidget(QLabel(t("img_dpi"))); exp_row1.addWidget(self.png_dpi)
        exp_form.addRow(exp_row1)
        
        exp_row2 = QHBoxLayout()
        self.bg_white = QCheckBox(t("bg_white")); self.bg_white.setChecked(True)
        self.bg_trans = QCheckBox(t("bg_trans"))
        self.bg_white.toggled.connect(lambda v: (self.bg_trans.setChecked(False) if v else None))
        self.bg_trans.toggled.connect(lambda v: (self.bg_white.setChecked(False) if v else None))
        self.chk_ray = QCheckBox(t("raytrace")); self.chk_ray.setChecked(True)
        self.btn_viewport = QPushButton(t("btn_viewport")); self.btn_viewport.setObjectName("refresh_btn"); self.btn_viewport.clicked.connect(self.fill_viewport_size)
        exp_row2.addWidget(self.bg_white); exp_row2.addWidget(self.bg_trans); exp_row2.addWidget(self.chk_ray); exp_row2.addStretch(1); exp_row2.addWidget(self.btn_viewport)
        exp_form.addRow(exp_row2)
        
        exp_row3 = QHBoxLayout()
        self.btn_export_png = QPushButton(t("btn_export_png")); self.btn_export_png.setObjectName("save_btn"); self.btn_export_png.clicked.connect(self.export_png)
        self.btn_export_dx = QPushButton(t("btn_export_dx")); self.btn_export_dx.setObjectName("save_btn"); self.btn_export_dx.clicked.connect(self.export_dx)
        exp_row3.addWidget(self.btn_export_png); exp_row3.addWidget(self.btn_export_dx); exp_row3.addStretch(1)
        exp_form.addRow(exp_row3)
        lay.addWidget(grp)
        lay.addLayout(btn_row)
        lay.addWidget(grp_exp)
        lay.addStretch(1)
        return w

    # --- 行为与回调 ---
    def refresh_objects(self):
        names = []
        try:
            from pymol import cmd
            # Prefer standard API: cmd.get_names("objects")
            # Fallback to legacy/custom cmd.get_object_list() if available
            if hasattr(cmd, "get_names"):
                names = cmd.get_names("objects")
            elif hasattr(cmd, "get_object_list"):
                names = cmd.get_object_list()
            else:
                names = []
        except Exception as e:
            self.log(f"Error refreshing objects: {e}")
            pass
        
        # Ensure names is a list
        if names is None: names = []
        
        if not names: names = [t("no_object")]
        for cb in (getattr(self, "obj_combo_gm", None),
                   getattr(self, "obj_combo_analysis", None),
                   getattr(self, "obj_combo_csv", None),
                   getattr(self, "obj_combo_apbs", None),
                   getattr(self, "pl_obj_combo", None),
                   getattr(self, "tc_obj_combo", None),
                   getattr(self, "ll_obj_combo", None),
                   getattr(self, "pn_obj_combo", None),
                   getattr(self, "ap_obj_combo", None),
                   getattr(self, "score_obj", None),
                   getattr(self, "pocket_obj_combo", None),
                   getattr(self, "pocket_comp_obj_a", None),
                   getattr(self, "pocket_comp_obj_b", None),
                   getattr(self, "pocket_interface_obj", None),
                   getattr(self, "gmotif_pocket_obj", None),
                   getattr(self, "glue_obj_combo", None),
                   getattr(self, "mut_obj_combo", None)):  # 添加突变分析对象下拉框
            if cb is not None:
                cb.blockSignals(True); cb.clear()
                for n in names: cb.addItem(n)
                cb.blockSignals(False)
        self.update_enablement()
        self.log(f"{t('refresh')} OK: {', '.join(names)}")

    def browse_csv(self):
        fn, _ = QFileDialog.getOpenFileName(self, t("select_csv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.csv_path.setText(fn); self.update_enablement()

    def browse_out_csv(self):
        fn, _ = QFileDialog.getSaveFileName(self, t("select_outcsv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.out_csv.setText(fn); self.update_enablement()

    def _browse_file(self, line_edit, filter_str):
        fn, _ = QFileDialog.getOpenFileName(self, "Select File", "", filter_str)
        if fn:
            line_edit.setText(fn)
            
    def _browse_save_file(self, line_edit, filter_str):
        fn, _ = QFileDialog.getSaveFileName(self, "Save File", "", filter_str)
        if fn:
            line_edit.setText(fn)

    def browse_gm_pdb(self):
        fn, _ = QFileDialog.getOpenFileName(self, t("select_pdb"), "", "PDB (*.pdb *.cif);;All Files (*)")
        if fn: self.gm_pdb.setText(fn); self.update_enablement()

    def browse_gm_out_csv(self):
        fn, _ = QFileDialog.getSaveFileName(self, t("select_outcsv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.gm_out_csv.setText(fn); self.update_enablement()

    def highlight_csv_clicked(self):
        try:
            obj = self.obj_combo_csv.currentText().strip()
            csv_file = self.csv_path.text().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object")); return
            if not csv_file or not os.path.exists(csv_file):
                QMessageBox.warning(self, t("title"), t("select_csv")); return

            self.log(f"Start highlighting: {os.path.basename(csv_file)}")
            result = highlight_csv_residues(csv_file, obj, show_labels=1, clear_old=1, debug=0, stick_by_element=1)

            if result and isinstance(result, dict):
                pairs = result.get("pairs", 0)
                unique = result.get("unique_residues", 0)
                self.log(f"Highlight done: {pairs} pairs, {unique} unique residues")
            else:
                self.log(t("log_highlight_ok"))
        except Exception as e:
            self.on_error(str(e))
            import traceback
            traceback.print_exc()

    def clear_highlight_clicked(self):
        try:
            self.log(t("log_clear_ok"))
        except Exception as e:
            self.on_error(str(e))

    def start_analysis(self):
        obj = self.obj_combo_analysis.currentText().strip()
        out_csv = self.out_csv.text().strip() or None
        if not obj or obj == t("no_object"):
            QMessageBox.warning(self, t("title"), t("no_object")); return
        self.analyze_btn.setEnabled(False)
        self.progress_bar.setVisible(True); self.progress_bar.setRange(0, 0)
        self.analysis_thread = AnalysisWorker(obj, None, out_csv)  # pdb_file=None，只从PyMOL加载
        self.analysis_thread.progress.connect(self.log)
        self.analysis_thread.error.connect(self.on_error)
        self.analysis_thread.finished.connect(self.on_finished_analysis)
        self.analysis_thread.start()

    def start_gmotif(self):
        obj = self.obj_combo_gm.currentText().strip()
        pdb = self.gm_pdb.text().strip() or None
        outcsv = self.gm_out_csv.text().strip() or None
        if not obj or obj == t("no_object"):
            QMessageBox.warning(self, t("title"), t("no_object")); return
        try:
            rmsd = float(self.gm_rmsd.text().strip() or "3.5")
        except Exception:
            rmsd = 3.5
        require_gly = self.gm_require_gly.isChecked()

        # Template source mapping
        idx = self.gm_template_mode.currentIndex()
        if idx == 0:
            template_mode, template_sel, template_builtin = "ideal", None, None
        elif idx in (1, 2, 3):
            template_mode = "builtin"
            template_sel = None
            template_builtin = [
                "GSPT1 (6H0G A:60-67)",
                "CK1α (3M51 A:36-43)",
                "VAV1 (2MC1 A:95-102)",
            ][idx - 1]
        else:
            template_mode, template_sel, template_builtin = "selection", (self.gm_template_sel.text().strip() or None), None

        self.gm_btn.setEnabled(False)
        self.progress_bar.setVisible(True); self.progress_bar.setRange(0, 0)
        self.gmotif_thread = GMotifWorker(obj, pdb, rmsd, require_gly, outcsv,
                                        template_mode, template_sel, template_builtin)
        self.gmotif_thread.progress.connect(self.log)
        self.gmotif_thread.error.connect(self.on_error)
        self.gmotif_thread.finished.connect(self.on_finished_gmotif)
        self.gmotif_thread.start()

    # ---- 一键渲染（G-Motif + 电势 + PNG）----
    def render_gmotif_with_esp(self):
        try:
            obj = self.obj_combo_gm.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object")); return

            from pymol import cmd

            # 1) 若有最近一次 G-Motif CSV，使用专用高亮函数
            if self._last_gmotif_csv and os.path.exists(self._last_gmotif_csv):
                try:
                    highlight_gmotif_loops(self._last_gmotif_csv, obj, color="yellow", show_labels=True)
                    self.log("Applied G-Motif highlight (latest CSV)")
                except Exception as e:
                    self.log(f"G-Motif highlight skipped: {e}")

            # 2) 生成电势并着色
            if getattr(self, "obj_combo_apbs", None):
                idx = self.obj_combo_apbs.findText(obj)
                if idx >= 0: self.obj_combo_apbs.setCurrentIndex(idx)

            try:
                grid = float(self.apbs_grid.text().strip() or "1.0")
            except Exception:
                grid = 1.0
            rng_text = (self.apbs_range.text().strip() or "-5,0,5")
            try:
                vmin, v0, vmax = [float(x) for x in rng_text.split(",")]
            except Exception:
                vmin, v0, vmax = -5.0, 0.0, 5.0

            # 生成电势图
            map_name = f"{obj}_esp_map"
            ramp_name = f"{obj}_esp_ramp"
            cmd.map_new(map_name, "coulomb", grid, obj)
            cmd.ramp_new(ramp_name, map_name, [vmin, v0, vmax], ["blue", "white", "red"])
            self._esp_maps[obj] = (map_name, ramp_name)

            # 3) 显示整体 cartoon（半透明）
            cmd.hide("everything", obj)
            cmd.show("cartoon", obj)
            cmd.set("cartoon_transparency", 0.3, obj)

            # 4) 只在 G-loop 周围显示静电势表面
            if self._last_gmotif_csv and os.path.exists(self._last_gmotif_csv):
                try:
                    import csv
                    # 读取所有 G-loop 区域
                    gloop_regions = []
                    with open(self._last_gmotif_csv, "r", encoding="utf-8") as f:
                        reader = csv.DictReader(f)
                        for r in reader:
                            chain = r.get("Chain", "")
                            start = r.get("Start", "")
                            end = r.get("End", "")
                            if chain and start and end:
                                gloop_regions.append((chain, start, end))

                    if gloop_regions:
                        # 创建一个包含所有 G-loop 的选择
                        gloop_selections = [f"(chain {ch} and resi {st}-{ed})" for ch, st, ed in gloop_regions]
                        all_gloop_sel = " or ".join(gloop_selections)

                        # 创建 G-loop 周围 10 埃的选择（用于显示静电势表面）
                        surface_sel_name = "gloop_surface_area"
                        cmd.select(surface_sel_name, f"byres ({obj} within 10 of ({all_gloop_sel}))")

                        # 只在这个区域显示表面并用电势着色
                        cmd.show("surface", surface_sel_name)
                        cmd.set("surface_quality", 1, surface_sel_name)
                        cmd.set("surface_color_smoothing", 1, surface_sel_name)
                        cmd.set("transparency", 0.2, surface_sel_name)
                        cmd.color(ramp_name, surface_sel_name)

                        # 在 G-loop 自身隐藏表面（只显示 sticks）
                        for ch, st, ed in gloop_regions:
                            gloop_sel = f"{obj} and chain {ch} and resi {st}-{ed}"
                            cmd.hide("surface", gloop_sel)

                        self.log("Electrostatic surface shown only around G-loop (10 Å region)")
                    else:
                        # 没有 G-loop 数据，显示整个蛋白表面
                        self._show_full_surface_esp(obj, ramp_name)
                except Exception as e:
                    self.log(f"ESP around G-loop failed; fallback to full surface: {e}")
                    self._show_full_surface_esp(obj, ramp_name)
            else:
                # No G-Motif CSV; show full protein surface
                self._show_full_surface_esp(obj, ramp_name)


            # 4) 设置光照参数
            cmd.set("ambient", 0.2)
            cmd.set("spec_power", 80)
            cmd.set("spec_reflect", 0.3)
            cmd.set("depth_cue", 1)
            cmd.set("fog_start", 0.45)
            cmd.orient(obj)

            self.log(f"Rendered: G-Motif + ESP (grid={grid} Å, range=({vmin},{v0},{vmax}))")

            # 5) 导出 PNG
            self._export_png_for_object(obj)

        except Exception as e:
            self.on_error(str(e))

    def _show_full_surface_esp(self, obj: str, ramp_name: str):
        """Show full-protein electrostatic surface (fallback)."""
        from pymol import cmd
        cmd.show("surface", obj)
        cmd.set("surface_quality", 1, obj)
        cmd.set("surface_color_smoothing", 1, obj)
        cmd.set("transparency", 0.2, obj)
        cmd.color(ramp_name, obj)
        self.log("Showing full-protein electrostatic surface (no G-loop data)")


    # --- 线程回调 ---
    def on_finished_analysis(self, interactions: List[Dict[str, Any]]):
        self._interactions = interactions or []
        self.log(f"✅ Analysis complete: {len(interactions)} interactions found")
        
        if interactions:
            self._last_pp_interactions = interactions
            msg = f"Found {len(interactions)} protein-protein interactions."
            QMessageBox.information(self, "Success", msg)

    def run_pp_heatmap(self):
        if not hasattr(self, '_last_pp_interactions') or not self._last_pp_interactions:
            QMessageBox.warning(self, "Data Missing", "Please run Protein-Protein analysis first.")
            return
            
        try:
            from .interaction_analyzer import generate_interaction_heatmap
            
            fn, _ = QFileDialog.getSaveFileName(self, "Save Heatmap", "ppi_heatmap.png", "PNG (*.png)")
            if fn:
                output_path = generate_interaction_heatmap(
                    interactions_result=self._last_pp_interactions,
                    output_path=fn,
                    show_plot=True
                )
                if output_path:
                    self.log(f"✅ Heatmap saved: {output_path}")
        except Exception as e:
            self.on_error(f"Heatmap failed: {e}")
        if self.out_csv.text().strip() and os.path.exists(self.out_csv.text().strip()):
            self.log(f"   Saved to: {os.path.basename(self.out_csv.text().strip())}")
        self.progress_bar.setVisible(False); self.progress_bar.setRange(0, 1)
        self.analyze_btn.setEnabled(True)

    def on_finished_gmotif(self, hits: List[Tuple], out_csv_path: str):
        self._gmotif_hits = hits or []
        self._last_gmotif_csv = out_csv_path
        self.log(f"✅ G-Motif detection complete: {len(hits)} hits found")
        if os.path.exists(out_csv_path):
            self.csv_path.setText(out_csv_path)
            self.log(f"   Saved to: {os.path.basename(out_csv_path)}")
        self.progress_bar.setVisible(False); self.progress_bar.setRange(0, 1)
        self.gm_btn.setEnabled(True)

    def on_error(self, msg: str):
        self.log(t("log_error").format(msg=msg))
        QMessageBox.critical(self, t("title"), msg)
        self.progress_bar.setVisible(False); self.progress_bar.setRange(0, 1)
        for b in (getattr(self, "analyze_btn", None), getattr(self, "gm_btn", None), getattr(self, "gm_btn_render", None), 
                  getattr(self, "pn_analyze_btn", None), getattr(self, "ll_analyze_btn", None), 
                  getattr(self, "tc_analyze_btn", None)):
            if b: b.setEnabled(True)


    # --- APBS/Quick & 导出 ---
    def apbs_run_quick(self):
        try:
            obj = self.obj_combo_apbs.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object")); return
            try:
                grid = float(self.apbs_grid.text().strip() or "1.0")
            except Exception:
                grid = 1.0
            rng_text = (self.apbs_range.text().strip() or "-5,0,5")
            try:
                vmin, v0, vmax = [float(x) for x in rng_text.split(",")]
            except Exception:
                vmin, v0, vmax = -5.0, 0.0, 5.0

            from pymol import cmd
            map_name = f"{obj}_esp_map"
            ramp_name = f"{obj}_esp_ramp"
            cmd.map_new(map_name, "coulomb", grid, obj)
            cmd.ramp_new(ramp_name, map_name, [vmin, v0, vmax], ["blue", "white", "red"])
            cmd.show("surface", obj)
            cmd.color(ramp_name, obj)
            cmd.set("surface_quality", 1, obj)
            cmd.set("surface_color_smoothing", 1, obj)
            self._esp_maps[obj] = (map_name, ramp_name)
            self.log(f"🔷 Quick ESP: map={map_name}, ramp={ramp_name}, range=({vmin},{v0},{vmax}), grid={grid} Å")
        except Exception as e:
            self.on_error(str(e))

    def apbs_run_true(self):
        try:
            obj = self.obj_combo_apbs.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object")); return

            apbs_tools = None
            try:
                from pymol.plugins import apbs_tools as _apbs
                apbs_tools = _apbs
            except Exception:
                try:
                    import apbs_tools as _apbs
                    apbs_tools = _apbs
                except Exception:
                    apbs_tools = None

            if apbs_tools is None or not hasattr(apbs_tools, "run_apbs"):
                self.log("APBS tools unavailable (apbs_tools.run_apbs not found). Falling back to Quick mode.")
                self.apbs_run_quick()
                return

            try:
                apbs_tools.run_apbs(selection=obj)
                self.log("APBS job submitted; if no visualization appears, check external paths in APBS Tools.")
            except Exception as ee:
                self.log(f"APBS call failed: {ee}. Falling back to Quick mode.")
                self.apbs_run_quick()
        except Exception as e:
            self.on_error(str(e))

    def fill_viewport_size(self):
        try:
            from pymol import cmd
            w, h = cmd.get_viewport()
            if w and h:
                self.png_w.setText(str(int(w)))
                self.png_h.setText(str(int(h)))
                self.log(f"Viewport: {w}x{h}px → filled into export settings")
            else:
                self.log("Failed to get viewport size; kept defaults")
        except Exception as e:
            self.on_error(str(e))

    def _export_png_for_object(self, obj: str):
        """按导出面板参数导出当前视图 PNG。"""
        from pymol import cmd
        fn, _ = QFileDialog.getSaveFileName(self, t("btn_export_png"), f"{obj}_gmotif_esp.png", "PNG (*.png)")
        if not fn: return
        if not fn.lower().endswith(".png"): fn += ".png"
        try:
            W = int(float(self.png_w.text().strip())); H = int(float(self.png_h.text().strip()))
        except Exception:
            W, H = 3000, 2000
        try:
            dpi = int(float(self.png_dpi.text().strip()))
        except Exception:
            dpi = 300
        ray = 1 if self.chk_ray.isChecked() else 0
        want_trans = self.bg_trans.isChecked()
        want_white = self.bg_white.isChecked() or not want_trans

        old_bg = cmd.get("bg_rgb")
        old_ray_bg = cmd.get("ray_opaque_background")
        if want_trans:
            cmd.bg_color("white"); cmd.set("ray_opaque_background", 0)
        elif want_white:
            cmd.bg_color("white"); cmd.set("ray_opaque_background", 1)

        cmd.png(fn, width=W, height=H, dpi=dpi, ray=ray)
        self.log(f"PNG export done: {os.path.basename(fn)} | {W}x{H}px @ {dpi} dpi | ray={ray} | bg={'transparent' if want_trans else 'white'}")

        # 恢复背景设置
        try:
            if isinstance(old_bg, (list, tuple)) and len(old_bg) == 3:
                r, g, b = [max(0.0, min(1.0, float(c))) for c in old_bg]
                cmd.set("bg_rgb", [r, g, b])
            cmd.set("ray_opaque_background", int(old_ray_bg))
        except Exception:
            pass

    def export_png(self):
        """从 APBS 页显式导出 PNG（与一键渲染共享同一套参数）。"""
        try:
            obj = self.obj_combo_apbs.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object")); return
            self._export_png_for_object(obj)
        except Exception as e:
            self.on_error(str(e))


    # ==============  New Workflow Tabs (Refactored) ==============\n    \n    def create_target_discovery_tab(self) -> QWidget:\n        \"\"\"Target Discovery Tab: G-Motif, Disease, Pocket\"\"\"\n        scroll_area = QScrollArea()\n        scroll_area.setWidgetResizable(True)\n        scroll_area.setFrameShape(QFrame.Shape.NoFrame)\n        \n        content_widget = QWidget()\n        content_widget.setObjectName(\"scroll_content\")\n        self._target_scroll_content = content_widget\n        bg_color = \"#0d1117\" if self._dark_mode else \"#f8fafc\"\n        content_widget.setStyleSheet(f\"#scroll_content {{ background-color: {bg_color}; }}\")\n        \n        layout = QVBoxLayout(content_widget)\n        layout.setSpacing(14)\n        layout.setContentsMargins(12, 12, 12, 12)\n        \n        # 1. Disease Analysis (if available)\n        if search_disease is not None:\n            try:\n                disease_card = self.create_disease_analysis_tab() \n                grp_disease = QGroupBox(\"Disease Target Analysis\")\n                d_layout = QVBoxLayout(grp_disease)\n                d_layout.addWidget(disease_card)\n                layout.addWidget(grp_disease)\n            except Exception as e:\n                self.log(f\"Failed to load Disease Analysis: {e}\")\n        \n        # 2. G-Motif Detection (Refactored from create_molecular_glue_tab)\n        grp_gm = QGroupBox(\"G-Motif (CRBN G-loop) Detection\")\n        gm_grid = QGridLayout(grp_gm)\n        gm_grid.setColumnStretch(1, 1); gm_grid.setColumnStretch(3, 1)\n        gm_grid.setHorizontalSpacing(8); gm_grid.setVerticalSpacing(10)\n        \n        # Row 0\n        gm_grid.addWidget(QLabel(\"Target Object:\"), 0, 0, Qt.AlignmentFlag.AlignRight)\n        self.obj_combo_gm = QComboBox(); self.obj_combo_gm.setMinimumHeight(32)\n        self.refresh_obj_gm = QPushButton(t(\"refresh\")); self.refresh_obj_gm.clicked.connect(self.refresh_objects)\n        r0 = QHBoxLayout(); r0.addWidget(self.obj_combo_gm, 1); r0.addWidget(self.refresh_obj_gm)\n        gm_grid.addLayout(r0, 0, 1)\n        \n        gm_grid.addWidget(QLabel(\"PDB File (opt):\"), 0, 2, Qt.AlignmentFlag.AlignRight)\n        self.gm_pdb = QLineEdit(); self.gm_pdb_browse = QPushButton(t(\"browse\"))\n        self.gm_pdb_browse.clicked.connect(self.browse_gm_pdb)\n        r0b = QHBoxLayout(); r0b.addWidget(self.gm_pdb, 1); r0b.addWidget(self.gm_pdb_browse)\n        gm_grid.addLayout(r0b, 0, 3)\n        \n        # Row 1\n        gm_grid.addWidget(QLabel(\"RMSD cutoff (Å):\"), 1, 0, Qt.AlignmentFlag.AlignRight)\n        self.gm_rmsd = QLineEdit(\"3.5\")\n        gm_grid.addWidget(self.gm_rmsd, 1, 1)\n        \n        self.gm_require_gly = QCheckBox(t(\"require_gly\")); self.gm_require_gly.setChecked(True)\n        gm_grid.addWidget(self.gm_require_gly, 1, 3)\n        \n        # Row 2\n        gm_grid.addWidget(QLabel(\"Template:\"), 2, 0, Qt.AlignmentFlag.AlignRight)\n        self.gm_template_mode = QComboBox()\n        self.gm_template_mode.addItems([\"Idealized (8×Cα)\", \"Built-in: GSPT1\", \"Built-in: CK1α\", \"Built-in: VAV1\", \"From Selection\"])\n        gm_grid.addWidget(self.gm_template_mode, 2, 1)\n        \n        gm_grid.addWidget(QLabel(\"Selection:\"), 2, 2, Qt.AlignmentFlag.AlignRight)\n        self.gm_template_sel = QLineEdit()\n        self.gm_template_pick = QPushButton(\"Pick (sele)\"); self.gm_template_pick.clicked.connect(lambda: self.gm_template_sel.setText(\"sele\"))\n        r2b = QHBoxLayout(); r2b.addWidget(self.gm_template_sel, 1); r2b.addWidget(self.gm_template_pick)\n        gm_grid.addLayout(r2b, 2, 3)\n        \n        def _toggle_template_inputs(idx):\n            use_sel = (idx == 4)\n            self.gm_template_sel.setEnabled(use_sel); self.gm_template_pick.setEnabled(use_sel)\n        self.gm_template_mode.currentIndexChanged.connect(_toggle_template_inputs)\n        _toggle_template_inputs(0)\n        \n        # Row 3\n        gm_grid.addWidget(QLabel(\"Output CSV:\"), 3, 0, Qt.AlignmentFlag.AlignRight)\n        self.gm_out_csv = QLineEdit()\n        self.gm_out_browse = QPushButton(t(\"browse\")); self.gm_out_browse.clicked.connect(self.browse_gm_out_csv)\n        r3 = QHBoxLayout(); r3.addWidget(self.gm_out_csv, 1); r3.addWidget(self.gm_out_browse)\n        gm_grid.addLayout(r3, 3, 1, 1, 3)\n        \n        # Buttons\n        gm_btn_row = QHBoxLayout()\n        self.gm_btn = QPushButton(\"Detect POI\"); self.gm_btn.setObjectName(\"highlight_btn\")\n        self.gm_btn.clicked.connect(self.start_gmotif)\n        self.gm_btn_render = QPushButton(\"Render All (POI + ESP + PNG)\"); self.gm_btn_render.setObjectName(\"highlight_btn\")\n        self.gm_btn_render.clicked.connect(self.render_gmotif_with_esp)\n        gm_btn_row.addWidget(self.gm_btn); gm_btn_row.addWidget(self.gm_btn_render); gm_btn_row.addStretch(1)\n        \n        layout.addWidget(grp_gm)\n        layout.addLayout(gm_btn_row)\n        \n        # 3. Pocket Detection (Existing logic)\n        layout.addWidget(self._create_pocket_detection_card())\n        \n        # 4. Advanced Pocket Analysis (Existing logic)\n        layout.addWidget(self._create_advanced_pocket_card())\n        \n        layout.addStretch(1)\n        scroll_area.setWidget(content_widget)\n        \n        wrapper = QWidget()\n        wl = QVBoxLayout(wrapper); wl.setContentsMargins(0,0,0,0); wl.addWidget(scroll_area)\n        return wrapper\n\n    def create_hit_identification_tab(self) -> QWidget:\n        \"\"\"Hit Identification Tab: Vina Docking\"\"\"\n        scroll_area = QScrollArea()\n        scroll_area.setWidgetResizable(True)\n        scroll_area.setFrameShape(QFrame.Shape.NoFrame)\n        \n        content_widget = QWidget()\n        content_widget.setObjectName(\"scroll_content\")\n        bg_color = \"#0d1117\" if self._dark_mode else \"#f8fafc\"\n        content_widget.setStyleSheet(f\"#scroll_content {{ background-color: {bg_color}; }}\")\n        \n        layout = QVBoxLayout(content_widget)\n        layout.setSpacing(14)\n        layout.setContentsMargins(12, 12, 12, 12)\n        \n        # Vina Docking Card\n        layout.addWidget(self._create_vina_docking_card())\n        \n        # Placeholder for future Virtual Screening\n        grp_vs = QGroupBox(\"Virtual Screening (Coming Soon)\")\n        vs_layout = QVBoxLayout(grp_vs)\n        vs_layout.addWidget(QLabel(\"Batch docking and scoring functionality will be available in future updates.\"))\n        layout.addWidget(grp_vs)\n        \n        layout.addStretch(1)\n        scroll_area.setWidget(content_widget)\n        \n        wrapper = QWidget()\n        wl = QVBoxLayout(wrapper); wl.setContentsMargins(0,0,0,0); wl.addWidget(scroll_area)\n        return wrapper\n        \n    def create_lead_optimization_tab(self) -> QWidget:\n        \"\"\"Lead Optimization: PPI, Glue, Ternary, Mutation\"\"\"\n        scroll_area = QScrollArea()\n        scroll_area.setWidgetResizable(True)\n        scroll_area.setFrameShape(QFrame.Shape.NoFrame)\n        \n        content_widget = QWidget()\n        content_widget.setObjectName(\"scroll_content\")\n        self._lead_scroll_content = content_widget\n        bg_color = \"#0d1117\" if self._dark_mode else \"#f8fafc\"\n        content_widget.setStyleSheet(f\"#scroll_content {{ background-color: {bg_color}; }}\")\n        \n        layout = QVBoxLayout(content_widget)\n        layout.setSpacing(14)\n        layout.setContentsMargins(12, 12, 12, 12)\n        \n        # 1. Molecular Glue Specifics\n        grp_glue = QGroupBox(\"Molecular Glue Analysis (PPI + Neo-Epitope)\")\n        glue_grid = QGridLayout(grp_glue)\n        glue_grid.setColumnStretch(1, 1); glue_grid.setColumnStretch(3, 1)\n        glue_grid.setHorizontalSpacing(8); glue_grid.setVerticalSpacing(10)\n        \n        glue_grid.addWidget(QLabel(\"Target Object:\"), 0, 0, Qt.AlignmentFlag.AlignRight)\n        self.glue_obj_combo = QComboBox(); self.glue_obj_combo.setMinimumHeight(36)\n        self.glue_refresh_btn = QPushButton(t(\"refresh\")); self.glue_refresh_btn.clicked.connect(self.refresh_objects)\n        r0 = QHBoxLayout(); r0.addWidget(self.glue_obj_combo, 1); r0.addWidget(self.glue_refresh_btn)\n        glue_grid.addLayout(r0, 0, 1)\n        \n        glue_grid.addWidget(QLabel(\"Glue Resname:\"), 0, 2, Qt.AlignmentFlag.AlignRight)\n        self.glue_resname = QLineEdit(); self.glue_resname.setPlaceholderText(\"e.g. CC885\")\n        glue_grid.addWidget(self.glue_resname, 0, 3)\n        \n        glue_grid.addWidget(QLabel(\"E3 Chains:\"), 1, 0, Qt.AlignmentFlag.AlignRight)\n        self.glue_e3_chains = QLineEdit(); self.glue_e3_chains.setPlaceholderText(\"e.g. A\")\n        glue_grid.addWidget(self.glue_e3_chains, 1, 1)\n        \n        glue_grid.addWidget(QLabel(\"Substrate Chains:\"), 1, 2, Qt.AlignmentFlag.AlignRight)\n        self.glue_sub_chains = QLineEdit(); self.glue_sub_chains.setPlaceholderText(\"e.g. B\")\n        glue_grid.addWidget(self.glue_sub_chains, 1, 3)\n        \n        glue_grid.addWidget(QLabel(\"Interface Dist (Å):\"), 2, 0, Qt.AlignmentFlag.AlignRight)\n        self.glue_interface_dist = QLineEdit(\"4.5\")\n        glue_grid.addWidget(self.glue_interface_dist, 2, 1)\n        \n        glue_grid.addWidget(QLabel(\"Neo-Epitope Dist (Å):\"), 2, 2, Qt.AlignmentFlag.AlignRight)\n        self.glue_neo_dist = QLineEdit(\"5.0\")\n        glue_grid.addWidget(self.glue_neo_dist, 2, 3)\n        \n        glue_btn_row = QHBoxLayout()\n        self.glue_full_btn = QPushButton(\"Full Glue Analysis\"); self.glue_full_btn.setObjectName(\"highlight_btn\")\n        self.glue_full_btn.clicked.connect(self.run_glue_full_analysis)\n        glue_btn_row.addWidget(self.glue_full_btn)\n        glue_btn_row.addStretch(1)\n        \n        layout.addWidget(grp_glue)\n        layout.addLayout(glue_btn_row)\n        \n        # 2. Ternary Complex\n        grp_ternary = QGroupBox(\"Ternary Complex Analysis\")\n        t_grid = QGridLayout(grp_ternary)\n        t_grid.setColumnStretch(1, 1); t_grid.setColumnStretch(3, 1)\n        \n        t_grid.addWidget(QLabel(\"Target Object:\"), 0, 0)\n        self.tc_obj_combo = QComboBox(); \n        t_grid.addWidget(self.tc_obj_combo, 0, 1)\n        \n        t_grid.addWidget(QLabel(\"Ligand:\"), 0, 2)\n        self.tc_ligand_name = QLineEdit(); self.tc_ligand_name.setPlaceholderText(\"Auto\")\n        t_grid.addWidget(self.tc_ligand_name, 0, 3)\n        \n        t_grid.addWidget(QLabel(\"E3 Chains:\"), 1, 0)\n        self.tc_protein1_chains = QLineEdit(); self.tc_protein1_chains.setPlaceholderText(\"e.g. A\")\n        t_grid.addWidget(self.tc_protein1_chains, 1, 1)\n        \n        t_grid.addWidget(QLabel(\"POI Chains:\"), 1, 2)\n        self.tc_protein2_chains = QLineEdit(); self.tc_protein2_chains.setPlaceholderText(\"e.g. B\")\n        t_grid.addWidget(self.tc_protein2_chains, 1, 3)\n        \n        t_btn_row = QHBoxLayout()\n        self.tc_analyze_btn = QPushButton(\"Analyze Complex\"); self.tc_analyze_btn.setObjectName(\"highlight_btn\")\n        self.tc_analyze_btn.clicked.connect(self.run_tc_analysis)\n        self.tc_render_btn = QPushButton(\"Render All\"); self.tc_render_btn.setObjectName(\"highlight_btn\")\n        self.tc_render_btn.clicked.connect(self.run_tc_render)\n        t_btn_row.addWidget(self.tc_analyze_btn); t_btn_row.addWidget(self.tc_render_btn); t_btn_row.addStretch(1)\n        \n        layout.addWidget(grp_ternary)\n        layout.addLayout(t_btn_row)\n        \n        # 3. General Interactions (PP, PL)\n        grp_pl = QGroupBox(\"Protein-Ligand Interactions\")\n        pl_layout = QHBoxLayout(grp_pl)\n        self.pl_obj_combo = QComboBox(); self.pl_obj_combo.setMinimumWidth(150)\n        self.pl_analyze_btn = QPushButton(\"Analyze PL\"); self.pl_analyze_btn.clicked.connect(self.run_pl_analysis)\n        pl_layout.addWidget(QLabel(\"Object:\")); pl_layout.addWidget(self.pl_obj_combo)\n        pl_layout.addWidget(self.pl_analyze_btn)\n        layout.addWidget(grp_pl)\n        \n        # 4. Mutation Analysis\n        layout.addWidget(self._create_mutation_analysis_card())\n        \n        layout.addStretch(1)\n        scroll_area.setWidget(content_widget)\n        \n        wrapper = QWidget()\n        wl = QVBoxLayout(wrapper); wl.setContentsMargins(0,0,0,0); wl.addWidget(scroll_area)\n        return wrapper\n        \n    def create_visualization_tab(self) -> QWidget:\n        \"\"\"Visualization: APBS, Export\"\"\"\n        return self.create_apbs_tab()\n\n

    # ==============  蛋白-配体分析标签页 ==============
    def create_disease_analysis_tab(self) -> QWidget:
        """Create Disease Analysis Tab (V3)"""
        if DiseaseAnalysisTab:
            return DiseaseAnalysisTab()
        else:
            w = QWidget()
            layout = QVBoxLayout(w)
            layout.addWidget(QLabel("⚠️ Disease Analysis module not available."))
            layout.addWidget(QLabel("Please ensure 'open_targets_api.py' is present."))
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
        self.mut_obj_combo = QComboBox()
        self.mut_obj_combo.setFixedHeight(28)
        self.mut_refresh_btn = QPushButton("Refresh")
        self.mut_refresh_btn.setObjectName("refresh_btn")
        self.mut_refresh_btn.setFixedHeight(28)
        self.mut_refresh_btn.clicked.connect(self.refresh_objects)
        obj_row = QHBoxLayout()
        obj_row.addWidget(self.mut_obj_combo, 1)
        obj_row.addWidget(self.mut_refresh_btn)
        input_grid.addLayout(obj_row, 0, 1)
        
        # 突变输入
        input_grid.addWidget(QLabel("Mutations:"), 1, 0)
        self.mut_input = QLineEdit()
        self.mut_input.setPlaceholderText("e.g., A:23:ALA, B:45:GLY")
        self.mut_input.setFixedHeight(28)
        input_grid.addWidget(self.mut_input, 1, 1)
        
        # 方法选择 (FoldX only)
        input_grid.addWidget(QLabel("Method:"), 2, 0)
        method_row = QHBoxLayout()
        self.mut_method_combo = QComboBox()
        self.mut_method_combo.addItems(["FoldX"])
        self.mut_method_combo.setFixedHeight(28)
        self.mut_method_combo.setFixedWidth(120)
        self.mut_method_combo.setToolTip("ΔΔG calculation requires FoldX\nDownload: https://foldxsuite.crg.eu/")
        method_row.addWidget(self.mut_method_combo)
        method_row.addStretch()
        input_grid.addLayout(method_row, 2, 1)
        
        layout.addLayout(input_grid)
        
        # 按钮行
        btn_row = QHBoxLayout()
        btn_row.setSpacing(6)
        
        self.mut_perform_btn = QPushButton("Perform Mutation")
        self.mut_perform_btn.setObjectName("primary_btn")
        self.mut_perform_btn.setFixedHeight(32)
        self.mut_perform_btn.clicked.connect(self.run_mutation)
        
        self.mut_minimize_btn = QPushButton("Minimize Energy")
        self.mut_minimize_btn.setObjectName("highlight_btn")
        self.mut_minimize_btn.setFixedHeight(32)
        self.mut_minimize_btn.clicked.connect(self.run_minimize)
        
        self.mut_analyze_btn = QPushButton("Full Analysis")
        self.mut_analyze_btn.setObjectName("highlight_btn")
        self.mut_analyze_btn.setFixedHeight(32)
        self.mut_analyze_btn.clicked.connect(self.run_mutation_analysis)
        
        btn_row.addWidget(self.mut_perform_btn)
        btn_row.addWidget(self.mut_minimize_btn)
        btn_row.addWidget(self.mut_analyze_btn)
        btn_row.addStretch()
        
        layout.addLayout(btn_row)
        
        return card

    def run_mutation(self):
        """执行突变"""
        obj_name = self.mut_obj_combo.currentText()
        mutations_str = self.mut_input.text().strip()
        
        if not obj_name or not mutations_str:
            self.log("❌ 请选择对象并输入突变")
            return
        
        try:
            # 解析突变字符串
            mutations = []
            for mut in mutations_str.split(','):
                mut = mut.strip()
                if ':' in mut:
                    parts = mut.split(':')
                    if len(parts) == 3:
                        mutations.append((parts[0], parts[1], parts[2]))
            
            if not mutations:
                self.log("❌ 突变格式错误，示例: A:23:ALA, B:45:GLY")
                return
            
            self.log(f"🧬 执行突变: {obj_name}")
            
            # 调用突变分析模块
            from pymol import cmd
            cmd.perform_mutation(obj_name, mutations, method='pymol')
            
            self.log("✅ 突变完成")
            
        except Exception as e:
            self.log(f"❌ 突变失败: {e}")
            import traceback
            traceback.print_exc()

    def run_minimize(self):
        """能量最小化"""
        obj_name = self.mut_obj_combo.currentText()
        
        if not obj_name:
            self.log("❌ 请选择对象")
            return
        
        try:
            self.log(f"⚡ 能量最小化: {obj_name}")
            
            from pymol import cmd
            cmd.minimize_energy(obj_name, cycles=100)
            
            self.log("✅ 最小化完成")
            
        except Exception as e:
            self.log(f"❌ 最小化失败: {e}")
            import traceback
            traceback.print_exc()

    def run_mutation_analysis(self):
        """完整突变分析"""
        obj_name = self.mut_obj_combo.currentText()
        mutations_str = self.mut_input.text().strip()
        method = self.mut_method_combo.currentText().lower()
        
        if not obj_name or not mutations_str:
            self.log("❌ 请选择对象并输入突变")
            return
        
        try:
            # 解析突变
            mutations = []
            for mut in mutations_str.split(','):
                mut = mut.strip()
                if ':' in mut:
                    parts = mut.split(':')
                    if len(parts) == 3:
                        mutations.append((parts[0], parts[1], parts[2]))
            
            if not mutations:
                self.log("❌ 突变格式错误")
                return
            
            self.log(f"🧬 开始完整分析: {obj_name}")
            self.log(f"突变数量: {len(mutations)}")
            self.log(f"方法: {method}")
            
            # 调用完整分析
            from pymol import cmd
            cmd.analyze_mutation_effects(obj_name, mutations, method=method)
            
            self.log("✅ 分析完成")
            
        except Exception as e:
            self.log(f"❌ 分析失败: {e}")
            import traceback
            traceback.print_exc()
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setSpacing(12)
        layout.setContentsMargins(10, 10, 10, 10)

        grp = QGroupBox("蛋白-配体相互作用分析" if get_lang() == "zh" else "Protein-Ligand Interactions")
        form = QFormLayout(grp)
        form.setVerticalSpacing(10)
        form.setSpacing(10)

        # 对象选择
        obj_row = QHBoxLayout()
        self.pl_obj_combo = QComboBox()
        self.pl_obj_combo.setMinimumHeight(26)
        self.pl_refresh_btn = QPushButton(t("refresh"))
        self.pl_refresh_btn.setObjectName("refresh_btn")
        self.pl_refresh_btn.setMinimumHeight(26)
        self.pl_refresh_btn.clicked.connect(self.refresh_objects)
        obj_row.addWidget(self.pl_obj_combo, 1)
        obj_row.addWidget(self.pl_refresh_btn)
        form.addRow(QLabel(t("target_obj")), obj_row)

        # 配体残基名
        self.pl_ligand_name = QLineEdit()
        self.pl_ligand_name.setMinimumHeight(26)
        self.pl_ligand_name.setPlaceholderText("留空自动检测，例如: LIG, ATP" if get_lang() == "zh" else "Auto-detect if blank, e.g.: LIG, ATP")
        form.addRow(QLabel("配体残基名:" if get_lang() == "zh" else "Ligand Resname:"), self.pl_ligand_name)

        # 蛋白链
        self.pl_protein_chains = QLineEdit()
        self.pl_protein_chains.setMinimumHeight(26)
        self.pl_protein_chains.setPlaceholderText("留空自动检测，例如: A,B" if get_lang() == "zh" else "Auto-detect if blank, e.g.: A,B")
        form.addRow(QLabel("蛋白质链:" if get_lang() == "zh" else "Protein Chains:"), self.pl_protein_chains)

        # 距离截断
        self.pl_distance = QLineEdit("4.5")
        self.pl_distance.setMinimumHeight(26)
        form.addRow(QLabel("距离截断 (Å):" if get_lang() == "zh" else "Distance cutoff (Å):"), self.pl_distance)

        # 输出CSV
        csv_row = QHBoxLayout()
        self.pl_csv = QLineEdit()
        self.pl_csv.setMinimumHeight(26)
        self.pl_csv.setPlaceholderText("可选，留空不保存" if get_lang() == "zh" else "Optional")
        self.pl_csv_btn = QPushButton(t("browse"))
        self.pl_csv_btn.setObjectName("browse_btn")
        self.pl_csv_btn.setMinimumHeight(32)
        self.pl_csv_btn.clicked.connect(lambda: self._browse_save_file(self.pl_csv, "CSV (*.csv)"))
        pl_csv_row_container = QWidget()
        pl_csv_row = QHBoxLayout(pl_csv_row_container)
        pl_csv_row.setContentsMargins(0, 0, 0, 0)
        pl_csv_row.addWidget(self.pl_csv, 1)
        pl_csv_row.addWidget(self.pl_csv_btn)
        pl_grid.addWidget(pl_csv_row_container, 2, 1, 1, 3)  # 跨三列
        layout.addWidget(grp)
        
        # 提示信息
        info_label = QLabel()
        info_label.setText(
            "自动检测并使用最优分析模式（RDKit高级分析 或 基础分析）"
            if get_lang() == "zh" else
            "Automatically uses best mode (RDKit advanced or basic analysis)"
        )
        info_label.setStyleSheet("color: #6b7280; padding: 8px; font-size: 12px;")
        layout.addWidget(info_label)

        # 按钮行（垂直布局，更醒目）
        btn_grp = QGroupBox("操作步骤" if get_lang() == "zh" else "Operations")
        btn_layout = QVBoxLayout(btn_grp)
        btn_layout.setSpacing(10)
        
        # 分析按钮（大按钮）
        self.pl_analyze_btn = QPushButton("1. Analyze Interactions" if get_lang() == "en" else "1. 分析相互作用")
        self.pl_analyze_btn.setObjectName("highlight_btn")
        self.pl_analyze_btn.setMinimumHeight(45)
        self.pl_analyze_btn.clicked.connect(self.run_pl_analysis)
        btn_layout.addWidget(self.pl_analyze_btn)

        # 可视化和网络图按钮（水平排列）
        sub_btn_row = QHBoxLayout()
        
        self.pl_visualize_btn = QPushButton("2. 3D Visualize" if get_lang() == "en" else "2. PyMOL 3D可视化")
        self.pl_visualize_btn.setObjectName("highlight_btn")
        self.pl_visualize_btn.setMinimumHeight(40)
        self.pl_visualize_btn.clicked.connect(self.run_pl_visualize)

        self.pl_network_btn = QPushButton("3. Network Plot" if get_lang() == "en" else "3. 生成网络图")
        self.pl_network_btn.setObjectName("highlight_btn")
        self.pl_network_btn.setMinimumHeight(40)
        self.pl_network_btn.clicked.connect(self.run_pl_network)

        sub_btn_row.addWidget(self.pl_visualize_btn)
        sub_btn_row.addWidget(self.pl_network_btn)
        
        btn_layout.addLayout(sub_btn_row)
        
        layout.addWidget(btn_grp)
        layout.addStretch()
        
        return w

    def create_ternary_tab(self) -> QWidget:
        """创建三元复合体分析标签页"""
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setSpacing(12)
        layout.setContentsMargins(10, 10, 10, 10)

        grp = QGroupBox("三元复合体分析 (蛋白-配体-蛋白)" if get_lang() == "zh" else "Ternary Complex (Prot-Lig-Prot)")
        form = QFormLayout(grp)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        form.setSpacing(10)

        # 对象选择
        obj_row = QHBoxLayout()
        self.tc_obj_combo = QComboBox()
        self.tc_obj_combo.setMinimumHeight(26)
        self.tc_refresh_btn = QPushButton(t("refresh"))
        self.tc_refresh_btn.setObjectName("refresh_btn")
        self.tc_refresh_btn.setMinimumHeight(26)
        self.tc_refresh_btn.clicked.connect(self.refresh_objects)
        obj_row.addWidget(self.tc_obj_combo, 1)
        obj_row.addWidget(self.tc_refresh_btn)
        form.addRow(QLabel(t("target_obj")), obj_row)

        # 配体残基名
        self.tc_ligand_name = QLineEdit()
        self.tc_ligand_name.setMinimumHeight(26)
        self.tc_ligand_name.setPlaceholderText("PROTAC分子名称，留空自动检测" if get_lang() == "zh" else "PROTAC name, auto-detect if blank")
        form.addRow(QLabel("配体残基名:" if get_lang() == "zh" else "Ligand Resname:"), self.tc_ligand_name)

        # 蛋白1链
        self.tc_protein1_chains = QLineEdit()
        self.tc_protein1_chains.setMinimumHeight(26)
        self.tc_protein1_chains.setPlaceholderText("例如: A" if get_lang() == "zh" else "e.g.: A")
        form.addRow(QLabel("蛋白质1链:" if get_lang() == "zh" else "Protein1 Chains:"), self.tc_protein1_chains)

        # 蛋白2链
        self.tc_protein2_chains = QLineEdit()
        self.tc_protein2_chains.setMinimumHeight(26)
        self.tc_protein2_chains.setPlaceholderText("例如: B" if get_lang() == "zh" else "e.g.: B")
        form.addRow(QLabel("蛋白质2链:" if get_lang() == "zh" else "Protein2 Chains:"), self.tc_protein2_chains)

        # 距离截断
        self.tc_distance = QLineEdit("4.5")
        self.tc_distance.setMinimumHeight(26)
        form.addRow(QLabel("距离截断 (Å):" if get_lang() == "zh" else "Distance cutoff (Å):"), self.tc_distance)

        # 输出CSV
        csv_row = QHBoxLayout()
        self.tc_csv = QLineEdit()
        self.tc_csv.setMinimumHeight(26)
        self.tc_csv.setPlaceholderText("可选" if get_lang() == "zh" else "Optional")
        self.tc_csv_btn = QPushButton(t("browse"))
        self.tc_csv_btn.setObjectName("browse_btn")
        self.tc_csv_btn.setMinimumHeight(32)
        self.tc_csv_btn.clicked.connect(lambda: self._browse_save_file(self.tc_csv, "CSV (*.csv)"))
        tc_csv_row_container = QWidget()
        tc_csv_row = QHBoxLayout(tc_csv_row_container)
        tc_csv_row.setContentsMargins(0, 0, 0, 0)
        tc_csv_row.addWidget(self.tc_csv, 1)
        tc_csv_row.addWidget(self.tc_csv_btn)
        ternary_grid.addWidget(tc_csv_row_container, 2, 3)
        layout.addWidget(grp)

        # 按钮组（醒目）
        btn_grp = QGroupBox("操作步骤" if get_lang() == "zh" else "Operations")
        btn_layout = QVBoxLayout(btn_grp)
        btn_layout.setSpacing(10)
        
        self.tc_analyze_btn = QPushButton("1. Analyze Ternary" if get_lang() == "en" else "1. 分析三元复合体")
        self.tc_analyze_btn.setObjectName("highlight_btn")
        self.tc_analyze_btn.setMinimumHeight(45)
        self.tc_analyze_btn.clicked.connect(self.run_tc_analysis)

        self.tc_network_btn = QPushButton("2. Network Plot" if get_lang() == "en" else "2. 生成网络图")
        self.tc_network_btn.setObjectName("highlight_btn")
        self.tc_network_btn.setMinimumHeight(40)
        self.tc_network_btn.clicked.connect(self.run_tc_network)

        btn_layout.addWidget(self.tc_analyze_btn)
        btn_layout.addWidget(self.tc_network_btn)

        layout.addWidget(btn_grp)
        layout.addStretch()
        
        return w

    def create_atom_pair_tab(self) -> QWidget:
        """创建原子对分析标签页"""
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setSpacing(12)
        layout.setContentsMargins(10, 10, 10, 10)

        grp = QGroupBox("原子对相互作用分析 (精确到原子)" if get_lang() == "zh" else "Atom Pair Analysis")
        form = QFormLayout(grp)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
        form.setSpacing(10)

        # 对象选择
        obj_row = QHBoxLayout()
        self.ap_obj_combo = QComboBox()
        self.ap_obj_combo.setMinimumHeight(26)
        self.ap_refresh_btn = QPushButton(t("refresh"))
        self.ap_refresh_btn.setObjectName("refresh_btn")
        self.ap_refresh_btn.setMinimumHeight(26)
        self.ap_refresh_btn.clicked.connect(self.refresh_objects)
        obj_row.addWidget(self.ap_obj_combo, 1)
        obj_row.addWidget(self.ap_refresh_btn)
        form.addRow(QLabel(t("target_obj")), obj_row)

        # 原子1选择
        self.ap_atom1 = QLineEdit()
        self.ap_atom1.setMinimumHeight(26)
        self.ap_atom1.setPlaceholderText('例如: "resn LIG and name N1" 或 "LIG/301/N1"' if get_lang() == "zh" else 'e.g.: "resn LIG and name N1"')
        form.addRow(QLabel("原子1选择:" if get_lang() == "zh" else "Atom1 Selection:"), self.ap_atom1)

        # 原子2选择
        self.ap_atom2 = QLineEdit()
        self.ap_atom2.setMinimumHeight(26)
        self.ap_atom2.setPlaceholderText('例如: "elem O" 或 "SER/50/OG"' if get_lang() == "zh" else 'e.g.: "elem O"')
        form.addRow(QLabel("原子2选择:" if get_lang() == "zh" else "Atom2 Selection:"), self.ap_atom2)

        # 距离截断
        self.ap_distance = QLineEdit("5.0")
        self.ap_distance.setMinimumHeight(26)
        form.addRow(QLabel("距离截断 (Å):" if get_lang() == "zh" else "Distance cutoff (Å):"), self.ap_distance)

        # 输出CSV
        csv_row = QHBoxLayout()
        self.ap_csv = QLineEdit()
        self.ap_csv.setMinimumHeight(26)
        self.ap_csv.setPlaceholderText("可选" if get_lang() == "zh" else "Optional")
        self.ap_csv_btn = QPushButton(t("browse"))
        self.ap_csv_btn.setObjectName("browse_btn")
        self.ap_csv_btn.setMinimumHeight(32)
        self.ap_csv_btn.clicked.connect(lambda: self._browse_save_file(self.ap_csv, "CSV (*.csv)"))
        ap_csv_row_container = QWidget()
        ap_csv_row = QHBoxLayout(ap_csv_row_container)
        ap_csv_row.setContentsMargins(0, 0, 0, 0)
        ap_csv_row.addWidget(self.ap_csv, 1)
        ap_csv_row.addWidget(self.ap_csv_btn)
        ap_grid.addWidget(ap_csv_row_container, 2, 1, 1, 3)  # 跨三列
        layout.addWidget(grp)

        # 快速模板
        template_grp = QGroupBox("快速模板" if get_lang() == "zh" else "Quick Templates")
        template_layout = QHBoxLayout(template_grp)
        template_layout.setSpacing(8)
        
        templates = [
            ("N-O氢键" if get_lang() == "zh" else "N-O H-bonds", '"elem N"', '"elem O"', "3.5"),
            ("配体-SER" if get_lang() == "zh" else "Lig-SER", '"resn LIG"', '"resn SER"', "4.5"),
            ("二硫键" if get_lang() == "zh" else "S-S", '"name SG"', '"name SG"', "2.5"),
            ("金属配位" if get_lang() == "zh" else "Metal", '"resn ZN"', '"elem N or elem O or elem S"', "3.0"),
        ]

        for label, atom1, atom2, dist in templates:
            btn = QPushButton(label)
            btn.setObjectName("browse_btn")
            btn.setMinimumHeight(30)
            btn.clicked.connect(lambda checked, a1=atom1, a2=atom2, d=dist: self.apply_ap_template(a1, a2, d))
            template_layout.addWidget(btn)

        layout.addWidget(template_grp)

        # 按钮组
        btn_grp = QGroupBox("操作步骤" if get_lang() == "zh" else "Operations")
        btn_layout = QHBoxLayout(btn_grp)
        btn_layout.setSpacing(10)
        
        self.ap_analyze_btn = QPushButton("1. Analyze Pairs" if get_lang() == "en" else "1. 分析原子对")
        self.ap_analyze_btn.setObjectName("highlight_btn")
        self.ap_analyze_btn.setMinimumHeight(45)
        self.ap_analyze_btn.clicked.connect(self.run_ap_analysis)

        self.ap_visualize_btn = QPushButton("2. Visualize" if get_lang() == "en" else "2. PyMOL可视化")
        self.ap_visualize_btn.setObjectName("highlight_btn")
        self.ap_visualize_btn.setMinimumHeight(45)
        self.ap_visualize_btn.clicked.connect(self.run_ap_visualize)

        btn_layout.addWidget(self.ap_analyze_btn)
        btn_layout.addWidget(self.ap_visualize_btn)

        layout.addWidget(btn_grp)
        layout.addStretch()
        
        return w


    # ========== V3: Disease Analysis Tab ==========
    def create_disease_analysis_tab(self) -> QWidget:
        """创建疾病分析标签页（V3 功能）"""
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.addWidget(QLabel("Disease Analysis module has been removed."))
        layout.addStretch()
        return w

    def create_readme_tab(self) -> QWidget:
        """创建README页面"""
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setSpacing(16)
        layout.setContentsMargins(20, 20, 20, 20)
        
        # 顶部：Check Environment 移至此处
        top_row = QHBoxLayout()
        check_btn = QPushButton("Check Environment")
        check_btn.setObjectName("bottom_nav_btn")
        check_btn.setFlat(True)
        check_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        check_btn.clicked.connect(self.check_environment)
        top_row.addWidget(check_btn)
        top_row.addStretch(1)
        layout.addLayout(top_row)
        
        # 标题
        title = QLabel("GLINT - README" if get_lang() == "en" else "GLINT - 说明文档")
        title.setStyleSheet("""
            font-size: 20px;
            font-weight: bold;
            color: #3b82f6;
            padding: 10px;
        """)
        layout.addWidget(title)
        
        # 内容文本框
        text = QTextEdit()
        text.setReadOnly(True)
        
        if get_lang() == "zh":
            readme_content = """
<h2>欢迎使用 GLINT</h2>

<p><b>GLINT</b> 是一个专为分子胶降解剂分析设计的 PyMOL 插件。</p>

<h3>主要功能</h3>
<ul>
<li><b>G-Motif 识别</b> - 识别 CRBN G-loop 结构域</li>
<li><b>相互作用分析</b> - 分析蛋白-蛋白相互作用</li>
<li><b>蛋白-配体</b> - 高级蛋白-配体相互作用分析</li>
<li><b>三元复合体</b> - 分析三组分相互作用</li>
<li><b>原子对分析</b> - 精确到原子级别的分析</li>
<li><b>静电势</b> - APBS/Coulomb 静电势计算和可视化</li>
</ul>

<h3>快速上手</h3>
<ol>
<li>在 PyMOL 中加载你的结构文件</li>
<li>从左侧导航栏选择功能模块</li>
<li>配置参数并执行分析</li>
<li>在右侧查看结果和日志</li>
</ol>

<h3>提示</h3>
<ul>
<li>所有功能都支持 CSV 导出</li>
<li>结果可以直接在 PyMOL 中可视化</li>
<li>支持高分辨率图片导出 (300 DPI)</li>
<li>部分功能需要安装 RDKit 和 SciPy</li>
</ul>

<h3>文档</h3>
<p>查看插件目录中的以下文档：</p>
<ul>
<li><code>QUICK_START.md</code> - 快速开始指南</li>
<li><code>MODERN_UI_GUIDE.md</code> - 界面设计指南</li>
<li><code>README_MODERN_UI.md</code> - 现代化界面总结</li>
</ul>

<h3>主题</h3>
<p>当前使用的是<b>现代化深色主题</b>，也可切换到浅色主题。</p>
<p>编辑 <code>unified_gui.py</code> 的第 2030 行进行切换。</p>

<p style="margin-top: 20px; color: #64748b;">版本: 1.0 | 更新: 2025-11-05</p>
            """
        else:
            readme_content = """
<h2>Welcome to GLINT</h2>

<p><b>GLINT</b> is a PyMOL plugin designed for molecular glue degrader analysis and classification.</p>

<h3>Main Features</h3>
<ul>
<li><b>G-Motif Detection</b> - Identify CRBN G-loop structural motifs</li>
<li><b>Interaction Analysis</b> - Analyze protein-protein interactions</li>
<li><b>Protein-Ligand</b> - Advanced protein-ligand interaction analysis</li>
<li><b>Ternary Complex</b> - Analyze three-component interactions</li>
<li><b>Atom Pairs</b> - Atom-level precision analysis</li>
<li><b>Electrostatics</b> - APBS/Coulomb electrostatic calculation and visualization</li>
</ul>

<h3>Quick Start</h3>
<ol>
<li>Load your structure file in PyMOL</li>
<li>Select a function module from the left navigation</li>
<li>Configure parameters and run analysis</li>
<li>View results and logs on the right</li>
</ol>

<h3>Tips</h3>
<ul>
<li>All functions support CSV export</li>
<li>Results can be directly visualized in PyMOL</li>
<li>High-resolution image export supported (300 DPI)</li>
<li>Some features require RDKit and SciPy</li>
</ul>

<h3>Documentation</h3>
<p>Check these documents in the plugin directory:</p>
<ul>
<li><code>QUICK_START.md</code> - Quick start guide</li>
<li><code>MODERN_UI_GUIDE.md</code> - UI design guide</li>
<li><code>README_MODERN_UI.md</code> - Modern UI summary</li>
</ul>

<h3>Theme</h3>
<p>Currently using <b>Modern Dark Theme</b>, switchable to Light Theme.</p>
<p>Edit line 2030 in <code>unified_gui.py</code> to switch themes.</p>

<p style="margin-top: 20px; color: #64748b;">Version: 1.0 | Updated: 2025-11-05</p>
            """
        
        text.setHtml(readme_content)
        layout.addWidget(text)
        
        return w
    
    def create_contact_tab(self) -> QWidget:
        """创建 Contact Us 页面"""
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.setSpacing(16)
        layout.setContentsMargins(20, 20, 20, 20)
        
        # 标题
        title = QLabel("Contact Us" if get_lang() == "en" else "联系我们")
        title.setObjectName("contact_title")
        layout.addWidget(title)
        
        # 内容文本框
        text = QTextEdit()
        text.setReadOnly(True)
        text.setObjectName("contact_body")

        # --- Dynamic style for privacy box ---
        if self._dark_mode:
            privacy_style = "margin-top: 30px; padding: 15px; background: #2d3748; border: 1px solid #4a5568; border-radius: 8px; color: #a0aec0;"
            footer_color = "#718096"
        else:
            privacy_style = "margin-top: 30px; padding: 15px; background: #e7f3ff; border: 1px solid #93c5fd; border-radius: 8px; color: #1e40af;"
            footer_color = "#64748b"

        if get_lang() == "zh":
            contact_content = f'''
<h3>👋 感谢使用 GLINT！</h3>

<p>如果你有任何问题、建议或反馈，欢迎联系我们！</p>

<h3>📩 联系方式</h3>

<p style='font-size: 14px; line-height: 1.8;'>
<b>开发者：</b> Roufen Chen<br>
<b>Email：</b> <a href="mailto:12319021@zju.edu.cn">12319021@zju.edu.cn</a><br>
<b>GitHub：</b> <a href="https://github.com/VesperChen01/glue-pymol">https://github.com/VesperChen01/glue-pymol</a><br>
</p>

<h3>问题报告</h3>
<p>发现 Bug？请在 GitHub 上提交 Issue，并包含：</p>
<ul><li>PyMOL 版本</li><li>Python 版本</li><li>操作系统</li><li>错误信息和日志</li><li>复现步骤</li></ul>

<h3>功能建议</h3>
<p>有新功能想法？欢迎在 GitHub Discussions 中分享！</p>

<h3>贡献</h3>
<p>欢迎提交 Pull Request！请阅读 <code>CONTRIBUTING.md</code> 了解贡献指南。</p>

<h3>支持项目</h3>
<p>如果 GLINT 对你的研究有帮助，请考虑：</p><ul><li>在 GitHub 上给我们一个 Star ⭐</li><li>在论文中引用 GLINT</li><li>分享给同事</li></ul>

<p style='{privacy_style}'>
<b>隐私声明：</b>GLINT 不会收集任何个人数据或结构信息。所有分析都在本地进行。
</p>

<p style='margin-top: 20px; color: {footer_color}; text-align: center;'>
感谢你的支持！🚀
</p>'''
        else:
            contact_content = f'''
<h3>👋 Thank you for using GLINT!</h3>

<p>If you have any questions, suggestions, or feedback, please don't hesitate to contact us!</p>

<h3>📩 Contact Information</h3>

<p style='font-size: 14px; line-height: 1.8;'>
<b>Developer:</b> Roufen Chen<br>
<b>Email:</b> <a href="mailto:12319021@zju.edu.cn">12319021@zju.edu.cn</a><br>
<b>GitHub:</b> <a href="https://github.com/VesperChen01/glue-pymol">https://github.com/VesperChen01/glue-pymol</a><br>
</p>

<h3>Bug Reports</h3>
<p>Found a bug? Please submit an Issue on GitHub with:</p>
<ul><li>PyMOL version</li><li>Python version</li><li>Operating system</li><li>Error messages and logs</li><li>Steps to reproduce</li></ul>

<h3>Feature Requests</h3>
<p>Have ideas for new features? Share them in GitHub Discussions!</p>

<h3>Contributing</h3>
<p>Pull requests are welcome! Please read <code>CONTRIBUTING.md</code> for contribution guidelines.</p>

<h3>Support the Project</h3>
<p>If GLINT helped your research, please consider:</p>
<ul><li>Giving us a Star on GitHub</li><li>Citing GLINT in your papers</li><li>Sharing with colleagues</li></ul>

<p style='{privacy_style}'>
<b>Privacy:</b> GLINT does not collect any personal data or structural information. All analyses are performed locally.
</p>

<p style='margin-top: 20px; color: {footer_color}; text-align: center;'>
Thank you for your support! 🚀
</p>'''
        
        text.setHtml(contact_content)
        layout.addWidget(text)
        
        return w

    # ==============  新增：蛋白-配体分析功能实现 ==============
    def _browse_save_file(self, line_edit, file_filter):
        """通用文件保存浏览"""
        fn, _ = QFileDialog.getSaveFileName(self, "保存文件" if get_lang() == "zh" else "Save File", "", file_filter)
        if fn:
            line_edit.setText(fn)

    # def apply_ap_template(self, atom1, atom2, dist):
    #     """应用原子对模板"""
    #     self.ap_atom1.setText(atom1)
    #     self.ap_atom2.setText(atom2)
    #     self.ap_distance.setText(dist)
    #     msg = f"已应用模板: {atom1} - {atom2}" if get_lang() == "zh" else f"Template applied: {atom1} - {atom2}"
    #     self.log(msg)

    def _export_pymol_to_pdb_sdf(self, obj_name, ligand_resname=None):
        """
        将 PyMOL 对象导出为 PDB 和 SDF 文件供高级分析使用

        返回:
            (protein_pdb_path, ligand_sdf_path) 或 (None, None)
        """
        import tempfile
        from pymol import cmd

        try:
            # 创建临时文件
            protein_fd, protein_pdb = tempfile.mkstemp(suffix=".pdb", prefix="protein_")
            ligand_fd, ligand_sdf = tempfile.mkstemp(suffix=".sdf", prefix="ligand_")
            os.close(protein_fd)
            os.close(ligand_fd)

            # 如果指定了配体残基名，分别保存
            if ligand_resname:
                # 保存蛋白（排除配体）
                cmd.save(protein_pdb, f"{obj_name} and not resn {ligand_resname}")
                # 保存配体
                cmd.save(ligand_sdf, f"{obj_name} and resn {ligand_resname}")
            else:
                # 自动检测配体（非标准残基）
                # 标准残基列表
                standard_residues = {
                    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                    'HOH', 'WAT', 'A', 'C', 'G', 'T', 'U', 'DA', 'DC', 'DG', 'DT'
                }

                # 获取所有残基名
                model = cmd.get_model(obj_name)
                all_resnames = set(atom.resn for atom in model.atom)
                ligand_resnames = [r for r in all_resnames if r not in standard_residues]

                if not ligand_resnames:
                    self.log("未检测到配体残基" if get_lang() == "zh" else "No ligand residue detected")
                    return None, None

                ligand_resname = ligand_resnames[0]
                self.log(f"📌 自动检测到配体: {ligand_resname}")

                # 保存文件
                cmd.save(protein_pdb, f"{obj_name} and not resn {ligand_resname}")
                cmd.save(ligand_sdf, f"{obj_name} and resn {ligand_resname}")

            # 验证文件
            if not os.path.exists(protein_pdb) or not os.path.exists(ligand_sdf):
                self.log("导出文件失败" if get_lang() == "zh" else "Export failed")
                return None, None

            if os.path.getsize(protein_pdb) == 0 or os.path.getsize(ligand_sdf) == 0:
                self.log("导出文件为空" if get_lang() == "zh" else "Exported files are empty")
                return None, None

            self.log(f"已导出: {os.path.basename(protein_pdb)}, {os.path.basename(ligand_sdf)}")
            return protein_pdb, ligand_sdf

        except Exception as e:
            self.log(f"导出失败: {e}")
            import traceback
            traceback.print_exc()
            return None, None

    def _check_and_install_dependencies(self):
        """
        检查高级分析依赖，如果缺失则提示用户安装

        返回:
            True: 依赖满足
            False: 用户取消或安装失败
        """
        # 检查各个包
        missing = []

        try:
            import rdkit
        except ImportError:
            missing.append("rdkit-pypi")

        try:
            import scipy
        except ImportError:
            missing.append("scipy")

        try:
            import numpy
        except ImportError:
            missing.append("numpy")

        if not missing:
            return True  # 所有依赖都满足

        # 有缺失的包，询问用户
        packages_str = " ".join(missing)
        msg = (
            f"高级分析需要以下 Python 包：\n\n{', '.join(missing)}\n\n"
            f"是否现在自动安装？\n\n"
            f"（将在后台运行: pip install {packages_str}）"
            if get_lang() == "zh" else
            f"Advanced analysis requires:\n\n{', '.join(missing)}\n\n"
            f"Install now?\n\n"
            f"(Will run: pip install {packages_str})"
        )

        reply = QMessageBox.question(
            self,
            "安装依赖？" if get_lang() == "zh" else "Install Dependencies?",
            msg,
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.Yes
        )

        if reply != QMessageBox.StandardButton.Yes:
            return False

        # 用户同意，开始安装
        self.log(f"\n{'开始安装依赖包...' if get_lang() == 'zh' else 'Installing dependencies...'}")
        self.log(f"   包: {', '.join(missing)}")

        import subprocess
        import sys

        try:
            # 使用 PyMOL 的 Python 解释器
            python_exe = sys.executable

            # 构建安装命令
            cmd = [python_exe, "-m", "pip", "install"] + missing

            self.log(f"   命令: {' '.join(cmd)}")
            self.log("   ⏳ 正在安装，请稍候...")

            # 运行安装（显示输出）
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )

            # 实时显示输出
            for line in process.stdout:
                line = line.strip()
                if line:
                    self.log(f"     {line}")

            process.wait()

            if process.returncode == 0:
                self.log("\n依赖安装成功！")
                QMessageBox.information(
                    self,
                    "成功" if get_lang() == "zh" else "Success",
                    "依赖包已成功安装！\n请重新运行分析。" if get_lang() == "zh" else "Dependencies installed!\nPlease re-run the analysis."
                )
                return True
            else:
                self.log(f"\n安装失败，返回码: {process.returncode}")
                QMessageBox.critical(
                    self,
                    "安装失败" if get_lang() == "zh" else "Installation Failed",
                    f"安装失败。\n\n请手动运行:\npip install {packages_str}\n\n或联系管理员。"
                    if get_lang() == "zh" else
                    f"Installation failed.\n\nPlease run manually:\npip install {packages_str}\n\nOr contact admin."
                )
                return False

        except Exception as e:
            self.log(f"\n安装异常: {e}")
            import traceback
            traceback.print_exc()

            QMessageBox.critical(
                self,
                "安装错误" if get_lang() == "zh" else "Installation Error",
                f"自动安装失败: {e}\n\n请手动安装:\npip install {packages_str}"
                if get_lang() == "zh" else
                f"Auto-install failed: {e}\n\nPlease install manually:\npip install {packages_str}"
            )
            return False

    def run_pl_analysis(self):
        """运行蛋白-配体分析"""
        try:
            obj_name = self.pl_obj_combo.currentText()
            if obj_name == t("no_object"):
                QMessageBox.warning(self, t("title"), "请先加载PDB结构" if get_lang() == "zh" else "Load PDB first")
                return

            ligand_resname = self.pl_ligand_name.text().strip() or None
            output_csv = self.pl_csv.text().strip() or None

            # ========== 使用严格标准分析 ==========
            from .interaction_analyzer import analyze_protein_ligand_interactions

            protein_chains_str = self.pl_protein_chains.text().strip()
            protein_chains = [c.strip() for c in protein_chains_str.split(",")] if protein_chains_str else None
            
            try:
                distance_cutoff = float(self.pl_distance.text())
            except (ValueError, TypeError):  # float() 转换可能失败
                distance_cutoff = 4.5
                self.log(f"距离参数无效，使用默认值 4.5 Å")

            lig_txt = ligand_resname or ('自动' if get_lang() == 'zh' else 'auto')
            self.log(f"\n▶ {obj_name} | Lig: {lig_txt} | Cutoff: {distance_cutoff}Å")

            result = analyze_protein_ligand_interactions(
                obj_name=obj_name,
                ligand_resname=ligand_resname,
                protein_chains=protein_chains,
                distance_cutoff=distance_cutoff,
                output_csv=output_csv
            )

            if result:
                mode = result.get("mode", "strict")
                standard = result.get("standard", "Publication")
                mode_text = f"{standard}标准" if get_lang() == "zh" else f"{standard} Standard"
                
                # 处理结果
                if mode == "advanced":
                    # 高级分析模式（如果有定制advanced模块）
                    advanced_results = result.get("advanced_results", {})
                    total = sum(len(v) for v in advanced_results.values())
                    details = ', '.join([f"{k}:{len(v)}" for k,v in advanced_results.items() if v])
                    self.log(f"✓ {total} total ({details})")
                    
                    # 保存结果供可视化使用
                    self.current_pl_result_advanced = advanced_results
                    self.current_pl_result = result
                    self.current_pl_ligand_sdf = result.get("ligand_sdf", None)
                    
                    result_msg = f"分析完成 ({mode_text})\n总计: {total} 个相互作用" if get_lang() == "zh" else f"Done ({mode_text})\nTotal: {total} interactions"
                else:
                    # 严格标准模式
                    n_interactions = len(result["interactions"])
                    self.log(f"✓ {n_interactions} key interactions")
                    self.current_pl_result = result
                    if hasattr(self, 'current_pl_result_advanced'):
                        delattr(self, 'current_pl_result_advanced')
                    
                    result_msg = f"发现 {n_interactions} 个关键相互作用 ({mode_text})" if get_lang() == "zh" else f"Found {n_interactions} key interactions ({mode_text})"
                
                # 显示CSV表格
                if output_csv and os.path.exists(output_csv):
                    try:
                        self.fill_table_from_csv(output_csv)
                    except Exception as e:
                        self.log(f"✗ Table error: {e}")
                
                # 显示结果对话框
                QMessageBox.information(self, "完成" if get_lang() == "zh" else "Done", result_msg)
            else:
                self.log("未找到相互作用" if get_lang() == "zh" else "No interactions found")

        except Exception as e:
            self.log(f"错误: {e}")
            import traceback; traceback.print_exc()

    def run_pl_visualize(self):
        """运行3D可视化"""
        try:
            if not hasattr(self, 'current_pl_result') and not hasattr(self, 'current_pl_result_advanced'):
                QMessageBox.warning(self, "警告" if get_lang() == "zh" else "Warning",
                    "请先运行分析" if get_lang() == "zh" else "Run analysis first")
                return

            obj_name = self.pl_obj_combo.currentText()
            ligand_resname = self.pl_ligand_name.text().strip() or None
            
            # 获取CSV路径（如果有的话）
            csv_path = self.pl_csv.text().strip() or None
            
            # 获取可视化选项
            show_hydrophobic = self.pl_show_hydrophobic.isChecked()
            
            # 解析置信度阈值
            conf_text = self.pl_min_confidence.currentText()
            try:
                min_confidence = float(conf_text.split()[0])  # 提取数字部分，如"0.8 (default)"→0.8
            except (ValueError, TypeError):  # float() 转换可能失败
                min_confidence = 0.8  # 默认值

            # 静默生成,仅显示结果
            
            # 检查是否使用了高级分析模式
            if hasattr(self, 'current_pl_result_advanced'):
                # 高级分析模式 - 使用专用的可视化函数
                from .interaction_analyzer_advanced import visualize_advanced_interactions
                visualize_advanced_interactions(obj_name=obj_name, csv_path=csv_path, ligand_resname=ligand_resname)
                self.log("✓")
            else:
                # 标准分析模式
                from .interaction_analyzer import visualize_protein_ligand_3d
                if csv_path and os.path.exists(csv_path):
                    visualize_protein_ligand_3d(
                        obj_name=obj_name,
                        csv_path=csv_path,
                        ligand_resname=ligand_resname,
                        show_hydrophobic=show_hydrophobic,
                        min_confidence=min_confidence
                    )
                elif hasattr(self, 'current_pl_result'):
                    visualize_protein_ligand_3d(
                        obj_name=obj_name,
                        interactions_result=self.current_pl_result,
                        ligand_resname=ligand_resname,
                        show_hydrophobic=show_hydrophobic,
                        min_confidence=min_confidence
                    )
                else:
                    self.log("没有可用的相互作用数据" if get_lang() == "zh" else "No interaction data available")
                    return
                self.log("✓")

        except Exception as e:
            self.log(f"错误: {e}")
            import traceback; traceback.print_exc()

    def run_pl_network(self):
        """生成网络图"""
        try:
            from .interaction_analyzer import generate_interaction_network_plot

            if not hasattr(self, 'current_pl_result') and not hasattr(self, 'current_pl_result_advanced'):
                QMessageBox.warning(self, "警告" if get_lang() == "zh" else "Warning",
                    "请先运行分析" if get_lang() == "zh" else "Run analysis first")
                return

            fn, _ = QFileDialog.getSaveFileName(self, "保存网络图" if get_lang() == "zh" else "Save Network",
                "interaction_network.png", "PNG (*.png)")
            if not fn:
                return

            # 静默生成
            
            # 优先使用CSV文件（支持两种模式）
            csv_path = self.pl_csv.text().strip() or None
            ligand_sdf = getattr(self, 'current_pl_ligand_sdf', None)
            obj_name = self.pl_obj_combo.currentText()
            ligand_resname = self.pl_ligand_name.text().strip() or None
            
            # 获取绘图样式
            plot_style = self.pl_plot_style.currentText().split(" (")[0].lower()
            
            if csv_path and os.path.exists(csv_path):
                output_path = generate_interaction_network_plot(
                    csv_path=csv_path,
                    output_path=fn,
                    show_plot=False,
                    ligand_sdf=ligand_sdf,
                    obj_name=obj_name,
                    ligand_resname=ligand_resname,
                    plot_style=plot_style
                )
            elif hasattr(self, 'current_pl_result'):
                output_path = generate_interaction_network_plot(
                    interactions_result=self.current_pl_result,
                    output_path=fn,
                    show_plot=False,
                    ligand_sdf=ligand_sdf,
                    obj_name=obj_name,
                    ligand_resname=ligand_resname,
                    plot_style=plot_style
                )
            else:
                self.log("没有可用的数据生成网络图" if get_lang() == "zh" else "No data for network plot")
                return

            if output_path:
                self.log(f"✓ {os.path.basename(output_path)}")
            else:
                self.log(f"✗ Failed" if get_lang() == "zh" else "✗ Failed")

        except Exception as e:
            self.log(f"错误: {e}")
            import traceback; traceback.print_exc()

    def run_ll_analysis(self):
        """运行小分子-小分子相互作用分析"""
        if analyze_ligand_ligand_interactions is None:
             QMessageBox.critical(self, "Error", "Ligand-Ligand analysis module not found.")
             return

        # Check for PDB file first
        pdb_file = self.ll_pdb.text().strip()
        obj = self.ll_obj_combo.currentText()
        
        if pdb_file and os.path.exists(pdb_file):
            try:
                # Load PDB file
                from pymol import cmd
                loaded_obj = os.path.basename(pdb_file).split('.')[0]
                # Ensure unique name
                loaded_obj = cmd.get_unused_name(loaded_obj)
                cmd.load(pdb_file, loaded_obj)
                obj = loaded_obj
                self.log(f"Loaded {pdb_file} as {obj}")
                self.refresh_objects() # Refresh combos
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to load PDB file: {e}")
                return
        
        sel1 = self.ll_sel1.text().strip()
        sel2 = self.ll_sel2.text().strip()
        dist_str = self.ll_dist.text().strip()
        csv_path = self.ll_csv.text().strip() or None
        
        if not obj:
            QMessageBox.warning(self, "Missing Input", "Please select a target object or load a PDB file.")
            return
        if not sel1 or not sel2:
            QMessageBox.warning(self, "Missing Input", "Please define both Selection 1 and Selection 2.")
            return
            
        try:
            dist = float(dist_str) if dist_str else 4.5
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Distance must be a number.")
            return
            
        try:
            # Import pymol to catch specific exceptions
            import pymol
            
            # 构建完整的 PyMOL selection
            # 如果用户只输入了 resn LIG, 我们需要将其限制在 obj 内吗?
            # 通常用户输入 selection string, 最好结合 obj
            # E.g. "obj and (sel1)"
            full_sel1 = f"({obj}) and ({sel1})"
            full_sel2 = f"({obj}) and ({sel2})"
            
            interactions = analyze_ligand_ligand_interactions(
                obj_name=obj,
                sel1=full_sel1,
                sel2=full_sel2,
                cutoff=dist,
                output_csv=csv_path,
                visualize=True
            )
            
            msg = f"Analysis complete. Found {len(interactions)} interactions."
            if csv_path:
                msg += f"\nSaved to: {csv_path}"
            
            QMessageBox.information(self, "Success", msg)
            
        except pymol.CmdException as e:
            msg = str(e)
            if "Invalid selection name" in msg:
                 QMessageBox.critical(self, "Selection Error", 
                     f"PyMOL could not understand your selection.\n\n"
                     f"Error: {msg}\n\n"
                     f"Tip: Please use valid PyMOL selection syntax.\n"
                     f"Examples:\n"
                     f"• resn LIG (by residue name)\n"
                     f"• resi 900 (by residue index)\n"
                     f"• chain A (by chain)\n\n"
                     f"You entered: '{sel1}' and '{sel2}'")
            else:
                 QMessageBox.critical(self, "PyMOL Error", str(e))
            return
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))
            import traceback
            traceback.print_exc()

    def run_pn_analysis(self):
        obj = self.pn_obj_combo.currentText().strip()
        pdb_file = self.pn_pdb.text().strip() or None
        
        if (not obj or obj == t("no_object")) and not pdb_file:
            QMessageBox.warning(self, "Missing Input", "Please select a target object or load a PDB file.")
            return
            
        nucleic_chains = self.pn_nucleic_chains.text().strip()
        protein_chains = self.pn_protein_chains.text().strip()
        dist_str = self.pn_distance.text().strip()
        csv_path = self.pn_csv.text().strip() or None
        
        try:
            dist = float(dist_str)
        except ValueError:
            QMessageBox.warning(self, "Invalid Input", "Distance cutoff must be a number.")
            return
            
        n_chains = [c.strip() for c in nucleic_chains.split(",")] if nucleic_chains else None
        p_chains = [c.strip() for c in protein_chains.split(",")] if protein_chains else None
        
        self.pn_analyze_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        
        self.log(f"Starting Protein-Nucleic Acid analysis on '{obj}'...")
        
        self.pn_worker = PNAnalysisWorker(obj, n_chains, p_chains, csv_path, dist, pdb_file)
        self.pn_worker.progress.connect(self.log)
        self.pn_worker.error.connect(self.on_pn_error)
        self.pn_worker.finished.connect(self.on_finished_pn_analysis)
        self.pn_worker.start()
        
    def on_finished_pn_analysis(self, interactions):
        self.pn_analyze_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 1)
        
        self._last_pn_interactions = interactions
        
        msg = f"Analysis complete. Found {len(interactions)} protein-nucleic acid interactions."
        csv_path = self.pn_csv.text().strip()
        if csv_path:
            msg += f"\nResults saved to: {csv_path}"
        
        self.log(f"✅ {msg}")
        
        if interactions:
            reply = QMessageBox.question(self, "Success", 
                                         f"{msg}\n\nDo you want to generate an interaction network plot?",
                                         QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            if reply == QMessageBox.StandardButton.Yes:
                self.run_pn_network()
        else:
            QMessageBox.information(self, "Success", msg)

    def run_pn_network(self):
        if not hasattr(self, '_last_pn_interactions') or not self._last_pn_interactions:
            return
            
        try:
            from .interaction_analyzer import generate_interaction_network_plot
            
            # 让用户选择保存路径（默认文件名：pn_network.png）
            fn, _ = QFileDialog.getSaveFileName(self, "Save Network Plot", "pn_network.png", "PNG (*.png)")
            if fn:
                output_path = generate_interaction_network_plot(
                    interactions_result=self._last_pn_interactions,
                    output_path=fn,
                    show_plot=True,
                    plot_style="professional"  # 指定默认样式
                )
                if output_path:
                    self.log(f"✅ Network plot saved: {output_path}")
        except Exception as e:
            self.on_pn_error(f"Plotting failed: {e}")

    def run_pn_heatmap(self):
        if not hasattr(self, '_last_pn_interactions') or not self._last_pn_interactions:
            QMessageBox.warning(self, "Data Missing", "Please run analysis first.")
            return
            
        try:
            from .interaction_analyzer import generate_interaction_heatmap
            
            fn, _ = QFileDialog.getSaveFileName(self, "Save Heatmap", "interaction_heatmap.png", "PNG (*.png)")
            if fn:
                output_path = generate_interaction_heatmap(
                    interactions_result=self._last_pn_interactions,
                    output_path=fn,
                    show_plot=True
                )
                if output_path:
                    self.log(f"✅ Heatmap saved: {output_path}")
        except Exception as e:
            self.on_pn_error(f"Heatmap failed: {e}")
    
    def on_pn_error(self, msg: str):
        """Handle protein-nucleic acid analysis errors"""
        self.pn_analyze_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 1)
        self.log(f"❌ Protein-Nucleic Analysis Error: {msg}")
        QMessageBox.critical(self, "Analysis Error", f"Protein-Nucleic Acid analysis failed:\n\n{msg}")

    def run_tc_analysis(self):
        """运行三元复合体分析"""
        try:
            from .interaction_analyzer import analyze_ternary_complex

            obj_name = self.tc_obj_combo.currentText()
            if obj_name == t("no_object"):
                QMessageBox.warning(self, t("title"), "请先加载PDB结构" if get_lang() == "zh" else "Load PDB first")
                return

            ligand_resname = self.tc_ligand_name.text().strip() or None
            p1_chains_str = self.tc_protein1_chains.text().strip()
            p2_chains_str = self.tc_protein2_chains.text().strip()
            protein1_chains = [c.strip() for c in p1_chains_str.split(",")] if p1_chains_str else None
            protein2_chains = [c.strip() for c in p2_chains_str.split(",")] if p2_chains_str else None
            distance_cutoff = float(self.tc_distance.text())
            output_csv = self.tc_csv.text().strip() or None

            self.log(f"\n▶ Ternary: {obj_name}")

            result = analyze_ternary_complex(
                obj_name=obj_name,
                ligand_resname=ligand_resname,
                protein1_chains=protein1_chains,
                protein2_chains=protein2_chains,
                distance_cutoff=distance_cutoff,
                output_csv=output_csv
            )

            if result:
                stats = result["bridging_analysis"]
                self.log(f"✓ P1: {stats['protein1_interactions_count']} | P2: {stats['protein2_interactions_count']}")
                self.current_tc_result = result
                QMessageBox.information(self, "✓" if get_lang() == "zh" else "✓",
                    f"P1: {stats['protein1_interactions_count']}  P2: {stats['protein2_interactions_count']}" if get_lang() == "zh" else
                    f"P1: {stats['protein1_interactions_count']}  P2: {stats['protein2_interactions_count']}")
            else:
                self.log("分析失败" if get_lang() == "zh" else "Analysis failed")

        except Exception as e:
            self.log(f"错误: {e}")
            import traceback; traceback.print_exc()

    def run_tc_network(self):
        """生成三元复合体网络图"""
        try:
            from .interaction_analyzer import generate_interaction_network_plot

            if not hasattr(self, 'current_tc_result'):
                QMessageBox.warning(self, "警告" if get_lang() == "zh" else "Warning",
                    "请先运行三元复合体分析" if get_lang() == "zh" else "Run ternary analysis first")
                return

            fn, _ = QFileDialog.getSaveFileName(self, "保存网络图" if get_lang() == "zh" else "Save Network",
                "ternary_network.png", "PNG (*.png)")
            if not fn:
                return

            # 静默生成
            output_path = generate_interaction_network_plot(
                interactions_result=self.current_tc_result,
                output_path=fn,
                show_plot=False
            )

            if output_path:
                self.log(f"✓ {os.path.basename(output_path)}")

        except Exception as e:
            self.log(f"错误: {e}")

    def run_tc_interface(self):
        """分析三元复合体蛋白-蛋白界面"""
        try:
            from pymol import cmd
            
            obj_name = self.tc_obj_combo.currentText()
            if obj_name == t("no_object"):
                QMessageBox.warning(self, t("title"), "请先加载PDB结构" if get_lang() == "zh" else "Load PDB first")
                return

            p1_chains = self.tc_protein1_chains.text().strip()
            p2_chains = self.tc_protein2_chains.text().strip()
            
            if not p1_chains or not p2_chains:
                QMessageBox.warning(self, "警告" if get_lang() == "zh" else "Warning",
                    "请指定E3和POI链" if get_lang() == "zh" else "Please specify E3 and POI chains")
                return
                
            cutoff = float(self.tc_distance.text())
            
            self.log(f"\n▶ PPI: E3[{p1_chains}] - POI[{p2_chains}]")
            
            # 创建界面选择
            interface_sel = f"interface_{obj_name}"
            e3_sel = f"{obj_name} and chain {p1_chains}"
            poi_sel = f"{obj_name} and chain {p2_chains}"
            
            cmd.select(interface_sel, f"(byres ({e3_sel} within {cutoff} of {poi_sel})) or (byres ({poi_sel} within {cutoff} of {e3_sel}))")
            
            # 高亮界面
            cmd.hide("everything", obj_name)
            cmd.show("cartoon", obj_name)
            cmd.show("sticks", interface_sel)
            cmd.color("gray80", obj_name)
            cmd.color("cyan", f"{interface_sel} and chain {p1_chains}")
            cmd.color("orange", f"{interface_sel} and chain {p2_chains}")
            
            # 计算界面统计
            n_e3_residues = cmd.count_atoms(f"{interface_sel} and chain {p1_chains} and name CA")
            n_poi_residues = cmd.count_atoms(f"{interface_sel} and chain {p2_chains} and name CA")
            
            self.log(f"✓ E3: {n_e3_residues} res | POI: {n_poi_residues} res")
            
            # 如果选中了计算ΔΔG
            if hasattr(self, 'tc_include_ddg') and self.tc_include_ddg.isChecked():
                self.log(f"计算界面ΔΔG..." if get_lang() == "zh" else "Calculating interface ΔΔG...")
                # 这里可以调用FoldX或其他方法计算ΔΔG
                # 暂时使用简化的ASA方法估算
                import numpy as np
                asa = cmd.get_area(interface_sel)
                ddg_estimate = asa * 0.01  # 简化估算: ~0.01 kcal/mol per Å²
                self.log(f"  估算ΔΔG: {ddg_estimate:.2f} kcal/mol (基于ASA)")
                
        except Exception as e:
            self.log(f"错误: {e}")
            import traceback; traceback.print_exc()

    def run_tc_render(self):
        """渲染三元复合体的完整展示"""
        try:
            from pymol import cmd
            
            obj_name = self.tc_obj_combo.currentText()
            if obj_name == t("no_object"):
                QMessageBox.warning(self, t("title"), "请先加载PDB结构" if get_lang() == "zh" else "Load PDB first")
                return
                
            self.log(f"\n▶ Rendering {obj_name}")
            
            # 先运行分析（如果还没有）
            if not hasattr(self, 'current_tc_result'):
                self.run_tc_analysis()
                if not hasattr(self, 'current_tc_result'):
                    return
            
            p1_chains = self.tc_protein1_chains.text().strip()
            p2_chains = self.tc_protein2_chains.text().strip()
            ligand_resname = self.tc_ligand_name.text().strip()
            
            # 设置显示样式
            cmd.hide("everything", obj_name)
            cmd.show("cartoon", obj_name)
            
            # E3 ligase (CRBN/VHL) - 蓝色系
            if p1_chains:
                cmd.color("slate", f"{obj_name} and chain {p1_chains}")
                cmd.set("cartoon_transparency", 0.2, f"{obj_name} and chain {p1_chains}")
            
            # POI - 橙色系
            if p2_chains:
                cmd.color("wheat", f"{obj_name} and chain {p2_chains}")
                cmd.set("cartoon_transparency", 0.2, f"{obj_name} and chain {p2_chains}")
            
            # PROTAC/分子胶 - 绿色球棍
            if ligand_resname:
                ligand_sel = f"{obj_name} and resn {ligand_resname}"
                cmd.show("sticks", ligand_sel)
                cmd.show("spheres", ligand_sel)
                cmd.color("forest", ligand_sel)
                cmd.set("sphere_scale", 0.3, ligand_sel)
                cmd.set("stick_radius", 0.2, ligand_sel)
            
            # 显示界面残基
            if self.tc_analyze_interface.isChecked():
                cutoff = float(self.tc_distance.text())
                interface_sel = f"interface_{obj_name}_render"
                e3_sel = f"{obj_name} and chain {p1_chains}"
                poi_sel = f"{obj_name} and chain {p2_chains}"
                
                cmd.select(interface_sel, f"(byres ({e3_sel} within {cutoff} of {poi_sel})) or (byres ({poi_sel} within {cutoff} of {e3_sel}))")
                cmd.show("sticks", f"{interface_sel} and sidechain")
                cmd.set("stick_radius", 0.15, interface_sel)
            
            # 设置视角和光照
            cmd.orient(obj_name)
            cmd.zoom(obj_name, 5)
            cmd.set("ambient", 0.3)
            cmd.set("spec_power", 200)
            cmd.set("spec_reflect", 0.2)
            
            # 导出PNG
            output_path = os.path.join(os.getcwd(), f"ternary_complex_{obj_name}.png")
            cmd.png(output_path, width=2400, height=2400, dpi=300, ray=1)
            
            self.log(f"✓ {os.path.basename(output_path)}")
            
        except Exception as e:
            self.log(f"错误: {e}")
            import traceback; traceback.print_exc()
    
    # ========== 新增: Molecular Glue 特异分析函数 ==========
    def run_glue_ppi_analysis(self):
        """分析蛋白-蛋白界面 (PPI Interface)"""
        try:
            from .ppi_analyzer import analyze_protein_protein_interface
            
            obj_name = self.glue_obj_combo.currentText()
            if obj_name == t("no_object") or not obj_name:
                QMessageBox.warning(self, "Warning", "Please select a structure object")
                return
            
            e3_chains_str = self.glue_e3_chains.text().strip()
            sub_chains_str = self.glue_sub_chains.text().strip()
            
            if not e3_chains_str or not sub_chains_str:
                QMessageBox.warning(self, "Warning", "Please specify both E3 and Substrate chains")
                return
            
            e3_chains = [c.strip() for c in e3_chains_str.split(",")]
            sub_chains = [c.strip() for c in sub_chains_str.split(",")]
            interface_dist = float(self.glue_interface_dist.text())
            ppi_csv = self.glue_ppi_csv.text().strip() or None
            
            self.log(f"\n▶ PPI: E3{e3_chains} - Sub{sub_chains}")
            
            result = analyze_protein_protein_interface(
                obj_name=obj_name,
                protein1_chains=e3_chains,
                protein2_chains=sub_chains,
                interface_distance=interface_dist,
                output_csv=ppi_csv
            )
            
            if result:
                # 保存结果
                self.current_glue_ppi_result = result
                
                # 显示结果
                contacts = result.get('interface_contacts', 0)
                bsa = result.get('bsa')
                is_strong = result.get('is_strong_interface', False)
                strength = result.get('interface_strength', 0)
                
                self.log(f"\nPPI Analysis Complete:")
                self.log(f"   Interface Contacts: {contacts}")
                if bsa:
                    self.log(f"   BSA: {bsa:.1f} Ų")
                self.log(f"   Interface Strength: {strength:.1f}/10")
                self.log(f"   Classification: {'Strong Interface' if is_strong else 'Weak Interface'}")
                
                # 判断机制
                if is_strong or contacts >= 10:
                    self.log(f"   Likely: Molecular Glue (strong PPI)")
                elif contacts < 3:
                    self.log(f"   Likely: PROTAC (weak PPI)")
                
                QMessageBox.information(self, "PPI Analysis Complete",
                    f"Interface Contacts: {contacts}\n"
                    f"{'BSA: ' + str(round(bsa, 1)) + ' Ų' if bsa else 'BSA: N/A'}\n"
                    f"Interface Strength: {strength:.1f}/10\n\n"
                    f"{'Strong Interface (likely Molecular Glue)' if is_strong else 'Weak Interface (likely PROTAC)'}")
            else:
                self.log("PPI analysis failed")
        
        except Exception as e:
            self.on_error(str(e))
            import traceback; traceback.print_exc()
    
    def run_glue_neo_epitope(self):
        """识别 Neo-表位 (Neo-Substrate Epitope)"""
        try:
            from .ppi_analyzer import identify_neo_epitope
            
            obj_name = self.glue_obj_combo.currentText()
            if obj_name == t("no_object") or not obj_name:
                QMessageBox.warning(self, "Warning", "Please select a structure object")
                return
            
            glue_resname = self.glue_resname.text().strip()
            e3_chains_str = self.glue_e3_chains.text().strip()
            sub_chains_str = self.glue_sub_chains.text().strip()
            
            if not glue_resname or not e3_chains_str or not sub_chains_str:
                QMessageBox.warning(self, "Warning", "Please specify Glue residue name, E3 and Substrate chains")
                return
            
            e3_chains = [c.strip() for c in e3_chains_str.split(",")]
            sub_chains = [c.strip() for c in sub_chains_str.split(",")]
            neo_dist = float(self.glue_neo_dist.text())
            neo_csv = self.glue_neo_csv.text().strip() or None
            
            self.log(f"\n▶ Neo-Epitope: Glue={glue_resname}")
            self.log(f"   E3: {e3_chains} → Substrate: {sub_chains}")
            
            result = identify_neo_epitope(
                obj_name=obj_name,
                e3_ligase_chains=e3_chains,
                substrate_chains=sub_chains,
                glue_resname=glue_resname,
                distance_threshold=neo_dist,
                output_csv=neo_csv
            )
            
            if result:
                # 保存结果
                self.current_glue_neo_result = result
                
                # 显示结果
                neo_count = result.get('neo_epitope_count', 0)
                bridging_atoms = result.get('bridging_glue_atoms', 0)
                is_glue = result.get('is_molecular_glue', False)
                confidence = result.get('confidence', 0)
                
                self.log(f"\nNeo-Epitope Detection Complete:")
                self.log(f"   Bridging Glue Atoms: {bridging_atoms}")
                self.log(f"   Neo-Epitope Residues: {neo_count}")
                self.log(f"   Confidence: {confidence:.2f}")
                self.log(f"   Classification: {'Molecular Glue' if is_glue else 'Not Typical Glue'}")
                
                # 列出Neo-表位残基
                if result.get('neo_substrate_residues'):
                    self.log(f"\n   Neo-Epitope Residues:")
                    for res in result['neo_substrate_residues'][:10]:  # 显示前10个
                        self.log(f"      {res['chain']}:{res['resname']} {res['resid']}")
                
                QMessageBox.information(self, "Neo-Epitope Detection Complete",
                    f"Bridging Glue Atoms: {bridging_atoms}\n"
                    f"Neo-Epitope Residues: {neo_count}\n"
                    f"Confidence: {confidence:.2f}\n\n"
                    f"{'Classified as Molecular Glue' if is_glue else 'Not typical Glue mechanism'}")
            else:
                self.log("Neo-epitope detection failed")
        
        except Exception as e:
            self.on_error(str(e))
            import traceback; traceback.print_exc()
    
    def run_glue_full_analysis(self):
        """完整分子胶分析: PPI + Neo-Epitope + Scoring + Classification"""
        try:
            obj_name = self.glue_obj_combo.currentText()
            if obj_name == t("no_object") or not obj_name:
                QMessageBox.warning(self, "Warning", "Please select a structure object")
                return
            
            glue_resname = self.glue_resname.text().strip()
            e3_chains_str = self.glue_e3_chains.text().strip()
            sub_chains_str = self.glue_sub_chains.text().strip()
            
            if not glue_resname or not e3_chains_str or not sub_chains_str:
                QMessageBox.warning(self, "Warning", "Please specify Glue residue name, E3 and Substrate chains")
                return
            
            self.log("\n" + "="*60)
            self.log("Starting Full Molecular Glue Analysis")
            self.log("="*60)
            
            # Step 1: PPI Analysis
            self.log("\n▶ [1/3] PPI Analysis...")
            self.run_glue_ppi_analysis()
            
            if not hasattr(self, 'current_glue_ppi_result'):
                self.log("PPI analysis failed, aborting")
                return
            
            # Step 2: Neo-Epitope Detection
            self.log("\n▶ [2/3] Neo-Epitope...")
            self.run_glue_neo_epitope()
            
            if not hasattr(self, 'current_glue_neo_result'):
                self.log("Neo-epitope detection failed, continuing...")
            
            # Step 3: Integrated Analysis & Classification
            self.log("\n[3/3] Integrated Mechanism Classification...")
            
            ppi_result = self.current_glue_ppi_result
            neo_result = getattr(self, 'current_glue_neo_result', None)
            
            # 提取关键指标
            ppi_contacts = ppi_result.get('interface_contacts', 0)
            bsa = ppi_result.get('bsa')
            is_strong_ppi = ppi_result.get('is_strong_interface', False)
            
            neo_count = neo_result.get('neo_epitope_count', 0) if neo_result else 0
            is_glue_neo = neo_result.get('is_molecular_glue', False) if neo_result else False
            
            # 综合判断
            final_mechanism = "Unknown"
            confidence = 0.0
            
            if is_glue_neo and is_strong_ppi:
                final_mechanism = "Molecular Glue"
                confidence = 0.95
            elif is_strong_ppi and neo_count >= 3:
                final_mechanism = "Molecular Glue"
                confidence = 0.85
            elif ppi_contacts >= 10 or (bsa and bsa > 800):
                final_mechanism = "Likely Molecular Glue"
                confidence = 0.75
            elif ppi_contacts < 3 and neo_count == 0:
                final_mechanism = "PROTAC (Linker-based)"
                confidence = 0.80
            else:
                final_mechanism = "Uncertain"
                confidence = 0.50
            
            # 输出最终报告
            self.log("\n" + "="*60)
            self.log("FINAL CLASSIFICATION")
            self.log("="*60)
            self.log(f"   Mechanism: {final_mechanism}")
            self.log(f"   Confidence: {confidence:.0%}")
            self.log(f"\n   Key Evidence:")
            self.log(f"      PPI Contacts: {ppi_contacts}")
            if bsa:
                self.log(f"      BSA: {bsa:.1f} Ų")
            self.log(f"      Neo-Epitope: {neo_count} residues")
            self.log(f"      Strong PPI: {'Yes' if is_strong_ppi else 'No'}")
            self.log("="*60)
            
            # 弹窗显示
            icon = QMessageBox.Icon.Information if "Glue" in final_mechanism else QMessageBox.Icon.Warning
            glue_msg = 'This complex exhibits Molecular Glue characteristics!'
            protac_msg = 'This complex likely uses a PROTAC/linker mechanism.'
            final_msg = glue_msg if 'Glue' in final_mechanism else protac_msg
            QMessageBox.information(self, "Molecular Glue Analysis Complete",
                f"Classification: {final_mechanism}\n"
                f"Confidence: {confidence:.0%}\n\n"
                f"Evidence:\n"
                f"  • PPI Contacts: {ppi_contacts}\n"
                f"  • BSA: {round(bsa, 1) if bsa else 'N/A'} Ų\n"
                f"  • Neo-Epitope: {neo_count} residues\n\n"
                f"{final_msg}")
        
        except Exception as e:
            self.on_error(str(e))
            import traceback; traceback.print_exc()

    # DISABLED - Atom Pair Analysis functions (not needed for general users)
    # def run_ap_analysis(self):
    #     """运行原子对分析"""
    #     pass
    # 
    # def run_ap_visualize(self):
    #     """可视化原子对"""
    #     pass

    # --- 杂项 ---
    def render_interactions_beautifully_clicked(self):
        """相互作用一键渲染：分析 → 高亮 → 美化渲染 → PNG"""
        try:
            obj = self.obj_combo_analysis.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object")); return

            # 步骤1: 如果没有CSV文件或相互作用数据，先进行分析
            csv_path = self.out_csv.text().strip()
            has_csv = csv_path and os.path.exists(csv_path)

            if not has_csv or not self._interactions:
                self.log("▶ Analyzing interactions...")
                pdb = self.pdb_path.text().strip() or None
                only_between = True  # 默认只分析链间相互作用

                # 直接调用分析函数（不使用线程，避免UI卡顿）
                interactions = analyze_pdb_interactions(
                    obj_name=obj,
                    pdb_file=pdb,
                    only_between_chains=only_between,
                    output_csv=csv_path if csv_path else None,
                    auto_highlight=False  # 不自动高亮，我们手动控制
                )

                self._interactions = interactions
                if csv_path:
                    has_csv = True
                    self.log(f"✓ {len(interactions)} interactions")
                    # 填充表格
                    self.fill_table_from_interactions(interactions)
                else:
                    # 创建临时CSV
                    import tempfile
                    fd, csv_path = tempfile.mkstemp(suffix="_interactions.csv")
                    os.close(fd)
                    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
                        import csv as csv_module
                        writer = csv_module.writer(f)
                        writer.writerow(["Chain1", "Residue1", "Chain2", "Residue2", "Distance", "Interaction"])
                        for inter in interactions:
                            writer.writerow([
                                inter["Chain1"], inter["Residue1"],
                                inter["Chain2"], inter["Residue2"],
                                inter["Distance"], inter["Interaction"]
                            ])
                    has_csv = True
                    self.log(f"✓ {len(interactions)} interactions")
                    self.fill_table_from_csv(csv_path)

            # 步骤2: 先高亮相互作用残基
            if has_csv:
                # 静默高亮
                try:
                    highlight_csv_residues(csv_path, obj=obj, show_labels=1,
                                         stick_by_element=1, clear_old=1)
                    # 静默完成
                except Exception as e:
                    self.log(f"高亮失败: {e}")

            # 步骤3: 美化渲染（不再重复高亮）
            # 静默渲染
            render_interactions_beautifully(obj, csv_path=None)  # 不传csv_path，避免重复高亮
            self.log("✓ Done")

            # 步骤4: 导出 PNG
            self._export_png_for_object(obj)

        except Exception as e:
            self.on_error(str(e))

    def load_csv_to_table(self):
        path = self.csv_path.text().strip()
        if path and os.path.exists(path):
            self.fill_table_from_csv(path)
        else:
            QMessageBox.warning(self, t("title"), t("select_csv"))

    def _browse_open_file(self, line_edit, file_filter):
        """通用文件打开浏览"""
        fn, _ = QFileDialog.getOpenFileName(self, "选择文件" if get_lang() == "zh" else "Select File", "", file_filter)
        if fn:
            line_edit.setText(fn)
    
    # Note: CRBN interface tools have been integrated into Ternary Complex analysis
    
    def check_environment(self):
        """检查环境和依赖状态"""
        try:
            from .env_checker import _quick_check_deps
            
            # 快速检查依赖
            all_ok, missing = _quick_check_deps()
            
            if all_ok:
                QMessageBox.information(
                    self,
                    "Environment OK",
                    "✅ All dependencies are installed!\n\nGLINT is ready to use."
                )
                self.log("✅ Environment check passed - all dependencies installed")
                return
            
            # 缺少依赖，提供安装选项
            reply = QMessageBox.question(
                self,
                "Missing Dependencies",
                f"❌ Missing: {', '.join(missing)}\n\n"
                "Would you like to run the installation script in Terminal?\n\n"
                "This will:\n"
                "  1. Create conda environment 'glint'\n"
                "  2. Install all dependencies\n"
                "  3. After completion, restart PyMOL with:\n"
                "     conda activate glint && pymol",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            
            if reply == QMessageBox.StandardButton.Yes:
                self._run_install_script()
            else:
                self.log(f"⚠️  Missing dependencies: {', '.join(missing)}")
                self.log("   Run install.sh manually or install with conda")
                
        except Exception as e:
            self.log(f"❌ Environment check failed: {e}")
            import traceback
            traceback.print_exc()
    
    def _run_install_script(self):
        """在终端中运行安装脚本"""
        import subprocess
        import sys
        
        # 多种方式查找 install.sh
        candidates = []
        
        # 方法1: 从 __file__ 推断
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            candidates.append(os.path.join(os.path.dirname(script_dir), "install.sh"))
            candidates.append(os.path.join(script_dir, "install.sh"))
        except Exception:  # __file__ 可能未定义
            pass
        
        # 方法2: 从 glint 模块位置推断
        try:
            import glint
            glint_dir = os.path.dirname(os.path.abspath(glint.__file__))
            candidates.append(os.path.join(os.path.dirname(glint_dir), "install.sh"))
        except (ImportError, AttributeError):  # 模块导入或属性访问可能失败
            pass
        
        # 方法3: 用户目录下的常见位置
        candidates.append(os.path.expanduser("~/git/GLINT/install.sh"))
        
        # 查找存在的脚本
        install_script = None
        for path in candidates:
            if os.path.exists(path):
                install_script = path
                break
        
        # 生成内联安装命令（当找不到脚本时使用）
        inline_cmd = '''echo "============================================================"
echo "🧬 GLINT - Installing Dependencies"
echo "============================================================"
conda create -n glint python=3.9 -y
conda install -n glint -c conda-forge rdkit scipy matplotlib pillow numpy pandas seaborn pyqt openbabel pymol-open-source boost boost-cpp swig -y
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate glint
pip install vina meeko || echo "Vina install failed (optional)"
echo ""
echo "✅ Installation complete!"
echo ""
echo "🚀 To start GLINT:"
echo "   conda activate glint && pymol"
echo "============================================================"'''
        
        if install_script:
            self.log(f"🚀 Opening Terminal to run: {install_script}")
            run_cmd = f"bash '{install_script}'"
        else:
            self.log("🚀 Opening Terminal to run inline installation...")
            run_cmd = inline_cmd
        
        try:
            if sys.platform == "darwin":  # macOS
                # 使用 osascript 打开 Terminal 并运行脚本
                apple_script = f'''
                tell application "Terminal"
                    activate
                    do script "{run_cmd}"
                end tell
                '''
                subprocess.run(["osascript", "-e", apple_script])
                
            elif sys.platform == "linux":
                # Linux: 尝试常见的终端
                terminals = ["gnome-terminal", "xterm", "konsole"]
                for term in terminals:
                    try:
                        if term == "gnome-terminal":
                            subprocess.Popen([term, "--", "bash", "-c", run_cmd])
                        else:
                            subprocess.Popen([term, "-e", f"bash -c '{run_cmd}'"])
                        break
                    except FileNotFoundError:
                        continue
                        
            elif sys.platform == "win32":  # Windows
                # Windows: 写入临时批处理文件
                import tempfile
                bat_file = os.path.join(tempfile.gettempdir(), "glint_install.bat")
                with open(bat_file, "w") as f:
                    f.write(run_cmd.replace("echo", "@echo"))
                subprocess.Popen(["cmd", "/c", "start", "cmd", "/k", bat_file])
            
            QMessageBox.information(
                self,
                "Installation Started",
                "Installation is running in Terminal.\n\n"
                "After it completes:\n"
                "  1. Close this PyMOL\n"
                "  2. Run: conda activate glint && pymol"
            )
            
        except Exception as e:
            QMessageBox.warning(
                self,
                "Error",
                f"Failed to open terminal: {e}\n\n"
                "Please run manually in terminal:\n"
                "  conda create -n glint python=3.9 -y\n"
                "  conda install -n glint -c conda-forge rdkit scipy "
                "matplotlib pillow numpy pandas seaborn pyqt autodock-vina openbabel -y"
            )
    
    def _check_environment_simple(self):
        """简单版环境检查（备用）"""
        self.log(f"\n{'='*50}")
        self.log("检查环境和依赖 (简单模式)" if get_lang() == "zh" else "Checking Environment (Simple Mode)")
        self.log(f"{'='*50}\n")
        
        # 检查 Python 依赖
        try:
            from .env_setup import get_dependency_status
            status = get_dependency_status()
            
            self.log("Python 依赖:" if get_lang() == "zh" else "Python Dependencies:")
            for pkg, available in status.items():
                symbol = "[OK]" if available else "[MISS]" 
                self.log(f"  {symbol} {pkg:20} {'已安装' if available else '未安装'}")
            
            missing = [pkg for pkg, avail in status.items() if not avail]
            if missing:
                self.log(f"\n缺失依赖: {', '.join(missing)}")
                self.log("安装命令: pip install " + " ".join(missing))
            else:
                self.log("\n所有 Python 依赖已满足" if get_lang() == "zh" else "\nAll Python dependencies satisfied")
        except Exception as e:
            self.log(f"无法检查 Python 依赖: {e}")
        
        self.log(f"\n{'='*50}")
        self.log("环境检查完成" if get_lang() == "zh" else "Environment check complete")
        self.log(f"{'='*50}\n")
    
    # ==========================
    # Scoring Callbacks
    # ==========================
    def run_quick_scoring(self):
        """运行快速经验评分 - DEPRECATED"""
        QMessageBox.information(
            self, 
            "Feature Removed", 
            "Empirical scoring has been removed.\n\n"
            "Please use:\n"
            "• Interaction Analysis tab for interaction analysis\n"
            "• Docking tab for AutoDock Vina scoring"
        )
        self.log("⚠️ Empirical scoring feature has been removed")
    
    def run_vina_scoring(self):
        """运行Vina评分 - DEPRECATED"""
        QMessageBox.information(
            self,
            "Feature Moved",
            "Vina scoring has been integrated into the Docking tab.\n\n"
            "Please use:\n"
            "• Docking tab → Pocket-Based Docking for full workflow"
        )
        self.log("⚠️ Please use Docking tab for Vina scoring")
    
    # ==========================
    # Pocket & Docking Callbacks
    # ==========================
    def run_pocket_detection(self):
        """运行口袋检测"""
        try:
            obj = self.pocket_obj_combo.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select an object")
                return
            
            # Get parameters
            try:
                grid_spacing = float(self.pocket_grid_spacing.text().strip() or "0.6")
                min_volume = float(self.pocket_min_volume.text().strip() or "20")
            except ValueError:
                QMessageBox.warning(self, "Warning", "Invalid parameter values")
                return
            
            self.log(f"\n🔍 Detecting pockets in {obj}...")
            self.log(f"   Grid spacing: {grid_spacing} Å")
            self.log(f"   Min volume: {min_volume} Ų")
            
            # Import pocket detector
            try:
                from .pocket_detector import detect_pockets
            except ImportError:
                from pocket_detector import detect_pockets
            
            # Detect pockets
            pockets = detect_pockets(
                obj_name=obj,
                grid_spacing=grid_spacing,
                min_volume=min_volume
            )
            
            if not pockets:
                self.log("⚠️  No pockets detected")
                QMessageBox.information(self, "Result", "No pockets detected with current parameters.\n\nTry adjusting grid spacing or min volume.")
                return
            
            # Store pockets for later use
            self._detected_pockets = pockets
            
            # Log results
            self.log(f"\n✅ Detected {len(pockets)} pockets:")
            for i, pocket in enumerate(pockets[:5], 1):  # Show top 5
                self.log(f"   {i}. Volume: {pocket['volume']:.1f} Ų, "
                        f"Druggability: {pocket['druggability_score']:.2f}, "
                        f"Center: ({pocket['center'][0]:.1f}, {pocket['center'][1]:.1f}, {pocket['center'][2]:.1f})")
            
            if len(pockets) > 5:
                self.log(f"   ... and {len(pockets) - 5} more")
            
            # Auto-visualize
            self.run_pocket_visualization()
            
        except Exception as e:
            self.on_error(str(e))
            import traceback
            traceback.print_exc()
    
    def run_pocket_visualization(self):
        """可视化检测到的口袋"""
        try:
            if not hasattr(self, '_detected_pockets') or not self._detected_pockets:
                QMessageBox.warning(self, "Warning", "Please run pocket detection first")
                return
            
            color_by_map = {
                "Volume": "volume",
                "Druggability": "druggability",
                "Hydrophobicity": "hydrophobicity",
                "Depth": "depth"
            }
            color_by = color_by_map.get(self.pocket_color_by.currentText(), "volume")
            
            self.log(f"\nVisualizing pockets (colored by {color_by})...")
            
            # Import visualizer
            try:
                from .pocket_visualizer import visualize_pockets
            except ImportError:
                from pocket_visualizer import visualize_pockets
            
            # Visualize
            visualize_pockets(
                self._detected_pockets,
                obj_name='glint_pockets',
                color_by=color_by,
                show_spheres=True,
                sphere_radius=1.5
            )
            
            self.log("Pockets visualized")
            
        except Exception as e:
            self.on_error(str(e))
            import traceback
            traceback.print_exc()
    
    def browse_vina_ligand(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Ligand", "", "MOL2/SDF/PDBQT (*.mol2 *.sdf *.pdbqt);;All Files (*)")
        if fn:
            self.vina_ligand.setText(fn)
    
    def run_vina_docking(self):
        """运行 Vina 对接"""
        try:
            from pymol import cmd
            
            # 获取受体对象（从 PyMOL）
            receptor_obj = self.pocket_obj_combo.currentText().strip()
            if not receptor_obj or receptor_obj == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select a receptor object from Pocket tab")
                return
            
            if receptor_obj not in cmd.get_object_list():
                QMessageBox.warning(self, "Warning", f"Object '{receptor_obj}' not found in PyMOL")
                return
            
            # 获取配体文件
            ligand = self.vina_ligand.text().strip()
            if not ligand or not os.path.exists(ligand):
                QMessageBox.warning(self, "Warning", "Please select a valid ligand file")
                return
            
            # 获取参数
            max_pockets = self.vina_max_pockets.value()
            exhaustiveness = self.vina_exhaustiveness.value()
            
            self.log(f"\nStarting Vina docking...")
            self.log(f"   Receptor: {receptor_obj}")
            self.log(f"   Ligand: {os.path.basename(ligand)}")
            self.log(f"   Exhaustiveness: {exhaustiveness}")
            
            # 导入 vina_integration
            try:
                from .vina_integration import pocket_based_docking, manual_box_docking, check_vina_available
            except ImportError:
                from vina_integration import pocket_based_docking, manual_box_docking, check_vina_available
            
            # 检查 Vina 可用性
            if not check_vina_available():
                QMessageBox.warning(
                    self, 
                    "Vina Not Found",
                    "AutoDock Vina is not installed.\n\n"
                    "Please run:\n"
                    "  python glint/env_checker.py --auto-install\n\n"
                    "Or install manually:\n"
                    "  conda install -c conda-forge vina"
                )
                return
            
            # 自定义盒子模式
            if self.vina_use_custom_box.isChecked():
                try:
                    cx = float(self.vina_cx.text().strip())
                    cy = float(self.vina_cy.text().strip())
                    cz = float(self.vina_cz.text().strip())
                    sx = float(self.vina_sx.text().strip())
                    sy = float(self.vina_sy.text().strip())
                    sz = float(self.vina_sz.text().strip())
                except ValueError:
                    QMessageBox.warning(self, "Warning", "Invalid custom box parameters")
                    return
                box = {
                    'center_x': cx, 'center_y': cy, 'center_z': cz,
                    'size_x': sx, 'size_y': sy, 'size_z': sz
                }
                self.log("   Using custom docking box (skip pocket detection)")
                result = manual_box_docking(
                    obj_name=receptor_obj,
                    ligand_file=ligand,
                    box_params=box,
                    exhaustiveness=exhaustiveness
                )
            else:
                self.log(f"   Max pockets: {max_pockets}")
                self.log("   Auto-detecting pockets for docking...")
                result = pocket_based_docking(
                    obj_name=receptor_obj,
                    ligand_file=ligand,
                    max_pockets=max_pockets,
                    exhaustiveness=exhaustiveness
                )
            
            if result.get('success'):
                self.log("Docking complete!")
                self.log(f"   Output directory: {result.get('output_dir')}")
                
                if 'results' in result:
                    self.log(f"\nTop docking results:")
                    sorted_results = sorted(
                        [r for r in result['results'] if r.get('success')],
                        key=lambda r: r.get('affinity', 0)
                    )
                    for i, r in enumerate(sorted_results[:3], 1):
                        self.log(f"   {i}. Pocket {r['pocket_id']}: {r['affinity']:.2f} kcal/mol")
                        self.log(f"      Output: {os.path.basename(r['output_pdbqt'])}")
                
                QMessageBox.information(
                    self,
                    "Success",
                    f"Docking completed successfully!\n\n"
                    f"Results saved to:\n{result.get('output_dir')}"
                )
            else:
                self.log(f"Docking failed: {result.get('error')}")
                QMessageBox.warning(self, "Error", f"Docking failed:\n{result.get('error')}")
            
        except Exception as e:
            self.on_error(str(e))
            import traceback
            traceback.print_exc()
    
    def load_vina_result(self):
        """加载 Vina 对接结果"""
        fn, _ = QFileDialog.getOpenFileName(self, "Select Docking Result", "", "PDBQT (*.pdbqt);;All Files (*)")
        if fn:
            try:
                from pymol import cmd
                obj_name = os.path.splitext(os.path.basename(fn))[0]
                cmd.load(fn, obj_name)
                self.log(f"Loaded docking result: {obj_name}")
                self.refresh_objects()
            except Exception as e:
                self.on_error(str(e))
    
    # ==========================
    # Advanced Pocket Analysis Callbacks
    # ==========================
    def run_pocket_comparison(self):
        """运行口袋对比分析"""
        try:
            obj_a = self.pocket_comp_obj_a.currentText().strip()
            obj_b = self.pocket_comp_obj_b.currentText().strip()
            
            if not obj_a or obj_a == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select Object A")
                return
            if not obj_b or obj_b == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select Object B")
                return
            
            align = self.pocket_comp_align.isChecked()
            
            self.log(f"\nComparing pockets: {obj_a} vs {obj_b}...")
            if align:
                self.log("   Aligning structures first...")
            
            # Import pocket modules
            try:
                from .pocket_detector import compare_pockets
            except ImportError:
                from pocket_detector import compare_pockets
            
            # Compare pockets
            pockets_a, pockets_b, comparison = compare_pockets(
                obj_a, obj_b, 
                align=align,
                grid_spacing=float(self.pocket_grid_spacing.text() or "0.6"),
                min_volume=float(self.pocket_min_volume.text() or "20")
            )
            
            # Store for visualization
            self._comparison_pockets_a = pockets_a
            self._comparison_pockets_b = pockets_b
            self._comparison_result = comparison
            
            # Log results
            self.log(f"\nComparison complete:")
            self.log(f"   {obj_a}: {len(pockets_a)} pockets")
            self.log(f"   {obj_b}: {len(pockets_b)} pockets")
            
            matched = sum(1 for c in comparison if c['match_type'] == 'matched')
            new = sum(1 for c in comparison if c['match_type'] == 'new')
            lost = sum(1 for c in comparison if c['match_type'] == 'lost')
            
            self.log(f"   Matched: {matched}")
            self.log(f"   New in {obj_b}: {new}")
            self.log(f"   Lost from {obj_a}: {lost}")
            
            # Show top changes
            if matched > 0:
                self.log(f"\nTop volume changes:")
                sorted_changes = sorted(
                    [c for c in comparison if c['match_type'] == 'matched'],
                    key=lambda c: abs(c['delta_volume']),
                    reverse=True
                )[:3]
                for c in sorted_changes:
                    delta = c['delta_volume']
                    sign = "+" if delta > 0 else ""
                    self.log(f"   Pocket {c['pocket_a_id']} → {c['pocket_b_id']}: {sign}{delta:.1f} Ų")
            
        except Exception as e:
            self.on_error(str(e))
            import traceback
            traceback.print_exc()
    
    def visualize_pocket_comparison(self):
        """可视化口袋对比结果"""
        try:
            if not hasattr(self, '_comparison_result') or not self._comparison_result:
                QMessageBox.warning(self, "Warning", "Please run pocket comparison first")
                return
            
            self.log(f"\nVisualizing pocket comparison...")
            
            # Import visualizer
            try:
                from .pocket_visualizer import visualize_pocket_comparison
            except ImportError:
                from pocket_visualizer import visualize_pocket_comparison
            
            # Visualize
            visualize_pocket_comparison(
                self._comparison_result,
                self._comparison_pockets_a,
                self._comparison_pockets_b,
                self.pocket_comp_obj_a.currentText(),
                self.pocket_comp_obj_b.currentText()
            )
            
            self.log("✅ Comparison visualized")
            self.log("   Colors: Blue=matched, Green=expanded, Yellow=shrank, Orange=new, Red=lost")
            
        except Exception as e:
            self.on_error(str(e))
            import traceback
            traceback.print_exc()
    
    def run_interface_pockets(self):
        """运行 PPI 界面口袋分析"""
        try:
            obj = self.pocket_interface_obj.currentText().strip()
            chain_a = self.pocket_interface_chain_a.text().strip()
            chain_b = self.pocket_interface_chain_b.text().strip()
            
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select an object")
                return
            if not chain_a or not chain_b:
                QMessageBox.warning(self, "Warning", "Please enter both chain IDs")
                return
            
            self.log(f"\n🔗 Analyzing PPI interface pockets...")
            self.log(f"   Object: {obj}")
            self.log(f"   Interface: Chain {chain_a} + Chain {chain_b}")
            
            # Import pocket integration module
            try:
                from .pocket_glue_integration import analyze_pockets_in_ppi_interface
            except ImportError:
                from pocket_glue_integration import analyze_pockets_in_ppi_interface
            
            # Analyze
            interface_pockets = analyze_pockets_in_ppi_interface(
                obj, chain_a, chain_b,
                visualize=True,
                grid_spacing=float(self.pocket_grid_spacing.text() or "0.6"),
                min_volume=float(self.pocket_min_volume.text() or "20")
            )
            
            # Store for later use
            self._interface_pockets = interface_pockets
            
            # Log results
            self.log(f"\n✅ Found {len(interface_pockets)} interface pockets")
            
            if len(interface_pockets) > 0:
                self.log(f"\n📊 Top interface pockets:")
                sorted_pockets = sorted(
                    interface_pockets,
                    key=lambda p: p['interface_overlap'],
                    reverse=True
                )[:5]
                
                for p in sorted_pockets:
                    self.log(f"   Pocket {p['id']}: "
                            f"Vol={p['volume']:.1f} Ų, "
                            f"Overlap={p['interface_overlap']:.1%}, "
                            f"Drug={p['druggability_score']:.2f}")
            
        except Exception as e:
            self.on_error(str(e))
            import traceback
            traceback.print_exc()
    
    def browse_pocket_corr_csv(self):
        fn, _ = QFileDialog.getOpenFileName(self, "Select Interaction CSV", "", "CSV (*.csv);;All Files (*)")
        if fn:
            self.pocket_corr_csv.setText(fn)
    
    def run_pocket_correlation(self):
        """运行口袋-相互作用关联分析"""
        try:
            csv_file = self.pocket_corr_csv.text().strip()
            
            if not csv_file or not os.path.exists(csv_file):
                QMessageBox.warning(self, "Warning", "Please select a valid interaction CSV file")
                return
            
            # Check if pockets are detected
            if not hasattr(self, '_detected_pockets') or not self._detected_pockets:
                QMessageBox.warning(self, "Warning", 
                    "Please detect pockets first using the 'Pocket Detection' section above")
                return
            
            self.log(f"\n🔍 Correlating pockets with interactions...")
            self.log(f"   CSV: {os.path.basename(csv_file)}")
            self.log(f"   Pockets: {len(self._detected_pockets)}")
            
            # Import correlation module
            try:
                from .pocket_glue_integration import correlate_pockets_with_interactions
            except ImportError:
                from pocket_glue_integration import correlate_pockets_with_interactions
            
            # Correlate
            correlations = correlate_pockets_with_interactions(
                self._detected_pockets,
                csv_file
            )
            
            # Log results
            self.log(f"\n✅ Correlation complete")
            
            if len(correlations) > 0:
                self.log(f"\n📈 Pocket-Interaction correlation:")
                sorted_corr = sorted(
                    correlations,
                    key=lambda c: c['num_interactions'],
                    reverse=True
                )[:5]
                
                for corr in sorted_corr:
                    types_str = ", ".join([f"{k}:{v}" for k, v in list(corr['interaction_types'].items())[:3]])
                    self.log(f"   Pocket {corr['pocket_id']}: "
                            f"{corr['num_interactions']} interactions "
                            f"(density={corr['interaction_density']:.4f}/Ų)")
                    self.log(f"      Types: {types_str}")
            else:
                self.log("No interactions found in pockets")
            
        except Exception as e:
            self.on_error(str(e))
            import traceback
            traceback.print_exc()
    
    def run_gmotif_pocket_analysis(self):
        """运行 G-motif 口袋综合分析"""
        try:
            obj = self.gmotif_pocket_obj.currentText().strip()
            e3_chain = self.gmotif_pocket_e3.text().strip()
            sub_chain = self.gmotif_pocket_sub.text().strip()
            glue_chain = self.gmotif_pocket_glue.text().strip() or None
            
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, "Warning", "Please select an object")
                return
            if not e3_chain or not sub_chain:
                QMessageBox.warning(self, "Warning", "Please enter E3 and Substrate chain IDs")
                return
            
            self.log(f"\nStarting comprehensive G-motif pocket analysis...")
            self.log(f"   Object: {obj}")
            self.log(f"   E3 chain: {e3_chain}")
            self.log(f"   Substrate chain: {sub_chain}")
            if glue_chain:
                self.log(f"   Glue chain: {glue_chain}")
            
            self.log("\nThis may take a few minutes...")
            
            # Import comprehensive analysis module
            try:
                from .pocket_glue_integration import comprehensive_gmotif_pocket_analysis
            except ImportError:
                from pocket_glue_integration import comprehensive_gmotif_pocket_analysis
            
            # Run comprehensive analysis
            results = comprehensive_gmotif_pocket_analysis(
                obj_name=obj,
                e3_chain=e3_chain,
                substrate_chain=sub_chain,
                glue_chain=glue_chain,
                output_dir='gmotif_pocket_analysis'
            )
            
            # Log summary
            self.log(f"\nComprehensive analysis complete!")
            self.log(f"   Output directory: gmotif_pocket_analysis/")
            
            if results.get('gmotif'):
                gmotif = results['gmotif'][0]
                self.log(f"\nG-motif detected:")
                self.log(f"   Position: {gmotif['start_resi']}-{gmotif['end_resi']}")
                self.log(f"   RMSD: {gmotif['rmsd']:.2f} Å")
            
            if results.get('pockets'):
                self.log(f"\nPockets around G-motif: {len(results['pockets'])}")
            
            if results.get('main_binding_pocket'):
                main = results['main_binding_pocket']
                self.log(f"\nMain binding pocket:")
                self.log(f"   Pocket ID: {main['pocket_id']}")
                self.log(f"   Interactions: {main['num_interactions']}")
                self.log(f"   Density: {main['interaction_density']:.4f}/Ų")
                self.log(f"   Druggability: {main['druggability_score']:.3f}")
            
            if results.get('gmotif_interface_pockets'):
                self.log(f"\nG-motif interface pockets: {len(results['gmotif_interface_pockets'])}")
            
            self.log("\nCheck the output directory for detailed reports and CSVs")
            
        except Exception as e:
            self.on_error(str(e))
            import traceback
            traceback.print_exc()
    
    # ==================== 突变分析事件处理 ====================
    
    def run_mutation(self):
        """执行突变"""
        obj_name = self.mut_obj_combo.currentText()
        mutations_str = self.mut_input.text().strip()
        
        if not obj_name or not mutations_str:
            self.log("❌ 请选择对象并输入突变")
            return
        
        try:
            # 解析突变字符串
            mutations = []
            for mut in mutations_str.split(','):
                mut = mut.strip()
                if ':' in mut:
                    parts = mut.split(':')
                    if len(parts) == 3:
                        mutations.append((parts[0], parts[1], parts[2]))
            
            if not mutations:
                self.log("❌ 突变格式错误，示例: A:23:ALA, B:45:GLY")
                return
            
            self.log(f"🧬 执行突变: {obj_name}")
            
            # 调用突变分析模块
            from pymol import cmd
            cmd.perform_mutation(obj_name, mutations, method='pymol')
            
            self.log("✅ 突变完成")
            
        except Exception as e:
            self.log(f"❌ 突变失败: {e}")
            import traceback
            traceback.print_exc()
    
    def run_minimize(self):
        """能量最小化"""
        obj_name = self.mut_obj_combo.currentText()
        
        if not obj_name:
            self.log("❌ 请选择对象")
            return
        
        try:
            self.log(f"⚡ 能量最小化: {obj_name}")
            
            from pymol import cmd
            cmd.minimize_energy(obj_name, cycles=100)
            
            self.log("✅ 最小化完成")
            
        except Exception as e:
            self.log(f"❌ 最小化失败: {e}")
            import traceback
            traceback.print_exc()
    
    def run_mutation_analysis(self):
        """完整突变分析"""
        obj_name = self.mut_obj_combo.currentText()
        mutations_str = self.mut_input.text().strip()
        method = self.mut_method_combo.currentText().lower()
        
        if not obj_name or not mutations_str:
            self.log("❌ 请选择对象并输入突变")
            return
        
        try:
            # 解析突变
            mutations = []
            for mut in mutations_str.split(','):
                mut = mut.strip()
                if ':' in mut:
                    parts = mut.split(':')
                    if len(parts) == 3:
                        mutations.append((parts[0], parts[1], parts[2]))
            
            if not mutations:
                self.log("❌ 突变格式错误")
                return
            
            self.log(f"🧬 开始完整分析: {obj_name}")
            self.log(f"突变数量: {len(mutations)}")
            self.log(f"方法: {method}")
            
            # 调用完整分析
            from pymol import cmd
            cmd.analyze_mutation_effects(obj_name, mutations, method=method)
            
            self.log("✅ 分析完成")
            
        except Exception as e:
            self.log(f"❌ 分析失败: {e}")
            import traceback
            traceback.print_exc()
    
    # ==================== 其他方法 ====================
    
    def run_compare_scoring(self):
        """对比评分方法 - DEPRECATED"""
        QMessageBox.information(
            self,
            "Feature Removed",
            "Scoring comparison has been removed.\n\n"
            "Please use:\n"
            "• Docking tab for AutoDock Vina scoring\n"
            "• Interaction Analysis tab for interaction details"
        )
        self.log("Scoring comparison feature has been removed")
    
    def run_ternary_scoring(self):
        """运行三元复合物评分 - DEPRECATED"""
        QMessageBox.information(
            self,
            "Feature Removed",
            "Ternary complex scoring has been removed.\n\n"
            "Please use:\n"
            "• Molecular Glue tab → Ternary Complex Analysis for interaction analysis"
        )
        self.log("Ternary scoring feature has been removed")
    
    def browse_heatmap_folder(self):
        """浏览选择热图 CSV 文件夹"""
        from PyQt5.QtWidgets import QFileDialog
        folder = QFileDialog.getExistingDirectory(self, "Select Folder Containing CSV Files", "")
        if folder:
            self.heatmap_folder.setText(folder)
    
    def run_generate_heatmap(self):
        """生成批量热图"""
        try:
            folder = self.heatmap_folder.text().strip()
            pattern = self.heatmap_pattern.text().strip() or "*_scores.csv"
            
            if not folder:
                QMessageBox.warning(self, "Warning", "Please select a folder")
                return
            
            if not os.path.exists(folder):
                QMessageBox.warning(self, "Warning", f"Folder does not exist: {folder}")
                return
            
            self.log(f"\nGenerating heatmap from: {folder}")
            self.log(f"   Pattern: {pattern}")
            self.score_result_text.clear()
            self.score_result_text.setPlainText("Generating heatmap...\nThis may take a few seconds...")
            
            # Import heatmap module
            try:
                from .binding_heatmap import generate_binding_heatmap
            except ImportError:
                from binding_heatmap import generate_binding_heatmap
            
            # Generate heatmap
            result = generate_binding_heatmap(folder, pattern=pattern)
            
            if result['success']:
                report = f"""
{'='*60}
Binding Energy Heatmap Generated
{'='*60}
Receptors: {result['n_receptors']}
Ligands:   {result['n_ligands']}
Output:    {result['output_path']}
{'='*60}
Heatmap saved successfully!
{'='*60}
                """
                self.score_result_text.setPlainText(report)
                self.log(f"Heatmap saved: {result['output_path']}")
                
                # 询问是否打开图片
                reply = QMessageBox.question(
                    self, 
                    "Success", 
                    f"Heatmap generated successfully!\n\nOpen the image?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
                )
                
                if reply == QMessageBox.StandardButton.Yes:
                    import subprocess
                    import sys
                    if sys.platform == 'darwin':  # macOS
                        subprocess.run(['open', result['output_path']])
                    elif sys.platform == 'win32':  # Windows
                        os.startfile(result['output_path'])
                    else:  # Linux
                        subprocess.run(['xdg-open', result['output_path']])
            else:
                error_msg = result.get('error', 'Unknown error')
                self.log(f"Heatmap generation failed: {error_msg}")
                self.score_result_text.setPlainText(f"Failed to generate heatmap:\n{error_msg}")
                QMessageBox.warning(self, "Error", f"Failed to generate heatmap:\n{error_msg}")
            
        except Exception as e:
            self.on_error(str(e))
            import traceback
            traceback.print_exc()
    
    def update_enablement(self):
        gm_ok = getattr(self, "obj_combo_gm", None) and self.obj_combo_gm.currentText().strip() not in ("", t("no_object"))
        if getattr(self, "gm_btn", None): self.gm_btn.setEnabled(bool(gm_ok))
        if getattr(self, "gm_btn_render", None): self.gm_btn_render.setEnabled(bool(gm_ok))

        ana_ok = getattr(self, "obj_combo_analysis", None) and self.obj_combo_analysis.currentText().strip() not in ("", t("no_object"))
        if getattr(self, "analyze_btn", None): self.analyze_btn.setEnabled(bool(ana_ok))

        obj_ok = getattr(self, "obj_combo_csv", None) and self.obj_combo_csv.currentText().strip() not in ("", t("no_object"))
        if getattr(self, "highlight_btn", None): self.highlight_btn.setEnabled(bool(obj_ok and self.csv_path.text().strip()))
        if getattr(self, "clear_btn", None): self.clear_btn.setEnabled(True)

        apbs_ok = getattr(self, "obj_combo_apbs", None) and self.obj_combo_apbs.currentText().strip() not in ("", t("no_object"))
        for b in (getattr(self, "btn_apbs_quick", None), getattr(self, "btn_apbs_true", None),
                  getattr(self, "btn_export_png", None), getattr(self, "btn_export_dx", None),
                  getattr(self, "btn_viewport", None)):
            if b: b.setEnabled(bool(apbs_ok))

    def log(self, text: str):
        self.log_edit.append(text)
        self.log_edit.moveCursor(self.log_edit.textCursor().End)
    
    def toggle_theme(self):
        """Toggle theme (dark/light)"""
        self._dark_mode = not self._dark_mode
        # 先清空样式再重新应用，避免残留规则
        try:
            self.setStyleSheet("")
        except Exception:
            pass
        self.setup_style()
        # Update theme icon
        # self.theme_toggle_btn.setText("☾" if self._dark_mode else "☀")
        # Update Modules button color (easter egg effect)
        self.update_modules_button_style()
        # Update separator color based on theme
        if hasattr(self, 'nav_separator'):
            sep_color = "#30363d" if self._dark_mode else "#e2e8f0"
            self.nav_separator.setStyleSheet(f"background-color: {sep_color}; margin: 8px 0;")
        # Update tab widget backgrounds (macOS compatibility)
        if hasattr(self, 'pocket_advanced_tabs'):
            bg_color = "#161b22" if self._dark_mode else "white"
            for i in range(self.pocket_advanced_tabs.count()):
                tab_widget = self.pocket_advanced_tabs.widget(i)
                if tab_widget:
                    tab_widget.setStyleSheet(f"#tab_content {{ background-color: {bg_color}; }}")
        # Update scroll area content backgrounds
        scroll_bg = "#0d1117" if self._dark_mode else "#f8fafc"
        for attr in ('_interaction_scroll_content', '_molecular_glue_scroll_content', '_docking_scroll_content'):
            if hasattr(self, attr):
                widget = getattr(self, attr)
                if widget:
                    widget.setStyleSheet(f"#scroll_content {{ background-color: {scroll_bg}; }}")
        # 主题变化后重新应用自动缩放，以确保字体与行距匹配
        self.apply_auto_scaling()
        try:
            self.update()
        except Exception:
            pass
        theme_name_en = "Dark Theme" if self._dark_mode else "Light Theme"
        self.log(f"Switched to {theme_name_en}")
    
    
    def update_modules_button_style(self):
        """Make Modules look like a centered, borderless label"""
        style = """
            QPushButton#modules_btn {
                background: transparent;
                color: #64748b;
                border: none;
                padding: 0px;
                font-weight: 600;
                font-size: 14px;
            }
            QPushButton#modules_btn:hover {
                color: #3b82f6;
            }
        """
        self.modules_btn.setStyleSheet(style)

    # --- 自动缩放（随屏幕/系统字体）---
    def _compute_ui_scale(self) -> float:
        try:
            screen = QApplication.primaryScreen()
            dpi = screen.logicalDotsPerInch()
            # 基准 DPI 96
            return max(1.0, dpi / 96.0)
        except Exception:
            return 1.0

    def apply_auto_scaling(self):
        """应用自动缩放"""
        # 暂时禁用，因为使用了布局调优
        pass

# ============================================================================
# Disease Analysis Components (Integrated from disease_analysis_gui.py)
# ============================================================================

class DiseaseQueryWorker(QThread):
    """
    疾病查询工作线程（避免 GUI 卡顿）
    """
    progress = pyqtSignal(str)
    finished = pyqtSignal(object)  # DataFrame
    error = pyqtSignal(str)
    
    def __init__(self, disease_name: str, top_n: int, include_e3: bool, e3_symbol: str):
        super().__init__()
        self.disease_name = disease_name
        self.top_n = top_n
        self.include_e3 = include_e3
        self.e3_symbol = e3_symbol
    
    def run(self):
        try:
            # 步骤 1: 搜索疾病
            self.progress.emit(f"Searching for disease: {self.disease_name}...")
            if search_disease is None:
                raise ImportError("Open Targets API modules not available")

            candidates = search_disease(self.disease_name, max_results=5)
            
            if not candidates:
                self.error.emit(f"No disease found for '{self.disease_name}'")
                return
            
            # 使用第一个候选
            disease_id = candidates[0]['id']
            disease_name = candidates[0]['name']
            
            # 步骤 2: 获取靶点
            self.progress.emit(f"Fetching targets for {disease_name}...")
            df = get_disease_targets(disease_id, top_n=self.top_n)
            
            if df is None or df.empty:
                self.error.emit(f"No targets found for {disease_name}")
                return
            
            # 步骤 3: 添加 E3 评分
            if self.include_e3:
                self.progress.emit(f"Calculating E3 compatibility scores ({self.e3_symbol})...")
                df = enrich_targets_with_e3_scores(df, e3_symbol=self.e3_symbol)
            
            # 添加疾病信息
            df['disease_name'] = disease_name
            df['disease_id'] = disease_id
            
            self.progress.emit(f"Query complete: {len(df)} targets found")
            self.finished.emit(df)
            
        except Exception as e:
            self.error.emit(f"Query failed: {str(e)}")


class DiseaseAnalysisTab(QWidget):
    """
    疾病分析主标签页
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_results = None  # 当前查询结果 DataFrame
        self.query_worker = None
        
        self.init_ui()
    
    def init_ui(self):
        """初始化界面"""
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(12)
        main_layout.setContentsMargins(12, 12, 12, 12)
        
        # ========== 疾病搜索区域 ==========
        search_group = QGroupBox("Disease Search")
        search_layout = QVBoxLayout(search_group)
        search_layout.setSpacing(10)
        
        # 疾病名称输入
        disease_row = QHBoxLayout()
        disease_row.addWidget(QLabel("Disease Name:"))
        self.disease_input = QLineEdit()
        self.disease_input.setPlaceholderText("e.g., multiple myeloma, breast cancer")
        self.disease_input.setMinimumHeight(32)
        disease_row.addWidget(self.disease_input, 1)
        search_layout.addLayout(disease_row)
        
        # 参数设置
        params_layout = QHBoxLayout()
        
        # Top N
        params_layout.addWidget(QLabel("Top N Targets:"))
        self.top_n_spin = QSpinBox()
        self.top_n_spin.setRange(5, 100)
        self.top_n_spin.setValue(30)
        self.top_n_spin.setMinimumHeight(32)
        params_layout.addWidget(self.top_n_spin)
        
        params_layout.addSpacing(20)
        
        # E3 连接酶
        params_layout.addWidget(QLabel("E3 Ligase:"))
        self.e3_combo = QComboBox()
        self.e3_combo.addItems(E3_LIGASES)
        self.e3_combo.setMinimumHeight(32)
        params_layout.addWidget(self.e3_combo)
        
        params_layout.addSpacing(20)
        
        # E3 评分选项
        self.e3_checkbox = QCheckBox("Include E3 Score")
        self.e3_checkbox.setChecked(True)
        params_layout.addWidget(self.e3_checkbox)
        
        params_layout.addStretch()
        search_layout.addLayout(params_layout)
        
        # 查询按钮
        btn_row = QHBoxLayout()
        self.query_btn = QPushButton("🔍 Query Targets")
        self.query_btn.setMinimumHeight(36)
        self.query_btn.setObjectName("primary_btn")
        self.query_btn.clicked.connect(self.on_query_clicked)
        btn_row.addWidget(self.query_btn)
        btn_row.addStretch()
        search_layout.addLayout(btn_row)
        
        # 进度条
        self.progress_bar = QProgressBar()
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setVisible(False)
        search_layout.addWidget(self.progress_bar)
        
        main_layout.addWidget(search_group)
        
        # ========== 结果展示区域 ==========
        results_group = QGroupBox("Results")
        results_layout = QVBoxLayout(results_group)
        results_layout.setSpacing(10)
        
        # 结果表格
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels([
            "Symbol", "Name", "Disease Score", "E3 Score", "Composite Score", "Actions"
        ])
        
        # 设置表格属性
        from PyQt5.QtWidgets import QAbstractItemView, QHeaderView
        self.results_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.results_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.results_table.setSortingEnabled(True)
        
        # 设置列宽
        header = self.results_table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)  # Symbol
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)  # Name
        header.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)  # Disease Score
        header.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)  # E3 Score
        header.setSectionResizeMode(4, QHeaderView.ResizeMode.ResizeToContents)  # Composite Score
        header.setSectionResizeMode(5, QHeaderView.ResizeMode.Fixed)  # Actions
        self.results_table.setColumnWidth(5, 120)
        
        results_layout.addWidget(self.results_table)
        
        # 导出按钮
        export_row = QHBoxLayout()
        self.export_csv_btn = QPushButton("📄 Export CSV")
        self.export_csv_btn.setMinimumHeight(32)
        self.export_csv_btn.clicked.connect(self.on_export_csv)
        self.export_csv_btn.setEnabled(False)
        export_row.addWidget(self.export_csv_btn)
        export_row.addStretch()
        results_layout.addLayout(export_row)
        
        main_layout.addWidget(results_group, 1)
        
        # ========== 状态栏 ==========
        self.status_label = QLabel("Ready")
        self.status_label.setStyleSheet("color: #666; font-size: 11px;")
        main_layout.addWidget(self.status_label)
    
    def on_query_clicked(self):
        """查询按钮点击事件"""
        disease_name = self.disease_input.text().strip()
        
        if not disease_name:
            QMessageBox.warning(self, "Input Required", "Please enter a disease name")
            return
        
        # 检查模块是否可用
        if search_disease is None or get_disease_targets is None:
            QMessageBox.critical(
                self,
                "Module Not Available",
                "Open Targets API modules are not available.\n"
                "Please ensure the required modules are installed."
            )
            return
        
        # 禁用查询按钮
        self.query_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)  # 不确定进度
        
        # 创建并启动 Worker 线程
        top_n = self.top_n_spin.value()
        include_e3 = self.e3_checkbox.isChecked()
        e3_symbol = self.e3_combo.currentText()
        
        self.query_worker = DiseaseQueryWorker(disease_name, top_n, include_e3, e3_symbol)
        self.query_worker.progress.connect(self.on_query_progress)
        self.query_worker.finished.connect(self.on_query_finished)
        self.query_worker.error.connect(self.on_query_error)
        self.query_worker.start()
    
    def on_query_progress(self, message: str):
        """查询进度更新"""
        self.status_label.setText(message)
    
    def on_query_finished(self, df):
        """查询完成"""
        self.current_results = df
        self.populate_results_table(df)
        
        # 恢复 UI
        self.query_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.export_csv_btn.setEnabled(True)
        
        disease_name = df['disease_name'].iloc[0] if not df.empty else "Unknown"
        self.status_label.setText(f"Query complete: {len(df)} targets found for {disease_name}")
    
    def on_query_error(self, error_msg: str):
        """查询错误"""
        self.query_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        
        QMessageBox.warning(
            self,
            "Query Error",
            f"Failed to query disease targets:\n\n{error_msg}\n\n"
            "Please check:\n"
            "• Network connection\n"
            "• Disease name spelling\n"
            "• API availability"
        )
        
        self.status_label.setText(f"Error: {error_msg}")
    
    def populate_results_table(self, df):
        """填充结果表格"""
        if df is None or df.empty:
            self.results_table.setRowCount(0)
            return
        
        self.results_table.setRowCount(len(df))
        self.results_table.setSortingEnabled(False)  # 填充时禁用排序
        
        # Check if pandas is available
        try:
            import pandas as pd
        except ImportError:
            return

        for i, row in df.iterrows():
            # Symbol
            self.results_table.setItem(i, 0, QTableWidgetItem(row['symbol']))
            
            # Name
            name = row.get('name', '')
            if len(name) > 50:
                name = name[:47] + "..."
            self.results_table.setItem(i, 1, QTableWidgetItem(name))
            
            # Disease Score
            disease_score = f"{row['score']:.4f}"
            self.results_table.setItem(i, 2, QTableWidgetItem(disease_score))
            
            # E3 Score
            if 'e3_score' in row and pd.notna(row['e3_score']):
                e3_score = f"{row['e3_score']:.4f}"
                self.results_table.setItem(i, 3, QTableWidgetItem(e3_score))
            else:
                self.results_table.setItem(i, 3, QTableWidgetItem("N/A"))
            
            # Composite Score
            if 'composite_score' in row and pd.notna(row['composite_score']):
                composite_score = f"{row['composite_score']:.4f}"
                self.results_table.setItem(i, 4, QTableWidgetItem(composite_score))
            else:
                self.results_table.setItem(i, 4, QTableWidgetItem("N/A"))
            
            # Actions - Load Structure 按钮
            load_btn = QPushButton("Load")
            load_btn.setMinimumHeight(28)
            load_btn.clicked.connect(lambda checked, symbol=row['symbol']: self.on_load_structure(symbol))
            self.results_table.setCellWidget(i, 5, load_btn)
        
        self.results_table.setSortingEnabled(True)  # 重新启用排序
    
    def on_load_structure(self, symbol: str):
        """加载蛋白结构到 PyMOL"""
        try:
            from pymol import cmd
            
            reply = QMessageBox.question(
                self,
                "Load Structure",
                f"Load structure for {symbol}?\n\n"
                "This will attempt to fetch the structure from PDB.\n"
                "You can also manually load from AlphaFold or other sources.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No
            )
            
            if reply == QMessageBox.StandardButton.Yes:
                # 尝试 fetch
                try:
                    # 先显示提示
                    QMessageBox.information(
                        self,
                        "Load Structure",
                        f"Instructions:\n\n"
                        f"1. Fetching {symbol} from PDB...\n"
                        f"2. If PDB not found, try loading AlphaFold model.\n"
                        f"3. Run: ot_glue_insight protein_obj=\"{symbol}\", gene_symbol=\"{symbol}\""
                    )
                    # 简单的 fetch (假设 symbol 也是 PDB ID，或者需要查询)
                    # 实际上 gene symbol != PDB ID. 
                    # 这里我们只是提供指导，或者如果 PyMOL 连接了网络，可以尝试 fetch alphafold
                    cmd.fetch(symbol, type="alphafold") # PyMOL 2.5+ supports fetch <uniprot> or alphafold
                except Exception as e:
                    QMessageBox.warning(self, "Fetch Failed", f"Could not fetch structure: {e}")

        except ImportError:
            QMessageBox.warning(
                self,
                "PyMOL Not Available",
                "PyMOL is not available in this environment"
            )
    
    def on_export_csv(self):
        """导出 CSV"""
        if self.current_results is None or self.current_results.empty:
            QMessageBox.warning(self, "No Data", "No results to export")
            return
        
        # 文件对话框
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export CSV",
            os.path.join(DEFAULT_OUTPUT_DIR, "disease_targets.csv"),
            "CSV Files (*.csv)"
        )
        
        if file_path:
            try:
                # 确保目录存在
                os.makedirs(os.path.dirname(file_path), exist_ok=True)
                
                # 导出
                self.current_results.to_csv(file_path, index=False)
                
                QMessageBox.information(
                    self,
                    "Export Successful",
                    f"Results exported to:\n{file_path}"
                )
                
                self.status_label.setText(f"Exported to: {file_path}")
                
            except Exception as e:
                QMessageBox.critical(
                    self,
                    "Export Failed",
                    f"Failed to export CSV:\n{str(e)}"
                )


# ============================================================================
# GLINTDialog continued - Style and Scaling Methods
# (These belong to GLINTDialog class, not DiseaseAnalysisTab)
# ============================================================================

# Monkey-patch these methods onto GLINTDialog
def _glint_setup_style(self):
    # 现代化样式 - 根据深色/浅色主题应用不同的样式
    # Mac兼容性优先
    if self._dark_mode:
        # 深色主题 - 优化版（更现代，更舒适，大字体）
        self.setStyleSheet("""
QDialog {
    background-color: #0d1117;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Microsoft YaHei", sans-serif;
    font-size: 14px;
    color: #e6edf3;
}

/* 导航区域 - 透明背景 */
QWidget#nav_widget {
    background-color: transparent;
}

/* 导航列表 - 更现代的深色 */
QListWidget {
    background-color: #161b22;
    border: 1px solid #30363d;
    border-radius: 10px;
    padding: 10px;
    outline: none;
}
QListWidget::item {
    padding: 14px 16px;
    margin: 4px 0;
    border-radius: 8px;
    color: #c9d1d9;
    font-size: 15px;
    font-weight: 500;
}
QListWidget::item:hover {
    background-color: #21262d;
    color: #e6edf3;
}
QListWidget::item:selected {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #3b82f6, stop:1 #2563eb);
    color: white;
    font-weight: 600;
    border: none;
}

/* 分组框 - 更清晰的边框和背景 */
QGroupBox {
    font-weight: 600;
    font-size: 16px;
    border: 1px solid #30363d;
    border-radius: 12px;
    margin-top: 14px;
    padding: 24px 16px 16px 16px;
    background-color: #161b22;
    color: #e6edf3;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 16px;
    padding: 0 8px;
    background-color: #161b22;
    color: #58a6ff;
}

/* 标签 - 更好的对比度 */
QLabel {
    color: #c9d1d9;
    font-size: 14px;
    padding: 4px 0;
    min-height: 24px;
}

/* 输入框和下拉框 - 更清晰的视觉效果 */
QLineEdit, QComboBox {
    border: 1px solid #30363d;
    border-radius: 8px;
    padding: 10px 14px;
    background-color: #0d1117;
    color: #e6edf3;
    font-size: 14px;
    selection-background-color: #3b82f6;
    min-height: 20px;
}
QLineEdit:focus, QComboBox:focus {
    border: 2px solid #58a6ff;
    background-color: #161b22;
}
QLineEdit:hover, QComboBox:hover {
    border-color: #58a6ff;
    background-color: #161b22;
}
QComboBox::drop-down {
    border: none;
    padding-right: 10px;
    width: 20px;
}
QComboBox::down-arrow {
    image: none;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 6px solid #8b949e;
    margin-right: 8px;
}

/* 数值输入框 */
QSpinBox {
    border: 1px solid #30363d;
    border-radius: 8px;
    padding: 10px 14px;
    background-color: #0d1117;
    color: #e6edf3;
    font-size: 14px;
    selection-background-color: #3b82f6;
    min-height: 20px;
}
QSpinBox:focus {
    border: 2px solid #58a6ff;
    background-color: #161b22;
}
QSpinBox:hover {
    border-color: #58a6ff;
    background-color: #161b22;
}
QSpinBox::up-button, QSpinBox::down-button {
    background-color: #21262d;
    border: none;
    width: 20px;
    border-radius: 4px;
    margin: 1px;
}
QSpinBox::up-button:hover, QSpinBox::down-button:hover {
    background-color: #30363d;
}
QSpinBox::up-arrow {
    image: none;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-bottom: 5px solid #8b949e;
    width: 0; height: 0;
}
QSpinBox::down-arrow {
    image: none;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-top: 5px solid #8b949e;
    width: 0; height: 0;
}

/* 普通按钮 - 更现代的样式 */
QPushButton {
    border: 1px solid #30363d;
    border-radius: 8px;
    padding: 10px 20px;
    background-color: #21262d;
    color: #c9d1d9;
    font-weight: 600;
    font-size: 14px;
}
QPushButton:hover {
    background-color: #30363d;
    border-color: #58a6ff;
    color: #e6edf3;
}
QPushButton:pressed {
    background-color: #161b22;
    transform: translateY(1px);
}
QPushButton:disabled {
    background-color: #161b22;
    color: #484f58;
    border-color: #21262d;
}

/* 主操作按钮 - 更鲜艳的蓝色渐变 */
QPushButton#primary_btn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #3b82f6, stop:1 #2563eb);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 700;
    font-size: 14px;
    padding: 12px 24px;
}
QPushButton#primary_btn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #60a5fa, stop:1 #3b82f6);
    box-shadow: 0 4px 12px rgba(59, 130, 246, 0.3);
}
QPushButton#primary_btn:pressed {
    background: #1d4ed8;
}
QPushButton#primary_btn:disabled {
    background-color: #30363d;
    color: #6e7681;
    background: #30363d;
}

/* 次要按钮 */
QPushButton#secondary_btn {
    background-color: #30363d;
    color: #c9d1d9;
    border: 1px solid #484f58;
    border-radius: 8px;
    font-weight: 600;
    font-size: 14px;
    padding: 12px 24px;
}
QPushButton#secondary_btn:hover {
    background-color: #484f58;
    border-color: #6e7681;
    color: #e6edf3;
}
QPushButton#secondary_btn:pressed {
    background-color: #21262d;
}

/* 图标按钮 */
QPushButton#icon_btn {
    background-color: #21262d;
    color: #8b949e;
    border: 1px solid #30363d;
    border-radius: 8px;
    font-size: 18px;
    padding: 8px;
}
QPushButton#icon_btn:hover {
    background-color: #30363d;
    border-color: #58a6ff;
    color: #c9d1d9;
}
QPushButton#icon_btn:pressed {
    background-color: #161b22;
}

/* Highlight 按钮 */
QPushButton#highlight_btn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #3b82f6, stop:1 #2563eb);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    font-size: 14px;
    padding: 10px 20px;
}
QPushButton#highlight_btn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #60a5fa, stop:1 #3b82f6);
}
QPushButton#highlight_btn:pressed {
    background: #1d4ed8;
}

/* 刷新按钮 */
QPushButton#refresh_btn {
    background-color: #21262d;
    color: #8b949e;
    border: 1px solid #30363d;
    border-radius: 8px;
    padding: 10px 14px;
    font-size: 12px;
    min-width: 70px;
}
QPushButton#refresh_btn:hover {
    background-color: #30363d;
    color: #c9d1d9;
    border-color: #58a6ff;
}

/* 主题切换按钮 */
QPushButton#theme_toggle_btn {
    background-color: #21262d;
    color: #f39c12;
    border: 1px solid #30363d;
    border-radius: 8px;
    padding: 6px;
    font-size: 20px;
    font-weight: normal;
}
QPushButton#theme_toggle_btn:hover {
    background-color: #30363d;
    border-color: #58a6ff;
    color: #f1c40f;
}
QPushButton#theme_toggle_btn:pressed {
    background-color: #161b22;
}

/* 底部导航按钮 */
QPushButton#bottom_nav_btn {
    background-color: transparent;
    color: #8b949e;
    border: none;
    text-align: left;
    padding: 12px 16px;
    font-size: 14px;
    font-weight: 500;
    border-radius: 8px;
}
QPushButton#bottom_nav_btn:hover {
    background-color: #21262d;
    color: #58a6ff;
}

/* 浏览按钮 */
QPushButton#browse_btn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #3b82f6, stop:1 #2563eb);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    min-width: 80px;
    padding: 10px 16px;
    font-size: 13px;
}
QPushButton#browse_btn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #60a5fa, stop:1 #3b82f6);
}

/* 保存按钮 - 绿色 */
QPushButton#save_btn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #10b981, stop:1 #059669);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    min-width: 80px;
    padding: 10px 16px;
    font-size: 13px;
}
QPushButton#save_btn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #34d399, stop:1 #10b981);
}

/* 清除按钮 */
QPushButton#clear_btn {
    background-color: #484f58;
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    padding: 10px 20px;
}
QPushButton#clear_btn:hover {
    background-color: #6e7681;
}

/* 关闭按钮 - 红色 */
QPushButton#close_btn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #ef4444, stop:1 #dc2626);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    min-width: 100px;
    padding: 10px 20px;
}
QPushButton#close_btn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #f87171, stop:1 #ef4444);
}

/* 文本编辑器 - 代码风格 */
QTextEdit {
    border: 1px solid #30363d;
    border-radius: 8px;
    padding: 12px;
    background-color: #0d1117;
    font-family: 'SF Mono', Menlo, Monaco, 'Courier New', monospace;
    font-size: 13px;
    color: #c9d1d9;
    selection-background-color: #3b82f6;
}

/* 表格 - 更现代的设计 */
QTableWidget {
    border: 1px solid #30363d;
    border-radius: 10px;
    background-color: #0d1117;
    gridline-color: #21262d;
    selection-background-color: #3b82f6;
    selection-color: white;
    font-size: 14px;
}
QHeaderView::section {
    background-color: #161b22;
    color: #e6edf3;
    padding: 10px 12px;
    border: none;
    border-bottom: 2px solid #30363d;
    font-weight: 600;
    font-size: 13px;
}
QTableWidget::item {
    padding: 12px 10px;
    min-height: 40px;
    border-bottom: 1px solid #21262d;
    color: #c9d1d9;
}
QTableWidget::item:hover {
    background-color: #161b22;
}
QTableWidget::item:selected {
    background-color: #3b82f6;
    color: white;
}
QTableWidget QTableCornerButton::section {
    background-color: #161b22;
    border: none;
    border-bottom: 2px solid #30363d;
}

/* 复选框 */
QCheckBox {
    color: #c9d1d9;
    spacing: 10px;
    font-size: 14px;
    padding: 6px 0;
}
QCheckBox::indicator {
    width: 20px;
    height: 20px;
    border: 1.5px solid #30363d;
    border-radius: 5px;
    background-color: #0d1117;
}
QCheckBox::indicator:hover {
    border-color: #58a6ff;
    background-color: #161b22;
}
QCheckBox::indicator:checked {
    background-color: #3b82f6;
    border-color: #3b82f6;
    image: none;
}
QCheckBox::indicator:checked:after {
    content: "✓";
    color: white;
    padding-left: 4px;
    font-weight: bold;
}

/* 进度条 */
QProgressBar {
    border: 1px solid #30363d;
    border-radius: 8px;
    text-align: center;
    background-color: #161b22;
    height: 28px;
    color: #c9d1d9;
    font-size: 12px;
    font-weight: 500;
}
QProgressBar::chunk {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #3b82f6, stop:1 #2563eb);
    border-radius: 6px;
}

/* 标签页 - QTabWidget */
QTabWidget::pane {
    border: 1px solid #30363d;
    border-radius: 8px;
    background-color: #161b22;
    top: -1px;
}
QTabWidget > QWidget {
    background-color: #161b22;
}
QTabBar::tab {
    background-color: #0d1117;
    color: #8b949e;
    border: 1px solid #30363d;
    border-bottom: none;
    border-top-left-radius: 8px;
    border-top-right-radius: 8px;
    padding: 10px 20px;
    margin-right: 2px;
    font-size: 13px;
    font-weight: 500;
}
QTabBar::tab:selected {
    background-color: #161b22;
    color: #e6edf3;
    font-weight: 600;
    border-bottom: 3px solid #3b82f6;
}
QTabBar::tab:hover {
    background-color: #21262d;
    color: #c9d1d9;
}

/* 滚动条 - 更精致 */
QScrollBar:vertical {
    border: none;
    background: #0d1117;
    width: 14px;
    border-radius: 7px;
    margin: 2px;
}
QScrollBar::handle:vertical {
    background: #30363d;
    border-radius: 7px;
    min-height: 40px;
}
QScrollBar::handle:vertical:hover {
    background: #484f58;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0;
}
QScrollBar:horizontal {
    border: none;
    background: #0d1117;
    height: 14px;
    border-radius: 7px;
    margin: 2px;
}
QScrollBar::handle:horizontal {
    background: #30363d;
    border-radius: 7px;
    min-width: 40px;
}
QScrollBar::handle:horizontal:hover {
    background: #484f58;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
    width: 0;
}
""")
    else:
        # 浅色主题 - 优化版（更现代，与深色主题一致，大字体）
        self.setStyleSheet("""
QDialog {
    background-color: #f8fafc;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Microsoft YaHei", sans-serif;
    font-size: 14px;
    color: #1e293b;
}

/* 导航区域 - 透明背景 */
QWidget#nav_widget {
    background-color: transparent;
}

/* 导航列表 - 更清爭的白色卡片 */
QListWidget {
    background-color: white;
    border: 1px solid #e2e8f0;
    border-radius: 10px;
    padding: 10px;
    outline: none;
}
QListWidget::item {
    padding: 14px 16px;
    margin: 4px 0;
    border-radius: 8px;
    color: #475569;
    font-size: 15px;
    font-weight: 500;
}
QListWidget::item:hover {
    background-color: #f1f5f9;
    color: #1e293b;
}
QListWidget::item:selected {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #3b82f6, stop:1 #2563eb);
    color: white;
    font-weight: 600;
    border: none;
}

/* 分组框 - 更清晰的卡片风格 */
QGroupBox {
    font-weight: 600;
    font-size: 16px;
    border: 1px solid #e2e8f0;
    border-radius: 12px;
    margin-top: 14px;
    padding: 24px 16px 16px 16px;
    background-color: white;
    color: #0f172a;
    box-shadow: 0 1px 3px rgba(0,0,0,0.05);
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 16px;
    padding: 0 8px;
    background-color: white;
    color: #3b82f6;
}

/* 标签 - 更好的对比度 */
QLabel {
    color: #475569;
    font-size: 14px;
    padding: 4px 0;
    min-height: 24px;
}

/* 输入框和下拉框 - 纯白背景 */
QLineEdit, QComboBox {
    border: 1px solid #cbd5e1;
    border-radius: 8px;
    padding: 10px 14px;
    background-color: white;
    color: #1e293b;
    font-size: 14px;
    selection-background-color: #3b82f6;
    min-height: 20px;
}
QLineEdit:focus, QComboBox:focus {
    border: 2px solid #3b82f6;
    background-color: white;
}
QLineEdit:hover, QComboBox:hover {
    border-color: #3b82f6;
    background-color: white;
}
QComboBox::drop-down {
    border: none;
    padding-right: 10px;
    width: 20px;
}
QComboBox::down-arrow {
    image: none;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 6px solid #64748b;
    margin-right: 8px;
}

/* 数值输入框 */
QSpinBox {
    border: 1px solid #cbd5e1;
    border-radius: 8px;
    padding: 10px 14px;
    background-color: white;
    color: #1e293b;
    font-size: 14px;
    selection-background-color: #3b82f6;
    min-height: 20px;
}
QSpinBox:focus {
    border: 2px solid #3b82f6;
    background-color: white;
}
QSpinBox:hover {
    border-color: #3b82f6;
    background-color: white;
}
QSpinBox::up-button, QSpinBox::down-button {
    background-color: #f1f5f9;
    border: none;
    width: 20px;
    border-radius: 4px;
    margin: 1px;
}
QSpinBox::up-button:hover, QSpinBox::down-button:hover {
    background-color: #e2e8f0;
}
QSpinBox::up-arrow {
    image: none;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-bottom: 5px solid #64748b;
    width: 0; height: 0;
}
QSpinBox::down-arrow {
    image: none;
    border-left: 4px solid transparent;
    border-right: 4px solid transparent;
    border-top: 5px solid #64748b;
    width: 0; height: 0;
}

/* 普通按钮 - 更现代的样式 */
QPushButton {
    border: 1px solid #cbd5e1;
    border-radius: 8px;
    padding: 10px 20px;
    background-color: white;
    color: #475569;
    font-weight: 600;
    font-size: 14px;
}
QPushButton:hover {
    background-color: #f1f5f9;
    border-color: #3b82f6;
    color: #1e293b;
}
QPushButton:pressed {
    background-color: #e2e8f0;
    transform: translateY(1px);
}
QPushButton:disabled {
    background-color: #f8fafc;
    color: #cbd5e1;
    border-color: #e2e8f0;
}

/* 主操作按钮 - 蓝色渐变 */
QPushButton#primary_btn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #3b82f6, stop:1 #2563eb);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 700;
    font-size: 14px;
    padding: 12px 24px;
}
QPushButton#primary_btn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #60a5fa, stop:1 #3b82f6);
    box-shadow: 0 4px 12px rgba(59, 130, 246, 0.3);
}
QPushButton#primary_btn:pressed {
    background: #1d4ed8;
}
QPushButton#primary_btn:disabled {
    background-color: #e2e8f0;
    color: #94a3b8;
    background: #e2e8f0;
}

/* 次要按钮 */
QPushButton#secondary_btn {
    background-color: #f1f5f9;
    color: #475569;
    border: 1px solid #cbd5e1;
    border-radius: 8px;
    font-weight: 600;
    font-size: 14px;
    padding: 12px 24px;
}
QPushButton#secondary_btn:hover {
    background-color: #e2e8f0;
    border-color: #94a3b8;
    color: #1e293b;
}
QPushButton#secondary_btn:pressed {
    background-color: #cbd5e1;
}

/* 图标按钮 */
QPushButton#icon_btn {
    background-color: white;
    color: #64748b;
    border: 1px solid #e2e8f0;
    border-radius: 8px;
    font-size: 18px;
    padding: 8px;
}
QPushButton#icon_btn:hover {
    background-color: #f1f5f9;
    border-color: #3b82f6;
    color: #1e293b;
}
QPushButton#icon_btn:pressed {
    background-color: #e2e8f0;
}

/* Highlight 按钮 */
QPushButton#highlight_btn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #3b82f6, stop:1 #2563eb);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    font-size: 14px;
    padding: 10px 20px;
}
QPushButton#highlight_btn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #60a5fa, stop:1 #3b82f6);
}
QPushButton#highlight_btn:pressed {
    background: #1d4ed8;
}

/* 刷新按钮 */
QPushButton#refresh_btn {
    background-color: white;
    color: #64748b;
    border: 1px solid #e2e8f0;
    border-radius: 8px;
    padding: 10px 14px;
    font-size: 12px;
    min-width: 70px;
}
QPushButton#refresh_btn:hover {
    background-color: #f1f5f9;
    color: #1e293b;
}

/* 主题切换按钮 */
QPushButton#theme_toggle_btn {
    background-color: white;
    color: #f39c12;
    border: 1px solid #e2e8f0;
    border-radius: 8px;
    padding: 6px;
    font-size: 20px;
    font-weight: normal;
}
QPushButton#theme_toggle_btn:hover {
    background-color: #f1f5f9;
    border-color: #3b82f6;
    color: #f1c40f;
}
QPushButton#theme_toggle_btn:pressed {
    background-color: #e2e8f0;
}

/* 底部导航按钮 */
QPushButton#bottom_nav_btn {
    background-color: transparent;
    color: #64748b;
    border: none;
    text-align: left;
    padding: 12px 16px;
    font-size: 14px;
    font-weight: 500;
    border-radius: 8px;
}
QPushButton#bottom_nav_btn:hover {
    background-color: #f1f5f9;
    color: #3b82f6;
}

/* 浏览按钮 */
QPushButton#browse_btn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #3b82f6, stop:1 #2563eb);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    min-width: 80px;
    padding: 10px 16px;
    font-size: 13px;
}
QPushButton#browse_btn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #60a5fa, stop:1 #3b82f6);
}

/* 保存按钮 - 绿色 */
QPushButton#save_btn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #10b981, stop:1 #059669);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    min-width: 80px;
    padding: 10px 16px;
    font-size: 13px;
}
QPushButton#save_btn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #34d399, stop:1 #10b981);
}

/* 清除按钮 */
QPushButton#clear_btn {
    background-color: #64748b;
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    padding: 10px 20px;
}
QPushButton#clear_btn:hover {
    background-color: #475569;
}

/* 关闭按钮 - 红色 */
QPushButton#close_btn {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #ef4444, stop:1 #dc2626);
    color: white;
    border: none;
    border-radius: 8px;
    font-weight: 600;
    min-width: 100px;
    padding: 10px 20px;
}
QPushButton#close_btn:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #f87171, stop:1 #ef4444);
}

/* 文本编辑器 - 纯白背景 */
QTextEdit {
    border: 1px solid #e2e8f0;
    border-radius: 8px;
    padding: 12px;
    background-color: white;
    font-family: 'SF Mono', Menlo, Monaco, 'Courier New', monospace;
    font-size: 13px;
    color: #1e293b;
    selection-background-color: #3b82f6;
}

/* 表格 - 更现代的设计 */
QTableWidget {
    border: 1px solid #e2e8f0;
    border-radius: 10px;
    background-color: white;
    gridline-color: #f1f5f9;
    selection-background-color: #3b82f6;
    selection-color: white;
    font-size: 14px;
}
QHeaderView::section {
    background-color: #f8fafc;
    color: #1e293b;
    padding: 10px 12px;
    border: none;
    border-bottom: 2px solid #e2e8f0;
    font-weight: 600;
    font-size: 13px;
}
QTableWidget::item {
    padding: 12px 10px;
    min-height: 40px;
    border-bottom: 1px solid #f1f5f9;
    color: #475569;
}
QTableWidget::item:hover {
    background-color: #f8fafc;
}
QTableWidget::item:selected {
    background-color: #3b82f6;
    color: white;
}
QTableWidget QTableCornerButton::section {
    background-color: #f8fafc;
    border: none;
    border-bottom: 2px solid #e2e8f0;
}

/* 复选框 */
QCheckBox {
    color: #475569;
    spacing: 10px;
    font-size: 14px;
    padding: 6px 0;
}
QCheckBox::indicator {
    width: 20px;
    height: 20px;
    border: 1.5px solid #cbd5e1;
    border-radius: 5px;
    background-color: white;
}
QCheckBox::indicator:hover {
    border-color: #3b82f6;
    background-color: #f8fafc;
}
QCheckBox::indicator:checked {
    background-color: #3b82f6;
    border-color: #3b82f6;
    image: none;
}
QCheckBox::indicator:checked:after {
    content: "✓";
    color: white;
    padding-left: 4px;
    font-weight: bold;
}

/* 进度条 */
QProgressBar {
    border: 1px solid #e2e8f0;
    border-radius: 8px;
    text-align: center;
    background-color: #f8fafc;
    height: 28px;
    color: #475569;
    font-size: 12px;
    font-weight: 500;
}
QProgressBar::chunk {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #3b82f6, stop:1 #2563eb);
    border-radius: 6px;
}

/* 标签页 - QTabWidget */
QTabWidget::pane {
    border: 1px solid #e2e8f0;
    border-radius: 8px;
    background-color: white;
    top: -1px;
}
QTabWidget > QWidget {
    background-color: white;
}
QTabBar::tab {
    background-color: #f8fafc;
    color: #64748b;
    border: 1px solid #e2e8f0;
    border-bottom: none;
    border-top-left-radius: 8px;
    border-top-right-radius: 8px;
    padding: 10px 20px;
    margin-right: 2px;
    font-size: 13px;
    font-weight: 500;
}
QTabBar::tab:selected {
    background-color: white;
    color: #1e293b;
    font-weight: 600;
    border-bottom: 3px solid #3b82f6;
}
QTabBar::tab:hover {
    background-color: #f1f5f9;
    color: #475569;
}

/* 滚动条 - 更精致 */
QScrollBar:vertical {
    border: none;
    background: #f8fafc;
    width: 14px;
    border-radius: 7px;
    margin: 2px;
}
QScrollBar::handle:vertical {
    background: #cbd5e1;
    border-radius: 7px;
    min-height: 40px;
}
QScrollBar::handle:vertical:hover {
    background: #94a3b8;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0;
}
QScrollBar:horizontal {
    border: none;
    background: #f8fafc;
    height: 14px;
    border-radius: 7px;
    margin: 2px;
}
QScrollBar::handle:horizontal {
    background: #cbd5e1;
    border-radius: 7px;
    min-width: 40px;
}
QScrollBar::handle:horizontal:hover {
    background: #94a3b8;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
    width: 0;
}
""")

    # ==============  New Workflow Tabs (Refactored) ==============\n    \n    def create_target_discovery_tab(self) -> QWidget:\n        \"\"\"Target Discovery Tab: G-Motif, Disease, Pocket\"\"\"\n        scroll_area = QScrollArea()\n        scroll_area.setWidgetResizable(True)\n        scroll_area.setFrameShape(QFrame.Shape.NoFrame)\n        \n        content_widget = QWidget()\n        content_widget.setObjectName(\"scroll_content\")\n        self._target_scroll_content = content_widget\n        bg_color = \"#0d1117\" if self._dark_mode else \"#f8fafc\"\n        content_widget.setStyleSheet(f\"#scroll_content {{ background-color: {bg_color}; }}\")\n        \n        layout = QVBoxLayout(content_widget)\n        layout.setSpacing(14)\n        layout.setContentsMargins(12, 12, 12, 12)\n        \n        # 1. Disease Analysis (if available)\n        if search_disease is not None:\n            try:\n                disease_card = self.create_disease_analysis_tab() \n                grp_disease = QGroupBox(\"Disease Target Analysis\")\n                d_layout = QVBoxLayout(grp_disease)\n                d_layout.addWidget(disease_card)\n                layout.addWidget(grp_disease)\n            except Exception as e:\n                self.log(f\"Failed to load Disease Analysis: {e}\")\n        \n        # 2. G-Motif Detection (Refactored from create_molecular_glue_tab)\n        grp_gm = QGroupBox(\"G-Motif (CRBN G-loop) Detection\")\n        gm_grid = QGridLayout(grp_gm)\n        gm_grid.setColumnStretch(1, 1); gm_grid.setColumnStretch(3, 1)\n        gm_grid.setHorizontalSpacing(8); gm_grid.setVerticalSpacing(10)\n        \n        # Row 0\n        gm_grid.addWidget(QLabel(\"Target Object:\"), 0, 0, Qt.AlignmentFlag.AlignRight)\n        self.obj_combo_gm = QComboBox(); self.obj_combo_gm.setMinimumHeight(32)\n        self.refresh_obj_gm = QPushButton(t(\"refresh\")); self.refresh_obj_gm.clicked.connect(self.refresh_objects)\n        r0 = QHBoxLayout(); r0.addWidget(self.obj_combo_gm, 1); r0.addWidget(self.refresh_obj_gm)\n        gm_grid.addLayout(r0, 0, 1)\n        \n        gm_grid.addWidget(QLabel(\"PDB File (opt):\"), 0, 2, Qt.AlignmentFlag.AlignRight)\n        self.gm_pdb = QLineEdit(); self.gm_pdb_browse = QPushButton(t(\"browse\"))\n        self.gm_pdb_browse.clicked.connect(self.browse_gm_pdb)\n        r0b = QHBoxLayout(); r0b.addWidget(self.gm_pdb, 1); r0b.addWidget(self.gm_pdb_browse)\n        gm_grid.addLayout(r0b, 0, 3)\n        \n        # Row 1\n        gm_grid.addWidget(QLabel(\"RMSD cutoff (Å):\"), 1, 0, Qt.AlignmentFlag.AlignRight)\n        self.gm_rmsd = QLineEdit(\"3.5\")\n        gm_grid.addWidget(self.gm_rmsd, 1, 1)\n        \n        self.gm_require_gly = QCheckBox(t(\"require_gly\")); self.gm_require_gly.setChecked(True)\n        gm_grid.addWidget(self.gm_require_gly, 1, 3)\n        \n        # Row 2\n        gm_grid.addWidget(QLabel(\"Template:\"), 2, 0, Qt.AlignmentFlag.AlignRight)\n        self.gm_template_mode = QComboBox()\n        self.gm_template_mode.addItems([\"Idealized (8×Cα)\", \"Built-in: GSPT1\", \"Built-in: CK1α\", \"Built-in: VAV1\", \"From Selection\"])\n        gm_grid.addWidget(self.gm_template_mode, 2, 1)\n        \n        gm_grid.addWidget(QLabel(\"Selection:\"), 2, 2, Qt.AlignmentFlag.AlignRight)\n        self.gm_template_sel = QLineEdit()\n        self.gm_template_pick = QPushButton(\"Pick (sele)\"); self.gm_template_pick.clicked.connect(lambda: self.gm_template_sel.setText(\"sele\"))\n        r2b = QHBoxLayout(); r2b.addWidget(self.gm_template_sel, 1); r2b.addWidget(self.gm_template_pick)\n        gm_grid.addLayout(r2b, 2, 3)\n        \n        def _toggle_template_inputs(idx):\n            use_sel = (idx == 4)\n            self.gm_template_sel.setEnabled(use_sel); self.gm_template_pick.setEnabled(use_sel)\n        self.gm_template_mode.currentIndexChanged.connect(_toggle_template_inputs)\n        _toggle_template_inputs(0)\n        \n        # Row 3\n        gm_grid.addWidget(QLabel(\"Output CSV:\"), 3, 0, Qt.AlignmentFlag.AlignRight)\n        self.gm_out_csv = QLineEdit()\n        self.gm_out_browse = QPushButton(t(\"browse\")); self.gm_out_browse.clicked.connect(self.browse_gm_out_csv)\n        r3 = QHBoxLayout(); r3.addWidget(self.gm_out_csv, 1); r3.addWidget(self.gm_out_browse)\n        gm_grid.addLayout(r3, 3, 1, 1, 3)\n        \n        # Buttons\n        gm_btn_row = QHBoxLayout()\n        self.gm_btn = QPushButton(\"Detect POI\"); self.gm_btn.setObjectName(\"highlight_btn\")\n        self.gm_btn.clicked.connect(self.start_gmotif)\n        self.gm_btn_render = QPushButton(\"Render All (POI + ESP + PNG)\"); self.gm_btn_render.setObjectName(\"highlight_btn\")\n        self.gm_btn_render.clicked.connect(self.render_gmotif_with_esp)\n        gm_btn_row.addWidget(self.gm_btn); gm_btn_row.addWidget(self.gm_btn_render); gm_btn_row.addStretch(1)\n        \n        layout.addWidget(grp_gm)\n        layout.addLayout(gm_btn_row)\n        \n        # 3. Pocket Detection (Existing logic)\n        layout.addWidget(self._create_pocket_detection_card())\n        \n        # 4. Advanced Pocket Analysis (Existing logic)\n        layout.addWidget(self._create_advanced_pocket_card())\n        \n        layout.addStretch(1)\n        scroll_area.setWidget(content_widget)\n        \n        wrapper = QWidget()\n        wl = QVBoxLayout(wrapper); wl.setContentsMargins(0,0,0,0); wl.addWidget(scroll_area)\n        return wrapper\n\n    def create_hit_identification_tab(self) -> QWidget:\n        \"\"\"Hit Identification Tab: Vina Docking\"\"\"\n        scroll_area = QScrollArea()\n        scroll_area.setWidgetResizable(True)\n        scroll_area.setFrameShape(QFrame.Shape.NoFrame)\n        \n        content_widget = QWidget()\n        content_widget.setObjectName(\"scroll_content\")\n        bg_color = \"#0d1117\" if self._dark_mode else \"#f8fafc\"\n        content_widget.setStyleSheet(f\"#scroll_content {{ background-color: {bg_color}; }}\")\n        \n        layout = QVBoxLayout(content_widget)\n        layout.setSpacing(14)\n        layout.setContentsMargins(12, 12, 12, 12)\n        \n        # Vina Docking Card\n        layout.addWidget(self._create_vina_docking_card())\n        \n        # Placeholder for future Virtual Screening\n        grp_vs = QGroupBox(\"Virtual Screening (Coming Soon)\")\n        vs_layout = QVBoxLayout(grp_vs)\n        vs_layout.addWidget(QLabel(\"Batch docking and scoring functionality will be available in future updates.\"))\n        layout.addWidget(grp_vs)\n        \n        layout.addStretch(1)\n        scroll_area.setWidget(content_widget)\n        \n        wrapper = QWidget()\n        wl = QVBoxLayout(wrapper); wl.setContentsMargins(0,0,0,0); wl.addWidget(scroll_area)\n        return wrapper\n        \n    def create_lead_optimization_tab(self) -> QWidget:\n        \"\"\"Lead Optimization: PPI, Glue, Ternary, Mutation\"\"\"\n        scroll_area = QScrollArea()\n        scroll_area.setWidgetResizable(True)\n        scroll_area.setFrameShape(QFrame.Shape.NoFrame)\n        \n        content_widget = QWidget()\n        content_widget.setObjectName(\"scroll_content\")\n        self._lead_scroll_content = content_widget\n        bg_color = \"#0d1117\" if self._dark_mode else \"#f8fafc\"\n        content_widget.setStyleSheet(f\"#scroll_content {{ background-color: {bg_color}; }}\")\n        \n        layout = QVBoxLayout(content_widget)\n        layout.setSpacing(14)\n        layout.setContentsMargins(12, 12, 12, 12)\n        \n        # 1. Molecular Glue Specifics\n        grp_glue = QGroupBox(\"Molecular Glue Analysis (PPI + Neo-Epitope)\")\n        glue_grid = QGridLayout(grp_glue)\n        glue_grid.setColumnStretch(1, 1); glue_grid.setColumnStretch(3, 1)\n        glue_grid.setHorizontalSpacing(8); glue_grid.setVerticalSpacing(10)\n        \n        glue_grid.addWidget(QLabel(\"Target Object:\"), 0, 0, Qt.AlignmentFlag.AlignRight)\n        self.glue_obj_combo = QComboBox(); self.glue_obj_combo.setMinimumHeight(36)\n        self.glue_refresh_btn = QPushButton(t(\"refresh\")); self.glue_refresh_btn.clicked.connect(self.refresh_objects)\n        r0 = QHBoxLayout(); r0.addWidget(self.glue_obj_combo, 1); r0.addWidget(self.glue_refresh_btn)\n        glue_grid.addLayout(r0, 0, 1)\n        \n        glue_grid.addWidget(QLabel(\"Glue Resname:\"), 0, 2, Qt.AlignmentFlag.AlignRight)\n        self.glue_resname = QLineEdit(); self.glue_resname.setPlaceholderText(\"e.g. CC885\")\n        glue_grid.addWidget(self.glue_resname, 0, 3)\n        \n        glue_grid.addWidget(QLabel(\"E3 Chains:\"), 1, 0, Qt.AlignmentFlag.AlignRight)\n        self.glue_e3_chains = QLineEdit(); self.glue_e3_chains.setPlaceholderText(\"e.g. A\")\n        glue_grid.addWidget(self.glue_e3_chains, 1, 1)\n        \n        glue_grid.addWidget(QLabel(\"Substrate Chains:\"), 1, 2, Qt.AlignmentFlag.AlignRight)\n        self.glue_sub_chains = QLineEdit(); self.glue_sub_chains.setPlaceholderText(\"e.g. B\")\n        glue_grid.addWidget(self.glue_sub_chains, 1, 3)\n        \n        glue_grid.addWidget(QLabel(\"Interface Dist (Å):\"), 2, 0, Qt.AlignmentFlag.AlignRight)\n        self.glue_interface_dist = QLineEdit(\"4.5\")\n        glue_grid.addWidget(self.glue_interface_dist, 2, 1)\n        \n        glue_grid.addWidget(QLabel(\"Neo-Epitope Dist (Å):\"), 2, 2, Qt.AlignmentFlag.AlignRight)\n        self.glue_neo_dist = QLineEdit(\"5.0\")\n        glue_grid.addWidget(self.glue_neo_dist, 2, 3)\n        \n        glue_btn_row = QHBoxLayout()\n        self.glue_full_btn = QPushButton(\"Full Glue Analysis\"); self.glue_full_btn.setObjectName(\"highlight_btn\")\n        self.glue_full_btn.clicked.connect(self.run_glue_full_analysis)\n        glue_btn_row.addWidget(self.glue_full_btn)\n        glue_btn_row.addStretch(1)\n        \n        layout.addWidget(grp_glue)\n        layout.addLayout(glue_btn_row)\n        \n        # 2. Ternary Complex\n        grp_ternary = QGroupBox(\"Ternary Complex Analysis\")\n        t_grid = QGridLayout(grp_ternary)\n        t_grid.setColumnStretch(1, 1); t_grid.setColumnStretch(3, 1)\n        \n        t_grid.addWidget(QLabel(\"Target Object:\"), 0, 0)\n        self.tc_obj_combo = QComboBox(); \n        t_grid.addWidget(self.tc_obj_combo, 0, 1)\n        \n        t_grid.addWidget(QLabel(\"Ligand:\"), 0, 2)\n        self.tc_ligand_name = QLineEdit(); self.tc_ligand_name.setPlaceholderText(\"Auto\")\n        t_grid.addWidget(self.tc_ligand_name, 0, 3)\n        \n        t_grid.addWidget(QLabel(\"E3 Chains:\"), 1, 0)\n        self.tc_protein1_chains = QLineEdit(); self.tc_protein1_chains.setPlaceholderText(\"e.g. A\")\n        t_grid.addWidget(self.tc_protein1_chains, 1, 1)\n        \n        t_grid.addWidget(QLabel(\"POI Chains:\"), 1, 2)\n        self.tc_protein2_chains = QLineEdit(); self.tc_protein2_chains.setPlaceholderText(\"e.g. B\")\n        t_grid.addWidget(self.tc_protein2_chains, 1, 3)\n        \n        t_btn_row = QHBoxLayout()\n        self.tc_analyze_btn = QPushButton(\"Analyze Complex\"); self.tc_analyze_btn.setObjectName(\"highlight_btn\")\n        self.tc_analyze_btn.clicked.connect(self.run_tc_analysis)\n        self.tc_render_btn = QPushButton(\"Render All\"); self.tc_render_btn.setObjectName(\"highlight_btn\")\n        self.tc_render_btn.clicked.connect(self.run_tc_render)\n        t_btn_row.addWidget(self.tc_analyze_btn); t_btn_row.addWidget(self.tc_render_btn); t_btn_row.addStretch(1)\n        \n        layout.addWidget(grp_ternary)\n        layout.addLayout(t_btn_row)\n        \n        # 3. General Interactions (PP, PL)\n        grp_pl = QGroupBox(\"Protein-Ligand Interactions\")\n        pl_layout = QHBoxLayout(grp_pl)\n        self.pl_obj_combo = QComboBox(); self.pl_obj_combo.setMinimumWidth(150)\n        self.pl_analyze_btn = QPushButton(\"Analyze PL\"); self.pl_analyze_btn.clicked.connect(self.run_pl_analysis)\n        pl_layout.addWidget(QLabel(\"Object:\")); pl_layout.addWidget(self.pl_obj_combo)\n        pl_layout.addWidget(self.pl_analyze_btn)\n        layout.addWidget(grp_pl)\n        \n        # 4. Mutation Analysis\n        layout.addWidget(self._create_mutation_analysis_card())\n        \n        layout.addStretch(1)\n        scroll_area.setWidget(content_widget)\n        \n        wrapper = QWidget()\n        wl = QVBoxLayout(wrapper); wl.setContentsMargins(0,0,0,0); wl.addWidget(scroll_area)\n        return wrapper\n        \n    def create_visualization_tab(self) -> QWidget:\n        \"\"\"Visualization: APBS, Export\"\"\"\n        return self.create_apbs_tab()\n\n# -------- 将 setup_style 绑定到 GLINTDialog --------
GLINTDialog.setup_style = _glint_setup_style
# -------- 独立运行入口（可选）--------
def launch_standalone():
    app = QApplication.instance() or QApplication(sys.argv)
    dlg = GLINTDialog()
    try:
        sys.exit(app.exec())
    except AttributeError:
        sys.exit(app.exec_())

if __name__ == "__main__":
    launch_standalone()
