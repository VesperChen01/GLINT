# -*- coding: utf-8 -*-
"""
MolStruct 统一 GUI
- 浅色主题 + 中英自动切换
- 表格化结果
- G-Motif 识别（置于首个标签）
- 相互作用分析
- 静电势（APBS/Quick）+ PNG(300dpi)/DX 导出
- G-Motif 一键渲染（电势+高亮+出图）

说明：
1) 本文件与 __init__.py 配合使用（GUI 非模态，避免“未响应”）
2) Quick 电势不依赖 APBS；APBS 真解算可选，需配置 apbs_tools
"""

from __future__ import annotations
import os, sys, csv, importlib.util
from typing import List, Dict, Any, Tuple

# -------- Qt 兼容（优先 PyQt5）--------
QT_LIB = None
try:
    from PyQt5.QtCore import Qt, QThread, pyqtSignal, QLocale, QTimer
    from PyQt5.QtWidgets import (
        QApplication, QDialog, QWidget, QVBoxLayout, QHBoxLayout, QTextEdit, QLabel,
        QLineEdit, QPushButton, QCheckBox, QComboBox, QFileDialog, QGroupBox,
        QFormLayout, QMessageBox, QTextEdit, QProgressBar, QFrame, QTabWidget,
        QTableWidget, QTableWidgetItem, QSizePolicy, QGridLayout, QListWidget, QStackedWidget
    )
    QT_LIB = "PyQt5"
except Exception:
    try:
        from PyQt6.QtCore import Qt, QThread, pyqtSignal, QLocale, QTimer
        from PyQt6.QtWidgets import (
            QApplication, QDialog, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
            QLineEdit, QPushButton, QCheckBox, QComboBox, QFileDialog, QGroupBox,
            QFormLayout, QMessageBox, QTextEdit, QProgressBar, QFrame, QTabWidget,
            QTableWidget, QTableWidgetItem, QSizePolicy, QGridLayout, QListWidget, QStackedWidget
        )
        QT_LIB = "PyQt6"
    except Exception as e:
        raise RuntimeError("需要安装 PyQt5 或 PyQt6") from e

# -------- 模块导入助手（相对→绝对→动态）--------
def _dynamic_load_by_filenames(names: list[str], symbol: str):
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
    """返回：highlight_csv_residues, highlight_gmotif_loops, analyze_pdb_interactions, find_crbn_g_motif, render_interactions_beautifully, generate_2d_interaction_diagram, analyze_protein_ligand_interactions, visualize_protein_ligand_3d, generate_interaction_network_plot, analyze_ternary_complex, analyze_atom_pair_interactions, visualize_atom_pairs"""
    # 包内尝试
    try:
        from .highlight_residues import highlight_csv_residues, highlight_gmotif_loops  # type: ignore
        from .interaction_analyzer import (  # type: ignore
            analyze_pdb_interactions, render_interactions_beautifully,
            analyze_protein_ligand_interactions, visualize_protein_ligand_3d,
            generate_interaction_network_plot, analyze_ternary_complex,
            analyze_atom_pair_interactions, visualize_atom_pairs
        )
        from .interaction_2d_plot import generate_2d_interaction_diagram  # type: ignore
        try:
            try:
                from .g_motif_analyzer import find_crbn_g_motif  # type: ignore
            except Exception:
                from .g_motif import find_crbn_g_motif  # type: ignore
        except Exception:
            find_crbn_g_motif = _dynamic_load_by_filenames(
                ["g_motif_analyzer.py", "g_motif.py", "g-motif.py"], "find_crbn_g_motif"
            )
        return highlight_csv_residues, highlight_gmotif_loops, analyze_pdb_interactions, find_crbn_g_motif, render_interactions_beautifully, generate_2d_interaction_diagram, analyze_protein_ligand_interactions, visualize_protein_ligand_3d, generate_interaction_network_plot, analyze_ternary_complex, analyze_atom_pair_interactions, visualize_atom_pairs
    except Exception:
        pass
    # 同目录绝对
    here = os.path.dirname(os.path.abspath(__file__))
    if here and here not in sys.path:
        sys.path.insert(0, here)
    try:
        from highlight_residues import highlight_csv_residues, highlight_gmotif_loops  # type: ignore
        from interaction_analyzer import (  # type: ignore
            analyze_pdb_interactions, render_interactions_beautifully,
            analyze_protein_ligand_interactions, visualize_protein_ligand_3d,
            generate_interaction_network_plot, analyze_ternary_complex,
            analyze_atom_pair_interactions, visualize_atom_pairs
        )
        from interaction_2d_plot import generate_2d_interaction_diagram  # type: ignore
    except Exception as e:
        raise ModuleNotFoundError(
            "未找到 highlight_residues / interaction_analyzer；请确认与 unified_gui.py 同目录或包内存在。"
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
    return highlight_csv_residues, highlight_gmotif_loops, analyze_pdb_interactions, find_crbn_g_motif, render_interactions_beautifully, generate_2d_interaction_diagram, analyze_protein_ligand_interactions, visualize_protein_ligand_3d, generate_interaction_network_plot, analyze_ternary_complex, analyze_atom_pair_interactions, visualize_atom_pairs

highlight_csv_residues, highlight_gmotif_loops, analyze_pdb_interactions, find_crbn_g_motif, render_interactions_beautifully, generate_2d_interaction_diagram, analyze_protein_ligand_interactions, visualize_protein_ligand_3d, generate_interaction_network_plot, analyze_ternary_complex, analyze_atom_pair_interactions, visualize_atom_pairs = _import_helpers()

# -------- 依赖检查 --------
def _check_and_install_deps():
    """检查并安装依赖，GUI启动时调用"""
    try:
        # 尝试导入env_setup模块
        here = os.path.dirname(os.path.abspath(__file__))
        sys.path.insert(0, here)
        try:
            from .env_setup import ensure_dependencies, get_dependency_status
        except ImportError:
            from env_setup import ensure_dependencies, get_dependency_status
        
        # 检查依赖状态
        status = get_dependency_status()
        missing = [pkg for pkg, avail in status.items() if not avail and pkg in ['rdkit', 'scipy', 'matplotlib', 'pillow', 'numpy']]
        
        if missing:
            print(f"[MolStruct GUI] 检测到缺失依赖: {', '.join(missing)}")
            print("[MolStruct GUI] 正在自动安装...")
            success = ensure_dependencies()
            if not success:
                print("[MolStruct GUI] ⚠️ 部分依赖安装失败，请手动安装")
                return False
        return True
    except Exception as e:
        print(f"[MolStruct GUI] 依赖检查失败: {e}")
        return False

# -------- Language & Text --------
LANG_FORCE = "en"  # Force English for all GUI
def get_lang() -> str:
    return "en"  # Always return English

T = {
    "title": {"zh": "MolStruct 统一 GUI", "en": "MolStruct Unified GUI"},
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
    "between_chains": {"zh": "仅不同链之间", "en": "Only between chains"},
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

# -------- 工作线程 --------
class AnalysisWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(list)
    error = pyqtSignal(str)
    def __init__(self, obj_name: str, pdb_file: str | None, only_between_chains: bool, output_csv: str | None):
        super().__init__()
        self.obj_name = obj_name
        self.pdb_file = pdb_file
        self.only_between_chains = only_between_chains
        self.output_csv = output_csv
    def run(self):
        try:
            self.progress.emit(t("log_start"))
            interactions = analyze_pdb_interactions(
                obj_name=self.obj_name,
                pdb_file=self.pdb_file,
                only_between_chains=self.only_between_chains,
                output_csv=self.output_csv,
                auto_highlight=False,
            )
            self.progress.emit(t("log_done").format(n=len(interactions)))
            self.finished.emit(interactions)
        except Exception as e:
            self.error.emit(str(e))

class GMotifWorker(QThread):
    progress = pyqtSignal(str)
    finished = pyqtSignal(list, str)
    error = pyqtSignal(str)
    def __init__(self, obj_name: str, pdb_file: str | None, rmsd: float, require_gly: bool,
                 out_csv: str | None, template_mode: str, template_sel: str | None, template_builtin: str | None):
        super().__init__()
        self.obj_name = obj_name; self.pdb_file = pdb_file
        self.rmsd = rmsd; self.require_gly = require_gly; self.out_csv = out_csv
        self.template_mode = template_mode; self.template_sel = template_sel; self.template_builtin = template_builtin
    def run(self):
        try:
            if find_crbn_g_motif is None:
                raise RuntimeError("未找到 find_crbn_g_motif，请确认 g_motif_analyzer.py 在插件目录中。")
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
                auto_highlight=0,
                require_gly_pos6=bool(self.require_gly),
            ) or []
            self.progress.emit("[G-Motif] " + (f"完成，命中 {len(hits)} 条" if get_lang()=="zh" else f"Done, {len(hits)} hits"))
            self.finished.emit(hits, out_csv_path)
        except Exception as e:
            self.error.emit(str(e))

# -------- 主对话框 --------
class MolStructDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(t("title"))
        # 固定窗口尺寸（锁定大小，不可调整）
        fixed_w = 1280
        fixed_h = 720
        self.setFixedSize(fixed_w, fixed_h)  # 锁定窗口大小
        
        # 居中显示
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
                x = (geom.width() - fixed_w) // 2
                y = (geom.height() - fixed_h) // 2
                self.move(x, y)

        self.analysis_thread: AnalysisWorker | None = None
        self.gmotif_thread: GMotifWorker | None = None

        self._interactions: List[Dict[str, Any]] = []
        self._gmotif_hits: List[Tuple] = []
        self._last_gmotif_csv: str | None = None
        self._esp_maps: Dict[str, Tuple[str, str]] = {}  # obj -> (map, ramp)
        self._dark_mode: bool = True  # 默认深色主题
        self._ui_scale: float = 1.0   # 自动缩放比例

        # 检查并安装依赖
        _check_and_install_deps()
        
        self.build_ui()
        self.setup_style()
        # 初始化主题图标
        self.theme_toggle_btn.setText("🌙" if self._dark_mode else "☀️")
        # 禁用自动缩放以保持固定高度
        # self.apply_auto_scaling()

        # ❗关键修复：延后首次 PyMOL 调用，避免构造期阻塞
        QTimer.singleShot(0, self.refresh_objects)

        self.update_enablement()
        self.update_modules_button_style()  # Initialize Modules button color
        self.log(t("log_ready"))

    # --- UI 结构 ---
    def build_ui(self):
        # ========== 主布局 ==========
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(15, 15, 15, 15)  # 恢复到舒适的外边距
        main_layout.setSpacing(12)  # 恢复到舒适的间距

        # 水平分割：导航 | 内容 | 结果
        content_row = QHBoxLayout()
        content_row.setSpacing(10)  # 减小间距

        # ========== Left: Navigation List ==========
        nav_widget = QWidget()
        nav_layout = QVBoxLayout(nav_widget)
        nav_layout.setSpacing(6)  # 减小间距
        nav_layout.setContentsMargins(0, 0, 0, 0)

        # Header: Theme icon + Modules label (centered)
        nav_header = QHBoxLayout()
        nav_header.setSpacing(4)
        
        # Theme toggle icon (sun/moon) - no border
        self.theme_toggle_btn = QPushButton("🌙")  # Moon icon for dark mode by default
        self.theme_toggle_btn.setObjectName("theme_icon_btn")
        self.theme_toggle_btn.setFlat(True)
        self.theme_toggle_btn.setFixedSize(28, 28)
        self.theme_toggle_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.theme_toggle_btn.setToolTip("Toggle theme")
        self.theme_toggle_btn.clicked.connect(self.toggle_theme)
        
        # Modules label (centered, borderless)
        self.modules_btn = QPushButton("Modules")
        self.modules_btn.setObjectName("modules_btn")
        self.modules_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self.modules_btn.setFlat(True)
        self.modules_btn.setToolTip("Modules")
        self.modules_btn.setMinimumWidth(120)  # 设置最小宽度以显示完整文本
        
        nav_header.addStretch()
        nav_header.addWidget(self.theme_toggle_btn)
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

        nav_items = [
            ("Welcome", "Welcome"),
            ("POI Discovery", "POI Discovery"),
            ("Interaction Analysis", "Interaction Analysis"),
            ("Electrostatics", "Electrostatics"),
        ]

        for zh_text, en_text in nav_items:
            text = zh_text if get_lang() == "zh" else en_text
            self.nav_list.addItem(text)

        self.nav_list.setCurrentRow(0)
        self.nav_list.currentRowChanged.connect(self.on_nav_changed)
        nav_layout.addWidget(self.nav_list, 1)
        
        # Bottom items (Check Environment, README and Contact)
        from PyQt5.QtWidgets import QFrame
        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setFrameShadow(QFrame.Shadow.Sunken)
        separator.setStyleSheet("background-color: #e2e8f0; margin: 8px 0;")
        nav_layout.addWidget(separator)
        
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
        readme_btn.clicked.connect(lambda: self.content_stack.setCurrentIndex(4))
        nav_layout.addWidget(readme_btn)
        
        # Contact button
        contact_btn = QPushButton("Contact Us")
        contact_btn.setObjectName("bottom_nav_btn")
        contact_btn.setFlat(True)
        contact_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        contact_btn.clicked.connect(lambda: self.content_stack.setCurrentIndex(5))
        nav_layout.addWidget(contact_btn)

        # ========== 中间：内容堆栈 ==========
        self.content_stack = QStackedWidget()

        # 创建各个页面
        self.content_stack.addWidget(self.create_welcome_tab())         # 0: Welcome
        self.content_stack.addWidget(self.create_molecular_glue_tab())  # 1: POI Discovery (G-Motif + Ternary Complex)
        self.content_stack.addWidget(self.create_interaction_tab())     # 2: Interaction Analysis (Prot-Prot + Prot-Lig + Atom)
        self.content_stack.addWidget(self.create_apbs_tab())            # 3: Electrostatics
        self.content_stack.addWidget(self.create_readme_tab())          # 4: README
        self.content_stack.addWidget(self.create_contact_tab())         # 5: Contact

        # ========== 右侧：结果与日志 ==========
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setSpacing(6)  # 减小间距
        right_layout.setContentsMargins(0, 0, 0, 0)

        grp = QGroupBox(t("right_results"))
        grp_layout = QVBoxLayout(grp)
        grp_layout.setSpacing(6)  # 减小间距

        self.table = QTableWidget()
        self.table.setColumnCount(6)
        self.table.setHorizontalHeaderLabels(T["table_header"][get_lang()])
        self.table.setSortingEnabled(True)
        self.table.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        # 行距与可读性修复（按比例）
        try:
            self.table.setAlternatingRowColors(True)
            vh = self.table.verticalHeader()
            base = 28
            scale = getattr(self, "_ui_scale", 1.0)
            vh.setDefaultSectionSize(int(base * scale))
            self.table.setWordWrap(False)
        except Exception:
            pass
        grp_layout.addWidget(self.table)

        self.log_edit = QTextEdit()
        self.log_edit.setReadOnly(True)
        self.log_edit.setMinimumHeight(100)  # 减小日志框高度
        grp_layout.addWidget(self.log_edit)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        grp_layout.addWidget(self.progress_bar)

        right_layout.addWidget(grp)
        
        # Store right widget reference for show/hide
        self.right_widget = right_widget

        # ========== Assemble horizontal layout (wider nav, balanced center, narrow right) ==========
        content_row.addWidget(nav_widget, 0)  # Left nav: auto-size to content
        nav_widget.setMaximumWidth(240)  # Further increased max width for nav panel to show full text
        nav_widget.setMinimumWidth(220)  # Further increased min width for nav panel
        content_row.addWidget(self.content_stack, 4)  # Center content: 4 parts (reduced from 6)
        content_row.addWidget(right_widget, 1)  # Right results: 1 part (窄结果栏)
        right_widget.setMaximumWidth(320)  # Slightly increased for better readability
        right_widget.setMinimumWidth(280)  # Adjusted min width

        # ========== Bottom: Empty (no status bar needed) ==========
        btn_layout = QHBoxLayout()
        btn_layout.setSpacing(8)  # 减小间距
        btn_layout.setContentsMargins(0, 6, 0, 0)  # 减小上边距
        btn_layout.addStretch(1)

        # ========== 添加到主布局 ==========
        main_layout.addLayout(content_row)
        main_layout.addLayout(btn_layout)

    def on_nav_changed(self, index):
        """导航切换"""
        self.content_stack.setCurrentIndex(index)
        # Hide right panel (CSV/log) when on Welcome page (index 0)
        if hasattr(self, 'right_widget'):
            self.right_widget.setVisible(index != 0)

    def create_interaction_tab(self) -> QWidget:
        """创建整合的相互作用分析标签页（Protein-Protein + Protein-Ligand + Atom Pairs）"""
        w = QWidget()
        main_layout = QVBoxLayout(w)
        main_layout.setSpacing(6)  # 大幅减小间距
        main_layout.setContentsMargins(6, 6, 6, 6)  # 减小边距
        
        # ========== 1. Protein-Protein Interaction ==========
        grp_pp = QGroupBox("Protein-Protein Interaction Analysis")
        # 使用网格布局实现两栏式
        pp_grid = QGridLayout(grp_pp)
        pp_grid.setColumnStretch(0, 0)
        pp_grid.setColumnStretch(1, 1)
        pp_grid.setColumnStretch(2, 0)
        pp_grid.setColumnStretch(3, 1)
        pp_grid.setHorizontalSpacing(8)  # 恢复水平间距
        pp_grid.setVerticalSpacing(8)  # 恢复垂直间距
        
        # 第一行: Target Object | [combo + refresh] | PDB File | [input + browse]
        pp_grid.addWidget(QLabel("Target Object"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.obj_combo_analysis = QComboBox()
        self.obj_combo_analysis.setMinimumHeight(26)
        self.refresh_obj_analysis = QPushButton(t("refresh"))
        self.refresh_obj_analysis.setObjectName("refresh_btn")
        self.refresh_obj_analysis.setMinimumHeight(26)
        self.refresh_obj_analysis.clicked.connect(self.refresh_objects)
        row0_container = QWidget()
        row0 = QHBoxLayout(row0_container)
        row0.setContentsMargins(0, 0, 0, 0)
        row0.addWidget(self.obj_combo_analysis, 1)
        row0.addWidget(self.refresh_obj_analysis)
        pp_grid.addWidget(row0_container, 0, 1)
        
        pp_grid.addWidget(QLabel("PDB File (optional)"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.pdb_path = QLineEdit()
        self.pdb_path.setMinimumHeight(26)
        self.pdb_browse = QPushButton(t("browse"))
        self.pdb_browse.setMinimumHeight(26)
        self.pdb_browse.setObjectName("browse_btn")
        self.pdb_browse.clicked.connect(self.browse_pdb)
        row1_container = QWidget()
        row1 = QHBoxLayout(row1_container)
        row1.setContentsMargins(0, 0, 0, 0)
        row1.addWidget(self.pdb_path, 1)
        row1.addWidget(self.pdb_browse)
        pp_grid.addWidget(row1_container, 0, 3)
        
        # 第二行: Only between chains | [checkbox] | Output CSV | [input + browse]
        pp_grid.addWidget(QLabel(""), 1, 0)  # 空位
        self.chk_between = QCheckBox(t("between_chains"))
        pp_grid.addWidget(self.chk_between, 1, 1)
        
        pp_grid.addWidget(QLabel("Output CSV (optional)"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.out_csv = QLineEdit()
        self.out_csv.setMinimumHeight(26)
        self.out_browse = QPushButton(t("browse"))
        self.out_browse.setMinimumHeight(26)
        self.out_browse.setObjectName("save_btn")
        self.out_browse.clicked.connect(self.browse_out_csv)
        row2_container = QWidget()
        row2 = QHBoxLayout(row2_container)
        row2.setContentsMargins(0, 0, 0, 0)
        row2.addWidget(self.out_csv, 1)
        row2.addWidget(self.out_browse)
        pp_grid.addWidget(row2_container, 1, 3)
        
        btn_pp_row = QHBoxLayout()
        self.analyze_btn = QPushButton("Analyze")
        self.analyze_btn.setObjectName("highlight_btn")
        self.analyze_btn.clicked.connect(self.start_analysis)
        
        self.render_interact_btn = QPushButton("Render (Analyze + Beautify + PNG)")
        self.render_interact_btn.setObjectName("highlight_btn")
        self.render_interact_btn.clicked.connect(self.render_interactions_beautifully_clicked)
        
        btn_pp_row.addWidget(self.analyze_btn)
        btn_pp_row.addWidget(self.render_interact_btn)
        btn_pp_row.addStretch(1)
        
        main_layout.addWidget(grp_pp)
        main_layout.addLayout(btn_pp_row)
        
        # ========== 2. Protein-Ligand Interaction ==========
        grp_pl = QGroupBox("Protein-Ligand Interaction Analysis")
        # 使用网格布局实现两栏式
        pl_grid = QGridLayout(grp_pl)
        pl_grid.setColumnStretch(0, 0)
        pl_grid.setColumnStretch(1, 1)
        pl_grid.setColumnStretch(2, 0)
        pl_grid.setColumnStretch(3, 1)
        pl_grid.setHorizontalSpacing(8)  # 恢复水平间距
        pl_grid.setVerticalSpacing(8)  # 恢复垂直间距
        
        # 第一行: Target Object | [combo + refresh] | Ligand Resname | [input]
        pl_grid.addWidget(QLabel("Target Object"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.pl_obj_combo = QComboBox()
        self.pl_obj_combo.setMinimumHeight(26)
        self.pl_refresh_btn = QPushButton(t("refresh"))
        self.pl_refresh_btn.setObjectName("refresh_btn")
        self.pl_refresh_btn.setMinimumHeight(26)
        self.pl_refresh_btn.clicked.connect(self.refresh_objects)
        pl_obj_row_container = QWidget()
        pl_obj_row = QHBoxLayout(pl_obj_row_container)
        pl_obj_row.setContentsMargins(0, 0, 0, 0)
        pl_obj_row.addWidget(self.pl_obj_combo, 1)
        pl_obj_row.addWidget(self.pl_refresh_btn)
        pl_grid.addWidget(pl_obj_row_container, 0, 1)
        
        pl_grid.addWidget(QLabel("Ligand Resname:"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.pl_ligand_name = QLineEdit()
        self.pl_ligand_name.setMinimumHeight(26)
        self.pl_ligand_name.setPlaceholderText("Auto-detect if blank")
        pl_grid.addWidget(self.pl_ligand_name, 0, 3)
        
        # 第二行: Protein Chains | [input] | Distance cutoff | [input]
        pl_grid.addWidget(QLabel("Protein Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.pl_protein_chains = QLineEdit()
        self.pl_protein_chains.setMinimumHeight(26)
        self.pl_protein_chains.setPlaceholderText("Auto-detect if blank")
        pl_grid.addWidget(self.pl_protein_chains, 1, 1)
        
        pl_grid.addWidget(QLabel("Distance cutoff (Å):"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.pl_distance = QLineEdit("4.5")
        self.pl_distance.setMinimumHeight(26)
        pl_grid.addWidget(self.pl_distance, 1, 3)
        
        # 第三行: Output CSV | [input + browse] (跨两列)
        pl_grid.addWidget(QLabel("Output CSV (optional)"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.pl_csv = QLineEdit()
        self.pl_csv.setMinimumHeight(26)
        self.pl_csv.setPlaceholderText("Optional")
        self.pl_csv_btn = QPushButton(t("browse"))
        self.pl_csv_btn.setMinimumHeight(26)
        self.pl_csv_btn.setObjectName("browse_btn")
        self.pl_csv_btn.clicked.connect(lambda: self._browse_save_file(self.pl_csv, "CSV (*.csv)"))
        pl_csv_row = QHBoxLayout()
        pl_csv_row.addWidget(self.pl_csv, 1)
        pl_csv_row.addWidget(self.pl_csv_btn)
        pl_grid.addLayout(pl_csv_row, 2, 1, 1, 3)  # 跨三列
        
        btn_pl_row = QHBoxLayout()
        self.pl_analyze_btn = QPushButton("Analyze")
        self.pl_analyze_btn.setObjectName("highlight_btn")
        self.pl_analyze_btn.clicked.connect(self.run_pl_analysis)
        
        self.pl_visualize_btn = QPushButton("3D Visualize")
        self.pl_visualize_btn.setObjectName("highlight_btn")
        self.pl_visualize_btn.clicked.connect(self.run_pl_visualize)
        
        self.pl_network_btn = QPushButton("Network Plot")
        self.pl_network_btn.setObjectName("highlight_btn")
        self.pl_network_btn.clicked.connect(self.run_pl_network)
        
        btn_pl_row.addWidget(self.pl_analyze_btn)
        btn_pl_row.addWidget(self.pl_visualize_btn)
        btn_pl_row.addWidget(self.pl_network_btn)
        btn_pl_row.addStretch(1)
        
        main_layout.addWidget(grp_pl)
        main_layout.addLayout(btn_pl_row)
        
        # ========== 3. Atom Pair Analysis ==========
        grp_ap = QGroupBox("Atom Pair Analysis (Atomic-Level Precision)")
        # 使用网格布局实现两栏式
        ap_grid = QGridLayout(grp_ap)
        ap_grid.setColumnStretch(0, 0)
        ap_grid.setColumnStretch(1, 1)
        ap_grid.setColumnStretch(2, 0)
        ap_grid.setColumnStretch(3, 1)
        ap_grid.setHorizontalSpacing(8)  # 恢复水平间距
        ap_grid.setVerticalSpacing(8)  # 恢复垂直间距
        
        # 第一行: Target Object | [combo + refresh] | Atom1 Selection | [input]
        ap_grid.addWidget(QLabel("Target Object"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.ap_obj_combo = QComboBox()
        self.ap_obj_combo.setMinimumHeight(26)
        self.ap_refresh_btn = QPushButton(t("refresh"))
        self.ap_refresh_btn.setObjectName("refresh_btn")
        self.ap_refresh_btn.setMinimumHeight(26)
        self.ap_refresh_btn.clicked.connect(self.refresh_objects)
        ap_obj_row_container = QWidget()
        ap_obj_row = QHBoxLayout(ap_obj_row_container)
        ap_obj_row.setContentsMargins(0, 0, 0, 0)
        ap_obj_row.addWidget(self.ap_obj_combo, 1)
        ap_obj_row.addWidget(self.ap_refresh_btn)
        ap_grid.addWidget(ap_obj_row_container, 0, 1)
        
        ap_grid.addWidget(QLabel("Atom1 Selection:"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.ap_atom1 = QLineEdit()
        self.ap_atom1.setMinimumHeight(26)
        self.ap_atom1.setPlaceholderText('e.g.: "resn LIG and name N1"')
        ap_grid.addWidget(self.ap_atom1, 0, 3)
        
        # 第二行: Atom2 Selection | [input] | Distance cutoff | [input]
        ap_grid.addWidget(QLabel("Atom2 Selection:"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.ap_atom2 = QLineEdit()
        self.ap_atom2.setMinimumHeight(26)
        self.ap_atom2.setPlaceholderText('e.g.: "elem O"')
        ap_grid.addWidget(self.ap_atom2, 1, 1)
        
        ap_grid.addWidget(QLabel("Distance cutoff (Å):"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.ap_distance = QLineEdit("5.0")
        self.ap_distance.setMinimumHeight(26)
        ap_grid.addWidget(self.ap_distance, 1, 3)
        
        # 第三行: Output CSV | [input + browse] (跨两列)
        ap_grid.addWidget(QLabel("Output CSV (optional)"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.ap_csv = QLineEdit()
        self.ap_csv.setMinimumHeight(26)
        self.ap_csv.setPlaceholderText("Optional")
        self.ap_csv_btn = QPushButton(t("browse"))
        self.ap_csv_btn.setMinimumHeight(26)
        self.ap_csv_btn.setObjectName("browse_btn")
        self.ap_csv_btn.clicked.connect(lambda: self._browse_save_file(self.ap_csv, "CSV (*.csv)"))
        ap_csv_row = QHBoxLayout()
        ap_csv_row.addWidget(self.ap_csv, 1)
        ap_csv_row.addWidget(self.ap_csv_btn)
        ap_grid.addLayout(ap_csv_row, 2, 1, 1, 3)  # 跨三列
        
        # Quick templates - 确保所有按钮在一行显示
        template_row = QHBoxLayout()
        templates = [
            ("N-O H-bonds", '"elem N"', '"elem O"', "3.5"),
            ("Lig-SER", '"resn LIG"', '"resn SER"', "4.5"),
            ("S-S", '"name SG"', '"name SG"', "2.5"),
            ("π-Stacking", '"aromatic"', '"aromatic"', "4.5"),
            ("π-Cation", '"aromatic"', '"basic"', "4.0"),
        ]
        for label, atom1, atom2, dist in templates:
            btn = QPushButton(label)
            btn.setObjectName("browse_btn") # 使用一个比较紧凑的样式
            btn.setMinimumHeight(28)
            btn.clicked.connect(lambda checked, a1=atom1, a2=atom2, d=dist: self.apply_ap_template(a1, a2, d))
            template_row.addWidget(btn)
        template_row.addStretch(1)
        
        btn_ap_row = QHBoxLayout()
        self.ap_analyze_btn = QPushButton("Analyze Pairs")
        self.ap_analyze_btn.setObjectName("highlight_btn")
        self.ap_analyze_btn.clicked.connect(self.run_ap_analysis)
        
        self.ap_visualize_btn = QPushButton("Visualize")
        self.ap_visualize_btn.setObjectName("highlight_btn")
        self.ap_visualize_btn.clicked.connect(self.run_ap_visualize)
        
        btn_ap_row.addWidget(self.ap_analyze_btn)
        btn_ap_row.addWidget(self.ap_visualize_btn)
        btn_ap_row.addStretch(1)
        
        main_layout.addWidget(grp_ap)
        main_layout.addWidget(QLabel("Quick Templates:"))
        main_layout.addLayout(template_row)
        main_layout.addLayout(btn_ap_row)
        
        main_layout.addStretch(1)
        return w
    
    def create_analysis_tab(self) -> QWidget:
        w = QWidget(); lay = QVBoxLayout(w); lay.setSpacing(8)

        # ===== 相互作用分析组 =====
        grp = QGroupBox(t("grp_analysis")); form = QFormLayout(grp); form.setLabelAlignment(Qt.AlignmentFlag.AlignRight); form.setSpacing(10)
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
        self.chk_between = QCheckBox(t("between_chains")); form.addRow(QLabel(""), self.chk_between)
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
        grp_csv = QGroupBox(t("grp_csv")); form_csv = QFormLayout(grp_csv); form_csv.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
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
        layout.setSpacing(20)
        layout.setContentsMargins(50, 30, 50, 30)
        
        # Welcome title
        title_label = QLabel("MolStruct Plugin for PyMOL")
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
        w = QWidget()
        main_layout = QVBoxLayout(w)
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
        gm_grid.setHorizontalSpacing(8)  # 恢复水平间距
        gm_grid.setVerticalSpacing(16)  # 再次增加此页面的垂直间距以解决重叠
        
        # 第一行： Target Object | [combo + refresh] | PDB File | [input + browse]
        gm_grid.addWidget(QLabel("Target Object"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.obj_combo_gm = QComboBox()
        self.obj_combo_gm.setMinimumHeight(24)
        self.refresh_obj_gm = QPushButton(t("refresh"))
        self.refresh_obj_gm.setObjectName("refresh_btn")
        self.refresh_obj_gm.setMinimumHeight(24)
        self.refresh_obj_gm.clicked.connect(self.refresh_objects)
        obj_row_container = QWidget()
        obj_row = QHBoxLayout(obj_row_container)
        obj_row.setContentsMargins(0, 0, 0, 0)
        obj_row.addWidget(self.obj_combo_gm, 1)
        obj_row.addWidget(self.refresh_obj_gm)
        gm_grid.addWidget(obj_row_container, 0, 1)
        
        gm_grid.addWidget(QLabel("PDB File (optional)"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.gm_pdb = QLineEdit()
        self.gm_pdb.setMinimumHeight(24)
        self.gm_pdb_browse = QPushButton(t("browse"))
        self.gm_pdb_browse.setMinimumHeight(24)
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
        self.gm_rmsd.setMinimumHeight(24)
        gm_grid.addWidget(self.gm_rmsd, 1, 1)
        
        gm_grid.addWidget(QLabel(""), 1, 2)  # 空位
        self.gm_require_gly = QCheckBox(t("require_gly"))
        self.gm_require_gly.setChecked(True)
        gm_grid.addWidget(self.gm_require_gly, 1, 3)
        
        # 第三行： Template Source | [combo] | Template Selection | [input + button]
        gm_grid.addWidget(QLabel("Template Source"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.gm_template_mode = QComboBox()
        self.gm_template_mode.setMinimumHeight(24)
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
        self.gm_template_sel.setMinimumHeight(24)
        self.gm_template_sel.setPlaceholderText("Enter selection, e.g., sele")
        self.gm_template_pick = QPushButton("Get Current (sele)")
        self.gm_template_pick.setMinimumHeight(24)
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
        self.gm_out_csv.setMinimumHeight(24)
        self.gm_out_csv.setPlaceholderText("Optional - leave blank for temp CSV")
        self.gm_out_browse = QPushButton(t("browse"))
        self.gm_out_browse.setMinimumHeight(24)
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
        
        # ========== Ternary Complex Analysis (with inte        grp_ternary = QGroupBox("Ternary Complex Analysis (E3-PROTAC-POI with Interface Tools)")
        # 使用网格布局实现两栏式
        ternary_grid = QGridLayout(grp_ternary)
        ternary_grid.setColumnStretch(0, 0)
        ternary_grid.setColumnStretch(1, 1)
        ternary_grid.setColumnStretch(2, 0)
        ternary_grid.setColumnStretch(3, 1)
        ternary_grid.setHorizontalSpacing(8)  # 恢复水平间距
        ternary_grid.setVerticalSpacing(16)  # 再次增加此页面的垂直间距以解决重叠
        
        # 第一行: Target Object | [combo + refresh] | Ligand Resname | [input]
        ternary_grid.addWidget(QLabel("Target Object"), 0, 0, Qt.AlignmentFlag.AlignRight)
        self.tc_obj_combo = QComboBox()
        self.tc_obj_combo.setMinimumHeight(24)
        self.tc_refresh_btn = QPushButton(t("refresh"))
        self.tc_refresh_btn.setObjectName("refresh_btn")
        self.tc_refresh_btn.setMinimumHeight(24)
        self.tc_refresh_btn.clicked.connect(self.refresh_objects)
        tc_obj_row_container = QWidget()
        tc_obj_row = QHBoxLayout(tc_obj_row_container)
        tc_obj_row.setContentsMargins(0, 0, 0, 0)
        tc_obj_row.addWidget(self.tc_obj_combo, 1)
        tc_obj_row.addWidget(self.tc_refresh_btn)
        ternary_grid.addWidget(tc_obj_row_container, 0, 1)
        
        ternary_grid.addWidget(QLabel("Ligand Resname:"), 0, 2, Qt.AlignmentFlag.AlignRight)
        self.tc_ligand_name = QLineEdit()
        self.tc_ligand_name.setMinimumHeight(24)
        self.tc_ligand_name.setPlaceholderText("PROTAC/Glue name, auto-detect if blank")
        ternary_grid.addWidget(self.tc_ligand_name, 0, 3)
        
        # 第二行: E3 Chains (CRBN/VHL) | [input] | POI Chains | [input]
        ternary_grid.addWidget(QLabel("E3 Ligase Chains:"), 1, 0, Qt.AlignmentFlag.AlignRight)
        self.tc_protein1_chains = QLineEdit()
        self.tc_protein1_chains.setMinimumHeight(24)
        self.tc_protein1_chains.setPlaceholderText("e.g.: A (CRBN/VHL/IAP)")
        ternary_grid.addWidget(self.tc_protein1_chains, 1, 1)
        
        ternary_grid.addWidget(QLabel("POI Chains:"), 1, 2, Qt.AlignmentFlag.AlignRight)
        self.tc_protein2_chains = QLineEdit()
        self.tc_protein2_chains.setMinimumHeight(24)
        self.tc_protein2_chains.setPlaceholderText("e.g.: B (Target Protein)")
        ternary_grid.addWidget(self.tc_protein2_chains, 1, 3)
        
        # 第三行: Interface cutoff | [input] | Output CSV | [input + browse]
        ternary_grid.addWidget(QLabel("Interface cutoff (Å):"), 2, 0, Qt.AlignmentFlag.AlignRight)
        self.tc_distance = QLineEdit("4.5")
        self.tc_distance.setMinimumHeight(24)
        ternary_grid.addWidget(self.tc_distance, 2, 1)
        
        ternary_grid.addWidget(QLabel("Output CSV (optional)"), 2, 2, Qt.AlignmentFlag.AlignRight)
        self.tc_csv = QLineEdit()
        self.tc_csv.setMinimumHeight(24)
        self.tc_csv.setPlaceholderText("Optional")
        self.tc_csv_btn = QPushButton(t("browse"))
        self.tc_csv_btn.setMinimumHeight(24)
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
        
        main_layout.addStretch(1)
        return w
    
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
            names = cmd.get_object_list() or []
        except Exception:
            pass
        if not names: names = [t("no_object")]
        for cb in (getattr(self, "obj_combo_gm", None),
                   getattr(self, "obj_combo_analysis", None),
                   getattr(self, "obj_combo_csv", None),
                   getattr(self, "obj_combo_apbs", None),
                   getattr(self, "pl_obj_combo", None),
                   getattr(self, "tc_obj_combo", None),
                   getattr(self, "ap_obj_combo", None)):
            if cb is not None:
                cb.blockSignals(True); cb.clear()
                for n in names: cb.addItem(n)
                cb.blockSignals(False)
        self.update_enablement()
        self.log(f"{t('refresh')} OK: {', '.join(names)}")

    def browse_csv(self):
        fn, _ = QFileDialog.getOpenFileName(self, t("select_csv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.csv_path.setText(fn); self.update_enablement()

    def browse_pdb(self):
        fn, _ = QFileDialog.getOpenFileName(self, t("select_pdb"), "", "PDB (*.pdb *.cif);;All Files (*)")
        if fn: self.pdb_path.setText(fn); self.update_enablement()

    def browse_out_csv(self):
        fn, _ = QFileDialog.getSaveFileName(self, t("select_outcsv"), "", "CSV (*.csv);;All Files (*)")
        if fn: self.out_csv.setText(fn); self.update_enablement()

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

            self.log(f"开始高亮显示：{os.path.basename(csv_file)}")
            result = highlight_csv_residues(csv_file, obj, show_labels=1, clear_old=1, debug=0, stick_by_element=1)

            if result and isinstance(result, dict):
                pairs = result.get("pairs", 0)
                unique = result.get("unique_residues", 0)
                self.log(f"高亮完成：{pairs} 对相互作用，{unique} 个唯一残基")
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
        pdb = self.pdb_path.text().strip() or None
        outcsv = self.out_csv.text().strip() or None
        only_between = self.chk_between.isChecked()
        if not obj or obj == t("no_object"):
            QMessageBox.warning(self, t("title"), t("no_object")); return
        self.analyze_btn.setEnabled(False)
        self.progress_bar.setVisible(True); self.progress_bar.setRange(0, 0)
        self.analysis_thread = AnalysisWorker(obj, pdb, only_between, outcsv)
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
                    self.log("已应用 G-Motif 高亮（基于最新 CSV）")
                except Exception as e:
                    self.log(f"G-Motif 高亮跳过：{e}")

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

                        self.log("静电势表面仅显示在 G-loop 周围 10Å 区域")
                    else:
                        # 没有 G-loop 数据，显示整个蛋白表面
                        self._show_full_surface_esp(obj, ramp_name)
                except Exception as e:
                    self.log(f"G-loop 区域静电势失败，回退到全表面: {e}")
                    self._show_full_surface_esp(obj, ramp_name)
            else:
                # 没有 G-Motif CSV，显示整个蛋白表面
                self._show_full_surface_esp(obj, ramp_name)


            # 4) 设置光照参数
            cmd.set("ambient", 0.2)
            cmd.set("spec_power", 80)
            cmd.set("spec_reflect", 0.3)
            cmd.set("depth_cue", 1)
            cmd.set("fog_start", 0.45)
            cmd.orient(obj)

            self.log(f"已渲染：G-Motif + ESP（grid={grid} Å, range=({vmin},{v0},{vmax})）")

            # 5) 导出 PNG
            self._export_png_for_object(obj)

        except Exception as e:
            self.on_error(str(e))

    def _show_full_surface_esp(self, obj: str, ramp_name: str):
        """显示整个蛋白的静电势表面（回退方案）"""
        from pymol import cmd
        cmd.show("surface", obj)
        cmd.set("surface_quality", 1, obj)
        cmd.set("surface_color_smoothing", 1, obj)
        cmd.set("transparency", 0.2, obj)
        cmd.color(ramp_name, obj)
        self.log("显示整个蛋白表面静电势（无 G-loop 数据）")


    # --- 线程回调 ---
    def on_finished_analysis(self, interactions: List[Dict[str, Any]]):
        self._interactions = interactions or []
        if self._interactions:
            self.fill_table_from_interactions(self._interactions)
        elif self.out_csv.text().strip() and os.path.exists(self.out_csv.text().strip()):
            self.fill_table_from_csv(self.out_csv.text().strip())
        self.progress_bar.setVisible(False); self.progress_bar.setRange(0, 1)
        self.analyze_btn.setEnabled(True)

    def on_finished_gmotif(self, hits: List[Tuple], out_csv_path: str):
        self._gmotif_hits = hits or []
        self._last_gmotif_csv = out_csv_path
        if os.path.exists(out_csv_path):
            self.fill_table_from_csv(out_csv_path)
            self.csv_path.setText(out_csv_path)
            self.log(f"G-Motif → CSV → Table: {os.path.basename(out_csv_path)}")
        else:
            self.fill_table_from_gmotif_hits(self._gmotif_hits)
        self.progress_bar.setVisible(False); self.progress_bar.setRange(0, 1)
        self.gm_btn.setEnabled(True)

    def on_error(self, msg: str):
        self.log(t("log_error").format(msg=msg))
        QMessageBox.critical(self, t("title"), msg)
        self.progress_bar.setVisible(False); self.progress_bar.setRange(0, 1)
        for b in (getattr(self, "analyze_btn", None), getattr(self, "gm_btn", None), getattr(self, "gm_btn_render", None)):
            if b: b.setEnabled(True)

    # --- 表格填充 ---
    def fill_table_from_interactions(self, rows: List[Dict[str, Any]]):
        headers = T["table_header"][get_lang()]
        self.table.clearContents(); self.table.setRowCount(len(rows))
        self.table.setHorizontalHeaderLabels(headers)
        keys = ["Chain1", "Residue1", "Chain2", "Residue2", "Distance", "Interaction"]
        for r, item in enumerate(rows):
            values = [item.get(k, "") for k in keys]
            for c, val in enumerate(values):
                it = QTableWidgetItem(str(val) if val is not None else "")
                if c == 4:
                    it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.table.setItem(r, c, it)
        self.table.resizeColumnsToContents()

    def fill_table_from_csv(self, csv_file: str):
        """智能识别CSV格式并填充表格"""
        data = []
        csv_type = "interaction"  # 默认为相互作用格式

        with open(csv_file, "r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)

            # 检测 CSV 类型 - 检查表头字段名
            fieldnames = reader.fieldnames or []
            
            if ("Chain" in fieldnames and "Sequence" in fieldnames and "Start" in fieldnames) or \
               ("链" in fieldnames and "序列" in fieldnames and "起始" in fieldnames):
                csv_type = "gmotif"
            elif "Protein_Atom" in fieldnames and "Ligand_Atom" in fieldnames:
                # 高级分析格式 (RDKit)
                csv_type = "advanced"

            # 读取所有数据行
            for row in reader:
                if csv_type == "gmotif":
                    data.append([
                        row.get("Chain", row.get("链", "")),
                        row.get("Sequence", row.get("序列", "")),
                        row.get("Start", row.get("起始", "")),
                        row.get("End", row.get("结束", "")),
                        row.get("RMSD", row.get("RMSD (Å)", "")),
                        row.get("Type", row.get("类型", "")),
                    ])
                elif csv_type == "advanced":
                    # 高级分析格式：Type, Protein_Atom, Ligand_Atom, Distance, Protein_Residue, Extra_Info
                    protein_residue = row.get("Protein_Residue", "")
                    if protein_residue:
                        protein_display = protein_residue
                    else:
                        protein_display = f"Atom {row.get('Protein_Atom', '')}"
                    
                    data.append([
                        row.get("Type", ""),
                        protein_display,
                        f"Atom {row.get('Ligand_Atom', '')}",
                        row.get("Distance", ""),
                        row.get("Extra_Info", ""),
                    ])
                else:
                    # 标准格式：Chain1, Residue1, Chain2, Residue2, Distance, Interaction
                    data.append([
                        row.get("Chain1", row.get("链1", "")),
                        row.get("Residue1", row.get("残基1", "")),
                        row.get("Chain2", row.get("链2", "")),
                        row.get("Residue2", row.get("残基2", "")),
                        row.get("Distance", row.get("距离", "")),
                        row.get("Interaction", row.get("相互作用", "")),
                    ])

        # 根据 CSV 类型设置表头
        if csv_type == "gmotif":
            headers = T["table_header_gmotif"][get_lang()]
        elif csv_type == "advanced":
            headers = ["类型", "蛋白残基", "配体原子", "距离(Å)", "详细信息"] if get_lang() == "zh" else ["Type", "Protein Residue", "Ligand Atom", "Distance(Å)", "Details"]
        else:
            headers = T["table_header"][get_lang()]

        # 设置表格列数和表头
        self.table.setColumnCount(len(headers))
        self.table.clearContents()
        self.table.setRowCount(len(data))
        self.table.setHorizontalHeaderLabels(headers)
        
        # 填充数据
        for r, values in enumerate(data):
            for c, val in enumerate(values):
                it = QTableWidgetItem(str(val))
                if c == 3 and csv_type == "advanced":  # Distance 列右对齐（高级格式）
                    it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                elif c == 4 and csv_type != "advanced":  # Distance 列右对齐（标准格式）
                    it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.table.setItem(r, c, it)
        self.table.resizeColumnsToContents()

    def fill_table_from_gmotif_hits(self, hits: List[Tuple]):
        # 使用 G-Motif 专用表头
        headers = T["table_header_gmotif"][get_lang()]
        self.table.clearContents(); self.table.setRowCount(len(hits))
        self.table.setHorizontalHeaderLabels(headers)
        # hits 格式: (chain, start, end, seq8, rmsd)
        for r, (ch, s, e, seq8, rmsd) in enumerate(hits):
            values = [ch, seq8, str(s), str(e), f"{rmsd:.2f}", "G-Motif"]
            for c, val in enumerate(values):
                it = QTableWidgetItem(str(val))
                if c == 4:  # RMSD 列右对齐
                    it.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.table.setItem(r, c, it)
        self.table.resizeColumnsToContents()

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
                self.log("APBS 工具不可用（未找到 apbs_tools.run_apbs）。已自动降级到 Quick 模式。")
                self.apbs_run_quick()
                return

            try:
                apbs_tools.run_apbs(selection=obj)
                self.log("APBS 计算已提交；若无可视化结果，请在 APBS Tools 中检查外部路径配置。")
            except Exception as ee:
                self.log(f"APBS 调用失败：{ee}。已自动降级到 Quick 模式。")
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
                self.log(f"视口尺寸: {w}x{h}px → 已填入导出设置")
            else:
                self.log("未能获取视口尺寸，已保持默认值")
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
        self.log(f"PNG 导出完成：{os.path.basename(fn)} | {W}x{H}px @ {dpi} dpi | ray={ray} | bg={'transparent' if want_trans else 'white'}")

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

    def export_dx(self):
        try:
            from pymol import cmd
            obj = self.obj_combo_apbs.currentText().strip()
            if not obj or obj == t("no_object"):
                QMessageBox.warning(self, t("title"), t("no_object")); return
            map_name = None
            if obj in self._esp_maps:
                map_name = self._esp_maps[obj][0]
            if not map_name:
                guess = f"{obj}_esp_map"
                map_name = guess
            fn, _ = QFileDialog.getSaveFileName(self, t("btn_export_dx"), f"{obj}_esp.dx", "DX (*.dx)")
            if not fn: return
            if not fn.lower().endswith(".dx"): fn += ".dx"
            cmd.save(fn, map_name)
            self.log(f"DX 导出完成：{os.path.basename(fn)} （源：{map_name}）")
        except Exception as e:
            self.on_error(str(e))

    # ==============  蛋白-配体分析标签页 ==============
    def create_prot_lig_tab(self) -> QWidget:
        """创建蛋白-配体分析标签页"""
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
        check_btn.clicked.connect(self.run_crbn_doctor)
        top_row.addWidget(check_btn)
        top_row.addStretch(1)
        layout.addLayout(top_row)
        
        # 标题
        title = QLabel("MolStruct - README" if get_lang() == "en" else "MolStruct - 说明文档")
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
<h2>欢迎使用 MolStruct</h2>

<p><b>MolStruct</b> 是一个功能强大的 PyMOL 插件，用于分子结构分析。</p>

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
<h2>Welcome to MolStruct</h2>

<p><b>MolStruct</b> is a powerful PyMOL plugin for molecular structure analysis.</p>

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
<h3>👋 感谢使用 MolStruct！</h3>

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
<p>如果 MolStruct 对你的研究有帮助，请考虑：</p><ul><li>在 GitHub 上给我们一个 Star ⭐</li><li>在论文中引用 MolStruct</li><li>分享给同事</li></ul>

<p style='{privacy_style}'>
<b>隐私声明：</b>MolStruct 不会收集任何个人数据或结构信息。所有分析都在本地进行。
</p>

<p style='margin-top: 20px; color: {footer_color}; text-align: center;'>
感谢你的支持！🚀
</p>'''
        else:
            contact_content = f'''
<h3>👋 Thank you for using MolStruct!</h3>

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
<p>If MolStruct helped your research, please consider:</p>
<ul><li>Giving us a Star on GitHub</li><li>Citing MolStruct in your papers</li><li>Sharing with colleagues</li></ul>

<p style='{privacy_style}'>
<b>Privacy:</b> MolStruct does not collect any personal data or structural information. All analyses are performed locally.
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

    def apply_ap_template(self, atom1, atom2, dist):
        """应用原子对模板"""
        self.ap_atom1.setText(atom1)
        self.ap_atom2.setText(atom2)
        self.ap_distance.setText(dist)
        msg = f"已应用模板: {atom1} - {atom2}" if get_lang() == "zh" else f"Template applied: {atom1} - {atom2}"
        self.log(msg)

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

            # ========== 统一分析（自动检测使用RDKit或基础模式） ==========
            from .interaction_analyzer import analyze_protein_ligand_interactions

            protein_chains_str = self.pl_protein_chains.text().strip()
            protein_chains = [c.strip() for c in protein_chains_str.split(",")] if protein_chains_str else None
            
            try:
                distance_cutoff = float(self.pl_distance.text())
            except:
                distance_cutoff = 4.5
                self.log(f"距离参数无效，使用默认值 4.5 Å")

            self.log(f"\n{'开始蛋白-配体分析...' if get_lang() == 'zh' else 'Starting Protein-Ligand analysis...'}")
            self.log(f"   对象: {obj_name}")
            self.log(f"   配体: {ligand_resname or ('自动检测' if get_lang() == 'zh' else 'auto-detect')}")
            self.log(f"   蛋白链: {protein_chains or '自动检测'}")
            self.log(f"   距离截断: {distance_cutoff} Å")

            result = analyze_protein_ligand_interactions(
                obj_name=obj_name,
                ligand_resname=ligand_resname,
                protein_chains=protein_chains,
                distance_cutoff=distance_cutoff,
                output_csv=output_csv
            )

            if result:
                mode = result.get("mode", "unknown")
                mode_text = "RDKit 高级模式" if mode == "advanced" else "基础模式"
                
                # 处理结果
                if mode == "advanced":
                    # 高级分析模式
                    advanced_results = result.get("advanced_results", {})
                    total = sum(len(v) for v in advanced_results.values())
                    self.log(f"\n分析完成 ({mode_text}): 总计 {total} 个相互作用")
                    for itype, interactions in advanced_results.items():
                        if interactions:
                            self.log(f"   - {itype}: {len(interactions)}")
                    
                    # 保存结果供可视化使用
                    self.current_pl_result_advanced = advanced_results
                    self.current_pl_result = result
                    self.current_pl_ligand_sdf = result.get("ligand_sdf", None)
                    
                    result_msg = f"分析完成 ({mode_text})\n总计: {total} 个相互作用" if get_lang() == "zh" else f"Done ({mode_text})\nTotal: {total} interactions"
                else:
                    # 基础分析模式
                    n_interactions = len(result["interactions"])
                    self.log(f"\n分析完成 ({mode_text}): 发现 {n_interactions} 个相互作用")
                    self.current_pl_result = result
                    if hasattr(self, 'current_pl_result_advanced'):
                        delattr(self, 'current_pl_result_advanced')
                    
                    result_msg = f"发现 {n_interactions} 个相互作用 ({mode_text})" if get_lang() == "zh" else f"Found {n_interactions} interactions ({mode_text})"
                
                # 显示CSV表格
                if output_csv and os.path.exists(output_csv):
                    try:
                        self.fill_table_from_csv(output_csv)
                        self.log(f"结果已显示在表格中")
                    except Exception as e:
                        self.log(f"无法显示表格: {e}")
                
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

            self.log(f"\n{'生成3D可视化...' if get_lang() == 'zh' else 'Generating 3D visualization...'}")
            
            # 检查是否使用了高级分析模式
            if hasattr(self, 'current_pl_result_advanced'):
                # 高级分析模式 - 使用专用的可视化函数
                from .interaction_analyzer_advanced import visualize_advanced_interactions
                visualize_advanced_interactions(obj_name=obj_name, csv_path=csv_path, ligand_resname=ligand_resname)
                self.log("3D可视化完成（高级分析模式）" if get_lang() == "zh" else "3D visualization done (advanced mode)")
            else:
                # 标准分析模式
                from .interaction_analyzer import visualize_protein_ligand_3d
                if csv_path and os.path.exists(csv_path):
                    visualize_protein_ligand_3d(
                        obj_name=obj_name,
                        csv_path=csv_path,
                        ligand_resname=ligand_resname
                    )
                elif hasattr(self, 'current_pl_result'):
                    visualize_protein_ligand_3d(
                        obj_name=obj_name,
                        interactions_result=self.current_pl_result,
                        ligand_resname=ligand_resname
                    )
                else:
                    self.log("没有可用的相互作用数据" if get_lang() == "zh" else "No interaction data available")
                    return
                self.log("3D可视化完成" if get_lang() == "zh" else "3D visualization done")

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

            self.log(f"\n{'生成交互网络图...' if get_lang() == 'zh' else 'Generating network plot...'}")
            
            # 优先使用CSV文件（支持两种模式）
            csv_path = self.pl_csv.text().strip() or None
            ligand_sdf = getattr(self, 'current_pl_ligand_sdf', None)
            
            if csv_path and os.path.exists(csv_path):
                output_path = generate_interaction_network_plot(
                    csv_path=csv_path,
                    output_path=fn,
                    show_plot=False,
                    ligand_sdf=ligand_sdf
                )
            elif hasattr(self, 'current_pl_result'):
                output_path = generate_interaction_network_plot(
                    interactions_result=self.current_pl_result,
                    output_path=fn,
                    show_plot=False,
                    ligand_sdf=ligand_sdf
                )
            else:
                self.log("没有可用的数据生成网络图" if get_lang() == "zh" else "No data for network plot")
                return

            if output_path:
                self.log(f"网络图已保存: {output_path}")
            else:
                self.log(f"网络图生成失败" if get_lang() == "zh" else "Network plot failed")

        except Exception as e:
            self.log(f"错误: {e}")
            import traceback; traceback.print_exc()

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

            self.log(f"\n{'开始三元复合体分析...' if get_lang() == 'zh' else 'Starting ternary complex analysis...'}")

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
                self.log(f"分析完成:" if get_lang() == "zh" else "Done:")
                self.log(f"   蛋白1: {stats['protein1_interactions_count']} 个相互作用")
                self.log(f"   蛋白2: {stats['protein2_interactions_count']} 个相互作用")
                self.current_tc_result = result
                QMessageBox.information(self, "完成" if get_lang() == "zh" else "Done",
                    f"三元复合体分析完成\n\n蛋白1: {stats['protein1_interactions_count']} 个相互作用\n蛋白2: {stats['protein2_interactions_count']} 个相互作用" if get_lang() == "zh" else
                    f"Ternary analysis done\n\nProtein1: {stats['protein1_interactions_count']}\nProtein2: {stats['protein2_interactions_count']}")
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

            self.log(f"\n{'生成三元复合体网络图...' if get_lang() == 'zh' else 'Generating ternary network...'}")
            output_path = generate_interaction_network_plot(
                interactions_result=self.current_tc_result,
                output_path=fn,
                show_plot=False
            )

            if output_path:
                self.log(f"网络图已保存: {output_path}")

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
            
            self.log(f"\n{'分析蛋白-蛋白界面...' if get_lang() == 'zh' else 'Analyzing protein-protein interface...'}")
            
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
            
            self.log(f"界面分析完成:")
            self.log(f"  E3侧: {n_e3_residues} 个残基")
            self.log(f"  POI侧: {n_poi_residues} 个残基")
            
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
                
            self.log(f"\n{'渲染三元复合体...' if get_lang() == 'zh' else 'Rendering ternary complex...'}")
            
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
            
            self.log(f"渲染完成，已保存至: {output_path}")
            
        except Exception as e:
            self.log(f"错误: {e}")
            import traceback; traceback.print_exc()

    def run_ap_analysis(self):
        """运行原子对分析"""
        try:
            from .interaction_analyzer import analyze_atom_pair_interactions

            obj_name = self.ap_obj_combo.currentText()
            if obj_name == t("no_object"):
                QMessageBox.warning(self, t("title"), "请先加载PDB结构" if get_lang() == "zh" else "Load PDB first")
                return

            atom1_sel = self.ap_atom1.text().strip()
            atom2_sel = self.ap_atom2.text().strip()

            if not atom1_sel or not atom2_sel:
                QMessageBox.warning(self, "警告" if get_lang() == "zh" else "Warning",
                    "请输入原子选择条件" if get_lang() == "zh" else "Enter atom selections")
                return

            distance_cutoff = float(self.ap_distance.text())
            output_csv = self.ap_csv.text().strip() or None

            self.log(f"\n{'开始原子对分析...' if get_lang() == 'zh' else 'Starting atom pair analysis...'}")
            self.log(f"   原子1: {atom1_sel}")
            self.log(f"   原子2: {atom2_sel}")

            result = analyze_atom_pair_interactions(
                obj_name=obj_name,
                atom1_selection=atom1_sel,
                atom2_selection=atom2_sel,
                distance_cutoff=distance_cutoff,
                output_csv=output_csv
            )

            if result:
                msg = f"分析完成: 发现 {len(result)} 个原子对相互作用" if get_lang() == "zh" else f"Done: {len(result)} atom pairs"
                self.log(msg)
                self.current_ap_result = result
                QMessageBox.information(self, "完成" if get_lang() == "zh" else "Done",
                    f"发现 {len(result)} 个原子对相互作用" if get_lang() == "zh" else f"Found {len(result)} atom pairs")
            else:
                self.log("未找到符合条件的原子对" if get_lang() == "zh" else "No atom pairs found")

        except Exception as e:
            self.log(f"错误: {e}")
            import traceback; traceback.print_exc()

    def run_ap_visualize(self):
        """可视化原子对"""
        try:
            from .interaction_analyzer import visualize_atom_pairs

            if not hasattr(self, 'current_ap_result'):
                QMessageBox.warning(self, "警告" if get_lang() == "zh" else "Warning",
                    "请先运行原子对分析" if get_lang() == "zh" else "Run atom pair analysis first")
                return

            obj_name = self.ap_obj_combo.currentText()

            self.log(f"\n{'可视化原子对...' if get_lang() == 'zh' else 'Visualizing atom pairs...'}")
            visualize_atom_pairs(
                obj_name=obj_name,
                interactions_result=self.current_ap_result
            )
            self.log("可视化完成" if get_lang() == "zh" else "Visualization done")

        except Exception as e:
            self.log(f"错误: {e}")

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
                self.log("📊 开始相互作用分析...")
                pdb = self.pdb_path.text().strip() or None
                only_between = self.chk_between.isChecked()

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
                    self.log(f"分析完成，发现 {len(interactions)} 个相互作用")
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
                    self.log(f"分析完成，发现 {len(interactions)} 个相互作用（临时CSV）")
                    self.fill_table_from_csv(csv_path)

            # 步骤2: 先高亮相互作用残基
            if has_csv:
                self.log("高亮相互作用残基...")
                try:
                    highlight_csv_residues(csv_path, obj=obj, show_labels=1,
                                         stick_by_element=1, clear_old=1)
                    self.log("残基高亮完成")
                except Exception as e:
                    self.log(f"高亮失败: {e}")

            # 步骤3: 美化渲染（不再重复高亮）
            self.log("应用美化渲染...")
            render_interactions_beautifully(obj, csv_path=None)  # 不传csv_path，避免重复高亮
            self.log("渲染完成")

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
            self.log(f"\n{'='*50}")
            self.log("⚡ 检查环境和依赖" if get_lang() == "zh" else "⚡ Checking Environment and Dependencies")
            self.log(f"{'='*50}\n")
            
            # 检查 Python 依赖
            try:
                from .env_setup import get_dependency_status
                status = get_dependency_status()
                
                self.log("📦 Python 依赖:" if get_lang() == "zh" else "📦 Python Dependencies:")
                for pkg, available in status.items():
                    symbol = "✅" if available else "❌"
                    self.log(f"  {symbol} {pkg:20} {'已安装' if available else '未安装'}")
                
                missing = [pkg for pkg, avail in status.items() if not avail]
                if missing:
                    self.log(f"\n⚠️  缺失依赖: {', '.join(missing)}")
                    self.log("💡 安装命令: pip install " + " ".join(missing))
                else:
                    self.log("\n✅ 所有 Python 依赖已满足" if get_lang() == "zh" else "\n✅ All Python dependencies satisfied")
            except Exception as e:
                self.log(f"⚠️  无法检查 Python 依赖: {e}")
            
            # 检查 FoldX
            self.log("\n🛠️  外部工具:" if get_lang() == "zh" else "\n🛠️  External Tools:")
            try:
                self._ensure_crbn_tools_loaded()
                from pymol import cmd
                ok = cmd.crbn_tools_doctor(verbose=0)
                if ok:
                    self.log("  ✅ FoldX: 已检测到")
                else:
                    self.log("  ❌ FoldX: 未检测到")
                    self.log("     💡 下载: https://foldxsuite.crg.eu/")
                    self.log("     💡 设置: export FOLDX=/path/to/foldx")
            except Exception as e:
                self.log(f"  ⚠️  FoldX 检查失败: {e}")
            
            # PyMOL 版本
            try:
                from pymol import cmd
                version = cmd.get_version()[0]
                self.log(f"\n🐍 PyMOL: {version}")
            except Exception:
                pass
            
            self.log(f"\n{'='*50}")
            self.log("✅ 环境检查完成" if get_lang() == "zh" else "✅ Environment check complete")
            self.log(f"{'='*50}\n")
            
        except Exception as e:
            self.log(f"❌ 检查失败: {e}")
            import traceback
            traceback.print_exc()
    
    def run_crbn_doctor(self):
        """调用 PyMOL CRBN crbn_tools_doctor 命令检查环境"""
        try:
            self._ensure_crbn_tools_loaded()
            
            self.log(f"\n{'检查环境...' if get_lang() == 'zh' else 'Checking environment...'}")
            
            from pymol import cmd
            ok = cmd.crbn_tools_doctor()
            
            if ok:
                self.log("环境检查通过！FoldX 可用。" if get_lang() == "zh" else "Environment OK! FoldX available.")
            else:
                self.log("FoldX 未检测到，ΔΔG 计算将使用 ASA 近似。" if get_lang() == "zh" else "FoldX not found, ΔΔG will use ASA proxy.")
            
        except Exception as e:
            self.log(f"错误: {e}")
            import traceback
            traceback.print_exc()
    
    def _ensure_crbn_tools_loaded(self):
        """确保 pymol_crbn_tools.py 已加载"""
        from pymol import cmd
        
        # 检查是否已加载（通过检查命令是否存在）
        if hasattr(cmd, 'interface_map'):
            return  # 已加载
        
        # 尝试加载
        here = os.path.dirname(os.path.abspath(__file__))
        parent = os.path.dirname(here)
        
        # 在多个位置查找
        search_paths = [
            os.path.join(parent, "pymol_crbn_tools.py"),
            os.path.join(here, "pymol_crbn_tools.py"),
            os.path.join(os.getcwd(), "pymol_crbn_tools.py"),
        ]
        
        for path in search_paths:
            if os.path.exists(path):
                self.log(f"加载 CRBN 工具: {os.path.basename(path)}")
                try:
                    cmd.do(f"run {path}")
                    self.log("CRBN 工具加载成功")
                    return
                except Exception as e:
                    self.log(f"加载失败: {e}")
        
        # 如果找不到，提示用户
        raise FileNotFoundError(
            "pymol_crbn_tools.py not found. Please ensure it's in the plugin directory."
        )
    
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
        self.setup_style()
        # Update theme icon (sun/moon)
        self.theme_toggle_btn.setText("🌙" if self._dark_mode else "☀️")
        # Update Modules button color (easter egg effect)
        self.update_modules_button_style()
        # 主题变化后重新应用自动缩放，以确保字体与行距匹配
        self.apply_auto_scaling()
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
            # 兼容 PyQt5/6
            from PyQt5.QtGui import QGuiApplication as _QGA  # type: ignore
        except Exception:
            try:
                from PyQt6.QtGui import QGuiApplication as _QGA  # type: ignore
            except Exception:
                _QGA = None  # type: ignore
        dpi = 96.0
        if _QGA is not None:
            scr = _QGA.primaryScreen()
            try:
                if scr is not None:
                    dpi = float(scr.logicalDotsPerInch())
            except Exception:
                pass
        scale = dpi / 96.0
        # 合理范围
        if scale < 0.9:
            scale = 0.9
        if scale > 1.6:
            scale = 1.6
        return scale

    def apply_auto_scaling(self):
        # 禁用自动缩放，保持固定高度
        # s = self._ui_scale = float(self._compute_ui_scale())
        self._ui_scale = 1.0  # 固定缩放比例为1.0
        # 只更新表格行高，不改变控件高度
        try:
            vh = self.table.verticalHeader()
            vh.setDefaultSectionSize(34)  # 固定行高
        except Exception:
            pass
        # 只缩放主题切换按钮
        try:
            self.theme_toggle_btn.setFixedSize(28, 28)  # 固定尺寸
        except Exception:
            pass

    def setup_style(self):
        # 现代化样式 - 科技感与专业性并重
        try:
            from .modern_style import get_modern_stylesheet
            # 使用实例变量控制主题
            self.setStyleSheet(get_modern_stylesheet(dark_mode=self._dark_mode))
            return
        except Exception as e:
            # 回退到内置样式
            print(f"Modern style not available: {e}")
        
        # 回退样式 - Mac兼容性优先
        self.setStyleSheet("""
QDialog {
    background-color: #f5f7fa;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Microsoft YaHei", sans-serif;
    font-size: 11px;
    color: #222;
}

QListWidget {
    background-color: white;
    border: 1px solid #ddd;
    border-radius: 6px;
    padding: 6px;
    outline: none;
}
QListWidget::item {
    padding: 10px 12px;
    margin: 2px 0;
    border-radius: 4px;
    color: #333;
}
QListWidget::item:hover {
    background-color: #f0f0f0;
}
QListWidget::item:selected {
    background-color: #2563eb;
    color: white;
    font-weight: bold;
}

QGroupBox {
    font-weight: bold;
    font-size: 12px;
    border: 1px solid #ddd;
    border-radius: 6px;
    margin-top: 10px;
    padding: 18px 12px 10px 12px;
    background-color: white;
    color: #111;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 4px;
    background-color: white;
}

QLabel {
    color: #333;
    font-size: 12px;     /* 适中的字体大小 */
    padding: 2px 0;      /* 适中的垂直间距 */
    min-height: 20px;    /* 适中的最小高度 */
}

QLineEdit, QComboBox {
    border: 1px solid #ccc;
    border-radius: 4px;
    padding: 6px 8px;      /* 更紧凑的垂直内边距 */
    background-color: white;
    min-height: 26px;      /* 更小的最小高度 */
    color: #222;
    font-size: 12px;       /* 适中的字体大小 */
}
QLineEdit:focus, QComboBox:focus {
    border: 2px solid #2563eb;
}
QLineEdit:hover, QComboBox:hover {
    border-color: #999;
}
QComboBox::drop-down {
    border: none;
    padding-right: 4px;
}

QPushButton {
    border: 1px solid #ccc;
    border-radius: 4px;
    padding: 8px 14px;
    background-color: #f5f5f5;
    color: #222;
    font-weight: bold;
    min-height: 28px;
    font-size: 11px;
}
QPushButton:hover {
    background-color: #e8e8e8;
    border-color: #999;
}
QPushButton:pressed {
    background-color: #d0d0d0;
}
QPushButton:disabled {
    background-color: #f5f5f5;
    color: #999;
    border-color: #ddd;
}

QPushButton#highlight_btn {
    background-color: #f97316;
    color: white;
    border: none;
    font-weight: bold;
    font-size: 11px;
    padding: 9px 16px;
    min-height: 32px;
}
QPushButton#highlight_btn:hover {
    background-color: #ea580c;
}
QPushButton#highlight_btn:pressed {
    background-color: #c2410c;
}

QPushButton#refresh_btn {
    background-color: #f5f5f5;
    color: #666;
    border: 1px solid #ddd;
    padding: 6px 10px;
    font-size: 10px;
    min-width: 60px;
}
QPushButton#refresh_btn:hover {
    background-color: #e8e8e8;
}

QPushButton#theme_toggle_btn_small {
    background-color: #f5f5f5;
    color: #666;
    border: 1px solid #ddd;
    border-radius: 4px;
    padding: 2px;
    font-size: 16px;
    font-weight: normal;
}
QPushButton#theme_toggle_btn_small:hover {
    background-color: #e8e8e8;
    border-color: #999;
}
QPushButton#theme_toggle_btn_small:pressed {
    background-color: #d0d0d0;
}

QPushButton#bottom_nav_btn {
    background-color: transparent;
    color: #666;
    border: none;
    text-align: left;
    padding: 8px 12px;
    font-size: 12px;
    font-weight: 500;
}
QPushButton#bottom_nav_btn:hover {
    background-color: #f1f5f9;
    color: #3b82f6;
}

QPushButton#browse_btn {
    background-color: #2563eb;
    color: white;
    border: none;
    font-weight: bold;
    min-width: 70px;
    padding: 6px 12px;
    font-size: 10px;
}
QPushButton#browse_btn:hover {
    background-color: #1d4ed8;
}

QPushButton#save_btn {
    background-color: #059669;
    color: white;
    border: none;
    font-weight: bold;
    min-width: 70px;
    padding: 6px 12px;
    font-size: 10px;
}
QPushButton#save_btn:hover {
    background-color: #047857;
}

QPushButton#clear_btn {
    background-color: #666;
    color: white;
    border: none;
    font-weight: bold;
    padding: 8px 14px;
}
QPushButton#clear_btn:hover {
    background-color: #555;
}

QPushButton#close_btn {
    background-color: #dc2626;
    color: white;
    border: none;
    font-weight: bold;
    min-width: 90px;
    padding: 8px 14px;
}
QPushButton#close_btn:hover {
    background-color: #b91c1c;
}

QTextEdit {
    border: 1px solid #ddd;
    border-radius: 4px;
    padding: 8px;
    background-color: #fafafa;
    font-family: Menlo, Monaco, 'Courier New', monospace;
    font-size: 10px;
    color: #222;
}

QTableWidget {
    border: 1px solid #ddd;
    border-radius: 4px;
    background-color: white;
    gridline-color: #eee;
    selection-background-color: #2563eb;
    selection-color: white;
}
QHeaderView::section {
    background-color: #f0f0f0;
    color: #222;
    padding: 6px;
    border: none;
    border-bottom: 1px solid #ddd;
    font-weight: bold;
    font-size: 11px;
}
QTableWidget::item {
    padding: 8px 6px;   /* 适中的垂直内边距 */
    min-height: 32px;   /* 适中的最小行高 */
    border-bottom: 1px solid #eee;
}
QTableWidget::item:selected {
    background-color: #2563eb;
    color: white;
}
QTableWidget QTableCornerButton::section {
    background-color: #f0f0f0;
    border: none;
    border-bottom: 1px solid #ddd;
}

QCheckBox {
    color: #333;
    spacing: 6px;
    font-size: 11px;
    padding: 3px 0;
}
QCheckBox::indicator {
    width: 14px;
    height: 14px;
    border: 1px solid #999;
    border-radius: 3px;
    background-color: white;
}
QCheckBox::indicator:hover {
    border-color: #2563eb;
}
QCheckBox::indicator:checked {
    background-color: #2563eb;
    border-color: #2563eb;
}

QProgressBar {
    border: 1px solid #ddd;
    border-radius: 4px;
    text-align: center;
    background-color: white;
    height: 22px;
    color: #222;
    font-size: 10px;
}
QProgressBar::chunk {
    background-color: #2563eb;
    border-radius: 2px;
}

QScrollBar:vertical {
    border: none;
    background: #f5f5f5;
    width: 10px;
    border-radius: 5px;
}
QScrollBar::handle:vertical {
    background: #999;
    border-radius: 5px;
    min-height: 20px;
}
QScrollBar::handle:vertical:hover {
    background: #666;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0;
}
QScrollBar:horizontal {
    border: none;
    background: #f5f5f5;
    height: 10px;
}
QScrollBar::handle:horizontal {
    background: #999;
    border-radius: 5px;
    min-width: 20px;
}
QScrollBar::handle:horizontal:hover {
    background: #666;
}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
    width: 0;
}
""")

# -------- 独立运行入口（可选）--------
def launch_standalone():
    app = QApplication.instance() or QApplication(sys.argv)
    dlg = MolStructDialog()
    try:
        sys.exit(app.exec())
    except AttributeError:
        sys.exit(app.exec_())

if __name__ == "__main__":
    launch_standalone()
