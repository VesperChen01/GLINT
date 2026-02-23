# -*- coding: utf-8 -*-
"""
Utility functions and common constants for GLINT GUI.

This module provides:
- Internationalization (i18n) support with translation functions
- Logo and asset path utilities
- Dependency checking and installation
- Dynamic module loading helpers
- Import helpers for analysis modules
"""
import os
import sys
import importlib.util
from typing import List, Optional, Any, Dict, Callable, Tuple, Union

# ============================================================================
# Language & Internationalization
# ============================================================================

LANG_FORCE = "en"  # Force English for all GUI


def get_lang() -> str:
    """
    Get the current language setting.
    
    Returns:
        Language code (currently always 'en' for English)
    """
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
    "btn_load_csv_to_table": {"zh": "载入到表格", "en": "Load to Table"},
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
    """
    Translate a key to the current language.
    
    Args:
        key: Translation key from the T dictionary
        
    Returns:
        Translated string, or the key itself if not found
    """
    lang = get_lang()
    d = T.get(key)
    if isinstance(d, dict):
        return d.get(lang, list(d.values())[0])
    if d is None:
        return key  # Return key if not found
    return str(d)


# ============================================================================
# QMessageBox Helpers (Fix display issues)
# ============================================================================

def show_message_box(parent, title: str, message: str, icon_type: str = "information",
                     min_width: int = 400, font_size: int = 13) -> Any:
    """
    显示一个格式良好的 QMessageBox，确保内容完整显示
    弹窗大小会根据文本内容自动调整

    Args:
        parent: 父窗口
        title: 对话框标题
        message: 消息内容
        icon_type: 图标类型 ("information", "warning", "critical", "question")
        min_width: 最小宽度（像素）
        font_size: 字体大小

    Returns:
        对话框的返回值

    Example:
        show_message_box(self, "Success", "Analysis complete!\nResults saved to:\n/path/to/file.csv")
    """
    try:
        from .qt_adapter import QMessageBox
    except ImportError:
        # Fallback if qt_adapter not available
        try:
            from PyQt5.QtWidgets import QMessageBox
        except ImportError:
            from PyQt6.QtWidgets import QMessageBox

    msg_box = QMessageBox(parent)
    msg_box.setWindowTitle(title)

    # 设置图标
    icon_map = {
        "information": QMessageBox.Icon.Information,
        "warning": QMessageBox.Icon.Warning,
        "critical": QMessageBox.Icon.Critical,
        "question": QMessageBox.Icon.Question
    }
    msg_box.setIcon(icon_map.get(icon_type, QMessageBox.Icon.Information))

    # 设置消息文本
    msg_box.setText(message)

    # 根据文本内容动态计算宽度
    lines = message.split('\n')
    max_line_len = max(len(line) for line in lines) if lines else 0
    num_lines = len(lines)

    # 估算宽度：每个字符约 10 像素，加上边距
    estimated_width = max(min_width, min(max_line_len * 10 + 100, 800))
    # 估算高度：每行约 30 像素，加上边距 (Base height 150)
    estimated_height = max(150, min(num_lines * 30 + 150, 600))

    # 设置最小尺寸和样式，确保内容完整显示
    msg_box.setMinimumWidth(estimated_width)
    msg_box.setMinimumHeight(estimated_height)
    msg_box.setStyleSheet(f"""
        QMessageBox {{
            min-width: {estimated_width}px;
            min-height: {estimated_height}px;
        }}
        QMessageBox QLabel {{
            min-width: {estimated_width - 80}px;
            font-size: {font_size}px;
        }}
    """)

    return msg_box.exec()


def show_question_box(parent, title: str, message: str,
                      min_width: int = 400, font_size: int = 13) -> bool:
    """
    显示一个 Yes/No 问题对话框，大小根据内容自动调整

    Args:
        parent: 父窗口
        title: 对话框标题
        message: 消息内容
        min_width: 最小宽度（像素）
        font_size: 字体大小

    Returns:
        True if Yes, False if No
    """
    try:
        from .qt_adapter import QMessageBox
    except ImportError:
        try:
            from PyQt5.QtWidgets import QMessageBox
        except ImportError:
            from PyQt6.QtWidgets import QMessageBox

    msg_box = QMessageBox(parent)
    msg_box.setWindowTitle(title)
    msg_box.setIcon(QMessageBox.Icon.Question)
    msg_box.setText(message)
    msg_box.setStandardButtons(QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)

    # 根据文本内容动态计算宽度
    lines = message.split('\n')
    max_line_len = max(len(line) for line in lines) if lines else 0
    estimated_width = max(min_width, min(max_line_len * 9 + 100, 800))

    msg_box.setMinimumWidth(estimated_width)
    msg_box.setStyleSheet(f"""
        QMessageBox {{
            min-width: {estimated_width}px;
        }}
        QMessageBox QLabel {{
            min-width: {estimated_width - 80}px;
            font-size: {font_size}px;
        }}
    """)

    result = msg_box.exec()
    return result == QMessageBox.StandardButton.Yes


# ============================================================================
# Asset Paths
# ============================================================================

def _get_logo_path() -> Optional[str]:
    """
    Get the path to the logo file.
    
    Returns:
        Path to logo.png if it exists, None otherwise
    """
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    logo_path = os.path.join(here, "assets", "logo.png")
    if os.path.exists(logo_path):
        return logo_path
    return None


def _get_asset_path(filename: str) -> Optional[str]:
    """
    Get the path to an asset file.
    
    Args:
        filename: Name of the asset file
        
    Returns:
        Full path to the asset if it exists, None otherwise
    """
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    asset_path = os.path.join(here, "assets", filename)
    if os.path.exists(asset_path):
        return asset_path
    return None


# ============================================================================
# Dependency Checking
# ============================================================================

def _check_and_install_deps() -> bool:
    """
    Check and install dependencies at GUI startup.
    
    Returns:
        True if all dependencies are satisfied, False otherwise
    """
    try:
        # Try to import env_checker module
        try:
            from ..env_checker import ensure_dependencies
        except (ImportError, ValueError):
            sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            from env_checker import ensure_dependencies
        
        return ensure_dependencies(silent=False)
    except Exception as e:
        print(f"[GLINT] Dependency check failed: {e}")
        return False


# ============================================================================
# Dynamic Module Loading
# ============================================================================

def _dynamic_load_by_filenames(names: List[str], symbol: str) -> Optional[Any]:
    """
    Dynamically load a symbol from one of several possible module files.
    
    Args:
        names: List of possible filenames to try
        symbol: Name of the symbol to import from the module
        
    Returns:
        The imported symbol, or None if not found
    """
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
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
    """Return helper functions, handling imports robustly."""
    try:
        from ..highlight_residues import highlight_csv_residues, highlight_gmotif_loops
        from ..interaction_analyzer import (
            analyze_pdb_interactions, render_interactions_beautifully,
            analyze_protein_ligand_interactions, visualize_protein_ligand_3d,
            generate_interaction_network_plot, analyze_ternary_complex,
            analyze_atom_pair_interactions, visualize_atom_pairs,
            analyze_protein_nucleic_interactions
        )
        from ..interaction_2d_plot import generate_2d_interaction_diagram
        try:
            from ..ligand_ligand_analyzer import analyze_ligand_ligand_interactions
        except ImportError:
            analyze_ligand_ligand_interactions = None

        find_crbn_g_motif = None
        try:
            try:
                from ..g_motif_analyzer import find_crbn_g_motif
            except Exception:
                from ..g_motif import find_crbn_g_motif
        except Exception:
            find_crbn_g_motif = _dynamic_load_by_filenames(
                ["g_motif_analyzer.py", "g_motif.py", "g-motif.py"], "find_crbn_g_motif"
            )
        

        
        return highlight_csv_residues, highlight_gmotif_loops, analyze_pdb_interactions, find_crbn_g_motif, render_interactions_beautifully, generate_2d_interaction_diagram, analyze_protein_ligand_interactions, visualize_protein_ligand_3d, generate_interaction_network_plot, analyze_ternary_complex, analyze_atom_pair_interactions, visualize_atom_pairs, analyze_ligand_ligand_interactions, analyze_protein_nucleic_interactions
    except (ImportError, ValueError):
        pass
        
    # Fallback: add parent dir to sys.path
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if here and here not in sys.path:
        sys.path.insert(0, here)
    
    try:
        import highlight_residues
        highlight_csv_residues = highlight_residues.highlight_csv_residues
        highlight_gmotif_loops = highlight_residues.highlight_gmotif_loops
        
        import interaction_analyzer
        analyze_pdb_interactions = interaction_analyzer.analyze_pdb_interactions
        render_interactions_beautifully = interaction_analyzer.render_interactions_beautifully
        analyze_protein_ligand_interactions = interaction_analyzer.analyze_protein_ligand_interactions
        visualize_protein_ligand_3d = interaction_analyzer.visualize_protein_ligand_3d
        generate_interaction_network_plot = interaction_analyzer.generate_interaction_network_plot
        analyze_ternary_complex = interaction_analyzer.analyze_ternary_complex
        analyze_atom_pair_interactions = interaction_analyzer.analyze_atom_pair_interactions
        visualize_atom_pairs = interaction_analyzer.visualize_atom_pairs
        analyze_protein_nucleic_interactions = interaction_analyzer.analyze_protein_nucleic_interactions
        
        import interaction_2d_plot
        generate_2d_interaction_diagram = interaction_2d_plot.generate_2d_interaction_diagram
        
        try:
            import ligand_ligand_analyzer
            analyze_ligand_ligand_interactions = ligand_ligand_analyzer.analyze_ligand_ligand_interactions
        except ImportError:
            analyze_ligand_ligand_interactions = None
            
    except Exception as e:
        # Raise a warning but return None placeholders if running outside env (e.g. CI)
        print(f"Warning: Helper modules not found: {e}")
        # Return Nones or mock functions to allow GUI to load
        return (None,) * 14  # 14 items
        
    find_crbn_g_motif = None
    try:
        import g_motif_analyzer
        find_crbn_g_motif = g_motif_analyzer.find_crbn_g_motif
    except Exception:
        try:
            import g_motif
            find_crbn_g_motif = g_motif.find_crbn_g_motif
        except Exception:
            find_crbn_g_motif = _dynamic_load_by_filenames(
                ["g_motif_analyzer.py", "g_motif.py", "g-motif.py"], "find_crbn_g_motif"
            )
    
    # Keep helper return tuple consistent across normal and fallback paths.
    return (highlight_csv_residues, highlight_gmotif_loops,
            analyze_pdb_interactions, find_crbn_g_motif,
            render_interactions_beautifully, generate_2d_interaction_diagram,
            analyze_protein_ligand_interactions, visualize_protein_ligand_3d,
            generate_interaction_network_plot, analyze_ternary_complex,
            analyze_atom_pair_interactions, visualize_atom_pairs,
            analyze_ligand_ligand_interactions, analyze_protein_nucleic_interactions)

# Load helpers immediately
(highlight_csv_residues, highlight_gmotif_loops, analyze_pdb_interactions, 
 find_crbn_g_motif, render_interactions_beautifully, generate_2d_interaction_diagram, 
 analyze_protein_ligand_interactions, visualize_protein_ligand_3d, 
 generate_interaction_network_plot, analyze_ternary_complex, 
 analyze_atom_pair_interactions, visualize_atom_pairs, 
 analyze_ligand_ligand_interactions, analyze_protein_nucleic_interactions) = _import_helpers()
