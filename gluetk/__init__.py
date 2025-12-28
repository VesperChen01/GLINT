# -*- coding: utf-8 -*-
"""
GlueTK - PyMOL Plugin for Molecular Glue Analysis
Molecular Glue vs PROTAC Classification Toolkit

Author: Vesper
Version: v0.1.14-beta-contact-immersive-bg-minimalist-contact-height-fix-contact-final-en-fix-contact-final-v2-hotfix-qcolor-contact-redesign-slogan-fix-final-visuals-polished-hotfix-v2-hotfix
"""

from __future__ import print_function
import locale

__version__ = "v0.1.14-beta"
__author__ = "Vesper"

# ---- 环境依赖检查 ----
# 延迟依赖检查，避免在导入时触发 PyQt 崩溃
_DEPS_OK = False
_DEPS_CHECKED = False

def _check_deps_safe():
    """安全地检查依赖（延迟到实际需要时）"""
    global _DEPS_OK, _DEPS_CHECKED
    if _DEPS_CHECKED:
        return _DEPS_OK
    
    try:
        from .env_checker import ensure_dependencies
        _DEPS_OK = ensure_dependencies(silent=True)  # 静默检查
        _DEPS_CHECKED = True
    except Exception as e:
        _DEPS_OK = False
        _DEPS_CHECKED = True
    
    return _DEPS_OK

# ---- 语言工具 ----
def _zh():
    try:
        lang = (locale.getdefaultlocale() or ["en"])[0] or "en"
        return str(lang).lower().startswith("zh")
    except Exception:
        return False

def _info(cn, en):
    # Force English output for all prompts
    print(en)

# Global flag for Vina availability
_vina_available = False

# ---- 命令注册（与 GUI 解耦）----
def _register_commands():
    global _vina_available
    try:
        from .highlight_residues import highlight_csv_residues
        from .interaction_analyzer import (
            analyze_pdb_interactions,
            analyze_protein_ligand_interactions,
            analyze_protein_nucleic_interactions,
            analyze_ternary_complex,
            analyze_atom_pair_interactions,
            visualize_protein_ligand_3d,
            generate_interaction_network_plot,
            toggle_interaction_lines
        )
        from .interaction_2d_plot import generate_2d_interaction_diagram
        from .binding_heatmap import plot_binding_heatmap
        
        # Ligand-Ligand Interaction
        try:
            from .ligand_ligand_analyzer import analyze_ligand_ligand_interactions
        except ImportError:
            analyze_ligand_ligand_interactions = None
        
        # 新增: 分子胶特异功能 (PPI 分析 & Neo-表位)
        from .ppi_analyzer import (
            analyze_protein_protein_interface,
            identify_neo_epitope,
            calculate_interface_bsa,
            visualize_ppi_interface,
            visualize_neo_epitope,
            ppi_analyze,
            neo_epitope_find
        )
        from .g_motif_analyzer import (
            find_crbn_g_motif,
            analyze_g_motif_glue_binding,
            validate_crbn_hbonds,
            validate_g_motif_geometry,
            degron_annotate
        )
        

        
        # 批量分析模块
        try:
            from .batch_analyzer import (
                batch_gmotif,
                batch_ppi,
                batch_pockets,
                batch_interactions,
                BatchAnalyzer
            )
            _batch_available = True
        except ImportError as e:
            print(f"⚠️ Batch analyzer not available: {e}")
            _batch_available = False
        
        # 分子胶数据库模块
        try:
            from .glue_database import (
                GlueDatabase,
                glue_db_search,
                glue_db_info,
                glue_db_fetch,
                glue_db_export
            )
            _glue_db_available = True
        except ImportError as e:
            print(f"⚠️ Glue database not available: {e}")
            _glue_db_available = False
        
        # 分子胶设计分析（Ternary complex 建模）
        from .glue_design_analyzer import (
            align_gloop_for_modeling,
            detect_clashes_at_interface,
            identify_exit_vectors,
            analyze_electrostatic_environment,
            comprehensive_glue_design_analysis
        )
        
        # 配体电性互补性分析 (Electrostatic Complementarity)
        try:
            from .ligand_ec_calculator import (
                calculate_ligand_ec,
                analyze_ternary_ec,
                compare_ligand_ec,
                analyze_multiconformer_ec,
                calculate_ec_hotspots,
                analyze_substituent_ec_effect,
                ECCalculator,
                DXGrid,
                LigandSurfaceSampler,
                GasteigerChargeCalculator
            )
            _ec_available = True
        except ImportError as e:
            print(f"⚠️ EC Calculator not available: {e}")
            _ec_available = False
        
        # 口袋检测与分析
        from .pocket_detector import detect_pockets, compare_pockets
        from .pocket_visualizer import (
            visualize_pockets,
            show_pocket_labels,
            visualize_pocket_comparison,
            overlay_pocket_electrostatics,
            visualize_pockets_with_interactions
        )
        from .pocket_glue_integration import (
            analyze_pockets_in_ppi_interface,
            analyze_pockets_with_glue,
            correlate_pockets_with_interactions,
            integrate_pockets_with_electrostatics,
            comprehensive_glue_pocket_analysis,
            comprehensive_gmotif_pocket_analysis
        )
        
        # 突变分析模块
        from .mutation_analyzer import (
            perform_mutation,
            minimize_energy,
            calculate_mutation_ddg,
            analyze_mutation_effects,
            ddg_heatmap
        )
        
        # 表面相似性与互补性分析模块
        try:
            from .surface_similarity import (
                analyze_surface_similarity,
                analyze_surface_complementarity,
                SurfaceSimilarityAnalyzer
            )
            _surface_similarity_available = True
        except ImportError as e:
            print(f"⚠️ Surface similarity analysis not available: {e}")
            _surface_similarity_available = False
        
        # Open Targets 疾病靶点分析
        try:
            from .disease_target_commands import (
                ot_disease_targets,
                ot_glue_insight
            )
            _disease_analysis_available = True
        except ImportError as e:
            print(f"⚠️ Disease analysis not available: {e}")
            _disease_analysis_available = False
        
        # Vina集成(可选,需要安装Vina)
        try:
            from .vina_integration import (
                vina_score_complex,
                compare_scoring_methods,
                pocket_based_docking
            )
            _vina_available = True
        except ImportError:
            _vina_available = False
    except Exception as e:
        _info(
            f"⚠️ 插件加载失败：{e}",
            f"⚠️ Plugin load failed: {e}"
        )
        return
    
    try:
        from pymol import cmd
        cmd.extend("highlight_csv_residues", highlight_csv_residues)
        cmd.extend("analyze_pdb_interactions", analyze_pdb_interactions)
        cmd.extend("analyze_protein_ligand_interactions", analyze_protein_ligand_interactions)
        cmd.extend("analyze_protein_nucleic_interactions", analyze_protein_nucleic_interactions)
        cmd.extend("analyze_ternary_complex", analyze_ternary_complex)
        cmd.extend("analyze_atom_pair_interactions", analyze_atom_pair_interactions)
        cmd.extend("visualize_protein_ligand_3d", visualize_protein_ligand_3d)
        cmd.extend("toggle_interaction_lines", toggle_interaction_lines)
        cmd.extend("generate_interaction_network_plot", generate_interaction_network_plot)
        cmd.extend("generate_2d_diagram", generate_2d_interaction_diagram)
        cmd.extend("plot_binding_heatmap", plot_binding_heatmap)
        if analyze_ligand_ligand_interactions:
            cmd.extend("analyze_ligand_ligand_interactions", analyze_ligand_ligand_interactions)
        
        # 分子胶特异命令
        cmd.extend("ppi_analyze", ppi_analyze)
        cmd.extend("neo_epitope_find", neo_epitope_find)
        cmd.extend("analyze_protein_protein_interface", analyze_protein_protein_interface)
        cmd.extend("identify_neo_epitope", identify_neo_epitope)
        cmd.extend("calculate_interface_bsa", calculate_interface_bsa)
        cmd.extend("visualize_ppi_interface", visualize_ppi_interface)
        cmd.extend("visualize_neo_epitope", visualize_neo_epitope)
        cmd.extend("find_crbn_g_motif", find_crbn_g_motif)
        cmd.extend("analyze_g_motif_glue_binding", analyze_g_motif_glue_binding)
        cmd.extend("validate_crbn_hbonds", validate_crbn_hbonds)
        cmd.extend("validate_g_motif_geometry", validate_g_motif_geometry)
        cmd.extend("degron_annotate", degron_annotate)
        

        
        # 分子胶设计分析命令
        cmd.extend("align_gloop_for_modeling", align_gloop_for_modeling)
        cmd.extend("detect_clashes_at_interface", detect_clashes_at_interface)
        cmd.extend("identify_exit_vectors", identify_exit_vectors)
        cmd.extend("analyze_electrostatic_environment", analyze_electrostatic_environment)
        cmd.extend("comprehensive_glue_design_analysis", comprehensive_glue_design_analysis)
        
        # 口袋检测命令
        cmd.extend("detect_pockets", detect_pockets)
        cmd.extend("compare_pockets", compare_pockets)
        cmd.extend("visualize_pockets", visualize_pockets)
        cmd.extend("show_pocket_labels", show_pocket_labels)
        cmd.extend("visualize_pocket_comparison", visualize_pocket_comparison)
        cmd.extend("overlay_pocket_electrostatics", overlay_pocket_electrostatics)
        cmd.extend("visualize_pockets_with_interactions", visualize_pockets_with_interactions)
        
        # 口袋-分子胶联动命令
        cmd.extend("analyze_pockets_in_ppi_interface", analyze_pockets_in_ppi_interface)
        cmd.extend("analyze_pockets_with_glue", analyze_pockets_with_glue)
        cmd.extend("correlate_pockets_with_interactions", correlate_pockets_with_interactions)
        cmd.extend("integrate_pockets_with_electrostatics", integrate_pockets_with_electrostatics)
        cmd.extend("comprehensive_glue_pocket_analysis", comprehensive_glue_pocket_analysis)
        cmd.extend("comprehensive_gmotif_pocket_analysis", comprehensive_gmotif_pocket_analysis)
        
        # 突变分析命令
        cmd.extend("perform_mutation", perform_mutation)
        cmd.extend("minimize_energy", minimize_energy)
        cmd.extend("calculate_mutation_ddg", calculate_mutation_ddg)
        cmd.extend("analyze_mutation_effects", analyze_mutation_effects)
        cmd.extend("ddg_heatmap", ddg_heatmap)
        
        # Open Targets 疾病靶点分析命令
        if _disease_analysis_available:
            cmd.extend("ot_disease_targets", ot_disease_targets)
            cmd.extend("ot_glue_insight", ot_glue_insight)
        
        # Vina集成命令(可选)
        if _vina_available:
            cmd.extend("vina_score_complex", vina_score_complex)
            cmd.extend("compare_scoring_methods", compare_scoring_methods)
            cmd.extend("pocket_based_docking", pocket_based_docking)
        
        # 批量分析命令
        if _batch_available:
            cmd.extend("batch_gmotif", batch_gmotif)
            cmd.extend("batch_ppi", batch_ppi)
            cmd.extend("batch_pockets", batch_pockets)
            cmd.extend("batch_interactions", batch_interactions)
        
        # 表面相似性与互补性分析命令
        if _surface_similarity_available:
            cmd.extend("analyze_surface_similarity", analyze_surface_similarity)
            cmd.extend("analyze_surface_complementarity", analyze_surface_complementarity)
        
        # 分子胶数据库命令
        if _glue_db_available:
            cmd.extend("glue_db_search", glue_db_search)
            cmd.extend("glue_db_info", glue_db_info)
            cmd.extend("glue_db_fetch", glue_db_fetch)
            cmd.extend("glue_db_export", glue_db_export)
        
        # 配体电性互补性分析命令
        if _ec_available:
            cmd.extend("calculate_ligand_ec", calculate_ligand_ec)
            cmd.extend("analyze_ternary_ec", analyze_ternary_ec)
            cmd.extend("compare_ligand_ec", compare_ligand_ec)
            cmd.extend("analyze_multiconformer_ec", analyze_multiconformer_ec)
            cmd.extend("calculate_ec_hotspots", calculate_ec_hotspots)
            cmd.extend("analyze_substituent_ec_effect", analyze_substituent_ec_effect)
        
        # 静默注册,避免终端输出过多
        # _info("✅ 已注册命令", "✅ Commands registered")
    except Exception as e:
        _info(f"⚠️ 无法注册到 PyMOL 命令空间：{e}",
              f"⚠️ Failed to register commands to PyMOL: {e}")

# ---- GUI 启动（非模态，防卡死）----
_dlg = None

def _import_gui_dialog():
    """尝试导入 GlueTKDialog (优先使用新的模块化 GUI)

    兼容两种加载方式:
    1) PyMOL 插件机制从 ~/.pymol/startup/gluetk 导入包
    2) 用户/Launcher 直接 `pymol gluetk/__init__.py` 作为脚本运行

    在某些情况下，__file__ 可能被解析到 PyMOL 自己的目录 (site-packages/pymol)，
    因此这里增加多重路径修正逻辑，尽量找到真正的 gluetk 根目录
    """
    import sys
    import os
    
    def _looks_like_plugin_root(path: str) -> bool:
        return (
            bool(path)
            and os.path.isdir(path)
            and os.path.exists(os.path.join(path, "gui"))
            and os.path.exists(os.path.join(path, "env_checker.py"))
        )
    
    # 1) 首选: 基于当前 __file__ 推断
    plugin_dir = os.path.dirname(os.path.abspath(__file__))
    if not _looks_like_plugin_root(plugin_dir):
        print(f"[GlueTK Debug] __file__ path suspicious: {plugin_dir}")
        # 2) 退而求其次: 使用 inspect 获取真实源文件路径
        try:
            import inspect
            frame = inspect.currentframe()
            if frame is not None:
                file_from_frame = inspect.getfile(frame)
                cand = os.path.dirname(os.path.abspath(file_from_frame))
                if _looks_like_plugin_root(cand):
                    plugin_dir = cand
                    print(f"[GlueTK Debug] Corrected plugin_dir via inspect: {plugin_dir}")
        except Exception as e:
            print(f"[GlueTK Debug] Inspect failed: {e}")
    
    # 3) 仍然不对: 搜索常见安装路径
    if not _looks_like_plugin_root(plugin_dir):
        home = os.path.expanduser("~")
        candidates = [
            os.path.join(home, ".pymol", "startup", "gluetk"),
            os.path.join(home, "pymol", "startup", "gluetk"),
        ]
        # 开发者环境: 当前工作目录下的 gluetk 目录
        cwd = os.getcwd()
        candidates.append(os.path.join(cwd, "gluetk"))
        
        for cand in candidates:
            if _looks_like_plugin_root(cand):
                plugin_dir = cand
                print(f"[GlueTK Debug] Found plugin root candidate: {plugin_dir}")
                break
    
    parent_dir = os.path.dirname(plugin_dir)
    
    # Debug: 打印路径信息
    print(f"[GlueTK Debug] plugin_dir = {plugin_dir}")
    print(f"[GlueTK Debug] parent_dir = {parent_dir}")
    print(f"[GlueTK Debug] parent_dir in sys.path? {parent_dir in sys.path}")
    
    if parent_dir and parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)
        print(f"[GlueTK Debug] Added {parent_dir} to sys.path")
    
    # 使用绝对导入 gluetk.gui.main_window
    try:
        import gluetk.gui.main_window as gui_module
        return gui_module.GlueTKDialog
    except ImportError as e:
        print(f"❌ Error importing modular GUI: {e}")
        print(f"[GlueTK Debug] sys.path = {sys.path[:5]}")
        import traceback
        traceback.print_exc()
        return None

def _check_qt_safe():
    """安全地检查 Qt 是否可用（使用子进程避免崩溃）"""
    import subprocess
    import sys
    import os
    
    # [macOS Fix] 设置环境变量以避免某些 Qt 绘图引起的崩溃
    if sys.platform == "darwin":
        os.environ["QT_MAC_WANTS_LAYER"] = "1"
        # 兼容性修复：避免 macOS 上的 OpenMP 冲突和多线程驱动问题
        os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
        if "OMP_NUM_THREADS" not in os.environ:
            os.environ["OMP_NUM_THREADS"] = "1"

    # 使用子进程测试常见的 Qt 绑定
    for binding in ["PyQt6", "PySide6", "PyQt5", "PySide2"]:
        test_code = f"import {binding}; print('OK')"
        try:
            result = subprocess.run(
                [sys.executable, "-c", test_code],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0 and "OK" in result.stdout:
                return True
        except:
            pass
    
    return False

def _get_qapp():
    """获取或创建 QApplication 实例 (针对 macOS/PyMOL 优化)"""
    try:
        from .gui.qt_adapter import QtWidgets
        if QtWidgets is None:
            return None
        app = QtWidgets.QApplication.instance()
        if app is None:
            # 如果 PyMOL 还没初始化 Qt 循环，我们也不强制创建，除非真的需要
            # 在某些 macOS 环境下，直接创建 QApplication 会导致 Segfault
            app = QtWidgets.QApplication([])
        return app
    except Exception as e:
        print(f"[GlueTK] Failed to get QApplication: {e}")
        return None

def gluetk_gui():
    """启动 GlueTK 统一 GUI 窗口（非模态，不阻塞事件循环）"""
    global _dlg
    
    # 安全地检查 PyQt 是否可用（避免在主进程中导入导致崩溃）
    print("[GlueTK] Checking Qt availability...")
    if not _check_qt_safe():
        _info(
            "Qt 绑定 (PyQt5/6/PySide2/6) 未安装或无法使用，无法启动 GUI",
            "Qt binding (PyQt5/6/PySide2/6) not installed or unavailable, cannot start GUI"
        )
        print("\n💡 Install PyQt5 or PyQt6 to use the GUI:")
        print("   conda install -c conda-forge pyqt -y")
        print("   # or")
        print("   pip install PyQt5")
        _print_cli_fallback()
        return
    
    print("[GlueTK] Qt is available, loading GUI...")
    GlueTKDialog = _import_gui_dialog()
    if GlueTKDialog is None:
        _info("GUI 导入失败", "Failed to import GUI module")
        import sys, os
        pkg_dir = os.path.dirname(os.path.realpath(__file__))
        print(f"  Package dir: {pkg_dir}")
        print(f"  gui/main_window.py exists: {os.path.exists(os.path.join(pkg_dir, 'gui', 'main_window.py'))}")
        print(f"  sys.path[0:3]: {sys.path[:3]}")
        _print_cli_fallback()
        return

    try:
        # [macOS Fix] 确保获取现有的 QApplication 实例，避免重复初始化导致的 Segfault
        app = _get_qapp()
        
        # 已有窗口则激活
        if _dlg is not None:
            try:
                _dlg.show(); _dlg.raise_(); _dlg.activateWindow()
                return
            except Exception:
                _dlg = None
        
        # 新建并非模态展示
        _dlg = GlueTKDialog()
        _dlg.setModal(False)
        _dlg.show()
        _dlg.raise_()
        _dlg.activateWindow()
    except Exception as e:
        _info(f"GUI 启动失败: {e}", f"GUI start failed: {e}")
        import traceback; traceback.print_exc()
        print("\n💡 If you see a segmentation fault:")
        print("   1. Make sure PyQt5 is properly installed in your conda environment")
        print("   2. Try: conda install -c conda-forge pyqt --force-reinstall")
        print("   3. Restart PyMOL after reinstalling PyQt5")
        _print_cli_fallback()

# Backward compatibility alias
def molstruct_gui():
    """Legacy alias for gluetk_gui() - for backward compatibility"""
    return gluetk_gui()

def _print_cli_fallback():
    _info("请使用命令行：", "Use CLI instead:")
    print("  highlight_csv_residues csv_path='file.csv', obj='object'")
    print("  analyze_pdb_interactions obj_name='object', output_csv='output.csv'")

# ---- 插件入口 ----
def __init_plugin__(app=None):
    """
    PyMOL 插件入口函数
    
    工作流程:
    1. 检查依赖 (已在模块加载时完成)
    2. 如果依赖正常，注册所有命令
    3. 总是注册 GUI 命令（用于显示错误信息）
    """
    # 只有依赖检查通过时才注册全部命令
    if _check_deps_safe():
        _register_commands()
    
    # 总是注册 GUI 命令
    try:
        from pymol import cmd
        cmd.extend("gluetk_gui", gluetk_gui)
        cmd.extend("molstruct_gui", molstruct_gui)  # 兼容别名
    except Exception as e:
        print(f"Warning: Failed to register GUI command: {e}")
    
    # 添加菜单项
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('GlueTK - Molecular Glue Analyzer', gluetk_gui)
    except Exception as e:
        pass  # 静默处理

    # 欢迎信息（延迟检查）
    if _check_deps_safe():
        print("\n🧬 GlueTK - Molecular Glue Analyzer v0.1.14-beta")
        print("┌" + "─" * 48 + "┐")
        print("│  Quick Start:                                   │")
        print("│    • gluetk_gui            - Launch GUI          │")
        print("│    • help(gluetk_gui)       - Show help         │")
        print("│    • Plugins → GlueTK       - Menu access       │")
        print("└" + "─" * 48 + "┘")
    else:
        print("\n🧬 GlueTK v0.1.14-beta - ⚠️  Some dependencies may be missing")
        print("💡 Most features are available. Use 'gluetk_gui' to launch GUI.\n")

# Auto-register if running within PyMOL environment (e.g. via 'run' command or import)
# IMPORTANT: Only auto-register in actual PyMOL environment, not in standalone Python
try:
    import pymol
    if hasattr(pymol, 'cmd') and hasattr(pymol, 'stored'):
        # Double-check we're in actual PyMOL environment
        # Avoid re-registering if already registered (check one key command)
        if 'gluetk_gui' not in pymol.cmd.keyword:
            __init_plugin__()
except Exception as e:
    # Silently fail if not in PyMOL environment
    pass
