# -*- coding: utf-8 -*-
"""
GlueTK - PyMOL Plugin for Molecular Glue Analysis
Molecular Glue vs PROTAC Classification Toolkit

Author: Vesper
Version: 1.0.0
"""

from __future__ import print_function
import locale

__version__ = "1.0.0"
__author__ = "Vesper"

# ---- 环境依赖检查与自动安装 ----
try:
    from .env_setup import ensure_dependencies
    _DEPS_OK = ensure_dependencies()
except Exception as e:
    print(f"[GlueTK] ⚠️ Dependency check failed: {e}")
    _DEPS_OK = False

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

# ---- 命令注册（与 GUI 解耦）----
def _register_commands():
    try:
        from .highlight_residues import highlight_csv_residues
        from .interaction_analyzer import (
            analyze_pdb_interactions,
            analyze_protein_ligand_interactions,
            analyze_ternary_complex,
            analyze_atom_pair_interactions,
            visualize_protein_ligand_3d,
            generate_interaction_network_plot
        )
        from .interaction_2d_plot import generate_2d_interaction_diagram
        from .binding_score import score_protein_ligand, score_ternary_complex
        from .binding_heatmap import plot_binding_heatmap
        
        # 新增: 分子胶特异功能 (PPI 分析 & Neo-表位)
        from .ppi_analyzer import (
            analyze_protein_protein_interface,
            identify_neo_epitope,
            calculate_interface_bsa,
            ppi_analyze,
            neo_epitope_find
        )
        from .g_motif_analyzer import (
            find_crbn_g_motif,
            analyze_g_motif_glue_binding
        )
        
        # Vina评分(可选,需要安装Vina)
        try:
            from .vina_scoring import vina_score_complex, compare_scoring_methods
            _vina_available = True
        except ImportError:
            _vina_available = False
    except Exception as e:
        _info(
            f"⚠️ 无法导入命令模块：{e}\n请确认 highlight_residues.py 与 interaction_analyzer.py 在插件目录。",
            f"⚠️ Failed to import command modules: {e}\nEnsure highlight_residues.py and interaction_analyzer.py exist."
        )
        return
    
    try:
        from pymol import cmd
        cmd.extend("highlight_csv_residues", highlight_csv_residues)
        cmd.extend("analyze_pdb_interactions", analyze_pdb_interactions)
        cmd.extend("analyze_protein_ligand_interactions", analyze_protein_ligand_interactions)
        cmd.extend("analyze_ternary_complex", analyze_ternary_complex)
        cmd.extend("analyze_atom_pair_interactions", analyze_atom_pair_interactions)
        cmd.extend("visualize_protein_ligand_3d", visualize_protein_ligand_3d)
        cmd.extend("generate_interaction_network_plot", generate_interaction_network_plot)
        cmd.extend("generate_2d_diagram", generate_2d_interaction_diagram)
        cmd.extend("score_protein_ligand", score_protein_ligand)
        cmd.extend("score_ternary_complex", score_ternary_complex)
        cmd.extend("plot_binding_heatmap", plot_binding_heatmap)
        
        # 分子胶特异命令
        cmd.extend("ppi_analyze", ppi_analyze)
        cmd.extend("neo_epitope_find", neo_epitope_find)
        cmd.extend("analyze_protein_protein_interface", analyze_protein_protein_interface)
        cmd.extend("identify_neo_epitope", identify_neo_epitope)
        cmd.extend("calculate_interface_bsa", calculate_interface_bsa)
        cmd.extend("find_crbn_g_motif", find_crbn_g_motif)
        cmd.extend("analyze_g_motif_glue_binding", analyze_g_motif_glue_binding)
        
        # Vina评分命令(可选)
        if _vina_available:
            cmd.extend("vina_score_complex", vina_score_complex)
            cmd.extend("compare_scoring_methods", compare_scoring_methods)
        
        _info("✅ 已注册命令：highlight_csv_residues / analyze_pdb_interactions / analyze_protein_ligand_interactions 等",
              "✅ Commands registered: highlight_csv_residues / analyze_pdb_interactions / analyze_protein_ligand_interactions etc.")
    except Exception as e:
        _info(f"⚠️ 无法注册到 PyMOL 命令空间：{e}",
              f"⚠️ Failed to register commands to PyMOL: {e}")

# ---- GUI 启动（非模态，防卡死）----
_dlg = None

def gluetk_gui():
    """启动 GlueTK 统一 GUI 窗口（非模态，不阻塞事件循环）"""
    global _dlg
    _info("[GlueTK] 正在启动统一GUI界面...", "[GlueTK] Launching unified GUI...")
    try:
        from .unified_gui import GlueTKDialog
    except Exception as e:
        _info(f"GUI 导入失败：{e}", f"Failed to import GUI: {e}")
        import traceback; traceback.print_exc()
        _print_cli_fallback()
        return

    try:
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
        _info("[GlueTK] GUI 已打开（非模态）", "[GlueTK] GUI opened (non-modal)")
    except Exception as e:
        _info(f"GUI 启动失败: {e}", f"GUI start failed: {e}")
        import traceback; traceback.print_exc()
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
    _info("🧬 GlueTK 插件已加载", "🧬 GlueTK plugin loaded")
    _register_commands()
    
    # Register GUI commands
    try:
        from pymol import cmd
        cmd.extend("gluetk_gui", gluetk_gui)
        cmd.extend("molstruct_gui", molstruct_gui)  # Legacy alias
    except Exception as e:
        print(f"Warning: Failed to register GUI command: {e}")
    
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('GlueTK - Molecular Glue Analyzer', gluetk_gui)
        _info("✅ GUI 菜单已启用：Plugins → GlueTK", "✅ GUI menu enabled: Plugins → GlueTK")
    except Exception as e:
        _info(f"⚠️ 无法添加 GUI 菜单（命令行仍可用）：{e}",
              f"⚠️ Unable to add GUI menu (CLI still available): {e}")

    print("📋 Commands:")
    print("")
    print("  🔬 Interaction Analysis:")
    print("    • analyze_pdb_interactions - Analyze interactions in a PDB structure (strict standards)")
    print("    • analyze_protein_ligand_interactions - Analyze protein-ligand interactions (strict standards)")
    print("    • analyze_ternary_complex - Analyze ternary complex")
    print("    • analyze_atom_pair_interactions - Analyze atom pair interactions")
    print("")
    print("  ✨ Molecular Glue Analysis (NEW):")
    print("    • ppi_analyze - Analyze protein-protein interface (PPI)")
    print("    • neo_epitope_find - Identify neo-substrate epitope")
    print("    • analyze_g_motif_glue_binding - Validate G-motif as glue substrate")
    print("    • find_crbn_g_motif - Find CRBN G-motif/G-loop")
    print("")
    print("  ⚡ Binding Energy Scoring:")
    print("    • score_protein_ligand - Calculate binding energy for protein-ligand complex (empirical)")
    print("    • score_ternary_complex - Calculate binding energy for molecular glue/PROTAC (with cooperativity)")
    print("    • plot_binding_heatmap - Generate batch binding energy heatmap from CSV files")
    try:
        from .vina_scoring import check_vina_available
        if check_vina_available():
            print("    • vina_score_complex - Calculate binding energy using AutoDock Vina (requires vina installed)")
            print("    • compare_scoring_methods - Compare empirical vs Vina scoring")
    except:
        pass
    print("")
    print("  🎨 Visualization:")
    print("    • highlight_csv_residues - Highlight residue interactions from CSV")
    print("    • visualize_protein_ligand_3d - 3D visualization of protein-ligand interactions")
    print("    • generate_interaction_network_plot - Generate interaction network plot")
    print("")
    print("  🖥️  GUI:")
    print("    • gluetk_gui - Open GlueTK unified analysis GUI (all features)")
    print("    • molstruct_gui - Legacy alias for gluetk_gui (for backward compatibility)")
    print("")
    print("💡 Use help(command_name) for details")
    print("⭐ Using strict standards: H-bond ≤2.8Å, Salt bridge ≤4.0Å, suitable for publication")
