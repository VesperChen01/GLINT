# -*- coding: utf-8 -*-
"""
GLINT - PyMOL Plugin for Molecular Glue Interface Analysis
Glue Interface Analyzer

Author: Roufen Chen
Version:
"""

from __future__ import print_function
import locale
import os
import sys
from glint.path_utils import get_conda_search_roots, get_conda_site_packages_patterns


# Fix relative import: ensure package path is correct when executed via `run` command
_this_dir = os.path.dirname(os.path.abspath(__file__))
_parent_dir = os.path.dirname(_this_dir)
if _parent_dir not in sys.path:
    sys.path.insert(0, _parent_dir)

# Inject conda environment site-packages into sys.path
# When launched via GLINT.app, conda is activated (CONDA_PREFIX is set),
# but PyMOL's embedded Python doesn't inherit conda's site-packages.
def _inject_conda_site_packages():
    import glob as _glob
    
    _debug = os.environ.get('GLINT_DEBUG', '')
    def _log(msg):
        if _debug:
            print(f"[GLINT Debug] {msg}")
    
    # 0. Check if PYTHONPATH already provides site-packages (e.g. set by launcher)
    _pythonpath = os.environ.get('PYTHONPATH', '')
    for _pp in _pythonpath.split(os.pathsep):
        if _pp.endswith('site-packages') and os.path.isdir(_pp) and _pp not in sys.path:
            sys.path.insert(0, _pp)
            _log(f"Injected from PYTHONPATH: {_pp}")
    
    # 1. Try CONDA_PREFIX (set by conda activate)
    _prefix = os.environ.get('CONDA_PREFIX', '')
    _log(f"CONDA_PREFIX={_prefix!r}")
    
    # 2. Fallback: infer from sys.prefix if it looks like a conda env
    if not _prefix:
        _sys_prefix = sys.prefix
        if os.path.isfile(os.path.join(_sys_prefix, 'conda-meta', 'history')):
            _prefix = _sys_prefix
            _log(f"Inferred prefix from sys.prefix: {_prefix}")
    
    # 3. Fallback: find glint env directly from conda base
    if not _prefix:
        _conda_base = os.environ.get('CONDA_EXE', '')
        if _conda_base:
            _conda_base = os.path.dirname(os.path.dirname(_conda_base))
        if not _conda_base:
            for _candidate in get_conda_search_roots():
                if os.path.isdir(_candidate):
                    _conda_base = _candidate
                    _log(f"Found conda base: {_conda_base}")
                    break
        if _conda_base:
            _env_path = os.path.join(_conda_base, 'envs', 'glint')
            if os.path.isdir(_env_path):
                _prefix = _env_path
                _log(f"Found glint env: {_prefix}")
    
    if not _prefix:
        _log("No conda prefix found, skipping injection")
        return
    
    _patterns = get_conda_site_packages_patterns(_prefix)
    for _pat in _patterns:
        for _sp in _glob.glob(_pat):
            if _sp not in sys.path:
                sys.path.insert(0, _sp)
                _log(f"Injected: {_sp}")

_inject_conda_site_packages()


# Import version from _version.py (supports two loading methods)
try:
    from ._version import __version__
except ImportError:
    try:
        from glint._version import __version__
    except ImportError:
        # Last resort: read version file directly
        _version_file = os.path.join(_this_dir, "_version.py")
        if os.path.exists(_version_file):
            __version__ = "unknown"
            with open(_version_file, "r") as f:
                for line in f:
                    if line.startswith("__version__"):
                        __version__ = line.split("=")[1].strip().strip('"\'')
                        break
        else:
            __version__ = "unknown"

__author__ = "Roufen Chen"

# ---- Dependency check ----
# Deferred dependency check to avoid triggering PyQt crash on import
_DEPS_OK = False
_DEPS_CHECKED = False

def _check_deps_safe():
    """
    Safely check dependencies (deferred until actually needed).
    
    Includes detailed error logging to help debug import issues.
    """
    global _DEPS_OK, _DEPS_CHECKED
    if _DEPS_CHECKED:
        return _DEPS_OK
    
    try:
        from .env_checker import ensure_dependencies
        _DEPS_OK = ensure_dependencies(silent=True)  # Silent check
        _DEPS_CHECKED = True
    except ImportError as e:
        # Detailed import error logging
        import traceback
        print(f"[GLINT] ⚠️ Failed to import dependency checker: {e}")
        print(f"[GLINT] Detailed error:")
        traceback.print_exc()
        _DEPS_OK = False
        _DEPS_CHECKED = True
    except Exception as e:
        # Other errors
        import traceback
        print(f"[GLINT] ⚠️ Error during dependency check: {e}")
        print(f"[GLINT] Detailed error:")
        traceback.print_exc()
        _DEPS_OK = False
        _DEPS_CHECKED = True
    
    return _DEPS_OK

# ---- Language utility ----
def _zh():
    """Detect whether the current environment is Chinese"""
    try:
        import os as _os
        # locale.getdefaultlocale() is deprecated in Python 3.11+
        # Check environment variables first, then fall back to locale.getlocale()
        for env_var in ('LANG', 'LANGUAGE', 'LC_ALL', 'LC_MESSAGES'):
            val = _os.environ.get(env_var, '')
            if val.lower().startswith('zh'):
                return True
        loc = locale.getlocale()
        if loc and loc[0] and str(loc[0]).lower().startswith('zh'):
            return True
        return False
    except Exception:
        return False

def _info(cn, en):
    # Force English output for all prompts
    print(en)

# Global flag for Vina availability
_vina_available = False

# ---- Command registration (decoupled from GUI) ----
def _register_commands():
    """
    Register all PyMOL commands.
    
    Includes detailed import error logging to help debug module loading issues.
    """
    global _vina_available
    import traceback
    
    # Track modules that failed to import
    _import_errors = []
    
    def _safe_import(module_name, items=None):
        """Safely import a module, logging errors without interrupting"""
        try:
            if items:
                module = __import__(module_name, globals(), locals(), items, 1)
                return tuple(getattr(module, item) for item in items)
            else:
                return __import__(module_name, globals(), locals(), [], 1)
        except ImportError as e:
            _import_errors.append((module_name, str(e)))
            print(f"[GLINT] ⚠️ Module import failed: {module_name}")
            print(f"[GLINT]   Error: {e}")
            return None if not items else tuple([None] * len(items))
        except Exception as e:
            _import_errors.append((module_name, str(e)))
            print(f"[GLINT] ❌ Module import exception: {module_name}")
            print(f"[GLINT]   Error: {e}")
            traceback.print_exc()
            return None if not items else tuple([None] * len(items))
    
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
        except ImportError as e:
            print(f"[GLINT] ℹ️ ligand_ligand_analyzer not available: {e}")
            analyze_ligand_ligand_interactions = None
        
        # Molecular glue specific features (PPI analysis & Neo-epitope)
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
        

        
        # Batch analysis module
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
        
        # Molecular glue design analysis (Ternary complex modeling)
        from .glue_design_analyzer import (
            align_gloop_for_modeling,
            detect_clashes_at_interface,
            identify_exit_vectors,
            analyze_electrostatic_environment,
            comprehensive_glue_design_analysis
        )
        
        # Ligand electrostatic complementarity analysis (Electrostatic Complementarity)
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
        
        # Pocket detection and analysis
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
        
        # Mutation analysis module
        from .mutation_analyzer import (
            perform_mutation,
            minimize_energy,
            calculate_mutation_ddg,
            analyze_mutation_effects,
            ddg_heatmap
        )
        
        # PyMOL visualization styles
        try:
            from .pymol_styles import (
                apply_professional_style,
                setup_protein_cartoon,
                setup_protein_surface,
                setup_ligand_sticks,
                setup_binding_site,
                apply_publication_figure_style,
                quick_protein_view,
                quick_ligand_view
            )
            _pymol_styles_available = True
        except ImportError as e:
            print(f"⚠️ PyMOL styles not available: {e}")
            _pymol_styles_available = False
        
        # Surface similarity and complementarity analysis module
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
        
        # Vina integration (optional, requires Vina installation)
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
            f"⚠️ Plugin load failed: {e}",
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
        
        # Molecular glue specific commands
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
        

        
        # Molecular glue design analysis commands
        cmd.extend("align_gloop_for_modeling", align_gloop_for_modeling)
        cmd.extend("detect_clashes_at_interface", detect_clashes_at_interface)
        cmd.extend("identify_exit_vectors", identify_exit_vectors)
        cmd.extend("analyze_electrostatic_environment", analyze_electrostatic_environment)
        cmd.extend("comprehensive_glue_design_analysis", comprehensive_glue_design_analysis)
        
        # Pocket detection commands
        cmd.extend("detect_pockets", detect_pockets)
        cmd.extend("compare_pockets", compare_pockets)
        cmd.extend("visualize_pockets", visualize_pockets)
        cmd.extend("show_pocket_labels", show_pocket_labels)
        cmd.extend("visualize_pocket_comparison", visualize_pocket_comparison)
        cmd.extend("overlay_pocket_electrostatics", overlay_pocket_electrostatics)
        cmd.extend("visualize_pockets_with_interactions", visualize_pockets_with_interactions)
        
        # Pocket-glue integration commands
        cmd.extend("analyze_pockets_in_ppi_interface", analyze_pockets_in_ppi_interface)
        cmd.extend("analyze_pockets_with_glue", analyze_pockets_with_glue)
        cmd.extend("correlate_pockets_with_interactions", correlate_pockets_with_interactions)
        cmd.extend("integrate_pockets_with_electrostatics", integrate_pockets_with_electrostatics)
        cmd.extend("comprehensive_glue_pocket_analysis", comprehensive_glue_pocket_analysis)
        cmd.extend("comprehensive_gmotif_pocket_analysis", comprehensive_gmotif_pocket_analysis)
        
        # Mutation analysis commands
        cmd.extend("perform_mutation", perform_mutation)
        cmd.extend("minimize_energy", minimize_energy)
        cmd.extend("calculate_mutation_ddg", calculate_mutation_ddg)
        cmd.extend("analyze_mutation_effects", analyze_mutation_effects)
        cmd.extend("ddg_heatmap", ddg_heatmap)

        # Vina integration commands (optional)
        if _vina_available:
            cmd.extend("vina_score_complex", vina_score_complex)
            cmd.extend("compare_scoring_methods", compare_scoring_methods)
            cmd.extend("pocket_based_docking", pocket_based_docking)
        
        # Batch analysis commands
        if _batch_available:
            cmd.extend("batch_gmotif", batch_gmotif)
            cmd.extend("batch_ppi", batch_ppi)
            cmd.extend("batch_pockets", batch_pockets)
            cmd.extend("batch_interactions", batch_interactions)
        
        # Surface similarity and complementarity analysis commands
        if _surface_similarity_available:
            cmd.extend("analyze_surface_similarity", analyze_surface_similarity)
            cmd.extend("analyze_surface_complementarity", analyze_surface_complementarity)

        # Ligand electrostatic complementarity analysis commands
        if _ec_available:
            cmd.extend("calculate_ligand_ec", calculate_ligand_ec)
            cmd.extend("analyze_ternary_ec", analyze_ternary_ec)
            cmd.extend("compare_ligand_ec", compare_ligand_ec)
            cmd.extend("analyze_multiconformer_ec", analyze_multiconformer_ec)
            cmd.extend("calculate_ec_hotspots", calculate_ec_hotspots)
            cmd.extend("analyze_substituent_ec_effect", analyze_substituent_ec_effect)
        
        # Silent registration, avoid excessive terminal output
        # _info("✅ Commands registered", "✅ Commands registered")
    except Exception as e:
        _info(f"⚠️ Failed to register commands to PyMOL: {e}",
              f"⚠️ Failed to register commands to PyMOL: {e}")

# ---- GUI launch (non-modal, prevent freezing) ----
_dlg = None

def _import_gui_dialog():
    """Try to import GLINTDialog (prefer the new modular GUI)

    Compatible with two loading methods:
    1) PyMOL plugin mechanism importing from ~/.pymol/startup/glint
    2) User/Launcher directly running `pymol glint/__init__.py` as a script

    In some cases, __file__ may resolve to PyMOL's own directory (site-packages/pymol),
    so multiple path correction logic is added here to find the real glint root directory.
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
    
    # 1) Preferred: infer from current __file__
    plugin_dir = os.path.dirname(os.path.abspath(__file__))
    if not _looks_like_plugin_root(plugin_dir):
        print(f"[GLINT Debug] __file__ path suspicious: {plugin_dir}")
        # 2) Fallback: use inspect to get the real source file path
        try:
            import inspect
            frame = inspect.currentframe()
            if frame is not None:
                file_from_frame = inspect.getfile(frame)
                cand = os.path.dirname(os.path.abspath(file_from_frame))
                if _looks_like_plugin_root(cand):
                    plugin_dir = cand
                    print(f"[GLINT Debug] Corrected plugin_dir via inspect: {plugin_dir}")
        except Exception as e:
            print(f"[GLINT Debug] Inspect failed: {e}")
    
    # 3) Still incorrect: search common installation paths
    if not _looks_like_plugin_root(plugin_dir):
        home = os.path.expanduser("~")
        candidates = [
            os.path.join(home, ".pymol", "startup", "glint"),
            os.path.join(home, "pymol", "startup", "glint"),
        ]
        # Developer environment: glint directory under current working directory
        cwd = os.getcwd()
        candidates.append(os.path.join(cwd, "glint"))
        
        for cand in candidates:
            if _looks_like_plugin_root(cand):
                plugin_dir = cand
                print(f"[GLINT Debug] Found plugin root candidate: {plugin_dir}")
                break
    
    parent_dir = os.path.dirname(plugin_dir)
    
    # Debug: print path information
    print(f"[GLINT Debug] plugin_dir = {plugin_dir}")
    print(f"[GLINT Debug] parent_dir = {parent_dir}")
    print(f"[GLINT Debug] parent_dir in sys.path? {parent_dir in sys.path}")
    
    if parent_dir and parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)
        print(f"[GLINT Debug] Added {parent_dir} to sys.path")
    
    # Use absolute import glint.gui.main_window
    try:
        import glint.gui.main_window as gui_module
        return gui_module.GLINTDialog
    except ImportError as e:
        print(f"❌ Error importing modular GUI: {e}")
        print(f"[GLINT Debug] sys.path = {sys.path[:5]}")
        import traceback
        traceback.print_exc()
        return None

def _check_qt_safe():
    """Safely check if Qt is available (using subprocess to avoid crashes)"""
    import subprocess
    import sys
    import os
    
    # [macOS Fix] Set environment variables to avoid Qt rendering crashes
    if sys.platform == "darwin":
        os.environ["QT_MAC_WANTS_LAYER"] = "1"
        # Compatibility fix: avoid OpenMP conflicts and multi-threaded driver issues on macOS
        os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
        if "OMP_NUM_THREADS" not in os.environ:
            os.environ["OMP_NUM_THREADS"] = "1"

    # Use subprocess to test common Qt bindings
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
        except (OSError, subprocess.SubprocessError):  # Subprocess call may fail
            pass
    
    return False

def _get_qapp():
    """Get or create QApplication instance (optimized for macOS/PyMOL)"""
    try:
        from .gui.qt_adapter import QtWidgets
        if QtWidgets is None:
            return None
        app = QtWidgets.QApplication.instance()
        if app is None:
            # If PyMOL hasn't initialized the Qt event loop, don't force creation unless needed
            # On some macOS environments, directly creating QApplication causes Segfault
            app = QtWidgets.QApplication([])
        return app
    except Exception as e:
        print(f"[GLINT] Failed to get QApplication: {e}")
        return None

def glint_gui():
    """Launch GLINT unified GUI window (non-modal, non-blocking event loop)"""
    global _dlg
    
    # Safely check if PyQt is available (avoid crash from importing in main process)
    print("[GLINT] Checking Qt availability...")
    if not _check_qt_safe():
        _info(
            "Qt binding (PyQt5/6/PySide2/6) not installed or unavailable, cannot start GUI"
        )
        print("\n💡 Install PyQt5 or PyQt6 to use the GUI:")
        print("   conda install -c conda-forge pyqt -y")
        print("   # or")
        print("   pip install PyQt5")
        _print_cli_fallback()
        return
    
    print("[GLINT] Qt is available, loading GUI...")
    GLINTDialog = _import_gui_dialog()
    if GLINTDialog is None:
        _info("Failed to import GUI module", "Failed to import GUI module")
        import sys, os
        pkg_dir = os.path.dirname(os.path.realpath(__file__))
        print(f"  Package dir: {pkg_dir}")
        print(f"  gui/main_window.py exists: {os.path.exists(os.path.join(pkg_dir, 'gui', 'main_window.py'))}")
        print(f"  sys.path[0:3]: {sys.path[:3]}")
        _print_cli_fallback()
        return

    try:
        # [macOS Fix] Ensure we get the existing QApplication instance to avoid Segfault from re-initialization
        app = _get_qapp()
        
        # If window already exists, activate it
        if _dlg is not None:
            try:
                _dlg.show(); _dlg.raise_(); _dlg.activateWindow()
                return
            except Exception:
                _dlg = None
        
        # Create new dialog and show non-modally
        _dlg = GLINTDialog()
        _dlg.setModal(False)
        _dlg.show()
        _dlg.raise_()
        _dlg.activateWindow()
    except Exception as e:
        _info(f"GUI start failed: {e}", f"GUI start failed: {e}")
        import traceback; traceback.print_exc()
        print("\n💡 If you see a segmentation fault:")
        print("   1. Make sure PyQt5 is properly installed in your conda environment")
        print("   2. Try: conda install -c conda-forge pyqt --force-reinstall")
        print("   3. Restart PyMOL after reinstalling PyQt5")
        _print_cli_fallback()

# Backward compatibility alias
def molstruct_gui():
    """Legacy alias for glint_gui() - for backward compatibility"""
    return glint_gui()

def _print_cli_fallback():
    _info("Use CLI instead:", "Use CLI instead:")
    print("  highlight_csv_residues csv_path='file.csv', obj='object'")
    print("  analyze_pdb_interactions obj_name='object', output_csv='output.csv'")

# ---- Plugin entry point ----
def __init_plugin__(app=None):
    """
    PyMOL plugin entry point function
    
    Workflow:
    1. Check dependencies (already done at module load time)
    2. If dependencies are OK, register all commands
    3. Always register GUI command (for displaying error messages)
    """
    # Only register all commands if dependency check passes
    if _check_deps_safe():
        _register_commands()
    
    # Always register GUI commands
    try:
        from pymol import cmd
        cmd.extend("glint_gui", glint_gui)
        cmd.extend("molstruct_gui", molstruct_gui)  # Backward compatibility alias
    except Exception as e:
        print(f"Warning: Failed to register GUI command: {e}")
    
    # Add menu item
    try:
        from pymol.plugins import addmenuitemqt
        addmenuitemqt('GLINT - Molecular Glue Analyzer', glint_gui)
    except Exception as e:
        pass  # Silently handle

    # Welcome message (deferred check)
    # Prevent double 'v' in version (e.g. __version__ = "v0.2.3" → "vv0.2.3")
    _display_ver = __version__.lstrip('v')
    if _check_deps_safe():
        print(f"\n🧬 GLINT - Molecular Glue Analyzer v{_display_ver}")
        print("┌" + "─" * 48 + "┐")
        print("│  Quick Start:                                   │")
        print("│    • glint_gui            - Launch GUI          │")
        print("│    • help(glint_gui)       - Show help         │")
        print("│    • Plugins → GLINT       - Menu access       │")
        print("└" + "─" * 48 + "┘")
    else:
        print(f"\n🧬 GLINT v{_display_ver} - ⚠️  Some dependencies may be missing")
        print("💡 Most features are available. Use 'glint_gui' to launch GUI.\n")

# Auto-register if running within PyMOL environment (e.g. via 'run' command or import)
# IMPORTANT: Only auto-register in actual PyMOL environment, not in standalone Python
try:
    import pymol
    if hasattr(pymol, 'cmd') and hasattr(pymol, 'stored'):
        # Double-check we're in actual PyMOL environment
        # Avoid re-registering if already registered (check one key command)
        if 'glint_gui' not in pymol.cmd.keyword:
            __init_plugin__()
except Exception as e:
    # Silently fail if not in PyMOL environment
    pass
