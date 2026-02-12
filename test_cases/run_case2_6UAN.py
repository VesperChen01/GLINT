# -*- coding: utf-8 -*-
"""
GLINT Test Case 2: DDB1-DCAF15-Indisulam-RBM39 (6UAN)
======================================================
非 CRBN 系统分子胶，验证 GLINT 通用性

结构说明:
  - Chain A: DDB1
  - Chain B: DCAF15 (E3 adapter)
  - Chain C: RBM39 (neo-substrate)
  - Ligand: 8BS (Indisulam/E7820 analog)

运行方式: 在 PyMOL 中执行  run test_cases/run_case2_6UAN.py
输出: test_cases/output/case2_6UAN/ 下的 CSV 和可视化图
"""

import os
import sys

_script_dir = os.path.dirname(os.path.abspath(__file__))
# 兼容 PyMOL headless (-cq) 和 GUI (run) 两种模式
if 'test_cases' not in _script_dir:
    for _candidate in sys.argv:
        if 'test_cases' in _candidate and os.path.exists(_candidate):
            _script_dir = os.path.dirname(os.path.abspath(_candidate))
            break
    else:
        _glint_root = os.environ.get('GLINT_ROOT', os.path.expanduser('~/Desktop/git/GLINT'))
        _script_dir = os.path.join(_glint_root, 'test_cases')
_project_dir = os.path.dirname(_script_dir)
if _project_dir not in sys.path:
    sys.path.insert(0, _project_dir)

from pymol import cmd

OUTPUT_DIR = os.path.join(_script_dir, "output", "case2_6UAN")
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("  GLINT Test Case 2: 6UAN (DDB1-DCAF15-Indisulam-RBM39)")
print("  ★ Non-CRBN molecular glue system")
print("=" * 60)

# ====== Step 0: 加载结构 ======
print("\n[Step 0] Fetching 6UAN...")
cmd.delete("all")
cmd.fetch("6UAN", async_=0)
cmd.remove("solvent")
cmd.remove("not alt ''+A")
cmd.alter("all", "alt=''")
print("  ✅ Structure loaded: 6UAN")

# ====== Step 1: PPI 界面分析 (DCAF15 ↔ RBM39) ======
# analyze_protein_protein_interface returns dict with keys:
#   "interface_residues", "interface_interactions", "interface_contacts",
#   "interface_strength", "bsa", "is_strong_interface"
# It also calls visualize_ppi_interface() internally when visualize=True (default).
print("\n[Step 1] PPI Interface Analysis (DCAF15-B ↔ RBM39-C)...")
ppi_result = None
try:
    from glint.ppi_analyzer import analyze_protein_protein_interface

    ppi_csv = os.path.join(OUTPUT_DIR, "ppi_interface.csv")
    ppi_result = analyze_protein_protein_interface(
        "6UAN", ["B"], ["C"], output_csv=ppi_csv
    )
    if ppi_result:
        # Fix: correct dict key is "interface_contacts", not "interactions"
        n_contacts = ppi_result.get("interface_contacts", 0)
        n_interactions = len(ppi_result.get("interface_interactions", []))
        print(f"  ✅ PPI analysis done: {n_contacts} contact pairs, {n_interactions} interactions")
        print(f"     Interface strength: {ppi_result.get('interface_strength', 'N/A')}")
        print(f"     Strong interface: {ppi_result.get('is_strong_interface', 'N/A')}")
        print(f"     CSV → {ppi_csv}")

        # Fix: BSA is already computed inside analyze_protein_protein_interface
        # calculate_interface_bsa() returns a float, not a dict
        bsa_value = ppi_result.get("bsa", None)
        if bsa_value is not None:
            print(f"  ✅ BSA = {bsa_value:.1f} Å²")
        else:
            print("  ⚠️ BSA not available in PPI result")
except Exception as e:
    print(f"  ⚠️ PPI analysis error: {e}")

# ====== Step 2: PPI 可视化（保存渲染图） ======
# Fix: visualize_ppi_interface is already called by analyze_protein_protein_interface
# (visualize=True by default). No need to call it again; just save the scene.
print("\n[Step 2] PPI Interface Visualization (saving renders)...")
try:
    cmd.set("ray_shadow", 0)
    cmd.bg_color("white")
    cmd.set("opaque_background", 1)
    png_path = os.path.join(OUTPUT_DIR, "ppi_interface_3d.png")
    cmd.png(png_path, width=1600, height=1200, ray=1, dpi=300)
    print(f"  ✅ 3D visualization saved → {png_path}")

    # 保存 PPI 可视化 PSE
    pse_ppi_path = os.path.join(OUTPUT_DIR, "ppi_interface.pse")
    cmd.save(pse_ppi_path)
    print(f"  ✅ PPI session saved → {pse_ppi_path}")
except Exception as e:
    print(f"  ⚠️ PPI visualization error: {e}")

# ====== Step 3: 蛋白-配体相互作用 ======
# Fix: correct argument order with keyword args
# Signature: analyze_protein_ligand_interactions(obj_name, ligand_resname, protein_chains, output_csv, ...)
print("\n[Step 3] Protein-Ligand Interaction Analysis...")
pli_result = None
pli_csv = os.path.join(OUTPUT_DIR, "protein_ligand_interactions.csv")
try:
    from glint.interaction_analyzer import analyze_protein_ligand_interactions

    pli_result = analyze_protein_ligand_interactions(
        "6UAN", ligand_resname="8BS", protein_chains=["B"], output_csv=pli_csv
    )
    if pli_result:
        print(f"  ✅ Protein-Ligand interactions analyzed")
        print(f"     CSV → {pli_csv}")
except Exception as e:
    print(f"  ⚠️ Protein-ligand analysis error: {e}")

# ====== Step 4: 口袋检测 ======
print("\n[Step 4] Pocket Detection...")
try:
    from glint.pocket_detector import detect_pockets
    from glint.pocket_visualizer import visualize_pockets

    pocket_csv = os.path.join(OUTPUT_DIR, "pockets.csv")
    pockets = detect_pockets("6UAN", output_csv=pocket_csv)
    if pockets:
        print(f"  ✅ Detected {len(pockets)} pocket(s)")
        print(f"     CSV → {pocket_csv}")

        # Fix: visualize_pockets(pockets, obj_name=...) — pockets is the first arg
        visualize_pockets(pockets, obj_name="6UAN")
        png_pocket = os.path.join(OUTPUT_DIR, "pockets_3d.png")
        cmd.png(png_pocket, width=1600, height=1200, ray=1, dpi=300)
        print(f"  ✅ Pocket visualization → {png_pocket}")
    else:
        print("  ℹ️ No pockets detected (all filtered out with default parameters)")
except Exception as e:
    print(f"  ⚠️ Pocket detection error: {e}")

# ====== Step 5: 口袋-PPI 界面联动 ======
print("\n[Step 5] Pocket-PPI Interface Integration...")
try:
    from glint.pocket_glue_integration import analyze_pockets_in_ppi_interface

    pocket_ppi_csv = os.path.join(OUTPUT_DIR, "pocket_ppi_interface.csv")
    pocket_ppi = analyze_pockets_in_ppi_interface(
        "6UAN", "B", "C", output_csv=pocket_ppi_csv
    )
    if pocket_ppi:
        print(f"  ✅ Interface pockets analyzed")
        print(f"     CSV → {pocket_ppi_csv}")
except Exception as e:
    print(f"  ⚠️ Pocket-PPI integration error: {e}")

# ====== Step 6: Binding Heatmap ======
# Fix: plot_binding_heatmap(input_folder, output_path=None) expects a folder of CSVs.
# It also requires pandas.
print("\n[Step 6] Binding Heatmap...")
try:
    from glint.binding_heatmap import plot_binding_heatmap

    heatmap_path = os.path.join(OUTPUT_DIR, "binding_heatmap.png")
    plot_binding_heatmap(OUTPUT_DIR, output_path=heatmap_path)
    print(f"  ✅ Binding heatmap → {heatmap_path}")
except ImportError as e:
    print(f"  ⚠️ Binding heatmap skipped (missing dependency, likely pandas): {e}")
except Exception as e:
    print(f"  ⚠️ Binding heatmap error: {e}")

# ====== Step 7: Interaction Network ======
# Fix: output → output_path, use keyword args
# Signature: generate_interaction_network_plot(interactions_result, csv_path, output_path, show_plot, obj_name, ligand_resname, ...)
print("\n[Step 7] Interaction Network Plot...")
try:
    from glint.interaction_analyzer import generate_interaction_network_plot

    network_path = os.path.join(OUTPUT_DIR, "interaction_network.png")
    generate_interaction_network_plot(
        interactions_result=pli_result,
        csv_path=pli_csv if os.path.exists(pli_csv) else None,
        obj_name="6UAN",
        ligand_resname="8BS",
        output_path=network_path,
        show_plot=False
    )
    print(f"  ✅ Interaction network → {network_path}")
except Exception as e:
    print(f"  ⚠️ Interaction network error: {e}")

# ====== Step 8: 保存 PyMOL 会话文件 ======
print("\n[Step 8] Saving PyMOL session...")
try:
    pse_path = os.path.join(OUTPUT_DIR, "case2_6UAN_final.pse")
    cmd.save(pse_path)
    print(f"  ✅ PyMOL session saved → {pse_path}")
except Exception as e:
    print(f"  ⚠️ Session save error: {e}")

# ====== 最终汇总 ======
print("\n" + "=" * 60)
print("  Case 2 Complete!")
print(f"  Output directory: {OUTPUT_DIR}")
print("  Generated files:")
for f in sorted(os.listdir(OUTPUT_DIR)):
    fpath = os.path.join(OUTPUT_DIR, f)
    size = os.path.getsize(fpath)
    print(f"    📄 {f}  ({size:,} bytes)")
print("=" * 60)
