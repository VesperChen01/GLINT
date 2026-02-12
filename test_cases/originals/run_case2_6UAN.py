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
print("\n[Step 1] PPI Interface Analysis (DCAF15-B ↔ RBM39-C)...")
try:
    from glint.ppi_analyzer import analyze_protein_protein_interface, calculate_interface_bsa
    
    ppi_csv = os.path.join(OUTPUT_DIR, "ppi_interface.csv")
    ppi_result = analyze_protein_protein_interface(
        "6UAN", ["B"], ["C"], output_csv=ppi_csv
    )
    if ppi_result:
        n_contacts = len(ppi_result.get("interactions", []))
        print(f"  ✅ PPI analysis done: {n_contacts} interface contacts")
        print(f"     CSV → {ppi_csv}")
    
    bsa = calculate_interface_bsa("6UAN", "B", "C")
    if bsa:
        print(f"  ✅ BSA = {bsa.get('bsa', 'N/A'):.1f} Å²")
except Exception as e:
    print(f"  ⚠️ PPI analysis error: {e}")

# ====== Step 2: PPI 可视化 ======
print("\n[Step 2] PPI Interface Visualization...")
try:
    from glint.ppi_analyzer import visualize_ppi_interface
    visualize_ppi_interface("6UAN", "B", "C")
    
    cmd.set("ray_shadow", 0)
    cmd.bg_color("white")
    cmd.set("opaque_background", 1)
    png_path = os.path.join(OUTPUT_DIR, "ppi_interface_3d.png")
    cmd.png(png_path, width=1600, height=1200, ray=1, dpi=300)
    print(f"  ✅ 3D visualization saved → {png_path}")
except Exception as e:
    print(f"  ⚠️ PPI visualization error: {e}")

# ====== Step 3: 蛋白-配体相互作用 ======
print("\n[Step 3] Protein-Ligand Interaction Analysis...")
try:
    from glint.interaction_analyzer import analyze_protein_ligand_interactions
    
    pli_csv = os.path.join(OUTPUT_DIR, "protein_ligand_interactions.csv")
    pli_result = analyze_protein_ligand_interactions("6UAN", "B", "8BS", output_csv=pli_csv)
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
        
        visualize_pockets("6UAN", pockets)
        png_pocket = os.path.join(OUTPUT_DIR, "pockets_3d.png")
        cmd.png(png_pocket, width=1600, height=1200, ray=1, dpi=300)
        print(f"  ✅ Pocket visualization → {png_pocket}")
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
print("\n[Step 6] Binding Heatmap...")
try:
    from glint.binding_heatmap import plot_binding_heatmap
    
    heatmap_path = os.path.join(OUTPUT_DIR, "binding_heatmap.png")
    plot_binding_heatmap("6UAN", "B", "C", output=heatmap_path)
    print(f"  ✅ Binding heatmap → {heatmap_path}")
except Exception as e:
    print(f"  ⚠️ Binding heatmap error: {e}")

# ====== Step 7: Interaction Network ======
print("\n[Step 7] Interaction Network Plot...")
try:
    from glint.interaction_analyzer import generate_interaction_network_plot
    
    network_path = os.path.join(OUTPUT_DIR, "interaction_network.png")
    generate_interaction_network_plot("6UAN", "B", "8BS", output=network_path)
    print(f"  ✅ Interaction network → {network_path}")
except Exception as e:
    print(f"  ⚠️ Interaction network error: {e}")

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
