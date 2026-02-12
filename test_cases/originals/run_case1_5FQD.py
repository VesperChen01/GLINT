# -*- coding: utf-8 -*-
"""
GLINT Test Case 1: CRBN-Lenalidomide-CK1α (5FQD)
==================================================
经典 IMiD 分子胶全流程测试

结构说明:
  - Chain A: CRBN (E3 ligase)
  - Chain C: CK1α (neo-substrate)
  - Ligand: LEN (Lenalidomide)

运行方式: 在 PyMOL 中执行  run test_cases/run_case1_5FQD.py
输出: test_cases/output/case1_5FQD/ 下的 CSV 和可视化图
"""

import os
import sys

# 确保 glint 可导入
# 兼容 PyMOL headless (-cq) 和 GUI (run) 两种模式
_script_dir = os.path.dirname(os.path.abspath(__file__))
# 如果检测到 __file__ 是相对路径且 _script_dir 不包含 test_cases，
# 则尝试从 sys.argv 或已知路径推断
if 'test_cases' not in _script_dir:
    # PyMOL -cq 模式: __file__ 可能只是文件名，需要手动寻找
    for _candidate in sys.argv:
        if 'test_cases' in _candidate and os.path.exists(_candidate):
            _script_dir = os.path.dirname(os.path.abspath(_candidate))
            break
    else:
        # 兜底: 使用固定路径或环境变量
        _glint_root = os.environ.get('GLINT_ROOT', os.path.expanduser('~/Desktop/git/GLINT'))
        _script_dir = os.path.join(_glint_root, 'test_cases')
_project_dir = os.path.dirname(_script_dir)
if _project_dir not in sys.path:
    sys.path.insert(0, _project_dir)

from pymol import cmd

# ====== 输出目录 ======
OUTPUT_DIR = os.path.join(_script_dir, "output", "case1_5FQD")
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("  GLINT Test Case 1: 5FQD (CRBN-Lenalidomide-CK1α)")
print("=" * 60)

# ====== Step 0: 加载结构 ======
print("\n[Step 0] Fetching 5FQD...")
cmd.delete("all")
cmd.fetch("5FQD", async_=0)
cmd.remove("solvent")
cmd.remove("not alt ''+A")  # 去掉 altloc
cmd.alter("all", "alt=''")
print("  ✅ Structure loaded: 5FQD")

# ====== Step 1: PPI 界面分析 ======
print("\n[Step 1] PPI Interface Analysis (CRBN A ↔ CK1α C)...")
try:
    from glint.ppi_analyzer import analyze_protein_protein_interface, calculate_interface_bsa
    
    ppi_csv = os.path.join(OUTPUT_DIR, "ppi_interface.csv")
    ppi_result = analyze_protein_protein_interface(
        "5FQD", ["A"], ["C"], output_csv=ppi_csv
    )
    if ppi_result:
        n_contacts = len(ppi_result.get("interactions", []))
        print(f"  ✅ PPI analysis done: {n_contacts} interface contacts")
        print(f"     CSV → {ppi_csv}")
    
    # BSA 计算
    bsa = calculate_interface_bsa("5FQD", "A", "C")
    if bsa:
        print(f"  ✅ BSA = {bsa.get('bsa', 'N/A'):.1f} Å²")
except Exception as e:
    print(f"  ⚠️ PPI analysis error: {e}")

# ====== Step 2: PPI 界面可视化 ======
print("\n[Step 2] PPI Interface Visualization...")
try:
    from glint.ppi_analyzer import visualize_ppi_interface
    visualize_ppi_interface("5FQD", "A", "C")
    
    # 保存 PyMOL 渲染图
    cmd.set("ray_shadow", 0)
    cmd.bg_color("white")
    cmd.set("opaque_background", 1)
    png_path = os.path.join(OUTPUT_DIR, "ppi_interface_3d.png")
    cmd.png(png_path, width=1600, height=1200, ray=1, dpi=300)
    print(f"  ✅ 3D visualization saved → {png_path}")
except Exception as e:
    print(f"  ⚠️ PPI visualization error: {e}")

# ====== Step 3: G-Motif 识别 ======
print("\n[Step 3] CRBN G-Motif (G-loop) Detection...")
try:
    from glint.g_motif_analyzer import find_crbn_g_motif
    
    gmotif_csv = os.path.join(OUTPUT_DIR, "g_motif.csv")
    gmotif_result = find_crbn_g_motif("5FQD", output_csv=gmotif_csv, template_mode="builtin")
    if gmotif_result:
        n_hits = len(gmotif_result) if isinstance(gmotif_result, list) else 1
        print(f"  ✅ G-motif found: {n_hits} hit(s)")
        print(f"     CSV → {gmotif_csv}")
except Exception as e:
    print(f"  ⚠️ G-motif detection error: {e}")

# ====== Step 4: Neo-Epitope 识别 ======
print("\n[Step 4] Neo-Epitope Identification...")
try:
    from glint.ppi_analyzer import identify_neo_epitope
    
    neo_csv = os.path.join(OUTPUT_DIR, "neo_epitope.csv")
    neo_result = identify_neo_epitope("5FQD", ["A"], ["C"], "LEN", output_csv=neo_csv)
    if neo_result:
        n_neo = len(neo_result.get("neo_contacts", []))
        print(f"  ✅ Neo-epitope contacts: {n_neo}")
        print(f"     CSV → {neo_csv}")
except Exception as e:
    print(f"  ⚠️ Neo-epitope error: {e}")

# ====== Step 5: 蛋白-配体相互作用分析 ======
print("\n[Step 5] Protein-Ligand Interaction Analysis...")
try:
    from glint.interaction_analyzer import analyze_protein_ligand_interactions
    
    pli_csv = os.path.join(OUTPUT_DIR, "protein_ligand_interactions.csv")
    pli_result = analyze_protein_ligand_interactions("5FQD", "A", "LEN", output_csv=pli_csv)
    if pli_result:
        print(f"  ✅ Protein-Ligand interactions analyzed")
        print(f"     CSV → {pli_csv}")
except Exception as e:
    print(f"  ⚠️ Protein-ligand analysis error: {e}")

# ====== Step 6: 口袋检测 ======
print("\n[Step 6] Pocket Detection...")
try:
    from glint.pocket_detector import detect_pockets
    from glint.pocket_visualizer import visualize_pockets
    
    pocket_csv = os.path.join(OUTPUT_DIR, "pockets.csv")
    pockets = detect_pockets("5FQD", output_csv=pocket_csv)
    if pockets:
        print(f"  ✅ Detected {len(pockets)} pocket(s)")
        print(f"     CSV → {pocket_csv}")
        
        # 可视化口袋
        visualize_pockets("5FQD", pockets)
        png_pocket = os.path.join(OUTPUT_DIR, "pockets_3d.png")
        cmd.png(png_pocket, width=1600, height=1200, ray=1, dpi=300)
        print(f"  ✅ Pocket visualization → {png_pocket}")
except Exception as e:
    print(f"  ⚠️ Pocket detection error: {e}")

# ====== Step 7: 口袋-PPI 界面联动 ======
print("\n[Step 7] Pocket-PPI Interface Integration...")
try:
    from glint.pocket_glue_integration import analyze_pockets_in_ppi_interface
    
    pocket_ppi_csv = os.path.join(OUTPUT_DIR, "pocket_ppi_interface.csv")
    pocket_ppi = analyze_pockets_in_ppi_interface(
        "5FQD", "A", "C", output_csv=pocket_ppi_csv
    )
    if pocket_ppi:
        print(f"  ✅ Interface pockets analyzed")
        print(f"     CSV → {pocket_ppi_csv}")
except Exception as e:
    print(f"  ⚠️ Pocket-PPI integration error: {e}")

# ====== Step 8: Binding Heatmap ======
print("\n[Step 8] Binding Heatmap...")
try:
    from glint.binding_heatmap import plot_binding_heatmap
    
    heatmap_path = os.path.join(OUTPUT_DIR, "binding_heatmap.png")
    plot_binding_heatmap("5FQD", "A", "C", output=heatmap_path)
    print(f"  ✅ Binding heatmap → {heatmap_path}")
except Exception as e:
    print(f"  ⚠️ Binding heatmap error: {e}")

# ====== Step 9: 综合分子胶设计分析 ======
print("\n[Step 9] Comprehensive Glue Design Analysis...")
try:
    from glint.glue_design_analyzer import comprehensive_glue_design_analysis
    
    glue_csv = os.path.join(OUTPUT_DIR, "glue_design_analysis.csv")
    glue_result = comprehensive_glue_design_analysis(
        "5FQD", "A", "C", "LEN", output_csv=glue_csv
    )
    if glue_result:
        print(f"  ✅ Comprehensive glue design analysis done")
        print(f"     CSV → {glue_csv}")
except Exception as e:
    print(f"  ⚠️ Glue design analysis error: {e}")

# ====== 最终汇总 ======
print("\n" + "=" * 60)
print("  Case 1 Complete!")
print(f"  Output directory: {OUTPUT_DIR}")
print("  Generated files:")
for f in sorted(os.listdir(OUTPUT_DIR)):
    fpath = os.path.join(OUTPUT_DIR, f)
    size = os.path.getsize(fpath)
    print(f"    📄 {f}  ({size:,} bytes)")
print("=" * 60)
