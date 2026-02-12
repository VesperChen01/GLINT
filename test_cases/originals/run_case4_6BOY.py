# -*- coding: utf-8 -*-
"""
GLINT Test Case 4: CRBN-CC-885-GSPT1 (6BOY)
=============================================
突变扫描 + 蛋白表面分析
Mutation scanning + protein surface analysis

结构说明 / Structure Info:
  - Chain B: CRBN (E3 ligase)
  - Chain D: GSPT1 (neo-substrate)
  - Ligand: 6GY (CC-885, molecular glue)

分析重点 / Focus:
  - G-loop 关键残基突变效应分析
  - 蛋白表面静电势渲染
  - 综合分子胶设计评估

运行方式: 在 PyMOL 中执行  run test_cases/run_case4_6BOY.py
输出: test_cases/output/case4_6BOY/ 下的 CSV 和可视化图
"""

import os
import sys

# 确保 glint 可导入 / Ensure glint is importable
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

# ====== 输出目录 / Output directory ======
OUTPUT_DIR = os.path.join(_script_dir, "output", "case4_6BOY")
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("  GLINT Test Case 4: 6BOY (CRBN-CC-885-GSPT1)")
print("  ★ Mutation Scanning + Surface Analysis")
print("=" * 60)

# ====== Step 0: 加载结构 / Load structure ======
print("\n[Step 0] Fetching 6BOY...")
cmd.delete("all")
cmd.fetch("6BOY", async_=0)
cmd.remove("solvent")
cmd.remove("not alt ''+A")  # 去掉 altloc / Remove alternate conformations
cmd.alter("all", "alt=''")
print("  ✅ Structure loaded: 6BOY")

# ====== Step 1: PPI 界面分析 / PPI Interface Analysis ======
print("\n[Step 1] PPI Interface Analysis (CRBN-B ↔ GSPT1-D)...")
try:
    from glint.ppi_analyzer import analyze_protein_protein_interface, calculate_interface_bsa

    ppi_csv = os.path.join(OUTPUT_DIR, "ppi_interface.csv")
    ppi_result = analyze_protein_protein_interface(
        "6BOY", ["B"], ["D"], output_csv=ppi_csv
    )
    if ppi_result:
        n_contacts = len(ppi_result.get("interactions", []))
        print(f"  ✅ PPI analysis done: {n_contacts} interface contacts")
        print(f"     CSV → {ppi_csv}")

    # BSA 计算 / BSA calculation
    bsa = calculate_interface_bsa("6BOY", "B", "D")
    if bsa:
        print(f"  ✅ BSA = {bsa.get('bsa', 'N/A'):.1f} Å²")
except Exception as e:
    print(f"  ⚠️ PPI analysis error: {e}")

# ====== Step 2: 突变分析 / Mutation Analysis ======
# 对 GSPT1 (Chain D) G-loop 关键残基做丙氨酸突变扫描
# Alanine scanning on key G-loop residues of GSPT1
print("\n[Step 2] Mutation Analysis (G-loop key residues)...")
try:
    from glint.mutation_analyzer import perform_mutation, analyze_mutation_effects

    # G-loop 关键残基: D571, D574 → ALA
    mutations = [("D", "571", "ALA"), ("D", "574", "ALA")]
    mutation_csv = os.path.join(OUTPUT_DIR, "mutation_analysis.csv")
    mut_result = analyze_mutation_effects(
        "6BOY", mutations,
        partner_sel="chain B",
        output_csv=mutation_csv,
        method="pymol"
    )
    if mut_result:
        print(f"  ✅ Mutation analysis done: {len(mutations)} mutations scanned")
        print(f"     CSV → {mutation_csv}")
except Exception as e:
    print(f"  ⚠️ Mutation analysis error: {e}")

# ====== Step 3: 蛋白表面分析 / Protein Surface Analysis ======
# 渲染蛋白表面静电势 / Render protein surface electrostatic potential
print("\n[Step 3] Protein Surface Rendering...")
try:
    from glint.protein_surface_analyzer import render_protein_surface

    render_protein_surface("6BOY", use_apbs=False, color_by='potential')

    # 保存 3D 截图 / Save 3D screenshot
    cmd.set("ray_shadow", 0)
    cmd.bg_color("white")
    cmd.set("opaque_background", 1)
    png_surface = os.path.join(OUTPUT_DIR, "protein_surface_3d.png")
    cmd.png(png_surface, width=1600, height=1200, ray=1, dpi=300)
    print(f"  ✅ Surface rendering saved → {png_surface}")
except Exception as e:
    print(f"  ⚠️ Surface rendering error: {e}")

# ====== Step 4: 蛋白-配体相互作用 / Protein-Ligand Interactions ======
print("\n[Step 4] Protein-Ligand Interaction Analysis...")
try:
    from glint.interaction_analyzer import analyze_protein_ligand_interactions

    pli_csv = os.path.join(OUTPUT_DIR, "protein_ligand_interactions.csv")
    pli_result = analyze_protein_ligand_interactions("6BOY", "6GY", output_csv=pli_csv)
    if pli_result:
        print(f"  ✅ Protein-Ligand interactions analyzed")
        print(f"     CSV → {pli_csv}")
except Exception as e:
    print(f"  ⚠️ Protein-ligand analysis error: {e}")

# ====== Step 5: G-Motif 检测 / G-Motif Detection ======
print("\n[Step 5] CRBN G-Motif (G-loop) Detection...")
try:
    from glint.g_motif_analyzer import find_crbn_g_motif

    gmotif_csv = os.path.join(OUTPUT_DIR, "g_motif.csv")
    gmotif_result = find_crbn_g_motif("6BOY", output_csv=gmotif_csv, template_mode="builtin")
    if gmotif_result:
        n_hits = len(gmotif_result) if isinstance(gmotif_result, list) else 1
        print(f"  ✅ G-motif found: {n_hits} hit(s)")
        print(f"     CSV → {gmotif_csv}")
except Exception as e:
    print(f"  ⚠️ G-motif detection error: {e}")

# ====== Step 6: Binding Heatmap ======
print("\n[Step 6] Binding Heatmap...")
try:
    from glint.binding_heatmap import plot_binding_heatmap

    heatmap_path = os.path.join(OUTPUT_DIR, "binding_heatmap.png")
    plot_binding_heatmap("6BOY", "B", "D", output=heatmap_path)
    print(f"  ✅ Binding heatmap → {heatmap_path}")
except Exception as e:
    print(f"  ⚠️ Binding heatmap error: {e}")

# ====== Step 7: 综合分子胶分析 / Comprehensive Glue Design Analysis ======
print("\n[Step 7] Comprehensive Glue Design Analysis...")
try:
    from glint.glue_design_analyzer import comprehensive_glue_design_analysis

    glue_csv = os.path.join(OUTPUT_DIR, "glue_design_analysis.csv")
    glue_result = comprehensive_glue_design_analysis(
        "6BOY", "B", "D", "6GY", output_csv=glue_csv
    )
    if glue_result:
        print(f"  ✅ Comprehensive glue design analysis done")
        print(f"     CSV → {glue_csv}")
except Exception as e:
    print(f"  ⚠️ Glue design analysis error: {e}")

# ====== 最终汇总 / Final Summary ======
print("\n" + "=" * 60)
print("  Case 4 Complete!")
print(f"  Output directory: {OUTPUT_DIR}")
print("  Generated files:")
for f in sorted(os.listdir(OUTPUT_DIR)):
    fpath = os.path.join(OUTPUT_DIR, f)
    size = os.path.getsize(fpath)
    print(f"    📄 {f}  ({size:,} bytes)")
print("=" * 60)
