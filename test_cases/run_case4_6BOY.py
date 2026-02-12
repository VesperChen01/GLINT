# -*- coding: utf-8 -*-
"""
GLINT Test Case 4: CRBN-CC-885-GSPT1 (6BOY)
=============================================
突变扫描 + 蛋白表面分析
Mutation scanning + protein surface analysis

结构说明 / Structure Info:
  - Chain A: DDB1
  - Chain B: CRBN (E3 ligase, with ligand RN6)
  - Chain C: GSPT1 (neo-substrate)
  - Ligand: RN6 (CC-885, molecular glue)

分析重点 / Focus:
  - G-loop 关键残基突变效应分析
  - 蛋白表面静电势渲染
  - PPI 界面分析

运行方式: 在 PyMOL 中执行  run test_cases/run_case4_6BOY.py
输出: test_cases/output/case4_6BOY/ 下的 CSV 和可视化图
"""

import os
import sys

_script_dir = os.path.dirname(os.path.abspath(__file__))
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

OUTPUT_DIR = os.path.join(_script_dir, "output", "case4_6BOY")
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("  GLINT Test Case 4: 6BOY (CRBN-CC-885-GSPT1)")
print("  ★ Mutation Scanning + Surface Analysis")
print("=" * 60)

# ====== Step 0: 加载结构 ======
print("\n[Step 0] Fetching 6BOY...")
cmd.delete("all")
cmd.fetch("6BOY", async_=0)
cmd.remove("solvent")
cmd.remove("not alt ''+A")
cmd.alter("all", "alt=''")
print("  ✅ Structure loaded: 6BOY")
print("  Chains present:", cmd.get_chains("6BOY"))

# ====== Step 1: PPI 界面分析 ======
print("\n[Step 1] PPI Interface Analysis (CRBN-B ↔ GSPT1-C)...")
ppi_result = None
try:
    from glint.ppi_analyzer import analyze_protein_protein_interface

    ppi_csv = os.path.join(OUTPUT_DIR, "ppi_interface.csv")
    ppi_result = analyze_protein_protein_interface(
        "6BOY", ["B"], ["C"], output_csv=ppi_csv
    )
    if ppi_result:
        n_contacts = ppi_result.get("interface_contacts", 0)
        n_interactions = len(ppi_result.get("interface_interactions", []))
        print(f"  ✅ PPI: {n_contacts} contact pairs, {n_interactions} interactions")
        print(f"     CSV → {ppi_csv}")
        bsa_value = ppi_result.get("bsa", None)
        if bsa_value is not None:
            print(f"  ✅ BSA = {bsa_value:.1f} Å²")

    cmd.set("ray_shadow", 0)
    cmd.bg_color("white")
    cmd.set("opaque_background", 1)
    png_ppi = os.path.join(OUTPUT_DIR, "ppi_interface_3d.png")
    cmd.png(png_ppi, width=1600, height=1200, ray=1, dpi=300)
    print(f"  ✅ PPI 3D → {png_ppi}")

    pse_ppi = os.path.join(OUTPUT_DIR, "ppi_interface.pse")
    cmd.save(pse_ppi)
    print(f"  ✅ PPI session → {pse_ppi}")
except Exception as e:
    print(f"  ⚠️ PPI analysis error: {e}")

# ====== Step 2: 突变分析 ======
print("\n[Step 2] Mutation Analysis (G-loop key residues)...")
try:
    from glint.mutation_analyzer import analyze_mutation_effects

    mutations = [("C", "571", "ALA"), ("C", "574", "ALA")]
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

# ====== Step 3: 蛋白表面分析 ======
print("\n[Step 3] Protein Surface Rendering...")
try:
    from glint.protein_surface_analyzer import render_protein_surface

    render_protein_surface("6BOY", use_apbs=False, color_by='potential')

    cmd.set("ray_shadow", 0)
    cmd.bg_color("white")
    cmd.set("opaque_background", 1)
    png_surface = os.path.join(OUTPUT_DIR, "protein_surface_3d.png")
    cmd.png(png_surface, width=1600, height=1200, ray=1, dpi=300)
    print(f"  ✅ Surface rendering saved → {png_surface}")
except Exception as e:
    print(f"  ⚠️ Surface rendering error: {e}")

# ====== Step 4: 蛋白-配体相互作用 ======
print("\n[Step 4] Protein-Ligand Interaction Analysis...")
try:
    from glint.interaction_analyzer import analyze_protein_ligand_interactions

    pli_csv = os.path.join(OUTPUT_DIR, "protein_ligand_interactions.csv")
    pli_result = analyze_protein_ligand_interactions(
        "6BOY", ligand_resname="RN6", protein_chains=["B"], output_csv=pli_csv
    )
    if pli_result:
        print(f"  ✅ Protein-Ligand interactions analyzed")
        print(f"     CSV → {pli_csv}")
except Exception as e:
    print(f"  ⚠️ Protein-ligand analysis error: {e}")

# ====== Step 5: G-Motif 检测 ======
print("\n[Step 5] CRBN G-Motif (G-loop) Detection...")
try:
    from glint.g_motif_analyzer import find_crbn_g_motif

    gmotif_csv = os.path.join(OUTPUT_DIR, "g_motif.csv")
    gmotif_result = find_crbn_g_motif("6BOY", out_csv=gmotif_csv, template_mode="builtin")
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
    plot_binding_heatmap(OUTPUT_DIR, output_path=heatmap_path)
    print(f"  ✅ Binding heatmap → {heatmap_path}")
except ImportError as e:
    print(f"  ⚠️ Binding heatmap skipped (missing dependency): {e}")
except Exception as e:
    print(f"  ⚠️ Binding heatmap error: {e}")

# ====== Step 7: 综合分子胶分析 ======
print("\n[Step 7] Comprehensive Glue Design Analysis...")
print("  ℹ️  Skipped: requires both template and target structures for G-loop alignment.")

# ====== Step 8: 保存 PyMOL 会话文件 ======
print("\n[Step 8] Saving PyMOL session...")
try:
    pse_path = os.path.join(OUTPUT_DIR, "case4_6BOY_final.pse")
    cmd.save(pse_path)
    print(f"  ✅ PyMOL session saved → {pse_path}")
except Exception as e:
    print(f"  ⚠️ Session save error: {e}")

# ====== 最终汇总 ======
print("\n" + "=" * 60)
print("  Case 4 Complete!")
print(f"  Output directory: {OUTPUT_DIR}")
print("  Generated files:")
for f in sorted(os.listdir(OUTPUT_DIR)):
    fpath = os.path.join(OUTPUT_DIR, f)
    size = os.path.getsize(fpath)
    print(f"    📄 {f}  ({size:,} bytes)")
print("=" * 60)
