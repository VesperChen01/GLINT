# -*- coding: utf-8 -*-
"""
GLINT Test Case 3: CRBN-Pomalidomide-IKZF1/GSPT1 (5HXB)
=========================================================
IMiD 分子胶 + 2D 配体图 + G-motif 对比

结构说明:
  - Chain A: GSPT1 (substrate)
  - Chain B: DDB1
  - Chain C: CRBN
  - Ligand: LEN (pomalidomide)

运行方式: 在 PyMOL 中执行  run test_cases/run_case3_5HXB.py
输出: test_cases/output/case3_5HXB/ 下的 CSV 和可视化图
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

OUTPUT_DIR = os.path.join(_script_dir, "output", "case3_5HXB")
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("=" * 60)
print("  GLINT Test Case 3: 5HXB (CRBN-Pomalidomide-GSPT1)")
print("  ★ 2D Diagram + G-motif + Ternary Analysis")
print("=" * 60)

# ====== Step 0: 加载 ======
print("\n[Step 0] Fetching 5HXB...")
cmd.delete("all")
cmd.fetch("5HXB", async_=0)
cmd.remove("solvent")
cmd.remove("not alt ''+A")
cmd.alter("all", "alt=''")
print("  ✅ Structure loaded: 5HXB")

# ====== Step 1: PPI 界面 (CRBN C ↔ GSPT1 A) ======
print("\n[Step 1] PPI Interface Analysis (CRBN-C ↔ GSPT1-A)...")
try:
    from glint.ppi_analyzer import analyze_protein_protein_interface, calculate_interface_bsa
    
    ppi_csv = os.path.join(OUTPUT_DIR, "ppi_interface.csv")
    ppi_result = analyze_protein_protein_interface("5HXB", ["C"], ["A"], output_csv=ppi_csv)
    if ppi_result:
        n = len(ppi_result.get("interactions", []))
        print(f"  ✅ PPI: {n} contacts → {ppi_csv}")
    
    bsa = calculate_interface_bsa("5HXB", "C", "A")
    if bsa:
        print(f"  ✅ BSA = {bsa.get('bsa', 'N/A'):.1f} Å²")
except Exception as e:
    print(f"  ⚠️ PPI error: {e}")

# ====== Step 2: G-Motif (GSPT1 模板) ======
print("\n[Step 2] G-Motif Detection (GSPT1 template)...")
try:
    from glint.g_motif_analyzer import find_crbn_g_motif, validate_g_motif_geometry
    
    gmotif_csv = os.path.join(OUTPUT_DIR, "g_motif.csv")
    gmotif = find_crbn_g_motif("5HXB", output_csv=gmotif_csv, template_mode="builtin")
    if gmotif:
        print(f"  ✅ G-motif found → {gmotif_csv}")
    
    # 几何验证
    geom = validate_g_motif_geometry("5HXB", "A", "570-577")
    if geom:
        print(f"  ✅ G-motif geometry validated")
except Exception as e:
    print(f"  ⚠️ G-motif error: {e}")

# ====== Step 3: CRBN 氢键验证 ======
print("\n[Step 3] CRBN H-bond Validation...")
try:
    from glint.g_motif_analyzer import validate_crbn_hbonds
    
    hbond_csv = os.path.join(OUTPUT_DIR, "crbn_hbonds.csv")
    hbonds = validate_crbn_hbonds("5HXB", output_csv=hbond_csv)
    if hbonds:
        print(f"  ✅ CRBN H-bonds validated → {hbond_csv}")
except Exception as e:
    print(f"  ⚠️ H-bond validation error: {e}")

# ====== Step 4: 2D 相互作用图 ======
print("\n[Step 4] 2D Interaction Diagram...")
try:
    from glint.interaction_2d_plot import generate_2d_interaction_diagram
    
    diagram_path = os.path.join(OUTPUT_DIR, "interaction_2d.png")
    generate_2d_interaction_diagram("5HXB", "C", "LEN", output=diagram_path)
    print(f"  ✅ 2D Diagram → {diagram_path}")
except Exception as e:
    print(f"  ⚠️ 2D diagram error: {e}")

# ====== Step 5: 蛋白-配体相互作用 ======
print("\n[Step 5] Protein-Ligand Interactions...")
try:
    from glint.interaction_analyzer import analyze_protein_ligand_interactions, visualize_protein_ligand_3d
    
    pli_csv = os.path.join(OUTPUT_DIR, "protein_ligand.csv")
    pli = analyze_protein_ligand_interactions("5HXB", "C", "LEN", output_csv=pli_csv)
    if pli:
        print(f"  ✅ Protein-ligand → {pli_csv}")
    
    # 3D 可视化
    visualize_protein_ligand_3d("5HXB", "C", "LEN")
    cmd.set("ray_shadow", 0)
    cmd.bg_color("white")
    png_3d = os.path.join(OUTPUT_DIR, "protein_ligand_3d.png")
    cmd.png(png_3d, width=1600, height=1200, ray=1, dpi=300)
    print(f"  ✅ 3D visualization → {png_3d}")
except Exception as e:
    print(f"  ⚠️ Protein-ligand error: {e}")

# ====== Step 6: 三元复合物分析 ======
print("\n[Step 6] Ternary Complex Analysis...")
try:
    from glint.interaction_analyzer import analyze_ternary_complex
    
    ternary_csv = os.path.join(OUTPUT_DIR, "ternary_complex.csv")
    ternary = analyze_ternary_complex("5HXB", "C", "A", "LEN", output_csv=ternary_csv)
    if ternary:
        print(f"  ✅ Ternary complex → {ternary_csv}")
except Exception as e:
    print(f"  ⚠️ Ternary complex error: {e}")

# ====== Step 7: Neo-Epitope ======
print("\n[Step 7] Neo-Epitope Identification...")
try:
    from glint.ppi_analyzer import identify_neo_epitope, visualize_neo_epitope
    
    neo_csv = os.path.join(OUTPUT_DIR, "neo_epitope.csv")
    neo = identify_neo_epitope("5HXB", ["C"], ["A"], "LEN", output_csv=neo_csv)
    if neo:
        n_neo = len(neo.get("neo_contacts", []))
        print(f"  ✅ Neo-epitope: {n_neo} contacts → {neo_csv}")
    
    visualize_neo_epitope("5HXB", ["C"], ["A"], "LEN")
    png_neo = os.path.join(OUTPUT_DIR, "neo_epitope_3d.png")
    cmd.png(png_neo, width=1600, height=1200, ray=1, dpi=300)
    print(f"  ✅ Neo-epitope visualization → {png_neo}")
except Exception as e:
    print(f"  ⚠️ Neo-epitope error: {e}")

# ====== 汇总 ======
print("\n" + "=" * 60)
print("  Case 3 Complete!")
print(f"  Output directory: {OUTPUT_DIR}")
print("  Generated files:")
for f in sorted(os.listdir(OUTPUT_DIR)):
    fpath = os.path.join(OUTPUT_DIR, f)
    size = os.path.getsize(fpath)
    print(f"    📄 {f}  ({size:,} bytes)")
print("=" * 60)
