# -*- coding: utf-8 -*-
"""
GLINT Test Case 3: CRBN-CC885-GSPT1 (5HXB)
=============================================
IMiD 分子胶 + 2D 配体图 + G-motif + 三元复合物分析

结构说明:
  - Chain A: GSPT1 (substrate)
  - Chain B: DDB1
  - Chain C: CRBN
  - Ligand: 85C (CC-885S)

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
print("  GLINT Test Case 3: 5HXB (CRBN-CC885-GSPT1)")
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
# analyze_protein_protein_interface returns dict with keys:
#   "interface_residues", "interface_interactions", "interface_contacts",
#   "interface_strength", "bsa", "is_strong_interface"
# It calls visualize_ppi_interface() internally when visualize=True (default).
print("\n[Step 1] PPI Interface Analysis (CRBN-C ↔ GSPT1-A)...")
ppi_result = None
try:
    from glint.ppi_analyzer import analyze_protein_protein_interface

    ppi_csv = os.path.join(OUTPUT_DIR, "ppi_interface.csv")
    ppi_result = analyze_protein_protein_interface("5HXB", ["C"], ["A"], output_csv=ppi_csv)
    if ppi_result:
        # Fix: correct dict keys
        n_contacts = ppi_result.get("interface_contacts", 0)
        n_interactions = len(ppi_result.get("interface_interactions", []))
        print(f"  ✅ PPI: {n_contacts} contact pairs, {n_interactions} interactions → {ppi_csv}")
        print(f"     Interface strength: {ppi_result.get('interface_strength', 'N/A')}")

        # Fix: BSA is already in ppi_result as a float
        bsa_value = ppi_result.get("bsa", None)
        if bsa_value is not None:
            print(f"  ✅ BSA = {bsa_value:.1f} Å²")

    # 保存 PPI 渲染图
    cmd.set("ray_shadow", 0)
    cmd.bg_color("white")
    cmd.set("opaque_background", 1)
    png_ppi = os.path.join(OUTPUT_DIR, "ppi_interface_3d.png")
    cmd.png(png_ppi, width=1600, height=1200, ray=1, dpi=300)
    print(f"  ✅ PPI 3D visualization → {png_ppi}")

    pse_ppi = os.path.join(OUTPUT_DIR, "ppi_interface.pse")
    cmd.save(pse_ppi)
    print(f"  ✅ PPI session → {pse_ppi}")
except Exception as e:
    print(f"  ⚠️ PPI error: {e}")

# ====== Step 2: G-Motif (GSPT1 模板) ======
# Fix: parameter is "out_csv", not "output_csv"
print("\n[Step 2] G-Motif Detection (GSPT1 template)...")
try:
    from glint.g_motif_analyzer import find_crbn_g_motif, validate_g_motif_geometry

    gmotif_csv = os.path.join(OUTPUT_DIR, "g_motif.csv")
    gmotif = find_crbn_g_motif("5HXB", out_csv=gmotif_csv, template_mode="builtin")
    if gmotif:
        print(f"  ✅ G-motif found → {gmotif_csv}")

    # 几何验证
    geom = validate_g_motif_geometry("5HXB", "A", "570-577")
    if geom:
        print(f"  ✅ G-motif geometry validated")
except Exception as e:
    print(f"  ⚠️ G-motif error: {e}")

# ====== Step 3: CRBN 氢键验证 ======
# Fix: validate_crbn_hbonds(obj_name, g_motif_chain, g_motif_resi_range, crbn_chain, ...)
# It does NOT have an output_csv parameter
print("\n[Step 3] CRBN H-bond Validation...")
try:
    from glint.g_motif_analyzer import validate_crbn_hbonds

    hbonds = validate_crbn_hbonds("5HXB", "A", "570-577", "C")
    if hbonds:
        print(f"  ✅ CRBN H-bonds validated")
        # Save result manually as CSV if needed
        hbond_csv = os.path.join(OUTPUT_DIR, "crbn_hbonds.csv")
        if isinstance(hbonds, dict) and hbonds.get("hbonds"):
            import csv
            with open(hbond_csv, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(["donor_chain", "donor_resi", "donor_resn", "donor_atom",
                                 "acceptor_chain", "acceptor_resi", "acceptor_resn", "acceptor_atom",
                                 "distance"])
                for hb in hbonds.get("hbonds", []):
                    writer.writerow([
                        hb.get("donor_chain", ""), hb.get("donor_resi", ""),
                        hb.get("donor_resn", ""), hb.get("donor_atom", ""),
                        hb.get("acceptor_chain", ""), hb.get("acceptor_resi", ""),
                        hb.get("acceptor_resn", ""), hb.get("acceptor_atom", ""),
                        hb.get("distance", "")
                    ])
            print(f"     CSV → {hbond_csv}")
except Exception as e:
    print(f"  ⚠️ H-bond validation error: {e}")

# ====== Step 4: 蛋白-配体相互作用 ======
# Fix: correct argument order with keyword args, ligand is "85C" not "LEN"
# Signature: analyze_protein_ligand_interactions(obj_name, ligand_resname, protein_chains, output_csv, ...)
print("\n[Step 4] Protein-Ligand Interactions...")
pli_result = None
pli_csv = os.path.join(OUTPUT_DIR, "protein_ligand.csv")
try:
    from glint.interaction_analyzer import analyze_protein_ligand_interactions

    pli_result = analyze_protein_ligand_interactions(
        "5HXB", ligand_resname="85C", protein_chains=["C"], output_csv=pli_csv
    )
    if pli_result:
        print(f"  ✅ Protein-ligand → {pli_csv}")
except Exception as e:
    print(f"  ⚠️ Protein-ligand error: {e}")

# ====== Step 5: 3D 蛋白-配体可视化 ======
# Fix: visualize_protein_ligand_3d(obj_name, interactions_result=None, ligand_resname=None, ...)
print("\n[Step 5] Protein-Ligand 3D Visualization...")
try:
    from glint.interaction_analyzer import visualize_protein_ligand_3d

    visualize_protein_ligand_3d("5HXB", interactions_result=pli_result, ligand_resname="85C")
    cmd.set("ray_shadow", 0)
    cmd.bg_color("white")
    png_3d = os.path.join(OUTPUT_DIR, "protein_ligand_3d.png")
    cmd.png(png_3d, width=1600, height=1200, ray=1, dpi=300)
    print(f"  ✅ 3D visualization → {png_3d}")
except Exception as e:
    print(f"  ⚠️ Protein-ligand 3D error: {e}")

# ====== Step 6: 2D 相互作用图 ======
# Fix: generate_2d_interaction_diagram(csv_path, ligand_resname, ..., output_path=None)
# First arg should be CSV path from PLI analysis, not the object name
print("\n[Step 6] 2D Interaction Diagram...")
try:
    from glint.interaction_2d_plot import generate_2d_interaction_diagram

    diagram_path = os.path.join(OUTPUT_DIR, "interaction_2d.png")
    if os.path.exists(pli_csv):
        generate_2d_interaction_diagram(
            pli_csv, "85C", obj_name="5HXB", output_path=diagram_path
        )
        print(f"  ✅ 2D Diagram → {diagram_path}")
    else:
        print("  ⚠️ 2D diagram skipped: PLI CSV not available")
except Exception as e:
    print(f"  ⚠️ 2D diagram error: {e}")

# ====== Step 7: 三元复合物分析 ======
# Fix: analyze_ternary_complex(obj_name, ligand_resname, protein1_chains, protein2_chains, output_csv, ...)
print("\n[Step 7] Ternary Complex Analysis...")
try:
    from glint.interaction_analyzer import analyze_ternary_complex

    ternary_csv = os.path.join(OUTPUT_DIR, "ternary_complex.csv")
    ternary = analyze_ternary_complex(
        "5HXB", ligand_resname="85C",
        protein1_chains=["C"], protein2_chains=["A"],
        output_csv=ternary_csv
    )
    if ternary:
        print(f"  ✅ Ternary complex → {ternary_csv}")
except Exception as e:
    print(f"  ⚠️ Ternary complex error: {e}")

# ====== Step 8: Neo-Epitope ======
# Fix: ligand is "85C" not "LEN"
# Fix: visualize_neo_epitope(obj_name, neo_result, ...) - second arg is result dict
print("\n[Step 8] Neo-Epitope Identification...")
try:
    from glint.ppi_analyzer import identify_neo_epitope, visualize_neo_epitope

    neo_csv = os.path.join(OUTPUT_DIR, "neo_epitope.csv")
    neo = identify_neo_epitope("5HXB", ["C"], ["A"], "85C", output_csv=neo_csv)
    if neo:
        n_neo = len(neo.get("neo_contacts", []))
        print(f"  ✅ Neo-epitope: {n_neo} contacts → {neo_csv}")

        # Fix: pass neo_result dict, not chain strings
        visualize_neo_epitope("5HXB", neo)
        png_neo = os.path.join(OUTPUT_DIR, "neo_epitope_3d.png")
        cmd.png(png_neo, width=1600, height=1200, ray=1, dpi=300)
        print(f"  ✅ Neo-epitope visualization → {png_neo}")
    else:
        print("  ⚠️ No neo-epitope contacts found")
except Exception as e:
    print(f"  ⚠️ Neo-epitope error: {e}")

# ====== Step 9: 保存 PyMOL 会话文件 ======
print("\n[Step 9] Saving PyMOL session...")
try:
    pse_path = os.path.join(OUTPUT_DIR, "case3_5HXB_final.pse")
    cmd.save(pse_path)
    print(f"  ✅ PyMOL session saved → {pse_path}")
except Exception as e:
    print(f"  ⚠️ Session save error: {e}")

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
