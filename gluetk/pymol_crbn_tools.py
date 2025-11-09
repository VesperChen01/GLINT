# -*- coding: utf-8 -*-
"""
PyMOL plugin: CRBN molecular glue toolkit (interface maps, ΔΔG heatmaps, score coloring, degron annotation)

Commands
- interface_map CRBN_sel, POI_sel [, name=iface, cutoff=4.0]
    Build interface contact visuals and per-residue contact counts on POI (stored in b-factor).

- ddg_heatmap CRBN_sel, POI_sel [, method=auto|foldx|asa, name=ddg, scale=0.025]
    Color POI by predicted ΔΔG of binding. If FoldX is detected, uses alanine scan; else ASA-based proxy.

- score_color POI_sel, scores_csv [, col_chain=chain, col_resi=resi, col_score=score, min=None, max=None]
    Read CSV with chain,resi,score and color POI by score; stores score into b-factor.

- degron_annotate POI_sel [, csv=None, auto=true, name=degron]
    Annotate degron-like regions. If csv provided, expects columns: chain,start,end,label (1-based resi).
    Auto mode marks simple C2H2-ZF-like motifs and β-hairpin heuristics.

Installation
- Place this file on PYTHONPATH or "run pymol_crbn_tools.py" inside PyMOL, or add to Plugin Manager.

Notes
- Requires only PyMOL. FoldX (optional) used if available in PATH or via env FOLDX.
- Colors and objects created:
    iface_hbond, iface_salt, iface_hydroph, iface_resi_map, ddg_heatmap

"""
from __future__ import annotations
import os
import csv
import math
import json
import tempfile
import subprocess
from collections import defaultdict
from typing import Dict, Tuple, List, Optional

try:
    from pymol import cmd, stored
except Exception as e:  # pragma: no cover
    raise RuntimeError("This module must be loaded inside PyMOL (pymol.cmd unavailable)") from e

# ------------- Utilities -------------

def _which(exe: str) -> Optional[str]:
    path = os.environ.get("FOLDX") if exe.lower() == "foldx" else None
    if path and os.path.isfile(path) and os.access(path, os.X_OK):
        return path
    for p in os.environ.get("PATH", "").split(os.pathsep):
        cand = os.path.join(p, exe)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
        # macOS app bundle
        if exe.lower() == "foldx" and os.path.isdir(os.path.join(p, "FoldX.app")):
            app_bin = os.path.join(p, "FoldX.app", "Contents", "MacOS", "FoldX")
            if os.path.isfile(app_bin) and os.access(app_bin, os.X_OK):
                return app_bin
    return None


def _normalize(v: List[float]) -> Tuple[float, float]:
    if not v:
        return (0.0, 1.0)
    vmin, vmax = min(v), max(v)
    if math.isclose(vmin, vmax):
        vmax = vmin + 1.0
    return vmin, vmax


def _set_b_factors(sel: str, resi_to_value: Dict[Tuple[str, str, str], float]):
    # resi key: (chain, resi, icode)
    def _alter(atom):
        key = (atom.chain, atom.resi, atom.q if hasattr(atom, 'q') else atom.icode)
        if key in resi_to_value:
            atom.b = float(resi_to_value[key])
        return atom
    cmd.alter(sel, "b=b", space={"b": 0.0})  # reset
    cmd.iterate_state(1, sel, "b=b", space={})  # no-op to force load
    cmd.alter(sel, "b=b", space={})  # no-op
    cmd.iterate(sel, "b=b", space={})  # no-op
    # The efficient way is using iterating with a callback; emulate with Python loop over model
    model = cmd.get_model(sel)
    for a in model.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        if key in resi_to_value:
            icode_part = getattr(a, 'q', a.icode)
            icode_str = icode_part if icode_part else '""'
            cmd.alter(f"{sel} and chain {a.chain} and resi {a.resi} and icode {icode_str}", f"b={float(resi_to_value[key])}")
    cmd.rebuild()


def _color_by_b(sel: str, palette: str = "blue_white_red", min_val: Optional[float] = None, max_val: Optional[float] = None, ramp_name: str = ""):
    # Use spectrum b
    if min_val is None or max_val is None:
        vals = []
        model = cmd.get_model(sel)
        for a in model.atom:
            vals.append(a.b)
        vmin, vmax = _normalize(vals)
    else:
        vmin, vmax = min_val, max_val
    cmd.spectrum("b", palette, selection=sel, minimum=vmin, maximum=vmax)
    if ramp_name:
        cmd.ramp_new(ramp_name, sel, [vmin, (vmin+vmax)/2.0, vmax], ["blue", "white", "red"])


# ------------- Interface Map -------------

def interface_map(CRBN_sel: str, POI_sel: str, name: str = "iface", cutoff: float = 4.0):
    """Build interface contacts and per-residue contact counts on POI selection.
    Produces distance objects and colors residues by contact count (stored in b-factor).
    """
    crbn = f"({CRBN_sel}) and not elem H"
    poi = f"({POI_sel}) and not elem H"

    # Distance objects for interaction classes
    # Hydrogen bonds (approx): N-O / O-N within 3.5 Å
    cmd.delete(f"{name}_hbond")
    cmd.distance(f"{name}_hbond", f"({crbn}) and (elem N+O)", f"({poi}) and (elem N+O)", cutoff=3.5)
    cmd.set("dash_color", "cyan", f"{name}_hbond")
    cmd.set("dash_width", 2.0, f"{name}_hbond")

    # Salt bridges: acidic O vs basic N within 4.0 Å
    acidic = "resn ASP+GLU and elem O"
    basic = "resn LYS+ARG+HIS and elem N"
    cmd.delete(f"{name}_salt")
    cmd.distance(f"{name}_salt", f"({crbn}) and ({acidic})", f"({poi}) and ({basic})", cutoff=cutoff)
    cmd.distance(f"{name}_salt", f"({crbn}) and ({basic})", f"({poi}) and ({acidic})", cutoff=cutoff)
    cmd.set("dash_color", "yellow", f"{name}_salt")

    # Hydrophobics: C-C contacts within 4.0 Å among hydrophobic residues
    hyd_res = "resn ALA+VAL+LEU+ILE+PHE+PRO+MET+TRP+TYR"
    cmd.delete(f"{name}_hydroph")
    cmd.distance(f"{name}_hydroph", f"({crbn}) and ({hyd_res}) and elem C", f"({poi}) and ({hyd_res}) and elem C", cutoff=cutoff)
    cmd.set("dash_color", "orange", f"{name}_hydroph")

    # Per-residue contact counts on POI
    # Count atoms within cutoff of CRBN heavy atoms
    contact_counts: Dict[Tuple[str, str, str], int] = defaultdict(int)
    crbn_model = cmd.get_model(crbn)
    poi_model = cmd.get_model(poi)

    # Build simple spatial index by binning (to avoid O(N^2) on big systems)
    bin_size = cutoff + 0.5
    bins: Dict[Tuple[int, int, int], List[Tuple[int, float, float, float]]] = defaultdict(list)
    for idx, a in enumerate(crbn_model.atom):
        key = (int(a.coord[0] // bin_size), int(a.coord[1] // bin_size), int(a.coord[2] // bin_size))
        bins[key].append((idx, a.coord[0], a.coord[1], a.coord[2]))

    def neighbor_keys(x, y, z):
        bx, by, bz = int(x // bin_size), int(y // bin_size), int(z // bin_size)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    yield (bx+dx, by+dy, bz+dz)

    c2 = cutoff * cutoff
    for a in poi_model.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        x, y, z = a.coord
        for k in neighbor_keys(x, y, z):
            for _, cx, cy, cz in bins.get(k, []):
                dx = x - cx; dy = y - cy; dz = z - cz
                if (dx*dx + dy*dy + dz*dz) <= c2:
                    contact_counts[key] += 1
                    break  # count per-atom at least one contact

    # Map counts to b-factors and color
    _set_b_factors(poi, contact_counts)
    cmd.show("surface", poi)
    _color_by_b(poi, palette="rainbow", ramp_name=f"{name}_resi_map")

    # Store a selection of interface residues for downstream use
    cmd.select(f"{name}_poi_resi", f"byres ({poi} within {cutoff} of {crbn})")
    cmd.select(f"{name}_crbn_resi", f"byres ({crbn} within {cutoff} of {poi})")

cmd.extend("interface_map", interface_map)


# ------------- ΔΔG Heatmap -------------

def _detect_foldx() -> Optional[str]:
    return _which("foldx")


def _prep_complex_tmp(crbn_sel: str, poi_sel: str) -> Tuple[str, Dict[Tuple[str,str,str], Tuple[str,str,str]]]:
    # Export CRBN+POI complex to PDB with consistent chain/resi mapping
    tmp = tempfile.mkdtemp(prefix="pymol_crbn_")
    pdb_path = os.path.join(tmp, "complex.pdb")
    # Create a temporary object combining selections
    tmp_obj = "__crbn_tmp_complex__"
    cmd.delete(tmp_obj)
    cmd.create(tmp_obj, f"({crbn_sel}) or ({poi_sel})")
    cmd.sort(tmp_obj)
    cmd.save(pdb_path, tmp_obj)
    # Build atom-to-resi mapping for later results
    mapping: Dict[Tuple[str,str,str], Tuple[str,str,str]] = {}
    model = cmd.get_model(poi_sel)
    for a in model.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        mapping[key] = key
    return pdb_path, mapping


def _foldx_alanine_scan(foldx_bin: str, pdb_path: str, poi_sel: str) -> Dict[Tuple[str,str,str], float]:
    # Interface residues on POI
    poi_iface = cmd.get_model(f"byres ({poi_sel} within 5.0 of not {poi_sel})")
    # Mutation list format for FoldX: "<PDB>; <chain><resi><icode><WT><Mut>"
    muts = []
    seen = set()
    for a in poi_iface.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        if key in seen:
            continue
        seen.add(key)
        wt = a.resn[:1]  # FoldX needs one-letter; map properly
        wt1 = _resn3_to_1(a.resn)
        if wt1 == "X":
            continue
        code = (a.chain or "A") + a.resi + (getattr(a, 'q', a.icode) or " ") + wt1 + "A"
        muts.append(code)
    if not muts:
        return {}

    tmpdir = os.path.dirname(pdb_path)
    mut_file = os.path.join(tmpdir, "mutations.txt")
    with open(mut_file, "w") as f:
        for m in muts:
            f.write(m + "\n")

    # Run FoldX BuildModel
    cmd_line = [foldx_bin, "--command=BuildModel", f"--pdb={os.path.basename(pdb_path)}", f"--pdb-dir={tmpdir}", f"--mutant-file={mut_file}"]
    try:
        subprocess.run(cmd_line, cwd=tmpdir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception:
        return {}

    # Parse DifferencesBetweenMutantAndWildType_fxout.csv if exists
    out_csv = os.path.join(tmpdir, "DifferencesBetweenMutantAndWildType_fxout.csv")
    ddg: Dict[Tuple[str,str,str], float] = {}
    if os.path.isfile(out_csv):
        with open(out_csv, newline='') as f:
            reader = csv.DictReader(f, delimiter=';')
            for row in reader:
                # Name like: <pdb>;<mutation> e.g., complex.pdb;A123A->A
                mut = row.get('Mutation', '') or row.get('mutation', '')
                # Extract chain/resi
                if len(mut) >= 5:
                    chain = mut[0]
                    resi = mut[1:4].strip()
                    icode = ""
                    key = (chain, resi, icode)
                    try:
                        ddg[key] = float(row.get('Total Energy', row.get('Total energy', '0')))
                    except Exception:
                        pass
    return ddg


def _resn3_to_1(resn: str) -> str:
    table = {
        'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
        'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    }
    return table.get(resn.upper(), 'X')


def _asa_ddg_proxy(poi_sel: str, complex_sel: str, scale: float = 0.025) -> Dict[Tuple[str,str,str], float]:
    # Use PyMOL get_area per residue in monomer vs complex; ΔASA scaled to ddG
    # Compute complex area and monomer area by hiding the partner temporarily
    # Prepare selections
    poi = f"({poi_sel})"
    comp = f"({complex_sel})"
    # Build list of residues
    model = cmd.get_model(poi)
    residues = []
    seen = set()
    for a in model.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        if key not in seen:
            seen.add(key)
            residues.append(key)

    ddg: Dict[Tuple[str,str,str], float] = {}
    # Ensure surface areas are computed with dot settings
    prev_dot_solvent = cmd.get("dot_solvent")
    prev_dot_density = cmd.get("dot_density")
    cmd.set("dot_solvent", 1)
    cmd.set("dot_density", 3)

    for (ch, resi, icode) in residues:
        sel_res = f"{poi} and chain {ch} and resi {resi}"
        asa_complex = cmd.get_area(sel_res, load_b=0, state=1)
        # To get monomer ASA, hide partner by excluding atoms within 100 Å (effectively isolate POI)
        # Simpler: duplicate POI only into temp object and measure
        tmp_obj = "__poi_tmp__"
        cmd.delete(tmp_obj)
        cmd.create(tmp_obj, sel_res)
        asa_monomer = cmd.get_area(tmp_obj, load_b=0, state=1)
        cmd.delete(tmp_obj)
        dASA = max(0.0, asa_monomer - asa_complex)
        ddg[(ch, resi, icode)] = -scale * dASA  # burial stabilizes binding (negative)

    # Restore settings
    cmd.set("dot_solvent", prev_dot_solvent)
    cmd.set("dot_density", prev_dot_density)
    return ddg


def ddg_heatmap(CRBN_sel: str, POI_sel: str, method: str = "auto", name: str = "ddg", scale: float = 0.025):
    """Color POI by ΔΔG. method: auto|foldx|asa. Stores value to b-factor and colors by spectrum.
    
    For publication-quality results, install FoldX:
    - Download from: https://foldxsuite.crg.eu/
    - Set environment variable: export FOLDX=/path/to/foldx
    - Or place foldx binary in your PATH
    
    If FoldX is not available, falls back to ASA-based proxy (less accurate).
    """
    poi = f"({POI_sel})"
    crbn = f"({CRBN_sel})"
    # Try FoldX if auto
    ddg: Dict[Tuple[str,str,str], float] = {}
    use_foldx = False
    if method in ("auto", "foldx"):
        fx = _detect_foldx()
        if fx:
            use_foldx = True
            print(f"[ddg_heatmap] ✅ Using FoldX at: {fx}")
            pdb_path, _ = _prep_complex_tmp(crbn, poi)
            ddg = _foldx_alanine_scan(fx, pdb_path, poi)
        elif method == "foldx":
            # User explicitly requested FoldX but it's not available
            print("[ddg_heatmap] ❌ FoldX not found in PATH or $FOLDX")
            print("[ddg_heatmap] 💡 Download FoldX from: https://foldxsuite.crg.eu/")
            print("[ddg_heatmap] 💡 Then set: export FOLDX=/path/to/foldx")
            if method == "foldx":  # Don't fallback if explicitly requested
                return
    
    if (not ddg) and method in ("auto", "asa"):
        # ASA proxy fallback
        print("[ddg_heatmap] ⚠️ Using ASA-based proxy (ΔΔG approximation, not publication quality)")
        print("[ddg_heatmap] 💡 For accurate ΔΔG values, install FoldX: https://foldxsuite.crg.eu/")
        ddg = _asa_ddg_proxy(poi, f"{crbn} or {poi}", scale=scale)

    if not ddg:
        cmd.feedback("pop", "all", "actions")
        print("[ddg_heatmap] ❌ No ΔΔG values computed. Check selections.")
        return

    _set_b_factors(poi, ddg)
    _color_by_b(poi, palette="blue_white_red", ramp_name=f"{name}_ramp")
    cmd.show("surface", poi)

cmd.extend("ddg_heatmap", ddg_heatmap)


# ------------- Score Color -------------

def score_color(POI_sel: str, scores_csv: str, col_chain: str = "chain", col_resi: str = "resi", col_score: str = "score", vmin: Optional[float] = None, vmax: Optional[float] = None):
    """Color POI by scores from CSV (chain,resi,score). Stores score into b-factor and colors by spectrum."""
    poi = f"({POI_sel})"
    if not os.path.isfile(scores_csv):
        raise FileNotFoundError(scores_csv)
    mapping: Dict[Tuple[str,str,str], float] = {}
    with open(scores_csv, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            chain = (row.get(col_chain, "")).strip() or "A"
            resi = str(row.get(col_resi, "")).strip()
            icode = ""
            try:
                val = float(row.get(col_score, 0.0))
            except Exception:
                continue
            mapping[(chain, resi, icode)] = val
    if not mapping:
        print("[score_color] No valid scores parsed.")
        return
    _set_b_factors(poi, mapping)
    _color_by_b(poi, palette="blue_white_red", min_val=vmin, max_val=vmax, ramp_name="score_ramp")
    cmd.show("surface", poi)

cmd.extend("score_color", score_color)


# ------------- Degron Annotation -------------

def _seq_of_selection(sel: str) -> Tuple[str, List[Tuple[str,str,str]]]:
    model = cmd.get_model(sel)
    residues: List[Tuple[str,str,str,str]] = []  # (chain,resi,icode,resn)
    seen = set()
    for a in model.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        if key not in seen:
            seen.add(key)
            residues.append((a.chain, a.resi, getattr(a, 'q', a.icode), a.resn))
    seq = ''.join(_resn3_to_1(r[3]) for r in residues)
    res_keys = [(r[0], r[1], r[2]) for r in residues]
    return seq, res_keys


def _find_c2h2(seq: str) -> List[Tuple[int,int]]:
    # Very rough C2H2 motif: C-X(2-4)-C-...-H-X(3-5)-H, window 20-40
    results = []
    n = len(seq)
    for i in range(n):
        if seq[i] != 'C':
            continue
        for j in range(i+2, min(i+5, n)):
            if seq[j] != 'C':
                continue
            for k in range(j+8, min(j+35, n)):
                if seq[k] != 'H':
                    continue
                for l in range(k+3, min(k+6, n)):
                    if seq[l] == 'H':
                        results.append((i, l))
                        break
    return results


def _find_beta_hairpins(seq: str) -> List[Tuple[int,int]]:
    # Heuristic: short segments (6-10) with alternating hydrophobicity signal
    hyd = set("VILMFYW")
    res = []
    n = len(seq)
    for i in range(n-6):
        window = seq[i:i+8]
        score = sum((1 if ((c in hyd) == (idx % 2 == 0)) else 0) for idx, c in enumerate(window))
        if score >= 6:
            res.append((i, i+7))
    return res


def degron_annotate(POI_sel: str, csv: Optional[str] = None, auto: bool = True, name: str = "degron"):
    poi = f"({POI_sel})"
    regions: List[Tuple[int,int,str]] = []  # start,end,label (0-based indices)

    # CSV regions (1-based resi numbers, chain-agnostic)
    if csv:
        if not os.path.isfile(csv):
            raise FileNotFoundError(csv)
        with open(csv, newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    s = int(str(row.get('start', '')).strip()) - 1
                    e = int(str(row.get('end', '')).strip()) - 1
                    label = (row.get('label', '') or 'degron').strip()
                    regions.append((s, e, label))
                except Exception:
                    continue

    if auto:
        seq, res_keys = _seq_of_selection(poi)
        for s, e in _find_c2h2(seq):
            regions.append((s, e, 'C2H2_like'))
        for s, e in _find_beta_hairpins(seq):
            regions.append((s, e, 'beta_hairpin_like'))

    if not regions:
        print("[degron_annotate] No regions to annotate.")
        return

    # Build selections and visuals
    # Map sequence index to (chain,resi,icode)
    _, res_keys = _seq_of_selection(poi)
    colors = {
        'C2H2_like': 'magenta',
        'beta_hairpin_like': 'tv_green',
        'degron': 'tv_red',
    }
    for idx, (s, e, label) in enumerate(regions, 1):
        s = max(0, s); e = max(s, e)
        keys = res_keys[s:e+1]
        if not keys:
            continue
        sel_parts = [f"(chain {ch} and resi {resi})" for (ch, resi, _ic) in keys]
        sel_expr = f"({poi}) and (" + " or ".join(sel_parts) + ")"
        sel_name = f"{name}_{label}_{idx}"
        cmd.select(sel_name, sel_expr)
        cmd.show("cartoon", sel_name)
        cmd.color(colors.get(label, 'yellow'), sel_name)
        # Label start-end
        cmd.label(f"first {sel_name}", f"\"{label}:{idx} start\"")
        cmd.label(f"last {sel_name}", f"\"{label}:{idx} end\"")

cmd.extend("degron_annotate", degron_annotate)


# ------------- Convenience: environment check -------------

def crbn_tools_doctor(verbose: int = 1):
    """Check environment and report optional components (FoldX)."""
    ok = True
    fx = _detect_foldx()
    if fx:
        print(f"[doctor] ✅ FoldX detected: {fx}")
    else:
        ok = False
        print("[doctor] ⚠️ FoldX not found")
        print("[doctor] 💡 ddg_heatmap will use ASA-based approximation (not publication quality)")
        print("[doctor] 💡 For publication-quality ΔΔG: download FoldX from https://foldxsuite.crg.eu/")
        print("[doctor] 💡 Then set environment variable: export FOLDX=/path/to/foldx")
    print(f"[doctor] PyMOL version: {cmd.get_version()[0]}")
    return ok

cmd.extend("crbn_tools_doctor", crbn_tools_doctor)


# ------------- Default visual settings -------------
cmd.set("dash_gap", 0.3)
cmd.set("dash_radius", 0.08)
cmd.set("cartoon_transparency", 0.2)
