# -*- coding: utf-8 -*-
"""
GlueTK benchmark runner (for JCIM manuscript).

Designed to be executed inside PyMOL:

    run validation/benchmark_analysis.py
    benchmark_all()

Writes:
    validation/results/benchmark_summary.csv
    validation/results/cases/<PDB>_ppi.csv
    validation/results/cases/<PDB>_neo.csv
"""

from __future__ import annotations

import csv
import os
import time
from dataclasses import dataclass
from typing import Iterable, Optional

from pymol import cmd

try:
    from gluetk.ppi_analyzer import ppi_analyze, neo_epitope_find
except Exception:  # pragma: no cover
    # If GlueTK isn't importable, the script still loads; user can fix sys.path / installation.
    ppi_analyze = None
    neo_epitope_find = None

try:
    from gluetk.interaction_analyzer import parse_pdb_structure, identify_molecule_type, distance
except Exception:  # pragma: no cover
    parse_pdb_structure = None
    identify_molecule_type = None
    distance = None


@dataclass(frozen=True)
class BenchmarkCase:
    pdb_id: str
    e3_chains: list[str]
    substrate_chains: list[str]
    expected_mechanism: str  # "Glue" or "PROTAC"
    glue_resname: Optional[str] = None  # if None, try to guess
    note: str = ""


# NOTE: Please adjust this list to match the exact PDB IDs / chain IDs used in the manuscript.
BENCHMARK_CASES: list[BenchmarkCase] = [
    # Molecular glues (examples)
    BenchmarkCase("6H0F", e3_chains=["A"], substrate_chains=["B"], expected_mechanism="Glue", glue_resname=None, note="CRBN ternary (glue)"),
    BenchmarkCase("6H0G", e3_chains=["A"], substrate_chains=["B"], expected_mechanism="Glue", glue_resname=None, note="CRBN ternary (glue)"),
    # PROTACs (examples; chain IDs vary by structure)
    BenchmarkCase("6BN7", e3_chains=["E"], substrate_chains=["A"], expected_mechanism="PROTAC", glue_resname=None, note="VHL-BRD4-dBET1"),
    BenchmarkCase("6SIS", e3_chains=["C"], substrate_chains=["A"], expected_mechanism="PROTAC", glue_resname=None, note="VHL-BRD4-MZ1"),
]


def _ensure_dirs(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _is_nontrivial_small_molecule(resn: str) -> bool:
    # Exclude common solvent/ions; keep this conservative.
    resn_u = (resn or "").upper()
    if resn_u in {"HOH", "WAT"}:
        return False
    if len(resn_u) <= 1:
        return False
    return True


def _guess_glue_resname(obj_name: str, e3_chains: Iterable[str], substrate_chains: Iterable[str], cutoff: float = 5.0) -> Optional[str]:
    """
    Heuristic: find a non-protein residue name whose atoms contact both E3 and substrate within cutoff.
    """
    if parse_pdb_structure is None or identify_molecule_type is None or distance is None:
        return None

    atoms = parse_pdb_structure(obj_name)
    if not atoms:
        return None

    e3_atoms = []
    sub_atoms = []
    ligand_atoms_by_resn: dict[str, list[tuple]] = {}

    for atom in atoms:
        chain, resn, _resi, _aname, coords = atom
        mol_type = identify_molecule_type(resn)
        if mol_type == "protein":
            if chain in set(e3_chains):
                e3_atoms.append(atom)
            elif chain in set(substrate_chains):
                sub_atoms.append(atom)
        else:
            if _is_nontrivial_small_molecule(resn):
                ligand_atoms_by_resn.setdefault(resn, []).append(atom)

    if not e3_atoms or not sub_atoms or not ligand_atoms_by_resn:
        return None

    best_resn = None
    best_bridge = -1
    for resn, lig_atoms in ligand_atoms_by_resn.items():
        bridge = 0
        for lig in lig_atoms:
            lig_xyz = lig[4]
            contacts_e3 = any(distance(lig_xyz, a[4]) <= cutoff for a in e3_atoms)
            if not contacts_e3:
                continue
            contacts_sub = any(distance(lig_xyz, a[4]) <= cutoff for a in sub_atoms)
            if contacts_sub:
                bridge += 1
        if bridge > best_bridge:
            best_bridge = bridge
            best_resn = resn

    if best_bridge <= 0:
        return None
    return best_resn


def _predict_mechanism(ppi_result: Optional[dict], neo_result: Optional[dict]) -> str:
    if neo_result and neo_result.get("is_molecular_glue") is True:
        return "Glue"
    if ppi_result:
        # Conservative fallback: very weak PPI suggests PROTAC-like linker mechanism.
        if (ppi_result.get("interface_contacts") or 0) < 3:
            return "PROTAC"
    return "Unknown"


def benchmark_all(
    out_csv: str = "validation/results/benchmark_summary.csv",
    fetch: bool = True,
    quiet: bool = False,
) -> str:
    """
    Run the default benchmark set and write a summary CSV.

    Returns the written CSV path.
    """
    if ppi_analyze is None or neo_epitope_find is None:
        raise RuntimeError("GlueTK functions not importable. Ensure GlueTK is on PyMOL's Python path.")

    out_dir = os.path.dirname(out_csv) or "."
    cases_dir = os.path.join(out_dir, "cases")
    _ensure_dirs(out_dir)
    _ensure_dirs(cases_dir)

    rows = []
    t0 = time.time()

    for case in BENCHMARK_CASES:
        obj = case.pdb_id
        if not quiet:
            print(f"\n[benchmark] === {case.pdb_id} ({case.expected_mechanism}) ===")

        cmd.reinitialize()
        if fetch:
            cmd.fetch(case.pdb_id, name=obj, type="pdb", async_=0)
        else:
            # User can pre-load structures manually with matching object names.
            if obj not in cmd.get_object_list():
                raise RuntimeError(f"Object '{obj}' not loaded and fetch=False")

        # Basic cleanup for reproducibility.
        cmd.remove("solvent")

        ppi_csv = os.path.join(cases_dir, f"{case.pdb_id}_ppi.csv")
        neo_csv = os.path.join(cases_dir, f"{case.pdb_id}_neo.csv")

        case_t0 = time.time()
        ppi = ppi_analyze(obj, case.e3_chains, case.substrate_chains, out_csv=ppi_csv)

        glue_resn = case.glue_resname or _guess_glue_resname(obj, case.e3_chains, case.substrate_chains)
        neo = None
        if glue_resn:
            neo = neo_epitope_find(obj, case.e3_chains, case.substrate_chains, glue_resn, out_csv=neo_csv)
        else:
            if not quiet:
                print("[benchmark] ⚠️ Could not guess glue resname; skipping neo-epitope for this case.")

        runtime_s = round(time.time() - case_t0, 2)
        predicted = _predict_mechanism(ppi, neo)
        correct = predicted == case.expected_mechanism

        rows.append(
            {
                "pdb_id": case.pdb_id,
                "expected": case.expected_mechanism,
                "predicted": predicted,
                "correct": int(bool(correct)),
                "e3_chains": "".join(case.e3_chains),
                "substrate_chains": "".join(case.substrate_chains),
                "glue_resname_used": glue_resn or "",
                "ppi_contacts": (ppi or {}).get("interface_contacts", ""),
                "ppi_strength": (ppi or {}).get("interface_strength", ""),
                "bsa": (ppi or {}).get("bsa", ""),
                "bridging_glue_atoms": (neo or {}).get("bridging_glue_atoms", ""),
                "neo_epitope_count": (neo or {}).get("neo_epitope_count", ""),
                "neo_confidence": (neo or {}).get("confidence", ""),
                "runtime_s": runtime_s,
                "note": case.note,
            }
        )

    fieldnames = list(rows[0].keys()) if rows else []
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    total = len(rows)
    acc = (sum(r["correct"] for r in rows) / total) if total else 0.0
    if not quiet:
        print(f"\n[benchmark] Wrote: {out_csv}")
        print(f"[benchmark] Accuracy: {acc:.3f} ({sum(r['correct'] for r in rows)}/{total})")
        print(f"[benchmark] Wall time: {round(time.time() - t0, 2)} s")
    return out_csv

