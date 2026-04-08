# -*- coding: utf-8 -*-
"""
CRBN G-motif (G-loop) detection.
Supports built-in structural templates and custom selections.

- template_mode: "builtin" (default, uses GSPT1) or "selection"
- In "builtin" mode, an internal PDB residue segment is used to extract 8×Cα atoms as the template;
  cmd.fetch(async_=0) may be invoked automatically when required.
- In "selection" mode, 8×Cα atoms are taken from template_sel (a selection expression or sele name).

Key behaviors:
- Handles insertion codes (for example 100A), altloc values (''/A), and residue sorting correctly.
- Default RMSD cutoff is 3.5 Å; pos6=Gly filtering can be disabled.
- Can print the top lowest RMSD values to help tune thresholds.
- CSV output is compatible with highlight_residues:
  Chain1,Residue1,Chain2,Residue2,Distance,Interaction
"""
from __future__ import print_function
import os, csv, re
from collections import defaultdict
from typing import Optional
from pymol import cmd

# ====== 1) Built-in template definitions (customize if needed) ======
# Template names must stay aligned with the GUI definitions in target_discovery.py
# Simplified protein names are used as the primary keys for GUI consumption

BUILTIN_TEMPLATES = {
    # ===== Primary templates (recommended) =====
    # GSPT1 G-loop from 5HXB (CRBN-pomalidomide-IKZF1 complex)
    # This is the most common molecular glue substrate template
    "GSPT1": {
        "pdb": "5HXB",
        "chain": "A",
        "resi_range": "570-577",
        "selection_fmt": "{obj} and chain A and resi 570+571+572+573+574+575+576+577 and name CA",
        "description": "GSPT1 G-loop from 5HXB (CRBN-pomalidomide-IKZF1 complex)",
    },
    # CK1α G-loop from 5FQD (CRBN-lenalidomide-CK1α complex)
    "CK1α": {
        "pdb": "5FQD",
        "chain": "C",
        "resi_range": "35-42",
        "selection_fmt": "{obj} and chain C and resi 35+36+37+38+39+40+41+42 and name CA",
        "description": "CK1α G-loop from 5FQD (CRBN-lenalidomide-CK1α complex)",
    },

    
    # ===== Legacy template aliases (kept for backward compatibility) =====
    "GSPT1 (5HXB A:570-577)": {
        "pdb": "5HXB",
        "chain": "A",
        "resi_range": "570-577",
        "selection_fmt": "{obj} and chain A and resi 570+571+572+573+574+575+576+577 and name CA",
        "description": "GSPT1 G-loop (legacy format)",
    },
    "CK1α (5FQD C:35-42)": {
        "pdb": "5FQD",
        "chain": "C",
        "resi_range": "35-42",
        "selection_fmt": "{obj} and chain C and resi 35+36+37+38+39+40+41+42 and name CA",
        "description": "CK1α G-loop from 5FQD (legacy format)",
    },
    "VAV1 (9NFR C:791-798)": {
        "pdb": "9NFR",
        "chain": "C",
        "resi_range": "791-798",
        "selection_fmt": "{obj} and chain C and resi 791+792+793+794+795+796+797+798 and name CA",
        "description": "VAV1 G-loop (legacy format)",
    },
}

# Template alias mapping (supports multiple naming styles)
TEMPLATE_ALIASES = {
    # Simplified name -> canonical name
    "gspt1": "GSPT1",
    "ck1a": "CK1α",
    "ck1alpha": "CK1α",
    "vav1": "VAV1",
    # GUI using的格式
    "GSPT1 (5HXB A:570-577)": "GSPT1",
    "CK1α (5FQD C:35-42)": "CK1α",
    "VAV1 (9NFR C:791-798)": "VAV1",
}

def get_builtin_template(name: str) -> Optional[dict]:
    """
    获取内置模板Configuration，支持多种命名格式
    
    Parameters:
        name: 模板Name（支持 "GSPT1", "gspt1", "GSPT1 (6H0G A:60-67)" 等格式）
    
    Return:
        dict: 模板Configuration，如果未找到Return None
    """
    # 直接匹配
    if name in BUILTIN_TEMPLATES:
        return BUILTIN_TEMPLATES[name]
    
    # 通过别名匹配
    normalized = name.strip()
    if normalized.lower() in TEMPLATE_ALIASES:
        canonical = TEMPLATE_ALIASES[normalized.lower()]
        return BUILTIN_TEMPLATES.get(canonical)
    
    # 尝试直接小写匹配
    for key in BUILTIN_TEMPLATES:
        if key.lower() == normalized.lower():
            return BUILTIN_TEMPLATES[key]
    
    return None

def list_builtin_templates() -> list:
    """
    列出所有可用的内置模板
    
    Return:
        list: [(name, description), ...]
    """
    result = []
    seen = set()
    for name, config in BUILTIN_TEMPLATES.items():
        # Skip legacy format entries to avoid duplicates
        if "legacy" in config.get("description", "").lower():
            continue
        if name not in seen:
            seen.add(name)
            result.append((name, config.get("description", "")))
    return result

# ====== CRBN key residue configuration ======
# Used for dynamic residue matching in validate_crbn_hbonds

# Default configuration (human CRBN, UniProt Q96SW2)
CRBN_KEY_RESIDUES_DEFAULT = {
    'N351': {
        'resn': 'ASN',           # Residue type
        'resi': '351',           # Residue number
        'atoms': ['ND2', 'OD1'], # Atoms participating in key H-bonds
        'description': 'Asn351 sidechain (NH2/O donor/acceptor)'
    },
    'H357': {
        'resn': 'HIS',
        'resi': '357',
        'atoms': ['ND1', 'NE2'],
        'description': 'His357 imidazole ring'
    },
    'W400': {
        'resn': 'TRP',
        'resi': '400',
        'atoms': ['NE1'],
        'description': 'Trp400 indole NH'
    },
}

# Common PDB presets for CRBN residue indices
# Format: 'PDB_CODE': {'N351': 'actual_resi', 'H357': 'actual_resi', 'W400': 'actual_resi'}
CRBN_PRESETS = {
    # CRBN-lenalidomide-GSPT1 ternary complex
    '6H0G': {'N351': '351', 'H357': '357', 'W400': '400'},
    # CRBN-pomalidomide-IKZF1 complex
    '5HXB': {'N351': '351', 'H357': '357', 'W400': '400'},
    # CRBN-lenalidomide-CK1α complex
    '5FQD': {'N351': '351', 'H357': '357', 'W400': '400'},
    # CRBN-thalidomide complex
    '4CI1': {'N351': '351', 'H357': '357', 'W400': '400'},
    # CRBN-CC-885-GSPT1 complex
    '5HXB': {'N351': '351', 'H357': '357', 'W400': '400'},
    # Mouse CRBN (add entries if residue numbering differs)
    # 'XXXX': {'N351': 'XXX', 'H357': 'XXX', 'W400': 'XXX'},
}

def get_crbn_key_residues(pdb_preset: Optional[str] = None,
                          custom_config: Optional[dict] = None) -> dict:
    """
    Return the CRBN key residue configuration.

    Priority: custom_config > pdb_preset > default

    Parameters:
        pdb_preset: PDB code (e.g. "6H0G"), uses a preset configuration
        custom_config: Custom configuration dictionary

    Return:
        dict: Key residue configuration
    """
    if custom_config:
        # Validate custom configuration
        required_keys = {'N351', 'H357', 'W400'}
        if not required_keys.issubset(custom_config.keys()):
            return CRBN_KEY_RESIDUES_DEFAULT.copy()

        # Build the full configuration
        result = {}
        for key in required_keys:
            if isinstance(custom_config[key], dict):
                result[key] = custom_config[key]
            else:
                # Simplified format: {'N351': '351'} -> full format
                base = CRBN_KEY_RESIDUES_DEFAULT[key].copy()
                base['resi'] = str(custom_config[key])
                result[key] = base
        return result

    if pdb_preset:
        preset_key = pdb_preset.upper()
        if preset_key in CRBN_PRESETS:
            # Build a full configuration from the preset
            preset = CRBN_PRESETS[preset_key]
            result = {}
            for key in ['N351', 'H357', 'W400']:
                base = CRBN_KEY_RESIDUES_DEFAULT[key].copy()
                base['resi'] = preset[key]
                result[key] = base
            return result

    return CRBN_KEY_RESIDUES_DEFAULT.copy()

# Amino-acid 3-letter to 1-letter mapping
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    # Non-standard or modified residues
    'MSE': 'M',  # Selenomethionine
    'SEC': 'U',  # Selenocysteine
    'PYL': 'O',  # Pyrrolysine
}

def _aa_3to1(resn: str) -> str:
    """Convert a 3-letter amino-acid code to a 1-letter code."""
    resn_upper = resn.strip().upper()
    return AA_3TO1.get(resn_upper, 'X')  # Unknown residues map to 'X'

def _parse_resi(resi_str: str):
    if resi_str is None:
        return (0, "")
    s = str(resi_str).strip()
    m = re.match(r'(-?\d+)\s*([A-Za-z]?)', s)
    if m:
        return (int(m.group(1)), m.group(2) or "")
    try:
        return (int(s), "")
    except Exception:
        return (0, s or "")

def _collect_ca_by_chain(obj_name: str):
    """Return {chain: [(sort_key, resi_label, resn, (x,y,z))...]}, altloc ''/A only."""
    by_chain = defaultdict(list)
    m = cmd.get_model(obj_name)
    for a in m.atom:
        if a.name.strip().upper() != 'CA':
            continue
        alt = (a.alt or '').strip().upper()
        if alt not in ("", "A"):
            continue
        chain = (a.chain or '').strip()
        resn = (a.resn or '').strip().upper()
        resi_label = (a.resi or '').strip()
        sort_key = _parse_resi(resi_label)
        by_chain[chain].append((sort_key, resi_label, resn, (a.coord[0], a.coord[1], a.coord[2])))
    for ch in by_chain:
        by_chain[ch].sort(key=lambda x: (x[0][0], x[0][1]))
    return by_chain

def _kabsch_rmsd(P, Q):
    import numpy as np
    P = np.asarray(P, float); Q = np.asarray(Q, float)
    Pc = P.mean(0); Qc = Q.mean(0)
    P0 = P - Pc; Q0 = Q - Qc
    C = P0.T @ Q0
    V, S, Wt = np.linalg.svd(C)
    if (np.linalg.det(V @ Wt) < 0.0):
        V[:, -1] *= -1.0
    R = V @ Wt
    P_aln = (R @ P0.T).T + Qc
    diff = P_aln - Q
    return float((diff**2).sum() / len(P))**0.5

# ====== 3) Template coordinate acquisition: idealized / built-in / selection ======
def _ideal_beta_hairpin_template():
    left = [(0.0,0.0,0.0),(3.8,0.2,0.0),(7.6,0.1,0.1),(11.4,0.0,0.0)]
    right = [(11.4,4.0,0.2),(7.6,4.2,0.0),(3.8,4.1,-0.1),(0.0,4.0,0.0)]
    return left + right  # 8×3

def _coords_from_selection(sel: str):
    """Collect 8×Cα coordinates from a selection (sorted by resi)."""
    if not sel:
        raise ValueError("template_sel is empty. Provide a selection containing 8 CA atoms.")
    sel_ca = f"({sel}) and name CA"
    m = cmd.get_model(sel_ca)
    items = []
    for a in m.atom:
        alt = (a.alt or '').strip().upper()
        if alt not in ("", "A"):
            continue
        resi_label = (a.resi or '').strip()
        sort_key = _parse_resi(resi_label)
        items.append((sort_key, (a.coord[0], a.coord[1], a.coord[2])))
    if len(items) < 8:
        raise ValueError(f"Insufficient CA atoms (<8) in template selection (got {len(items)}): {sel_ca}")
    items.sort(key=lambda x: (x[0][0], x[0][1]))
    return [xyz for _, xyz in items[:8]]

def _find_loaded_object_contains(code: str):
    code = (code or "").lower()
    try:
        objs = cmd.get_names("objects")
    except AttributeError:
        objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
        
    for obj in (objs or []):
        if code in obj.lower():
            return obj
    return None

def _coords_from_builtin(name: str):
    """
    Fetch 8×Cα coordinates from a built-in template.

    Parameters:
        name: Template name in one of the supported formats:
              - Canonical: "GSPT1", "CK1α", "VAV1"
              - Lowercase: "gspt1", "ck1a"
              - Legacy: "GSPT1 (6H0G A:60-67)"

    Return:
        list: Coordinates for 8 Cα atoms as (x, y, z)
    """
    # Use the unified template lookup helper
    info = get_builtin_template(name)
    if not info:
        available = [t[0] for t in list_builtin_templates()]
        raise ValueError(f"Unknown builtin template: '{name}'. Available templates: {available}")

    pdb_code = info["pdb"]

    # Auto-fetch the PDB if it is not already loaded
    obj = _find_loaded_object_contains(pdb_code)
    if obj is None:
        try:
            cmd.fetch(pdb_code, async_=0)  # Use async_ keyword for compatibility
            obj = _find_loaded_object_contains(pdb_code) or pdb_code
        except Exception as e:
            raise RuntimeError(
                f"Failed to fetch {pdb_code} ({e}). Load the PDB manually or use the 'selection' template."
            )

    sel = info["selection_fmt"].format(obj=obj)
    return _coords_from_selection(sel)

def _get_template_coords(template_mode: str, template_sel: Optional[str], template_builtin: Optional[str]):
    """
    Return template coordinates.

    Parameters:
        template_mode: "builtin" (default) or "selection"
        template_sel: Selection expression when template_mode == "selection"
        template_builtin: Built-in template name (default "GSPT1")

    Return:
        list: Coordinates for 8 Cα atoms as (x, y, z)
    """
    mode = (template_mode or "builtin").lower()
    if mode == "selection":
        return _coords_from_selection(template_sel or "")
    # Default to builtin mode with GSPT1
    return _coords_from_builtin(template_builtin or "GSPT1")

def _calculate_sasa_for_window(obj_name, chain, window):
    """
    Compute SASA for an 8-residue window.

    Parameters:
        obj_name: PyMOL object name
        chain: Chain ID
        window: [(resi_label, resn, xyz), ...] list of 8 residues

    Return:
        float: Mean SASA (Å²/residue), or None if the calculation fails
    """
    try:
        # Build residue list
        resi_list = '+'.join([w[0] for w in window])
        selection = f"{obj_name} and chain {chain} and resi {resi_list}"

        # Compute SASA
        total_sasa = cmd.get_area(selection, state=1)

        # Compute mean SASA
        avg_sasa = total_sasa / 8.0
        return avg_sasa
    except Exception:
        # Return None on failure (e.g. object missing)
        return None

# ====== 4) Main function ======
def find_crbn_g_motif(obj_name=None, pdb_file=None,
                       template_mode="builtin", template_sel=None, template_builtin="GSPT1",
                       rmsd_cutoff=1.0, out_csv=None, auto_highlight=0, require_gly_pos="6,3",
                       exclude_proline=False, check_surface_exposure=True,
                       min_sasa_per_residue=15.0, topk_debug=10,
                       highlight_surface=False, export_coords=False, coords_csv=None):
    """
    Main entry point for G-loop motif discovery.

    Parameters:
        obj_name: PyMOL object name
        pdb_file: PDB file path
        template_mode: "builtin" (default) or "selection"
        template_sel: Selection expression when template_mode == "selection"
        template_builtin: Built-in template name (default "GSPT1", options "CK1α", "VAV1")
        rmsd_cutoff: RMSD cutoff (default 1.0 Å)
        out_csv: Output CSV path
        auto_highlight: Auto-highlight toggle (0/1)
        require_gly_pos: Required glycine positions:
            - "6,3" or "3,6": glycine at position 6 or 3 (default)
            - "6": glycine at position 6 only
            - "3": glycine at position 3 only
            - None or "": no glycine requirement
        exclude_proline: Exclude windows containing proline (default False)
        check_surface_exposure: Check surface exposure (default True)
        min_sasa_per_residue: Minimum SASA threshold (default 15.0 Å²/residue)
        topk_debug: Show top N lowest RMSD values for debugging
        highlight_surface: Highlight G-loop surface (default False)
        export_coords: Export coordinates (default False)
        coords_csv: Optional coordinates CSV output path

    Return:
        dict: {
            'hits': [(chain, start_resi_label, end_resi_label, seq8, rmsd), ...],
            'csv_path': str,  # Output CSV path
            'coordinates': {...} or None,  # Coordinate data when export_coords=True
            'surface_info': {...} or None,  # Surface info when highlight_surface=True
        }

        For backward compatibility, you may iterate directly over hits-only results.
    """
    # Load object
    tmp_obj = None
    if pdb_file:
        tmp_obj = "_gmotif_tmp_obj"
        cmd.load(pdb_file, tmp_obj, quiet=1)
        obj = tmp_obj
    else:
        try:
            objs = cmd.get_names("objects")
        except AttributeError:
            objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
            
        obj = obj_name or (objs[0] if objs else None)
    if not obj:
        return []

    # Template coordinates
    try:
        tmpl = _get_template_coords(template_mode, template_sel, template_builtin)
    except Exception:
        tmpl = _coords_from_builtin("GSPT1")

    hits = []
    best_rmsd_pool = []
    by_chain = _collect_ca_by_chain(obj)
    total_windows = 0
    gly_pos_windows = 0
    proline_excluded_windows = 0
    buried_windows = 0

    # Parse require_gly_pos values
    gly_positions = []
    if require_gly_pos:
        if isinstance(require_gly_pos, str):
            # Supports "6,3", "3,6", "6", or "3"
            for pos in require_gly_pos.replace(" ", "").split(","):
                try:
                    gly_positions.append(int(pos))
                except ValueError:
                    pass
        elif isinstance(require_gly_pos, (int, float)):
            gly_positions.append(int(require_gly_pos))
        elif isinstance(require_gly_pos, (list, tuple)):
            gly_positions = [int(p) for p in require_gly_pos]

    for ch, rows in by_chain.items():
        seq = [(r[1], r[2], r[3]) for r in rows]  # (resi_label, resn, xyz)
        if len(seq) < 8:
            continue
        for i in range(0, len(seq) - 7):
            window = seq[i:i+8]
            total_windows += 1
            # Check required glycine positions (pos 6 or 3)
            if gly_positions:
                has_gly_at_required_pos = False
                for pos in gly_positions:
                    # Positions 1-8 map to indices 0-7
                    idx = pos - 1
                    if 0 <= idx < 8 and window[idx][1] in ("GLY", "G"):
                        has_gly_at_required_pos = True
                        break
                if not has_gly_at_required_pos:
                    continue
                gly_pos_windows += 1
            # Exclude windows containing proline
            if exclude_proline:
                has_proline = any(w[1] in ("PRO", "P") for w in window)
                if has_proline:
                    proline_excluded_windows += 1
                    continue
            P = [w[2] for w in window]
            try:
                rmsd = _kabsch_rmsd(P, tmpl)
            except Exception:
                continue
            best_rmsd_pool.append((rmsd, ch, window))
            if rmsd <= float(rmsd_cutoff):
                # Surface exposure check
                if check_surface_exposure:
                    avg_sasa = _calculate_sasa_for_window(obj, ch, window)
                    if avg_sasa is None or avg_sasa < min_sasa_per_residue:
                        buried_windows += 1
                        continue

                resi_s = window[0][0]; resi_e = window[-1][0]
                seq8 = ''.join(_aa_3to1(aa) if aa else 'X' for aa in [w[1] for w in window])
                hits.append((ch, resi_s, resi_e, seq8, rmsd))

    # Write CSV - output all scanned windows (not just hits)
    if out_csv is None:
        import tempfile
        fd, out_csv = tempfile.mkstemp(suffix="_gmotif.csv"); os.close(fd)

    # Sort best_rmsd_pool by RMSD
    best_rmsd_pool.sort(key=lambda x: x[0])

    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        # Add Status column
        w.writerow(["Chain", "Sequence", "Start", "End", "RMSD", "Status", "Type"])

        # Write all windows
        for rmsd, ch, window in best_rmsd_pool:
            resi_s = window[0][0]
            resi_e = window[-1][0]
            seq8 = ''.join(_aa_3to1(aa) if aa else 'X' for aa in [w[1] for w in window])

            # Determine pass/fail by RMSD cutoff
            status = "Pass" if rmsd <= float(rmsd_cutoff) else "Fail"

            w.writerow([ch, seq8, resi_s, resi_e, f"{rmsd:.2f}", status, "G-Motif"])

    # Sort hits by RMSD (ascending) so hits[0] is the best one
    hits.sort(key=lambda x: x[4])

    # Prepare return payload
    result = {
        'hits': hits,
        'csv_path': out_csv,
        'coordinates': None,
        'surface_info': None,
        'all_windows': best_rmsd_pool  # All scanned windows (for advanced analysis)
    }

    # Auto-highlight
    if auto_highlight and hits:
        try:
            try:
                from .highlight_residues import highlight_gmotif_loops
            except Exception:
                from highlight_residues import highlight_gmotif_loops
            # Use the dedicated G-loop highlight helper
            highlight_gmotif_loops(out_csv, obj=obj, color="yellow", show_labels=True, clear_old=True)
        except Exception:
            pass

    # Highlight surface (optional)
    if highlight_surface and hits:
        try:
            try:
                from .highlight_residues import highlight_gloop_surface
            except Exception:
                from highlight_residues import highlight_gloop_surface

            # Highlight the first hit surface
            first_hit = hits[0]
            ch, resi_s, resi_e, seq8, rmsd = first_hit
            surface_result = highlight_gloop_surface(
                obj=obj, chain=ch, start_resi=resi_s, end_resi=resi_e,
                surface_color="yellow", surface_transparency=0.3
            )
            result['surface_info'] = surface_result
        except Exception:
            pass

    # Export coordinates (optional)
    if export_coords and hits:
        try:
            try:
                from .highlight_residues import get_gloop_coordinates
            except Exception:
                from highlight_residues import get_gloop_coordinates

            # Export coordinates for the first hit
            first_hit = hits[0]
            ch, resi_s, resi_e, seq8, rmsd = first_hit

            # Auto-generate coordinates CSV if none was provided
            if coords_csv is None:
                import tempfile
                fd, coords_csv = tempfile.mkstemp(suffix="_gloop_coords.csv")
                os.close(fd)

            coords_result = get_gloop_coordinates(
                obj=obj, chain=ch, start_resi=resi_s, end_resi=resi_e,
                atom_types=["all"],  # Export all atoms
                output_csv=coords_csv
            )
            result['coordinates'] = coords_result
        except Exception:
            pass

    if tmp_obj:
        try:
            cmd.delete(tmp_obj)
        except Exception:
            pass

    # Return a list-iterable wrapper for backward compatibility
    return GMotifResult(result)


class GMotifResult:
    """
    G-motif result wrapper with backward-compatible iteration and dict-like access.
    """
    def __init__(self, data):
        self._data = data
        self._hits = data.get('hits', [])
    
    def __iter__(self):
        """Iterate over hits (for 'for hit in result')."""
        return iter(self._hits)

    def __len__(self):
        """Return the number of hits."""
        return len(self._hits)

    def __getitem__(self, key):
        """Support result[0] and result['hits'] access."""
        if isinstance(key, int):
            return self._hits[key]
        return self._data.get(key)

    def __bool__(self):
        """Allow truthiness checks like 'if result:'"""
        return len(self._hits) > 0

    @property
    def hits(self):
        """Return the hits list."""
        return self._hits

    @property
    def csv_path(self):
        """Return the CSV path."""
        return self._data.get('csv_path')

    @property
    def coordinates(self):
        """Return coordinate data."""
        return self._data.get('coordinates')

    @property
    def surface_info(self):
        """Return surface information."""
        return self._data.get('surface_info')

    @property
    def all_windows(self):
        """Return all scanned windows."""
        return self._data.get('all_windows', [])

    def get(self, key, default=None):
        """Dictionary-style getter."""
        return self._data.get(key, default)

    def to_dict(self):
        """Convert to a plain dictionary."""
        return self._data.copy()

    def __repr__(self):
        return f"GMotifResult(hits={len(self._hits)}, csv='{self.csv_path}')"


# ====== 5) G-motif internal geometry validation ======
def validate_g_motif_geometry(obj_name, g_motif_chain, g_motif_resi_range,
                               max_internal_hbond=3.5, reference_rmsd_threshold=1.0,
                               template_name: Optional[str] = None):
    """
    Validate internal G-motif geometry (α-turn).

    Supports two G-loop types:
    - Type A (Gly at position 6): G₋₄ backbone O → G₀ backbone N (indices 1 → 5)
    - Type B (Gly at position 3): G₀ backbone O → G₊₄ backbone N (indices 2 → 6)

    Parameters:
        obj_name: PyMOL object name
        g_motif_chain: Chain ID of the G-motif
        g_motif_resi_range: Residue range (start, end) or "start-end"
        max_internal_hbond: Max hydrogen bond distance (default 3.5 Å)
        reference_rmsd_threshold: RMSD threshold (default 1.0 Å)
        template_name: Template name for RMSD comparison (default None, auto-detect)

    Return:
        dict: {
            'gloop_type': str,  # "Type_A" (Gly@6) or "Type_B" (Gly@3)
            'gly_position': int,  # Glycine position (3 or 6)
            'has_alpha_turn_hbond': bool,
            'alpha_turn_distance': float,
            'alpha_turn_description': str,  # H-bond description
            'has_secondary_hbond': bool,
            'secondary_hbond_distance': float,
            'secondary_hbond_description': str,
            'gly_confirmed': bool,
            'rmsd_to_template': float,
            'template_used': str,
            'is_valid_geometry': bool
        }
    """
    import numpy as np

    # Parse range
    if isinstance(g_motif_resi_range, str):
        parts = g_motif_resi_range.split('-')
        start_resi = parts[0].strip()
        end_resi = parts[1].strip() if len(parts) > 1 else start_resi
    else:
        start_resi, end_resi = g_motif_resi_range
    
    g_sel = f"{obj_name} and chain {g_motif_chain} and resi {start_resi}-{end_resi}"
    
    # Collect backbone atoms for 8 residues
    g_model = cmd.get_model(g_sel)
    g_residues = {}
    for a in g_model.atom:
        resi = (a.resi or '').strip()
        sort_key = _parse_resi(resi)
        if resi not in g_residues:
            g_residues[resi] = {'sort': sort_key, 'atoms': {}, 'resn': (a.resn or '').strip().upper()}
        aname = (a.name or '').strip().upper()
        if aname in ('N', 'CA', 'C', 'O'):
            g_residues[resi]['atoms'][aname] = np.array([a.coord[0], a.coord[1], a.coord[2]])
    
    sorted_resis = sorted(g_residues.keys(), key=lambda r: (g_residues[r]['sort'][0], g_residues[r]['sort'][1]))
    if len(sorted_resis) < 8:
        print(f"[G-motif Geometry] ⚠️ Only {len(sorted_resis)} residues, expected 8")
        return None
    
    # 检测 G-loop Type：检查第 3 位和第 6 位是否为甘氨酸
    # 位置 1-8 对应索引 0-7
    gly_at_pos3 = g_residues[sorted_resis[2]]['resn'] in ('GLY', 'G')  # 索引 2 = 位置 3
    gly_at_pos6 = g_residues[sorted_resis[5]]['resn'] in ('GLY', 'G')  # 索引 5 = 位置 6
    
    # 确定 G-loop Type
    if gly_at_pos6:
        gloop_type = "Type_A"
        gly_position = 6
        gly_confirmed = True
        gly_resn = g_residues[sorted_resis[5]]['resn']
        # Type A: G₋₄ O → G₀ N (索引 1 → 5)
        hbond_donor_idx = 1   # G₋₄
        hbond_acceptor_idx = 5  # G₀
        secondary_acceptor_idx = 6  # G₊₁
        alpha_turn_desc = "G₋₄ O → G₀ N (索引 1 → 5)"
        secondary_desc = "G₋₄ O → G₊₁ N (索引 1 → 6)"
    elif gly_at_pos3:
        gloop_type = "Type_B"
        gly_position = 3
        gly_confirmed = True
        gly_resn = g_residues[sorted_resis[2]]['resn']
        # Type B: G₀ O → G₊₄ N (索引 2 → 6)
        hbond_donor_idx = 2   # G₀ (位置 3)
        hbond_acceptor_idx = 6  # G₊₄ (位置 7)
        secondary_acceptor_idx = 7  # G₊₅ (位置 8)
        alpha_turn_desc = "G₀ O → G₊₄ N (索引 2 → 6)"
        secondary_desc = "G₀ O → G₊₅ N (索引 2 → 7)"
    else:
        gloop_type = "Unknown"
        gly_position = 0
        gly_confirmed = False
        gly_resn = "N/A"
        # 默认using Type A 逻辑
        hbond_donor_idx = 1
        hbond_acceptor_idx = 5
        secondary_acceptor_idx = 6
        alpha_turn_desc = "G₋₄ O → G₀ N (default)"
        secondary_desc = "G₋₄ O → G₊₁ N (default)"
    
    result = {
        'gloop_type': gloop_type,
        'gly_position': gly_position,
        'has_alpha_turn_hbond': False,
        'alpha_turn_distance': None,
        'alpha_turn_description': alpha_turn_desc,
        'has_secondary_hbond': False,
        'secondary_hbond_distance': None,
        'secondary_hbond_description': secondary_desc,
        'gly_confirmed': gly_confirmed,
        'rmsd_to_template': None,
        'template_used': template_name or "auto",
        'is_valid_geometry': False
    }
    
    # 计算主要 α-turn 氢Key
    donor_resi = sorted_resis[hbond_donor_idx]
    acceptor_resi = sorted_resis[hbond_acceptor_idx]
    if 'O' in g_residues[donor_resi]['atoms'] and 'N' in g_residues[acceptor_resi]['atoms']:
        o_donor = g_residues[donor_resi]['atoms']['O']
        n_acceptor = g_residues[acceptor_resi]['atoms']['N']
        dist = np.linalg.norm(o_donor - n_acceptor)
        result['alpha_turn_distance'] = round(float(dist), 2)
        result['has_alpha_turn_hbond'] = (dist <= max_internal_hbond)
    
    # 计算次级氢Key
    if secondary_acceptor_idx < len(sorted_resis):
        secondary_resi = sorted_resis[secondary_acceptor_idx]
        if 'O' in g_residues[donor_resi]['atoms'] and 'N' in g_residues[secondary_resi]['atoms']:
            o_donor = g_residues[donor_resi]['atoms']['O']
            n_secondary = g_residues[secondary_resi]['atoms']['N']
            dist = np.linalg.norm(o_donor - n_secondary)
            result['secondary_hbond_distance'] = round(float(dist), 2)
            result['has_secondary_hbond'] = (dist <= max_internal_hbond)
    
    # RMSD 与模板比较（using Cα）
    ca_coords = []
    for resi in sorted_resis[:8]:
        if 'CA' in g_residues[resi]['atoms']:
            ca_coords.append(g_residues[resi]['atoms']['CA'])
    
    if len(ca_coords) == 8:
        # 获取模板坐标
        template_coords = None
        template_used = "idealized"
        
        if template_name:
            # 尝试using指定的模板
            try:
                template_coords = _coords_from_builtin(template_name)
                template_used = template_name
            except Exception as e:
                print(f"[G-motif Geometry] ⚠️ 无法Load模板 '{template_name}': {e}")
                template_coords = None
        
        if template_coords is None:
            # using理想化模板
            template_coords = _ideal_beta_hairpin_template()
            template_used = "idealized"
        
        result['template_used'] = template_used
        
        try:
            rmsd = _kabsch_rmsd(np.array(ca_coords), np.array(template_coords))
            result['rmsd_to_template'] = round(rmsd, 3)
        except Exception:
            pass
    
    # 综合判定
    result['is_valid_geometry'] = (
        result['gly_confirmed'] and
        result['has_alpha_turn_hbond'] and
        (result['rmsd_to_template'] is None or result['rmsd_to_template'] < reference_rmsd_threshold)
    )
    
    # PrintResults
    print("\n" + "=" * 70)
    print("G-motif 内部几何验证 (G-motif Internal Geometry Validation)")
    print("=" * 70)
    print(f"G-loop Type: {gloop_type} (Gly @ 位置 {gly_position})")
    print("-" * 70)
    
    gly_status = "✅" if result['gly_confirmed'] else "❌"
    print(f"{gly_status} 甘氨酸Confirm:  {gly_resn} @ 位置 {gly_position}")
    
    alpha_status = "✅" if result['has_alpha_turn_hbond'] else "❌"
    alpha_dist = f"{result['alpha_turn_distance']} Å" if result['alpha_turn_distance'] else "N/A"
    print(f"{alpha_status} α-turn 氢Key:  {alpha_dist}")
    print(f"   {alpha_turn_desc}")
    
    if result['secondary_hbond_distance']:
        sec_status = "✅" if result['has_secondary_hbond'] else "➖"
        sec_dist = f"{result['secondary_hbond_distance']} Å"
        print(f"{sec_status} 次级氢Key:  {sec_dist}")
        print(f"   {secondary_desc}")
    
    if result['rmsd_to_template']:
        rmsd_status = "✅" if result['rmsd_to_template'] < reference_rmsd_threshold else "⚠️"
        print(f"{rmsd_status} RMSD vs {template_used}:  {result['rmsd_to_template']} Å")
    
    valid_status = "✅ 是" if result['is_valid_geometry'] else "❌ 否"
    print(f"\n有效 α-turn 几何: {valid_status}")
    print("=" * 70 + "\n")
    
    return result


# ====== 6) CRBN关Key残基氢Key验证 ======
def validate_crbn_hbonds(obj_name, g_motif_chain, g_motif_resi_range, crbn_chain,
                         max_hbond_dist=3.5, min_donor_angle=120.0,
                         pdb_preset: Optional[str] = None,
                         crbn_key_residues: Optional[dict] = None):
    """
    验证 G-loop 与 CRBN 的 3 个关Key氢Key（基于 Schrödinger 标准）：
    - G-3 backbone O → CRBN Asn351 (sidechain NH2)
    - G-2 backbone O → CRBN His357 (sidechain)
    - G-1 backbone O → CRBN Trp400 (sidechain NH)
    
    Parameters:
        obj_name: PyMOL 对象名
        g_motif_chain: G-motif 所在链 ID
        g_motif_resi_range: G-motif 残基范围 (start, end) 或 "start-end"
        crbn_chain: CRBN 链 ID
        max_hbond_dist: 氢Key最大距离阈Value（默认 3.5 Å）
        min_donor_angle: 最小供体角degrees（默认 120°，当前未using）
        pdb_preset: PDB 代码预设（如 '6H0G'），自动using对应的残基Number
        crbn_key_residues: 自定义 CRBN 关Key残基Configuration
            简化格式: {'N351': '351', 'H357': '357', 'W400': '400'}
            完整格式: {'N351': {'resn': 'ASN', 'resi': '351', 'atoms': ['ND2', 'OD1']}, ...}
    
    Return:
        dict: {
            'N351_hbond': {'found': bool, 'distance': float, 'g_pos': int},
            'H357_hbond': {'found': bool, 'distance': float, 'g_pos': int},
            'W400_hbond': {'found': bool, 'distance': float, 'g_pos': int},
            'total_hbonds': int,
            'is_canonical_gloop': bool,
            'config_used': dict  # using的ConfigurationInformation
        }
    """
    import numpy as np
    
    # 解析 G-motif 范围
    if isinstance(g_motif_resi_range, str):
        parts = g_motif_resi_range.split('-')
        start_resi = parts[0].strip()
        end_resi = parts[1].strip() if len(parts) > 1 else start_resi
    else:
        start_resi, end_resi = g_motif_resi_range
    
    g_sel = f"{obj_name} and chain {g_motif_chain} and resi {start_resi}-{end_resi}"
    crbn_sel = f"{obj_name} and chain {crbn_chain} and polymer.protein"
    
    # 获取 G-loop 8 个残基的 backbone O
    g_model = cmd.get_model(g_sel)
    g_residues = {}
    for a in g_model.atom:
        resi = (a.resi or '').strip()
        sort_key = _parse_resi(resi)
        if resi not in g_residues:
            g_residues[resi] = {'sort': sort_key, 'atoms': []}
        g_residues[resi]['atoms'].append(a)
    
    sorted_resis = sorted(g_residues.keys(), key=lambda r: (g_residues[r]['sort'][0], g_residues[r]['sort'][1]))
    if len(sorted_resis) < 8:
        print(f"[CRBN H-bond] ⚠️ G-loop has only {len(sorted_resis)} residues, expected 8")
        return None
    
    # Locate G-3, G-2, G-1 的 backbone O (对应位置 5, 6, 7 in 0-indexed)
    g_positions = {}
    for i, resi in enumerate(sorted_resis[:8]):
        backbone_o = None
        for atom in g_residues[resi]['atoms']:
            if atom.name.strip().upper() == 'O':  # backbone carbonyl
                backbone_o = np.array([atom.coord[0], atom.coord[1], atom.coord[2]])
                break
        g_positions[i] = {'resi': resi, 'O': backbone_o}
    
    # 获取 CRBN 关Key残基Configuration（支持动态Configuration）
    key_residue_config = get_crbn_key_residues(pdb_preset=pdb_preset, custom_config=crbn_key_residues)
    
    # 获取 CRBN 关Key残基的侧链原子
    crbn_model = cmd.get_model(crbn_sel)
    crbn_key_atoms = {'N351': [], 'H357': [], 'W400': []}
    
    for a in crbn_model.atom:
        resi_str = (a.resi or '').strip()
        resn = (a.resn or '').strip().upper()
        aname = (a.name or '').strip().upper()
        coord = np.array([a.coord[0], a.coord[1], a.coord[2]])
        
        # 动态匹配 CRBN 关Key残基
        for key in ['N351', 'H357', 'W400']:
            config = key_residue_config[key]
            expected_resn = config['resn']
            expected_resi = config['resi']
            expected_atoms = config['atoms']
            
            # 匹配残基Type和Number
            if resn == expected_resn and expected_resi in resi_str:
                if aname in expected_atoms:
                    crbn_key_atoms[key].append((f"{resn}{resi_str}", aname, coord))
    
    # 检查 3 个氢Key
    result = {
        'N351_hbond': {'found': False, 'distance': None, 'g_pos': -3},
        'H357_hbond': {'found': False, 'distance': None, 'g_pos': -2},
        'W400_hbond': {'found': False, 'distance': None, 'g_pos': -1},
        'total_hbonds': 0,
        'is_canonical_gloop': False,
        'config_used': {
            'N351_resi': key_residue_config['N351']['resi'],
            'H357_resi': key_residue_config['H357']['resi'],
            'W400_resi': key_residue_config['W400']['resi'],
            'preset': pdb_preset or 'default'
        }
    }
    
    # G-3 (index 5) → N351
    if 5 in g_positions and g_positions[5]['O'] is not None:
        for crbn_atom in crbn_key_atoms['N351']:
            dist = np.linalg.norm(g_positions[5]['O'] - crbn_atom[2])
            if dist <= max_hbond_dist:
                result['N351_hbond'] = {'found': True, 'distance': round(float(dist), 2), 'g_pos': -3}
                result['total_hbonds'] += 1
                break
    
    # G-2 (index 6) → H357
    if 6 in g_positions and g_positions[6]['O'] is not None:
        for crbn_atom in crbn_key_atoms['H357']:
            dist = np.linalg.norm(g_positions[6]['O'] - crbn_atom[2])
            if dist <= max_hbond_dist:
                result['H357_hbond'] = {'found': True, 'distance': round(float(dist), 2), 'g_pos': -2}
                result['total_hbonds'] += 1
                break
    
    # G-1 (index 7) → W400
    if 7 in g_positions and g_positions[7]['O'] is not None:
        for crbn_atom in crbn_key_atoms['W400']:
            dist = np.linalg.norm(g_positions[7]['O'] - crbn_atom[2])
            if dist <= max_hbond_dist:
                result['W400_hbond'] = {'found': True, 'distance': round(float(dist), 2), 'g_pos': -1}
                result['total_hbonds'] += 1
                break
    
    # 判定是否为标准 G-loop (至少 2/3 氢Key)
    result['is_canonical_gloop'] = result['total_hbonds'] >= 2
    
    # PrintResults
    print("\n" + "=" * 70)
    print("CRBN G-loop 氢Key验证 (CRBN-G-loop H-bond Validation)")
    print("=" * 70)
    print(f"Configuration: {result['config_used']['preset']} "
          f"(N351={result['config_used']['N351_resi']}, "
          f"H357={result['config_used']['H357_resi']}, "
          f"W400={result['config_used']['W400_resi']})")
    print("-" * 70)
    
    # 动态生成Label
    n351_label = f"G-3 → Asn{key_residue_config['N351']['resi']}"
    h357_label = f"G-2 → His{key_residue_config['H357']['resi']}"
    w400_label = f"G-1 → Trp{key_residue_config['W400']['resi']}"
    
    for key, label in [('N351_hbond', n351_label),
                        ('H357_hbond', h357_label),
                        ('W400_hbond', w400_label)]:
        hb = result[key]
        status = "✅" if hb['found'] else "❌"
        dist_str = f"{hb['distance']} Å" if hb['distance'] else "N/A"
        print(f"{status} {label:20s}  距离: {dist_str}")
    print(f"\n总氢Key数: {result['total_hbonds']}/3")
    canonical = "✅ 是" if result['is_canonical_gloop'] else "⚠️ 否"
    print(f"标准 G-loop: {canonical}")
    print("=" * 70 + "\n")
    
    return result


# ====== 7) Glue结合验证功能（改进版）======
def analyze_g_motif_glue_binding(obj_name,
                                  g_motif_chain,
                                  g_motif_resi_range,  # tuple: (start, end) or "start-end"
                                  glue_resname,
                                  crbn_chain,
                                  distance_threshold=5.0,
                                  validate_hbonds=True,
                                  validate_geometry=True):
    """
    验证检测到的 G-motif 是否真的是分子胶底物（基于 Annual Review 2023）
    
    检测逻辑:
    1. 内部几何验证：α-turn 氢Key (G₋₄→G₀)
    2. CRBN 关Key氢Key验证 (N351, H357, W400)
    3. Glue 是否插入 G-motif 和 CRBN 之间
    4. 关Key残基（如 Gly₀）与 MGD 的 vdW 接触
    
    Parameters:
        obj_name: PyMOL对象Name
        g_motif_chain: G-motif所在链 ID
        g_motif_resi_range: G-motif残基范围 (start, end) 或 "start-end"
        glue_resname: 分子胶Name
        crbn_chain: CRBN链 ID
        distance_threshold: 接触距离阈Value（埃）
    
    Return:
        dict: {
            'is_glue_substrate': bool,  # 是否为分子胶底物
            'binding_mode': str,  # "glue-induced" / "direct" / "no-binding"
            'g_motif_crbn_distance': float,  # G-motif到CRBN最短距离
            'glue_contacts_g_motif': int,  # Glue与G-motifcontact count
            'glue_contacts_crbn': int,  # Glue与CRBNcontact count
            'key_residues_engaged': list,  # 参与相互作用的关Key残基
            'confidence': float  # 置信degrees [0-1]
        }
    """
    try:
        objs = cmd.get_names("objects")
    except AttributeError:
        objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
        
    if obj_name not in objs:
        print(f"[analyze_g_motif_glue_binding] Object '{obj_name}' not found")
        return None
    
    # 解析残基范围
    if isinstance(g_motif_resi_range, str):
        parts = g_motif_resi_range.split('-')
        start_resi = parts[0].strip()
        end_resi = parts[1].strip() if len(parts) > 1 else start_resi
    else:
        start_resi, end_resi = g_motif_resi_range
    
    print(f"[G-motif Glue Binding] Analyzing: {g_motif_chain}:{start_resi}-{end_resi} with {glue_resname} and CRBN {crbn_chain}")
    
    # 1. 验证内部几何（如果Enable）
    geometry_result = None
    if validate_geometry:
        geometry_result = validate_g_motif_geometry(obj_name, g_motif_chain, g_motif_resi_range)
        if geometry_result and not geometry_result['is_valid_geometry']:
            print("⚠️ Warning: α-turn 几何不符合，可能不是有效 G-motif\n")
    
    # 2. 验证 CRBN 氢Key（如果Enable）
    hbond_result = None
    if validate_hbonds:
        hbond_result = validate_crbn_hbonds(obj_name, g_motif_chain, g_motif_resi_range, crbn_chain)
        if hbond_result and not hbond_result['is_canonical_gloop']:
            print("⚠️ Warning: CRBN 关Key氢Key不足，可能不是标准 G-loop\n")
    
    # 获取原子坐标
    g_motif_sel = f"{obj_name} and chain {g_motif_chain} and resi {start_resi}-{end_resi}"
    crbn_sel = f"{obj_name} and chain {crbn_chain} and polymer.protein"
    glue_sel = f"{obj_name} and resn {glue_resname}"
    
    # 检查Glue是否存在
    glue_model = cmd.get_model(glue_sel)
    if len(glue_model.atom) == 0:
        print(f"[G-motif Glue Binding] ⚠️ Glue molecule '{glue_resname}' not found")
        return None
    
    # 获取原子
    g_motif_atoms = [(a.coord[0], a.coord[1], a.coord[2], a.resn, a.resi, a.name) for a in cmd.get_model(g_motif_sel).atom]
    crbn_atoms = [(a.coord[0], a.coord[1], a.coord[2], a.resn, a.resi, a.name) for a in cmd.get_model(crbn_sel).atom]
    glue_atoms = [(a.coord[0], a.coord[1], a.coord[2], a.resn, a.resi, a.name) for a in glue_model.atom]
    
    if not g_motif_atoms or not crbn_atoms:
        print("[G-motif Glue Binding] ⚠️ Missing atoms for analysis")
        return None
    
    # 1. 计算 G-motif 与 CRBN 的最短距离
    min_g_crbn_dist = float('inf')
    for g_atom in g_motif_atoms:
        for c_atom in crbn_atoms:
            d = _distance_3d(g_atom[:3], c_atom[:3])
            min_g_crbn_dist = min(min_g_crbn_dist, d)
    
    # 2. 计算 Glue 与 G-motif 和 CRBN 的接触
    glue_contacts_g_motif = 0
    glue_contacts_crbn = 0
    glue_g_contacts_list = []  # 记录具体接触
    glue_c_contacts_list = []
    
    for glue_atom in glue_atoms:
        # 与 G-motif 接触
        for g_atom in g_motif_atoms:
            d = _distance_3d(glue_atom[:3], g_atom[:3])
            if d <= distance_threshold:
                glue_contacts_g_motif += 1
                glue_g_contacts_list.append((g_atom[3], g_atom[4], d))  # (resn, resi, dist)
                break  # 每个Glue原子只计数一次
        
        # 与 CRBN 接触
        for c_atom in crbn_atoms:
            d = _distance_3d(glue_atom[:3], c_atom[:3])
            if d <= distance_threshold:
                glue_contacts_crbn += 1
                glue_c_contacts_list.append((c_atom[3], c_atom[4], d))
                break
    
    # 3. 分析关Key残基（Gly6 及侧链位点 G-4, G-3, G-2, G+1）
    key_residues_engaged = []
    # 获取 G-motif 序列位置
    g_motif_residues = defaultdict(list)
    for atom in g_motif_atoms:
        res_key = (atom[3], atom[4])  # (resn, resi)
        g_motif_residues[res_key].append(atom)
    
    # 按resiSort
    sorted_residues = sorted(g_motif_residues.keys(), key=lambda x: _parse_resi(x[1]))
    
    # 检查第6位（Gly，G位）与 MGD 的 vdW 接触
    vdw_threshold = 4.5  # vdW 接触阈Value
    if len(sorted_residues) >= 6:
        gly6_key = sorted_residues[5]  # 0-indexed, 对应 G 位
        gly6_atoms = g_motif_residues[gly6_key]
        
        # 检查Gly（G位）是否与Glue有 vdW 接触（文献指出保守Gly面向MGD）
        gly_vdw_glue = False
        gly_close_glue = False
        gly_min_dist = float('inf')
        
        for gly_atom in gly6_atoms:
            for glue_atom in glue_atoms:
                d = _distance_3d(gly_atom[:3], glue_atom[:3])
                gly_min_dist = min(gly_min_dist, d)
                if d <= vdw_threshold:
                    gly_vdw_glue = True
                if d <= distance_threshold:
                    gly_close_glue = True
        
        if gly_vdw_glue:
            key_residues_engaged.append({
                'position': 'G (pos 6)',
                'residue': f"{gly6_key[0]} {gly6_key[1]}",
                'contacts_glue_vdw': True,
                'min_distance': round(gly_min_dist, 2),
                'note': '保守Gly面向MGD'
            })
    
    # 检查侧链位点 (G-4, G-3, G-2, G+1) 与 CRBN 的接触
    sidechain_positions = [1, 2, 3, 6]  # 0-indexed: G-4, G-3, G-2, G+1
    for idx in sidechain_positions:
        if idx >= len(sorted_residues):
            continue
        res_key = sorted_residues[idx]
        res_atoms = g_motif_residues[res_key]
        
        contacts_crbn = False
        for res_atom in res_atoms:
            # 排除 backbone atoms
            if res_atom[5].strip().upper() in ('N', 'CA', 'C', 'O'):
                continue
            for c_atom in crbn_atoms:
                if _distance_3d(res_atom[:3], c_atom[:3]) <= distance_threshold:
                    contacts_crbn = True
                    break
            if contacts_crbn:
                break
        
        if contacts_crbn:
            g_pos_label = ['G-4', 'G-3', 'G-2', 'G+1'][sidechain_positions.index(idx)]
            key_residues_engaged.append({
                'position': g_pos_label,
                'residue': f"{res_key[0]} {res_key[1]}",
                'contacts_crbn_sidechain': True,
                'note': '侧链与CRBN接触'
            })
    
    # 4. 判定结合模式
    binding_mode = "no-binding"
    is_glue_substrate = False
    confidence = 0.0
    
    # 检测是否为Glue诱导结合
    if glue_contacts_g_motif >= 2 and glue_contacts_crbn >= 2:
        # Glue同时与G-motif和CRBN接触 → 分子胶机制
        binding_mode = "glue-induced"
        is_glue_substrate = True
        confidence = min(1.0, (glue_contacts_g_motif + glue_contacts_crbn) / 10.0)
        
        # 如果关Key残基参与，提高置信degrees
        if key_residues_engaged:
            confidence = min(1.0, confidence + 0.2)
    
    elif min_g_crbn_dist <= distance_threshold and glue_contacts_g_motif == 0:
        # G-motif直接与CRBN接触，但Glue不参与 → 直接结合（非分子胶）
        binding_mode = "direct"
        is_glue_substrate = False
        confidence = 0.3
    
    elif glue_contacts_g_motif >= 2 or glue_contacts_crbn >= 2:
        # 部分接触
        binding_mode = "partial"
        is_glue_substrate = (glue_contacts_g_motif >= 2 and glue_contacts_crbn >= 1)
        confidence = 0.5 if is_glue_substrate else 0.2
    
    result = {
        'is_glue_substrate': is_glue_substrate,
        'binding_mode': binding_mode,
        'g_motif_crbn_distance': round(min_g_crbn_dist, 2) if min_g_crbn_dist != float('inf') else None,
        'glue_contacts_g_motif': glue_contacts_g_motif,
        'glue_contacts_crbn': glue_contacts_crbn,
        'key_residues_engaged': key_residues_engaged,
        'confidence': round(confidence, 2),
        'geometry_validation': geometry_result,  # Add几何验证Results
        'crbn_hbond_validation': hbond_result     # Add氢Key验证Results
    }
    
    # PrintResults
    print("\n" + "=" * 60)
    print("G-motif Glue结合验证 (G-motif Glue Binding Analysis)")
    print("=" * 60)
    print(f"G-motif 到 CRBN 距离: {result['g_motif_crbn_distance']} Å")
    print(f"Glue - G-motif contact count: {glue_contacts_g_motif}")
    print(f"Glue - CRBN contact count: {glue_contacts_crbn}")
    print(f"结合模式: {binding_mode}")
    status_text = "✨ 是" if is_glue_substrate else "⚠️ 否"
    print(f"是否为分子胶底物: {status_text}")
    print(f"置信degrees: {confidence:.2f}")
    
    if key_residues_engaged:
        print(f"\n关Key残基参与:")
        for res in key_residues_engaged:
            pos = res['position']
            residue = res['residue']
            if 'contacts_glue_vdw' in res:
                print(f"  {pos}: {residue} - vdW接触Glue (距离: {res['min_distance']} Å) - {res['note']}")
            elif 'contacts_crbn_sidechain' in res:
                print(f"  {pos}: {residue} - {res['note']}")
    
    if geometry_result:
        print(f"\nα-turn 几何: {'✅ 有效' if geometry_result['is_valid_geometry'] else '❌ 无效'}")
        if geometry_result['alpha_turn_distance']:
            print(f"  G₋₄→G₀: {geometry_result['alpha_turn_distance']} Å")
    
    if hbond_result:
        print(f"CRBN氢Key: {hbond_result['total_hbonds']}/3 (标准G-loop: {'是' if hbond_result['is_canonical_gloop'] else '否'})")
    print("=" * 60)
    
    return result


def _distance_3d(coord1, coord2):
    """计算3D距离"""
    return ((coord1[0] - coord2[0])**2 + 
            (coord1[1] - coord2[1])**2 + 
            (coord1[2] - coord2[2])**2)**0.5

# ===========================================================================
# Degron Annotation & Helpers (Moved from pymol_crbn_tools.py)
# ===========================================================================

def _resn3_to_1(resn):
    table = {
        'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
        'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V',
    }
    return table.get(resn.upper(), 'X')

def _seq_of_selection(sel):
    model = cmd.get_model(sel)
    residues = []  # (chain,resi,icode,resn)
    seen = set()
    for a in model.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        if key not in seen:
            seen.add(key)
            residues.append((a.chain, a.resi, getattr(a, 'q', a.icode), a.resn))
    seq = ''.join(_resn3_to_1(r[3]) for r in residues)
    res_keys = [(r[0], r[1], r[2]) for r in residues]
    return seq, res_keys

def _find_beta_hairpins(seq):
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

def degron_annotate(POI_sel: str, csv_path: str = None, auto: bool = True, name: str = "degron"):
    """
    Annotate degron-like regions. If csv provided, expects columns: chain,start,end,label.
    Auto mode marks beta-hairpin-like regions using sequence heuristics.
    """
    poi = f"({POI_sel})"
    regions = []  # start,end,label (0-based indices)

    # CSV regions (1-based resi numbers, chain-agnostic)
    if csv_path:
        if not os.path.isfile(csv_path):
            print(f"[degron_annotate] File not found: {csv_path}")
            return
        with open(csv_path, newline='') as f:
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
        for s, e in _find_beta_hairpins(seq):
            regions.append((s, e, 'beta_hairpin_like'))

    if not regions:
        print("[degron_annotate] No regions to annotate.")
        return

    # Build selections and visuals
    # Map sequence index to (chain,resi,icode)
    _, res_keys = _seq_of_selection(poi)
    colors = {
        'beta_hairpin_like': 'tv_green',
        'degron': 'tv_red',
    }
    for idx, (s, e, label) in enumerate(regions, 1):
        s = max(0, s); e = max(s, e)
        if s >= len(res_keys) or e >= len(res_keys):
            continue
            
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
