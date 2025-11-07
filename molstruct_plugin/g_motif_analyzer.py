# -*- coding: utf-8 -*-
"""
CRBN G-MOTIF（G-loop）识别（支持：理想化 / 内置真实模板 / 自定义选择）
- template_mode: "ideal"（默认） / "builtin" / "selection"
- 当 "builtin" 时，用内置的 PDB+残基段选择取 8×Cα 作为模板；必要时自动 cmd.fetch(async_=0)
- 当 "selection" 时，从 template_sel（选择表达式/sele 名）中取 8×Cα

修复/特性：
- 正确处理插入码（如 100A）、altloc（''/A）、残基排序
- 默认 RMSD_cutoff = 3.5 Å；可选 pos6=Gly（可关闭）
- 控制台输出前若干最小 RMSD，便于调阈值
- 输出 CSV 兼容 highlight_residues：Chain1,Residue1,Chain2,Residue2,Distance,Interaction
"""
from __future__ import print_function
import os, csv, re
from collections import defaultdict
from pymol import cmd

# ====== 1) 内置模板定义（可按需自改）======
BUILTIN_TEMPLATES = {
    "GSPT1 (6H0G A:60-67)": {
        "pdb": "6H0G",
        "selection_fmt": "{obj} and chain A and resi 60+61+62+63+64+65+66+67 and name CA",
    },
    "CK1α (3M51 A:36-43)": {
        "pdb": "3M51",
        "selection_fmt": "{obj} and chain A and resi 36+37+38+39+40+41+42+43 and name CA",
    },
    "VAV1 RT-loop (2MC1 A:95-102)": {
        "pdb": "2MC1",
        "selection_fmt": "{obj} and chain A and resi 95+96+97+98+99+100+101+102 and name CA",
    },
}

# ====== 2) 基础工具：resi 解析/排序、Kabsch RMSD ======
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
    """返回 {chain: [(sort_key, resi_label, resn, (x,y,z))...按 sort_key 排序]}（altloc 取 ''/A）"""
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

# ====== 3) 模板坐标获取：理想化 / 内置 / 选择 ======
def _ideal_beta_hairpin_template():
    left = [(0.0,0.0,0.0),(3.8,0.2,0.0),(7.6,0.1,0.1),(11.4,0.0,0.0)]
    right = [(11.4,4.0,0.2),(7.6,4.2,0.0),(3.8,4.1,-0.1),(0.0,4.0,0.0)]
    return left + right  # 8×3

def _coords_from_selection(sel: str):
    """从选择中抓 8×Cα（按 resi 排序取前 8 个）。"""
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
    for obj in (cmd.get_object_list() or []):
        if code in obj.lower():
            return obj
    return None

def _coords_from_builtin(name: str):
    info = BUILTIN_TEMPLATES.get(name)
    if not info:
        raise ValueError(f"Unknown builtin template: {name}")
    pdb_code = info["pdb"]
    # 若会话中不存在，自动 fetch（注意 async_）
    obj = _find_loaded_object_contains(pdb_code)
    if obj is None:
        try:
            cmd.fetch(pdb_code, async_=0)  # ✅ 修复：用 async_ 避免语法错误
            obj = _find_loaded_object_contains(pdb_code) or pdb_code
        except Exception as e:
            raise RuntimeError(f"Failed to fetch {pdb_code} ({e}). Load the PDB manually or use the 'selection' template.")
    sel = info["selection_fmt"].format(obj=obj)
    return _coords_from_selection(sel)

def _get_template_coords(template_mode: str, template_sel: str | None, template_builtin: str | None):
    mode = (template_mode or "ideal").lower()
    if mode == "selection":
        return _coords_from_selection(template_sel or "")
    if mode == "builtin":
        return _coords_from_builtin(template_builtin or "")
    return _ideal_beta_hairpin_template()

# ====== 4) 主函数 ======
def find_crbn_g_motif(obj_name=None, pdb_file=None,
                       template_mode="ideal", template_sel=None, template_builtin=None,
                       rmsd_cutoff=3.5, out_csv=None, auto_highlight=0, require_gly_pos6=True,
                       topk_debug=10):
    """
    - template_mode: "ideal" / "builtin" / "selection"
    - template_sel: 选择表达式（当 template_mode="selection"）
    - template_builtin: 内置模板名（在 BUILTIN_TEMPLATES 的 key 中任选）
    返回：[(chain, start_resi_label, end_resi_label, seq8, rmsd), ...]
    """
    # 载入对象
    tmp_obj = None
    if pdb_file:
        tmp_obj = "_gmotif_tmp_obj"
        cmd.load(pdb_file, tmp_obj, quiet=1)
        obj = tmp_obj
    else:
        obj = obj_name or (cmd.get_object_list()[0] if cmd.get_object_list() else None)
    if not obj:
        print("[G-MOTIF] No object available; load a structure or provide pdb_file")
        return []

    # 模板坐标
    try:
        tmpl = _get_template_coords(template_mode, template_sel, template_builtin)
    except Exception as e:
        print(f"[G-MOTIF] Template error: {e}; falling back to idealized template.")
        tmpl = _ideal_beta_hairpin_template()

    hits = []
    best_rmsd_pool = []
    by_chain = _collect_ca_by_chain(obj)
    total_windows = 0
    gly_pos6_windows = 0

    for ch, rows in by_chain.items():
        seq = [(r[1], r[2], r[3]) for r in rows]  # (resi_label, resn, xyz)
        if len(seq) < 8:
            continue
        for i in range(0, len(seq) - 7):
            window = seq[i:i+8]
            total_windows += 1
            if require_gly_pos6:
                if window[5][1] not in ("GLY", "G"):
                    continue
                gly_pos6_windows += 1
            P = [w[2] for w in window]
            try:
                rmsd = _kabsch_rmsd(P, tmpl)
            except Exception:
                continue
            best_rmsd_pool.append((rmsd, ch, window))
            if rmsd <= float(rmsd_cutoff):
                resi_s = window[0][0]; resi_e = window[-1][0]
                seq8 = ''.join((aa[:1] if aa else 'X') for aa in [w[1] for w in window])
                hits.append((ch, resi_s, resi_e, seq8, rmsd))

    # 调试输出
    print(f"[G-MOTIF] Total windows: {total_windows}")
    if require_gly_pos6:
        print(f"[G-MOTIF] Windows with pos6=Gly: {gly_pos6_windows}")
    if best_rmsd_pool:
        best_rmsd_pool.sort(key=lambda x: x[0])
        nshow = min(topk_debug, len(best_rmsd_pool))
        print(f"[G-MOTIF] Smallest RMSD (top {nshow}):")
        for k in range(nshow):
            r, ch, window = best_rmsd_pool[k]
            seq8 = ''.join((aa[:1] if aa else 'X') for aa in [w[1] for w in window])
            resi_s = window[0][0]; resi_e = window[-1][0]
            print("  #{:02d} chain={} {:>6s}-{:>6s}  seq={}  RMSD={:.2f} Å".format(
                k+1, ch or '.', str(resi_s), str(resi_e), seq8, r
            ))

    # 写 CSV
    if out_csv is None:
        import tempfile
        fd, out_csv = tempfile.mkstemp(suffix="_gmotif.csv"); os.close(fd)
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        # 使用 G-Motif 专用表头
        w.writerow(["Chain", "Sequence", "Start", "End", "RMSD", "Type"])
        for (ch, s, e, seq8, rmsd) in hits:
            w.writerow([ch, seq8, s, e, f"{rmsd:.2f}", "G-Motif"])

    print(f"[G-MOTIF] Hits: {len(hits)}; output: {out_csv}")

    # 自动高亮
    if auto_highlight and hits:
        try:
            try:
                from .highlight_residues import highlight_csv_residues
            except Exception:
                from highlight_residues import highlight_csv_residues
            highlight_csv_residues(out_csv, obj=obj, show_labels=1, stick_by_element=1)
        except Exception as e:
            print(f"[G-MOTIF] Auto highlight failed: {e}")

    if tmp_obj:
        try: cmd.delete(tmp_obj)
        except Exception: pass

    return hits
