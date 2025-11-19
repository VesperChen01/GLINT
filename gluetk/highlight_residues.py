# -*- coding: utf-8 -*-
"""
highlight_residues.py
从CSV文件高亮显示残基相互作用的功能模块

基于原始的highlight_csv_residues_chain_or_element_sticks.py改进
"""

from __future__ import print_function
import csv
import os
import re
from pymol import cmd

# --- Robust CSV helpers (encoding, delimiter, header normalization) ---
def _read_csv_robust(csv_path):
    """
    Return: list of dict rows (normalized lowercase keys), and fmap.
    Rules:
      - Detect encoding among utf-8, utf-8-sig, gb18030, shift_jis, cp1252.
      - Detect delimiter via csv.Sniffer; fallback to [',',';','\t','|'].
      - Normalize headers: lower(), strip(), remove spaces/underscores.
    """
    # 1) read raw bytes
    with open(csv_path, "rb") as fb:
        raw = fb.read()
    # 2) try encodings
    encodings = ["utf-8", "utf-8-sig", "gb18030", "shift_jis", "cp1252"]
    last_err = None
    text = None
    for enc in encodings:
        try:
            text = raw.decode(enc)
            break
        except Exception as e:
            last_err = e
            text = None
    if text is None:
        raise UnicodeDecodeError("csv", raw, 0, 1, "Cannot decode with tried encodings")
    # 3) detect delimiter
    import csv
    sample = text[:4096]
    delimiter = None
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",",";","\t","|"])
        delimiter = dialect.delimiter
    except Exception:
        for d in [",",";","\t","|"]:
            if d in sample:
                delimiter = d
                break
        if delimiter is None:
            delimiter = ","
    # 4) parse with DictReader
    from io import StringIO
    sio = StringIO(text)
    reader = csv.DictReader(sio, delimiter=delimiter)
    if not reader.fieldnames:
        raise ValueError("CSV file has no header row.")
    # 5) normalize header mapping
    import re as _re
    def _norm(s):
        return _re.sub(r"[_\s]+", "", (s or "").strip().lower())
    fmap = {_norm(k): k for k in reader.fieldnames}
    rows = []
    for r in reader:
        r2 = {}
        for nk, ok in fmap.items():
            r2[nk] = (r.get(ok, "") or "").strip()
        rows.append(r2)
    return rows, fmap

def _pick_value(row, fmap, *candidates):
    """Try multiple header candidates (normalized) and return the first non-empty string."""
    import re as _re
    def _norm(s):
        return _re.sub(r"[_\s]+", "", (s or "").strip().lower())
    for cand in candidates:
        key = _norm(cand)
        if key in row and row[key]:
            return row[key]
    return ""

def _split_residue_tag(tag):
    """Accept 'A:ARG:123' or 'ARG 123' or 'ARG:123' optionally with chain prefix.
       Returns (chain, resn, resi) strings (may be empty)."""
    tag = (tag or "").strip()
    if not tag:
        return "","", ""
    import re as _re
    m = _re.match(r'(?:(?P<chain>[A-Za-z0-9]):)?(?P<resn>[A-Za-z0-9\*]+)[:\s]+(?P<resi>-?\d+)', tag)
    if m:
        d = m.groupdict()
        return d.get("chain",""), d.get("resn",""), d.get("resi","")
    toks = _re.split(r'[:\s]+', tag)
    if len(toks) == 3:
        return toks[0], toks[1], toks[2]
    if len(toks) == 2:
        return "", toks[0], toks[1]
    return "","", tag

def _parse_residue_tag(tag):
    """解析残基标签，支持多种格式"""
    tag = tag.strip()
    # 尝试匹配 "CHAIN:RES NUM" 或 "RES NUM" 格式
    m = re.match(r'(?:(?P<chain>[A-Za-z0-9]):)?(?P<resn>[A-Za-z0-9\*]+)\s+(?P<resi>-?\d+)', tag)
    if m:
        d = m.groupdict()
        return d.get("chain") or "", d["resn"], d["resi"]
    
    # 尝试匹配 "CHAIN:RES:NUM" 格式
    m = re.match(r'(?:(?P<chain>[A-Za-z0-9]):)?(?P<resn>[A-Za-z0-9\*]+):(?P<resi>-?\d+)', tag)
    if m:
        d = m.groupdict()
        return d.get("chain") or "", d["resn"], d["resi"]
    
    # 如果都不匹配，使用简单的分割方法
    toks = re.split(r'[:\s]+', tag)
    chain = ""
    resn = toks[0] if len(toks) >= 1 else ""
    resi = toks[1] if len(toks) >= 2 else ""
    return chain, resn, resi

def _res_sel(obj, chain, resn, resi):
    """构建残基选择表达式"""
    parts = []
    if obj:   parts.append(f"model {obj}")
    if chain: parts.append(f"chain {chain}")
    if resi:  parts.append(f"resi {resi}")
    if resn and resn != "*":
        parts.append(f"resn {resn}")
    return " and ".join(parts) if parts else "all"

def _setup_view(obj, protein_chain, partner_chain, colorA, colorB):
    """设置基本视图和链颜色"""
    cmd.hide("everything", obj)
    if protein_chain:
        cmd.show("cartoon", f"model {obj} and chain {protein_chain}")
        cmd.color(colorA, f"model {obj} and chain {protein_chain}")
    if partner_chain:
        cmd.show("cartoon", f"model {obj} and chain {partner_chain}")
        cmd.color(colorB, f"model {obj} and chain {partner_chain}")
    cmd.set("cartoon_transparency", 0.15)
    cmd.bg_color("white")

def _label_for_residue(resn, resi):
    """生成残基标签，使用一字母氨基酸代码 + 残基号，例如 "T100"""    
    resn = (resn or "").strip().upper()
    resi_str = str(resi).strip()

    # 3-letter 到 1-letter 的氨基酸映射
    aa_map = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
        "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    }

    if resn in aa_map:
        prefix = aa_map[resn]
    elif len(resn) == 1:
        prefix = resn
    else:
        prefix = resn[:1] if resn else ""

    return f"{prefix}{resi_str}" if resi_str else prefix

def _interaction_to_abbr(interaction_text):
    """
    将中文相互作用类型转换为英文缩写
    PyMOL标签不支持中文，需要转换为ASCII字符
    """
    abbr_map = {
        "氢键": "HB",           # Hydrogen Bond
        "盐桥": "SB",           # Salt Bridge
        "疏水相互作用": "HP",   # Hydrophobic
        "π–π 堆积": "Pi-Pi",   # Pi-Pi stacking
        "π–阳离子相互作用": "Pi-Cat",  # Pi-Cation
        "范德华力": "VDW",      # van der Waals
    }

    # 尝试匹配已知的相互作用类型
    interaction_text = interaction_text.strip()
    for cn, en in abbr_map.items():
        if cn in interaction_text:
            return en

    # 如果没有匹配到，返回前10个ASCII字符（避免中文乱码）
    return ''.join(c for c in interaction_text if ord(c) < 128)[:10] or "INT"

def _place_label_pseudoatom(sel_name, label_text, idx):
    """在残基位置放置标签（Times-like 字体，黑色）"""
    # 尝试使用CA原子位置
    for name_sel in [f"({sel_name}) and name CA", f"({sel_name}) and name C1'", f"({sel_name}) and name P"]:
        try:
            coords = cmd.get_atom_coords(name_sel)
            obj_name = f"rlab_{idx}"
            cmd.pseudoatom(obj_name, pos=coords, label=label_text)
            cmd.set("label_size", 16, obj_name)
            cmd.set("label_color", "black", obj_name)
            cmd.set("label_font_id", 5, obj_name)  # 使用接近 Times Roman 的矢量字体
            return
        except Exception:
            pass
    
    # 如果没有CA原子，使用质心
    try:
        model = cmd.get_model(sel_name)
        if model.atom:
            xs = [a.coord[0] for a in model.atom]
            ys = [a.coord[1] for a in model.atom]
            zs = [a.coord[2] for a in model.atom]
            cx, cy, cz = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
            obj_name = f"rlab_{idx}"
            cmd.pseudoatom(obj_name, pos=(cx,cy,cz), label=label_text)
            cmd.set("label_size", 16, obj_name)
            cmd.set("label_color", "black", obj_name)
            cmd.set("label_font_id", 5, obj_name)
    except Exception:
        pass

def _color_sticks_by_element(selection_expr):
    """按元素类型给棍状模型着色"""
    cmd.color("grey70", f"({selection_expr}) and elem C")
    cmd.color("blue",   f"({selection_expr}) and elem N")
    cmd.color("red",    f"({selection_expr}) and elem O")
    cmd.color("yellow", f"({selection_expr}) and elem S")

def highlight_csv_residues(csv_path, obj=None,
                           protein_chain="A", partner_chain="B",
                           colorA="lightblue", colorB="lightorange",
                           show_labels=1, clear_old=1, debug=0, stick_by_element=1,
                           show_interaction_type=0):
    """
    从CSV文件高亮显示残基相互作用

    参数:
        csv_path: CSV文件路径
        obj: PyMOL对象名称，如果为None则使用第一个加载的对象
        protein_chain: 蛋白质链ID (默认"A")
        partner_chain: 伙伴链ID (默认"B")
        colorA: 链A颜色 (默认"lightblue")
        colorB: 链B颜色 (默认"lightorange")
        show_labels: 是否显示标签 (默认1)
        clear_old: 是否清除旧的选择 (默认1)
        debug: 调试模式 (默认0)
        stick_by_element: 是否按元素着色棍状模型 (默认1)
        show_interaction_type: 是否在标签中显示相互作用类型 (默认0，不显示)

    CSV格式要求:
        必须包含列: Chain1, Residue1, Chain2, Residue2
        可选列: Interaction, Distance
    """
    csv_path = os.path.abspath(csv_path)
    if not os.path.exists(csv_path):
        print(f"[highlight_csv_residues] CSV not found: {csv_path}")
        return

    # 获取对象
    if obj is None:
        objs = cmd.get_object_list()
        if not objs:
            print("[highlight_csv_residues] No objects loaded. Please load a structure first.")
            return
        obj = objs[0]
        print(f"[highlight_csv_residues] Using first object: {obj}")

    if obj not in cmd.get_object_list():
        print(f"[highlight_csv_residues] Object '{obj}' not found. Available: {cmd.get_object_list()}")
        return

    # 清除旧的选择和标签
    if clear_old:
        for name in cmd.get_names("objects"):
            if name.startswith("intsel_") or name.startswith("rlab_"):
                cmd.delete(name)

    # 设置基本视图
    _setup_view(obj, protein_chain, partner_chain, colorA, colorB)

    # 读取CSV文件
    rows = []
    try:
        with open(csv_path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            # 创建字段名映射（不区分大小写）
            fmap = {k.strip().lower(): k for k in reader.fieldnames}
            
            def has(key): 
                return key in fmap
            
            def get(r, key): 
                return r[fmap[key]].strip() if key in fmap else ""
            
            for r in reader:
                try:
                    res1 = get(r, "residue1") or get(r, "res1")
                    res2 = get(r, "residue2") or get(r, "res2")
                    ch1 = get(r, "chain1") if has("chain1") else ""
                    ch2 = get(r, "chain2") if has("chain2") else ""
                    interaction = get(r, "interaction") or get(r, "type") or ""
                    
                    if res1 and res2:
                        rows.append((res1, res2, ch1, ch2, interaction))
                except Exception as e:
                    print(f"[highlight_csv_residues] Skipping row: {r}, Error: {e}")
                    continue
    except Exception as e:
        print(f"[highlight_csv_residues] Failed to read CSV: {e}")
        return

    # 处理每一行数据
    built = []
    labeled_residues = set()  # 跟踪已标记的残基，避免重复标签

    for idx, (res1, res2, ch1_csv, ch2_csv, interaction) in enumerate(rows, start=1):
        # 解析残基信息
        c1, n1, i1 = _parse_residue_tag(res1)
        c2, n2, i2 = _parse_residue_tag(res2)

        # 如果解析不到链信息，使用CSV中的链信息
        if not c1 and ch1_csv:
            c1 = ch1_csv
        if not c2 and ch2_csv:
            c2 = ch2_csv

        # 创建选择
        s1 = f"intsel_1_{idx}"
        s2 = f"intsel_2_{idx}"
        sel1 = _res_sel(obj, c1, n1, i1)
        sel2 = _res_sel(obj, c2, n2, i2)

        if debug and idx <= 5:
            print(f"[debug] Residue1 selection -> {sel1}")
            print(f"[debug] Residue2 selection -> {sel2}")

        # 创建选择并显示为棍状模型
        cmd.select(s1, sel1)
        cmd.select(s2, sel2)
        cmd.show("sticks", s1)
        cmd.show("sticks", s2)
        built.extend([s1, s2])

        # 按元素着色（如果启用）
        if stick_by_element:
            _color_sticks_by_element(sel1)
            _color_sticks_by_element(sel2)

        # 添加标签（如果启用）- 只为每个唯一残基创建一次标签
        if show_labels:
            # 为残基1创建唯一标识符
            res1_key = (c1, n1, i1)
            if res1_key not in labeled_residues:
                label1 = _label_for_residue(n1, i1)
                _place_label_pseudoatom(s1, label1, len(labeled_residues) + 1)
                labeled_residues.add(res1_key)

            # 为残基2创建唯一标识符
            res2_key = (c2, n2, i2)
            if res2_key not in labeled_residues:
                label2 = _label_for_residue(n2, i2)
                _place_label_pseudoatom(s2, label2, len(labeled_residues) + 1)
                labeled_residues.add(res2_key)

    # 缩放到显示的残基
    if built:
        cmd.zoom(" or ".join(built), buffer=5.0, complete=1)

    # 强制刷新PyMOL视图
    cmd.refresh()
    cmd.rebuild()

    print(f"[highlight_csv_residues] Highlighted {len(rows)} interaction residue pairs, {len(labeled_residues)} unique residues")

    # 返回统计信息，供GUI使用
    return {
        "pairs": len(rows),
        "unique_residues": len(labeled_residues),
        "selections": built
    }

# 注册命令到PyMOL
cmd.extend("highlight_csv_residues", highlight_csv_residues)

def highlight_gmotif_loops(csv_path, obj=None, color="yellow", show_labels=True, clear_old=True):
    """
    专门用于高亮 G-Motif G-loop 区域

    参数:
        csv_path: G-Motif CSV文件路径
        obj: PyMOL对象名称
        color: G-loop 高亮颜色 (默认"yellow")
        show_labels: 是否显示标签
        clear_old: 是否清除旧的高亮
    """
    csv_path = os.path.abspath(csv_path)
    if not os.path.exists(csv_path):
        print(f"[highlight_gmotif_loops] CSV not found: {csv_path}")
        return

    # 获取对象
    if obj is None:
        objs = cmd.get_object_list()
        if not objs:
            print("[highlight_gmotif_loops] No objects loaded")
            return
        obj = objs[0]

    if obj not in cmd.get_object_list():
        print(f"[highlight_gmotif_loops] Object '{obj}' not found")
        return

    # 清除旧的高亮
    if clear_old:
        for name in cmd.get_names("objects"):
            if name.startswith("gloop_") or name.startswith("glab_"):
                cmd.delete(name)

    # 读取 G-Motif CSV 文件
    loops = []
    try:
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            fmap = {k.strip().lower(): k for k in reader.fieldnames}

            def get(r, key):
                k = key.strip().lower()
                return r[fmap[k]].strip() if k in fmap else ""

            for r in reader:
                try:
                    chain = get(r, "chain")
                    sequence = get(r, "sequence")
                    start = get(r, "start")
                    end = get(r, "end")
                    rmsd = get(r, "rmsd")

                    if chain and start and end:
                        loops.append((chain, sequence, start, end, rmsd))
                except Exception as e:
                    print(f"[highlight_gmotif_loops] Skipping row: {r}, Error: {e}")
                    continue
    except Exception as e:
        print(f"[highlight_gmotif_loops] Failed to read CSV: {e}")
        return

    if not loops:
        print("[highlight_gmotif_loops] No G-loop data found")
        return

    # 高亮每个 G-loop
    for idx, (chain, seq, start, end, rmsd) in enumerate(loops, start=1):
        sel_name = f"gloop_{idx}"
        sel_expr = f"model {obj} and chain {chain} and resi {start}-{end}"

        # 创建选择并显示
        cmd.select(sel_name, sel_expr)
        cmd.show("sticks", sel_name)
        cmd.show("cartoon", sel_name)
        cmd.set("cartoon_thickness", 0.4, sel_name)

        # 按元素着色棍状模型
        _color_sticks_by_element(sel_expr)

        # 高亮 cartoon (使用指定颜色)
        cmd.color(color, f"{sel_expr} and backbone")

        # 添加标签
        if show_labels:
            try:
                # 使用 G-loop 中间残基的 CA 原子位置
                mid_resi = (int(start) + int(end)) // 2
                ca_sel = f"model {obj} and chain {chain} and resi {mid_resi} and name CA"
                coords = cmd.get_atom_coords(ca_sel)
                label_text = f"G-loop {seq} (RMSD={rmsd}Å)"
                obj_name = f"glab_{idx}"
                cmd.pseudoatom(obj_name, pos=coords, label=label_text)
                cmd.set("label_size", 18, obj_name)
                cmd.set("label_color", "black", obj_name)
                cmd.set("label_bg_color", "yellow", obj_name)
                cmd.set("label_bg_transparency", 0.3, obj_name)
            except Exception as e:
                print(f"[highlight_gmotif_loops] Failed to create label: {e}")

    # 缩放到 G-loops
    all_loops = " or ".join([f"gloop_{i}" for i in range(1, len(loops)+1)])
    if all_loops:
        cmd.zoom(all_loops, buffer=8.0, complete=1)

    print(f"[highlight_gmotif_loops] Highlighted {len(loops)} G-loop regions")

cmd.extend("highlight_gmotif_loops", highlight_gmotif_loops)
