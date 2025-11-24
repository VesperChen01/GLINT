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
import math
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
    
    # 自动检测所有链
    all_chains = set()
    model = cmd.get_model(obj)
    for atom in model.atom:
        if atom.chain:
            all_chains.add(atom.chain.strip())
    
    # 设置所有蛋白链的 cartoon 显示
    for chain in all_chains:
        cmd.show("cartoon", f"model {obj} and chain {chain} and polymer.protein")
        
        # 根据链设置颜色
        if chain == protein_chain:
            cmd.color(colorA, f"model {obj} and chain {chain}")
        elif chain == partner_chain:
            cmd.color(colorB, f"model {obj} and chain {chain}")
        else:
            # 其他链使用默认颜色
            cmd.color("gray70", f"model {obj} and chain {chain}")
    
    # 配体显示为 sticks 并使用独特颜色
    ligand_sel = f"model {obj} and organic and not polymer"
    if cmd.count_atoms(ligand_sel) > 0:
        cmd.show("sticks", ligand_sel)
        cmd.color("tv_orange", f"{ligand_sel} and elem C")
        cmd.color("blue", f"{ligand_sel} and elem N")
        cmd.color("red", f"{ligand_sel} and elem O")
        cmd.color("yellow", f"{ligand_sel} and elem S")
    
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

def _color_sticks_by_chain(selection_expr, chain, colorA, colorB):
    """根据链设置碳原子颜色，其他元素保持标准颜色"""
    # 根据链选择颜色
    if chain == "A":
        carbon_color = colorA
    elif chain == "B":
        carbon_color = colorB
    else:
        carbon_color = "grey70"
    
    # 设置碳原子颜色继承链的颜色
    cmd.color(carbon_color, f"({selection_expr}) and elem C")
    # 其他元素保持标准颜色
    cmd.color("blue",   f"({selection_expr}) and elem N")
    cmd.color("red",    f"({selection_expr}) and elem O")
    cmd.color("yellow", f"({selection_expr}) and elem S")

def draw_atom_interaction_lines(obj, chain1, resid1, atom1,
                                chain2, resid2, atom2,
                                interaction_type, idx):
    """
    在两个原子之间绘制相互作用连接线
    
    参数:
        obj: PyMOL对象名称
        chain1, resid1, atom1: 第一个原子的链、残基号、原子名
        chain2, resid2, atom2: 第二个原子的链、残基号、原子名
        interaction_type: 相互作用类型(氢键、盐桥等)
        idx: 索引编号,用于生成唯一的距离对象名称
    
    返回:
        distance_obj_name: 创建的距离对象名称,如果失败则返回None
    """
    # 相互作用类型到颜色的映射
    color_map = {
        "氢键": "yellow",
        "盐桥": "magenta",
        "疏水相互作用": "green",
        "π–π 堆积": "orange",
        "π–阳离子相互作用": "tv_orange",
        "卤素键": "cyan",
        "水桥": "lightblue",
        "金属配位": "purple",
    }
    
    # 获取颜色,默认为灰色
    color = color_map.get(interaction_type, "grey70")
    
    # 构建原子选择表达式
    # 处理原子名称中的特殊字符(如单引号)
    atom1_clean = atom1.replace("'", "\\'") if atom1 else "*"
    atom2_clean = atom2.replace("'", "\\'") if atom2 else "*"
    
    # 构建选择表达式
    parts1 = [f"model {obj}"]
    if chain1:
        parts1.append(f"chain {chain1}")
    if resid1:
        parts1.append(f"resi {resid1}")
    if atom1 and atom1 != "ring" and atom1 != "ring/cation":
        parts1.append(f"name {atom1_clean}")
    sel1 = " and ".join(parts1)
    
    parts2 = [f"model {obj}"]
    if chain2:
        parts2.append(f"chain {chain2}")
    if resid2:
        parts2.append(f"resi {resid2}")
    if atom2 and atom2 != "ring" and atom2 != "ring/cation":
        parts2.append(f"name {atom2_clean}")
    sel2 = " and ".join(parts2)
    
    # 对于π相互作用(ring/ring),使用残基质心
    if atom1 in ["ring", "ring/cation"] or atom2 in ["ring", "ring/cation"]:
        # 使用残基的所有原子质心
        pass  # PyMOL的distance命令会自动处理
    
    # 创建距离对象
    dist_name = f"dist_{idx}"
    
    try:
        # 使用distance命令创建连接线
        cmd.distance(dist_name, sel1, sel2)
        
        # 设置距离对象的显示样式
        cmd.hide("labels", dist_name)  # 隐藏距离标签
        cmd.color(color, dist_name)     # 设置颜色
        cmd.set("dash_width", 2.5, dist_name)  # 设置线条粗细
        cmd.set("dash_gap", 0.3, dist_name)    # 设置虚线间隙(0.3 = 较密集的虚线)
        cmd.set("dash_length", 0.2, dist_name) # 设置虚线长度
        
        # 确保distance对象可见
        cmd.show("dashes", dist_name)
        
        return dist_name
    except Exception as e:
        # 输出错误信息便于调试
        print(f"[draw_atom_interaction_lines] ⚠️  Failed to create distance '{dist_name}':")
        print(f"    Selection 1: {sel1}")
        print(f"    Selection 2: {sel2}")
        print(f"    Error: {e}")
        return None

def highlight_csv_residues(csv_path, obj=None,
                           protein_chain="A", partner_chain="B",
                           colorA="lightblue", colorB="lightorange",
                           show_labels=1, clear_old=1, debug=0, stick_by_element=1,
                           show_interaction_type=0, show_atom_lines=True,
                           show_only_interactions=True):
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
        show_atom_lines: 是否显示原子间相互作用连接线 (默认True)

    CSV格式要求:
        必须包含列: Chain1, Residue1, Chain2, Residue2
        可选列: Interaction, Distance, Ligand_Atom, Protein_Atom (或其他原子列)
    """
    csv_path = os.path.abspath(csv_path)
    if not os.path.exists(csv_path):
        print(f"[highlight_csv_residues] CSV not found: {csv_path}")
        return

    # 获取对象
    if obj is None:
        try:
            objs = cmd.get_names("objects")
        except AttributeError:
            objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
            
        if not objs:
            print("[highlight_csv_residues] No objects loaded. Please load a structure first.")
            return
        obj = objs[0]
        print(f"[highlight_csv_residues] Using first object: {obj}")

    # Check existence
    current_objs = []
    try:
        current_objs = cmd.get_names("objects")
    except AttributeError:
        current_objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []

    if obj not in current_objs:
        print(f"[highlight_csv_residues] Object '{obj}' not found. Available: {current_objs}")
        return

    # 清除旧的选择和标签
    if clear_old:
        for name in cmd.get_names("objects"):
            if name.startswith("intsel_") or name.startswith("rlab_") or name.startswith("dist_"):
                cmd.delete(name)

    # 设置基本视图（如果只显示相互作用，稍后设置）
    if not show_only_interactions:
        _setup_view(obj, protein_chain, partner_chain, colorA, colorB)
    else:
        # 简单的基础设置
        cmd.hide("everything", obj)
        cmd.bg_color("white")

    # 读取CSV文件
    rows = []
    try:
        with open(csv_path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            # 创建字段名映射（不区分大小写）
            fmap = {k.strip().lower(): k for k in reader.fieldnames}
            
            def has(key): 
                key_lower = key.strip().lower()
                return key_lower in fmap
            
            def get(r, key):
                key_lower = key.strip().lower()
                if key_lower in fmap:
                    original_key = fmap[key_lower]
                    return r.get(original_key, "").strip()
                return ""
            
            for r in reader:
                try:
                    # 支持标准格式
                    res1 = get(r, "residue1") or get(r, "res1")
                    res2 = get(r, "residue2") or get(r, "res2")
                    ch1 = get(r, "chain1") if has("chain1") else ""
                    ch2 = get(r, "chain2") if has("chain2") else ""
                    
                    # 支持核酸格式 (Nucleic -> 1, Protein -> 2)
                    if not res1:
                        res1 = get(r, "nucleic_residue")
                        ch1 = get(r, "nucleic_chain")
                    if not res2:
                        res2 = get(r, "protein_residue")
                        ch2 = get(r, "protein_chain")
                        
                    # 支持配体格式 (Ligand -> 1, Protein -> 2)
                    if not res1:
                        res1 = get(r, "ligand_residue")
                        ch1 = get(r, "ligand_chain")
                    
                    interaction = get(r, "interaction") or get(r, "type") or ""
                    
                    # 读取原子信息(如果存在) - 支持大小写
                    # 注意: get 函数会将键转为小写，所以要用小写版本
                    atom1 = get(r, "ligand_atom") or get(r, "atom1") or get(r, "nucleic_atom") or ""
                    atom2 = get(r, "protein_atom") or get(r, "atom2") or ""
                    
                    # 直接检查原始大写版本（如果上面没找到）
                    if not atom1:
                        for key in r.keys():
                            if key in ['Atom1', 'Ligand_Atom']:
                                atom1 = r[key].strip() if r[key] else ""
                                break
                    
                    if not atom2:
                        for key in r.keys():
                            if key in ['Atom2', 'Protein_Atom']:
                                atom2 = r[key].strip() if r[key] else ""
                                break
                    
                    if res1 and res2:
                        rows.append((res1, res2, ch1, ch2, interaction, atom1, atom2))
                except Exception as e:
                    print(f"[highlight_csv_residues] Skipping row: {r}, Error: {e}")
                    continue
    except Exception as e:
        print(f"[highlight_csv_residues] Failed to read CSV: {e}")
        return

    # 处理每一行数据
    built = []
    labeled_residues = set()  # 跟踪已标记的残基，避免重复标签
    distance_objects = []  # 跟踪创建的距离对象

    for idx, (res1, res2, ch1_csv, ch2_csv, interaction, atom1, atom2) in enumerate(rows, start=1):
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

        # 颜色设置：如果是配体则用配体颜色，否则继承链的颜色
        # 检测是否是配体（非标准残基）
        standard_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        is_ligand1 = n1 not in standard_residues
        is_ligand2 = n2 not in standard_residues
        
        # 设置颜色
        if stick_by_element:
            # 配体残基：橙色碳原子
            if is_ligand1:
                cmd.color("tv_orange", f"({sel1}) and elem C")
                cmd.color("blue", f"({sel1}) and elem N")
                cmd.color("red", f"({sel1}) and elem O")
                cmd.color("yellow", f"({sel1}) and elem S")
            else:
                # 蛋白残基：继承链的颜色的碳原子
                _color_sticks_by_chain(sel1, c1, colorA, colorB)
            
            if is_ligand2:
                cmd.color("tv_orange", f"({sel2}) and elem C")
                cmd.color("blue", f"({sel2}) and elem N")
                cmd.color("red", f"({sel2}) and elem O")
                cmd.color("yellow", f"({sel2}) and elem S")
            else:
                # 蛋白残基：继承链的颜色的碳原子
                _color_sticks_by_chain(sel2, c2, colorA, colorB)

        # 绘制原子级相互作用连接线（如果启用且有原子信息）
        if show_atom_lines:
            if atom1 and atom2:
                # 从残基标签中提取残基号
                # res1/res2格式: "LIG 1" 或 "ARG 123"
                resid1 = i1  # 已经从_parse_residue_tag解析出来
                resid2 = i2
                
                if debug:
                    print(f"[debug] Drawing atom line {idx}: {c1}:{resid1}:{atom1} <-> {c2}:{resid2}:{atom2} ({interaction})")
                
                dist_obj = draw_atom_interaction_lines(
                    obj, c1, resid1, atom1,
                    c2, resid2, atom2,
                    interaction, idx
                )
                if dist_obj:
                    distance_objects.append(dist_obj)
                    if debug:
                        print(f"[debug]   ✓ Created distance object: {dist_obj}")
            elif debug and idx <= 5:
                print(f"[debug] No atom info for row {idx}: atom1='{atom1}', atom2='{atom2}'")

        # 添加标签（如果启用）- 只为每个唯一残基创建一次标签
        if show_labels:
            # 为残基1创建唯一标识符
            res1_key = (c1, n1, i1)
            if res1_key not in labeled_residues:
                # 生成标签文本
                label1 = _label_for_residue(n1, i1)
                # 如果需要显示相互作用类型
                if show_interaction_type and interaction:
                    abbr = _interaction_to_abbr(interaction)
                    label1 = f"{label1}({abbr})"
                _place_label_pseudoatom(s1, label1, len(labeled_residues) + 1)
                labeled_residues.add(res1_key)

            # 为残基2创建唯一标识符
            res2_key = (c2, n2, i2)
            if res2_key not in labeled_residues:
                label2 = _label_for_residue(n2, i2)
                if show_interaction_type and interaction:
                    abbr = _interaction_to_abbr(interaction)
                    label2 = f"{label2}({abbr})"
                _place_label_pseudoatom(s2, label2, len(labeled_residues) + 1)
                labeled_residues.add(res2_key)

    # 如果启用了只显示相互作用，现在设置蛋白背景
    if show_only_interactions and built:
        # 显示蛋白链为半透明 cartoon 背景
        cmd.show("cartoon", f"{obj} and polymer.protein")
        cmd.set("cartoon_transparency", 0.6, obj)
        cmd.color("gray80", f"{obj} and polymer.protein")
        
        # 强调显示相互作用的残基
        for sel in built:
            # 相互作用残基显示为不透明
            cmd.set("stick_transparency", 0.0, sel)
    
    # 缩放到显示的残基
    if built:
        cmd.zoom(" or ".join(built), buffer=5.0, complete=1)

    # 确保所有distance对象可见
    if show_atom_lines and distance_objects:
        cmd.show("dashes")  # 全局显示所有虚线对象
        print(f"[highlight_csv_residues] 💡 Tip: If you don't see dashes, run: show dashes")

    # 强制刷新PyMOL视图
    cmd.refresh()
    cmd.rebuild()

    if show_atom_lines and distance_objects:
        print(f"[highlight_csv_residues] Highlighted {len(rows)} interaction residue pairs, {len(labeled_residues)} unique residues, {len(distance_objects)} atom-level connections")
    else:
        print(f"[highlight_csv_residues] Highlighted {len(rows)} interaction residue pairs, {len(labeled_residues)} unique residues")

    # 返回统计信息，供GUI使用
    return {
        "pairs": len(rows),
        "unique_residues": len(labeled_residues),
        "selections": built,
        "distance_objects": distance_objects
    }

# 注册命令到PyMOL
cmd.extend("highlight_csv_residues", highlight_csv_residues)

# ===========================================================================
# Score Coloring & B-factor Utilities (Moved from pymol_crbn_tools.py)
# ===========================================================================

def normalize_values(v):
    """Normalize a list of values to [0, 1]"""
    if not v:
        return (0.0, 1.0)
    vmin, vmax = min(v), max(v)
    if math.isclose(vmin, vmax):
        vmax = vmin + 1.0
    return vmin, vmax

def set_b_factors(sel, resi_to_value):
    """
    Set B-factors for residues based on a dictionary {(chain, resi, icode): value}
    resi key: (chain, resi, icode)
    """
    # Reset B-factors
    cmd.alter(sel, "b=b", space={"b": 0.0})
    cmd.iterate_state(1, sel, "b=b", space={})  # Force update
    
    # Helper for updating atoms
    model = cmd.get_model(sel)
    for a in model.atom:
        key = (a.chain, a.resi, getattr(a, 'q', a.icode))
        if key in resi_to_value:
            icode_part = getattr(a, 'q', a.icode)
            icode_str = icode_part if icode_part else '""'
            cmd.alter(f"{sel} and chain {a.chain} and resi {a.resi} and icode {icode_str}", 
                     f"b={float(resi_to_value[key])}")
    cmd.rebuild()

def color_by_b(sel, palette="blue_white_red", min_val=None, max_val=None, ramp_name=""):
    """
    Color selection by B-factor spectrum
    """
    if min_val is None or max_val is None:
        vals = []
        model = cmd.get_model(sel)
        for a in model.atom:
            vals.append(a.b)
        vmin, vmax = normalize_values(vals)
    else:
        vmin, vmax = min_val, max_val
        
    cmd.spectrum("b", palette, selection=sel, minimum=vmin, maximum=vmax)
    
    if ramp_name:
        cmd.ramp_new(ramp_name, sel, [vmin, (vmin+vmax)/2.0, vmax], ["blue", "white", "red"])

def score_color(POI_sel, scores_csv, col_chain="chain", col_resi="resi", col_score="score", 
                vmin=None, vmax=None):
    """
    Color POI by scores from CSV (chain,resi,score). 
    Stores score into b-factor and colors by spectrum.
    """
    poi = f"({POI_sel})"
    if not os.path.isfile(scores_csv):
        print(f"[score_color] File not found: {scores_csv}")
        return
        
    mapping = {}
    try:
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
    except Exception as e:
        print(f"[score_color] Failed to read CSV: {e}")
        return

    if not mapping:
        print("[score_color] No valid scores parsed.")
        return
        
    set_b_factors(poi, mapping)
    color_by_b(poi, palette="blue_white_red", min_val=vmin, max_val=vmax, ramp_name="score_ramp")
    cmd.show("surface", poi)

cmd.extend("score_color", score_color)

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
        try:
            objs = cmd.get_names("objects")
        except AttributeError:
            objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
            
        if not objs:
            print("[highlight_gmotif_loops] No objects loaded")
            return
        obj = objs[0]

    # Check existence
    current_objs = []
    try:
        current_objs = cmd.get_names("objects")
    except AttributeError:
        current_objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
        
    if obj not in current_objs:
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
