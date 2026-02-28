# -*- coding: utf-8 -*-
"""
highlight_residues.py
从CSVFile高亮Display残基相互作用的功能Module

基于原始的highlight_csv_residues_chain_or_element_sticks.py改进
"""

from __future__ import print_function
import csv
import os
import re
import math
from pymol import cmd

# Import unified color scheme
try:
    from .color_scheme import INTERACTION_COLORS_HEX
except ImportError:
    from color_scheme import INTERACTION_COLORS_HEX

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
    """解析残基Label，支持多种格式"""
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
    
    # 如果都不匹配，using简单的分割Method
    toks = re.split(r'[:\s]+', tag)
    chain = ""
    resn = toks[0] if len(toks) >= 1 else ""
    resi = toks[1] if len(toks) >= 2 else ""
    return chain, resn, resi

def _res_sel(obj, chain, resn, resi):
    """构建残基Select表达式"""
    parts = []
    if obj:   parts.append(f"model {obj}")
    if chain: parts.append(f"chain {chain}")
    if resi:  parts.append(f"resi {resi}")
    if resn and resn != "*":
        parts.append(f"resn {resn}")
    return " and ".join(parts) if parts else "all"

def _setup_view(obj, protein_chain, partner_chain, colorA, colorB):
    """Settings基本视图和链颜色"""
    cmd.hide("everything", obj)
    
    # 自动检测所有链
    all_chains = set()
    model = cmd.get_model(obj)
    for atom in model.atom:
        if atom.chain:
            all_chains.add(atom.chain.strip())
    
    # Settings所有蛋白链的 cartoon Display
    for chain in all_chains:
        cmd.show("cartoon", f"model {obj} and chain {chain} and polymer.protein")
        
        # 根据链Settings颜色
        if chain == protein_chain:
            cmd.color(colorA, f"model {obj} and chain {chain}")
        elif chain == partner_chain:
            cmd.color(colorB, f"model {obj} and chain {chain}")
        else:
            # 其他链using默认颜色
            cmd.color("gray70", f"model {obj} and chain {chain}")
    
    # 配体Display为 sticks 并using独特颜色
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
    """生成残基Label，using一字母氨基酸代码 + 残基号，例如 "T100"""    
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
    将中文相互作用Type转换为英文缩写
    PyMOLLabel不支持中文，需要转换为ASCII字符
    """
    abbr_map = {
        "Hydrogen Bond": "HB",
        "Salt Bridge": "SB",
        "Hydrophobic": "HP",
        "Pi-Pi Stacking": "Pi-Pi",
        "Pi-Cation": "Pi-Cat",
        "van der Waals": "VDW",
    }

    # 尝试匹配已知的相互作用Type
    interaction_text = interaction_text.strip()
    for cn, en in abbr_map.items():
        if cn in interaction_text:
            return en

    # 如果没有匹配到，Return前10个ASCII字符（避免中文乱码）
    return ''.join(c for c in interaction_text if ord(c) < 128)[:10] or "INT"

def _place_label_pseudoatom(sel_name, label_text, idx):
    """在残基位置放置Label（Times-like 字体，黑色）"""
    # 尝试usingCA原子位置
    for name_sel in [f"({sel_name}) and name CA", f"({sel_name}) and name C1'", f"({sel_name}) and name P"]:
        try:
            coords = cmd.get_atom_coords(name_sel)
            obj_name = f"rlab_{idx}"
            cmd.pseudoatom(obj_name, pos=coords, label=label_text)
            cmd.set("label_size", 16, obj_name)
            cmd.set("label_color", "black", obj_name)
            cmd.set("label_font_id", 5, obj_name)  # using接近 Times Roman 的矢量字体
            return
        except Exception:
            pass
    
    # 如果没有CA原子，using质心
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
    """按元素Type给棍状模型着色"""
    cmd.color("grey70", f"({selection_expr}) and elem C")
    cmd.color("blue",   f"({selection_expr}) and elem N")
    cmd.color("red",    f"({selection_expr}) and elem O")
    cmd.color("yellow", f"({selection_expr}) and elem S")

def _color_sticks_by_chain(selection_expr, chain, colorA, colorB):
    """根据链Settings碳原子颜色，其他元素保持标准颜色"""
    # 根据链Select颜色
    if chain == "A":
        carbon_color = colorA
    elif chain == "B":
        carbon_color = colorB
    else:
        carbon_color = "grey70"
    
    # Settings碳原子颜色继承链的颜色
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
    
    Parameters:
        obj: PyMOL对象Name
        chain1, resid1, atom1: First原子的链、残基号、原子名
        chain2, resid2, atom2: 第二个原子的链、残基号、原子名
        interaction_type: Interaction type (Hydrogen Bond, Salt Bridge, etc.)
        idx: 索引Number,用于生成唯一的距离对象Name
    
    Return:
        distance_obj_name: Create的距离对象Name,如果Failed则ReturnNone
    """
    # using统一配色方案
    # 将 hex 颜色转换为 PyMOL 颜色Name（或直接using相近的 PyMOL 内置颜色）
    color_map = {
        "Hydrogen Bond": "marine",           # Blue (#2196F3)
        "Salt Bridge": "tv_red",          # Orange-red (#FF5722)
        "Hydrophobic": "green",    # Green (#4CAF50)
        "Pi-Pi Stacking": "purple",      # Purple (#9C27B0)
        "Pi-Cation": "magenta",  # Pink (#E91E63)
        "Halogen Bond": "tv_orange",     # Orange (#FF9800)
        "Water Bridge": "cyan",            # Cyan (#00BCD4)
        "Metal Coordination": "violet",      # Deep purple (#673AB7)
    }
    
    # 获取颜色,默认为灰色
    color = color_map.get(interaction_type, "grey70")
    
    # 构建原子Select表达式
    # 处理原子Name中的特殊字符(如单引号)
    atom1_clean = atom1.replace("'", "\\'") if atom1 else "*"
    atom2_clean = atom2.replace("'", "\\'") if atom2 else "*"
    
    # 构建Select表达式
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
    
    # 对于π相互作用(ring/ring),using残基质心
    if atom1 in ["ring", "ring/cation"] or atom2 in ["ring", "ring/cation"]:
        # using残基的所有原子质心
        pass  # PyMOL的distance命令会自动处理
    
    # Create距离对象
    dist_name = f"dist_{idx}"
    
    # Determine cutoff based on interaction type
    # Since interactions are already filtered by the analyzer, we use a generous cutoff
    # to ensure they are drawn in PyMOL regardless of minor coordinate discrepancies.
    cutoff = 10.0 
    
    # if "氢Key" in interaction_type or "Hydrogen" in interaction_type:
    #     cutoff = 3.8
    # elif "盐桥" in interaction_type or "Salt" in interaction_type:
    #     cutoff = 5.0
    # elif "疏水" in interaction_type or "Hydrophobic" in interaction_type:
    #     cutoff = 5.5
    # elif "Pi" in interaction_type or "π" in interaction_type:
    #     cutoff = 6.5
    # elif "金属" in interaction_type or "Metal" in interaction_type:
    #     cutoff = 4.0
    
    try:
        # usingdistance命令Create连接线
        # using cutoff ParametersFilter掉距离过远的Error连接（例如 > 3.8A 的氢Key）
        cmd.distance(dist_name, sel1, sel2, cutoff=cutoff)
        
        # Settings距离对象的Display样式
        cmd.hide("labels", dist_name)  # Hide距离Label
        cmd.color(color, dist_name)     # Settings颜色
        cmd.set("dash_width", 2.5, dist_name)  # Settings线条粗细
        cmd.set("dash_gap", 0.3, dist_name)    # Settings虚线间隙(0.3 = 较密集的虚线)
        cmd.set("dash_length", 0.2, dist_name) # Settings虚线长degrees
        
        # 确保distance对象可见
        cmd.show("dashes", dist_name)
        
        return dist_name
    except Exception as e:
        # 输出ErrorInformation便于调试
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
    从CSVFile高亮Display残基相互作用

    Parameters:
        csv_path: CSVFilePath
        obj: PyMOL对象Name，如果为None则usingFirstLoad的对象
        protein_chain: 蛋白质链ID (默认"A")
        partner_chain: 伙伴链ID (默认"B")
        colorA: 链A颜色 (默认"lightblue")
        colorB: 链B颜色 (默认"lightorange")
        show_labels: 是否DisplayLabel (默认1)
        clear_old: 是否清除旧的Select (默认1)
        debug: 调试模式 (默认0)
        stick_by_element: 是否按元素着色棍状模型 (默认1)
        show_interaction_type: 是否在Label中Display相互作用Type (默认0，不Display)
        show_atom_lines: 是否Display原子间相互作用连接线 (默认True)

    CSV格式要求:
        必须Package含列: Chain1, Residue1, Chain2, Residue2
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

    # 清除旧的Select和Label
    if clear_old:
        for name in cmd.get_names("objects"):
            if name.startswith("intsel_") or name.startswith("rlab_") or name.startswith("dist_"):
                cmd.delete(name)

    # Settings基本视图（如果只Display相互作用，稍后Settings）
    if not show_only_interactions:
        _setup_view(obj, protein_chain, partner_chain, colorA, colorB)
    else:
        # 简单的基础Settings
        cmd.hide("everything", obj)
        cmd.bg_color("white")

    # 读取CSVFile
    rows = []
    try:
        with open(csv_path, "r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            # Create字段名映射（不区分Size写）
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
                    
                    # 读取Atom information(如果存在) - 支持Size写
                    # 注意: get Function会将Key转为小写，所以要用小写Version
                    atom1 = get(r, "ligand_atom") or get(r, "atom1") or get(r, "nucleic_atom") or ""
                    atom2 = get(r, "protein_atom") or get(r, "atom2") or ""
                    
                    # 直接检查原始大写Version（如果上面没找到）
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
    labeled_residues = set()  # 跟踪已标记的残基，避免重复Label
    distance_objects = []  # 跟踪Create的距离对象

    for idx, (res1, res2, ch1_csv, ch2_csv, interaction, atom1, atom2) in enumerate(rows, start=1):
        # 解析残基Information
        c1, n1, i1 = _parse_residue_tag(res1)
        c2, n2, i2 = _parse_residue_tag(res2)

        # 如果解析不到Chain information，usingCSV中的Chain information
        if not c1 and ch1_csv:
            c1 = ch1_csv
        if not c2 and ch2_csv:
            c2 = ch2_csv

        # CreateSelect
        s1 = f"intsel_1_{idx}"
        s2 = f"intsel_2_{idx}"
        sel1 = _res_sel(obj, c1, n1, i1)
        sel2 = _res_sel(obj, c2, n2, i2)

        if debug and idx <= 5:
            print(f"[debug] Residue1 selection -> {sel1}")
            print(f"[debug] Residue2 selection -> {sel2}")

        # CreateSelect并Display为棍状模型
        cmd.select(s1, sel1)
        cmd.select(s2, sel2)
        cmd.show("sticks", s1)
        cmd.show("sticks", s2)
        built.extend([s1, s2])

        # 颜色Settings：如果是配体则用配体颜色，否则继承链的颜色
        # 检测是否是配体（非标准残基）
        standard_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        is_ligand1 = n1 not in standard_residues
        is_ligand2 = n2 not in standard_residues
        
        # Settings颜色
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

        # 绘制原子级相互作用连接线（如果Enable且有Atom information）
        if show_atom_lines:
            if atom1 and atom2:
                # 从残基Label中提取残基号
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

        # AddLabel（如果Enable）- 只为每个唯一残基Create一次Label
        if show_labels:
            # 为残基1Create唯一ID符
            res1_key = (c1, n1, i1)
            if res1_key not in labeled_residues:
                # 生成Label文本
                label1 = _label_for_residue(n1, i1)
                # 如果需要Display相互作用Type
                if show_interaction_type and interaction:
                    abbr = _interaction_to_abbr(interaction)
                    label1 = f"{label1}({abbr})"
                _place_label_pseudoatom(s1, label1, len(labeled_residues) + 1)
                labeled_residues.add(res1_key)

            # 为残基2Create唯一ID符
            res2_key = (c2, n2, i2)
            if res2_key not in labeled_residues:
                label2 = _label_for_residue(n2, i2)
                if show_interaction_type and interaction:
                    abbr = _interaction_to_abbr(interaction)
                    label2 = f"{label2}({abbr})"
                _place_label_pseudoatom(s2, label2, len(labeled_residues) + 1)
                labeled_residues.add(res2_key)

    # 如果Enable了只Display相互作用，现在Settings蛋白背景
    if show_only_interactions and built:
        # Display蛋白链为半透明 cartoon 背景
        cmd.show("cartoon", f"{obj} and polymer.protein")
        cmd.set("cartoon_transparency", 0.6, obj)
        cmd.color("gray80", f"{obj} and polymer.protein")
        
        # 强调Display相互作用的残基
        for sel in built:
            # 相互作用残基Display为不透明
            cmd.set("stick_transparency", 0.0, sel)
    
    # 缩放到Display的残基
    if built:
        cmd.zoom(" or ".join(built), buffer=5.0, complete=1)

    # 确保所有distance对象可见
    if show_atom_lines and distance_objects:
        cmd.show("dashes")  # 全局Display所有虚线对象
        print(f"[highlight_csv_residues] 💡 Tip: If you don't see dashes, run: show dashes")

    # 强制RefreshPyMOL视图
    cmd.refresh()
    cmd.rebuild()

    if show_atom_lines and distance_objects:
        print(f"[highlight_csv_residues] Highlighted {len(rows)} interaction residue pairs, {len(labeled_residues)} unique residues, {len(distance_objects)} atom-level connections")
    else:
        print(f"[highlight_csv_residues] Highlighted {len(rows)} interaction residue pairs, {len(labeled_residues)} unique residues")

    # Return统计Information，供GUIusing
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

    Parameters:
        csv_path: G-Motif CSVFilePath
        obj: PyMOL对象Name
        color: G-loop 高亮颜色 (默认"yellow")
        show_labels: 是否DisplayLabel
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

    # 读取 G-Motif CSV File
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
        
        # 处理残基范围Select - 支持带插入码的残基（如 100A）
        # PyMOL 的 resi Select器对于范围可能不支持插入码，所以我们需要特殊处理
        try:
            start_int = int(start)
            end_int = int(end)
            # 如果 start 和 end 都是纯数字，using范围Select
            sel_expr = f"model {obj} and chain {chain} and resi {start_int}-{end_int}"
        except ValueError:
            # 如果Package含插入码，using逐个残基Select
            # 尝试解析残基号和插入码
            import re
            start_match = re.match(r'(-?\d+)([A-Za-z]?)', str(start))
            end_match = re.match(r'(-?\d+)([A-Za-z]?)', str(end))
            
            if start_match and end_match:
                start_num = int(start_match.group(1))
                end_num = int(end_match.group(1))
                # 生成残基列表（简化处理：只using数字范围）
                resi_list = '+'.join(str(i) for i in range(start_num, end_num + 1))
                sel_expr = f"model {obj} and chain {chain} and resi {resi_list}"
            else:
                # 回退到原始方式
                sel_expr = f"model {obj} and chain {chain} and resi {start}-{end}"

        # 检查是否有原子
        if cmd.count_atoms(sel_expr) == 0:
            print(f"[highlight_gmotif_loops] Warning: No atoms selected for G-loop {idx} ({chain}:{start}-{end})")
            continue

        # Create新对象
        cmd.create(sel_name, sel_expr)
        
        # SettingsDisplay
        cmd.hide("everything", sel_name)
        cmd.show("sticks", sel_name)
        cmd.show("cartoon", sel_name)
        cmd.set("cartoon_thickness", 0.4, sel_name)

        # 按元素着色棍状模型 (对新对象)
        _color_sticks_by_element(sel_name)

        # 高亮 cartoon (using指定颜色)
        cmd.color(color, f"{sel_name} and backbone")

        # AddLabel
        if show_labels:
            try:
                # using G-loop 中间残基的 CA 原子位置
                # 首先尝试从Select中获取 CA 原子
                ca_sel = f"{sel_name} and name CA"
                ca_count = cmd.count_atoms(ca_sel)
                
                if ca_count > 0:
                    # 获取所有 CA 原子坐标，取中间一个
                    model = cmd.get_model(ca_sel)
                    if model.atom:
                        mid_idx = len(model.atom) // 2
                        mid_atom = model.atom[mid_idx]
                        coords = (mid_atom.coord[0], mid_atom.coord[1], mid_atom.coord[2])
                        
                        label_text = f"G-loop {seq} (RMSD={rmsd}A)"
                        label_obj_name = f"glab_{idx}"
                        cmd.pseudoatom(label_obj_name, pos=coords, label=label_text)
                        cmd.set("label_size", 18, label_obj_name)
                        cmd.set("label_color", "black", label_obj_name)
                        cmd.set("label_bg_color", "yellow", label_obj_name)
                        cmd.set("label_bg_transparency", 0.3, label_obj_name)
                else:
                    print(f"[highlight_gmotif_loops] No CA atoms found for label in G-loop {idx}")
            except Exception as e:
                print(f"[highlight_gmotif_loops] Failed to create label for G-loop {idx}: {e}")

    # 缩放到 G-loops
    all_loops = " or ".join([f"gloop_{i}" for i in range(1, len(loops)+1)])
    if all_loops:
        cmd.zoom(all_loops, buffer=8.0, complete=1)

    print(f"[highlight_gmotif_loops] Highlighted {len(loops)} G-loop regions")

cmd.extend("highlight_gmotif_loops", highlight_gmotif_loops)


def highlight_gloop_surface(obj=None, chain=None, start_resi=None, end_resi=None,
                            surface_color="yellow", surface_transparency=0.3,
                            show_cartoon=True, cartoon_color="tv_yellow",
                            selection_name=None, clear_old=True):
    """
    高亮 G-loop 区域的分子表面
    
    Parameters:
        obj: PyMOL 对象Name
        chain: 链 ID
        start_resi: 起始残基Number
        end_resi: 结束残基Number
        surface_color: 表面颜色 (默认 "yellow")
        surface_transparency: 表面透明degrees (默认 0.3)
        show_cartoon: 是否同时Display cartoon (默认 True)
        cartoon_color: cartoon 颜色 (默认 "tv_yellow")
        selection_name: 自定义SelectName (默认自动生成)
        clear_old: 是否清除旧的 G-loop 表面 (默认 True)
    
    Return:
        dict: {
            'selection_name': str,  # Create的SelectName
            'atom_count': int,      # 选中的原子数
            'surface_area': float,  # 表面积 (Å²)
        }
    """
    # 获取对象
    if obj is None:
        try:
            objs = cmd.get_names("objects")
        except AttributeError:
            objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
        if not objs:
            print("[highlight_gloop_surface] No objects loaded")
            return None
        obj = objs[0]
    
    # 验证Parameters
    if chain is None or start_resi is None or end_resi is None:
        print("[highlight_gloop_surface] Error: chain, start_resi, end_resi are required")
        return None
    
    # 清除旧的 G-loop 表面对象
    if clear_old:
        for name in cmd.get_names("objects"):
            if name.startswith("gloop_surf_"):
                cmd.delete(name)
    
    # CreateSelectName (作为新的对象名)
    if selection_name is None:
        selection_name = f"gloop_surf_{chain}_{start_resi}_{end_resi}"
    
    # 构建源Select表达式
    try:
        start_int = int(start_resi)
        end_int = int(end_resi)
        sel_expr = f"model {obj} and chain {chain} and resi {start_int}-{end_int}"
    except ValueError:
        # 处理带插入码的残基
        resi_list = '+'.join(str(i) for i in range(int(start_resi.rstrip('ABCDEFGHIJ')),
                                                    int(end_resi.rstrip('ABCDEFGHIJ')) + 1))
        sel_expr = f"model {obj} and chain {chain} and resi {resi_list}"
    
    # 检查是否有原子
    if cmd.count_atoms(sel_expr) == 0:
        print(f"[highlight_gloop_surface] Warning: No atoms selected for {chain}:{start_resi}-{end_resi}")
        return None

    # Create新对象 (extract/create)
    cmd.create(selection_name, sel_expr)
    
    # 获取新对象的原子数
    atom_count = cmd.count_atoms(selection_name)
    
    # Settings新对象的Display方式
    cmd.hide("everything", selection_name) # Hide所有默认表示
    cmd.show("surface", selection_name)     # 只Display表面
    
    # Settings颜色和透明degrees
    cmd.color(surface_color, selection_name)
    cmd.set("surface_transparency", surface_transparency, selection_name)
    
    # 可选：Display cartoon (在新对象上)
    if show_cartoon:
        cmd.show("cartoon", selection_name)
        cmd.color(cartoon_color, selection_name)
    
    # 计算表面积
    try:
        surface_area = cmd.get_area(selection_name, state=1)
    except Exception:
        surface_area = None
    
    # 缩放到Select区域
    cmd.zoom(selection_name, buffer=5.0, complete=1)
    
    print(f"[highlight_gloop_surface] ✅ Highlighted G-loop surface: {chain}:{start_resi}-{end_resi}")
    print(f"    Atoms: {atom_count}, Surface area: {surface_area:.1f} Å²" if surface_area else f"    Atoms: {atom_count}")
    
    return {
        'selection_name': selection_name,
        'atom_count': atom_count,
        'surface_area': surface_area
    }

cmd.extend("highlight_gloop_surface", highlight_gloop_surface)


def get_gloop_coordinates(obj=None, chain=None, start_resi=None, end_resi=None,
                          atom_types=None, output_csv=None):
    """
    获取 G-loop 区域的原子坐标
    
    Parameters:
        obj: PyMOL 对象Name
        chain: 链 ID
        start_resi: 起始残基Number
        end_resi: 结束残基Number
        atom_types: 要提取的原子Type列表 (默认 ["CA"] 只提取 Cα)
                    可选: ["CA"], ["CA", "CB"], ["all"], None (等同于 ["CA"])
        output_csv: 输出 CSV FilePath (可选)
    
    Return:
        dict: {
            'coordinates': [
                {
                    'chain': str,
                    'resi': str,
                    'resn': str,
                    'atom': str,
                    'x': float,
                    'y': float,
                    'z': float,
                    'b_factor': float
                },
                ...
            ],
            'centroid': (x, y, z),  # 质心坐标
            'ca_coords': [(x, y, z), ...],  # Cα 坐标列表 (用于 RMSD 计算)
            'csv_path': str or None
        }
    """
    # 获取对象
    if obj is None:
        try:
            objs = cmd.get_names("objects")
        except AttributeError:
            objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []
        if not objs:
            print("[get_gloop_coordinates] No objects loaded")
            return None
        obj = objs[0]
    
    # 验证Parameters
    if chain is None or start_resi is None or end_resi is None:
        print("[get_gloop_coordinates] Error: chain, start_resi, end_resi are required")
        return None
    
    # 默认只提取 Cα
    if atom_types is None:
        atom_types = ["CA"]
    
    # 构建Select表达式
    try:
        start_int = int(start_resi)
        end_int = int(end_resi)
        sel_expr = f"model {obj} and chain {chain} and resi {start_int}-{end_int}"
    except ValueError:
        resi_list = '+'.join(str(i) for i in range(int(start_resi.rstrip('ABCDEFGHIJ')),
                                                    int(end_resi.rstrip('ABCDEFGHIJ')) + 1))
        sel_expr = f"model {obj} and chain {chain} and resi {resi_list}"
    
    # 如果不是 "all"，Add原子TypeFilter
    if atom_types != ["all"]:
        atom_filter = "+".join(atom_types)
        sel_expr = f"({sel_expr}) and name {atom_filter}"
    
    # 获取原子模型
    model = cmd.get_model(sel_expr)
    
    if len(model.atom) == 0:
        print(f"[get_gloop_coordinates] Warning: No atoms found for {chain}:{start_resi}-{end_resi}")
        return None
    
    # 提取坐标
    coordinates = []
    ca_coords = []
    sum_x, sum_y, sum_z = 0.0, 0.0, 0.0
    
    for atom in model.atom:
        coord_dict = {
            'chain': atom.chain,
            'resi': atom.resi,
            'resn': atom.resn,
            'atom': atom.name,
            'x': round(atom.coord[0], 3),
            'y': round(atom.coord[1], 3),
            'z': round(atom.coord[2], 3),
            'b_factor': round(atom.b, 2)
        }
        coordinates.append(coord_dict)
        
        # 累加用于Calculate center of mass
        sum_x += atom.coord[0]
        sum_y += atom.coord[1]
        sum_z += atom.coord[2]
        
        # 收集 Cα 坐标
        if atom.name.strip().upper() == 'CA':
            ca_coords.append((atom.coord[0], atom.coord[1], atom.coord[2]))
    
    # Calculate center of mass
    n_atoms = len(coordinates)
    centroid = (
        round(sum_x / n_atoms, 3),
        round(sum_y / n_atoms, 3),
        round(sum_z / n_atoms, 3)
    ) if n_atoms > 0 else (0.0, 0.0, 0.0)
    
    # 输出到 CSV
    csv_path = None
    if output_csv:
        import csv as csv_module
        csv_path = os.path.abspath(output_csv)
        with open(csv_path, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv_module.writer(f)
            writer.writerow(['Chain', 'Resi', 'Resn', 'Atom', 'X', 'Y', 'Z', 'B_factor'])
            for coord in coordinates:
                writer.writerow([
                    coord['chain'], coord['resi'], coord['resn'], coord['atom'],
                    coord['x'], coord['y'], coord['z'], coord['b_factor']
                ])
            # Add质心行
            writer.writerow(['', '', 'CENTROID', '', centroid[0], centroid[1], centroid[2], ''])
        print(f"[get_gloop_coordinates] Coordinates saved to: {csv_path}")
    
    print(f"[get_gloop_coordinates] ✅ Extracted {len(coordinates)} atoms, {len(ca_coords)} Cα atoms")
    print(f"    Centroid: ({centroid[0]:.2f}, {centroid[1]:.2f}, {centroid[2]:.2f})")
    
    return {
        'coordinates': coordinates,
        'centroid': centroid,
        'ca_coords': ca_coords,
        'csv_path': csv_path
    }

cmd.extend("get_gloop_coordinates", get_gloop_coordinates)


def highlight_gloop_with_coords(csv_path=None, obj=None, chain=None, start_resi=None, end_resi=None,
                                 surface_color="yellow", show_surface=True, show_coords=True,
                                 output_coords_csv=None):
    """
    综合功能：高亮 G-loop 表面并Export坐标
    
    可以从 CSV File读取 G-loop Information，或直接指定Parameters
    
    Parameters:
        csv_path: G-Motif CSV FilePath (可选，如果提供则从中读取First hit)
        obj: PyMOL 对象Name
        chain: 链 ID (如果 csv_path 提供则可选)
        start_resi: 起始残基Number
        end_resi: 结束残基Number
        surface_color: 表面颜色
        show_surface: 是否Display表面
        show_coords: 是否提取坐标
        output_coords_csv: 坐标输出 CSV Path
    
    Return:
        dict: {
            'surface_result': {...},  # highlight_gloop_surface 的ReturnValue
            'coords_result': {...},   # get_gloop_coordinates 的ReturnValue
        }
    """
    # 如果提供了 CSV，从中读取 G-loop Information
    if csv_path and os.path.exists(csv_path):
        import csv as csv_module
        with open(csv_path, 'r', encoding='utf-8-sig') as f:
            reader = csv_module.DictReader(f)
            for row in reader:
                # 只处理 Pass Status的 hit
                if row.get('Status', '').strip() == 'Pass':
                    chain = row.get('Chain', chain)
                    start_resi = row.get('Start', start_resi)
                    end_resi = row.get('End', end_resi)
                    break
    
    result = {
        'surface_result': None,
        'coords_result': None
    }
    
    # 高亮表面
    if show_surface:
        result['surface_result'] = highlight_gloop_surface(
            obj=obj, chain=chain, start_resi=start_resi, end_resi=end_resi,
            surface_color=surface_color
        )
    
    # 提取坐标
    if show_coords:
        result['coords_result'] = get_gloop_coordinates(
            obj=obj, chain=chain, start_resi=start_resi, end_resi=end_resi,
            atom_types=["all"],  # 提取所有原子
            output_csv=output_coords_csv
        )
    
    return result

cmd.extend("highlight_gloop_with_coords", highlight_gloop_with_coords)
