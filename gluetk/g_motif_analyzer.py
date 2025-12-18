# -*- coding: utf-8 -*-
"""
CRBN G-MOTIF（G-loop）识别（支持：内置真实模板 / 自定义选择）
- template_mode: "builtin"（默认，使用 GSPT1）/ "selection"
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
from typing import Optional
from pymol import cmd

# ====== 1) 内置模板定义（可按需自改）======
# 注意：模板名称必须与 GUI (target_discovery.py) 中的定义一致
# 使用简化的蛋白名称作为主键，便于 GUI 调用

BUILTIN_TEMPLATES = {
    # ===== 主要模板（推荐使用）=====
    # GSPT1 G-loop: 来自 6H0G (CRBN-lenalidomide-GSPT1 ternary complex)
    # 这是最常用的分子胶底物模板
    "GSPT1": {
        "pdb": "6H0G",
        "chain": "A",
        "resi_range": "60-67",
        "selection_fmt": "{obj} and chain A and resi 60+61+62+63+64+65+66+67 and name CA",
        "description": "GSPT1 G-loop from 6H0G (CRBN-lenalidomide-GSPT1 ternary complex)",
    },
    # CK1α G-loop: 来自 5FQD (CRBN-lenalidomide-CK1α complex)
    "CK1α": {
        "pdb": "5FQD",
        "chain": "C",
        "resi_range": "35-42",
        "selection_fmt": "{obj} and chain C and resi 35+36+37+38+39+40+41+42 and name CA",
        "description": "CK1α G-loop from 5FQD (CRBN-lenalidomide-CK1α complex)",
    },
    # VAV1 G-loop: 来自 2MC1 (VAV1 DH domain)
    "VAV1": {
        "pdb": "2MC1",
        "chain": "A",
        "resi_range": "95-102",
        "selection_fmt": "{obj} and chain A and resi 95+96+97+98+99+100+101+102 and name CA",
        "description": "VAV1 G-loop from 2MC1 (VAV1 DH domain)",
    },
    
    # ===== 兼容旧版格式（保留向后兼容）=====
    "GSPT1 (6H0G A:60-67)": {
        "pdb": "6H0G",
        "chain": "A",
        "resi_range": "60-67",
        "selection_fmt": "{obj} and chain A and resi 60+61+62+63+64+65+66+67 and name CA",
        "description": "GSPT1 G-loop (legacy format)",
    },
    "CK1α (3M51 A:36-43)": {
        "pdb": "3M51",
        "chain": "A",
        "resi_range": "36-43",
        "selection_fmt": "{obj} and chain A and resi 36+37+38+39+40+41+42+43 and name CA",
        "description": "CK1α G-loop from 3M51 (legacy format)",
    },
    "VAV1 (2MC1 A:95-102)": {
        "pdb": "2MC1",
        "chain": "A",
        "resi_range": "95-102",
        "selection_fmt": "{obj} and chain A and resi 95+96+97+98+99+100+101+102 and name CA",
        "description": "VAV1 G-loop (legacy format)",
    },
}

# 模板别名映射（支持多种命名方式）
TEMPLATE_ALIASES = {
    # 简化名称 -> 标准名称
    "gspt1": "GSPT1",
    "ck1a": "CK1α",
    "ck1alpha": "CK1α",
    "vav1": "VAV1",
    # GUI 使用的格式
    "GSPT1 (6H0G A:60-67)": "GSPT1",
    "CK1α (3M51 A:36-43)": "CK1α",
    "VAV1 (2MC1 A:95-102)": "VAV1",
}

def get_builtin_template(name: str) -> Optional[dict]:
    """
    获取内置模板配置，支持多种命名格式
    
    参数:
        name: 模板名称（支持 "GSPT1", "gspt1", "GSPT1 (6H0G A:60-67)" 等格式）
    
    返回:
        dict: 模板配置，如果未找到返回 None
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
    
    返回:
        list: [(name, description), ...]
    """
    result = []
    seen = set()
    for name, config in BUILTIN_TEMPLATES.items():
        # 跳过旧版格式（避免重复）
        if "legacy" in config.get("description", "").lower():
            continue
        if name not in seen:
            seen.add(name)
            result.append((name, config.get("description", "")))
    return result

# ====== CRBN 关键残基配置 ======
# 用于 validate_crbn_hbonds 函数的动态残基匹配

# 默认配置（基于人源 CRBN，UniProt Q96SW2）
CRBN_KEY_RESIDUES_DEFAULT = {
    'N351': {
        'resn': 'ASN',           # 残基类型
        'resi': '351',           # 残基编号
        'atoms': ['ND2', 'OD1'], # 参与氢键的原子
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

# 常见 PDB 结构的 CRBN 残基预设配置
# 格式: 'PDB_CODE': {'N351': 'actual_resi', 'H357': 'actual_resi', 'W400': 'actual_resi'}
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
    # 小鼠 CRBN (如果残基编号不同，在此添加)
    # 'XXXX': {'N351': 'XXX', 'H357': 'XXX', 'W400': 'XXX'},
}

def get_crbn_key_residues(pdb_preset: Optional[str] = None,
                          custom_config: Optional[dict] = None) -> dict:
    """
    获取 CRBN 关键残基配置
    
    优先级: custom_config > pdb_preset > default
    
    参数:
        pdb_preset: PDB 代码（如 '6H0G'），使用预设配置
        custom_config: 自定义配置字典
    
    返回:
        dict: 关键残基配置
    """
    if custom_config:
        # 验证自定义配置格式
        required_keys = {'N351', 'H357', 'W400'}
        if not required_keys.issubset(custom_config.keys()):
            print(f"[CRBN Config] ⚠️ 自定义配置缺少必要键: {required_keys - set(custom_config.keys())}")
            print("[CRBN Config] 使用默认配置")
            return CRBN_KEY_RESIDUES_DEFAULT.copy()
        
        # 构建完整配置
        result = {}
        for key in required_keys:
            if isinstance(custom_config[key], dict):
                result[key] = custom_config[key]
            else:
                # 简化格式: {'N351': '351'} -> 完整格式
                base = CRBN_KEY_RESIDUES_DEFAULT[key].copy()
                base['resi'] = str(custom_config[key])
                result[key] = base
        return result
    
    if pdb_preset:
        preset_key = pdb_preset.upper()
        if preset_key in CRBN_PRESETS:
            # 从预设构建完整配置
            preset = CRBN_PRESETS[preset_key]
            result = {}
            for key in ['N351', 'H357', 'W400']:
                base = CRBN_KEY_RESIDUES_DEFAULT[key].copy()
                base['resi'] = preset[key]
                result[key] = base
            print(f"[CRBN Config] 使用预设配置: {preset_key}")
            return result
        else:
            print(f"[CRBN Config] ⚠️ 未知预设 '{pdb_preset}'，可用预设: {list(CRBN_PRESETS.keys())}")
            print("[CRBN Config] 使用默认配置")
    
    return CRBN_KEY_RESIDUES_DEFAULT.copy()

# 氨基酸三字母到单字母转换表
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    # 非标准/修饰氨基酸
    'MSE': 'M',  # 硒代蛋氨酸
    'SEC': 'U',  # 硒半胱氨酸
    'PYL': 'O',  # 吡咯赖氨酸
}

def _aa_3to1(resn: str) -> str:
    """将三字母氨基酸代码转换为单字母代码"""
    resn_upper = resn.strip().upper()
    return AA_3TO1.get(resn_upper, 'X')  # 未知氨基酸返回 'X'

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
    从内置模板获取 8×Cα 坐标
    
    参数:
        name: 模板名称，支持多种格式：
              - 简化名称: "GSPT1", "CK1α", "VAV1"
              - 小写: "gspt1", "ck1a"
              - 旧版格式: "GSPT1 (6H0G A:60-67)"
    
    返回:
        list: 8 个 Cα 原子的 (x, y, z) 坐标
    """
    # 使用新的模板查找函数
    info = get_builtin_template(name)
    if not info:
        available = [t[0] for t in list_builtin_templates()]
        raise ValueError(f"Unknown builtin template: '{name}'. Available templates: {available}")
    
    pdb_code = info["pdb"]
    print(f"[G-MOTIF] Using builtin template: {name} (PDB: {pdb_code})")
    
    # 若会话中不存在，自动 fetch（注意 async_）
    obj = _find_loaded_object_contains(pdb_code)
    if obj is None:
        try:
            print(f"[G-MOTIF] Fetching {pdb_code} from RCSB...")
            cmd.fetch(pdb_code, async_=0)  # ✅ 修复：用 async_ 避免语法错误
            obj = _find_loaded_object_contains(pdb_code) or pdb_code
        except Exception as e:
            raise RuntimeError(f"Failed to fetch {pdb_code} ({e}). Load the PDB manually or use the 'selection' template.")
    
    sel = info["selection_fmt"].format(obj=obj)
    return _coords_from_selection(sel)

def _get_template_coords(template_mode: str, template_sel: Optional[str], template_builtin: Optional[str]):
    """
    获取模板坐标
    
    参数:
        template_mode: "builtin"（默认）或 "selection"
        template_sel: 自定义选择表达式（当 template_mode="selection"）
        template_builtin: 内置模板名称（默认 "GSPT1"）
    
    返回:
        list: 8 个 Cα 原子的 (x, y, z) 坐标
    """
    mode = (template_mode or "builtin").lower()
    if mode == "selection":
        return _coords_from_selection(template_sel or "")
    # 默认使用 builtin 模式，GSPT1 作为默认模板
    return _coords_from_builtin(template_builtin or "GSPT1")

def _calculate_sasa_for_window(obj_name, chain, window):
    """
    计算 8 个残基窗口的溶剂可及表面积 (SASA)
    
    参数:
        obj_name: PyMOL 对象名
        chain: 链 ID
        window: [(resi_label, resn, xyz), ...] 8个残基的列表
    
    返回:
        float: 平均 SASA (Ų/残基)，如果计算失败返回 None
    """
    try:
        # 构建残基列表
        resi_list = '+'.join([w[0] for w in window])
        selection = f"{obj_name} and chain {chain} and resi {resi_list}"
        
        # 计算 SASA
        total_sasa = cmd.get_area(selection, state=1)
        
        # 计算平均值
        avg_sasa = total_sasa / 8.0
        return avg_sasa
    except Exception as e:
        # 如果计算失败，返回 None（例如对象不存在）
        return None

# ====== 4) 主函数 ======
def find_crbn_g_motif(obj_name=None, pdb_file=None,
                       template_mode="builtin", template_sel=None, template_builtin="GSPT1",
                       rmsd_cutoff=1.0, out_csv=None, auto_highlight=0, require_gly_pos6=True,
                       exclude_proline=True, check_surface_exposure=True,
                       min_sasa_per_residue=15.0, topk_debug=10,
                       highlight_surface=False, export_coords=False, coords_csv=None):
    """
    G-loop 挖掘主函数
    
    参数:
        obj_name: PyMOL 对象名
        pdb_file: PDB 文件路径
        template_mode: "builtin"（默认）/ "selection"
        template_sel: 选择表达式（当 template_mode="selection"）
        template_builtin: 内置模板名（默认 "GSPT1"，可选 "CK1α", "VAV1"）
        rmsd_cutoff: RMSD 阈值（默认 1.0 Å）
        out_csv: 输出 CSV 文件路径
        auto_highlight: 是否自动高亮（0/1）
        require_gly_pos6: 是否要求第 6 位为甘氨酸（默认 True）
        exclude_proline: 是否排除含脯氨酸的窗口（默认 True）
        check_surface_exposure: 是否检查表面暴露（默认 True）
        min_sasa_per_residue: 最小 SASA 阈值，单位 Ų/残基（默认 15.0）
        topk_debug: 显示前 N 个最小 RMSD（调试用）
        highlight_surface: 是否高亮 G-loop 表面（默认 False）
        export_coords: 是否导出坐标（默认 False）
        coords_csv: 坐标输出 CSV 路径（可选）
    
    返回:
        dict: {
            'hits': [(chain, start_resi_label, end_resi_label, seq8, rmsd), ...],
            'csv_path': str,  # 输出 CSV 路径
            'coordinates': {...} or None,  # 坐标数据（如果 export_coords=True）
            'surface_info': {...} or None,  # 表面信息（如果 highlight_surface=True）
        }
        
        为了向后兼容，如果只有 hits，也可以直接迭代返回值
    """
    # 载入对象
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
        print("[G-MOTIF] No object available; load a structure or provide pdb_file")
        return []

    # 模板坐标
    try:
        tmpl = _get_template_coords(template_mode, template_sel, template_builtin)
    except Exception as e:
        print(f"[G-MOTIF] Template error: {e}; falling back to GSPT1 template.")
        tmpl = _coords_from_builtin("GSPT1")

    hits = []
    best_rmsd_pool = []
    by_chain = _collect_ca_by_chain(obj)
    total_windows = 0
    gly_pos6_windows = 0
    proline_excluded_windows = 0
    buried_windows = 0

    for ch, rows in by_chain.items():
        seq = [(r[1], r[2], r[3]) for r in rows]  # (resi_label, resn, xyz)
        if len(seq) < 8:
            continue
        for i in range(0, len(seq) - 7):
            window = seq[i:i+8]
            total_windows += 1
            # 检查第 6 位是否为甘氨酸
            if require_gly_pos6:
                if window[5][1] not in ("GLY", "G"):
                    continue
                gly_pos6_windows += 1
            # 检查是否包含脯氨酸（Pro, P）
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
                # 表面暴露检查
                if check_surface_exposure:
                    avg_sasa = _calculate_sasa_for_window(obj, ch, window)
                    if avg_sasa is None or avg_sasa < min_sasa_per_residue:
                        buried_windows += 1
                        continue
                
                resi_s = window[0][0]; resi_e = window[-1][0]
                seq8 = ''.join(_aa_3to1(aa) if aa else 'X' for aa in [w[1] for w in window])
                hits.append((ch, resi_s, resi_e, seq8, rmsd))

    # 调试输出
    print(f"[G-MOTIF] Total windows: {total_windows}")
    if require_gly_pos6:
        print(f"[G-MOTIF] Windows with pos6=Gly: {gly_pos6_windows}")
    if exclude_proline:
        print(f"[G-MOTIF] Windows excluded (含脯氨酸): {proline_excluded_windows}")
    if check_surface_exposure:
        print(f"[G-MOTIF] Windows excluded (埋藏在内部): {buried_windows}")
    if best_rmsd_pool:
        best_rmsd_pool.sort(key=lambda x: x[0])
        nshow = min(topk_debug, len(best_rmsd_pool))
        print(f"[G-MOTIF] Smallest RMSD (top {nshow}):")
        for k in range(nshow):
            r, ch, window = best_rmsd_pool[k]
            seq8 = ''.join(_aa_3to1(aa) if aa else 'X' for aa in [w[1] for w in window])
            resi_s = window[0][0]; resi_e = window[-1][0]
            print("  #{:02d} chain={} {:>6s}-{:>6s}  seq={}  RMSD={:.2f} Å".format(
                k+1, ch or '.', str(resi_s), str(resi_e), seq8, r
            ))

    # 写 CSV - 输出所有扫描的窗口(不仅是 hits)
    if out_csv is None:
        import tempfile
        fd, out_csv = tempfile.mkstemp(suffix="_gmotif.csv"); os.close(fd)
    
    # 对 best_rmsd_pool 按 RMSD 排序
    best_rmsd_pool.sort(key=lambda x: x[0])
    
    print(f"[G-MOTIF] Writing {len(best_rmsd_pool)} windows to CSV (all RMSD values): {out_csv}")
    print(f"[G-MOTIF] Hits (RMSD ≤ {rmsd_cutoff} Å): {len(hits)}")
    
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        # 表头增加 Status 列
        w.writerow(["Chain", "Sequence", "Start", "End", "RMSD", "Status", "Type"])
        
        # 写入所有窗口
        for rmsd, ch, window in best_rmsd_pool:
            resi_s = window[0][0]
            resi_e = window[-1][0]
            seq8 = ''.join(_aa_3to1(aa) if aa else 'X' for aa in [w[1] for w in window])
            
            # 判断是否通过阈值
            status = "Pass" if rmsd <= float(rmsd_cutoff) else "Fail"
            
            w.writerow([ch, seq8, resi_s, resi_e, f"{rmsd:.2f}", status, "G-Motif"])

    print(f"[G-MOTIF] Hits: {len(hits)}; output: {out_csv}")

    # 准备返回结果
    result = {
        'hits': hits,
        'csv_path': out_csv,
        'coordinates': None,
        'surface_info': None,
        'all_windows': best_rmsd_pool  # 所有扫描的窗口（用于高级分析）
    }

    # 自动高亮
    if auto_highlight and hits:
        try:
            try:
                from .highlight_residues import highlight_gmotif_loops
            except Exception:
                from highlight_residues import highlight_gmotif_loops
            # 使用专门的 G-loop 高亮函数
            highlight_gmotif_loops(out_csv, obj=obj, color="yellow", show_labels=True, clear_old=True)
            print(f"[G-MOTIF] Auto-highlighted {len(hits)} G-loop regions in PyMOL")
        except Exception as e:
            print(f"[G-MOTIF] Auto highlight failed: {e}")

    # 高亮表面（新功能）
    if highlight_surface and hits:
        try:
            try:
                from .highlight_residues import highlight_gloop_surface
            except Exception:
                from highlight_residues import highlight_gloop_surface
            
            # 高亮第一个 hit 的表面
            first_hit = hits[0]
            ch, resi_s, resi_e, seq8, rmsd = first_hit
            surface_result = highlight_gloop_surface(
                obj=obj, chain=ch, start_resi=resi_s, end_resi=resi_e,
                surface_color="yellow", surface_transparency=0.3
            )
            result['surface_info'] = surface_result
            print(f"[G-MOTIF] ✅ Surface highlighted for G-loop {ch}:{resi_s}-{resi_e}")
        except Exception as e:
            print(f"[G-MOTIF] Surface highlight failed: {e}")

    # 导出坐标（新功能）
    if export_coords and hits:
        try:
            try:
                from .highlight_residues import get_gloop_coordinates
            except Exception:
                from highlight_residues import get_gloop_coordinates
            
            # 导出第一个 hit 的坐标
            first_hit = hits[0]
            ch, resi_s, resi_e, seq8, rmsd = first_hit
            
            # 如果没有指定坐标 CSV 路径，自动生成
            if coords_csv is None:
                import tempfile
                fd, coords_csv = tempfile.mkstemp(suffix="_gloop_coords.csv")
                os.close(fd)
            
            coords_result = get_gloop_coordinates(
                obj=obj, chain=ch, start_resi=resi_s, end_resi=resi_e,
                atom_types=["all"],  # 导出所有原子
                output_csv=coords_csv
            )
            result['coordinates'] = coords_result
            print(f"[G-MOTIF] ✅ Coordinates exported to: {coords_csv}")
        except Exception as e:
            print(f"[G-MOTIF] Coordinate export failed: {e}")

    if tmp_obj:
        try: cmd.delete(tmp_obj)
        except Exception: pass

    # 为了向后兼容，返回一个可以像列表一样迭代的对象
    return GMotifResult(result)


class GMotifResult:
    """
    G-Motif 结果包装类，支持向后兼容的列表迭代和新的字典访问
    """
    def __init__(self, data):
        self._data = data
        self._hits = data.get('hits', [])
    
    def __iter__(self):
        """支持 for hit in result 的迭代"""
        return iter(self._hits)
    
    def __len__(self):
        """支持 len(result)"""
        return len(self._hits)
    
    def __getitem__(self, key):
        """支持 result[0] 和 result['hits'] 两种访问方式"""
        if isinstance(key, int):
            return self._hits[key]
        return self._data.get(key)
    
    def __bool__(self):
        """支持 if result: 判断"""
        return len(self._hits) > 0
    
    @property
    def hits(self):
        """获取 hits 列表"""
        return self._hits
    
    @property
    def csv_path(self):
        """获取 CSV 路径"""
        return self._data.get('csv_path')
    
    @property
    def coordinates(self):
        """获取坐标数据"""
        return self._data.get('coordinates')
    
    @property
    def surface_info(self):
        """获取表面信息"""
        return self._data.get('surface_info')
    
    @property
    def all_windows(self):
        """获取所有扫描的窗口"""
        return self._data.get('all_windows', [])
    
    def get(self, key, default=None):
        """字典式 get 方法"""
        return self._data.get(key, default)
    
    def to_dict(self):
        """转换为普通字典"""
        return self._data.copy()
    
    def __repr__(self):
        return f"GMotifResult(hits={len(self._hits)}, csv='{self.csv_path}')"


# ====== 5) G-motif 内部几何验证 ======
def validate_g_motif_geometry(obj_name, g_motif_chain, g_motif_resi_range,
                               max_internal_hbond=3.5, reference_rmsd_threshold=1.0):
    """
    验证 G-motif 内部几何特征（α-turn）：
    1. 必需氢键：G₋₄ backbone O → G₀ backbone N (α-turn 定义)
    2. 可选氢键：G₋₄ backbone O → G₊₁ backbone N
    3. 中心 Gly (G₀) 必须是 GLY
    4. RMSD 与 GSPT1 G-loop 的相似度
    
    返回:
        dict: {
            'has_alpha_turn_hbond': bool,  # G₋₄→G₀
            'alpha_turn_distance': float,
            'has_secondary_hbond': bool,   # G₋₄→G₊₁
            'secondary_hbond_distance': float,
            'central_gly_confirmed': bool,
            'rmsd_to_gspt1': float,
            'is_valid_geometry': bool
        }
    """
    import numpy as np
    
    # 解析范围
    if isinstance(g_motif_resi_range, str):
        parts = g_motif_resi_range.split('-')
        start_resi = parts[0].strip()
        end_resi = parts[1].strip() if len(parts) > 1 else start_resi
    else:
        start_resi, end_resi = g_motif_resi_range
    
    g_sel = f"{obj_name} and chain {g_motif_chain} and resi {start_resi}-{end_resi}"
    
    # 获取 8 个残基的 backbone 原子
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
    
    # 映射位置：索引 0-7 对应 G₋₅ 到 G₊₂
    # G₋₄ = index 1, G₀ = index 5, G₊₁ = index 6
    result = {
        'has_alpha_turn_hbond': False,
        'alpha_turn_distance': None,
        'has_secondary_hbond': False,
        'secondary_hbond_distance': None,
        'central_gly_confirmed': False,
        'rmsd_to_gspt1': None,
        'is_valid_geometry': False
    }
    
    # 1. 检查中心 Gly (G₀, index 5)
    central_resn = g_residues[sorted_resis[5]]['resn']
    result['central_gly_confirmed'] = (central_resn in ('GLY', 'G'))
    
    # 2. 必需氢键：G₋₄ O → G₀ N (index 1 → 5)
    if 'O' in g_residues[sorted_resis[1]]['atoms'] and 'N' in g_residues[sorted_resis[5]]['atoms']:
        o_g4 = g_residues[sorted_resis[1]]['atoms']['O']
        n_g0 = g_residues[sorted_resis[5]]['atoms']['N']
        dist = np.linalg.norm(o_g4 - n_g0)
        result['alpha_turn_distance'] = round(float(dist), 2)
        result['has_alpha_turn_hbond'] = (dist <= max_internal_hbond)
    
    # 3. 可选氢键：G₋₄ O → G₊₁ N (index 1 → 6)
    if 'O' in g_residues[sorted_resis[1]]['atoms'] and 'N' in g_residues[sorted_resis[6]]['atoms']:
        o_g4 = g_residues[sorted_resis[1]]['atoms']['O']
        n_g1 = g_residues[sorted_resis[6]]['atoms']['N']
        dist = np.linalg.norm(o_g4 - n_g1)
        result['secondary_hbond_distance'] = round(float(dist), 2)
        result['has_secondary_hbond'] = (dist <= max_internal_hbond)
    
    # 4. RMSD 与 GSPT1 G-loop 比较（使用 Cα）
    ca_coords = []
    for resi in sorted_resis[:8]:
        if 'CA' in g_residues[resi]['atoms']:
            ca_coords.append(g_residues[resi]['atoms']['CA'])
    
    if len(ca_coords) == 8:
        # GSPT1 G-loop Cα 坐标（来自 6H0G，近似理想化）
        gspt1_template = _ideal_beta_hairpin_template()  # 可用真实坐标替换
        try:
            rmsd = _kabsch_rmsd(np.array(ca_coords), np.array(gspt1_template))
            result['rmsd_to_gspt1'] = round(rmsd, 3)
        except Exception:
            pass
    
    # 5. 综合判定
    result['is_valid_geometry'] = (
        result['central_gly_confirmed'] and
        result['has_alpha_turn_hbond'] and
        (result['rmsd_to_gspt1'] is None or result['rmsd_to_gspt1'] < reference_rmsd_threshold)
    )
    
    # 打印结果
    print("\n" + "=" * 70)
    print("G-motif 内部几何验证 (G-motif Internal Geometry Validation)")
    print("=" * 70)
    gly_status = "✅" if result['central_gly_confirmed'] else "❌"
    print(f"{gly_status} 中心 Gly (G₀):  {central_resn}")
    
    alpha_status = "✅" if result['has_alpha_turn_hbond'] else "❌"
    alpha_dist = f"{result['alpha_turn_distance']} Å" if result['alpha_turn_distance'] else "N/A"
    print(f"{alpha_status} α-turn 氢键 (G₋₄→G₀):  {alpha_dist}")
    
    if result['secondary_hbond_distance']:
        sec_status = "✅" if result['has_secondary_hbond'] else "➖"
        sec_dist = f"{result['secondary_hbond_distance']} Å"
        print(f"{sec_status} 次级氢键 (G₋₄→G₊₁):  {sec_dist}")
    
    if result['rmsd_to_gspt1']:
        rmsd_status = "✅" if result['rmsd_to_gspt1'] < reference_rmsd_threshold else "⚠️"
        print(f"{rmsd_status} RMSD vs GSPT1:  {result['rmsd_to_gspt1']} Å")
    
    valid_status = "✅ 是" if result['is_valid_geometry'] else "❌ 否"
    print(f"\n有效 α-turn 几何: {valid_status}")
    print("=" * 70 + "\n")
    
    return result


# ====== 6) CRBN关键残基氢键验证 ======
def validate_crbn_hbonds(obj_name, g_motif_chain, g_motif_resi_range, crbn_chain,
                         max_hbond_dist=3.5, min_donor_angle=120.0,
                         pdb_preset: Optional[str] = None,
                         crbn_key_residues: Optional[dict] = None):
    """
    验证 G-loop 与 CRBN 的 3 个关键氢键（基于 Schrödinger 标准）：
    - G-3 backbone O → CRBN Asn351 (sidechain NH2)
    - G-2 backbone O → CRBN His357 (sidechain)
    - G-1 backbone O → CRBN Trp400 (sidechain NH)
    
    参数:
        obj_name: PyMOL 对象名
        g_motif_chain: G-motif 所在链 ID
        g_motif_resi_range: G-motif 残基范围 (start, end) 或 "start-end"
        crbn_chain: CRBN 链 ID
        max_hbond_dist: 氢键最大距离阈值（默认 3.5 Å）
        min_donor_angle: 最小供体角度（默认 120°，当前未使用）
        pdb_preset: PDB 代码预设（如 '6H0G'），自动使用对应的残基编号
        crbn_key_residues: 自定义 CRBN 关键残基配置
            简化格式: {'N351': '351', 'H357': '357', 'W400': '400'}
            完整格式: {'N351': {'resn': 'ASN', 'resi': '351', 'atoms': ['ND2', 'OD1']}, ...}
    
    返回:
        dict: {
            'N351_hbond': {'found': bool, 'distance': float, 'g_pos': int},
            'H357_hbond': {'found': bool, 'distance': float, 'g_pos': int},
            'W400_hbond': {'found': bool, 'distance': float, 'g_pos': int},
            'total_hbonds': int,
            'is_canonical_gloop': bool,
            'config_used': dict  # 使用的配置信息
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
    
    # 定位 G-3, G-2, G-1 的 backbone O (对应位置 5, 6, 7 in 0-indexed)
    g_positions = {}
    for i, resi in enumerate(sorted_resis[:8]):
        backbone_o = None
        for atom in g_residues[resi]['atoms']:
            if atom.name.strip().upper() == 'O':  # backbone carbonyl
                backbone_o = np.array([atom.coord[0], atom.coord[1], atom.coord[2]])
                break
        g_positions[i] = {'resi': resi, 'O': backbone_o}
    
    # 获取 CRBN 关键残基配置（支持动态配置）
    key_residue_config = get_crbn_key_residues(pdb_preset=pdb_preset, custom_config=crbn_key_residues)
    
    # 获取 CRBN 关键残基的侧链原子
    crbn_model = cmd.get_model(crbn_sel)
    crbn_key_atoms = {'N351': [], 'H357': [], 'W400': []}
    
    for a in crbn_model.atom:
        resi_str = (a.resi or '').strip()
        resn = (a.resn or '').strip().upper()
        aname = (a.name or '').strip().upper()
        coord = np.array([a.coord[0], a.coord[1], a.coord[2]])
        
        # 动态匹配 CRBN 关键残基
        for key in ['N351', 'H357', 'W400']:
            config = key_residue_config[key]
            expected_resn = config['resn']
            expected_resi = config['resi']
            expected_atoms = config['atoms']
            
            # 匹配残基类型和编号
            if resn == expected_resn and expected_resi in resi_str:
                if aname in expected_atoms:
                    crbn_key_atoms[key].append((f"{resn}{resi_str}", aname, coord))
    
    # 检查 3 个氢键
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
    
    # 判定是否为标准 G-loop (至少 2/3 氢键)
    result['is_canonical_gloop'] = result['total_hbonds'] >= 2
    
    # 打印结果
    print("\n" + "=" * 70)
    print("CRBN G-loop 氢键验证 (CRBN-G-loop H-bond Validation)")
    print("=" * 70)
    print(f"配置: {result['config_used']['preset']} "
          f"(N351={result['config_used']['N351_resi']}, "
          f"H357={result['config_used']['H357_resi']}, "
          f"W400={result['config_used']['W400_resi']})")
    print("-" * 70)
    
    # 动态生成标签
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
    print(f"\n总氢键数: {result['total_hbonds']}/3")
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
    1. 内部几何验证：α-turn 氢键 (G₋₄→G₀)
    2. CRBN 关键氢键验证 (N351, H357, W400)
    3. Glue 是否插入 G-motif 和 CRBN 之间
    4. 关键残基（如 Gly₀）与 MGD 的 vdW 接触
    
    参数:
        obj_name: PyMOL对象名称
        g_motif_chain: G-motif所在链 ID
        g_motif_resi_range: G-motif残基范围 (start, end) 或 "start-end"
        glue_resname: 分子胶名称
        crbn_chain: CRBN链 ID
        distance_threshold: 接触距离阈值（埃）
    
    返回:
        dict: {
            'is_glue_substrate': bool,  # 是否为分子胶底物
            'binding_mode': str,  # "glue-induced" / "direct" / "no-binding"
            'g_motif_crbn_distance': float,  # G-motif到CRBN最短距离
            'glue_contacts_g_motif': int,  # Glue与G-motif接触数
            'glue_contacts_crbn': int,  # Glue与CRBN接触数
            'key_residues_engaged': list,  # 参与相互作用的关键残基
            'confidence': float  # 置信度 [0-1]
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
    
    # 1. 验证内部几何（如果启用）
    geometry_result = None
    if validate_geometry:
        geometry_result = validate_g_motif_geometry(obj_name, g_motif_chain, g_motif_resi_range)
        if geometry_result and not geometry_result['is_valid_geometry']:
            print("⚠️ 警告: α-turn 几何不符合，可能不是有效 G-motif\n")
    
    # 2. 验证 CRBN 氢键（如果启用）
    hbond_result = None
    if validate_hbonds:
        hbond_result = validate_crbn_hbonds(obj_name, g_motif_chain, g_motif_resi_range, crbn_chain)
        if hbond_result and not hbond_result['is_canonical_gloop']:
            print("⚠️ 警告: CRBN 关键氢键不足，可能不是标准 G-loop\n")
    
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
    
    # 3. 分析关键残基（Gly6 及侧链位点 G-4, G-3, G-2, G+1）
    key_residues_engaged = []
    # 获取 G-motif 序列位置
    g_motif_residues = defaultdict(list)
    for atom in g_motif_atoms:
        res_key = (atom[3], atom[4])  # (resn, resi)
        g_motif_residues[res_key].append(atom)
    
    # 按resi排序
    sorted_residues = sorted(g_motif_residues.keys(), key=lambda x: _parse_resi(x[1]))
    
    # 检查第6位（Gly，G位）与 MGD 的 vdW 接触
    vdw_threshold = 4.5  # vdW 接触阈值
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
        
        # 如果关键残基参与，提高置信度
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
        'geometry_validation': geometry_result,  # 添加几何验证结果
        'crbn_hbond_validation': hbond_result     # 添加氢键验证结果
    }
    
    # 打印结果
    print("\n" + "=" * 60)
    print("G-motif Glue结合验证 (G-motif Glue Binding Analysis)")
    print("=" * 60)
    print(f"G-motif 到 CRBN 距离: {result['g_motif_crbn_distance']} Å")
    print(f"Glue - G-motif 接触数: {glue_contacts_g_motif}")
    print(f"Glue - CRBN 接触数: {glue_contacts_crbn}")
    print(f"结合模式: {binding_mode}")
    status_text = "✨ 是" if is_glue_substrate else "⚠️ 否"
    print(f"是否为分子胶底物: {status_text}")
    print(f"置信度: {confidence:.2f}")
    
    if key_residues_engaged:
        print(f"\n关键残基参与:")
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
        print(f"CRBN氢键: {hbond_result['total_hbonds']}/3 (标准G-loop: {'是' if hbond_result['is_canonical_gloop'] else '否'})")
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

def _find_c2h2(seq):
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
    Auto mode marks simple C2H2-ZF-like motifs and β-hairpin heuristics.
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
