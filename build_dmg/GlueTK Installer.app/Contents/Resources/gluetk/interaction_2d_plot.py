# -*- coding: utf-8 -*-
"""
interaction_2d_plot.py
生成蛋白质-配体2D相互作用图 (Schrödinger Style v2.0)

特性：
- **Sector-Based Layout**: 基于扇区分配的智能布局，彻底解决重叠问题
- **Professional Aesthetics**: 使用圆角矩形 Badge 和专业字体，媲美商业软件
- **High Resolution**: 默认 300 DPI 输出
- **Visual Clarity**: 样条曲线(Spline)处理疏水作用，虚线处理氢键
"""

from __future__ import print_function
import os
import csv
import tempfile
import math
import numpy as np

# 氨基酸单字母对照表
AA_THREE_TO_ONE = {
    'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'PHE': 'F', 'MET': 'M', 'PRO': 'P', 'TRP': 'W',
    'GLY': 'G', 'SER': 'S', 'THR': 'T', 'ASN': 'N', 'GLN': 'Q', 'TYR': 'Y', 'CYS': 'C',
    'ASP': 'D', 'GLU': 'E',
    'HIS': 'H', 'LYS': 'K', 'ARG': 'R'
}

# 氨基酸分类颜色方案
AA_CATEGORY_COLORS = {
    'nonpolar': '#FFB74D',      # 🌀 橙色 - 非极性
    'polar': '#64B5F6',         # 💧 蓝色 - 极性不带电
    'negative': '#81C784',      # 🌿 绿色 - 负电荷
    'positive': '#E57373',      # 🌱 红色 - 正电荷
    'unknown': '#BDBDBD'        # 灰色 - 未知
}

def get_aa_category(aa_code):
    """根据氨基酸单字母代码返回类别"""
    nonpolar = set('AVLIFMPW')     # 非极性
    polar = set('GSTNQYC')         # 极性不带电
    negative = set('DE')           # 负电荷
    positive = set('HKR')          # 正电荷
    
    if aa_code in nonpolar:
        return 'nonpolar'
    elif aa_code in polar:
        return 'polar'
    elif aa_code in negative:
        return 'negative'
    elif aa_code in positive:
        return 'positive'
    else:
        return 'unknown'

def parse_residue_label(res_str):
    """解析残基标签并转换为单字母格式
    输入: 'LYS 383' 或 'K383'
    输出: ('K', '383', 'positive')
    """
    parts = res_str.strip().split()
    if len(parts) >= 2:
        # 三字母格式
        aa_three = parts[0].upper()
        aa_one = AA_THREE_TO_ONE.get(aa_three, aa_three[0] if aa_three else 'X')
        num = parts[1]
    elif len(parts) == 1 and any(c.isdigit() for c in parts[0]):
        # 已经是单字母格式 'K383'
        import re
        match = re.match(r'([A-Z]+)(\d+)', parts[0])
        if match:
            aa_one = match.group(1)
            num = match.group(2)
        else:
            aa_one = 'X'
            num = '0'
    else:
        aa_one = 'X'
        num = '0'
    
    category = get_aa_category(aa_one)
    return aa_one, num, category

# GlueTK 统一配色方案 (与 3D 视图保持一致)
SCHRODINGER_STYLE = {
    "Hbond": {
        "color": "#2196F3",      # Blue
        "label": "Hydrogen Bond",
        "style": "-",           
        "arrow": True,
        "linewidth": 1.5,
        "show": True
    },
    "Salt": {
        "color": "#FF5722",      # Orange
        "label": "Salt Bridge",
        "style": "-",
        "arrow": True,
        "linewidth": 1.5,
        "show": True
    },
    "PiPi": {
        "color": "#9C27B0",      # Purple
        "label": "π-π Stack",
        "style": "--",
        "arrow": False,
        "linewidth": 1.5,
        "show": True
    },
    "Pi": {
        "color": "#E91E63",      # Pink
        "label": "π-Cation",
        "style": "--",
        "arrow": False,
        "linewidth": 1.5,
        "show": True
    },
    "Hydrophobic": {
        "color": "#4CAF50",      # Green
        "label": "Hydrophobic",
        "style": "arc",         # Special handling
        "arrow": False,
        "linewidth": 1.0,
        "show": True
    },
    "Halogen": {
        "color": "#FF9800",      # Orange-Yellow
        "label": "Halogen",
        "style": "-",
        "arrow": True,
        "linewidth": 1.5,
        "show": True
    },
    "Metal": {
        "color": "#673AB7",      # Deep Purple
        "label": "Metal",
        "style": "-",
        "arrow": False,
        "linewidth": 1.5,
        "show": True
    },
    "Water": {
        "color": "#00BCD4",      # Cyan
        "label": "Water Bridge",
        "style": ":",
        "arrow": False,
        "linewidth": 1.5,
        "show": True
    },
    "VDW": {
        "color": "#BDBDBD",      # Light Gray
        "label": "vdW",
        "style": ":",
        "arrow": False,
        "linewidth": 1.0,
        "show": False
    }
}

def get_interaction_style(itype):
    """获取相互作用样式"""
    itype = str(itype).lower()
    if "氢键" in itype or "hbond" in itype or "hydrogen" in itype:
        return SCHRODINGER_STYLE["Hbond"]
    elif "盐桥" in itype or "salt" in itype or "ionic" in itype:
        return SCHRODINGER_STYLE["Salt"]
    elif "π-π" in itype or "pipi" in itype or "stacking" in itype:
        return SCHRODINGER_STYLE["PiPi"]
    elif "π" in itype or "pi" in itype or "cation" in itype:
        return SCHRODINGER_STYLE["Pi"]
    elif "卤素" in itype or "halogen" in itype:
        return SCHRODINGER_STYLE["Halogen"]
    elif "金属" in itype or "metal" in itype or "coord" in itype:
        return SCHRODINGER_STYLE["Metal"]
    elif "水桥" in itype or "water" in itype:
        return SCHRODINGER_STYLE["Water"]
    elif "疏水" in itype or "hydrophobic" in itype or "alkyl" in itype:
        return SCHRODINGER_STYLE["Hydrophobic"]
    else:
        return SCHRODINGER_STYLE["VDW"]

def generate_2d_interaction_diagram(csv_path, ligand_resname, pdb_file=None, obj_name=None,
                                     output_path=None, width=2400, height=2000, dpi=300,
                                     show_distance=False, show_vdw=False,
                                     protein_name=None, compact=True):
    """
    生成高精度且清晰的 2D 相互作用图 (Schrödinger 风格 v2.0)
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        print("[2D Diagram] RDKit is required: pip install rdkit")
        return None

    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib.patches import FancyBboxPatch, PathPatch
        from matplotlib.path import Path
        from matplotlib.lines import Line2D
        from PIL import Image
    except ImportError:
        print("[2D Diagram] matplotlib and Pillow are required")
        return None

    print(f"[2D Diagram] Starting Pro Generation for {ligand_resname}...")

    # 1. 提取配体结构 (PDB格式)
    mol_map = None
    mol_draw = None
    temp_pdb = None
    
    if obj_name:
        try:
            from pymol import cmd
            temp_pdb = tempfile.mktemp(suffix=".pdb")
            cmd.save(temp_pdb, f"{obj_name} and resn {ligand_resname}", format="pdb")
            
            if os.path.exists(temp_pdb):
                mol_map = Chem.MolFromPDBFile(temp_pdb, removeHs=False, sanitize=False)
                mol_raw = Chem.MolFromPDBFile(temp_pdb, removeHs=False, sanitize=False)
                try:
                    Chem.SanitizeMol(mol_map)
                    Chem.SanitizeMol(mol_raw)
                    # 强制移除氢原子用于绘图
                    mol_draw = Chem.RemoveHs(mol_raw, implicitOnly=False, updateExplicitCount=True)
                except Exception as e:
                    print(f"[2D Diagram] Sanitization issue: {e}")
                    mol_draw = Chem.RemoveHs(mol_raw) if mol_raw else None
        except Exception as e:
            print(f"[2D Diagram] Failed to extract ligand: {e}")
            if temp_pdb and os.path.exists(temp_pdb): os.remove(temp_pdb)
            return None
    
    if temp_pdb and os.path.exists(temp_pdb):
        os.remove(temp_pdb)

    if mol_map is None or mol_draw is None:
        print(f"[2D Diagram] ❌ Failed to load ligand molecule.")
        return None
    
    # 2. 建立原子映射 (增强版)
    # 目的: 将CSV中的原子名映射到 mol_draw 的原子索引
    # 风险: PDB原子名可能不唯一，或者RemoveHs导致丢失
    
    # 3. 读取相互作用并映射 (基于 3D 坐标匹配 - v3.0 Robust Mapping)
    # 策略: Name -> 3D Coords (in raw PDB) -> Nearest Atom Index (in clean Mol)
    # 优势: 自动处理 RemoveHs 导致的原子名丢失问题；自动处理氢键映射到重原子
    
    # A. 建立原始 PDB 的 Name -> Atom 映射
    map_name_to_atom = {}
    for atom in mol_map.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info:
            name = info.GetName().strip().upper()
            map_name_to_atom[name] = atom 
            
    # B. 检查 3D 构象
    if mol_draw.GetNumConformers() == 0 or mol_map.GetNumConformers() == 0:
        print("[2D Diagram] ❌ Error: Missing 3D coordinates for mapping.")
        return None
        
    draw_conf = mol_draw.GetConformer()
    raw_conf = mol_map.GetConformer()

    interactions = []
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lig_atom_name = row.get("Ligand_Atom", "").strip().upper()
            prot_res = row.get("Protein_Residue", "").strip()
            itype = row.get("Interaction", "").strip()
            dist = row.get("Distance", "").strip()
            if not prot_res: continue

            target_idx = None
            
            # 尝试通过坐标匹配
            if lig_atom_name in map_name_to_atom:
                raw_atom = map_name_to_atom[lig_atom_name]
                raw_pos = raw_conf.GetAtomPosition(raw_atom.GetIdx())
                
                # 在 mol_draw 中寻找最近的原子
                min_d2 = 9999.0
                best_idx = None
                
                for i in range(mol_draw.GetNumAtoms()):
                    pos = draw_conf.GetAtomPosition(i)
                    d2 = (pos.x - raw_pos.x)**2 + (pos.y - raw_pos.y)**2 + (pos.z - raw_pos.z)**2
                    if d2 < min_d2:
                        min_d2 = d2
                        best_idx = i
                
                # 判定阈值: 1.6埃 (覆盖 C-H 键长约 1.09A, N-H 约 1.0A)
                # 如果是重原子对应重原子，距离应接近0
                # 如果是 H 对应重原子，距离约 1.0-1.1
                if best_idx is not None and min_d2 < 2.56: # 1.6^2
                    target_idx = best_idx
                else:
                    print(f"[2D Diagram] ⚠️ Mapping Warning: {lig_atom_name} closest atom dist {math.sqrt(min_d2):.2f}A > 1.6A")
            
            if target_idx is not None:
                interactions.append({
                    "protein": prot_res,
                    "target_idx": target_idx,
                    "type": itype,
                    "distance": dist,
                    "lig_atom_name": lig_atom_name
                })
            else:
                 # Fallback: 如果名字直接匹配成功且坐标匹配失败 (极罕见), 尝试直接名字匹配
                 pass

    # 4. 生成 2D 坐标 (Clean Molecule)
    try:
        from rdkit.Chem import rdCoordGen
        rdCoordGen.AddCoords(mol_draw)
    except:
        AllChem.Compute2DCoords(mol_draw)
            
    num_atoms = mol_draw.GetNumAtoms()
    is_large_molecule = num_atoms > 40
    
    # 绘图设置
    drawer_w, drawer_h = (1200, 1000)
    drawer = rdMolDraw2D.MolDraw2DCairo(drawer_w, drawer_h)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.padding = 0.25
    opts.annotationFontScale = 0.8
    opts.bondLineWidth = 3.5 if not is_large_molecule else 2.5
    opts.comicMode = False # Professional mode

    drawer.DrawMolecule(mol_draw)
    drawer.FinishDrawing()
    png_data = drawer.GetDrawingText()
    
    # 获取坐标
    atom_coords = {}
    for i in range(mol_draw.GetNumAtoms()):
        pt = drawer.GetDrawCoords(i)
        atom_coords[i] = (pt.x, pt.y)

    xs = [p[0] for p in atom_coords.values()]
    ys = [p[1] for p in atom_coords.values()]
    center_x = sum(xs) / len(xs) if xs else drawer_w/2
    center_y = sum(ys) / len(ys) if ys else drawer_h/2

    # 5. Matplotlib 组装 (High Res)
    fig_w, fig_h = width / dpi, height / dpi
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    ax = fig.add_subplot(111)
    
    tmp_png = tempfile.mktemp(suffix=".png")
    with open(tmp_png, "wb") as f:
        f.write(png_data)
    lig_img = Image.open(tmp_png)
    
    ax.imshow(lig_img)
    # 移除 extent 参数，避免坐标系映射错误
    # ax.imshow(lig_img, origin='upper', extent=[0, drawer_w, 0, drawer_h])
    ax.imshow(lig_img, zorder=5) # 提升分子图层级，遮盖背后的连接线
    
    # 显式设置坐标轴范围，确保与 RDKit 像素坐标一致
    # 并在顶部预留空间给标题 (Y轴负方向是向上)
    top_margin = 20 # 缩小顶部留白 (无标题模式)
    ax.set_xlim(0, drawer_w)
    ax.set_ylim(drawer_h, -top_margin)
    ax.axis('off')
    
    os.remove(tmp_png)
    
    # Debug: 绘制原子名 (已禁用，避免视觉混乱)
    # for i, (ax_x, ax_y) in atom_coords.items():
    #     atom = mol_draw.GetAtomWithIdx(i)
    #     info = atom.GetPDBResidueInfo()
    #     if info:
    #         aname = info.GetName().strip()
    #         ax.text(ax_x, ax_y, aname, ha='center', va='center', 
    #                 fontsize=5, color='gray', alpha=0.5, zorder=1)

    # === Layout Engine ===
    residue_data = {}
    for inter in interactions:
        p = inter["protein"]
        if p not in residue_data: residue_data[p] = []
        residue_data[p].append(inter)
    
    # Base layout parameters
    boundary_dist = 100.0 # 距离分子边缘的基础距离
    circle_radius = 32.0 # 氨基酸标签圆形半径 (稍调大)
    if compact: 
        boundary_dist = 90.0
        circle_radius = 25.0
    if is_large_molecule: boundary_dist += 25
    
    # Calculate geometric center of interactions for each residue
    layout_items = []
    for prot_res, inter_list in residue_data.items():
        # Find centroid of interacting atoms
        indices = [i["target_idx"] for i in inter_list if i["target_idx"] in atom_coords]
        if not indices: continue
        
        cx = sum(atom_coords[i][0] for i in indices) / len(indices)
        cy = sum(atom_coords[i][1] for i in indices) / len(indices)
        
        # Vector from molecule center
        vx, vy = cx - center_x, cy - center_y
        raw_angle = math.atan2(vy, vx)
        
        layout_items.append({
            "res": prot_res,
            "inters": inter_list,
            "raw_angle": raw_angle,
            "centroid": (cx, cy),
            "indices": indices
        })
    
    # Sort by angle
    layout_items.sort(key=lambda x: x["raw_angle"])
    
    # Sector Allocation (Fan out)
    min_sep_angle = math.radians(20) # 最小间隔角度
    if len(layout_items) > 15: min_sep_angle = math.radians(15)
    
    adjusted_items = []
    if layout_items:
        # Simple collision resolution
        current_angle = layout_items[0]["raw_angle"]
        adjusted_items.append({**layout_items[0], "final_angle": current_angle})
        
        for i in range(1, len(layout_items)):
            next_angle = layout_items[i]["raw_angle"]
            diff = next_angle - current_angle
            
            # Handle wrap around pi/-pi
            while diff < 0: diff += 2*math.pi
            
            if diff < min_sep_angle:
                current_angle += min_sep_angle
            else:
                current_angle = next_angle
                
            adjusted_items.append({**layout_items[i], "final_angle": current_angle})

    # Draw Items
    drawn_boxes = [] # Keep track for overlap check if needed
    placed_bubbles = [] # list of (x, y, radius)
    
    for idx, item in enumerate(adjusted_items):
        prot_res = item["res"]
        angle = item["final_angle"]
        cx, cy = item["centroid"] # Interaction centroid on ligand
        
        # Calculate badge position
        # Project out from ligand center, but respecting the interaction centroid
        # A mix of (Center->Out) and (Centroid->Out)
        
        # Base vector from center
        dx, dy = math.cos(angle), math.sin(angle)
        
        # Determine "Surface" distance approximation
        # Find furthest atom in this direction
        max_r = 0
        for i in range(num_atoms):
            ax_x, ax_y = atom_coords[i]
            # Project onto direction
            proj = (ax_x - center_x)*dx + (ax_y - center_y)*dy
            if proj > max_r: max_r = proj
            
        radius = max_r + boundary_dist
        
        # Smart collision avoidance
        current_r = radius
        valid_position = False
        
        # Try to place, push out if collision detected
        # Max 5 layers push
        for attempt in range(5):
             px = center_x + dx * current_r
             py = center_y + dy * current_r
             
             collision = False
             # Check against already placed bubbles
             for (bx, by, br) in placed_bubbles:
                 dist_sq = (px - bx)**2 + (py - by)**2
                 min_dist = circle_radius + br + 10.0 # 10px padding
                 if dist_sq < min_dist**2:
                     collision = True
                     break
             
             if not collision:
                 valid_position = True
                 break
             else:
                 # Push out one diameter + padding
                 current_r += (circle_radius * 2 + 15)
        
        radius = current_r # Update radius for subsequent logic if needed
        placed_bubbles.append((px, py, circle_radius))
        
        # Draw Residue Badge (Rounded Rectangle)
        # Determine interaction type color
        color_weights = {"hydrogen": 10, "salt": 9, "pipi": 8, "pi": 7, "metal": 6, "halogen": 5, "hydrophobic": 2}
        best_style = SCHRODINGER_STYLE["VDW"]
        max_w = 0
        
        for inter in item["inters"]:
            itype = inter["type"]
            w = 0
            for k, v in color_weights.items():
                if k in release_itype(itype): w = max(w, v)
            if w > max_w:
                max_w = w
                best_style = get_interaction_style(itype)

        # 解析氨基酸信息
        aa_letter, aa_num, aa_category = parse_residue_label(prot_res)
        badge_color = AA_CATEGORY_COLORS[aa_category]
        
        # 气泡感设计 (Bubble Style)
        import matplotlib.colors as mcolors
        import matplotlib.patches as patches
        
        # 1. 阴影 (Drop Shadow)
        shadow_offset = 5
        shadow = patches.Circle(
            (px + shadow_offset, py + shadow_offset), circle_radius,
            fc='#CFD8DC', ec='none', alpha=0.5, zorder=9
        )
        ax.add_patch(shadow)

        # 2. 主体圆圈 (Main Bubble)
        # 计算极淡的填充色
        base_rgb = mcolors.to_rgb(badge_color)
        light_fill = (*base_rgb, 0.3) # 8%透明度的填充
        
        circle = patches.Circle(
            (px, py), circle_radius,
            fc='white', ec=badge_color, lw=2.0, zorder=10
        )
        # 内层填充
        circle_fill = patches.Circle(
            (px, py), circle_radius - 2,
            fc=light_fill, ec='none', zorder=10
        )
        
        ax.add_patch(circle)
        ax.add_patch(circle_fill)
        
        # 3. 高光 (Highlight)
        highlight = patches.Ellipse(
            (px - circle_radius*0.35, py - circle_radius*0.35), 
            circle_radius*0.7, circle_radius*0.35, angle=140,
            fc='white', alpha=0.5, zorder=11
        )
        ax.add_patch(highlight)
        
        # 在圆形内部显示单字母+数字 (自适应字体大小)
        label_text = f"{aa_letter}{aa_num}"
        # 根据字符长度自动调整字体大小
        if len(label_text) <= 3:
            font_size = 9
        elif len(label_text) == 4:
            font_size = 8
        else:
            font_size = 7
        
        ax.text(px, py, label_text, 
                ha='center', va='center', 
                fontsize=font_size, fontweight='bold', 
                color='#37474F', zorder=11, 
                fontfamily='sans-serif')
        
        # Draw Interactions (Improved Overlap Handling)
        # Group by target atom to separate overlapping lines
        target_groups = {}
        for inter in item["inters"]:
            tid = inter["target_idx"]
            if tid not in atom_coords: continue
            if tid not in target_groups: target_groups[tid] = []
            target_groups[tid].append(inter)
            
        for tid, inters in target_groups.items():
            tx, ty = atom_coords[tid]
            
            # Separate interactions into Lines vs Halos
            line_inters = []
            halo_inters = []
            
            for inter in inters:
                style = get_interaction_style(inter["type"])
                if not show_vdw and not style["show"]: continue
                if style["label"] == "Hydrophobic":
                    halo_inters.append(inter)
                else:
                    line_inters.append(inter)
            
            # 1. Draw Hydrophobic Halos (Stacking)
            for inter in halo_inters:
                style = get_interaction_style(inter["type"])
                import matplotlib.patches as patches
                hydro_radius = 14  # 14px radius
                hydro_circle = patches.Circle(
                    (tx, ty), hydro_radius,
                    fc=style["color"], ec='none', 
                    alpha=0.25, zorder=3 
                )
                ax.add_patch(hydro_circle)
                
            # 2. Draw Lines with Offset (Separating overlaps)
            n_lines = len(line_inters)
            for i, inter in enumerate(line_inters):
                style = get_interaction_style(inter["type"])
                
                # Calculate Offset to prevent overlap
                off_x, off_y = 0, 0
                if n_lines > 1:
                    vx, vy = tx - px, py - py
                    mag = math.sqrt(vx*vx + vy*vy)
                    if mag > 0.1:
                        ux, uy = -vy/mag, vx/mag # Perpendicular unit vector
                        spacing = 5.0 # 5 pixels separation
                        shift = (i - (n_lines - 1) / 2.0) * spacing
                        off_x = ux * shift
                        off_y = uy * shift
                
                # Apply Offset
                start_x, start_y = px + off_x, py + off_y
                end_x, end_y = tx + off_x, ty + off_y
                
                ls = style["style"]
                lw = style["linewidth"]
                lc = style["color"]
                line_alpha = 0.6 # Ensure visibility
                
                if style["arrow"]:
                     # Arrow with offset
                     ax.annotate("", xy=(end_x, end_y), xytext=(start_x, start_y),
                                arrowprops=dict(arrowstyle="->", color=lc, lw=lw, ls=ls, 
                                              shrinkA=5, shrinkB=5, alpha=line_alpha),
                                zorder=2) 
                else:
                    # Line with offset
                    ax.plot([start_x, end_x], [start_y, end_y], color=lc, lw=lw, ls=ls, zorder=2, alpha=line_alpha)
                    
    # Legend - 相互作用类型
    interaction_handles = []
    seen_interactions = set()
    order = ["Hbond", "Salt", "PiPi", "Pi", "Metal", "Halogen", "Water", "Hydrophobic"]
    
    for inter in interactions:
        s = get_interaction_style(inter["type"])
        l = s["label"]
        if not show_vdw and not s["show"]: continue
        if l not in seen_interactions:
            seen_interactions.add(l)
    
    for k in order:
        s = SCHRODINGER_STYLE[k]
        if s["label"] in seen_interactions:
            if s["label"] == "Hydrophobic":
                # 疏水作用用圆形标识
                interaction_handles.append(Line2D([0], [0], 
                                                   marker='o', color='w',
                                                   markerfacecolor=s["color"], 
                                                   markeredgecolor=s["color"],
                                                   markersize=8, alpha=0.5,
                                                   linestyle='None',
                                                   label=s["label"]))
            else:
                interaction_handles.append(Line2D([0], [0], color=s["color"], lw=2, 
                                     ls=s["style"] if s["style"] != "arc" else "-", 
                                     label=s["label"]))
    
    # Legend - 氨基酸分类 (英文标签)
    aa_legend_labels = {
        'nonpolar': 'Nonpolar',
        'polar': 'Polar',
        'negative': 'Negative',
        'positive': 'Positive'
    }
    
    # 检测实际出现的氨基酸类型
    seen_categories = set()
    for prot_res in residue_data.keys():
        _, _, cat = parse_residue_label(prot_res)
        if cat != 'unknown':
            seen_categories.add(cat)
    
    aa_handles = []
    for cat in ['nonpolar', 'polar', 'negative', 'positive']:
        if cat in seen_categories:
            aa_handles.append(mpatches.Patch(
                facecolor='white', 
                edgecolor=AA_CATEGORY_COLORS[cat], 
                linewidth=2,
                label=aa_legend_labels[cat]
            ))
    
    # 组合图例
    all_handles = interaction_handles + aa_handles
    if all_handles:
        ax.legend(handles=all_handles, loc='lower right', 
                  frameon=True, fancybox=True, framealpha=0.95, 
                  fontsize=9, ncol=1)

    # Title (已禁用)
    # t_txt = f"{protein_name} : {ligand_resname}" if protein_name else ligand_resname
    # ax.text(drawer_w/2, -80, t_txt, ha='center', fontsize=20, fontweight='bold', color='#263238', fontfamily='sans-serif')

    plt.tight_layout()
    if output_path is None: 
        output_path = f"{ligand_resname}_2d.png"
    
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight", facecolor='white')
    plt.close()
    
    print(f"[2D Diagram] ✅ Saved v2.0 Pro Diagram to {output_path}")
    return output_path

def release_itype(t):
    return str(t).lower()

try:
    from pymol import cmd
    cmd.extend("generate_2d_diagram", generate_2d_interaction_diagram)
except:
    pass
