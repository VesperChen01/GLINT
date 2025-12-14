# -*- coding: utf-8 -*-
"""
interaction_2d_plot.py
生成蛋白质-配体2D相互作用图 (Schrödinger Style)

特性：
- Schrödinger Maestro 风格的现代简洁设计
- 双分子策略：使用全原子模型映射相互作用，使用重原子模型进行绘图
- 自动将 H 原子的相互作用重映射到相邻的重原子，避免"线团"混乱
- 准确的 PDB原子名映射
- 环形紧凑布局，智能避免重叠
- 氢键箭头指示供受体方向
- 疏水接触使用灰色弧形表示
"""

from __future__ import print_function
import os
import csv
import tempfile
import math
import numpy as np

# GlueTK 统一配色方案 (与 3D 视图保持一致)
SCHRODINGER_STYLE = {
    "Hbond": {
        "color": "#2196F3",      # 蓝色 - 氢键
        "label": "H-bond",
        "style": "-",           
        "arrow": True,
        "linewidth": 2.5,
        "show": True
    },
    "Salt": {
        "color": "#FF5722",      # 橙色 - 盐桥
        "label": "Salt Bridge",
        "style": "-",
        "arrow": True,
        "linewidth": 2.5,
        "show": True
    },
    "PiPi": {
        "color": "#9C27B0",      # 紫色 - π-π 堆积
        "label": "π-π Stack",
        "style": "-",
        "arrow": False,
        "linewidth": 2.5,
        "show": True
    },
    "Pi": {
        "color": "#E91E63",      # 粉红 - π-阳离子
        "label": "π-Cation",
        "style": "-",
        "arrow": False,
        "linewidth": 2.5,
        "show": True
    },
    "Hydrophobic": {
        "color": "#4CAF50",      # 绿色 - 疏水相互作用
        "label": "Hydrophobic",
        "style": "arc",
        "arrow": False,
        "linewidth": 1.5,
        "show": True
    },
    "Halogen": {
        "color": "#FF9800",      # 橙黄 - 卤素键
        "label": "Halogen",
        "style": "-",
        "arrow": True,
        "linewidth": 2.5,
        "show": True
    },
    "Metal": {
        "color": "#673AB7",      # 深紫 - 金属配位
        "label": "Metal",
        "style": "-",
        "arrow": False,
        "linewidth": 2.5,
        "show": True
    },
    "Water": {
        "color": "#00BCD4",      # 青色 - 水桥
        "label": "Water Bridge",
        "style": "--",
        "arrow": False,
        "linewidth": 2.0,
        "show": True
    },
    "VDW": {
        "color": "#BDBDBD",      # 浅灰 - 范德华力
        "label": "vdW",
        "style": ":",
        "arrow": False,
        "linewidth": 1.0,
        "show": False
    }
}

# 保留旧的配色方案别名以兼容
DS_STYLE = SCHRODINGER_STYLE

def get_interaction_style(itype):
    """获取相互作用样式"""
    itype = itype.lower()
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
                                     output_path=None, width=1600, height=1200, dpi=150,
                                     style="schrodinger", show_distance=False, show_vdw=False,
                                     protein_name=None, compact=True):
    """
    生成高精度且清晰的 2D 相互作用图 (Schrödinger 风格)
    
    参数:
        csv_path: 相互作用 CSV 文件路径
        ligand_resname: 配体残基名称
        pdb_file: PDB 文件路径（可选）
        obj_name: PyMOL 对象名（可选）
        output_path: 输出图片路径
        width, height, dpi: 图像尺寸参数
        style: 可选 "schrodinger" (默认，简洁现代) 或 "discovery" (经典)
        show_distance: 是否显示距离标注 (默认 False - Schrödinger 风格)
        show_vdw: 是否显示范德华相互作用 (默认 False)
        protein_name: 蛋白名称（可选，用于标题）
        compact: 是否使用紧凑布局 (默认 True)
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
        from matplotlib.patches import Circle
        from matplotlib.lines import Line2D
        from PIL import Image
    except ImportError:
        print("[2D Diagram] matplotlib and Pillow are required")
        return None

    print(f"[2D Diagram] Starting Clean Generation for {ligand_resname}...")

    # 1. 提取配体结构 (PDB格式)
    # 我们需要加载两次:
    # mol_map: 带 H (removeHs=False)，用于通过 PDB 原子名查找原子
    # mol_draw: 不带 H (removeHs=True)，用于清晰绘图
    
    mol_map = None
    mol_draw = None
    temp_pdb = None
    
    if obj_name:
        try:
            from pymol import cmd
            temp_pdb = tempfile.mktemp(suffix=".pdb")
            # 通过 PyMOL 保存 PDB
            cmd.save(temp_pdb, f"{obj_name} and resn {ligand_resname}", format="pdb")
            
            if os.path.exists(temp_pdb):
                # 1. Mapping Molecule (Explicit H for mapping)
                mol_map = Chem.MolFromPDBFile(temp_pdb, removeHs=False, sanitize=False)
                
                # 2. Drawing Molecule (Implicit H for clean look)
                mol_draw = Chem.MolFromPDBFile(temp_pdb, removeHs=True, sanitize=False)
                
                try:
                    Chem.SanitizeMol(mol_map)
                    Chem.SanitizeMol(mol_draw)
                except:
                    pass
        except Exception as e:
            print(f"[2D Diagram] Failed to extract ligand from PyMOL: {e}")
            if temp_pdb and os.path.exists(temp_pdb): os.remove(temp_pdb)
            return None
    
    # 清理临时文件
    if temp_pdb and os.path.exists(temp_pdb):
        os.remove(temp_pdb)

    if mol_map is None or mol_draw is None:
        print(f"[2D Diagram] ❌ Failed to load ligand molecule.")
        return None
    
    # 2. 建立原子映射
    # Step A: 建立 Name -> Atom Object (in mol_map)
    map_name_to_atom = {}
    for atom in mol_map.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info:
            name = info.GetName().strip()
            map_name_to_atom[name] = atom

    # Step B: 建立 Name -> Index (in mol_draw)
    # 我们的目标是找到 mol_draw 中的对应原子索引
    # mol_draw 中的原子都是重原子。
    draw_name_to_idx = {}
    for atom in mol_draw.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info:
            name = info.GetName().strip()
            draw_name_to_idx[name] = atom.GetIdx()
            
    print(f"[2D Diagram] Mapping prepared. Heavy atoms in drawing: {len(draw_name_to_idx)}")

    # 3. 读取相互作用数据并关联到 mol_draw 索引
    interactions = []
    
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            lig_atom_name = row.get("Ligand_Atom", "").strip()
            prot_res = row.get("Protein_Residue", "").strip()
            itype = row.get("Interaction", "").strip()
            dist = row.get("Distance", "").strip()
            
            if not prot_res: continue

            # --- Resolving Target Atom Index in Clean Molecule ---
            target_idx = None
            
            # Case 1: 直接是重原子，且存在于 mol_draw
            if lig_atom_name in draw_name_to_idx:
                target_idx = draw_name_to_idx[lig_atom_name]
            else:
                # Case 2: 可能是 H 原子，或者名字没匹配上
                # 在 mol_map 中找这个原子
                if lig_atom_name in map_name_to_atom:
                    map_atom = map_name_to_atom[lig_atom_name]
                    
                    # 检查是否为氢原子 (AtomicNum = 1)
                    if map_atom.GetAtomicNum() == 1:
                        # 找到相连的重原子 neighbor
                        neighbors = map_atom.GetNeighbors()
                        if neighbors:
                            neighbor = neighbors[0] # H 只有一个邻居
                            neighbor_info = neighbor.GetPDBResidueInfo()
                            if neighbor_info:
                                neighbor_name = neighbor_info.GetName().strip()
                                # 尝试在 draw mol 中找这个重原子
                                if neighbor_name in draw_name_to_idx:
                                    target_idx = draw_name_to_idx[neighbor_name]
            
            # 如果还是找不到 (e.g. Ring interaction 标记为 Center)，暂忽略或映射到中心
            
            if target_idx is not None:
                interactions.append({
                    "protein": prot_res,
                    "target_idx": target_idx,
                    "type": itype,
                    "distance": dist
                })

    print(f"[2D Diagram] Loaded {len(interactions)} interactions mapped to heavy atoms.")

    # 4. 生成 2D 坐标 (Clean Molecule)
    try:
        # 使用 rdCoordGen 获得更好的 2D 拓扑结构 (特别是对于大环)
        from rdkit.Chem import rdCoordGen
        rdCoordGen.AddCoords(mol_draw)
    except:
        try:
            AllChem.Compute2DCoords(mol_draw)
        except:
            pass
            
    # 计算分子大小，自适应调整参数
    num_atoms = mol_draw.GetNumAtoms()
    is_large_molecule = num_atoms > 50
    
    # 绘图设置 - 动态调整
    drawer_w, drawer_h = (1200, 900) if is_large_molecule else (1000, 750)
    drawer = rdMolDraw2D.MolDraw2DCairo(drawer_w, drawer_h)
    
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.padding = 0.20 if is_large_molecule else 0.25
    
    # 大分子参数调整：线条更细，字体更小
    if is_large_molecule:
        opts.bondLineWidth = 2.5
        opts.baseFontSize = 0.6
        opts.multipleBondOffset = 0.15
        opts.fixedBondLength = 25
        opts.minFontSize = 10
        opts.annotationFontScale = 0.6
    else:
        opts.bondLineWidth = 4.5
        opts.baseFontSize = 0.9
        opts.multipleBondOffset = 0.20
        opts.fixedBondLength = 32
        opts.minFontSize = 14
        opts.annotationFontScale = 0.8

    # 绘制
    drawer.DrawMolecule(mol_draw)
    drawer.FinishDrawing()
    
    png_data = drawer.GetDrawingText()
    
    # 获取坐标
    atom_coords = {}
    for i in range(mol_draw.GetNumAtoms()):
        pt = drawer.GetDrawCoords(i)
        atom_coords[i] = (pt.x, pt.y)

    # 计算中心
    xs = [p[0] for p in atom_coords.values()]
    ys = [p[1] for p in atom_coords.values()]
    center_x = sum(xs) / len(xs) if xs else drawer_w/2
    center_y = sum(ys) / len(ys) if ys else drawer_h/2

    # 辅助函数：绘制水滴形 (Teardrop) 残基
    def create_teardrop_patch(center, radius, angle, facecolor, edgecolor, alpha=1.0):
        cx, cy = center
        tip_dist = radius * 1.35  # 尖端轻微突出的距离
        
        # 路径顶点生成
        verts = []
        codes = []
        
        # 1. 尖端 (Tip)
        tx = cx + tip_dist * math.cos(angle)
        ty = cy + tip_dist * math.sin(angle)
        verts.append((tx, ty))
        codes.append(Path.MOVETO)
        
        # 2. 侧翼角度 (Wing Angle)
        # 从尖端向两侧展开，连接到圆的切点
        # 我们保留背部约 220 度的圆弧
        wing_angle_offset = math.radians(65) # 偏离中心轴的角度
        
        # 生成圆弧点 (背部)
        # 从 angle + wing_offset 开始，绕一圈到 angle - wing_offset
        start_angle = angle + wing_angle_offset
        end_angle = angle - wing_angle_offset + 2 * math.pi
        
        num_arc_points = 24
        for i in range(num_arc_points + 1):
            theta = start_angle + (end_angle - start_angle) * i / num_arc_points
            px = cx + radius * math.cos(theta)
            py = cy + radius * math.sin(theta)
            verts.append((px, py))
            codes.append(Path.LINETO)
            
        # 闭合回到尖端
        verts.append((tx, ty))
        codes.append(Path.LINETO)
        
        path = Path(verts, codes)
        return mpatches.PathPatch(path, facecolor=facecolor, edgecolor=edgecolor, 
                                 linewidth=2.0, alpha=alpha, zorder=10, capstyle='round')

    # 5. Matplotlib 组装
    fig = plt.figure(figsize=(width/100, height/100), dpi=dpi)
    ax = fig.add_subplot(111)
    
    # 必须导入 Path
    from matplotlib.path import Path
    
    tmp_png = tempfile.mktemp(suffix=".png")
    with open(tmp_png, "wb") as f:
        f.write(png_data)
    lig_img = Image.open(tmp_png)
    
    ax.imshow(lig_img, origin='upper', extent=[0, drawer_w, 0, drawer_h])
    os.remove(tmp_png)
    
    ax.set_xlim(0, drawer_w)
    ax.set_ylim(drawer_h, 0)
    ax.axis('off')
    
    # 分组相互作用
    residue_data = {}
    for inter in interactions:
        p = inter["protein"]
        if p not in residue_data: residue_data[p] = []
        residue_data[p].append(inter)

    # 布局参数 - Schrödinger 风格：紧凑、环形
    # 大分子需要更大的布局半径
    base_dist = 110.0 if is_large_molecule else 90.0
    label_dist = base_dist if compact else base_dist + 20
    
    # 大分子使用更小的气泡
    base_r = 24.0 if is_large_molecule else 28.0
    bubble_r = base_r if compact else base_r + 4.0
    
    drawn_bubbles = []

    print(f"[2D Diagram] Placing {len(residue_data)} residue bubbles...")
    
    # 【优化】使用极坐标均匀环形布局 - Schrödinger 风格
    # 先按角度排序残基，避免连线交叉
    residue_angles = []
    for prot_res, inter_list in residue_data.items():
        coords = [atom_coords[i["target_idx"]] for i in inter_list if i["target_idx"] in atom_coords]
        if coords:
            anc_x = sum(c[0] for c in coords) / len(coords)
            anc_y = sum(c[1] for c in coords) / len(coords)
            angle = math.atan2(anc_y - center_y, anc_x - center_x)
        else:
            angle = 0
            anc_x, anc_y = center_x, center_y # fallback
            
        residue_angles.append({
            "res": prot_res, 
            "inters": inter_list, 
            "angle": angle,
            "anchor": (anc_x, anc_y)
        })
    
    # 按角度排序
    residue_angles.sort(key=lambda x: x["angle"])
    num_residues = len(residue_angles)
    
    for idx, item in enumerate(residue_angles):
        prot_res = item["res"]
        inter_list = item["inters"]
        base_angle = item["angle"]
        anchor_x, anchor_y = item["anchor"]
        
        # 向量方向（从中心指向锚点）
        vx, vy = anchor_x - center_x, anchor_y - center_y
        norm = math.sqrt(vx**2 + vy**2)
        if norm < 0.1: vx, vy = 1, 0; norm = 1
        
        # 均匀分布角度调整
        if num_residues > 1:
            # 计算相邻残基的角度，避免挤在一起
            angle_step = 2 * math.pi / max(num_residues, 10)
            min_angle_diff = angle_step * 0.5
            
            # 检查与已放置气泡的角度冲突
            for (ox, oy, _) in drawn_bubbles:
                other_angle = math.atan2(oy - center_y, ox - center_x)
                angle_diff = abs(base_angle - other_angle)
                if angle_diff > math.pi:
                    angle_diff = 2 * math.pi - angle_diff
                if angle_diff < min_angle_diff:
                    # 轻微偏移，交替方向
                    shift = min_angle_diff * (1 if idx % 2 == 0 else -1)
                    base_angle += shift
            
            dx = math.cos(base_angle)
            dy = math.sin(base_angle)
        else:
            dx, dy = vx/norm, vy/norm
        
        # 初始位置 - 基于锚点距离动态调整
        # 大分子的标签距离需要更远一点，避免压住分子
        dist_factor = 0.35 if is_large_molecule else 0.25
        effective_dist = label_dist + min(norm * dist_factor, 60)
        px = center_x + dx * effective_dist
        py = center_y + dy * effective_dist
        
        # 简化的斥力迭代
        for iteration in range(8): # 增加迭代次数
            moved = False
            
            # Bubble-Bubble 斥力
            for (ox, oy, _) in drawn_bubbles:
                d = math.sqrt((px-ox)**2 + (py-oy)**2)
                min_d = bubble_r * 2.5 # 增加间距
                if d < min_d and d > 0.1:
                    rx, ry = px-ox, py-oy
                    rn = math.sqrt(rx**2+ry**2)
                    push = (min_d - d) * 0.5
                    px += (rx/rn) * push
                    py += (ry/rn) * push
                    moved = True
            
            # 避免与配体重叠（配体区域保护）
            d_center = math.sqrt((px-center_x)**2 + (py-center_y)**2)
            # 保护半径动态计算
            min_center_dist = (drawer_w * 0.22) if is_large_molecule else 180
            if compact: min_center_dist *= 0.9
            
            if d_center < min_center_dist:
                cx, cy = px-center_x, py-center_y
                cn = math.sqrt(cx**2+cy**2)
                if cn > 0.1:
                    px += (cx/cn) * 15
                    py += (cy/cn) * 15
                    moved = True
            
            if not moved:
                break
        
        drawn_bubbles.append((px, py, bubble_r))
        
        # 确定优先级最高的颜色
        best_style = SCHRODINGER_STYLE["VDW"]
        w_max = 0
        priority_map = {
            "氢键": 10, "hbond": 10, "hydrogen": 10,
            "盐桥": 9, "salt": 9, "ionic": 9,
            "π-π": 8, "pipi": 8, "stacking": 8,
            "π": 7, "pi": 7, "cation": 7,
            "金属": 6, "metal": 6,
            "卤素": 5, "halogen": 5,
            "水桥": 4, "water": 4,
            "疏水": 3, "hydrophobic": 3,
        }
        
        for i in inter_list:
            s_name = i["type"].lower()
            w = 0
            for key, priority in priority_map.items():
                if key in s_name:
                    w = max(w, priority)
                    break
            
            if w > w_max:
                w_max = w
                best_style = get_interaction_style(i["type"])
        
        # 名字简化：只显示三字母代码 + 残基编号
        parts = prot_res.split()
        if len(parts) >= 2:
            res_name = parts[0][:3].upper()
            res_num = parts[1]
            txt = f"{res_name}{res_num}"
        else:
            txt = prot_res[:7]
        
        # 使用淡化的相互作用颜色填充气泡
        import matplotlib.colors as mcolors
        fill_color = mcolors.to_rgba(best_style["color"], alpha=0.15)
        edge_color = mcolors.to_rgba(best_style["color"], alpha=0.95)
        
        # 【新】绘制 Teardrop (吉他拨片) 形状
        # 计算拨片指向：从气泡中心指向锚点（或者配体中心）
        # 这里指向锚点 (anchor_x, anchor_y) 会更准确
        aim_dx = anchor_x - px
        aim_dy = anchor_y - py
        aim_angle = math.atan2(aim_dy, aim_dx)
        
        patch = create_teardrop_patch(
            (px, py), bubble_r, aim_angle, 
            facecolor=fill_color, edgecolor=edge_color
        )
        ax.add_patch(patch)
        
        # 标签字体优化
        text_color = '#263238' # 深蓝灰
        font_weight = 'bold'
        font_size = 9 if is_large_molecule else 9.5
        
        ax.text(px, py, txt, ha='center', va='center', 
               fontsize=font_size, fontweight=font_weight, 
               color=text_color, zorder=11, family='sans-serif')
        
        # 【改进】绘制连线 - 支持箭头、弧形
        for inter in inter_list:
            if inter["target_idx"] not in atom_coords: continue
            
            line_style = get_interaction_style(inter["type"])
            
            # 跳过不显示的相互作用类型（如 VDW）
            if not show_vdw and not line_style.get("show", True):
                continue
            
            tx, ty = atom_coords[inter["target_idx"]]
            
            # 计算连线切点（从原子到气泡边缘）
            lx, ly = px-tx, py-ty
            ld = math.sqrt(lx**2+ly**2)
            if ld > bubble_r:
                ex = tx + (lx/ld)*(ld-bubble_r)
                ey = ty + (ly/ld)*(ld-bubble_r)
            else:
                ex, ey = px, py
            
            linewidth = line_style.get("linewidth", 2.0)
            
            # 疏水相互作用使用灰色弧形
            if line_style.get("style") == "arc":
                import matplotlib.patches as patches
                # 绘制简洁的扇形弧（表示疏水接触区域）
                arc_angle = 20  # 弧的角度范围（更窄更简洁）
                center_angle = math.degrees(math.atan2(ey-ty, ex-tx))
                
                wedge = patches.Wedge(
                    (tx, ty), ld*0.55, 
                    center_angle - arc_angle/2, 
                    center_angle + arc_angle/2,
                    facecolor=line_style["color"], 
                    alpha=0.12, 
                    edgecolor=line_style["color"],
                    linewidth=linewidth,
                    zorder=4
                )
                ax.add_patch(wedge)
            else:
                # 常规连线
                ls = line_style.get("style", "-")
                
                # 氢键和盐桥添加箭头（从配体指向残基）
                if line_style.get("arrow", False):
                    from matplotlib.patches import FancyArrowPatch
                    arrow = FancyArrowPatch(
                        (tx, ty), (ex, ey),
                        arrowstyle='->,head_width=0.35,head_length=0.5',
                        color=line_style["color"],
                        linewidth=linewidth,
                        linestyle=ls,
                        alpha=0.9,
                        zorder=5,
                        mutation_scale=12
                    )
                    ax.add_patch(arrow)
                else:
                    # 普通线条
                    line = Line2D([tx, ex], [ty, ey], 
                                 color=line_style["color"], 
                                 linestyle=ls, 
                                 linewidth=linewidth, 
                                 alpha=0.85, 
                                 zorder=5,
                                 solid_capstyle='round')
                    ax.add_line(line)
            
            # 【简化】移除距离标注（Schrödinger 风格默认不显示）
            if show_distance:
                d_val = inter["distance"]
                if d_val and d_val not in ["-", ""]:
                    mx, my = (tx+ex)/2, (ty+ey)/2
                    ax.text(mx, my, f"{d_val}Å", fontsize=7, color='#424242', 
                           bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='#BDBDBD', alpha=0.92, linewidth=0.5),
                           ha='center', va='center', zorder=6)

    # 【改进】简化图例 - 只显示实际出现的相互作用类型
    handles = []
    seen_types = set()
    
    # 收集实际使用的相互作用类型（保持一致的顺序）
    type_order = ["Hbond", "Salt", "PiPi", "Pi", "Metal", "Halogen", "Water", "Hydrophobic", "VDW"]
    type_handles = {}
    
    for inter in interactions:
        itype = inter["type"]
        istyle = get_interaction_style(itype)
        label = istyle["label"]
        
        # 跳过 VDW（除非明确要求显示）
        if not show_vdw and label == "vdW":
            continue
            
        if label not in seen_types:
            seen_types.add(label)
            type_handles[label] = mpatches.Patch(color=istyle["color"], label=label)
    
    # 按预定义顺序排列图例
    for t in type_order:
        if t in SCHRODINGER_STYLE:
            label = SCHRODINGER_STYLE[t]["label"]
            if label in type_handles:
                handles.append(type_handles[label])
    
    # 图例放在右下角，更紧凑
    if handles:
        legend = ax.legend(handles=handles, loc='lower right', 
                          bbox_to_anchor=(1.0, 0.0), 
                          fontsize=8, 
                          framealpha=0.92,
                          edgecolor='#BDBDBD',
                          fancybox=True,
                          shadow=False,
                          borderpad=0.6,
                          handlelength=1.2,
                          handleheight=0.8)
        legend.get_frame().set_linewidth(0.8)
    
    # 标题样式优化
    title_text = ligand_resname
    if protein_name:
        title_text = f"{protein_name} : {ligand_resname}"
    
    ax.text(drawer_w/2, 25, title_text, 
           ha='center', fontsize=13, fontweight='bold', 
           color='#212121', family='sans-serif')

    plt.tight_layout()
    if output_path is None: 
        output_path = f"{ligand_resname}_2d.png"
    
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight", 
                facecolor='white', edgecolor='none')
    plt.close()
    
    print(f"[2D Diagram] ✅ Saved Schrödinger-style diagram to {output_path}")
    return output_path

try:
    from pymol import cmd
    cmd.extend("generate_2d_diagram", generate_2d_interaction_diagram)
except:
    pass
