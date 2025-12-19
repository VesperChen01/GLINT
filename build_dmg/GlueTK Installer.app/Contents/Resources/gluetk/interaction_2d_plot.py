# -*- coding: utf-8 -*-
"""
interaction_2d_plot.py
生成蛋白质-配体2D相互作用图 (Clean Visualization Version)

特性：
- 双分子策略：使用全原子模型映射相互作用，使用重原子模型进行绘图
- 自动将 H 原子的相互作用重映射到相邻的重原子，避免"线团"混乱
- 准确的 PDB原子名映射
- 优化的视觉样式
"""

from __future__ import print_function
import os
import csv
import tempfile
import math
import numpy as np

# 相互作用颜色方案 (Discovery Studio 风格 + GlueTK Consistency)
DS_STYLE = {
    "Hbond": {
        "color": "#43A047",      # Green (Classic DS style for 2D)
        "label": "Hydrogen Bond",
        "style": "--"
    },
    "Salt": {
        "color": "#F4511E",      # Deep Orange
        "label": "Salt Bridge",
        "style": "--"
    },
    "Pi": {
        "color": "#FB8C00",      # Orange
        "label": "Pi-Interaction",
        "style": "--"
    },
    "PiPi": {
        "color": "#FDD835",      # Yellow (Match 3D GlueTK)
        "label": "Pi-Pi Stacking",
        "style": "--"
    },
    "Hydrophobic": {
        "color": "#F06292",      # Pink
        "label": "Hydrophobic",
        "style": "--"
    },
    "Halogen": {
        "color": "#00ACC1",      # Cyan
        "label": "Halogen Bond",
        "style": "--"
    },
    "VDW": {
        "color": "#C5E1A5",      # Light Green
        "label": "van der Waals",
        "style": ":"
    }
}

def get_interaction_style(itype):
    """获取相互作用样式"""
    itype = itype.lower()
    if "氢键" in itype or "hbond" in itype or "hydrogen" in itype:
        return DS_STYLE["Hbond"]
    elif "盐桥" in itype or "salt" in itype or "ionic" in itype:
        return DS_STYLE["Salt"]
    elif "π-π" in itype or "pipi" in itype or "stacking" in itype:
        return DS_STYLE["PiPi"]
    elif "π" in itype or "pi" in itype or "cation" in itype:
        return DS_STYLE["Pi"]
    elif "卤素" in itype or "halogen" in itype:
        return DS_STYLE["Halogen"]
    elif "疏水" in itype or "hydrophobic" in itype or "alkyl" in itype:
        return DS_STYLE["Hydrophobic"]
    else:
        return DS_STYLE["VDW"]

def generate_2d_interaction_diagram(csv_path, ligand_resname, pdb_file=None, obj_name=None,
                                     output_path=None, width=1600, height=1200, dpi=150):
    """
    生成高精度且清晰的2D相互作用图
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
    
    # Force Explicit Hydrogen Removal for Drawing (Robust Method)
    try:
        # Try full sanitization first
        Chem.SanitizeMol(mol_draw)
        mol_draw = Chem.RemoveHs(mol_draw)
    except Exception as e:
        print(f"[2D Diagram] Warning: Full sanitization failed ({e}), attempting partial sanitization for drawing...")
        try:
             # Partial sanitization (skip properties that might fail)
             mol_draw.UpdatePropertyCache(strict=False)
             Chem.SanitizeMol(mol_draw, Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS, catchErrors=True)
             mol_draw = Chem.RemoveHs(mol_draw, implicitOnly=False)
        except Exception as e2:
             print(f"[2D Diagram] Critical: Failed to remove Hs even with partial sanitization: {e2}")
             # Last resort: Manually delete H atoms
             mw = Chem.RWMol(mol_draw)
             atoms_to_remove = [a.GetIdx() for a in mw.GetAtoms() if a.GetAtomicNum() == 1]
             atoms_to_remove.sort(reverse=True)
             for idx in atoms_to_remove:
                 mw.RemoveAtom(idx)
             mol_draw = mw.GetMol()

    
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
        AllChem.Compute2DCoords(mol_draw)
    except:
        pass

    # 绘图设置
    drawer_w, drawer_h = 1000, 750
    drawer = rdMolDraw2D.MolDraw2DCairo(drawer_w, drawer_h)
    
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.padding = 0.15      # 增加内边距
    opts.bondLineWidth = 3   # 加粗化学键
    opts.baseFontSize = 0.7  # 更大的原子标签 (相对缩放)
    opts.multipleBondOffset = 0.15
    
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

    # 5. Matplotlib 组装
    fig = plt.figure(figsize=(width/100, height/100), dpi=dpi)
    ax = fig.add_subplot(111)
    
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

    # 布局参数
    label_dist = 140.0 # 距离原子的长度
    bubble_r = 45.0
    drawn_bubbles = []

    print(f"[2D Diagram] Placing {len(residue_data)} residue bubbles...")

    for prot_res, inter_list in residue_data.items():
        # 计算锚点 (所有关联重原子的平均位置)
        coords = [atom_coords[i["target_idx"]] for i in inter_list if i["target_idx"] in atom_coords]
        if not coords:
            anchor_x, anchor_y = center_x, center_y # Fallback
        else:
            anchor_x = sum(c[0] for c in coords)/len(coords)
            anchor_y = sum(c[1] for c in coords)/len(coords)
        
        # 向量方向
        vx, vy = anchor_x - center_x, anchor_y - center_y
        norm = math.sqrt(vx**2 + vy**2)
        if norm < 0.1: vx, vy = 1, 0; norm=1
        dx, dy = vx/norm, vy/norm
        
        # 初始位置
        px = anchor_x + dx * label_dist
        py = anchor_y + dy * label_dist
        
        # 简单的斥力迭代 (防重叠)
        for _ in range(5):
            moved = False
            # 1. Bubble-Bubble Repulsion
            for (ox, oy, _) in drawn_bubbles:
                d = math.sqrt((px-ox)**2 + (py-oy)**2)
                min_d = bubble_r * 2.5
                if d < min_d:
                    rx, ry = px-ox, py-oy
                    rn = math.sqrt(rx**2+ry**2)
                    if rn < 0.1: rx, ry = 1,0; rn=1
                    px += (rx/rn)*40
                    py += (ry/rn)*40
                    moved = True
            
            # 2. Bubble-Atom Repulsion (Optional, avoid overlapping ligand)
            # 简单检查离中心太近
            d_center = math.sqrt((px-center_x)**2 + (py-center_y)**2)
            # 如果进入了配体区域 (假设半径 250)
            if d_center < 250:
                 cx, cy = px-center_x, py-center_y
                 cn = math.sqrt(cx**2+cy**2)
                 if cn<0.1: cx,cy=1,0; cn=1
                 px += (cx/cn)*20
                 py += (cy/cn)*20
                 moved=True
                 
            if not moved: break
        
        drawn_bubbles.append((px, py, bubble_r))
        
        # 确定优先级最高的颜色
        best_style = DS_STYLE["VDW"]
        w_max = 0
        for i in inter_list:
            s_name = i["type"]
            w = 0
            if "氢键" in s_name or "Hbond" in s_name: w=10
            elif "盐桥" in s_name or "Salt" in s_name: w=9
            elif "π-π" in s_name or "PiPi" in s_name: w=8
            elif "π" in s_name or "Pi" in s_name: w=7
            
            if w > w_max:
                w_max = w
                best_style = get_interaction_style(s_name)
        
        # 绘制气泡
        # 名字分行
        parts = prot_res.split()
        if len(parts) >= 2: txt = f"{parts[0][:3]}\n{parts[1]}"
        else: txt = prot_res[:3]
        
        circle = Circle((px, py), bubble_r, facecolor='white', edgecolor=best_style["color"], linewidth=2, zorder=10)
        ax.add_patch(circle)
        ax.text(px, py, txt, ha='center', va='center', fontsize=10, fontweight='bold', zorder=11)
        
        # 绘制连线
        for inter in inter_list:
            if inter["target_idx"] not in atom_coords: continue
            tx, ty = atom_coords[inter["target_idx"]]
            
            # 计算连线切点
            lx, ly = px-tx, py-ty
            ld = math.sqrt(lx**2+ly**2)
            if ld > bubble_r:
                ex = tx + (lx/ld)*(ld-bubble_r)
                ey = ty + (ly/ld)*(ld-bubble_r)
            else:
                ex, ey = px, py
            
            style = get_interaction_style(inter["type"])
            line = Line2D([tx, ex], [ty, ey], 
                         color=style["color"], 
                         linestyle=style["style"], 
                         linewidth=1.8, # Thicker lines
                         alpha=0.85, 
                         zorder=5)
            ax.add_line(line)
            
            # Distance label
            d_val = inter["distance"]
            if d_val and d_val not in ["-", ""]:
                mx, my = (tx+ex)/2, (ty+ey)/2
                ax.text(mx, my, str(d_val), fontsize=8, color='black', 
                       bbox=dict(boxstyle='round,pad=0.1', fc='white', ec='none', alpha=0.9),
                       ha='center', va='center', zorder=6)

    # Legend
    handles = []
    seen = set()
    sort_order = ["Hbond", "Salt", "PiPi", "Pi", "Halogen", "Hydrophobic", "VDW"]
    
    for k in sort_order:
        s = DS_STYLE[k]
        if s["label"] not in seen:
            handles.append(mpatches.Patch(color=s["color"], label=s["label"]))
            seen.add(s["label"])
            
    ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.15, 1), fontsize=9, title="Interaction Types")
    ax.text(drawer_w/2, 40, f"{ligand_resname} Interactions", ha='center', fontsize=16, fontweight='bold')

    plt.tight_layout()
    if output_path is None: output_path = f"{ligand_resname}_2d.png"
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    
    print(f"[2D Diagram] ✅ Saved standard clean diagram to {output_path}")
    return output_path

try:
    from pymol import cmd
    cmd.extend("generate_2d_diagram", generate_2d_interaction_diagram)
except:
    pass
