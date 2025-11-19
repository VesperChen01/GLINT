# -*- coding: utf-8 -*-
"""
interaction_2d_plot.py
生成蛋白质-配体2D相互作用图

增强版特性：
- 配体化学结构式
- 蛋白质残基结构式（氨基酸侧链）
- 按类别着色（疏水、极性、带电等）
- 标注具体相互作用原子
"""

from __future__ import print_function
import os
import csv
import tempfile
import math


# 相互作用颜色方案 (Discovery Studio 风格)
DS_STYLE = {
    "Hbond": {
        "color": "#43A047",      # Green
        "label": "Hydrogen Bond",
        "style": "--"
    },
    "Salt": {
        "color": "#F4511E",      # Deep Orange
        "label": "Salt Bridge",
        "style": "--"
    },
    "Pi": {
        "color": "#FB8C00",      # Orange (Pi-Cation)
        "label": "Pi-Interaction",
        "style": "--"
    },
    "PiPi": {
        "color": "#8E24AA",      # Purple
        "label": "Pi-Pi Stacking",
        "style": "--"
    },
    "Hydrophobic": {
        "color": "#F06292",      # Pink (Alkyl/Pi-Alkyl)
        "label": "Hydrophobic/Alkyl",
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
    生成2D相互作用图（Discovery Studio 风格）
    
    特性：
    - 残基显示为圆形气泡
    - 颜色根据相互作用类型编码（H键绿色，烷基/疏水粉色，盐桥橙色）
    - 虚线连接
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
        from matplotlib.patches import Circle, FancyArrowPatch, Rectangle
        from matplotlib.lines import Line2D
        from PIL import Image
    except ImportError:
        print("[2D Diagram] matplotlib and Pillow are required")
        return None

    print(f"[2D Diagram] 2.0 Loaded - Discovery Studio Style")
    print(f"[2D Diagram] Reading interaction data: {csv_path}")
    
    # 读取相互作用数据
    interactions = []
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if "Ligand_Residue" in row:
                lig_res = row.get("Ligand_Residue", "").strip()
                prot_res = row.get("Protein_Residue", "").strip()
                lig_atom = row.get("Ligand_Atom", "").strip()
                prot_atom = row.get("Protein_Atom", "").strip()
            else:
                res1 = row.get("Residue1", "").strip()
                res2 = row.get("Residue2", "").strip()
                if ligand_resname.upper() in res1.upper():
                    lig_res, prot_res = res1, res2
                    lig_atom, prot_atom = "", ""
                elif ligand_resname.upper() in res2.upper():
                    lig_res, prot_res = res2, res1
                    lig_atom, prot_atom = "", ""
                else:
                    continue
            
            interactions.append({
                "ligand": lig_res,
                "protein": prot_res,
                "ligand_atom": lig_atom,
                "protein_atom": prot_atom,
                "type": row.get("Interaction", "").strip(),
                "distance": row.get("Distance", "").strip()
            })

    if not interactions:
        print(f"[2D Diagram] No interactions found for ligand '{ligand_resname}'")
        return None

    print(f"[2D Diagram] Found {len(interactions)} interactions")

    # 提取配体结构
    mol = None
    if obj_name:
        try:
            from pymol import cmd
            temp_sdf = tempfile.mktemp(suffix=".sdf")
            cmd.save(temp_sdf, f"{obj_name} and resn {ligand_resname}", format="sdf")
            if os.path.exists(temp_sdf):
                mol = Chem.SDMolSupplier(temp_sdf, removeHs=True)[0] # 移除氢原子，图更清晰
                os.remove(temp_sdf)
        except Exception as e:
            print(f"[2D Diagram] Failed to extract ligand from PyMOL: {e}")

    if mol is None:
        print(f"[2D Diagram] Unable to extract ligand structure")
        return None

    # 生成2D坐标
    if not mol.GetNumConformers():
        AllChem.Compute2DCoords(mol)
    else:
        # 尝试保留原有构象投影到2D，或者重新生成
        try:
            AllChem.GenerateDepictionMatching3DStructure(mol, mol)
        except:
            AllChem.Compute2DCoords(mol)

    # 创建画布
    fig = plt.figure(figsize=(width/100, height/100), dpi=dpi)
    ax = fig.add_subplot(111)
    # 扩大视野以容纳周围的残基
    ax.set_xlim(-10, 10)
    ax.set_ylim(-8, 8)
    ax.axis('off')

    # 绘制配体结构
    print(f"[2D Diagram] Drawing ligand structure...")
    # 使用RDKit绘制配体，背景透明
    dopts = rdMolDraw2D.MolDrawOptions()
    dopts.clearBackground = False
    dopts.fixedBondLength = 40
    
    drawer = rdMolDraw2D.MolDraw2DCairo(800, 600)
    drawer.SetDrawOptions(dopts)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    temp_img = tempfile.mktemp(suffix=".png")
    with open(temp_img, "wb") as f:
        f.write(drawer.GetDrawingText())

    ligand_img = Image.open(temp_img)
    # 调整配体显示大小和位置
    ax.imshow(ligand_img, extent=[-4, 4, -3, 3], zorder=10)
    os.remove(temp_img)

    # 整理残基相互作用
    residue_interactions = {}
    residue_styles = {} # 记录每个残基的主要样式（颜色）
    
    # 优先级：Hbond > Salt > Pi > Hydrophobic
    style_priority = {
        "Hbond": 10, "Salt": 9, "PiPi": 8, "Pi": 7, "Halogen": 6, "Hydrophobic": 5, "VDW": 1
    }

    for inter in interactions:
        prot_res = inter["protein"].strip()
        if not prot_res: continue
        
        # 提取残基名和编号 (e.g. "VAL 103")
        # 假设格式 "RES ID" 或 "RES:ID"
        parts = prot_res.split()
        if len(parts) >= 2:
            res_name = parts[0]
            res_id = parts[1]
        else:
            res_name = prot_res
            res_id = "?"
            
        unique_key = prot_res # 使用完整字符串作为键
        
        if unique_key not in residue_interactions:
            residue_interactions[unique_key] = []
            residue_styles[unique_key] = {"priority": 0, "style": DS_STYLE["VDW"]}
            
        residue_interactions[unique_key].append(inter)
        
        # 更新残基颜色样式
        style = get_interaction_style(inter["type"])
        # 反向查找 key
        style_key = "VDW"
        for k, v in DS_STYLE.items():
            if v == style:
                style_key = k
                break
        
        prio = style_priority.get(style_key, 0)
        if prio > residue_styles[unique_key]["priority"]:
            residue_styles[unique_key] = {"priority": prio, "style": style}

    residues = list(residue_interactions.keys())
    n_residues = len(residues)
    radius = 6.0 # 布局半径

    print(f"[2D Diagram] Drawing {n_residues} interacting residues (Discovery Studio style)...")

    # 绘制残基（圆形布局）
    for i, res_key in enumerate(residues):
        angle = 2 * math.pi * i / n_residues - math.pi/2
        # 稍微错开半径，避免太整齐
        r_offset = 0.2 * (i % 2)
        x = (radius + r_offset) * math.cos(angle)
        y = (radius + r_offset) * math.sin(angle) * 0.8 # 压扁一点椭圆

        # 获取样式
        style = residue_styles[res_key]["style"]
        bubble_color = style["color"]
        
        # 解析显示名称
        parts = res_key.split()
        if len(parts) >= 2:
            res_n = parts[0][:3].upper()
            res_i = parts[1]
            display_text = f"{res_n}\n{res_i}"
        else:
            display_text = res_key[:3]

        # 1. 绘制圆形气泡
        circle = Circle((x, y), 0.9, facecolor='white',
                      edgecolor=bubble_color, linewidth=2.5, alpha=1.0, zorder=20)
        ax.add_patch(circle)
        
        # 2. 内部填充（淡色）
        circle_fill = Circle((x, y), 0.9, facecolor=bubble_color,
                           alpha=0.15, zorder=19)
        ax.add_patch(circle_fill)

        # 3. 文字标签
        ax.text(x, y, display_text, fontsize=11, fontweight='bold',
               ha='center', va='center', color='black', zorder=21,
               multialignment='center')

        # 4. 绘制虚线连接
        interactions_list = residue_interactions[res_key]
        
        # 计算连线起始点（配体侧）
        # 稍微向外一点，避免穿过配体主体
        lig_x = 3.0 * math.cos(angle)
        lig_y = 2.4 * math.sin(angle)
        
        for inter in interactions_list:
            istyle = get_interaction_style(inter["type"])
            line_color = istyle["color"]
            dash_style = istyle["style"] # '--' or ':'
            
            # 绘制虚线
            # 终点是圆圈边缘而不是中心
            # 简单的向量计算
            dx, dy = x - lig_x, y - lig_y
            dist_len = math.sqrt(dx*dx + dy*dy)
            if dist_len > 0:
                # 缩短终点，停在气泡边缘 (半径 ~0.9)
                end_x = x - (dx / dist_len) * 0.9
                end_y = y - (dy / dist_len) * 0.9
            else:
                end_x, end_y = x, y

            line = Line2D([lig_x, end_x], [lig_y, end_y],
                         color=line_color,
                         linewidth=1.8,
                         linestyle=dash_style,
                         alpha=0.8,
                         zorder=5)
            ax.add_line(line)
            
            # 距离标签
            dist = inter["distance"]
            if dist and dist != "-" and dist != "0.0":
                mid_x = (lig_x + end_x) / 2
                mid_y = (lig_y + end_y) / 2
                ax.text(mid_x, mid_y, f"{dist}",
                       fontsize=9, ha='center', va='center',
                       bbox=dict(boxstyle='round,pad=0.15',
                               facecolor='white', edgecolor=line_color, linewidth=1, alpha=0.9),
                       color=line_color, fontweight='bold', zorder=6)

    # 添加图例
    legend_handles = []
    seen_labels = set()
    
    # 按照优先级排序图例
    sorted_styles = sorted(DS_STYLE.items(), key=lambda x: style_priority.get(x[0], 0), reverse=True)
    
    for key, style in sorted_styles:
        label = style["label"]
        if label not in seen_labels:
            # 组合图例：圆形 + 颜色 + 虚线
            # 这里简单用Patch代表颜色
            patch = mpatches.Patch(color=style["color"], label=label, alpha=0.6)
            legend_handles.append(patch)
            seen_labels.add(label)

    ax.legend(handles=legend_handles, loc='upper right',
              title="Interaction Types", fontsize=9,
              framealpha=0.9, edgecolor='gray')

    # 添加标题
    ax.text(0, 7.5, f"{ligand_resname} Interaction Diagram",
           fontsize=16, fontweight='bold', ha='center',
           bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='none', alpha=0.9))

    # 保存
    if output_path is None:
        output_path = f"{ligand_resname}_interaction_2d.png"

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"[2D Diagram] Saved to: {output_path}")
    return output_path


# PyMOL命令注册
try:
    from pymol import cmd
    cmd.extend("generate_2d_diagram", generate_2d_interaction_diagram)
except Exception:
    pass
