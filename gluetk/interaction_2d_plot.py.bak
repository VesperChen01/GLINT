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

# 氨基酸分类和颜色
AA_CATEGORIES = {
    "hydrophobic": {
        "residues": ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "PRO", "TRP"],
        "color": "#4dd0e1",
        "name_zh": "疏水",
        "name_en": "Hydrophobic"
    },
    "polar": {
        "residues": ["SER", "THR", "CYS", "TYR", "ASN", "GLN"],
        "color": "#81c784",
        "name_zh": "极性",
        "name_en": "Polar"
    },
    "positive": {
        "residues": ["ARG", "LYS", "HIS"],
        "color": "#64b5f6",
        "name_zh": "正电",
        "name_en": "Positive"
    },
    "negative": {
        "residues": ["ASP", "GLU"],
        "color": "#e57373",
        "name_zh": "负电",
        "name_en": "Negative"
    },
    "special": {
        "residues": ["GLY"],
        "color": "#bdbdbd",
        "name_zh": "特殊",
        "name_en": "Special"
    }
}

# 氨基酸侧链SMILES
AA_SMILES = {
    "ALA": "CC([NH3+])C([O-])=O",
    "ARG": "NCCCC([NH3+])C([O-])=O",
    "ASN": "NC(=O)CC([NH3+])C([O-])=O",
    "ASP": "[O-]C(=O)CC([NH3+])C([O-])=O",
    "CYS": "SCC([NH3+])C([O-])=O",
    "GLN": "NC(=O)CCC([NH3+])C([O-])=O",
    "GLU": "[O-]C(=O)CCC([NH3+])C([O-])=O",
    "GLY": "[NH3+]CC([O-])=O",
    "HIS": "c1c[nH]cn1CC([NH3+])C([O-])=O",
    "ILE": "CCC(C)C([NH3+])C([O-])=O",
    "LEU": "CC(C)CC([NH3+])C([O-])=O",
    "LYS": "[NH3+]CCCCC([NH3+])C([O-])=O",
    "MET": "CSCCC([NH3+])C([O-])=O",
    "PHE": "c1ccc(CC([NH3+])C([O-])=O)cc1",
    "PRO": "C1CC([NH2+]C1)C([O-])=O",
    "SER": "OCC([NH3+])C([O-])=O",
    "THR": "CC(O)C([NH3+])C([O-])=O",
    "TRP": "c1ccc2c(c1)c(c[nH]2)CC([NH3+])C([O-])=O",
    "TYR": "Oc1ccc(CC([NH3+])C([O-])=O)cc1",
    "VAL": "CC(C)C([NH3+])C([O-])=O"
}

def get_residue_category(resname):
    """获取残基类别和颜色"""
    res = resname[:3].upper()
    for cat, info in AA_CATEGORIES.items():
        if res in info["residues"]:
            return cat, info
    return "other", {"color": "#9e9e9e", "name_zh": "其他", "name_en": "Other"}


def generate_2d_interaction_diagram(csv_path, ligand_resname, pdb_file=None, obj_name=None,
                                     output_path=None, width=1600, height=1200, dpi=150):
    """
    生成2D相互作用图（显示残基结构式和颜色分类）
    
    参数:
        csv_path: 相互作用CSV文件路径
        ligand_resname: 配体残基名称
        pdb_file: PDB文件路径
        obj_name: PyMOL对象名称
        output_path: 输出路径
        width: 图片宽度
        height: 图片高度
        dpi: 分辨率
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
                mol = Chem.SDMolSupplier(temp_sdf, removeHs=False)[0]
                os.remove(temp_sdf)
        except Exception as e:
            print(f"[2D Diagram] Failed to extract ligand from PyMOL: {e}")

    if mol is None:
        print(f"[2D Diagram] Unable to extract ligand structure")
        return None

    if not mol.GetNumConformers():
        AllChem.Compute2DCoords(mol)

    # 创建画布
    fig = plt.figure(figsize=(width/100, height/100), dpi=dpi)
    ax = fig.add_subplot(111)
    ax.set_xlim(-8, 8)
    ax.set_ylim(-6, 6)
    ax.axis('off')

    # 绘制配体结构
    print(f"[2D Diagram] Drawing ligand structure...")
    drawer = rdMolDraw2D.MolDraw2DCairo(600, 450)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    temp_img = tempfile.mktemp(suffix=".png")
    with open(temp_img, "wb") as f:
        f.write(drawer.GetDrawingText())

    ligand_img = Image.open(temp_img)
    ax.imshow(ligand_img, extent=[-3, 3, -2.2, 2.2], zorder=10)
    os.remove(temp_img)

    # 收集残基信息
    residue_interactions = {}
    for inter in interactions:
        prot_res = inter["protein"].strip()
        res_name = prot_res.split()[0] if prot_res else "UNK"
        
        if res_name not in residue_interactions:
            residue_interactions[res_name] = []
        residue_interactions[res_name].append(inter)

    residues = list(residue_interactions.keys())
    n_residues = len(residues)
    radius = 5.5

    print(f"[2D Diagram] Drawing {n_residues} interacting residues...")

    # 绘制残基（圆形布局）
    for i, res_name in enumerate(residues):
        angle = 2 * math.pi * i / n_residues - math.pi/2
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)

        cat, cat_info = get_residue_category(res_name)
        res_color = cat_info["color"]
        cat_name = cat_info["name_en"]

        # 尝试绘制残基结构
        res_3letter = res_name[:3].upper()
        if res_3letter in AA_SMILES:
            try:
                aa_mol = Chem.MolFromSmiles(AA_SMILES[res_3letter])
                if aa_mol:
                    AllChem.Compute2DCoords(aa_mol)
                    
                    aa_drawer = rdMolDraw2D.MolDraw2DCairo(180, 180)
                    aa_drawer.DrawMolecule(aa_mol)
                    aa_drawer.FinishDrawing()
                    
                    temp_aa_img = tempfile.mktemp(suffix=".png")
                    with open(temp_aa_img, "wb") as f:
                        f.write(aa_drawer.GetDrawingText())
                    
                    aa_img = Image.open(temp_aa_img)
                    
                    rect = Rectangle((x-0.9, y-0.9), 1.8, 1.8,
                                    facecolor='white', edgecolor=res_color,
                                    linewidth=4, zorder=8, alpha=0.95)
                    ax.add_patch(rect)
                    
                    ax.imshow(aa_img, extent=[x-0.85, x+0.85, y-0.85, y+0.85], zorder=9)
                    os.remove(temp_aa_img)
                    
                    ax.text(x, y-1.15, f"{res_name}",
                           fontsize=9, fontweight='bold',
                           ha='center', va='top',
                           bbox=dict(boxstyle='round,pad=0.3', 
                                   facecolor=res_color, edgecolor='none', alpha=0.9),
                           color='white', zorder=11)
                    
                    ax.text(x, y+1.15, cat_name,
                           fontsize=7, ha='center', va='bottom',
                           bbox=dict(boxstyle='round,pad=0.2',
                                   facecolor=res_color, edgecolor='none', alpha=0.7),
                           color='white', zorder=11)
                    
            except Exception:
                circle = Circle((x, y), 0.8, facecolor=res_color,
                              edgecolor='white', linewidth=3, alpha=0.85, zorder=9)
                ax.add_patch(circle)
                ax.text(x, y, res_name[:3], fontsize=11, fontweight='bold',
                       ha='center', va='center', color='white', zorder=10)
        else:
            circle = Circle((x, y), 0.8, facecolor=res_color,
                          edgecolor='white', linewidth=3, alpha=0.85, zorder=9)
            ax.add_patch(circle)
            ax.text(x, y, res_name[:3], fontsize=11, fontweight='bold',
                   ha='center', va='center', color='white', zorder=10)

        # 绘制连线
        interactions_list = residue_interactions[res_name]
        
        interaction_colors = {
            "氢键": "#2196F3",
            "Hbond": "#2196F3",
            "盐桥": "#FF5722",
            "SaltBridge": "#FF5722",
            "疏水": "#4CAF50",
            "Hydrophobic": "#4CAF50",
            "π-π": "#9C27B0",
            "PiPi": "#9C27B0",
            "π-阳": "#E91E63",
            "PiCation": "#E91E63",
        }
        
        for inter in interactions_list:
            itype = inter["type"]
            distance = inter["distance"]
            
            line_color = "#666666"
            for key, color in interaction_colors.items():
                if key in itype:
                    line_color = color
                    break
            
            lig_x = 2.5 * math.cos(angle)
            lig_y = 1.8 * math.sin(angle)
            
            arrow = FancyArrowPatch((lig_x, lig_y), (x, y),
                                  arrowstyle='-',
                                  color=line_color,
                                  linewidth=2.5,
                                  alpha=0.6,
                                  zorder=2)
            ax.add_patch(arrow)
            
            if distance and distance != "-":
                mid_x = (lig_x + x) / 2
                mid_y = (lig_y + y) / 2
                ax.text(mid_x, mid_y, f"{distance}Å",
                       fontsize=7, ha='center', va='center',
                       bbox=dict(boxstyle='round,pad=0.2',
                               facecolor='white', edgecolor=line_color, alpha=0.9),
                       color=line_color, fontweight='bold', zorder=3)

    # 添加图例
    legend_elements = []
    for cat, info in AA_CATEGORIES.items():
        legend_elements.append(mpatches.Patch(color=info["color"], 
                                             label=info["name_en"],
                                             alpha=0.85))
    
    interaction_legend = [
        Line2D([0], [0], color="#2196F3", linewidth=3, label="Hydrogen bond"),
        Line2D([0], [0], color="#FF5722", linewidth=3, label="Salt bridge"),
        Line2D([0], [0], color="#4CAF50", linewidth=3, label="Hydrophobic"),
        Line2D([0], [0], color="#9C27B0", linewidth=3, label="Pi-Pi"),
        Line2D([0], [0], color="#E91E63", linewidth=3, label="Pi-Cation"),
    ]
    
    leg1 = ax.legend(handles=legend_elements, loc='upper left', 
                    title="Residue categories", fontsize=10, title_fontsize=11,
                    framealpha=0.95, edgecolor='black')
    ax.add_artist(leg1)
    
    ax.legend(handles=interaction_legend, loc='upper right',
                    title="Interactions", fontsize=10, title_fontsize=11,
                    framealpha=0.95, edgecolor='black')

    # 添加标题
    ax.text(0, 5.8, f"{ligand_resname} Protein-Ligand 2D Interaction Diagram",
           fontsize=18, fontweight='bold', ha='center',
           bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                   edgecolor='#333', linewidth=2, alpha=0.95))
    
    ax.text(0, -5.6, f"Total: {len(residue_interactions)} interacting residues",
           fontsize=12, ha='center', style='italic', color='#666')

    # 保存
    if output_path is None:
        output_path = f"{ligand_resname}_interaction_2d.png"

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"[2D Diagram] 2D interaction diagram saved: {output_path}")
    print(f"[2D Diagram] Residues colored by category: Hydrophobic (cyan), Polar (green), Positive (blue), Negative (red)")
    
    return output_path


# PyMOL命令注册
try:
    from pymol import cmd
    cmd.extend("generate_2d_diagram", generate_2d_interaction_diagram)
except Exception:
    pass
