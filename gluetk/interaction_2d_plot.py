# -*- coding: utf-8 -*-
"""
interaction_2d_plot.py
Generate protein-ligand 2D interaction diagrams (Schrödinger Style v4.2)

Features:
- **Sector-Based Layout**: Smart layout based on sector allocation to prevent overlaps
- **Professional Aesthetics**: Rounded badge design with professional fonts
- **High Resolution**: Default 300 DPI output
- **Visual Clarity**: Spline curves for hydrophobic interactions, dashed lines for H-bonds
- **PIL Overlay**: Uses PIL for reliable badge and line rendering

Author: GlueTK Team
"""

from __future__ import print_function
import os
import csv
import tempfile
import math
import re
import io
from typing import Dict, List, Tuple, Optional, Any

import numpy as np

# Amino acid three-letter to one-letter code mapping
AA_THREE_TO_ONE: Dict[str, str] = {
    'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'PHE': 'F', 'MET': 'M', 'PRO': 'P', 'TRP': 'W',
    'GLY': 'G', 'SER': 'S', 'THR': 'T', 'ASN': 'N', 'GLN': 'Q', 'TYR': 'Y', 'CYS': 'C',
    'ASP': 'D', 'GLU': 'E',
    'HIS': 'H', 'LYS': 'K', 'ARG': 'R'
}

# Amino acid category color scheme
AA_CATEGORY_COLORS: Dict[str, str] = {
    'nonpolar': '#FFB74D',      # Orange - Nonpolar residues
    'polar': '#64B5F6',         # Blue - Polar uncharged residues
    'negative': '#81C784',      # Green - Negatively charged residues
    'positive': '#E57373',      # Red - Positively charged residues
    'unknown': '#BDBDBD'        # Gray - Unknown residues
}

def get_aa_category(aa_code: str) -> str:
    """Get amino acid category based on single-letter code.
    
    Args:
        aa_code: Single-letter amino acid code (e.g., 'K', 'D', 'A')
        
    Returns:
        Category string: 'nonpolar', 'polar', 'negative', 'positive', or 'unknown'
    """
    nonpolar = set('AVLIFMPW')     # Nonpolar (hydrophobic)
    polar = set('GSTNQYC')         # Polar uncharged
    negative = set('DE')           # Negatively charged (acidic)
    positive = set('HKR')          # Positively charged (basic)
    
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

def parse_residue_label(res_str: str) -> Tuple[str, str, str]:
    """Parse residue label and convert to single-letter format.
    
    Args:
        res_str: Residue string like 'LYS 383' or 'K383'
        
    Returns:
        Tuple of (single_letter_code, residue_number, category)
        Example: ('K', '383', 'positive')
    """
    parts = res_str.strip().split()
    if len(parts) >= 2:
        # Three-letter format
        aa_three = parts[0].upper()
        aa_one = AA_THREE_TO_ONE.get(aa_three, aa_three[0] if aa_three else 'X')
        num = parts[1]
    elif len(parts) == 1 and any(c.isdigit() for c in parts[0]):
        # Already single-letter format 'K383'
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

# GlueTK unified color scheme (consistent with 3D view)
SCHRODINGER_STYLE: Dict[str, Dict[str, Any]] = {
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
        "show": True # Kept True to generate badge, but will skip line drawing
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

def get_interaction_style(itype: str) -> Dict[str, Any]:
    """Get interaction style configuration based on interaction type.
    
    Args:
        itype: Interaction type string (e.g., 'hydrogen bond', 'salt bridge')
        
    Returns:
        Dictionary containing color, label, style, arrow, linewidth, and show settings
    """
    itype = str(itype).lower()
    if "氢键" in itype or "hbond" in itype or "hydrogen" in itype:
        return SCHRODINGER_STYLE["Hbond"]
    elif "盐桥" in itype or "salt" in itype or "ionic" in itype:
        return SCHRODINGER_STYLE["Salt"]
    elif "cation" in itype:  # Check cation first (catch pi-cation)
        return SCHRODINGER_STYLE["Pi"]
    elif "π" in itype or "pipi" in itype or "pi-pi" in itype or "pi–pi" in itype or "stacking" in itype or "stack" in itype or "堆积" in itype:
        return SCHRODINGER_STYLE["PiPi"]
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

def generate_2d_interaction_diagram(
    csv_path: str,
    ligand_resname: str,
    pdb_file: Optional[str] = None,
    obj_name: Optional[str] = None,
    output_path: Optional[str] = None,
    width: int = 2400,
    height: int = 2000,
    dpi: int = 300,
    show_distance: bool = False,
    show_vdw: bool = False,
    protein_name: Optional[str] = None,
    compact: bool = True
) -> Optional[str]:
    """Generate high-quality 2D interaction diagram (Schrödinger Style v2.1).
    
    Args:
        csv_path: Path to CSV file containing interaction data
        ligand_resname: Ligand residue name (e.g., 'LIG', 'ATP')
        pdb_file: Optional path to PDB file containing ligand structure
        obj_name: Optional PyMOL object name to extract ligand from
        output_path: Output path for the generated image
        width: Image width in pixels (default: 2400)
        height: Image height in pixels (default: 2000)
        dpi: Image resolution (default: 300)
        show_distance: Whether to show interaction distances
        show_vdw: Whether to show van der Waals interactions
        protein_name: Optional protein name for title
        compact: Use compact layout mode (default: True)
        
    Returns:
        Path to the generated image file, or None if generation failed
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
        import matplotlib.patches as patches
        from matplotlib.patches import FancyBboxPatch
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
    elif pdb_file and os.path.exists(pdb_file):
        try:
            # Handle loading from direct PDB file
            mol_map = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)
            mol_raw = Chem.MolFromPDBFile(pdb_file, removeHs=False, sanitize=False)
            try:
                Chem.SanitizeMol(mol_map)
                Chem.SanitizeMol(mol_raw)
                mol_draw = Chem.RemoveHs(mol_raw, implicitOnly=False, updateExplicitCount=True)
            except Exception as e:
                print(f"[2D Diagram] Sanitization issue (from file): {e}")
                mol_draw = Chem.RemoveHs(mol_raw) if mol_raw else None     
        except Exception as e:
            print(f"[2D Diagram] Failed to load ligand from file: {e}")
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
    
    # 调试：打印原子映射表
    print(f"[2D Diagram] DEBUG: map_name_to_atom contains {len(map_name_to_atom)} entries")
    if len(map_name_to_atom) <= 20:
        print(f"[2D Diagram] DEBUG: Available atom names: {list(map_name_to_atom.keys())}")
    else:
        print(f"[2D Diagram] DEBUG: First 20 atom names: {list(map_name_to_atom.keys())[:20]}") 
            
    # B. 检查 3D 构象
    if mol_draw.GetNumConformers() == 0 or mol_map.GetNumConformers() == 0:
        print("[2D Diagram] ❌ Error: Missing 3D coordinates for mapping.")
        return None
        
    draw_conf = mol_draw.GetConformer()
    raw_conf = mol_map.GetConformer()

    interactions = []
    unmapped_atoms = []  # 收集映射失败的原子名
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        # 调试：打印CSV字段名
        print(f"[2D Diagram] DEBUG: CSV fieldnames = {reader.fieldnames}")
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
            
            # 容错处理：处理 "RING", "CATION" 等特殊原子名
            elif lig_atom_name in ["RING", "AROMATIC", "PI"] and mol_draw:
                print(f"[2D Diagram] ℹ️ Mapping Ring/Aromatic interaction for {prot_res}")
                # 1. Try Aromatic Atoms located in rings
                # Filter to ensure they are actually in a ring (paranoia check) and prefer Carbons
                candidates = [a.GetIdx() for a in mol_draw.GetAtoms() 
                              if a.GetIsAromatic() or (a.IsInRing() and a.GetSymbol() == 'C')]
                
                # 2. If no aromatic, just try any Ring Atom
                if not candidates:
                    candidates = [a.GetIdx() for a in mol_draw.GetAtoms() if a.IsInRing()]

                if candidates:
                    # 取中间的一个
                    target_idx = candidates[len(candidates)//2]
                else:
                     # 3. Fallback: use the first atom BUT only if no ring found
                     target_idx = 0
                    
            elif lig_atom_name in ["CATION", "QUATERNARY NITROGEN", "POS"] and mol_draw:
                print(f"[2D Diagram] ℹ️ Mapping Cation interaction for {prot_res}")
                # 寻找带正电原子 (N+)
                pos_indices = [a.GetIdx() for a in mol_draw.GetAtoms() 
                              if a.GetFormalCharge() > 0 or (a.GetSymbol() == 'N' and a.GetExplicitValence() == 4)]
                if pos_indices:
                    target_idx = pos_indices[0]
                elif "RING" not in lig_atom_name: 
                     # Fallback to aromatic
                     candidates = [a.GetIdx() for a in mol_draw.GetAtoms() if a.GetIsAromatic()]
                     if candidates:
                         target_idx = candidates[len(candidates)//2]
                     else:
                         target_idx = 0 # Ultimate fallback
            else:
                # 未映射成功
                if lig_atom_name and lig_atom_name not in unmapped_atoms:
                    unmapped_atoms.append(lig_atom_name)

            if target_idx is not None:
                interactions.append({
                    "protein": prot_res,
                    "target_idx": target_idx,
                    "type": itype,
                    "distance": dist,
                    "lig_atom_name": lig_atom_name
                })
    
    # 调试：打印未映射的原子
    if unmapped_atoms:
        print(f"[2D Diagram] ⚠️ DEBUG: {len(unmapped_atoms)} unique atom names not mapped: {unmapped_atoms[:10]}{'...' if len(unmapped_atoms) > 10 else ''}")

    # 4. 生成 2D 坐标 (Clean Molecule)
    try:
        from rdkit.Chem import rdCoordGen
        rdCoordGen.AddCoords(mol_draw)
    except ImportError:
        AllChem.Compute2DCoords(mol_draw)
            
    num_atoms = mol_draw.GetNumAtoms()
    is_large_molecule = num_atoms > 40
    
    # Dynamic Canvas Sizing
    # User feedback: "Squeezed in a clump" for large molecules. 
    # Solution: Increase canvas resolution for large molecules to give more layout space.
    if num_atoms > 100:
        drawer_w, drawer_h = (2400, 2000)
        scale_factor = 2.0
    elif num_atoms > 50:
        drawer_w, drawer_h = (1800, 1500)
        scale_factor = 1.5
    else:
        drawer_w, drawer_h = (1200, 1000)
        scale_factor = 1.0

    drawer = rdMolDraw2D.MolDraw2DCairo(drawer_w, drawer_h)
    opts = drawer.drawOptions()
    opts.clearBackground = True
    # Set white background color for the molecule drawing area
    opts.setBackgroundColour((1.0, 1.0, 1.0, 1.0))  # White background
    opts.padding = 0.15  # Slightly more padding to leave room for badges
    opts.annotationFontScale = 0.8 * scale_factor # Scale font inside molecule
    opts.bondLineWidth = 3.5 if scale_factor < 1.5 else 4.0
    opts.comicMode = False # Professional mode

    drawer.DrawMolecule(mol_draw)
    # [Critical Note] DO NOT call FinishDrawing here, otherwise overlays won't render
    
    # 获取坐标
    atom_coords = {}
    for i in range(mol_draw.GetNumAtoms()):
        pt = drawer.GetDrawCoords(i)
        atom_coords[i] = (pt.x, pt.y)

    xs = [p[0] for p in atom_coords.values()]
    ys = [p[1] for p in atom_coords.values()]
    center_x = sum(xs) / len(xs) if xs else drawer_w/2
    center_y = sum(ys) / len(ys) if ys else drawer_h/2

    print(f"[2D Diagram] Captured {len(atom_coords)} atom positions. Processing {len(interactions)} interactions.")

    # === New Layout Engine: Local Surface Normal + Force-Directed Placement ===
    
    # 1. Group interactions by Residue
    residue_data = {}
    for inter in interactions:
        p = inter["protein"]
        if p not in residue_data: residue_data[p] = []
        residue_data[p].append(inter)
        
    # 2. Prepare Items for Layout
    layout_items = []
    
    # Helper: calculate outward vector for an atom
    def get_outward_vector(atom_idx: int, all_coords: Dict[int, Tuple[float, float]],
                           center_p: Tuple[float, float]) -> Tuple[float, float]:
        ax_local, ay_local = all_coords[atom_idx]
        mx, my = 0.0, 0.0
        
        has_neighbors = False
        try:
            atom = mol_draw.GetAtomWithIdx(atom_idx)
            neighbors = atom.GetNeighbors()
            
            vecs_x, vecs_y = 0.0, 0.0
            count = 0
            for nb in neighbors:
                nb_idx = nb.GetIdx()
                if nb_idx in all_coords:
                    nx, ny = all_coords[nb_idx]
                    dx, dy = ax_local - nx, ay_local - ny
                    d = math.sqrt(dx*dx + dy*dy)
                    if d > 1e-4:
                        vecs_x += dx/d
                        vecs_y += dy/d
                        count += 1
            if count > 0:
                mx, my = vecs_x, vecs_y
                has_neighbors = True
        except:
            pass
            
        cx, cy = center_p
        gx, gy = ax_local - cx, ay_local - cy
        g_len = math.sqrt(gx*gx + gy*gy)
        if g_len > 1e-4:
            gx /= g_len
            gy /= g_len
            if has_neighbors:
                mx = mx * 0.7 + gx * 0.3
                my = my * 0.7 + gy * 0.3
            else:
                mx, my = gx, gy
        
        m_len = math.sqrt(mx*mx + my*my)
        if m_len < 1e-6: return (1.0, 0.0) 
        return (mx/m_len, my/m_len)

    badge_radius = 45.0 * (1.0 + (scale_factor - 1.0) * 0.5)
    # Increase base distance significantly to prevent overlap with molecule
    # For large molecules, we need even more distance
    base_distance = 200.0 * scale_factor
    
    # Calculate molecule bounding box for collision avoidance
    mol_min_x = min(xs) - 30 * scale_factor
    mol_max_x = max(xs) + 30 * scale_factor
    mol_min_y = min(ys) - 30 * scale_factor
    mol_max_y = max(ys) + 30 * scale_factor
    
    for prot_res, inter_list in residue_data.items():
        # Find the best anchor - use the first non-hydrophobic interaction if available
        # Otherwise use the first interaction
        anchor_idx = None
        primary_inter = None
        for inter in inter_list:
            style = get_interaction_style(inter["type"])
            if style["label"] != "Hydrophobic":
                anchor_idx = inter["target_idx"]
                primary_inter = inter
                break
        if anchor_idx is None:
            anchor_idx = inter_list[0]["target_idx"]
            primary_inter = inter_list[0]
            
        if anchor_idx not in atom_coords: continue
        
        cx, cy = atom_coords[anchor_idx]
        vx, vy = get_outward_vector(anchor_idx, atom_coords, (center_x, center_y))
        
        layout_items.append({
            "id": prot_res,
            "x": cx + vx * base_distance,
            "y": cy + vy * base_distance,
            "r": badge_radius,
            "vx": vx, "vy": vy,
            "anchor": (cx, cy),
            "anchor_idx": anchor_idx,
            "inters": inter_list,
            "primary_inter": primary_inter,
            "res_data": parse_residue_label(prot_res)
        })

    # 调试：打印layout_items数量
    print(f"[2D Diagram] DEBUG: residue_data contains {len(residue_data)} residues")
    print(f"[2D Diagram] DEBUG: layout_items contains {len(layout_items)} items")
    print(f"[2D Diagram] DEBUG: Canvas size = {drawer_w} x {drawer_h}")
    if len(layout_items) == 0 and len(residue_data) > 0:
        # 检查为什么没有生成layout_items
        for prot_res, inter_list in list(residue_data.items())[:3]:
            anchor_idx = inter_list[0]["target_idx"]
            print(f"[2D Diagram] DEBUG: Residue {prot_res}: anchor_idx={anchor_idx}, in atom_coords={anchor_idx in atom_coords}")
    elif len(layout_items) > 0:
        # 打印前3个气泡的坐标
        for item in layout_items[:3]:
            print(f"[2D Diagram] DEBUG: Badge '{item['id']}' pos=({item['x']:.1f}, {item['y']:.1f}), anchor=({item['anchor'][0]:.1f}, {item['anchor'][1]:.1f})")

    # Collision Resolution - Badge vs Badge AND Badge vs Molecule
    for iteration in range(120):
        for i in range(len(layout_items)):
            item = layout_items[i]
            
            # A. Badge vs Badge collision
            for j in range(i + 1, len(layout_items)):
                it_i, it_j = layout_items[i], layout_items[j]
                dx, dy = it_i["x"] - it_j["x"], it_i["y"] - it_j["y"]
                d = math.sqrt(dx*dx + dy*dy)
                min_d = it_i["r"] + it_j["r"] + 25 # Increased spacing
                if d < min_d:
                    if d < 1e-4: dx, dy, d = 1.0, 0.0, 1.0
                    push = (min_d - d) * 0.55
                    it_i["x"] += (dx/d)*push; it_i["y"] += (dy/d)*push
                    it_j["x"] -= (dx/d)*push; it_j["y"] -= (dy/d)*push
            
            # B. Badge vs Molecule collision - push badge away from molecule center
            bx, by, br = item["x"], item["y"], item["r"]
            # Check if badge overlaps with molecule bounding box (with margin)
            margin_mol = br + 15 * scale_factor
            if (bx + br > mol_min_x - margin_mol and bx - br < mol_max_x + margin_mol and
                by + br > mol_min_y - margin_mol and by - br < mol_max_y + margin_mol):
                # Badge is too close to molecule - push it outward
                ax, ay = item["anchor"]
                dx, dy = bx - ax, by - ay
                d = math.sqrt(dx*dx + dy*dy)
                if d < 1e-4:
                    dx, dy, d = item["vx"], item["vy"], 1.0
                # Push further out along the anchor->badge direction
                push_dist = 15 * scale_factor
                item["x"] += (dx/d) * push_dist
                item["y"] += (dy/d) * push_dist

    # Boundary Clamping (Prevent clipping)
    margin = 25 * scale_factor
    for item in layout_items:
        r = item["r"]
        # Clamp X
        if item["x"] < r + margin: item["x"] = r + margin
        if item["x"] > drawer_w - r - margin: item["x"] = drawer_w - r - margin
        # Clamp Y
        if item["y"] < r + margin: item["y"] = r + margin
        if item["y"] > drawer_h - r - margin: item["y"] = drawer_h - r - margin

    # 5. Finish RDKit drawing and get base image
    drawer.FinishDrawing()
    
    # Convert RDKit drawing to PIL Image for overlay
    png_data = drawer.GetDrawingText()
    base_img = Image.open(io.BytesIO(png_data)).convert("RGBA")
    
    # Replace gray/transparent background with white
    # Create a white background image
    white_bg = Image.new("RGBA", base_img.size, (255, 255, 255, 255))
    # Composite the molecule image onto white background
    base_img = Image.alpha_composite(white_bg, base_img)
    
    # Create overlay layer for badges and lines
    overlay = Image.new("RGBA", base_img.size, (255, 255, 255, 0))
    
    # Import PIL drawing tools
    from PIL import ImageDraw, ImageFont
    draw = ImageDraw.Draw(overlay)
    
    def hex_to_rgb(h):
        h = h.lstrip('#')
        return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
    
    # Try to load a nice font, fallback to default
    try:
        # Try common system fonts
        font_size = int(badge_radius * 0.55)
        font_paths = [
            "/System/Library/Fonts/Helvetica.ttc",
            "/System/Library/Fonts/SFNSText.ttf",
            "/Library/Fonts/Arial.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
            "C:/Windows/Fonts/arial.ttf"
        ]
        font = None
        for fp in font_paths:
            if os.path.exists(fp):
                try:
                    font = ImageFont.truetype(fp, font_size)
                    break
                except:
                    continue
        if font is None:
            font = ImageFont.load_default()
    except:
        font = ImageFont.load_default()
    
    # A. Draw interaction lines
    # Strategy: Draw ONE line per residue badge to its anchor atom
    # For hydrophobic: draw ellipses on atoms but NO line to badge
    # For other types: draw ONE line from badge to anchor
    lines_drawn = 0
    hydro_drawn = 0
    hydro_atoms_drawn = set()  # Track which atoms already have hydrophobic ellipses
    
    for item in layout_items:
        bx, by, br = item["x"], item["y"], item["r"]
        anchor_idx = item["anchor_idx"]
        
        # Collect interaction types for this residue
        has_non_hydrophobic = False
        interaction_types = set()
        
        for inter in item["inters"]:
            itype = inter["type"]
            style = get_interaction_style(itype)
            interaction_types.add(style["label"])
            
            if style["label"] != "Hydrophobic":
                has_non_hydrophobic = True
            
            # Draw hydrophobic ellipses on atoms (deduplicated)
            if style["label"] == "Hydrophobic":
                target_idx = inter["target_idx"]
                if target_idx not in atom_coords:
                    continue
                if target_idx in hydro_atoms_drawn:
                    continue  # Skip duplicate ellipses
                    
                tx, ty = atom_coords[target_idx]
                radius_h = 16 * scale_factor
                color_rgb = hex_to_rgb(style["color"])
                # Semi-transparent fill
                draw.ellipse(
                    [tx - radius_h, ty - radius_h, tx + radius_h, ty + radius_h],
                    fill=color_rgb + (50,),  # Alpha = 50 (more transparent)
                    outline=color_rgb + (120,),
                    width=2
                )
                hydro_atoms_drawn.add(target_idx)
                hydro_drawn += 1
        
        # Draw ONE line from badge to anchor for non-hydrophobic interactions
        if has_non_hydrophobic and anchor_idx in atom_coords:
            tx, ty = atom_coords[anchor_idx]
            
            # Determine line style based on primary interaction type
            primary_inter = item.get("primary_inter", item["inters"][0])
            style = get_interaction_style(primary_inter["type"])
            
            # Calculate line start point (from badge edge toward target atom)
            dx, dy = tx - bx, ty - by
            dist = math.sqrt(dx*dx + dy*dy)
            if dist > 1e-5:
                # Start from badge rim
                rim_x = bx + (dx/dist) * br
                rim_y = by + (dy/dist) * br
                
                # Line color and width
                color_rgb = hex_to_rgb(style["color"])
                line_width = int(2.5 * scale_factor)
                
                # Draw line based on style
                if style["style"] == "--":
                    # Dashed line
                    _draw_dashed_line(draw, rim_x, rim_y, tx, ty, color_rgb + (220,), line_width, dash_length=12)
                elif style["style"] == ":":
                    # Dotted line
                    _draw_dashed_line(draw, rim_x, rim_y, tx, ty, color_rgb + (220,), line_width, dash_length=6)
                else:
                    # Solid line
                    draw.line([(rim_x, rim_y), (tx, ty)], fill=color_rgb + (220,), width=line_width)
                
                lines_drawn += 1

    print(f"[2D Diagram] DEBUG: Drew {lines_drawn} interaction lines, {hydro_drawn} hydrophobic ellipses")

    # B. Draw residue badges
    print(f"[2D Diagram] DEBUG: Drawing {len(layout_items)} badges...")
    for item in layout_items:
        bx, by, br = item["x"], item["y"], item["r"]
        aa_letter, aa_num, aa_cat = item["res_data"]
        badge_color = AA_CATEGORY_COLORS.get(aa_cat, '#9E9E9E')
        color_rgb = hex_to_rgb(badge_color)
        
        # Badge background (white fill)
        draw.ellipse(
            [bx - br, by - br, bx + br, by + br],
            fill=(255, 255, 255, 240),
            outline=color_rgb + (255,),
            width=3
        )
        
        # Badge text
        label = f"{aa_letter}{aa_num}"
        # Get text bounding box for centering
        try:
            bbox = draw.textbbox((0, 0), label, font=font)
            tw = bbox[2] - bbox[0]
            th = bbox[3] - bbox[1]
        except:
            # Fallback for older PIL versions
            tw, th = draw.textsize(label, font=font) if hasattr(draw, 'textsize') else (len(label)*10, 20)
        
        text_x = bx - tw / 2
        text_y = by - th / 2
        draw.text((text_x, text_y), label, fill=(40, 40, 40, 255), font=font)

    # C. Draw Legend
    _draw_legend(draw, drawer_w, drawer_h, scale_factor, font)

    # Composite overlay onto base image
    result = Image.alpha_composite(base_img, overlay)
    
    # Convert to RGB for saving (PNG with transparency or JPEG)
    if output_path is None:
        output_path = f"{ligand_resname}_2d.png"
    
    # Save as PNG
    result.save(output_path, "PNG", dpi=(dpi, dpi))
        
    print(f"[2D Diagram] ✅ v4.2 PIL-Overlay Diagram successfully saved to {output_path}")
    return output_path


def _draw_dashed_line(draw, x1, y1, x2, y2, color, width, dash_length=10):
    """Draw a dashed line using PIL.
    
    Args:
        draw: PIL ImageDraw object
        x1, y1: Start coordinates
        x2, y2: End coordinates
        color: Line color (RGBA tuple)
        width: Line width
        dash_length: Length of each dash segment
    """
    dx = x2 - x1
    dy = y2 - y1
    length = math.sqrt(dx*dx + dy*dy)
    if length < 1:
        return
    
    dx /= length
    dy /= length
    
    pos = 0
    drawing = True
    while pos < length:
        seg_len = min(dash_length, length - pos)
        if drawing:
            sx = x1 + dx * pos
            sy = y1 + dy * pos
            ex = x1 + dx * (pos + seg_len)
            ey = y1 + dy * (pos + seg_len)
            draw.line([(sx, sy), (ex, ey)], fill=color, width=width)
        pos += dash_length
        drawing = not drawing

def _draw_legend(draw, width, height, scale_factor, font):
    """Draw legend for interaction types."""
    # Legend settings
    box_width = 300 * scale_factor
    line_height = 40 * scale_factor
    padding = 20 * scale_factor
    
    # Filter active styles
    active_styles = [s for k, s in SCHRODINGER_STYLE.items() if s["show"] and k != "VDW"]
    
    box_height = len(active_styles) * line_height + padding * 2
    
    # Position: Bottom Right
    start_x = width - box_width - padding
    start_y = height - box_height - padding
    
    # Draw Background
    draw.rectangle(
        [start_x, start_y, start_x + box_width, start_y + box_height],
        fill=(255, 255, 255, 200),
        outline=(200, 200, 200, 255),
        width=1
    )
    
    cursor_y = start_y + padding
    
    for style in active_styles:
        # Draw Line/Icon
        icon_x = start_x + padding
        icon_y = cursor_y + line_height / 2
        icon_w = 40 * scale_factor
        
        color_rgb = tuple(int(style["color"].lstrip('#')[i:i+2], 16) for i in (0, 2, 4))
        
        if style["label"] == "Hydrophobic":
            # Draw green sphere/ellipse
            r = 8 * scale_factor
            draw.ellipse(
                [icon_x, icon_y - r, icon_x + 2*r, icon_y + r],
                fill=color_rgb + (100,),
                outline=color_rgb + (255,),
                width=2
            )
        else:
            # Draw Line
            line_w = int(3 * scale_factor)
            if style["style"] == "--":
                 _draw_dashed_line(draw, icon_x, icon_y, icon_x + icon_w, icon_y, color_rgb + (255,), line_w, dash_length=8)
            elif style["style"] == ":":
                 _draw_dashed_line(draw, icon_x, icon_y, icon_x + icon_w, icon_y, color_rgb + (255,), line_w, dash_length=4)
            else:
                 draw.line([(icon_x, icon_y), (icon_x + icon_w, icon_y)], fill=color_rgb + (255,), width=line_w)
                 
        # Draw Label
        text_x = icon_x + icon_w + padding
        text_y = cursor_y + (line_height - 20) / 2 # Approx vert center
        draw.text((text_x, text_y), style["label"], fill=(50, 50, 50, 255), font=font)
        
        cursor_y += line_height
        
def _normalize_interaction_type(itype: str) -> str:
    """Normalize interaction type string to lowercase.
    
    Args:
        itype: Raw interaction type string
        
    Returns:
        Lowercase normalized string
    """
    return str(itype).lower()


# Register PyMOL command if available
try:
    from pymol import cmd
    cmd.extend("generate_2d_diagram", generate_2d_interaction_diagram)
except ImportError:
    pass
