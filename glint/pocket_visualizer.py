# -*- coding: utf-8 -*-
"""
GLINT Pocket Visualizer
=========================
口袋可视化模块

功能：
- PyMOL CGO 对象生成
- 按性质着色（体积、疏水性、可成药性、静电势）
- 球体、网格、表面多种显示模式
- 与相互作用网络联动
- 与 APBS 静电势叠加
"""

import numpy as np
import os

try:
    from pymol import cmd
    from pymol.cgo import COLOR, SPHERE, CYLINDER, BEGIN, END, TRIANGLES, NORMAL, VERTEX
except ImportError:
    cmd = None
    COLOR = SPHERE = CYLINDER = BEGIN = END = TRIANGLES = NORMAL = VERTEX = None


def _zh():
    """检测中文环境"""
    import locale
    try:
        lang = locale.getdefaultlocale()[0]
        return lang and lang.startswith('zh')
    except:
        return False


def _info(cn, en):
    """双语信息输出"""
    print(cn if _zh() else en)


def _get_color_gradient(value, vmin, vmax, colormap='blue_white_red'):
    """
    根据数值返回颜色
    
    colormap:
        'blue_white_red': 蓝-白-红（负-中-正）
        'white_orange': 白-橙（疏水性）
        'gray_green': 灰-绿（可成药性）
        'rainbow': 彩虹色谱
    """
    # 归一化到 [0, 1]
    if vmax - vmin == 0:
        t = 0.5
    else:
        t = (value - vmin) / (vmax - vmin)
    t = max(0.0, min(1.0, t))
    
    if colormap == 'blue_white_red':
        if t < 0.5:
            # 蓝 -> 白
            s = t * 2
            return (s, s, 1.0)
        else:
            # 白 -> 红
            s = (t - 0.5) * 2
            return (1.0, 1.0 - s, 1.0 - s)
    
    elif colormap == 'white_orange':
        # 白 -> 橙
        return (1.0, 1.0 - t * 0.5, 1.0 - t)
    
    elif colormap == 'gray_green':
        # 灰 -> 绿
        return (0.5 - t * 0.3, 0.5 + t * 0.5, 0.5 - t * 0.3)
    
    elif colormap == 'rainbow':
        # 彩虹：蓝-青-绿-黄-红
        if t < 0.25:
            s = t * 4
            return (0.0, s, 1.0)
        elif t < 0.5:
            s = (t - 0.25) * 4
            return (0.0, 1.0, 1.0 - s)
        elif t < 0.75:
            s = (t - 0.5) * 4
            return (s, 1.0, 0.0)
        else:
            s = (t - 0.75) * 4
            return (1.0, 1.0 - s, 0.0)
    
    else:
        return (0.8, 0.8, 0.8)


def visualize_pockets(pockets, obj_name='pockets', color_by='volume',
                     show_spheres=True, show_surface=False, show_mesh=False,
                     sphere_radius=1.0, transparency=0.5):
    """
    可视化口袋
    
    参数：
        pockets: detect_pockets 返回的口袋列表
        obj_name: PyMOL 对象名前缀
        color_by: 着色依据 ('volume', 'hydrophobicity', 'druggability', 'charge', 'depth')
        show_spheres: 显示球体（口袋中心）
        show_surface: 显示表面
        show_mesh: 显示网格
        sphere_radius: 球体半径倍数
        transparency: 透明度
    """
    if not cmd:
        _info("错误：需要 PyMOL 环境", "Error: PyMOL required")
        return
    
    if len(pockets) == 0:
        _info("没有口袋可显示", "No pockets to visualize")
        return
    
    # 确定颜色范围
    if color_by == 'volume':
        values = [p['volume'] for p in pockets]
        colormap = 'rainbow'
    elif color_by == 'hydrophobicity':
        values = [p['hydrophobicity'] for p in pockets]
        colormap = 'white_orange'
    elif color_by == 'druggability':
        values = [p['druggability_score'] for p in pockets]
        colormap = 'gray_green'
    elif color_by == 'charge':
        values = [p['net_charge'] for p in pockets]
        colormap = 'blue_white_red'
    elif color_by == 'depth':
        values = [p['depth'] for p in pockets]
        colormap = 'rainbow'
    else:
        values = [0.5] * len(pockets)
        colormap = 'gray_green'
    
    vmin, vmax = min(values), max(values)
    
    _info(f"可视化 {len(pockets)} 个口袋（按 {color_by} 着色）",
          f"Visualizing {len(pockets)} pockets (colored by {color_by})")
    
    # 显示球体（口袋中心）
    if show_spheres:
        for i, pocket in enumerate(pockets):
            color = _get_color_gradient(values[i], vmin, vmax, colormap)
            center = pocket['center']
            radius = sphere_radius * (pocket['volume'] / 100.0) ** (1/3)
            
            # 创建 CGO 球体
            cgo_obj = [
                COLOR, color[0], color[1], color[2],
                SPHERE, center[0], center[1], center[2], radius
            ]
            
            pocket_obj_name = f"{obj_name}_sphere_{pocket['id']}"
            cmd.load_cgo(cgo_obj, pocket_obj_name)
            cmd.set('cgo_transparency', transparency, pocket_obj_name)
    
    # 显示网格点
    if show_mesh:
        for i, pocket in enumerate(pockets):
            color = _get_color_gradient(values[i], vmin, vmax, colormap)
            grid_coords = np.array(pocket['grid_coords'])
            
            # 采样（避免太多点）
            if len(grid_coords) > 1000:
                indices = np.random.choice(len(grid_coords), 1000, replace=False)
                grid_coords = grid_coords[indices]
            
            # 转换为实际坐标（从 pocket_detector 获取 origin 和 spacing）
            # 这里简化：直接使用 center 附近的点
            cgo_obj = [COLOR, color[0], color[1], color[2]]
            
            for coord in grid_coords[:200]:  # 限制数量
                # 这里需要从 detector 传入 origin 和 spacing，暂时跳过精确实现
                pass
            
            if len(cgo_obj) > 3:
                pocket_obj_name = f"{obj_name}_mesh_{pocket['id']}"
                cmd.load_cgo(cgo_obj, pocket_obj_name)
    
    # 高亮口袋残基
    for pocket in pockets:
        residues = pocket['residues']
        if len(residues) > 0:
            selection = ' or '.join([
                f"(chain {r['chain']} and resi {r['resi']})"
                for r in residues
            ])
            pocket_sel_name = f"{obj_name}_residues_{pocket['id']}"
            cmd.select(pocket_sel_name, selection)
            
            # 显示为 sticks
            cmd.show('sticks', pocket_sel_name)
            
            # 按口袋着色
            color = _get_color_gradient(values[pockets.index(pocket)], vmin, vmax, colormap)
            cmd.set_color(f'pocket_color_{pocket["id"]}', list(color))
            cmd.color(f'pocket_color_{pocket["id"]}', pocket_sel_name)
    
    _info(f"可视化完成：{obj_name}_*", f"Visualization complete: {obj_name}_*")


def show_pocket_labels(pockets, obj_name='pocket_labels'):
    """显示口袋标签"""
    if not cmd:
        return
    
    for pocket in pockets:
        center = pocket['center']
        label_text = f"P{pocket['id']}\\nV:{pocket['volume']:.0f}\\nD:{pocket['druggability_score']:.2f}"
        
        # 创建伪原子用于标签
        pseudo_name = f"{obj_name}_{pocket['id']}"
        cmd.pseudoatom(pseudo_name, pos=center)
        cmd.label(pseudo_name, f'"{label_text}"')


def visualize_pocket_comparison(comparison, pockets_a, pockets_b, 
                                obj_a_name='obj_a', obj_b_name='obj_b'):
    """
    可视化口袋对比结果
    
    参数：
        comparison: compare_pockets 返回的对比结果
        pockets_a, pockets_b: 两组口袋
        obj_a_name, obj_b_name: 对象名
    """
    if not cmd:
        return
    
    _info("可视化口袋对比...", "Visualizing pocket comparison...")
    
    # 颜色方案
    COLOR_MATCHED = [0.5, 0.5, 1.0]  # 蓝色 - 匹配
    COLOR_EXPANDED = [0.0, 1.0, 0.0]  # 绿色 - 扩大
    COLOR_SHRANK = [1.0, 1.0, 0.0]    # 黄色 - 缩小
    COLOR_NEW = [1.0, 0.5, 0.0]       # 橙色 - 新增
    COLOR_LOST = [1.0, 0.0, 0.0]      # 红色 - 消失
    
    for comp in comparison:
        if comp['match_type'] == 'matched':
            pocket_a = next(p for p in pockets_a if p['id'] == comp['pocket_a_id'])
            pocket_b = next(p for p in pockets_b if p['id'] == comp['pocket_b_id'])
            
            # 根据体积变化选择颜色
            if comp['delta_volume'] > 10:
                color = COLOR_EXPANDED
            elif comp['delta_volume'] < -10:
                color = COLOR_SHRANK
            else:
                color = COLOR_MATCHED
            
            # 绘制连线
            cgo_obj = [
                COLOR, color[0], color[1], color[2],
                CYLINDER,
                pocket_a['center'][0], pocket_a['center'][1], pocket_a['center'][2],
                pocket_b['center'][0], pocket_b['center'][1], pocket_b['center'][2],
                0.2, color[0], color[1], color[2], color[0], color[1], color[2]
            ]
            cmd.load_cgo(cgo_obj, f"pocket_link_{comp['pocket_a_id']}_{comp['pocket_b_id']}")
        
        elif comp['match_type'] == 'new':
            pocket_b = next(p for p in pockets_b if p['id'] == comp['pocket_b_id'])
            cgo_obj = [
                COLOR, COLOR_NEW[0], COLOR_NEW[1], COLOR_NEW[2],
                SPHERE, pocket_b['center'][0], pocket_b['center'][1], pocket_b['center'][2], 2.0
            ]
            cmd.load_cgo(cgo_obj, f"pocket_new_{comp['pocket_b_id']}")
        
        elif comp['match_type'] == 'lost':
            pocket_a = next(p for p in pockets_a if p['id'] == comp['pocket_a_id'])
            cgo_obj = [
                COLOR, COLOR_LOST[0], COLOR_LOST[1], COLOR_LOST[2],
                SPHERE, pocket_a['center'][0], pocket_a['center'][1], pocket_a['center'][2], 2.0
            ]
            cmd.load_cgo(cgo_obj, f"pocket_lost_{comp['pocket_a_id']}")
    
    _info("对比可视化完成", "Comparison visualization complete")


def overlay_pocket_electrostatics(pocket, apbs_map_file=None, obj_name='protein'):
    """
    将口袋与 APBS 静电势叠加显示
    
    参数：
        pocket: 单个口袋字典
        apbs_map_file: APBS 输出的 .dx 文件
        obj_name: 蛋白对象名
    """
    if not cmd:
        return
    
    if apbs_map_file and os.path.exists(apbs_map_file):
        # 加载静电势图
        map_name = f"pocket_{pocket['id']}_electro"
        cmd.load(apbs_map_file, map_name)
        
        # 在口袋位置显示等势面
        center = pocket['center']
        cmd.isomesh(f"{map_name}_pos", map_name, 1.0, 
                   f"({center[0]},{center[1]},{center[2]})", carve=5.0)
        cmd.isomesh(f"{map_name}_neg", map_name, -1.0, 
                   f"({center[0]},{center[1]},{center[2]})", carve=5.0)
        
        cmd.color('red', f"{map_name}_pos")
        cmd.color('blue', f"{map_name}_neg")
        
        _info(f"已叠加口袋 {pocket['id']} 的静电势", 
              f"Overlayed electrostatics for pocket {pocket['id']}")
    else:
        _info("未找到 APBS 地图文件", "APBS map file not found")


def highlight_pocket_interactions(pocket, interaction_csv=None):
    """
    高亮口袋内的相互作用
    
    参数：
        pocket: 单个口袋字典
        interaction_csv: 相互作用 CSV 文件
    """
    if not cmd:
        return
    
    # 获取口袋残基
    pocket_residues = set((r['chain'], r['resi']) for r in pocket['residues'])
    
    # 读取相互作用 CSV
    if interaction_csv and os.path.exists(interaction_csv):
        import csv
        interactions_in_pocket = []
        
        with open(interaction_csv, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # 检查是否在口袋内
                res_key = None
                if 'Protein_Chain' in row and 'Protein_Residue' in row:
                    chain = row['Protein_Chain']
                    resi = row['Protein_Residue'].split()[-1]  # 提取残基号
                    res_key = (chain, resi)
                elif 'Chain1' in row and 'Residue1' in row:
                    chain = row['Chain1']
                    resi = row['Residue1'].split()[-1]
                    res_key = (chain, resi)
                
                if res_key in pocket_residues:
                    interactions_in_pocket.append(row)
        
        if len(interactions_in_pocket) > 0:
            _info(f"口袋 {pocket['id']} 包含 {len(interactions_in_pocket)} 个相互作用",
                  f"Pocket {pocket['id']} contains {len(interactions_in_pocket)} interactions")
            
            # TODO: 可视化这些相互作用（虚线、标签等）
    
    return pocket_residues


def create_pocket_surface(pocket, obj_name='pocket_surface', grid_data=None):
    """
    创建口袋表面（marching cubes）
    
    参数：
        pocket: 单个口袋字典
        obj_name: 对象名
        grid_data: (grid, origin) 元组（从 detector 传入）
    """
    # 需要实现 marching cubes 算法或使用 PyMOL 的内置功能
    # 这里简化：使用球体近似
    if not cmd:
        return
    
    center = pocket['center']
    radius = (3 * pocket['volume'] / (4 * np.pi)) ** (1/3)
    
    # 创建球形表面
    cmd.pseudoatom(obj_name, pos=center, vdw=radius)
    cmd.show('surface', obj_name)
    cmd.set('transparency', 0.5, obj_name)


def visualize_pockets_with_interactions(pockets, interaction_csv, 
                                        obj_name='pocket_int'):
    """
    联合显示口袋与相互作用
    
    参数：
        pockets: 口袋列表
        interaction_csv: 相互作用 CSV
        obj_name: 对象名前缀
    """
    if not cmd:
        return
    
    # 可视化口袋
    visualize_pockets(pockets, obj_name=obj_name, color_by='druggability')
    
    # 为每个口袋高亮相互作用
    for pocket in pockets:
        highlight_pocket_interactions(pocket, interaction_csv)


def export_pocket_to_pdb(pocket, output_pdb, grid_data=None):
    """
    将口袋导出为 PDB 伪原子（方便其他软件分析）
    
    参数：
        pocket: 单个口袋字典
        output_pdb: 输出 PDB 文件
        grid_data: (grid, origin) 元组
    """
    with open(output_pdb, 'w') as f:
        f.write(f"REMARK Pocket {pocket['id']}\n")
        f.write(f"REMARK Volume: {pocket['volume']:.2f} A^3\n")
        f.write(f"REMARK Druggability: {pocket['druggability_score']:.3f}\n")
        
        # 中心点
        center = pocket['center']
        f.write(f"HETATM    1  CTR POC A   1    "
                f"{center[0]:8.3f}{center[1]:8.3f}{center[2]:8.3f}"
                f"  1.00  0.00           C\n")
        
        # 网格点（采样）
        if 'grid_coords' in pocket and len(pocket['grid_coords']) > 0:
            grid_coords = np.array(pocket['grid_coords'])
            # 采样 100 个点
            if len(grid_coords) > 100:
                indices = np.random.choice(len(grid_coords), 100, replace=False)
                grid_coords = grid_coords[indices]
            
            for i, coord in enumerate(grid_coords, start=2):
                # 需要转换为实际坐标
                # 这里简化：使用相对坐标
                f.write(f"HETATM{i:5d}  DOT POC A   1    "
                        f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                        f"  1.00  0.00           H\n")
        
        f.write("END\n")
    
    _info(f"口袋已导出到 {output_pdb}", f"Pocket exported to {output_pdb}")


if __name__ == '__main__':
    print("GLINT Pocket Visualizer - 请在 PyMOL 中使用")
    print("示例：visualize_pockets(pockets, color_by='druggability')")
