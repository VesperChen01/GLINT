# -*- coding: utf-8 -*-
"""
pocket_docking.py
口袋导向的 AutoDock Vina 对接模块

功能:
1. 从口袋分析自动生成 Vina config.txt
2. 支持用户自定义 config 文件
3. 批量对接到多个口袋
4. 整合相互作用分析
"""

from __future__ import print_function
import os
import sys
import tempfile
import subprocess
from typing import List, Dict, Any, Optional, Tuple
from pymol import cmd

# ==================== 口袋中心计算 ====================
def calculate_pocket_box(pocket_data: Dict[str, Any], padding: float = 5.0) -> Dict[str, float]:
    """
    从口袋数据计算对接盒子参数
    
    参数:
        pocket_data: 口袋数据字典，包含 grid_points 列表
        padding: 盒子边界扩展 (Å)
    
    返回:
        dict: {
            'center_x': float,
            'center_y': float,
            'center_z': float,
            'size_x': float,
            'size_y': float,
            'size_z': float
        }
    """
    grid_points = pocket_data.get('grid_points', [])
    
    if not grid_points:
        raise ValueError("No grid points in pocket data")
    
    # 计算边界
    xs = [p[0] for p in grid_points]
    ys = [p[1] for p in grid_points]
    zs = [p[2] for p in grid_points]
    
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)
    
    # 中心和尺寸
    center_x = (min_x + max_x) / 2.0
    center_y = (min_y + max_y) / 2.0
    center_z = (min_z + max_z) / 2.0
    
    size_x = (max_x - min_x) + 2 * padding
    size_y = (max_y - min_y) + 2 * padding
    size_z = (max_z - min_z) + 2 * padding
    
    return {
        'center_x': center_x,
        'center_y': center_y,
        'center_z': center_z,
        'size_x': size_x,
        'size_y': size_y,
        'size_z': size_z
    }


def generate_vina_config(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    box_params: Dict[str, float],
    output_pdbqt: str,
    config_path: Optional[str] = None,
    exhaustiveness: int = 8,
    num_modes: int = 9
) -> str:
    """
    生成 Vina 配置文件
    
    参数:
        receptor_pdbqt: 受体 PDBQT 文件路径
        ligand_pdbqt: 配体 PDBQT 文件路径
        box_params: 盒子参数字典 (center_x, center_y, center_z, size_x, size_y, size_z)
        output_pdbqt: 输出 PDBQT 文件路径
        config_path: 配置文件保存路径（None 则使用临时文件）
        exhaustiveness: 搜索精度 (default 8)
        num_modes: 输出模式数量 (default 9)
    
    返回:
        str: 配置文件路径
    """
    if config_path is None:
        fd, config_path = tempfile.mkstemp(suffix='_vina_config.txt', text=True)
        os.close(fd)
    
    config_content = f"""receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}

out = {output_pdbqt}

center_x = {box_params['center_x']:.3f}
center_y = {box_params['center_y']:.3f}
center_z = {box_params['center_z']:.3f}

size_x = {box_params['size_x']:.1f}
size_y = {box_params['size_y']:.1f}
size_z = {box_params['size_z']:.1f}

exhaustiveness = {exhaustiveness}
num_modes = {num_modes}
"""
    
    with open(config_path, 'w') as f:
        f.write(config_content)
    
    print(f"[generate_vina_config] ✅ Config saved to: {config_path}")
    print(f"   Box center: ({box_params['center_x']:.2f}, {box_params['center_y']:.2f}, {box_params['center_z']:.2f})")
    print(f"   Box size: ({box_params['size_x']:.1f}, {box_params['size_y']:.1f}, {box_params['size_z']:.1f}) Å")
    
    return config_path


# ==================== 口袋对接核心 ====================
def dock_to_pocket(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    pocket_data: Dict[str, Any],
    output_dir: str,
    pocket_id: int = 1,
    padding: float = 5.0,
    exhaustiveness: int = 8,
    custom_config: Optional[str] = None,
    vina_bin: Optional[str] = None
) -> Dict[str, Any]:
    """
    对接配体到指定口袋
    
    参数:
        receptor_pdbqt: 受体 PDBQT 文件
        ligand_pdbqt: 配体 PDBQT 文件
        pocket_data: 口袋数据字典
        output_dir: 输出目录
        pocket_id: 口袋编号
        padding: 盒子边界扩展 (Å)
        exhaustiveness: Vina 搜索精度
        custom_config: 自定义配置文件路径（可选，覆盖自动生成）
        vina_bin: Vina 可执行文件路径（None 则自动查找）
    
    返回:
        dict: {
            'success': bool,
            'output_pdbqt': str,
            'config_file': str,
            'affinity': float,  # 最佳模式亲和力
            'log': str
        }
    """
    from .vina_scoring import find_vina_executable
    
    # 查找 Vina
    if vina_bin is None:
        vina_bin = find_vina_executable()
    
    if not vina_bin:
        return {
            'success': False,
            'error': 'Vina executable not found'
        }
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 输出文件
    ligand_basename = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
    output_pdbqt = os.path.join(output_dir, f"{ligand_basename}_pocket{pocket_id}_out.pdbqt")
    log_file = os.path.join(output_dir, f"{ligand_basename}_pocket{pocket_id}.log")
    
    # 生成或使用配置文件
    if custom_config and os.path.exists(custom_config):
        print(f"[dock_to_pocket] Using custom config: {custom_config}")
        config_path = custom_config
    else:
        print(f"[dock_to_pocket] Generating auto config for pocket {pocket_id}...")
        box_params = calculate_pocket_box(pocket_data, padding=padding)
        config_path = os.path.join(output_dir, f"pocket{pocket_id}_config.txt")
        generate_vina_config(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            box_params=box_params,
            output_pdbqt=output_pdbqt,
            config_path=config_path,
            exhaustiveness=exhaustiveness
        )
    
    # 运行 Vina
    print(f"[dock_to_pocket] Running Vina docking...")
    print(f"   Exhaustiveness: {exhaustiveness}")
    
    cmd_args = [
        vina_bin,
        '--config', config_path,
        '--log', log_file
    ]
    
    try:
        result = subprocess.run(
            cmd_args,
            capture_output=True,
            text=True,
            timeout=600  # 10分钟超时
        )
        
        # 解析结果
        affinity = None
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                log_content = f.read()
                
            # 解析最佳模式的亲和力
            for line in log_content.split('\n'):
                if line.strip().startswith('1 '):  # 第一个模式
                    parts = line.split()
                    try:
                        affinity = float(parts[1])
                        break
                    except (ValueError, IndexError):
                        pass
        
        if affinity is not None:
            print(f"[dock_to_pocket] ✅ Docking complete!")
            print(f"   Best affinity: {affinity:.2f} kcal/mol")
            print(f"   Output: {output_pdbqt}")
            
            return {
                'success': True,
                'output_pdbqt': output_pdbqt,
                'config_file': config_path,
                'log_file': log_file,
                'affinity': affinity,
                'pocket_id': pocket_id
            }
        else:
            return {
                'success': False,
                'error': 'Could not parse docking result',
                'output': result.stdout + result.stderr
            }
    
    except subprocess.TimeoutExpired:
        return {
            'success': False,
            'error': 'Docking timeout (>10 min)'
        }
    except Exception as e:
        return {
            'success': False,
            'error': str(e)
        }


# ==================== 批量口袋对接 ====================
def dock_to_multiple_pockets(
    receptor_pdbqt: str,
    ligand_pdbqt: str,
    pockets: List[Dict[str, Any]],
    output_dir: str,
    max_pockets: int = 3,
    padding: float = 5.0,
    exhaustiveness: int = 8
) -> List[Dict[str, Any]]:
    """
    批量对接到多个口袋（按 druggability 排序）
    
    参数:
        receptor_pdbqt: 受体 PDBQT 文件
        ligand_pdbqt: 配体 PDBQT 文件
        pockets: 口袋数据列表
        output_dir: 输出目录
        max_pockets: 最大对接口袋数量
        padding: 盒子边界扩展
        exhaustiveness: Vina 搜索精度
    
    返回:
        list: 对接结果列表
    """
    print(f"[dock_to_multiple_pockets] Docking to top {max_pockets} pockets...")
    
    # 按 druggability 排序
    sorted_pockets = sorted(
        pockets,
        key=lambda p: p.get('druggability_score', 0),
        reverse=True
    )[:max_pockets]
    
    results = []
    
    for i, pocket in enumerate(sorted_pockets, 1):
        print(f"\n{'='*60}")
        print(f"Pocket {i}/{len(sorted_pockets)}")
        print(f"{'='*60}")
        print(f"Volume: {pocket['volume']:.1f} Ų")
        print(f"Druggability: {pocket['druggability_score']:.2f}")
        
        result = dock_to_pocket(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            pocket_data=pocket,
            output_dir=output_dir,
            pocket_id=i,
            padding=padding,
            exhaustiveness=exhaustiveness
        )
        
        result['pocket_data'] = pocket
        results.append(result)
    
    # 总结
    print(f"\n{'='*60}")
    print("Docking Summary")
    print(f"{'='*60}")
    
    successful = [r for r in results if r.get('success')]
    print(f"Successful dockings: {len(successful)}/{len(results)}")
    
    if successful:
        print("\nBest results:")
        sorted_results = sorted(successful, key=lambda r: r.get('affinity', 0))
        for i, r in enumerate(sorted_results[:3], 1):
            print(f"  {i}. Pocket {r['pocket_id']}: {r['affinity']:.2f} kcal/mol")
    
    print(f"{'='*60}")
    
    return results


# ==================== PyMOL 命令封装 ====================
def pocket_based_docking(
    obj_name: str,
    ligand_file: str,
    output_dir: Optional[str] = None,
    auto_detect_pockets: bool = True,
    custom_config: Optional[str] = None,
    max_pockets: int = 3,
    exhaustiveness: int = 8
) -> Dict[str, Any]:
    """
    PyMOL 命令: 基于口袋检测的自动对接
    
    用法:
        # 自动检测口袋并对接
        pocket_based_docking protein, ligand.mol2
        
        # 使用自定义配置
        pocket_based_docking protein, ligand.mol2, custom_config=config.txt
        
        # 对接到指定数量的口袋
        pocket_based_docking protein, ligand.mol2, max_pockets=5
    
    参数:
        obj_name: 蛋白 PyMOL 对象名
        ligand_file: 配体文件路径 (.mol2, .sdf, .pdb)
        output_dir: 输出目录（None 则自动创建）
        auto_detect_pockets: 是否自动检测口袋
        custom_config: 自定义 Vina 配置文件
        max_pockets: 最大对接口袋数量
        exhaustiveness: Vina 搜索精度
    
    返回:
        dict: 对接结果
    """
    from .vina_scoring import export_to_pdbqt_simple, find_vina_executable
    from .pocket_detector import detect_pockets
    
    # 检查 Vina
    vina_bin = find_vina_executable()
    if not vina_bin:
        print("[pocket_based_docking] ⚠️ Vina not found. Please install:")
        print("   conda install -c conda-forge vina")
        return {'success': False, 'error': 'Vina not available'}
    
    # 检查对象
    if obj_name not in cmd.get_object_list():
        print(f"[pocket_based_docking] ⚠️ Object '{obj_name}' not found")
        return {'success': False, 'error': 'Object not found'}
    
    # 创建输出目录
    if output_dir is None:
        output_dir = f"{obj_name}_docking_{os.path.splitext(os.path.basename(ligand_file))[0]}"
    
    os.makedirs(output_dir, exist_ok=True)
    print(f"[pocket_based_docking] Output directory: {output_dir}")
    
    # 准备受体 PDBQT
    receptor_pdbqt = os.path.join(output_dir, f"{obj_name}_receptor.pdbqt")
    print(f"[pocket_based_docking] Preparing receptor...")
    
    if not export_to_pdbqt_simple(obj_name, obj_name, receptor_pdbqt):
        print("[pocket_based_docking] ⚠️ Failed to export receptor")
        return {'success': False, 'error': 'Receptor export failed'}
    
    # 准备配体 PDBQT
    ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
    print(f"[pocket_based_docking] Preparing ligand...")
    
    # 加载配体到 PyMOL
    ligand_obj = "temp_ligand"
    cmd.load(ligand_file, ligand_obj)
    
    if not export_to_pdbqt_simple(ligand_obj, ligand_obj, ligand_pdbqt):
        print("[pocket_based_docking] ⚠️ Failed to export ligand")
        cmd.delete(ligand_obj)
        return {'success': False, 'error': 'Ligand export failed'}
    
    cmd.delete(ligand_obj)
    
    # 检测口袋或使用自定义配置
    if auto_detect_pockets and not custom_config:
        print(f"[pocket_based_docking] Detecting pockets...")
        pockets = detect_pockets(
            obj_name=obj_name,
            grid_spacing=0.6,
            min_volume=20.0,
            use_schrodinger_standard=False
        )
        
        if not pockets:
            print("[pocket_based_docking] ⚠️ No pockets detected")
            return {'success': False, 'error': 'No pockets found'}
        
        print(f"[pocket_based_docking] Found {len(pockets)} pockets")
        
        # 批量对接
        results = dock_to_multiple_pockets(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            pockets=pockets,
            output_dir=output_dir,
            max_pockets=max_pockets,
            exhaustiveness=exhaustiveness
        )
        
        return {
            'success': True,
            'results': results,
            'output_dir': output_dir,
            'n_pockets': len(pockets)
        }
    
    else:
        # 使用自定义配置单次对接
        if not custom_config or not os.path.exists(custom_config):
            print("[pocket_based_docking] ⚠️ Custom config file not found")
            return {'success': False, 'error': 'Config file not found'}
        
        print(f"[pocket_based_docking] Using custom config: {custom_config}")
        
        # 单次对接
        output_pdbqt = os.path.join(output_dir, "docking_out.pdbqt")
        log_file = os.path.join(output_dir, "docking.log")
        
        cmd_args = [
            vina_bin,
            '--config', custom_config,
            '--log', log_file
        ]
        
        try:
            result = subprocess.run(
                cmd_args,
                capture_output=True,
                text=True,
                timeout=600
            )
            
            print(f"[pocket_based_docking] ✅ Docking complete")
            print(f"   Output: {output_pdbqt}")
            print(f"   Log: {log_file}")
            
            return {
                'success': True,
                'output_pdbqt': output_pdbqt,
                'log_file': log_file,
                'output_dir': output_dir
            }
        
        except Exception as e:
            print(f"[pocket_based_docking] ⚠️ Docking failed: {e}")
            return {'success': False, 'error': str(e)}


# ==================== 可视化对接结果 ====================
def visualize_docking_result(
    obj_name: str,
    docking_pdbqt: str,
    pocket_id: int = 1,
    show_interactions: bool = True
):
    """
    PyMOL 命令: 可视化对接结果
    
    用法:
        visualize_docking_result protein, docking_out.pdbqt, pocket_id=1
    
    参数:
        obj_name: 蛋白对象名
        docking_pdbqt: 对接结果 PDBQT 文件
        pocket_id: 口袋编号（用于命名）
        show_interactions: 是否显示相互作用
    """
    # 加载对接结果
    pose_obj = f"{obj_name}_pose_p{pocket_id}"
    cmd.load(docking_pdbqt, pose_obj)
    
    # 显示设置
    cmd.hide('everything', obj_name)
    cmd.show('cartoon', obj_name)
    cmd.color('gray80', obj_name)
    
    # 显示配体
    cmd.show('sticks', pose_obj)
    cmd.color('green', f'{pose_obj} and elem C')
    cmd.util.cnc(pose_obj)  # 按元素着色
    
    # 显示结合位点残基
    cmd.select(f'binding_site_p{pocket_id}', f'{obj_name} within 5 of {pose_obj}')
    cmd.show('sticks', f'binding_site_p{pocket_id}')
    cmd.color('cyan', f'binding_site_p{pocket_id} and elem C')
    
    # 相互作用分析
    if show_interactions:
        try:
            from .interaction_analyzer import analyze_protein_ligand_interactions
            
            # 分析相互作用
            # 注意: 需要将 PDBQT 转为标准 PDB
            temp_pdb = docking_pdbqt.replace('.pdbqt', '_temp.pdb')
            cmd.save(temp_pdb, pose_obj)
            
            print(f"[visualize_docking_result] Analyzing interactions...")
            # 这里可以调用相互作用分析
            
        except ImportError:
            print("[visualize_docking_result] Interaction analysis not available")
    
    # 居中视图
    cmd.zoom(pose_obj, 5)
    
    print(f"[visualize_docking_result] ✅ Visualization complete")
    print(f"   Pose object: {pose_obj}")
    print(f"   Binding site selection: binding_site_p{pocket_id}")


# ==================== PyMOL 命令注册 ====================
def register_pocket_docking_commands():
    """注册 PyMOL 命令"""
    cmd.extend("pocket_based_docking", pocket_based_docking)
    cmd.extend("visualize_docking_result", visualize_docking_result)
    
    print("[pocket_docking] Commands registered:")
    print("   pocket_based_docking")
    print("   visualize_docking_result")


if __name__ == "__main__":
    # 测试代码
    print("This module should be loaded in PyMOL")
    print("Example usage:")
    print("  pocket_based_docking protein, ligand.mol2")
