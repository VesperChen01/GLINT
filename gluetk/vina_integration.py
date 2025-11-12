# -*- coding: utf-8 -*-
"""
vina_integration.py
GlueTK - AutoDock Vina 完整集成模块

整合功能：
1. Vina 可执行文件检测与 PDBQT 导出
2. 分子对接评分（score_only, local_only）
3. 基于口袋的自动对接
4. 批量对接与评分对比
5. 与经验评分函数对比

依赖：
- AutoDock Vina (conda install -c conda-forge vina)
- Open Babel 或 MGLTools (PDBQT 转换)
"""

from __future__ import print_function
import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import List, Dict, Any, Optional
from pymol import cmd


# ==================== 工具检测 ====================
def find_vina_executable():
    """
    自动查找 Vina 可执行文件
    
    返回:
        str: Vina 可执行文件路径，未找到返回 None
    """
    vina = shutil.which('vina') or shutil.which('vina.exe')
    if vina:
        return vina
    
    if sys.platform == 'win32':
        common_paths = [
            r'C:\Program Files\vina\vina.exe',
            r'C:\Program Files (x86)\vina\vina.exe',
            os.path.join(os.environ.get('LOCALAPPDATA', ''), 'Programs', 'VinaTools', 'vina', 'vina.exe'),
        ]
        for path in common_paths:
            if os.path.exists(path):
                return path
    
    if 'CONDA_PREFIX' in os.environ:
        conda_vina = os.path.join(os.environ['CONDA_PREFIX'], 'bin', 'vina')
        if os.path.exists(conda_vina):
            return conda_vina
    
    return None


def find_obabel_executable():
    """查找 Open Babel 可执行文件"""
    return shutil.which('obabel') or shutil.which('obabel.exe')


def find_mgltools_scripts():
    """查找 MGLTools 脚本路径"""
    # 常见 MGLTools 安装位置
    possible_paths = [
        '/usr/local/MGLTools',
        '/opt/mgltools',
        os.path.join(os.path.expanduser('~'), 'MGLTools'),
    ]
    
    if 'MGLTOOLS_HOME' in os.environ:
        possible_paths.insert(0, os.environ['MGLTOOLS_HOME'])
    
    for base_path in possible_paths:
        prepare_ligand = os.path.join(base_path, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_ligand4.py')
        prepare_receptor = os.path.join(base_path, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_receptor4.py')
        
        if os.path.exists(prepare_ligand) and os.path.exists(prepare_receptor):
            return {
                'prepare_ligand': prepare_ligand,
                'prepare_receptor': prepare_receptor,
                'python': os.path.join(base_path, 'bin', 'pythonsh')
            }
    
    return None


def check_vina_available():
    """检查 Vina 是否可用"""
    vina = find_vina_executable()
    if vina:
        print(f"[vina_integration] ✅ Found Vina: {vina}")
        return True
    else:
        print("[vina_integration] ⚠️ Vina not found")
        print("[vina_integration] Install: conda install -c conda-forge vina")
        return False


# ==================== PDBQT 转换（标准流程）====================
def export_to_pdbqt(obj_name, selection, output_pdbqt, is_receptor=True):
    """
    标准 PDBQT 导出（使用 Open Babel 或 MGLTools）
    
    优先级：
    1. Open Babel (obabel) - 最推荐
    2. MGLTools (prepare_ligand4.py / prepare_receptor4.py)
    
    参数:
        obj_name: PyMOL 对象名
        selection: PyMOL 选择表达式
        output_pdbqt: 输出 PDBQT 文件路径
        is_receptor: 是否为受体（True=蛋白，False=配体）
    
    返回:
        bool: 是否成功
    """
    # 先导出为 PDB
    temp_pdb = output_pdbqt.replace('.pdbqt', '_temp.pdb')
    cmd.save(temp_pdb, selection)
    
    if not os.path.exists(temp_pdb):
        print(f"[export_to_pdbqt] ⚠️ Failed to save PDB: {temp_pdb}")
        return False
    
    # 方法 1: 使用 Open Babel
    obabel = find_obabel_executable()
    if obabel:
        try:
            # 为配体添加氢原子，为受体移除水分子
            cmd_args = [obabel, temp_pdb, '-O', output_pdbqt, '-xh']
            
            if is_receptor:
                # 受体：移除氢原子，保留极性氢
                cmd_args.extend(['-d'])  # 删除氢原子
            else:
                # 配体：添加氢原子
                cmd_args.extend(['-h'])  # 添加氢原子
            
            result = subprocess.run(cmd_args, capture_output=True, text=True, timeout=60)
            
            if result.returncode == 0 and os.path.exists(output_pdbqt):
                print(f"[export_to_pdbqt] ✅ Converted using Open Babel")
                os.unlink(temp_pdb)
                return True
            else:
                print(f"[export_to_pdbqt] ⚠️ Open Babel conversion failed: {result.stderr}")
        except Exception as e:
            print(f"[export_to_pdbqt] ⚠️ Open Babel error: {e}")
    
    # 方法 2: 使用 MGLTools
    mgltools = find_mgltools_scripts()
    if mgltools:
        try:
            script = mgltools['prepare_receptor'] if is_receptor else mgltools['prepare_ligand']
            python_exe = mgltools['python']
            
            cmd_args = [python_exe, script, '-r' if is_receptor else '-l', temp_pdb, '-o', output_pdbqt]
            
            result = subprocess.run(cmd_args, capture_output=True, text=True, timeout=60)
            
            if result.returncode == 0 and os.path.exists(output_pdbqt):
                print(f"[export_to_pdbqt] ✅ Converted using MGLTools")
                os.unlink(temp_pdb)
                return True
            else:
                print(f"[export_to_pdbqt] ⚠️ MGLTools conversion failed: {result.stderr}")
        except Exception as e:
            print(f"[export_to_pdbqt] ⚠️ MGLTools error: {e}")
    
    # 工具缺失：提示使用环境检查
    print("\n" + "=" * 60)
    print("⚠️  PDBQT Conversion Tools Not Found")
    print("=" * 60)
    print("PDBQT conversion requires Open Babel or MGLTools.")
    print("\nPlease run the environment checker:")
    print("  python gluetk/env_checker.py --auto-install")
    print("\nOr install manually:")
    print("  conda install -c conda-forge openbabel")
    print("=" * 60)
    
    if os.path.exists(temp_pdb):
        os.unlink(temp_pdb)
    
    return False


# ==================== Vina 评分核心 ====================
def score_with_vina(protein_pdbqt, ligand_pdbqt, vina_bin=None, mode='score_only'):
    """
    使用 Vina 评分复合物
    
    参数:
        protein_pdbqt: 蛋白 PDBQT 文件路径
        ligand_pdbqt: 配体 PDBQT 文件路径
        vina_bin: Vina 可执行文件路径（None 则自动查找）
        mode: 'score_only' 或 'local_only'
    
    返回:
        dict: {'affinity': float, 'success': bool, 'output': str}
    """
    if vina_bin is None:
        vina_bin = find_vina_executable()
    
    if not vina_bin:
        return {'success': False, 'error': 'Vina not found', 'affinity': None}
    
    cmd_args = [
        vina_bin,
        '--receptor', protein_pdbqt,
        '--ligand', ligand_pdbqt,
    ]
    
    if mode == 'score_only':
        cmd_args.append('--score_only')
    elif mode == 'local_only':
        cmd_args.append('--local_only')
    
    try:
        result = subprocess.run(
            cmd_args,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        output = result.stdout + result.stderr
        
        affinity = None
        for line in output.split('\n'):
            if 'Affinity:' in line or 'kcal/mol' in line:
                parts = line.split()
                for part in parts:
                    try:
                        val = float(part)
                        affinity = val
                        break
                    except ValueError:
                        continue
                if affinity is not None:
                    break
        
        if affinity is None:
            return {'success': False, 'error': 'Could not parse affinity', 'output': output, 'affinity': None}
        
        return {'success': True, 'affinity': affinity, 'output': output}
    
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Timeout', 'affinity': None}
    except Exception as e:
        return {'success': False, 'error': str(e), 'affinity': None}


def vina_score_complex(protein_obj, ligand_selection, mode='score_only', show_report=True):
    """
    PyMOL 命令：使用 Vina 评分蛋白-配体复合物
    
    用法:
        vina_score_complex protein, resn LIG
        vina_score_complex protein, resn LIG, mode=local_only
    """
    if not check_vina_available():
        print("[vina_score_complex] ⚠️ Vina not available")
        return None
    
    if protein_obj not in cmd.get_object_list():
        print(f"[vina_score_complex] ⚠️ Object '{protein_obj}' not found")
        return None
    
    with tempfile.TemporaryDirectory() as tmpdir:
        protein_pdbqt = os.path.join(tmpdir, 'protein.pdbqt')
        ligand_pdbqt = os.path.join(tmpdir, 'ligand.pdbqt')
        
        print(f"[vina_score_complex] Exporting protein: {protein_obj}")
        if not export_to_pdbqt(protein_obj, protein_obj, protein_pdbqt, is_receptor=True):
            print("[vina_score_complex] ⚠️ Failed to export protein")
            return None
        
        print(f"[vina_score_complex] Exporting ligand: {ligand_selection}")
        if not export_to_pdbqt(protein_obj, ligand_selection, ligand_pdbqt, is_receptor=False):
            print("[vina_score_complex] ⚠️ Failed to export ligand")
            return None
        
        print(f"[vina_score_complex] Running Vina ({mode})...")
        result = score_with_vina(protein_pdbqt, ligand_pdbqt, mode=mode)
    
    if show_report and result['success']:
        print("=" * 60)
        print("Vina Scoring Report")
        print("=" * 60)
        print(f"Protein: {protein_obj}")
        print(f"Ligand:  {ligand_selection}")
        print(f"Mode:    {mode}")
        print(f"\nAffinity: {result['affinity']:.2f} kcal/mol")
        print("=" * 60)
    elif not result['success']:
        print(f"[vina_score_complex] ⚠️ Vina scoring failed: {result.get('error', 'Unknown error')}")
    
    return result


def compare_scoring_methods(protein_obj, ligand_resname):
    """
    PyMOL 命令：对比经验评分 vs Vina 评分 - DEPRECATED
    
    用法:
        compare_scoring_methods protein, LIG
    """
    print("=" * 60)
    print("Scoring Method Comparison - DEPRECATED")
    print("=" * 60)
    print("\nThis function has been deprecated.")
    print("Empirical scoring (binding_score.py) has been removed.")
    print("\nPlease use:")
    print("  • vina_score_complex() for Vina scoring")
    print("  • analyze_protein_ligand_interactions() for interaction analysis")
    print("=" * 60)
    return None


# ==================== 口袋对接功能 ====================
def calculate_pocket_box(pocket_data: Dict[str, Any], padding: float = 5.0) -> Dict[str, float]:
    """
    从口袋数据计算对接盒子参数
    
    参数:
        pocket_data: 口袋数据字典，包含 grid_points 或 center
        padding: 盒子边界扩展 (Å)
    
    返回:
        dict: {center_x, center_y, center_z, size_x, size_y, size_z}
    """
    import numpy as np
    
    if 'grid_coords' in pocket_data:
        grid_points = np.array(pocket_data['grid_coords'])
        
        xs = grid_points[:, 0]
        ys = grid_points[:, 1]
        zs = grid_points[:, 2]
        
        min_x, max_x = xs.min(), xs.max()
        min_y, max_y = ys.min(), ys.max()
        min_z, max_z = zs.min(), zs.max()
        
        center_x = (min_x + max_x) / 2.0
        center_y = (min_y + max_y) / 2.0
        center_z = (min_z + max_z) / 2.0
        
        size_x = (max_x - min_x) + 2 * padding
        size_y = (max_y - min_y) + 2 * padding
        size_z = (max_z - min_z) + 2 * padding
    
    elif 'center' in pocket_data:
        center = pocket_data['center']
        center_x, center_y, center_z = center
        
        volume = pocket_data.get('volume', 100.0)
        approx_radius = (3 * volume / (4 * 3.14159)) ** (1/3)
        
        size_x = size_y = size_z = 2 * approx_radius + 2 * padding
    
    else:
        raise ValueError("Pocket data must contain 'grid_coords' or 'center'")
    
    return {
        'center_x': center_x,
        'center_y': center_y,
        'center_z': center_z,
        'size_x': size_x,
        'size_y': size_y,
        'size_z': size_z
    }


def generate_vina_config(receptor_pdbqt: str, ligand_pdbqt: str, box_params: Dict[str, float],
                        output_pdbqt: str, config_path: Optional[str] = None,
                        exhaustiveness: int = 8, num_modes: int = 9) -> str:
    """
    生成 Vina 配置文件
    
    参数:
        receptor_pdbqt: 受体 PDBQT 文件路径
        ligand_pdbqt: 配体 PDBQT 文件路径
        box_params: 盒子参数字典
        output_pdbqt: 输出 PDBQT 文件路径
        config_path: 配置文件保存路径
        exhaustiveness: 搜索精度
        num_modes: 输出模式数量
    
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


def dock_to_pocket(receptor_pdbqt: str, ligand_pdbqt: str, pocket_data: Dict[str, Any],
                   output_dir: str, pocket_id: int = 1, padding: float = 5.0,
                   exhaustiveness: int = 8, vina_bin: Optional[str] = None) -> Dict[str, Any]:
    """
    对接配体到指定口袋
    
    参数:
        receptor_pdbqt: 受体 PDBQT 文件
        ligand_pdbqt: 配体 PDBQT 文件
        pocket_data: 口袋数据字典
        output_dir: 输出目录
        pocket_id: 口袋编号
        padding: 盒子边界扩展
        exhaustiveness: Vina 搜索精度
        vina_bin: Vina 可执行文件路径
    
    返回:
        dict: 对接结果
    """
    if vina_bin is None:
        vina_bin = find_vina_executable()
    
    if not vina_bin:
        return {'success': False, 'error': 'Vina executable not found'}
    
    os.makedirs(output_dir, exist_ok=True)
    
    ligand_basename = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
    output_pdbqt = os.path.join(output_dir, f"{ligand_basename}_pocket{pocket_id}_out.pdbqt")
    log_file = os.path.join(output_dir, f"{ligand_basename}_pocket{pocket_id}.log")
    
    print(f"[dock_to_pocket] Generating config for pocket {pocket_id}...")
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
    
    print(f"[dock_to_pocket] Running Vina docking...")
    cmd_args = [vina_bin, '--config', config_path, '--log', log_file]
    
    try:
        result = subprocess.run(cmd_args, capture_output=True, text=True, timeout=600)
        
        affinity = None
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                for line in f:
                    if line.strip().startswith('1 '):
                        parts = line.split()
                        try:
                            affinity = float(parts[1])
                            break
                        except (ValueError, IndexError):
                            pass
        
        if affinity is not None:
            print(f"[dock_to_pocket] ✅ Docking complete!")
            print(f"   Best affinity: {affinity:.2f} kcal/mol")
            
            return {
                'success': True,
                'output_pdbqt': output_pdbqt,
                'config_file': config_path,
                'log_file': log_file,
                'affinity': affinity,
                'pocket_id': pocket_id
            }
        else:
            return {'success': False, 'error': 'Could not parse docking result'}
    
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Docking timeout (>10 min)'}
    except Exception as e:
        return {'success': False, 'error': str(e)}


def pocket_based_docking(obj_name: str, ligand_file: str, output_dir: Optional[str] = None,
                        max_pockets: int = 3, exhaustiveness: int = 8) -> Dict[str, Any]:
    """
    PyMOL 命令：基于口袋检测的自动对接
    
    用法:
        pocket_based_docking protein, ligand.mol2
        pocket_based_docking protein, ligand.mol2, max_pockets=5
    """
    from .pocket_detector import detect_pockets
    
    if not check_vina_available():
        print("[pocket_based_docking] ⚠️ Vina not available")
        return {'success': False, 'error': 'Vina not available'}
    
    if obj_name not in cmd.get_object_list():
        print(f"[pocket_based_docking] ⚠️ Object '{obj_name}' not found")
        return {'success': False, 'error': 'Object not found'}
    
    if output_dir is None:
        output_dir = f"{obj_name}_docking_{os.path.splitext(os.path.basename(ligand_file))[0]}"
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 准备受体
    receptor_pdbqt = os.path.join(output_dir, f"{obj_name}_receptor.pdbqt")
    if not export_to_pdbqt(obj_name, obj_name, receptor_pdbqt, is_receptor=True):
        return {'success': False, 'error': 'Receptor export failed'}
    
    # 准备配体
    ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
    ligand_obj = "temp_ligand"
    cmd.load(ligand_file, ligand_obj)
    
    if not export_to_pdbqt(ligand_obj, ligand_obj, ligand_pdbqt, is_receptor=False):
        cmd.delete(ligand_obj)
        return {'success': False, 'error': 'Ligand export failed'}
    
    cmd.delete(ligand_obj)
    
    # 检测口袋
    print(f"[pocket_based_docking] Detecting pockets...")
    pockets = detect_pockets(obj_name=obj_name)
    
    if not pockets:
        return {'success': False, 'error': 'No pockets found'}
    
    # 按 druggability 排序
    sorted_pockets = sorted(pockets, key=lambda p: p.get('druggability_score', 0), reverse=True)[:max_pockets]
    
    results = []
    for i, pocket in enumerate(sorted_pockets, 1):
        print(f"\n{'='*60}")
        print(f"Docking to Pocket {i}/{len(sorted_pockets)}")
        print(f"Volume: {pocket['volume']:.1f} Ų, Druggability: {pocket['druggability_score']:.2f}")
        print(f"{'='*60}")
        
        result = dock_to_pocket(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            pocket_data=pocket,
            output_dir=output_dir,
            pocket_id=i,
            padding=5.0,
            exhaustiveness=exhaustiveness
        )
        
        result['pocket_data'] = pocket
        results.append(result)
    
    successful = [r for r in results if r.get('success')]
    print(f"\n{'='*60}")
    print(f"Docking Summary: {len(successful)}/{len(results)} successful")
    if successful:
        best = min(successful, key=lambda r: r['affinity'])
        print(f"Best result: Pocket {best['pocket_id']}, {best['affinity']:.2f} kcal/mol")
    print(f"{'='*60}")
    
    return {'success': True, 'results': results, 'output_dir': output_dir}


def manual_box_docking(obj_name: str, ligand_file: str, box_params: Dict[str, float],
                       output_dir: Optional[str] = None, exhaustiveness: int = 8) -> Dict[str, Any]:
    """
    使用自定义对接盒参数进行 Vina 对接（不做口袋检测）

    参数:
        obj_name: 受体的 PyMOL 对象名
        ligand_file: 配体文件 (MOL2/SDF/PDBQT)
        box_params: {'center_x','center_y','center_z','size_x','size_y','size_z'}
        output_dir: 输出目录
        exhaustiveness: Vina 搜索精度
    返回:
        {'success': bool, 'output_dir': str, 'result': dict}
    """
    if not check_vina_available():
        return {'success': False, 'error': 'Vina not available'}

    if obj_name not in cmd.get_object_list():
        return {'success': False, 'error': f"Object '{obj_name}' not found"}

    if output_dir is None:
        output_dir = f"{obj_name}_manual_docking_{os.path.splitext(os.path.basename(ligand_file))[0]}"
    os.makedirs(output_dir, exist_ok=True)

    # 导出受体 PDBQT
    receptor_pdbqt = os.path.join(output_dir, f"{obj_name}_receptor.pdbqt")
    if not export_to_pdbqt(obj_name, obj_name, receptor_pdbqt, is_receptor=True):
        return {'success': False, 'error': 'Receptor export failed'}

    # 导出配体 PDBQT
    ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
    ligand_obj = "temp_ligand"
    cmd.load(ligand_file, ligand_obj)
    if not export_to_pdbqt(ligand_obj, ligand_obj, ligand_pdbqt, is_receptor=False):
        cmd.delete(ligand_obj)
        return {'success': False, 'error': 'Ligand export failed'}
    cmd.delete(ligand_obj)

    # 生成配置并运行对接
    ligand_basename = os.path.splitext(os.path.basename(ligand_pdbqt))[0]
    output_pdbqt = os.path.join(output_dir, f"{ligand_basename}_manual_out.pdbqt")
    log_file = os.path.join(output_dir, f"{ligand_basename}_manual.log")

    config_path = os.path.join(output_dir, f"manual_config.txt")
    generate_vina_config(
        receptor_pdbqt=receptor_pdbqt,
        ligand_pdbqt=ligand_pdbqt,
        box_params=box_params,
        output_pdbqt=output_pdbqt,
        config_path=config_path,
        exhaustiveness=exhaustiveness
    )

    vina_bin = find_vina_executable()
    cmd_args = [vina_bin, '--config', config_path, '--log', log_file]
    try:
        subprocess.run(cmd_args, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Docking timeout (>10 min)'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

    # 解析打分
    affinity = None
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            for line in f:
                if line.strip().startswith('1 '):
                    parts = line.split()
                    try:
                        affinity = float(parts[1])
                        break
                    except (ValueError, IndexError):
                        pass

    if affinity is None:
        return {'success': False, 'error': 'Could not parse docking result'}

    return {
        'success': True,
        'output_dir': output_dir,
        'results': [{
            'success': True,
            'output_pdbqt': output_pdbqt,
            'config_file': config_path,
            'log_file': log_file,
            'affinity': affinity,
            'pocket_id': 'manual'
        }]
    }


if __name__ == "__main__":
    print("[vina_integration] This is a PyMOL plugin module")
    print("Commands: vina_score_complex, compare_scoring_methods, pocket_based_docking, manual_box_docking")
