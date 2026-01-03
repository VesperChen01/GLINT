# -*- coding: utf-8 -*-
"""
vina_integration.py
GlueTK - AutoDock Vina 完整集成模块
"""

from __future__ import print_function
import os
import sys
import subprocess
import tempfile
import shutil
import glob
from pathlib import Path
from typing import List, Dict, Any, Optional, Callable
from pymol import cmd

# 版本标记 - 用于确认代码是否被正确加载
_VINA_MODULE_VERSION = "2026-01-03-v2"
print(f"[vina_integration] 模块版本: {_VINA_MODULE_VERSION}")


def find_vina_executable():
    """自动查找 Vina 可执行文件"""
    vina = shutil.which('vina') or shutil.which('vina.exe')
    if vina:
        return vina
    
    if 'CONDA_PREFIX' in os.environ:
        conda_vina = os.path.join(os.environ['CONDA_PREFIX'], 'bin', 'vina')
        if os.path.exists(conda_vina):
            return conda_vina
    
    home = os.path.expanduser('~')
    user_paths = [
        os.path.join(home, 'bin', 'vina'),
        os.path.join(home, '.local', 'bin', 'vina'),
        '/usr/local/bin/vina',
        '/opt/homebrew/bin/vina',
    ]
    
    for path in user_paths:
        if '*' in path:
            matches = glob.glob(path)
            if matches:
                return matches[0]
        elif os.path.exists(path):
            return path
    return None


def find_obabel_executable():
    """查找 Open Babel 可执行文件"""
    if os.environ.get("OBABEL_BINARY"):
        return os.environ.get("OBABEL_BINARY")
    
    obabel = shutil.which("obabel") or shutil.which("obabel.exe")
    if obabel:
        return obabel
    
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        conda_obabel = os.path.join(conda_prefix, "bin", "obabel")
        if os.path.exists(conda_obabel):
            return conda_obabel
    
    home = os.path.expanduser("~")
    for env_name in ["gluetk", "base"]:
        for base in [f"{home}/miniconda3", f"{home}/anaconda3", "/opt/homebrew/Caskroom/miniconda/base"]:
            obabel_path = f"{base}/envs/{env_name}/bin/obabel"
            if os.path.exists(obabel_path):
                return obabel_path
    return None


def find_mgltools_scripts():
    """查找 MGLTools 脚本路径"""
    possible_paths = ['/usr/local/MGLTools', '/opt/mgltools', os.path.join(os.path.expanduser('~'), 'MGLTools')]
    if 'MGLTOOLS_HOME' in os.environ:
        possible_paths.insert(0, os.environ['MGLTOOLS_HOME'])
    
    for base_path in possible_paths:
        prepare_ligand = os.path.join(base_path, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_ligand4.py')
        prepare_receptor = os.path.join(base_path, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24', 'prepare_receptor4.py')
        if os.path.exists(prepare_ligand) and os.path.exists(prepare_receptor):
            return {'prepare_ligand': prepare_ligand, 'prepare_receptor': prepare_receptor,
                    'python': os.path.join(base_path, 'bin', 'pythonsh')}
    return None


def check_vina_available():
    """检查 Vina 是否可用"""
    vina = find_vina_executable()
    if vina:
        print(f"[vina_integration] ✅ Found Vina: {vina}")
        return True
    print("[vina_integration] ⚠️ Vina not found. Install: conda install -c conda-forge vina")
    return False


def export_to_pdbqt(obj_name, selection, output_pdbqt, is_receptor=True):
    """标准 PDBQT 导出
    
    对于受体：使用 -xr 参数生成刚性受体格式（无 ROOT/ENDROOT 标签）
    对于配体：使用 -h 参数添加氢原子
    """
    print(f"[vina_integration] export_to_pdbqt: is_receptor={is_receptor}, output={output_pdbqt}")
    
    temp_pdb = output_pdbqt.replace('.pdbqt', '_temp.pdb')
    cmd.save(temp_pdb, selection)
    
    if not os.path.exists(temp_pdb):
        print(f"[vina_integration] ❌ 临时 PDB 文件创建失败")
        return False
    
    obabel = find_obabel_executable()
    if obabel:
        try:
            if is_receptor:
                # 受体：使用 -xr 生成刚性受体格式，不包含 ROOT 标签
                cmd_args = [obabel, temp_pdb, '-O', output_pdbqt, '-xr']
                print(f"[vina_integration] 使用 -xr 参数（刚性受体格式）")
            else:
                # 配体：添加氢原子
                cmd_args = [obabel, temp_pdb, '-O', output_pdbqt, '-h']
                print(f"[vina_integration] 使用 -h 参数（配体格式）")
            
            result = subprocess.run(cmd_args, capture_output=True, text=True, timeout=60)
            if result.returncode == 0 and os.path.exists(output_pdbqt):
                os.unlink(temp_pdb)
                print(f"[vina_integration] ✅ PDBQT 导出成功")
                return True
            else:
                print(f"[vina_integration] ❌ obabel 失败: {result.stderr}")
        except Exception as e:
            print(f"[vina_integration] ❌ 异常: {e}")
    
    if os.path.exists(temp_pdb):
        os.unlink(temp_pdb)
    return False


def convert_ligand_to_pdbqt(ligand_file: str, output_pdbqt: str) -> bool:
    """将配体文件转换为 PDBQT（不通过 PyMOL）"""
    # 如果源文件和目标文件相同，直接返回成功
    if os.path.abspath(ligand_file) == os.path.abspath(output_pdbqt):
        return True
    obabel = find_obabel_executable()
    if not obabel:
        return False
    try:
        result = subprocess.run([obabel, ligand_file, '-O', output_pdbqt, '-h'], 
                               capture_output=True, text=True, timeout=60)
        return result.returncode == 0 and os.path.exists(output_pdbqt)
    except Exception:
        return False


def generate_vina_config(receptor_pdbqt: str, ligand_pdbqt: str, box_params: Dict[str, float],
                        output_pdbqt: str, config_path: Optional[str] = None,
                        exhaustiveness: int = 8, num_modes: int = 9) -> str:
    """生成 Vina 配置文件"""
    if config_path is None:
        fd, config_path = tempfile.mkstemp(suffix='_vina_config.txt', text=True)
        os.close(fd)
    
    config_content = f"""receptor = {os.path.abspath(receptor_pdbqt)}
ligand = {os.path.abspath(ligand_pdbqt)}
out = {os.path.abspath(output_pdbqt)}
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
    return config_path


def run_vina_docking(receptor_pdbqt: str, ligand_pdbqt: str, box_params: Dict[str, float],
                     output_dir: str, ligand_name: str, exhaustiveness: int = 8,
                     num_modes: int = 9) -> Dict[str, Any]:
    """运行单个 Vina 对接"""
    vina_bin = find_vina_executable()
    if not vina_bin:
        return {'success': False, 'error': 'Vina not found'}
    
    output_pdbqt = os.path.join(output_dir, f"{ligand_name}_out.pdbqt")
    log_file = os.path.join(output_dir, f"{ligand_name}.log")
    config_path = os.path.join(output_dir, f"{ligand_name}_config.txt")
    
    generate_vina_config(receptor_pdbqt, ligand_pdbqt, box_params, output_pdbqt,
                        config_path, exhaustiveness, num_modes)
    
    try:
        print(f"[vina_integration] 运行 Vina: {vina_bin}")
        # 新版 Vina 不支持 --log 参数，输出到 stdout
        result = subprocess.run([vina_bin, '--config', config_path],
                      capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            error_msg = result.stderr.strip() if result.stderr else 'Vina execution failed'
            print(f"[vina_integration] ❌ Vina 失败: {error_msg}")
            return {'success': False, 'error': error_msg, 'ligand': ligand_name}

        # 保存输出到日志文件
        with open(log_file, 'w') as f:
            f.write(result.stdout)

    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Timeout', 'ligand': ligand_name}
    except Exception as e:
        return {'success': False, 'error': str(e), 'ligand': ligand_name}

    # 从 stdout 或日志文件中解析结果
    affinity = None
    output_text = result.stdout if result.stdout else ""
    if os.path.exists(log_file):
        with open(log_file, 'r') as f:
            output_text = f.read()

    for line in output_text.split('\n'):
        if line.strip().startswith('1 '):
            parts = line.split()
            try:
                affinity = float(parts[1])
                break
            except (ValueError, IndexError):
                pass

    if affinity is None:
        return {'success': False, 'error': 'Parse failed', 'ligand': ligand_name}
    
    print(f"[vina_integration] ✅ 对接成功: {affinity:.2f} kcal/mol")
    return {'success': True, 'affinity': affinity, 'output_pdbqt': output_pdbqt,
            'log_file': log_file, 'ligand': ligand_name}


def manual_box_docking(obj_name: str, ligand_file: str, box_params: Dict[str, float],
                       output_dir: Optional[str] = None, exhaustiveness: int = 8,
                       num_modes: int = 9, remove_selection: Optional[str] = None) -> Dict[str, Any]:
    """使用自定义对接盒参数进行 Vina 对接"""
    if not check_vina_available():
        return {'success': False, 'error': 'Vina not available'}

    try:
        objs = cmd.get_names("objects")
    except AttributeError:
        objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []

    if obj_name not in objs:
        return {'success': False, 'error': f"Object '{obj_name}' not found"}

    if output_dir is None:
        output_dir = f"{obj_name}_docking"
    os.makedirs(output_dir, exist_ok=True)

    # 导出受体 PDBQT
    receptor_pdbqt = os.path.join(output_dir, f"{obj_name}_receptor.pdbqt")
    receptor_sel = f"({obj_name}) and not ({remove_selection})" if remove_selection else obj_name
    
    if not export_to_pdbqt(obj_name, receptor_sel, receptor_pdbqt, is_receptor=True):
        return {'success': False, 'error': 'Receptor export failed'}

    # 导出配体 PDBQT
    ligand_pdbqt = os.path.join(output_dir, "ligand.pdbqt")
    ligand_obj = "temp_ligand"
    cmd.load(ligand_file, ligand_obj)
    if not export_to_pdbqt(ligand_obj, ligand_obj, ligand_pdbqt, is_receptor=False):
        cmd.delete(ligand_obj)
        return {'success': False, 'error': 'Ligand export failed'}
    cmd.delete(ligand_obj)

    ligand_name = os.path.splitext(os.path.basename(ligand_file))[0]
    result = run_vina_docking(receptor_pdbqt, ligand_pdbqt, box_params, output_dir,
                              ligand_name, exhaustiveness, num_modes)
    
    if result['success']:
        return {'success': True, 'output_dir': output_dir, 'results': [result]}
    return result


# ==================== 批量对接功能 ====================
def get_ligand_files(ligand_path: str) -> List[str]:
    """获取配体文件列表（支持单文件或文件夹）"""
    SUPPORTED_EXTS = {'.mol2', '.sdf', '.pdbqt', '.mol', '.pdb'}
    
    if os.path.isfile(ligand_path):
        return [ligand_path]
    elif os.path.isdir(ligand_path):
        files = []
        for f in os.listdir(ligand_path):
            if os.path.splitext(f)[1].lower() in SUPPORTED_EXTS:
                files.append(os.path.join(ligand_path, f))
        return sorted(files)
    return []


def batch_docking(obj_name: str, ligand_path: str, box_params: Dict[str, float],
                  output_dir: Optional[str] = None, exhaustiveness: int = 8,
                  num_modes: int = 9, remove_selection: Optional[str] = None,
                  progress_callback: Optional[Callable[[int, int, str], None]] = None) -> Dict[str, Any]:
    """
    批量对接：支持单个配体文件或配体文件夹
    
    参数:
        obj_name: 受体的 PyMOL 对象名
        ligand_path: 配体文件或文件夹路径
        box_params: 对接盒参数
        output_dir: 输出目录
        exhaustiveness: Vina 搜索精度
        num_modes: 输出构象数量
        remove_selection: 要移除的共晶配体 selection
        progress_callback: 进度回调函数 (current, total, ligand_name)
    """
    print(f"[vina_integration] batch_docking 开始，模块版本: {_VINA_MODULE_VERSION}")
    
    if not check_vina_available():
        return {'success': False, 'error': 'Vina not available'}

    ligand_files = get_ligand_files(ligand_path)
    if not ligand_files:
        return {'success': False, 'error': 'No ligand files found'}

    try:
        objs = cmd.get_names("objects")
    except AttributeError:
        objs = cmd.get_object_list() if hasattr(cmd, "get_object_list") else []

    if obj_name not in objs:
        return {'success': False, 'error': f"Object '{obj_name}' not found"}

    if output_dir is None:
        output_dir = f"{obj_name}_batch_docking"
    os.makedirs(output_dir, exist_ok=True)

    # 导出受体 PDBQT（只做一次）
    receptor_pdbqt = os.path.join(output_dir, f"{obj_name}_receptor.pdbqt")
    receptor_sel = f"({obj_name}) and not ({remove_selection})" if remove_selection else obj_name
    
    print(f"[vina_integration] 导出受体: {receptor_pdbqt}")
    if not export_to_pdbqt(obj_name, receptor_sel, receptor_pdbqt, is_receptor=True):
        return {'success': False, 'error': 'Receptor export failed'}

    results = []
    total = len(ligand_files)
    
    for i, lig_file in enumerate(ligand_files):
        ligand_name = os.path.splitext(os.path.basename(lig_file))[0]
        
        if progress_callback:
            progress_callback(i + 1, total, ligand_name)
        
        # 转换配体为 PDBQT
        ligand_pdbqt = os.path.join(output_dir, f"{ligand_name}.pdbqt")
        # 如果已经是 pdbqt 格式
        if lig_file.lower().endswith(".pdbqt"):
            # 检查源文件和目标文件是否相同，避免 shutil.copy 报错
            if os.path.abspath(lig_file) != os.path.abspath(ligand_pdbqt):
                shutil.copy(lig_file, ligand_pdbqt)
            # 如果相同，直接使用原文件路径
            else:
                ligand_pdbqt = lig_file
        elif not convert_ligand_to_pdbqt(lig_file, ligand_pdbqt):
            # 尝试通过 PyMOL 转换
            ligand_obj = f"temp_lig_{i}"
            try:
                cmd.load(lig_file, ligand_obj)
                if not export_to_pdbqt(ligand_obj, ligand_obj, ligand_pdbqt, is_receptor=False):
                    cmd.delete(ligand_obj)
                    results.append({'success': False, 'error': 'Conversion failed', 'ligand': ligand_name})
                    continue
                cmd.delete(ligand_obj)
            except Exception as e:
                results.append({'success': False, 'error': str(e), 'ligand': ligand_name})
                continue
        
        # 运行对接
        result = run_vina_docking(receptor_pdbqt, ligand_pdbqt, box_params, output_dir,
                                  ligand_name, exhaustiveness, num_modes)
        results.append(result)
    
    # 按亲和力排序
    successful = [r for r in results if r.get('success')]
    successful.sort(key=lambda x: x.get('affinity', 0))
    
    # 生成汇总 CSV
    csv_path = os.path.join(output_dir, "docking_results.csv")
    with open(csv_path, 'w') as f:
        f.write("Ligand,Affinity (kcal/mol),Output File\n")
        for r in successful:
            f.write(f"{r['ligand']},{r['affinity']:.2f},{r.get('output_pdbqt', '')}\n")
    
    return {
        'success': True,
        'output_dir': output_dir,
        'results': results,
        'csv_path': csv_path,
        'total': total,
        'successful': len(successful),
        'failed': total - len(successful)
    }


if __name__ == "__main__":
    print("[vina_integration] This is a PyMOL plugin module")
