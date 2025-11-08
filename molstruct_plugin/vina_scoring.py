# -*- coding: utf-8 -*-
"""
vina_scoring.py
MolStruct插件的AutoDock Vina评分集成

将Vina评分功能整合到PyMOL环境,提供更精确的结合能估算
- 自动检测Vina安装
- 支持快速评分模式(无重对接)
- 可与经验评分函数对比
"""

from __future__ import print_function
import os
import sys
import subprocess
import tempfile
import shutil
from pathlib import Path
from pymol import cmd

# ==================== Vina检测 ====================
def find_vina_executable():
    """
    自动查找Vina可执行文件
    
    返回:
        str: Vina可执行文件路径,未找到返回None
    """
    # 1. 检查PATH
    vina = shutil.which('vina') or shutil.which('vina.exe')
    if vina:
        return vina
    
    # 2. 检查常见安装位置(Windows)
    if sys.platform == 'win32':
        common_paths = [
            r'C:\Program Files\vina\vina.exe',
            r'C:\Program Files (x86)\vina\vina.exe',
            os.path.join(os.environ.get('LOCALAPPDATA', ''), 'Programs', 'VinaTools', 'vina', 'vina.exe'),
        ]
        for path in common_paths:
            if os.path.exists(path):
                return path
    
    # 3. 检查conda环境
    if 'CONDA_PREFIX' in os.environ:
        conda_vina = os.path.join(os.environ['CONDA_PREFIX'], 'bin', 'vina')
        if os.path.exists(conda_vina):
            return conda_vina
    
    return None


def check_vina_available():
    """检查Vina是否可用"""
    vina = find_vina_executable()
    if vina:
        print(f"[vina_scoring] ✅ Found Vina: {vina}")
        return True
    else:
        print("[vina_scoring] ⚠️ Vina not found in PATH")
        print("[vina_scoring] 💡 Install: conda install -c conda-forge vina")
        print("[vina_scoring]    Or download: https://github.com/ccsb-scripps/AutoDock-Vina/releases")
        return False


# ==================== PyMOL导出PDBQT ====================
def export_to_pdbqt_simple(obj_name, selection, output_pdbqt):
    """
    简化版PDBQT导出(不依赖MGLTools)
    
    注意: 这是简化版本,仅适用于快速评分
    对于严格对接,建议使用MGLTools的prepare_ligand4.py
    
    参数:
        obj_name: PyMOL对象名
        selection: PyMOL选择表达式
        output_pdbqt: 输出PDBQT文件路径
    """
    # 先导出PDB
    temp_pdb = output_pdbqt.replace('.pdbqt', '_temp.pdb')
    cmd.save(temp_pdb, selection)
    
    # 简单转换: PDB → PDBQT
    # 这里只是添加基本的AD4字段,不处理氢原子和电荷
    # 仅用于Vina的--score_only模式
    try:
        with open(temp_pdb, 'r') as f_in, open(output_pdbqt, 'w') as f_out:
            for line in f_in:
                if line.startswith(('ATOM', 'HETATM')):
                    # 添加简单的AD4类型(基于元素)
                    atom_name = line[12:16].strip()
                    element = line[76:78].strip() if len(line) > 77 else atom_name[0]
                    
                    # 简化的AD4类型映射
                    ad4_type = element[0]
                    if element in ['C', 'CA']:
                        ad4_type = 'C'
                    elif element in ['N', 'NA', 'NS']:
                        ad4_type = 'N'
                    elif element in ['O', 'OA', 'OS']:
                        ad4_type = 'O'
                    elif element in ['S', 'SA']:
                        ad4_type = 'S'
                    elif element in ['H', 'HD']:
                        ad4_type = 'H'
                    
                    # 写入PDBQT格式(添加电荷和AD4类型列)
                    # PDBQT格式: PDB格式 + charge + AD4_type
                    pdbqt_line = line.rstrip() + f"  0.00 {ad4_type:>2s}\n"
                    f_out.write(pdbqt_line)
                elif line.startswith(('MODEL', 'ENDMDL', 'TER', 'END')):
                    f_out.write(line)
        
        os.unlink(temp_pdb)
        return True
    
    except Exception as e:
        print(f"[export_to_pdbqt_simple] ⚠️ Export failed: {e}")
        if os.path.exists(temp_pdb):
            os.unlink(temp_pdb)
        return False


# ==================== Vina评分核心 ====================
def score_with_vina(protein_pdbqt, ligand_pdbqt, vina_bin=None, mode='score_only'):
    """
    使用Vina评分复合物
    
    参数:
        protein_pdbqt: 蛋白PDBQT文件路径
        ligand_pdbqt: 配体PDBQT文件路径
        vina_bin: Vina可执行文件路径(None则自动查找)
        mode: 'score_only' 或 'local_only'
            - score_only: 仅评分,不优化(最快,<1秒)
            - local_only: 局部优化后评分(较慢,约5-10秒)
    
    返回:
        dict: {
            'affinity': float,  # kcal/mol
            'success': bool,
            'output': str
        }
    """
    if vina_bin is None:
        vina_bin = find_vina_executable()
    
    if not vina_bin:
        return {'success': False, 'error': 'Vina not found', 'affinity': None}
    
    # 构建Vina命令
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
        
        # 解析评分
        affinity = None
        for line in output.split('\n'):
            if 'Affinity:' in line or 'kcal/mol' in line:
                # 查找数字
                parts = line.split()
                for i, part in enumerate(parts):
                    try:
                        val = float(part)
                        affinity = val
                        break
                    except ValueError:
                        continue
                if affinity is not None:
                    break
        
        if affinity is None:
            return {
                'success': False,
                'error': 'Could not parse affinity',
                'output': output,
                'affinity': None
            }
        
        return {
            'success': True,
            'affinity': affinity,
            'output': output
        }
    
    except subprocess.TimeoutExpired:
        return {'success': False, 'error': 'Timeout', 'affinity': None}
    except Exception as e:
        return {'success': False, 'error': str(e), 'affinity': None}


# ==================== PyMOL命令封装 ====================
def vina_score_complex(protein_obj, ligand_selection, mode='score_only', show_report=True):
    """
    PyMOL命令: 使用Vina评分蛋白-配体复合物
    
    用法:
        vina_score_complex protein, resn LIG
        vina_score_complex protein, resn LIG, mode=local_only
    
    参数:
        protein_obj: 蛋白PyMOL对象名
        ligand_selection: 配体选择表达式(如 "resn LIG" 或 "chain L and resn MK1")
        mode: 'score_only' 或 'local_only'
        show_report: 是否打印报告
    
    返回:
        dict: Vina评分结果
    """
    # 检查Vina
    if not check_vina_available():
        print("[vina_score_complex] ⚠️ Vina not available, cannot proceed")
        return None
    
    # 检查对象
    if protein_obj not in cmd.get_object_list():
        print(f"[vina_score_complex] ⚠️ Object '{protein_obj}' not found")
        return None
    
    # 创建临时文件
    with tempfile.TemporaryDirectory() as tmpdir:
        protein_pdbqt = os.path.join(tmpdir, 'protein.pdbqt')
        ligand_pdbqt = os.path.join(tmpdir, 'ligand.pdbqt')
        
        # 导出蛋白
        print(f"[vina_score_complex] Exporting protein: {protein_obj}")
        if not export_to_pdbqt_simple(protein_obj, protein_obj, protein_pdbqt):
            print("[vina_score_complex] ⚠️ Failed to export protein")
            return None
        
        # 导出配体
        print(f"[vina_score_complex] Exporting ligand: {ligand_selection}")
        if not export_to_pdbqt_simple(protein_obj, ligand_selection, ligand_pdbqt):
            print("[vina_score_complex] ⚠️ Failed to export ligand")
            return None
        
        # Vina评分
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
        print("⚠️  Note: This is Vina's empirical scoring function")
        print("    For publication, use full flexible docking workflow")
        print("=" * 60)
    elif not result['success']:
        print(f"[vina_score_complex] ⚠️ Vina scoring failed: {result.get('error', 'Unknown error')}")
    
    return result


def compare_scoring_methods(protein_obj, ligand_resname):
    """
    PyMOL命令: 对比经验评分 vs Vina评分
    
    用法:
        compare_scoring_methods protein, LIG
    """
    print("=" * 60)
    print("Scoring Method Comparison")
    print("=" * 60)
    
    # 1. 经验评分(MolStruct内置)
    try:
        from .binding_score import calculate_binary_score
        from .interaction_analyzer import analyze_protein_ligand_interactions
        
        print("\n[1/2] Running MolStruct empirical scoring...")
        result = analyze_protein_ligand_interactions(protein_obj, ligand_resname)
        empirical_score = calculate_binary_score(result)
        
        if empirical_score:
            print(f"✅ MolStruct Score: {empirical_score['total']:.2f} kcal/mol")
            print(f"   Components: H-bond={empirical_score['components']['hbond']:.2f}, "
                  f"Ionic={empirical_score['components']['ionic']:.2f}, "
                  f"Hydrophobic={empirical_score['components']['hydrophobic']:.2f}")
        else:
            print("⚠️  MolStruct scoring failed")
    except ImportError:
        print("⚠️  MolStruct scoring module not available")
        empirical_score = None
    
    # 2. Vina评分
    print("\n[2/2] Running Vina scoring...")
    vina_result = vina_score_complex(
        protein_obj,
        f"resn {ligand_resname}",
        mode='score_only',
        show_report=False
    )
    
    if vina_result and vina_result['success']:
        vina_score = vina_result['affinity']
        print(f"✅ Vina Score: {vina_score:.2f} kcal/mol")
    else:
        print("⚠️  Vina scoring failed")
        vina_score = None
    
    # 对比总结
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    
    if empirical_score and vina_score:
        print(f"MolStruct Empirical: {empirical_score['total']:>7.2f} kcal/mol  (Fast, ~0.1s)")
        print(f"Vina Score:          {vina_score:>7.2f} kcal/mol  (Medium, ~1s)")
        print(f"\nDifference:          {abs(empirical_score['total'] - vina_score):>7.2f} kcal/mol")
        
        print("\n💡 Recommendations:")
        print("   • Use MolStruct for quick screening (>1000 compounds)")
        print("   • Use Vina for validation (top 10-50 hits)")
        print("   • Use full Vina docking for final candidates")
    
    print("=" * 60)
    
    return {
        'empirical': empirical_score,
        'vina': vina_result
    }


# ==================== 批量评分(与现有Vina项目整合) ====================
def batch_score_vina_results(results_dir):
    """
    批量评分Vina对接结果,添加经验评分作为补充
    
    用法(在PyMOL外):
        python -c "from molstruct_plugin.vina_scoring import batch_score_vina_results; batch_score_vina_results('./results')"
    
    参数:
        results_dir: Vina对接结果目录(包含*_out子目录)
    """
    import glob
    import csv
    from pymol import cmd
    
    print(f"[batch_score_vina_results] Scanning: {results_dir}")
    
    # 查找所有*_out目录
    out_dirs = glob.glob(os.path.join(results_dir, '*_out'))
    
    if not out_dirs:
        print("[batch_score_vina_results] ⚠️ No *_out directories found")
        return
    
    print(f"[batch_score_vina_results] Found {len(out_dirs)} result directories")
    
    # 尝试导入评分模块
    try:
        from .binding_score import calculate_binary_score
        from .interaction_analyzer import analyze_protein_ligand_interactions
        has_empirical = True
    except ImportError:
        print("[batch_score_vina_results] ⚠️ Empirical scoring not available")
        has_empirical = False
    
    for out_dir in out_dirs:
        receptor_name = os.path.basename(out_dir).replace('_out', '')
        print(f"\n[batch_score_vina_results] Processing: {receptor_name}")
        
        # 查找所有PDBQT结果
        pdbqt_files = glob.glob(os.path.join(out_dir, '*.pdbqt'))
        
        if not pdbqt_files:
            print(f"[batch_score_vina_results]   ⚠️ No PDBQT files in {out_dir}")
            continue
        
        # 创建增强评分CSV
        enhanced_csv = os.path.join(results_dir, f"{receptor_name}_scores_enhanced.csv")
        
        with open(enhanced_csv, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Ligand', 'Vina_Affinity', 'Empirical_Score', 
                'H-bonds', 'Salt_Bridges', 'Hydrophobic_Contacts'
            ])
            
            for pdbqt_file in pdbqt_files:
                lig_name = os.path.splitext(os.path.basename(pdbqt_file))[0]
                
                # 读取Vina评分
                vina_affinity = 'NA'
                try:
                    with open(pdbqt_file, 'r') as pf:
                        for line in pf:
                            if line.startswith('REMARK VINA RESULT'):
                                parts = line.split()
                                if len(parts) >= 4:
                                    vina_affinity = parts[3]
                                    break
                except:
                    pass
                
                # 计算经验评分
                empirical_score = 'NA'
                hbonds = ionic = hydrophobic = 0
                
                if has_empirical:
                    try:
                        # 加载到PyMOL
                        obj_name = f"temp_{lig_name}"
                        cmd.load(pdbqt_file, obj_name)
                        
                        # 分析相互作用
                        result = analyze_protein_ligand_interactions(obj_name, lig_name)
                        if result:
                            score_dict = calculate_binary_score(result)
                            if score_dict:
                                empirical_score = f"{score_dict['total']:.2f}"
                                hbonds = score_dict['interaction_counts'].get('hbond', 0)
                                ionic = score_dict['interaction_counts'].get('ionic', 0)
                                hydrophobic = score_dict['interaction_counts'].get('hydrophobic', 0)
                        
                        cmd.delete(obj_name)
                    except Exception as e:
                        print(f"      ⚠️ {lig_name}: {e}")
                
                writer.writerow([
                    lig_name, vina_affinity, empirical_score,
                    hbonds, ionic, hydrophobic
                ])
        
        print(f"[batch_score_vina_results]   ✅ Enhanced CSV: {enhanced_csv}")


if __name__ == "__main__":
    print("[vina_scoring] This is a PyMOL plugin module")
    print("[vina_scoring] Commands:")
    print("  - vina_score_complex")
    print("  - compare_scoring_methods")
