# -*- coding: utf-8 -*-
"""
binding_heatmap.py
GLINTPlugin的批量结合能热图生成Module

功能:
- 自动扫描多个 *_scores.csv File
- 合并数据生成受体-配体结合能矩阵
- 生成蓝色渐变热图
"""

from __future__ import print_function
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # 非交互式后端，避免GUI冲突
import matplotlib.pyplot as plt
import seaborn as sns


def load_binding_csv(csv_path):
    """
    Load单个结合能CSVFile
    
    自动识别列名（支持多种格式）:
    - Ligand/Compound/Name → 配体Name
    - BindingEnergy/Score/Energy/Affinity → 结合能
    
    Parameters:
        csv_path: CSVFilePath
    
    Return:
        DataFrame: Package含 Ligand, BindingEnergy 两列
    """
    df = pd.read_csv(csv_path)
    
    # 自动识别列名（不区分Size写）
    cols_lower = [c.lower() for c in df.columns]
    lig_col, val_col = None, None
    
    for i, c in enumerate(cols_lower):
        if 'ligand' in c or 'compound' in c or 'name' in c:
            lig_col = df.columns[i]
        if 'binding' in c or 'score' in c or 'energy' in c or 'affinity' in c:
            val_col = df.columns[i]
    
    # 如果没找到，using前两列
    if lig_col is None:
        lig_col = df.columns[0]
    if val_col is None:
        val_candidates = [c for c in df.columns if c != lig_col]
        val_col = val_candidates[0] if val_candidates else df.columns[1]
    
    # 提取并重命名
    out = df[[lig_col, val_col]].copy()
    out.columns = ["Ligand", "BindingEnergy"]
    
    # 清洗数据
    out["Ligand"] = out["Ligand"].astype(str).str.strip()
    out["BindingEnergy"] = pd.to_numeric(out["BindingEnergy"], errors="coerce")
    
    return out.dropna(subset=["Ligand", "BindingEnergy"])


def generate_binding_heatmap(input_folder, output_path=None, pattern="*_scores.csv"):
    """
    批量生成结合能热图
    
    Parameters:
        input_folder: Package含多个 CSV File的Directory
        output_path: 输出PNGPath（可选，默认自动生成）
        pattern: File匹配模式（默认 *_scores.csv）
    
    Return:
        dict: {
            'success': bool,
            'output_path': str,
            'n_receptors': int,
            'n_ligands': int,
            'error': str (if failed)
        }
    """
    try:
        # 1. 扫描CSVFile
        import glob
        files = glob.glob(os.path.join(input_folder, pattern))
        
        if not files:
            return {
                'success': False,
                'error': f'No files matching "{pattern}" found in {input_folder}'
            }
        
        print(f"[binding_heatmap] Found {len(files)} score files:")
        for f in files:
            print(f"  - {os.path.basename(f)}")
        
        # 2. Load并合并所有CSV
        merged = None
        for fpath in files:
            # 受体Name：从File名提取（去掉 _scores.csv）
            fname = os.path.basename(fpath)
            receptor = os.path.splitext(fname)[0].replace("_scores", "")
            
            df = load_binding_csv(fpath)
            df = df.rename(columns={"BindingEnergy": receptor})
            
            if merged is None:
                merged = df
            else:
                merged = pd.merge(merged, df, on="Ligand", how="outer")
        
        # 3. 准备矩阵
        matrix = merged.set_index("Ligand").sort_index()
        
        # 按平均结合能Sort（最强结合在上）
        matrix = matrix.reindex(matrix.mean(axis=1).sort_values(ascending=True).index)
        
        # 4. 生成热图
        plt.figure(figsize=(max(11, len(matrix.columns) * 1.2), max(6, len(matrix) * 0.4)), dpi=180)
        
        ax = sns.heatmap(
            matrix,
            cmap="Blues",              # 蓝色渐变
            annot=True,                # Display数Value
            fmt=".2f",                 # 保留2位小数
            linewidths=1.0,
            linecolor="white",
            cbar_kws={"label": "Binding Energy (kcal/mol)"},
            vmin=-10, vmax=-4,         # 动态范围
            center=-7,                 # 中心点
            square=False
        )
        
        ax.set_xlabel("Receptor (Gene/Protein)", fontsize=12, fontweight='bold')
        ax.set_ylabel("Ligand", fontsize=12, fontweight='bold')
        plt.xticks(rotation=-30, ha="left")
        plt.yticks(rotation=0)
        
        plt.title(
            "Receptor–Ligand Binding Energy Heatmap\n"
            "(Darker blue = lower binding energy → stronger affinity)",
            fontsize=14, pad=16, fontweight='bold'
        )
        
        plt.tight_layout()
        
        # 5. Save
        if output_path is None:
            output_path = os.path.join(input_folder, "binding_energy_heatmap.png")
        
        plt.savefig(output_path, bbox_inches="tight", dpi=300)
        plt.close()
        
        print(f"[binding_heatmap] ✅ Heatmap saved to: {output_path}")
        
        return {
            'success': True,
            'output_path': output_path,
            'n_receptors': len(matrix.columns),
            'n_ligands': len(matrix),
            'matrix': matrix
        }
    
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {
            'success': False,
            'error': str(e)
        }


def generate_heatmap_from_files(csv_files, output_path=None):
    """
    从指定的CSVFile列表生成热图
    
    Parameters:
        csv_files: CSVFilePath列表
        output_path: 输出PNGPath（可选）
    
    Return:
        dict: 同 generate_binding_heatmap
    """
    try:
        if not csv_files:
            return {'success': False, 'error': 'No CSV files provided'}
        
        print(f"[binding_heatmap] Loading {len(csv_files)} score files...")
        
        # Load并合并
        merged = None
        for fpath in csv_files:
            fname = os.path.basename(fpath)
            receptor = os.path.splitext(fname)[0].replace("_scores", "")
            
            df = load_binding_csv(fpath)
            df = df.rename(columns={"BindingEnergy": receptor})
            
            if merged is None:
                merged = df
            else:
                merged = pd.merge(merged, df, on="Ligand", how="outer")
        
        # 准备矩阵
        matrix = merged.set_index("Ligand").sort_index()
        matrix = matrix.reindex(matrix.mean(axis=1).sort_values(ascending=True).index)
        
        # 生成热图
        plt.figure(figsize=(max(11, len(matrix.columns) * 1.2), max(6, len(matrix) * 0.4)), dpi=180)
        
        ax = sns.heatmap(
            matrix,
            cmap="Blues",
            annot=True,
            fmt=".2f",
            linewidths=1.0,
            linecolor="white",
            cbar_kws={"label": "Binding Energy (kcal/mol)"},
            vmin=-10, vmax=-4,
            center=-7,
            square=False
        )
        
        ax.set_xlabel("Receptor (Gene/Protein)", fontsize=12, fontweight='bold')
        ax.set_ylabel("Ligand", fontsize=12, fontweight='bold')
        plt.xticks(rotation=-30, ha="left")
        plt.yticks(rotation=0)
        
        plt.title(
            "Receptor–Ligand Binding Energy Heatmap\n"
            "(Darker blue = lower binding energy → stronger affinity)",
            fontsize=14, pad=16, fontweight='bold'
        )
        
        plt.tight_layout()
        
        # Save
        if output_path is None:
            output_dir = os.path.dirname(csv_files[0]) if csv_files else os.getcwd()
            output_path = os.path.join(output_dir, "binding_energy_heatmap.png")
        
        plt.savefig(output_path, bbox_inches="tight", dpi=300)
        plt.close()
        
        print(f"[binding_heatmap] ✅ Heatmap saved to: {output_path}")
        
        return {
            'success': True,
            'output_path': output_path,
            'n_receptors': len(matrix.columns),
            'n_ligands': len(matrix),
            'matrix': matrix
        }
    
    except Exception as e:
        import traceback
        traceback.print_exc()
        return {
            'success': False,
            'error': str(e)
        }


# PyMOL命令封装
def plot_binding_heatmap(input_folder, output_path=None):
    """
    PyMOL命令: 生成批量结合能热图
    
    用法:
        plot_binding_heatmap /path/to/csv_folder
        plot_binding_heatmap /path/to/csv_folder, /path/to/output.png
    """
    result = generate_binding_heatmap(input_folder, output_path)
    
    if result['success']:
        print("=" * 60)
        print("Binding Energy Heatmap Generated")
        print("=" * 60)
        print(f"Receptors: {result['n_receptors']}")
        print(f"Ligands:   {result['n_ligands']}")
        print(f"Output:    {result['output_path']}")
        print("=" * 60)
    else:
        print(f"[plot_binding_heatmap] ⚠️  Failed: {result.get('error', 'Unknown error')}")
    
    return result
