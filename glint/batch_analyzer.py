# -*- coding: utf-8 -*-
"""
GLINT Batch Analyzer
=====================
批量 PDB 分析模块 - 一次性分析多个结构

功能：
1. 批量 G-motif 检测
2. 批量 PPI 界面分析
3. 批量口袋检测
4. 批量相互作用分析
5. 结果汇总与导出

Author: Roufen Chen
Date: 2025-12
"""

from __future__ import print_function
import os
import csv
import json
import tempfile
from datetime import datetime
from typing import List, Dict, Any, Optional, Callable
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import traceback

try:
    from pymol import cmd
except ImportError:
    cmd = None


class BatchAnalyzer:
    """批量分析器"""
    
    def __init__(self, output_dir: Optional[str] = None, max_workers: int = 4):
        """
        初始化批量分析器
        
        参数:
            output_dir: 输出目录（默认当前目录下的 batch_results）
            max_workers: 最大并行工作线程数
        """
        self.output_dir = Path(output_dir) if output_dir else Path("batch_results")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers
        self.results = []
        self.errors = []
        
    def _load_pdb(self, pdb_path: str, obj_name: Optional[str] = None) -> Optional[str]:
        """加载 PDB 文件到 PyMOL"""
        if not cmd:
            print("❌ PyMOL 不可用")
            return None
            
        if obj_name is None:
            obj_name = Path(pdb_path).stem
            
        # 确保对象名唯一
        base_name = obj_name
        counter = 1
        while obj_name in cmd.get_names("objects"):
            obj_name = f"{base_name}_{counter}"
            counter += 1
            
        try:
            cmd.load(pdb_path, obj_name)
            return obj_name
        except Exception as e:
            print(f"❌ 加载失败 {pdb_path}: {e}")
            return None
    
    def _fetch_pdb(self, pdb_id: str) -> Optional[str]:
        """从 PDB 数据库获取结构"""
        if not cmd:
            print("❌ PyMOL 不可用")
            return None
            
        obj_name = pdb_id.lower()
        try:
            cmd.fetch(pdb_id, obj_name)
            return obj_name
        except Exception as e:
            print(f"❌ 获取失败 {pdb_id}: {e}")
            return None
    
    def batch_gmotif_detection(self, 
                               pdb_sources: List[str],
                               rmsd_cutoff: float = 3.5,
                               require_gly: bool = True,
                               template: str = "GSPT1",
                               output_csv: Optional[str] = None) -> Dict[str, Any]:
        """
        批量 G-motif 检测
        
        参数:
            pdb_sources: PDB 文件路径或 PDB ID 列表
            rmsd_cutoff: RMSD 阈值
            require_gly: 是否要求中心 Gly
            template: 模板名称 (GSPT1, CK1α, VAV1)
            output_csv: 输出 CSV 路径
            
        返回:
            汇总结果字典
        """
        print(f"\n{'='*60}")
        print(f"🔬 批量 G-motif 检测")
        print(f"{'='*60}")
        print(f"输入: {len(pdb_sources)} 个结构")
        print(f"参数: RMSD≤{rmsd_cutoff}Å, 模板={template}")
        print(f"{'='*60}\n")
        
        try:
            from .g_motif_analyzer import find_crbn_g_motif
        except ImportError:
            try:
                from g_motif_analyzer import find_crbn_g_motif
            except ImportError:
                print("❌ G-motif 分析模块不可用")
                return {"success": False, "error": "Module not available"}
        
        results = []
        total_hits = 0
        
        for i, source in enumerate(pdb_sources, 1):
            print(f"[{i}/{len(pdb_sources)}] 分析: {source}")
            
            # 判断是文件路径还是 PDB ID
            if os.path.exists(source):
                obj_name = self._load_pdb(source)
            elif len(source) == 4 and source.isalnum():
                obj_name = self._fetch_pdb(source)
            else:
                print(f"  ⚠️ 无效输入: {source}")
                results.append({
                    "source": source,
                    "success": False,
                    "error": "Invalid input",
                    "hits": []
                })
                continue
            
            if not obj_name:
                results.append({
                    "source": source,
                    "success": False,
                    "error": "Load failed",
                    "hits": []
                })
                continue
            
            try:
                # 运行 G-motif 检测
                hits = find_crbn_g_motif(
                    obj_name=obj_name,
                    rmsd_cutoff=rmsd_cutoff,
                    require_gly=require_gly,
                    template_name=template
                )
                
                hit_count = len(hits) if hits else 0
                total_hits += hit_count
                
                results.append({
                    "source": source,
                    "obj_name": obj_name,
                    "success": True,
                    "hit_count": hit_count,
                    "hits": hits or []
                })
                
                print(f"  ✅ 发现 {hit_count} 个 G-motif")
                
            except Exception as e:
                print(f"  ❌ 分析失败: {e}")
                results.append({
                    "source": source,
                    "obj_name": obj_name,
                    "success": False,
                    "error": str(e),
                    "hits": []
                })
            
            # 清理对象（可选）
            # cmd.delete(obj_name)
        
        # 汇总
        summary = {
            "analysis_type": "G-motif Detection",
            "timestamp": datetime.now().isoformat(),
            "parameters": {
                "rmsd_cutoff": rmsd_cutoff,
                "require_gly": require_gly,
                "template": template
            },
            "total_structures": len(pdb_sources),
            "successful": sum(1 for r in results if r["success"]),
            "total_hits": total_hits,
            "results": results
        }
        
        # 导出 CSV
        if output_csv:
            self._export_gmotif_csv(results, output_csv)
        
        # 打印汇总
        print(f"\n{'='*60}")
        print(f"📊 G-motif 检测汇总")
        print(f"{'='*60}")
        print(f"总结构数: {len(pdb_sources)}")
        print(f"成功分析: {summary['successful']}")
        print(f"总 G-motif 数: {total_hits}")
        print(f"{'='*60}\n")
        
        return summary
    
    def batch_ppi_analysis(self,
                          pdb_sources: List[str],
                          chain_pairs: List[tuple],
                          interface_distance: float = 4.5,
                          output_csv: Optional[str] = None) -> Dict[str, Any]:
        """
        批量 PPI 界面分析
        
        参数:
            pdb_sources: PDB 文件路径或 PDB ID 列表
            chain_pairs: 链对列表 [(['A'], ['B']), (['C'], ['D']), ...]
                        如果只有一个元组，则应用于所有结构
            interface_distance: 界面距离阈值
            output_csv: 输出 CSV 路径
            
        返回:
            汇总结果字典
        """
        print(f"\n{'='*60}")
        print(f"🔬 批量 PPI 界面分析")
        print(f"{'='*60}")
        print(f"输入: {len(pdb_sources)} 个结构")
        print(f"参数: 界面距离≤{interface_distance}Å")
        print(f"{'='*60}\n")
        
        try:
            from .ppi_analyzer import analyze_protein_protein_interface
        except ImportError:
            try:
                from ppi_analyzer import analyze_protein_protein_interface
            except ImportError:
                print("❌ PPI 分析模块不可用")
                return {"success": False, "error": "Module not available"}
        
        results = []
        
        # 如果只提供一个链对，应用于所有结构
        if len(chain_pairs) == 1:
            chain_pairs = chain_pairs * len(pdb_sources)
        elif len(chain_pairs) != len(pdb_sources):
            print("⚠️ 链对数量与结构数量不匹配，使用第一个链对")
            chain_pairs = [chain_pairs[0]] * len(pdb_sources)
        
        for i, (source, (chains1, chains2)) in enumerate(zip(pdb_sources, chain_pairs), 1):
            print(f"[{i}/{len(pdb_sources)}] 分析: {source}")
            print(f"  链对: {chains1} vs {chains2}")
            
            # 加载结构
            if os.path.exists(source):
                obj_name = self._load_pdb(source)
            elif len(source) == 4 and source.isalnum():
                obj_name = self._fetch_pdb(source)
            else:
                print(f"  ⚠️ 无效输入: {source}")
                results.append({
                    "source": source,
                    "success": False,
                    "error": "Invalid input"
                })
                continue
            
            if not obj_name:
                results.append({
                    "source": source,
                    "success": False,
                    "error": "Load failed"
                })
                continue
            
            try:
                # 运行 PPI 分析
                ppi_result = analyze_protein_protein_interface(
                    obj_name=obj_name,
                    protein1_chains=chains1,
                    protein2_chains=chains2,
                    interface_distance=interface_distance,
                    visualize=False  # 批量模式不可视化
                )
                
                if ppi_result:
                    results.append({
                        "source": source,
                        "obj_name": obj_name,
                        "success": True,
                        "chains1": chains1,
                        "chains2": chains2,
                        "interface_contacts": ppi_result.get("interface_contacts", 0),
                        "interface_strength": ppi_result.get("interface_strength", 0),
                        "bsa": ppi_result.get("bsa"),
                        "is_strong_interface": ppi_result.get("is_strong_interface", False)
                    })
                    print(f"  ✅ 接触数: {ppi_result.get('interface_contacts', 0)}, "
                          f"强度: {ppi_result.get('interface_strength', 0):.1f}")
                else:
                    results.append({
                        "source": source,
                        "obj_name": obj_name,
                        "success": False,
                        "error": "No result"
                    })
                    
            except Exception as e:
                print(f"  ❌ 分析失败: {e}")
                results.append({
                    "source": source,
                    "obj_name": obj_name,
                    "success": False,
                    "error": str(e)
                })
        
        # 汇总
        successful_results = [r for r in results if r["success"]]
        summary = {
            "analysis_type": "PPI Interface Analysis",
            "timestamp": datetime.now().isoformat(),
            "parameters": {
                "interface_distance": interface_distance
            },
            "total_structures": len(pdb_sources),
            "successful": len(successful_results),
            "avg_contacts": sum(r.get("interface_contacts", 0) for r in successful_results) / len(successful_results) if successful_results else 0,
            "avg_strength": sum(r.get("interface_strength", 0) for r in successful_results) / len(successful_results) if successful_results else 0,
            "strong_interfaces": sum(1 for r in successful_results if r.get("is_strong_interface")),
            "results": results
        }
        
        # 导出 CSV
        if output_csv:
            self._export_ppi_csv(results, output_csv)
        
        # 打印汇总
        print(f"\n{'='*60}")
        print(f"📊 PPI 分析汇总")
        print(f"{'='*60}")
        print(f"总结构数: {len(pdb_sources)}")
        print(f"成功分析: {summary['successful']}")
        print(f"平均接触数: {summary['avg_contacts']:.1f}")
        print(f"平均强度: {summary['avg_strength']:.1f}")
        print(f"强界面数: {summary['strong_interfaces']}")
        print(f"{'='*60}\n")
        
        return summary
    
    def batch_pocket_detection(self,
                              pdb_sources: List[str],
                              min_volume: float = 30.0,
                              min_depth: float = 2.5,
                              output_csv: Optional[str] = None) -> Dict[str, Any]:
        """
        批量口袋检测
        
        参数:
            pdb_sources: PDB 文件路径或 PDB ID 列表
            min_volume: 最小口袋体积
            min_depth: 最小口袋深度
            output_csv: 输出 CSV 路径
            
        返回:
            汇总结果字典
        """
        print(f"\n{'='*60}")
        print(f"🔬 批量口袋检测")
        print(f"{'='*60}")
        print(f"输入: {len(pdb_sources)} 个结构")
        print(f"参数: 体积≥{min_volume}ų, 深度≥{min_depth}Å")
        print(f"{'='*60}\n")
        
        try:
            from .pocket_detector import detect_pockets
        except ImportError:
            try:
                from pocket_detector import detect_pockets
            except ImportError:
                print("❌ 口袋检测模块不可用")
                return {"success": False, "error": "Module not available"}
        
        results = []
        total_pockets = 0
        
        for i, source in enumerate(pdb_sources, 1):
            print(f"[{i}/{len(pdb_sources)}] 分析: {source}")
            
            # 加载结构
            if os.path.exists(source):
                obj_name = self._load_pdb(source)
            elif len(source) == 4 and source.isalnum():
                obj_name = self._fetch_pdb(source)
            else:
                print(f"  ⚠️ 无效输入: {source}")
                results.append({
                    "source": source,
                    "success": False,
                    "error": "Invalid input",
                    "pockets": []
                })
                continue
            
            if not obj_name:
                results.append({
                    "source": source,
                    "success": False,
                    "error": "Load failed",
                    "pockets": []
                })
                continue
            
            try:
                # 运行口袋检测
                pockets = detect_pockets(
                    obj_name=obj_name,
                    min_volume=min_volume,
                    min_depth=min_depth
                )
                
                pocket_count = len(pockets) if pockets else 0
                total_pockets += pocket_count
                
                # 简化口袋信息
                pocket_summary = []
                for p in (pockets or []):
                    pocket_summary.append({
                        "id": p.get("id"),
                        "volume": p.get("volume"),
                        "depth": p.get("depth"),
                        "druggability": p.get("druggability_score")
                    })
                
                results.append({
                    "source": source,
                    "obj_name": obj_name,
                    "success": True,
                    "pocket_count": pocket_count,
                    "pockets": pocket_summary
                })
                
                print(f"  ✅ 发现 {pocket_count} 个口袋")
                
            except Exception as e:
                print(f"  ❌ 分析失败: {e}")
                results.append({
                    "source": source,
                    "obj_name": obj_name,
                    "success": False,
                    "error": str(e),
                    "pockets": []
                })
        
        # 汇总
        summary = {
            "analysis_type": "Pocket Detection",
            "timestamp": datetime.now().isoformat(),
            "parameters": {
                "min_volume": min_volume,
                "min_depth": min_depth
            },
            "total_structures": len(pdb_sources),
            "successful": sum(1 for r in results if r["success"]),
            "total_pockets": total_pockets,
            "results": results
        }
        
        # 导出 CSV
        if output_csv:
            self._export_pocket_csv(results, output_csv)
        
        # 打印汇总
        print(f"\n{'='*60}")
        print(f"📊 口袋检测汇总")
        print(f"{'='*60}")
        print(f"总结构数: {len(pdb_sources)}")
        print(f"成功分析: {summary['successful']}")
        print(f"总口袋数: {total_pockets}")
        print(f"{'='*60}\n")
        
        return summary
    
    def batch_interaction_analysis(self,
                                   pdb_sources: List[str],
                                   ligand_resnames: Optional[List[str]] = None,
                                   distance_cutoff: float = 4.5,
                                   output_csv: Optional[str] = None) -> Dict[str, Any]:
        """
        批量蛋白-配体相互作用分析
        
        参数:
            pdb_sources: PDB 文件路径或 PDB ID 列表
            ligand_resnames: 配体残基名列表（None 则自动检测）
            distance_cutoff: 距离阈值
            output_csv: 输出 CSV 路径
            
        返回:
            汇总结果字典
        """
        print(f"\n{'='*60}")
        print(f"🔬 批量蛋白-配体相互作用分析")
        print(f"{'='*60}")
        print(f"输入: {len(pdb_sources)} 个结构")
        print(f"参数: 距离≤{distance_cutoff}Å")
        print(f"{'='*60}\n")
        
        try:
            from .interaction_analyzer import analyze_protein_ligand_interactions
        except ImportError:
            try:
                from interaction_analyzer import analyze_protein_ligand_interactions
            except ImportError:
                print("❌ 相互作用分析模块不可用")
                return {"success": False, "error": "Module not available"}
        
        results = []
        
        # 处理配体名称列表
        if ligand_resnames is None:
            ligand_resnames = [None] * len(pdb_sources)
        elif len(ligand_resnames) == 1:
            ligand_resnames = ligand_resnames * len(pdb_sources)
        elif len(ligand_resnames) != len(pdb_sources):
            print("⚠️ 配体名称数量与结构数量不匹配，使用自动检测")
            ligand_resnames = [None] * len(pdb_sources)
        
        for i, (source, lig_name) in enumerate(zip(pdb_sources, ligand_resnames), 1):
            print(f"[{i}/{len(pdb_sources)}] 分析: {source}")
            if lig_name:
                print(f"  配体: {lig_name}")
            
            # 加载结构
            if os.path.exists(source):
                obj_name = self._load_pdb(source)
            elif len(source) == 4 and source.isalnum():
                obj_name = self._fetch_pdb(source)
            else:
                print(f"  ⚠️ 无效输入: {source}")
                results.append({
                    "source": source,
                    "success": False,
                    "error": "Invalid input"
                })
                continue
            
            if not obj_name:
                results.append({
                    "source": source,
                    "success": False,
                    "error": "Load failed"
                })
                continue
            
            try:
                # 运行相互作用分析
                ia_result = analyze_protein_ligand_interactions(
                    obj_name=obj_name,
                    ligand_resname=lig_name,
                    distance_cutoff=distance_cutoff
                )
                
                if ia_result:
                    interactions = ia_result.get("interactions", [])
                    
                    # 统计相互作用类型
                    type_counts = {}
                    for inter in interactions:
                        itype = inter.get("type", "Unknown")
                        type_counts[itype] = type_counts.get(itype, 0) + 1
                    
                    results.append({
                        "source": source,
                        "obj_name": obj_name,
                        "success": True,
                        "ligand": ia_result.get("ligand_resname", lig_name),
                        "interaction_count": len(interactions),
                        "type_counts": type_counts
                    })
                    
                    print(f"  ✅ 发现 {len(interactions)} 个相互作用")
                else:
                    results.append({
                        "source": source,
                        "obj_name": obj_name,
                        "success": False,
                        "error": "No result"
                    })
                    
            except Exception as e:
                print(f"  ❌ 分析失败: {e}")
                results.append({
                    "source": source,
                    "obj_name": obj_name,
                    "success": False,
                    "error": str(e)
                })
        
        # 汇总
        successful_results = [r for r in results if r["success"]]
        summary = {
            "analysis_type": "Protein-Ligand Interaction Analysis",
            "timestamp": datetime.now().isoformat(),
            "parameters": {
                "distance_cutoff": distance_cutoff
            },
            "total_structures": len(pdb_sources),
            "successful": len(successful_results),
            "avg_interactions": sum(r.get("interaction_count", 0) for r in successful_results) / len(successful_results) if successful_results else 0,
            "results": results
        }
        
        # 导出 CSV
        if output_csv:
            self._export_interaction_csv(results, output_csv)
        
        # 打印汇总
        print(f"\n{'='*60}")
        print(f"📊 相互作用分析汇总")
        print(f"{'='*60}")
        print(f"总结构数: {len(pdb_sources)}")
        print(f"成功分析: {summary['successful']}")
        print(f"平均相互作用数: {summary['avg_interactions']:.1f}")
        print(f"{'='*60}\n")
        
        return summary
    
    def comprehensive_batch_analysis(self,
                                     pdb_sources: List[str],
                                     chain_pairs: Optional[List[tuple]] = None,
                                     ligand_resnames: Optional[List[str]] = None,
                                     output_dir: Optional[str] = None) -> Dict[str, Any]:
        """
        综合批量分析 - 运行所有分析类型
        
        参数:
            pdb_sources: PDB 文件路径或 PDB ID 列表
            chain_pairs: PPI 分析的链对列表
            ligand_resnames: 配体名称列表
            output_dir: 输出目录
            
        返回:
            所有分析结果的汇总
        """
        if output_dir:
            self.output_dir = Path(output_dir)
            self.output_dir.mkdir(parents=True, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        print(f"\n{'🔬'*30}")
        print(f"GLINT 综合批量分析")
        print(f"{'🔬'*30}")
        print(f"输入: {len(pdb_sources)} 个结构")
        print(f"输出目录: {self.output_dir}")
        print(f"{'🔬'*30}\n")
        
        all_results = {
            "timestamp": timestamp,
            "input_count": len(pdb_sources),
            "analyses": {}
        }
        
        # 1. G-motif 检测
        print("\n" + "="*60)
        print("📍 步骤 1/4: G-motif 检测")
        print("="*60)
        gmotif_csv = self.output_dir / f"gmotif_{timestamp}.csv"
        gmotif_results = self.batch_gmotif_detection(
            pdb_sources, output_csv=str(gmotif_csv)
        )
        all_results["analyses"]["gmotif"] = gmotif_results
        
        # 2. PPI 分析（如果提供了链对）
        if chain_pairs:
            print("\n" + "="*60)
            print("📍 步骤 2/4: PPI 界面分析")
            print("="*60)
            ppi_csv = self.output_dir / f"ppi_{timestamp}.csv"
            ppi_results = self.batch_ppi_analysis(
                pdb_sources, chain_pairs, output_csv=str(ppi_csv)
            )
            all_results["analyses"]["ppi"] = ppi_results
        else:
            print("\n⏭️ 跳过 PPI 分析（未提供链对）")
        
        # 3. 口袋检测
        print("\n" + "="*60)
        print("📍 步骤 3/4: 口袋检测")
        print("="*60)
        pocket_csv = self.output_dir / f"pockets_{timestamp}.csv"
        pocket_results = self.batch_pocket_detection(
            pdb_sources, output_csv=str(pocket_csv)
        )
        all_results["analyses"]["pockets"] = pocket_results
        
        # 4. 相互作用分析
        print("\n" + "="*60)
        print("📍 步骤 4/4: 蛋白-配体相互作用分析")
        print("="*60)
        interaction_csv = self.output_dir / f"interactions_{timestamp}.csv"
        interaction_results = self.batch_interaction_analysis(
            pdb_sources, ligand_resnames, output_csv=str(interaction_csv)
        )
        all_results["analyses"]["interactions"] = interaction_results
        
        # 保存汇总 JSON
        summary_json = self.output_dir / f"batch_summary_{timestamp}.json"
        with open(summary_json, 'w', encoding='utf-8') as f:
            json.dump(all_results, f, indent=2, ensure_ascii=False, default=str)
        
        print(f"\n{'🔬'*30}")
        print(f"✅ 综合批量分析完成")
        print(f"{'🔬'*30}")
        print(f"结果保存至: {self.output_dir}")
        print(f"汇总文件: {summary_json.name}")
        print(f"{'🔬'*30}\n")
        
        return all_results
    
    # ==================== CSV 导出方法 ====================
    
    def _export_gmotif_csv(self, results: List[Dict], output_csv: str):
        """导出 G-motif 结果到 CSV"""
        with open(output_csv, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Source', 'Object', 'Success', 'Hit_Count',
                'Chain', 'Start', 'End', 'Sequence', 'RMSD'
            ])
            
            for r in results:
                if r["success"] and r.get("hits"):
                    for hit in r["hits"]:
                        writer.writerow([
                            r["source"],
                            r.get("obj_name", ""),
                            "Yes",
                            r.get("hit_count", 0),
                            hit[0] if isinstance(hit, tuple) else hit.get("chain", ""),
                            hit[1] if isinstance(hit, tuple) else hit.get("start", ""),
                            hit[2] if isinstance(hit, tuple) else hit.get("end", ""),
                            hit[3] if isinstance(hit, tuple) else hit.get("sequence", ""),
                            hit[4] if isinstance(hit, tuple) else hit.get("rmsd", "")
                        ])
                else:
                    writer.writerow([
                        r["source"],
                        r.get("obj_name", ""),
                        "No" if not r["success"] else "Yes",
                        0,
                        "", "", "", "", ""
                    ])
        
        print(f"📄 G-motif 结果已保存: {output_csv}")
    
    def _export_ppi_csv(self, results: List[Dict], output_csv: str):
        """导出 PPI 结果到 CSV"""
        with open(output_csv, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Source', 'Object', 'Success', 'Chains1', 'Chains2',
                'Interface_Contacts', 'Interface_Strength', 'BSA', 'Is_Strong'
            ])
            
            for r in results:
                writer.writerow([
                    r["source"],
                    r.get("obj_name", ""),
                    "Yes" if r["success"] else "No",
                    ",".join(r.get("chains1", [])) if r.get("chains1") else "",
                    ",".join(r.get("chains2", [])) if r.get("chains2") else "",
                    r.get("interface_contacts", ""),
                    f"{r.get('interface_strength', 0):.2f}" if r.get("interface_strength") else "",
                    f"{r.get('bsa', 0):.1f}" if r.get("bsa") else "",
                    "Yes" if r.get("is_strong_interface") else "No"
                ])
        
        print(f"📄 PPI 结果已保存: {output_csv}")
    
    def _export_pocket_csv(self, results: List[Dict], output_csv: str):
        """导出口袋检测结果到 CSV"""
        with open(output_csv, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Source', 'Object', 'Success', 'Pocket_Count',
                'Pocket_ID', 'Volume', 'Depth', 'Druggability'
            ])
            
            for r in results:
                if r["success"] and r.get("pockets"):
                    for p in r["pockets"]:
                        writer.writerow([
                            r["source"],
                            r.get("obj_name", ""),
                            "Yes",
                            r.get("pocket_count", 0),
                            p.get("id", ""),
                            f"{p.get('volume', 0):.1f}" if p.get("volume") else "",
                            f"{p.get('depth', 0):.1f}" if p.get("depth") else "",
                            f"{p.get('druggability', 0):.3f}" if p.get("druggability") else ""
                        ])
                else:
                    writer.writerow([
                        r["source"],
                        r.get("obj_name", ""),
                        "No" if not r["success"] else "Yes",
                        0,
                        "", "", "", ""
                    ])
        
        print(f"📄 口袋检测结果已保存: {output_csv}")
    
    def _export_interaction_csv(self, results: List[Dict], output_csv: str):
        """导出相互作用结果到 CSV"""
        with open(output_csv, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Source', 'Object', 'Success', 'Ligand',
                'Total_Interactions', 'HBonds', 'Salt_Bridges',
                'Hydrophobic', 'Pi_Pi', 'Cation_Pi'
            ])
            
            for r in results:
                type_counts = r.get("type_counts", {})
                writer.writerow([
                    r["source"],
                    r.get("obj_name", ""),
                    "Yes" if r["success"] else "No",
                    r.get("ligand", ""),
                    r.get("interaction_count", 0),
                    type_counts.get("Hydrogen Bond", 0),
                    type_counts.get("Salt Bridge", 0),
                    type_counts.get("Hydrophobic", 0),
                    type_counts.get("Pi-Pi Stacking", 0),
                    type_counts.get("Pi-Cation", 0)
                ])
        
        print(f"📄 相互作用结果已保存: {output_csv}")


# ==================== PyMOL 命令封装 ====================

def batch_gmotif(pdb_list: str, output_csv: Optional[str] = None, 
                 rmsd_cutoff: float = 3.5, template: str = "GSPT1"):
    """
    PyMOL 命令: 批量 G-motif 检测
    
    用法:
        batch_gmotif "6H0G,6H0F,5FQD"
        batch_gmotif "6H0G,6H0F,5FQD", output_csv="gmotif_results.csv"
        batch_gmotif "/path/to/file1.pdb,/path/to/file2.pdb"
    """
    pdb_sources = [p.strip() for p in pdb_list.split(",")]
    analyzer = BatchAnalyzer()
    return analyzer.batch_gmotif_detection(
        pdb_sources, 
        rmsd_cutoff=float(rmsd_cutoff),
        template=template,
        output_csv=output_csv
    )


def batch_ppi(pdb_list: str, chains1: str, chains2: str, 
              output_csv: Optional[str] = None, interface_dist: float = 4.5):
    """
    PyMOL 命令: 批量 PPI 分析
    
    用法:
        batch_ppi "6H0G,6H0F", "A", "B"
        batch_ppi "6H0G,6H0F", "A,C", "B,D", output_csv="ppi_results.csv"
    """
    pdb_sources = [p.strip() for p in pdb_list.split(",")]
    c1 = [c.strip() for c in chains1.split(",")]
    c2 = [c.strip() for c in chains2.split(",")]
    
    analyzer = BatchAnalyzer()
    return analyzer.batch_ppi_analysis(
        pdb_sources,
        chain_pairs=[(c1, c2)],
        interface_distance=float(interface_dist),
        output_csv=output_csv
    )


def batch_pockets(pdb_list: str, output_csv: Optional[str] = None,
                  min_volume: float = 30.0, min_depth: float = 2.5):
    """
    PyMOL 命令: 批量口袋检测
    
    用法:
        batch_pockets "6H0G,6H0F,5FQD"
        batch_pockets "6H0G,6H0F", output_csv="pocket_results.csv"
    """
    pdb_sources = [p.strip() for p in pdb_list.split(",")]
    analyzer = BatchAnalyzer()
    return analyzer.batch_pocket_detection(
        pdb_sources,
        min_volume=float(min_volume),
        min_depth=float(min_depth),
        output_csv=output_csv
    )


def batch_interactions(pdb_list: str, ligands: Optional[str] = None,
                       output_csv: Optional[str] = None, distance: float = 4.5):
    """
    PyMOL 命令: 批量相互作用分析
    
    用法:
        batch_interactions "6H0G,6H0F"
        batch_interactions "6H0G,6H0F", ligands="CC9,LEN"
    """
    pdb_sources = [p.strip() for p in pdb_list.split(",")]
    lig_list = [l.strip() for l in ligands.split(",")] if ligands else None
    
    analyzer = BatchAnalyzer()
    return analyzer.batch_interaction_analysis(
        pdb_sources,
        ligand_resnames=lig_list,
        distance_cutoff=float(distance),
        output_csv=output_csv
    )


# 注册 PyMOL 命令
if cmd:
    cmd.extend("batch_gmotif", batch_gmotif)
    cmd.extend("batch_ppi", batch_ppi)
    cmd.extend("batch_pockets", batch_pockets)
    cmd.extend("batch_interactions", batch_interactions)


if __name__ == "__main__":
    print("GLINT Batch Analyzer")
    print("=====================")
    print("PyMOL 命令:")
    print("  batch_gmotif '6H0G,6H0F,5FQD'")
    print("  batch_ppi '6H0G,6H0F', 'A', 'B'")
    print("  batch_pockets '6H0G,6H0F'")
    print("  batch_interactions '6H0G,6H0F'")