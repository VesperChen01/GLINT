# -*- coding: utf-8 -*-
"""
Ternary Complex Evaluation Worker
三元复合物评估的后台工作线程
"""
import traceback
from typing import Optional, Dict, Any

from .qt_adapter import QThread, Signal as pyqtSignal


class TernaryEvaluationWorker(QThread):
    """三元复合物评估Worker"""
    progress = pyqtSignal(str)
    finished = pyqtSignal(dict, str)  # (result_dict, out_csv)
    error = pyqtSignal(str)
    
    def __init__(self, pdb_path: str, e3_chain: str, poi_chain: str,
                 lig_resn: str, smiles: Optional[str] = None,
                 out_csv: Optional[str] = None, obj_name: Optional[str] = None):
        super().__init__()
        self.pdb_path = pdb_path
        self.e3_chain = e3_chain
        self.poi_chain = poi_chain
        self.lig_resn = lig_resn
        self.smiles = smiles
        self.out_csv = out_csv
        self.obj_name = obj_name
        
    def run(self):
        try:
            self.progress.emit("Loading ternary complex evaluator...")
            
            from ..ternary_complex_evaluator import TernaryComplexEvaluator
            
            # 创建评估器（不传参数）
            evaluator = TernaryComplexEvaluator()
            
            self.progress.emit("Evaluating ternary complex...")
            
            # 使用evaluate方法进行完整评估
            features = evaluator.evaluate(
                pdb_path=self.pdb_path,
                e3_chain=self.e3_chain,
                poi_chain=self.poi_chain,
                ligand_resn=self.lig_resn,
                ligand_smiles=self.smiles,
                obj_name=self.obj_name
            )
            
            # 转换为字典
            result = evaluator.to_dict(features)
            
            # 保存CSV
            if self.out_csv:
                import csv
                with open(self.out_csv, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Feature', 'Value'])
                    for k, v in result.items():
                        writer.writerow([k, v])
                self.progress.emit(f"Results saved to {self.out_csv}")
            
            self.finished.emit(result, self.out_csv or "")
            
        except Exception as e:
            self.error.emit(f"{str(e)}\n{traceback.format_exc()}")