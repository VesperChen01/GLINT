# -*- coding: utf-8 -*-
"""
改进的 EC 分析模块
添加了更好的错误处理、诊断信息和 fallback 机制

主要改进：
1. 详细的依赖检查和诊断
2. 多层 fallback 机制
3. 更好的错误消息和建议
4. 性能监控和进度显示
5. 自动降级策略
Author: GLINT Team
"""

from __future__ import print_function
import os
import sys
import time
import tempfile
from typing import Dict, List, Tuple, Optional, Any

# ========== 依赖检查 ==========

class DependencyChecker:
    """依赖检查和诊断工具"""
    
    def __init__(self):
        self.status = {}
        self.warnings = []
        self.errors = []
    
    def check_core_dependencies(self) -> bool:
        """检查核心依赖"""
        print("[EC] 检查核心依赖...")
        
        core_deps = {
            'numpy': 'NumPy',
            'scipy': 'SciPy',
            'rdkit': 'RDKit',
        }
        
        all_available = True
        for import_name, display_name in core_deps.items():
            try:
                __import__(import_name)
                print(f"  ✅ {display_name}")
                self.status[import_name] = True
            except ImportError:
                print(f"  ❌ {display_name} 未安装")
                self.status[import_name] = False
                self.errors.append(f"{display_name} 是必需的")
                all_available = False
        
        return all_available
    
    def check_ec_dependencies(self) -> Dict[str, bool]:
        """检查 EC 专用依赖"""
        print("[EC] 检查 EC 分析依赖...")
        
        ec_deps = {
            'pdb2pqr': 'PDB2PQR',
            'apbs': 'APBS',
        }
        
        results = {}
        for import_name, display_name in ec_deps.items():
            try:
                __import__(import_name)
                print(f"  ✅ {display_name} Python API")
                results[import_name] = True
                self.status[import_name] = True
            except ImportError:
                print(f"  ⚠️  {display_name} Python API 未安装")
                results[import_name] = False
                self.status[import_name] = False
                self.warnings.append(f"{display_name} 可通过命令行工具使用")
        
        return results
    
    def check_external_tools(self) -> Dict[str, bool]:
        """检查外部命令行工具"""
        print("[EC] 检查外部工具...")
        
        import shutil
        
        tools = {
            'pdb2pqr': 'PDB2PQR',
            'apbs': 'APBS',
        }
        
        results = {}
        for cmd, display_name in tools.items():
            path = shutil.which(cmd)
            if path:
                print(f"  ✅ {display_name} CLI: {path}")
                results[cmd] = True
                self.status[f"{cmd}_cli"] = True
            else:
                print(f"  ⚠️  {display_name} CLI 未找到")
                results[cmd] = False
                self.status[f"{cmd}_cli"] = False
        
        return results
    
    def get_status_report(self) -> str:
        """获取状态报告"""
        report = []
        report.append("\n" + "="*60)
        report.append("EC 分析依赖状态报告")
        report.append("="*60)
        
        for dep, available in self.status.items():
            status_str = "✅ 可用" if available else "❌ 不可用"
            report.append(f"  {dep}: {status_str}")
        
        if self.warnings:
            report.append("\n⚠️  警告:")
            for warning in self.warnings:
                report.append(f"  - {warning}")
        
        if self.errors:
            report.append("\n❌ 错误:")
            for error in self.errors:
                report.append(f"  - {error}")
        
        report.append("="*60 + "\n")
        return "\n".join(report)
    
    def print_installation_guide(self):
        """打印安装指南"""
        print("\n" + "="*60)
        print("EC 分析依赖安装指南")
        print("="*60)
        
        print("\n1. 安装 PDB2PQR:")
        print("   pip3 install pdb2pqr")
        
        print("\n2. 安装 APBS:")
        print("   pip3 install apbs")
        print("   或 (macOS): brew install brewsci/bio/apbs")
        
        print("\n3. 验证安装:")
        print("   python3 ec_dependency_checker.py")
        
        print("\n" + "="*60 + "\n")


class ECAnalysisManager:
    """EC 分析管理器，处理 fallback 和错误恢复"""
    
    def __init__(self):
        self.checker = DependencyChecker()
        self.has_pdb2pqr = False
        self.has_apbs = False
        self.has_pdb2pqr_cli = False
        self.has_apbs_cli = False
        self.fallback_mode = False
    
    def initialize(self) -> bool:
        """初始化 EC 分析环境"""
        print("\n" + "="*60)
        print("初始化 EC 分析环境")
        print("="*60 + "\n")
        
        # 检查核心依赖
        if not self.checker.check_core_dependencies():
            print("\n❌ 核心依赖缺失，无法继续")
            print(self.checker.get_status_report())
            return False
        
        # 检查 EC 依赖
        ec_deps = self.checker.check_ec_dependencies()
        self.has_pdb2pqr = ec_deps.get('pdb2pqr', False)
        self.has_apbs = ec_deps.get('apbs', False)
        
        # 检查外部工具
        ext_tools = self.checker.check_external_tools()
        self.has_pdb2pqr_cli = ext_tools.get('pdb2pqr', False)
        self.has_apbs_cli = ext_tools.get('apbs', False)
        
        # 判断是否需要 fallback
        if not (self.has_pdb2pqr or self.has_pdb2pqr_cli):
            print("\n⚠️  PDB2PQR 不可用，将使用简化的 PQR 转换")
            self.fallback_mode = True
        
        if not (self.has_apbs or self.has_apbs_cli):
            print("\n⚠️  APBS 不可用，将使用 Coulomb 近似计算静电势")
            self.fallback_mode = True
        
        print(self.checker.get_status_report())
        
        if self.fallback_mode:
            print("💡 已启用 Fallback 模式")
            print("   功能将受到限制，但基本 EC 分析仍可进行")
        
        return True
    
    def get_pdb2pqr_runner(self):
        """获取 PDB2PQR 运行器"""
        from glint.ligand_ec_calculator import PDB2PQRRunner
        return PDB2PQRRunner()
    
    def get_apbs_runner(self):
        """获取 APBS 运行器"""
        from glint.ligand_ec_calculator import APBSRunner
        return APBSRunner()
    
    def run_ec_analysis(self, obj_name: str, ligand_resname: str,
                       output_dir: str = None, **kwargs) -> Optional[Dict[str, Any]]:
        """运行 EC 分析，带有完整的错误处理"""
        
        if not self.initialize():
            return None
        
        try:
            from glint.ligand_ec_calculator import calculate_ligand_ec
            
            print("\n" + "="*60)
            print("开始 EC 分析")
            print("="*60 + "\n")
            
            result = calculate_ligand_ec(
                obj_name=obj_name,
                ligand_resname=ligand_resname,
                output_dir=output_dir,
                **kwargs
            )
            
            if result:
                print("\n✅ EC 分析完成")
                return result
            else:
                print("\n❌ EC 分析失败")
                return None
        
        except Exception as e:
            print(f"\n❌ EC 分析出错: {e}")
            print("\n建议:")
            print("  1. 检查依赖: python3 ec_dependency_checker.py")
            print("  2. 安装缺失的工具: bash install_ec_dependencies.sh")
            print("  3. 查看详细日志")
            return None


# ========== 快速诊断函数 ==========

def quick_diagnose() -> bool:
    """快速诊断 EC 分析环境"""
    manager = ECAnalysisManager()
    return manager.initialize()


def print_help():
    """打印帮助信息"""
    print("""
GLINT EC 分析 - 快速开始指南
================================

1. 检查依赖:
   python3 -c "from glint.ec_analysis_improved import quick_diagnose; quick_diagnose()"

2. 安装缺失的依赖:
   bash install_ec_dependencies.sh

3. 在 PyMOL 中运行 EC 分析:
   calculate_ligand_ec 'complex', 'LIG', output_dir='./ec_output'

4. 进行三元复合物分析:
   analyze_ternary_ec 'complex', 'GLUE', ['A'], ['B']

常见问题:
---------

Q: PDB2PQR 未找到
A: 运行 pip3 install pdb2pqr

Q: APBS 未找到
A: 运行 pip3 install apbs 或 brew install brewsci/bio/apbs

Q: 分析失败
A: 运行诊断工具检查依赖状态

更多信息: https://github.com/VesperChen01/GLINT
    """)


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--help':
        print_help()
    else:
        quick_diagnose()
