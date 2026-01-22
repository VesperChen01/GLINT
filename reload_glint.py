# -*- coding: utf-8 -*-
"""
重新加载 GLINT 模块的脚本
在 PyMOL 中运行此脚本以重新加载修改后的代码
"""

import sys
import os
import importlib

# 添加路径
sys.path.insert(0, "/Users/vesper/Desktop/git/GlueTK")
sys.path.insert(0, os.path.expanduser("~/.pymol/startup"))

print("=" * 60)
print("重新加载 GLINT 模块...")
print("=" * 60)

# 重新加载 interaction_analyzer
try:
    from glint import interaction_analyzer
    importlib.reload(interaction_analyzer)
    print("✅ interaction_analyzer 重新加载成功")
    
    # 验证是否加载了新版本(检查是否有最终清理步骤)
    import inspect
    source = inspect.getsource(interaction_analyzer.visualize_protein_ligand_3d)
    if "最终清理" in source:
        print("✅ 确认加载了修复后的版本")
    else:
        print("⚠️  警告: 可能仍在使用旧版本,建议重启 PyMOL")
except Exception as e:
    print(f"❌ interaction_analyzer 重新加载失败: {e}")

# 重新加载 GUI utils (如果通过 GUI 调用)
try:
    from glint.gui import utils
    importlib.reload(utils)
    print("✅ GUI utils 重新加载成功")
except Exception as e:
    print(f"⚠️  GUI utils 重新加载失败 (如果不使用 GUI 可以忽略): {e}")

print("=" * 60)
print("重新加载完成!")
print("=" * 60)
print()
print("现在可以重新运行分析:")
print("  result = interaction_analyzer.analyze_protein_ligand_interactions('complex', '***')")
print("  interaction_analyzer.visualize_protein_ligand_3d('complex', result, '***')")
print()
print("验证修复成功的标志:")
print("  - 日志中应该看到 '🧹 最终清理: 确保只显示配体和相互作用残基'")
print("  - 不应该看到 'Selector-Error: Invalid selection name \"lig_pocket\"'")
print("  - 只有相互作用的残基显示为 licorice")
print()

