# -*- coding: utf-8 -*-
"""
重新加载 GLINT 所有模块的脚本
在 PyMOL 中运行: run /path/to/reload_glint.py
"""

import sys
import os
import importlib
import inspect

# 动态检测项目路径
_script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else os.getcwd()
sys.path.insert(0, _script_dir)
sys.path.insert(0, os.path.expanduser("~/.pymol/startup"))

print("=" * 60)
print("🔄 重新加载 GLINT 所有模块...")
print(f"   脚本路径: {_script_dir}")
print("=" * 60)

# ============================================================
# 第一步：从 sys.modules 中清除所有 glint 模块缓存
# 这是确保重新加载生效的关键步骤
# ============================================================
glint_modules = [key for key in sys.modules if key == 'glint' or key.startswith('glint.')]
for mod_name in sorted(glint_modules):
    del sys.modules[mod_name]
print(f"🗑️  已清除 {len(glint_modules)} 个 glint 模块缓存")

# ============================================================
# 第二步：重新导入所有核心模块
# ============================================================
reload_targets = [
    ("glint", "glint 包初始化"),
    ("glint.color_scheme", "统一颜色方案"),
    ("glint.interaction_types", "相互作用类型定义"),
    ("glint.interaction_2d_plot", "2D 相互作用图"),
    ("glint.interaction_analyzer", "相互作用分析器"),
    ("glint.ppi_analyzer", "蛋白-蛋白相互作用分析"),
    ("glint.ligand_ligand_analyzer", "配体-配体相互作用分析"),
    ("glint.ec_visualization", "EC 可视化"),
    ("glint.ligand_ec_calculator", "配体 EC 计算器"),
    ("glint.gui.utils", "GUI 工具函数"),
]

success_count = 0
fail_count = 0

for mod_name, description in reload_targets:
    try:
        mod = importlib.import_module(mod_name)
        print(f"  ✅ {mod_name} ({description})")
        success_count += 1
    except Exception as e:
        print(f"  ❌ {mod_name} ({description}): {e}")
        fail_count += 1

# ============================================================
# 第三步：验证关键模块是否为最新版本
# ============================================================
print()
print("-" * 60)
print("🔍 版本验证:")

# 验证 2D 布局算法（6阶段 Discovery Studio 风格）
try:
    from glint import interaction_2d_plot
    source = inspect.getsource(interaction_2d_plot)
    if "Phase 1" in source and "辐射扫描" in source:
        print("  ✅ interaction_2d_plot: 6阶段 Discovery Studio 布局算法 (最新)")
    else:
        print("  ⚠️  interaction_2d_plot: 可能是旧版本布局算法")
except Exception as e:
    print(f"  ❌ interaction_2d_plot 验证失败: {e}")

# 验证 interaction_analyzer
try:
    from glint import interaction_analyzer
    source = inspect.getsource(interaction_analyzer)
    if "最终清理" in source:
        print("  ✅ interaction_analyzer: 包含最终清理逻辑 (最新)")
    else:
        print("  ⚠️  interaction_analyzer: 可能缺少最终清理逻辑")
except Exception as e:
    print(f"  ❌ interaction_analyzer 验证失败: {e}")

# 验证版本号
try:
    from glint._version import __version__
    print(f"  📌 GLINT 版本: {__version__}")
except Exception:
    print("  ⚠️  无法读取版本号")

print()
print("=" * 60)
print(f"🎉 重新加载完成! 成功: {success_count}, 失败: {fail_count}")
if fail_count > 0:
    print("   ⚠️ 部分模块加载失败，建议完全退出并重启 PyMOL")
else:
    print("   所有模块已更新为最新版本，可以直接使用")
print("=" * 60)

