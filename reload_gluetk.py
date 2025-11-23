# -*- coding: utf-8 -*-
"""
GlueTK 重新加载脚本
在 PyMOL 中运行此脚本以重新加载 GlueTK
"""

import sys
import importlib

# 清除 GlueTK 相关模块的缓存
modules_to_reload = []
for module_name in list(sys.modules.keys()):
    if 'gluetk' in module_name.lower():
        modules_to_reload.append(module_name)

print(f"清除 {len(modules_to_reload)} 个 GlueTK 模块缓存...")
for module_name in modules_to_reload:
    del sys.modules[module_name]

# 重新导入 GlueTK
print("重新加载 GlueTK...")
try:
    import gluetk
    importlib.reload(gluetk)
    print("✅ GlueTK 重新加载成功！")
    
    # 启动 GUI
    print("启动 GUI...")
    gluetk.gluetk_gui()
    print("✅ GUI 已打开！")
    
except Exception as e:
    print(f"❌ 重新加载失败: {e}")
    import traceback
    traceback.print_exc()
