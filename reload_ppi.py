"""
GlueTK PPI 模块重新加载脚本

使用方法：
在PyMOL中运行：
run /Users/vesper/Desktop/git/GlueTK/GlueTK/reload_ppi.py
"""

import sys
import importlib

# 重新加载ppi_analyzer模块
if 'gluetk.ppi_analyzer' in sys.modules:
    print("[Reload] Reloading gluetk.ppi_analyzer...")
    importlib.reload(sys.modules['gluetk.ppi_analyzer'])
    print("[Reload] ✅ gluetk.ppi_analyzer reloaded")
elif 'ppi_analyzer' in sys.modules:
    print("[Reload] Reloading ppi_analyzer...")
    importlib.reload(sys.modules['ppi_analyzer'])
    print("[Reload] ✅ ppi_analyzer reloaded")
else:
    print("[Reload] ⚠️ ppi_analyzer not loaded yet")

print("[Reload] 请重新运行PPI分析以查看更新的颜色")
