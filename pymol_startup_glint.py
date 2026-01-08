# PyMOL 启动脚本 - GLINT 开发环境
# 保存为: ~/.pymol/startup/01_glint_dev.py

import sys
import os

# 添加开发路径（优先级最高）
dev_path = "/Users/vesper/Desktop/git/GlueTK"
if dev_path not in sys.path:
    sys.path.insert(0, dev_path)

# 添加 PyMOL 启动目录
startup_path = os.path.expanduser("~/.pymol/startup")
if startup_path not in sys.path:
    sys.path.insert(0, startup_path)

print("✅ GLINT dev paths loaded")
print(f"   - Dev: {dev_path}")
print(f"   - Startup: {startup_path}")

# 自动导入 GLINT
try:
    import glint
    print(f"✅ GLINT v{glint.__version__} loaded")
    print("   Type 'glint_gui' to launch GUI")
except ImportError as e:
    print(f"⚠️  GLINT import failed: {e}")
    print("   Check if 'glint' directory exists in:")
    print(f"   - {dev_path}/glint")
    print(f"   - {startup_path}/glint")

