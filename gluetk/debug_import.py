
import sys
import os

# 模拟 PyMOL 环境路径
package_dir = "/Volumes/data/git/GlueTK"
if package_dir not in sys.path:
    sys.path.insert(0, package_dir)

print(f"Checking import of gluetk.gui.main_window from {package_dir}...")

try:
    import gluetk.gui.main_window
    print("✅ Import successful!")
except Exception as e:
    print(f"❌ Import failed: {e}")
    import traceback
    traceback.print_exc()
