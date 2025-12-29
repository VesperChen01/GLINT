#!/bin/bash
# Linux PyInstaller Build Script for GlueTK Installer

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Building GlueTK Linux Installer with PyInstaller..."

# 清理旧的构建
rm -rf build dist

# 确保 PyInstaller 已安装
if ! command -v pyinstaller &> /dev/null; then
    echo "PyInstaller not found. Installing PyInstaller..."
    python3 -m pip install pyinstaller
fi

# 运行 PyInstaller
# --noconfirm: 不询问确认
# --onefile: 打包成一个文件
# --windowed: 无控制台窗口 (对于 GUI 应用)
# --name: 可执行文件名
# --add-data: 添加数据文件 (gluetk 源码, assets)
# --icon: 图标文件 (Linux .desktop 文件直接使用 PNG，这里可以省略或指向一个 PNG)
# --hidden-import: 隐藏导入，解决 PyInstaller 无法自动检测的模块
pyinstaller --noconfirm \
            --onefile \
            --windowed \
            --name "GlueTK_Installer" \
            --add-data "../gluetk:gluetk" \
            --add-data "../gluetk/assets:gluetk/assets" \
            --hidden-import "tkinter" \
            --hidden-import "tkinter.ttk" \
            --hidden-import "tkinter.filedialog" \
            --hidden-import "tkinter.messagebox" \
            "GlueTK_Installer.py"

if [ $? -eq 0 ]; then
    echo "✅ PyInstaller build successful!"
    echo "Installer created at: dist/GlueTK_Installer"
else
    echo "❌ PyInstaller build failed!"
    exit 1
fi

# 清理 PyInstaller 产生的临时文件
rm -rf build
rm -rf GlueTK_Installer.spec

echo "Build process finished."