#!/bin/bash
# 同步 GLINT 修改到 PyMOL 启动目录

set -e

PYMOL_STARTUP="$HOME/.pymol/startup/glint"
SOURCE_DIR="./glint"

echo "🔄 同步 GLINT 文件到 PyMOL 启动目录..."
echo "   源目录: $SOURCE_DIR"
echo "   目标目录: $PYMOL_STARTUP"
echo ""

# 检查目标目录是否存在
if [ ! -d "$PYMOL_STARTUP" ]; then
    echo "❌ PyMOL 启动目录不存在: $PYMOL_STARTUP"
    echo "   请先运行 install_glint.sh 安装 GLINT"
    exit 1
fi

# 检查源目录是否存在
if [ ! -d "$SOURCE_DIR" ]; then
    echo "❌ 源目录不存在: $SOURCE_DIR"
    exit 1
fi

# 同步文件（排除缓存和临时文件）
echo "📦 复制文件..."
rsync -av --delete \
      --exclude='__pycache__' \
      --exclude='*.pyc' \
      --exclude='.DS_Store' \
      --exclude='*.sh' \
      "$SOURCE_DIR/" "$PYMOL_STARTUP/"

echo ""
echo "✅ 同步完成！"
echo ""
echo "💡 提示："
echo "   1. 如果 PyMOL 正在运行，请重启 PyMOL 以加载更新"
echo "   2. 或在 PyMOL 中运行: reinitialize"
echo "   3. 然后重新导入: import glint; glint.glint_gui()"
echo ""

