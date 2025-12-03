#!/bin/bash
# GlueTK 快速修复脚本
# 用于修复已安装但无法启动的 GlueTK

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

echo -e "${BLUE}🔧 GlueTK 修复工具${NC}"
echo "================================"

# 查找 conda
CONDA_EXE=""
if command -v conda &> /dev/null; then
    CONDA_EXE=$(command -v conda)
else
    for p in "$HOME/miniconda3/bin/conda" \
             "$HOME/anaconda3/bin/conda" \
             "/opt/homebrew/Caskroom/miniconda/base/bin/conda" \
             "/opt/homebrew/bin/conda"; do
        if [ -x "$p" ]; then
            CONDA_EXE="$p"
            break
        fi
    done
fi

if [ -z "$CONDA_EXE" ]; then
    echo -e "${RED}❌ 未找到 Conda${NC}"
    exit 1
fi

echo -e "${GREEN}✅ 找到 Conda: $CONDA_EXE${NC}"
eval "$($CONDA_EXE shell.bash hook)"

# 检查 gluetk 环境
if ! conda env list | grep -q "^gluetk "; then
    echo -e "${RED}❌ gluetk 环境不存在${NC}"
    echo "请运行完整安装: bash install_gluetk.sh"
    exit 1
fi

echo -e "${GREEN}✅ gluetk 环境存在${NC}"

# 检查 PyMOL 安装方式
echo ""
echo -e "${BLUE}检测 PyMOL 安装...${NC}"

PYMOL_TYPE=""
PYMOL_PATH=""

# 1. 检查系统 PyMOL.app
if [ -d "/Applications/PyMOL.app" ]; then
    PYMOL_APP_BIN="/Applications/PyMOL.app/Contents/MacOS/PyMOL"
    if [ -x "$PYMOL_APP_BIN" ]; then
        PYMOL_TYPE="app"
        PYMOL_PATH="$PYMOL_APP_BIN"
        echo -e "${GREEN}✅ 找到 PyMOL.app${NC}"
    fi
fi

# 2. 检查 conda 环境中的 PyMOL
if [ -z "$PYMOL_TYPE" ]; then
    echo -e "${YELLOW}未找到 PyMOL.app，检查 conda 环境...${NC}"
    if conda run -n gluetk which pymol &> /dev/null; then
        PYMOL_TYPE="conda"
        echo -e "${GREEN}✅ conda 环境中有 PyMOL${NC}"
    else
        echo -e "${YELLOW}⚠️  conda 环境中没有 PyMOL，正在安装...${NC}"
        if conda install -n gluetk -c conda-forge pymol-open-source -y; then
            PYMOL_TYPE="conda"
            echo -e "${GREEN}✅ PyMOL 安装成功${NC}"
        else
            echo -e "${RED}❌ PyMOL 安装失败${NC}"
            echo ""
            echo "建议方案："
            echo "  1. 手动安装 PyMOL.app: https://pymol.org/"
            echo "  2. 然后重新运行此脚本"
            exit 1
        fi
    fi
fi

# 3. 重建桌面 App
echo ""
echo -e "${BLUE}重建桌面应用...${NC}"

DESKTOP_APP="$HOME/Desktop/GlueTK.app"
if [ -d "$DESKTOP_APP" ]; then
    rm -rf "$DESKTOP_APP"
fi

CONTENTS="$DESKTOP_APP/Contents"
MACOS="$CONTENTS/MacOS"
RESOURCES="$CONTENTS/Resources"

mkdir -p "$MACOS"
mkdir -p "$RESOURCES"

# 查找并复制图标
ICON_FOUND=false
for icon_path in \
    "$HOME/.pymol/startup/gluetk/assets/AppIcon.icns" \
    "$(pwd)/gluetk/assets/AppIcon.icns"; do
    if [ -f "$icon_path" ]; then
        cp "$icon_path" "$RESOURCES/AppIcon.icns"
        ICON_FOUND=true
        break
    fi
done

if [ "$ICON_FOUND" = false ]; then
    echo -e "${YELLOW}⚠️  未找到图标文件${NC}"
fi

# 创建 Info.plist
cat > "$CONTENTS/Info.plist" << 'EOF'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleExecutable</key>
    <string>launcher</string>
    <key>CFBundleIconFile</key>
    <string>AppIcon</string>
    <key>CFBundleIdentifier</key>
    <string>com.vesper.gluetk</string>
    <key>CFBundleName</key>
    <string>GlueTK</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleShortVersionString</key>
    <string>1.0.0</string>
</dict>
</plist>
EOF

# 创建启动脚本
if [ "$PYMOL_TYPE" = "app" ]; then
    # 使用 PyMOL.app
    cat > "$MACOS/launcher" << LAUNCHER_APP_EOF
#!/bin/bash
# GlueTK Launcher (PyMOL.app)

# 激活 conda 环境以加载依赖
export CONDA_EXE="$CONDA_EXE"
eval "\$(\$CONDA_EXE shell.bash hook)" 2>/dev/null
conda activate gluetk 2>/dev/null || true

# 使用系统 PyMOL
"$PYMOL_PATH" -d "import sys, os; sys.path.insert(0, os.path.expanduser('~/.pymol/startup')); import gluetk; gluetk.gluetk_gui()"
LAUNCHER_APP_EOF
else
    # 使用 conda PyMOL
    cat > "$MACOS/launcher" << LAUNCHER_CONDA_EOF
#!/bin/bash
# GlueTK Launcher (conda pymol)

CONDA_EXE="$CONDA_EXE"
"\$CONDA_EXE" run -n gluetk pymol -d "import sys, os; sys.path.insert(0, os.path.expanduser('~/.pymol/startup')); import gluetk; gluetk.gluetk_gui()"
LAUNCHER_CONDA_EOF
fi

chmod +x "$MACOS/launcher"

echo -e "${GREEN}✅ 桌面应用已重建${NC}"
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}✅ 修复完成！${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "现在可以双击桌面上的 GlueTK.app 启动"
echo ""
if [ "$PYMOL_TYPE" = "app" ]; then
    echo -e "${BLUE}使用的 PyMOL: 系统 PyMOL.app${NC}"
else
    echo -e "${YELLOW}使用的 PyMOL: conda 环境${NC}"
    echo -e "${YELLOW}💡 建议安装 PyMOL.app 以获得更好性能${NC}"
fi
echo ""
