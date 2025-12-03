#!/bin/bash
# GlueTK 一键安装脚本 (改进版)
# 支持多种 PyMOL 安装方式，更健壮的依赖检查

set -e

# 配置
ENV_NAME="gluetk"
PYTHON_VERSION="3.9"
INSTALL_DIR="$HOME/.pymol/startup/gluetk"
DESKTOP_APP="$HOME/Desktop/GlueTK.app"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🧬 GlueTK 安装脚本${NC}"
echo "================================"

# 1. 检查 Conda
echo -e "\n${BLUE}[1/6]${NC} 检查 Conda 环境..."
CONDA_EXE=""
if command -v conda &> /dev/null; then
    CONDA_EXE=$(command -v conda)
else
    # 尝试常见路径
    for p in "$HOME/miniconda3/bin/conda" \
             "$HOME/anaconda3/bin/conda" \
             "/opt/miniconda3/bin/conda" \
             "/opt/anaconda3/bin/conda" \
             "/usr/local/bin/conda" \
             "/opt/homebrew/bin/conda" \
             "/opt/homebrew/Caskroom/miniconda/base/bin/conda" \
             "$HOME/opt/miniconda3/bin/conda"; do
        if [ -x "$p" ]; then
            CONDA_EXE="$p"
            break
        fi
    done
fi

if [ -z "$CONDA_EXE" ]; then
    echo -e "${RED}❌ 未找到 Conda${NC}"
    echo "请先安装 Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo -e "${GREEN}✅ 找到 Conda: $CONDA_EXE${NC}"

# 初始化 conda
eval "$($CONDA_EXE shell.bash hook)"

# 2. 检查/创建 conda 环境
echo -e "\n${BLUE}[2/6]${NC} 检查 Conda 环境 '$ENV_NAME'..."
if conda env list | grep -q "^${ENV_NAME} "; then
    echo -e "${YELLOW}⚠️  环境已存在，将更新依赖${NC}"
else
    echo "   创建新环境..."
    conda create -n "$ENV_NAME" python=$PYTHON_VERSION -y
    echo -e "${GREEN}✅ 环境创建成功${NC}"
fi

# 3. 检测 PyMOL 安装方式
echo -e "\n${BLUE}[3/6]${NC} 检测 PyMOL..."
PYMOL_TYPE=""
PYMOL_PATH=""

# 3a. 检查系统 PyMOL.app (macOS 商业版/开源版)
if [ -d "/Applications/PyMOL.app" ]; then
    PYMOL_APP_BIN="/Applications/PyMOL.app/Contents/MacOS/PyMOL"
    if [ -x "$PYMOL_APP_BIN" ]; then
        PYMOL_TYPE="app"
        PYMOL_PATH="$PYMOL_APP_BIN"
        echo -e "${GREEN}✅ 找到 PyMOL.app: $PYMOL_PATH${NC}"
    fi
fi

# 3b. 检查 conda 环境中的 pymol
if [ -z "$PYMOL_TYPE" ]; then
    echo "   系统中未找到 PyMOL.app，将在 conda 环境中安装..."
    echo "   这可能需要几分钟..."
    
    # 尝试安装 pymol-open-source
    if conda install -n "$ENV_NAME" -c conda-forge pymol-open-source -y; then
        # 验证安装
        if conda run -n "$ENV_NAME" which pymol &> /dev/null; then
            PYMOL_TYPE="conda"
            PYMOL_PATH="conda"
            echo -e "${GREEN}✅ PyMOL 已安装到 conda 环境${NC}"
        else
            echo -e "${RED}❌ PyMOL 安装失败${NC}"
            echo "建议手动安装 PyMOL.app: https://pymol.org/"
            exit 1
        fi
    else
        echo -e "${RED}❌ 无法通过 conda 安装 PyMOL${NC}"
        echo "请手动安装 PyMOL.app 到 /Applications/ 目录"
        exit 1
    fi
fi

# 4. 安装依赖包
echo -e "\n${BLUE}[4/6]${NC} 安装依赖包..."
echo "   这可能需要几分钟，请耐心等待..."

PACKAGES=(
    "rdkit"
    "scipy"
    "matplotlib"
    "pillow"
    "numpy"
    "pandas"
    "seaborn"
    "pyqt"
    "openbabel"
    "requests"
)

conda install -n "$ENV_NAME" -c conda-forge "${PACKAGES[@]}" -y

echo -e "${GREEN}✅ 依赖包安装完成${NC}"

# 5. 复制 GlueTK 文件
echo -e "\n${BLUE}[5/6]${NC} 安装 GlueTK 插件..."
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SOURCE_DIR="$SCRIPT_DIR/gluetk"

if [ ! -d "$SOURCE_DIR" ]; then
    # 可能脚本在 gluetk 目录内
    SOURCE_DIR="$SCRIPT_DIR"
fi

if [ ! -f "$SOURCE_DIR/__init__.py" ]; then
    echo -e "${RED}❌ 找不到 GlueTK 源码${NC}"
    exit 1
fi

# 创建安装目录
mkdir -p "$(dirname "$INSTALL_DIR")"

# 复制文件
echo "   复制文件到 $INSTALL_DIR ..."
if [ -d "$INSTALL_DIR" ]; then
    rm -rf "$INSTALL_DIR"
fi

mkdir -p "$INSTALL_DIR"
rsync -av --exclude='__pycache__' \
          --exclude='*.pyc' \
          --exclude='.DS_Store' \
          --exclude='.git' \
          --exclude='*.sh' \
          "$SOURCE_DIR/" "$INSTALL_DIR/"

echo -e "${GREEN}✅ GlueTK 文件安装完成${NC}"

# 6. 创建启动器 (根据 PyMOL 类型)
echo -e "\n${BLUE}[6/6]${NC} 创建桌面应用..."

if [ -d "$DESKTOP_APP" ]; then
    rm -rf "$DESKTOP_APP"
fi

CONTENTS="$DESKTOP_APP/Contents"
MACOS="$CONTENTS/MacOS"
RESOURCES="$CONTENTS/Resources"

mkdir -p "$MACOS"
mkdir -p "$RESOURCES"

# 复制图标
if [ -f "$SOURCE_DIR/assets/AppIcon.icns" ]; then
    cp "$SOURCE_DIR/assets/AppIcon.icns" "$RESOURCES/AppIcon.icns"
fi

# 创建 Info.plist
cat > "$CONTENTS/Info.plist" << 'PLIST_EOF'
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
PLIST_EOF

# 创建启动脚本（根据 PyMOL 类型）
if [ "$PYMOL_TYPE" = "app" ]; then
    # 使用系统 PyMOL.app + conda 环境的依赖
    cat > "$MACOS/launcher" << LAUNCHER_EOF
#!/bin/bash
# GlueTK Launcher (PyMOL.app + conda dependencies)

# 初始化 conda 环境变量
export CONDA_EXE="$CONDA_EXE"
eval "\$(\$CONDA_EXE shell.bash hook)" 2>/dev/null
conda activate $ENV_NAME 2>/dev/null || true

# 启动 PyMOL
"$PYMOL_PATH" -d "import sys, os; sys.path.insert(0, os.path.expanduser('~/.pymol/startup')); import gluetk; gluetk.gluetk_gui()"
LAUNCHER_EOF
else
    # 使用 conda 环境中的 PyMOL
    cat > "$MACOS/launcher" << LAUNCHER_EOF
#!/bin/bash
# GlueTK Launcher (conda pymol)

CONDA_EXE="$CONDA_EXE"
if [ -n "\$CONDA_EXE" ] && [ -x "\$CONDA_EXE" ]; then
    echo "Starting GlueTK via conda env $ENV_NAME..."
    "\$CONDA_EXE" run -n $ENV_NAME pymol -d "import sys, os; sys.path.insert(0, os.path.expanduser('~/.pymol/startup')); import gluetk; gluetk.gluetk_gui()"
else
    osascript -e 'display alert "Error" message "Conda not found. Please reinstall GlueTK."'
    exit 1
fi
LAUNCHER_EOF
fi

chmod +x "$MACOS/launcher"

echo -e "${GREEN}✅ 桌面应用创建完成: $DESKTOP_APP${NC}"

# 完成
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}🎉 GlueTK 安装成功！${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "启动方式："
echo "  1. 双击桌面上的 GlueTK.app"
echo "  2. 或在 PyMOL 中运行:"
echo "     run ~/.pymol/startup/gluetk/__init__.py"
echo "     gluetk_gui"
echo ""
echo "PyMOL 类型: $PYMOL_TYPE"
if [ "$PYMOL_TYPE" = "app" ]; then
    echo -e "${BLUE}💡 提示: 使用的是系统 PyMOL.app，性能更好${NC}"
else
    echo -e "${YELLOW}💡 提示: 使用的是 conda PyMOL，建议安装 PyMOL.app 以获得更好体验${NC}"
fi
echo ""
