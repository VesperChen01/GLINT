#!/bin/bash
# GLINT 一键安装脚本
# 支持多种 PyMOL 安装方式，更健壮的依赖检查

set -e

# 配置
ENV_NAME="glint"
PYTHON_VERSION="3.10"
INSTALL_DIR="$HOME/.pymol/startup/glint"
DESKTOP_APP="$HOME/Desktop/GLINT.app"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🧬 GLINT 安装脚本${NC}"
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
HAS_SYSTEM_PYMOL=false

# 3a. 检查系统 PyMOL.app (macOS 商业版/开源版)
if [ -d "/Applications/PyMOL.app" ]; then
    PYMOL_APP_BIN="/Applications/PyMOL.app/Contents/MacOS/PyMOL"
    if [ -x "$PYMOL_APP_BIN" ]; then
        HAS_SYSTEM_PYMOL=true
        echo -e "${GREEN}✅ 检测到系统 PyMOL.app: $PYMOL_APP_BIN${NC}"
    fi
fi

# 3b. 优先在 conda 环境中安装 PyMOL（推荐方式，避免依赖冲突）
echo "   在 conda 环境中安装 PyMOL（推荐方式）..."
echo "   这可能需要几分钟..."

# 尝试安装 pymol-open-source
if conda install -n "$ENV_NAME" -c conda-forge pymol-open-source -y; then
    # 验证安装
    if conda run -n "$ENV_NAME" which pymol &> /dev/null; then
        PYMOL_TYPE="conda"
        PYMOL_PATH="conda"
        echo -e "${GREEN}✅ PyMOL 已安装到 conda 环境${NC}"

        if [ "$HAS_SYSTEM_PYMOL" = true ]; then
            echo -e "${YELLOW}💡 提示: 检测到系统已有 PyMOL.app，但 GLINT 将使用 conda 环境中的 PyMOL${NC}"
            echo -e "${YELLOW}   这样可以避免依赖冲突，确保所有功能正常工作${NC}"
        fi
    else
        echo -e "${RED}❌ PyMOL 安装验证失败${NC}"

        if [ "$HAS_SYSTEM_PYMOL" = true ]; then
            echo -e "${YELLOW}⚠️  将尝试使用系统 PyMOL.app（可能存在依赖冲突）${NC}"
            PYMOL_TYPE="app"
            PYMOL_PATH="$PYMOL_APP_BIN"
        else
            echo -e "${RED}❌ 无法安装 PyMOL${NC}"
            echo "建议手动安装 PyMOL.app: https://pymol.org/"
            exit 1
        fi
    fi
else
    echo -e "${YELLOW}⚠️  无法通过 conda 安装 PyMOL${NC}"

    if [ "$HAS_SYSTEM_PYMOL" = true ]; then
        echo -e "${YELLOW}   将使用系统 PyMOL.app（可能存在依赖冲突）${NC}"
        PYMOL_TYPE="app"
        PYMOL_PATH="$PYMOL_APP_BIN"
    else
        echo -e "${RED}❌ 未找到可用的 PyMOL${NC}"
        echo "请手动安装 PyMOL.app 到 /Applications/ 目录"
        exit 1
    fi
fi

# 4. 安装依赖包
echo -e "\n${BLUE}[4/6]${NC} 安装依赖包..."
echo "   这可能需要几分钟，请耐心等待..."

# 核心依赖（通过 conda 安装）
CONDA_PACKAGES=(
    "rdkit"
    "scipy"
    "matplotlib"
    "pillow"
    "numpy=1.26.4" # 指定兼容 PyMOL 的版本
    "pandas"
    "seaborn"
    "pyqt"
    "openbabel"
    "requests"
    "vina"
    "pdb2pqr" # 通过 conda 安装
)

# 临时禁用 InsecureRequestWarning
export PYTHONWARNINGS="ignore::urllib3.exceptions.InsecureRequestWarning"
conda install -n "$ENV_NAME" -c conda-forge "${CONDA_PACKAGES[@]}" -y
unset PYTHONWARNINGS # 安装完成后取消设置

# 表面分析依赖（通过 pip 安装，因为 open3d 在 conda 上不稳定）
echo "   安装表面分析依赖 (Open3D, scikit-image)..."
conda run -n "$ENV_NAME" python -m pip install open3d scikit-image --quiet --disable-pip-version-check || {
    echo -e "${YELLOW}⚠️  Open3D 安装失败，表面分析将使用内置回退方案${NC}"
}
unset PYTHONWARNINGS # 安装完成后取消设置

# 安装 APBS（EC 静电互补性分析必需）
echo "   安装 APBS 求解器..."
conda run -n "$ENV_NAME" python -m pip install apbs-binary --quiet --disable-pip-version-check || {
    echo -e "${YELLOW}⚠️  apbs-binary 安装失败，尝试通过 conda 安装...${NC}"
    conda install -n "$ENV_NAME" -c conda-forge apbs -y || {
        echo -e "${YELLOW}⚠️  APBS 安装失败，EC分析将无法运行${NC}"
    }
}

echo -e "${GREEN}✅ 依赖包安装完成${NC}"

# 4b. 检查并安装 GCC (HADDOCK3/CNS 依赖)
echo -e "\n${BLUE}[4b]${NC} 检查 GCC 依赖..."
if ! brew list gcc &> /dev/null; then
    echo "   安装 GCC (HADDOCK3 需要)..."
    echo "   ${YELLOW}这可能需要 5-10 分钟，请耐心等待...${NC}"
    if brew install gcc; then
        echo -e "${GREEN}✅ GCC 安装成功${NC}"
    else
        echo -e "${YELLOW}⚠️  GCC 安装失败，HADDOCK3 可能无法正常工作${NC}"
    fi
else
    echo -e "${GREEN}✅ GCC 已安装${NC}"
fi

# 4c. 使用 pip 安装 HADDOCK3（蛋白-蛋白对接引擎）
echo -e "\n${BLUE}[4c]${NC} 安装 HADDOCK3..."
echo "   ${YELLOW}使用 pip 安装 HADDOCK3...${NC}"
if conda run -n "$ENV_NAME" python -m pip install -U haddock3 --quiet --disable-pip-version-check; then
    echo -e "${GREEN}✅ HADDOCK3 (pip) 安装成功${NC}"
else
    echo -e "${YELLOW}⚠️  HADDOCK3 安装失败，蛋白-蛋白对接功能将不可用${NC}"
    echo "   可稍后手动执行: conda run -n $ENV_NAME python -m pip install -U haddock3"
fi

# 4d. 验证 EC 分析依赖
echo -e "\n${BLUE}[4d]${NC} 验证 EC 分析依赖..."
EC_DEPS_OK=true

# 检查 PDB2PQR
if conda run -n "$ENV_NAME" python -c "import pdb2pqr" 2>/dev/null; then
    echo -e "${GREEN}✅ PDB2PQR Python API 可用${NC}"
else
    echo -e "${YELLOW}⚠️  PDB2PQR Python API 不可用${NC}"
    EC_DEPS_OK=false
fi

if command -v pdb2pqr &> /dev/null; then
    echo -e "${GREEN}✅ PDB2PQR CLI 可用${NC}"
else
    echo -e "${YELLOW}⚠️  PDB2PQR CLI 不可用${NC}"
fi

# 检查 APBS
if conda run -n "$ENV_NAME" python -c "from apbs_binary import run_apbs" 2>/dev/null; then
    echo -e "${GREEN}✅ APBS (apbs-binary) 可用${NC}"
elif conda run -n "$ENV_NAME" python -c "import apbs" 2>/dev/null; then
    echo -e "${GREEN}✅ APBS (Python API) 可用${NC}"
elif command -v apbs &> /dev/null; then
    echo -e "${GREEN}✅ APBS CLI 可用${NC}"
else
    echo -e "${YELLOW}⚠️  APBS 不可用，EC分析将无法运行${NC}"
    EC_DEPS_OK=false
fi

if [ "$EC_DEPS_OK" = false ]; then
    echo -e "${YELLOW}💡 提示: 可稍后运行 bash install_ec_dependencies.sh 来完成 EC 分析环境设置${NC}"
fi

# 5. 复制 GLINT 文件
echo -e "\n${BLUE}[5/6]${NC} 安装 GLINT 插件..."
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SOURCE_DIR="$SCRIPT_DIR/glint"

if [ ! -d "$SOURCE_DIR" ]; then
    echo -e "${RED}❌ 找不到 GLINT 源码目录: $SOURCE_DIR${NC}"
    exit 1
fi

if [ ! -f "$SOURCE_DIR/__init__.py" ]; then
    echo -e "${RED}❌ 找不到 GLINT 源码${NC}"
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

echo -e "${GREEN}✅ GLINT 文件安装完成${NC}"

# 5b. 创建 PyMOL 启动脚本
echo -e "\n${BLUE}[5b]${NC} 配置 PyMOL 启动脚本..."
PYMOL_STARTUP_DIR="$HOME/.pymol/startup"
mkdir -p "$PYMOL_STARTUP_DIR"

# 复制启动脚本
cp "$SCRIPT_DIR/pymol_startup_glint.py" "$PYMOL_STARTUP_DIR/01_glint.py"
echo -e "${GREEN}✅ PyMOL 启动脚本已配置${NC}"

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
    <string>com.vesper.glint</string>
    <key>CFBundleName</key>
    <string>GLINT</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleShortVersionString</key>
    <string>1.0.0</string>
</dict>
</plist>
PLIST_EOF


# 创建启动脚本（根据 PyMOL 类型）
cat > "$MACOS/launcher" << 'LAUNCHER_EOF'
#!/bin/bash
# GLINT Launcher (macOS .app)

# 日志文件（用于调试）
LOG_FILE="$HOME/.glint_launch.log"
echo "=== GLINT Launch $(date) ===" >> "$LOG_FILE"

# 查找 conda
CONDA_EXE=""
if command -v conda &> /dev/null; then
    CONDA_EXE=$(command -v conda)
else
    for p in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" "/opt/miniconda3/bin/conda" "/opt/anaconda3/bin/conda" "/usr/local/bin/conda" "/opt/homebrew/bin/conda" "/opt/homebrew/Caskroom/miniconda/base/bin/conda" "$HOME/opt/miniconda3/bin/conda"; do
        if [ -x "$p" ]; then
            CONDA_EXE="$p"
            echo "Found conda at: $CONDA_EXE" >> "$LOG_FILE"
            break
        fi
    done
fi

if [ -z "$CONDA_EXE" ]; then
    echo "ERROR: Conda not found" >> "$LOG_FILE"
    osascript -e 'display alert "GLINT Error" message "Conda not found.\n\nPlease install Miniconda first:\nhttps://docs.conda.io/en/latest/miniconda.html"'
    exit 1
fi

echo "Using conda: $CONDA_EXE" >> "$LOG_FILE"

# 查找 glint 环境
CONDA_BASE=$("$CONDA_EXE" info --base 2>/dev/null)
ENV_PATH="$CONDA_BASE/envs/glint"

if [ ! -d "$ENV_PATH" ]; then
    echo "ERROR: glint env not found at $ENV_PATH" >> "$LOG_FILE"
    osascript -e 'display alert "GLINT Error" message "Conda environment glint not found.\n\nPlease run the GLINT installer first."'
    exit 1
fi

echo "Using env: $ENV_PATH" >> "$LOG_FILE"

# 设置环境变量（避免 OpenMP 冲突）
export KMP_DUPLICATE_LIB_OK=TRUE
export OMP_NUM_THREADS=1
export QT_MAC_WANTS_LAYER=1

# 初始化 conda
eval "$($CONDA_EXE shell.bash hook)" 2>/dev/null

# 激活环境
conda activate glint 2>/dev/null || {
    echo "ERROR: Failed to activate conda env" >> "$LOG_FILE"
    osascript -e 'display alert "GLINT Error" message "Failed to activate conda environment glint.\n\nPlease check your conda installation."'
    exit 1
}

echo "Conda env activated, starting PyMOL..." >> "$LOG_FILE"

LAUNCHER_EOF

# 根据 PyMOL 类型添加启动命令
if [ "$PYMOL_TYPE" = "app" ]; then
    # 使用系统 PyMOL.app（需要注入 conda 环境路径）
    cat >> "$MACOS/launcher" << 'LAUNCHER_EOF'
# 启动系统 PyMOL.app 并加载 GLINT
# 注意：需要将 conda 环境的 site-packages 添加到 PYTHONPATH
CONDA_ENV_PYTHON="$ENV_PATH/lib/python3.10/site-packages"
if [ -d "$CONDA_ENV_PYTHON" ]; then
    export PYTHONPATH="$CONDA_ENV_PYTHON:$PYTHONPATH"
    echo "Added conda env to PYTHONPATH: $CONDA_ENV_PYTHON" >> "$LOG_FILE"
fi

/Applications/PyMOL.app/Contents/MacOS/PyMOL -d "
import sys, os
# 添加 conda 环境路径（确保能找到依赖）
conda_site_packages = os.path.expanduser('~/miniconda3/envs/glint/lib/python3.10/site-packages')
if os.path.exists(conda_site_packages) and conda_site_packages not in sys.path:
    sys.path.insert(0, conda_site_packages)
    print(f'Added conda env to sys.path: {conda_site_packages}')

# 添加 GLINT 启动路径
startup_path = os.path.expanduser('~/.pymol/startup')
if startup_path not in sys.path:
    sys.path.insert(0, startup_path)

try:
    import glint
    glint.glint_gui()
except Exception as e:
    print(f'GLINT Error: {e}')
    import traceback
    traceback.print_exc()
"

echo "PyMOL exited with code: $?" >> "$LOG_FILE"
LAUNCHER_EOF
else
    # 使用 conda 环境中的 PyMOL（推荐方式）
    cat >> "$MACOS/launcher" << 'LAUNCHER_EOF'
# 启动 conda PyMOL 并加载 GLINT
pymol -d "
import sys, os
startup_path = os.path.expanduser('~/.pymol/startup')
if startup_path not in sys.path:
    sys.path.insert(0, startup_path)
try:
    import glint
    glint.glint_gui()
except Exception as e:
    print(f'GLINT Error: {e}')
    import traceback
    traceback.print_exc()
"

echo "PyMOL exited with code: $?" >> "$LOG_FILE"
LAUNCHER_EOF
fi

chmod +x "$MACOS/launcher"

echo -e "${GREEN}✅ 桌面应用创建完成: $DESKTOP_APP${NC}"

# 完成
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}🎉 GLINT 安装成功！${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "启动方式："
echo "  1. 双击桌面上的 GLINT.app"
echo "  2. 或在 PyMOL 中运行:"
echo "     import glint"
echo "     glint.glint_gui()"
echo ""
echo "PyMOL 类型: $PYMOL_TYPE"
if [ "$PYMOL_TYPE" = "conda" ]; then
    echo -e "${GREEN}✅ 使用 conda 环境中的 PyMOL（推荐）${NC}"
    echo -e "${GREEN}   所有依赖已隔离，避免与系统 PyMOL 冲突${NC}"
    if [ "$HAS_SYSTEM_PYMOL" = true ]; then
        echo -e "${BLUE}💡 提示: 检测到系统 PyMOL.app，但 GLINT 使用独立的 conda PyMOL${NC}"
        echo -e "${BLUE}   这样可以确保依赖兼容性，您的系统 PyMOL 不会受影响${NC}"
    fi
else
    echo -e "${YELLOW}⚠️  使用系统 PyMOL.app（可能存在依赖冲突）${NC}"
    echo -e "${YELLOW}   如果遇到问题，建议重新安装并使用 conda PyMOL${NC}"
fi
echo ""
echo -e "${BLUE}📝 注意: 项目已从 GlueTK 重命名为 GLINT${NC}"
echo ""


