#!/bin/bash
# GlueTK 一键安装脚本
# 用法: bash install.sh

set -e

ENV_NAME="gluetk"
PYTHON_VERSION="3.9"

echo "============================================================"
echo "🧬 GlueTK - One-Click Installation"
echo "============================================================"

# 检查 conda
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install Miniconda first:"
    echo "   https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "✓ Conda found"

# 检查环境是否已存在
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "⚠️  Environment '${ENV_NAME}' already exists."
    read -p "   Remove and recreate? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n ${ENV_NAME} -y
    else
        echo "   Skipping environment creation..."
    fi
fi

# 创建环境（如果不存在）
if ! conda env list | grep -q "^${ENV_NAME} "; then
    echo ""
    echo "📦 Creating conda environment '${ENV_NAME}'..."
    conda create -n ${ENV_NAME} python=${PYTHON_VERSION} -y
fi

# 安装 conda 依赖
echo ""
echo "📦 Installing conda dependencies..."
echo "   (This might take a few minutes)"

# 1. 尝试安装核心计算库 + pymol-open-source
conda install -n ${ENV_NAME} -c conda-forge \
    rdkit scipy matplotlib pillow numpy pandas seaborn \
    pyqt openbabel pymol-open-source -y

# 2. 检查 PyMOL 是否安装成功
PYMOL_EXECUTABLE="pymol"
USE_SYSTEM_PYMOL=false

if ! conda run -n ${ENV_NAME} command -v pymol &> /dev/null; then
    echo "⚠️  Conda 'pymol' not found. Checking for system PyMOL..."
    
    if [ -d "/Applications/PyMOL.app" ]; then
        echo "✓ Found system PyMOL at /Applications/PyMOL.app"
        PYMOL_EXECUTABLE="/Applications/PyMOL.app/Contents/bin/pymol"
        USE_SYSTEM_PYMOL=true
    else
        echo "❌ Error: No PyMOL found (neither in conda nor /Applications)."
        echo "   Please install PyMOL manually."
        # 不退出，允许只安装环境
    fi
else
    echo "✓ Conda PyMOL installed successfully."
fi

# 激活环境
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${ENV_NAME}

# 使用 pip 尝试安装 Vina + Meeko（可选，如失败仅给出警告，不中断安装）
echo ""
echo "📦 Installing optional docking stack via pip (vina + meeko)..."
set +e
pip install vina meeko --quiet --disable-pip-version-check
PIP_STATUS=$?
set -e
if [[ ${PIP_STATUS} -ne 0 ]]; then
  echo "⚠️  pip install vina/meeko failed (optional). Docking features may be unavailable."
  echo "   You can try manually inside the environment:"
  echo "   conda activate ${ENV_NAME} && pip install vina meeko"
else
  echo "✅ Vina + Meeko installed via pip"
fi

# 创建插件符号链接 (作为备份加载方式)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p ~/.pymol/startup
ln -sf "${SCRIPT_DIR}" ~/.pymol/startup/gluetk
echo "✅ Plugin installed to ~/.pymol/startup/gluetk"

# 创建 MacOS APP 启动器 (带图标)
APP_NAME="GlueTK"
APP_DIR="$HOME/Desktop/${APP_NAME}.app"
CONTENTS_DIR="${APP_DIR}/Contents"
MACOS_DIR="${CONTENTS_DIR}/MacOS"
RESOURCES_DIR="${CONTENTS_DIR}/Resources"

echo ""
echo "📦 Creating ${APP_NAME}.app on Desktop..."

# 清理旧版
rm -rf "${APP_DIR}" "$HOME/Desktop/Start_GlueTK.command"

# 创建目录结构
mkdir -p "${MACOS_DIR}"
mkdir -p "${RESOURCES_DIR}"

# 复制图标
ICON_SOURCE="${SCRIPT_DIR}/assets/AppIcon.icns"
if [ -f "${ICON_SOURCE}" ]; then
    cp "${ICON_SOURCE}" "${RESOURCES_DIR}/AppIcon.icns"
else
    echo "⚠️  Icon file not found at ${ICON_SOURCE}"
fi

# 创建 Info.plist
cat > "${CONTENTS_DIR}/Info.plist" << EOF
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
    <string>1.0</string>
</dict>
</plist>
EOF

# 获取 Conda Base 路径 (用于硬编码到启动脚本)
CONDA_BASE_PATH=$(conda info --base)

# 创建启动脚本
LAUNCHER_SCRIPT="${MACOS_DIR}/launcher"
cat > "${LAUNCHER_SCRIPT}" << EOF
#!/bin/bash
# GlueTK Launcher

# 1. Source Conda to get dependencies (rdkit, etc)
source "${CONDA_BASE_PATH}/etc/profile.d/conda.sh"
conda activate ${ENV_NAME}

# 2. Define Paths
PLUGIN_PATH="${SCRIPT_DIR}/__init__.py"
PYMOL_CMD="${PYMOL_EXECUTABLE}"

echo "Starting GlueTK..."
echo "PyMOL: \${PYMOL_CMD}"
echo "Plugin: \${PLUGIN_PATH}"

# 3. Launch
if [[ "\${PYMOL_CMD}" == *"/Applications/PyMOL.app"* ]]; then
    # 如果使用的是系统 PyMOL，尝试注入 PYTHONPATH 以加载 Conda 里的库
    # 注意：如果 Python 版本差异过大，这可能会失败，但在 Mac 上这是唯一让系统 PyMOL 用 Conda 库的办法
    export PYTHONPATH="\$CONDA_PREFIX/lib/python${PYTHON_VERSION}/site-packages:\$PYTHONPATH"
    "\${PYMOL_CMD}" "\${PLUGIN_PATH}"
else
    # 使用 Conda 内部 PyMOL
    pymol "\${PLUGIN_PATH}"
fi
EOF

chmod +x "${LAUNCHER_SCRIPT}"

echo ""
echo "============================================================"
echo "✅ Installation complete!"
echo "============================================================"
echo ""
echo "📌 Created App: ${APP_NAME}.app"
echo "   (Double-click it on Desktop to start GlueTK)"
echo ""

# 询问是否立即启动
read -p "🚀 Launch GlueTK now? [Y/n] " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Nn]$ ]]; then
    open "${APP_DIR}"
fi
echo "============================================================"
