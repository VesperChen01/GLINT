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

# 安装 conda 依赖（不包含 Vina，本身可选，且在 osx-arm64 上可能没有官方包）
echo ""
echo "📦 Installing conda dependencies (core stack, without Vina)..."
conda install -n ${ENV_NAME} -c conda-forge \
    rdkit scipy matplotlib pillow numpy pandas seaborn \
    pyqt openbabel pymol-open-source -y

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

# Source Conda
source "${CONDA_BASE_PATH}/etc/profile.d/conda.sh"

# 激活环境
conda activate ${ENV_NAME}

# 插件路径
PLUGIN_PATH="${SCRIPT_DIR}/__init__.py"

echo "Starting GlueTK..."
# 启动 PyMOL 并加载插件
# -q: quiet launch (optional)
pymol "\${PLUGIN_PATH}"
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
