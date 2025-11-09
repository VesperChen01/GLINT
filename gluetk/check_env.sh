#!/bin/bash
# -*- coding: utf-8 -*-
# GlueTK 环境检测与自动配置脚本
# 支持 Windows (Git Bash/WSL)、macOS、Linux
# 检测系统 → 检测/安装 conda → 创建/检测环境 → 安装/检测依赖 → 验证 GUI 功能

# 不使用 set -e，手动处理错误以便更好的交互

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 日志函数
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# 环境名称和 Python 版本
ENV_NAME="gluetk"
PYTHON_VERSION="3.9"

# 必要的包列表
REQUIRED_PACKAGES=(
    "rdkit:RDKit (化学信息学)"
    "scipy:SciPy (科学计算)"
    "matplotlib:Matplotlib (绘图)"
    "PIL:Pillow (图像处理)"
    "numpy:NumPy (数值计算)"
    "PyQt5:PyQt5 (GUI框架)"
)

REQUIRED_COMMANDS=(
    "vina:AutoDock Vina (分子对接)"
    "obabel:Open Babel (格式转换)"
)

# 检测操作系统
detect_os() {
    log_info "检测操作系统..."
    
    # 检测 Windows
    if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" || "$OSTYPE" == "cygwin" ]]; then
        OS="Windows"
        ARCH="x86_64"
        CONDA_INSTALLER="Miniconda3-latest-Windows-x86_64.exe"
        IS_WINDOWS=true
        log_warning "检测到 Windows 系统，建议在 WSL 或 Git Bash 中运行"
    # 检测 macOS
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        OS="macOS"
        ARCH=$(uname -m)
        IS_WINDOWS=false
        if [[ "$ARCH" == "arm64" ]]; then
            CONDA_INSTALLER="Miniconda3-latest-MacOSX-arm64.sh"
        else
            CONDA_INSTALLER="Miniconda3-latest-MacOSX-x86_64.sh"
        fi
    # 检测 Linux
    elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
        OS="Linux"
        ARCH=$(uname -m)
        IS_WINDOWS=false
        if [[ "$ARCH" == "x86_64" ]]; then
            CONDA_INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
        elif [[ "$ARCH" == "aarch64" ]]; then
            CONDA_INSTALLER="Miniconda3-latest-Linux-aarch64.sh"
        else
            log_error "不支持的 Linux 架构: $ARCH"
            exit 1
        fi
    else
        log_error "不支持的操作系统: $OSTYPE"
        exit 1
    fi
    
    log_success "系统: $OS ($ARCH)"
}

# 检测 conda 是否已安装
check_conda() {
    log_info "检查 conda 是否已安装..."
    
    if command -v conda &> /dev/null; then
        CONDA_PATH=$(which conda)
        CONDA_VERSION=$(conda --version 2>&1)
        log_success "conda 已安装: $CONDA_VERSION"
        log_info "路径: $CONDA_PATH"
        return 0
    else
        log_warning "未检测到 conda"
        return 1
    fi
}

# 安装 Miniconda
install_miniconda() {
    log_info "开始安装 Miniconda..."
    
    INSTALL_DIR="$HOME/miniconda3"
    
    # 下载 Miniconda
    DOWNLOAD_URL="https://repo.anaconda.com/miniconda/$CONDA_INSTALLER"
    log_info "下载 Miniconda: $DOWNLOAD_URL"
    
    cd /tmp
    curl -O "$DOWNLOAD_URL"
    
    # 安装
    log_info "安装到: $INSTALL_DIR"
    bash "$CONDA_INSTALLER" -b -p "$INSTALL_DIR"
    
    # 清理安装文件
    rm "$CONDA_INSTALLER"
    
    # 初始化 conda
    log_info "初始化 conda..."
    source "$INSTALL_DIR/etc/profile.d/conda.sh"
    
    # 检测 shell 类型并初始化
    if [[ -n "$ZSH_VERSION" ]]; then
        SHELL_RC="$HOME/.zshrc"
        SHELL_NAME="zsh"
    elif [[ -n "$BASH_VERSION" ]]; then
        SHELL_RC="$HOME/.bashrc"
        SHELL_NAME="bash"
    else
        SHELL_RC="$HOME/.profile"
        SHELL_NAME="shell"
    fi
    
    "$INSTALL_DIR/bin/conda" init "$SHELL_NAME"
    
    log_success "Miniconda 安装完成！"
    log_warning "请运行以下命令激活 conda:"
    echo -e "${YELLOW}source $SHELL_RC${NC}"
    echo -e "或重启终端"
    
    # 临时激活当前会话
    export PATH="$INSTALL_DIR/bin:$PATH"
}

# 检测 conda 环境是否存在
check_conda_env() {
    log_info "检查 conda 环境 '$ENV_NAME'..."
    
    if conda env list | grep -q "^${ENV_NAME} \|^${ENV_NAME}$"; then
        log_success "环境 '$ENV_NAME' 已存在"
        return 0
    else
        log_warning "环境 '$ENV_NAME' 不存在"
        return 1
    fi
}

# 创建 GlueTK 环境
create_gluetk_env() {
    log_info "创建 GlueTK conda 环境..."
    
    # 检查环境是否已存在
    if check_conda_env; then
        read -p "是否删除并重新创建? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            log_info "删除旧环境..."
            conda env remove -n "$ENV_NAME" -y
        else
            log_info "使用现有环境"
            return 0
        fi
    fi
    
    log_info "创建环境: $ENV_NAME (Python $PYTHON_VERSION)"
    conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y
    
    log_success "环境创建完成"
}

# 激活 conda 环境
activate_env() {
    log_info "激活环境: $ENV_NAME"
    source "$(conda info --base)/etc/profile.d/conda.sh" 2>/dev/null || true
    conda activate "$ENV_NAME" 2>/dev/null || {
        log_error "无法激活环境 '$ENV_NAME'"
        return 1
    }
    log_success "环境已激活: $(python --version 2>&1)"
}

# 检测 Python 包是否安装
check_package() {
    local pkg_import="$1"
    local pkg_name="$2"
    
    if python -c "import $pkg_import" 2>/dev/null; then
        log_success "✓ $pkg_name"
        return 0
    else
        log_error "✗ $pkg_name (未安装)"
        return 1
    fi
}

# 检测命令行工具是否安装
check_command() {
    local cmd="$1"
    local name="$2"
    
    if command -v "$cmd" &> /dev/null; then
        local version=$("$cmd" --version 2>&1 | head -n 1 || echo "unknown")
        log_success "✓ $name: $version"
        return 0
    else
        log_error "✗ $name (未找到)"
        return 1
    fi
}

# 检测所有依赖
check_all_dependencies() {
    log_info "检测已安装的依赖..."
    
    activate_env || return 1
    
    local all_ok=true
    local missing_packages=()
    
    # 检测 Python 包
    log_info "
=== Python 包 ==="
    for pkg_info in "${REQUIRED_PACKAGES[@]}"; do
        IFS=':' read -r pkg_import pkg_name <<< "$pkg_info"
        if ! check_package "$pkg_import" "$pkg_name"; then
            all_ok=false
            missing_packages+=("$pkg_import")
        fi
    done
    
    # 检测命令行工具
    log_info "
=== 命令行工具 ==="
    for cmd_info in "${REQUIRED_COMMANDS[@]}"; do
        IFS=':' read -r cmd name <<< "$cmd_info"
        if ! check_command "$cmd" "$name"; then
            all_ok=false
        fi
    done
    
    if $all_ok; then
        log_success "
✅ 所有依赖已安装！"
        return 0
    else
        log_warning "
⚠️  部分依赖未安装"
        return 1
    fi
}

# 安装依赖包
install_dependencies() {
    log_info "安装依赖包..."
    
    activate_env || return 1
    
    # 更新 conda (静默模式)
    log_info "更新 conda..."
    conda update -n base -c defaults conda -y -q 2>/dev/null || true
    
    # 安装 conda 包
    log_info "安装 conda 包 (rdkit, scipy, matplotlib, pillow, numpy, pyqt, vina, openbabel)..."
    conda install -c conda-forge \
        rdkit \
        scipy \
        matplotlib \
        pillow \
        numpy \
        pyqt \
        autodock-vina \
        openbabel \
        -y -q
    
    # 安装 pip 包
    log_info "安装 pip 包..."
    pip install --upgrade pip -q
    pip install pymol-open-source -q 2>/dev/null || log_warning "PyMOL 安装失败 (可选依赖)"
    
    log_success "依赖包安装完成！"
}

# 测试 GUI 功能
test_gui_functionality() {
    log_info "测试 GUI 功能..."
    
    activate_env || return 1
    
    # 创建测试脚本
    local test_script="/tmp/molstruct_gui_test.py"
    cat > "$test_script" << 'EOF'
# -*- coding: utf-8 -*-
import sys

def test_gui():
    errors = []
    
    # 测试 Qt 框架
    try:
        from PyQt5.QtWidgets import QApplication, QMainWindow
        from PyQt5.QtCore import Qt
        print("✓ PyQt5 基础功能")
    except ImportError as e:
        errors.append(f"PyQt5: {e}")
    
    # 测试绘图功能
    try:
        import matplotlib
        matplotlib.use('Agg')  # 非交互式后端
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])
        plt.close(fig)
        print("✓ Matplotlib 绘图功能")
    except Exception as e:
        errors.append(f"Matplotlib: {e}")
    
    # 测试化学信息学功能
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
        mol = Chem.MolFromSmiles('CCO')
        if mol:
            AllChem.Compute2DCoords(mol)
            mw = Descriptors.MolWt(mol)
            print(f"✓ RDKit 化学计算 (乙醇 MW={mw:.2f})")
        else:
            errors.append("RDKit: 无法创建分子")
    except Exception as e:
        errors.append(f"RDKit: {e}")
    
    # 测试数值计算
    try:
        import numpy as np
        from scipy import spatial
        points = np.random.rand(10, 3)
        dist = spatial.distance.cdist(points, points)
        print(f"✓ NumPy/SciPy 距离计算 (10x10 矩阵)")
    except Exception as e:
        errors.append(f"NumPy/SciPy: {e}")
    
    # 测试图像处理
    try:
        from PIL import Image, ImageDraw
        img = Image.new('RGB', (100, 100), color='white')
        draw = ImageDraw.Draw(img)
        draw.rectangle([10, 10, 90, 90], outline='black')
        print("✓ Pillow 图像处理")
    except Exception as e:
        errors.append(f"Pillow: {e}")
    
    if errors:
        print("\n❌ 发现错误:")
        for err in errors:
            print(f"  - {err}")
        return False
    else:
        print("\n✅ 所有 GUI 功能测试通过！")
        return True

if __name__ == '__main__':
    success = test_gui()
    sys.exit(0 if success else 1)
EOF
    
    # 运行测试
    if python "$test_script"; then
        log_success "GUI 功能测试通过"
        rm -f "$test_script"
        return 0
    else
        log_error "GUI 功能测试失败"
        rm -f "$test_script"
        return 1
    fi
}

# 验证安装
verify_installation() {
    log_info "
======================================"
    log_info "验证环境配置"
    log_info "======================================"
    
    # 检测所有依赖
    if ! check_all_dependencies; then
        log_error "依赖检测失败，请先安装缺失的包"
        return 1
    fi
    
    # 测试 GUI 功能
    echo ""
    if ! test_gui_functionality; then
        log_warning "GUI 功能测试未完全通过"
    fi
    
    echo ""
    log_success "======================================"
    log_success "环境配置完成！"
    log_success "======================================"
    echo ""
    log_info "使用方法:"
    echo -e "  1. 激活环境: ${GREEN}conda activate $ENV_NAME${NC}"
    echo -e "  2. 运行插件: 在 PyMOL 中载入 gluetk"
    echo -e "  3. 打开 GUI: ${GREEN}molstruct_gui${NC} (在 PyMOL 命令行)"
    echo ""
    log_info "测试命令:"
    echo -e "  ${GREEN}vina --version${NC}"
    echo -e "  ${GREEN}python -c 'from PyQt5.QtWidgets import QApplication; print(\"Qt OK\")'${NC}"
}

# 主函数
main() {
    echo ""
    log_info "======================================"
    log_info "   GlueTK 环境检测与配置脚本"
    log_info "======================================"
    echo ""
    
    # 步骤 1: 检测操作系统
    log_info "[步骤 1/6] 检测操作系统"
    detect_os
    echo ""
    
    # 步骤 2: 检查/安装 conda
    log_info "[步骤 2/6] 检查 Conda 环境"
    if ! check_conda; then
        echo ""
        log_warning "未检测到 Conda，需要安装 Miniconda"
        read -p "是否立即安装? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            install_miniconda
            echo ""
            log_info "重新检测 conda..."
            if ! check_conda; then
                log_error "conda 安装后仍无法检测到，请手动执行:"
                echo "  source ~/.zshrc  # 或 source ~/.bashrc"
                echo "然后重新运行此脚本"
                exit 1
            fi
        else
            log_error "需要 conda 才能继续，退出"
            exit 1
        fi
    fi
    echo ""
    
    # 步骤 3: 检查/创建 conda 环境
    log_info "[步骤 3/6] 检查 Conda 环境 '$ENV_NAME'"
    if check_conda_env; then
        log_info "环境已存在，跳过创建"
    else
        log_warning "环境不存在，即将创建"
        create_gluetk_env
    fi
    echo ""
    
    # 步骤 4: 检查/安装依赖
    log_info "[步骤 4/6] 检查依赖包"
    if check_all_dependencies; then
        log_info "所有依赖已安装，跳过安装步骤"
    else
        echo ""
        log_warning "检测到缺失的依赖包"
        read -p "是否立即安装? (y/n): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            install_dependencies
        else
            log_warning "跳过依赖安装"
        fi
    fi
    echo ""
    
    # 步骤 5: 测试 GUI 功能
    log_info "[步骤 5/6] 测试 GUI 功能"
    if ! test_gui_functionality; then
        log_warning "GUI 功能测试未完全通过，但不影响基本使用"
    fi
    echo ""
    
    # 步骤 6: 最终验证
    log_info "[步骤 6/6] 最终验证"
    verify_installation
}

# 运行主函数
main
