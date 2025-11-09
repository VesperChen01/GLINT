# -*- coding: utf-8 -*-
"""
env_setup.py
环境依赖自动检测与静默安装模块

用于确保 GlueTK 插件所需的所有依赖包在运行前已正确安装
"""

from __future__ import print_function
import sys
import subprocess
import importlib
from typing import List, Tuple

# 必需依赖包（用于高质量分析模式）
REQUIRED_PACKAGES = [
    ("rdkit", "rdkit"),  # (import_name, pip_name)
    ("scipy", "scipy"),
    ("matplotlib", "matplotlib"),
    ("PIL", "pillow"),
    ("numpy", "numpy"),
]

# 可选依赖（增强功能）
OPTIONAL_PACKAGES = [
    ("openbabel", "openbabel-wheel"),
]


def _check_package(import_name: str) -> bool:
    """检查包是否可导入"""
    try:
        importlib.import_module(import_name)
        return True
    except ImportError:
        return False


def _install_package(pip_name: str, silent: bool = True) -> bool:
    """
    静默安装包
    
    参数:
        pip_name: pip 包名
        silent: 是否静默安装（不显示安装输出）
    
    返回:
        bool: 安装是否成功
    """
    try:
        # 使用当前 Python 解释器的 pip 进行安装
        cmd = [sys.executable, "-m", "pip", "install", pip_name, "--quiet", "--disable-pip-version-check"]
        
        if silent:
            result = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=300  # 5分钟超时
            )
        else:
            result = subprocess.run(cmd, timeout=300)
        
        return result.returncode == 0
    except Exception as e:
        print(f"[env_setup] Failed to install {pip_name}: {e}")
        return False


def check_and_install_dependencies(silent: bool = True, required_only: bool = True) -> Tuple[bool, List[str]]:
    """
    检查并安装所有必需依赖
    
    参数:
        silent: 是否静默安装（不显示详细输出）
        required_only: 是否只安装必需依赖（不安装可选依赖）
    
    返回:
        (all_ok, missing_packages): 所有依赖是否满足，缺失的包列表
    """
    missing_packages = []
    packages_to_check = REQUIRED_PACKAGES if required_only else REQUIRED_PACKAGES + OPTIONAL_PACKAGES
    
    print("[env_setup] 🔍 Checking dependencies...")
    
    for import_name, pip_name in packages_to_check:
        if not _check_package(import_name):
            missing_packages.append(pip_name)
            print(f"[env_setup] ⚠️  Missing dependency: {pip_name}")
    
    if not missing_packages:
        print("[env_setup] ✅ All dependencies satisfied")
        return True, []
    
    # 开始安装缺失的包
    print(f"[env_setup] 📦 Installing {len(missing_packages)} missing package(s)...")
    
    failed_packages = []
    for pip_name in missing_packages:
        print(f"[env_setup] ⏳ Installing {pip_name}...", end=" ")
        if _install_package(pip_name, silent=silent):
            print("✓")
        else:
            print("✗")
            failed_packages.append(pip_name)
    
    if failed_packages:
        print(f"[env_setup] ⚠️  Failed to install: {', '.join(failed_packages)}")
        print("[env_setup] 💡 Please install manually: pip install " + " ".join(failed_packages))
        return False, failed_packages
    
    print("[env_setup] ✅ All dependencies installed")
    return True, []


def ensure_dependencies() -> bool:
    """
    确保所有依赖已安装（入口函数）
    
    返回:
        bool: 是否所有依赖都已满足
    """
    success, failed = check_and_install_dependencies(silent=True, required_only=True)
    return success


def get_dependency_status() -> dict:
    """
    获取依赖包状态
    
    返回:
        dict: {package_name: is_available}
    """
    status = {}
    all_packages = REQUIRED_PACKAGES + OPTIONAL_PACKAGES
    
    for import_name, pip_name in all_packages:
        status[pip_name] = _check_package(import_name)
    
    return status


def print_dependency_report():
    """Print dependency status report"""
    print("\n" + "="*50)
    print("GlueTK Dependency Status Report")
    print("="*50)
    
    status = get_dependency_status()
    
    print("\n[Required]")
    for import_name, pip_name in REQUIRED_PACKAGES:
        available = status.get(pip_name, False)
        symbol = "✓" if available else "✗"
        print(f"  {symbol} {pip_name:20} {'installed' if available else 'not installed'}")
    
    print("\n[Optional]")
    for import_name, pip_name in OPTIONAL_PACKAGES:
        available = status.get(pip_name, False)
        symbol = "✓" if available else "✗"
        print(f"  {symbol} {pip_name:20} {'installed' if available else 'not installed'}")
    
    print("="*50 + "\n")


if __name__ == "__main__":
    # 命令行测试
    print_dependency_report()
    ensure_dependencies()
    print("\nStatus after check:")
    print_dependency_report()
