# -*- coding: utf-8 -*-
"""
GlueTK 环境检测和配置模块
支持 Windows/macOS/Linux 环境检测、Conda 安装、依赖安装、GUI 功能测试
"""

import os
import sys
import platform
import subprocess
import tempfile
from typing import Dict, List, Tuple, Optional

# 环境配置
ENV_NAME = "gluetk"
PYTHON_VERSION = "3.9"

# 必需的 Python 包 (import_name, display_name, pip_name)
REQUIRED_PACKAGES = [
    ("rdkit", "RDKit", "rdkit"),
    ("scipy", "SciPy", "scipy"),
    ("matplotlib", "Matplotlib", "matplotlib"),
    ("PIL", "Pillow", "pillow"),
    ("numpy", "NumPy", "numpy"),
    ("PyQt5", "PyQt5", "pyqt5"),
]

# 必需的命令行工具
REQUIRED_COMMANDS = [
    ("vina", "AutoDock Vina"),
    ("obabel", "Open Babel"),
]


class EnvironmentChecker:
    """环境检测和配置类"""
    
    def __init__(self, log_callback=None, auto_install=False):
        """
        Args:
            log_callback: 日志回调函数，用于输出消息到 GUI
            auto_install: 是否自动安装缺失的依赖
        """
        self.log_callback = log_callback or print
        self.auto_install = auto_install
        self.os_type = None
        self.os_arch = None
        self.conda_path = None
        
    def log(self, msg: str):
        """输出日志"""
        if self.log_callback:
            self.log_callback(msg)
    
    # ========== 系统检测 ==========
    
    def detect_os(self) -> Tuple[str, str]:
        """
        检测操作系统和架构
        
        Returns:
            (os_type, arch): 例如 ("macOS", "arm64")
        """
        system = platform.system()
        machine = platform.machine()
        
        if system == "Darwin":
            self.os_type = "macOS"
        elif system == "Linux":
            self.os_type = "Linux"
        elif system == "Windows":
            self.os_type = "Windows"
        else:
            self.os_type = "Unknown"
        
        self.os_arch = machine
        
        self.log(f"✓ 系统: {self.os_type} ({self.os_arch})")
        return self.os_type, self.os_arch
    
    # ========== Conda 检测 ==========
    
    def check_conda(self) -> bool:
        """
        检查 conda 是否已安装
        
        Returns:
            bool: True 如果 conda 可用
        """
        try:
            result = subprocess.run(
                ["conda", "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                version = result.stdout.strip()
                self.log(f"✓ Conda: {version}")
                self.conda_path = subprocess.run(
                    ["which", "conda"],
                    capture_output=True,
                    text=True
                ).stdout.strip()
                return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
        
        self.log("✗ Conda: 未检测到")
        return False
    
    def check_conda_env(self, env_name: str = ENV_NAME) -> bool:
        """
        检查指定的 conda 环境是否存在
        
        Args:
            env_name: 环境名称
            
        Returns:
            bool: True 如果环境存在
        """
        if not self.check_conda():
            return False
        
        try:
            result = subprocess.run(
                ["conda", "env", "list"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                for line in result.stdout.split('\n'):
                    if line.strip().startswith(env_name + " ") or line.strip().startswith(env_name + "\t"):
                        self.log(f"✓ 环境 '{env_name}' 已存在")
                        return True
        except subprocess.TimeoutExpired:
            pass
        
        self.log(f"✗ 环境 '{env_name}' 不存在")
        return False
    
    # ========== 依赖检测 ==========
    
    def check_package(self, import_name: str, display_name: str) -> bool:
        """
        检查 Python 包是否已安装
        
        Args:
            import_name: 导入名称
            display_name: 显示名称
            
        Returns:
            bool: True 如果包已安装
        """
        try:
            __import__(import_name)
            self.log(f"  ✓ {display_name}")
            return True
        except ImportError:
            self.log(f"  ✗ {display_name} (未安装)")
            return False
    
    def check_command(self, cmd: str, name: str) -> bool:
        """
        检查命令行工具是否可用
        
        Args:
            cmd: 命令名称
            name: 显示名称
            
        Returns:
            bool: True 如果命令可用
        """
        try:
            result = subprocess.run(
                [cmd, "--version"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                version = result.stdout.split('\n')[0]
                self.log(f"  ✓ {name}: {version}")
                return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            pass
        
        self.log(f"  ✗ {name} (未找到)")
        return False
    
    def check_all_dependencies(self) -> Dict[str, bool]:
        """
        检查所有依赖项
        
        Returns:
            Dict[str, bool]: {依赖名: 是否已安装}
        """
        self.log("\n📦 检查 Python 包:")
        
        status = {}
        all_ok = True
        
        for import_name, display_name, _ in REQUIRED_PACKAGES:
            available = self.check_package(import_name, display_name)
            status[display_name] = available
            if not available:
                all_ok = False
        
        self.log("\n🛠️  检查命令行工具:")
        
        for cmd, name in REQUIRED_COMMANDS:
            available = self.check_command(cmd, name)
            status[name] = available
            if not available:
                all_ok = False
        
        if all_ok:
            self.log("\n✅ 所有依赖已安装")
        else:
            missing = [k for k, v in status.items() if not v]
            self.log(f"\n⚠️  缺失依赖: {', '.join(missing)}")
        
        return status
    
    # ========== GUI 功能测试 ==========
    
    def test_gui_functionality(self) -> bool:
        """
        测试 GUI 相关功能是否正常
        
        Returns:
            bool: True 如果所有测试通过
        """
        self.log("\n🧪 测试 GUI 功能:")
        
        errors = []
        
        # 测试 PyQt5
        try:
            from PyQt5.QtWidgets import QApplication
            from PyQt5.QtCore import Qt
            self.log("  ✓ PyQt5 基础功能")
        except ImportError as e:
            errors.append(f"PyQt5: {e}")
            self.log(f"  ✗ PyQt5: {e}")
        
        # 测试 Matplotlib
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.plot([1, 2, 3], [1, 2, 3])
            plt.close(fig)
            self.log("  ✓ Matplotlib 绘图功能")
        except Exception as e:
            errors.append(f"Matplotlib: {e}")
            self.log(f"  ✗ Matplotlib: {e}")
        
        # 测试 RDKit
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, Descriptors
            mol = Chem.MolFromSmiles('CCO')
            if mol:
                AllChem.Compute2DCoords(mol)
                mw = Descriptors.MolWt(mol)
                self.log(f"  ✓ RDKit 化学计算 (乙醇 MW={mw:.2f})")
            else:
                errors.append("RDKit: 无法创建分子")
                self.log("  ✗ RDKit: 无法创建分子")
        except Exception as e:
            errors.append(f"RDKit: {e}")
            self.log(f"  ✗ RDKit: {e}")
        
        # 测试 NumPy/SciPy
        try:
            import numpy as np
            from scipy import spatial
            points = np.random.rand(10, 3)
            dist = spatial.distance.cdist(points, points)
            self.log(f"  ✓ NumPy/SciPy 距离计算")
        except Exception as e:
            errors.append(f"NumPy/SciPy: {e}")
            self.log(f"  ✗ NumPy/SciPy: {e}")
        
        # 测试 Pillow
        try:
            from PIL import Image, ImageDraw
            img = Image.new('RGB', (100, 100), color='white')
            draw = ImageDraw.Draw(img)
            draw.rectangle([10, 10, 90, 90], outline='black')
            self.log("  ✓ Pillow 图像处理")
        except Exception as e:
            errors.append(f"Pillow: {e}")
            self.log(f"  ✗ Pillow: {e}")
        
        if errors:
            self.log(f"\n⚠️  部分功能测试失败")
            return False
        else:
            self.log(f"\n✅ 所有 GUI 功能测试通过")
            return True
    
    # ========== 安装指南 ==========
    
    def get_conda_install_url(self) -> str:
        """
        获取 Miniconda 下载链接
        
        Returns:
            str: 下载 URL
        """
        if not self.os_type or not self.os_arch:
            self.detect_os()
        
        base_url = "https://repo.anaconda.com/miniconda/"
        
        if self.os_type == "macOS":
            if self.os_arch == "arm64":
                return f"{base_url}Miniconda3-latest-MacOSX-arm64.sh"
            else:
                return f"{base_url}Miniconda3-latest-MacOSX-x86_64.sh"
        elif self.os_type == "Linux":
            if self.os_arch == "x86_64":
                return f"{base_url}Miniconda3-latest-Linux-x86_64.sh"
            elif self.os_arch == "aarch64":
                return f"{base_url}Miniconda3-latest-Linux-aarch64.sh"
        elif self.os_type == "Windows":
            return f"{base_url}Miniconda3-latest-Windows-x86_64.exe"
        
        return base_url
    
    # ========== 自动安装 ==========
    
    def auto_install_conda_env(self) -> bool:
        """
        自动创建 conda 环境
        
        Returns:
            bool: True 如果成功
        """
        if not self.check_conda():
            self.log("✗ 无法自动创建环境: Conda 未安装")
            return False
        
        if self.check_conda_env(ENV_NAME):
            self.log(f"✓ 环境 '{ENV_NAME}' 已存在，跳过创建")
            return True
        
        self.log(f"\n🔧 正在创建 Conda 环境 '{ENV_NAME}'...")
        
        try:
            result = subprocess.run(
                ["conda", "create", "-n", ENV_NAME, f"python={PYTHON_VERSION}", "-y"],
                capture_output=True,
                text=True,
                timeout=300
            )
            
            if result.returncode == 0:
                self.log(f"✅ 环境 '{ENV_NAME}' 创建成功")
                return True
            else:
                self.log(f"✗ 环境创建失败: {result.stderr}")
                return False
        except subprocess.TimeoutExpired:
            self.log("✗ 环境创建超时")
            return False
        except Exception as e:
            self.log(f"✗ 环境创建出错: {e}")
            return False
    
    def auto_install_dependencies(self) -> Dict[str, bool]:
        """
        自动安装缺失的依赖
        
        Returns:
            Dict[str, bool]: {依赖名: 是否安装成功}
        """
        if not self.check_conda():
            self.log("✗ 无法自动安装: Conda 未安装")
            return {}
        
        status = self.check_all_dependencies()
        missing = [name for name, avail in status.items() if not avail]
        
        if not missing:
            self.log("✅ 所有依赖已安装，无需操作")
            return status
        
        self.log(f"\n📦 开始自动安装缺失的依赖...")
        
        install_results = {}
        
        # 区分 Python 包和命令行工具
        py_packages = [p for p in missing if p not in ["AutoDock Vina", "Open Babel"]]
        cmd_tools = [p for p in missing if p in ["AutoDock Vina", "Open Babel"]]
        
        # 批量安装 Python 包
        if py_packages:
            self.log(f"\n  安装 Python 包: {', '.join(py_packages)}")
            packages_to_install = []
            
            for _, display_name, pip_name in REQUIRED_PACKAGES:
                if display_name in py_packages:
                    packages_to_install.append(pip_name)
            
            if packages_to_install:
                try:
                    cmd = [
                        "conda", "install", "-c", "conda-forge",
                        "-y"
                    ] + packages_to_install
                    
                    self.log(f"  执行: {' '.join(cmd)}")
                    
                    result = subprocess.run(
                        cmd,
                        capture_output=True,
                        text=True,
                        timeout=600
                    )
                    
                    if result.returncode == 0:
                        for pkg in py_packages:
                            install_results[pkg] = True
                            self.log(f"  ✅ {pkg} 安装成功")
                    else:
                        for pkg in py_packages:
                            install_results[pkg] = False
                            self.log(f"  ✗ {pkg} 安装失败")
                        self.log(f"  错误信息: {result.stderr[:200]}")
                except subprocess.TimeoutExpired:
                    for pkg in py_packages:
                        install_results[pkg] = False
                        self.log(f"  ✗ {pkg} 安装超时")
                except Exception as e:
                    for pkg in py_packages:
                        install_results[pkg] = False
                        self.log(f"  ✗ {pkg} 安装出错: {e}")
        
        # 安装命令行工具
        if "AutoDock Vina" in cmd_tools:
            self.log(f"\n  安装 AutoDock Vina...")
            success = self._install_conda_package("autodock-vina")
            install_results["AutoDock Vina"] = success
            if success:
                self.log("  ✅ AutoDock Vina 安装成功")
            else:
                self.log("  ✗ AutoDock Vina 安装失败")
        
        if "Open Babel" in cmd_tools:
            self.log(f"\n  安装 Open Babel...")
            success = self._install_conda_package("openbabel")
            install_results["Open Babel"] = success
            if success:
                self.log("  ✅ Open Babel 安装成功")
            else:
                self.log("  ✗ Open Babel 安装失败")
        
        # 总结
        success_count = sum(1 for v in install_results.values() if v)
        total_count = len(install_results)
        
        if success_count == total_count:
            self.log(f"\n✅ 所有依赖安装完成 ({success_count}/{total_count})")
        else:
            self.log(f"\n⚠️  部分依赖安装失败 ({success_count}/{total_count})")
        
        return install_results
    
    def _install_conda_package(self, package_name: str) -> bool:
        """
        安装单个 conda 包
        
        Args:
            package_name: 包名
            
        Returns:
            bool: True 如果成功
        """
        try:
            result = subprocess.run(
                ["conda", "install", "-c", "conda-forge", package_name, "-y"],
                capture_output=True,
                text=True,
                timeout=300
            )
            return result.returncode == 0
        except:
            return False
    
    def get_install_instructions(self) -> Dict[str, str]:
        """
        获取安装说明
        
        Returns:
            Dict[str, str]: 安装说明字典
        """
        status = self.check_all_dependencies()
        missing_packages = [name for name, avail in status.items() if not avail]
        
        instructions = {}
        
        # Conda 安装
        if not self.check_conda():
            instructions["conda"] = f"""
📥 安装 Miniconda:

1. 下载: {self.get_conda_install_url()}
2. 安装后运行: source ~/.zshrc  (或 ~/.bashrc)
3. 验证: conda --version
"""
        
        # 环境创建
        if not self.check_conda_env():
            instructions["environment"] = f"""
🔧 创建 Conda 环境:

conda create -n {ENV_NAME} python={PYTHON_VERSION} -y
conda activate {ENV_NAME}
"""
        
        # 依赖安装
        if missing_packages:
            # 区分 Python 包和命令行工具
            py_packages = [p for p in missing_packages if p not in ["AutoDock Vina", "Open Babel"]]
            cmd_tools = [p for p in missing_packages if p in ["AutoDock Vina", "Open Babel"]]
            
            install_cmd = []
            
            if py_packages:
                install_cmd.append("# 安装 Python 包")
                install_cmd.append("conda install -c conda-forge rdkit scipy matplotlib pillow numpy pyqt -y")
            
            if cmd_tools:
                install_cmd.append("\n# 安装命令行工具")
                if "AutoDock Vina" in cmd_tools:
                    install_cmd.append("conda install -c conda-forge autodock-vina -y")
                if "Open Babel" in cmd_tools:
                    install_cmd.append("conda install -c conda-forge openbabel -y")
            
            instructions["dependencies"] = f"""
📦 安装依赖包:

{chr(10).join(install_cmd)}
"""
        
        # 使用说明
        instructions["usage"] = f"""
✅ 使用方法:

1. 激活环境:
   conda activate {ENV_NAME}

2. 启动 PyMOL (在激活的环境中)

3. 加载插件:
   run /path/to/gluetk/__init__.py

4. 打开 GUI:
   molstruct_gui
"""
        
        return instructions
    
    # ========== 完整检查 ==========
    
    def run_full_check(self) -> Dict[str, any]:
        """
        运行完整的环境检查
        
        Returns:
            Dict: 检查结果摘要
        """
        self.log("=" * 60)
        self.log("⚡ GlueTK 环境检查")
        self.log("=" * 60)
        
        # 1. 检测操作系统
        self.log("\n🖥️  系统信息:")
        os_type, os_arch = self.detect_os()
        
        # 2. 检查 Conda
        self.log("\n🐍 Conda:")
        has_conda = self.check_conda()
        
        # 3. 检查环境
        if has_conda:
            self.log(f"\n📁 Conda 环境:")
            has_env = self.check_conda_env(ENV_NAME)
        else:
            has_env = False
        
        # 4. 检查依赖
        dep_status = self.check_all_dependencies()
        
        # 5. 自动安装（如果启用）
        if self.auto_install and has_conda:
            self.log("\n" + "=" * 60)
            self.log("🚀 自动安装模式")
            self.log("=" * 60)
            
            # 创建环境（如果不存在）
            if not has_env:
                has_env = self.auto_install_conda_env()
            
            # 安装依赖
            install_results = self.auto_install_dependencies()
            
            # 重新检查依赖状态
            self.log("\n🔍 重新检查依赖状态...")
            dep_status = self.check_all_dependencies()
        
        # 6. 测试 GUI 功能
        gui_ok = self.test_gui_functionality()
        
        # 7. 生成安装说明或总结
        self.log("\n" + "=" * 60)
        
        all_ok = has_conda and has_env and all(dep_status.values()) and gui_ok
        
        if all_ok:
            self.log("✅ 环境配置完美！所有功能可用")
        else:
            if self.auto_install:
                self.log("⚠️  部分依赖安装失败或需要手动配置")
            else:
                self.log("⚠️  环境需要配置")
            instructions = self.get_install_instructions()
            for key, text in instructions.items():
                self.log(f"\n{text}")
        
        self.log("=" * 60)
        
        return {
            "os_type": os_type,
            "os_arch": os_arch,
            "has_conda": has_conda,
            "has_env": has_env,
            "dependencies": dep_status,
            "gui_ok": gui_ok,
            "all_ok": all_ok
        }


# ========== 便捷函数 ==========

def check_environment(log_callback=None, auto_install=False) -> Dict[str, any]:
    """
    检查环境（便捷函数）
    
    Args:
        log_callback: 日志回调函数
        auto_install: 是否自动安装缺失的依赖
        
    Returns:
        Dict: 检查结果
    """
    checker = EnvironmentChecker(log_callback, auto_install=auto_install)
    return checker.run_full_check()


def get_dependency_status() -> Dict[str, bool]:
    """
    获取依赖状态（便捷函数）
    
    Returns:
        Dict[str, bool]: {依赖名: 是否可用}
    """
    checker = EnvironmentChecker(log_callback=None)
    return checker.check_all_dependencies()


# ========== 命令行接口 ==========

if __name__ == "__main__":
    # 命令行模式
    import argparse
    
    parser = argparse.ArgumentParser(description="GlueTK 环境检查和自动配置工具")
    parser.add_argument(
        "--auto-install",
        action="store_true",
        help="自动安装缺失的依赖（需要 Conda）"
    )
    
    args = parser.parse_args()
    
    checker = EnvironmentChecker(auto_install=args.auto_install)
    result = checker.run_full_check()
    
    sys.exit(0 if result["all_ok"] else 1)
