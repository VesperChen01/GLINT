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

# 第一次运行标记文件
_FIRST_RUN_MARKER = os.path.join(os.path.expanduser("~"), ".gluetk_initialized")

# 必需的 Python 包 (import_name, display_name, pip_name)
REQUIRED_PACKAGES = [
    ("rdkit", "RDKit", "rdkit"),
    ("scipy", "SciPy", "scipy"),
    ("matplotlib", "Matplotlib", "matplotlib"),
    ("PIL", "Pillow", "pillow"),
    ("numpy", "NumPy", "numpy"),
    ("pandas", "Pandas", "pandas"),
    ("seaborn", "Seaborn", "seaborn"),
    ("PyQt5", "PyQt5", "pyqt5"),
    ("open3d", "Open3D", "open3d"),  # 表面分析必需
    ("skimage", "scikit-image", "scikit-image"),  # Marching Cubes 算法
]

# 可选但推荐的 Python 包 (import_name, display_name, pip_name, description)
OPTIONAL_PACKAGES = [
    ("pyhmmer", "pyhmmer", "pyhmmer", "C2H2 锌指 HMM 检测（提高精度）"),
    ("haddock", "HADDOCK3", "haddock3", "Protein-Protein Docking Engine"),
    ("trimesh", "trimesh", "trimesh", "轻量级网格处理库"),
]

# EC 分析必需的 Python 包 (import_name, display_name, pip_name, description)
EC_REQUIRED_PACKAGES = [
    ("pdb2pqr", "PDB2PQR", "pdb2pqr", "蛋白结构准备（EC分析必需）"),
    ("apbs", "APBS Python", "apbs", "泊松-玻尔兹曼求解器 Python API"),
]

# 外部工具配置 (cmd_name, display_name, description, install_info)
EXTERNAL_TOOLS = {
    "msms": {
        "display_name": "MSMS",
        "description": "分子表面生成工具（最精确）",
        "install_info": {
            "macOS": "brew install brewsci/bio/msms 或从 https://ccsb.scripps.edu/msms/ 下载",
            "Linux": "从 https://ccsb.scripps.edu/msms/ 下载并添加到 PATH",
            "Windows": "从 https://ccsb.scripps.edu/msms/ 下载 Windows 版本",
        },
        "search_paths": {
            "macOS": ["/usr/local/bin/msms", "/opt/homebrew/bin/msms", "~/bin/msms"],
            "Linux": ["/usr/bin/msms", "/usr/local/bin/msms", "~/bin/msms"],
            "Windows": [r"C:\Program Files\MSMS\msms.exe", r"C:\msms\msms.exe"],
        }
    },
    "apbs": {
        "display_name": "APBS",
        "description": "自适应泊松-玻尔兹曼求解器（精确静电势，EC分析必需）",
        "install_info": {
            "macOS": "pip install apbs 或 brew install brewsci/bio/apbs",
            "Linux": "pip install apbs 或 apt install apbs",
            "Windows": "pip install apbs 或从 https://www.poissonboltzmann.org/ 下载",
        },
        "search_paths": {
            "macOS": ["/usr/local/bin/apbs", "/opt/homebrew/bin/apbs", "~/bin/apbs"],
            "Linux": ["/usr/bin/apbs", "/usr/local/bin/apbs", "~/bin/apbs"],
            "Windows": [r"C:\Program Files\APBS\apbs.exe", r"C:\APBS\apbs.exe"],
        },
        "python_module": "apbs",  # 也可以通过 Python API 使用
    },
    "pdb2pqr": {
        "display_name": "PDB2PQR",
        "description": "PDB 到 PQR 格式转换（EC分析必需）",
        "install_info": {
            "macOS": "pip install pdb2pqr",
            "Linux": "pip install pdb2pqr",
            "Windows": "pip install pdb2pqr",
        },
        "search_paths": {
            "macOS": ["/usr/local/bin/pdb2pqr", "/opt/homebrew/bin/pdb2pqr"],
            "Linux": ["/usr/bin/pdb2pqr", "/usr/local/bin/pdb2pqr"],
            "Windows": [],
        },
        "python_module": "pdb2pqr",  # 推荐通过 Python API 使用
    },
}

# Pfam HMM 文件配置
HMM_FILES = {
    "zf-C2H2": {
        "pfam_id": "PF00096",
        "url": "https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/PF00096?annotation=hmm",
        "local_name": "zf-C2H2.hmm",
        "description": "C2H2 zinc finger domain HMM",
    },
}

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
    
    def _find_command_path(self, cmd: str) -> Optional[str]:
        """
        查找命令的完整路径，支持用户目录下的工具
        
        GUI 应用（如 PyMOL.app）可能不继承终端的 shell PATH，
        因此需要额外搜索用户常用路径。
        
        Args:
            cmd: 命令名称
            
        Returns:
            str: 命令路径，未找到返回 None
        """
        import shutil
        
        # 1. 系统 PATH
        path = shutil.which(cmd)
        if path:
            return path
        
        # 2. Conda 环境
        if 'CONDA_PREFIX' in os.environ:
            conda_path = os.path.join(os.environ['CONDA_PREFIX'], 'bin', cmd)
            if os.path.exists(conda_path):
                return conda_path
        
        # 3. 用户常用路径
        home = os.path.expanduser('~')
        user_paths = [
            os.path.join(home, 'bin', cmd),
            os.path.join(home, '.local', 'bin', cmd),
            os.path.join(home, 'local', 'bin', cmd),
            f'/usr/local/bin/{cmd}',
        ]
        
        # 4. 平台特定路径
        if sys.platform == 'darwin':  # macOS
            user_paths.extend([
                f'/opt/homebrew/bin/{cmd}',  # Apple Silicon Homebrew
            ])
        
        for path in user_paths:
            if os.path.exists(path):
                return path
        
        return None
    
    def check_command(self, cmd: str, name: str) -> bool:
        """
        检查命令行工具是否可用
        
        Args:
            cmd: 命令名称
            name: 显示名称
            
        Returns:
            bool: True 如果命令可用
        """
        # 使用增强的路径查找
        cmd_path = self._find_command_path(cmd)
        if not cmd_path:
            self.log(f"  ✗ {name} (未找到)")
            return False
        
        try:
            result = subprocess.run(
                [cmd_path, "--version"],
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
        
        # 测试 PyQt5 - 使用子进程避免崩溃 (非常重要，在 macOS 上尤其如此)
        try:
            import subprocess
            import os
            
            # [macOS Fix] 设置环境变量以避免某些 Qt 绘图引起的崩溃
            env = os.environ.copy()
            if sys.platform == "darwin":
                env["QT_MAC_WANTS_LAYER"] = "1"

            test_code = "from PyQt5.QtWidgets import QApplication; from PyQt5.QtCore import Qt; print('OK')"
            result = subprocess.run(
                [sys.executable, "-c", test_code],
                capture_output=True,
                text=True,
                timeout=5,
                env=env
            )
            if result.returncode == 0 and "OK" in result.stdout:
                self.log("  ✓ PyQt5 基础功能")
            else:
                errors.append(f"PyQt5: Subprocess test failed")
                self.log(f"  ✗ PyQt5: Subprocess test failed")
                if result.stderr:
                    self.log(f"    Error: {result.stderr.strip()}")
        except Exception as e:
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
    
    # ========== 可选依赖检测与安装 ==========
    
    def check_optional_packages(self) -> Dict[str, bool]:
        """
        检查可选 Python 包
        
        Returns:
            Dict[str, bool]: {包名: 是否已安装}
        """
        self.log("\n📦 检查可选 Python 包:")
        status = {}
        
        for import_name, display_name, _, description in OPTIONAL_PACKAGES:
            try:
                __import__(import_name)
                self.log(f"  ✓ {display_name} - {description}")
                status[display_name] = True
            except ImportError:
                self.log(f"  ✗ {display_name} (未安装) - {description}")
                status[display_name] = False
        
        return status
    
    def check_external_tools(self) -> Dict[str, Dict]:
        """
        检查外部工具（MSMS, APBS, PDB2PQR 等）
        
        Returns:
            Dict[str, Dict]: {工具名: {available: bool, path: str, python_api: bool, install_info: str}}
        """
        self.log("\n🔧 检查外部工具:")
        
        if not self.os_type:
            self.detect_os()
        
        results = {}
        
        for tool_name, tool_info in EXTERNAL_TOOLS.items():
            display_name = tool_info["display_name"]
            description = tool_info["description"]
            
            # 首先检查 Python API 是否可用
            python_api_available = False
            python_module = tool_info.get("python_module")
            if python_module:
                try:
                    __import__(python_module)
                    python_api_available = True
                except ImportError:
                    pass
            
            # 查找命令行工具路径
            tool_path = self._find_external_tool(tool_name, tool_info)
            
            if python_api_available or tool_path:
                if python_api_available and tool_path:
                    self.log(f"  ✓ {display_name} - Python API + CLI ({tool_path})")
                elif python_api_available:
                    self.log(f"  ✓ {display_name} - Python API 可用")
                else:
                    self.log(f"  ✓ {display_name} - {tool_path}")
                
                results[tool_name] = {
                    "available": True,
                    "path": tool_path,
                    "python_api": python_api_available,
                    "description": description,
                    "install_info": None
                }
            else:
                install_info = tool_info["install_info"].get(self.os_type, "请参考官方文档安装")
                self.log(f"  ✗ {display_name} (未找到) - {description}")
                self.log(f"      安装方法: {install_info}")
                results[tool_name] = {
                    "available": False,
                    "path": None,
                    "python_api": False,
                    "description": description,
                    "install_info": install_info
                }
        
        return results
    
    def _find_external_tool(self, tool_name: str, tool_info: Dict) -> Optional[str]:
        """
        查找外部工具路径
        
        Args:
            tool_name: 工具名称
            tool_info: 工具配置信息
            
        Returns:
            str: 工具路径，未找到返回 None
        """
        import shutil
        
        # 1. 系统 PATH
        path = shutil.which(tool_name)
        if path:
            return path
        
        # 2. 平台特定搜索路径
        if self.os_type and self.os_type in tool_info.get("search_paths", {}):
            for search_path in tool_info["search_paths"][self.os_type]:
                expanded_path = os.path.expanduser(search_path)
                if os.path.isfile(expanded_path) and os.access(expanded_path, os.X_OK):
                    return expanded_path
        
        # 3. Conda 环境
        if 'CONDA_PREFIX' in os.environ:
            conda_path = os.path.join(os.environ['CONDA_PREFIX'], 'bin', tool_name)
            if os.path.exists(conda_path):
                return conda_path
        
        return None
    
    def setup_surface_analysis(self) -> bool:
        """
        一键设置表面分析环境
        
        包括：
        1. 安装 Open3D
        2. 安装 scikit-image
        3. 检查 MSMS/APBS（提供安装指南）
        
        Returns:
            bool: 是否全部成功
        """
        self.log("\n" + "=" * 60)
        self.log("🔧 设置表面分析环境 (MaSIF-style)")
        self.log("=" * 60)
        
        success = True
        
        # 1. 安装 Open3D
        try:
            import open3d
            self.log("  ✓ Open3D 已安装")
        except ImportError:
            self.log("  安装 Open3D...")
            if self._install_pip_package("open3d"):
                self.log("  ✅ Open3D 安装成功")
            else:
                self.log("  ⚠️ Open3D 安装失败，将使用内置回退方案")
                success = False
        
        # 2. 安装 scikit-image
        try:
            import skimage
            self.log("  ✓ scikit-image 已安装")
        except ImportError:
            self.log("  安装 scikit-image...")
            if self._install_pip_package("scikit-image"):
                self.log("  ✅ scikit-image 安装成功")
            else:
                self.log("  ⚠️ scikit-image 安装失败，将使用简化算法")
        
        # 3. 检查外部工具
        self.log("\n  检查外部工具...")
        tools = self.check_external_tools()
        
        msms_available = tools.get("msms", {}).get("available", False)
        apbs_available = tools.get("apbs", {}).get("available", False)
        
        if not msms_available:
            self.log("\n  💡 MSMS 未安装，表面生成将使用 Open3D 或内置方法")
            self.log(f"     推荐安装: {tools.get('msms', {}).get('install_info', '参考官方文档')}")
        
        if not apbs_available:
            self.log("\n  💡 APBS 未安装，静电势将使用 Coulomb 近似")
            self.log(f"     推荐安装: {tools.get('apbs', {}).get('install_info', '参考官方文档')}")
        
        if success:
            self.log("\n✅ 表面分析环境设置完成")
        else:
            self.log("\n⚠️ 表面分析可用（基础模式），完整功能需要手动配置")
        
        return success
    
    def install_optional_package(self, package_name: str) -> bool:
        """
        安装可选包
        
        Args:
            package_name: 包的显示名称
            
        Returns:
            bool: 是否安装成功
        """
        pip_name = None
        for _, display_name, pip, _ in OPTIONAL_PACKAGES:
            if display_name == package_name:
                pip_name = pip
                break
        
        if not pip_name:
            self.log(f"  ✗ 未知的可选包: {package_name}")
            return False
        
        return self._install_pip_package(pip_name)
    
    def install_all_optional_packages(self) -> Dict[str, bool]:
        """
        安装所有可选包
        
        Returns:
            Dict[str, bool]: {包名: 是否安装成功}
        """
        self.log("\n📦 安装可选 Python 包...")
        results = {}
        
        for import_name, display_name, pip_name, description in OPTIONAL_PACKAGES:
            try:
                __import__(import_name)
                self.log(f"  ✓ {display_name} 已安装")
                results[display_name] = True
            except ImportError:
                self.log(f"  安装 {display_name} ({description})...")
                if self._install_pip_package(pip_name):
                    self.log(f"  ✅ {display_name} 安装成功")
                    results[display_name] = True
                else:
                    self.log(f"  ✗ {display_name} 安装失败")
                    results[display_name] = False
        
        return results
    
    # ========== HMM 文件管理 ==========
    
    def get_hmm_data_dir(self) -> str:
        """
        获取 HMM 数据目录路径
        
        Returns:
            str: 数据目录路径
        """
        # 优先使用插件目录下的 data 文件夹
        plugin_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(plugin_dir, "data")
        
        if not os.path.exists(data_dir):
            try:
                os.makedirs(data_dir, exist_ok=True)
            except Exception:
                # 回退到用户目录
                data_dir = os.path.join(os.path.expanduser("~"), ".gluetk", "data")
                os.makedirs(data_dir, exist_ok=True)
        
        return data_dir
    
    def check_hmm_files(self) -> Dict[str, bool]:
        """
        检查 HMM 文件是否存在
        
        Returns:
            Dict[str, bool]: {HMM名称: 是否存在}
        """
        self.log("\n🔬 检查 HMM Profile 文件:")
        data_dir = self.get_hmm_data_dir()
        status = {}
        
        for name, info in HMM_FILES.items():
            local_path = os.path.join(data_dir, info["local_name"])
            if os.path.exists(local_path):
                self.log(f"  ✓ {name} ({info['pfam_id']}) - {info['description']}")
                status[name] = True
            else:
                self.log(f"  ✗ {name} ({info['pfam_id']}) - 未下载")
                status[name] = False
        
        return status
    
    def download_hmm_file(self, hmm_name: str) -> bool:
        """
        下载单个 HMM 文件
        
        Args:
            hmm_name: HMM 名称（如 "zf-C2H2"）
            
        Returns:
            bool: 是否下载成功
        """
        if hmm_name not in HMM_FILES:
            self.log(f"  ✗ 未知的 HMM: {hmm_name}")
            return False
        
        info = HMM_FILES[hmm_name]
        data_dir = self.get_hmm_data_dir()
        local_path = os.path.join(data_dir, info["local_name"])
        
        if os.path.exists(local_path):
            self.log(f"  ✓ {hmm_name} 已存在")
            return True
        
        self.log(f"  下载 {hmm_name} ({info['pfam_id']})...")
        
        try:
            import urllib.request
            
            # 下载 HMM 文件
            urllib.request.urlretrieve(info["url"], local_path)
            
            # 验证文件
            if os.path.exists(local_path) and os.path.getsize(local_path) > 100:
                self.log(f"  ✅ {hmm_name} 下载成功: {local_path}")
                return True
            else:
                self.log(f"  ✗ {hmm_name} 下载失败: 文件无效")
                if os.path.exists(local_path):
                    os.remove(local_path)
                return False
                
        except Exception as e:
            self.log(f"  ✗ {hmm_name} 下载失败: {e}")
            return False
    
    def download_all_hmm_files(self) -> Dict[str, bool]:
        """
        下载所有 HMM 文件
        
        Returns:
            Dict[str, bool]: {HMM名称: 是否下载成功}
        """
        self.log("\n🔬 下载 HMM Profile 文件...")
        results = {}
        
        for name in HMM_FILES:
            results[name] = self.download_hmm_file(name)
        
        return results
    
    def setup_c2h2_detection(self) -> bool:
        """
        一键设置 C2H2 锌指检测环境
        
        包括：
        1. 安装 pyhmmer
        2. 下载 zf-C2H2 HMM 文件
        
        Returns:
            bool: 是否全部成功
        """
        self.log("\n" + "=" * 60)
        self.log("🔧 设置 C2H2 锌指检测环境")
        self.log("=" * 60)
        
        success = True
        
        # 1. 安装 pyhmmer
        try:
            import pyhmmer
            self.log("  ✓ pyhmmer 已安装")
        except ImportError:
            self.log("  安装 pyhmmer...")
            if self._install_pip_package("pyhmmer"):
                self.log("  ✅ pyhmmer 安装成功")
            else:
                self.log("  ⚠️ pyhmmer 安装失败，将使用 regex 回退方案")
                success = False
        
        # 2. 下载 HMM 文件
        if not self.download_hmm_file("zf-C2H2"):
            self.log("  ⚠️ HMM 文件下载失败，将使用 regex 回退方案")
            success = False
        
        if success:
            self.log("\n✅ C2H2 检测环境设置完成（HMM 模式）")
        else:
            self.log("\n⚠️ C2H2 检测可用（regex 模式），HMM 模式需要手动配置")
        
        return success
    
    def setup_ec_analysis(self) -> bool:
        """
        一键设置电性互补性（EC）分析环境
        
        包括：
        1. 安装 PDB2PQR (Python API)
        2. 检查/安装 APBS
        3. 确保 NumPy, SciPy, RDKit 可用
        
        Returns:
            bool: 是否全部成功
        """
        self.log("\n" + "=" * 60)
        self.log("🔧 设置电性互补性（EC）分析环境")
        self.log("=" * 60)
        
        success = True
        
        # 1. 检查核心依赖
        self.log("\n  检查核心依赖...")
        core_deps = [("numpy", "NumPy"), ("scipy", "SciPy"), ("rdkit", "RDKit")]
        for import_name, display_name in core_deps:
            try:
                __import__(import_name)
                self.log(f"  ✓ {display_name} 已安装")
            except ImportError:
                self.log(f"  ✗ {display_name} 未安装 - EC分析需要此依赖")
                success = False
        
        # 2. 安装 PDB2PQR
        self.log("\n  检查 PDB2PQR...")
        try:
            import pdb2pqr
            self.log("  ✓ PDB2PQR Python API 已安装")
        except ImportError:
            self.log("  安装 PDB2PQR...")
            if self._install_pip_package("pdb2pqr"):
                self.log("  ✅ PDB2PQR 安装成功")
            else:
                self.log("  ⚠️ PDB2PQR 安装失败，将使用简化的PQR转换")
        
        # 3. 检查 APBS
        self.log("\n  检查 APBS...")
        apbs_available = False
        
        # 检查 Python API
        try:
            import apbs
            self.log("  ✓ APBS Python API 已安装")
            apbs_available = True
        except ImportError:
            self.log("  ✗ APBS Python API 未安装")
        
        # 检查命令行工具
        import shutil
        apbs_path = shutil.which('apbs')
        if apbs_path:
            self.log(f"  ✓ APBS CLI 可用: {apbs_path}")
            apbs_available = True
        else:
            self.log("  ✗ APBS CLI 未找到")
        
        if not apbs_available:
            self.log("\n  ⚠️ APBS 未安装，尝试自动安装...")
            if self._install_pip_package("apbs"):
                self.log("  ✅ APBS (pip) 安装成功")
                apbs_available = True
            else:
                self.log("  ⚠️ APBS pip 安装失败，EC分析将无法运行")
                self.log("     推荐安装方法:")
                if self.os_type == "macOS":
                    self.log("     - pip install apbs")
                    self.log("     - 或 brew install brewsci/bio/apbs")
                elif self.os_type == "Linux":
                    self.log("     - pip install apbs")
                    self.log("     - 或 apt install apbs")
                else:
                    self.log("     - pip install apbs")
                    self.log("     - 或从 https://www.poissonboltzmann.org/ 下载")
            success = False
        
        # 4. 总结
        if success:
            self.log("\n✅ EC分析环境设置完成")
            self.log("   可以使用以下功能:")
            self.log("   - calculate_ligand_ec: 蛋白-配体EC分析")
            self.log("   - analyze_ternary_ec: 三元复合物（分子胶）EC分析")
            self.log("   - compare_ligand_ec: 多配体EC比较（SAR分析）")
            self.log("   - calculate_ec_hotspots: EC热点识别")
        else:
            self.log("\n⚠️ EC分析环境部分可用，完整功能需要手动配置")
        
        return success
    
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
    def _install_pip_package(self, package_name: str) -> bool:
        """
        通过 pip 安装包 (fallback)
        """
        try:
            # 使用当前 python 环境的 pip
            cmd = [sys.executable, "-m", "pip", "install", package_name, "--quiet", "--disable-pip-version-check"]
            self.log(f"  执行 pip: {' '.join(cmd)}")
            subprocess.check_call(cmd)
            return True
        except Exception as e:
            self.log(f"  ✗ Pip 安装失败: {e}")
            return False

    def auto_install_dependencies(self) -> Dict[str, bool]:
        """
        自动安装缺失的依赖
        
        Returns:
            Dict[str, bool]: {依赖名: 是否安装成功}
        """
        use_conda = self.check_conda()
        if not use_conda:
             self.log("⚠️ Conda 未安装，将尝试使用 pip 安装")
        
        status = self.check_all_dependencies()
        missing = [name for name, avail in status.items() if not avail]
        
        if not missing:
            self.log("✅ 所有依赖已安装，无需操作")
            return status
        
        self.log(f"\\n📦 开始自动安装缺失的依赖...")
        
        install_results = {}
        
        # 批量安装 Python 包
        py_packages = missing
        if py_packages:
            self.log(f"\\n  安装 Python 包: {', '.join(py_packages)}")
            
            # 尝试 Conda 安装（仅限真正有 Conda 包的 Python 库，不包含 Vina）
            if use_conda:
                packages_to_install = []
                for import_name, display_name, pip_name in REQUIRED_PACKAGES:
                    if display_name in py_packages:
                        # Meeko 在 conda-forge 上有包，但 Vina 仅通过 pip 提供
                        if display_name == "Meeko":
                            packages_to_install.append("meeko")
                        else:
                            packages_to_install.append(pip_name)
                
                if packages_to_install:
                    try:
                        cmd = ["conda", "install", "-c", "conda-forge", "-y"] + packages_to_install
                        self.log(f"  执行: {' '.join(cmd)}")
                        
                        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
                        
                        if result.returncode == 0:
                            for pkg in py_packages:
                                install_results[pkg] = True
                                self.log(f"  ✅ {pkg} 安装成功 (Conda)")
                        else:
                            self.log(f"  ⚠️ Conda 安装失败，尝试 Pip...")
                    except Exception as e:
                        self.log(f"  ⚠️ Conda 出错 ({e})，尝试 Pip...")
            
            # 无论 Conda 成功与否，对所有缺失包再尝试一次 pip（避免架构差异问题）
            for pkg in py_packages:
                pip_name = next((p[2] for p in REQUIRED_PACKAGES if p[1] == pkg), None)
                if not pip_name:
                    continue
                if self._install_pip_package(pip_name):
                    install_results[pkg] = True
                    self.log(f"  ✅ {pkg} 安装成功 (Pip)")
                else:
                    install_results[pkg] = False
                    self.log(f"  ✗ {pkg} 安装失败")
            else:
                # No Conda, use Pip directly
                for pkg in py_packages:
                    pip_name = next((p[2] for p in REQUIRED_PACKAGES if p[1] == pkg), None)
                    if pip_name:
                        if self._install_pip_package(pip_name):
                            install_results[pkg] = True
                            self.log(f"  ✅ {pkg} 安装成功 (Pip)")
                        else:
                            install_results[pkg] = False
                            self.log(f"  ✗ {pkg} 安装失败")

        # 总结
        success_count = sum(1 for v in install_results.values() if v)
        total_count = len(install_results)
        
        if success_count == total_count:
            self.log(f"\\n✅ 所有依赖安装完成 ({success_count}/{total_count})")
        else:
            self.log(f"\\n⚠️  部分依赖安装失败 ({success_count}/{total_count})")
        
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
            instructions["dependencies"] = f"""
📦 安装所有依赖（Vina 可选，推荐用 pip 安装）:

# 核心栈（不含 Vina）
conda install -c conda-forge rdkit scipy matplotlib pillow \\
    numpy pandas seaborn pyqt openbabel -y

# 可选：在激活环境中安装 Vina + Meeko
pip install vina meeko

# 可选：C2H2 锌指检测 HMM 支持
pip install pyhmmer
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
        
        # 4. 检查所有依赖
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


def setup_c2h2_detection(log_callback=None) -> bool:
    """
    一键设置 C2H2 锌指检测环境（便捷函数）
    
    包括：
    1. 安装 pyhmmer
    2. 下载 zf-C2H2 HMM 文件
    
    Args:
        log_callback: 日志回调函数
        
    Returns:
        bool: 是否全部成功
    """
    checker = EnvironmentChecker(log_callback=log_callback)
    return checker.setup_c2h2_detection()


def setup_surface_analysis(log_callback=None) -> bool:
    """
    一键设置表面分析环境（便捷函数）
    
    包括：
    1. 安装 Open3D
    2. 安装 scikit-image
    3. 检查 MSMS/APBS
    
    Args:
        log_callback: 日志回调函数
        
    Returns:
        bool: 是否全部成功
    """
    checker = EnvironmentChecker(log_callback=log_callback)
    return checker.setup_surface_analysis()


def setup_ec_analysis(log_callback=None) -> bool:
    """
    一键设置电性互补性（EC）分析环境（便捷函数）
    
    包括：
    1. 安装 PDB2PQR (Python API)
    2. 检查/安装 APBS
    3. 确保 NumPy, SciPy, RDKit 可用
    
    Args:
        log_callback: 日志回调函数
        
    Returns:
        bool: 是否全部成功
    """
    checker = EnvironmentChecker(log_callback=log_callback)
    return checker.setup_ec_analysis()


def check_ec_analysis_deps(log_callback=None) -> Dict[str, bool]:
    """
    检查EC分析依赖状态（便捷函数）
    
    Args:
        log_callback: 日志回调函数
        
    Returns:
        Dict[str, bool]: {依赖名: 是否可用}
    """
    checker = EnvironmentChecker(log_callback=log_callback)
    
    status = {}
    
    # 核心 Python 包
    core_packages = [
        ("numpy", "NumPy"),
        ("scipy", "SciPy"),
        ("rdkit", "RDKit"),
    ]
    
    for import_name, display_name in core_packages:
        try:
            __import__(import_name)
            status[display_name] = True
        except ImportError:
            status[display_name] = False
    
    # EC 专用包
    ec_packages = [
        ("pdb2pqr", "PDB2PQR"),
        ("apbs", "APBS Python"),
    ]
    
    for import_name, display_name in ec_packages:
        try:
            __import__(import_name)
            status[display_name] = True
        except ImportError:
            status[display_name] = False
    
    # 外部工具（命令行）
    tools = checker.check_external_tools()
    for tool_name in ["apbs", "pdb2pqr"]:
        tool_info = tools.get(tool_name, {})
        cli_name = f"{tool_info.get('display_name', tool_name)} CLI"
        status[cli_name] = tool_info.get("path") is not None
    
    return status


def check_surface_analysis_deps(log_callback=None) -> Dict[str, bool]:
    """
    检查表面分析依赖状态（便捷函数）
    
    Args:
        log_callback: 日志回调函数
        
    Returns:
        Dict[str, bool]: {依赖名: 是否可用}
    """
    checker = EnvironmentChecker(log_callback=log_callback)
    
    status = {}
    
    # Python 包
    for import_name, display_name, _, _ in OPTIONAL_PACKAGES:
        if display_name in ["Open3D", "scikit-image", "trimesh"]:
            try:
                __import__(import_name)
                status[display_name] = True
            except ImportError:
                status[display_name] = False
    
    # 外部工具
    tools = checker.check_external_tools()
    for tool_name, tool_info in tools.items():
        status[tool_info.get("display_name", tool_name)] = tool_info.get("available", False)
    
    return status


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

def is_in_conda_env() -> bool:
    """
    检查当前是否在 conda 环境中运行
    
    Returns:
        bool: True 如果在 conda 环境中
    """
    return bool(os.environ.get('CONDA_PREFIX', ''))


def is_first_run() -> bool:
    """
    检查是否是第一次运行
    
    Returns:
        bool: True 如果是第一次运行
    """
    return not os.path.exists(_FIRST_RUN_MARKER)


def mark_initialized():
    """标记已初始化"""
    try:
        with open(_FIRST_RUN_MARKER, 'w') as f:
            f.write('initialized')
    except:
        pass


def _quick_check_deps(skip_pyqt=True) -> Tuple[bool, List[str]]:
    """
    快速检查依赖（不输出日志）
    
    Args:
        skip_pyqt: 是否跳过 PyQt5 检查（避免导入时崩溃）
    
    Returns:
        (all_ok, missing_list)
    """
    import shutil
    missing = []
    
    # 检查 Python 包（跳过 PyQt5 以避免崩溃）
    for import_name, display_name, _ in REQUIRED_PACKAGES:
        # 跳过 PyQt5 检查，因为导入它可能导致崩溃
        if skip_pyqt and import_name == "PyQt5":
            continue
        try:
            __import__(import_name)
        except ImportError:
            missing.append(display_name)
    
    # 检查命令行工具
    for cmd, name in REQUIRED_COMMANDS:
        if not shutil.which(cmd):
            missing.append(name)
    
    return len(missing) == 0, missing


def _print_setup_guide(missing: List[str]):
    """打印完整的环境设置指南"""
    print("\n" + "=" * 60)
    print("🚨 GlueTK - Environment Setup Required")
    print("=" * 60)
    
    if missing:
        print(f"\n❌ Missing: {', '.join(missing)}")
    
    print("\n📖 One-click setup (run in terminal):")
    print("─" * 60)
    print("   bash /path/to/GlueTK/install.sh")
    print("─" * 60)
    
    print("\n📖 Or manual setup:")
    print(f"""
  conda create -n {ENV_NAME} python={PYTHON_VERSION} -y
  conda activate {ENV_NAME}
  conda install -c conda-forge rdkit scipy matplotlib pillow \\
      numpy pandas seaborn pyqt autodock-vina openbabel -y
""")
    print("💡 After setup, always start PyMOL from activated environment:")
    print(f"   conda activate {ENV_NAME} && pymol")
    print("=" * 60 + "\n")


def _print_missing_deps_hint(missing: List[str]):
    """打印缺失依赖的简短提示"""
    print("\n" + "=" * 60)
    print("⚠️  GlueTK - Missing Dependencies")
    print("=" * 60)
    print(f"\n❌ Missing: {', '.join(missing)}")
    
    # 构建包名映射
    pkg_map = {p[1]: p[2] for p in REQUIRED_PACKAGES}
    pkg_map["AutoDock Vina"] = "autodock-vina"
    pkg_map["Open Babel"] = "openbabel"
    
    conda_names = [pkg_map.get(m, m.lower()) for m in missing]
    
    print(f"\n💡 Install missing packages:")
    print(f"   conda activate {ENV_NAME}")
    print(f"   conda install -c conda-forge {' '.join(conda_names)} -y")
    print("=" * 60 + "\n")


def ensure_dependencies(silent=False) -> bool:
    """
    检查并确保所有依赖可用。
    
    第一次运行时显示完整的环境设置指南。
    之后运行时，如果依赖缺失，显示简短提示。
    
    Args:
        silent: 如果为 True，则不输出任何信息
    
    Returns:
        bool: True 如果所有依赖都可用
    """
    # 快速检查依赖
    all_ok, missing = _quick_check_deps()
    
    if all_ok:
        # 所有依赖都存在，标记已初始化
        mark_initialized()
        return True
    
    # 依赖缺失
    if silent:
        return False
    
    # 检查是否是第一次运行
    if is_first_run():
        # 第一次运行，显示完整设置指南
        _print_setup_guide(missing)
    else:
        # 非第一次运行，显示简短提示
        _print_missing_deps_hint(missing)
    
    return False
