# -*- coding: utf-8 -*-
"""
GLINT Environment Detection and Configuration Module
Supports Windows/macOS/Linux environment detection, Conda installation, dependency installation, GUI functionality testing
"""

import os
import sys
import platform
import subprocess
import tempfile
from typing import Dict, List, Tuple, Optional

# Environment configuration
ENV_NAME = "glint"
PYTHON_VERSION = "3.9"

# First-run marker file
_FIRST_RUN_MARKER = os.path.join(os.path.expanduser("~"), ".glint_initialized")

# Required Python packages (import_name, display_name, pip_name)
# All installed during setup, checked at startup
REQUIRED_PACKAGES = [
    ("rdkit", "RDKit", "rdkit"),
    ("scipy", "SciPy", "scipy"),
    ("matplotlib", "Matplotlib", "matplotlib"),
    ("PIL", "Pillow", "pillow"),
    ("numpy", "NumPy", "numpy"),
    ("pandas", "Pandas", "pandas"),
    ("seaborn", "Seaborn", "seaborn"),
    ("PyQt5", "PyQt5", "pyqt5"),
]

# Optional advanced feature packages - not checked at startup, prompted on demand
# (import_name, display_name, pip_name, description)
OPTIONAL_PACKAGES = [
    ("open3d", "Open3D", "open3d", "Surface analysis"),
    ("skimage", "scikit-image", "scikit-image", "Marching Cubes algorithm"),
    ("haddock", "HADDOCK3", "haddock3", "Protein-protein docking"),
    ("trimesh", "trimesh", "trimesh", "Mesh processing"),
]

# EC analysis required Python packages (import_name, display_name, pip_name, description)
EC_REQUIRED_PACKAGES = [
    ("pdb2pqr", "PDB2PQR", "pdb2pqr", "Protein structure preparation (required for EC analysis)"),
    ("apbs", "APBS Python", "apbs", "Poisson-Boltzmann solver Python API"),
]

# External tool configuration (cmd_name, display_name, description, install_info)
EXTERNAL_TOOLS = {
    "apbs": {
        "display_name": "APBS",
        "description": "Adaptive Poisson-Boltzmann solver (accurate electrostatic potential, required for EC and surface analysis)",
        "install_info": {
            "macOS": "pip install apbs or brew install brewsci/bio/apbs",
            "Linux": "pip install apbs or apt install apbs",
            "Windows": "pip install apbs or download from https://www.poissonboltzmann.org/",
        },
        "search_paths": {
            "macOS": ["/usr/local/bin/apbs", "/opt/homebrew/bin/apbs", "~/bin/apbs"],
            "Linux": ["/usr/bin/apbs", "/usr/local/bin/apbs", "~/bin/apbs"],
            "Windows": [r"C:\Program Files\APBS\apbs.exe", r"C:\APBS\apbs.exe"],
        },
        "python_module": "apbs",  # Can also be used via Python API
    },
    "pdb2pqr": {
        "display_name": "PDB2PQR",
        "description": "PDB to PQR format conversion (required for EC analysis)",
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
        "python_module": "pdb2pqr",  # Recommended to use via Python API
    },
}

# Core required command-line tools (currently empty, core features don't depend on external commands)
REQUIRED_COMMANDS = []

# Optional command-line tools - only checked when specific features are needed
# (cmd_name, display_name, description)
OPTIONAL_COMMANDS = [
    ("vina", "AutoDock Vina", "Molecular docking"),
    ("obabel", "Open Babel", "Molecular format conversion"),
]


class EnvironmentChecker:
    """Environment detection and configuration class"""

    def __init__(self, log_callback=None, auto_install=False):
        """
        Args:
            log_callback: Log callback function for outputting messages to GUI
            auto_install: Whether to automatically install missing dependencies
        """
        self.log_callback = log_callback or print
        self.auto_install = auto_install
        self.os_type = None
        self.os_arch = None
        self.conda_path = None

    def log(self, msg: str):
        """Output a log message."""
        if self.log_callback:
            self.log_callback(msg)

    # ========== System Detection ==========

    def detect_os(self) -> Tuple[str, str]:
        """
        Detect operating system and architecture.

        Returns:
            (os_type, arch): e.g. ("macOS", "arm64")
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

        self.log(f"✓ System: {self.os_type} ({self.os_arch})")
        return self.os_type, self.os_arch

    # ========== Conda Detection ==========

    def check_conda(self) -> bool:
        """
        Check if conda is installed.

        Returns:
            bool: True if conda is available
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

        self.log("✗ Conda: not detected")
        return False

    def check_conda_env(self, env_name: str = ENV_NAME) -> bool:
        """
        Check if the specified conda environment exists.

        Args:
            env_name: Environment name

        Returns:
            bool: True if the environment exists
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
                        self.log(f"✓ Environment '{env_name}' exists")
                        return True
        except subprocess.TimeoutExpired:
            pass

        self.log(f"✗ Environment '{env_name}' does not exist")
        return False

    # ========== Dependency Detection ==========

    def check_package(self, import_name: str, display_name: str) -> bool:
        """
        Check if a Python package is installed.

        Args:
            import_name: Import name
            display_name: Display name

        Returns:
            bool: True if the package is installed
        """
        try:
            __import__(import_name)
            self.log(f"  ✓ {display_name}")
            return True
        except ImportError:
            self.log(f"  ✗ {display_name} (not installed)")
            return False

    def _find_command_path(self, cmd: str) -> Optional[str]:
        """
        Find the full path of a command, including tools under user directories.

        GUI launches (e.g., PyMOL.app) may not inherit the terminal shell PATH,
        so we additionally search common user locations.

        Args:
            cmd: Command name

        Returns:
            str: Command path, or None if not found
        """
        import shutil

        # 1. System PATH
        path = shutil.which(cmd)
        if path:
            return path

        # 2. Conda environment
        if 'CONDA_PREFIX' in os.environ:
            conda_path = os.path.join(os.environ['CONDA_PREFIX'], 'bin', cmd)
            if os.path.exists(conda_path):
                return conda_path

        # 3. Common user paths
        home = os.path.expanduser('~')
        user_paths = [
            os.path.join(home, 'bin', cmd),
            os.path.join(home, '.local', 'bin', cmd),
            os.path.join(home, 'local', 'bin', cmd),
            f'/usr/local/bin/{cmd}',
        ]

        # 4. Platform-specific paths
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
        Check if a command-line tool is available.

        Args:
            cmd: Command name
            name: Display name

        Returns:
            bool: True if the command is available
        """
        # Use enhanced path lookup
        cmd_path = self._find_command_path(cmd)
        if not cmd_path:
            self.log(f"  ✗ {name} (not found)")
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

        self.log(f"  ✗ {name} (not found)")
        return False

    def check_all_dependencies(self) -> Dict[str, bool]:
        """
        Check all dependency items.

        Returns:
            Dict[str, bool]: {Dependency name: whether installed}
        """
        self.log("\n📦 Checking Python packages:")

        status = {}
        all_ok = True

        # Required packages
        for import_name, display_name, _ in REQUIRED_PACKAGES:
            available = self.check_package(import_name, display_name)
            status[display_name] = available
            if not available:
                all_ok = False

        # Command-line tools (if any are required)
        if REQUIRED_COMMANDS:
            self.log("\n🛠️  Checking command-line tools:")
            for cmd, name in REQUIRED_COMMANDS:
                available = self.check_command(cmd, name)
                status[name] = available
                if not available:
                    all_ok = False

        if all_ok:
            self.log("\n✅ All dependencies are installed")
        else:
            missing = [k for k, v in status.items() if not v]
            self.log(f"\n⚠️  Missing dependencies: {', '.join(missing)}")

        return status

    # ========== GUI Functionality Tests ==========

    def test_gui_functionality(self) -> bool:
        """
        Test whether GUI-related functionality works properly.

        Returns:
            bool: True if all tests pass
        """
        self.log("\n🧪 Testing GUI functionality:")

        errors = []

        # Test PyQt5 - use a subprocess to avoid crashing (very important, especially on macOS)
        try:
            import subprocess
            import os

            # [macOS Fix] Set environment variables to avoid some Qt rendering crashes
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
                self.log("  ✓ PyQt5 basic functionality")
            else:
                errors.append("PyQt5: Subprocess test failed")
                self.log("  ✗ PyQt5: Subprocess test failed")
                if result.stderr:
                    self.log(f"    Error: {result.stderr.strip()}")
        except Exception as e:
            errors.append(f"PyQt5: {e}")
            self.log(f"  ✗ PyQt5: {e}")

        # Test Matplotlib
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.plot([1, 2, 3], [1, 2, 3])
            plt.close(fig)
            self.log("  ✓ Matplotlib plotting")
        except Exception as e:
            errors.append(f"Matplotlib: {e}")
            self.log(f"  ✗ Matplotlib: {e}")

        # Test RDKit
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem, Descriptors
            mol = Chem.MolFromSmiles('CCO')
            if mol:
                AllChem.Compute2DCoords(mol)
                mw = Descriptors.MolWt(mol)
                self.log(f"  ✓ RDKit chemistry calculation (ethanol MW={mw:.2f})")
            else:
                errors.append("RDKit: failed to create molecule")
                self.log("  ✗ RDKit: failed to create molecule")
        except Exception as e:
            errors.append(f"RDKit: {e}")
            self.log(f"  ✗ RDKit: {e}")

        # Test NumPy/SciPy
        try:
            import numpy as np
            from scipy import spatial
            points = np.random.rand(10, 3)
            dist = spatial.distance.cdist(points, points)
            self.log("  ✓ NumPy/SciPy distance computation")
        except Exception as e:
            errors.append(f"NumPy/SciPy: {e}")
            self.log(f"  ✗ NumPy/SciPy: {e}")

        # Test Pillow
        try:
            from PIL import Image, ImageDraw
            img = Image.new('RGB', (100, 100), color='white')
            draw = ImageDraw.Draw(img)
            draw.rectangle([10, 10, 90, 90], outline='black')
            self.log("  ✓ Pillow image processing")
        except Exception as e:
            errors.append(f"Pillow: {e}")
            self.log(f"  ✗ Pillow: {e}")

        if errors:
            self.log("\n⚠️  Some functionality tests failed")
            return False
        else:
            self.log("\n✅ All GUI functionality tests passed")
            return True

    # ========== Optional Dependency Checks and Installation ==========

    def check_optional_packages(self) -> Dict[str, bool]:
        """
        Check optional Python packages.

        Returns:
            Dict[str, bool]: {Package name: whether installed}
        """
        self.log("\n📦 Checking optional Python packages:")
        status = {}

        for import_name, display_name, _, description in OPTIONAL_PACKAGES:
            try:
                __import__(import_name)
                self.log(f"  ✓ {display_name} - {description}")
                status[display_name] = True
            except ImportError:
                self.log(f"  ✗ {display_name} (not installed) - {description}")
                status[display_name] = False

        return status

    def check_external_tools(self) -> Dict[str, Dict]:
        """
        Check external tools (MSMS, APBS, PDB2PQR, etc.).

        Returns:
            Dict[str, Dict]: {tool: {available: bool, path: str, python_api: bool, install_info: str}}
        """
        self.log("\n🔧 Checking external tools:")

        if not self.os_type:
            self.detect_os()

        results = {}

        for tool_name, tool_info in EXTERNAL_TOOLS.items():
            display_name = tool_info["display_name"]
            description = tool_info["description"]

            # First, check whether the Python API is available
            python_api_available = False
            python_module = tool_info.get("python_module")
            if python_module:
                try:
                    __import__(python_module)
                    python_api_available = True
                except ImportError:
                    pass

            # Find CLI tool path
            tool_path = self._find_external_tool(tool_name, tool_info)

            if python_api_available or tool_path:
                if python_api_available and tool_path:
                    self.log(f"  ✓ {display_name} - Python API + CLI ({tool_path})")
                elif python_api_available:
                    self.log(f"  ✓ {display_name} - Python API available")
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
                install_info = tool_info["install_info"].get(self.os_type, "Please refer to the official documentation for installation")
                self.log(f"  ✗ {display_name} (not found) - {description}")
                self.log(f"      Install method: {install_info}")
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
        Find an external tool path.

        Args:
            tool_name: Tool name
            tool_info: Tool configuration info

        Returns:
            str: Tool path, or None if not found
        """
        import shutil

        # 1. System PATH
        path = shutil.which(tool_name)
        if path:
            return path

        # 2. Platform-specific search paths
        if self.os_type and self.os_type in tool_info.get("search_paths", {}):
            for search_path in tool_info["search_paths"][self.os_type]:
                expanded_path = os.path.expanduser(search_path)
                if os.path.isfile(expanded_path) and os.access(expanded_path, os.X_OK):
                    return expanded_path

        # 3. Conda environment
        if 'CONDA_PREFIX' in os.environ:
            conda_path = os.path.join(os.environ['CONDA_PREFIX'], 'bin', tool_name)
            if os.path.exists(conda_path):
                return conda_path

        return None

    def setup_surface_analysis(self) -> bool:
        """
        One-click setup for the surface analysis environment.

        Includes:
        1. Install Open3D
        2. Install scikit-image
        3. Check APBS (provide installation guidance)

        Returns:
            bool: Whether everything succeeded
        """
        self.log("\n" + "=" * 60)
        self.log("🔧 Setting up surface analysis environment (MaSIF-style)")
        self.log("=" * 60)

        success = True

        # 1. Install Open3D
        try:
            import open3d
            self.log("  ✓ Open3D is installed")
        except ImportError:
            self.log("  Installing Open3D...")
            if self._install_pip_package("open3d"):
                self.log("  ✅ Open3D installation succeeded")
            else:
                self.log("  ⚠️ Open3D installation failed; the built-in fallback will be used")
                success = False

        # 2. Install scikit-image
        try:
            import skimage
            self.log("  ✓ scikit-image is installed")
        except ImportError:
            self.log("  Installing scikit-image...")
            if self._install_pip_package("scikit-image"):
                self.log("  ✅ scikit-image installation succeeded")
            else:
                self.log("  ⚠️ scikit-image installation failed; a simplified algorithm will be used")

        # 3. Check external tools
        self.log("\n  Checking external tools...")
        tools = self.check_external_tools()
        apbs_available = tools.get("apbs", {}).get("available", False)

        if not apbs_available:
            self.log("\n  💡 APBS is not installed; electrostatic potential will use a Coulomb approximation")
            self.log(f"     Recommended install: {tools.get('apbs', {}).get('install_info', 'See the official documentation')}")

        if success:
            self.log("\n✅ Surface analysis environment setup completed")
        else:
            self.log("\n⚠️ Surface analysis is available (fallback mode); install Open3D/APBS for preferred behavior")

        return success

    def install_optional_package(self, package_name: str) -> bool:
        """
        Install an optional package.

        Args:
            package_name: Package display name

        Returns:
            bool: Whether installation succeeded
        """
        pip_name = None
        for _, display_name, pip, _ in OPTIONAL_PACKAGES:
            if display_name == package_name:
                pip_name = pip
                break

        if not pip_name:
            self.log(f"  ✗ Unknown optional package: {package_name}")
            return False

        return self._install_pip_package(pip_name)

    def install_all_optional_packages(self) -> Dict[str, bool]:
        """
        Install all optional packages.

        Returns:
            Dict[str, bool]: {Package name: whether installation succeeded}
        """
        self.log("\n📦 Installing optional Python packages...")
        results = {}

        for import_name, display_name, pip_name, description in OPTIONAL_PACKAGES:
            try:
                __import__(import_name)
                self.log(f"  ✓ {display_name} is installed")
                results[display_name] = True
            except ImportError:
                self.log(f"  Installing {display_name} ({description})...")
                if self._install_pip_package(pip_name):
                    self.log(f"  ✅ {display_name} installation succeeded")
                    results[display_name] = True
                else:
                    self.log(f"  ✗ {display_name} installation failed")
                    results[display_name] = False

        return results

    # ========== HMM File Management ==========

    def get_hmm_data_dir(self) -> str:
        """
        Get the HMM data directory path (created lazily, with permission error handling).

        Returns:
            str: Data directory path. If it cannot be created, a temporary directory is returned.
        """
        # Prefer using the plugin directory's "data" folder
        plugin_dir = os.path.dirname(os.path.abspath(__file__))
        data_dir = os.path.join(plugin_dir, "data")

        # Try to create the plugin directory's "data" folder
        if not os.path.exists(data_dir):
            try:
                os.makedirs(data_dir, exist_ok=True)
                return data_dir
            except PermissionError:
                self.log(f"  ⚠️ Cannot create plugin data directory {data_dir}: insufficient permissions")
            except OSError as e:
                self.log(f"  ⚠️ Cannot create plugin data directory {data_dir}: {e}")
        else:
            return data_dir

        # Fall back to user directory ~/.glint/data
        user_data_dir = os.path.join(os.path.expanduser("~"), ".glint", "data")
        try:
            os.makedirs(user_data_dir, exist_ok=True)
            return user_data_dir
        except PermissionError:
            self.log(f"  ⚠️ Cannot create user data directory {user_data_dir}: insufficient permissions")
        except OSError as e:
            self.log(f"  ⚠️ Cannot create user data directory {user_data_dir}: {e}")

        # Final fallback to a temporary directory
        temp_data_dir = os.path.join(tempfile.gettempdir(), "glint_data")
        try:
            os.makedirs(temp_data_dir, exist_ok=True)
            self.log(f"  ℹ️ Using temporary data directory: {temp_data_dir}")
            return temp_data_dir
        except Exception as e:
            self.log(f"  ❌ Unable to create any data directory: {e}")
            # Return the temp directory path; even if creation failed, let the caller handle it
            return temp_data_dir

    def setup_ec_analysis(self) -> bool:
        """
        One-click setup for the electrostatic complementarity (EC) analysis environment.

        Includes:
        1. Install PDB2PQR (Python API)
        2. Check/install APBS
        3. Ensure NumPy, SciPy, RDKit are available

        Returns:
            bool: Whether everything succeeded
        """
        self.log("\n" + "=" * 60)
        self.log("🔧 Setting up electrostatic complementarity (EC) analysis environment")
        self.log("=" * 60)

        success = True

        # 1. Check core dependencies
        self.log("\n  Checking core dependencies...")
        core_deps = [("numpy", "NumPy"), ("scipy", "SciPy"), ("rdkit", "RDKit")]
        for import_name, display_name in core_deps:
            try:
                __import__(import_name)
                self.log(f"  ✓ {display_name} is installed")
            except ImportError:
                self.log(f"  ✗ {display_name} is not installed - EC analysis requires this dependency")
                success = False

        # 2. Install PDB2PQR
        self.log("\n  Checking PDB2PQR...")
        try:
            import pdb2pqr
            self.log("  ✓ PDB2PQR Python API is installed")
        except ImportError:
            self.log("  Installing PDB2PQR...")
            if self._install_pip_package("pdb2pqr"):
                self.log("  ✅ PDB2PQR installation succeeded")
            else:
                self.log("  ⚠️ PDB2PQR installation failed; a simplified PQR conversion will be used")

        # 3. Check APBS
        self.log("\n  Checking APBS...")
        apbs_available = False

        # Check Python API
        try:
            import apbs
            self.log("  ✓ APBS Python API is installed")
            apbs_available = True
        except ImportError:
            self.log("  ✗ APBS Python API is not installed")

        # Check CLI tool
        import shutil
        apbs_path = shutil.which('apbs')
        if apbs_path:
            self.log(f"  ✓ APBS CLI is available: {apbs_path}")
            apbs_available = True
        else:
            self.log("  ✗ APBS CLI not found")

        if not apbs_available:
            self.log("\n  ⚠️ APBS is not installed; attempting automatic installation...")
            if self._install_pip_package("apbs"):
                self.log("  ✅ APBS (pip) installation succeeded")
                apbs_available = True
            else:
                self.log("  ⚠️ APBS pip installation failed; EC analysis will not be able to run")
                self.log("     Recommended installation methods:")
                if self.os_type == "macOS":
                    self.log("     - pip install apbs")
                    self.log("     - or brew install brewsci/bio/apbs")
                elif self.os_type == "Linux":
                    self.log("     - pip install apbs")
                    self.log("     - or apt install apbs")
                else:
                    self.log("     - pip install apbs")
                    self.log("     - or download from https://www.poissonboltzmann.org/")
            success = False

        # 4. Summary
        if success:
            self.log("\n✅ EC analysis environment setup completed")
            self.log("   You can use the following features:")
            self.log("   - calculate_ligand_ec: Protein-ligand EC analysis")
            self.log("   - analyze_ternary_ec: Ternary complex (molecular glue) EC analysis")
            self.log("   - compare_ligand_ec: Multi-ligand EC comparison (SAR analysis)")
            self.log("   - calculate_ec_hotspots: EC hotspot identification")
        else:
            self.log("\n⚠️ EC analysis environment is partially available; full functionality requires manual configuration")

        return success

    # ========== Installation Guide ==========

    def get_conda_install_url(self) -> str:
        """
        Get the Miniconda download URL.

        Returns:
            str: Download URL
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

    # ========== Automatic Installation ==========

    def auto_install_conda_env(self) -> bool:
        """
        Automatically create the conda environment.

        Returns:
            bool: True if succeeded
        """
        if not self.check_conda():
            self.log("✗ Cannot auto-create environment: Conda is not installed")
            return False

        if self.check_conda_env(ENV_NAME):
            self.log(f"✓ Environment '{ENV_NAME}' already exists, skipping creation")
            return True

        self.log(f"\n🔧 Creating Conda environment '{ENV_NAME}'...")

        try:
            result = subprocess.run(
                ["conda", "create", "-n", ENV_NAME, f"python={PYTHON_VERSION}", "-y"],
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode == 0:
                self.log(f"✅ Environment '{ENV_NAME}' created successfully")
                return True
            else:
                self.log(f"✗ Environment creation failed: {result.stderr}")
                return False
        except subprocess.TimeoutExpired:
            self.log("✗ Environment creation timed out")
            return False
        except Exception as e:
            self.log(f"✗ Environment creation error: {e}")
            return False

    def _install_pip_package(self, package_name: str) -> bool:
        """
        Install a package via pip (fallback).
        """
        try:
            # Use pip from the current Python environment
            cmd = [sys.executable, "-m", "pip", "install", package_name, "--quiet", "--disable-pip-version-check"]
            self.log(f"  Running pip: {' '.join(cmd)}")
            subprocess.check_call(cmd)
            return True
        except Exception as e:
            self.log(f"  ✗ Pip install failed: {e}")
            return False

    def auto_install_dependencies(self) -> Dict[str, bool]:
        """
        Automatically install missing dependencies.

        Returns:
            Dict[str, bool]: {Dependency name: whether installation succeeded}
        """
        use_conda = self.check_conda()
        if not use_conda:
            self.log("⚠️ Conda is not installed; will try installing via pip")

        status = self.check_all_dependencies()
        missing = [name for name, avail in status.items() if not avail]

        if not missing:
            self.log("✅ All dependencies are installed; nothing to do")
            return status

        self.log("\\n📦 Starting automatic installation of missing dependencies...")

        install_results = {}

        # Batch install Python packages
        py_packages = missing
        if py_packages:
            self.log(f"\\n  Installing Python packages: {', '.join(py_packages)}")

            # Try conda install first (only for Python libraries that actually exist on conda)
            if use_conda:
                packages_to_install = []
                for import_name, display_name, pip_name in REQUIRED_PACKAGES:
                    if display_name in py_packages:
                        # Meeko exists on conda-forge, but Vina is only provided via pip
                        if display_name == "Meeko":
                            packages_to_install.append("meeko")
                        else:
                            packages_to_install.append(pip_name)

                if packages_to_install:
                    try:
                        cmd = ["conda", "install", "-c", "conda-forge", "-y"] + packages_to_install
                        self.log(f"  Running: {' '.join(cmd)}")

                        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

                        if result.returncode == 0:
                            for pkg in py_packages:
                                install_results[pkg] = True
                                self.log(f"  ✅ {pkg} installed successfully (Conda)")
                        else:
                            self.log("  ⚠️ Conda install failed; trying pip...")
                    except Exception as e:
                        self.log(f"  ⚠️ Conda error ({e}); trying pip...")

            # Regardless of conda success or failure, try pip for all missing packages (to avoid arch/solver issues)
            for pkg in py_packages:
                pip_name = next((p[2] for p in REQUIRED_PACKAGES if p[1] == pkg), None)
                if not pip_name:
                    continue
                if self._install_pip_package(pip_name):
                    install_results[pkg] = True
                    self.log(f"  ✅ {pkg} installed successfully (Pip)")
                else:
                    install_results[pkg] = False
                    self.log(f"  ✗ {pkg} installation failed")
            else:
                # No Conda, use Pip directly
                for pkg in py_packages:
                    pip_name = next((p[2] for p in REQUIRED_PACKAGES if p[1] == pkg), None)
                    if pip_name:
                        if self._install_pip_package(pip_name):
                            install_results[pkg] = True
                            self.log(f"  ✅ {pkg} installed successfully (Pip)")
                        else:
                            install_results[pkg] = False
                            self.log(f"  ✗ {pkg} installation failed")

        # Summary
        success_count = sum(1 for v in install_results.values() if v)
        total_count = len(install_results)

        if success_count == total_count:
            self.log(f"\\n✅ All dependency installations completed ({success_count}/{total_count})")
        else:
            self.log(f"\\n⚠️  Some dependency installations failed ({success_count}/{total_count})")

        return install_results

    def _install_conda_package(self, package_name: str) -> bool:
        """
        Install a single conda package.

        Args:
            package_name: Package name

        Returns:
            bool: True if succeeded
        """
        try:
            result = subprocess.run(
                ["conda", "install", "-c", "conda-forge", package_name, "-y"],
                capture_output=True,
                text=True,
                timeout=300
            )
            return result.returncode == 0
        except (OSError, subprocess.SubprocessError):  # Subprocess invocation may fail
            return False

    def get_install_instructions(self) -> Dict[str, str]:
        """
        Get installation instructions.

        Returns:
            Dict[str, str]: Instructions dictionary
        """
        status = self.check_all_dependencies()
        missing_packages = [name for name, avail in status.items() if not avail]

        instructions = {}

        # Conda installation
        if not self.check_conda():
            instructions["conda"] = f"""
📥 Install Miniconda:

1. Download: {self.get_conda_install_url()}
2. After installation, run: source ~/.zshrc  (or ~/.bashrc)
3. Verify: conda --version
"""

        # Environment creation
        if not self.check_conda_env():
            instructions["environment"] = f"""
🔧 Create Conda environment:

conda create -n {ENV_NAME} python={PYTHON_VERSION} -y
conda activate {ENV_NAME}
"""

        # Dependencies installation
        if missing_packages:
            instructions["dependencies"] = """
📦 Install all dependencies (Vina is optional; pip is recommended):

# Core stack (without Vina)
conda install -c conda-forge rdkit scipy matplotlib pillow \\
    numpy pandas seaborn pyqt openbabel -y

# Optional: install Vina + Meeko within the activated environment
pip install vina meeko

"""

        # Usage instructions
        instructions["usage"] = f"""
✅ Usage:

1. Activate environment:
   conda activate {ENV_NAME}

2. Start PyMOL (from the activated environment)

3. Load plugin:
   run /path/to/glint/__init__.py

4. Open GUI:
   molstruct_gui
"""

        return instructions

    # ========== Full Check ==========

    def run_full_check(self) -> Dict[str, any]:
        """
        Run the full environment check.

        Returns:
            Dict: Summary of check results
        """
        self.log("=" * 60)
        self.log("⚡ GLINT Environment Check")
        self.log("=" * 60)

        # 1. Detect OS
        self.log("\n🖥️  System information:")
        os_type, os_arch = self.detect_os()

        # 2. Check Conda
        self.log("\n🐍 Conda:")
        has_conda = self.check_conda()

        # 3. Check environment
        if has_conda:
            self.log("\n📁 Conda environment:")
            has_env = self.check_conda_env(ENV_NAME)
        else:
            has_env = False

        # 4. Check all dependencies
        dep_status = self.check_all_dependencies()

        # 5. Auto install (if enabled)
        if self.auto_install and has_conda:
            self.log("\n" + "=" * 60)
            self.log("🚀 Auto-install mode")
            self.log("=" * 60)

            # Create environment (if it does not exist)
            if not has_env:
                has_env = self.auto_install_conda_env()

            # Install dependencies
            install_results = self.auto_install_dependencies()

            # Re-check dependency status
            self.log("\n🔍 Re-checking dependency status...")
            dep_status = self.check_all_dependencies()

        # 6. Test GUI functionality
        gui_ok = self.test_gui_functionality()

        # 7. Generate installation instructions or summary
        self.log("\n" + "=" * 60)

        all_ok = has_conda and has_env and all(dep_status.values()) and gui_ok

        if all_ok:
            self.log("✅ Environment is perfectly configured! All features are available.")
        else:
            if self.auto_install:
                self.log("⚠️  Some dependencies failed to install or require manual configuration")
            else:
                self.log("⚠️  Environment setup is required")
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


# ========== Convenience Functions ==========

def check_environment(log_callback=None, auto_install=False) -> Dict[str, any]:
    """
    Check the environment (convenience function).

    Args:
        log_callback: Log callback function
        auto_install: Whether to auto-install missing dependencies

    Returns:
        Dict: Check results
    """
    checker = EnvironmentChecker(log_callback, auto_install=auto_install)
    return checker.run_full_check()


def get_dependency_status() -> Dict[str, bool]:
    """
    Get dependency status (convenience function).

    Returns:
        Dict[str, bool]: {Dependency name: available or not}
    """
    checker = EnvironmentChecker(log_callback=None)
    return checker.check_all_dependencies()


def setup_surface_analysis(log_callback=None) -> bool:
    """
    One-click setup for surface analysis (convenience function).

    Includes:
    1. Install Open3D
    2. Install scikit-image
    3. Check MSMS/APBS

    Args:
        log_callback: Log callback function

    Returns:
        bool: Whether everything succeeded
    """
    checker = EnvironmentChecker(log_callback=log_callback)
    return checker.setup_surface_analysis()


def setup_ec_analysis(log_callback=None) -> bool:
    """
    One-click setup for electrostatic complementarity (EC) analysis (convenience function).

    Includes:
    1. Install PDB2PQR (Python API)
    2. Check/install APBS
    3. Ensure NumPy, SciPy, RDKit are available

    Args:
        log_callback: Log callback function

    Returns:
        bool: Whether everything succeeded
    """
    checker = EnvironmentChecker(log_callback=log_callback)
    return checker.setup_ec_analysis()


def check_ec_analysis_deps(log_callback=None) -> Dict[str, bool]:
    """
    Check dependency status for EC analysis (convenience function).

    Args:
        log_callback: Log callback function

    Returns:
        Dict[str, bool]: {Dependency name: available or not}
    """
    checker = EnvironmentChecker(log_callback=log_callback)

    status = {}

    # Core Python packages
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

    # EC-specific packages
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

    # External tools (CLI)
    tools = checker.check_external_tools()
    for tool_name in ["apbs", "pdb2pqr"]:
        tool_info = tools.get(tool_name, {})
        cli_name = f"{tool_info.get('display_name', tool_name)} CLI"
        status[cli_name] = tool_info.get("path") is not None

    return status


def check_surface_analysis_deps(log_callback=None) -> Dict[str, bool]:
    """
    Check dependency status for surface analysis (convenience function).

    Args:
        log_callback: Log callback function

    Returns:
        Dict[str, bool]: {Dependency name: available or not}
    """
    checker = EnvironmentChecker(log_callback=log_callback)

    status = {}

    # Python packages
    for import_name, display_name, _, _ in OPTIONAL_PACKAGES:
        if display_name in ["Open3D", "scikit-image", "trimesh"]:
            try:
                __import__(import_name)
                status[display_name] = True
            except ImportError:
                status[display_name] = False

    # External tools
    tools = checker.check_external_tools()
    for tool_name, tool_info in tools.items():
        status[tool_info.get("display_name", tool_name)] = tool_info.get("available", False)

    return status


# ========== Command Line Interface ==========

if __name__ == "__main__":
    # Command-line mode
    import argparse

    parser = argparse.ArgumentParser(description="GLINT environment check and auto-configuration tool")
    parser.add_argument(
        "--auto-install",
        action="store_true",
        help="Automatically install missing dependencies (requires Conda)"
    )

    args = parser.parse_args()

    checker = EnvironmentChecker(auto_install=args.auto_install)
    result = checker.run_full_check()

    sys.exit(0 if result["all_ok"] else 1)


def is_in_conda_env() -> bool:
    """
    Check whether the current process is running inside a conda environment.

    Returns:
        bool: True if running in a conda environment
    """
    return bool(os.environ.get('CONDA_PREFIX', ''))


def is_first_run() -> bool:
    """
    Check whether this is the first run.

    Returns:
        bool: True if this is the first run
    """
    return not os.path.exists(_FIRST_RUN_MARKER)


def mark_initialized():
    """
    Mark as initialized (with permission error handling).

    Attempts to create an initialization marker file in the user's home directory.
    If it fails (permission issues, etc.), ignore silently and do not affect normal plugin usage.
    """
    try:
        # Ensure the parent directory exists
        marker_dir = os.path.dirname(_FIRST_RUN_MARKER)
        if marker_dir and not os.path.exists(marker_dir):
            try:
                os.makedirs(marker_dir, exist_ok=True)
            except (PermissionError, OSError):
                # Cannot create directory; ignore silently
                return

        with open(_FIRST_RUN_MARKER, 'w') as f:
            f.write('initialized')
    except PermissionError:
        # Insufficient permissions; ignore silently
        pass
    except OSError:
        # Other OS errors; ignore silently
        pass
    except Exception:
        # Any other error; ignore silently
        pass


def _quick_check_deps(skip_pyqt=True) -> Tuple[bool, List[str]]:
    """
    Quickly check dependencies (no log output).

    Args:
        skip_pyqt: Whether to skip the PyQt5 check (to avoid potential crashes during import)

    Returns:
        (all_ok, missing_list)
    """
    import shutil
    missing = []

    # Check required Python packages (skip PyQt5 to avoid potential crashes)
    for import_name, display_name, _ in REQUIRED_PACKAGES:
        if skip_pyqt and import_name == "PyQt5":
            continue
        try:
            __import__(import_name)
        except ImportError:
            missing.append(display_name)

    # Check command-line tools (currently empty)
    for cmd, name in REQUIRED_COMMANDS:
        if not shutil.which(cmd):
            missing.append(name)

    return len(missing) == 0, missing


def _print_setup_guide(missing: List[str]):
    """Print the full environment setup guide."""
    print("\n" + "=" * 60)
    print("🚨 GLINT - Environment Setup Required")
    print("=" * 60)

    if missing:
        print(f"\n❌ Missing: {', '.join(missing)}")

    print("\n📖 One-click setup (run in terminal):")
    print("─" * 60)
    print("   bash /path/to/GLINT/install.sh")
    print("─" * 60)

    print("\n📖 Or manual setup:")
    print(f"""
  conda create -n {ENV_NAME} python={PYTHON_VERSION} -y
  conda activate {ENV_NAME}
  conda install -c conda-forge rdkit scipy matplotlib pillow \\
      numpy pandas seaborn pyqt -y
""")
    print("💡 After setup, always start PyMOL from the activated environment:")
    print(f"   conda activate {ENV_NAME} && pymol")
    print("=" * 60 + "\n")


def _print_missing_deps_hint(missing: List[str]):
    """Print a short hint for missing dependencies."""
    print("\n" + "=" * 60)
    print("⚠️  GLINT - Missing Dependencies")
    print("=" * 60)

    if missing:
        print(f"\n❌ Missing: {', '.join(missing)}")

    # Build package name mapping
    pkg_map = {p[1]: p[2] for p in REQUIRED_PACKAGES}

    conda_names = [pkg_map.get(m, m.lower()) for m in missing]

    if conda_names:
        print("\n💡 Install missing packages:")
        print(f"   conda activate {ENV_NAME}")
        print(f"   conda install -c conda-forge {' '.join(conda_names)} -y")
    print("=" * 60 + "\n")


def ensure_dependencies(silent=False) -> bool:
    """
    Check and ensure all dependencies are available.

    Args:
        silent: If True, do not print any information

    Returns:
        bool: True if all dependencies are available
    """
    # Quick dependency check
    all_ok, missing = _quick_check_deps()

    if all_ok:
        # All dependencies are present; mark initialized
        mark_initialized()
        return True

    # Dependencies are missing
    if silent:
        return False

    # Check whether this is the first run
    if is_first_run():
        # First run: show full setup guide
        _print_setup_guide(missing)
    else:
        # Not the first run: show short hint
        _print_missing_deps_hint(missing)

    return False
