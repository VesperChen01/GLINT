#!/bin/bash
# release.sh - GlueTK Release Script (Merged Version)
# Combines version update, rebuild installer, and DMG packaging
#
# Usage: ./release.sh <new_version>
# Example: ./release.sh v0.1.6-beta

set -e

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
APP_NAME="GlueTK Installer.app"
CONTENTS="${APP_NAME}/Contents"
MACOS="${CONTENTS}/MacOS"
RESOURCES="${CONTENTS}/Resources"
BUILD_DIR="build_dmg"
GLUETK_LIB="gluetk"

# Check conda environment
echo -e "${BLUE}[0/6] Checking conda environment...${NC}"
if [ -z "$CONDA_PREFIX" ]; then
    echo -e "${YELLOW}⚠️  Conda environment not activated${NC}"
    echo "Please activate the gluetk environment:"
    echo "  conda activate gluetk"
    exit 1
fi
echo -e "${GREEN}✅ Conda environment: $CONDA_PREFIX${NC}"

if [ -z "$1" ]; then
    CURRENT_VERSION=$(python3 -c "from gluetk._version import __version__; print(__version__)")
    echo ""
    echo -e "${BLUE}GlueTK Release Script (Merged)${NC}"
    echo "================================"
    echo ""
    echo "Current version: $CURRENT_VERSION"
    echo ""
    echo "Usage: $0 <new_version>"
    echo "Example: $0 v0.1.6-beta"
    echo ""
    echo "This script will:"
    echo "  1. Update version in gluetk/_version.py"
    echo "  2. Rebuild GlueTK Installer.app"
    echo "  3. Create DMG: GlueTK_Installer_<version>.dmg"
    echo "  4. Copy to PyMOL startup directory"
    echo "  5. Clean up build artifacts"
    echo ""
    exit 0
fi

NEW_VERSION="$1"
OLD_VERSION=$(python3 -c "from gluetk._version import __version__; print(__version__)")

echo ""
echo -e "${BLUE}🚀 GlueTK Release: $OLD_VERSION → $NEW_VERSION${NC}"
echo ""

# ============================================================================
# STEP 1: Update version in _version.py (from update_version.sh)
# ============================================================================
echo -e "${BLUE}[1/6] Updating version in source files...${NC}"
sed -i '' "s/__version__ = \".*\"/__version__ = \"$NEW_VERSION\"/" gluetk/_version.py
echo -e "${GREEN}   ✅ Updated gluetk/_version.py${NC}"

# ============================================================================
# STEP 2: Rebuild GlueTK Installer.app (from rebuild_installer.sh)
# ============================================================================
echo -e "${BLUE}[2/6] Rebuilding GlueTK Installer.app...${NC}"

# 2.1 Create directories
rm -rf "${APP_NAME}"
mkdir -p "${MACOS}"
mkdir -p "${RESOURCES}"
echo "   Creating app structure..."

# 2.2 Copy icon if exists
if [ -f "gluetk/assets/AppIcon.icns" ]; then
    cp "gluetk/assets/AppIcon.icns" "${RESOURCES}/AppIcon.icns"
    cp "gluetk/assets/logo.png" "${RESOURCES}/AppIcon.png" 2>/dev/null || true
    echo -e "${GREEN}   ✅ Icon copied${NC}"
else
    echo -e "${YELLOW}   ⚠️  Icon not found in gluetk/assets/${NC}"
fi

# 2.3 Create Info.plist
cat > "${CONTENTS}/Info.plist" << EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleExecutable</key>
    <string>launcher</string>
    <key>CFBundleIconFile</key>
    <string>AppIcon</string>
    <key>CFBundleIdentifier</key>
    <string>com.vesper.gluetk.installer</string>
    <key>CFBundleName</key>
    <string>GlueTK Installer</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleShortVersionString</key>
    <string>$NEW_VERSION</string>
</dict>
</plist>
EOF
echo "   Info.plist created"

# 2.4 Create launcher script
cat > "${MACOS}/launcher" << 'LAUNCHER_EOF'
#!/bin/bash
# GlueTK Installer Launcher (Tk GUI via conda python)

CONDA_EXE=""
if command -v conda &> /dev/null; then
    CONDA_EXE=$(command -v conda)
else
    for p in "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" "/opt/miniconda3/bin/conda" "/opt/anaconda3/bin/conda" "/usr/local/bin/conda" "/opt/homebrew/bin/conda" "$HOME/opt/miniconda3/bin/conda"; do
        if [ -x "$p" ]; then
            CONDA_EXE="$p"
            break
        fi
    done
fi

if [ -z "$CONDA_EXE" ]; then
  osascript -e 'display alert "Conda Not Found" message "Please install Miniconda (conda) first, then re-run GlueTK Installer."'
  exit 1
fi

eval "$($CONDA_EXE shell.bash hook)"

CONDA_BASE="$(conda info --base 2>/dev/null)"
if [ ! -d "$CONDA_BASE" ]; then
  osascript -e 'display alert "Conda Base Not Found" message "Conda is installed but CONDA_BASE is invalid. Please check your conda installation."'
  exit 1
fi

PYTHON="${CONDA_BASE}/bin/python"
if [ ! -x "$PYTHON" ]; then
  if command -v python &> /dev/null; then
    PYTHON=$(command -v python)
  else
    osascript -e 'display alert "Python Not Found" message "Cannot find a suitable python interpreter."'
    exit 1
  fi
fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RESOURCES_DIR="$(dirname "$DIR")/Resources"
INSTALLER_SCRIPT="${RESOURCES_DIR}/GlueTK_Installer.py"

"$PYTHON" "$INSTALLER_SCRIPT"
LAUNCHER_EOF

chmod +x "${MACOS}/launcher"
echo -e "${GREEN}   ✅ Launcher script created${NC}"

# 2.5 Copy Tk GUI installer (from rebuild_installer.sh lines 108-747)
cat > "${RESOURCES}/GlueTK_Installer.py" << 'PYTHON_INSTALLER_EOF'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GlueTK macOS Installer - 图形化安装程序
Professional GUI installer for macOS, similar to Windows version
"""

import os
import sys
import subprocess
import shutil
import threading
import webbrowser
from pathlib import Path

try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
except ImportError:
    print("Error: tkinter not found")
    sys.exit(1)

# 配置
ENV_NAME = "gluetk"
PYTHON_VERSION = "3.10"
DEFAULT_INSTALL_PATH = os.path.join(os.path.expanduser("~"), ".pymol", "startup", "gluetk")

# Conda 依赖包
CONDA_PACKAGES = [
    "rdkit", "scipy", "matplotlib", "pillow", "numpy=1.26.4",
    "pandas", "seaborn", "pyqt", "openbabel", "pymol-open-source",
    "meeko", "vina", "haddock_biobb", "scikit-image",
    "pdb2pqr",
]

# Pip 包
PIP_PACKAGES = ["requests", "open3d"]


def get_gluetk_source_dir():
    """获取 GlueTK 源码目录"""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    bundled_dir = os.path.join(script_dir, "gluetk")
    if os.path.isdir(bundled_dir):
        return bundled_dir

    try:
        repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(script_dir))))
        gluetk_dir = os.path.join(repo_dir, "gluetk")
        if os.path.isdir(gluetk_dir):
            return gluetk_dir
    except:
        pass

    return None


class InstallerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GlueTK Installer")
        self.root.geometry("750x780")
        self.root.minsize(700, 750)
        self.root.resizable(True, True)

        self._configure_styles()
        self._set_icon()

        self.install_path = tk.StringVar(value=DEFAULT_INSTALL_PATH)
        self.create_shortcut = tk.BooleanVar(value=True)
        self.install_deps = tk.BooleanVar(value=True)
        self.conda_ok = False
        self.env_ok = False
        self.conda_exe = None
        self.env_path = None

        self._build_ui()
        self._center_window()

        self.root.update()
        self.root.lift()
        self.root.attributes('-topmost',True)
        self.root.after_idle(self.root.attributes,'-topmost',False)

        self.root.after(500, lambda: threading.Thread(target=self._check_environment, daemon=True).start())

    def _configure_styles(self):
        """Configure better fonts and styles for macOS"""
        style = ttk.Style()

        available_themes = style.theme_names()
        if 'aqua' in available_themes:
            style.theme_use('aqua')
        elif 'clam' in available_themes:
            style.theme_use('clam')

        self.title_font = ("Helvetica Neue", 22, "bold")
        self.subtitle_font = ("Helvetica Neue", 11)
        self.normal_font = ("Helvetica Neue", 10)
        self.small_font = ("Helvetica Neue", 9)
        self.mono_font = ("Menlo", 9)

        try:
            import tkinter.font as tkfont
            available_fonts = tkfont.families()

            if "Helvetica Neue" not in available_fonts:
                if "Arial" in available_fonts:
                    self.title_font = ("Arial", 22, "bold")
                    self.subtitle_font = ("Arial", 11)
                    self.normal_font = ("Arial", 10)
                    self.small_font = ("Arial", 9)

            if "Menlo" not in available_fonts:
                if "Courier New" in available_fonts:
                    self.mono_font = ("Courier New", 9)
        except:
            pass

        style.configure("TLabel", font=self.normal_font)
        style.configure("TButton", font=self.normal_font, padding=6)
        style.configure("TCheckbutton", font=self.normal_font)
        style.configure("TEntry", font=self.normal_font)
        style.configure("TLabelframe", font=self.normal_font)
        style.configure("TLabelframe.Label", font=self.normal_font)

        self.root.configure(bg='#f0f0f0')

    def _set_icon(self):
        """设置窗口图标"""
        try:
            icon_path = os.path.join(os.path.dirname(__file__), "gluetk", "assets", "logo.png")
            if os.path.exists(icon_path):
                img = tk.PhotoImage(file=icon_path)
                self.root.iconphoto(True, img)
        except:
            pass

    def _center_window(self):
        """居中窗口"""
        self.root.update_idletasks()
        w = self.root.winfo_width()
        h = self.root.winfo_height()
        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        x = (sw - w) // 2
        y = (sh - h) // 2
        self.root.geometry(f"{w}x{h}+{x}+{y}")

    def _build_ui(self):
        """Build user interface with improved fonts"""
        main = ttk.Frame(self.root, padding=20)
        main.pack(fill=tk.BOTH, expand=True)

        title_frame = ttk.Frame(main)
        title_frame.pack(fill=tk.X, pady=(0, 15))

        title = ttk.Label(title_frame, text="🧬 GlueTK Installer",
                         font=self.title_font)
        title.pack()

        subtitle = ttk.Label(title_frame,
                            text="PyMOL Plugin for Molecular Glue Analysis",
                            font=self.subtitle_font, foreground="#666666")
        subtitle.pack(pady=(5, 0))

        version = ttk.Label(title_frame, text="Version: v0.1.6-beta",
                           font=self.small_font, foreground="#888888")
        version.pack(pady=(3, 0))

        status_frame = ttk.LabelFrame(main, text=" Environment Status ", padding=12)
        status_frame.pack(fill=tk.X, pady=(0, 15))

        conda_row = ttk.Frame(status_frame)
        conda_row.pack(fill=tk.X, pady=5)
        ttk.Label(conda_row, text="Conda:", font=self.normal_font, width=22).pack(side=tk.LEFT)
        self.conda_status = ttk.Label(conda_row, text="⏳ Checking...",
                                      font=self.normal_font, foreground="#E67E22")
        self.conda_status.pack(side=tk.LEFT, padx=10)
        self.conda_install_btn = ttk.Button(conda_row, text="Download Miniconda",
                                            command=self._install_conda, state=tk.DISABLED)
        self.conda_install_btn.pack(side=tk.RIGHT)

        env_row = ttk.Frame(status_frame)
        env_row.pack(fill=tk.X, pady=5)
        ttk.Label(env_row, text=f"Conda Env '{ENV_NAME}':",
                  font=self.normal_font, width=22).pack(side=tk.LEFT)
        self.env_status = ttk.Label(env_row, text="⏳ Checking...",
                                    font=self.normal_font, foreground="#E67E22")
        self.env_status.pack(side=tk.LEFT, padx=10)

        path_frame = ttk.LabelFrame(main, text=" Installation Path ", padding=12)
        path_frame.pack(fill=tk.X, pady=(0, 15))

        path_row = ttk.Frame(path_frame)
        path_row.pack(fill=tk.X)
        self.path_entry = ttk.Entry(path_row, textvariable=self.install_path,
                                    width=60, font=self.normal_font)
        self.path_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(path_row, text="Browse...", command=self._browse_path).pack(side=tk.RIGHT, padx=(10, 0))

        path_note = ttk.Label(path_frame,
                             text="📌 GlueTK will be installed here. PyMOL loads plugins from ~/.pymol/startup/",
                             font=self.small_font, foreground="#888888")
        path_note.pack(anchor=tk.W, pady=(8, 0))

        opts_frame = ttk.LabelFrame(main, text=" Options ", padding=12)
        opts_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Checkbutton(opts_frame, text="Install/Update dependencies (conda packages + PyMOL)",
                       variable=self.install_deps).pack(anchor=tk.W, pady=3)
        ttk.Checkbutton(opts_frame, text="Create desktop app (GlueTK.app)",
                       variable=self.create_shortcut).pack(anchor=tk.W, pady=3)

        progress_frame = ttk.LabelFrame(main, text=" Progress ", padding=12)
        progress_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))

        self.progress = ttk.Progressbar(progress_frame, mode="determinate", length=500)
        self.progress.pack(fill=tk.X, pady=(0, 10))

        log_frame = ttk.Frame(progress_frame)
        log_frame.pack(fill=tk.BOTH, expand=True)

        self.log_text = tk.Text(log_frame, height=12, state=tk.DISABLED,
                               font=self.mono_font, bg="#1a1a2e", fg="#eaeaea",
                               insertbackground="#ffffff", selectbackground="#3d5a80",
                               relief=tk.FLAT, padx=10, pady=8)
        scrollbar = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        btn_frame = ttk.Frame(main)
        btn_frame.pack(fill=tk.X, pady=(5, 0))

        self.install_btn = ttk.Button(btn_frame, text="🚀 Install GlueTK",
                                      command=self._start_install, width=20)
        self.install_btn.pack(side=tk.RIGHT, padx=(10, 0))

        ttk.Button(btn_frame, text="Cancel", command=self.root.quit, width=12).pack(side=tk.RIGHT)

    def _log(self, msg):
        """写入日志"""
        self.log_text.configure(state=tk.NORMAL)
        self.log_text.insert(tk.END, msg + "\n")
        self.log_text.see(tk.END)
        self.log_text.configure(state=tk.DISABLED)
        self.root.update()

    def _browse_path(self):
        """浏览安装路径"""
        path = filedialog.askdirectory(title="Select Installation Directory",
                                       initialdir=os.path.dirname(self.install_path.get()))
        if path:
            self.install_path.set(path)

    def _check_environment(self):
        """检查环境"""
        self._log("Checking environment...")

        self.conda_exe = self._find_conda()

        if self.conda_exe:
            try:
                result = subprocess.run([self.conda_exe, "--version"],
                                       capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    version = result.stdout.strip()
                    self.root.after(0, lambda: self.conda_status.configure(
                        text=f"✅ {version}", foreground="green"))
                    self.conda_ok = True
                    self.root.after(0, lambda: self.conda_install_btn.configure(state=tk.DISABLED))
                    self._log(f"  Conda found: {self.conda_exe}")
                    self._log(f"  Version: {version}")
                else:
                    raise Exception("conda command failed")
            except Exception as e:
                self.root.after(0, lambda: self.conda_status.configure(
                    text="❌ Error", foreground="red"))
                self._log(f"  Conda error: {e}")
        else:
            self.root.after(0, lambda: self.conda_status.configure(
                text="❌ Not found", foreground="red"))
            self.root.after(0, lambda: self.conda_install_btn.configure(state=tk.NORMAL))
            self._log("  Conda: Not found")
            self._log("  Please install Miniconda first!")
            return

        try:
            base_result = subprocess.run([self.conda_exe, "info", "--base"],
                                        capture_output=True, text=True, timeout=10)
            conda_base = base_result.stdout.strip() if base_result.returncode == 0 else ""

            if not conda_base:
                raise Exception("Cannot determine conda base")

            self.env_path = os.path.join(conda_base, "envs", ENV_NAME)

            if os.path.isdir(self.env_path):
                self.root.after(0, lambda: self.env_status.configure(
                    text="✅ Exists", foreground="green"))
                self.env_ok = True
                self._log(f"  Environment '{ENV_NAME}': {self.env_path}")
            else:
                self.root.after(0, lambda: self.env_status.configure(
                    text="⚠️ Will create", foreground="orange"))
                self._log(f"  Environment '{ENV_NAME}': Will be created")
        except Exception as e:
            self.root.after(0, lambda: self.env_status.configure(
                text="❓ Unknown", foreground="gray"))
            self._log(f"  Environment check failed: {e}")

        self._log("Ready to install.")

    def _find_conda(self):
        """查找 Conda 安装路径 (macOS)"""
        try:
            result = subprocess.run(["which", "conda"], capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                conda_path = result.stdout.strip().split("\n")
                if os.path.exists(conda_path):
                    return conda_path
        except:
            pass

        home = os.path.expanduser("~")
        candidates = [
            os.path.join(home, "miniconda3", "bin", "conda"),
            os.path.join(home, "anaconda3", "bin", "conda"),
            "/opt/miniconda3/bin/conda",
            "/opt/anaconda3/bin/conda",
            "/usr/local/bin/conda",
            "/opt/homebrew/bin/conda",
            "/opt/homebrew/Caskroom/miniconda/base/bin/conda"
        ]

        for path in candidates:
            if os.path.exists(path):
                return path

        return None

    def _install_conda(self):
        """打开 Miniconda 下载页面"""
        webbrowser.open("https://docs.conda.io/en/latest/miniconda.html")
        messagebox.showinfo("Install Miniconda",
                           "Please download and install Miniconda for macOS.\n\n"
                           "After installation, restart this installer.")

    def _start_install(self):
        """开始安装"""
        if not self.conda_ok:
            messagebox.showerror("Error", "Please install Miniconda first!")
            return

        self.install_btn.configure(state=tk.DISABLED)
        threading.Thread(target=self._do_install, daemon=True).start()

    def _do_install(self):
        """执行安装"""
        try:
            self.progress["value"] = 0

            self._log("\n" + "=" * 50)
            self._log("[1/5] Checking conda environment...")
            self._log("=" * 50)
            self.progress["value"] = 10

            if not self.env_path:
                try:
                    base_result = subprocess.run([self.conda_exe, "info", "--base"],
                                                capture_output=True, text=True, timeout=10)
                    conda_base = base_result.stdout.strip()
                    if conda_base:
                        self.env_path = os.path.join(conda_base, "envs", ENV_NAME)
                    else:
                        raise Exception("Cannot determine conda base")
                except Exception as e:
                    raise Exception(f"Failed to determine environment path: {e}")

            if not self.env_ok:
                self._log(f"  Creating environment at {self.env_path}...")
                result = subprocess.run(
                    [self.conda_exe, "create", "-p", self.env_path, f"python={PYTHON_VERSION}", "-y"],
                    capture_output=True, text=True, timeout=600
                )
                if result.returncode != 0:
                    self._log(f"  Error: {result.stderr}")
                    raise Exception("Failed to create conda environment")
                self._log("  ✅ Environment created")
            else:
                self._log(f"  Environment exists: {self.env_path}")

            self.progress["value"] = 20

            if self.install_deps.get():
                self._log("\n" + "=" * 50)
                self._log("[2/5] Installing dependencies...")
                self._log("=" * 50)
                self._log("  This may take 10-20 minutes, please wait...")

                pkg_str = " ".join(CONDA_PACKAGES)
                self._log(f"  Installing: {pkg_str}")

                cmd = [
                    self.conda_exe, "install",
                    "-p", self.env_path,
                    "-c", "conda-forge",
                    "-c", "schrodinger",
                    "-c", "haddocking",
                    "-y"
                ] + CONDA_PACKAGES

                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

                if result.returncode == 0:
                    self._log("  ✅ Conda packages installed")
                else:
                    self._log(f"  ⚠️ Some packages may have failed")
                    if result.stderr:
                        self._log(f"  {result.stderr[:500]}")

                self.progress["value"] = 50

                if PIP_PACKAGES:
                    self._log("\n  Installing pip packages...")
                    pip_cmd = [
                        self.conda_exe, "run", "-p", self.env_path,
                        "python", "-m", "pip", "install", "--quiet"
                    ] + PIP_PACKAGES

                    pip_result = subprocess.run(pip_cmd, capture_output=True, text=True, timeout=600)
                    if pip_result.returncode == 0:
                        self._log("  ✅ Pip packages installed")
                    else:
                        self._log("  ⚠️ Some pip packages may have failed")
            else:
                self._log("\n[2/5] Skipping dependencies (unchecked)")

            self.progress["value"] = 60

            self._log("\n" + "=" * 50)
            self._log("[3/5] Installing plugin files...")
            self._log("=" * 50)

            install_path = self.install_path.get()
            source_dir = get_gluetk_source_dir()

            if source_dir and os.path.isdir(source_dir):
                os.makedirs(install_path, exist_ok=True)

                if os.path.exists(install_path):
                    for item in os.listdir(install_path):
                        item_path = os.path.join(install_path, item)
                        try:
                            if os.path.isdir(item_path):
                                shutil.rmtree(item_path)
                            else:
                                os.remove(item_path)
                        except:
                            pass

                IGNORED = {'__pycache__', '.DS_Store', '.git', '.gitignore', '*.pyc'}
                for item in os.listdir(source_dir):
                    if item in IGNORED:
                        continue
                    src = os.path.join(source_dir, item)
                    dst = os.path.join(install_path, item)
                    try:
                        if os.path.isdir(src):
                            shutil.copytree(src, dst,
                                          ignore=shutil.ignore_patterns(*IGNORED))
                        else:
                            shutil.copy2(src, dst)
                    except Exception as e:
                        self._log(f"  ⚠️ Failed to copy {item}: {e}")

                self._log(f"  ✅ Copied to {install_path}")
            else:
                self._log(f"  ⚠️ Source directory not found: {source_dir}")

            self.progress["value"] = 80

            if self.create_shortcut.get():
                self._log("\n" + "=" * 50)
                self._log("[4/5] Creating desktop App...")
                self._log("=" * 50)
                self._create_desktop_app()
            else:
                self._log("\n[4/5] Skipping app creation")

            self.progress["value"] = 100

            self._log("\n" + "=" * 50)
            self._log("[5/5] ✅ Installation complete!")
            self._log("=" * 50)
            self._log("\nHow to use GlueTK:")
            self._log("  1. Double-click 'GlueTK.app' on Desktop")
            self._log("  2. Or run in PyMOL: gluetk_gui")

            self.root.after(0, lambda: messagebox.showinfo(
                "Success",
                "GlueTK installed successfully!\n\n"
                "You can now launch GlueTK from the Desktop shortcut."))

        except Exception as e:
            self._log(f"\n❌ Error: {e}")
            import traceback
            self._log(traceback.format_exc())
            self.root.after(0, lambda: messagebox.showerror("Error", str(e)))
        finally:
            self.root.after(0, lambda: self.install_btn.configure(state=tk.NORMAL))

    def _create_desktop_app(self):
        """Create macOS .app bundle"""
        desktop_app = os.path.expanduser("~/Desktop/GlueTK.app")
        if os.path.exists(desktop_app):
            shutil.rmtree(desktop_app)

        contents = os.path.join(desktop_app, "Contents")
        macos = os.path.join(contents, "MacOS")
        resources = os.path.join(contents, "Resources")
        os.makedirs(macos, exist_ok=True)
        os.makedirs(resources, exist_ok=True)

        icon_src = None
        c1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "AppIcon.icns")
        source_dir = get_gluetk_source_dir()
        if source_dir:
            c2 = os.path.join(source_dir, "assets", "AppIcon.icns")
            if os.path.exists(c1): icon_src = c1
            elif os.path.exists(c2): icon_src = c2

        if icon_src and os.path.exists(icon_src):
            shutil.copy2(icon_src, os.path.join(resources, "AppIcon.icns"))

        with open(os.path.join(contents, "Info.plist"), "w") as f:
            f.write('''<?xml version="1.0" encoding="UTF-8"?>
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
</plist>''')

        conda_exe_str = self.conda_exe or ""
        env_path_str = self.env_path or ""
        launcher_script = os.path.join(macos, "launcher")

        with open(launcher_script, "w") as f:
            f.write(f'''#!/bin/bash
# GlueTK Launcher (macOS .app)

CONDA_EXE="{conda_exe_str}"
ENV_PATH="{env_path_str}"

if [ -n "$CONDA_EXE" ] && [ -x "$CONDA_EXE" ] && [ -n "$ENV_PATH" ] && [ -d "$ENV_PATH" ]; then
  echo "Starting GlueTK via conda env: $ENV_PATH"
  "$CONDA_EXE" run -p "$ENV_PATH" bash -c "export KMP_DUPLICATE_LIB_OK=TRUE && export OMP_NUM_THREADS=1 && pymol -d 'import sys, os; sys.path.insert(0, os.path.expanduser('\"'\"'~/.pymol/startup'\"'\"')); import gluetk; gluetk.gluetk_gui()'"
else
  osascript -e 'display alert "Error" message "Conda env gluetk not found or invalid. Please re-run GlueTK Installer to create it."'
  exit 1
fi
''')

        os.chmod(launcher_script, 0o755)
        self._log(f"  ✅ Created {desktop_app}")


def main():
    root = tk.Tk()
    app = InstallerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
PYTHON_INSTALLER_EOF

echo -e "${GREEN}   ✅ Installer.app rebuilt${NC}"

# ============================================================================
# STEP 3: Package DMG (from package_dmg.sh)
# ============================================================================
echo -e "${BLUE}[3/6] Packaging GlueTK into DMG...${NC}"

DMG_NAME="GlueTK_Installer_${NEW_VERSION}"

# 3.1 Clean and create build directory
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"

# 3.2 Copy Installer.app
cp -R "${APP_NAME}" "${BUILD_DIR}/${APP_NAME}"

# 3.3 Bundle gluetk source code
mkdir -p "${BUILD_DIR}/${APP_NAME}/Contents/Resources/gluetk"
rsync -av --exclude='__pycache__' --exclude='*.pyc' --exclude='.DS_Store' \
      "${GLUETK_LIB}/" "${BUILD_DIR}/${APP_NAME}/Contents/Resources/gluetk/"

# 3.4 Create DMG
rm -f "${DMG_NAME}.dmg"
hdiutil create -volname "GlueTK Installer" \
    -srcfolder "${BUILD_DIR}" \
    -ov -format UDZO \
    "${DMG_NAME}.dmg"

echo -e "${GREEN}   ✅ DMG created: ${DMG_NAME}.dmg${NC}"

# ============================================================================
# STEP 4: Update PyMOL startup directory
# ============================================================================
echo -e "${BLUE}[4/6] Updating PyMOL startup directory...${NC}"
if [ -d ~/.pymol/startup/gluetk ]; then
    rsync -av --exclude='__pycache__' --exclude='*.pyc' --exclude='.DS_Store' \
          gluetk/ ~/.pymol/startup/gluetk/
    echo -e "${GREEN}   ✅ Updated ~/.pymol/startup/gluetk/${NC}"
else
    echo -e "${YELLOW}   ⚠️  PyMOL startup directory not found (will be created on first install)${NC}"
fi

# ============================================================================
# STEP 5: Clean up
# ============================================================================
echo -e "${BLUE}[5/6] Cleaning up...${NC}"
rm -rf build_dmg/GlueTK\ Installer.app
echo -e "${GREEN}   ✅ Cleaned up build artifacts${NC}"

# ============================================================================
# STEP 6: Summary
# ============================================================================
echo ""
echo -e "${BLUE}═══════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}✅ Release $NEW_VERSION complete!${NC}"
echo ""
echo -e "${BLUE}📦 DMG: ${DMG_NAME}.dmg${NC}"
echo ""
echo -e "${BLUE}Next steps:${NC}"
echo "  1. Test the DMG on a clean system"
echo "  2. Upload to GitHub releases"
echo "  3. Update README if needed"
echo -e "${BLUE}═══════════════════════════════════════════════════════════${NC}"

