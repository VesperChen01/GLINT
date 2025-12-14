#!/bin/bash
# rebuild_installer.sh
# 重新创建 GlueTK Installer.app

set -e

APP_NAME="GlueTK Installer.app"
CONTENTS="${APP_NAME}/Contents"
MACOS="${CONTENTS}/MacOS"
RESOURCES="${CONTENTS}/Resources"

echo "🏗️  Rebuilding ${APP_NAME}..."

# 1. 创建目录
rm -rf "${APP_NAME}"
mkdir -p "${MACOS}"
mkdir -p "${RESOURCES}"

# 2. 复制图标 (如果存在)
if [ -f "gluetk/assets/AppIcon.icns" ]; then
    cp "gluetk/assets/AppIcon.icns" "${RESOURCES}/AppIcon.icns"
    # 同时生成一个 png 用于资源
    cp "gluetk/assets/logo.png" "${RESOURCES}/AppIcon.png" 2>/dev/null || true
else
    echo "⚠️  Icon not found in gluetk/assets/"
fi

# 3. 创建 Info.plist
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
    <string>0.1.4-beta</string>
</dict>
</plist>
EOF

# 4. 创建启动脚本 (launcher)
# 使用 Conda 的 Python 启动 Tk GUI 安装器（与老版本界面一致）
cat > "${MACOS}/launcher" << 'EOF'
#!/bin/bash
# GlueTK Installer Launcher (Tk GUI via conda python)

# 查找 conda
# 查找 conda
CONDA_EXE=""
if command -v conda &> /dev/null; then
    CONDA_EXE=$(command -v conda)
else
    # 尝试常见路径
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

# 初始化 conda 环境
eval "$($CONDA_EXE shell.bash hook)"

CONDA_BASE="$(conda info --base 2>/dev/null)"
if [ ! -d "$CONDA_BASE" ]; then
  osascript -e 'display alert "Conda Base Not Found" message "Conda is installed but CONDA_BASE is invalid. Please check your conda installation."'
  exit 1
fi

PYTHON="${CONDA_BASE}/bin/python"
if [ ! -x "$PYTHON" ]; then
  # 退而求其次，使用当前 shell 中的 python
  if command -v python &> /dev/null; then
    PYTHON=$(command -v python)
  else
    osascript -e 'display alert "Python Not Found" message "Cannot find a suitable python interpreter."'
    exit 1
  fi
fi

# 获取资源目录
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
RESOURCES_DIR="$(dirname "$DIR")/Resources"
INSTALLER_SCRIPT="${RESOURCES_DIR}/GlueTK_Installer.py"

"$PYTHON" "$INSTALLER_SCRIPT"
EOF

chmod +x "${MACOS}/launcher"

# 5. （可选）保留旧的 Tk 安装器脚本（目前不再使用，仅作为参考）
#    主要安装逻辑改为 Run_Install.command + gluetk/install.sh
cat > "${RESOURCES}/GlueTK_Installer.py" << 'EOF'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GlueTK Installer - 图形化安装程序
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
PYTHON_VERSION = "3.9"
DEFAULT_INSTALL_PATH = os.path.expanduser("~/.pymol/startup/gluetk")

# 包含 Vina / Meeko / PyMOL / HADDOCK3 的完整依赖列表
# 说明：
# - HADDOCK3 官方 PyPI（pip install haddock3）在 macOS 上可能触发本地编译并失败（例如 -march=native）。
# - 因此这里优先走 conda。
# - HADDOCK3 的 conda 构建目前在 bioconda 提供（包名：haddock_biobb），可提供 haddock3 CLI/模块。
# - autodock-vina 在 conda 频道中不稳定，改用 pip 安装 vina-split
CONDA_PACKAGES = [
    "rdkit", "scipy", "matplotlib", "pillow", "numpy",
    "pandas", "seaborn", "pyqt", "openbabel", "pymol-open-source",
    "meeko",
    "haddock_biobb"
]

# Pip 包：requests + vina-split（AutoDock Vina 的 pip 封装）
PIP_PACKAGES = ["requests", "vina-split"]

def get_gluetk_source_dir():
    """获取 GlueTK 源码目录"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 1. 优先查找 Bundle 内部 (Resources/gluetk) - 用于打包后的 App
    bundled_dir = os.path.join(script_dir, "gluetk")
    if os.path.isdir(bundled_dir):
        return bundled_dir
        
    # 2. 开发模式：从 app bundle 向上查找 repo 目录
    # Resources -> Contents -> app -> repo
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
        self.root.geometry("700x650")
        self.root.resizable(True, True)
        
        self.install_path = tk.StringVar(value=DEFAULT_INSTALL_PATH)
        self.create_launcher = tk.BooleanVar(value=True)
        self.install_deps = tk.BooleanVar(value=True)
        self.conda_ok = False
        self.env_ok = False
        self.conda_exe = None  # Store detected conda path
        self.env_path = None   # Store full path to conda env (force standard conda base)
        
        self._build_ui()
        # Fix for macOS Dark Mode / Blank Screen
        self.root.update()
        self.root.lift()
        self.root.attributes('-topmost',True)
        self.root.after_idle(self.root.attributes,'-topmost',False)
        
        self._center_window()
        
        # 只有在 UI 完全加载后才启动检查
        self.root.after(1000, lambda: threading.Thread(target=self._check_environment, daemon=True).start())
    
    def _center_window(self):
        self.root.update_idletasks()
        w = self.root.winfo_width()
        h = self.root.winfo_height()
        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        x = (sw - w) // 2
        y = (sh - h) // 2
        self.root.geometry(f"{w}x{h}+{x}+{y}")
    
    def _build_ui(self):
        main = ttk.Frame(self.root, padding=20)
        main.pack(fill=tk.BOTH, expand=True)
        
        title = ttk.Label(main, text="🧬 GlueTK Installer", font=("Helvetica", 24, "bold"))
        title.pack(pady=(0, 5))
        subtitle = ttk.Label(main, text="PyMOL Plugin for Molecular Glue Analysis", font=("Helvetica", 12), foreground="gray")
        subtitle.pack(pady=(0, 20))
        
        # Status
        status_frame = ttk.LabelFrame(main, text="Environment Status", padding=10)
        status_frame.pack(fill=tk.X, pady=(0, 15))
        
        conda_row = ttk.Frame(status_frame)
        conda_row.pack(fill=tk.X, pady=2)
        ttk.Label(conda_row, text="Conda:").pack(side=tk.LEFT)
        self.conda_status = ttk.Label(conda_row, text="Checking...", foreground="orange")
        self.conda_status.pack(side=tk.LEFT, padx=10)
        self.conda_install_btn = ttk.Button(conda_row, text="Install Miniconda", command=self._install_conda, state=tk.DISABLED)
        self.conda_install_btn.pack(side=tk.RIGHT)
        
        env_row = ttk.Frame(status_frame)
        env_row.pack(fill=tk.X, pady=2)
        ttk.Label(env_row, text=f"Conda Env '{ENV_NAME}':").pack(side=tk.LEFT)
        self.env_status = ttk.Label(env_row, text="Checking...", foreground="orange")
        self.env_status.pack(side=tk.LEFT, padx=10)
        
        # Path
        path_frame = ttk.LabelFrame(main, text="Installation Path", padding=10)
        path_frame.pack(fill=tk.X, pady=(0, 15))
        path_row = ttk.Frame(path_frame)
        path_row.pack(fill=tk.X)
        self.path_entry = ttk.Entry(path_row, textvariable=self.install_path, width=50)
        self.path_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(path_row, text="Browse...", command=self._browse_path).pack(side=tk.RIGHT, padx=(10, 0))
        path_note = ttk.Label(path_frame, text="📌 GlueTK will be installed here. PyMOL loads plugins from ~/.pymol/startup/", font=("Helvetica", 10), foreground="gray")
        path_note.pack(anchor=tk.W, pady=(5, 0))
        
        # Options
        opts_frame = ttk.LabelFrame(main, text="Options", padding=10)
        opts_frame.pack(fill=tk.X, pady=(0, 15))
        ttk.Checkbutton(opts_frame, text="Install/Update dependencies (conda packages)", variable=self.install_deps).pack(anchor=tk.W)
        ttk.Checkbutton(opts_frame, text="Create desktop app (GlueTK.app)", variable=self.create_launcher).pack(anchor=tk.W)
        
        # Progress
        progress_frame = ttk.LabelFrame(main, text="Progress", padding=10)
        progress_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))
        self.progress = ttk.Progressbar(progress_frame, mode="determinate", length=400)
        self.progress.pack(fill=tk.X, pady=(0, 10))
        log_frame = ttk.Frame(progress_frame)
        log_frame.pack(fill=tk.BOTH, expand=True)
        self.log_text = tk.Text(log_frame, height=8, state=tk.DISABLED, font=("Courier", 10), bg="#1e1e1e", fg="#d4d4d4")
        scrollbar = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Buttons
        btn_frame = ttk.Frame(main)
        btn_frame.pack(fill=tk.X)
        self.install_btn = ttk.Button(btn_frame, text="🚀 Install GlueTK", command=self._start_install)
        self.install_btn.pack(side=tk.RIGHT, padx=(10, 0))
        ttk.Button(btn_frame, text="Cancel", command=self.root.quit).pack(side=tk.RIGHT)
        
    def _log(self, msg):
        self.log_text.configure(state=tk.NORMAL)
        self.log_text.insert(tk.END, msg + "\n")
        self.log_text.see(tk.END)
        self.log_text.configure(state=tk.DISABLED)
        self.root.update()
    
    def _browse_path(self):
        path = filedialog.askdirectory(title="Select Installation Directory", initialdir=os.path.dirname(self.install_path.get()))
        if path:
            self.install_path.set(path)
            
    def _check_environment(self):
        self._log("Checking environment...")
        
        # 1. Find Conda Executable
        self.conda_exe = shutil.which("conda")
        if not self.conda_exe:
            # Try standard paths
            possible_paths = [
                os.path.expanduser("~/miniconda3/bin/conda"),
                os.path.expanduser("~/opt/miniconda3/bin/conda"),
                os.path.expanduser("~/anaconda3/bin/conda"),
                "/usr/local/bin/conda",
                "/opt/homebrew/bin/conda",
                "/opt/homebrew/Caskroom/miniconda/base/bin/conda"
            ]
            for p in possible_paths:
                if os.path.exists(p) and os.access(p, os.X_OK):
                    self.conda_exe = p
                    break
        
        if self.conda_exe:
            try:
                result = subprocess.run([self.conda_exe, "--version"], capture_output=True, text=True, timeout=5)
                if result.returncode == 0:
                    version = result.stdout.strip()
                    self.conda_status.configure(text=f"✅ {version}", foreground="green")
                    self.conda_ok = True
                    self.conda_install_btn.configure(state=tk.DISABLED)
                    self._log(f"  Conda found: {self.conda_exe}")
                    self._log(f"  Version: {version}")
                else:
                    raise Exception("conda command failed")
            except Exception as e:
                self.conda_status.configure(text="❌ Error", foreground="red")
                self._log(f"  Conda error: {e}")
        else:
            self.conda_status.configure(text="❌ Not found", foreground="red")
            self.conda_install_btn.configure(state=tk.NORMAL)
            self._log("  Conda: Not found in PATH or standard locations")
            return
            
        # 2. Check Environment（强制使用 conda base/envs 下的环境，避免误用 PyMOL.app 内的 gluetk）
        try:
            base_result = subprocess.run([self.conda_exe, "info", "--base"], capture_output=True, text=True, timeout=10)
            conda_base = base_result.stdout.strip() if base_result.returncode == 0 else ""
            if not conda_base:
                raise Exception("Cannot determine conda base")

            self.env_path = os.path.join(conda_base, "envs", ENV_NAME)

            if os.path.isdir(self.env_path):
                self.env_status.configure(text="✅ Exists", foreground="green")
                self.env_ok = True
                self._log(f"  Environment '{ENV_NAME}': {self.env_path}")
            else:
                self.env_status.configure(text="⚠️ Will create", foreground="orange")
                self._log(f"  Environment '{ENV_NAME}': Will be created at {self.env_path}")
        except Exception as e:
            self.env_status.configure(text="❓ Unknown", foreground="gray")
            self._log(f"  Environment check failed: {e}")
        self._log("Ready to install.")
        
    def _install_conda(self):
        webbrowser.open("https://docs.conda.io/en/latest/miniconda.html")
        messagebox.showinfo("Install Miniconda", "Please download and install Miniconda, then restart this installer.")
        
    def _start_install(self):
        self.install_btn.configure(state=tk.DISABLED)
        threading.Thread(target=self._do_install, daemon=True).start()
        
    def _do_install(self):
        try:
            self.progress["value"] = 0
            
            # 1. Conda Env
            self._log("\n[1/5] Checking conda environment...")
            self.progress["value"] = 20
            
            # 确保 env_path 已设置（防止 None 错误）
            if not self.env_path:
                try:
                    base_result = subprocess.run([self.conda_exe, "info", "--base"], capture_output=True, text=True, timeout=10)
                    conda_base = base_result.stdout.strip() if base_result.returncode == 0 else ""
                    if conda_base:
                        self.env_path = os.path.join(conda_base, "envs", ENV_NAME)
                    else:
                        raise Exception("Cannot determine conda base directory")
                except Exception as e:
                    raise Exception(f"Failed to determine environment path: {e}")
            
            if not self.env_ok:
                self._log(f"  Creating environment at {self.env_path}...")
                # 强制使用 -p 指定路径，避免与 PyMOL.app 内同名环境冲突
                result = subprocess.run([self.conda_exe, "create", "-p", self.env_path, f"python={PYTHON_VERSION}", "-y"], capture_output=True, text=True)
                if result.returncode != 0:
                    self._log(f"  Error: {result.stderr}")
                    raise Exception(result.stderr.strip() or "conda create failed")
                else:
                    self._log("  ✅ Environment created")
            else:
                self._log(f"  Environment already exists: {self.env_path}")
                
            # 2. Dependencies
            if self.install_deps.get():
                self._log("\\n[2/5] Installing dependencies (conda)...")
                self.progress["value"] = 40
                pkg_str = " ".join(CONDA_PACKAGES)
                # 使用 conda-forge + bioconda，并启用 libmamba solver（更稳定）
                self._log(f"  Running: conda install -p {self.env_path} -c conda-forge -c bioconda {pkg_str}")
                cmd = [
                    self.conda_exe, "install",
                    "-p", self.env_path,
                    "-c", "conda-forge",
                    "-c", "bioconda",
                    "--solver=libmamba",
                    "-y",
                ] + CONDA_PACKAGES
                
                # Run without shell=True
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    self._log("  ✅ Conda packages installed")
                else:
                    # 输出更可诊断的信息，并中止（否则后续一定检测不到）
                    if result.stdout:
                        self._log("  ---- conda stdout ----")
                        self._log(result.stdout.strip()[-4000:])
                    if result.stderr:
                        self._log("  ---- conda stderr ----")
                        self._log(result.stderr.strip()[-4000:])
                    raise Exception("Conda dependency installation failed")
                
                # 验证 PyMOL 是否安装成功（我们强制使用 conda PyMOL）
                pymol_check = subprocess.run([self.conda_exe, "run", "-p", self.env_path, "which", "pymol"],
                                            capture_output=True, text=True)
                if pymol_check.returncode != 0:
                    raise Exception("Conda PyMOL (pymol-open-source) not found in gluetk env")
                
                # 2b. Pip dependencies inside the environment
                if PIP_PACKAGES:
                    self._log("\\n[2b/5] Installing pip packages inside environment...")
                    env_pip = os.path.join(self.env_path, "bin", "pip") if self.env_path else None
                    if env_pip and os.path.exists(env_pip):
                        pip_cmd = [env_pip, "install", "--upgrade", "--disable-pip-version-check"] + PIP_PACKAGES
                    else:
                        pip_cmd = [self.conda_exe, "run", "-p", self.env_path, "python", "-m", "pip", "install", "--upgrade", "--disable-pip-version-check"] + PIP_PACKAGES
                    self._log(f"  Running: {' '.join(pip_cmd)}")
                    pip_result = subprocess.run(pip_cmd, capture_output=True, text=True)
                    if pip_result.returncode == 0:
                        self._log("  ✅ Pip packages installed")
                    else:
                        self._log(f"  ⚠️ Pip install issues: {pip_result.stderr}")
            else:
                self._log("\\n[2/5] Skipping dependencies (unchecked)")
                
            self.progress["value"] = 60
            
            # 3. Copy Plugin Files
            self._log("\n[3/5] Installing plugin files...")
            install_path = self.install_path.get()
            source_dir = get_gluetk_source_dir()
            
            if source_dir and os.path.isdir(source_dir):
                if os.path.exists(install_path) and os.path.samefile(source_dir, install_path):
                    self._log(f"  ⚠️ Destination same as source. Skipping copy.")
                else:
                    os.makedirs(install_path, exist_ok=True)
                    IGNORED_FILES = {'.DS_Store', '__pycache__', '.git', '.gitignore'}
                    for item in os.listdir(source_dir):
                        if item in IGNORED_FILES: continue
                        src = os.path.join(source_dir, item)
                        dst = os.path.join(install_path, item)
                        try:
                            if os.path.isdir(src):
                                if os.path.exists(dst):
                                    if os.path.samefile(src, dst): continue
                                    shutil.rmtree(dst)
                                shutil.copytree(src, dst, ignore=shutil.ignore_patterns(*IGNORED_FILES))
                            else:
                                if os.path.exists(dst) and os.path.samefile(src, dst): continue
                                shutil.copy2(src, dst)
                        except Exception as e:
                            self._log(f"  ⚠️ Failed to copy {item}: {e}")
                    self._log(f"  ✅ Copied to {install_path}")
            else:
                self._log(f"  ⚠️ Source directory not found: {source_dir}")
                
            self.progress["value"] = 80
            
            # 4. Create Desktop App
            if self.create_launcher.get():
                self._log("\n[4/5] Creating desktop App...")
                desktop_app = os.path.expanduser("~/Desktop/GlueTK.app")
                if os.path.exists(desktop_app):
                    shutil.rmtree(desktop_app)
                
                contents = os.path.join(desktop_app, "Contents")
                macos = os.path.join(contents, "MacOS")
                resources = os.path.join(contents, "Resources")
                os.makedirs(macos, exist_ok=True)
                os.makedirs(resources, exist_ok=True)
                
                # Icon
                icon_src = None
                if source_dir:
                    c1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "AppIcon.icns")
                    c2 = os.path.join(source_dir, "assets", "AppIcon.icns")
                    if os.path.exists(c1): icon_src = c1
                    elif os.path.exists(c2): icon_src = c2
                
                if icon_src:
                    shutil.copy2(icon_src, os.path.join(resources, "AppIcon.icns"))
                    
                # Info.plist
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
</dict>
</plist>''')

                # Launcher（智能检测 PyMOL 类型）
                conda_exe_str = self.conda_exe or ""
                launcher_script = os.path.join(macos, "launcher")
                
                # 强制统一：永远使用 conda 环境里的 PyMOL（pymol-open-source）
                has_conda_pymol = False
                try:
                    check_result = subprocess.run([self.conda_exe, "run", "-p", self.env_path, "which", "pymol"],
                                                capture_output=True, text=True, timeout=5)
                    has_conda_pymol = (check_result.returncode == 0)
                except:
                    pass

                with open(launcher_script, "w") as f:
                    if has_conda_pymol:
                        self._log("  Using conda PyMOL (forced)")
                        env_path_str = self.env_path or ""
                        f.write(f'''#!/bin/bash
# GlueTK Launcher (FORCED conda pymol-open-source)

CONDA_EXE="{conda_exe_str}"
ENV_PATH="{env_path_str}"
if [ -n "$CONDA_EXE" ] && [ -x "$CONDA_EXE" ] && [ -n "$ENV_PATH" ] && [ -d "$ENV_PATH" ]; then
  echo "Starting GlueTK via conda env: $ENV_PATH"
  "$CONDA_EXE" run -p "$ENV_PATH" pymol -d "import sys, os; sys.path.insert(0, os.path.expanduser('~/.pymol/startup')); import gluetk; gluetk.gluetk_gui()"
else
  osascript -e 'display alert "Error" message "Conda env gluetk not found. Please re-run GlueTK Installer to create it."'
  exit 1
fi
''')
                    else:
                        self._log("  ⚠️ Conda PyMOL not found (pymol-open-source install may have failed)")
                        f.write(f'''#!/bin/bash
# GlueTK Launcher (conda pymol missing)

osascript -e 'display alert "PyMOL Missing" message "Conda PyMOL (pymol-open-source) was not found in the gluetk environment. Please re-run the installer with dependency installation enabled." buttons {"OK"}'
exit 1
''')
                
                os.chmod(launcher_script, 0o755)
                self._log(f"  ✅ Created {desktop_app}")
            else:
                self._log("\n[4/5] Skipping app creation")
                
            self.progress["value"] = 100
            self._log("\n[5/5] ✅ Installation complete!")
            self._log("\nStart GlueTK by double-clicking GlueTK.app on your Desktop.")
            
            self.root.after(0, lambda: messagebox.showinfo("Success", "GlueTK installed successfully!"))
            
        except Exception as e:
            self._log(f"\n❌ Error: {e}")
            self.root.after(0, lambda: messagebox.showerror("Error", str(e)))
        finally:
            self.root.after(0, lambda: self.install_btn.configure(state=tk.NORMAL))

def main():
    root = tk.Tk()
    app = InstallerApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
EOF

# 6. （保留 Terminal 版入口作为高级用户备用，不在 launcher 中默认调用）
cat > "${RESOURCES}/Run_Install.command" << 'EOF'
#!/bin/bash
# Optional: Terminal-based GlueTK installer entry point (fallback)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GLUETK_DIR="${SCRIPT_DIR}/gluetk"

if [ ! -d "$GLUETK_DIR" ]; then
  echo "[GlueTK Installer] ERROR: GlueTK sources not found in: $GLUETK_DIR"
  echo "If you are running from a source checkout, please run: bash gluetk/install.sh"
  read -n1 -p "Press any key to exit..." _
  exit 1
fi

cd "$GLUETK_DIR"

if [ -x "install.sh" ]; then
  bash install.sh
else
  echo "[GlueTK Installer] ERROR: install.sh not found or not executable."
  read -n1 -p "Press any key to exit..." _
  exit 1
fi

read -n1 -p "Press any key to close this window..." _
EOF

chmod +x "${RESOURCES}/Run_Install.command"

echo "✅ Rebuilt GlueTK Installer.app"
echo "   Now run: bash package_dmg.sh"
