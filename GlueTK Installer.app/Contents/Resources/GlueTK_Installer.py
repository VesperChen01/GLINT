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

# 包含 Vina 和 Meeko 的完整依赖列表
CONDA_PACKAGES = [
    "rdkit", "scipy", "matplotlib", "pillow", "numpy", 
    "pandas", "seaborn", "pyqt", "openbabel", "pymol-open-source",
    "autodock-vina", "meeko"
]

PIP_PACKAGES = ["requests"]

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
            
        # 2. Check Environment
        try:
            result = subprocess.run([self.conda_exe, "env", "list"], capture_output=True, text=True, timeout=10)
            if ENV_NAME in result.stdout:
                self.env_status.configure(text="✅ Exists", foreground="green")
                self.env_ok = True
                self._log(f"  Environment '{ENV_NAME}': Exists")
            else:
                self.env_status.configure(text="⚠️ Will create", foreground="orange")
                self._log(f"  Environment '{ENV_NAME}': Will be created")
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
            if not self.env_ok:
                self._log(f"  Creating environment '{ENV_NAME}'...")
                # Use full conda path
                result = subprocess.run([self.conda_exe, "create", "-n", ENV_NAME, f"python={PYTHON_VERSION}", "-y"], capture_output=True, text=True)
                if result.returncode != 0:
                    self._log(f"  Error: {result.stderr}")
                else:
                    self._log("  ✅ Environment created")
            else:
                self._log("  Environment already exists")
                
            # 2. Dependencies
            if self.install_deps.get():
                self._log("\\n[2/5] Installing dependencies (conda)...")
                self.progress["value"] = 40
                pkg_str = " ".join(CONDA_PACKAGES)
                # Use full conda path and list args instead of shell=True for better safety/stability
                self._log(f"  Running: conda install -n {ENV_NAME} -c conda-forge {pkg_str}")
                cmd = [self.conda_exe, "install", "-n", ENV_NAME, "-c", "conda-forge", "-y"] + CONDA_PACKAGES
                
                # Run without shell=True
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    self._log("  ✅ Conda packages installed")
                else:
                    self._log(f"  ⚠️ Some packages may have failed")
                
                # 验证 PyMOL 是否安装成功
                pymol_check = subprocess.run([self.conda_exe, "run", "-n", ENV_NAME, "which", "pymol"], 
                                            capture_output=True, text=True)
                if pymol_check.returncode != 0:
                    self._log("  ⚠️ PyMOL not installed in conda, will check system PyMOL.app")
                
                # 2b. Pip dependencies inside the environment
                if PIP_PACKAGES:
                    self._log("\\n[2b/5] Installing pip packages inside environment...")
                    pip_cmd = [self.conda_exe, "run", "-n", ENV_NAME, "python", "-m", "pip", "install", "--upgrade"] + PIP_PACKAGES
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
                
                # 检测 PyMOL 安装方式
                pymol_app_path = "/Applications/PyMOL.app/Contents/MacOS/PyMOL"
                has_pymol_app = os.path.exists(pymol_app_path) and os.access(pymol_app_path, os.X_OK)
                
                # 检查 conda 环境中是否有 pymol
                has_conda_pymol = False
                try:
                    check_result = subprocess.run([self.conda_exe, "run", "-n", ENV_NAME, "which", "pymol"],
                                                capture_output=True, text=True, timeout=5)
                    has_conda_pymol = (check_result.returncode == 0)
                except:
                    pass
                
                with open(launcher_script, "w") as f:
                    if has_pymol_app:
                        # 使用系统 PyMOL.app
                        self._log("  Using system PyMOL.app")
                        f.write(f'''#!/bin/bash
# GlueTK Launcher (PyMOL.app + conda dependencies)

# 激活 conda 环境以加载依赖
export CONDA_EXE="{conda_exe_str}"
eval "$(\$CONDA_EXE shell.bash hook)" 2>/dev/null
conda activate {ENV_NAME} 2>/dev/null || true

# 使用系统 PyMOL
"{pymol_app_path}" -d "import sys, os; sys.path.insert(0, os.path.expanduser('~/.pymol/startup')); import gluetk; gluetk.gluetk_gui()"
''')
                    elif has_conda_pymol:
                        # 使用 conda PyMOL
                        self._log("  Using conda PyMOL")
                        f.write(f'''#!/bin/bash
# GlueTK Launcher (conda pymol)

CONDA_EXE="{conda_exe_str}"
if [ -n "$CONDA_EXE" ] && [ -x "$CONDA_EXE" ]; then
  echo "Starting GlueTK via conda env {ENV_NAME}..."
  "$CONDA_EXE" run -n {ENV_NAME} pymol -d "import sys, os; sys.path.insert(0, os.path.expanduser('~/.pymol/startup')); import gluetk; gluetk.gluetk_gui()"
else
  osascript -e 'display alert "Error" message "Conda not found. Please reinstall GlueTK."'
  exit 1
fi
''')
                    else:
                        # 没有找到 PyMOL，创建错误提示启动器
                        self._log("  ⚠️ PyMOL not found! Creating error handler launcher")
                        f.write(f'''#!/bin/bash
# GlueTK Launcher (PyMOL not found)

osascript -e 'display alert "PyMOL Not Found" message "Please install PyMOL.app from https://pymol.org/ or run the installer again to install conda PyMOL." buttons {{"OK"}}'
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
