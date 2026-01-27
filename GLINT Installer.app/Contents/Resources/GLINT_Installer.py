#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GLINT macOS Installer - 图形化安装程序
Professional GUI installer for macOS
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
    print("Error: tkinter not found. Please install it (e.g., sudo apt-get install python3-tk).")
    sys.exit(1)

# 配置
ENV_NAME = "glint"
PYTHON_VERSION = "3.10"
DEFAULT_INSTALL_PATH = os.path.join(os.path.expanduser("~"), ".pymol", "startup", "glint")

# Conda 依赖包 (与 macOS 版本保持一致)
# 注意: haddock_biobb 已移除，因为 haddocking channel 不可用 (HTTP 404)
CONDA_PACKAGES = [
    "rdkit", "scipy", "matplotlib", "pillow", "numpy=1.26.4", # 指定 numpy 版本
    "pandas", "seaborn", "pyqt", "openbabel", "pymol-open-source",
    "meeko", "vina", "scikit-image",
    "pdb2pqr",
]

# Pip 包 (open3d 在 conda 上不稳定, haddock3 因 haddocking channel 不可用改用 pip)
PIP_PACKAGES = ["requests", "open3d", "haddock3"]


def get_glint_source_dir():
    """获取 GLINT 源码目录"""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. macOS .app bundle Resources 目录 (优先)
    resources_dir = os.path.join(os.path.dirname(script_dir), "Resources", "glint")
    if os.path.isdir(resources_dir):
        return resources_dir

    # 2. 同级目录的 glint
    bundled_dir = os.path.join(script_dir, "glint")
    if os.path.isdir(bundled_dir):
        return bundled_dir

    # 3. PyInstaller 打包后的临时目录
    if hasattr(sys, '_MEIPASS'):
        bundled_dir = os.path.join(sys._MEIPASS, "glint")
        if os.path.isdir(bundled_dir):
            return bundled_dir

    return None


class InstallerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GLINT Installer")
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
        
        self.root.after(500, lambda: threading.Thread(target=self._check_environment, daemon=True).start())
    
    def _configure_styles(self):
        """Configure better fonts and styles for macOS"""
        style = ttk.Style()

        # Use aqua theme on macOS
        available_themes = style.theme_names()
        if 'aqua' in available_themes:
            style.theme_use('aqua')

        # macOS system fonts
        self.title_font = ("SF Pro Display", 22, "bold")
        self.subtitle_font = ("SF Pro Text", 11)
        self.normal_font = ("SF Pro Text", 10)
        self.small_font = ("SF Pro Text", 9)
        self.mono_font = ("Menlo", 9)

        # Fallback to Helvetica if SF Pro is not available
        try:
            import tkinter.font as tkfont
            available_fonts = tkfont.families()

            if "SF Pro Display" not in available_fonts:
                self.title_font = ("Helvetica Neue", 22, "bold")
                self.subtitle_font = ("Helvetica Neue", 11)
                self.normal_font = ("Helvetica Neue", 10)
                self.small_font = ("Helvetica Neue", 9)

            if "Menlo" not in available_fonts:
                self.mono_font = ("Monaco", 9)
        except:
            pass

        # Configure ttk styles with better fonts
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
            script_dir = os.path.dirname(os.path.abspath(__file__))

            # Try macOS .app bundle Resources directory first
            resources_dir = os.path.join(os.path.dirname(script_dir), "Resources")
            icon_path = os.path.join(resources_dir, "AppIcon.png")

            if not os.path.exists(icon_path):
                icon_path = os.path.join(script_dir, "glint", "assets", "logo.png")

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
        # Main frame
        main = ttk.Frame(self.root, padding=20)
        main.pack(fill=tk.BOTH, expand=True)
        
        # Title section
        title_frame = ttk.Frame(main)
        title_frame.pack(fill=tk.X, pady=(0, 15))
        
        title = ttk.Label(title_frame, text="🧬 GLINT Installer",
                         font=self.title_font)
        title.pack()
        
        subtitle = ttk.Label(title_frame,
                            text="PyMOL Plugin for Molecular Glue Analysis",
                            font=self.subtitle_font, foreground="#666666")
        subtitle.pack(pady=(5, 0))
        
        # 从 _version.py 动态获取版本
        try:
            from glint._version import __version__
        except ImportError:
            __version__ = "Unknown"
        version = ttk.Label(title_frame, text=f"Version: v{__version__}",
                           font=self.small_font, foreground="#888888")
        version.pack(pady=(3, 0))
        
        # Environment Status section
        status_frame = ttk.LabelFrame(main, text=" Environment Status ", padding=12)
        status_frame.pack(fill=tk.X, pady=(0, 15))
        
        # Conda status row
        conda_row = ttk.Frame(status_frame)
        conda_row.pack(fill=tk.X, pady=5)
        ttk.Label(conda_row, text="Conda:", font=self.normal_font, width=22).pack(side=tk.LEFT)
        self.conda_status = ttk.Label(conda_row, text="⏳ Checking...",
                                      font=self.normal_font, foreground="#E67E22")
        self.conda_status.pack(side=tk.LEFT, padx=10)
        self.conda_install_btn = ttk.Button(conda_row, text="Download Miniconda",
                                            command=self._install_conda, state=tk.DISABLED)
        self.conda_install_btn.pack(side=tk.RIGHT)
        
        # Environment status row
        env_row = ttk.Frame(status_frame)
        env_row.pack(fill=tk.X, pady=5)
        ttk.Label(env_row, text=f"Conda Env '{ENV_NAME}':",
                  font=self.normal_font, width=22).pack(side=tk.LEFT)
        self.env_status = ttk.Label(env_row, text="⏳ Checking...",
                                    font=self.normal_font, foreground="#E67E22")
        self.env_status.pack(side=tk.LEFT, padx=10)
        
        # Installation Path section
        path_frame = ttk.LabelFrame(main, text=" Installation Path ", padding=12)
        path_frame.pack(fill=tk.X, pady=(0, 15))
        
        path_row = ttk.Frame(path_frame)
        path_row.pack(fill=tk.X)
        self.path_entry = ttk.Entry(path_row, textvariable=self.install_path,
                                    width=60, font=self.normal_font)
        self.path_entry.pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(path_row, text="Browse...", command=self._browse_path).pack(side=tk.RIGHT, padx=(10, 0))
        
        path_note = ttk.Label(path_frame,
                             text="📌 GLINT will be installed here. PyMOL loads plugins from ~/.pymol/startup/",
                             font=self.small_font, foreground="#888888")
        path_note.pack(anchor=tk.W, pady=(8, 0))
        
        # Options section
        opts_frame = ttk.LabelFrame(main, text=" Options ", padding=12)
        opts_frame.pack(fill=tk.X, pady=(0, 15))
        
        ttk.Checkbutton(opts_frame, text="Install/Update dependencies (conda packages + PyMOL)",
                       variable=self.install_deps).pack(anchor=tk.W, pady=3)
        ttk.Checkbutton(opts_frame, text="Create desktop shortcut",
                       variable=self.create_shortcut).pack(anchor=tk.W, pady=3)
        
        # Progress section
        progress_frame = ttk.LabelFrame(main, text=" Progress ", padding=12)
        progress_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 15))
        
        self.progress = ttk.Progressbar(progress_frame, mode="determinate", length=500)
        self.progress.pack(fill=tk.X, pady=(0, 10))
        
        # Log text area with better styling
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
        
        # Button section
        btn_frame = ttk.Frame(main)
        btn_frame.pack(fill=tk.X, pady=(5, 0))
        
        self.install_btn = ttk.Button(btn_frame, text="🚀 Install GLINT",
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
        
        # 1. 查找 Conda
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
        
        # 2. 检查环境
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
        # 检查 PATH
        try:
            result = subprocess.run(["which", "conda"], capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                conda_path = result.stdout.strip()
                if os.path.exists(conda_path):
                    return conda_path
        except:
            pass

        # 常见路径 (macOS)
        home = os.path.expanduser("~")
        candidates = [
            os.path.join(home, "miniconda3", "bin", "conda"),
            os.path.join(home, "anaconda3", "bin", "conda"),
            "/opt/miniconda3/bin/conda",
            "/opt/anaconda3/bin/conda",
            "/opt/homebrew/Caskroom/miniconda/base/bin/conda",
            "/opt/homebrew/Caskroom/miniconda/base/condabin/conda",
            "/usr/local/Caskroom/miniconda/base/bin/conda",
            "/usr/local/miniconda3/bin/conda",
            "/usr/local/anaconda3/bin/conda",
        ]

        for path in candidates:
            if os.path.exists(path):
                return path

        return None
    
    def _install_conda(self):
        """打开 Miniconda 下载页面"""
        webbrowser.open("https://docs.conda.io/en/latest/miniconda.html")
        messagebox.showinfo("Install Miniconda", 
                           "Please download and install Miniconda for Linux.\n\n"
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
            
            # 1. 创建/检查 Conda 环境
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
            
            # 2. 安装依赖
            if self.install_deps.get():
                self._log("\n" + "=" * 50)
                self._log("[2/5] Installing dependencies...")
                self._log("=" * 50)
                self._log("  This may take 10-20 minutes, please wait...")
                
                pkg_str = " ".join(CONDA_PACKAGES)
                self._log(f"  Installing: {pkg_str}")
                
                # 注意: haddocking channel 已不可用 (HTTP 404)，移除该 channel
                cmd = [
                    self.conda_exe, "install",
                    "-p", self.env_path,
                    "-c", "conda-forge",
                    "-c", "schrodinger", # Add schrodinger channel for pymol-open-source
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
                
                # 安装 pip 包
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
            
            # 3. 复制插件文件
            self._log("\n" + "=" * 50)
            self._log("[3/5] Installing plugin files...")
            self._log("=" * 50)
            
            install_path = self.install_path.get()
            source_dir = get_glint_source_dir()
            
            if source_dir and os.path.isdir(source_dir):
                # 创建目标目录
                os.makedirs(install_path, exist_ok=True)
                
                # 删除旧安装
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
                
                # 复制文件
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
            
            # 4. 创建快捷方式
            if self.create_shortcut.get():
                self._log("\n" + "=" * 50)
                self._log("[4/5] Creating desktop shortcut...")
                self._log("=" * 50)
                self._create_shortcut()
            else:
                self._log("\n[4/5] Skipping shortcut creation")
            
            self.progress["value"] = 100
            
            # 5. 完成
            self._log("\n" + "=" * 50)
            self._log("[5/5] ✅ Installation complete!")
            self._log("=" * 50)
            self._log("\nHow to use GLINT:")
            self._log("  1. Open GLINT from ~/Applications/GLINT.app")
            self._log("  2. Or run in PyMOL: glint_gui")

            self.root.after(0, lambda: messagebox.showinfo(
                "Success",
                "GLINT installed successfully!\n\n"
                "You can now launch GLINT from ~/Applications/GLINT.app"))
            
        except Exception as e:
            self._log(f"\n❌ Error: {e}")
            import traceback
            self._log(traceback.format_exc())
            self.root.after(0, lambda: messagebox.showerror("Error", str(e)))
        finally:
            self.root.after(0, lambda: self.install_btn.configure(state=tk.NORMAL))
    
    def _create_shortcut(self):
        """Create macOS application shortcut"""
        home = os.path.expanduser("~")
        app_dir = os.path.join(home, "Applications")
        os.makedirs(app_dir, exist_ok=True)

        app_name = "GLINT.app"
        app_path = os.path.join(app_dir, app_name)

        # Create .app bundle structure
        contents_dir = os.path.join(app_path, "Contents")
        macos_dir = os.path.join(contents_dir, "MacOS")
        resources_dir = os.path.join(contents_dir, "Resources")

        os.makedirs(macos_dir, exist_ok=True)
        os.makedirs(resources_dir, exist_ok=True)

        # Create launcher script
        launcher_script = os.path.join(macos_dir, "GLINT")
        with open(launcher_script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f'eval "$({self.conda_exe} shell.bash hook)"\n')
            f.write(f'conda activate {ENV_NAME}\n')
            f.write('export KMP_DUPLICATE_LIB_OK=TRUE\n')
            f.write('export OMP_NUM_THREADS=1\n')
            f.write('pymol -d "import sys, os; sys.path.insert(0, os.path.expanduser(\'~/.pymol/startup\')); import glint; glint.glint_gui()"\n')
        os.chmod(launcher_script, 0o755)

        # Create Info.plist
        plist_path = os.path.join(contents_dir, "Info.plist")
        with open(plist_path, "w") as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write('<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">\n')
            f.write('<plist version="1.0">\n')
            f.write('<dict>\n')
            f.write('    <key>CFBundleExecutable</key>\n')
            f.write('    <string>GLINT</string>\n')
            f.write('    <key>CFBundleIconFile</key>\n')
            f.write('    <string>AppIcon</string>\n')
            f.write('    <key>CFBundleIdentifier</key>\n')
            f.write('    <string>com.glint.app</string>\n')
            f.write('    <key>CFBundleName</key>\n')
            f.write('    <string>GLINT</string>\n')
            f.write('    <key>CFBundlePackageType</key>\n')
            f.write('    <string>APPL</string>\n')
            f.write('    <key>CFBundleShortVersionString</key>\n')
            f.write('    <string>1.0</string>\n')
            f.write('</dict>\n')
            f.write('</plist>\n')

        # Copy icon if available
        icon_src = os.path.join(self.install_path.get(), "assets", "logo.png")
        if os.path.exists(icon_src):
            icon_dst = os.path.join(resources_dir, "AppIcon.png")
            shutil.copy2(icon_src, icon_dst)

        self._log(f"  ✅ Created application: {app_path}")
        self._log(f"  You can find GLINT in ~/Applications/")


def main():
    root = tk.Tk()
    app = InstallerApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()