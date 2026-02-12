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

HAS_TK = False
_tk_version = 0.0
try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    HAS_TK = True
    _tk_version = float(tk.TkVersion)
except ImportError:
    tk = None

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
    "metis",  # apbs-binary (pip) 依赖 libmetis.dylib
]

# Pip 包 (open3d 在 conda 上不稳定, haddock3 因 haddocking channel 不可用改用 pip)
# apbs-binary: 通过 pip 安装 macOS APBS 二进制，用于 EC 静电互补性分析
PIP_PACKAGES = ["requests", "open3d", "haddock3", "apbs-binary"]

# Conda Terms of Service channels (newer conda requires explicit acceptance in non-interactive mode)
CONDA_TOS_CHANNELS = [
    "https://repo.anaconda.com/pkgs/main",
    "https://repo.anaconda.com/pkgs/r",
]



def find_conda():
    """查找 Conda 安装路径 (macOS) — shared by GUI and CLI installers."""
    # 检查 PATH
    try:
        result = subprocess.run(["which", "conda"], capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            conda_path = result.stdout.strip()
            if os.path.exists(conda_path):
                return conda_path
    except Exception:
        pass

    # 常见路径 (macOS)
    home = os.path.expanduser("~")
    candidates = [
        os.path.join(home, "miniforge3", "bin", "conda"),
        os.path.join(home, "mambaforge", "bin", "conda"),
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

        # macOS: bring window to front
        if sys.platform == 'darwin':
            self.root.lift()
            self.root.focus_force()

        self.root.after(500, lambda: threading.Thread(target=self._check_environment, daemon=True).start())

    def _configure_styles(self):
        """Configure ttk styles."""
        self.title_font = ("Helvetica", 22, "bold")
        self.subtitle_font = ("Helvetica", 11)
        self.normal_font = ("Helvetica", 10)
        self.small_font = ("Helvetica", 9)
        self.mono_font = ("Menlo", 9)

        style = ttk.Style()
        style.configure("TLabel", font=self.normal_font)
        style.configure("TButton", font=self.normal_font, padding=6)
        style.configure("TCheckbutton", font=self.normal_font)
        style.configure("TEntry", font=self.normal_font)
        style.configure("TLabelframe", font=self.normal_font)
        style.configure("TLabelframe.Label", font=self.normal_font)

        try:
            import tkinter.font as tkfont
            available_fonts = tkfont.families()
            if "Menlo" not in available_fonts:
                self.mono_font = ("Monaco", 9)
        except Exception:
            pass

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
        version = ttk.Label(title_frame, text=f"Version: {__version__}",
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
                               font=self.mono_font, bg="#ffffff", fg="#222222",
                               insertbackground="#222222", selectbackground="#cfe8ff",
                               relief=tk.SOLID, bd=1, padx=10, pady=8)
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
        """写入日志 — thread-safe."""
        def _do_log():
            try:
                self.log_text.configure(state=tk.NORMAL)
                self.log_text.insert(tk.END, msg + "\n")
                self.log_text.see(tk.END)
                self.log_text.configure(state=tk.DISABLED)
            except Exception:
                pass

        if threading.current_thread() is threading.main_thread():
            _do_log()
        else:
            self.root.after(0, _do_log)

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
        return find_conda()

    def _install_conda(self):
        """打开 Miniconda 下载页面（macOS）"""
        webbrowser.open("https://docs.conda.io/en/latest/miniconda.html")
        messagebox.showinfo(
            "Install Miniconda",
            "Please download and install Miniconda for macOS.\n\n"
            "Recommended: choose Miniconda3 macOS installer matching your CPU (Apple Silicon or Intel).\n"
            "After installation, restart this installer."
        )


    def _ensure_conda_tos_accepted(self):
        """Ensure required conda ToS channels are accepted for non-interactive install."""
        self._log("  Checking conda Terms of Service...")
        for ch in CONDA_TOS_CHANNELS:
            cmd = [
                self.conda_exe, "tos", "accept",
                "--override-channels",
                "--channel", ch,
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                self._log(f"  ✅ ToS accepted: {ch}")
            else:
                stderr = (result.stderr or "").strip()
                stdout = (result.stdout or "").strip()
                # If already accepted, continue
                merged = f"{stdout}\n{stderr}".lower()
                if "already" in merged and "accept" in merged:
                    self._log(f"  ℹ️ ToS already accepted: {ch}")
                    continue
                self._log(f"  ❌ Failed accepting ToS for: {ch}")
                if stderr:
                    self._log(f"  {stderr[:400]}")
                raise Exception(f"Conda ToS acceptance failed for channel: {ch}")




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

            # Ensure conda ToS is accepted before any non-interactive conda action
            self._ensure_conda_tos_accepted()

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

                # 验证 APBS 是否安装成功（EC 分析的核心依赖）
                apbs_check = subprocess.run(
                    [self.conda_exe, "run", "-p", self.env_path, "apbs", "--version"],
                    capture_output=True, text=True, timeout=30
                )
                if apbs_check.returncode == 0:
                    self._log("  ✅ APBS verified")
                else:
                    self._log("  ⚠️ APBS not found — EC analysis may not work")
                    self._log("  💡 Try: pip install apbs-binary (in glint env)")
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


class NativeOSXInstaller:
    """macOS native installer using osascript dialogs.

    When the system Python only has Tk 8.5, tkinter cannot render any
    widgets on modern macOS (Monterey+).  This class bypasses tkinter
    entirely and drives the installation through native macOS dialogs
    via ``osascript`` (AppleScript).

    All installation logic (conda env, packages, file copy, shortcut)
    reuses the same constants and helper functions as InstallerApp.
    """

    def __init__(self):
        self.conda_exe = find_conda()
        self.env_path = None

    # ── osascript helpers ─────────────────────────────────────────

    @staticmethod
    def _osascript_dialog(title, message, buttons=None, default_button=None,
                          icon="note"):
        """Show a native macOS dialog and return the button text clicked.

        Parameters
        ----------
        title : str
            Window title (used as the alert heading).
        message : str
            Body text.
        buttons : list[str] or None
            Button labels.  Defaults to ``["OK"]``.
        default_button : str or None
            Which button is the default (Enter key).
        icon : str
            "note", "caution", or "stop".

        Returns
        -------
        str or None
            The button label that was clicked, or *None* if the user
            cancelled / osascript failed.
        """
        if buttons is None:
            buttons = ["OK"]
        # Build AppleScript button list:  {"Cancel", "Install"}
        btn_list = "{" + ", ".join(f'"{b}"' for b in buttons) + "}"
        default_part = ""
        if default_button:
            default_part = f' default button "{default_button}"'

        # Escape backslashes and double-quotes inside message / title
        safe_title = title.replace("\\", "\\\\").replace('"', '\\"')
        safe_msg = message.replace("\\", "\\\\").replace('"', '\\"')
        # Convert Python newlines to AppleScript linefeed
        safe_msg = safe_msg.replace("\n", "\\n")

        script = (
            f'display dialog "{safe_msg}" '
            f'with title "{safe_title}" '
            f'buttons {btn_list}{default_part} '
            f'with icon {icon}'
        )
        try:
            result = subprocess.run(
                ["osascript", "-e", script],
                capture_output=True, text=True, timeout=600,
            )
            if result.returncode != 0:
                return None
            # stdout looks like: "button returned:Install\n"
            out = result.stdout.strip()
            if "button returned:" in out:
                return out.split("button returned:")[-1].strip()
            return out or buttons[-1]
        except Exception:
            return None

    @staticmethod
    def _osascript_choose_folder(prompt, default_path=None):
        """Show a native folder-chooser dialog.

        Returns the chosen POSIX path string, or *None* on cancel.
        """
        safe_prompt = prompt.replace("\\", "\\\\").replace('"', '\\"')
        default_part = ""
        if default_path:
            safe_default = default_path.replace("\\", "\\\\").replace('"', '\\"')
            default_part = f' default location POSIX file "{safe_default}"'
        script = f'POSIX path of (choose folder with prompt "{safe_prompt}"{default_part})'
        try:
            result = subprocess.run(
                ["osascript", "-e", script],
                capture_output=True, text=True, timeout=600,
            )
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout.strip().rstrip("/")
            return None
        except Exception:
            return None

    @staticmethod
    def _osascript_input(title, message, default_answer=""):
        """Show a dialog with a text input field.

        Returns the entered string, or *None* on cancel.
        """
        safe_title = title.replace("\\", "\\\\").replace('"', '\\"')
        safe_msg = message.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")
        safe_default = default_answer.replace("\\", "\\\\").replace('"', '\\"')
        script = (
            f'text returned of (display dialog "{safe_msg}" '
            f'with title "{safe_title}" '
            f'default answer "{safe_default}" '
            f'buttons {{"Cancel", "OK"}} default button "OK")'
        )
        try:
            result = subprocess.run(
                ["osascript", "-e", script],
                capture_output=True, text=True, timeout=600,
            )
            if result.returncode == 0:
                return result.stdout.strip()
            return None
        except Exception:
            return None

    @staticmethod
    def _osascript_notification(message, title="GLINT Installer"):
        """Post a macOS notification (non-blocking)."""
        safe_title = title.replace("\\", "\\\\").replace('"', '\\"')
        safe_msg = message.replace("\\", "\\\\").replace('"', '\\"')
        script = (
            f'display notification "{safe_msg}" with title "{safe_title}"'
        )
        try:
            subprocess.Popen(
                ["osascript", "-e", script],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            )
        except Exception:
            pass

    # ── Main flow ─────────────────────────────────────────────────

    def run(self):
        """Entry point — drives the full install flow via native dialogs."""
        print("GLINT Installer — Native macOS mode (Tk < 8.6)")

        # ① Welcome
        btn = self._osascript_dialog(
            "GLINT Installer",
            "Welcome to GLINT Installer\n"
            "PyMOL Plugin for Molecular Glue Analysis\n\n"
            "This will install GLINT and its dependencies.\n"
            "The process may take 10–20 minutes.",
            buttons=["Cancel", "Install"],
            default_button="Install",
        )
        if btn != "Install":
            print("User cancelled.")
            return

        # ② Check conda
        if not self.conda_exe:
            btn = self._osascript_dialog(
                "Conda Not Found",
                "GLINT requires Conda (Miniforge / Miniconda) to manage "
                "its Python environment and dependencies.\n\n"
                "Please install Miniforge first, then relaunch this installer.\n\n"
                "Click \"Download\" to open the Miniforge download page.",
                buttons=["Cancel", "Download"],
                default_button="Download",
                icon="caution",
            )
            if btn == "Download":
                webbrowser.open("https://docs.conda.io/en/latest/miniconda.html")
            print("❌ Conda not found. Cannot proceed.")
            return

        # Report conda version
        try:
            ver_result = subprocess.run(
                [self.conda_exe, "--version"],
                capture_output=True, text=True, timeout=10,
            )
            conda_ver = ver_result.stdout.strip() if ver_result.returncode == 0 else "unknown"
        except Exception:
            conda_ver = "unknown"
        print(f"✅ Conda found: {self.conda_exe} ({conda_ver})")

        # Determine env path
        try:
            base_result = subprocess.run(
                [self.conda_exe, "info", "--base"],
                capture_output=True, text=True, timeout=10,
            )
            conda_base = base_result.stdout.strip()
            if not conda_base:
                raise Exception("Cannot determine conda base")
            self.env_path = os.path.join(conda_base, "envs", ENV_NAME)
        except Exception as exc:
            self._osascript_dialog(
                "Error",
                f"Cannot determine conda base directory:\n{exc}",
                icon="stop",
            )
            return

        env_exists = os.path.isdir(self.env_path)

        # ③ Installation path
        install_path = DEFAULT_INSTALL_PATH
        path_answer = self._osascript_input(
            "Installation Path",
            "Where should GLINT be installed?\n"
            "(PyMOL loads plugins from ~/.pymol/startup/)\n\n"
            "Edit the path below or accept the default:",
            default_answer=install_path,
        )
        if path_answer is None:
            print("User cancelled path selection.")
            return
        install_path = path_answer.strip() or install_path

        # ④ Confirm
        source_dir = get_glint_source_dir()
        summary = (
            f"Conda: {conda_ver}\n"
            f"Environment: {ENV_NAME} ({'exists' if env_exists else 'will be created'})\n"
            f"Install path: {install_path}\n"
            f"Source: {source_dir or 'bundled'}\n\n"
            "Click \"Install\" to start. This may take 10–20 minutes."
        )
        btn = self._osascript_dialog(
            "Confirm Installation",
            summary,
            buttons=["Cancel", "Install"],
            default_button="Install",
        )
        if btn != "Install":
            print("User cancelled.")
            return

        # ⑤ Run installation
        self._do_install(install_path, source_dir, env_exists)

    def _ensure_conda_tos_accepted(self):
        """Accept conda Terms of Service for non-interactive install."""
        print("  Checking conda Terms of Service...")
        for ch in CONDA_TOS_CHANNELS:
            cmd = [
                self.conda_exe, "tos", "accept",
                "--override-channels",
                "--channel", ch,
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            if result.returncode == 0:
                print(f"  ✅ ToS accepted: {ch}")
            else:
                stderr = (result.stderr or "").strip()
                stdout = (result.stdout or "").strip()
                merged = f"{stdout}\n{stderr}".lower()
                if "already" in merged and "accept" in merged:
                    print(f"  ℹ️ ToS already accepted: {ch}")
                    continue
                print(f"  ❌ Failed accepting ToS for: {ch}")
                if stderr:
                    print(f"  {stderr[:400]}")
                raise Exception(f"Conda ToS acceptance failed for channel: {ch}")

    def _do_install(self, install_path, source_dir, env_exists):
        """Execute the full installation (blocking)."""
        try:
            self._osascript_notification("Starting installation…")

            # 1. Conda environment
            print("\n" + "=" * 50)
            print("[1/5] Checking conda environment...")
            print("=" * 50)

            self._ensure_conda_tos_accepted()

            if not env_exists:
                print(f"  Creating environment at {self.env_path}...")
                self._osascript_notification("Creating conda environment…")
                result = subprocess.run(
                    [self.conda_exe, "create", "-p", self.env_path,
                     f"python={PYTHON_VERSION}", "-y"],
                    capture_output=True, text=True, timeout=600,
                )
                if result.returncode != 0:
                    raise Exception(
                        f"Failed to create conda environment:\n{result.stderr[:500]}"
                    )
                print("  ✅ Environment created")
            else:
                print(f"  Environment exists: {self.env_path}")

            # 2. Install dependencies
            print("\n" + "=" * 50)
            print("[2/5] Installing dependencies...")
            print("=" * 50)
            print("  This may take 10-20 minutes, please wait...")
            self._osascript_notification("Installing conda packages (this takes a while)…")

            pkg_str = " ".join(CONDA_PACKAGES)
            print(f"  Installing: {pkg_str}")

            cmd = [
                self.conda_exe, "install",
                "-p", self.env_path,
                "-c", "conda-forge",
                "-c", "schrodinger",
                "-y",
            ] + CONDA_PACKAGES

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            if result.returncode == 0:
                print("  ✅ Conda packages installed")
            else:
                print("  ⚠️ Some packages may have failed")
                if result.stderr:
                    print(f"  {result.stderr[:500]}")

            # Pip packages
            if PIP_PACKAGES:
                print("\n  Installing pip packages...")
                self._osascript_notification("Installing pip packages…")
                pip_cmd = [
                    self.conda_exe, "run", "-p", self.env_path,
                    "python", "-m", "pip", "install", "--quiet",
                ] + PIP_PACKAGES
                pip_result = subprocess.run(
                    pip_cmd, capture_output=True, text=True, timeout=600,
                )
                if pip_result.returncode == 0:
                    print("  ✅ Pip packages installed")
                else:
                    print("  ⚠️ Some pip packages may have failed")

            # Verify APBS
            apbs_check = subprocess.run(
                [self.conda_exe, "run", "-p", self.env_path, "apbs", "--version"],
                capture_output=True, text=True, timeout=30,
            )
            if apbs_check.returncode == 0:
                print("  ✅ APBS verified")
            else:
                print("  ⚠️ APBS not found — EC analysis may not work")

            # 3. Copy plugin files
            print("\n" + "=" * 50)
            print("[3/5] Installing plugin files...")
            print("=" * 50)
            self._osascript_notification("Copying plugin files…")

            if source_dir and os.path.isdir(source_dir):
                os.makedirs(install_path, exist_ok=True)

                # Remove old installation
                if os.path.exists(install_path):
                    for item in os.listdir(install_path):
                        item_path = os.path.join(install_path, item)
                        try:
                            if os.path.isdir(item_path):
                                shutil.rmtree(item_path)
                            else:
                                os.remove(item_path)
                        except Exception:
                            pass

                IGNORED = {'__pycache__', '.DS_Store', '.git', '.gitignore', '*.pyc'}
                for item in os.listdir(source_dir):
                    if item in IGNORED:
                        continue
                    src = os.path.join(source_dir, item)
                    dst = os.path.join(install_path, item)
                    try:
                        if os.path.isdir(src):
                            shutil.copytree(
                                src, dst,
                                ignore=shutil.ignore_patterns(*IGNORED),
                            )
                        else:
                            shutil.copy2(src, dst)
                    except Exception as exc:
                        print(f"  ⚠️ Failed to copy {item}: {exc}")

                print(f"  ✅ Copied to {install_path}")
            else:
                print(f"  ⚠️ Source directory not found: {source_dir}")

            # 4. Create shortcut
            print("\n" + "=" * 50)
            print("[4/5] Creating desktop shortcut...")
            print("=" * 50)
            self._create_shortcut(install_path)

            # 5. Done
            print("\n" + "=" * 50)
            print("[5/5] ✅ Installation complete!")
            print("=" * 50)

            self._osascript_dialog(
                "Installation Complete",
                "GLINT has been installed successfully!\n\n"
                "How to use GLINT:\n"
                "  1. Open GLINT from ~/Applications/GLINT.app\n"
                "  2. Or run in PyMOL: glint_gui",
                buttons=["OK"],
                icon="note",
            )

        except Exception as exc:
            import traceback
            print(f"\n❌ Error: {exc}")
            traceback.print_exc()
            self._osascript_dialog(
                "Installation Failed",
                f"An error occurred during installation:\n\n{exc}",
                buttons=["OK"],
                icon="stop",
            )

    def _create_shortcut(self, install_path):
        """Create macOS .app bundle shortcut (same logic as InstallerApp)."""
        home = os.path.expanduser("~")
        app_dir = os.path.join(home, "Applications")
        os.makedirs(app_dir, exist_ok=True)

        app_name = "GLINT.app"
        app_path = os.path.join(app_dir, app_name)

        contents_dir = os.path.join(app_path, "Contents")
        macos_dir = os.path.join(contents_dir, "MacOS")
        resources_dir = os.path.join(contents_dir, "Resources")

        os.makedirs(macos_dir, exist_ok=True)
        os.makedirs(resources_dir, exist_ok=True)

        launcher_script = os.path.join(macos_dir, "GLINT")
        with open(launcher_script, "w") as f:
            f.write("#!/bin/bash\n")
            f.write(f'eval "$({self.conda_exe} shell.bash hook)"\n')
            f.write(f'conda activate {ENV_NAME}\n')
            f.write('export KMP_DUPLICATE_LIB_OK=TRUE\n')
            f.write('export OMP_NUM_THREADS=1\n')
            f.write(
                'pymol -d "import sys, os; sys.path.insert(0, '
                "os.path.expanduser('~/.pymol/startup')); "
                'import glint; glint.glint_gui()"\n'
            )
        os.chmod(launcher_script, 0o755)

        plist_path = os.path.join(contents_dir, "Info.plist")
        with open(plist_path, "w") as f:
            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write('<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" '
                    '"http://www.apple.com/DTDs/PropertyList-1.0.dtd">\n')
            f.write('<plist version="1.0">\n<dict>\n')
            f.write('    <key>CFBundleExecutable</key>\n    <string>GLINT</string>\n')
            f.write('    <key>CFBundleIconFile</key>\n    <string>AppIcon</string>\n')
            f.write('    <key>CFBundleIdentifier</key>\n    <string>com.glint.app</string>\n')
            f.write('    <key>CFBundleName</key>\n    <string>GLINT</string>\n')
            f.write('    <key>CFBundlePackageType</key>\n    <string>APPL</string>\n')
            f.write('    <key>CFBundleShortVersionString</key>\n    <string>1.0</string>\n')
            f.write('</dict>\n</plist>\n')

        icon_src = os.path.join(install_path, "assets", "logo.png")
        if os.path.exists(icon_src):
            shutil.copy2(icon_src, os.path.join(resources_dir, "AppIcon.png"))

        print(f"  ✅ Created application: {app_path}")
        print(f"  You can find GLINT in ~/Applications/")


def _activate_foreground():
    """Bring the tkinter window to the foreground on macOS.
    When launched from a .app bundle via 'open', the process may start
    in the background. This uses osascript to bring it forward.
    Must be called AFTER tk.Tk() to avoid NSException conflicts."""
    try:
        subprocess.Popen([
            'osascript', '-e',
            'tell application "System Events" to set frontmost of '
            'first process whose unix id is ' + str(os.getpid()) + ' to true'
        ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        pass


class CLIInstaller:
    """Minimal CLI fallback installer for when the GUI cannot render."""

    def __init__(self):
        self.conda_exe = find_conda()
        self.env_path = None

    def run(self):
        print("=" * 55)
        print("  GLINT Installer (CLI mode)")
        print("  PyMOL Plugin for Molecular Glue Analysis")
        print("=" * 55)

        # Check conda
        if not self.conda_exe:
            print("\n❌ Conda not found!")
            print("Please install Miniconda first:")
            print("  https://docs.conda.io/en/latest/miniconda.html")
            sys.exit(1)

        try:
            result = subprocess.run([self.conda_exe, "--version"],
                                    capture_output=True, text=True, timeout=10)
            print(f"\n✅ Conda found: {result.stdout.strip()}")
        except Exception as e:
            print(f"\n❌ Conda error: {e}")
            sys.exit(1)

        # Determine env path
        try:
            base_result = subprocess.run([self.conda_exe, "info", "--base"],
                                         capture_output=True, text=True, timeout=10)
            conda_base = base_result.stdout.strip()
            self.env_path = os.path.join(conda_base, "envs", ENV_NAME)
        except Exception as e:
            print(f"❌ Cannot determine conda base: {e}")
            sys.exit(1)

        env_exists = os.path.isdir(self.env_path)
        print(f"  Environment '{ENV_NAME}': {'exists' if env_exists else 'will be created'}")

        install_path = DEFAULT_INSTALL_PATH
        source_dir = get_glint_source_dir()

        print(f"\n📁 Install path: {install_path}")
        print(f"📦 Source: {source_dir or 'not found'}")

        # Confirm
        try:
            answer = input("\nProceed with installation? [Y/n] ").strip().lower()
        except EOFError:
            answer = "y"
        if answer and answer != "y":
            print("Cancelled.")
            return

        # Create environment
        if not env_exists:
            print(f"\n[1/4] Creating conda environment...")
            result = subprocess.run(
                [self.conda_exe, "create", "-p", self.env_path,
                 f"python={PYTHON_VERSION}", "-y"],
                timeout=600
            )
            if result.returncode != 0:
                print("❌ Failed to create environment")
                sys.exit(1)
            print("✅ Environment created")
        else:
            print(f"\n[1/4] Environment exists: {self.env_path}")

        # Install packages
        print("\n[2/4] Installing dependencies (this may take 10-20 min)...")
        cmd = [self.conda_exe, "install", "-p", self.env_path,
               "-c", "conda-forge", "-c", "schrodinger", "-y"] + CONDA_PACKAGES
        result = subprocess.run(cmd, timeout=3600)
        if result.returncode == 0:
            print("✅ Conda packages installed")
        else:
            print("⚠️ Some conda packages may have failed")

        if PIP_PACKAGES:
            print("  Installing pip packages...")
            pip_cmd = [self.conda_exe, "run", "-p", self.env_path,
                       "python", "-m", "pip", "install", "--quiet"] + PIP_PACKAGES
            subprocess.run(pip_cmd, timeout=600)

        # Copy files
        print("\n[3/4] Installing plugin files...")
        if source_dir and os.path.isdir(source_dir):
            os.makedirs(install_path, exist_ok=True)
            IGNORED = {'__pycache__', '.DS_Store', '.git', '.gitignore', '*.pyc'}
            for item in os.listdir(source_dir):
                if item in IGNORED:
                    continue
                src = os.path.join(source_dir, item)
                dst = os.path.join(install_path, item)
                try:
                    if os.path.isdir(src):
                        if os.path.exists(dst):
                            shutil.rmtree(dst)
                        shutil.copytree(src, dst,
                                        ignore=shutil.ignore_patterns(*IGNORED))
                    else:
                        shutil.copy2(src, dst)
                except Exception as e:
                    print(f"  ⚠️ Failed to copy {item}: {e}")
            print(f"✅ Copied to {install_path}")
        else:
            print(f"⚠️ Source directory not found: {source_dir}")

        print("\n[4/4] ✅ Installation complete!")
        print("\nHow to use GLINT:")
        print("  1. Open PyMOL from the glint conda environment")
        print("  2. Run: glint_gui")


def main():
    # When launched from .app bundle, remove bundle identifier
    # to prevent tkinter NSApplication conflicts
    if sys.platform == 'darwin':
        os.environ.pop('__CFBundleIdentifier', None)
        os.environ['TK_SILENCE_DEPRECATION'] = '1'

    # CLI mode: explicit flag or no tkinter available
    if '--cli' in sys.argv or not HAS_TK:
        if not HAS_TK:
            print("tkinter not available.")
        # On macOS without tkinter, try native osascript dialogs first
        if sys.platform == 'darwin' and not ('--cli' in sys.argv):
            print("Using native macOS installer (osascript dialogs).")
            native = NativeOSXInstaller()
            native.run()
            return
        print("Using CLI installer.")
        cli = CLIInstaller()
        cli.run()
        return

    # macOS + Tk < 8.6: use native osascript dialogs instead of tkinter.
    # Tk 8.5 (shipped with macOS system Python) cannot render any widgets
    # on modern macOS (Monterey+), even plain tk.* widgets.
    if sys.platform == 'darwin' and _tk_version < 8.6:
        print(f"Tk {_tk_version} detected — using native macOS installer.")
        native = NativeOSXInstaller()
        native.run()
        return

    # Try the full tkinter GUI (requires Tk >= 8.6)
    try:
        root = tk.Tk()

        # Quick rendering sanity check on macOS: create a test label
        # and verify the window can be drawn. Some Tk builds report
        # version 8.6 but still fail to render on modern macOS.
        if sys.platform == 'darwin':
            try:
                _test = tk.Label(root, text="test")
                _test.pack()
                root.update_idletasks()
                _test.destroy()
            except Exception:
                root.destroy()
                raise RuntimeError("Tk rendering test failed")

        # On macOS, bring the window to the foreground after creation
        if sys.platform == 'darwin':
            try:
                root.attributes('-topmost', True)
                root.after(100, lambda: root.attributes('-topmost', False))
            except Exception:
                pass
            root.after(200, _activate_foreground)
            root.lift()
            root.focus_force()

        app = InstallerApp(root)
        root.mainloop()
    except Exception as e:
        print(f"\nGUI failed to start: {e}")
        # On macOS, fall back to native osascript dialogs (better UX than CLI)
        if sys.platform == 'darwin':
            print("Falling back to native macOS installer...\n")
            try:
                native = NativeOSXInstaller()
                native.run()
                return
            except Exception as e2:
                print(f"Native installer also failed: {e2}")
        print("Falling back to CLI installer...\n")
        cli = CLIInstaller()
        cli.run()


if __name__ == "__main__":
    main()