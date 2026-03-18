#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GLINT Windows Installer - 图形化安装程序
Professional GUI installer for Windows, similar to macOS version
"""

import os
import sys
import subprocess
import shutil
import threading
import webbrowser
from pathlib import Path
from datetime import datetime
import traceback

try:
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
except ImportError:
    print("Error: tkinter not found")
    sys.exit(1)

# 配置
ENV_NAME = "glint"
PYTHON_VERSION = "3.10"  # 更新为 3.10
DEFAULT_INSTALL_PATH = os.path.join(os.path.expanduser("~"), ".pymol", "startup", "glint")

# Conda 依赖包 (与 macOS 版本保持一致)
# 注意: haddock_biobb 已移除，因为 haddocking channel 不可用 (HTTP 404)
CONDA_PACKAGES = [
    "rdkit", "scipy", "matplotlib", "pillow", "numpy=1.26.4",  # 指定 numpy 版本以兼容 PyMOL
    "pandas", "seaborn", "pyqt", "openbabel", "pymol-open-source",
    "meeko", "scikit-image",
    "pdb2pqr",
]

# Pip 包 (open3d 在 conda 上不稳定, haddock3 因 haddocking channel 不可用改用 pip)
# vina: conda-forge 当前无 win-64 构建，因此改用 pip 安装
PIP_PACKAGES = ["requests", "open3d", "haddock3", "vina"]

# APBS 1.5 预编译二进制文件下载地址
APBS_DOWNLOAD_URLS = {
    "darwin_x86_64": "https://github.com/Electrostatics/apbs-pdb2pqr/releases/download/vAPBS-1.5.0/APBS-1.5-osx.zip",
    "darwin_arm64": "https://github.com/Electrostatics/apbs-pdb2pqr/releases/download/vAPBS-1.5.0/APBS-1.5-osx.zip",
    "linux_x86_64": "https://github.com/Electrostatics/apbs-pdb2pqr/releases/download/vAPBS-1.5.0/APBS-1.5-linux64.tar.gz",
    "windows_x86_64": "https://github.com/Electrostatics/apbs-pdb2pqr/releases/download/vAPBS-1.5.0/APBS-1.5-win64.zip",
}


def get_glint_source_dir():
    """获取 GLINT 源码目录"""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. 优先查找同级目录的 glint
    bundled_dir = os.path.join(script_dir, "glint")
    if os.path.isdir(bundled_dir):
        return bundled_dir

    # 2. PyInstaller 打包后的临时目录
    if hasattr(sys, '_MEIPASS'):
        bundled_dir = os.path.join(sys._MEIPASS, "glint")
        if os.path.isdir(bundled_dir):
            return bundled_dir

    return None


class InstallerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("GLINT Installer")
        self.root.geometry("800x900")
        self.root.minsize(750, 850)
        self.root.resizable(True, True)

        # Configure better fonts and styling
        self._configure_styles()

        # Set window icon
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

        # Delayed environment check
        self.root.after(500, lambda: threading.Thread(target=self._check_environment, daemon=True).start())

    def _configure_styles(self):
        """Configure better fonts and styles for Windows"""
        style = ttk.Style()

        # Try to use a modern theme
        available_themes = style.theme_names()
        if 'vista' in available_themes:
            style.theme_use('vista')
        elif 'winnative' in available_themes:
            style.theme_use('winnative')
        elif 'clam' in available_themes:
            style.theme_use('clam')

        # Define better fonts - use Microsoft YaHei for Chinese support, fallback to Segoe UI
        # These fonts look much better on Windows
        self.title_font = ("Microsoft YaHei UI", 22, "bold")
        self.subtitle_font = ("Microsoft YaHei UI", 11)
        self.normal_font = ("Microsoft YaHei UI", 10)
        self.small_font = ("Microsoft YaHei UI", 9)
        self.mono_font = ("Cascadia Code", 9)  # Better monospace font

        # Fallback fonts if Microsoft YaHei is not available
        try:
            import tkinter.font as tkfont
            available_fonts = tkfont.families()

            if "Microsoft YaHei UI" not in available_fonts:
                if "Microsoft YaHei" in available_fonts:
                    self.title_font = ("Microsoft YaHei", 22, "bold")
                    self.subtitle_font = ("Microsoft YaHei", 11)
                    self.normal_font = ("Microsoft YaHei", 10)
                    self.small_font = ("Microsoft YaHei", 9)
                elif "Segoe UI" in available_fonts:
                    self.title_font = ("Segoe UI", 22, "bold")
                    self.subtitle_font = ("Segoe UI", 11)
                    self.normal_font = ("Segoe UI", 10)
                    self.small_font = ("Segoe UI", 9)

            if "Cascadia Code" not in available_fonts:
                if "Consolas" in available_fonts:
                    self.mono_font = ("Consolas", 9)
                else:
                    self.mono_font = ("Courier New", 9)
        except Exception:
            # 字体枚举失败时使用安全默认字体（避免在部分打包/Windows 环境中导致 GUI 渲染异常）
            self.title_font = ("Segoe UI", 22, "bold")
            self.subtitle_font = ("Segoe UI", 11)
            self.normal_font = ("Segoe UI", 10)
            self.small_font = ("Segoe UI", 9)
            self.mono_font = ("Consolas", 9)

        # Configure ttk styles with better fonts
        style.configure("TLabel", font=self.normal_font)
        style.configure("TButton", font=self.normal_font, padding=6)
        style.configure("TCheckbutton", font=self.normal_font)
        style.configure("TEntry", font=self.normal_font)
        style.configure("TLabelframe", font=self.normal_font)
        style.configure("TLabelframe.Label", font=self.normal_font)

        # Configure root window background
        self.root.configure(bg='#f0f0f0')

    def _set_icon(self):
        """设置窗口图标"""
        try:
            icon_path = os.path.join(os.path.dirname(__file__), "glint", "assets", "logo.png")
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
        ttk.Button(path_row, text="Browse...", command=self._browse_path, width=10).pack(side=tk.RIGHT, padx=(10, 0))

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

        # Button section (pack before progress so buttons are always visible)
        btn_frame = ttk.Frame(main)
        btn_frame.pack(fill=tk.X, side=tk.BOTTOM, pady=(5, 0))

        self.install_btn = ttk.Button(btn_frame, text="🚀 Install GLINT",
                                      command=self._start_install, width=20)
        self.install_btn.pack(side=tk.RIGHT, padx=(10, 0))

        ttk.Button(btn_frame, text="Cancel", command=self.root.quit, width=12).pack(side=tk.RIGHT)

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

    def _install_apbs_binary(self):
        """下载并安装 APBS 1.5 预编译二进制文件到 conda 环境"""
        import platform
        import urllib.request
        import zipfile
        import tarfile
        import tempfile

        self._log("  Installing APBS 1.5 binary...")

        # 确定平台
        system = platform.system().lower()
        machine = platform.machine().lower()

        if system == "darwin":
            if machine == "arm64":
                key = "darwin_arm64"
            else:
                key = "darwin_x86_64"
        elif system == "linux":
            key = "linux_x86_64"
        elif system == "windows":
            key = "windows_x86_64"
        else:
            self._log(f"  ⚠️ Unsupported platform: {system} {machine}")
            return False

        url = APBS_DOWNLOAD_URLS.get(key)
        if not url:
            self._log(f"  ⚠️ No APBS binary available for {key}")
            return False

        try:
            # 下载到临时目录
            temp_dir = tempfile.mkdtemp()
            filename = os.path.basename(url)
            download_path = os.path.join(temp_dir, filename)

            self._log(f"  Downloading from {url}...")
            urllib.request.urlretrieve(url, download_path)

            # 解压
            extract_dir = os.path.join(temp_dir, "apbs_extract")
            os.makedirs(extract_dir, exist_ok=True)

            if filename.endswith('.zip'):
                with zipfile.ZipFile(download_path, 'r') as zf:
                    zf.extractall(extract_dir)
            elif filename.endswith('.tar.gz'):
                with tarfile.open(download_path, 'r:gz') as tf:
                    tf.extractall(extract_dir)

            # 找到 apbs 可执行文件
            apbs_exe = None
            for root, dirs, files in os.walk(extract_dir):
                for f in files:
                    if f == 'apbs' or f == 'apbs.exe':
                        apbs_exe = os.path.join(root, f)
                        break
                if apbs_exe:
                    break

            if not apbs_exe:
                self._log("  ⚠️ APBS executable not found in archive")
                return False

            # 复制到 conda 环境的 bin 目录 (Windows 用 Scripts)
            if system == "windows":
                env_bin = os.path.join(self.env_path, "Scripts")
                dest_apbs = os.path.join(env_bin, "apbs.exe")
            else:
                env_bin = os.path.join(self.env_path, "bin")
                dest_apbs = os.path.join(env_bin, "apbs")
            os.makedirs(env_bin, exist_ok=True)

            shutil.copy2(apbs_exe, dest_apbs)
            if system != "windows":
                os.chmod(dest_apbs, 0o755)

            # 同时复制相关的库文件（如果有）
            apbs_dir = os.path.dirname(apbs_exe)
            lib_dir = os.path.join(apbs_dir, "..", "lib")
            if os.path.exists(lib_dir):
                env_lib = os.path.join(self.env_path, "lib")
                os.makedirs(env_lib, exist_ok=True)
                for item in os.listdir(lib_dir):
                    src = os.path.join(lib_dir, item)
                    dst = os.path.join(env_lib, item)
                    if os.path.isfile(src):
                        shutil.copy2(src, dst)

            # 清理临时文件
            shutil.rmtree(temp_dir, ignore_errors=True)

            self._log(f"  ✅ APBS 1.5 installed to {dest_apbs}")
            return True

        except Exception as e:
            self._log(f"  ⚠️ Failed to install APBS: {e}")
            return False

    def _log(self, msg):
        """线程安全的日志输出，通过 root.after 调度到主线程"""
        def _write():
            self.log_text.configure(state=tk.NORMAL)
            self.log_text.insert(tk.END, msg + "\n")
            self.log_text.see(tk.END)
            self.log_text.configure(state=tk.DISABLED)
        try:
            self.root.after(0, _write)
        except Exception:
            pass


    def _run_cmd_stream(self, cmd, timeout=3600):
        """流式执行命令，逐行输出到 GUI 日志区域，返回 returncode"""
        import time
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            creationflags=subprocess.CREATE_NO_WINDOW if hasattr(subprocess, 'CREATE_NO_WINDOW') else 0
        )
        start = time.time()
        try:
            for line in proc.stdout:
                stripped = line.rstrip('\n\r')
                if stripped:
                    self._log(f"    {stripped}")
                if time.time() - start > timeout:
                    proc.kill()
                    raise subprocess.TimeoutExpired(cmd, timeout)
            proc.wait()
        except subprocess.TimeoutExpired:
            proc.kill()
            raise
        finally:
            if proc.stdout:
                proc.stdout.close()
        return proc.returncode

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
        """查找 Conda 安装路径"""
        # Windows 常见路径
        home = os.path.expanduser("~")
        candidates = [
            os.path.join(home, "miniconda3", "Scripts", "conda.exe"),
            os.path.join(home, "anaconda3", "Scripts", "conda.exe"),
            os.path.join(home, "Miniconda3", "Scripts", "conda.exe"),
            os.path.join(home, "Anaconda3", "Scripts", "conda.exe"),
            os.path.join(home, "AppData", "Local", "miniconda3", "Scripts", "conda.exe"),
            os.path.join(home, "AppData", "Local", "anaconda3", "Scripts", "conda.exe"),
            r"C:\ProgramData\miniconda3\Scripts\conda.exe",
            r"C:\ProgramData\Anaconda3\Scripts\conda.exe",
            r"C:\miniconda3\Scripts\conda.exe",
            r"C:\Anaconda3\Scripts\conda.exe",
        ]

        # 检查 PATH
        try:
            result = subprocess.run(["where", "conda"], capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                conda_path = result.stdout.strip().split("\n")[0]
                if os.path.exists(conda_path):
                    return conda_path
        except:
            pass

        # 检查候选路径
        for path in candidates:
            if os.path.exists(path):
                return path

        return None

    def _install_conda(self):
        """打开 Miniconda 下载页面"""
        webbrowser.open("https://docs.conda.io/en/latest/miniconda.html")
        messagebox.showinfo("Install Miniconda",
                           "Please download and install Miniconda for Windows.\n\n"
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
                returncode = self._run_cmd_stream(
                    [self.conda_exe, "create", "-p", self.env_path, f"python={PYTHON_VERSION}", "-y"],
                    timeout=600
                )
                if returncode != 0:
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

                cmd = [
                    self.conda_exe, "install",
                    "-p", self.env_path,
                    "-c", "conda-forge",
                    "-y"
                ] + CONDA_PACKAGES

                returncode = self._run_cmd_stream(cmd, timeout=3600)

                if returncode == 0:
                    self._log("  ✅ Conda packages installed")
                else:
                    self._log("  ⚠️ Some packages may have failed")

                self.progress["value"] = 50

                # 安装 pip 包
                if PIP_PACKAGES:
                    self._log("\n  Installing pip packages...")
                    pip_cmd = [
                        self.conda_exe, "run", "-p", self.env_path,
                        "python", "-m", "pip", "install"
                    ] + PIP_PACKAGES

                    pip_rc = self._run_cmd_stream(pip_cmd, timeout=600)
                    if pip_rc == 0:
                        self._log("  ✅ Pip packages installed")
                    else:
                        self._log("  ⚠️ Some pip packages may have failed")

                # 安装 APBS 1.5 预编译二进制文件（EC 分析的核心依赖）
                self._log("\n  Installing APBS 1.5...")
                if self._install_apbs_binary():
                    # 验证安装
                    apbs_check = subprocess.run(
                        [self.conda_exe, "run", "-p", self.env_path, "apbs", "--version"],
                        capture_output=True, text=True, timeout=30
                    )
                    if apbs_check.returncode == 0:
                        self._log("  ✅ APBS verified")
                    else:
                        self._log("  ⚠️ APBS installed but verification failed")
                else:
                    self._log("  ⚠️ APBS installation failed — EC analysis may not work")
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
            self._log("  1. Double-click 'GLINT' shortcut on Desktop")
            self._log("  2. Or run in PyMOL: glint_gui")

            self.root.after(0, lambda: messagebox.showinfo(
                "Success",
                "GLINT installed successfully!\n\n"
                "You can now launch GLINT from the Desktop shortcut."))

        except Exception as e:
            self._log(f"\n❌ Error: {e}")
            self._log(traceback.format_exc())
            self.root.after(0, lambda: messagebox.showerror("Error", str(e)))
        finally:
            self.root.after(0, lambda: self.install_btn.configure(state=tk.NORMAL))

    def _create_shortcut(self):
        """Create desktop shortcut with proper icon"""
        home = os.path.expanduser("~")
        desktop = os.path.join(home, "Desktop")

        # Create launcher directory
        launcher_dir = os.path.join(home, ".glint")
        os.makedirs(launcher_dir, exist_ok=True)

        # Get conda base path
        conda_base = os.path.dirname(os.path.dirname(self.conda_exe))
        activate_bat = os.path.join(conda_base, "Scripts", "activate.bat")

        # Create icon file from PNG
        icon_path = self._create_icon(launcher_dir)

        # Create launcher batch file
        launcher_bat = os.path.join(launcher_dir, "launch_glint.bat")
        with open(launcher_bat, "w", encoding="utf-8") as f:
            f.write("@echo off\n")
            f.write(f'call "{activate_bat}" "{self.env_path}"\n')
            # 添加环境变量以解决 OpenMP 错误
            f.write('set KMP_DUPLICATE_LIB_OK=TRUE\n')
            f.write('set OMP_NUM_THREADS=1\n')
            f.write('pymol -d "import sys, os; sys.path.insert(0, os.path.expanduser(\'~/.pymol/startup\')); import glint; glint.glint_gui()"\n')

        # Create VBS script to hide command window
        vbs_file = os.path.join(launcher_dir, "launch_glint.vbs")
        with open(vbs_file, "w", encoding="utf-8") as f:
            f.write('Set WshShell = CreateObject("WScript.Shell")\n')
            f.write(f'WshShell.Run chr(34) & "{launcher_bat}" & chr(34), 0\n')
            f.write('Set WshShell = Nothing\n')

        # Create shortcut using PowerShell
        shortcut_path = os.path.join(desktop, "GLINT.lnk")

        try:
            # Build PowerShell script with icon if available
            icon_line = ""
            if icon_path and os.path.exists(icon_path):
                icon_line = f'$Shortcut.IconLocation = "{icon_path}"'

            ps_script = f'''
$WshShell = New-Object -ComObject WScript.Shell
$Shortcut = $WshShell.CreateShortcut("{shortcut_path}")
$Shortcut.TargetPath = "wscript.exe"
$Shortcut.Arguments = '"{vbs_file}"'
$Shortcut.WorkingDirectory = "{launcher_dir}"
$Shortcut.Description = "GLINT - Molecular Glue Analyzer"
{icon_line}
$Shortcut.Save()
'''
            result = subprocess.run(["powershell", "-Command", ps_script],
                                   capture_output=True, timeout=30)
            if result.returncode == 0:
                self._log(f"  ✅ Created shortcut: {shortcut_path}")
                if icon_path:
                    self._log(f"  ✅ Icon set: {icon_path}")
            else:
                raise Exception(result.stderr.decode() if result.stderr else "Unknown error")
        except Exception as e:
            self._log(f"  ⚠️ Shortcut creation failed: {e}")
            self._log(f"  You can manually run: {launcher_bat}")

    def _create_icon(self, launcher_dir):
        """Create .ico file from PNG logo"""
        try:
            # Find the logo PNG
            source_dir = get_glint_source_dir()
            if not source_dir:
                source_dir = self.install_path.get()

            logo_png = os.path.join(source_dir, "assets", "logo.png")

            # Also check installed path
            if not os.path.exists(logo_png):
                logo_png = os.path.join(self.install_path.get(), "assets", "logo.png")

            if not os.path.exists(logo_png):
                self._log("  ⚠️ Logo PNG not found, using default icon")
                return None

            icon_path = os.path.join(launcher_dir, "glint.ico")

            # Try to use PIL to convert PNG to ICO
            try:
                from PIL import Image
                img = Image.open(logo_png)
                # Create multiple sizes for better quality
                sizes = [(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)]
                img.save(icon_path, format='ICO', sizes=sizes)
                self._log(f"  ✅ Created icon: {icon_path}")
                return icon_path
            except ImportError:
                self._log("  ⚠️ PIL not available, trying alternative method...")

            # Alternative: Use conda environment's PIL
            try:
                convert_script = f'''
import sys
sys.path.insert(0, r"{self.env_path}\\Lib\\site-packages")
from PIL import Image
img = Image.open(r"{logo_png}")
sizes = [(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)]
img.save(r"{icon_path}", format='ICO', sizes=sizes)
print("OK")
'''
                result = subprocess.run(
                    [self.conda_exe, "run", "-p", self.env_path, "python", "-c", convert_script],
                    capture_output=True, text=True, timeout=30
                )
                if result.returncode == 0 and os.path.exists(icon_path):
                    self._log(f"  ✅ Created icon: {icon_path}")
                    return icon_path
            except Exception as e:
                self._log(f"  ⚠️ Icon conversion failed: {e}")

            return None

        except Exception as e:
            self._log(f"  ⚠️ Icon creation error: {e}")
            return None


def main():
    try:
        # 设置 DPI 感知 (Windows 10+)
        try:
            from ctypes import windll
            windll.shcore.SetProcessDpiAwareness(1)
        except:
            pass

        root = tk.Tk()
        app = InstallerApp(root)
        root.mainloop()
    except Exception as e:
        try:
            log_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "installer_error.log")
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(f"[{datetime.now().isoformat()}] Unhandled exception: {repr(e)}\n")
                f.write(traceback.format_exc())
                f.write("\n" + ("-" * 80) + "\n")
        except Exception:
            pass


if __name__ == "__main__":
    main()