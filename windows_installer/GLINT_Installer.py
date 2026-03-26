#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GLINT Windows Installer - GUI Setup Program
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

# Configuration
ENV_NAME = "glint"
PYTHON_VERSION = "3.10"  # Use Python 3.10
DEFAULT_INSTALL_PATH = os.path.join(os.path.expanduser("~"), ".pymol", "startup", "glint")

# Conda packages (aligned with macOS version)
# Note: haddock_biobb removed because haddocking channel is unavailable (HTTP 404)
CONDA_PACKAGES = [
    "rdkit", "scipy", "matplotlib", "pillow", "numpy=1.26.4",  # Pin numpy version for PyMOL compatibility
    "pandas", "seaborn", "pyqt", "openbabel", "pymol-open-source",
    "scikit-image",
    "pdb2pqr",
]

# 镜像源配置
TUNA_CONDA_CHANNELS = [
    "https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main",
    "https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r",
    "https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2",
]
# 清华 PyPI 镜像地址
TUNA_PIP_INDEX_URL = "https://pypi.tuna.tsinghua.edu.cn/simple"
# Pip packages (Windows keeps a smaller, stable set)
PIP_PACKAGES = ["requests", "open3d", "meeko"]


def get_glint_source_dir():
    """Locate the bundled GLINT source directory"""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. Check sibling 'glint' directory
    bundled_dir = os.path.join(script_dir, "glint")
    if os.path.isdir(bundled_dir):
        return bundled_dir

    # 2. PyInstaller bundled temp directory
    if hasattr(sys, '_MEIPASS'):
        bundled_dir = os.path.join(sys._MEIPASS, "glint")
        if os.path.isdir(bundled_dir):
            return bundled_dir

    return None


def get_external_dir():
    """Locate the bundled external directory"""
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 1. Check sibling 'external' directory
    bundled_dir = os.path.join(script_dir, "external")
    if os.path.isdir(bundled_dir):
        return bundled_dir

    # 2. PyInstaller bundled temp directory
    if hasattr(sys, '_MEIPASS'):
        bundled_dir = os.path.join(sys._MEIPASS, "external")
        if os.path.isdir(bundled_dir):
            return bundled_dir

    return None


def _verify_apbs_binary(env_path: str) -> bool:
    """Verify APBS binary in conda environment using file existence and execute permission."""
    dest_apbs = os.path.join(env_path, "Scripts", "apbs.exe")
    return os.path.isfile(dest_apbs) and os.access(dest_apbs, os.X_OK)



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
        self.use_tuna_mirror = tk.BooleanVar(value=False)
        self.conda_ok = False
        self.env_ok = False
        self.conda_exe = None
        self.env_path = None

        self._build_ui()
        self._center_window()

        # Delayed environment check
        self.root.after(500, lambda: threading.Thread(target=self._check_environment, daemon=True).start())

    def _configure_styles(self):
        """Configure fonts and styles for Windows"""
        style = ttk.Style()

        # Try to use a modern theme
        available_themes = style.theme_names()
        if 'vista' in available_themes:
            style.theme_use('vista')
        elif 'winnative' in available_themes:
            style.theme_use('winnative')
        elif 'clam' in available_themes:
            style.theme_use('clam')

        # Default fonts: Segoe UI for UI, Consolas for monospace
        self.title_font = ("Segoe UI", 22, "bold")
        self.subtitle_font = ("Segoe UI", 11)
        self.normal_font = ("Segoe UI", 10)
        self.small_font = ("Segoe UI", 9)
        self.mono_font = ("Consolas", 9)

        # Upgrade monospace font if Cascadia Code is available
        try:
            import tkinter.font as tkfont
            available_fonts = tkfont.families()
            if "Cascadia Code" in available_fonts:
                self.mono_font = ("Cascadia Code", 9)
        except Exception:
            # Font enumeration failed; keep safe defaults
            pass

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
        """Set window icon"""
        try:
            icon_path = os.path.join(os.path.dirname(__file__), "glint", "assets", "logo.png")
            if os.path.exists(icon_path):
                img = tk.PhotoImage(file=icon_path)
                self.root.iconphoto(True, img)
        except:
            pass

    def _center_window(self):
        """Center the window on screen"""
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

        # Dynamically read version from _version.py
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
                             text="GLINT will be installed here. PyMOL loads plugins from ~/.pymol/startup/",
                             font=self.small_font, foreground="#888888")
        path_note.pack(anchor=tk.W, pady=(8, 0))

        # Options section
        opts_frame = ttk.LabelFrame(main, text=" Options ", padding=12)
        opts_frame.pack(fill=tk.X, pady=(0, 15))

        ttk.Checkbutton(opts_frame, text="Install/Update dependencies (conda packages + PyMOL)",
                       variable=self.install_deps).pack(anchor=tk.W, pady=3)
        ttk.Checkbutton(opts_frame, text="Use Tsinghua mirror for conda/pip",
                       variable=self.use_tuna_mirror).pack(anchor=tk.W, pady=3)

        # Button section (pack before progress so buttons are always visible)
        btn_frame = ttk.Frame(main)
        btn_frame.pack(fill=tk.X, side=tk.BOTTOM, pady=(5, 0))

        self.install_btn = ttk.Button(btn_frame, text="Install GLINT",
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
        """Install APBS 1.5 from bundled MSI resources into the conda environment"""
        import tempfile

        self._log("  Installing APBS 1.5 binary...")

        external_dir = get_external_dir()
        if not external_dir:
            self._log("External directory not found, cannot install APBS binary")
            return False

        source_dir = os.path.join(external_dir, "apbs", "apbs-win")
        if not os.path.isdir(source_dir):
            self._log(f"APBS resource directory not found: {source_dir}")
            return False

        msi_path = os.path.join(source_dir, "apb1.5.msi")
        if not os.path.isfile(msi_path):
            self._log(f"APBS MSI not found: {msi_path}")
            return False

        temp_dir = tempfile.mkdtemp(prefix="glint-apbs-")
        try:
            extract_dir = os.path.join(temp_dir, "msi_extract")
            os.makedirs(extract_dir, exist_ok=True)

            creationflags = subprocess.CREATE_NO_WINDOW if hasattr(subprocess, "CREATE_NO_WINDOW") else 0
            result = subprocess.run(
                ["msiexec", "/a", msi_path, f"TARGETDIR={extract_dir}", "/qn"],
                capture_output=True,
                text=True,
                timeout=300,
                creationflags=creationflags,
            )
            if result.returncode != 0:
                self._log(f"Failed to extract APBS MSI (exit code {result.returncode})")
                details = "\n".join(part.strip() for part in [result.stdout or "", result.stderr or ""] if part.strip())
                if details:
                    self._log(details)
                return False

            apbs_exe = None
            for root, _, files in os.walk(extract_dir):
                if "apbs.exe" in files:
                    apbs_exe = os.path.join(root, "apbs.exe")
                    break

            if not apbs_exe:
                self._log("APBS executable not found in extracted MSI contents")
                return False

            source_bin_dir = os.path.dirname(apbs_exe)
            dest_dir = os.path.join(self.env_path, "Scripts")
            os.makedirs(dest_dir, exist_ok=True)

            for item in os.listdir(source_bin_dir):
                source = os.path.join(source_bin_dir, item)
                dest = os.path.join(dest_dir, item)
                if os.path.isdir(source):
                    if os.path.exists(dest):
                        shutil.rmtree(dest)
                    shutil.copytree(source, dest)
                else:
                    shutil.copy2(source, dest)

            dest_apbs = os.path.join(dest_dir, "apbs.exe")
            if not os.path.isfile(dest_apbs):
                self._log("APBS binary not found after installation")
                return False

            self._log(f"APBS installed to {dest_apbs}")
            return True

        except subprocess.TimeoutExpired:
            self._log("APBS MSI extraction timed out")
            return False
        except Exception as e:
            self._log(f"Failed to install APBS: {e}")
            return False
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

    def _log(self, msg):
        """Thread-safe log output, dispatched to main thread via root.after"""
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
        """Run a command with streaming stdout, log each line to GUI, return returncode"""
        import time
        import os

        env = os.environ.copy()
        env['PYTHONIOENCODING'] = 'utf-8'
        env['PYTHONUTF8'] = '1'

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            encoding='utf-8',
            errors='replace',
            bufsize=1,
            creationflags=subprocess.CREATE_NO_WINDOW if hasattr(subprocess, 'CREATE_NO_WINDOW') else 0,
            env=env
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
        """Browse for installation path"""
        path = filedialog.askdirectory(title="Select Installation Directory",
                                       initialdir=os.path.dirname(self.install_path.get()))
        if path:
            self.install_path.set(path)

    def _check_environment(self):
        """Check environment prerequisites"""
        self._log("Checking environment...")

        # 1. Find Conda
        self.conda_exe = self._find_conda()

        if self.conda_exe:
            try:
                result = subprocess.run([self.conda_exe, "--version"],
                                       capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    version = result.stdout.strip()
                    self.root.after(0, lambda: self.conda_status.configure(
                        text=f"{version}", foreground="green"))
                    self.conda_ok = True
                    self.root.after(0, lambda: self.conda_install_btn.configure(state=tk.DISABLED))
                    self._log(f"  Conda found: {self.conda_exe}")
                    self._log(f"  Version: {version}")
                else:
                    raise Exception("conda command failed")
            except Exception as e:
                self.root.after(0, lambda: self.conda_status.configure(
                    text="Error", foreground="red"))
                self._log(f"  Conda error: {e}")
        else:
            self.root.after(0, lambda: self.conda_status.configure(
                text="Not found", foreground="red"))
            self.root.after(0, lambda: self.conda_install_btn.configure(state=tk.NORMAL))
            self._log("  Conda: Not found")
            self._log("  Please install Miniconda first!")
            return

        # 2. Check conda environment
        try:
            base_result = subprocess.run([self.conda_exe, "info", "--base"],
                                        capture_output=True, text=True, timeout=10)
            conda_base = base_result.stdout.strip() if base_result.returncode == 0 else ""

            if not conda_base:
                raise Exception("Cannot determine conda base")

            self.env_path = os.path.join(conda_base, "envs", ENV_NAME)

            if os.path.isdir(self.env_path):
                self.root.after(0, lambda: self.env_status.configure(
                    text="Exists", foreground="green"))
                self.env_ok = True
                self._log(f"  Environment '{ENV_NAME}': {self.env_path}")
            else:
                self.root.after(0, lambda: self.env_status.configure(
                    text="Will create", foreground="orange"))
                self._log(f"  Environment '{ENV_NAME}': Will be created")
        except Exception as e:
            self.root.after(0, lambda: self.env_status.configure(
                text="❓ Unknown", foreground="gray"))
            self._log(f"  Environment check failed: {e}")

        self._log("Ready to install.")

    def _find_conda(self):
        """Find conda executable path"""
        # Common Windows installation paths
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

        # Check PATH
        try:
            result = subprocess.run(["where", "conda"], capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                conda_path = result.stdout.strip().split("\n")[0]
                if os.path.exists(conda_path):
                    return conda_path
        except:
            pass

        # Check candidate paths
        for path in candidates:
            if os.path.exists(path):
                return path

        return None

    def _install_conda(self):
        """Open Miniconda download page"""
        webbrowser.open("https://docs.conda.io/en/latest/miniconda.html")
        messagebox.showinfo("Install Miniconda",
                           "Please download and install Miniconda for Windows.\n\n"
                           "After installation, restart this installer.")

    def _install_vina_binary(self):
        """Install Vina binary into the conda environment"""
        external_dir = get_external_dir()
        if not external_dir:
            self._log("External directory not found, cannot install Vina binary")
            return False

        source = os.path.join(external_dir, "vina", "vina_1.2.7_win.exe")
        if not os.path.exists(source):
            self._log(f"Vina binary not found: {source}")
            return False

        try:
            dest_dir = os.path.join(self.env_path, "Scripts")
            os.makedirs(dest_dir, exist_ok=True)
            dest = os.path.join(dest_dir, "vina.exe")
            shutil.copy2(source, dest)
            self._log(f"Vina installed to {dest}")
            return True
        except Exception as e:
            self._log(f"Failed to install Vina: {e}")
            return False

    def _start_install(self):
        """Start installation"""
        if not self.conda_ok:
            messagebox.showerror("Error", "Please install Miniconda first!")
            return

        self.install_btn.configure(state=tk.DISABLED)
        threading.Thread(target=self._do_install, daemon=True).start()

    def _do_install(self):
        """Perform the installation"""
        try:
            self.progress["value"] = 0

            # 1. Create / check conda environment
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


            # Windows 下 conda 版本差异较大，旧版不一定支持自动 ToS 配置
            # 这里直接跳过，避免 'auto_accept_default_terms' 报错
            if not self.env_ok:
                self._log(f"  Creating environment at {self.env_path}...")
                returncode = self._run_cmd_stream(
                    [self.conda_exe, "create", "-p", self.env_path, f"python={PYTHON_VERSION}", "-y"],
                    timeout=600
                )
                if returncode != 0:
                    raise Exception("Failed to create conda environment")
                self._log(" Environment created")
            else:
                self._log(f"  Environment exists: {self.env_path}")

            self.progress["value"] = 20

            # 2. Install dependencies
            if self.install_deps.get():
                self._log("\n" + "=" * 50)
                self._log("[2/5] Installing dependencies...")
                self._log("=" * 50)

                # 使用 globals() 兜底，兼容旧打包产物或异常全局状态，避免 NameError 中断安装
                pip_packages = globals().get("PIP_PACKAGES", ["requests", "open3d", "meeko"])
                tuna_pip_index_url = globals().get("TUNA_PIP_INDEX_URL", "https://pypi.tuna.tsinghua.edu.cn/simple")
                total_pkgs = len(CONDA_PACKAGES) + len(pip_packages)
                installed_count = 0
                failed_pkgs = []


                # Progress range for dependency step: 20 -> 55
                progress_start = 20
                progress_end = 55

                # --- Conda packages (one by one) ---
                self._log(f"\n  Conda packages to install: {len(CONDA_PACKAGES)}")
                for idx, pkg in enumerate(CONDA_PACKAGES, 1):
                    overall_idx = installed_count + 1
                    self._log(
                        f"\n  Installing (overall {overall_idx}/{total_pkgs}) via conda: {pkg} ({idx}/{len(CONDA_PACKAGES)})..."
                    )
                    conda_channels = ["-c", "conda-forge"]
                    if self.use_tuna_mirror.get():
                        conda_channels = [
                            "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge",
                            "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main",
                            "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r",
                            "-c", "https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2",
                        ]
                    cmd = [
                        self.conda_exe, "install",
                        "-p", self.env_path,
                    ] + conda_channels + ["-y", pkg]
                    try:
                        rc = self._run_cmd_stream(cmd, timeout=600)
                        if rc == 0:
                            self._log(f" {pkg} installed")
                        else:
                            self._log(f" {pkg} failed (exit code {rc})")
                            failed_pkgs.append(pkg)
                    except subprocess.TimeoutExpired:
                        self._log(f" {pkg} timed out")
                        failed_pkgs.append(pkg)
                    except Exception as e:
                        self._log(f" {pkg} error: {e}")
                        failed_pkgs.append(pkg)

                    installed_count += 1
                    pct = progress_start + (progress_end - progress_start) * installed_count / total_pkgs
                    self.progress["value"] = pct

                # --- Pip packages (one by one) ---
                if pip_packages:
                    self._log(f"\n  Pip packages to install: {len(pip_packages)}")
                    for idx, pkg in enumerate(pip_packages, 1):
                        overall_idx = installed_count + 1
                        self._log(
                            f"\n  Installing (overall {overall_idx}/{total_pkgs}) via pip: {pkg} ({idx}/{len(pip_packages)})..."
                        )
                        pip_cmd = [
                            self.conda_exe, "run", "-p", self.env_path,
                            "python", "-m", "pip", "install"
                        ]
                        if self.use_tuna_mirror.get():
                            pip_cmd.extend(["-i", tuna_pip_index_url, "--trusted-host", "pypi.tuna.tsinghua.edu.cn"])
                        pip_cmd.append(pkg)
                        if pkg == "haddock3":
                            pip_cmd.extend(["--only-binary", ":all:"])
                        try:
                            rc = self._run_cmd_stream(pip_cmd, timeout=300)
                            if rc == 0:
                                self._log(f" {pkg} installed")
                            else:
                                self._log(f" {pkg} failed (exit code {rc})")
                                failed_pkgs.append(pkg)
                        except subprocess.TimeoutExpired:
                            self._log(f" {pkg} timed out")
                            failed_pkgs.append(pkg)
                        except Exception as e:
                            self._log(f" {pkg} error: {e}")
                            failed_pkgs.append(pkg)

                        installed_count += 1
                        pct = progress_start + (progress_end - progress_start) * installed_count / total_pkgs
                        self.progress["value"] = pct

                # Summary
                if failed_pkgs:
                    self._log(f"\n {len(failed_pkgs)} package(s) failed: {', '.join(failed_pkgs)}")
                else:
                    self._log("\n All packages installed successfully")

                # 3. Copy plugin files
                self._log("\n" + "=" * 50)
                self._log("[3/5] Installing plugin files...")
                self._log("=" * 50)

                install_path = self.install_path.get()
                source_dir = get_glint_source_dir()

                if source_dir and os.path.isdir(source_dir):
                    # Create destination directory
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
                            except:
                                pass

                    # Copy files
                    IGNORED = {'__pycache__', '.DS_Store', '.git', '.gitignore', '*.pyc'}
                    for item in os.listdir(source_dir):
                        if item in IGNORED:
                            continue
                        src = os.path.join(source_dir, item)
                        dst = os.path.join(install_path, item)
                        try:
                            if os.path.isdir(src):
                                shutil.copytree(src, dst, ignore=shutil.ignore_patterns(*IGNORED))
                            else:
                                shutil.copy2(src, dst)
                        except Exception as e:
                            self._log(f"Failed to copy {item}: {e}")

                    self._log(f"Copied to {install_path}")
                else:
                    self._log(f"Source directory not found: {source_dir}")

                self.progress["value"] = 80



            self.root.after(0, lambda: messagebox.showinfo(
                "Success",
                "GLINT installed successfully!\n\n"
                "You can now launch GLINT from the Desktop shortcut."))

        except Exception as e:
            self._log(f"\n Error: {e}")
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

        pymol_exe = os.path.join(self.env_path, "Scripts", "pymol.exe")
        if not os.path.exists(pymol_exe):
            self._log(f"PyMOL executable not found: {pymol_exe}")
            self._log("GLINT will not open until PyMOL is installed correctly")
            return False





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
                self._log(f"Created shortcut: {shortcut_path}")
                if icon_path:
                    self._log(f"Icon set: {icon_path}")
            else:
                raise Exception(result.stderr.decode() if result.stderr else "Unknown error")
        except Exception as e:
            self._log(f"Shortcut creation failed: {e}")
            self._log(f"You can manually run: {launcher_bat}")

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
                self._log("Logo PNG not found, using default icon")
                return None

            icon_path = os.path.join(launcher_dir, "glint.ico")

            # Try to use PIL to convert PNG to ICO
            try:
                from PIL import Image
                img = Image.open(logo_png)
                # Create multiple sizes for better quality
                sizes = [(16, 16), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)]
                img.save(icon_path, format='ICO', sizes=sizes)
                self._log(f"Created icon: {icon_path}")
                return icon_path
            except ImportError:
                self._log("PIL not available, trying alternative method...")

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
                    self._log(f"Created icon: {icon_path}")
                    return icon_path
            except Exception as e:
                self._log(f"Icon conversion failed: {e}")

            return None

        except Exception as e:
            self._log(f"Icon creation error: {e}")
            return None


def main():
    try:
        # Set DPI awareness (Windows 10+)
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
