"""
glint/path_utils.py — 跨平台路径工具模块

集中管理 Windows / macOS / Linux 三平台下的可执行文件查找逻辑，
消除各模块中分散的 _which()、find_*_executable() 等重复实现。

仅使用标准库：os, sys, shutil, platform, glob, pathlib
"""

import os
import sys
import shutil
import platform
import glob
from pathlib import Path
from typing import Optional, List, Sequence


# ==================== 平台判断 ====================

def is_windows() -> bool:
    """当前是否为 Windows 平台。"""
    return sys.platform == "win32"


def is_macos() -> bool:
    """当前是否为 macOS 平台。"""
    return sys.platform == "darwin"


def _exe_name(name: str) -> str:
    """Windows 下自动追加 .exe 后缀（若尚未包含）。"""
    if is_windows() and not name.lower().endswith(".exe"):
        return name + ".exe"
    return name


# ==================== Conda 路径工具 ====================

def get_conda_bin_dirs(prefix: str) -> List[str]:
    """
    给定 conda 环境前缀，返回该环境下的可执行文件目录列表。
    - Unix:    [prefix/bin]
    - Windows: [prefix/Scripts, prefix/Library/bin]
    """
    if is_windows():
        return [
            os.path.join(prefix, "Scripts"),
            os.path.join(prefix, "Library", "bin"),
        ]
    return [os.path.join(prefix, "bin")]


def get_conda_bin_dir(prefix: str) -> str:
    """
    给定 conda 环境前缀，返回主可执行文件目录。
    - Unix:    prefix/bin
    - Windows: prefix/Scripts
    """
    if is_windows():
        return os.path.join(prefix, "Scripts")
    return os.path.join(prefix, "bin")


def get_conda_search_roots() -> List[str]:
    """
    返回当前平台常见的 conda 安装根目录列表。
    用于在 CONDA_PREFIX 不可用时遍历查找。
    """
    home = os.path.expanduser("~")
    roots = []

    if is_windows():
        roots = [
            os.path.join(home, "miniconda3"),
            os.path.join(home, "miniforge3"),
            os.path.join(home, "anaconda3"),
            os.path.join(home, "Miniconda3"),
            os.path.join(home, "Anaconda3"),
            os.path.join(os.environ.get("LOCALAPPDATA", ""), "miniconda3"),
            os.path.join(os.environ.get("LOCALAPPDATA", ""), "Continuum", "anaconda3"),
            os.path.join("C:\\", "ProgramData", "miniconda3"),
            os.path.join("C:\\", "ProgramData", "Anaconda3"),
        ]
    elif is_macos():
        roots = [
            os.path.join(home, "miniconda3"),
            os.path.join(home, "miniforge3"),
            os.path.join(home, "anaconda3"),
            os.path.join(home, "opt", "miniconda3"),
            os.path.join("/opt", "miniconda3"),
            os.path.join("/opt", "anaconda3"),
            os.path.join("/opt", "homebrew", "Caskroom", "miniconda", "base"),
            os.path.join("/usr", "local", "Caskroom", "miniconda", "base"),
        ]
    else:  # Linux
        roots = [
            os.path.join(home, "miniconda3"),
            os.path.join(home, "miniforge3"),
            os.path.join(home, "anaconda3"),
            os.path.join(home, "opt", "miniconda3"),
            os.path.join("/opt", "miniconda3"),
            os.path.join("/opt", "anaconda3"),
        ]

    # 从 CONDA_EXE 推断根目录
    conda_exe = os.environ.get("CONDA_EXE", "")
    if conda_exe:
        inferred = os.path.dirname(os.path.dirname(conda_exe))
        if inferred and inferred not in roots:
            roots.insert(0, inferred)

    return [r for r in roots if r]  # 过滤空字符串


def get_conda_base_dirs() -> List[str]:
    """get_conda_search_roots 的别名，供 __init__.py 使用。"""
    return get_conda_search_roots()


def get_conda_site_packages_patterns(prefix: str) -> List[str]:
    """
    给定 conda 环境前缀，返回 site-packages glob 模式列表。
    - Unix:    prefix/lib/python*/site-packages, prefix/lib/site-packages
    - Windows: prefix/Lib/site-packages
    """
    if is_windows():
        return [os.path.join(prefix, "Lib", "site-packages")]
    return [
        os.path.join(prefix, "lib", "python*", "site-packages"),
        os.path.join(prefix, "lib", "site-packages"),
    ]




# ==================== 通用路径搜索 ====================

def get_common_bin_paths() -> List[str]:
    """
    返回平台相关的常见可执行文件目录列表（不含 conda）。
    """
    home = os.path.expanduser("~")

    if is_windows():
        program_files = os.environ.get("ProgramFiles", r"C:\Program Files")
        program_files_x86 = os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")
        return [
            program_files,
            program_files_x86,
            os.path.join(home, "AppData", "Local", "Programs"),
        ]
    elif is_macos():
        return [
            os.path.join(home, "bin"),
            os.path.join(home, ".local", "bin"),
            os.path.join(home, "local", "bin"),
            os.path.join("/usr", "local", "bin"),
            os.path.join("/opt", "homebrew", "bin"),
        ]
    else:  # Linux
        return [
            os.path.join(home, "bin"),
            os.path.join(home, ".local", "bin"),
            os.path.join(home, "local", "bin"),
            os.path.join("/usr", "local", "bin"),
            os.path.join("/usr", "bin"),
        ]


def _is_executable(path: str) -> bool:
    """判断路径是否为可执行文件。"""
    return os.path.isfile(path) and os.access(path, os.X_OK)



def find_in_conda_envs(
    exe_name: str,
    env_names: Sequence[str] = ("glint", "base"),
) -> Optional[str]:
    """
    遍历 conda 根目录和指定环境名，查找可执行文件。
    返回第一个找到的路径，或 None。
    """
    exe = _exe_name(exe_name)
    for root in get_conda_search_roots():
        if not os.path.isdir(root):
            continue
        for env in env_names:
            if env == "base":
                for bin_dir in get_conda_bin_dirs(root):
                    candidate = os.path.join(bin_dir, exe)
                    if _is_executable(candidate):
                        return candidate
            else:
                env_prefix = os.path.join(root, "envs", env)
                if os.path.isdir(env_prefix):
                    for bin_dir in get_conda_bin_dirs(env_prefix):
                        candidate = os.path.join(bin_dir, exe)
                        if _is_executable(candidate):
                            return candidate
    return None


def find_executable(
    exe_name: str,
    extra_names: Optional[Sequence[str]] = None,
    env_var: Optional[str] = None,
    search_conda_envs: bool = True,
    extra_paths: Optional[Sequence[str]] = None,
) -> Optional[str]:
    """
    终极可执行文件查找函数，组合以下策略（按优先级）：
    1. 环境变量 env_var 指定的路径
    2. shutil.which（系统 PATH，Windows 自动追加 .exe）
    3. 当前 CONDA_PREFIX 环境
    4. 通过 sys.executable 推断的 conda 环境
    5. 遍历常见 conda 安装根目录和环境
    6. 平台常见 bin 目录
    7. 用户提供的额外搜索路径
    """
    all_names = [exe_name]
    if extra_names:
        all_names.extend(extra_names)

    # 1. 环境变量
    if env_var:
        env_path = os.environ.get(env_var, "")
        if env_path and _is_executable(env_path):
            return env_path

    for name in all_names:
        exe = _exe_name(name)
        # 2. shutil.which
        found = shutil.which(exe)
        if found:
            return found
        if is_windows() and name != exe:
            found = shutil.which(name)
            if found:
                return found

    # 3. 当前 CONDA_PREFIX
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        for name in all_names:
            exe = _exe_name(name)
            for bin_dir in get_conda_bin_dirs(conda_prefix):
                candidate = os.path.join(bin_dir, exe)
                if _is_executable(candidate):
                    return candidate

    # 4. 从 sys.executable 推断 conda 环境 bin 目录
    try:
        python_dir = os.path.dirname(os.path.realpath(sys.executable))
        for name in all_names:
            exe = _exe_name(name)
            candidate = os.path.join(python_dir, exe)
            if _is_executable(candidate):
                return candidate
    except Exception:
        pass

    # 5. 遍历 conda 环境
    if search_conda_envs:
        for name in all_names:
            result = find_in_conda_envs(name)
            if result:
                return result

    # 6. 平台常见 bin 目录
    for bin_dir in get_common_bin_paths():
        for name in all_names:
            exe = _exe_name(name)
            candidate = os.path.join(bin_dir, exe)
            if _is_executable(candidate):
                return candidate

    # 7. 额外路径
    if extra_paths:
        for path in extra_paths:
            expanded = os.path.expanduser(path)
            if _is_executable(expanded):
                return expanded

    return None


# 向后兼容别名：各模块原 _which() 的直接替代
def which_tool(
    exe_name: str,
    extra_search_paths: Optional[Sequence[str]] = None,
    env_var: Optional[str] = None,
) -> Optional[str]:
    """which_tool — _which() 的统一替代。"""
    return find_executable(
        exe_name,
        env_var=env_var,
        extra_paths=extra_search_paths,
    )


# ==================== 平台搜索路径 ====================

def get_platform_search_paths() -> List[str]:
    """
    返回平台相关的通用搜索路径列表（含 conda + 常见 bin）。
    """
    paths = []
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        paths.extend(get_conda_bin_dirs(conda_prefix))
    try:
        python_dir = os.path.dirname(os.path.realpath(sys.executable))
        if python_dir not in paths:
            paths.append(python_dir)
    except Exception:
        pass
    paths.extend(get_common_bin_paths())
    return paths
