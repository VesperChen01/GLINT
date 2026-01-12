# -*- mode: python ; coding: utf-8 -*-

block_cipher = None

# 获取当前脚本的目录
import os
import sys
script_dir = os.path.dirname(os.path.abspath(__file__))

# 假设 glint 源码在项目根目录下的 glint 文件夹中
# PyInstaller 运行时，sys._MEIPASS 会指向临时解压目录
# 所以这里需要确保 PyInstaller 能够找到 glint 模块
glint_source_path = os.path.abspath(os.path.join(script_dir, "..", "glint"))

a = Analysis(
    ['GLINT_Installer.py'],
    pathex=[script_dir],
    binaries=[],
    datas=[
        (os.path.join(glint_source_path, "**", "*"), "glint"), # 递归添加整个 glint 目录
    ],
    hiddenimports=[
        'tkinter',
        'tkinter.ttk',
        'tkinter.filedialog',
        'tkinter.messagebox',
        'PIL', # Pillow 库
        'PIL.Image',
        'PIL.ImageTk',
        'webbrowser',
        'subprocess',
        'shutil',
        'threading',
        'pathlib',
        'os',
        'sys',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='GLINT_Installer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_info_entries=[('PyInstaller', 'PyInstaller')],
    console=False, # GUI 应用，不需要控制台
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    # Linux 不使用 .ico 文件，图标在 .desktop 文件中指定 PNG
)