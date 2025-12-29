# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller spec file for GlueTK Windows Installer
Creates a single .exe file with embedded gluetk source code
"""

import os
import sys

block_cipher = None

# Get the directory containing this spec file
spec_dir = os.path.dirname(os.path.abspath(SPEC))

# Path to gluetk source
gluetk_src = os.path.join(spec_dir, 'gluetk')

# Collect all gluetk files
gluetk_datas = []
if os.path.exists(gluetk_src):
    for root, dirs, files in os.walk(gluetk_src):
        # Skip __pycache__ directories
        dirs[:] = [d for d in dirs if d != '__pycache__']
        for file in files:
            if file.endswith('.pyc') or file == '.DS_Store':
                continue
            src_path = os.path.join(root, file)
            # Calculate relative path from gluetk_src
            rel_dir = os.path.relpath(root, spec_dir)
            gluetk_datas.append((src_path, rel_dir))

a = Analysis(
    ['GlueTK_Installer.py'],
    pathex=[spec_dir],
    binaries=[],
    datas=gluetk_datas,
    hiddenimports=[
        'tkinter',
        'tkinter.ttk',
        'tkinter.filedialog',
        'tkinter.messagebox',
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
    name='GlueTK_Installer',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,  # No console window
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,  # Add icon path here if available: icon='path/to/icon.ico'
)