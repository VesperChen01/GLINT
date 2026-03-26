# -*- mode: python ; coding: utf-8 -*-
"""
PyInstaller spec file for GLINT Windows Installer
Creates a single .exe file with embedded glint source code
"""

import os

spec_dir = os.path.dirname(os.path.abspath(SPEC))
glint_src = os.path.join(spec_dir, 'glint')
external_src = os.path.join(spec_dir, 'external')

glint_datas = []
if os.path.exists(glint_src):
    for root, dirs, files in os.walk(glint_src):
        dirs[:] = [d for d in dirs if d != '__pycache__']
        for file in files:
            if file.endswith('.pyc') or file == '.DS_Store':
                continue
            glint_datas.append((os.path.join(root, file), os.path.relpath(root, spec_dir)))

external_datas = []
if os.path.exists(external_src):
    for root, dirs, files in os.walk(external_src):
        dirs[:] = [d for d in dirs if d != '__pycache__']
        for file in files:
            if file.endswith('.pyc') or file == '.DS_Store':
                continue
            external_datas.append((os.path.join(root, file), os.path.relpath(root, spec_dir)))

a = Analysis(
    ['GLINT_Installer.py'],
    pathex=[spec_dir],
    binaries=[],
    datas=glint_datas + external_datas,
    hiddenimports=['tkinter', 'tkinter.ttk', 'tkinter.filedialog', 'tkinter.messagebox'],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data)

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
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='glint/assets/logo.ico',
)
