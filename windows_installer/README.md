# GlueTK Windows Installer

This directory contains files to build a professional Windows installer for GlueTK.

## Files

- `GlueTK_Installer.py` - GUI installer application (similar to macOS version)
- `build_exe.bat` - Batch script to build the .exe installer
- `build_exe.ps1` - PowerShell script (alternative, more reliable)
- `GlueTK_Installer.spec` - PyInstaller specification file for advanced builds

## Quick Start

### Prerequisites

1. **Windows 10/11** - Must be run on Windows
2. **Python 3.9+** - Install from https://www.python.org/ (check "Add to PATH" during install)
3. **PyInstaller** - Will be installed automatically by the build script

### Build Steps

#### Step 1: Prepare the folder structure

Copy the `gluetk` source folder into `windows_installer/`:

```
windows_installer/
├── GlueTK_Installer.py
├── build_exe.bat
├── build_exe.ps1
├── GlueTK_Installer.spec
├── README.md
└── gluetk/          <-- Copy this from the project root
    ├── __init__.py
    ├── gui/
    ├── assets/
    └── ...
```

#### Step 2: Run the build script

**Option A: Use PowerShell (Recommended)**
```powershell
# Right-click on build_exe.ps1 and select "Run with PowerShell"
# Or run in PowerShell:
cd windows_installer
powershell -ExecutionPolicy Bypass -File build_exe.ps1
```

**Option B: Use Command Prompt**
```cmd
cd windows_installer
build_exe.bat
```

#### Step 3: Find the output

The installer will be created at:
```
windows_installer/dist/GlueTK_Installer.exe
```

### Manual Build (Alternative)

```cmd
cd windows_installer

REM Install PyInstaller
pip install pyinstaller

REM Build using spec file
pyinstaller GlueTK_Installer.spec

REM Or build directly
pyinstaller --onefile --windowed --name GlueTK_Installer --add-data "gluetk;gluetk" GlueTK_Installer.py
```

## Distribution

Distribute the single `GlueTK_Installer.exe` file to Windows users. They can:

1. Double-click the .exe to launch the installer
2. The installer will:
   - Check for Miniconda installation
   - Create a `gluetk` conda environment
   - Install all dependencies (PyMOL, RDKit, etc.)
   - Copy GlueTK plugin files
   - Create a desktop shortcut

## End User Requirements

- **Miniconda** or **Anaconda** - Required for dependency management
  - Download from: https://docs.conda.io/en/latest/miniconda.html
- **~5GB disk space** - For conda environment and dependencies

## Troubleshooting

### "Python not found"
- Install Python 3.9+ from https://www.python.org/
- **Important**: Check "Add Python to PATH" during installation
- Restart Command Prompt/PowerShell after installation

### Batch file shows garbled characters
- Use the PowerShell script instead: `build_exe.ps1`
- Or run: `chcp 65001` before running the batch file

### "gluetk not found"
Make sure you copied the `gluetk` folder into the `windows_installer` directory.

### "Conda not found" (for end users)
Install Miniconda from https://docs.conda.io/en/latest/miniconda.html

### Antivirus blocks the .exe
PyInstaller-generated executables may trigger false positives. You may need to:
- Add an exception in your antivirus software
- Sign the executable with a code signing certificate (for production)

### Build fails with permission error
- Try running as Administrator
- Make sure no antivirus is blocking the build process