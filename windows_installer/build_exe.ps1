# build_exe.ps1
# PowerShell script to build GlueTK Windows Installer
# Run: powershell -ExecutionPolicy Bypass -File build_exe.ps1

Write-Host "========================================"
Write-Host "  GlueTK Windows EXE Builder"
Write-Host "========================================"
Write-Host ""

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path

# Check Python
Write-Host "[1/4] Checking Python..."
try {
    $pythonVersion = python --version 2>&1
    Write-Host "[OK] $pythonVersion"
} catch {
    Write-Host "[ERROR] Python not found!"
    Write-Host "Please install Python 3.9+ from https://www.python.org/"
    Read-Host "Press Enter to exit"
    exit 1
}
Write-Host ""

# Install PyInstaller
Write-Host "[2/4] Installing PyInstaller..."
python -m pip install pyinstaller --quiet --disable-pip-version-check
if ($LASTEXITCODE -ne 0) {
    Write-Host "[ERROR] Failed to install PyInstaller"
    Read-Host "Press Enter to exit"
    exit 1
}
Write-Host "[OK] PyInstaller ready"
Write-Host ""

# Prepare files
Write-Host "[3/4] Preparing files..."

$GluetkDir = Join-Path $ScriptDir "gluetk"
$ParentGluetkDir = Join-Path (Split-Path -Parent $ScriptDir) "gluetk"

if (-not (Test-Path $GluetkDir)) {
    Write-Host "[INFO] gluetk folder not found in current directory"
    
    if (Test-Path $ParentGluetkDir) {
        Write-Host "[INFO] Copying gluetk from parent directory..."
        Copy-Item -Path $ParentGluetkDir -Destination $GluetkDir -Recurse -Force
        Write-Host "[OK] gluetk copied successfully"
    } else {
        Write-Host "[ERROR] gluetk folder not found!"
        Write-Host ""
        Write-Host "Please copy the gluetk folder into this directory:"
        Write-Host "  $ScriptDir"
        Write-Host ""
        Write-Host "Expected structure:"
        Write-Host "  windows_installer\"
        Write-Host "  +-- build_exe.ps1"
        Write-Host "  +-- GlueTK_Installer.py"
        Write-Host "  +-- gluetk\"
        Write-Host "      +-- __init__.py"
        Write-Host "      +-- gui\"
        Write-Host "      +-- assets\"
        Write-Host ""
        Read-Host "Press Enter to exit"
        exit 1
    }
} else {
    Write-Host "[OK] gluetk folder found"
}

# Remove __pycache__ directories
Write-Host "[INFO] Cleaning __pycache__ directories..."
Get-ChildItem -Path $GluetkDir -Recurse -Directory -Filter "__pycache__" | Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
Write-Host "[OK] Files prepared"
Write-Host ""

# Run PyInstaller
Write-Host "[4/4] Building EXE with PyInstaller..."
Write-Host "This may take a few minutes..."
Write-Host ""

Set-Location $ScriptDir

pyinstaller --noconfirm `
    --onefile `
    --windowed `
    --name "GlueTK_Installer" `
    --add-data "gluetk;gluetk" `
    --hidden-import tkinter `
    --hidden-import tkinter.ttk `
    --hidden-import tkinter.filedialog `
    --hidden-import tkinter.messagebox `
    GlueTK_Installer.py

if ($LASTEXITCODE -ne 0) {
    Write-Host ""
    Write-Host "[ERROR] PyInstaller build failed!"
    Write-Host ""
    Write-Host "Common solutions:"
    Write-Host "  1. Make sure Python is in PATH"
    Write-Host "  2. Try running as Administrator"
    Write-Host "  3. Check if antivirus is blocking"
    Write-Host ""
    Read-Host "Press Enter to exit"
    exit 1
}

Write-Host ""
Write-Host "========================================"
Write-Host "  Build Complete!"
Write-Host "========================================"
Write-Host ""
Write-Host "Output: $ScriptDir\dist\GlueTK_Installer.exe"
Write-Host ""
Write-Host "You can distribute this EXE file to Windows users."
Write-Host ""

Read-Host "Press Enter to exit"