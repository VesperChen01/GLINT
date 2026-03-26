# build_exe.ps1
# PowerShell script to build GLINT Windows Installer
# Run: powershell -ExecutionPolicy Bypass -File build_exe.ps1

Write-Host "========================================"
Write-Host "  GLINT Windows EXE Builder"
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

$GlintDir = Join-Path $ScriptDir "glint"
$ParentGlintDir = Join-Path (Split-Path -Parent $ScriptDir) "glint"
$ExternalDir = Join-Path $ScriptDir "external"
$ParentExternalDir = Join-Path (Split-Path -Parent $ScriptDir) "external"

# Ensure glint source is present
if (-not (Test-Path $GlintDir)) {
    Write-Host "[INFO] glint folder not found in current directory"
    if (Test-Path $ParentGlintDir) {
        Write-Host "[INFO] Copying glint from parent directory..."
        Copy-Item -Path $ParentGlintDir -Destination $GlintDir -Recurse -Force
        Write-Host "[OK] glint copied successfully"
    } else {
        Write-Host "[ERROR] glint folder not found!"
        Write-Host "Please copy the glint folder into this directory or keep it in the parent folder."
        Read-Host "Press Enter to exit"
        exit 1
    }
} else {
    Write-Host "[OK] glint folder found"
}

# Ensure external resources are present (optional, but needed for APBS/Vina)
if (-not (Test-Path $ExternalDir)) {
    Write-Host "[INFO] external folder not found in current directory"
    if (Test-Path $ParentExternalDir) {
        Write-Host "[INFO] Copying external from parent directory..."
        Copy-Item -Path $ParentExternalDir -Destination $ExternalDir -Recurse -Force
        Write-Host "[OK] external copied successfully"
    } else {
        Write-Host "[WARNING] external folder not found; APBS/Vina resources will not be bundled"
    }
} else {
    Write-Host "[OK] external folder found"
}

# Remove __pycache__ directories
Write-Host "[INFO] Cleaning __pycache__ directories..."
Get-ChildItem -Path $GlintDir -Recurse -Directory -Filter "__pycache__" | Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
if (Test-Path $ExternalDir) {
    Get-ChildItem -Path $ExternalDir -Recurse -Directory -Filter "__pycache__" | Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
}
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
    --name "GLINT_Installer" `
    --add-data "glint;glint" `
    --add-data "glint\assets;glint\assets" `
    --add-data "external;external" `
    --hidden-import tkinter `
    --hidden-import tkinter.ttk `
    --hidden-import tkinter.filedialog `
    --hidden-import tkinter.messagebox `
    GLINT_Installer.py

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
Write-Host "Output: $ScriptDir\dist\GLINT_Installer.exe"
Write-Host ""
Write-Host "You can distribute this EXE file to Windows users."
Write-Host ""

Read-Host "Press Enter to exit"