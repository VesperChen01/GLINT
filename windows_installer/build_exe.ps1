# build_exe.ps1
# PowerShell script to build GLINT Windows Installer
# Run: powershell -ExecutionPolicy Bypass -File build_exe.ps1

Write-Host "========================================"
Write-Host "  GLINT Windows EXE Builder"
Write-Host "========================================"
Write-Host ""

$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ParentDir = Split-Path -Parent $ScriptDir
$GlintSrc = Join-Path $ParentDir "glint"
$ExternalSrc = Join-Path $ParentDir "external"
$GlintDst = Join-Path $ScriptDir "glint"
$ExternalDst = Join-Path $ScriptDir "external"

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
if (-not (Test-Path $GlintSrc)) {
    Write-Host "[ERROR] Source glint folder not found: $GlintSrc"
    Read-Host "Press Enter to exit"
    exit 1
}

if (Test-Path $GlintDst) {
    Write-Host "[INFO] Removing stale bundled glint directory..."
    Remove-Item -Path $GlintDst -Recurse -Force
}
Write-Host "[INFO] Copying latest glint source..."
Copy-Item -Path $GlintSrc -Destination $GlintDst -Recurse -Force

if (Test-Path $ExternalDst) {
    Write-Host "[INFO] Removing stale bundled external directory..."
    Remove-Item -Path $ExternalDst -Recurse -Force
}
if (Test-Path $ExternalSrc) {
    Write-Host "[INFO] Copying latest external resources..."
    Copy-Item -Path $ExternalSrc -Destination $ExternalDst -Recurse -Force
    Write-Host "[OK] external copied successfully"
} else {
    Write-Host "[WARNING] external folder not found; APBS/Vina resources will not be bundled"
}

Write-Host "[INFO] Cleaning old build artifacts..."
if (Test-Path (Join-Path $ScriptDir "build")) { Remove-Item -Path (Join-Path $ScriptDir "build") -Recurse -Force }
if (Test-Path (Join-Path $ScriptDir "dist")) { Remove-Item -Path (Join-Path $ScriptDir "dist") -Recurse -Force }
Get-ChildItem -Path $GlintDst -Recurse -Directory -Filter "__pycache__" -ErrorAction SilentlyContinue | Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
if (Test-Path $ExternalDst) {
    Get-ChildItem -Path $ExternalDst -Recurse -Directory -Filter "__pycache__" -ErrorAction SilentlyContinue | Remove-Item -Recurse -Force -ErrorAction SilentlyContinue
}
Write-Host "[OK] Files prepared"
Write-Host ""

# Run PyInstaller
Write-Host "[4/4] Building EXE with PyInstaller..."
Write-Host "This may take a few minutes..."
Write-Host ""

Set-Location $ScriptDir
pyinstaller --clean --noconfirm GLINT_Installer.spec

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