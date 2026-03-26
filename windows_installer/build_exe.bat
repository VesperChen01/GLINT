@echo off
chcp 65001 >nul 2>&1
REM build_exe.bat
REM Build GLINT Installer as .exe using PyInstaller
REM Must be run on Windows

setlocal enabledelayedexpansion

echo ========================================
echo   GLINT Windows EXE Builder
echo ========================================
echo.

set "SCRIPT_DIR=%~dp0"
set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"
set "PARENT_DIR=%SCRIPT_DIR%\.."
set "GLINT_SRC=%PARENT_DIR%\glint"
set "EXTERNAL_SRC=%PARENT_DIR%\external"
set "GLINT_DST=%SCRIPT_DIR%\glint"
set "EXTERNAL_DST=%SCRIPT_DIR%\external"

REM Check Python
echo [1/4] Checking Python...
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found!
    echo Please install Python 3.9+ from https://www.python.org/
    pause
    exit /b 1
)
echo [OK] Python found
echo.

REM Install PyInstaller
echo [2/4] Installing PyInstaller...
python -m pip install pyinstaller --quiet --disable-pip-version-check
if errorlevel 1 (
    echo [ERROR] Failed to install PyInstaller
    pause
    exit /b 1
)
echo [OK] PyInstaller ready
echo.

REM Prepare files
echo [3/4] Preparing files...

if not exist "%GLINT_SRC%" (
    echo [ERROR] Source glint folder not found: %GLINT_SRC%
    pause
    exit /b 1
)

if exist "%GLINT_DST%" (
    echo [INFO] Removing stale bundled glint directory...
    rmdir /S /Q "%GLINT_DST%"
)
echo [INFO] Copying latest glint source...
xcopy /E /I /Y /Q "%GLINT_SRC%" "%GLINT_DST%"
if errorlevel 1 (
    echo [ERROR] Failed to copy glint source
    pause
    exit /b 1
)

if exist "%EXTERNAL_DST%" (
    echo [INFO] Removing stale bundled external directory...
    rmdir /S /Q "%EXTERNAL_DST%"
)
if exist "%EXTERNAL_SRC%" (
    echo [INFO] Copying latest external resources...
    xcopy /E /I /Y /Q "%EXTERNAL_SRC%" "%EXTERNAL_DST%"
    if errorlevel 1 (
        echo [WARNING] Failed to copy external resources
    ) else (
        echo [OK] external copied successfully
    )
) else (
    echo [WARNING] external folder not found; APBS/Vina resources will not be bundled
)

echo [INFO] Cleaning old build artifacts...
if exist "%SCRIPT_DIR%\build" rmdir /S /Q "%SCRIPT_DIR%\build"
if exist "%SCRIPT_DIR%\dist" rmdir /S /Q "%SCRIPT_DIR%\dist"
for /d /r "%GLINT_DST%" %%d in (__pycache__) do @if exist "%%d" rmdir /S /Q "%%d"
for /d /r "%EXTERNAL_DST%" %%d in (__pycache__) do @if exist "%%d" rmdir /S /Q "%%d"
echo [OK] Files prepared
echo.

REM Run PyInstaller
echo [4/4] Building EXE with PyInstaller...
echo This may take a few minutes...
echo.

cd /d "%SCRIPT_DIR%"
pyinstaller --clean --noconfirm GLINT_Installer.spec
if errorlevel 1 (
    echo.
    echo [ERROR] PyInstaller build failed!
    echo.
    echo Common solutions:
    echo   1. Make sure Python is in PATH
    echo   2. Try running as Administrator
    echo   3. Check if antivirus is blocking
    echo.
    pause
    exit /b 1
)

echo.
echo ========================================
echo   Build Complete!
echo ========================================
echo.
echo Output: %SCRIPT_DIR%\dist\GLINT_Installer.exe
echo.
echo You can distribute this EXE file to Windows users.
echo.
pause
