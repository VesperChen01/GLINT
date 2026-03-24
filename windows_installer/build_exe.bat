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

REM Get current directory
set "SCRIPT_DIR=%~dp0"
set "SCRIPT_DIR=%SCRIPT_DIR:~0,-1%"

REM Check if glint folder exists
if not exist "%SCRIPT_DIR%\glint" (
    echo [INFO] glint folder not found in current directory
    echo [INFO] Checking parent directory...

    if exist "%SCRIPT_DIR%\..\glint" (
        echo [INFO] Copying glint from parent directory...
        xcopy /E /I /Y /Q "%SCRIPT_DIR%\..\glint" "%SCRIPT_DIR%\glint"
        if errorlevel 1 (
            echo [ERROR] Failed to copy glint source
            echo.
            echo Please manually copy the glint folder into this directory:
            echo   %SCRIPT_DIR%
            echo.
            pause
            exit /b 1
        )
        echo [OK] glint copied successfully
    ) else (
        echo [ERROR] glint folder not found!
        echo.
        echo Please copy the glint folder into this directory:
        echo   %SCRIPT_DIR%
        echo.
        echo Expected structure:
        echo   windows_installer\
        echo   +-- build_exe.bat
        echo   +-- GLINT_Installer.py
        echo   +-- glint\
        echo       +-- __init__.py
        echo       +-- gui\
        echo       +-- assets\
        echo       +-- ...
        echo.
        pause
        exit /b 1
    )
) else (
    echo [OK] glint folder found
)

REM Remove __pycache__ directories
echo [INFO] Cleaning __pycache__ directories...
for /d /r "%SCRIPT_DIR%\glint" %%d in (__pycache__) do (
    if exist "%%d" (
        rmdir /s /q "%%d" 2>nul
    )
)
echo [OK] Files prepared
echo.

REM Run PyInstaller
echo [4/4] Building EXE with PyInstaller...
echo This may take a few minutes...
echo.

cd /d "%SCRIPT_DIR%"

pyinstaller --noconfirm ^
    --onefile ^
    --windowed ^
    --name "GLINT_Installer" ^
    --add-data "glint;glint" ^
    --add-data "glint\assets;glint\assets" ^
    --hidden-import tkinter ^
    --hidden-import tkinter.ttk ^
    --hidden-import tkinter.filedialog ^
    --hidden-import tkinter.messagebox ^
    GLINT_Installer.py

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