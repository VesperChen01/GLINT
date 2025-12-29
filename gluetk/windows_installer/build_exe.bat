@echo off
chcp 65001 >nul 2>&1
REM build_exe.bat
REM Build GlueTK Installer as .exe using PyInstaller
REM Must be run on Windows

setlocal enabledelayedexpansion

echo ========================================
echo   GlueTK Windows EXE Builder
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

REM Check if gluetk folder exists
if not exist "%SCRIPT_DIR%\gluetk" (
    echo [INFO] gluetk folder not found in current directory
    echo [INFO] Checking parent directory...
    
    if exist "%SCRIPT_DIR%\..\gluetk" (
        echo [INFO] Copying gluetk from parent directory...
        xcopy /E /I /Y /Q "%SCRIPT_DIR%\..\gluetk" "%SCRIPT_DIR%\gluetk"
        if errorlevel 1 (
            echo [ERROR] Failed to copy gluetk source
            echo.
            echo Please manually copy the gluetk folder into this directory:
            echo   %SCRIPT_DIR%
            echo.
            pause
            exit /b 1
        )
        echo [OK] gluetk copied successfully
    ) else (
        echo [ERROR] gluetk folder not found!
        echo.
        echo Please copy the gluetk folder into this directory:
        echo   %SCRIPT_DIR%
        echo.
        echo Expected structure:
        echo   windows_installer\
        echo   +-- build_exe.bat
        echo   +-- GlueTK_Installer.py
        echo   +-- gluetk\
        echo       +-- __init__.py
        echo       +-- gui\
        echo       +-- assets\
        echo       +-- ...
        echo.
        pause
        exit /b 1
    )
) else (
    echo [OK] gluetk folder found
)

REM Remove __pycache__ directories
echo [INFO] Cleaning __pycache__ directories...
for /d /r "%SCRIPT_DIR%\gluetk" %%d in (__pycache__) do (
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
    --name "GlueTK_Installer" ^
    --add-data "gluetk;gluetk" ^
    --hidden-import tkinter ^
    --hidden-import tkinter.ttk ^
    --hidden-import tkinter.filedialog ^
    --hidden-import tkinter.messagebox ^
    GlueTK_Installer.py

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
echo Output: %SCRIPT_DIR%\dist\GlueTK_Installer.exe
echo.
echo You can distribute this EXE file to Windows users.
echo.

pause