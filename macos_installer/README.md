# GLINT macOS Installer

This directory contains the macOS installer for GLINT.

## Building the Installer

### 1. Build the .app bundle

```bash
cd macos_installer
./build_app.sh
```

This will create `GLINT Installer.app` in the parent directory.

### 2. Create the DMG

```bash
./create_dmg.sh
```

This will create a DMG file in the parent directory.

## Structure

- `GLINT_Installer.py` - Main installer GUI script
- `build_app.sh` - Script to build the .app bundle
- `create_dmg.sh` - Script to create the DMG installer

## How it Works

1. The `.app` bundle contains:
   - `Contents/MacOS/launcher` - Bash script that finds conda and launches the Python installer
   - `Contents/Resources/GLINT_Installer.py` - The actual installer GUI
   - `Contents/Resources/glint/` - The GLINT source code to be installed
   - `Contents/Resources/AppIcon.icns` - Application icon

2. When the user opens the app:
   - The launcher script finds the conda installation
   - Launches the Python installer GUI
   - The GUI allows the user to:
     - Check conda environment
     - Install dependencies
     - Copy GLINT files to ~/.pymol/startup/glint
     - Create a GLINT.app in ~/Applications

## Requirements

- macOS 10.13 or later
- Miniconda or Anaconda installed

## Distribution

The DMG file can be distributed to users. They can:
1. Open the DMG
2. Double-click "GLINT Installer.app"
3. Follow the installation wizard

