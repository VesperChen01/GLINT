# GLINT macOS Installer - Quick Start Guide

## For End Users

### Installation Steps

1. **Download the installer**
   - Get `GLINT_Installer_v0.1.28-beta.dmg`

2. **Open the DMG**
   - Double-click the downloaded DMG file
   - A new window will open showing "GLINT Installer.app"

3. **Run the installer**
   - Double-click "GLINT Installer.app"
   - If macOS shows a security warning, go to System Preferences → Security & Privacy and click "Open Anyway"

4. **Follow the installation wizard**
   - The installer will check for conda/miniconda
   - If conda is not found, it will prompt you to install Miniconda first
   - Choose installation options:
     - ✅ Install/Update dependencies (recommended)
     - ✅ Create application shortcut (recommended)
   - Click "Install GLINT"

5. **Wait for installation**
   - The installer will:
     - Create a conda environment named "glint"
     - Install PyMOL and all dependencies
     - Copy GLINT files to ~/.pymol/startup/glint
     - Create GLINT.app in ~/Applications

6. **Launch GLINT**
   - Open GLINT from ~/Applications/GLINT.app
   - Or launch PyMOL and run: `glint_gui`

### Requirements

- macOS 10.13 or later
- Miniconda or Anaconda (will prompt to install if not found)
- ~2GB free disk space for dependencies

### Troubleshooting

**"Cannot open because it is from an unidentified developer"**
- Go to System Preferences → Security & Privacy
- Click "Open Anyway" next to the GLINT Installer message

**"Conda not found"**
- Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html
- Restart the installer after installing Miniconda

**Installation fails**
- Check that you have internet connection
- Ensure you have enough disk space (~2GB)
- Try running the installer again

---

## For Developers

### Building the Installer

```bash
cd macos_installer
./build_release.sh
```

This will:
1. Build the .app bundle
2. Create the DMG
3. Place it in the parent directory

### Verifying the Installer

```bash
cd macos_installer
./verify_installer.sh
```

### Individual Build Steps

**Build only the .app:**
```bash
cd macos_installer
./build_app.sh
```

**Create only the DMG:**
```bash
cd macos_installer
./create_dmg.sh
```

### File Structure

```
macos_installer/
├── GLINT_Installer.py      # Main installer GUI
├── build_app.sh            # Build .app bundle
├── create_dmg.sh           # Create DMG
├── build_release.sh        # Complete build pipeline
├── verify_installer.sh     # Verification script
├── README.md               # Documentation
├── FIX_SUMMARY.md          # Fix details
└── COMPLETE_FIX_REPORT.md  # Comprehensive report
```

### Customization

**Change version:**
Edit `create_dmg.sh` and update the version detection logic.

**Change install path:**
Edit `GLINT_Installer.py` and modify `DEFAULT_INSTALL_PATH`.

**Change dependencies:**
Edit `GLINT_Installer.py` and modify `CONDA_PACKAGES` or `PIP_PACKAGES`.

### Distribution

The DMG file can be distributed via:
- Direct download
- GitHub releases
- Website hosting
- Cloud storage

Users simply download and run - no additional setup required.

