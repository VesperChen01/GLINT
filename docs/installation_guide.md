# GLINT Installation Guide

## Chapter 1: Installation

Welcome to the official GLINT installation guide. This chapter covers the system requirements, pre-requisites, and step-by-step installation procedures for macOS, Windows, and Linux.

---

## 1. System Requirements

| Component | Requirement |
|-----------|-------------|
| **Operating System** | macOS 10.14+, Windows 10/11, Linux (Ubuntu 20.04+) |
| **Processor** | 64-bit processor |
| **RAM** | 8 GB minimum (16 GB recommended) |
| **Disk Space** | 5 GB free space |
| **Software** | Miniconda or Anaconda |

---

## 2. Pre-requisites: Installing Conda

GLINT relies on **Conda** for dependency management.

1. Download Miniconda from the official website:
   [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
2. Run the installer and follow the on-screen prompts.
3. Restart your terminal to ensure Conda is in your PATH.

---

## 3. Installation Instructions

### 3.1 macOS Installation

**Method 1: GUI Installer (Recommended)**

1. Download the latest `GLINT_Installer_X.X.X.dmg` from the [GitHub Releases](https://github.com/VesperChen01/GLINT/releases) page.

    ![Download DMG file from GitHub Releases](./images/ch1/macos_01_github_release.png)
    *Figure 1: Download the latest macOS DMG installer from the GitHub releases page.*

2. Double-click the `.dmg` file to mount the disk image.

    ![Finder window showing DMG mounted and GLINT Installer icon visible](./images/ch1/macos_02_dmg_mounted.png)
    *Figure 2: DMG file opened in Finder showing GLINT Installer.app icon. This window confirms the image has been mounted successfully.*

3. Drag the **GLINT Installer** icon into your Applications folder.

    ![Drag and drop window showing GLINT Installer being dragged to Applications folder](./images/ch1/macos_03_drag_to_apps.png)
    *Figure 3: Standard macOS drag-and-drop window. Drag the installer icon from the left to the Applications folder on the right.*

4. Launch **GLINT Installer** from your Applications folder.
5. The installer will perform a **System Check**.

    ![Main window of GLINT Installer on macOS showing Environment Status and Installation options](./images/ch1/macos_04_installer_main.png)
    *Figure 4: GLINT Installer main window showing welcome screen and installation options. Verify that both "Conda" and "Environment 'glint'" status indicators show green checkmarks before proceeding.*

6. Click **Install** to begin.

    ![Progress bar filling up in the GLINT Installer window with conda setup logs visible](./images/ch1/macos_05_installation_progress.png)
    *Figure 5: Installation progress bar with conda environment setup. The progress bar indicates the installation status, and the log text area shows the conda environment creation steps.*

**Method 2: One-Shot Script**

```bash
bash install_glint.sh
```

---

### 3.2 Windows Installation

1. Download the latest `GLINT_Installer_X.X.X.exe` from the [GitHub Releases](https://github.com/VesperChen01/GLINT/releases) page.

    ![GitHub releases page highlighting the Windows .exe file](./images/ch1/windows_01_github_release.png)
    *Figure 6: GitHub Releases page highlighting the Windows .exe installer.*

2. Double-click the `.exe` file to run the installer. You may see a **User Account Control (UAC)** prompt.

    ![Windows UAC prompt asking for permission to run the installer](./images/ch1/windows_02_uac_prompt.png)
    *Figure 7: UAC prompt for administrator permission. Click "Yes" to allow the installer to make changes to your device.*

3. The **GLINT Installer** window will open.

    ![Main window of the Windows GLINT Installer showing Installation Path and options](./images/ch1/windows_03_installer_main.png)
    *Figure 8: Windows installer GUI main window. Note the status indicators, installation path field (defaulting to `C:\Users\<Username>\.pymol\startup\glint`), and the main action buttons.*

4. Click **Install GLINT** to proceed.

    ![Windows Installer showing progress bar and installation logs in real-time](./images/ch1/windows_04_installation_progress.png)
    *Figure 9: Windows installation progress. The progress bar is filling, and the log text area shows the real-time installation steps being performed.*

---

### 3.3 Linux Installation

1. Ensure **Miniconda** or **Anaconda** is installed.
2. Make the installer script executable and run it:

```bash
chmod +x install_glint.sh
./install_glint.sh
```

    ![Linux terminal showing the successful execution of the GLINT installation script](./images/ch1/linux_01_terminal_success.png)
    *Figure 10: Linux Terminal - Running Installation Script. Shows green checkmarks and progress indicators as the script executes successfully.*

---

## 4. Installation Paths

| Platform | Path |
|----------|------|
| **macOS / Linux** | `~/.pymol/startup/glint` |
| **Windows** | `C:\Users\<Username>\.pymol\startup\glint` |

---

## 5. Optional Dependencies

### 5.1 HADDOCK3 (Protein-Protein Docking)

Required for ternary complex modeling.
*   **Automatic:** Installed by the installer/script if build tools are available.
*   **Manual:** `pip install haddock3`

### 5.2 FoldX (Mutation Analysis)

Used for mutational binding energy prediction.
*   **macOS/Linux:** Requires a commercial license.
*   **Windows:** Available via official installation packages.

### 5.3 APBS (Electrostatic Analysis)

Used for surface electrostatic calculations.
*   **macOS:** Bundled with the GUI installer. Otherwise: `brew install brewsci/bio/apbs`
*   **Windows:** Download from [poissonboltzmann.org](https://www.poissonboltzmann.org/).

    ![GLINT Settings or Environment Checker interface showing APBS configuration](./images/ch1/apbs_configuration.png)
    *Figure 11: APBS Configuration (Optional). If APBS is not found automatically, users can manually specify the path in the GLINT Environment Checker or Settings menu.*

---

## 6. Verifying Installation

1. **Launch PyMOL**: Run `pymol` in your terminal or click the PyMOL icon.
2. **Check Menu**: Look for a **"GLINT"** menu in the PyMOL menu bar.

    ![PyMOL main interface with the GLINT menu visible in the top menu bar](./images/ch1/pymol_01_glint_menu.png)
    *Figure 12: GLINT menu visible in PyMOL GUI. The GLINT menu appears between the "Plugin" and "Wizard" menus in the top menu bar.*

3. **Run Environment Checker**: Click **GLINT -> Utils -> Environment Checker**.

    ![PyMOL GUI showing the GLINT Environment Checker window with all checks passing](./images/ch1/pymol_02_env_checker.png)
    *Figure 13: PyMOL command line showing successful glint import and Environment Checker displaying all dependencies with green status indicators.*

---

## 7. Troubleshooting

### Conda not found

*   Ensure Conda is in your system PATH. Restart your terminal or run `source ~/.bash_profile` (macOS/Linux).

### GLINT menu not appearing

*   Check if the startup script exists: `ls ~/.pymol/startup/01_glint.py`.
*   Check the PyMOL console for errors when launching.

### Missing dependencies

*   Run the installer again or manually install missing packages via `conda install <package_name>`.

---

## Next Steps

Once GLINT is installed, proceed to **Chapter 2: Quick Start** to learn how to load structures and analyze molecular glue complexes.

---

*Documentation Version: v0.3.0*

