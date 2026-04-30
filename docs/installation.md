---
layout: default
title: Installation
---

# Installation Guide

Complete guide to installing GLINT on your system.

## System Requirements

| Component | Requirement |
|-----------|-------------|
| **Operating System** | macOS 10.14+, Windows 10/11, Linux (Ubuntu 20.04+) |
| **Processor** | 64-bit processor |
| **RAM** | 8 GB minimum (16 GB recommended) |
| **Disk Space** | 5 GB free space |
| **Software** | Miniconda or Anaconda |

## Prerequisites: Installing Conda

GLINT relies on **Conda** for dependency management.

1. Download Miniconda from the official website:
   [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
2. Run the installer and follow the on-screen prompts.
3. Restart your terminal to ensure Conda is in your PATH.

## Installation Methods

### macOS Installation

#### Method 1: GUI Installer (Recommended)

1. Download the latest `GLINT_Installer_X.X.X.dmg` from the [GitHub Releases](https://github.com/VesperChen01/GLINT/releases) page.
2. Double-click the `.dmg` file to mount the disk image.
3. Drag the **GLINT Installer** icon into your Applications folder.
4. Launch **GLINT Installer** from your Applications folder.
5. Click **Install** to begin the installation process.

#### Method 2: One-Shot Script

```bash
bash install_glint.sh
```

### Windows Installation

1. Download the latest `GLINT_Installer_X.X.X.exe` from the [GitHub Releases](https://github.com/VesperChen01/GLINT/releases) page.
2. Double-click the `.exe` file to run the installer.
3. Follow the on-screen instructions to complete the installation.

### Linux Installation

1. Ensure **Miniconda** or **Anaconda** is installed.
2. Make the installer script executable and run it:

```bash
chmod +x install_glint.sh
./install_glint.sh
```

## Installation Paths

| Platform | Path |
|----------|------|
| **macOS / Linux** | `~/.pymol/startup/glint` |
| **Windows** | `C:\Users\<Username>\.pymol\startup\glint` |

## Optional Dependencies

### HADDOCK3 (Protein-Protein Docking)

Required for ternary complex modeling.
- **Automatic:** Installed by the installer/script if build tools are available.
- **Manual:** `pip install haddock3`

### FoldX (Mutation Analysis)

Used for mutational binding energy prediction.
- **macOS/Linux:** Requires a commercial license.
- **Windows:** Available via official installation packages.

### APBS (Electrostatic Analysis)

Used for surface electrostatic calculations.
- **macOS:** `brew install brewsci/bio/apbs`
- **Windows:** Download from [poissonboltzmann.org](https://www.poissonboltzmann.org/).

## Verifying Installation

1. **Launch PyMOL**: Run `pymol` in your terminal or click the PyMOL icon.
2. **Check Menu**: Look for a **"GLINT"** menu in the PyMOL menu bar.
3. **Run Environment Checker**: Click **GLINT -> Utils -> Environment Checker**.

## Troubleshooting

### Conda not found

Ensure Conda is in your system PATH. Restart your terminal or run `source ~/.bash_profile` (macOS/Linux).

### GLINT menu not appearing

Check if the startup script exists: `ls ~/.pymol/startup/01_glint.py`

### Missing dependencies

Run the installer again or manually install missing packages via `conda install <package_name>`.
