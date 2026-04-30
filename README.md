<div align="center">

<img src="glint/assets/logo.png" alt="GLINT Logo" width="180"/>


### 🧬 PyMOL Plugin for Molecular Glue Discovery & Analysis

[![GitHub](https://img.shields.io/badge/GitHub-VesperChen01/GLINT-181717?style=flat-square&logo=github)](https://github.com/VesperChen01/GLINT)
[![License](https://img.shields.io/badge/License-MIT-blue?style=flat-square)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=flat-square&logo=python&logoColor=white)](https://python.org)
[![PyMOL](https://img.shields.io/badge/PyMOL-Plugin-green?style=flat-square)](https://pymol.org)

[Features](#-features) • [Installation](#-installation) • [Quick Start](#-quick-start) • [Tutorials](#-tutorials) • [Documentation](#-documentation)

---

</div>

## 📖 Overview

**GLINT** is a comprehensive PyMOL plugin designed for molecular glue discovery and protein-ligand interaction analysis. It provides specialized tools for analyzing ternary complexes, detecting G-motifs, identifying neo-epitopes, and structure-guided optimization workflows.

<details>
<summary><b>🔬 What are Molecular Glues?</b></summary>

Molecular glues are small molecules that induce or stabilize protein-protein interactions, typically between an E3 ubiquitin ligase and a target protein (neo-substrate), leading to targeted protein degradation.

</details>

---

## ✨ Features

<table>
<tr>
<td width="50%">

### 🎯 Target Discovery
- **G-Motif Detection** - CRBN G-loop binding motifs via RMSD matching
- **Protein Surface Analysis** - Electrostatic and hydrophobic patch analysis
- **Surface Similarity & Complementarity** - Surface-based comparison with APBS-supported electrostatics

</td>
<td width="50%">

### 🔍 Hit Identification
- **Vina Docking** - Integrated AutoDock Vina
- **HADDOCK3 Integration** - Ternary complex modeling
- **FoldX Mutation Analysis** - Mutational binding energy prediction

</td>
</tr>
<tr>
<td width="50%">

### 🧩 Ternary Complex Evaluation
- **Automated Ternary Complex Analysis** - One-click full evaluation with automatic ligand property calculation
- **Comprehensive Interface Analysis** - BSA calculation using ternary complex formula, contact statistics
- **Ligand Property Extraction** - Automatic SMILES extraction and molecular property computation (MW, LogP, TPSA, etc.)
- **Geometric Metrics** - COG shift, E3-MG-POI angle, distance measurements, duality index
- **Enhanced CSV Export** - Consistent with GUI display, includes timestamps and detailed results

</td>
<td width="50%">

### ⚡ Lead Optimization
- **Electrostatic Complementarity Analysis** - Ligand and ternary EC analysis workflows
- **PPI Interface Analysis** - Interaction type classification
- **Protein-Ligand Interactions** - H-bonds, salt bridges, π-π, etc.

</td>
</tr>
</table>

---

## 📦 Installation

<details open>
<summary><b>🍎 macOS Installer (Recommended)</b></summary>

1. Download `GLINT_Installer_0.3.0.dmg` from [Releases](https://github.com/VesperChen01/GLINT/releases)
2. Open the DMG and run the installer
3. Follow the on-screen instructions

</details>

<details>
<summary><b>🚀 One-Shot Setup Script</b></summary>

```bash
bash install_glint.sh
```

</details>

<details>
<summary><b>🪟 Windows Installation</b></summary>

1. Download `GLINT_Installer_X.X.X.exe` from [Releases](https://github.com/VesperChen01/GLINT/releases)
2. Double-click the `.exe` file and follow the installation wizard
3. Restart PyMOL to load the plugin

</details>

<details>
<summary><b>🐧 Linux Installation</b></summary>

```bash
chmod +x install_glint.sh
./install_glint.sh
```

</details>

---

## 🚀 Quick Start

1. **Launch PyMOL** and load your protein structure
2. **Access GLINT Menu** - Look for the "GLINT" menu in the menu bar
3. **Run Analysis** - Navigate to the desired analysis tool
4. **View Results** - Results are displayed in real-time with visualization options

### Example: Ternary Complex Evaluation

```python
# Load your ternary complex structure
load ternary_complex.pdb

# Access GLINT -> Ternary Complex Evaluation
# Set parameters:
# - E3 Chain: A
# - POI Chain: B  
# - Ligand Resn: UNL

# Click "Run Full Evaluation"
# GLINT will automatically:
# - Calculate interface BSA using ternary complex formula
# - Extract ligand SMILES and compute molecular properties
# - Analyze geometric metrics (COG shift, angles, distances)
# - Calculate duality index for binding balance

# Export results via "Export" button
```

---

## 📚 Documentation

- **[Installation Guide](docs/installation.md)** - Detailed installation instructions
- **[Tutorial](docs/tutorial.md)** - Step-by-step tutorials for all features
- **[API Reference](docs/api.md)** - Complete API documentation

---

## 🎯 Key Features in v0.3.0

### ✨ New: Automated Ternary Complex Analysis
- **One-click full evaluation** - Run complete analysis with single button click
- **Automatic ligand property calculation** - SMILES extraction and molecular computation
- **Enhanced BSA calculation** - Using ternary complex formula with intermediate values
- **Improved CSV export** - Consistent with GUI display, includes timestamps

### 🔧 Improvements
- Better error handling and user feedback
- Real-time progress updates
- Enhanced visualization options
- Improved documentation and tutorials

---

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- PyMOL community for the excellent molecular visualization platform
- RDKit and OpenBabel for chemical informatics tools
- BioPython for structural biology utilities
- All contributors and users of GLINT

---

<div align="center">

**[⬆ Back to Top](#-overview)**

Made with ❤️ by the GLINT Team

</div>
