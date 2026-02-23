<div align="center">

<img src="glint/assets/logo.png" alt="GLINT Logo" width="180"/>

# GLINT

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
- **Surface Analysis** - Electrostatic & hydrophobic patches
- **Surface Similarity** - MaSIF-style binding site comparison

</td>
<td width="50%">

### 🔍 Hit Identification
- **Pocket Detection** - Druggable pockets with scoring
- **Vina Docking** - Integrated AutoDock Vina
- **HADDOCK3 Integration** - Ternary complex modeling

</td>
</tr>
<tr>
<td width="50%">

### ⚡ Lead Optimization
- **PPI Interface Analysis** - Interaction type classification
- **Neo-Epitope Identification** - Glue-induced epitopes
- **Protein-Ligand Interactions** - H-bonds, salt bridges, π-π, etc.
- **Mutation Analysis** - ΔΔG prediction (FoldX/PyRosetta)

</td>
<td width="50%">

### 🎨 Visualization
- **Publication-Quality Rendering** - Customizable 3D styles
- **2D Interaction Diagrams** - RDKit-powered diagrams
- **Electrostatic Surface** - Coulomb & APBS electrostatics

</td>
</tr>
</table>

---

## 📦 Installation

<details open>
<summary><b>🍎 macOS Installer (Recommended)</b></summary>

1. Download `GLINT_Installer_v0.2.3.dmg` from [Releases](https://github.com/VesperChen01/GLINT/releases)
2. Open the DMG and run the installer
3. Follow the on-screen instructions

</details>

<details>
<summary><b>🚀 One-Shot Setup Script</b></summary>

```bash
# Clone or download GLINT
cd /path/to/glint

# Run the setup script
bash install_glint.sh
```

This creates a Conda environment `glint` with all dependencies.

</details>

<details>
<summary><b>🔧 Manual Installation</b></summary>

```bash
# Create conda environment
conda create -n glint python=3.9 -y
conda activate glint

# Install dependencies
conda install -c conda-forge rdkit scipy matplotlib pillow numpy pandas seaborn pyqt -y

# Optional: Install PyMOL
conda install -c conda-forge pymol-open-source -y

# Copy glint folder to PyMOL startup directory
cp -r glint ~/.pymol/startup/
```

</details>

<details>
<summary><b>📦 PyMOL Plugin Manager</b></summary>

1. Download `glint.zip` from releases
2. In PyMOL: `Plugin → Plugin Manager → Install New Plugin`
3. Select the ZIP file

</details>

---

## 🚀 Quick Start

### Launch GUI

```python
# In PyMOL command line
glint_gui
```

Or via menu: `Plugin → GLINT - Molecular Glue Analyzer`

### Command Line Examples

```python
# Load a structure
fetch 6H0G

# Detect G-motif
find_crbn_g_motif 6h0g, rmsd_cutoff=3.5

# Analyze PPI interface
ppi_analyze 6h0g, ['A'], ['B']

# Identify neo-epitope
neo_epitope_find 6h0g, ['A'], ['B'], CC9

# Generate 2D diagram
generate_2d_diagram 6h0g, CC9
```

---

## 🖥️ GUI Interface Guide

GLINT provides a modern, intuitive graphical interface with five main modules:

### Welcome Page

The welcome page provides quick access to all modules and features classic molecular glue case studies:

| Case | PDB | Description |
|------|-----|-------------|
| Thalidomide | 4CIW | CRBN degrader prototype |
| Indisulam | 5S9M | DCAF15-RBM39 splicing factor degrader |
| QS1 | 8P1D | Next-gen GSPT1 degrader |

Click any case to automatically load the structure into PyMOL.

---

### 🎯 Target Discovery Tab

Identify potential molecular glue targets and binding sites.

<details>
<summary><b>G-Motif (CRBN G-loop) Detection</b></summary>

Detect CRBN G-loop binding motifs using RMSD-based template matching.

| Parameter | Description | Default |
|-----------|-------------|---------|
| Target Object | PyMOL object to analyze | - |
| RMSD Cutoff | Maximum RMSD for match (Å) | 3.5 |
| Template | Reference structure (GSPT1/CK1α/VAV1) | GSPT1 |
| Require Central Gly | Filter for glycine at position 4 | ✓ |
| Highlight Surface | Show molecular surface of G-loop | ✓ |
| Export Coordinates | Save G-loop atom coordinates | ✗ |

**Buttons:**
- **Detect POI** - Run G-motif detection
- **Show Surface** - Highlight G-loop surface
- **Export Coords** - Save coordinates to CSV
- **Surface Analysis** - Analyze electrostatic/hydrophobic patches

</details>

<details>
<summary><b>Pocket Detection & Interface Pockets</b></summary>

Detect druggable pockets and prioritize pockets located at protein-protein interfaces.

| Parameter | Description | Default |
|-----------|-------------|---------|
| Detection Scope | Whole structure or selected region | Whole structure |
| Druggability Filter | Minimum pocket score threshold | 0.0 |
| Interface Distance | Contact cutoff for interface pockets (Å) | 5.0 |
| Output | Export pocket summary CSV | Optional |

**Buttons:**
- **Detect Pockets** - Run global pocket detection
- **Interface Pockets** - Focus on pockets at protein interfaces
- **Visualize Pockets** - Color-code by score and show labels

</details>

<details>
<summary><b>Protein Surface Analysis</b></summary>

Detect electrostatic and hydrophobic patches on protein surfaces.

**Buttons:**
- **Analyze Surface** - Run patch detection
- **Visualize Patches** - Color patches (Red=positive, Blue=negative, Green=hydrophobic)

</details>

<details>
<summary><b>Surface Similarity & Complementarity (MaSIF-style)</b></summary>

Compare binding site features using geometric and chemical descriptors.

| Parameter | Description | Default |
|-----------|-------------|---------|
| Analysis Type | Similarity or Complementarity (PPI) | Similarity |
| Surface Method | auto/msms/open3d/edtsurf | auto |
| Patch Radius | Radius for comparison (Å) | 12.0 |
| Interface Distance | Contact threshold (Å) | 4.0 |
| Use APBS | Accurate electrostatics (slower) | ✗ |

**Output Metrics:**
- Geometric complementarity
- Electrostatic complementarity
- Hydrophobic complementarity
- Interface area

</details>

---

### 🔍 Hit Identification Tab

Identify binding sites and perform molecular docking.

<details>
<summary><b>Binding Site Detection</b></summary>

Detect druggable binding pockets with volume, depth, and druggability scoring.

| Parameter | Description | Default |
|-----------|-------------|---------|
| Grid Spacing | Detection resolution (Å) | 0.6 |
| Min Volume | Minimum pocket volume (ų) | 20 |
| Color By | Volume/Druggability/Hydrophobicity/Depth | Volume |

</details>

<details>
<summary><b>AutoDock Vina Docking</b></summary>

Integrated molecular docking to detected pockets.

| Parameter | Description | Default |
|-----------|-------------|---------|
| Ligand File | MOL2/SDF/PDBQT file | - |
| Max Pockets | Number of pockets to dock | 3 |
| Exhaustiveness | Search thoroughness | 8 |
| Custom Box | Manual docking box coordinates | ✗ |

</details>

<details>
<summary><b>HADDOCK3 Protein-Protein Docking</b></summary>

Model ternary complexes with protein-protein docking.

| Mode | Description |
|------|-------------|
| Blind (Random AIR) | No restraints, random sampling |
| Blind (Centroid) | Center-of-mass guided |
| Blind (Surface) | Surface-based sampling |
| Pocket-constrained | Residue-based AIR restraints |

</details>

---

### ⚡ Lead Optimization Tab

Analyze and optimize molecular interactions.

<details>
<summary><b>Protein-Protein Interface (PPI) Analysis</b></summary>

Analyze protein-protein interfaces with interaction classification.

| Parameter | Description | Default |
|-----------|-------------|---------|
| Protein1 Chains | First protein chain IDs | - |
| Protein2 Chains | Second protein chain IDs | - |
| Interface Distance | Contact threshold (Å) | 4.5 |
| Display Mode | Cartoon+Surface or Surface only | Cartoon+Surface |
| Show Labels | Display distance labels | ✓ |

**Output:**
- Interface contacts count
- Buried Surface Area (BSA)
- Interface strength score (0-10)
- Strong/Weak classification

</details>

<details>
<summary><b>Protein-Ligand Interactions</b></summary>

Comprehensive interaction analysis with visualization.

| Parameter | Description | Default |
|-----------|-------------|---------|
| Ligand Name | Residue name (auto-detect if blank) | - |
| Protein Chains | Specific chains (optional) | - |
| Distance | Interaction cutoff (Å) | 4.5 |

**Interaction Types Detected:**
- Hydrogen bonds
- Salt bridges
- Hydrophobic contacts
- π-π stacking
- Cation-π interactions
- Halogen bonds

**Buttons:**
- **Analyze Protein-Ligand** - Run analysis with 3D visualization
- **Generate 2D Diagram** - RDKit-powered 2D interaction diagram

</details>

<details>
<summary><b>Ligand-Ligand Interactions</b></summary>

Analyze interactions between two small molecules.

| Parameter | Description |
|-----------|-------------|
| Selection 1 | PyMOL selection for first ligand |
| Selection 2 | PyMOL selection for second ligand |
| Distance | Interaction cutoff (Å) |

</details>

<details>
<summary><b>Mutation Analysis</b></summary>

Predict mutation effects on binding.

- **Perform Mutation** - Apply mutations via PyMOL wizard
- **Minimize** - Energy minimization with sculpting
- **Full Analysis** - ΔΔG prediction (requires FoldX/PyRosetta)

</details>

### 🎨 Visualization Tab

Generate publication-quality images.

<details>
<summary><b>Electrostatic Surface (APBS)</b></summary>

| Parameter | Description | Default |
|-----------|-------------|---------|
| Grid Spacing | Resolution (Å) | 1.0 |
| Range | Color scale (min,mid,max) | -5,0,5 |

**Buttons:**
- **Quick ESP** - Fast Coulomb-based calculation
- **True APBS** - Full APBS calculation (requires APBS tools)

</details>

<details>
<summary><b>Image Export</b></summary>

| Parameter | Description | Default |
|-----------|-------------|---------|
| Width | Image width (px) | 3000 |
| Height | Image height (px) | 2000 |
| DPI | Resolution | 300 |
| Background | White or Transparent | White |
| Ray Trace | High-quality rendering | ✓ |

**Buttons:**
- **Export PNG** - Save current view
- **Export DX** - Save electrostatic map
- **Fill Viewport** - Use current viewport size

</details>

---

## 📚 Tutorials

<details>
<summary><b>Tutorial 1: Analyzing a CRBN-Glue-Substrate Complex</b></summary>

Analyze the CRBN-CC885-GSPT1 complex (PDB: 6H0G).

#### Step 1: Load Structure
```python
fetch 6H0G
```

#### Step 2: Launch GUI
```python
glint_gui
```

#### Step 3: G-Motif Detection
1. Go to **Target Discovery** tab
2. Select `6h0g` from dropdown
3. Set RMSD cutoff: `3.5` Å
4. Select template: `GSPT1 (6H0G)`
5. Click **Detect POI**

```
✅ G-Motif detection complete: 1 hits found
   Chain B, Residues 568-575, Sequence: GSGKGSSF, RMSD: 0.00 Å
```

#### Step 4: PPI Interface Analysis
1. Go to **Lead Optimization** tab
2. Set Protein1 Chains: `A` (CRBN)
3. Set Protein2 Chains: `B` (GSPT1)
4. Click **Analyze PPI Interface**

```
PPI Analysis Complete:
  Interface Contacts: 45
  BSA: 1250.3 Å²
  Interface Strength: 7.5/10
```

#### Step 5: Visualization
Go to **Visualization** tab → Click **Render All**

</details>

<details>
<summary><b>Tutorial 3: Molecular Glue Case Comparison</b></summary>

| Metric | Case A (6H0G) | Case B (5T35) |
|--------|---------------|---------------|
| Interface contacts | Compare interaction counts | Compare interaction counts |
| Neo-epitope pattern | Evaluate detected residues | Evaluate detected residues |
| Pocket profile | Compare pocket size/score | Compare pocket size/score |

```python
# Case A analysis
fetch 6H0G
ppi_analyze 6h0g, ['A'], ['B']
# → Inspect interface contacts

# Case B analysis
fetch 5T35
ppi_analyze 5t35, ['A'], ['B']
# → Compare with Case A
```

</details>

<details>
<summary><b>Tutorial 4: Pocket-Based Virtual Screening</b></summary>

```python
# Load target
fetch 6H0G

# Detect pockets
detect_pockets 6h0g, min_volume=50, min_depth=3.0

# Run docking (requires Vina)
pocket_based_docking 6h0g, "ligand.mol2", max_pockets=3

# Analyze results
load docking_results/ligand_pocket1_out.pdbqt, docked_pose
analyze_protein_ligand_interactions 6h0g, docked_pose
```

</details>

<details>
<summary><b>Tutorial 5: Surface Complementarity Analysis</b></summary>

1. Go to **Target Discovery** → **Surface Similarity**
2. Set Selection 1: `chain A` (CRBN)
3. Set Selection 2: `chain B` (GSPT1)
4. Select: `Complementarity (PPI)`
5. Click **Compare**

```
✅ Complementarity analysis complete:
   Overall Score: 0.75
   Geometric: 0.82
   Electrostatic: 0.68
   Hydrophobic: 0.71
```

</details>

---

## 📋 Command Reference

<details>
<summary><b>Core Analysis Commands</b></summary>

| Command | Description |
|---------|-------------|
| `analyze_protein_ligand_interactions` | Protein-ligand interaction detection |
| `analyze_pdb_interactions` | Generic PDB interaction analysis |
| `analyze_ternary_complex` | Ternary complex analysis |

</details>

<details>
<summary><b>Molecular Glue Commands</b></summary>

| Command | Description |
|---------|-------------|
| `find_crbn_g_motif` | G-motif detection |
| `ppi_analyze` | PPI interface analysis |
| `neo_epitope_find` | Neo-epitope detection |
| `detect_pockets` | Druggable pocket detection |

</details>

<details>
<summary><b>Visualization Commands</b></summary>

| Command | Description |
|---------|-------------|
| `visualize_protein_ligand_3d` | 3D visualization |
| `generate_2d_diagram` | 2D interaction diagram |
| `visualize_pockets` | Pocket visualization |

</details>

<details>
<summary><b>Docking Commands</b></summary>

| Command | Description |
|---------|-------------|
| `pocket_based_docking` | Vina docking to detected pockets |
| `vina_score_complex` | Score existing complex |

</details>

---

## ⚙️ Requirements

| Requirement | Version | Note |
|-------------|---------|------|
| PyMOL | Open-Source or Incentive | Required |
| Python | 3.8+ | Required |
| PyQt | 5 or 6 | Required |
| NumPy, SciPy, Pandas | Latest | Required |
| RDKit | Latest | Optional (2D diagrams) |
| AutoDock Vina | Latest | Optional (docking) |

---

## 📝 Citation

If you use GLINT in your research, please cite:

```bibtex
@software{glint2024,
  title = {GLINT: A PyMOL Plugin for Molecular Glue Discovery and Analysis},
  author = {Chen, Roufen},
  year = {2024},
  url = {https://github.com/VesperChen01/GLINT}
}
```

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).

---

## 📬 Contact

<div align="center">

**Roufen Chen**

[![Email](https://img.shields.io/badge/Email-12319021@zju.edu.cn-EA4335?style=flat-square&logo=gmail&logoColor=white)](mailto:12319021@zju.edu.cn)
[![GitHub](https://img.shields.io/badge/GitHub-VesperChen01-181717?style=flat-square&logo=github)](https://github.com/VesperChen01)

</div>

---

## 🙏 Acknowledgments

- [PyMOL](https://pymol.org) - Visualization framework
- [RDKit](https://rdkit.org) - Cheminformatics tools
- The molecular glue research community

---

<div align="center">

Made with ❤️ for structural biology and drug discovery research

**⭐ Star this repo if you find it useful!**

</div>
