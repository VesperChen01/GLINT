<p align="center">
  <img src="gluetk/assets/logo.png" alt="GlueTK Logo" width="200"/>
</p>

<h1 align="center">GlueTK</h1>

<p align="center">
  <strong>PyMOL Plugin for Molecular Glue Analysis</strong>
</p>

<p align="center">
  <a href="#features">Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#documentation">Documentation</a> •
  <a href="#license">License</a>
</p>

---

## Overview

**GlueTK** is a comprehensive PyMOL plugin for protein-ligand interaction analysis and visualization, with specialized support for molecular glue analysis. It implements strict Schrödinger standards suitable for publication-quality results.

## Features

- 🔬 **Interaction Analysis**: Hydrogen bonds, salt bridges, π-π stacking, π-cation, hydrophobic contacts, metal coordination
- 🧬 **G-Motif Detection**: CRBN G-loop/G-motif recognition with RMSD-based template matching
- 🔗 **PPI Analysis**: Protein-protein interface detection and neo-epitope identification
- ⚡ **Electrostatics**: Quick Coulomb electrostatics and APBS integration
- 📊 **2D Diagrams**: RDKit-powered 2D interaction diagrams
- 🎨 **Beautiful Visualization**: Publication-ready 3D rendering with customizable styles
- 🎯 **Dual Standards**: Default (fast screening) and Schrödinger-compatible (publication) modes

## Installation

### Option 1: One-Shot Setup (Recommended)

```bash
bash gluetk/check_env.sh
```

This creates a Conda environment `gluetk` with all dependencies.

### Option 2: Manual Installation

```bash
conda create -n gluetk python=3.9 -y
conda activate gluetk
conda install -c conda-forge rdkit scipy matplotlib pillow numpy pandas seaborn pyqt -y
```

### Option 3: PyMOL Plugin Manager

1. Download `gluetk.zip` from releases
2. In PyMOL: `Plugin → Plugin Manager → Install New Plugin`
3. Select the ZIP file

## Quick Start

### Launch GUI

```python
# In PyMOL command line:
run /path/to/gluetk/__init__.py
gluetk_gui
```

### Command Line Usage

```python
# Load a structure
fetch 1hsg

# Analyze protein-ligand interactions
analyze_protein_ligand_interactions('1hsg', 'MK1', use_schrodinger_standard=True)

# Detect G-motif
find_crbn_g_motif('structure_name', rmsd_cutoff=3.5)

# Generate 2D diagram
generate_2d_diagram('structure_name', 'LIG')
```

## Available Commands

### Core Analysis
- `analyze_protein_ligand_interactions` - Protein-ligand interaction detection
- `analyze_pdb_interactions` - Generic PDB interaction analysis
- `analyze_ternary_complex` - PROTAC/ternary complex analysis
- `analyze_atom_pair_interactions` - Atom-level interaction analysis

### Molecular Glue
- `find_crbn_g_motif` - G-motif detection
- `analyze_g_motif_glue_binding` - Glue binding analysis
- `validate_crbn_hbonds` - CRBN H-bond validation
- `ppi_analyze` - PPI interface analysis
- `identify_neo_epitope` - Neo-epitope detection

### Visualization
- `visualize_protein_ligand_3d` - 3D visualization
- `generate_2d_diagram` - 2D interaction diagram
- `highlight_csv_residues` - CSV-based highlighting
- `render_interactions_beautifully` - Publication-quality rendering

### Scoring
- `score_protein_ligand` - Empirical scoring
- `score_ternary_complex` - Ternary complex scoring
- `plot_binding_heatmap` - Batch heatmap generation

## Documentation

See [WARP.md](WARP.md) for detailed development documentation.

## Requirements

- PyMOL (Open-Source or Incentive)
- Python 3.8+
- PyQt5 or PyQt6
- Optional: RDKit, SciPy, Matplotlib, NumPy, Pandas

## License

This project is licensed under the MIT License.

## Author

**Vesper** - *Initial work and maintenance*

---

<p align="center">
  Made with ❤️ for structural biology research
</p>
