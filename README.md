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
- **Ternary Complex Evaluation** - BSA/contact analysis, molecular properties, and ternary complex metrics

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
