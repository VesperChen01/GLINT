#!/bin/bash
# GitHub Pages Setup Script for GLINT Documentation
# This script sets up a Jekyll-based documentation website

echo "=== GitHub Pages Documentation Setup ==="
echo ""

# Create necessary directories
mkdir -p docs/_includes
mkdir -p docs/_layouts
mkdir -p docs/_posts
mkdir -p docs/assets/css
mkdir -p docs/images

echo "✅ Created directory structure"

# Create Jekyll configuration file
cat > docs/_config.yml << 'EOF'
# GLINT Documentation Configuration

title: GLINT - PyMOL Plugin for Molecular Glue Discovery
description: Comprehensive PyMOL plugin for molecular glue discovery and analysis
theme: minima
markdown: kramdown
highlighter: rouge
url: "https://vesperchen01.github.io/GLINT"
baseurl: "/GLINT"

# Navigation
nav:
  - title: Home
    url: /
  - title: Installation
    url: /installation
  - title: Tutorial
    url: /tutorial
  - title: API Reference
    url: /api

# Plugins
plugins:
  - jekyll-seo-tag
  - jekyll-sitemap

# Exclude from processing
exclude:
  - Gemfile
  - Gemfile.lock
  - node_modules
  - vendor/bundle/
  - vendor/cache/
  - vendor/gems/
  - vendor/ruby/
  - "*.docx"
EOF

echo "✅ Created Jekyll configuration"

# Create custom CSS
cat > docs/assets/css/custom.css << 'EOF'
/* Custom styles for GLINT documentation */

:root {
  --primary-color: #3b82f6;
  --secondary-color: #10b981;
  --accent-color: #f59e0b;
  --text-color: #1e293b;
  --bg-color: #ffffff;
  --border-color: #e2e8f0;
}

body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  line-height: 1.6;
  color: var(--text-color);
}

/* Header styles */
.site-header {
  background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
  color: white;
  padding: 1rem 0;
}

.site-title {
  color: white !important;
  font-weight: 700;
}

/* Navigation */
.nav-links {
  display: flex;
  gap: 1.5rem;
  align-content: center;
}

.nav-links a {
  color: rgba(255, 255, 255, 0.9);
  text-decoration: none;
  transition: color 0.3s;
}

.nav-links a:hover {
  color: white;
}

/* Content styles */
.post-content {
  max-width: 800px;
  margin: 0 auto;
  padding: 2rem 1rem;
}

h1, h2, h3 {
  color: var(--text-color);
  margin-top: 2rem;
}

h1 {
  border-bottom: 3px solid var(--primary-color);
  padding-bottom: 0.5rem;
}

h2 {
  border-bottom: 2px solid var(--border-color);
  padding-bottom: 0.3rem;
}

code {
  background: #f1f5f9;
  padding: 0.2rem 0.4rem;
  border-radius: 4px;
  font-family: 'Courier New', monospace;
  font-size: 0.9em;
}

pre {
  background: #1e293b;
  color: #f8fafc;
  padding: 1rem;
  border-radius: 8px;
  overflow-x: auto;
}

pre code {
  background: none;
  padding: 0;
  color: inherit;
}

/* Table styles */
table {
  width: 100%;
  border-collapse: collapse;
  margin: 1.5rem 0;
}

th, td {
  padding: 0.75rem;
  text-align: left;
  border: 1px solid var(--border-color);
}

th {
  background: #f8fafc;
  font-weight: 600;
}

/* Alert boxes */
.alert {
  padding: 1rem;
  border-radius: 8px;
  margin: 1rem 0;
  border-left: 4px solid;
}

.alert-info {
  background: #eff6ff;
  border-color: var(--primary-color);
}

.alert-warning {
  background: #fffbeb;
  border-color: var(--accent-color);
}

.alert-success {
  background: #f0fdf4;
  border-color: var(--secondary-color);
}

/* Responsive design */
@media (max-width: 768px) {
  .nav-links {
    flex-direction: column;
    gap: 0.5rem;
  }
}
EOF

echo "✅ Created custom CSS"

# Create main page layout
cat > docs/_layouts/default.html << 'EOF'
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{% if page.title %}{{ page.title }} - {% endif %}{{ site.title }}</title>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css">
  <link rel="stylesheet" href="{{ '/assets/css/custom.css' | relative_url }}">
  {% seo %}
</head>
<body>
  <header class="site-header">
    <div class="container">
      <div class="d-flex justify-content-between align-items-center">
        <a class="site-title" href="{{ '/' | relative_url }}">{{ site.title }}</a>
        <nav class="nav-links">
          {% for item in site.nav %}
          <a href="{{ item.url | relative_url }}">{{ item.title }}</a>
          {% endfor %}
        </nav>
      </div>
    </div>
  </header>

  <main class="post-content">
    {{ content }}
  </main>

  <footer class="site-footer" style="background: #1e293b; color: white; padding: 2rem 0; margin-top: 3rem;">
    <div class="container">
      <div class="row">
        <div class="col-md-6">
          <h5>GLINT</h5>
          <p class="text-muted">PyMOL Plugin for Molecular Glue Discovery & Analysis</p>
        </div>
        <div class="col-md-6 text-md-end">
          <p class="text-muted small">
            © {{ site.time | date: '%Y' }} GLINT Project.<br>
            <a href="https://github.com/VesperChen01/GLINT" style="color: #94a3b8;">GitHub Repository</a>
          </p>
        </div>
      </div>
    </div>
  </footer>

  <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
EOF

echo "✅ Created default layout"

# Create index page
cat > docs/index.md << 'EOF'
---
layout: default
title: Home
---

# Welcome to GLINT Documentation

**GLINT** is a comprehensive PyMOL plugin designed for molecular glue discovery and protein-ligand interaction analysis.

## Quick Links

- [Installation Guide](installation) - Get started with GLINT
- [Tutorial](tutorial) - Step-by-step tutorials
- [API Reference](api) - Complete API documentation

## Features

### 🎯 Target Discovery
- G-Motif Detection via RMSD matching
- Protein Surface Analysis
- Surface Similarity & Complementarity

### 🔍 Hit Identification
- Integrated AutoDock Vina Docking
- HADDOCK3 Integration for ternary complex modeling
- FoldX Mutation Analysis

### 🧩 Ternary Complex Evaluation
- Automated ternary complex analysis with one-click evaluation
- Comprehensive interface analysis with BSA calculation
- Automatic ligand property extraction and computation
- Geometric metrics and duality index
- Enhanced CSV export consistent with GUI display

### ⚡ Lead Optimization
- Electrostatic Complementarity Analysis
- PPI Interface Analysis
- Protein-Ligand Interaction analysis

## Getting Started

1. **Install GLINT** - Follow the [Installation Guide](installation)
2. **Load Your Structure** - Open PyMOL and load your protein structure
3. **Run Analysis** - Use the GLINT menu to access analysis tools
4. **Explore Results** - Visualize and export your findings

## Support

- 📧 Email: support@glint-project.org
- 🐛 Issues: [GitHub Issues](https://github.com/VesperChen01/GLINT/issues)
- 💬 Discussions: [GitHub Discussions](https://github.com/VesperChen01/GLINT/discussions)

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/VesperChen01/GLINT/blob/main/LICENSE) file for details.
EOF

echo "✅ Created index page"

# Create installation page
cat > docs/installation.md << 'EOF'
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
EOF

echo "✅ Created installation page"

# Create API reference page
cat > docs/api.md << 'EOF'
---
layout: default
title: API Reference
---

# API Reference

Complete API documentation for GLINT modules.

## Ternary Complex Evaluator

### `TernaryComplexEvaluator`

Main class for evaluating ternary complexes.

```python
from glint.ternary_complex_evaluator import TernaryComplexEvaluator

evaluator = TernaryComplexEvaluator()
features = evaluator.evaluate(
    pdb_path="complex.pdb",
    e3_chain="A",
    poi_chain="B",
    ligand_resn="UNL",
    ligand_smiles="CC(=O)Nc1ccc(O)cc1",
    obj_name="my_complex"
)
```

### Methods

#### `evaluate(pdb_path, e3_chain, poi_chain, ligand_resn, ligand_smiles=None, obj_name=None)`

Evaluate a ternary complex and return comprehensive features.

**Parameters:**
- `pdb_path` (str): Path to PDB file
- `e3_chain` (str): E3 ligase chain ID
- `poi_chain` (str): POI chain ID
- `ligand_resn` (str): Ligand residue name
- `ligand_smiles` (str, optional): Ligand SMILES string
- `obj_name` (str, optional): PyMOL object name

**Returns:**
- `TernaryComplexFeatures`: Comprehensive feature object

## Ligand Calculator

### `LigandCalculator`

Calculate molecular properties for small molecules.

```python
from glint.ternary_complex_evaluator import LigandCalculator

calc = LigandCalculator()
props = calc.calculate(smiles="CC(=O)Nc1ccc(O)cc1")
```

### Properties

- `molecular_weight` (float): Molecular weight in Daltons
- `logp` (float): LogP value
- `tpsa` (float): Topological polar surface area
- `rotatable_bonds` (int): Number of rotatable bonds
- `hbd_count` (int): Hydrogen bond donor count
- `hba_count` (int): Hydrogen bond acceptor count
- `fsp3` (float): Fraction of sp3 carbons
- `num_rings` (int): Number of rings

## BSA Calculator

### `BSACalculator`

Calculate buried surface area for interfaces.

```python
from glint.ternary_complex_evaluator import BSACalculator

bsa_calc = BSACalculator(obj_name="my_complex", ligand_resn="UNL")
bsa_results = bsa_calc.calculate_bsa_ternary(e3_chain="A", poi_chain="B")
```

## GUI Components

### TernaryEvaluationTab

Main GUI tab for ternary complex evaluation.

```python
from glint.gui.tabs.ternary_evaluation import TernaryEvaluationTab

tab = TernaryEvaluationTab(parent_window)
```

### Methods

- `run_full_evaluation()`: Run complete ternary complex evaluation
- `run_ligand_only()`: Calculate ligand properties only
- `export_results()`: Export results to CSV
- `visualize_geometry()`: Visualize geometric features in PyMOL

## Data Structures

### `TernaryComplexFeatures`

```python
@dataclass
class TernaryComplexFeatures:
    interface: InterfaceFeatures
    ligand: LigandFeatures
    distances: DistanceFeatures
    geometry: GeometryFeatures
    balance_index: Optional[float]
```

### `InterfaceFeatures`

```python
@dataclass
class InterfaceFeatures:
    bsa_total: float
    bsa_mg_e3: float
    bsa_mg_poi: float
    bsa_e3_poi: float
    contact_count_45: int
    contact_count_50: int
    min_inter_chain_dist: float
```
EOF

echo "✅ Created API reference page"

# Create a placeholder for tutorial (will be replaced by converted content)
cat > docs/tutorial.md << 'EOF'
---
layout: default
title: Tutorial
---

# GLINT Tutorial

*This page will be updated with the converted content from the Word document.*

## Tutorial Overview

This tutorial will guide you through the main features of GLINT:

1. **Installation and Setup**
2. **Loading Structures**
3. **Target Discovery**
4. **Ternary Complex Evaluation**
5. **Hit Identification**
6. **Lead Optimization**

## Getting Started

To begin using GLINT, make sure you have:

- PyMOL installed and running
- GLINT plugin installed (see [Installation Guide](installation))
- A protein structure file (PDB format)

## Next Steps

*Full tutorial content coming soon after Word document conversion.*
EOF

echo "✅ Created tutorial placeholder"

# Create .nojekyll file (important for GitHub Pages)
touch docs/.nojekyll
echo "✅ Created .nojekyll file"

# Create GitHub Actions workflow for auto-deployment
mkdir -p .github/workflows
cat > .github/workflows/deploy-docs.yml << 'EOF'
name: Deploy Documentation

on:
  push:
    branches: [ main ]
    paths:
      - 'docs/**'
      - '.github/workflows/deploy-docs.yml'

permissions:
  contents: read
  pages: write
  id-token: write

concurrency:
  group: "pages"
  cancel-in-progress: false

jobs:
  deploy:
    environment:
      name: github-pages
      url: ${{ steps.deployment.outputs.page_url }}
    runs-on: ubuntu-latest
    steps:
      - name: Checkout
        uses: actions/checkout@v4

      - name: Setup Ruby
        uses: ruby/setup-ruby@v1
        with:
          ruby-version: '3.2'
          bundler-cache: true

      - name: Install dependencies
        run: |
          gem install jekyll jekyll-seo-tag jekyll-sitemap

      - name: Build site
        run: |
          cd docs
          jekyll build

      - name: Upload artifact
        uses: actions/upload-pages-artifact@v3
        with:
          path: docs/_site

      - name: Deploy to GitHub Pages
        id: deployment
        uses: actions/deploy-pages@v4
EOF

echo "✅ Created GitHub Actions workflow"

echo ""
echo "=== GitHub Pages Setup Complete! ==="
echo ""
echo "📁 Created files:"
echo "  - docs/_config.yml (Jekyll configuration)"
echo "  - docs/_layouts/default.html (Page layout)"
echo "  - docs/assets/css/custom.css (Custom styling)"
echo "  - docs/index.md (Home page)"
echo "  - docs/installation.md (Installation guide)"
echo "  - docs/tutorial.md (Tutorial - placeholder)"
echo "  - docs/api.md (API reference)"
echo "  - .github/workflows/deploy-docs.yml (Auto-deployment)"
echo ""
echo "🚀 Next steps:"
echo "1. Convert your Word document: bash convert_docx.sh"
echo "2. Review and edit the converted Markdown file"
echo "3. Test locally: cd docs && jekyll serve"
echo "4. Commit and push to GitHub"
echo "5. Enable GitHub Pages in repository settings"
echo ""
echo "📖 Documentation will be available at:"
echo "https://<username>.github.io/GLINT/"
