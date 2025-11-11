# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Project Overview

**GlueTK - PyMOL Plugin for Molecular Glue Analysis** - A PyMOL plugin for protein-ligand interaction analysis and visualization, with support for strict Schrödinger standards suitable for publication.

**Language**: Python (PyMOL plugin)  
**Version**: 1.0.0  
**Author**: Vesper

## Development Commands

### Installation
This is a PyMOL plugin. Installation is manual:
1. Copy `gluetk/` directory to PyMOL's plugin directory
2. Or use PyMOL's Plugin Manager to install the ZIP file

### Environment setup (recommended)
- One-shot guided setup (creates Conda env `gluetk`, installs deps, verifies GUI):
  ```bash
  bash gluetk/check_env.sh
  ```
- Manual setup with Conda:
  ```bash
  conda create -n gluetk python=3.9 -y
  conda activate gluetk
  conda install -c conda-forge rdkit scipy matplotlib pillow numpy pyqt autodock-vina openbabel -y
  # Optional (headless dev without system PyMOL):
  pip install pymol-open-source
  ```
- Diagnose/auto-fix in current shell:
  ```bash
  python gluetk/env_checker.py            # report only
  python gluetk/env_checker.py --auto-install  # try to install missing pieces
  ```

### Running in PyMOL
```python
# In PyMOL command line:
run /path/to/gluetk/__init__.py

# Or load via Plugin Manager (GUI):
# Plugin → Plugin Manager → Install New Plugin
```

### Testing and validation
- There are no automated unit tests; testing is performed inside PyMOL.
- Quick smoke test:
  ```python
  fetch 1hsg
  analyze_protein_ligand_interactions('1hsg', 'MK1', use_schrodinger_standard=True)
  gluetk_gui
  ```
- Benchmark set used for manuscript validation:
  ```python
  run validation/benchmark_analysis.py
  benchmark_all()
  ```

### Packaging for distribution
```bash
zip -r gluetk.zip gluetk
```

### Hot-reload during development (in PyMOL)
```python
import sys
mods = [m for m in list(sys.modules) if 'gluetk' in m.lower()]
for m in mods:
    sys.modules.pop(m)
run /path/to/gluetk/__init__.py
```

## Architecture Overview

### Module Structure

The plugin follows a **modular architecture** with clear separation of concerns:

```
gluetk/
├── __init__.py              # Plugin entry, command registration, GUI hooks
├── interaction_analyzer.py  # Core interaction engine (strict criteria)
├── highlight_residues.py    # Visualization and CSV parsing utilities
├── g_motif_analyzer.py      # CRBN G-motif/G-loop detection
├── ppi_analyzer.py          # PPI interface + neo-epitope (Glue-specific)
├── binding_score.py         # Empirical binding-energy scoring (binary/ternary)
├── binding_heatmap.py       # Batch heatmap from *_scores.csv
├── vina_scoring.py          # AutoDock Vina integration (optional)
├── interaction_2d_plot.py   # 2D interaction diagrams
├── unified_gui.py           # Qt GUI (non-modal)
├── modern_style.py          # GUI styling
├── env_setup.py, env_checker.py  # Dependency detection/auto-install
└── README.md                # User documentation (Chinese)
```

### Key Design Patterns

1. **Plugin System Integration**
   - `__init__.py` uses `cmd.extend()` to register Python functions as PyMOL commands
   - GUI is launched non-modally to prevent blocking PyMOL's event loop
   - Global `_dlg` variable maintains single GUI instance

2. **Dual Standards System**
   - **Default Standard**: Fast screening (H-bond ≤3.5Å)
   - **Schrödinger Standard**: Publication-quality (H-bond ≤2.8Å, matching Maestro)
   - Controlled via `INTERACTION_PARAMS` dict and `use_schrodinger_standard` flag

3. **Flexible Input Handling**
   - Can work with PyMOL objects (`obj_name`) OR PDB files (`pdb_file`)
   - Automatic ligand detection using `organic` selection or explicit residue name
   - Robust CSV parsing supporting multiple encodings (UTF-8, GB18030, Shift_JIS, CP1252)

4. **Optional Dependencies**
   - RDKit detection with graceful fallback: `RDKIT_AVAILABLE` flag
   - If RDKit present: advanced analysis (2D diagrams, SMILES, precise geometry)
   - If absent: basic analysis still functional

5. **Locale-Aware I18N**
   - `_zh()` helper detects Chinese locale
   - `_info(cn, en)` prints messages in appropriate language
   - GUI uses translation dictionary `T` keyed by language code

### Core Interaction Detection

**Location**: `interaction_analyzer.py`

The plugin implements **strict, publication-quality interaction criteria**:

- **Hydrogen Bonds**: Distance + angle validation (D-H-A geometry)
- **Salt Bridges**: Charged atom proximity with H-bond exclusion logic
- **Hydrophobic**: Carbon-carbon contacts between hydrophobic residues
- **π-π Stacking**: Ring plane geometry and distance
- **π-Cation**: Aromatic ring to cationic residue
- **Metal Coordination**: Metal ion to heteroatom coordination

**Key Functions**:
- `analyze_protein_ligand_interactions()` - Main analysis entry point
- `analyze_pdb_interactions()` - Generic PDB interaction analysis
- `analyze_ternary_complex()` - PROTAC/ternary complex analysis
- `visualize_protein_ligand_3d()` - 3D PyMOL visualization

### G-Motif Detection

**Location**: `g_motif_analyzer.py`

Specialized module for CRBN G-loop/G-motif identification with **publication-quality validation** based on Annual Review of Pharmacology and Toxicology 2023:

#### Core Detection
- **Template Modes**: 
  - `ideal`: Theoretical β-hairpin geometry
  - `builtin`: Real PDB templates (GSPT1, CK1α, VAV1)
  - `selection`: User-provided custom template
- **RMSD-based matching** with Kabsch alignment
- Supports insertion codes and altlocs correctly
- Default RMSD cutoff: 3.5Å, optional Gly requirement at position 6

#### CRBN H-bond Validation (NEW)
**Function**: `validate_crbn_hbonds()`

Verifies the **3 canonical backbone H-bonds** between G-loop and CRBN:
1. **G-3 carbonyl O** ↔ **CRBN Asn351** sidechain (NH2)
2. **G-2 carbonyl O** ↔ **CRBN His357** sidechain (ND1/NE2)
3. **G-1 carbonyl O** ↔ **CRBN Trp400** sidechain (NE1)

**Criteria for canonical G-loop**: ≥2/3 H-bonds (default threshold 3.5Å)

**Literature basis**:
- Deep mutational scanning identified N351, H357, W400 as resistance hotspots
- CRBN N351D mutant (loss of H-bond donor) abolishes GSPT1 degradation
- Matches Schrödinger Maestro H-bond standards

#### Enhanced Glue Binding Analysis
**Function**: `analyze_g_motif_glue_binding()` - now includes:
- CRBN H-bond validation (via `validate_hbonds=True`)
- **Gly vdW contact** with MGD (conserved Gly faces drug, threshold 4.5Å)
- **Sidechain contacts** at G-4, G-3, G-2, G+1 positions with CRBN
- Neosubstrate classification (glue-induced vs direct binding)

**Validation**: Tested on known structures (GSPT1/6H0G, CK1α/5FQD, IKZF3/6H0F)

See `G_LOOP_VALIDATION_GUIDE.md` for detailed usage and `test_gloop_validation.py` for test suite.

### GUI Architecture

**Location**: `unified_gui.py`

- **Qt Framework**: PyQt5/PyQt6 compatibility layer
- **Non-Modal Design**: Prevents PyMOL UI freezing
- **Tabbed Interface**:
  1. G-Motif Detection
  2. Interaction Analysis & Highlight
  3. Electrostatics (APBS/Quick)
- **Worker Threads**: Background processing with `QThread` for long operations
- **Light Theme**: Professional appearance with custom styling

## Important Conventions

### Command Registration
All PyMOL commands are registered in `__init__.py`:
```python
cmd.extend("command_name", function_name)
```

Commands available in PyMOL:
- Core analysis and viz:
  - `analyze_protein_ligand_interactions`, `analyze_pdb_interactions`, `analyze_ternary_complex`, `analyze_atom_pair_interactions`
  - `visualize_protein_ligand_3d`, `generate_interaction_network_plot`, `generate_2d_diagram`, `highlight_csv_residues`
- Molecular glue–specific:
  - `ppi_analyze`, `analyze_protein_protein_interface`, `identify_neo_epitope`, `calculate_interface_bsa`
  - `find_crbn_g_motif`, `analyze_g_motif_glue_binding`, `validate_crbn_hbonds`, `validate_g_motif_geometry`
- Scoring and batch plots:
  - `score_protein_ligand`, `score_ternary_complex`, `plot_binding_heatmap`
  - If Vina available: `vina_score_complex`, `compare_scoring_methods`
- GUI:
  - `gluetk_gui` (preferred), `molstruct_gui` (legacy alias)

### CSV Format Standard
All interaction CSVs follow this schema:
```csv
Chain1,Residue1,Chain2,Residue2,Distance,Interaction
A,SER 99,L,LIG 301,2.75,氢键
```

Or for protein-ligand specific:
```csv
Ligand_Chain,Ligand_Residue,Ligand_Atom,Protein_Chain,Protein_Residue,Protein_Atom,Distance,Interaction
```

### Interaction Parameter Configuration

Edit `INTERACTION_PARAMS` in `interaction_analyzer.py` to customize:
```python
INTERACTION_PARAMS = {
    "hbond": {"max_DA_dist": 2.5, "min_donor_angle": 120, "min_acceptor_angle": 90},
    "ionic": {"max_dist": 3.7, "exclude_if_hbond": True},
    "hydrophobic": {"pi_cation_max": 4.5, "other_max": 3.6},
    "metal_coord": {"max_dist": 3.4},
    # ... etc
}
```

### Residue Parsing
The codebase handles complex residue identifiers:
- Insertion codes: `100A`
- Chain prefixes: `A:ARG:123` or `ARG 123`
- Altlocs: Preferentially uses `''` or `A`

Helper: `_parse_resi()`, `_split_residue_tag()`

## Code Style & Patterns

1. **UTF-8 Source Encoding**: All files declare `# -*- coding: utf-8 -*-`
2. **Print for User Communication**: Use `print()` for status messages (PyMOL has no logging framework)
3. **Docstrings**: Chinese + English mixed documentation
4. **Error Handling**: Try-except with user-friendly Chinese/English messages via `_info()`
5. **PyMOL Object Safety**: Always check `cmd.get_object_list()` before operations
6. **Temp Files**: Use `tempfile` module, clean up with `os.unlink()`

## Common Development Patterns

### Adding a New Interaction Type
1. Add parameters to `INTERACTION_PARAMS` in `interaction_analyzer.py`
2. Implement detection function `is_newtype()`
3. Call in `analyze_interactions()` loop
4. Update CSV export in `analyze_pdb_interactions()`
5. Update README.md with Chinese/English descriptions

### Adding a New GUI Tab
1. Extend `GlueTKDialog` in `unified_gui.py`
2. Add translations to `T` dictionary
3. Create tab with `QWidget` + layout
4. Implement worker thread if needed (inherit `QThread`)
5. Connect signals: `progress`, `finished`, `error`

### Working with PyMOL Objects
```python
# Always check existence first
if obj_name not in cmd.get_object_list():
    print(f"Object '{obj_name}' not found")
    return

# Get atom data
model = cmd.get_model(obj_name)
for atom in model.atom:
    chain = atom.chain.strip()
    resn = atom.resn.strip()
    coord = (atom.coord[0], atom.coord[1], atom.coord[2])
```

### CSV Robustness Pattern
Use `_read_csv_robust()` from `highlight_residues.py`:
- Auto-detects encoding
- Auto-detects delimiter
- Normalizes headers (lowercase, no spaces)
- Handles BOM markers

## Plugin-Specific Considerations

1. **PyMOL Environment**: Code runs inside PyMOL's Python interpreter
2. **No `__main__`**: Entry point is `__init_plugin__(app=None)`
3. **GUI Lifecycle**: Must handle PyMOL restart (reset `_dlg = None`)
4. **Command Namespace**: Avoid name conflicts with PyMOL built-ins
5. **Chinese Characters**: PyMOL labels don't support Unicode well - use `_interaction_to_abbr()`

## Useful PyMOL Commands for Development

```python
# Reload plugin after changes
run /path/to/__init__.py

# Debug: print all registered commands
cmd.keyword

# Check if command exists
'command_name' in dir(cmd)

# Get selection info
cmd.iterate('selection', 'print(chain, resn, resi, name)')

# Clear all highlights
cmd.delete('highlight_*')
```

## Publication Standards

This plugin implements **strict interaction criteria suitable for academic publication**:

- Hydrogen bonds: Maestro-compatible ≤2.8Å
- Salt bridges: ≤4.0Å with angle validation
- Results comparable to Schrödinger Maestro
- Enable with `use_schrodinger_standard=True`

Reference: Schrödinger Maestro documentation standards

---

**Last Updated**: 2025-11-10  
**Maintained by**: Vesper
