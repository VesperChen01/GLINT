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

### Dependencies
```bash
# Required: PyMOL (host application)
# Optional but recommended for advanced features:
pip install rdkit scipy matplotlib pillow numpy
```

### Testing
**Note**: This codebase does not have automated tests. All testing is done manually within PyMOL.

To test the plugin:
1. Launch PyMOL
2. Load a PDB structure: `fetch 1hsg`
3. Test basic commands:
   ```python
   analyze_protein_ligand_interactions('1hsg', 'MK1')
   gluetk_gui
   ```

### Running in PyMOL
```python
# In PyMOL command line:
run /path/to/gluetk/__init__.py

# Or load via Plugin Manager (GUI):
# Plugin → Plugin Manager → Install New Plugin
```

## Architecture Overview

### Module Structure

The plugin follows a **modular architecture** with clear separation of concerns:

```
gluetk/
├── __init__.py              # Plugin entry point, command registration
├── interaction_analyzer.py   # Core interaction detection engine
├── highlight_residues.py     # Visualization and highlighting
├── g_motif_analyzer.py       # CRBN G-motif/G-loop detection
├── interaction_2d_plot.py    # 2D interaction diagram generation
├── unified_gui.py            # Qt-based GUI (non-modal)
├── modern_style.py           # GUI styling
├── pymol_crbn_tools.py       # CRBN-specific analysis tools
└── README.md                 # User documentation (Chinese)
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

Specialized module for CRBN G-loop/G-motif identification:

- **Template Modes**: 
  - `ideal`: Theoretical β-hairpin geometry
  - `builtin`: Real PDB templates (GSPT1, CK1α, VAV1)
  - `selection`: User-provided custom template
- **RMSD-based matching** with Kabsch alignment
- Supports insertion codes and altlocs correctly
- Default RMSD cutoff: 3.5Å, optional Gly requirement at position 6

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
- `analyze_protein_ligand_interactions`
- `analyze_pdb_interactions`
- `highlight_csv_residues`
- `gluetk_gui`
- `visualize_protein_ligand_3d`
- `generate_interaction_network_plot`
- `generate_2d_diagram`
- `analyze_ternary_complex`
- `analyze_atom_pair_interactions`

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

**Last Updated**: 2025-11-06  
**Maintained by**: Vesper
