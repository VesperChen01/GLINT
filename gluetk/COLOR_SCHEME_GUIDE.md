# GlueTK Unified Color Scheme Guide

## Overview

All GlueTK modules now use a **unified Schrödinger-inspired color scheme** defined in [`color_scheme.py`](color_scheme.py). This ensures visual consistency across:

- 2D interaction diagrams
- 3D PyMOL visualizations  
- PPI (protein-protein interface) analysis
- GUI components

## Color Definitions

### Interaction Types

| Interaction | Hex Color | RGB (0-1) | PyMOL Name | Description |
|-------------|-----------|-----------|------------|-------------|
| **Hydrogen Bond** | `#2196F3` | `(0.129, 0.588, 0.953)` | `glue_hbond` | Blue - H-bonds |
| **Salt Bridge** | `#FF5722` | `(1.000, 0.341, 0.133)` | `glue_salt` | Orange-red - Ionic |
| **Pi-Pi Stacking** | `#9C27B0` | `(0.612, 0.153, 0.690)` | `glue_pipi` | Purple - π-π |
| **Pi-Cation** | `#E91E63` | `(0.914, 0.118, 0.388)` | `glue_pication` | Pink - π-cation |
| **Hydrophobic** | `#4CAF50` | `(0.298, 0.686, 0.314)` | `glue_hydrophobic` | Green - Hydrophobic |
| **Halogen Bond** | `#FF9800` | `(1.000, 0.596, 0.000)` | `glue_halogen` | Orange - Halogen |
| **Metal Coordination** | `#673AB7` | `(0.404, 0.227, 0.718)` | `glue_metal` | Deep purple - Metal |
| **Water Bridge** | `#00BCD4` | `(0.000, 0.737, 0.831)` | `glue_water` | Cyan - Water |
| **Other/Unknown** | `#9E9E9E` | `(0.620, 0.620, 0.620)` | `glue_other` | Gray - Other |

### Residue Types (Discovery Studio Style)

| Type | Face Color | Edge Color | Text Color |
|------|------------|------------|------------|
| **Hydrophobic** | `#C5E1A5` | `#4CAF50` | `#2E7D32` |
| **Nonpolar** | `#FFE0B2` | `#FF9800` | `#E65100` |
| **Polar** | `#BBDEFB` | `#2196F3` | `#1565C0` |
| **Negative** | `#A5D6A7` | `#388E3C` | `#1B5E20` |
| **Positive** | `#FFCDD2` | `#E57373` | `#C62828` |

## Module Integration

### Updated Modules

All modules now import from the unified color scheme:

1. **[`interaction_2d_plot.py`](interaction_2d_plot.py)** ✅
   - Uses `INTERACTION_COLORS_HEX` for matplotlib plots
   - Uses `RESIDUE_COLORS_HEX` for residue bubbles
   - Uses `INTERACTION_LINE_STYLES` for line rendering

2. **[`interaction_analyzer.py`](interaction_analyzer.py)** ✅
   - Uses `INTERACTION_COLORS_PYMOL` for 3D visualization
   - Uses `register_pymol_colors()` to register custom colors
   - Uses `PYMOL_COLOR_NAMES` for color references

3. **[`ppi_analyzer.py`](ppi_analyzer.py)** ✅
   - Uses `PYMOL_COLOR_NAMES` for PPI visualization
   - Registers colors via `register_pymol_colors()`
   - Consistent with interaction_analyzer.py

4. **[`highlight_residues.py`](highlight_residues.py)** ✅
   - Imports `INTERACTION_COLORS_HEX` for reference
   - Maps to nearest PyMOL built-in colors for compatibility

## Usage Examples

### In Python Modules

```python
# Import the unified color scheme
from gluetk.color_scheme import (
    INTERACTION_COLORS_HEX,
    INTERACTION_COLORS_PYMOL,
    PYMOL_COLOR_NAMES,
    get_interaction_color,
    register_pymol_colors,
)

# Get color in different formats
hbond_hex = get_interaction_color('hbond', format='hex')      # '#2196F3'
hbond_pymol = get_interaction_color('hbond', format='pymol')  # [0.129, 0.588, 0.953]
hbond_name = get_interaction_color('hbond', format='name')    # 'glue_hbond'

# Register all colors in PyMOL
from pymol import cmd
register_pymol_colors(cmd)

# Use the registered color
cmd.color('glue_hbond', 'my_selection')
```

### In Matplotlib Plots

```python
from gluetk.color_scheme import INTERACTION_COLORS_HEX, get_interaction_style

# Get color for plotting
color = INTERACTION_COLORS_HEX['hbond']
plt.plot(x, y, color=color)

# Get complete style (color, linewidth, linestyle)
style = get_interaction_style('hbond')
plt.plot(x, y, color=style['color'], linewidth=style['linewidth'], 
         linestyle=style['linestyle'], label=style['label'])
```

### In PyMOL Scripts

```python
from pymol import cmd
from gluetk.color_scheme import register_pymol_colors, PYMOL_COLOR_NAMES

# Register all custom colors
register_pymol_colors(cmd)

# Use them in your visualizations
cmd.color(PYMOL_COLOR_NAMES['hbond'], 'hbond_selection')
cmd.color(PYMOL_COLOR_NAMES['salt'], 'saltbridge_selection')
```

## Line Styles

Each interaction type has a defined line style for consistent visualization:

```python
from gluetk.color_scheme import INTERACTION_LINE_STYLES

# Example: Hydrogen bond style
hbond_style = INTERACTION_LINE_STYLES['hbond']
# {
#     'color': '#2196F3',
#     'linewidth': 2.0,
#     'linestyle': '--',
#     'label': 'Hydrogen Bond',
#     'dash_width': 2.0,
#     'dash_gap': 0.3,
# }
```

## Benefits of Unified Colors

✅ **Consistency** - Same colors across 2D diagrams, 3D views, and GUI  
✅ **Professional** - Schrödinger-inspired color palette  
✅ **Maintainability** - Single source of truth for all colors  
✅ **Flexibility** - Multiple formats (hex, RGB, PyMOL names) for different contexts  
✅ **Documentation** - Clear mapping of interaction types to colors

## Migration Notes

All color definitions have been migrated from individual modules to `color_scheme.py`. The old hardcoded colors have been replaced with imports from the unified scheme. This ensures:

- No color conflicts between modules
- Easy updates (change once, applies everywhere)
- Consistent user experience across all GlueTK features

## Future Extensions

To add new interaction types or modify colors:

1. Update [`color_scheme.py`](color_scheme.py) with new definitions
2. Colors automatically propagate to all modules
3. No need to modify individual visualization modules