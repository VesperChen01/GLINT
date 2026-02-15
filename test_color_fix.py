#!/usr/bin/env python3
"""
Test script to verify color registration in PyMOL
"""

# This script can be run in PyMOL to test the color fix
print("=" * 60)
print("Testing GLINT Color Registration")
print("=" * 60)

# Reload the module
import sys
import importlib

# Reload pymol_styles
if 'glint.pymol_styles' in sys.modules:
    importlib.reload(sys.modules['glint.pymol_styles'])
    print("✓ Reloaded glint.pymol_styles")
else:
    import glint.pymol_styles
    print("✓ Imported glint.pymol_styles")

# Test color registration
from pymol import cmd

test_colors = ['deeporange', 'wheat', 'tan', 'limon']
print("\nTesting custom colors:")
for color in test_colors:
    try:
        # Try to use the color
        cmd.set_color(f'test_{color}', cmd.get_color_tuple(color))
        print(f"  ✓ {color}: OK")
    except Exception as e:
        print(f"  ✗ {color}: FAILED - {e}")

print("\n" + "=" * 60)
print("Test complete. You can now try running the analysis again.")
print("=" * 60)
