#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Quick reload script for GlueTK plugin
Run this in PyMOL command line: run /Users/vesper/Desktop/git/glue-pymol/reload_plugin.py
"""

import sys
import os

# Add plugin directory to path
plugin_dir = "/Users/vesper/Desktop/git/glue-pymol/gluetk"
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

# Remove old modules from cache
modules_to_remove = [k for k in sys.modules.keys() if 'gluetk' in k.lower() or 'molstruct' in k.lower() or 'unified_gui' in k.lower() or 'modern_style' in k.lower()]
for mod in modules_to_remove:
    del sys.modules[mod]

print(f"[Reload] Removed {len(modules_to_remove)} cached modules")

# Reload the plugin
try:
    from pymol import cmd
    
    # Remove existing dialog if present
    try:
        import unified_gui
        if hasattr(unified_gui, '_dlg') and unified_gui._dlg is not None:
            unified_gui._dlg.close()
            unified_gui._dlg = None
            print("[Reload] Closed existing dialog")
    except:
        pass
    
    # Re-run the init file
    init_file = os.path.join(plugin_dir, "__init__.py")
    with open(init_file, 'r', encoding='utf-8') as f:
        exec(f.read(), {'__file__': init_file})
    
    print("[Reload] ✓ Plugin reloaded successfully!")
    print("[Reload] Run 'gluetk_gui' to open the interface")
    
except Exception as e:
    print(f"[Reload] ✗ Error: {e}")
    import traceback
    traceback.print_exc()
