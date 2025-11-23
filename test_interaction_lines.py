#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test script to verify protein-ligand interaction lines are visible
测试蛋白-配体相互作用线条是否可见

Usage in PyMOL:
    run /Users/vesper/Desktop/git/GlueTK/GlueTK/test_interaction_lines.py
    test_interaction_lines()
"""

def test_interaction_lines():
    """
    Test that interaction lines are visible in protein-ligand visualization
    测试蛋白-配体可视化中相互作用线条是否可见
    """
    from pymol import cmd
    
    print("\n" + "="*60)
    print("Testing GlueTK Interaction Line Visualization")
    print("测试 GlueTK 相互作用线条可视化")
    print("="*60 + "\n")
    
    # Step 1: Fetch a test structure
    print("Step 1: Fetching test structure (1hsg)...")
    cmd.delete("all")
    cmd.fetch("1hsg")
    
    # Step 2: Run the interaction analysis
    print("\nStep 2: Analyzing protein-ligand interactions...")
    cmd.do("run /Users/vesper/Desktop/git/GlueTK/GlueTK/gluetk/__init__.py")
    result = cmd.do("analyze_protein_ligand_interactions('1hsg', 'MK1', use_schrodinger_standard=True)")
    
    # Step 3: Visualize with interaction lines
    print("\nStep 3: Visualizing with interaction lines...")
    cmd.do("visualize_protein_ligand_3d('1hsg', ligand_resname='MK1')")
    
    # Step 4: Check if interaction objects exist
    print("\nStep 4: Checking interaction objects...")
    all_objects = cmd.get_names("objects")
    interaction_objects = [obj for obj in all_objects if "interact_" in obj or "hbonds_" in obj]
    
    if interaction_objects:
        print(f"✅ Found {len(interaction_objects)} interaction objects:")
        for obj in interaction_objects[:5]:  # Show first 5
            print(f"   - {obj}")
        if len(interaction_objects) > 5:
            print(f"   ... and {len(interaction_objects)-5} more")
    else:
        print("❌ No interaction objects found!")
    
    # Step 5: Toggle visibility to test
    print("\nStep 5: Testing toggle function...")
    print("   Hiding lines...")
    cmd.do("toggle_interaction_lines(False)")
    
    import time
    time.sleep(1)
    
    print("   Showing lines with thicker width...")
    cmd.do("toggle_interaction_lines(True, 3.5)")
    
    # Step 6: Summary
    print("\n" + "="*60)
    print("Test Complete! 测试完成！")
    print("="*60)
    print("\nYou should now see:")
    print("  1. Protein-ligand complex with MK1 highlighted")
    print("  2. Colored dashed lines showing interactions:")
    print("     - Blue: Hydrogen bonds")
    print("     - Orange: Salt bridges")
    print("     - Green: π-π stacking")
    print("     - Magenta: π-cation interactions")
    print("\nUseful commands:")
    print("  toggle_interaction_lines(False)    # Hide all lines")
    print("  toggle_interaction_lines(True, 4)  # Show with width=4")
    print("  show labels, interact_*            # Show distance labels")
    print("  set dash_color, red, interact_*    # Change all to red")
    
    return interaction_objects

# If run directly in PyMOL
if __name__ == "pymol":
    test_interaction_lines()