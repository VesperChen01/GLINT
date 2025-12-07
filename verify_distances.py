"""
Check if distance objects exist in PyMOL

Run in PyMOL:
run /Users/vesper/Desktop/git/GlueTK/GlueTK/verify_distances.py
"""

from pymol import cmd

print("\n" + "="*60)
print("Verifying Distance Objects")
print("="*60)

# Get all objects
all_objects = cmd.get_names("all")
print(f"\nTotal objects: {len(all_objects)}")

# Filter distance objects
distance_objects = [obj for obj in all_objects if obj.startswith("ppi_") and not obj.startswith("ppi_res_") and not obj.startswith("ppi_legend_")]

print(f"Distance objects found: {len(distance_objects)}\n")

if distance_objects:
    # Group by type
    by_type = {}
    for obj in distance_objects:
        parts = obj.split("_")
        if len(parts) >= 2:
            obj_type = parts[1]
            if obj_type not in by_type:
                by_type[obj_type] = []
            by_type[obj_type].append(obj)
    
    print("Distance objects by type:")
    for obj_type in sorted(by_type.keys()):
        objs = by_type[obj_type]
        print(f"\n  {obj_type}: {len(objs)} objects")
        for obj in objs[:3]:
            # Try to get distance value
            try:
                dist = cmd.get_distance(obj)
                print(f"    - {obj}: {dist:.2f} Å")
            except:
                print(f"    - {obj}: (cannot get distance)")
        if len(objs) > 3:
            print(f"    ... and {len(objs)-3} more")
else:
    print("⚠️  No distance objects found!")
    print("\nAll ppi_* objects:")
    ppi_objects = [obj for obj in all_objects if obj.startswith("ppi_")]
    for obj in ppi_objects[:10]:
        print(f"  - {obj}")
    if len(ppi_objects) > 10:
        print(f"  ... and {len(ppi_objects)-10} more")

print("\n" + "="*60 + "\n")
