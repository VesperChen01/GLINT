"""
List all distance objects and try to enable them

Run in PyMOL:
run /Users/vesper/Desktop/git/GlueTK/GlueTK/show_distances.py
"""

from pymol import cmd

print("\n" + "="*60)
print("Showing Distance Objects")
print("="*60)

# Get all objects
all_objects = cmd.get_names("all")
distance_objects = [obj for obj in all_objects if obj.startswith("ppi_") and 
                    ("hbond" in obj or "saltbridge" in obj or "hydrophobic" in obj or 
                     "pipi" in obj or "cationpi" in obj)]

print(f"\nFound {len(distance_objects)} distance objects\n")

if distance_objects:
    # Try to enable/show them
    for obj in distance_objects[:5]:  # Show first 5
        print(f"Enabling: {obj}")
        try:
            cmd.enable(obj)
            cmd.show("dashes", obj)
        except Exception as e:
            print(f"  Error: {e}")
    
    print(f"\n... and {len(distance_objects)-5} more" if len(distance_objects) > 5 else "")
    
    # Try to get object info
    print("\nObject info:")
    for obj in distance_objects[:3]:
        try:
            info = cmd.get_object_list(obj)
            print(f"  {obj}: {info}")
        except:
            print(f"  {obj}: (no info)")
    
    print("\n✅ Distance objects exist!")
    print("   If you don't see dashes, try:")
    print("   1. enable ppi_*")
    print("   2. show dashes")
    print("   3. zoom")
else:
    print("❌ No distance objects found")

print("="*60 + "\n")
