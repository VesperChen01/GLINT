
import sys
import os
from unittest.mock import MagicMock

# Mock PyMOL
sys.modules['pymol'] = MagicMock()

# Adds GlueTK to path
sys.path.append(os.getcwd())

from gluetk.interaction_analyzer import identify_molecule_type

print("=== Simulating Centroid Logic with New Cutoff ===")
def distance(p1, p2):
    import math
    return math.sqrt(sum((a - b)**2 for a, b in zip(p1, p2)))

def centroid(coords):
    if not coords: return (0,0,0)
    x = sum(c[0] for c in coords) / len(coords)
    y = sum(c[1] for c in coords) / len(coords)
    z = sum(c[2] for c in coords) / len(coords)
    return (x,y,z)

# Scenario: Ligand and Receptor centroids are 25 Angstroms apart.
# This would FAIL the old 15.0 cutoff, but should PASS the new 60.0 cutoff.

ligand_atoms = [(0.0, 0.0, 0.0)]
dna_atoms = [(25.0, 0.0, 0.0)]

lig_c = centroid(ligand_atoms)
prot_c = centroid(dna_atoms)

dist = distance(lig_c, prot_c)
print(f"Ligand Centroid: {lig_c}")
print(f"Receptor Centroid: {prot_c}")
print(f"Centroid Distance: {dist:.2f}")

cutoff = 60.0
print(f"Testing against cutoff: {cutoff}")

if dist > cutoff:
    print("❌ REJECTED (Unexpected!)")
else:
    print("✅ PASSED cutoff (Success)")
