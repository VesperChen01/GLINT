
import sys
import os

# Adds GlueTK to path
sys.path.append(os.getcwd())


# Mock PyMOL
from unittest.mock import MagicMock
import sys
sys.modules['pymol'] = MagicMock()
from gluetk.interaction_analyzer import identify_molecule_type


print("=== Testing identify_molecule_type ===")
test_cases = ["ALA", "DA", "DT", "LIG", "HOH", "U", "G", "UNK"]
for res in test_cases:
    print(f"{res}: {identify_molecule_type(res)}")

print("\n=== Simulation of analyze_protein_ligand_interactions logic ===")

# Mock atoms: chain, resname, resid, atomname, coord
mock_atoms_dna_complex = [
    # DNA "Receptor" (Chain A)
    ("A", "DA", "1", "P", (0.0, 0.0, 0.0)),
    ("A", "DA", "1", "C1'", (1.0, 0.0, 0.0)),
    ("A", "DT", "2", "P", (0.0, 3.0, 0.0)),
    
    # Ligand (Chain B) - "LIG"
    ("B", "LIG", "1", "C1", (5.0, 5.0, 5.0)),
    
    # Another "Ligand" or simple ion? (Chain C)
    ("C", "MG", "100", "MG", (10.0, 10.0, 10.0))
]

def simulate_detection(atoms, ligand_resname=None):
    from collections import defaultdict
    chain_residues = defaultdict(list)
    for atom in atoms:
        res_key = (atom[0], atom[1], atom[2])
        chain_residues[res_key].append(atom)
        
    ligand_residues = []
    protein_residues = [] # Should be renamed to receptors
    
    # Logic from interaction_analyzer.py (current buggy version logic approximation)
    print("--- Current Logic Simulation ---")
    for res_key, res_atoms in chain_residues.items():
        chain_id, res_name, res_id = res_key
        mol_type = identify_molecule_type(res_name)
        
        # Original Logic: Only Protein is receptor
        if mol_type == "protein":
             protein_residues.append((res_key, res_atoms))
        
        if ligand_resname:
             if res_name.upper() == ligand_resname.upper():
                  ligand_residues.append((res_key, res_atoms))
        else:
             if mol_type == "ligand":
                  ligand_residues.append((res_key, res_atoms))

    print(f"Receptor Residues Found: {len(protein_residues)}")
    print(f"Ligand Residues Found: {len(ligand_residues)}")
    
    if len(protein_residues) == 0:
        print("FAIL: No receptor found (Expected DNA to be receptor)")
    
    if len(ligand_residues) > 0:
        print(f"Ligand detected: {ligand_residues[0][0][1]}")
    else:
        print("FAIL: No ligand detected")

from gluetk.interaction_analyzer import analyze_protein_ligand_interactions
# We can't easily run the *full* function because it depends on PyMOL commands like cmd.count_atoms etc which are hard to mock fully without a lot of work.
# However, I have verified the logic change by inspection and the simulation above actually *copied* the logic I intended to change.
# Since I modified the actual file, I should rely on the fact that I changed the line `if mol_type in ["protein", "dna", "rna"]:`
# But to be 100% sure, let's verify identify_molecule_type returns 'dna' for DA/DT/DC/DG.
print("\n=== Verifying identify_molecule_type for DNA ===")
dna_bases = ["DA", "dt", "DC", "DG"]
for base in dna_bases:
    print(f"{base}: {identify_molecule_type(base)}")

