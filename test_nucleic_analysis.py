import sys
import os
import unittest
from unittest.mock import MagicMock

# Mock pymol module before importing gluetk
sys.modules['pymol'] = MagicMock()
sys.modules['pymol.cmd'] = MagicMock()

# Add the current directory to sys.path so we can import gluetk
sys.path.append(os.getcwd())

from gluetk.interaction_analyzer import analyze_protein_nucleic_interactions

class TestNucleicAnalysis(unittest.TestCase):
    def setUp(self):
        # Mock atoms: 
        # Chain A: Protein (ARG 10)
        # Chain B: DNA (DG 1)
        
        # Coordinates setup for a Salt Bridge and H-bond
        # ARG 10 NH1 at (0, 0, 0)
        # DG 1 OP1 at (0, 0, 3.0) -> Salt Bridge (dist 3.0)
        
        self.mock_atoms = [
            # Chain A (Protein) - ARG 10
            ("A", "ARG", "10", "N", (2.0, 0.0, 0.0)),
            ("A", "ARG", "10", "CA", (2.5, 0.0, 0.0)),
            ("A", "ARG", "10", "C", (3.0, 0.0, 0.0)),
            ("A", "ARG", "10", "O", (3.0, 1.0, 0.0)),
            ("A", "ARG", "10", "CB", (2.0, 1.0, 0.0)),
            ("A", "ARG", "10", "CG", (1.5, 2.0, 0.0)),
            ("A", "ARG", "10", "CD", (1.0, 3.0, 0.0)),
            ("A", "ARG", "10", "NE", (0.5, 4.0, 0.0)),
            ("A", "ARG", "10", "CZ", (0.0, 5.0, 0.0)),
            ("A", "ARG", "10", "NH1", (0.0, 0.0, 0.0)), # Positive charge center
            ("A", "ARG", "10", "NH2", (0.5, 5.0, 0.0)),
            
            # Chain B (DNA) - DG 1
            ("B", "DG", "1", "P", (0.0, 0.0, 4.0)),
            ("B", "DG", "1", "OP1", (0.0, 0.0, 3.0)), # Negative charge center, dist=3.0 to NH1
            ("B", "DG", "1", "OP2", (0.0, 1.0, 4.0)),
            ("B", "DG", "1", "O5'", (1.0, 0.0, 4.0)),
            ("B", "DG", "1", "C5'", (2.0, 0.0, 4.0)),
            # ... simplified
        ]

    def test_analysis(self):
        # Mock parse_pdb_structure to return our mock atoms
        import gluetk.interaction_analyzer
        gluetk.interaction_analyzer.parse_pdb_structure = MagicMock(return_value=self.mock_atoms)
        
        # Run analysis
        result = analyze_protein_nucleic_interactions(obj_name="mock_obj")
        
        # Verify results
        self.assertIsNotNone(result)
        self.assertIn("B", result["nucleic_chains"])
        self.assertIn("A", result["protein_chains"])
        
        interactions = result["interactions"]
        self.assertTrue(len(interactions) > 0)
        
        # Check for Salt Bridge
        salt_bridges = [i for i in interactions if i["Interaction"] == "盐桥"]
        self.assertTrue(len(salt_bridges) > 0)
        sb = salt_bridges[0]
        self.assertEqual(sb["Protein_Residue"], "ARG 10")
        self.assertEqual(sb["Nucleic_Residue"], "DG 1")
        self.assertAlmostEqual(sb["Distance"], 3.0)
        
        print("\nTest Passed! Found interactions:")
        for i in interactions:
            print(f"  - {i['Interaction']}: {i['Protein_Residue']} - {i['Nucleic_Residue']} (d={i['Distance']})")

if __name__ == '__main__':
    unittest.main()
