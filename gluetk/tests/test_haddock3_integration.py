
import os
import sys
import unittest
from unittest.mock import MagicMock, patch

# Ensure we can import the module
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from gluetk.haddock3_integration import Haddock3Runner

class TestHaddock3Runner(unittest.TestCase):
    def setUp(self):
        # Mock availability check
        self.patcher = patch('gluetk.haddock3_integration.check_haddock3_available')
        self.mock_check = self.patcher.start()
        self.mock_check.return_value = {
            'available': True,
            'via': 'cli',
            'detail': '/mock/haddock3'
        }
        self.runner = Haddock3Runner(use_cli=True)

    def tearDown(self):
        self.patcher.stop()

    @patch('gluetk.haddock3_integration.subprocess.run')
    @patch('gluetk.haddock3_integration.Haddock3Runner._gather_pdbs')
    @patch('pathlib.Path.mkdir')
    @patch('pathlib.Path.write_text')
    @patch('os.path.exists')
    def test_run_docking_no_pdbs(self, mock_exists, mock_write, mock_mkdir, mock_gather, mock_run):
        # Setup mocks
        mock_exists.return_value = True # Simulate files exist
        mock_run.return_value = MagicMock(
            stdout="Running HADDOCK3...\nStep 1: topoaa\nStep 2: rigidbody\nFinished with errors.",
            returncode=0
        )
        mock_gather.return_value = [] # Simulate no PDBs found
        
        # Call run_docking
        result = self.runner.run_docking(
            receptor_pdb="rec.pdb",
            ligand_pdb="lig.pdb",
            output_dir="/tmp/haddock_out"
        )
        
        # Assertions
        self.assertFalse(result['success'])
        self.assertIn("未在 rigidbody 结果中找到 PDB 模型", result['error'])
        self.assertIn("Run Directory: /tmp/haddock_out/run1", result['error'])
        self.assertIn("Log Tail:", result['error'])
        self.assertIn("Finished with errors.", result['error'])

if __name__ == "__main__":
    unittest.main()
