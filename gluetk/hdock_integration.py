# -*- coding: utf-8 -*-
"""
HDOCK Integration Module
Wrapper for HDOCKlite protein-protein docking.
"""
import os
import subprocess
import shutil
import tempfile
from typing import Optional, List, Dict

class HDockRunner:
    def __init__(self, hdock_dir: str):
        """
        Initialize with path to HDOCKlite directory
        """
        self.hdock_dir = hdock_dir
        self.hdock_exe = os.path.join(hdock_dir, "hdock")
        self.createpl_exe = os.path.join(hdock_dir, "createpl")
        
        # Validate executables
        self._validate_executables()

    def _validate_executables(self):
        if not os.path.exists(self.hdock_exe):
            raise FileNotFoundError(f"HDOCK executable not found at {self.hdock_exe}")
        if not os.path.exists(self.createpl_exe):
            # Typo in README says 'creapl' or 'createpl', check actual file
            # The ls output showed 'createpl*'
            if os.path.exists(os.path.join(self.hdock_dir, "creapl")):
                self.createpl_exe = os.path.join(self.hdock_dir, "creapl")
            elif not os.path.exists(self.createpl_exe):
                 raise FileNotFoundError(f"createpl executable not found at {self.createpl_exe}")

    def run_docking(self, receptor_pdb: str, ligand_pdb: str, output_dir: str = None) -> Dict:
        """
        Run HDOCK docking
        
        Args:
            receptor_pdb: Path to receptor PDB (e.g. E3 + Glue)
            ligand_pdb: Path to ligand PDB (e.g. POI)
            output_dir: Directory to save results
            
        Returns:
            dict: Result info {'success': bool, 'output_file': str, 'models_file': str}
        """
        if output_dir is None:
            output_dir = os.path.join(os.getcwd(), "hdock_results")
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Filenames
        out_name = "Hdock.out"
        out_path = os.path.join(output_dir, out_name)
        models_name = "model_top10.pdb" # Generate top 10 models
        models_path = os.path.join(output_dir, models_name)
        
        # 1. Run HDOCK
        # Usage: hdock receptor.pdb ligand.pdb -out Hdock.out
        cmd_dock = [
            self.hdock_exe,
            receptor_pdb,
            ligand_pdb,
            "-out", out_name
        ]
        
        print(f"[HDOCK] Running docking: {' '.join(cmd_dock)}")
        try:
            # Run in output dir to avoid clutter
            process = subprocess.run(
                cmd_dock,
                cwd=output_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            return {'success': False, 'error': f"HDOCK failed: {e.stdout}"}
            
        if not os.path.exists(out_path):
            return {'success': False, 'error': "HDOCK did not generate output file."}
            
        # 2. Generate Models
        # Usage: createpl Hdock.out top100.pdb -nmax 100 -complex -models
        # Let's generate top 10 for visualization
        cmd_create = [
            self.createpl_exe,
            out_name,
            models_name,
            "-nmax", "10",
            "-complex",
            "-models"
        ]
        
        print(f"[HDOCK] Generating models: {' '.join(cmd_create)}")
        try:
            process = subprocess.run(
                cmd_create,
                cwd=output_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                check=True
            )
        except subprocess.CalledProcessError as e:
            return {'success': False, 'error': f"createpl failed: {e.stdout}"}
            
        if not os.path.exists(models_path):
             return {'success': False, 'error': "Failed to generate models PDB."}
             
        return {
            'success': True,
            'output_dir': output_dir,
            'docking_out': out_path,
            'models_pdb': models_path
        }

def check_hdock_available() -> Optional[str]:
    """
    Check if HDOCK is available in standard locations.
    Returns path to HDOCKlite directory if found, else None.
    """
    # Check relative to this file (dev env)
    # gluetk/hdock_integration.py -> gluetk/ -> .../GlueTK/HDOCKlite-v1.1
    # Current: /Volumes/data/git/GlueTK/gluetk/hdock_integration.py
    # Root: /Volumes/data/git/GlueTK
    
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here) # gluetk root
    project_root = os.path.dirname(root) # git root
    
    potential_paths = [
        os.path.join(project_root, "HDOCKlite-v1.1"),
        os.path.join(root, "HDOCKlite-v1.1"),
        os.path.expanduser("~/HDOCKlite-v1.1"),
        "/usr/local/bin/HDOCKlite-v1.1"
    ]
    
    for p in potential_paths:
        if os.path.exists(os.path.join(p, "hdock")):
            return p
            
    return None
