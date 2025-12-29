# -*- coding: utf-8 -*-
"""
Ligand Feature Extraction Module for Surface Analysis
======================================================
MaSIF-neosurf inspired small molecule feature extraction.

Provides:
- Atom-type based hydrophobicity for small molecules
- Partial charge estimation for ligand atoms
- Ligand detection and classification
"""

from typing import Dict, List, Optional, Tuple
import numpy as np

# =============================================================================
# Atom-Type Based Hydrophobicity (Wildman-Crippen LogP contributions)
# =============================================================================

# Simplified atom-type hydrophobicity based on Wildman-Crippen LogP
# Normalized to [0, 1] range for consistency with protein scale
ATOM_HYDROPHOBICITY = {
    # Carbon types
    'C.3': 0.65,    # sp3 carbon (alkyl)
    'C.2': 0.55,    # sp2 carbon (alkene, carbonyl)
    'C.ar': 0.60,   # aromatic carbon
    'C.1': 0.50,    # sp carbon (alkyne)
    'C.cat': 0.40,  # carbocation
    
    # Nitrogen types
    'N.3': 0.25,    # sp3 nitrogen (amine)
    'N.2': 0.30,    # sp2 nitrogen (imine)
    'N.ar': 0.35,   # aromatic nitrogen
    'N.am': 0.20,   # amide nitrogen
    'N.pl3': 0.25,  # planar nitrogen
    'N.4': 0.10,    # quaternary nitrogen (charged)
    
    # Oxygen types
    'O.3': 0.15,    # sp3 oxygen (ether, hydroxyl)
    'O.2': 0.10,    # sp2 oxygen (carbonyl)
    'O.co2': 0.05,  # carboxylate oxygen
    'O.spc': 0.10,  # water oxygen
    
    # Sulfur types
    'S.3': 0.70,    # sp3 sulfur (thiol, thioether)
    'S.2': 0.65,    # sp2 sulfur
    'S.O': 0.40,    # sulfoxide
    'S.O2': 0.30,   # sulfone
    
    # Phosphorus
    'P.3': 0.45,    # sp3 phosphorus
    
    # Halogens (very hydrophobic)
    'F': 0.75,
    'Cl': 0.85,
    'Br': 0.90,
    'I': 0.95,
    
    # Default by element
    'C': 0.60,
    'N': 0.25,
    'O': 0.12,
    'S': 0.68,
    'P': 0.45,
    'H': 0.50,
}

# =============================================================================
# Partial Charges for Ligand Atoms
# =============================================================================

# Gasteiger-like partial charges by atom type
ATOM_PARTIAL_CHARGES = {
    # Carbon types
    'C.3': 0.0,      # sp3 carbon
    'C.2': 0.1,      # sp2 carbon (slightly positive due to electronegativity)
    'C.ar': 0.0,     # aromatic carbon
    'C.1': 0.1,      # sp carbon
    'C.cat': 0.5,    # carbocation
    
    # Nitrogen types
    'N.3': -0.3,     # sp3 nitrogen (amine)
    'N.2': -0.2,     # sp2 nitrogen
    'N.ar': -0.15,   # aromatic nitrogen
    'N.am': -0.4,    # amide nitrogen
    'N.pl3': -0.3,   # planar nitrogen
    'N.4': 1.0,      # quaternary nitrogen (charged)
    
    # Oxygen types
    'O.3': -0.4,     # sp3 oxygen
    'O.2': -0.5,     # sp2 oxygen (carbonyl)
    'O.co2': -0.8,   # carboxylate oxygen (charged)
    
    # Sulfur types
    'S.3': -0.2,     # sp3 sulfur
    'S.2': -0.15,    # sp2 sulfur
    'S.O': 0.3,      # sulfoxide sulfur
    'S.O2': 0.5,     # sulfone sulfur
    
    # Phosphorus
    'P.3': 0.5,      # phosphate phosphorus
    
    # Halogens
    'F': -0.2,
    'Cl': -0.1,
    'Br': -0.05,
    'I': 0.0,
    
    # Default by element
    'C': 0.0,
    'N': -0.3,
    'O': -0.4,
    'S': -0.2,
    'P': 0.5,
    'H': 0.1,
}

# =============================================================================
# Common Ligand Codes (from MaSIF-neosurf)
# =============================================================================

COMMON_COFACTORS = {
    'ADP': 'Adenosine diphosphate',
    'ATP': 'Adenosine triphosphate',
    'COA': 'Coenzyme A',
    'FAD': 'Flavin adenine dinucleotide',
    'FMN': 'Flavin mononucleotide',
    'HEM': 'Heme',
    'NAD': 'Nicotinamide adenine dinucleotide',
    'NAP': 'NADP',
    'SAM': 'S-adenosyl methionine',
    'GDP': 'Guanosine diphosphate',
    'GTP': 'Guanosine triphosphate',
}

# Molecular glue related ligands
MOLECULAR_GLUE_LIGANDS = {
    'LEN': 'Lenalidomide',
    'POM': 'Pomalidomide',
    'THL': 'Thalidomide',
    'CC2': 'CC-220 (Iberdomide)',
    'DGX': 'Degronimid',
}

# =============================================================================
# Feature Extraction Functions
# =============================================================================

def is_ligand(resn: str) -> bool:
    """Check if residue is a ligand (not standard amino acid or water)."""
    standard_aa = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
        'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
    }
    water = {'HOH', 'WAT', 'H2O', 'DOD'}
    ions = {'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'CO', 'NI'}
    
    resn_upper = resn.upper().strip()
    return resn_upper not in standard_aa and resn_upper not in water and resn_upper not in ions


def guess_atom_type(atom_name: str, element: str, resn: str) -> str:
    """
    Guess SYBYL-like atom type from atom name and element.
    
    Args:
        atom_name: PDB atom name (e.g., 'CA', 'N1', 'O2')
        element: Element symbol (e.g., 'C', 'N', 'O')
        resn: Residue name
    
    Returns:
        SYBYL-like atom type string
    """
    elem = element.upper().strip()
    name = atom_name.upper().strip()
    
    # Halogens - return as-is
    if elem in ('F', 'CL', 'BR', 'I'):
        return elem
    
    # Carbon types
    if elem == 'C':
        # Aromatic carbons often have names like CA, CB in rings
        if resn in MOLECULAR_GLUE_LIGANDS or 'AR' in name:
            return 'C.ar'
        # Carbonyl carbon
        if name in ('C', 'C1', 'CO') or 'O' in name:
            return 'C.2'
        return 'C.3'  # Default sp3
    
    # Nitrogen types
    if elem == 'N':
        # Aromatic nitrogen
        if resn in MOLECULAR_GLUE_LIGANDS:
            return 'N.ar'
        # Amide nitrogen
        if 'AM' in name or name == 'N':
            return 'N.am'
        # Charged nitrogen (quaternary)
        if name.endswith('+') or 'NZ' in name:
            return 'N.4'
        return 'N.3'  # Default sp3
    
    # Oxygen types
    if elem == 'O':
        # Carboxylate
        if 'OXT' in name or name in ('OD1', 'OD2', 'OE1', 'OE2'):
            return 'O.co2'
        # Carbonyl
        if name in ('O', 'O1', 'O2') and 'C' not in name:
            return 'O.2'
        return 'O.3'  # Default sp3
    
    # Sulfur types
    if elem == 'S':
        return 'S.3'
    
    # Phosphorus
    if elem == 'P':
        return 'P.3'
    
    # Default: return element
    return elem


def get_ligand_hydrophobicity(atom_name: str, element: str, resn: str) -> float:
    """
    Get hydrophobicity value for a ligand atom.
    
    Args:
        atom_name: PDB atom name
        element: Element symbol
        resn: Residue name
    
    Returns:
        Hydrophobicity value [0, 1]
    """
    atom_type = guess_atom_type(atom_name, element, resn)
    
    # Try specific atom type first
    if atom_type in ATOM_HYDROPHOBICITY:
        return ATOM_HYDROPHOBICITY[atom_type]
    
    # Fall back to element
    elem = element.upper().strip()
    return ATOM_HYDROPHOBICITY.get(elem, 0.5)


def get_ligand_charge(atom_name: str, element: str, resn: str) -> float:
    """
    Get partial charge for a ligand atom.
    
    Args:
        atom_name: PDB atom name
        element: Element symbol
        resn: Residue name
    
    Returns:
        Partial charge value
    """
    atom_type = guess_atom_type(atom_name, element, resn)
    
    # Try specific atom type first
    if atom_type in ATOM_PARTIAL_CHARGES:
        return ATOM_PARTIAL_CHARGES[atom_type]
    
    # Fall back to element
    elem = element.upper().strip()
    return ATOM_PARTIAL_CHARGES.get(elem, 0.0)


def extract_ligand_features(atoms: List[Dict]) -> Dict[str, np.ndarray]:
    """
    Extract features for ligand atoms.
    
    Args:
        atoms: List of atom dicts with 'coord', 'element', 'resn', 'name' keys
    
    Returns:
        Dict with 'hydrophobicity', 'charges', 'coords', 'is_ligand' arrays
    """
    n_atoms = len(atoms)
    
    hydrophobicity = np.zeros(n_atoms)
    charges = np.zeros(n_atoms)
    is_lig = np.zeros(n_atoms, dtype=bool)
    coords = np.zeros((n_atoms, 3))
    
    for i, atom in enumerate(atoms):
        resn = atom.get('resn', '')
        name = atom.get('name', '')
        elem = atom.get('element', 'C')
        coord = atom.get('coord', [0, 0, 0])
        
        coords[i] = coord
        is_lig[i] = is_ligand(resn)
        
        if is_lig[i]:
            # Use ligand-specific features
            hydrophobicity[i] = get_ligand_hydrophobicity(name, elem, resn)
            charges[i] = get_ligand_charge(name, elem, resn)
        else:
            # Use default values (will be overwritten by protein features)
            hydrophobicity[i] = 0.5
            charges[i] = 0.0
    
    return {
        'hydrophobicity': hydrophobicity,
        'charges': charges,
        'coords': coords,
        'is_ligand': is_lig
    }


def classify_ligand(resn: str) -> str:
    """
    Classify ligand type.
    
    Args:
        resn: Residue name (3-letter code)
    
    Returns:
        Classification string
    """
    resn_upper = resn.upper().strip()
    
    if resn_upper in COMMON_COFACTORS:
        return f"cofactor:{COMMON_COFACTORS[resn_upper]}"
    
    if resn_upper in MOLECULAR_GLUE_LIGANDS:
        return f"molecular_glue:{MOLECULAR_GLUE_LIGANDS[resn_upper]}"
    
    return "small_molecule"


# =============================================================================
# Integration with Surface Analysis
# =============================================================================

def enhance_surface_features_with_ligand(
    surface_points: List,
    atoms: List[Dict],
    ligand_features: Dict[str, np.ndarray],
    distance_threshold: float = 4.0
) -> None:
    """
    Enhance surface point features with ligand-specific properties.
    
    Modifies surface_points in place to update hydrophobicity and
    electrostatic values for points near ligand atoms.
    
    Args:
        surface_points: List of SurfacePoint objects
        atoms: List of atom dicts
        ligand_features: Output from extract_ligand_features()
        distance_threshold: Distance to consider ligand influence (Å)
    """
    from scipy.spatial import cKDTree
    
    # Build KD-tree for ligand atoms only
    ligand_mask = ligand_features['is_ligand']
    if not np.any(ligand_mask):
        return  # No ligands
    
    ligand_coords = ligand_features['coords'][ligand_mask]
    ligand_hydro = ligand_features['hydrophobicity'][ligand_mask]
    ligand_charges = ligand_features['charges'][ligand_mask]
    
    ligand_tree = cKDTree(ligand_coords)
    
    # Update surface points near ligands
    for point in surface_points:
        # Find nearby ligand atoms
        distances, indices = ligand_tree.query(
            point.coord, k=min(5, len(ligand_coords)), distance_upper_bound=distance_threshold
        )
        
        # Filter valid results
        valid = distances < distance_threshold
        if not np.any(valid):
            continue
        
        valid_distances = distances[valid]
        valid_indices = indices[valid]
        
        # Distance-weighted average
        weights = 1.0 / (valid_distances + 0.1)
        weight_sum = np.sum(weights)
        
        # Update hydrophobicity (blend with existing)
        ligand_hydro_avg = np.sum(ligand_hydro[valid_indices] * weights) / weight_sum
        point.hydrophobicity = 0.5 * point.hydrophobicity + 0.5 * ligand_hydro_avg
        
        # Update electrostatic potential (add ligand contribution)
        ligand_charge_contribution = np.sum(ligand_charges[valid_indices] * weights) / weight_sum
        point.electrostatic_potential += ligand_charge_contribution * 5.0  # Scale factor