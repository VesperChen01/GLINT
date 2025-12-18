"""
Molecular Glue Database Module for GlueTK

This module provides a comprehensive database of known molecular glues,
their targets, E3 ligases, and structural information.

Author: Roufen Chen
Email: 12319021@zju.edu.cn
"""

import json
import os
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, asdict
from enum import Enum

try:
    from pymol import cmd
except ImportError:
    cmd = None


class E3LigaseType(Enum):
    """E3 ligase types used in molecular glue degradation."""
    CRBN = "CRBN"
    DCAF15 = "DCAF15"
    DCAF16 = "DCAF16"
    CDK_CYCLIN = "CDK/Cyclin"
    VHL = "VHL"
    MDM2 = "MDM2"
    IAP = "IAP"
    UNKNOWN = "Unknown"


class MechanismType(Enum):
    """Degradation mechanism types."""
    MOLECULAR_GLUE = "Molecular Glue"
    PROTAC = "PROTAC"
    HYBRID = "Hybrid"
    UNKNOWN = "Unknown"


@dataclass
class MolecularGlue:
    """Data class representing a molecular glue compound."""
    name: str
    aliases: List[str]
    smiles: str
    molecular_weight: float
    e3_ligase: str
    substrates: List[str]
    pdb_ids: List[str]
    clinical_status: str
    company: str
    mechanism_notes: str
    references: List[str]
    
    def to_dict(self) -> Dict:
        return asdict(self)


@dataclass
class Substrate:
    """Data class representing a molecular glue substrate."""
    name: str
    uniprot_id: str
    gene_name: str
    e3_ligase: str
    glues: List[str]
    recognition_motif: str
    key_residues: List[str]
    pdb_ids: List[str]
    disease_relevance: List[str]
    
    def to_dict(self) -> Dict:
        return asdict(self)


# ============================================================================
# MOLECULAR GLUE DATABASE
# ============================================================================

MOLECULAR_GLUES_DB = {
    # CRBN-based IMiDs
    'Thalidomide': MolecularGlue(
        name='Thalidomide',
        aliases=['Thalomid', 'α-(N-phthalimido)glutarimide'],
        smiles='O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1',
        molecular_weight=258.23,
        e3_ligase='CRBN',
        substrates=['IKZF1', 'IKZF3', 'SALL4', 'CK1α'],
        pdb_ids=['4CI1', '4CI2', '4CI3'],
        clinical_status='FDA Approved (Multiple Myeloma)',
        company='Celgene/BMS',
        mechanism_notes='First-generation IMiD, binds CRBN and recruits zinc finger proteins',
        references=['Fischer et al. Nature 2014', 'Krönke et al. Science 2014'],
    ),
    
    'Lenalidomide': MolecularGlue(
        name='Lenalidomide',
        aliases=['Revlimid', 'CC-5013'],
        smiles='Nc1cccc2C(=O)N(C3CCC(=O)NC3=O)C(=O)c12',
        molecular_weight=259.26,
        e3_ligase='CRBN',
        substrates=['IKZF1', 'IKZF3', 'CK1α', 'ZFP91'],
        pdb_ids=['4TZ4', '5FQD', '6H0F'],
        clinical_status='FDA Approved (Multiple Myeloma, MDS)',
        company='Celgene/BMS',
        mechanism_notes='Second-generation IMiD with improved potency',
        references=['Krönke et al. Science 2014', 'Petzold et al. Nature 2016'],
    ),
    
    'Pomalidomide': MolecularGlue(
        name='Pomalidomide',
        aliases=['Pomalyst', 'CC-4047'],
        smiles='Nc1ccc2C(=O)N(C3CCC(=O)NC3=O)C(=O)c2c1',
        molecular_weight=273.24,
        e3_ligase='CRBN',
        substrates=['IKZF1', 'IKZF3', 'SALL4', 'PLZF'],
        pdb_ids=['5FQD'],
        clinical_status='FDA Approved (Multiple Myeloma)',
        company='Celgene/BMS',
        mechanism_notes='Third-generation IMiD',
        references=['Matyskiela et al. Nature 2016'],
    ),
    
    'CC-885': MolecularGlue(
        name='CC-885',
        aliases=['CC885'],
        smiles='',  # Complex structure
        molecular_weight=0,
        e3_ligase='CRBN',
        substrates=['GSPT1', 'ZFP91'],
        pdb_ids=['6H0F', '6H0G'],
        clinical_status='Preclinical',
        company='Celgene/BMS',
        mechanism_notes='Selective GSPT1 degrader via G-loop recognition',
        references=['Matyskiela et al. Nature 2018'],
    ),
    
    'CC-90009': MolecularGlue(
        name='CC-90009',
        aliases=['CC90009'],
        smiles='',
        molecular_weight=0,
        e3_ligase='CRBN',
        substrates=['GSPT1'],
        pdb_ids=['6BOY'],
        clinical_status='Phase 1 (AML)',
        company='Celgene/BMS',
        mechanism_notes='Clinical GSPT1 degrader',
        references=['Hansen et al. 2020'],
    ),
    
    'CC-92480': MolecularGlue(
        name='CC-92480',
        aliases=['Mezigdomide', 'CC92480'],
        smiles='',
        molecular_weight=0,
        e3_ligase='CRBN',
        substrates=['IKZF1', 'IKZF3'],
        pdb_ids=['7ZNA'],
        clinical_status='Phase 3 (Multiple Myeloma)',
        company='BMS',
        mechanism_notes='Next-generation CELMoD with enhanced IKZF degradation',
        references=[''],
    ),
    
    'Iberdomide': MolecularGlue(
        name='Iberdomide',
        aliases=['CC-220'],
        smiles='',
        molecular_weight=0,
        e3_ligase='CRBN',
        substrates=['IKZF1', 'IKZF3', 'Aiolos'],
        pdb_ids=[],
        clinical_status='Phase 3 (SLE, Multiple Myeloma)',
        company='BMS',
        mechanism_notes='CELMoD for autoimmune diseases',
        references=[''],
    ),
    
    # DCAF15-based sulfonamides
    'Indisulam': MolecularGlue(
        name='Indisulam',
        aliases=['E7070'],
        smiles='Cc1ccc(NS(=O)(=O)c2ccc3[nH]c(=O)cc(Cl)c3c2)cc1',
        molecular_weight=396.85,
        e3_ligase='DCAF15',
        substrates=['RBM39', 'RBM23'],
        pdb_ids=['6SJ7'],
        clinical_status='Phase 2 (Solid Tumors)',
        company='Eisai',
        mechanism_notes='Aryl sulfonamide that recruits RRM domain proteins',
        references=['Han et al. Science 2017', 'Uehara et al. Nat Chem Biol 2017'],
    ),
    
    'E7820': MolecularGlue(
        name='E7820',
        aliases=[],
        smiles='',
        molecular_weight=0,
        e3_ligase='DCAF15',
        substrates=['RBM39'],
        pdb_ids=['6UE5'],
        clinical_status='Preclinical',
        company='Eisai',
        mechanism_notes='Sulfonamide molecular glue',
        references=['Faust et al. Nat Chem Biol 2020'],
    ),
    
    'Tasisulam': MolecularGlue(
        name='Tasisulam',
        aliases=['LY573636'],
        smiles='',
        molecular_weight=0,
        e3_ligase='DCAF15',
        substrates=['RBM39'],
        pdb_ids=[],
        clinical_status='Phase 2 (Discontinued)',
        company='Eli Lilly',
        mechanism_notes='Sulfonamide molecular glue',
        references=[''],
    ),
    
    # CDK-based glues
    'CR8': MolecularGlue(
        name='CR8',
        aliases=[],
        smiles='',
        molecular_weight=0,
        e3_ligase='CDK/Cyclin',
        substrates=['Cyclin K'],
        pdb_ids=['6TD3'],
        clinical_status='Tool Compound',
        company='Academic',
        mechanism_notes='CDK12/13 inhibitor that induces Cyclin K degradation',
        references=['Słabicki et al. Nature 2020'],
    ),
    
    # Other glues
    'BI-3802': MolecularGlue(
        name='BI-3802',
        aliases=[],
        smiles='',
        molecular_weight=0,
        e3_ligase='CRBN',
        substrates=['BCL6'],
        pdb_ids=[],
        clinical_status='Preclinical',
        company='Boehringer Ingelheim',
        mechanism_notes='BCL6 degrader via polymerization',
        references=['Kerres et al. Cell Rep 2017'],
    ),
    
    'NVP-DKY709': MolecularGlue(
        name='NVP-DKY709',
        aliases=['DKY709'],
        smiles='',
        molecular_weight=0,
        e3_ligase='CRBN',
        substrates=['IKZF2'],
        pdb_ids=[],
        clinical_status='Phase 1',
        company='Novartis',
        mechanism_notes='Selective IKZF2 (Helios) degrader',
        references=[''],
    ),
    
    'MRT-2359': MolecularGlue(
        name='MRT-2359',
        aliases=[],
        smiles='',
        molecular_weight=0,
        e3_ligase='CRBN',
        substrates=['GSPT1'],
        pdb_ids=[],
        clinical_status='Phase 1/2',
        company='Monte Rosa',
        mechanism_notes='GSPT1 degrader for MYC-driven cancers',
        references=[''],
    ),
}


# ============================================================================
# SUBSTRATE DATABASE
# ============================================================================

SUBSTRATES_DB = {
    'IKZF1': Substrate(
        name='IKZF1',
        uniprot_id='Q13422',
        gene_name='IKZF1',
        e3_ligase='CRBN',
        glues=['Lenalidomide', 'Pomalidomide', 'Thalidomide', 'CC-92480', 'Iberdomide'],
        recognition_motif='C2H2 Zinc Finger with G-loop',
        key_residues=['G151', 'G152', 'Q147'],
        pdb_ids=['4TZ4', '5FQD'],
        disease_relevance=['Multiple Myeloma', 'ALL', 'Autoimmune'],
    ),
    
    'IKZF3': Substrate(
        name='IKZF3',
        uniprot_id='Q9UKT9',
        gene_name='IKZF3',
        e3_ligase='CRBN',
        glues=['Lenalidomide', 'Pomalidomide', 'Thalidomide', 'CC-92480'],
        recognition_motif='C2H2 Zinc Finger with G-loop',
        key_residues=['G146', 'G147'],
        pdb_ids=['6H0F'],
        disease_relevance=['Multiple Myeloma', 'CLL'],
    ),
    
    'CK1α': Substrate(
        name='CK1α',
        uniprot_id='P48729',
        gene_name='CSNK1A1',
        e3_ligase='CRBN',
        glues=['Lenalidomide', 'Thalidomide'],
        recognition_motif='β-hairpin G-loop',
        key_residues=['G40'],
        pdb_ids=['4CI3', '5FQD'],
        disease_relevance=['MDS with del(5q)', 'AML'],
    ),
    
    'GSPT1': Substrate(
        name='GSPT1',
        uniprot_id='P15170',
        gene_name='GSPT1',
        e3_ligase='CRBN',
        glues=['CC-885', 'CC-90009', 'MRT-2359'],
        recognition_motif='β-hairpin G-loop',
        key_residues=['G575', 'G577'],
        pdb_ids=['6H0F', '6H0G', '6BOY'],
        disease_relevance=['AML', 'MYC-driven cancers'],
    ),
    
    'SALL4': Substrate(
        name='SALL4',
        uniprot_id='Q9UJQ4',
        gene_name='SALL4',
        e3_ligase='CRBN',
        glues=['Thalidomide', 'Pomalidomide'],
        recognition_motif='C2H2 Zinc Finger',
        key_residues=['G416', 'G432'],
        pdb_ids=['7SQE'],
        disease_relevance=['Teratogenicity', 'AML'],
    ),
    
    'RBM39': Substrate(
        name='RBM39',
        uniprot_id='Q14498',
        gene_name='RBM39',
        e3_ligase='DCAF15',
        glues=['Indisulam', 'E7820', 'Tasisulam'],
        recognition_motif='RRM domain',
        key_residues=['G268', 'R267'],
        pdb_ids=['6SJ7', '6UE5'],
        disease_relevance=['AML', 'Solid Tumors'],
    ),
    
    'Cyclin K': Substrate(
        name='Cyclin K',
        uniprot_id='O75461',
        gene_name='CCNK',
        e3_ligase='CDK/Cyclin',
        glues=['CR8'],
        recognition_motif='CDK interface',
        key_residues=[],
        pdb_ids=['6TD3'],
        disease_relevance=['Cancer'],
    ),
    
    'BCL6': Substrate(
        name='BCL6',
        uniprot_id='P41182',
        gene_name='BCL6',
        e3_ligase='CRBN',
        glues=['BI-3802'],
        recognition_motif='BTB domain',
        key_residues=[],
        pdb_ids=[],
        disease_relevance=['DLBCL', 'Lymphoma'],
    ),
}


# ============================================================================
# E3 LIGASE DATABASE
# ============================================================================

E3_LIGASES_DB = {
    'CRBN': {
        'full_name': 'Cereblon',
        'uniprot_id': 'Q96SW2',
        'complex': 'CRL4-CRBN',
        'recognition_features': ['β-hairpin G-loop', 'C2H2 Zinc Finger'],
        'binding_site_residues': ['W380', 'W386', 'W400', 'H378'],
        'pdb_examples': ['4CI1', '4TZ4', '6H0F'],
        'known_substrates': 50,
        'druggability': 'High',
        'notes': 'Most validated E3 for molecular glues, IMiD binding site well characterized',
    },
    
    'DCAF15': {
        'full_name': 'DDB1 and CUL4 Associated Factor 15',
        'uniprot_id': 'Q9NZQ3',
        'complex': 'CRL4-DCAF15',
        'recognition_features': ['RRM domain'],
        'binding_site_residues': [],
        'pdb_examples': ['6SJ7', '6UE5'],
        'known_substrates': 5,
        'druggability': 'Medium',
        'notes': 'Sulfonamide binding site, recruits RNA-binding proteins',
    },
    
    'DCAF16': {
        'full_name': 'DDB1 and CUL4 Associated Factor 16',
        'uniprot_id': 'Q9NXF7',
        'complex': 'CRL4-DCAF16',
        'recognition_features': ['Unknown'],
        'binding_site_residues': [],
        'pdb_examples': [],
        'known_substrates': 2,
        'druggability': 'Low',
        'notes': 'Emerging E3 for molecular glues',
    },
    
    'VHL': {
        'full_name': 'Von Hippel-Lindau',
        'uniprot_id': 'P40337',
        'complex': 'CRL2-VHL',
        'recognition_features': ['Hydroxyproline recognition'],
        'binding_site_residues': ['H110', 'S111', 'H115', 'W117'],
        'pdb_examples': ['1LM8', '5T35'],
        'known_substrates': 10,
        'druggability': 'High',
        'notes': 'Primarily used for PROTACs, limited molecular glue applications',
    },
}


# ============================================================================
# DATABASE CLASS
# ============================================================================

class GlueDatabase:
    """
    Database interface for molecular glue information.
    """
    
    def __init__(self):
        self.glues = MOLECULAR_GLUES_DB
        self.substrates = SUBSTRATES_DB
        self.e3_ligases = E3_LIGASES_DB
    
    def search_glue(self, query: str) -> List[MolecularGlue]:
        """Search for molecular glues by name or alias."""
        results = []
        query_lower = query.lower()
        
        for name, glue in self.glues.items():
            if query_lower in name.lower():
                results.append(glue)
            elif any(query_lower in alias.lower() for alias in glue.aliases):
                results.append(glue)
        
        return results
    
    def search_by_substrate(self, substrate: str) -> List[MolecularGlue]:
        """Find molecular glues that degrade a specific substrate."""
        results = []
        substrate_lower = substrate.lower()
        
        for name, glue in self.glues.items():
            if any(substrate_lower in s.lower() for s in glue.substrates):
                results.append(glue)
        
        return results
    
    def search_by_e3(self, e3_type: str) -> List[MolecularGlue]:
        """Find molecular glues that use a specific E3 ligase."""
        results = []
        
        for name, glue in self.glues.items():
            if e3_type.upper() in glue.e3_ligase.upper():
                results.append(glue)
        
        return results
    
    def get_substrate_info(self, substrate_name: str) -> Optional[Substrate]:
        """Get detailed information about a substrate."""
        return self.substrates.get(substrate_name)
    
    def get_e3_info(self, e3_name: str) -> Optional[Dict]:
        """Get detailed information about an E3 ligase."""
        return self.e3_ligases.get(e3_name)
    
    def get_pdb_structures(self, glue_name: str = None, 
                          substrate_name: str = None) -> List[str]:
        """Get PDB IDs for a glue or substrate."""
        pdb_ids = []
        
        if glue_name and glue_name in self.glues:
            pdb_ids.extend(self.glues[glue_name].pdb_ids)
        
        if substrate_name and substrate_name in self.substrates:
            pdb_ids.extend(self.substrates[substrate_name].pdb_ids)
        
        return list(set(pdb_ids))
    
    def get_clinical_glues(self, phase: str = None) -> List[MolecularGlue]:
        """Get molecular glues in clinical development."""
        results = []
        
        for name, glue in self.glues.items():
            status = glue.clinical_status.lower()
            if phase:
                if phase.lower() in status:
                    results.append(glue)
            else:
                if 'phase' in status or 'approved' in status:
                    results.append(glue)
        
        return results
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics."""
        stats = {
            'total_glues': len(self.glues),
            'total_substrates': len(self.substrates),
            'total_e3_ligases': len(self.e3_ligases),
            'glues_by_e3': {},
            'clinical_status': {
                'approved': 0,
                'phase_3': 0,
                'phase_2': 0,
                'phase_1': 0,
                'preclinical': 0,
            },
            'total_pdb_structures': 0,
        }
        
        pdb_set = set()
        
        for name, glue in self.glues.items():
            # Count by E3
            e3 = glue.e3_ligase
            stats['glues_by_e3'][e3] = stats['glues_by_e3'].get(e3, 0) + 1
            
            # Count by clinical status
            status = glue.clinical_status.lower()
            if 'approved' in status:
                stats['clinical_status']['approved'] += 1
            elif 'phase 3' in status:
                stats['clinical_status']['phase_3'] += 1
            elif 'phase 2' in status:
                stats['clinical_status']['phase_2'] += 1
            elif 'phase 1' in status:
                stats['clinical_status']['phase_1'] += 1
            else:
                stats['clinical_status']['preclinical'] += 1
            
            # Collect PDB IDs
            pdb_set.update(glue.pdb_ids)
        
        stats['total_pdb_structures'] = len(pdb_set)
        
        return stats
    
    def generate_report(self, glue_name: str = None) -> str:
        """Generate a formatted report."""
        report = []
        
        if glue_name:
            # Single glue report
            if glue_name not in self.glues:
                return f"Molecular glue '{glue_name}' not found in database."
            
            glue = self.glues[glue_name]
            report.append("=" * 60)
            report.append(f"MOLECULAR GLUE: {glue.name}")
            report.append("=" * 60)
            report.append(f"\nAliases: {', '.join(glue.aliases) if glue.aliases else 'None'}")
            report.append(f"E3 Ligase: {glue.e3_ligase}")
            report.append(f"Molecular Weight: {glue.molecular_weight if glue.molecular_weight else 'N/A'}")
            report.append(f"Clinical Status: {glue.clinical_status}")
            report.append(f"Company: {glue.company}")
            report.append(f"\nSubstrates: {', '.join(glue.substrates)}")
            report.append(f"PDB Structures: {', '.join(glue.pdb_ids) if glue.pdb_ids else 'None'}")
            report.append(f"\nMechanism: {glue.mechanism_notes}")
            report.append(f"\nReferences: {'; '.join(glue.references)}")
            
        else:
            # Database overview
            stats = self.get_statistics()
            
            report.append("=" * 60)
            report.append("MOLECULAR GLUE DATABASE OVERVIEW")
            report.append("=" * 60)
            report.append(f"\nTotal Molecular Glues: {stats['total_glues']}")
            report.append(f"Total Substrates: {stats['total_substrates']}")
            report.append(f"Total E3 Ligases: {stats['total_e3_ligases']}")
            report.append(f"Total PDB Structures: {stats['total_pdb_structures']}")
            
            report.append(f"\n{'─' * 40}")
            report.append("GLUES BY E3 LIGASE")
            report.append(f"{'─' * 40}")
            for e3, count in stats['glues_by_e3'].items():
                report.append(f"  {e3}: {count}")
            
            report.append(f"\n{'─' * 40}")
            report.append("CLINICAL STATUS")
            report.append(f"{'─' * 40}")
            report.append(f"  FDA Approved: {stats['clinical_status']['approved']}")
            report.append(f"  Phase 3: {stats['clinical_status']['phase_3']}")
            report.append(f"  Phase 2: {stats['clinical_status']['phase_2']}")
            report.append(f"  Phase 1: {stats['clinical_status']['phase_1']}")
            report.append(f"  Preclinical: {stats['clinical_status']['preclinical']}")
            
            report.append(f"\n{'─' * 40}")
            report.append("ALL MOLECULAR GLUES")
            report.append(f"{'─' * 40}")
            for name, glue in self.glues.items():
                report.append(f"  • {name} ({glue.e3_ligase}) - {glue.clinical_status}")
        
        report.append("\n" + "=" * 60)
        return '\n'.join(report)
    
    def export_to_json(self, output_path: str):
        """Export database to JSON file."""
        data = {
            'molecular_glues': {k: v.to_dict() for k, v in self.glues.items()},
            'substrates': {k: v.to_dict() for k, v in self.substrates.items()},
            'e3_ligases': self.e3_ligases,
            'statistics': self.get_statistics(),
        }
        
        with open(output_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"Database exported to: {output_path}")


# ============================================================================
# PYMOL COMMANDS
# ============================================================================

def glue_db_search(query: str, search_type: str = 'name'):
    """
    Search the molecular glue database.
    
    Usage:
        glue_db_search Lenalidomide
        glue_db_search IKZF1, search_type=substrate
        glue_db_search CRBN, search_type=e3
    
    Parameters:
    -----------
    query : str
        Search query
    search_type : str
        Type of search: 'name', 'substrate', or 'e3'
    """
    db = GlueDatabase()
    
    if search_type == 'name':
        results = db.search_glue(query)
    elif search_type == 'substrate':
        results = db.search_by_substrate(query)
    elif search_type == 'e3':
        results = db.search_by_e3(query)
    else:
        print(f"Unknown search type: {search_type}")
        return
    
    if results:
        print(f"\nFound {len(results)} result(s):\n")
        for glue in results:
            print(f"  • {glue.name}")
            print(f"    E3: {glue.e3_ligase}")
            print(f"    Substrates: {', '.join(glue.substrates)}")
            print(f"    Status: {glue.clinical_status}")
            print(f"    PDB: {', '.join(glue.pdb_ids) if glue.pdb_ids else 'None'}")
            print()
    else:
        print(f"No results found for '{query}'")


def glue_db_info(glue_name: str = None):
    """
    Display molecular glue database information.
    
    Usage:
        glue_db_info              # Show database overview
        glue_db_info Lenalidomide # Show specific glue info
    """
    db = GlueDatabase()
    report = db.generate_report(glue_name)
    print(report)


def glue_db_fetch(glue_name: str):
    """
    Fetch PDB structures for a molecular glue.
    
    Usage:
        glue_db_fetch CC-885
    """
    if cmd is None:
        print("PyMOL not available")
        return
    
    db = GlueDatabase()
    pdb_ids = db.get_pdb_structures(glue_name=glue_name)
    
    if not pdb_ids:
        print(f"No PDB structures found for '{glue_name}'")
        return
    
    print(f"Fetching {len(pdb_ids)} structure(s) for {glue_name}...")
    
    for pdb_id in pdb_ids:
        try:
            cmd.fetch(pdb_id)
            print(f"  ✓ Loaded {pdb_id}")
        except Exception as e:
            print(f"  ✗ Failed to load {pdb_id}: {e}")


def glue_db_export(output_path: str = 'glue_database.json'):
    """
    Export the molecular glue database to JSON.
    
    Usage:
        glue_db_export
        glue_db_export my_database.json
    """
    db = GlueDatabase()
    db.export_to_json(output_path)


# Register PyMOL commands
if cmd is not None:
    cmd.extend('glue_db_search', glue_db_search)
    cmd.extend('glue_db_info', glue_db_info)
    cmd.extend('glue_db_fetch', glue_db_fetch)
    cmd.extend('glue_db_export', glue_db_export)