# -*- coding: utf-8 -*-
"""
Benchmark Analysis for JCIM Manuscript
======================================

Automated analysis of 10 known structures (6 glues + 4 PROTACs)
to validate GlueTK's classification accuracy.

Usage in PyMOL:
    run validation/benchmark_analysis.py
    benchmark_all()  # Run all 10 structures
    
Or individual:
    analyze_cc885()  # Analyze CC-885 (6H0F)
"""

from __future__ import print_function
import os
import csv
from pymol import cmd

# Ensure commands are available
try:
    from gluetk.ppi_analyzer import analyze_protein_protein_interface, identify_neo_epitope
    # Note: binding_score module has been removed
    # from gluetk.binding_score import calculate_ternary_score, format_score_report
    from gluetk.interaction_analyzer import analyze_ternary_complex
except ImportError:
    print("[Benchmark] Import from package failed, trying direct import...")
    import sys
    plugin_dir = os.path.join(os.path.dirname(__file__), '..', 'gluetk')
    sys.path.insert(0, plugin_dir)
    from ppi_analyzer import analyze_protein_protein_interface, identify_neo_epitope
    # from binding_score import calculate_ternary_score, format_score_report
    from interaction_analyzer import analyze_ternary_complex

# ========== Benchmark Dataset ==========
BENCHMARK_STRUCTURES = {
    # Molecular Glues (Expected: Glue)
    "CC-885_6H0F": {
        "pdb": "6H0F",
        "glue_name": "CC885",  # Check actual residue name in PDB
        "e3_chain": "A",
        "substrate_chain": "B",
        "expected_mechanism": "Molecular Glue",
        "reference": "Matyskiela et al. (2018) Nature",
        "notes": "CRBN-CC-885-GSPT1, canonical glue"
    },
    "CC-90009_6BOY": {
        "pdb": "6BOY",
        "glue_name": "CC90009",
        "e3_chain": "A",
        "substrate_chain": "B",
        "expected_mechanism": "Molecular Glue",
        "reference": "Hansen et al. (2020)",
        "notes": "CRBN-CC90009-GSPT1"
    },
    "Lenalidomide_4TZ4": {
        "pdb": "4TZ4",
        "glue_name": "1N6",  # Lenalidomide residue name
        "e3_chain": "A",
        "substrate_chain": "B",
        "expected_mechanism": "Molecular Glue",
        "reference": "Kronke et al. (2014) Science",
        "notes": "CRBN-Lenalidomide-IKZF1, FDA approved"
    },
    "Thalidomide_4CI3": {
        "pdb": "4CI3",
        "glue_name": "3B9",  # Thalidomide analog
        "e3_chain": "A",
        "substrate_chain": "B",
        "expected_mechanism": "Molecular Glue",
        "reference": "Fischer et al. (2014)",
        "notes": "CRBN-Thalidomide-CK1α"
    },
    
    # PROTACs (Expected: PROTAC)
    "dBET1_6BN7": {
        "pdb": "6BN7",
        "glue_name": "QXQ",  # dBET1 residue name
        "e3_chain": "E",  # VHL
        "substrate_chain": "A",  # BRD4
        "expected_mechanism": "PROTAC",
        "reference": "Gadd et al. (2017) Nat Chem Biol",
        "notes": "BRD4-dBET1-VHL, PROTAC with linker"
    },
    "MZ1_6SIS": {
        "pdb": "6SIS",
        "glue_name": "P96",  # MZ1 residue name
        "e3_chain": "C",  # VHL
        "substrate_chain": "A",  # BRD4
        "expected_mechanism": "PROTAC",
        "reference": "Gadd et al. (2017)",
        "notes": "BRD4-MZ1-VHL, optimized PROTAC"
    },
}


def analyze_structure(name, info, output_dir="validation/results"):
    """
    Analyze a single structure from the benchmark dataset
    
    Returns:
        dict: {
            'name': str,
            'predicted_mechanism': str,
            'expected_mechanism': str,
            'correct': bool,
            'ppi_contacts': int,
            'bsa': float,
            'neo_epitope_count': int,
            'confidence': float,
            'cooperativity': float
        }
    """
    print("\n" + "="*80)
    print(f"Analyzing: {name}")
    print("="*80)
    
    pdb = info["pdb"]
    glue = info["glue_name"]
    e3_chain = info["e3_chain"]
    sub_chain = info["substrate_chain"]
    expected = info["expected_mechanism"]
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Fetch structure
    print(f"[1/4] Fetching {pdb}...")
    try:
        cmd.fetch(pdb, async_=0)
    except Exception as e:
        print(f"⚠️  Failed to fetch {pdb}: {e}")
        return None
    
    # Step 2: PPI Analysis
    print(f"[2/4] Analyzing protein-protein interface ({e3_chain} vs {sub_chain})...")
    ppi_csv = os.path.join(output_dir, f"{name}_ppi.csv")
    
    try:
        ppi_result = analyze_protein_protein_interface(
            obj_name=pdb,
            protein1_chains=[e3_chain],
            protein2_chains=[sub_chain],
            output_csv=ppi_csv
        )
    except Exception as e:
        print(f"⚠️  PPI analysis failed: {e}")
        ppi_result = None
    
    # Step 3: Neo-Epitope Detection
    print(f"[3/4] Detecting neo-epitope (Glue: {glue})...")
    neo_csv = os.path.join(output_dir, f"{name}_neo_epitope.csv")
    
    try:
        neo_result = identify_neo_epitope(
            obj_name=pdb,
            e3_ligase_chains=[e3_chain],
            substrate_chains=[sub_chain],
            glue_resname=glue,
            output_csv=neo_csv
        )
    except Exception as e:
        print(f"⚠️  Neo-epitope detection failed: {e}")
        neo_result = None
    
    # Step 4: Extract Results
    predicted_mechanism = "Unknown"
    ppi_contacts = 0
    bsa = None
    neo_epitope_count = 0
    confidence = 0.0
    cooperativity = 0.0
    
    if ppi_result:
        ppi_contacts = ppi_result.get('interface_contacts', 0)
        bsa = ppi_result.get('bsa')
        is_strong = ppi_result.get('is_strong_interface', False)
    
    if neo_result:
        neo_epitope_count = neo_result.get('neo_epitope_count', 0)
        confidence = neo_result.get('confidence', 0.0)
        predicted_mechanism = "Molecular Glue" if neo_result.get('is_molecular_glue', False) else "PROTAC"
    
    # Fallback classification if neo-epitope detection failed
    if predicted_mechanism == "Unknown" and ppi_result:
        if ppi_contacts >= 10 or (bsa and bsa > 800):
            predicted_mechanism = "Molecular Glue"
        elif ppi_contacts < 3:
            predicted_mechanism = "PROTAC"
    
    # Check correctness
    correct = (predicted_mechanism == expected)
    
    # Print Summary
    print("\n" + "-"*80)
    print(f"Results for {name}:")
    print("-"*80)
    print(f"  Expected:  {expected}")
    print(f"  Predicted: {predicted_mechanism}")
    print(f"  Status:    {'✅ CORRECT' if correct else '❌ WRONG'}")
    print(f"\n  PPI Contacts: {ppi_contacts}")
    if bsa:
        print(f"  BSA:          {bsa:.1f} Ų")
    print(f"  Neo-Epitope:  {neo_epitope_count} residues")
    print(f"  Confidence:   {confidence:.2f}")
    print("-"*80)
    
    return {
        'name': name,
        'pdb': pdb,
        'predicted_mechanism': predicted_mechanism,
        'expected_mechanism': expected,
        'correct': correct,
        'ppi_contacts': ppi_contacts,
        'bsa': bsa,
        'neo_epitope_count': neo_epitope_count,
        'confidence': confidence,
        'cooperativity': cooperativity,
        'reference': info.get('reference', ''),
        'notes': info.get('notes', '')
    }


def benchmark_all():
    """
    Run benchmark analysis on all structures
    """
    print("\n" + "="*80)
    print("GlueTK Benchmark Analysis - JCIM Manuscript")
    print("="*80)
    print(f"Total structures: {len(BENCHMARK_STRUCTURES)}")
    print("="*80)
    
    results = []
    
    for name, info in BENCHMARK_STRUCTURES.items():
        result = analyze_structure(name, info)
        if result:
            results.append(result)
    
    # ========== Summary Statistics ==========
    print("\n\n" + "="*80)
    print("BENCHMARK SUMMARY")
    print("="*80)
    
    total = len(results)
    correct = sum(1 for r in results if r['correct'])
    accuracy = (correct / total * 100) if total > 0 else 0
    
    # Calculate True Positives, True Negatives, etc.
    tp = sum(1 for r in results if r['predicted_mechanism'] == 'Molecular Glue' and r['expected_mechanism'] == 'Molecular Glue')
    tn = sum(1 for r in results if r['predicted_mechanism'] == 'PROTAC' and r['expected_mechanism'] == 'PROTAC')
    fp = sum(1 for r in results if r['predicted_mechanism'] == 'Molecular Glue' and r['expected_mechanism'] == 'PROTAC')
    fn = sum(1 for r in results if r['predicted_mechanism'] == 'PROTAC' and r['expected_mechanism'] == 'Molecular Glue')
    
    sensitivity = (tp / (tp + fn) * 100) if (tp + fn) > 0 else 0
    specificity = (tn / (tn + fp) * 100) if (tn + fp) > 0 else 0
    
    print(f"\nOverall Accuracy: {correct}/{total} ({accuracy:.1f}%)")
    print(f"\nConfusion Matrix:")
    print(f"  True Positives (Glue):  {tp}")
    print(f"  True Negatives (PROTAC): {tn}")
    print(f"  False Positives:        {fp}")
    print(f"  False Negatives:        {fn}")
    print(f"\nSensitivity: {sensitivity:.1f}%")
    print(f"Specificity: {specificity:.1f}%")
    
    # Detailed Results Table
    print("\n" + "-"*80)
    print("Detailed Results:")
    print("-"*80)
    print(f"{'Name':<25} {'Expected':<15} {'Predicted':<15} {'PPI':<5} {'BSA':<8} {'Neo':<5} {'Status':<10}")
    print("-"*80)
    
    for r in results:
        bsa_str = f"{r['bsa']:.0f}" if r['bsa'] else "N/A"
        status = "✅" if r['correct'] else "❌"
        print(f"{r['name']:<25} {r['expected_mechanism']:<15} {r['predicted_mechanism']:<15} "
              f"{r['ppi_contacts']:<5} {bsa_str:<8} {r['neo_epitope_count']:<5} {status:<10}")
    
    print("-"*80)
    
    # Save to CSV
    csv_path = "validation/results/benchmark_summary.csv"
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    
    with open(csv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Name', 'PDB', 'Expected', 'Predicted', 'Correct', 'PPI_Contacts', 
                        'BSA', 'Neo_Epitope_Count', 'Confidence', 'Reference', 'Notes'])
        for r in results:
            writer.writerow([
                r['name'], r['pdb'], r['expected_mechanism'], r['predicted_mechanism'],
                r['correct'], r['ppi_contacts'], r['bsa'] or '', r['neo_epitope_count'],
                r['confidence'], r['reference'], r['notes']
            ])
    
    print(f"\n✅ Results saved to: {csv_path}")
    print("="*80)
    
    return results


# ========== Individual Analysis Functions ==========
def analyze_cc885():
    """Analyze CC-885 (canonical molecular glue)"""
    return analyze_structure("CC-885_6H0F", BENCHMARK_STRUCTURES["CC-885_6H0F"])

def analyze_dbet1():
    """Analyze dBET1 (PROTAC negative control)"""
    return analyze_structure("dBET1_6BN7", BENCHMARK_STRUCTURES["dBET1_6BN7"])


# Register commands to PyMOL
cmd.extend("benchmark_all", benchmark_all)
cmd.extend("analyze_cc885", analyze_cc885)
cmd.extend("analyze_dbet1", analyze_dbet1)

print("[Benchmark] Commands registered:")
print("  benchmark_all()  - Run all 10 structures")
print("  analyze_cc885()  - Analyze CC-885 only")
print("  analyze_dbet1()  - Analyze dBET1 only")
