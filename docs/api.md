---
layout: default
title: API Reference
---

# API Reference

Complete API documentation for GLINT modules.

## Ternary Complex Evaluator

### `TernaryComplexEvaluator`

Main class for evaluating ternary complexes.

```python
from glint.ternary_complex_evaluator import TernaryComplexEvaluator

evaluator = TernaryComplexEvaluator()
features = evaluator.evaluate(
    pdb_path="complex.pdb",
    e3_chain="A",
    poi_chain="B",
    ligand_resn="UNL",
    ligand_smiles="CC(=O)Nc1ccc(O)cc1",
    obj_name="my_complex"
)
```

### Methods

#### `evaluate(pdb_path, e3_chain, poi_chain, ligand_resn, ligand_smiles=None, obj_name=None)`

Evaluate a ternary complex and return comprehensive features.

**Parameters:**
- `pdb_path` (str): Path to PDB file
- `e3_chain` (str): E3 ligase chain ID
- `poi_chain` (str): POI chain ID
- `ligand_resn` (str): Ligand residue name
- `ligand_smiles` (str, optional): Ligand SMILES string
- `obj_name` (str, optional): PyMOL object name

**Returns:**
- `TernaryComplexFeatures`: Comprehensive feature object

## Ligand Calculator

### `LigandCalculator`

Calculate molecular properties for small molecules.

```python
from glint.ternary_complex_evaluator import LigandCalculator

calc = LigandCalculator()
props = calc.calculate(smiles="CC(=O)Nc1ccc(O)cc1")
```

### Properties

- `molecular_weight` (float): Molecular weight in Daltons
- `logp` (float): LogP value
- `tpsa` (float): Topological polar surface area
- `rotatable_bonds` (int): Number of rotatable bonds
- `hbd_count` (int): Hydrogen bond donor count
- `hba_count` (int): Hydrogen bond acceptor count
- `fsp3` (float): Fraction of sp3 carbons
- `num_rings` (int): Number of rings

## BSA Calculator

### `BSACalculator`

Calculate buried surface area for interfaces.

```python
from glint.ternary_complex_evaluator import BSACalculator

bsa_calc = BSACalculator(obj_name="my_complex", ligand_resn="UNL")
bsa_results = bsa_calc.calculate_bsa_ternary(e3_chain="A", poi_chain="B")
```

## GUI Components

### TernaryEvaluationTab

Main GUI tab for ternary complex evaluation.

```python
from glint.gui.tabs.ternary_evaluation import TernaryEvaluationTab

tab = TernaryEvaluationTab(parent_window)
```

### Methods

- `run_full_evaluation()`: Run complete ternary complex evaluation
- `run_ligand_only()`: Calculate ligand properties only
- `export_results()`: Export results to CSV
- `visualize_geometry()`: Visualize geometric features in PyMOL

## Data Structures

### `TernaryComplexFeatures`

```python
@dataclass
class TernaryComplexFeatures:
    interface: InterfaceFeatures
    ligand: LigandFeatures
    distances: DistanceFeatures
    geometry: GeometryFeatures
    balance_index: Optional[float]
```

### `InterfaceFeatures`

```python
@dataclass
class InterfaceFeatures:
    bsa_total: float
    bsa_mg_e3: float
    bsa_mg_poi: float
    bsa_e3_poi: float
    contact_count_45: int
    contact_count_50: int
    min_inter_chain_dist: float
```
