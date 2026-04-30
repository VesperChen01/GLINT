---
layout: default
title: GLINT Tutorial
---

![](/GLINT/images/media/image1.png)

**Welcome to the GLINT tutorial!**

As molecular glues become increasingly important in targeted protein
degradation and the regulation of protein--protein interactions,
efficiently identifying potential targets, characterizing interface
features, and optimizing candidate molecules has emerged as a key
challenge in computational biology and drug design. GLINT was developed
in this context as an integrated analysis toolkit, aiming to provide
researchers with a systematic solution spanning from target discovery to
molecular optimization.

GLINT integrates several core functional modules, including:

- **Target Discovery**: for identifying potential binding sites and
  detecting G-motif features

- **Hit Identification**: for pocket detection and docking analysis

- **Ternary Evaluation**: for interface scoring and neo-epitope mapping

- ![](/GLINT/images/media/image2.png)**Lead Optimization**: for interaction
  profiling and mutation
  analysis![](/GLINT/images/media/image3.png)

> This tutorial starts with basic operations and gradually guides you
> through the workflow of using GLINT. Through step-by-step examples, it
> demonstrates how to complete a full molecular glue design and analysis
> task. Whether you are new to this field or looking to improve your
> analysis efficiency, this guide will help you get started quickly.

# **Chapter 1 Installation guide**

Welcome to the official GLINT installation guide. This chapter covers
the system requirements, pre-requisites, and step-by-step installation
procedures for macOS, Windows, and Linux.

## **System Requirements**

> Table 1. Minimum System Requirements

  -----------------------------------------------------------------------
               Component                          Requirement
  ----------------------------------- -----------------------------------
           Operating System           macOS 10.14+, Windows 10/11, Linux
                                                (Ubuntu 20.04+)

               Processor                       64-bit processor

                  RAM                   8 GB minimum (16 GB recommended

              Disk Space                        5 GB free space

               Software                 Miniconda or Anaconda, Python3
  -----------------------------------------------------------------------

## **1.2 Pre-requisites: Installing Conda**

![](/GLINT/images/media/image5.svg)

Figure 1. GLINT Installer interface with key steps highlighted

  -----------------------------------------------------------------------------------
              Platform                                    Path
  --------------------------------- -------------------------------------------------
            MacOS / Linux                        \~/.pymol/startup/glint

               Windows               C:\\Users\\\<Username\>\\.pymol\\startup\\glint
  -----------------------------------------------------------------------------------

  : Table 2. Installation Paths

## **1.3 Optional Dependencies**

**1.3.1 HADDOCK3 (Protein-Protein Docking)**

HADDOCK3 is required for ternary complex modeling. If the system has the
necessary build tools, it is typically installed automatically by the
installer or setup script. If not, it can be installed manually using
the command pip install haddock3. On Windows, it is recommended to use
the *HADDOCK Web Server* instead of local installation:
https://wenmr.science.uu.nl/haddock2.4/

**1.3.2 APBS (Electrostatic Analysis)**

APBS is used for calculating surface electrostatics. On macOS, it is
usually bundled with the GUI installer; otherwise, it can be installed
via Homebrew using brew install brewsci/bio/apbs. On Windows, it should
be downloaded from Poisson-Boltzmann.org.

**1.3.3 Missing Dependencies**

If any dependencies are missing, you can rerun the installer to attempt
automatic resolution. Alternatively, required packages can be installed
manually using Conda with the command conda install \<package_name\>.

# **Chapter 2 Target Discovery**

This chapter introduces the *Target Discovery* functionality of GLINT,
designed for molecular glue target identification. It includes three
core modules:

- G-Motif Detection

- Protein Surface Analysis

- Surface Similarity & Complementarity Analysis

![Figure 2.　 Target Discovery functionality of
GLINT.](/GLINT/images/media/image7.svg)

## **2.1 G-Motif Detection**

*G-Motif Detection* is designed to identify potential G-loop domains in
CRBN molecular glue substrate proteins. The G-loop is an 8-residue
structural motif whose conformation enables binding to the CRBN E3
ligase, forming a key interface in molecular glue--induced ternary
complexes.

To address the reliance of CRBN molecular glues on local loop
conformations of substrates, we developed a structure-based G-loop
template matching module. This module extracts Cα coordinates of
consecutive 8-residue segments from the target protein and aligns them
to known G-loop templates derived from validated substrates. Structural
similarity is evaluated *via* geometric superposition and RMSD
calculation. Compared to sequence motif searches, this approach
emphasizes structural equivalence, making it more suitable for
identifying candidates with low sequence homology but conserved
conformations.

The implementation supports multiple templates and alias mappings, and
is robust to insertion codes, alternate conformations, and
inconsistencies in PDB residue numbering for CRBN interface regions. A
default RMSD cutoff of 3.5 Å is applied, with an optional constraint
requiring glycine at position 6. The output is fully compatible with
downstream residue-highlighting tools and can be directly visualized in
PyMOL, enabling an integrated "search--locate--interpret" workflow.

**Core Algorithm**\
RMSD-based structural matching: align the Cα atoms of 8-residue segments
in the target protein to known G-loop templates and evaluate
conformational similarity using RMSD.

Table 3. Templates of G-Motif

  -----------------------------------------------------------------------
           Name                   PDB ID               Residue Range
  ----------------------- ----------------------- -----------------------
           GSPT1                   5HXB                 A:570--577

           CK1α                    5FQD                  C:35--42
  -----------------------------------------------------------------------

**2.1.1 Workflow**

**Load Structure:**\
Load the target protein structure in PyMOL (PDB file or existing
object).

**Open Target Discovery Panel:**\
Click **GLINT → Target Discovery** tab.

**Select Target Object:**

- Choose the object from the "Target Object" dropdown

- Click "refresh" to update the list

**Select Template:**

- **GSPT1 (5HXB):** for GSPT1-type G-loop detection

- **CK1α (5FQD):** for CK1α-type G-loop detection

- **From Selection:** use currently selected residues in PyMOL as a
  custom template

![Figure 3. G-Motif (CRBN G-loop) detection workflow and results in
GLINT.](/GLINT/images/media/image9.svg)

The top panel shows the Target Discovery interface with parameters set
for G-Motif detection (template: GSPT1, RMSD cutoff: 3.5 Å, glycine
requirement enabled). The middle panel displays the detection log in
PyMOL, listing identified hits with their chain IDs, residue ranges,
sequences, and RMSD values. The bottom panel visualizes one
representative hit (Hit 1), where the detected G-loop surface patch is
highlighted in yellow on the target protein structure (green cartoon),
indicating a potential molecular glue binding region.

## **2.2 Protein Surface Analysis**

Protein Surface Analysis identifies and classifies functional regions on
protein surfaces by characterizing their physicochemical properties.
Specifically, it detects electrostatic patches, including positively
charged regions (blue) and negatively charged regions (red), as well as
hydrophobic patches (green), which often correspond to potential
ligand-binding sites.

**2.2.1 workflow**

**Select Target Object:** Choose the object from the dropdown.

**Set Parameters:**

  -------------------------------------------------------------------------
   **Parameter**     **Default**                **Description**
  --------------- ----------------- ---------------------------------------
        pH               7.4          pH for APBS electrostatics (affects
                                              protonation states)

      Surface       Electrostatic                Coloring mode
     Property         Potential     

    Output CSV           --                Optional output file path
  -------------------------------------------------------------------------

Run Analysis:　Click Analyze Surface.**\
2.2.2 Color Scheme**

**Electrostatic Potential:**

🔴 Red: Negative (acidic)

⚪ White: Neutral

🔵 Blue: Positive (basic)

**Hydrophobicity**:

⚪ White: Hydrophilic

🟢 Green: Hydrophobic

![](/GLINT/images/media/image11.svg)

Figure 4. Protein surface hydrophobicity analysis using GLINT.

The top panel shows the Protein Surface Analysis interface with
hydrophobicity selected as the surface property (pH = 7.4). The middle
panel displays the PyMOL log output confirming successful surface
rendering using the Kyte--Doolittle scale. The bottom panel presents the
protein surface visualization, where hydrophobic regions are highlighted
in green and hydrophilic regions in white, revealing potential
ligand-binding sites and functional surface patches.

## **2.3 Surface Similarity & Complementarity**

This module provides a comprehensive framework for comparing geometric
and chemical features across two protein surfaces. It supports both
similarity and complementarity analyses. In similarity search, the
module identifies regions on a target surface that closely resemble a
given template, enabling the detection of structurally or chemically
conserved patterns. In complementarity search, it focuses on finding
regions that are structurally and chemically compatible with the
template surface, which is particularly useful for predicting
protein--protein interaction (PPI) interfaces.

**2.3.1 Workflow**

**Select Objects:**

- Object 1: template protein (with known functional surface)

- Object 2: target protein

**Define Template Region:**

  --------------- -------------------------------------------------------
    **Source**                        **Description**

   Ligand Pocket       Automatically constructs pocket around ligand
                                       (recommended)

   Custom Region                  Manually defined region
  --------------- -------------------------------------------------------

**Ligand Pocket Parameters:**

  -------------------------------------------------------------------------
  **Parameter**            **Default**   **Description**
  ------------------------ ------------- ----------------------------------
  Pocket Ligand            organic       Ligand selection expression

  Pocket Radius (Å)        6.0           Radius around ligand
  -------------------------------------------------------------------------

**Analysis Parameters:**

  ------------------------------------------------------------------------
  **Parameter**            **Default**   **Description**
  ------------------------ ------------- ---------------------------------
  Surface Method           auto          auto / open3d / edtsurf

  Patch Radius (Å)         12.0          Surface patch size

  Interface Dist (Å)       4.0           Contact distance threshold
  ------------------------------------------------------------------------

**Run Analysis:**

- Run Similarity Search

- Run Complementarity Search

**Visualize Results:　**Click Show Matched Regions

![Figure 5. Surface complementarity analysis using
GLINT.](/GLINT/images/media/image13.svg)

The top panel shows the Surface Complementarity interface with
parameters configured for ligand pocket--based analysis (Object 1: 5FQD,
Object 2: 8ARJ). The left panel illustrates the template binding pocket
automatically defined around the ligand (cyan surface) on the reference
protein (PDB ID: 5FQD, chain B). The right panel shows matched surface
patches identified on the target protein (PDB ID: 8ARJ), highlighting
structurally similar regions. The bottom panel presents electrostatic
surface representations of the matched regions, where red and blue
indicate negative and positive potentials, respectively, demonstrating
physicochemical complementarity between the template and target
surfaces.

# **Chapter 3 Hit Identification**

**The Hit Identification module** is a core functional component of the
GLINT framework for candidate molecule modeling and screening. It
integrates two key tools---**AutoDock Vina** for small-molecule docking
and **HADDOCK3** for protein--protein docking---forming a continuous
workflow from ligand pose prediction to ternary complex modeling.

## **3.1 AutoDock Vina Small-Molecule Docking** Evaluate the binding affinity between candidate small-molecule ligands and target proteins, and predict binding poses.

**Workflow**

Receptor Selection:\
Load the target protein structure in PyMOL and select the object to be
used as the receptor.

Ligand Input:

- Single-ligand docking:\
  Select a single ligand file (supported formats: .pdbqt, .sdf, .mol2)

- Batch docking:\
  Select a folder containing multiple ligand files; the system will
  automatically iterate through all ligands and perform batch docking

Docking Box Definition:

- Center selection:\
  Automatically calculate the geometric center based on a PyMOL
  selection

- Manual input:\
  Directly specify the center coordinates (X, Y, Z) and box dimensions

Parameter Settings:

- exhaustiveness (search thoroughness):\
  Range: 1--32 (default: 8). Higher values increase search accuracy but
  require more computation time

- num_modes (number of output poses):\
  Range: 1--20. Controls the number of top-ranked binding poses returned

Run Docking:\
Click "Run Docking". The system will automatically perform format
conversion, docking calculations, and result parsing.

## **3.2 HADDOCK3 Protein--Protein Docking**

**Use Case:**\
Build structural models of E3 ligase--molecular glue--substrate ternary
complexes based on known ligand binding poses.

Workflow

Receptor and Ligand Selection:\
Select two protein objects in PyMOL as the receptor and ligand for
docking.

Residue Restraint Definition:

- Active residues:\
  Key residues involved in molecular glue--mediated interface formation

- Passive residues:\
  Residues allowed to participate in contacts but not driving the
  interaction

Restraint File Generation:\
The system automatically generates an AMBIG.tbl restraint file, using
either the *haddock3-restraints* CLI or a built-in fallback method.

Configuration Generation and Execution:\
A minimal HADDOCK3 configuration file is automatically generated,
including the following stages:

\[topoaa\] \# Topology building

\[rigidbody\] \# Rigid-body sampling (default: 60 structures)

\[flexref\] \# Flexible refinement

\[emref\] \# Energy minimization refinement

Model Collection and Loading:\
After computation, output models are collected based on priority (emref
\> flexref \> rigidbody), merged into a multi-model PDB file, and loaded
into PyMOL.

Scoring Extraction:\
The system automatically parses capri_ss.tsv or capri_clt.tsv files to
extract HADDOCK scores for downstream analysis.

Alternative Options:\
If HADDOCK3 is not installed locally:

- Click "Open HADDOCK Server" to access the HADDOCK web server

- Use "Export PDB for Server" to export structure files for online
  submission

## **3.3 Typical Workflow Example**

Ligand Library → Vina Batch Docking → Top-N Pose Selection → HADDOCK3
Ternary Complex Modeling → Candidate Optimization

Step-by-Step Description

1.  Batch Docking with Vina:\
    Perform batch docking on a candidate ligand library using AutoDock
    Vina.\
    Select hit molecules with binding affinity \< −7 kcal/mol (E3--MG
    complexes).

2.  Prepare Docking Inputs for HADDOCK3:\
    Import the selected E3--molecular glue (E3--MG) complexes together
    with the target protein structure into HADDOCK3.\
    Define residues surrounding the molecular glue as active restraints.

3.  Ternary Complex Modeling:\
    Perform restraint-driven docking to construct E3--MG--target protein
    ternary complexes, generating multiple candidate conformations for
    further evaluation.

> ![一張含有 文字, 螢幕擷取畫面, 收據 的圖片 AI
> 產生的內容可能不正確。](/GLINT/images/media/image14.png)

Figure 6. Hit Identification module in GLINT.\
The interface integrates two core components for candidate modeling and
screening. The top panel shows the AutoDock Vina workflow for
small-molecule docking, including receptor selection, ligand input,
docking box definition, and parameter settings (center, size,
exhaustiveness, and number of modes). The bottom panel presents the
HADDOCK3 protein--protein docking setup, where receptor and ligand
proteins are selected and active/passive residue restraints are defined.
Together, these tools enable a streamlined workflow from ligand docking
to ternary complex modeling.

# **Chapter 4 Ternary Evaluation**

The **Ternary Complex Evaluation** module is one of the core components
of GLINT, specifically designed to analyze interactions between E3
ubiquitin ligases and target proteins (POIs) mediated by molecular
glues. This module integrates three major computational
components---**Interface Analysis**, **Ligand Properties**, and
**Ternary Geometry**---to provide a comprehensive structural assessment
of ternary complexes.

Interface Module

The **Interface Module** focuses on calculating **Buried Surface Area
(BSA)** and contact atom statistics. A dedicated formulation for ternary
complexes is used:

$$\text{BSA}_{MG} = \left( SA_{E3 - MG} + SA_{POI - MG} \right) - SA_{E3 - POI - MG}$$

where SA*~E3-MG~* and SA*~POI-MG~* represent the surface areas of the
E3--ligand and POI--ligand binary complexes, respectively, and
SA*~E3-POI-MG~* corresponds to the total surface area of the full
ternary complex.

Surface area calculations are performed using PyMOL's get_area()
function, with **FreeSASA** as an alternative for ligand surface
calculations. In addition to BSA, the module computes:

- Contact atom counts within **4.5 Å and 5.0 Å** cutoffs

- Minimum interatomic distances between components

These metrics provide quantitative indicators for evaluating binding
strength and specificity of molecular glues.

Ligand Module

The **Ligand Module** evaluates physicochemical properties of
small-molecule ligands using the **RDKit** library. Based on the
ligand's SMILES representation, the module calculates:

- Molecular weight (MW)

- LogP

- Topological polar surface area (TPSA)

- Hydrogen bond donors and acceptors (HBD/HBA)

- Number of rotatable bonds

- Fsp³ ratio

- Ring count

These descriptors are essential for assessing drug-likeness and
understanding the ligand's role within the ternary complex.

If a SMILES string cannot be directly extracted, the system attempts to
reconstruct it from the ligand structure using **Open Babel** or RDKit.

The **Ternary Geometry Module** analyzes spatial relationships within
the complex. Two key geometric descriptors are computed:

- **COG Shift:**\
  The perpendicular distance of the ligand's center of geometry (COG)
  from the E3--POI axis, reflecting how far the ligand deviates from the
  interface axis

- **E3--MG--POI Angle:**\
  The angle formed by the three components, describing their relative
  spatial arrangement

Additionally, the module calculates:

- Centers of geometry for each component

- Minimum interatomic distances between components

These parameters provide detailed insights into the structural
organization of the ternary complex.

The module introduces a **Balance Index** (also referred to as the
**Duality Index**) to evaluate binding symmetry:

$$\text{Balance Index} = \frac{\min\left( BSA_{MG - E3},BSA_{MG - POI} \right)}{\max\left( BSA_{MG - E3},BSA_{MG - POI} \right)}$$

This metric ranges from **0 to 1**:

- Values close to **1.0** indicate balanced binding between E3 and POI

- Values close to **0.0** indicate asymmetric binding toward one side

This index is particularly useful for assessing the "molecular glue"
behavior, where ideal candidates exhibit balanced interactions with both
proteins.

The module features a modern, card-based interface, including:

- **Input panel:** selection of PyMOL object, E3 chain, POI chain, and
  ligand residue

- **Three analysis panels:** Interface, Ligand, and Geometry

- **Summary panel and detailed results view**

Users can:

- Run a full evaluation via **"Run Full Evaluation"**

- Execute individual modules separately

- Export results as CSV

- Visualize geometric features directly in PyMOL

The system performs computations using multithreading, with results
displayed upon completion.

![](/GLINT/images/media/image15.png)

Figure 6. **Ternary Evaluation module in GLINT.**

The interface provides integrated analysis of E3 ligase--molecular
glue--target protein ternary complexes. The top panel allows selection
of the PyMOL object, E3 chain, POI chain, and ligand residue. The lower
panels display three analysis components: **Interface** (buried surface
area and contact analysis), **Ligand** (molecular properties such as
molecular weight, LogP, TPSA, and hydrogen bonding features), and
**Geometry** (ternary complex metrics including distances, angles, and
shape complementarity). Results can be computed and visualized, enabling
quantitative evaluation of ternary complex stability and quality.

# **Chapter 5 Lead optimization**

The **Lead Optimization** module is a core component of the GLINT
toolkit, designed to provide medicinal chemists and structural
biologists with a comprehensive suite of molecular interaction analysis
tools. This module integrates four key analysis
approaches---**electrostatic complementarity**, **protein--protein
interface analysis**, **protein--ligand interaction analysis**, and
**ligand--ligand interaction analysis**---offering multi-perspective
guidance for structure-based drug design and optimization.

### **5.1 Electrostatic Complementarity (EC) Analysis**

The **Electrostatic Complementarity (EC)** module is one of the central
tools for lead optimization, used to evaluate the electrostatic matching
between a ligand and its protein binding site. This module is based on
electrostatic potentials computed using **APBS (Adaptive
Poisson--Boltzmann Solver)**.

The workflow follows a rigorous computational pipeline:

1.  **Structure preprocessing:**\
    Protein structures are processed using **PDB2PQR**, including
    hydrogen addition, protonation state optimization, and charge
    assignment

2.  **Electrostatic calculation:**\
    APBS is used to compute the 3D electrostatic potential field of the
    protein, generating OpenDX grid data

3.  **Ligand surface sampling:**\
    The ligand solvent-accessible surface is generated using van der
    Waals radii, followed by dense surface sampling

4.  **Complementarity calculation:**\
    At each sampled surface point, the EC index is computed as the
    product of protein and ligand electrostatic potentials

- **Positive EC values:** electrostatic complementarity (attractive
  interactions)

- **Negative EC values:** electrostatic repulsion

Advanced features include:

- σ-hole virtual site generation (for halogen bonding)

- Lone pair virtual sites (e.g., carbonyl oxygen)

- Bridging water molecule filtering

These enhancements improve both the accuracy and biological relevance of
EC analysis.

### **5.2 Protein--Protein Interface (PPI) Analysis**

The **Protein--Protein Interface (PPI)** module is designed to
characterize interactions within protein complexes. It identifies and
analyzes interface regions between two proteins, including:

- Interface residue identification

- Contact residue counts and types

- Interface surface area and interaction strength

Users can specify:

- Chain IDs (e.g., A and B)

- Distance cutoff (default: 4.5 Å)

The module supports multiple visualization modes:

- Surface + interaction representation

- Cartoon + interaction representation

Optional display features include residue labels, distance labels, and
hydrophobic interactions.

Interface strength is evaluated based on contact counts, interface area,
and interaction diversity, providing quantitative insight into complex
stability.

### **5.3 Protein--Ligand Interaction Analysis**

The **Protein--Ligand Interaction** module is one of the most widely
used tools for lead optimization. It provides detailed characterization
of interactions between small molecules and protein targets.

The system supports automatic ligand detection or manual specification
and identifies multiple interaction types:

- Hydrogen bonds

- Salt bridges

- π--π stacking

- Cation--π interactions

- Hydrophobic interactions

- Halogen bonds

- Metal coordination

Interactions are detected using geometry-based rules combined with
confidence scoring to ensure reliability.

Users can:

- Set a minimum confidence threshold (default: 0.8)

- Toggle hydrophobic interactions and distance labels

Visualization options include:

- **3D visualization** in PyMOL

- **2D interaction diagrams**, generated with ligand structures and
  annotated interactions for reporting and presentation

### **5.4 Ligand--Ligand Interaction Analysis**

The **Ligand--Ligand Interaction** module analyzes interactions between
two small molecules. It is particularly useful for:

- Drug--drug interaction analysis

- Covalent inhibitor reaction modeling

- Ligand dimerization studies

Users define ligand selections using PyMOL syntax (e.g., resn LIG1 and
resn LIG2). The system then:

- Counts contacting atoms

- Measures interatomic distances

- Identifies potential hydrogen bonds and π--π interactions

In molecular glue--based degrader design, this module helps evaluate
interaction balance between:

- Molecular glue and E3 ligase

- Molecular glue and target protein
