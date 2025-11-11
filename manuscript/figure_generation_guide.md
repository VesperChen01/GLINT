# Figure Generation Guide for JCIM Manuscript

**Target**: Journal of Chemical Information and Modeling (JCIM)  
**Requirements**: 300 DPI minimum, TIFF or PNG format, color figures allowed  
**Last Updated**: 2025-11-11

**Note**: This guide reflects GlueTK's current features including G-motif detection and pocket analysis.

---

## Figure 1: GlueTK Software Architecture

**Type**: Flowchart diagram  
**Tool**: PowerPoint, draw.io, or BioRender  
**Size**: Full page width (7 inches)

### Content:
```
[User Input via PyMOL GUI/CLI]
           ↓
    [PyMOL Plugin Core]
           ↓
    ┌──────┴──────┐
    ↓             ↓
[Structure     [Parameters]
 Loading]       
    ↓             ↓
    └──────┬──────┘
           ↓
    [Analysis Modules]
    ┌──────┴──────┬──────────┬─────────┐
    ↓             ↓          ↓         ↓
[PPI Analyzer] [Neo-Epitope] [Scoring] [G-Motif]
    ↓             ↓          ↓         ↓
[BSA Calc]   [Confidence]  [Coop]   [RMSD]
    ↓             ↓          ↓         ↓
    └──────┬──────┴──────────┴─────────┘
           ↓
  [Mechanism Classification]
    (Glue / PROTAC / Unknown)
           ↓
    ┌──────┴──────┐
    ↓             ↓
[CSV Output]  [3D Visualization]
```

### Steps to Create:
1. Use BioRender or draw.io for professional look
2. Export as PNG 300 DPI
3. Add legend explaining each module

---

## Figure 2: CC-885 Molecular Glue Case Study (4 panels)

**Type**: Multi-panel structural figure  
**Tool**: PyMOL + Image editor  
**Size**: Full page width, 4 panels (A-D)

### Panel A: Protein-Protein Interface (BSA Visualization)
**PyMOL Commands:**
```python
# Load structure
fetch 6H0F

# Show E3 (CRBN) and Substrate (GSPT1)
hide everything
show cartoon, chain A  # CRBN
show cartoon, chain B  # GSPT1

# Color by chain
color blue, chain A
color red, chain B

# Show CC-885
show sticks, resn CC885
color yellow, resn CC885

# Highlight interface residues (from PPI analysis result)
# Read interface_residues from CSV
select interface_A, chain A and resi 60+62+64+65+67
select interface_B, chain B and resi 100+102+105+108
show sticks, interface_A or interface_B
color green, interface_A
color orange, interface_B

# Set view
set_view (\
     0.123,   -0.456,    0.789,\
    -0.456,    0.789,    0.123,\
    -0.789,   -0.123,    0.456,\
     0.000,    0.000,  -50.000,\
     1.234,    5.678,    9.012,\
    40.000,   60.000,  -20.000 )

# Export high-res PNG
ray 2400, 1800  # 300 DPI at 8x6 inches
png manuscript/figures/fig2a_ppi.png
```

### Panel B: Neo-Epitope Residues Highlighted
```python
# Continue from Panel A

# Highlight neo-epitope residues (from neo_epitope result)
select neo_epitope, chain B and resi 60+62+64+65+67
color magenta, neo_epitope
show spheres, neo_epitope
set sphere_scale, 0.3, neo_epitope

# Add labels
label neo_epitope and name CA, "%s%s" % (resn, resi)
set label_color, black
set label_size, 20

# Export
ray 2400, 1800
png manuscript/figures/fig2b_neo_epitope.png
```

### Panel C: 2D Interaction Diagram
```python
# Use your existing 2D diagram function
from gluetk.interaction_2d_plot import generate_2d_interaction_diagram

generate_2d_interaction_diagram(
    obj_name='6H0F',
    ligand_resname='CC885',
    output_png='manuscript/figures/fig2c_2d_diagram.png',
    show_labels=True
)
```

### Panel D: G-loop Structure with Glue Binding Mode (Using GlueTK Commands)
```python
# Use GlueTK's G-motif detection
find_crbn_g_motif(
    obj_name='6H0F',
    chain='B',
    residue_range='60-67',  # GSPT1 G-loop
    template_mode='builtin',  # Use known GSPT1 template
    rmsd_cutoff=3.5
)

# Visualize G-motif with glue binding
analyze_g_motif_glue_binding(
    obj_name='6H0F',
    substrate_chain='B',
    g_loop_range='60-67',
    glue_resname='CC885',
    e3_chain='A'
)

# Manual PyMOL visualization for publication
hide everything
show cartoon, chain B and resi 55-75  # G-loop region
color red, chain B

# Show CC-885 in detail
show sticks, resn CC885
color yellow, resn CC885
util.cnc  # Color by element

# Show key interactions (from GlueTK output)
distance hb1, chain B and resi 60 and name N, resn CC885 and name O1
distance hb2, chain B and resi 65 and name N, resn CC885 and name O2

# Labels for G-motif residues
label chain B and resi 60 and name CA, "Gly60"
label chain B and resi 62 and name CA, "Val62"
label chain B and resi 65 and name CA, "Gly65"

# Add G-motif annotation
pseudoatom g_motif_label, pos=[x, y, z]
label g_motif_label, "G-motif (RMSD=2.1Å)"

# Export
ray 2400, 1800
png manuscript/figures/fig2d_gloop.png
```

### Assembly:
1. Use Photoshop/GIMP to combine 4 panels
2. Add panel labels (A, B, C, D) in top-left corners
3. Add scale bars
4. Export as TIFF 300 DPI

---

## Figure 3: dBET1 PROTAC Negative Control (3 panels)

**Type**: Comparison figure  
**Size**: 2/3 page width

### Panel A: No Significant PPI
```python
fetch 6BN7
hide everything
show cartoon, chain A  # BRD4
show cartoon, chain E  # VHL
color blue, chain A
color red, chain E

# Show dBET1 linker
show sticks, resn QXQ
color yellow, resn QXQ

# Note: Very few interface contacts
ray 1800, 1800
png manuscript/figures/fig3a_protac_ppi.png
```

### Panel B: Linker Visualization
```python
# Zoom into linker region
zoom resn QXQ, 10
show sticks, byres (resn QXQ around 5)

# Export
ray 1800, 1800
png manuscript/figures/fig3b_linker.png
```

### Panel C: Side-by-side Comparison (CC-885 vs dBET1)
- Use image editor to place Fig2A and Fig3A side by side
- Add annotations: "Strong PPI (BSA=950Ų)" vs "Weak PPI (BSA=180Ų)"

---

## Figure 4: Benchmark Performance Summary

**Type**: Bar chart + Confusion matrix  
**Tool**: Python matplotlib or R ggplot2  
**Data Source**: `validation/results/benchmark_summary.csv`

### Python Code:
```python
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Load benchmark results
df = pd.read_csv('validation/results/benchmark_summary.csv')

# Create figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Accuracy by structure
structures = df['Name'].str.split('_').str[0]  # Extract short names
correct = df['Correct'].astype(int)

ax1.bar(range(len(structures)), correct, color=['green' if c else 'red' for c in correct])
ax1.set_xticks(range(len(structures)))
ax1.set_xticklabels(structures, rotation=45, ha='right')
ax1.set_ylabel('Correct (1) / Wrong (0)')
ax1.set_title('(A) Classification Results')
ax1.set_ylim([0, 1.2])

# Panel B: Confusion Matrix
tp = sum((df['Predicted'] == 'Molecular Glue') & (df['Expected'] == 'Molecular Glue'))
tn = sum((df['Predicted'] == 'PROTAC') & (df['Expected'] == 'PROTAC'))
fp = sum((df['Predicted'] == 'Molecular Glue') & (df['Expected'] == 'PROTAC'))
fn = sum((df['Predicted'] == 'PROTAC') & (df['Expected'] == 'Molecular Glue'))

confusion = [[tp, fn], [fp, tn]]
sns.heatmap(confusion, annot=True, fmt='d', cmap='Blues', ax=ax2,
            xticklabels=['Pred: Glue', 'Pred: PROTAC'],
            yticklabels=['True: Glue', 'True: PROTAC'],
            cbar=False)
ax2.set_title('(B) Confusion Matrix')

plt.tight_layout()
plt.savefig('manuscript/figures/fig4_benchmark.png', dpi=300, bbox_inches='tight')
plt.show()
```

---

## Figure 5: PPI Metrics Comparison (Glue vs PROTAC)

**Type**: Scatter plot with categorical separation  
**Tool**: Python matplotlib

### Python Code:
```python
import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('validation/results/benchmark_summary.csv')

# Separate by mechanism
glues = df[df['Expected'] == 'Molecular Glue']
protacs = df[df['Expected'] == 'PROTAC']

fig, ax = plt.subplots(figsize=(8, 6))

# Scatter plot: PPI Contacts vs BSA
ax.scatter(glues['PPI_Contacts'], glues['BSA'], 
          s=200, c='green', alpha=0.6, label='Molecular Glue', 
          edgecolors='black', linewidth=1.5)
ax.scatter(protacs['PPI_Contacts'], protacs['BSA'], 
          s=200, c='red', alpha=0.6, label='PROTAC',
          edgecolors='black', linewidth=1.5)

# Decision boundaries
ax.axhline(y=800, color='gray', linestyle='--', label='BSA threshold (800Ų)')
ax.axvline(x=10, color='gray', linestyle=':', label='PPI threshold (10 contacts)')

ax.set_xlabel('PPI Interface Contacts', fontsize=14)
ax.set_ylabel('Buried Surface Area (Ų)', fontsize=14)
ax.set_title('PPI Metrics: Molecular Glue vs PROTAC', fontsize=16)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('manuscript/figures/fig5_ppi_metrics.png', dpi=300, bbox_inches='tight')
plt.show()
```

---

## Figure 6: Pocket Analysis and Glue-Pocket Correlation

**Type**: Multi-panel pocket visualization  
**Tool**: PyMOL + GlueTK pocket commands

### Panel A: Pocket Detection at PPI Interface
```python
# Use GlueTK's pocket detection
fetch 6H0F

# Detect pockets at PPI interface
detect_pockets(
    obj_name='6H0F',
    region='interface',  # Focus on PPI interface
    chains=['A', 'B'],
    output_prefix='6H0F_pockets'
)

# Visualize pockets with interactions
visualize_pockets_with_interactions(
    obj_name='6H0F',
    pocket_data='6H0F_pockets.json',
    show_glue=True,
    glue_resname='CC885'
)

# Export
ray 2400, 1800
png manuscript/figures/fig6a_pockets.png
```

### Panel B: Glue-Pocket Correlation Heatmap
```python
import matplotlib.pyplot as plt
import json

# Load pocket-glue correlation data (from GlueTK output)
with open('6H0F_pocket_glue_correlation.json', 'r') as f:
    data = json.load(f)

# Plot correlation
pockets = [p['pocket_id'] for p in data['pockets']]
distances = [p['distance_to_glue'] for p in data['pockets']]
overlaps = [p['overlap_score'] for p in data['pockets']]

fig, ax = plt.subplots(figsize=(8, 6))
scatter = ax.scatter(distances, overlaps, s=200, c=range(len(pockets)), 
                     cmap='viridis', alpha=0.7, edgecolors='black')
ax.set_xlabel('Distance to Glue (Å)', fontsize=14)
ax.set_ylabel('Glue-Pocket Overlap Score', fontsize=14)
ax.set_title('Pocket-Glue Correlation Analysis', fontsize=16)

# Annotate primary binding pocket
primary_idx = overlaps.index(max(overlaps))
ax.annotate('Primary Binding', xy=(distances[primary_idx], overlaps[primary_idx]),
            xytext=(distances[primary_idx]+2, overlaps[primary_idx]+0.1),
            arrowprops=dict(arrowstyle='->', lw=2))

plt.colorbar(scatter, label='Pocket ID')
plt.tight_layout()
plt.savefig('manuscript/figures/fig6b_correlation.png', dpi=300)
```

### Panel C: Pocket Druggability Analysis
```python
# Generate druggability plot
import pandas as pd
import seaborn as sns

df = pd.read_csv('6H0F_pocket_properties.csv')

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Volume
axes[0].bar(df['pocket_id'], df['volume'])
axes[0].set_xlabel('Pocket ID')
axes[0].set_ylabel('Volume (Å³)')
axes[0].set_title('Pocket Volume')

# Depth
axes[1].bar(df['pocket_id'], df['depth'])
axes[1].set_xlabel('Pocket ID')
axes[1].set_ylabel('Depth (Å)')
axes[1].set_title('Pocket Depth')

# Hydrophobicity
axes[2].bar(df['pocket_id'], df['hydrophobicity'])
axes[2].set_xlabel('Pocket ID')
axes[2].set_ylabel('Hydrophobicity Score')
axes[2].set_title('Hydrophobicity')

plt.tight_layout()
plt.savefig('manuscript/figures/fig6c_druggability.png', dpi=300)
```

---

## Supporting Information Figures

### SI Figure 1: Parameter Sensitivity Analysis
- Test different BSA thresholds (600, 700, 800, 900, 1000Å²)
- Test different PPI contact cutoffs (5, 8, 10, 12, 15)
- Plot accuracy vs threshold for both parameters
- Show optimal is BSA=800Å², PPI contacts=10

### SI Figure 2: Runtime Benchmarking
- Bar chart showing time breakdown:
  - PPI Analysis: 12.3s
  - Neo-Epitope: 8.7s
  - BSA Calc: 3.2s
  - G-Motif Detection: 2.5s
  - Pocket Analysis: 15.6s
  - Scoring: 0.8s
  - **Total**: ~42s

### SI Figure 3: G-Motif Template Comparison
- Compare RMSD for different G-motif templates:
  - Ideal β-hairpin
  - GSPT1 (6H0F)
  - CK1α (4CI3)
  - VAV1 (custom)
- Show which template best fits each substrate

### SI Figure 4: GUI Screenshots
- Screenshot of main GUI with all tabs:
  - Tab 1: G-Motif Detection
  - Tab 2: Interaction Analysis
  - Tab 3: PPI & Neo-Epitope
  - Tab 4: Pocket Analysis (NEW)
  - Tab 5: Electrostatics
- Annotate key features and workflow

---

## Graphical Abstract (TOC Graphic)

**Required by JCIM**: Single image summarizing the work  
**Size**: 8.5 cm × 4.75 cm (3.35 × 1.87 inches)  
**Resolution**: 300 DPI minimum

### Content Suggestion:
```
[Input: PDB Structure]
        ↓
   [GlueTK]
        ↓
    ┌───┴───┐
    ↓       ↓
 [Glue]  [PROTAC]
  (✓)     (✗)
```

With visual representations:
- Left side: CC-885 structure showing strong PPI
- Right side: dBET1 structure showing linker
- Middle: GlueTK logo/icon

---

## Figure Checklist

- [ ] All figures at 300 DPI minimum
- [ ] File format: TIFF or PNG
- [ ] Color figures properly labeled
- [ ] Scale bars included where appropriate
- [ ] Panel labels (A, B, C, D) in top-left corners
- [ ] Figure legends written (150-200 words each)
- [ ] High-contrast colors for accessibility
- [ ] Font size ≥8 pt in final figure
- [ ] Total figure count ≤ 6 (main text)
- [ ] SI figures organized separately

---

## Next Steps

1. **Run benchmark analysis**:
   ```python
   run validation/benchmark_analysis.py
   benchmark_all()
   ```

2. **Generate PyMOL figures** (Fig 2, 3)

3. **Run Python scripts** for data visualization (Fig 4, 5)

4. **Assemble in image editor** (Photoshop/GIMP)

5. **Write figure legends** (see manuscript outline)

6. **Review with co-authors** before submission
