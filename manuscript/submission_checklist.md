# JCIM Submission Checklist for GlueTK

**Target Journal**: Journal of Chemical Information and Modeling (JCIM)  
**Submission Deadline**: TBD  
**Expected Timeline**: 3-4 weeks preparation → Submit

---

## 📝 **Manuscript Components** (6 weeks)

### Week 1-2: Data Collection & Analysis

- [ ] **Run Benchmark Analysis**
  ```bash
  # In PyMOL
  run validation/benchmark_analysis.py
  benchmark_all()
  ```
  - [ ] Verify all 6-10 structures analyzed successfully
  - [ ] Check accuracy > 85% (target: 90%)
  - [ ] Save results CSV
  - [ ] Screenshot terminal output for SI

- [ ] **Collect Performance Metrics**
  - [ ] Runtime for each analysis step
  - [ ] Memory usage monitoring
  - [ ] Test on different hardware (if possible)

- [ ] **Parameter Sensitivity Analysis**
  - [ ] Test BSA thresholds: 600, 700, 800, 900, 1000 Ų
  - [ ] Test PPI contact thresholds: 5, 8, 10, 12 contacts
  - [ ] Plot accuracy vs threshold

### Week 3-4: Figure Generation

- [ ] **Figure 1: Architecture** (Draw.io/BioRender)
  - [ ] Create flowchart
  - [ ] Export 300 DPI PNG
  - [ ] Write legend (100 words)

- [ ] **Figure 2: CC-885 Case Study** (PyMOL)
  - [ ] Panel A: PPI interface
  - [ ] Panel B: Neo-epitope highlighting
  - [ ] Panel C: 2D interaction diagram
  - [ ] Panel D: G-loop detail
  - [ ] Assemble in Photoshop/GIMP
  - [ ] Add scale bars and labels
  - [ ] Write legend (200 words)

- [ ] **Figure 3: dBET1 PROTAC** (PyMOL)
  - [ ] Panel A: Weak PPI
  - [ ] Panel B: Linker visualization
  - [ ] Panel C: Comparison with CC-885
  - [ ] Write legend (150 words)

- [ ] **Figure 4: Benchmark Results** (Python)
  - [ ] Accuracy bar chart
  - [ ] Confusion matrix heatmap
  - [ ] Write legend (150 words)

- [ ] **Figure 5: PPI Metrics** (Python)
  - [ ] Scatter plot (PPI contacts vs BSA)
  - [ ] Decision boundaries
  - [ ] Write legend (100 words)

- [ ] **Graphical Abstract (TOC)**
  - [ ] Design simple, clear summary graphic
  - [ ] 8.5 cm × 4.75 cm @ 300 DPI
  - [ ] Test in grayscale for accessibility

### Week 5: Writing

- [ ] **Title** (finalize)
  - Recommended: "GlueTK: A PyMOL Plugin for Molecular Glue Mechanism Analysis via Protein-Protein Interface and Neo-Epitope Detection"

- [ ] **Abstract** (250 words)
  - [ ] Background (2-3 sentences)
  - [ ] Methods (3-4 sentences)
  - [ ] Results (2-3 sentences)
  - [ ] Conclusion (1-2 sentences)

- [ ] **Introduction** (800-1000 words)
  - [ ] TPD background
  - [ ] PROTAC vs Glue distinction (Table 1)
  - [ ] Computational challenges
  - [ ] Our contribution

- [ ] **Methods** (1500-2000 words)
  - [ ] Algorithm descriptions (with pseudocode)
  - [ ] Implementation details
  - [ ] Interaction parameters (Table 2)
  - [ ] Validation dataset (Table 3)

- [ ] **Results** (1000-1500 words)
  - [ ] Case Study 1: CC-885
  - [ ] Negative Control: dBET1
  - [ ] Benchmark performance (Table 4)
  - [ ] Computational performance (Table 5)

- [ ] **Discussion** (800-1000 words)
  - [ ] Advantages
  - [ ] Limitations
  - [ ] Tool comparison (Table 6)
  - [ ] Experimental validation opportunities

- [ ] **Conclusions** (200-300 words)

- [ ] **References** (40-60 refs)
  - [ ] Glue degraders (Krönke, Matyskiela, Schreiber)
  - [ ] PROTACs (Sakamoto, Gadd)
  - [ ] Computational methods (Böhm, Barlow)
  - [ ] Existing tools (PLIP, ProLIF)

### Week 6: Supporting Information

- [ ] **SI-1**: Detailed Algorithm Pseudocode
- [ ] **SI-2**: Complete Benchmark Dataset Table
- [ ] **SI-3**: Parameter Sensitivity Analysis (figures)
- [ ] **SI-4**: Installation Guide
  ```markdown
  # Installation
  1. Install PyMOL (≥2.5)
  2. Download GlueTK from GitHub
  3. Install via Plugin Manager...
  ```
- [ ] **SI-5**: Tutorial with Example
  ```python
  # Example: Analyze CC-885
  fetch 6H0F
  ppi_analyze 6H0F, [A], [B]
  neo_epitope_find 6H0F, [A], [B], CC885
  ```
- [ ] **SI-6**: Example Output Files (CSV screenshots)
- [ ] **SI-7**: Video Tutorial (optional, YouTube)
  - [ ] Record 5-min demo
  - [ ] Upload to YouTube
  - [ ] Add link to SI

---

## 💻 **Code & Data Preparation**

- [ ] **GitHub Repository**
  - [ ] Clean up code
  - [ ] Add docstrings to all functions
  - [ ] Remove debug print statements
  - [ ] Add LICENSE file (MIT recommended)
  - [ ] Update README.md with:
    - [ ] Installation instructions
    - [ ] Quick start guide
    - [ ] Citation information
    - [ ] Link to manuscript (when published)

- [ ] **Zenodo Data Deposit**
  - [ ] Create Zenodo account
  - [ ] Upload benchmark dataset
    - [ ] PDB IDs list
    - [ ] Analysis results CSV
    - [ ] Figure source data
  - [ ] Get DOI
  - [ ] Add DOI to manuscript

- [ ] **Documentation**
  - [ ] ReadTheDocs setup (optional but recommended)
  - [ ] API documentation
  - [ ] Tutorials page

---

## 📋 **JCIM-Specific Requirements**

### Manuscript Format

- [ ] **File Format**: Word (.docx) or LaTeX
- [ ] **Font**: Times New Roman, 12 pt
- [ ] **Line Spacing**: Double-spaced
- [ ] **Margins**: 1 inch all sides
- [ ] **Page Numbers**: Bottom right

### Figures

- [ ] **Resolution**: 300 DPI minimum (600 DPI for line art)
- [ ] **Format**: TIFF or EPS (PNG acceptable for web)
- [ ] **Size**: Width ≤ 7 inches (full page) or 3.5 inches (half page)
- [ ] **Color**: RGB mode
- [ ] **File Naming**: Figure1.tif, Figure2.tif, etc.

### Supplementary Material

- [ ] **SI Document**: Single PDF
- [ ] **Data Files**: Separate ZIP archive
- [ ] **Code**: GitHub link (stable release)

### Cover Letter

- [ ] **Addressed to**: Editor-in-Chief, JCIM
- [ ] **Content**:
  - [ ] Brief summary (3-4 sentences)
  - [ ] Significance statement
  - [ ] Confirmation of no conflicts of interest
  - [ ] Suggested reviewers (3-5 names)
    - Expert in TPD (e.g., from Novartis, C4 Therapeutics)
    - Computational chemist (e.g., from academic labs)
    - PyMOL/structural biology tool developer

---

## ✅ **Pre-Submission Checks**

### Technical Validation

- [ ] **Code Works on Fresh Install**
  - [ ] Test on clean PyMOL installation
  - [ ] Test on Windows, macOS, Linux (if possible)
  - [ ] All dependencies install automatically

- [ ] **Benchmark Reproduces**
  - [ ] Re-run `benchmark_all()` from scratch
  - [ ] Verify results match reported values
  - [ ] Check runtime is consistent

- [ ] **Figures Match Data**
  - [ ] Every number in figures traceable to source data
  - [ ] No discrepancies between text and figures
  - [ ] All error bars/statistics correct

### Writing Quality

- [ ] **Spell Check** (US English)
- [ ] **Grammar Check** (Grammarly or equivalent)
- [ ] **Reference Check**
  - [ ] All citations in text appear in bibliography
  - [ ] All references cited in text
  - [ ] DOIs included where available

- [ ] **Consistency**
  - [ ] Terminology consistent throughout
    - "Molecular glue" vs "molecular glue degrader"
    - "PROTAC" capitalization
  - [ ] Abbreviations defined on first use
  - [ ] Units consistent (Å, kcal/mol, etc.)

### Legal & Ethical

- [ ] **Authorship**
  - [ ] All authors approve final version
  - [ ] Contributions clearly stated
  - [ ] No ghost authors

- [ ] **Conflicts of Interest**
  - [ ] Declared (or "none")

- [ ] **Funding**
  - [ ] Acknowledged (if any)
  - [ ] Grant numbers included

- [ ] **Data Availability**
  - [ ] GitHub link active
  - [ ] Zenodo DOI valid
  - [ ] All data publicly accessible

---

## 📤 **Submission Process**

1. **Create JCIM Account**
   - URL: https://pubs.acs.org/journal/jcisd8
   - Register as author

2. **Prepare Submission Files**
   ```
   submission/
   ├── manuscript.docx
   ├── figures/
   │   ├── Figure1.tif
   │   ├── Figure2.tif
   │   ├── Figure3.tif
   │   ├── Figure4.tif
   │   ├── Figure5.tif
   │   └── TOC_graphic.tif
   ├── supporting_info.pdf
   ├── data_files.zip
   └── cover_letter.pdf
   ```

3. **Online Submission**
   - [ ] Upload manuscript
   - [ ] Upload figures separately
   - [ ] Upload SI
   - [ ] Suggest reviewers
   - [ ] Fill out metadata (keywords, etc.)

4. **Review & Revise**
   - Expected review time: 4-8 weeks
   - Be prepared for revisions
   - Address all reviewer comments

---

## 🎯 **Success Criteria**

**Manuscript Acceptance** if:
- ✅ Technical novelty (PPI + Neo-epitope integration)
- ✅ Validation on ≥6 structures with >85% accuracy
- ✅ Clear biological relevance (glue vs PROTAC distinction)
- ✅ Open-source and user-friendly
- ✅ Publication-quality figures and writing

**Impact Factor**: JCIM ~5-6 (solid computational chemistry journal)

---

## 📊 **Timeline Summary**

| Week | Task | Deliverable |
|------|------|-------------|
| 1-2 | Data collection | Benchmark results CSV |
| 3-4 | Figure generation | 5 main figures + TOC |
| 5 | Writing | Complete manuscript draft |
| 6 | SI preparation | SI document + data files |
| 7 | Internal review | Revised draft |
| 8 | Final checks | Submit! |

---

## 🔗 **Useful Links**

- **JCIM Homepage**: https://pubs.acs.org/journal/jcisd8
- **Author Guidelines**: https://publish.acs.org/publish/author_guidelines?coden=jcisd8
- **ACS Paragon Plus**: https://paragonplus.acs.org/
- **Zenodo**: https://zenodo.org/
- **GitHub**: https://github.com/yourname/glue-pymol

---

**Last Updated**: 2025-01-09  
**Status**: ✅ Ready to start Week 1 (Data Collection)
