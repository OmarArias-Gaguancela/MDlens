# MDlens

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](COLAB_LINK_HERE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![MDAnalysis](https://img.shields.io/badge/MDAnalysis-2.x-orange)](https://www.mdanalysis.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)

**Multi-Complex, Multi-Replicate Molecular Dynamics Comparative Analysis**

Built by **Omar Arias-Gaguancela, PhD** | SciLearningWorkshops LLC

---

## Overview

MDlens is a Google Colab-based Python notebook for **comparative analysis of multiple protein–ligand MD simulations**. It handles multiple complexes simultaneously, supports multiple replicates per complex (with automatic averaging and error-band visualization), and produces publication-ready plots alongside quantitative tables.

MDlens is built on [MDAnalysis](https://www.mdanalysis.org/) and follows a clean two-class architecture (`MDComplex` + `MDComparator`) that keeps data loading and computation separate from all plotting logic. Every analysis block is wrapped in robust error handling so a failure in one metric never crashes the full notebook.

> **Looking for the single-complex version?**
> Check out [MD_quick_plot](https://github.com/OmarArias-Gaguancela/MD_quick_plot),
> the predecessor to MDlens for single-complex, single-replicate MD analysis.

---

## Architecture

### `MDComplex`
Represents **one protein–ligand system**. Accepts one topology and one or more trajectory files (replicates). Responsibilities:
- Load each replicate as a separate `mda.Universe`
- Compute all metrics per replicate, then average across replicates with standard deviation (error bands)
- Handle unequal replicate lengths via linear interpolation before averaging
- Store results internally so `MDComparator` can request them on demand

### `MDComparator`
Accepts a list of `MDComplex` objects. Responsible for **all plotting** — both overlay (multiple complexes on one figure) and per-complex (one figure per complex). Runs a three-level validation check at initialization and routes around incompatible analyses automatically.

---

## Features

| Feature | Description | Output |
|---|---|---|
| **Backbone RMSD** | vs time, normalized to frame 0, with replicate error bands | `MDlens_rmsd.png` |
| **Ligand RMSD** | Align protein backbone first, then RMSD of ligand heavy atoms | `MDlens_ligand_rmsd.png` |
| **Radius of Gyration** | Protein Rg vs time | `MDlens_rg.png` |
| **H-bond Count** | Protein–ligand H-bonds vs time (bidirectional) | `MDlens_hbonds.png` |
| **P–L Distance** | Minimum protein–ligand distance with 4 Å contact threshold | `MDlens_pl_distance.png` |
| **MM Binding Energy** | Simplified Coulomb + Lennard-Jones estimate in kJ/mol | `MDlens_binding_energy.png` |
| **RMSF Overlay** | Per-residue Cα RMSF, mean ± SD across replicates | `MDlens_rmsf.png` |
| **Segmented RMSF Heatmap** | Trajectory split into N segments, RMSF per segment as heatmap | `MDlens_rmsf_heatmap_[name].png` |
| **ΔRMSF Table** | Pairwise residue flexibility comparison, ranked and exported | `MDlens_delta_rmsf_table.csv` |
| **Contact Residues** | Residues within 5 Å of ligand, ranked by contact frequency | `MDlens_contact_residues_[name].csv` |
| **Binding Site PCA** | Joint PCA on contact-residue Cα positions, 95% ellipses | `MDlens_binding_site_pca.png` |
| **Free Energy Landscape** | RMSD vs Rg, Gaussian-smoothed, global minimum reported | `MDlens_fel_[name].png` |
| **Frame Extraction** | Save PDB snapshot at any requested simulation time | user-specified |

---

## Quick Start

```python
from MDlens import MDComplex, MDComparator  # or run Cell 3 in Colab

# Complex with three replicates
complex_a = MDComplex(
    name="Compound A",
    topology="compA.pdb",
    trajectories=["rep1.dcd", "rep2.dcd", "rep3.dcd"],
    protein_selection="protein",
    ligand_selection="resname LIG",
    dt_in_ps=None,       # reads timestep from trajectory; override if needed
    resid_offset=0,
)

# Complex with a single replicate
complex_b = MDComplex(
    name="Compound B",
    topology="compB.pdb",
    trajectories="single.dcd",    # string accepted for single replicate
    protein_selection="protein",
    ligand_selection="resname LIG",
)

# Compute all metrics
complex_a.compute(n_segments=10)
complex_b.compute(n_segments=10)

# Compare
comparator = MDComparator(
    complexes=[complex_a, complex_b],
    n_segments=10,
    delta_rmsf_threshold=0.0,   # Å; 0 = show all residues
    contact_cutoff=5.0,         # Å
    contact_min_freq=0.30,      # fraction of frames
)

comparator.plot_all_overlay()          # Cell 5
comparator.run_per_complex_analyses()  # Cell 6 — FEL + heatmaps
comparator.run_contact_and_delta_rmsf()# Cell 7
comparator.plot_binding_site_pca()     # Cell 8
```

---

## Supported Input Formats

| File type | Formats accepted |
|---|---|
| **Topology** | PDB, GRO, PSF, TOP (any MDAnalysis-supported format) |
| **Trajectory** | XTC, TRR, DCD, NetCDF (.nc), and any MDAnalysis-supported format |

---

## Validation and Warning System

MDlens runs a three-level validation check automatically when `MDComparator` is initialized:

**Level 1 — Hard stop (incompatible proteins):**
```
❌ CRITICAL: Protein size mismatch detected.
   Complex1 (Compound A): 487 residues
   Complex2 (Compound B): 312 residues

   Binding site PCA and ΔRMSF comparison require the same protein target.
   These analyses will be SKIPPED.
   All other analyses (RMSD, RMSF, Rg, FEL, etc.) will still run independently.
```

**Level 2 — Warning (low contact residue overlap):**
```
⚠️  WARNING: Contact residue overlap is 67% between Compound A and Compound B.
   Compound A contact residues: [87, 142, 156, 201, 203, 211]
   Compound B contact residues: [142, 156, 189, 201, 203, 244]
   Shared residues used for PCA: [142, 156, 201, 203]

   Tip: If structures have different residue numbering, use resid_offset parameter.
   e.g. complex2 = MDComplex(..., resid_offset=10)
```

**Level 3 — Passive advisory (always shown before PCA):**
```
ℹ️  PCA NOTE: Binding site PCA comparison assumes the same protein target
    across all complexes. If comparing different proteins, results are not valid.
```

---

## Notebook Cell Structure

| Cell | Contents |
|---|---|
| Cell 1 | Install Required Packages |
| Cell 2 | Mount Google Drive & Set Working Directory |
| Cell 3 | `MDComplex` and `MDComparator` Class Definitions |
| Cell 4 | Define Your Complexes (user edits file paths here) |
| Cell 5 | Run Comparator: All Overlay Plots |
| Cell 6 | Run Per-Complex Analyses (FEL, RMSF Heatmaps) |
| Cell 7 | Contact Residue Detection & ΔRMSF Table |
| Cell 8 | Binding Site PCA |
| Cell 9 | Individual Analysis Options (single metrics on demand) |
| Cell 10 | Frame Extraction Utility |

---

## Dependencies

```
MDAnalysis
matplotlib
seaborn
numpy
scipy
pandas
scikit-learn
ipywidgets
tqdm
```

All packages are auto-installed in Cell 1.

---

## Output Files

```
MDlens_rmsd.png
MDlens_rmsf.png
MDlens_rmsf_heatmap_[complexname].png
MDlens_rg.png
MDlens_ligand_rmsd.png
MDlens_hbonds.png
MDlens_binding_energy.png
MDlens_pl_distance.png
MDlens_fel_[complexname].png
MDlens_binding_site_pca.png
MDlens_delta_rmsf_table.csv
MDlens_contact_residues_[complexname].csv
```

---

## Citations

If you use MDlens in your research, please cite:

- **MDlens:** Arias-Gaguancela, O. (2025). MDlens: Multi-complex, multi-replicate MD simulation comparative analysis. SciLearningWorkshops LLC. https://github.com/OmarArias-Gaguancela/MDlens

- **MDAnalysis:** Michaud-Agrawal, N. et al. (2011) MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations. *J. Comput. Chem.* 32, 2319–2327. DOI: [10.1002/jcc.21787](https://doi.org/10.1002/jcc.21787)

- **Segmented RMSF heatmap approach inspired by eRMSF:**  
  Arantes, P.R. et al. (2025) eRMSF: An Enhanced Root Mean Square Fluctuation Analysis for Exploring Protein Flexibility from MD Simulations. *J. Chem. Inf. Model.* **65**(23), 12648–12654. DOI: [10.1021/acs.jcim.5c00906](https://doi.org/10.1021/acs.jcim.5c00906)

---

## License

MIT License — free to use for research and education with attribution required.  
Copyright (c) 2025 Omar Arias-Gaguancela, SciLearningWorkshops LLC.  
See [LICENSE](LICENSE) for full text.

---

## Author

**Omar Arias-Gaguancela, PhD**  
SciLearningWorkshops LLC  
[GitHub](https://github.com/OmarArias-Gaguancela)
