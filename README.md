# Computational Screening Reveals Molecular and Colloidal Determinants of Antibiotic Capture by Alginate–Chitosan Nanoparticles

Analysis code, data and Supporting Information accompanying the study.

**Authors:** Antonia Flores-Brizuela†, Gabriela Urra†, Javier Farías-Abarca, David Alarcón-Muñoz, Martina Bravo-Sandoval, Melissa Alegría-Arcos, Jonathan Palma, Luis Morales-Quintana\*, Daniel Bustos\*

† Equal contribution. \* Corresponding authors: dbustos@ucm.cl (D.B.), luis.morales@uautonoma.cl (L.M-Q.)

Laboratorio de Bioinformática y Química Computacional (LBQC), Departamento de Medicina Traslacional, Facultad de Medicina, Universidad Católica del Maule, Talca, Chile.

> **Link to the study:** _DOI to be added upon publication (manuscript under review at the Journal of Chemical Information and Modeling)._ <!-- Replace with the DOI/journal URL once available. -->

## Overview

The study screens 161 environmentally relevant antibiotics against three
calcium-cross-linked alginate–chitosan nanoparticle compositions
(alginate:chitosan 1:1, 1:2 and 2:1) by ensemble molecular docking, then
characterises six prioritized antibiotic–nanoparticle combinations by 100 ns
explicit-solvent molecular dynamics in triplicate. The analysis has two stages,
and this repository reproduces both:

| Stage | Language | Script | Produces |
|-------|----------|--------|----------|
| **Docking** — family performance, descriptor/ligand-efficiency analysis, consensus ranking, composition selectivity | R | `scripts/docking_analysis.R` | Figures 1–4, S1–S3; docking summary tables |
| **Molecular dynamics** — nanoparticle stability, adsorption by direct contact, interaction decomposition, statistics | Python | `scripts/01`–`05` | Figures 5, 6, S4, S5; MD tables |

It does **not** re-run the docking engine or the MD simulations themselves; it
analyses their stored outputs (see *Scope*).

## Quick start

### Docking stage (R ≥ 4.2)

```bash
# install the packages listed in environment_R.txt, then from the repository root:
Rscript scripts/docking_analysis.R
```

Reads the raw Vina scores in `data/raw/docking/` and writes Figures 1–4, S1–S3
to `results/figures/` and the docking summary tables to
`results/tables/docking/`.

### Molecular-dynamics stage (Python 3.11)

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

python scripts/01_verify_inputs.py          # checksum the raw inputs
python scripts/02_nanoparticle_stability.py # stability table
python scripts/03_adsorption_metrics.py     # adsorption tables
python scripts/04_statistical_tests.py      # significance tests
python scripts/05_figures.py                # Figures 5, 6, S4, S5
```

Python scripts 02–05 must run in order: 04 and 05 read tables written by 02 and
03. Total runtime is under three minutes on a laptop; nothing here needs a GPU.

## Repository layout

```
data/raw/docking/     raw AutoDock Vina scores (VINA_11/12/21.xlsx) and
                      physicochemical descriptors (PK.xlsx) — inputs to the R stage
data/raw/md/          consolidated per-frame trajectory analyses, one workbook per observable
data/processed/       input manifest with checksums (written by script 01)
src/algchi_screen/    Python analysis library imported by the MD scripts
scripts/              docking_analysis.R (R stage); numbered Python pipeline (MD stage);
                      dynamic_classifier.tcl (the VMD script used on the raw trajectories)
results/figures/      all manuscript figures: Fig1–Fig6, FigS1–FigS5
results/tables/       MD result tables; results/tables/docking/ holds the docking summaries
supporting_information/  Supporting Information of the article (Figures S1–S5, Tables S1–S5), PDF
requirements.txt      pinned Python dependencies
environment_R.txt     R package list for the docking stage
```

Every path, threshold and system definition for the Python stage lives in
`src/algchi_screen/config.py`.

## Figures

| File | Manuscript | Content | Stage |
|------|-----------|---------|-------|
| `Fig1.png` | Figure 1 | Mean docking energy by antibiotic family | docking (R) |
| `Fig2.png` | Figure 2 | Multi-criteria prioritization (energy vs ligand efficiency) | docking (R) |
| `Fig3.png` | Figure 3 | Consensus ranking of top candidates | docking (R) |
| `Fig4.png` | Figure 4 | Composition-dependent selectivity | docking (R) |
| `Fig5.png` | Figure 5 | Nanoparticle stability over 100 ns | MD (Python) |
| `Fig6.png` | Figure 6 | Antibiotic adsorption by direct contact | MD (Python) |
| `FigS1.png` | Figure S1 | Heavy-atom count vs mean docking energy | docking (R) |
| `FigS2.png` | Figure S1 | `results/tables/docking/dock_summary.xlsx` (plotted by `scripts/docking_analysis.R`) |
| Figure S2 | Ligand-efficiency distribution by family | docking (R) |
| `FigS3.png` | Figure S3 | Consensus-score heatmap, top 30 | docking (R) |
| `FigS4.png` | Figure S4 | One-dimensional density profiles | MD (Python) |
| `FigS5.png` | Figure S5 | Radius-classifier artifact | MD (Python) |

## Supporting Information

| File | Content |
|------|---------|
| `supporting_information/Supporting_Information_Figures.pdf` | Figures S1–S5 with captions |
| `supporting_information/Supporting_Information_Tables.pdf` | Tables S1–S5 |

Data behind the figures and Supporting Information tables:

| Item | Source file(s) |
|------|----------------|
| Figure 1 | `results/tables/docking/Table_Fig1.xlsx` |
| Figure 2 | `results/tables/docking/Table_Fig2.xlsx` |
| Figure 3 | `results/tables/docking/Table_Fig3.xlsx` |
| Figure 4 | `results/tables/docking/Table_Fig4.xlsx` |
| Figure 5 | `data/raw/md/NP_geometry_consolidado.xlsx`, `data/raw/md/Rg_consolidado.xlsx`, `results/tables/polymer_network.csv` (plotted by `scripts/05_figures.py`) |
| Figure 6 | `results/tables/adsorption_per_molecule.csv`, `adsorption_by_replica.csv`, `contact_events.csv` |
| Figure S1 | `results/tables/docking/dock_summary.xlsx` (plotted by `scripts/docking_analysis.R`) |
| Figure S2 | `results/tables/docking/Table_FigS2.xlsx` |
| Figure S3 | `results/tables/docking/Table_FigS3.xlsx` |
| Figure S4 | `data/raw/md/Density_1-1.xlsx`, `Density_1-2.xlsx`, `Density_2-1.xlsx` (plotted by `scripts/05_figures.py`) |
| Figure S5 | `results/tables/classifier_vs_contact.csv`, `data/raw/md/Classifier_positions_consolidado.xlsx` (plotted by `scripts/05_figures.py`) |
| Table S2 | `results/tables/docking/Table_S1_spearman_descriptors_vs_mean_energy.xlsx` |
| Table S3 | `results/tables/stability_summary.csv` |
| Table S4 | `results/tables/adsorption_summary.csv` |
| Table S5 | `results/tables/classifier_vs_contact.csv` |

## Data provenance

`data/raw/docking/` holds the AutoDock Vina ensemble docking scores aggregated
over 11 nanoparticle conformers × 5 docking grids per composition
(`VINA_11/12/21.xlsx`) and the physicochemical descriptors of the library
(`PK.xlsx`). The lowest-level per-frame docking engine output (per-grid logs and
poses) is bulky and is available from the corresponding authors on request.

`data/raw/md/` holds the per-frame analyses of the production trajectories:
molecule positions and contact states, nanoparticle geometry, radius of
gyration, solvent-accessible surface area, interaction counts, density profiles
and mean-squared displacements. Each workbook stores one sheet per system (or
per system–antibiotic) with the three replicas side by side as `Rep1_`, `Rep2_`
and `Rep3_` column prefixes. These were produced from the Desmond trajectories
with `scripts/dynamic_classifier.tcl` under VMD, included for completeness; it
requires the trajectory files, which are too large to distribute here and are
available from the corresponding authors on request.

## Two points of method worth knowing before reading the code

**Adsorption is measured by direct interatomic contact, not by radial
position.** A molecule counts as adsorbed in a frame when at least one of its
heavy atoms lies within 3.5 Å of a polymer heavy atom (`CONTACT_CUTOFF_A` in
`config.py`). The alternative — classifying molecules by whether their centre of
mass falls inside a nanoparticle radius — fails when the particle deforms,
because the classification boundary can expand past molecules that never touched
the polymer. Script 03 quantifies that failure
(`results/tables/classifier_vs_contact.csv`): in the 1:2 system 61.8% of frames
labelled internalized by the radial criterion record no polymer contact at all.
Both measures are reported; only the contact measure is used for any claim.

**Interaction counts are decomposed before use.** All antibiotics within a
system share one trajectory, so the interaction totals recorded against each of
them contain the same large polymer-internal contribution. Taking the
frame-wise minimum across the antibiotics of a system estimates that shared
term: subtracting it gives the antibiotic-specific excess, while the minimum
itself, normalized per chain, describes the polymer network
(`src/algchi_screen/interactions.py`).

## Reproduction status

The Python stage was executed end to end against the data in `data/raw/md/` on
Python 3.11, and the MD tables and figures in `results/` are that run's output.
The pipeline is deterministic: it reads stored per-frame analyses and computes
summary statistics, with no random sampling or model fitting, so re-running on
the same inputs reproduces the tables bit for bit. The only stochastic element
is the cosmetic horizontal jitter in Figure 6a, which uses a seeded generator
and carries no information. The R stage is likewise deterministic given the
stored Vina scores.

The molecular dynamics simulations themselves are stochastic and are not re-run
here; independent simulations with different initial velocities would not
reproduce these trajectories frame for frame. The findings that depend on
replicate agreement — particle stability and the ordering of formulations by
adsorption — were consistent across all three replicas, and replicate spread is
reported as a standard deviation throughout.

## Scope

This repository covers the analysis of docking scores and trajectories that
already exist. It does not include system construction (Packmol), simulation
setup or the Desmond production runs, nor the AutoDock Vina docking runs
themselves; those are described in the manuscript's Materials and Methods, and
their aggregated outputs are provided here as inputs.

## Citation

If you use this code or data, please cite the accompanying article (reference to
be completed on publication); citation metadata are in `CITATION.cff`. Contact
the corresponding authors for the trajectory and raw docking files.

## License

Code (`scripts/`, `src/`) is released under the MIT License. Data, results tables,
figures and Supporting Information are released under CC BY 4.0. See `LICENSE`.
