# Incomplete Resolution of Senescence-Associated Programs After Traumatic Brain Injury

This repository contains the complete analysis code for:
**Sheng J, et al. "A transient p21 wave and a persistent microglial p16 signature define incomplete resolution of senescence-associated programs after traumatic brain injury"** (Aging Cell, 2026, in submission)

## Overview

Three public mouse TBI single-cell/single-nucleus RNA-seq datasets (6 h to 6 months post-injury) were reanalyzed to track Cdkn1a (p21) and Cdkn2a (p16) detection per cell type and per mouse, after ambient-RNA removal, doublet exclusion and depth matching.

## Data

All data are from GEO:
- **GSE269748** (Jha et al., Neuron 2024): 262,226 cells; CCI; 24 h / 7 d / 6 mo
- **GSE277487** (Ji et al., Nat Commun 2026): 162,890 nuclei; CCI; 6 h / 2 d / 4 d
- **GSE271769** (Gupta et al., bioRxiv 2024): 139,219 cells; closed-skull; D1 / D7 / D28
- **GSE319409**: Visium spatial transcriptomics (7 d)
- **GSE290150**: SCENIC contrast (24 h)

## Pipeline

| Script | Step |
|--------|------|
| `s1_01`–`s1_02` | Ingest + QC + decontX + scDblFinder |
| `s1_03` | nUMI downsampling |
| `s1_04` | Cell-type label transfer |
| `s1_05a,b` | Cell table + pre-registered readouts |
| `s1_07`–`s1_09` | Identity audits |
| `s1_12`–`s1_13` | Dose gradient + ipsi/contra |
| `s1_14` | 7d/6mo arm + p21/p16 dual-gene readout |
| `s1_16a,b` | GSE277487 decontamination panel + stats |
| `s1_17a,b2` | Phenotyping (depth-controlled) |
| `s1_18d` | Optimal transport (origin analysis) |
| `s1_19` | Spatial transcriptomics (Visium) |
| `s1_20` | Translational candidates |
| `s1_21`–`s1_21d` | GSE271769 pipeline + floor sensitivity + depth stratification |
| `s1_22`–`s1_24` | Atlas UMAP + gap-fill analyses |

## Requirements

- R 4.5.3: Seurat 5.5.0, decontX (celda) 1.8.0, scDblFinder 1.24.10, nichenetr 2.2.1.1
- Python 3.14: scanpy 1.12.1, scipy 1.16.3, POT (optimal transport)

## Key findings

1. A glial p21 (Cdkn1a) wave peaks within 24 h and resolves in a cell-type-staggered manner (astrocytes by 7 d; microglia between D28 and 6 mo)
2. A microglial p16 (Cdkn2a) program persists to D28 (two datasets) and 6 months (one dataset)
3. Acute p21+ microglia are Mcl1-high; chronic p16+ microglia are Cst7+/Itgax+
4. Composite senescence scores (SenMayo, CellAge) miss the persistent p16 arm

## License

MIT

## Citation

If you use this code, please cite the original data sources (GSE269748, GSE277487, GSE271769) and our manuscript.
