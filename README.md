# Colorectal Cancer TLS Microenvironment and Prognosis Analysis Project

## Project Overview

This project performs a multi-omics integrative analysis to study prognosis in colorectal cancer. We use MOVICS (Multi-Omics Integration and Visualization in Cancer Subtyping) to combine transcriptomic, methylation, somatic mutation, and lncRNA data for molecular subtyping and prognosis association analysis. A particular focus is placed on the role of tertiary lymphoid structures (TLS) in colorectal cancer prognosis.

## Data Sources

- TCGA-COAD and TCGA-READ cohort
- GEO validation cohort (GSE39582)
- MSK-IMPACT colorectal cancer cohort

## Analyses Performed

### 1. Multi-omics Integration
- mRNA expression analysis
- DNA methylation analysis
- Somatic mutation analysis
- lncRNA expression analysis

### 2. Molecular Subtyping
- Apply 10 clustering algorithms for subtype discovery
- Derive consensus subtypes from combined results
- Visualize and validate subtype assignments

### 3. Survival Analysis
- Compare survival between molecular subtypes
- Build a 6-gene prognostic model composed of:
  - PDZD4
  - PPP1R1A
  - PCOLCE2
  - ACSL6
  - CALB2
  - PTH1R

### 4. Functional Analysis
- Differential expression analysis
- KEGG pathway enrichment
- Immune infiltration analyses
- Drug sensitivity prediction

### 5. Clinical Correlation
- Association with clinicopathologic features
- Multivariable Cox regression analyses
- Nomogram construction for individualized prediction

## Key Findings

1. Two molecular subtypes related to TLS features were identified.
2. A 6-gene prognostic signature was developed and evaluated.
3. Several potential therapeutic targets and biomarkers were highlighted.
4. The model was validated in independent cohorts.

## Code Structure

- `Colorectal_MOVICS_TLS.R`: Main analysis workflow
- `supplementary_code.R`: Supplementary analyses and plotting
- `data/`: Raw and processed data files

## Usage

### Environment requirements
```R
# R >= 4.0
# MOVICS
# survival
# survminer
# ggplot2
# dplyr
# Other R packages used in the scripts (see script headers for full list)
```

## References
Workflow and functions primarily follow the MOVICS vignette: https://xlucpu.github.io/MOVICS/MOVICS-VIGNETTE.html

## Contact
18135079495@163.com
