# ASAP_ULK1_sleep_cognition_MACE

This repository contains the R code used for the analysis of the study:

**"Associations of sleep parameters and autophagy biomarker ULK1 with cognitive function and cardiovascular outcomes in the Akershus Sleep Apnea Project (ASAP)"**

---

## Overview

This project investigates the associations of sleep-related hypoxia (SpO₂) and the autophagy biomarker ULK1 with:

- Cognitive function (cross-sectional analysis)
- Major adverse cardiovascular events (MACE) (survival analysis)

The analytical workflow includes:

- Data preparation and variable construction
- Cross-sectional linear regression analyses (Models 1–5)
- Survival analyses using Cox proportional hazards models
- Multiple imputation for missing covariates
- Sensitivity analyses
- Figure generation (publication-ready forest plots and diagnostics)

---

## Repository structure

```text
scripts/
├── 01_prepare_analysis_dataset.R
├── 02_cross_sectional_analysis.R
├── 03_survival_analysis.R
├── 04_multiple_imputation_sensitivity.R

README.md
```

---

## Data availability

The data used in this study are derived from the **Akershus Sleep Apnea Project (ASAP)** and are not publicly available due to ethical and data protection restrictions.

Access to the data may be granted upon reasonable request to the data custodians, subject to appropriate approvals.

This repository therefore does **not** include raw data or site-specific data import paths.

---

## Reproducibility

All analyses were conducted in R.

The scripts are organized according to the analytical workflow and can be run sequentially:

- `01_prepare_analysis_dataset.R`  
  Constructs analysis-ready datasets and defines exposures, outcomes, covariates, and survival variables  

- `02_cross_sectional_analysis.R`  
  Performs cross-sectional linear regression analyses for cognitive outcomes  

- `03_survival_analysis.R`  
  Performs survival analyses for MACE, including Cox models and Kaplan–Meier curves  

- `04_multiple_imputation_sensitivity.R`  
  Conducts multiple imputation for missing covariates and repeats analyses as sensitivity analyses  

Users must first load their approved dataset into R as:

```r
asap_raw
```
before running the scripts.

---

## Software and packages

All analyses were conducted using R (version 4.3.1).

Key packages include:

- tidyverse
- survival
- mice
- broom
- ggplot2
- patchwork


---

## Notes
Only covariates were imputed in the multiple imputation analyses
Exposures (SpO₂, ULK1) and outcomes were not imputed
AHI was treated as an adjustment covariate (categorical) rather than an exposure

---

## Code availability

The full analysis code is provided in this repository.

---

## Contact

For questions regarding the code or analysis, please contact the corresponding author.
