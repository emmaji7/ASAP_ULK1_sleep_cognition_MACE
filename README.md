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
session_info.txt
