############################################
# Pros-RODAM CKD EWAS
# 07_Sensitivity_analysis
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
############################################

## Overview

This folder contains all sensitivity analyses conducted for the CKD epigenome-wide association study (EWAS) in the Pros-RODAM cohort.

The sensitivity analysis was performed to evaluate the robustness of the main EWAS findings by additionally adjusting for:
- Hypertension (HTN)
- Diabetes Mellitus (DM)

All analyses were conducted at both:
- Baseline (Timepoint 1; T1)
- Follow-up (Timepoint 2; T2)

---

## Objectives

- Assess whether EWAS signals remain consistent after adjusting for major cardiometabolic risk factors
- Evaluate potential confounding by hypertension and diabetes
- Confirm robustness of identified CpGs across models

---

## Folder Structure

### 1. Differential Methylation Analyses (DMPs)

- `07a_DMP_Sensitivity_CKD_T1_T2.R`  
- `07b_DMP_Sensitivity_DeltaeGFR_T1_T2.R`  
- `07c_DMP_Sensitivity_DeltaACR_T1_T2.R`  

**Description:**
- Linear models (limma) including HTN and DM as covariates
- Outcomes:
  - Incident CKD (binary)
  - ΔeGFR (continuous)
  - ΔACR (continuous)
- Models adjusted for:
  - Age, sex, site
  - Estimated blood cell proportions
  - Plate
  - Diabetes (DM)
  - Hypertension (HTN)

**Output:**
- TopTables with:
  - Raw p-values
  - FDR-adjusted p-values
  - Bacon-corrected p-values

---

### 2. Quality Control and Visualization

#### QQ Plots
- `07d_QQplots_Sensitivity_allEWAS.R`

**Description:**
- Generates QQ plots for all traits (T1 and T2)
- Computes genomic inflation factor (λ)

---

#### Manhattan Plots
- `07e_Manhattan_Sensitivity_allEWAS.R`

**Description:**
- Genome-wide visualization of CpG associations
- Displays Bonferroni significance threshold

---

#### Volcano Plots
- `07f_Volcano_Sensitivity_allEWAS.R`

**Description:**
- Visualizes effect size vs significance
- Highlights significant CpGs based on predefined thresholds

---

### 3. CpG-Level Visualization (Biological Interpretation)

#### Violin Plots – CKD
- `07g_ViolinPlots_Sensitivity_CKD_T1_T2.R`

**Description:**
- Methylation levels (%) for significant CKD CpGs
- Stratified by:
  - Site
  - CKD status (case vs control)

---

#### Violin Plots – ΔeGFR
- `07h_ViolinPlots_Sensitivity_DeltaeGFR_T1_T2.R`

**Description:**
- Methylation levels (%) for ΔeGFR-associated CpGs
- Stratified by site

---

#### Violin Plots – ΔACR
- `07i_ViolinPlots_Sensitivity_DeltaACR_T1_T2.R`

**Description:**
- Methylation levels (%) for ΔACR-associated CpGs
- Stratified by site

---

## Statistical Analyses

- Linear modeling: `limma`
- Inflation correction: `bacon`
- Multiple testing correction: FDR
- Non-parametric testing:
  - Kruskal–Wallis tests for site differences in methylation levels

---

## Input Data

- DNA methylation data:
  - Beta values
  - M-values (functional normalization)
- Phenotype data:
  - CKD status
  - eGFR and ACR measures
  - Covariates (age, sex, site, cell counts)
  - Hypertension and diabetes status

---

## Output Files

- TopTables (sensitivity models)
- QQ plots (PDF)
- Manhattan plots (PNG)
- Volcano plots (JPEG)
- Violin plots (TIFF; publication quality)

---

## Interpretation Notes

- Sensitivity analyses assess robustness rather than discovery
- Consistent CpGs across main and sensitivity models strengthen causal inference
- Differences may indicate:
  - Confounding by cardiometabolic factors
  - Mediation pathways (e.g., via HTN or DM)

---

## Notes

- All file paths in scripts are generalized for reproducibility
- Users must define local paths to:
  - Methylation matrices
  - Phenotype datasets
- Scripts assume EPIC array annotation (hg19)

---

## End of README
############################################