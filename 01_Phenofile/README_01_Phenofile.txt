# Pros-RODAM CKD EWAS

## Step 01: Phenotype Preparation and Preprocessing

---

## Overview

This folder contains all scripts required to generate the phenotype files used in the epigenome-wide association study (EWAS) of chronic kidney disease (CKD) in the Pros-RODAM cohort.

The workflow includes:

* Selection of incident CKD cases
* Propensity score matching of controls (1:2 ratio)
* Construction of baseline (T1) and follow-up (T2) phenotype datasets
* Estimation of blood cell-type proportions
* DNA methylation preprocessing, normalization, and probe filtering

All steps are designed to ensure a consistent longitudinal analytic sample across timepoints.

---

## Study Design Summary

The analysis is based on a nested case–control design within the Pros-RODAM cohort:

* Individuals with CKD at baseline were excluded
* Incident CKD cases were identified at follow-up (~6 years)
* Of these, **92 cases had DNA methylation data available**
* A pool of eligible controls (n ≈ 736) was available
* Controls were matched to cases using **1:2 propensity score matching** based on:

  * Age
  * Sex
  * Study site

This resulted in:

* **92 CKD cases**
* **182 matched controls**

The same individuals were retained at both timepoints:

* **Baseline (T1): n = 274**
* **Follow-up (T2): n = 274**

Ensuring within-individual longitudinal comparability.

---

## Scripts Included

### 01a_propensity_matching_and_phenofile_creation.R

* Identifies incident CKD cases
* Performs 1:2 propensity score matching
* Generates baseline (T1) and follow-up (T2) phenotype files
* Ensures identical individuals across timepoints via RodamID

**Output:**

* CKD_EWAS_T1_phenofile.csv
* CKD_EWAS_T2_phenofile.csv

---

### 01b_cell_type_estimation_T1_T2.R

* Estimates blood cell-type proportions using the Houseman method
* Uses IDAT files and reference-based deconvolution
* Adds cell proportions to phenotype files

**Cell types estimated:**

* CD8T
* CD4T
* NK cells
* B cells
* Monocytes
* Neutrophils

**Output:**

* CKD_EWAS_T1_phenofile_wCells.csv
* CKD_EWAS_T2_phenofile_wCells.csv

---

### 01c_normalisation_and_filtering_T1_T2.R

* Performs functional normalization (funnorm)
* Removes:

  * Sex chromosome probes
  * Cross-hybridizing probes
  * Probes with infinite M-values
* Generates:

  * Raw beta values
  * Normalized beta values
  * M-values

**Output:**

* Beta and M-value matrices (T1 and T2)
* Annotation files
* Quality control plots

---

## Key Methodological Notes

* DNA methylation measured using Illumina EPIC array (~850K CpGs)
* Functional normalization performed using **minfi**
* M-values used for statistical analyses
* β-values retained for interpretation and visualization
* Blood cell composition estimated using the Houseman reference method

---

## Reproducibility Notes

* Scripts require access to:

  * Individual-level phenotype data
  * IDAT methylation files
* Data are not publicly shared due to ethical restrictions
* File paths must be adapted to local computing environments

---

## Output of This Step

This step produces the final EWAS-ready dataset:

* Matched phenotype files (T1 and T2)
* Cell-type adjusted covariates
* Normalized methylation matrices

These outputs serve as input for downstream analyses:

* EWAS (DMP analysis)
* DMR analysis
* DEG analysis
* eQTM analysis

---

## Relationship to Manuscript

This folder corresponds to the following sections:

* Methods: Study population and design
* Methods: DNA methylation processing
* Methods: Statistical analysis (preprocessing steps)

---

## Author

Muhulo Muhau Mungamba et al.
Pros-RODAM CKD EWAS Project
