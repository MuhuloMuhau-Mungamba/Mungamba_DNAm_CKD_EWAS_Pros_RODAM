Pros-RODAM CKD EWAS
Folder: 02_EWAS

Overview
This folder contains all scripts related to epigenome-wide association analyses (EWAS) performed in the Prospective RODAM cohort.

The EWAS analyses were conducted to identify differentially methylated positions (DMPs) associated with:

- Incident chronic kidney disease (CKD)
- Change in estimated glomerular filtration rate (ΔeGFR)
- Change in urinary albumin-to-creatinine ratio (ΔACR)

Analyses were performed at:
- Baseline (T1)
- Follow-up (T2)

All models were adjusted for demographic, biological, and technical covariates, including:
- Age, sex, and study site
- Estimated blood cell proportions (CD8T, CD4T, NK, Bcell, Mono, Neu)
- Technical variables (plate/batch)

---

Folder Structure

02a_EWAS_CKD_incidence_T1_T2.R
EWAS models for incident CKD at baseline and follow-up.

02b_EWAS_delta_eGFR_T1_T2.R
EWAS models for change in eGFR (ΔeGFR).

02c_EWAS_delta_ACR_T1_T2.R
EWAS models for change in albuminuria (ΔACR).

02d_QQplots_allEWAS.R
Generation of QQ plots for all EWAS results, including calculation of genomic inflation factor (lambda).

02e_Manhattan_allEWAS.R
Generation of Manhattan plots for all EWAS analyses.

02f_Volcano_allEWAS.R
Generation of volcano plots visualising effect sizes and statistical significance.

---

Inputs

- Normalised beta values and M-values (from 01_Phenofile/normalisation)
- Phenotype files with matched cases and controls and estimated cell proportions

---

Outputs

- EWAS summary statistics (DMP tables)
- QQ plots (with lambda estimates)
- Manhattan plots
- Volcano plots

---

Reproducibility Notes

- Scripts assume access to individual-level DNA methylation and phenotype data.
- File paths should be adapted to the local computing environment.
- Due to ethical and data protection constraints, raw data are not publicly available.

---

Relationship to Manuscript

These scripts correspond to:
- Differentially Methylated Positions (DMPs)
- Primary EWAS analyses
- Supplementary QC figures (QQ, Manhattan, volcano plots)

All analyses are aligned with the Methods and Results sections of the manuscript.