# Differentially Methylated Regions (DMR) Analysis

## Overview

This folder contains scripts and outputs for the identification of **Differentially Methylated Regions (DMRs)** in the Pros-RODAM cohort using the **bumphunter** algorithm.

Analyses were conducted for:
- Incident CKD (binary outcome)
- Change in kidney function:
  - ΔeGFR (delta eGFR)
  - ΔACR (delta log-transformed ACR)

All analyses were performed at:
- Baseline (T1)
- Follow-up (T2)

---

## Study Design

Participants were selected from a prospective cohort free of CKD at baseline.

- Incident CKD cases were identified at follow-up
- DNA methylation data were available for:
  - 92 cases
  - 736 controls

A **1:2 matched nested case–control design** was applied:
- Matching variables:
  - Age
  - Sex
  - Study site

Final analytical sample:
- 92 cases
- 182 controls
- Total per timepoint: 274 participants

The same individuals were used at baseline (T1) and follow-up (T2).

---

## Methods

### DNA Methylation Processing

- Platform: Illumina EPIC array
- Functional normalisation (minfi::preprocessFunnorm)
- Probe filtering:
  - Removal of sex chromosome probes (chrX, chrY)
  - Removal of cross-hybridising probes (population-specific SNP filtering)

---

### DMR Identification

DMRs were identified using the **bumphunter** algorithm.

For each outcome:
1. Initial scan without bootstrapping to explore candidate regions
2. Final analysis with **bootstrap (B = 500)**

Criteria:
- Minimum probes per region: **> 2**
- Regions annotated to nearest gene using UCSC annotation

---

### Model Specification

All models were adjusted for:

- Age  
- Sex  
- Study site  
- Blood cell proportions:
  - CD8T, CD4T, NK, B cells, Monocytes, Neutrophils  
- Technical covariate:
  - Plate  

---

## Cutoff Selection

Different cutoffs were used based on outcome scale:

| Outcome        | Cutoff (β difference) |
|----------------|----------------------|
| CKD (binary)   | 0.015                |
| ΔeGFR          | 0.0006               |
| ΔACR (log)     | 0.00035              |

---

## Scripts

| Script | Description |
|--------|------------|
| `03a_DMR_CKD_incident_T1_T2.R` | DMR analysis for incident CKD |
| `03b_DMR_delta_eGFR_T1_T2.R` | DMR analysis for ΔeGFR |
| `03c_DMR_delta_ACR_T1_T2.R` | DMR analysis for ΔACR |

Each script:
- Runs baseline (T1) and follow-up (T2) analyses
- Performs filtering, modelling, and DMR detection
- Outputs annotated DMR tables

---

## Outputs

Output files contain:

- Genomic region (chr, start, end)
- Number of CpGs (L)
- Effect size estimate
- Statistical significance
- Gene annotation (nearest gene)

Example outputs:
- `DMRS_CKD_T1_cutoff0.015_boot500.txt`
- `DMRS_deltaeGFR_T2_cutoff0.0006_boot500.txt`
- `DMRS_deltaACR_T2_cutoff0.00035_boot500.txt`

---

## Reporting

Results were summarised in tables for manuscript reporting.

No graphical visualisation was performed for DMRs. Instead:
- Significant regions were prioritised based on:
  - Effect size
  - Statistical significance
  - Biological relevance

---

## Notes

- Baseline and follow-up analyses use identical participants
- Missing values were removed for analyses prior to modelling
- Scripts require user-defined paths to:
  - IDAT files
  - Phenotype files
  - Probe filtering lists

---

## End