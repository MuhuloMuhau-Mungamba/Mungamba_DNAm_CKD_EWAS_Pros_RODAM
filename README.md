Integrative analysis of genome-wide DNA methylation and Gene Expression profiles associated with incident chronic kidney disease: The longitudinal Pros-RODAM Cohort Study.

Project Overview

This repository contains the full analysis pipeline supporting the manuscript:

"Integrative analysis of genome-wide DNA methylation and Gene Expression profiles associated with incident chronic kidney disease: The longitudinal Pros-RODAM Cohort Study."

AUTHORS
Muhulo Muhau Mungamba1,2#, Andrea Venema3,4, Manasa Kalya Puroshothama3,4,5, Karlijn A.C. Meeks1,6, Felix. P Chilunga1, Eva L. van der Linden1, Charles F. Hayfron-Benjamin1,7,8, Constance R. Sewani-Rusike2, Liffert Vogt9, Bert-Jan van den Born1,7,Samuel N. Darko10,Ellis Owusu-Dabo11, Benedicta N. Nkeh-Chungag12, Charles Agyemang1,13* and Peter Henneman3,4,5*

AFFILIATIONS
1.	Department of Public and Occupational Health, Amsterdam University Medical Centers, University of Amsterdam, Amsterdam Public Health Research Institute, Amsterdam, The Netherlands.
2.	Department of Human Biology, Faculty of Medicine and Health Sciences, Walter Sisulu University, Mthatha, South Africa.
3.	Department of Human Genetics, Epigenetics, Amsterdam Reproduction and Development, research Institute, Amsterdam University Medical Centers, Amsterdam, The Netherlands.
4.	Amsterdam Reproduction & Development, Amsterdam University Medical Centers, University of Amsterdam, Amsterdam, The Netherlands
5.	Emma Center for Personalized Medicine, Amsterdam UMC, Amsterdam, The Netherlands.
6.	Division of Endocrinology, Diabetes and Nutrition, Department of Medicine, University of Maryland School of Medicine, Baltimore, MD, USA.
7.	Department of Vascular Medicine, Amsterdam University Medical Centers, University of Amsterdam, Amsterdam Cardiovascular Sciences, Amsterdam, The Netherlands.
8.	Departments of Physiology, Anaesthesia and Critical Care, University of Ghana Medical School / Korle Bu Teaching Hospital, Ghana.
9.	Department of Internal Medicine, Section Nephrology, Amsterdam Cardiovascular Sciences, Amsterdam UMC, University of Amsterdam, Amsterdam, The Netherlands.
10.	Department of Molecular Medicine, Kwame Nkrumah University of Science and Technology, Ghana
11.	Department of Global and International Health, Kwame Nkrumah University of Science and Technology, Ghana
12.	Department of Biological and Environmental Sciences, Faculty of Natural Sciences, Walter Sisulu University, Mthatha, South Africa.
13.	Division of Endocrinology, Diabetes and Metabolism, Department of Medicine, The Johns Hopkins University School of Medicine, Baltimore, MD, United States.

* Shared authorship
  
#Corresponding author
Muhau M. Mungamba, Department of Public and Occupational Health, Amsterdam Public Health Research institute, Amsterdam University Medical Centers, University of Amsterdam, Amsterdam, The Netherlands. Email: m.m.mungamba@amsterdamumc.nl/mmungamba@wsu.ac.za. 

Abstract

Background:Chronic kidney disease (CKD) is increasing rapidly in African populations, yet its molecular determinants remain poorly characterized. While epigenome-wide studies have been largely confined to high-income settings, comparable longitudinal epigenetic data from African populations are scarce, despite distinct genetic and environmental contexts. We conducted the first longitudinal epigenome-wide study of CKD in an African population, integrating DNA methylation and gene expression to identify epigenetic changes associated with CKD and related kidney function decline.

Methods: DNA methylation was profiled using the Illumina EPIC (850K) array in 548 Ghanaian participants from the RODAM-Prospective cohort at baseline and ~6-year follow-up. Epigenome-wide association studies (EWAS) were performed for incident CKD, change in estimated glomerular filtration rate (ΔeGFR), and change in urinary albumin-to-creatinine ratio (ΔACR), adjusting for demographic, clinical, and technical covariates. Differentially methylated regions (DMRs) were identified using bumphunter with family-wise error rate control. Differential gene expression (RNA-seq) analyses were conducted, followed by expression quantitative trait methylation (eQTM) analyses to link methylation with gene expression.

Results: At baseline, significant differentially methylated positions (DMPs) were detected at INF2 (incident CKD), CLIP4 (ΔeGFR), and TMEM109 (ΔACR), as well as DMRs annotated to in SH2D1B (incident CKD) and a hypomethylated region in TMEM109 were also observed. At follow-up, a progression-related DMR in C4orf39/TRIM61 (ΔACR) was identified, with recurrent signals at hypermethylated HLA-DPB1 and hypomethylated ZFP57 and OR2T11 across timepoints. Integrative eQTM analysis revealed strong inverse methylation-expression relationships for ZFP57 and HLA-DRB5, highlighting metabolic and immune regulatory pathways as key mechanisms.

Conclusion: We identified epigenetic changes linked to CKD development in an African population. These methylation signals, supported by gene expression patterns, highlight immune and kidney-related pathways associated with CKD risk and progression. Our findings provide a first foundation for future studies investigating the potential of epigenetic markers in CKD risk stratification and prevention in African populations.

## Study Cohort

RODAM-Prospective cohort: Ghanaian migrants and non-migrants  
- Rural Ghana  
- Urban Ghana  
- Amsterdam (The Netherlands)  

---

## Software Environment

All analyses were performed in the R statistical computing environment.

Key R packages include:  
`minfi`, `limma`, `bacon`, `missMethyl`, `clusterProfiler`, `bumphunter`, `enrichR`, `ggplot2`, `dplyr`

---

## Repository Structure

The analysis pipeline is organized into sequential modules:

### 01_Preprocessing
- DNA methylation normalization (Funnorm)
- Probe filtering (sex chromosomes, cross-hybridization, SNPs)
- Generation of beta and M-values

---

### 02_EWAS
- Differentially methylated positions (DMPs)
- Outcomes:
  - Incident CKD
  - ΔeGFR
  - ΔACR
- Includes:
  - QQ plots
  - Manhattan plots
  - Volcano plots

---

### 03_DMR
- Differentially methylated regions (DMRs)
- Identified using `bumphunter`
- Multiple cutoffs explored
- Annotated to nearest genes

---

### 04_Enrichment
- Functional annotation of EWAS signals:
  - Gene Ontology (GO)
  - KEGG pathways
  - Transcription factor enrichment
  - Regulatory enrichment (eFORGE)

---

### 05_DEG
- Differential gene expression analysis (RNA-seq)  
- Conducted by collaborator

---

### 06_eQTM
- Integration of DNA methylation and gene expression  
- Identification of methylation–expression relationships  

---

### 07_Sensitivity_analysis
- EWAS repeated with additional adjustment for:
  - Hypertension
  - Diabetes
- Includes:
  - DMP analyses
  - QQ, Manhattan, Volcano plots
  - CpG-level violin plots

---

### 08_Additional_analyses
- Supplementary and reviewer-requested analyses:
  - Correlation analyses
  - Power calculations
  - Post hoc robustness checks

---

## Reproducibility Notes

- File paths in scripts are generalized for data security  
- Users must define local paths for:
  - Methylation data
  - Phenotype data  
- Scripts assume Illumina EPIC array annotation (hg19)

---

## End of Document
############################################
