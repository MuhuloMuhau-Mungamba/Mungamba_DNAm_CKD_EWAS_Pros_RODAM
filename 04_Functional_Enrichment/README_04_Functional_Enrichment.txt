Functional Enrichment Analysis (04_Functional_Enrichment)
Overview

This folder contains functional annotation analyses of DNA methylation signals associated with chronic kidney disease (CKD) identified in the EWAS (02_EWAS) and DMR (03_DMR) steps.

The purpose of these analyses is to provide biological interpretation of CKD-associated CpGs by identifying enriched:

Biological pathways
Gene ontologies
Transcription factor binding profiles
Regulatory genomic regions

All analyses were conducted separately for:

Baseline (Timepoint 1)
Follow-up (Timepoint 2)
Input Data

The enrichment analyses are based on EWAS-derived CpGs (DMPs), including:

CKD incidence
ΔeGFR (change in kidney function)
ΔACR (change in albuminuria)

CpGs were:

Ranked by bacon-corrected p-values
Combined across traits within each timepoint
Filtered and deduplicated prior to enrichment analyses
Analyses Performed
1. Gene-Based Pathway Enrichment (clusterProfiler)

Script:
04a_clusterProfiler_GO_KEGG.R

Description:

CpGs mapped to genes using annotation (UCSC_RefGene_Name)
Gene symbols converted to Entrez IDs
Enrichment performed for:
Gene Ontology (GO; Biological Processes)
KEGG pathways

Purpose:

Provide exploratory biological interpretation of associated genes

Output:

Results/clusterProfiler/GO_baseline.txt
Results/clusterProfiler/KEGG_baseline.txt
Results/clusterProfiler/GO_followup.txt
Results/clusterProfiler/KEGG_followup.txt
2. Bias-Corrected Pathway Enrichment (missMethyl)

Script:
04b_missMethyl_GO_KEGG.R

Description:

Uses gometh() to correct for unequal CpG representation per gene
CpGs ranked by significance
Enrichment performed using:
Top 1000 CpGs
Top 5000 CpGs
Background set: all tested CpGs

Purpose:

Primary pathway enrichment analysis for interpretation

Output:

Results/missMethyl/GO_top1000_baseline.txt
Results/missMethyl/GO_top5000_baseline.txt
Results/missMethyl/GO_top1000_followup.txt
Results/missMethyl/GO_top5000_followup.txt
3. Transcription Factor Enrichment (enrichR)

Script:
04c_TF_enrichment_enrichR.R

Description:

Gene sets derived from CpGs (p < 1×10⁻³)
Enrichment performed using:
ChEA 2016
ENCODE + ChEA consensus TF datasets

Purpose:

Identify regulatory transcription factors linked to methylation signals

Output:

Results/TF_enrichment/TF_baseline.txt
Results/TF_enrichment/TF_followup.txt
4. Regulatory Element Enrichment Preparation (eFORGE)

Script:
04d_eFORGE_preparation.R

Description:

Preparation of CpG lists for external eFORGE analysis
CpGs ranked by significance
Top CpGs selected:
Top 1000 CpGs (primary)
Top 5000 CpGs (optional)

Purpose:

Assess enrichment in:
DNase I hypersensitive sites
Tissue-specific regulatory regions

Output:

Results/eFORGE/top1000_baseline.txt
Results/eFORGE/top1000_followup.txt

Note:

eFORGE analysis is performed externally using these input files
Output Structure

All results are organized into subfolders:

Results/clusterProfiler/ → gene-based enrichment
Results/missMethyl/ → bias-corrected enrichment (primary results)
Results/TF_enrichment/ → transcription factor enrichment
Results/eFORGE/ → CpG lists for regulatory enrichment
Interpretation

Functional enrichment analyses indicate that CKD-associated methylation signals are enriched in:

Biological processes related to:
Cell signaling
Intercellular communication
Tissue remodeling
Immune and inflammatory pathways
Regulatory mechanisms involving:
Epigenetic repression (e.g., Polycomb complex)
Fibrotic signaling pathways (e.g., TGF-β)
Oxidative stress response
Genomic regions associated with:
Open chromatin (DNase I hypersensitive sites)
Epithelial and stromal cell types
Notes
clusterProfiler results are exploratory and gene-based
missMethyl results account for probe-number bias and are prioritized for interpretation
Transcription factor enrichment provides regulatory context
eFORGE results are generated externally using prepared CpG lists