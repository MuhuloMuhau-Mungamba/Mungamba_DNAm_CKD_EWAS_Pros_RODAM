############################################
# CKD EWAS – Functional Enrichment
# Analysis: GO and KEGG enrichment (clusterProfiler)
#
# Author: Muhulo Muhau Mungamba
# Date: 2026-04-16
#
# Description:
# Gene-based pathway enrichment analysis.
# CpGs mapped to genes → Entrez IDs → GO + KEGG.
# Baseline and follow-up analysed separately.
############################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

############################################
# Load libraries
############################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)

############################################
# Input files (relative paths)
############################################

ckd_base_file   <- "Results/Timepoint_1/DMPs/CKD/TopTable_CKD_T1.txt"
egfr_base_file  <- "Results/Timepoint_1/DMPs/Delta_eGFR/TopTable_Delta_eGFR_T1.txt"
acr_base_file   <- "Results/Timepoint_1/DMPs/Delta_ACR/TopTable_Delta_ACR_T1.txt"

ckd_follow_file <- "Results/Timepoint_2/DMPs/CKD/TopTable_CKD_T2.txt"
egfr_follow_file<- "Results/Timepoint_2/DMPs/Delta_eGFR/TopTable_Delta_eGFR_T2.txt"
acr_follow_file <- "Results/Timepoint_2/DMPs/Delta_ACR/TopTable_Delta_ACR_T2.txt"

############################################
# Function: Extract genes
############################################

extract_genes <- function(df, p_cutoff = 1e-3){
  df %>%
    filter(BaconPvalue_Mvalues < p_cutoff,
           !is.na(UCSC_RefGene_Name),
           UCSC_RefGene_Name != "") %>%
    pull(UCSC_RefGene_Name) %>%
    str_split(";") %>%
    unlist() %>%
    unique()
}

############################################
# -------- BASELINE --------
############################################

ckd_base  <- read.delim(ckd_base_file)
egfr_base <- read.delim(egfr_base_file)
acr_base  <- read.delim(acr_base_file)

combined_base <- bind_rows(ckd_base, egfr_base, acr_base)

genes_base <- extract_genes(combined_base)

gene_ids_base <- bitr(genes_base,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = org.Hs.eg.db) %>%
  pull(ENTREZID) %>%
  unique()

go_base <- enrichGO(gene_ids_base, OrgDb=org.Hs.eg.db, ont="BP", readable=TRUE)
kegg_base <- enrichKEGG(gene_ids_base, organism="hsa")

############################################
# -------- FOLLOW-UP --------
############################################

ckd_follow  <- read.delim(ckd_follow_file)
egfr_follow <- read.delim(egfr_follow_file)
acr_follow  <- read.delim(acr_follow_file)

combined_follow <- bind_rows(ckd_follow, egfr_follow, acr_follow)

genes_follow <- extract_genes(combined_follow)

gene_ids_follow <- bitr(genes_follow,
                        fromType="SYMBOL",
                        toType="ENTREZID",
                        OrgDb=org.Hs.eg.db) %>%
  pull(ENTREZID) %>%
  unique()

go_follow <- enrichGO(gene_ids_follow, OrgDb=org.Hs.eg.db, ont="BP", readable=TRUE)
kegg_follow <- enrichKEGG(gene_ids_follow, organism="hsa")

############################################
# Save results
############################################

dir.create("Results/clusterProfiler", recursive=TRUE, showWarnings=FALSE)

write.table(as.data.frame(go_base),   "Results/clusterProfiler/GO_baseline.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(kegg_base), "Results/clusterProfiler/KEGG_baseline.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(as.data.frame(go_follow),   "Results/clusterProfiler/GO_followup.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(kegg_follow), "Results/clusterProfiler/KEGG_followup.txt", sep="\t", quote=FALSE, row.names=FALSE)

############################################
# END OF SCRIPT
############################################
