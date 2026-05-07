############################################
# CKD EWAS – Transcription Factor Enrichment
#
# Author: Muhulo Muhau Mungamba
# Date: 2026-04-16
############################################

rm(list = ls())
gc()

library(enrichR)
library(dplyr)
library(stringr)

############################################
# Load combined CpG tables
############################################

combined_base   <- read.delim("Results/missMethyl/combined_base.txt")
combined_follow <- read.delim("Results/missMethyl/combined_follow.txt")

############################################
# Extract genes
############################################

extract_genes <- function(df){
  df %>%
    filter(BaconPvalue_Mvalues < 1e-3,
           !is.na(UCSC_RefGene_Name),
           UCSC_RefGene_Name != "") %>%
    pull(UCSC_RefGene_Name) %>%
    str_split(";") %>%
    unlist() %>%
    unique()
}

genes_base   <- extract_genes(combined_base)
genes_follow <- extract_genes(combined_follow)

############################################
# Run enrichment
############################################

dbs <- c("ChEA_2016","ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")

tf_base   <- enrichr(genes_base, dbs)
tf_follow <- enrichr(genes_follow, dbs)

############################################
# Save results
############################################

dir.create("Results/TF_enrichment", recursive=TRUE, showWarnings=FALSE)

write.table(tf_base[[1]], "Results/TF_enrichment/TF_baseline.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(tf_follow[[1]], "Results/TF_enrichment/TF_followup.txt", sep="\t", quote=FALSE, row.names=FALSE)

############################################
# END OF SCRIPT
############################################