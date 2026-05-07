############################################
# CKD EWAS – Functional Enrichment
# Analysis: GO and KEGG enrichment (missMethyl)
#
# Author: Muhulo Muhau Mungamba
# Date: 2026-04-15
#
# Description:
# Bias-corrected pathway enrichment accounting for CpG probe density.
# Primary enrichment analysis.
############################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

############################################
# Load libraries
############################################

library(missMethyl)
library(dplyr)

############################################
# Load DMP results
############################################

ckd_base  <- read.delim("Results/Timepoint_1/DMPs/CKD/TopTable_CKD_T1.txt")
egfr_base <- read.delim("Results/Timepoint_1/DMPs/Delta_eGFR/TopTable_Delta_eGFR_T1.txt")
acr_base  <- read.delim("Results/Timepoint_1/DMPs/Delta_ACR/TopTable_Delta_ACR_T1.txt")

ckd_follow  <- read.delim("Results/Timepoint_2/DMPs/CKD/TopTable_CKD_T2.txt")
egfr_follow <- read.delim("Results/Timepoint_2/DMPs/Delta_eGFR/TopTable_Delta_eGFR_T2.txt")
acr_follow  <- read.delim("Results/Timepoint_2/DMPs/Delta_ACR/TopTable_Delta_ACR_T2.txt")

############################################
# Combine CpGs
############################################

combine_cpgs <- function(...) {
  bind_rows(...) %>%
    filter(!is.na(CpG)) %>%
    distinct(CpG, .keep_all=TRUE) %>%
    arrange(BaconPvalue_Mvalues)
}

combined_base   <- combine_cpgs(ckd_base, egfr_base, acr_base)
combined_follow <- combine_cpgs(ckd_follow, egfr_follow, acr_follow)

############################################
# Define CpG sets
############################################

all_base  <- combined_base$CpG
top1000_b <- head(all_base, 1000)
top5000_b <- head(all_base, 5000)

all_follow  <- combined_follow$CpG
top1000_f <- head(all_follow, 1000)
top5000_f <- head(all_follow, 5000)

############################################
# Run enrichment
############################################

go_1000_base <- gometh(top1000_b, all_base, "GO", array.type="EPIC")
go_5000_base <- gometh(top5000_b, all_base, "GO", array.type="EPIC")

kegg_1000_base <- gometh(top1000_b, all_base, "KEGG", array.type="EPIC")
kegg_5000_base <- gometh(top5000_b, all_base, "KEGG", array.type="EPIC")

go_1000_follow <- gometh(top1000_f, all_follow, "GO", array.type="EPIC")
go_5000_follow <- gometh(top5000_f, all_follow, "GO", array.type="EPIC")

kegg_1000_follow <- gometh(top1000_f, all_follow, "KEGG", array.type="EPIC")
kegg_5000_follow <- gometh(top5000_f, all_follow, "KEGG", array.type="EPIC")

############################################
# Save results
############################################

dir.create("Results/missMethyl", recursive=TRUE, showWarnings=FALSE)

write.table(go_1000_base, "Results/missMethyl/GO_top1000_baseline.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(go_5000_base, "Results/missMethyl/GO_top5000_baseline.txt", sep="\t", quote=FALSE, row.names=FALSE)

write.table(go_1000_follow, "Results/missMethyl/GO_top1000_followup.txt", sep="\t", quote=FALSE, row.names=FALSE)
write.table(go_5000_follow, "Results/missMethyl/GO_top5000_followup.txt", sep="\t", quote=FALSE, row.names=FALSE)

############################################
# END OF SCRIPT
############################################