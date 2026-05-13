############################################
# CKD EWAS – eFORGE Input Preparation
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
# Date: 2026-04-16
############################################

rm(list = ls())
gc()

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
# Extract top CpGs
############################################

top1000_base   <- head(combined_base$CpG, 1000)
top1000_follow <- head(combined_follow$CpG, 1000)

############################################
# Save files
############################################

dir.create("Results/eFORGE", recursive=TRUE, showWarnings=FALSE)

write.table(top1000_base,
            "Results/eFORGE/top1000_baseline.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(top1000_follow,
            "Results/eFORGE/top1000_followup.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

############################################
# END OF SCRIPT
############################################
