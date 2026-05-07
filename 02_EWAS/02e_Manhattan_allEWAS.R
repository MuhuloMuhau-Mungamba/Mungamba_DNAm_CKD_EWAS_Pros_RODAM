############################################
# Pros-RODAM CKD EWAS
# 02e_Manhattan_allEWAS.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
# Date: 2025-08-12
#
# Purpose:
# - Generate Manhattan plots for all EWAS results
############################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

############################
# Load packages
############################

library(qqman)

############################
# Input files
############################

files <- list(
  "CKD_T1"  = "CKD_EWAS_T1_results.txt",
  "CKD_T2"  = "CKD_EWAS_T2_results.txt",
  "eGFR_T1" = "CKD_EWAS_T1_delta_eGFR_results.txt",
  "eGFR_T2" = "CKD_EWAS_T2_delta_eGFR_results.txt",
  "ACR_T1"  = "CKD_EWAS_T1_delta_ACR_results.txt",
  "ACR_T2"  = "CKD_EWAS_T2_delta_ACR_results.txt"
)

############################
# Manhattan plotting
############################

for (name in names(files)) {
  
  df <- read.table(files[[name]], header = TRUE, sep = "\t")
  
  # Ensure required columns exist
  df$CHR <- as.numeric(gsub("chr", "", df$chr))
  df$BP  <- as.numeric(df$pos)
  df$P   <- df$P.Value
  df$SNP <- rownames(df)
  
  png(paste0("Manhattan_", name, ".png"), width = 1600, height = 600)
  
  manhattan(df,
            chr = "CHR",
            bp  = "BP",
            p   = "P",
            snp = "SNP",
            main = paste("Manhattan plot:", name),
            genomewideline = -log10(1e-7),
            suggestiveline = -log10(1e-5))
  
  dev.off()
}