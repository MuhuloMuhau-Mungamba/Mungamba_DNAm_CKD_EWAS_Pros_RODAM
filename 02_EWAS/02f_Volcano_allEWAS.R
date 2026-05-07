############################################
# Pros-RODAM CKD EWAS
# 02f_Volcano_allEWAS.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
# Date: 2025-08-12
#
# Purpose:
# - Generate volcano plots for all EWAS results
############################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

############################
# Load packages
############################

library(ggplot2)

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
# Volcano plots
############################

for (name in names(files)) {
  
  df <- read.table(files[[name]], header = TRUE, sep = "\t")
  
  volplot <- data.frame(
    logFC = df$DeltaBeta,
    pval  = df$P.Value
  )
  
  volplot$threshold <- ifelse(
    volplot$pval < 1e-6 & abs(volplot$logFC) > 0.05,
    "Significant", "Not significant"
  )
  
  volplot$threshold <- factor(volplot$threshold,
                              levels = c("Significant", "Not significant"))
  
  p <- ggplot(volplot, aes(x = logFC, y = -log10(pval), color = threshold)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_hline(yintercept = -log10(1e-6), linetype = "dashed") +
    geom_vline(xintercept = c(-0.05, 0.05), linetype = "dashed") +
    theme_bw() +
    labs(
      title = paste("Volcano plot:", name),
      x = "Effect size (Delta Beta)",
      y = "-log10(p-value)"
    ) +
    theme(legend.title = element_blank())
  
  ggsave(
    filename = paste0("Volcano_", name, ".png"),
    plot = p,
    width = 6,
    height = 5
  )
}