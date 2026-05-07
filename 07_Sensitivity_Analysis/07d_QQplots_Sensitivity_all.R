############################################
# Pros-RODAM CKD EWAS – Sensitivity Analysis
# 07d_QQplots_Sensitivity_all.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
#
# Purpose:
# - Generate QQ plots for sensitivity EWAS results
# - Compute genomic inflation factor (lambda)
#
# Models:
# - Adjusted for DM + HTN
#
# Outcomes:
# - Incident CKD
# - ΔeGFR
# - ΔACR
#
# Timepoints:
# - Baseline (T1)
# - Follow-up (T2)
############################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

############################
# Load packages
############################

library(dplyr)

############################
# QQ plot function
############################

qq_plot_lambda <- function(pvalues, title, outfile) {
  
  # Clean p-values
  pvalues <- pvalues[is.finite(pvalues) & pvalues > 0 & pvalues < 1]
  
  # Observed vs expected
  observed <- -log10(sort(pvalues))
  expected <- -log10(ppoints(length(pvalues)))
  
  ############################
  # Lambda calculation
  ############################
  
  chisq <- qchisq(pvalues, df = 1, lower.tail = FALSE)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
  
  cat("Lambda for", title, "=", round(lambda, 3), "\n")
  
  ############################
  # Plot
  ############################
  
  pdf(outfile, width = 6, height = 6)
  
  plot(expected, observed,
       pch = 20,
       main = paste0(title, "\nLambda = ", round(lambda, 3)),
       xlab = expression(Expected~-log[10](italic(p))),
       ylab = expression(Observed~-log[10](italic(p))))
  
  abline(0, 1, col = "red", lwd = 2)
  
  dev.off()
}

############################
# Input files (EDIT PATHS)
############################

files <- list(
  "CKD_T1"  = "path/to/CKD_Sensitivity_T1.txt",
  "CKD_T2"  = "path/to/CKD_Sensitivity_T2.txt",
  "eGFR_T1" = "path/to/Delta_eGFR_Sensitivity_T1.txt",
  "eGFR_T2" = "path/to/Delta_eGFR_Sensitivity_T2.txt",
  "ACR_T1"  = "path/to/Delta_ACR_Sensitivity_T1.txt",
  "ACR_T2"  = "path/to/Delta_ACR_Sensitivity_T2.txt"
)

############################
# Run QQ plots (raw p-values)
############################

for (name in names(files)) {
  
  df <- read.table(files[[name]], header = TRUE, sep = "\t")
  
  qq_plot_lambda(
    pvalues = df$P.Value,
    title   = paste("QQ plot (Sensitivity):", name),
    outfile = paste0("QQ_", name, "_Sensitivity.pdf")
  )
}

############################
# Run QQ plots (Bacon-corrected)
############################

for (name in names(files)) {
  
  df <- read.table(files[[name]], header = TRUE, sep = "\t")
  
  qq_plot_lambda(
    pvalues = df$BaconPvalue_Mvalues,
    title   = paste("QQ plot Bacon (Sensitivity):", name),
    outfile = paste0("QQ_Bacon_", name, "_Sensitivity.pdf")
  )
}

############################################
# END OF SCRIPT
############################################