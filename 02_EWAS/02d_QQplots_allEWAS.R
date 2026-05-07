############################################
# Pros-RODAM CKD EWAS
# 02d_QQplots_allEWAS.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
#
#
# Purpose:
# - Generate QQ plots for EWAS results
# - Compute genomic inflation factor (lambda)
#
# Outcomes:
# - Incident CKD
# - ΔeGFR
# - ΔACR
# At:
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
  
  print(paste("Lambda for", title, "=", round(lambda, 3)))
  
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
# Run QQ plots
############################

for (name in names(files)) {
  
  df <- read.table(files[[name]], header = TRUE, sep = "\t")
  
  qq_plot_lambda(
    pvalues = df$P.Value,
    title   = paste("QQ plot:", name),
    outfile = paste0("QQ_", name, ".pdf")
  )
}