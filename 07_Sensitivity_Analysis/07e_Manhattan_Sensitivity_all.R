############################################
# Pros-RODAM CKD EWAS – Sensitivity Analysis
# 07e_Manhattan_Sensitivity_all.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
#
# Purpose:
# - Generate Manhattan plots for sensitivity EWAS results
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
# Manhattan plot function
############################

manhattan_plot <- function(df, p_col, title, outfile, p_cutoff = 1e-7) {
  
  ############################
  # Prepare data
  ############################
  
  df <- df %>%
    mutate(
      chr = as.numeric(gsub("chr", "", chr)),
      pos = as.numeric(pos)
    ) %>%
    filter(!is.na(chr), !is.na(pos), chr %in% 1:22, pos > 0)
  
  ############################
  # Create coordinates
  ############################
  
  cols <- c("seagreen3", "skyblue3")
  coor <- c()
  col_vec <- c()
  ticks <- c()
  pvals <- c()
  
  last_base <- 0
  
  for (i in 1:22) {
    chr_data <- df[df$chr == i, ]
    
    if (nrow(chr_data) == 0) next
    
    pvals <- c(pvals, chr_data[[p_col]])
    
    coor_chr <- chr_data$pos + last_base
    coor <- c(coor, coor_chr)
    
    col_vec <- c(col_vec, rep(cols[i %% 2 + 1], nrow(chr_data)))
    
    ticks <- c(ticks, mean(range(coor_chr)))
    
    last_base <- max(coor_chr)
  }
  
  ############################
  # Plot
  ############################
  
  png(outfile, width = 1600, height = 600, res = 120)
  
  plot(coor,
       -log10(pvals),
       col = col_vec,
       pch = 20,
       xlab = "Chromosome",
       ylab = "-log10(P-value)",
       main = title,
       xaxt = "n")
  
  axis(1, ticks, 1:22)
  
  abline(h = -log10(p_cutoff), col = "red", lty = 2)
  
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
# Run Manhattan plots (raw)
############################

for (name in names(files)) {
  
  df <- read.table(files[[name]], header = TRUE, sep = "\t")
  
  manhattan_plot(
    df      = df,
    p_col   = "P.Value",
    title   = paste("Manhattan (Sensitivity):", name),
    outfile = paste0("Manhattan_", name, "_Sensitivity.png")
  )
}

############################
# Run Manhattan plots (Bacon)
############################

for (name in names(files)) {
  
  df <- read.table(files[[name]], header = TRUE, sep = "\t")
  
  manhattan_plot(
    df      = df,
    p_col   = "BaconPvalue_Mvalues",
    title   = paste("Manhattan Bacon (Sensitivity):", name),
    outfile = paste0("Manhattan_Bacon_", name, "_Sensitivity.png")
  )
}

############################################
# END OF SCRIPT
############################################