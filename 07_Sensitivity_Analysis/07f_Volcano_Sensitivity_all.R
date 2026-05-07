############################################
# Pros-RODAM CKD EWAS – Sensitivity Analysis
# 07f_Volcano_Sensitivity_all.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
#
# Purpose:
# - Generate Volcano plots for sensitivity EWAS results
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
library(ggplot2)

############################
# Volcano plot function
############################

volcano_plot <- function(df, title, outfile,
                         p_thresh = 1e-6,
                         effect_thresh = 0.05) {
  
  ############################
  # Prepare data
  ############################
  
  df <- df %>%
    filter(is.finite(P.Value),
           is.finite(logFCBeta)) %>%
    mutate(
      logP = -log10(P.Value),
      threshold = ifelse(
        P.Value < p_thresh & abs(logFCBeta) > effect_thresh,
        "Significant",
        "Not significant"
      )
    )
  
  df$threshold <- factor(df$threshold,
                         levels = c("Significant", "Not significant"))
  
  ############################
  # Plot
  ############################
  
  jpeg(outfile, width = 1600, height = 1200, res = 150)
  
  print(
    ggplot(df, aes(x = logFCBeta, y = logP, color = threshold)) +
      geom_point(alpha = 0.5, size = 1.2) +
      
      geom_vline(xintercept = c(-effect_thresh, effect_thresh),
                 linetype = "dashed") +
      
      geom_hline(yintercept = -log10(p_thresh),
                 linetype = "dashed") +
      
      labs(
        title = title,
        x = "Effect size (Beta)",
        y = expression(-log[10](p))
      ) +
      
      scale_color_manual(values = c(
        "Significant" = "#B41A21",
        "Not significant" = "grey70"
      )) +
      
      theme_bw() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.title = element_blank(),
        legend.position = "bottom"
      )
  )
  
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
# Run Volcano plots
############################

for (name in names(files)) {
  
  df <- read.table(files[[name]], header = TRUE, sep = "\t")
  
  volcano_plot(
    df      = df,
    title   = paste("Volcano (Sensitivity):", name),
    outfile = paste0("Volcano_", name, "_Sensitivity.jpg")
  )
}

############################################
# END OF SCRIPT
############################################