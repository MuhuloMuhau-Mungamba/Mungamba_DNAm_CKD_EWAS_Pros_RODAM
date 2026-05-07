############################################
# CKD EWAS – Sensitivity Analysis (Violin Plots ΔACR)
# 07i_ViolinPlots_Sensitivity_DeltaACR_T1_T2.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
# Date: 2025-09-30
#
# Purpose:
# - Visualize methylation levels for ΔACR-associated CpGs
# - Stratified by site
# - For baseline (T1) and follow-up (T2)
############################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

############################
# Load packages
############################

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

############################
# USER INPUT (EDIT)
############################

# Significant CpGs (ΔACR)
cpgs <- c("cg05215481", "cg04646451", "cg15646763")

# Input files (generalised paths)
beta_file <- "path/to/beta_matrix.txt"
pheno_T1  <- "path/to/phenotype_T1.csv"
pheno_T2  <- "path/to/phenotype_T2.csv"

############################
# FUNCTION: Prepare data
############################

prepare_data <- function(beta, pheno, cpgs, timepoint) {
  
  # Subset CpGs
  beta_sub <- beta[cpgs, ]
  
  # Reshape to long → wide
  beta_tidy <- beta_sub %>%
    rownames_to_column("CpG") %>%
    pivot_longer(-CpG, names_to = "Basename", values_to = "Beta") %>%
    pivot_wider(names_from = CpG, values_from = Beta)
  
  # Match Basename format
  pheno$Basename <- paste0("X", pheno$Basename)
  
  # Merge
  df <- merge(pheno, beta_tidy, by = "Basename", all.x = TRUE)
  
  # Site variable (T1 vs T2)
  if (timepoint == "T1") {
    df$Site <- factor(df$R1_Site,
                      levels = c(1,2,3),
                      labels = c("Rural Ghana", "Urban Ghana", "Amsterdam Ghana"))
  } else {
    df$Site <- factor(df$R2_Site,
                      levels = c(1,2,3),
                      labels = c("Rural Ghana", "Urban Ghana", "Amsterdam Ghana"))
  }
  
  return(df)
}

############################
# FUNCTION: Plot violin
############################

plot_violin <- function(df, cpg, outfile) {
  
  df$Beta_percent <- df[[cpg]] * 100
  
  p <- ggplot(df, aes(x = Site, y = Beta_percent, fill = Site)) +
    geom_violin(trim = FALSE) +
    stat_summary(fun = median, geom = "point",
                 shape = 16, size = 3, color = "black") +
    theme_minimal() +
    labs(
      x = "Site",
      y = paste0("Beta-value (%) - ", cpg)
    ) +
    theme(
      text = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      legend.position = "none"
    )
  
  ggsave(outfile, p, width = 10, height = 6, dpi = 300)
}

############################
# FUNCTION: Kruskal test
############################

run_kruskal <- function(df, cpg) {
  
  df$Beta_percent <- df[[cpg]] * 100
  
  test <- kruskal.test(Beta_percent ~ Site, data = df)
  
  cat("\nKruskal-Wallis for", cpg, "\n")
  print(test)
}

############################
# LOAD DATA
############################

beta <- read.table(beta_file, header = TRUE, sep = "\t", check.names = FALSE)

############################
# =========================
# BASELINE (T1)
# =========================
############################

pheno1 <- read.csv(pheno_T1)

df_T1 <- prepare_data(beta, pheno1, cpgs, "T1")

for (cpg in cpgs) {
  
  plot_violin(
    df_T1,
    cpg,
    paste0("Violin_", cpg, "_DeltaACR_T1_Sensitivity.tiff")
  )
  
  run_kruskal(df_T1, cpg)
}

############################
# =========================
# FOLLOW-UP (T2)
# =========================
############################

pheno2 <- read.csv(pheno_T2)

df_T2 <- prepare_data(beta, pheno2, cpgs, "T2")

for (cpg in cpgs) {
  
  plot_violin(
    df_T2,
    cpg,
    paste0("Violin_", cpg, "_DeltaACR_T2_Sensitivity.tiff")
  )
  
  run_kruskal(df_T2, cpg)
}

############################################
# END OF SCRIPT
############################################