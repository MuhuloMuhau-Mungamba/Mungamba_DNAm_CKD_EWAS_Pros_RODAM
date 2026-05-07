############################################
# CKD EWAS – Sensitivity Analysis (Violin Plots)
# 07g_ViolinPlots_Sensitivity_CKD_T1_T2.R
#
# Author: Muhulo Muhau Mungamba & Venema Andrea
# Date: 2025-09-30
#
# Purpose:
# - Visualize methylation levels for significant CpGs
# - Stratified by site and CKD status
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

# CpGs of interest
cpgs <- c("cg26122413")

# Input files
beta_file <- "path/to/beta_matrix.txt"
pheno_T1  <- "path/to/phenotype_T1.csv"
pheno_T2  <- "path/to/phenotype_T2.csv"

############################
# FUNCTION: Prepare data
############################

prepare_violin_data <- function(beta, pheno, cpgs, timepoint) {
  
  # Subset CpGs
  beta_sub <- beta[cpgs, ]
  
  # Reshape
  beta_tidy <- beta_sub %>%
    rownames_to_column("CpG") %>%
    pivot_longer(-CpG, names_to = "Basename", values_to = "Beta") %>%
    pivot_wider(names_from = CpG, values_from = Beta)
  
  # Match Basename format
  pheno$Basename <- paste0("X", pheno$Basename)
  
  # Merge
  df <- merge(pheno, beta_tidy, by = "Basename", all.x = TRUE)
  
  ############################
  # Create variables
  ############################
  
  if (timepoint == "T1") {
    df$Site <- factor(df$R1_Site,
                      levels = c(1, 2, 3),
                      labels = c("Rural Ghana", "Urban Ghana", "Amsterdam Ghana"))
    df$CKD  <- df$defCKD_R1
  } else {
    df$Site <- factor(df$R2_Site,
                      levels = c(1, 2, 3),
                      labels = c("Rural Ghana", "Urban Ghana", "Amsterdam Ghana"))
    df$CKD  <- df$defCKD_R2
  }
  
  df$Group <- paste(df$Site,
                    ifelse(df$CKD == 1, "Case", "Control"),
                    sep = " - ")
  
  return(df)
}

############################
# FUNCTION: Plot violin
############################

plot_violin <- function(df, cpg, outfile) {
  
  df$Beta_percent <- df[[cpg]] * 100
  
  p <- ggplot(df, aes(x = Group, y = Beta_percent, fill = Group)) +
    geom_violin(trim = FALSE) +
    stat_summary(fun = median, geom = "point",
                 shape = 16, size = 3, color = "black") +
    theme_minimal() +
    labs(
      x = NULL,
      y = paste0("Beta-value (%) - ", cpg)
    ) +
    theme(
      text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  ggsave(outfile, p, width = 10, height = 6, dpi = 300)
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

df_T1 <- prepare_violin_data(beta, pheno1, cpgs, "T1")

# Plot all CpGs
for (cpg in cpgs) {
  
  plot_violin(
    df = df_T1,
    cpg = cpg,
    outfile = paste0("Violin_", cpg, "_T1_Sensitivity.tiff")
  )
}

############################
# =========================
# FOLLOW-UP (T2)
# =========================
############################

pheno2 <- read.csv(pheno_T2)

df_T2 <- prepare_violin_data(beta, pheno2, cpgs, "T2")

# Plot all CpGs
for (cpg in cpgs) {
  
  plot_violin(
    df = df_T2,
    cpg = cpg,
    outfile = paste0("Violin_", cpg, "_T2_Sensitivity.tiff")
  )
}

############################
# OPTIONAL: Kruskal-Wallis
############################

run_kruskal <- function(df, cpg) {
  
  df$Beta_percent <- df[[cpg]] * 100
  
  df$Group <- factor(df$Group)
  
  test <- kruskal.test(Beta_percent ~ Group, data = df)
  
  cat("\nKruskal-Wallis for", cpg, "\n")
  print(test)
}

# Example
for (cpg in cpgs) {
  run_kruskal(df_T2, cpg)
}

############################################
# END OF SCRIPT
############################################