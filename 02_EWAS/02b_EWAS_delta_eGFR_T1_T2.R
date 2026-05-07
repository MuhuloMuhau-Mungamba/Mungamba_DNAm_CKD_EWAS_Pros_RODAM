############################################
# Pros-RODAM CKD EWAS
# 02b_EWAS_delta_eGFR_T1_T2.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
#
# EWAS of change in kidney function (ΔeGFR)
# - Baseline (T1)
# - Follow-up (T2)
#
# Study design:
# Nested case-control (n = 274)
#
# Model:
# ΔeGFR ~ age + sex + site + cell counts + plate
#
# Notes:
# - Uses normalised methylation data from 01c
# - Uses phenotype + cell counts from 01b
############################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

############################
# Load required packages
############################

library(limma)
library(dplyr)
library(bacon)

############################################################
######################## BASELINE (T1) ######################
############################################################

############################
# File paths
############################

beta_file  <- "CKD_EWAS_T1_Beta.txt"
mval_file  <- "CKD_EWAS_T1_Mval.txt"
annot_file <- "CKD_annotation.txt"
pheno_file <- "CKD_EWAS_T1_phenofile_wCells.csv"

############################
# Load data
############################

beta  <- read.table(beta_file, sep = "\t", header = TRUE, check.names = FALSE)
mval  <- read.table(mval_file, sep = "\t", header = TRUE, check.names = FALSE)
annot <- read.table(annot_file, sep = "\t", header = TRUE)

targets <- read.csv(pheno_file)
head(targets)

############################
# Define variables
############################

delta_eGFR <- targets$delta_eGFR

sex  <- factor(targets$R1_Sex)
age  <- targets$R1_Age
site <- factor(targets$R1_Site)

plate <- factor(targets$plate)

CD8T <- targets$CD8T
CD4T <- targets$CD4T
NK   <- targets$NK
Bcell<- targets$Bcell
Mono <- targets$Mono
Neu  <- targets$Neu

############################
# Design matrices
############################

design_full <- model.matrix(
  ~ delta_eGFR + sex + age +
    CD8T + CD4T + NK + Bcell + Mono + Neu +
    site + plate
)

design_base <- model.matrix(~ delta_eGFR)

############################
# EWAS: M-values
############################

fit_m <- lmFit(mval, design_full)
fit_m <- eBayes(fit_m)

fit_m0 <- lmFit(mval, design_base)
fit_m0 <- eBayes(fit_m0)

top_m_fdr <- topTable(
  fit_m,
  coef = "delta_eGFR",
  number = nrow(mval),
  adjust = "fdr"
)

top_m_raw <- topTable(
  fit_m0,
  coef = "delta_eGFR",
  number = nrow(mval),
  adjust = "none"
)

############################
# Bacon correction
############################

tstats <- fit_m$t[, "delta_eGFR"]
bc <- bacon(tstats)

bc_p   <- pval(bc)
bc_fdr <- p.adjust(bc_p, method = "fdr")

############################
# Add annotation
############################

annot_match <- annot[match(rownames(top_m_fdr), rownames(annot)), ]

top_m_out <- cbind(
  top_m_fdr,
  DeltaMval = top_m_raw$logFC,
  BaconPvalue = bc_p,
  BaconFDR = bc_fdr
)

############################
# EWAS: Beta values
############################

fit_b <- lmFit(beta, design_full)
fit_b <- eBayes(fit_b)

fit_b0 <- lmFit(beta, design_base)
fit_b0 <- eBayes(fit_b0)

top_b_fdr <- topTable(
  fit_b,
  coef = "delta_eGFR",
  number = nrow(beta),
  adjust = "fdr"
)

top_b_raw <- topTable(
  fit_b0,
  coef = "delta_eGFR",
  number = nrow(beta),
  adjust = "none"
)

top_b_out <- cbind(
  top_b_fdr,
  DeltaBeta = top_b_raw$logFC
)

############################
# Combine results
############################

top_b_out <- top_b_out[match(rownames(top_m_out), rownames(top_b_out)), ]

final_T1 <- cbind(
  top_m_out,
  top_b_out,
  annot_match
)

############################
# Save output
############################

write.table(
  final_T1,
  "CKD_EWAS_T1_delta_eGFR_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)

############################################################
######################## FOLLOW-UP (T2) #####################
############################################################

############################
# File paths
############################

beta_file  <- "CKD_EWAS_T2_Beta.txt"
mval_file  <- "CKD_EWAS_T2_Mval.txt"
pheno_file <- "CKD_EWAS_T2_phenofile_wCells.csv"

############################
# Load data
############################

beta  <- read.table(beta_file, sep = "\t", header = TRUE, check.names = FALSE)
mval  <- read.table(mval_file, sep = "\t", header = TRUE, check.names = FALSE)

targets <- read.csv(pheno_file)
head(targets)

############################
# Define variables
############################

delta_eGFR <- targets$delta_eGFR

sex  <- factor(targets$R2_Sex)
age  <- targets$R2_Age
site <- factor(targets$R2_Site)

plate <- factor(targets$plate)

CD8T <- targets$CD8T
CD4T <- targets$CD4T
NK   <- targets$NK
Bcell<- targets$Bcell
Mono <- targets$Mono
Neu  <- targets$Neu

############################
# Design matrices
############################

design_full <- model.matrix(
  ~ delta_eGFR + sex + age +
    CD8T + CD4T + NK + Bcell + Mono + Neu +
    site + plate
)

design_base <- model.matrix(~ delta_eGFR)

############################
# EWAS: M-values
############################

fit_m <- lmFit(mval, design_full)
fit_m <- eBayes(fit_m)

fit_m0 <- lmFit(mval, design_base)
fit_m0 <- eBayes(fit_m0)

top_m_fdr <- topTable(
  fit_m,
  coef = "delta_eGFR",
  number = nrow(mval),
  adjust = "fdr"
)

top_m_raw <- topTable(
  fit_m0,
  coef = "delta_eGFR",
  number = nrow(mval),
  adjust = "none"
)

############################
# Bacon correction
############################

tstats <- fit_m$t[, "delta_eGFR"]
bc <- bacon(tstats)

bc_p   <- pval(bc)
bc_fdr <- p.adjust(bc_p, method = "fdr")

top_m_out <- cbind(
  top_m_fdr,
  DeltaMval = top_m_raw$logFC,
  BaconPvalue = bc_p,
  BaconFDR = bc_fdr
)

############################
# EWAS: Beta values
############################

fit_b <- lmFit(beta, design_full)
fit_b <- eBayes(fit_b)

fit_b0 <- lmFit(beta, design_base)
fit_b0 <- eBayes(fit_b0)

top_b_fdr <- topTable(
  fit_b,
  coef = "delta_eGFR",
  number = nrow(beta),
  adjust = "fdr"
)

top_b_raw <- topTable(
  fit_b0,
  coef = "delta_eGFR",
  number = nrow(beta),
  adjust = "none"
)

top_b_out <- cbind(
  top_b_fdr,
  DeltaBeta = top_b_raw$logFC
)

############################
# Combine results
############################

top_b_out <- top_b_out[match(rownames(top_m_out), rownames(top_b_out)), ]

final_T2 <- cbind(
  top_m_out,
  top_b_out,
  annot_match
)

############################
# Save output
############################

write.table(
  final_T2,
  "CKD_EWAS_T2_delta_eGFR_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)

############################################################
# END OF SCRIPT
############################################################