############################################
# CKD EWAS – Sensitivity Analysis (ΔACR)
# 07c_Sensitivity_EWAS_Delta_ACR_T1_T2.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
# Date: 2025-09-26
#
# Description:
# Sensitivity EWAS for ΔACR (log-transformed)
# Adjusted for:
# age, sex, site, blood cell composition, plate,
# hypertension (HTN), diabetes (DM)
#
# Structure:
# 1. Baseline (T1)
# 2. Follow-up (T2)
############################################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

############################
# Libraries
############################

library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(minfiData)
library(limma)
library(affy)
library(dplyr)
library(bacon)

############################
# Paths (EDIT BEFORE RUNNING)
############################

idatspath  <- "path/to/idat_files/"
cross_file <- "path/to/cross_reactive_probes.csv"

############################
# =========================
# BASELINE (T1)
# =========================
############################

pheno_T1 <- "path/to/phenotype_T1_withcells.csv"
out_T1   <- "2025.09.26_Sensitivity_EWAS_Delta_ACR_T1.txt"

targets <- read.csv(pheno_T1)

############################
# Methylation preprocessing
############################

RGset <- read.metharray.exp(base = idatspath, targets = targets, force = TRUE)

GMset <- preprocessFunnorm(RGset, nPCs = 2, bgCorr = TRUE, dyeCorr = TRUE)

annotation <- getAnnotation(GMset)

# Remove sex chromosomes
sex_idx <- which(annotation$chr %in% c("chrX", "chrY"))
GMset <- GMset[-sex_idx, ]

# Remove cross-reactive probes
probes_remove <- read.csv(cross_file, header = FALSE)[,1]
GMset <- GMset[!featureNames(GMset) %in% probes_remove, ]

# Extract M-values and Beta-values
mval <- getM(GMset)
beta <- getBeta(GMset)

# Remove infinite values
inf_idx <- which(mval == -Inf, arr.ind = TRUE)
bad <- unique(rownames(inf_idx))

mval <- mval[!rownames(mval) %in% bad, ]
beta <- beta[!rownames(beta) %in% bad, ]

############################
# Variables (T1)
############################

deltaacr <- targets$delta_LogACR

DM  <- factor(targets$R1_DM_Dichot)
HTN <- factor(targets$R1_HTN_dichot)

sex  <- factor(targets$R1_Sex)
site <- factor(targets$R1_Site)
age  <- targets$R1_Age

CD8T  <- targets$CD8T
CD4T  <- targets$CD4T
NK    <- targets$NK
Bcell <- targets$Bcell
Mono  <- targets$Mono
Neu   <- targets$Neu
plate <- targets$plate

############################
# Model (T1)
############################

design <- model.matrix(
  ~ deltaacr + sex + age +
    CD8T + CD4T + NK + Bcell + Mono + Neu +
    site + plate + DM + HTN
)

############################
# EWAS M-values (T1)
############################

fit <- lmFit(mval, design)
fit <- eBayes(fit)

top_m <- topTable(fit, coef = "deltaacr", number = nrow(mval), adjust = "fdr")

############################
# Bacon correction (T1)
############################

tstats <- fit$t[, "deltaacr"]
tstats_bc <- bacon(tstats)

bc_p <- pval(tstats_bc)
bc_fdr <- p.adjust(bc_p, method = "fdr")

detach("package:bacon", unload = TRUE)

############################
# EWAS Beta-values (T1)
############################

fit_b <- lmFit(beta, design)
fit_b <- eBayes(fit_b)

top_b <- topTable(fit_b, coef = "deltaacr", number = nrow(beta), adjust = "fdr")

############################
# Combine results (T1)
############################

top_b <- top_b[match(rownames(top_m), rownames(top_b)), ]

final_T1 <- cbind(
  top_m,
  DeltaMval = top_m$logFC,
  BaconPvalue = bc_p[match(rownames(top_m), names(bc_p))],
  BaconFDR = bc_fdr[match(rownames(top_m), names(bc_fdr))],
  top_b,
  DeltaBeta = top_b$logFC
)

write.table(final_T1, out_T1, sep = "\t", quote = FALSE)

############################
# =========================
# FOLLOW-UP (T2)
# =========================
############################

pheno_T2 <- "path/to/phenotype_T2_withcells.csv"
out_T2   <- "2025.09.26_Sensitivity_EWAS_Delta_ACR_T2.txt"

targets <- read.csv(pheno_T2)

############################
# Preprocessing (T2)
############################

RGset <- read.metharray.exp(base = idatspath, targets = targets, force = TRUE)

GMset <- preprocessFunnorm(RGset, nPCs = 2)

annotation <- getAnnotation(GMset)

sex_idx <- which(annotation$chr %in% c("chrX", "chrY"))
GMset <- GMset[-sex_idx, ]

probes_remove <- read.csv(cross_file, header = FALSE)[,1]
GMset <- GMset[!featureNames(GMset) %in% probes_remove, ]

mval <- getM(GMset)
beta <- getBeta(GMset)

inf_idx <- which(mval == -Inf, arr.ind = TRUE)
bad <- unique(rownames(inf_idx))

mval <- mval[!rownames(mval) %in% bad, ]
beta <- beta[!rownames(beta) %in% bad, ]

############################
# Variables (T2)
############################

deltaacr <- targets$delta_LogACR

DM  <- factor(targets$R2_DM_Dichot)
HTN <- factor(targets$R2_HTN_dichot)

sex  <- factor(targets$R2_Sex)
site <- factor(targets$R2_Site)
age  <- targets$R2_Age

CD8T  <- targets$CD8T
CD4T  <- targets$CD4T
NK    <- targets$NK
Bcell <- targets$Bcell
Mono  <- targets$Mono
Neu   <- targets$Neu
plate <- targets$plate

############################
# Model + EWAS (T2)
############################

design <- model.matrix(
  ~ deltaacr + sex + age +
    CD8T + CD4T + NK + Bcell + Mono + Neu +
    site + plate + DM + HTN
)

fit <- lmFit(mval, design)
fit <- eBayes(fit)

top_m <- topTable(fit, coef = "deltaacr", number = nrow(mval), adjust = "fdr")

############################
# Bacon correction (T2)
############################

tstats <- fit$t[, "deltaacr"]
tstats_bc <- bacon(tstats)

bc_p <- pval(tstats_bc)
bc_fdr <- p.adjust(bc_p, method = "fdr")

detach("package:bacon", unload = TRUE)

############################
# Beta-values (T2)
############################

fit_b <- lmFit(beta, design)
fit_b <- eBayes(fit_b)

top_b <- topTable(fit_b, coef = "deltaacr", number = nrow(beta), adjust = "fdr")

############################
# Combine results (T2)
############################

top_b <- top_b[match(rownames(top_m), rownames(top_b)), ]

final_T2 <- cbind(
  top_m,
  DeltaMval = top_m$logFC,
  BaconPvalue = bc_p[match(rownames(top_m), names(bc_p))],
  BaconFDR = bc_fdr[match(rownames(top_m), names(bc_fdr))],
  top_b,
  DeltaBeta = top_b$logFC
)

write.table(final_T2, out_T2, sep = "\t", quote = FALSE)

############################################
# END OF SCRIPT
############################################