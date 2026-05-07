############################################
# CKD EWAS – Sensitivity Analysis
# 07a_Sensitivity_EWAS_CKD_T1_T2.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
# Date: 2025-09-26
#
# Description:
# Sensitivity EWAS for incident CKD
# Additional adjustment for:
# - Hypertension (HTN)
# - Diabetes (DM)
#
# Structure:
# 1. Baseline (T1)
# 2. Follow-up (T2)
############################################

############################
# 0. Housekeeping
############################

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

############################
# 1. Load libraries
############################

library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(minfiData)
library(limma)
library(affy)
library(dplyr)
library(bacon)

############################
# 2. Define paths (EDIT)
############################

idatspath  <- "path/to/idat_files/"
cross_file <- "path/to/cross_reactive_probes.csv"

############################
# =========================
# 3. BASELINE (T1)
# =========================
############################

pheno_T1 <- "path/to/phenotype_T1_withcells.csv"
out_T1   <- "2025.09.26_Sensitivity_EWAS_CKD_T1.txt"

targets <- read.csv(pheno_T1)

############################
# Load + preprocess methylation
############################

RGset <- read.metharray.exp(base = idatspath, targets = targets, force = TRUE)

GMset <- preprocessFunnorm(RGset, nPCs = 2, bgCorr = TRUE, dyeCorr = TRUE)

annotation <- getAnnotation(GMset)

# Remove sex chromosomes
sex_idx <- which(annotation$chr %in% c("chrX", "chrY"))
GMset_noSex <- GMset[-sex_idx, ]

# Remove cross-reactive probes
probes_remove <- read.csv(cross_file, header = FALSE)[,1]
GMset_culled <- GMset_noSex[
  !featureNames(GMset_noSex) %in% probes_remove,
]

# Extract values
mval <- getM(GMset_culled)
beta <- getBeta(GMset_culled)

# Remove infinite probes
inf_idx <- which(mval == -Inf, arr.ind = TRUE)
bad_probes <- unique(rownames(inf_idx))

mval <- mval[!rownames(mval) %in% bad_probes, ]
beta <- beta[!rownames(beta) %in% bad_probes, ]

############################
# Define variables (T1)
############################

CKD <- factor(targets$defCKD_R2)
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
  ~0 + CKD + sex + age +
    CD8T + CD4T + NK + Bcell + Mono + Neu +
    site + plate + DM + HTN
)

contrast <- makeContrasts(CKD1 - CKD0, levels = design)

############################
# EWAS M-values (T1)
############################

fit <- lmFit(mval, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

top_m <- topTable(fit, number = nrow(mval), adjust = "fdr")

############################
# Bacon (T1)
############################

tstats <- fit$t
tstats_bc <- bacon(tstats)

bc_p <- pval(tstats_bc)
bc_fdr <- p.adjust(bc_p, "fdr")

detach("package:bacon", unload = TRUE)

############################
# EWAS Betas (T1)
############################

fit_b <- lmFit(beta, design)
fit_b <- contrasts.fit(fit_b, contrast)
fit_b <- eBayes(fit_b)

top_b <- topTable(fit_b, number = nrow(beta), adjust = "fdr")

############################
# Combine (T1)
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
# 4. FOLLOW-UP (T2)
# =========================
############################

pheno_T2 <- "path/to/phenotype_T2_withcells.csv"
out_T2   <- "2025.09.26_Sensitivity_EWAS_CKD_T2.txt"

targets <- read.csv(pheno_T2)

############################
# Reload methylation data (T2)
############################

RGset <- read.metharray.exp(base = idatspath, targets = targets, force = TRUE)

GMset <- preprocessFunnorm(RGset, nPCs = 2, bgCorr = TRUE, dyeCorr = TRUE)

annotation <- getAnnotation(GMset)

sex_idx <- which(annotation$chr %in% c("chrX", "chrY"))
GMset_noSex <- GMset[-sex_idx, ]

probes_remove <- read.csv(cross_file, header = FALSE)[,1]
GMset_culled <- GMset_noSex[
  !featureNames(GMset_noSex) %in% probes_remove,
]

mval <- getM(GMset_culled)
beta <- getBeta(GMset_culled)

inf_idx <- which(mval == -Inf, arr.ind = TRUE)
bad_probes <- unique(rownames(inf_idx))

mval <- mval[!rownames(mval) %in% bad_probes, ]
beta <- beta[!rownames(beta) %in% bad_probes, ]

############################
# Variables (T2)
############################

CKD <- factor(targets$defCKD_R2)
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
# Model (T2)
############################

design <- model.matrix(
  ~0 + CKD + sex + age +
    CD8T + CD4T + NK + Bcell + Mono + Neu +
    site + plate + DM + HTN
)

contrast <- makeContrasts(CKD1 - CKD0, levels = design)

############################
# EWAS (T2)
############################

fit <- lmFit(mval, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

top_m <- topTable(fit, number = nrow(mval), adjust = "fdr")

############################
# Bacon (T2)
############################

tstats <- fit$t
tstats_bc <- bacon(tstats)

bc_p <- pval(tstats_bc)
bc_fdr <- p.adjust(bc_p, "fdr")

detach("package:bacon", unload = TRUE)

############################
# Betas (T2)
############################

fit_b <- lmFit(beta, design)
fit_b <- contrasts.fit(fit_b, contrast)
fit_b <- eBayes(fit_b)

top_b <- topTable(fit_b, number = nrow(beta), adjust = "fdr")

############################
# Combine (T2)
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

############################
# END OF SCRIPT
############################