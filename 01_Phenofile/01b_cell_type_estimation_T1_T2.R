############################################################
# Pros-RODAM CKD EWAS
# Step 1b: Blood cell-type proportion estimation (T1 & T2)
#
# Purpose:
# - Estimate blood cell-type proportions using Houseman method
# - Apply to both baseline (T1) and follow-up (T2)
# - Merge cell counts with phenotype files
#
# Input:
# - CKD_EWAS_T1_phenofile.csv
# - CKD_EWAS_T2_phenofile.csv
#
# Output:
# - CKD_EWAS_T1_phenofile_wCells.csv
# - CKD_EWAS_T2_phenofile_wCells.csv
#
# Final sample:
# - 274 individuals per timepoint
#
# Author: Muhulo M. Mungamba & Andrea Venema
############################################################


## ---------------------------------------------------------
## 0. Housekeeping
## ---------------------------------------------------------

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)


## ---------------------------------------------------------
## 1. Load packages
## ---------------------------------------------------------

library(minfi)
library(IlluminaHumanMethylationEPICmanifest)
library(FlowSorted.Blood.EPIC)


## ---------------------------------------------------------
## 2. Define IDAT path
## ---------------------------------------------------------

idatspath <- "path/to/IDAT_files/"
BasePath <- idatspath


## =========================================================
## PART A: BASELINE (T1)
## =========================================================

## ---------------------------------------------------------
## 3A. Load phenotype (T1)
## ---------------------------------------------------------

targets_T1 <- read.csv("CKD_EWAS_T1_phenofile.csv")

## ---------------------------------------------------------
## 4A. Read IDAT files
## ---------------------------------------------------------

RGset_T1 <- read.metharray.exp(
  base = BasePath,
  targets = targets_T1
)

## ---------------------------------------------------------
## 5A. Estimate cell counts
## ---------------------------------------------------------

cell_est_T1 <- estimateCellCounts(
  RGset_T1,
  compositeCellType = "Blood",
  cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),
  returnAll = TRUE,
  meanPlot = FALSE,
  verbose = FALSE
)

## ---------------------------------------------------------
## 6A. Prepare cell count table
## ---------------------------------------------------------

cell_df_T1 <- as.data.frame(cell_est_T1$counts)
cell_df_T1$Basename <- rownames(cell_df_T1)

cell_df_T1 <- cell_df_T1[, c("Basename", setdiff(names(cell_df_T1), "Basename"))]

## ---------------------------------------------------------
## 7A. Merge with phenotype
## ---------------------------------------------------------

pheno_T1_wCells <- merge(
  targets_T1,
  cell_df_T1,
  by = "Basename",
  all.x = TRUE
)

## ---------------------------------------------------------
## 8A. Save T1 file
## ---------------------------------------------------------

write.csv(
  pheno_T1_wCells,
  "CKD_EWAS_T1_phenofile_wCells.csv",
  row.names = FALSE
)



## =========================================================
## PART B: FOLLOW-UP (T2)
## =========================================================

## ---------------------------------------------------------
## 3B. Load phenotype (T2)
## ---------------------------------------------------------

targets_T2 <- read.csv("CKD_EWAS_T2_phenofile.csv")

## ---------------------------------------------------------
## 4B. Read IDAT files
## ---------------------------------------------------------

RGset_T2 <- read.metharray.exp(
  base = BasePath,
  targets = targets_T2
)

## ---------------------------------------------------------
## 5B. Estimate cell counts
## ---------------------------------------------------------

cell_est_T2 <- estimateCellCounts(
  RGset_T2,
  compositeCellType = "Blood",
  cellTypes = c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran"),
  returnAll = TRUE,
  meanPlot = FALSE,
  verbose = FALSE
)

## ---------------------------------------------------------
## 6B. Prepare cell count table
## ---------------------------------------------------------

cell_df_T2 <- as.data.frame(cell_est_T2$counts)
cell_df_T2$Basename <- rownames(cell_df_T2)

cell_df_T2 <- cell_df_T2[, c("Basename", setdiff(names(cell_df_T2), "Basename"))]

## ---------------------------------------------------------
## 7B. Merge with phenotype
## ---------------------------------------------------------

pheno_T2_wCells <- merge(
  targets_T2,
  cell_df_T2,
  by = "Basename",
  all.x = TRUE
)

## ---------------------------------------------------------
## 8B. Save T2 file
## ---------------------------------------------------------

write.csv(
  pheno_T2_wCells,
  "CKD_EWAS_T2_phenofile_wCells.csv",
  row.names = FALSE
)


## ---------------------------------------------------------
## End of script
## ---------------------------------------------------------