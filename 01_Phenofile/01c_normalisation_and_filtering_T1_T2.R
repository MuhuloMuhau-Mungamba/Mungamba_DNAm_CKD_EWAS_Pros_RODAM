############################################################
# Pros-RODAM CKD EWAS
# Step 1c: Normalisation and probe filtering (T1 & T2)
#
# Purpose:
# - Functional normalisation of DNAm data (funnorm)
# - Remove sex chromosome probes
# - Remove SNPs and cross-hybridisation probes
# - Generate beta values and M-values
# - Perform QC checks
#
# Input:
# - CKD_EWAS_T1_phenofile_wCells.csv
# - CKD_EWAS_T2_phenofile_wCells.csv
# - IDAT files
#
# Output:
# - Normalised beta values (T1 & T2)
# - Normalised M-values (T1 & T2)
# - Annotation files
# - QC plots
#
# Final sample:
# - n = 274 (T1), n = 274 (T2)
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
library(minfiData)


## ---------------------------------------------------------
## 2. Define paths
## ---------------------------------------------------------

idatspath <- "path/to/IDAT_files/"
cross_file <- "path/to/cross_reactive_probes.csv"


## ---------------------------------------------------------
## 3. Function: Normalisation pipeline
## ---------------------------------------------------------

run_normalisation <- function(pheno_file, output_prefix) {
  
  ## Load phenotype
  targets <- read.csv(pheno_file)
  
  ## Read IDATs
  RGset <- read.metharray.exp(
    base = idatspath,
    targets = targets,
    force = TRUE
  )
  
  ## Raw beta
  beta_raw <- getBeta(RGset)
  
  ## Functional normalisation
  GMset <- preprocessFunnorm(
    RGset,
    nPCs = 2,
    bgCorr = TRUE,
    dyeCorr = TRUE,
    verbose = TRUE
  )
  
  ## Annotation
  annotation <- getAnnotation(GMset)
  
  ## Remove sex chromosomes
  sex_idx <- which(annotation$chr %in% c("chrX", "chrY"))
  GMset_noSex <- GMset[-sex_idx, ]
  
  ## Remove cross-hybridising probes
  probes4removal <- read.csv(
    cross_file,
    stringsAsFactors = FALSE,
    col.names = 1
  )[, 1]
  
  GMset_culled <- GMset_noSex[
    !featureNames(GMset_noSex) %in% probes4removal,
  ]
  
  ## Extract M-values and beta values
  mval <- getM(GMset_culled)
  beta <- getBeta(GMset_culled)
  
  ## Remove infinite M-values
  inf_idx <- which(mval == -Inf, arr.ind = TRUE)
  bad_probes <- unique(rownames(inf_idx))
  
  mval_clean <- mval[!rownames(mval) %in% bad_probes, ]
  beta_clean <- beta[!rownames(beta) %in% bad_probes, ]
  
  ## Save outputs
  write.table(
    beta_raw,
    paste0(output_prefix, "_betaraw.txt"),
    sep = "\t",
    quote = FALSE
  )
  
  write.table(
    mval_clean,
    paste0(output_prefix, "_Mval.txt"),
    sep = "\t",
    quote = FALSE
  )
  
  write.table(
    beta_clean,
    paste0(output_prefix, "_Beta.txt"),
    sep = "\t",
    quote = FALSE
  )
  
  ## QC plots
  pdf(paste0(output_prefix, "_QC_density_raw.pdf"))
  densityPlot(RGset, main = "Raw beta", xlab = "Beta")
  dev.off()
  
  pdf(paste0(output_prefix, "_QC_density_normalised.pdf"))
  densityPlot(beta_clean, main = "Normalised beta", xlab = "Beta")
  dev.off()
  
  return("Done")
}


## =========================================================
## PART A: BASELINE (T1)
## =========================================================

run_normalisation(
  pheno_file = "CKD_EWAS_T1_phenofile_wCells.csv",
  output_prefix = "CKD_EWAS_T1"
)


## =========================================================
## PART B: FOLLOW-UP (T2)
## =========================================================

run_normalisation(
  pheno_file = "CKD_EWAS_T2_phenofile_wCells.csv",
  output_prefix = "CKD_EWAS_T2"
)


## ---------------------------------------------------------
## End of script
## ---------------------------------------------------------