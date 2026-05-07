############################################
# Pros-RODAM CKD EWAS
# 03b_DMR_delta_eGFR_T1_T2.R
#
# Author: Muhulo Muhau Mungamba & Andrea Venema
#
# DMR analysis using bumphunter
# Outcome: Delta eGFR
# Timepoints: Baseline (T1) and Follow-up (T2)
############################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

############################
# Load required packages
############################

library(IlluminaHumanMethylationEPICmanifest)
library(minfi)
library(minfiData)
library(limma)
library(affy)
library(GenomicRanges)

############################
# Define paths (EDIT HERE)
############################

idatspath <- "path/to/idat_files/"
cross_file <- "path/to/crosshyb_probe_list.csv"

############################
# ==========================
# T1: BASELINE
# ==========================
############################

pheno_file_t1 <- "path/to/T1_phenofile_with_cells.csv"

targets <- read.table(pheno_file_t1, sep = ",", header = TRUE)

RGset <- read.metharray.exp(base = idatspath, targets = targets)

GMset <- preprocessFunnorm(
  RGset,
  nPCs = 2,
  bgCorr = TRUE,
  dyeCorr = TRUE,
  verbose = TRUE
)

annotation <- getAnnotation(GMset)

############################
# Probe filtering
############################

# Remove sex chromosomes
sex_idx <- which(annotation$chr %in% c("chrX", "chrY"))
GMset_noSex <- GMset[-sex_idx, ]

# Remove cross-hybridisation probes
probes4removal <- read.csv(
  cross_file,
  stringsAsFactors = FALSE,
  col.names = 1
)[, 1]

GMset_culled <- GMset_noSex[
  !featureNames(GMset_noSex) %in% probes4removal,
]

############################
# Design matrix (T1)
############################

deltaegfr <- targets$delta_eGFR

sex  <- factor(targets$R1_Sex)
site <- factor(targets$R1_Site)
age  <- targets$R1_Age

CD8T <- targets$CD8T
CD4T <- targets$CD4T
NK   <- targets$NK
Bcell<- targets$Bcell
Mono <- targets$Mono
Neu  <- targets$Neu

plate <- targets$plate

design <- model.matrix(
  ~ deltaegfr + sex + age +
    CD8T + CD4T + NK + Bcell + Mono + Neu +
    site + plate
)

############################
# Annotation
############################

annot <- getAnnotation(RGset)

annot_gr <- makeGRangesFromDataFrame(
  annot,
  keep.extra.columns = TRUE,
  start.field = "pos",
  end.field = "pos"
)

############################
# Bumphunter (cutoff 0.0006)
############################

tab <- bumphunter(GMset_culled, design, coef = 2, cutoff = 0.0006)
dmr <- subset(tab$table, L > 2)

tab <- bumphunter(
  GMset_culled,
  design,
  coef = 2,
  cutoff = 0.0006,
  B = 500,
  nullMethod = "bootstrap"
)

dmr <- subset(tab$table, L > 2)

dmr_gr <- makeGRangesFromDataFrame(dmr, keep.extra.columns = TRUE)
dmr_annot <- annot_gr[nearest(dmr_gr, annot_gr)]

dmr_out <- cbind(
  dmr,
  Gene_start = start(dmr_annot),
  Gene_end   = end(dmr_annot),
  Gene_Name  = dmr_annot$UCSC_RefGene_Name
)

write.table(
  dmr_out,
  "DMRS_deltaeGFR_T1_cutoff0.0006_boot500.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################
# ==========================
# T2: FOLLOW-UP
# ==========================
############################

pheno_file_t2 <- "path/to/T2_phenofile_with_cells.csv"

targets <- read.table(pheno_file_t2, sep = ",", header = TRUE)

RGset <- read.metharray.exp(base = idatspath, targets = targets)

GMset <- preprocessFunnorm(
  RGset,
  nPCs = 2,
  bgCorr = TRUE,
  dyeCorr = TRUE,
  verbose = TRUE
)

annotation <- getAnnotation(GMset)

############################
# Probe filtering
############################

sex_idx <- which(annotation$chr %in% c("chrX", "chrY"))
GMset_noSex <- GMset[-sex_idx, ]

probes4removal <- read.csv(
  cross_file,
  stringsAsFactors = FALSE,
  col.names = 1
)[, 1]

GMset_culled <- GMset_noSex[
  !featureNames(GMset_noSex) %in% probes4removal,
]

############################
# Design matrix (T2)
############################

deltaegfr <- targets$delta_eGFR

sex  <- factor(targets$R2_Sex)
site <- factor(targets$R2_Site)
age  <- targets$R2_Age

CD8T <- targets$CD8T
CD4T <- targets$CD4T
NK   <- targets$NK
Bcell<- targets$Bcell
Mono <- targets$Mono
Neu  <- targets$Neu

plate <- targets$plate

design <- model.matrix(
  ~ deltaegfr + sex + age +
    CD8T + CD4T + NK + Bcell + Mono + Neu +
    site + plate
)

############################
# Bumphunter (T2)
############################

tab <- bumphunter(GMset_culled, design, coef = 2, cutoff = 0.0006)
dmr <- subset(tab$table, L > 2)

tab <- bumphunter(
  GMset_culled,
  design,
  coef = 2,
  cutoff = 0.0006,
  B = 500,
  nullMethod = "bootstrap"
)

dmr <- subset(tab$table, L > 2)

dmr_gr <- makeGRangesFromDataFrame(dmr, keep.extra.columns = TRUE)
dmr_annot <- annot_gr[nearest(dmr_gr, annot_gr)]

dmr_out <- cbind(
  dmr,
  Gene_start = start(dmr_annot),
  Gene_end   = end(dmr_annot),
  Gene_Name  = dmr_annot$UCSC_RefGene_Name
)

write.table(
  dmr_out,
  "DMRS_deltaeGFR_T2_cutoff0.0006_boot500.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################
# END SCRIPT
############################