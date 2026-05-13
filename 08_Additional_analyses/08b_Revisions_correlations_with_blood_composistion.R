##Addressing reviewers comments for CKD EWAS##
#Muhau Mungamba: 2026.04.15
##Analysis Revisions for eBioMedicine##

#Load library pathway#

.libPaths("/appdata/users/P083854/R_packages/")

#set working directory to the L-Drive#

setwd("/net/cifs/prod1.umcinfra.nl/data/Basic/")## original working directory to L drive, run it first


##Reviewer 1 comment with Analysis##
#2. The authors have already adjusted for estimated blood cell fractions 
#using the Houseman method. It would be very informative if the authors 
#could examine correlations between methylation levels at the major hits 
#and estimated proportions of neutrophils, monocytes, CD4-positive T cells, 
#CD8-positive T cells, B cells, and NK cells, or assess interaction effects 
#between outcome and blood cell composition. Such supplementary analyses 
#would help readers more concretely understand whether the observed signals 
#may primarily arise from specific blood cell populations.

save.image()
savehistory()


#setworking directory for the revisions#

setwd("/net/cifs/prod1.umcinfra.nl/data/Basic/divjk/sg/ProsRODAM/MUHAU_prosRODAM/CKD_DNA_Methylation/Manuscript/Publication/eBioMedicine/Revisions/Analysis/")

##clear environment##

rm(list = ls())

library(minfi)
library(Hmisc)
library(dplyr)

#Load pheno with cells#

phenofile = "/net/cifs/prod1.umcinfra.nl/data/Basic/divjk/sg/ProsRODAM/MUHAU_prosRODAM/CKD_DNA_Methylation/Results/Timepoint_1/Repeat_Analysis/Pheno/ 2025.08.21_Repeat_Timepoint1pheno_withcells.csv"

targets <- read.csv(phenofile)

##Load methylation data##

##Path to idat data
idatspath= "/net/cifs/prod1.umcinfra.nl/data/Basic/divjk/sg/ProsRODAM/EpigeneticData_IlluminaEpic/RawData/idats/"


BasePath=idatspath

## this file has the probe names of probes that crosshybidize, the ones already filtered for MAF 0.05
## they need to be deleted from the data

cross="/net/cifs/prod1.umcinfra.nl/data/Basic/divjk/sg/ProsRODAM/Filtering/AFR_EUR_0.05_CrosHYB_SNPSlist.csv"

### load raw data from idats
RGset=read.metharray.exp(BasePath,targets)

### get unnormalized beta values
betaraw=getBeta(RGset)

#### Funnorm normalization of the data
GMset=preprocessFunnorm(RGset, nPCs=2, bgCorr = TRUE,
                        dyeCorr = TRUE, verbose = TRUE)
pd <- pData(RGset)

## get annotation to filter probes
annotation <- getAnnotation(GMset)

## make a matrix from the annotation
probe.features=as.matrix(annotation)

###Get X and yProbes
indices <- which(annotation$chr == "chrX"| annotation$chr == "chrY")

##filter out X and Y probes
GMset.nosex <- GMset[-indices,]

## filterout probes with snps (This step is not necessary here)
#GMset.culled1 <-dropLociWithSnps( GMset.nosex, snps = c("CpG", "SBE"), maf = 0.01, snpAnno = NULL)

## Read crosshybridization probe file
probes4removal <- read.csv(cross, stringsAsFactors  = F, col.names = 1)[,1]

## filter uit cross hybridization probes
GMset.culled <- GMset.nosex [which(!featureNames(GMset.nosex) %in% probes4removal), ]

## some M values give invinite values that give errors
## these probes will be removed from both the betas and the M values
## so both list are the same length

mval = getM(GMset.culled)
beta2=getBeta(GMset.culled)
x=which(mval == -Inf, arr.ind = T)
xx=sort(rownames(x))
mval2=mval[!rownames(mval) %in% xx, ]
beta=beta2[!rownames(beta2) %in% xx, ]
remove("x","xx")

##extract top CpGs for Timepoint 1##

top_cpgs <- c(
  "cg26122413","cg07333667","cg23714355","cg06034807",
  "cg05215481","cg04646451","cg15646763","cg05160409"
)

# Keep only CpGs that exist
top_cpgs <- top_cpgs[top_cpgs %in% rownames(mval)]

cpg_matrix <- t(mval[top_cpgs, ])   # samples × CpGs

##merge with pheno

data_full <- cbind(targets, cpg_matrix)

##run correlation##

cell_types <- c("CD8T","CD4T","NK","Bcell","Mono","Neu")

results_list <- list()

for(cpg in top_cpgs){
  for(cell in cell_types){
    
    test <- cor.test(data_full[[cpg]], data_full[[cell]], method="pearson")
    
    results_list[[paste(cpg, cell, sep="_")]] <- data.frame(
      CpG = cpg,
      Cell = cell,
      Correlation = test$estimate,
      Pvalue = test$p.value
    )
  }
}

results <- do.call(rbind, results_list)

# Adjust p-values
results$FDR <- p.adjust(results$Pvalue, method="fdr")

write.csv(results,
          "Reviewer1_CellCorrelation_Results_Timepoint1.csv",
          row.names=FALSE)


##Timepoint 2 (Follow-up correlations)

##clear environment##

rm(list = ls())

library(minfi)
library(Hmisc)
library(dplyr)

#Load pheno with cells#

phenofile = "/net/cifs/prod1.umcinfra.nl/data/Basic/divjk/sg/ProsRODAM/MUHAU_prosRODAM/CKD_DNA_Methylation/Results/Timepoint_2/Repeat_Analysis/Pheno/2025.08.23_Matchedtimepoint2_pheno_withcells.csv"

targets <- read.csv(phenofile)

##Load methylation data##

##Path to idat data
idatspath= "/net/cifs/prod1.umcinfra.nl/data/Basic/divjk/sg/ProsRODAM/EpigeneticData_IlluminaEpic/RawData/idats/"


BasePath=idatspath

## this file has the probe names of probes that crosshybidize, the ones already filtered for MAF 0.05
## they need to be deleted from the data

cross="/net/cifs/prod1.umcinfra.nl/data/Basic/divjk/sg/ProsRODAM/Filtering/AFR_EUR_0.05_CrosHYB_SNPSlist.csv"

### load raw data from idats
RGset=read.metharray.exp(BasePath,targets)

### get unnormalized beta values
betaraw=getBeta(RGset)

#### Funnorm normalization of the data
GMset=preprocessFunnorm(RGset, nPCs=2, bgCorr = TRUE,
                        dyeCorr = TRUE, verbose = TRUE)
pd <- pData(RGset)

## get annotation to filter probes
annotation <- getAnnotation(GMset)

## make a matrix from the annotation
probe.features=as.matrix(annotation)

###Get X and yProbes
indices <- which(annotation$chr == "chrX"| annotation$chr == "chrY")

##filter out X and Y probes
GMset.nosex <- GMset[-indices,]

## filterout probes with snps (This step is not necessary here)
#GMset.culled1 <-dropLociWithSnps( GMset.nosex, snps = c("CpG", "SBE"), maf = 0.01, snpAnno = NULL)

## Read crosshybridization probe file
probes4removal <- read.csv(cross, stringsAsFactors  = F, col.names = 1)[,1]

## filter uit cross hybridization probes
GMset.culled <- GMset.nosex [which(!featureNames(GMset.nosex) %in% probes4removal), ]

## some M values give invinite values that give errors
## these probes will be removed from both the betas and the M values
## so both list are the same length

mval = getM(GMset.culled)
beta2=getBeta(GMset.culled)
x=which(mval == -Inf, arr.ind = T)
xx=sort(rownames(x))
mval2=mval[!rownames(mval) %in% xx, ]
beta=beta2[!rownames(beta2) %in% xx, ]
remove("x","xx")

##extract top CpGs for Timepoint 2##

top_cpgs <- c(
  "cg26122413",
  "cg23714355",
  "cg06034807",
  "cg04646451",
  "cg15646763",
  "cg05215481"
)

# Keep only CpGs that exist
top_cpgs <- top_cpgs[top_cpgs %in% rownames(mval)]

cpg_matrix <- t(mval[top_cpgs, ])   # samples × CpGs

##merge with pheno

data_full <- cbind(targets, cpg_matrix)

##run correlation##

cell_types <- c("CD8T","CD4T","NK","Bcell","Mono","Neu")

results_list <- list()

for(cpg in top_cpgs){
  for(cell in cell_types){
    
    test <- cor.test(data_full[[cpg]], data_full[[cell]], method="pearson")
    
    results_list[[paste(cpg, cell, sep="_")]] <- data.frame(
      CpG = cpg,
      Cell = cell,
      Correlation = test$estimate,
      Pvalue = test$p.value
    )
  }
}

results <- do.call(rbind, results_list)

# Adjust p-values
results$FDR <- p.adjust(results$Pvalue, method="fdr")

write.csv(results,
          "Reviewer1_CellCorrelation_Results_Timepoint2.csv",
          row.names=FALSE)


