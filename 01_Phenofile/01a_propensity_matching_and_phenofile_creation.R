############################################################
# Pros-RODAM CKD EWAS
# Step 1a: Propensity score matching and phenotype creation
#
# Purpose:
# - Select incident CKD cases with DNAm data (n = 92)
# - Select controls with DNAm data (n ≈ 736)
# - Perform 1:2 matching (controls to cases) based on:
#     age, sex, and study site
# - Create matched cohort (92 cases, 182 controls)
# - Generate phenotype files for:
#     (1) Baseline (T1)
#     (2) Follow-up (T2)
#
# Final dataset:
# - 274 individuals (92 cases, 182 controls)
# - Baseline (T1): n = 274
# - Follow-up (T2): n = 274
# - Total observations: n = 548
#
# Author: Muhulo M. Mungamba
############################################################


## ---------------------------------------------------------
## 0. Housekeeping
## ---------------------------------------------------------

rm(list = ls())
gc()


## ---------------------------------------------------------
## 1. Load packages
## ---------------------------------------------------------

library(dplyr)


## ---------------------------------------------------------
## 2. Load full phenotype dataset
## ---------------------------------------------------------

pheno <- read.csv("path/to/full_pros_rodam_dataset.csv")


## ---------------------------------------------------------
## 3. Restrict to DNAm subset
## ---------------------------------------------------------

pheno_dnAm <- pheno %>%
  filter(!is.na(DNAm_ID))


## ---------------------------------------------------------
## 4. Define cases and controls (based on follow-up CKD)
## ---------------------------------------------------------

cases <- pheno_dnAm %>%
  filter(defCKD_R2 == 1)

controls <- pheno_dnAm %>%
  filter(defCKD_R2 == 0)

## Expected:
nrow(cases)     # 92
nrow(controls)  # ~736


## ---------------------------------------------------------
## 5. Estimate propensity scores (cases + controls)
## ---------------------------------------------------------

ps_model <- glm(defCKD_R2 ~ R1_Age + R1_Sex + R1_Site,
                data = rbind(cases, controls),
                family = binomial)

cases$ps <- predict(ps_model, newdata = cases, type = "response")
controls$ps <- predict(ps_model, newdata = controls, type = "response")


## ---------------------------------------------------------
## 6. Match controls to cases (1:2, without replacement)
## ---------------------------------------------------------

match_1to2 <- function(cases, controls) {
  
  matched_controls <- data.frame()
  remaining_controls <- controls
  
  for (i in 1:nrow(cases)) {
    
    case_ps <- cases$ps[i]
    
    # Distance to all controls
    dist <- abs(remaining_controls$ps - case_ps)
    
    # Select 2 closest controls
    idx <- order(dist)[1:2]
    
    matched_controls <- rbind(
      matched_controls,
      remaining_controls[idx, ]
    )
    
    # Remove selected controls
    remaining_controls <- remaining_controls[-idx, ]
  }
  
  matched <- rbind(cases, matched_controls)
  return(matched)
}

matched_data <- match_1to2(cases, controls)


## ---------------------------------------------------------
## 7. Check final matching
## ---------------------------------------------------------

table(matched_data$defCKD_R2)
# Expected:
# 0 = 182 controls
# 1 = 92 cases


## ---------------------------------------------------------
## 8. Create baseline (T1) phenotype file
## ---------------------------------------------------------

pheno_T1 <- matched_data %>%
  filter(DNA_BaselineFollowUp == 1)

nrow(pheno_T1)  # 274


## ---------------------------------------------------------
## 9. Create follow-up (T2) phenotype file
## ---------------------------------------------------------

# Use same individuals via RodamID
pheno_T2 <- pheno %>%
  filter(DNA_BaselineFollowUp == 2) %>%
  filter(RodamID %in% pheno_T1$RodamID)

nrow(pheno_T2)  # 274


## ---------------------------------------------------------
## 10. Consistency checks
## ---------------------------------------------------------

length(unique(pheno_T1$RodamID))  # 274
length(unique(pheno_T2$RodamID))  # 274


## ---------------------------------------------------------
## 11. Save phenotype files
## ---------------------------------------------------------

write.csv(pheno_T1,
          "CKD_EWAS_T1_phenofile.csv",
          row.names = FALSE)

write.csv(pheno_T2,
          "CKD_EWAS_T2_phenofile.csv",
          row.names = FALSE)


## ---------------------------------------------------------
## End of script
## ---------------------------------------------------------
