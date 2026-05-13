# Author: Muhulo Muhau Mungamba
# Date: 2026-04-16
# Purpose: Extract and validate literature CpGs in baseline and follow-up EWAS

############################################################
# LOAD LIBRARIES
############################################################

library(dplyr)
library(stringr)

############################################################
# SANITY CHECK
############################################################

exists("combined_base")
exists("combined_follow")

############################################################
# 1. DEFINE LITERATURE CpGs
############################################################

cpg_list <- c(
  "cg02304370","cg04460609","cg00501876","cg04864179",
  "cg23597162","cg26099045","cg12644285",
  "cg17944885","cg06158227","cg20777437",
  "cg18181703","cg02711608","cg23570810",
  "cg18194850","cg07242931","cg20146909",
  "cg16618493","cg03297731",
  "cg01068906","cg10632966"
)

############################################################
# 2. CpG ANNOTATION
############################################################

cpg_annotation <- data.frame(
  CpG = cpg_list,
  Gene = c(
    "PHRF1","LDB2","CSRNP1","IRF5",
    "JAZF1","PELI1","CHD2",
    "ZNF788","ZSCAN29","CDCP2",
    "SOCS3","SLC1A5","IFITM1",
    "SUCLG2","MAN1C1","LRRC8D",
    "ZBTB7B","GDPD3",
    "NOD2","RPEL1"
  ),
  Category = c(
    "Causal","Causal","Causal","Causal",
    "Kidney_tissue","Kidney_tissue","Kidney_tissue",
    "Top_hit","Top_hit","Top_hit",
    "UACR","UACR","UACR",
    "Replication","Replication","Replication",
    "Replication","Replication",
    "Immune","Immune"
  ),
  stringsAsFactors = FALSE
)

############################################################
# 3. EXTRACTION FUNCTION (FIXED FOR YOUR DATA)
############################################################

extract_cpgs <- function(df, label){
  
  out <- df %>%
    filter(CpG %in% cpg_list) %>%
    select(
      CpG,
      logFCMval,
      BaconPvalue_Mvalues,
      P.Value,
      adj.P.Val
    ) %>%
    mutate(
      Analysis = label,
      Direction = ifelse(logFCMval > 0, "Positive", "Negative"),
      Significant = ifelse(BaconPvalue_Mvalues < 0.05, "Yes", "No")
    ) %>%
    rename(
      Effect = logFCMval,
      Pvalue = BaconPvalue_Mvalues
    )
  
  return(out)
}

############################################################
# 4. RUN FOR BASELINE + FOLLOW-UP
############################################################

cpg_base   <- extract_cpgs(combined_base, "Baseline")
cpg_follow <- extract_cpgs(combined_follow, "Follow-up")

############################################################
# 5. COMBINE RESULTS
############################################################

cpg_all <- bind_rows(cpg_base, cpg_follow)

############################################################
# 6. MERGE WITH ANNOTATION
############################################################

cpg_final <- cpg_all %>%
  left_join(cpg_annotation, by = "CpG") %>%
  arrange(Analysis, Pvalue)

############################################################
# 7. QUICK VIEW
############################################################

cat("\n--- BASELINE CpGs ---\n")
print(cpg_final %>% filter(Analysis == "Baseline"))

cat("\n--- FOLLOW-UP CpGs ---\n")
print(cpg_final %>% filter(Analysis == "Follow-up"))

############################################################
# 8. SUMMARY (FOR RESULTS TEXT)
############################################################

summary_table <- cpg_final %>%
  group_by(Analysis) %>%
  summarise(
    Total_CpGs = n(),
    Significant_CpGs = sum(Significant == "Yes")
  )

print(summary_table)

############################################################
# 9. SAVE RESULTS
############################################################

outdir <- "Results/CKD_CpG_validation"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.table(
  cpg_final,
  file.path(outdir, "Literature_CpG_overlap.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
# END
############################################################
