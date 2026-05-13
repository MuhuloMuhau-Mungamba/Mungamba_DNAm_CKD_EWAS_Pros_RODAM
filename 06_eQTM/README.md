# CKD eQTM analysis scripts

This repository contains GitHub-ready R Markdown scripts for the CKD DNA methylation / gene expression eQTM workflow. Local absolute paths were removed and replaced with project-relative folders.

## Folder structure

Create the following folders before running the workflow:

```text
.
├── input/
├── output/
│   └── plots/
├── functions/
└── scripts/
```

The scripts assume they are run from the repository root. Input data and annotation files should be placed in `input/`, custom functions in `functions/`, and generated files will be written to `output/` and `output/plots/`.

## Scripts

| Script | Purpose |
|---|---|
| `ckd_preprocessing_GEX_liftover_hg38tohg19.Rmd` | Converts/lifts gene expression coordinates from hg38 to hg19 and prepares the filtered normalized gene expression matrix. |
| `CKD_preprocessing_DMRS_geneoverlaps.Rmd` | Extends DMR regions and maps overlapping hg19 genes using a GTF annotation file. |
| `PHSugg_eQTM_CKD.Rmd` | Runs the main binary CKD eQTM analysis. |
| `PHSugg_intraclass_eQTM_CKD.Rmd` | Runs intraclass eQTM analyses for controls and samples. |
| `PHSugg_eQTMoutput_Annotation.Rmd` | Annotates eQTM output with gene/region metadata. |
| `CorrelationPlot_eQTM.Rmd` | Generates selected eQTM correlation plots. |
| `Violinplot_DMR_betaMatrix.Rmd` | Generates DMR methylation beta-value violin plots. |

## Suggested run order

1. `ckd_preprocessing_GEX_liftover_hg38tohg19.Rmd`
2. `CKD_preprocessing_DMRS_geneoverlaps.Rmd`
3. `PHSugg_eQTM_CKD.Rmd`
4. `PHSugg_intraclass_eQTM_CKD.Rmd` if intraclass analyses are required
5. `PHSugg_eQTMoutput_Annotation.Rmd`
6. `CorrelationPlot_eQTM.Rmd`
7. `Violinplot_DMR_betaMatrix.Rmd`

## Required inputs

Place required input files in `input/`. Based on the scripts, expected inputs include phenotype metadata, normalized gene expression matrices, methylation beta matrices, filtered annotation files, DMR tables, and the hg19 GTF annotation file. File names are referenced directly in each script through `file.path(input_dir, "...")`.

Large or sensitive research data files should not be uploaded to GitHub. Keep those files outside the public repository and provide instructions for authorized users to obtain them.

## Notes for GitHub release

- Local absolute paths and active `setwd()` calls were removed.
- Scripts now use `input_dir`, `output_dir`, `plots_dir`, and `functions_dir` variables.
- Sample-level data and private file-system details should remain excluded from the repository.
- Add data files to `.gitignore` if they contain sensitive or controlled-access information.

## Minimal `.gitignore`

```gitignore
input/
output/
*.RData
*.rds
*.RDS
*.csv
*.txt
*.tsv
*.xlsx
.DS_Store
.Rhistory
.RData
```
