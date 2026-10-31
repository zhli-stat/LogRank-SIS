# Log-rank SIS: Feature Screening for Ultrahigh-Dimensional Interval-Censored Data

This repository contains the source code used to reproduce the numerical results in

> **An Efficient Feature Screening Procedure for Ultrahigh-Dimensional Interval-Censored Response**
> Zhuohang Li, Xuewei Cheng (2026).

The proposed **Log-rank SIS** procedure ranks binary predictors by the absolute
value of a variance-adjusted generalized log-rank score. It requires only a single
pooled NPMLE fit of the response distribution. In the ADNI application
(*n* = 278, *p* = 530,437 SNPs) genome-wide screening takes about 4 seconds,
compared with about 945 seconds for ADD-SIS.

## Repository structure

```
LogRank-SIS/
├── README.md
├── Plots/
│   ├── plot_screening_performance.R
│   ├── plot_gof_sensitivity.py
│   ├── plot_salary_cdf.R
│   └── plot_adni_manhattan.R
├── Simulation/
│   ├── Example1-3/Code/
│   │   ├── LM-Log-rank.R
│   │   ├── GLM-Log-rank.R
│   │   ├── COX-Log-rank.R
│   │   └── Proportion-Tables.R
│   └── Outlier_Case/Code/
│       └── Outlier_Case.R
└── CaseStudy/
    ├── Code/
    │   ├── ADNI_Data_Preprocessing.R
    │   ├── CaseStudyCode_Log-rank_ADNI.Rmd
    │   ├── CaseStudyCode_Log-rank_Salary.Rmd
    │   ├── get_top50_genomewide.R
    │   └── snp_annotation/
    │       ├── config.py
    │       ├── fetch_dbsnp.py
    │       ├── fetch_gwas_catalog.py
    │       ├── compare_and_summarize.py
    │       └── enrich_and_nearest.py
    └── Data/                 # not shipped; see "Data preparation"
```

`plot_adni_manhattan.R` is in `Plots/` and runs from the repository root. It
temporarily switches the working directory to `CaseStudy/Code/` to source
`get_top50_genomewide.R`, so the ADNI screening script's own relative paths
are unchanged.

## Code and paper correspondence

| Script | Paper output |
|---|---|
| `Simulation/Example1-3/Code/LM-Log-rank.R` | Example 1 (linear model), Table 1 |
| `Simulation/Example1-3/Code/GLM-Log-rank.R` | Example 2 (Poisson regression), Table 2 |
| `Simulation/Example1-3/Code/COX-Log-rank.R` | Example 3 (Cox model), Table 3 |
| `Simulation/Example1-3/Code/Proportion-Tables.R` | Tables 4-6 (selection proportions) |
| `Simulation/Outlier_Case/Code/Outlier_Case.R` | Table 7 (interval-expansion and upper-bound perturbations) |
| `Plots/plot_screening_performance.R` | Figure 1 |
| `Plots/plot_gof_sensitivity.py` | Figure 2 |
| `Plots/plot_adni_manhattan.R` | Figure 3 (screening utilities by chromosomal position) |
| `CaseStudy/Code/CaseStudyCode_Log-rank_ADNI.Rmd` | ADNI analysis (Section 4), Tables 8-11, inputs for Figures 2-3 |
| `CaseStudy/Code/CaseStudyCode_Log-rank_Salary.Rmd` | Salary analysis (Appendix B) |
| `Plots/plot_salary_cdf.R` | Figure B.1 |
| `CaseStudy/Code/get_top50_genomewide.R` | Genome-wide Log-rank SIS top-50 SNPs (input to Table C.1) |
| `CaseStudy/Code/snp_annotation/*.py` | Appendix C: annotation of the 50 selected SNPs (Table C.1) |
| `CaseStudy/Code/ADNI_Data_Preprocessing.R` | Prepares `ADNI_Y_Final.rds` and `ADNI_X_Final.rds` |

Table 11's ADD-SIS screening time and goodness-of-fit are computed inside
`CaseStudyCode_Log-rank_ADNI.Rmd` (the `ADD-SIS [Phase II]` chunk), so no
separate ADD-SIS script is shipped.

### The SNP annotation pipeline

The five scripts in `CaseStudy/Code/snp_annotation/` convert the 50
Log-rank SIS rsIDs into Appendix C, Table C.1 (dbSNP + Ensembl VEP
annotation, grouped into 33 loci, compared with GWAS Catalog Alzheimer's
disease associations). They operate only on the list of rsIDs and on public
databases; they do not read or write individual-level ADNI data.

| Script | Step | Output |
|---|---|---|
| `config.py` | Shared paths and parsing helpers | - |
| `fetch_dbsnp.py` | Fetch dbSNP annotations via NCBI E-utilities | `dbsnp_annotations.csv` |
| `fetch_gwas_catalog.py` | Download GWAS Catalog AD associations and match them to the 50 SNPs | `gwas_ad_matches.csv` |
| `compare_and_summarize.py` | Compare VEP vs dbSNP gene calls; recompute 33 loci and summary | `compare_vep_dbsnp.csv`, `summary.txt` |
| `enrich_and_nearest.py` | Enrich matches and find the nearest AD-associated locus per SNP/locus | `nearest_ad_hit.csv`, `pmid_39998322_rows.csv` |

## Requirements

R packages:

```r
install.packages(c("ADNIMERGE2", "data.table", "dplyr", "tidyr",
                   "doParallel", "parallel", "Icens", "interval",
                   "icenReg", "survival", "energy", "MASS",
                   "ggplot2", "scales", "patchwork", "readr"))
```

`ADNIMERGE2` is the ADNI R package used by `ADNI_Data_Preprocessing.R`.
`plot_adni_manhattan.R` optionally uses `ggrepel` if it is installed.

Python (only for `Plots/plot_gof_sensitivity.py`):

```bash
pip install matplotlib numpy pandas
```

The `snp_annotation` scripts use only the Python standard library.

## Data preparation

### Salary data (Appendix B)

Downloaded from the repository of Zhong et al. (2023):

- https://github.com/tsienchen/ADD-SIS/tree/main/CaseStudy

Place the following files in `CaseStudy/Data/`:

| File | Description |
|---|---|
| `REAL_DATA.csv` | Raw job-posting text data |
| `Words.csv` | Words generated by word segmentation |
| `Word_Filter.csv` | Words kept after low-frequency filtering |
| `WordMatrix_Filter.csv` | Binary design matrix (postings x words) |
| `SelectedWord.csv` | Top selected words with translations |

### ADNI data (Section 4)

ADNI data cannot be redistributed and requires registration at
[adni.loni.usc.edu](https://adni.loni.usc.edu). `ADNI_Data_Preprocessing.R`
builds the two inputs below from the `ADNIMERGE2` tables and an interactively
selected PLINK recode file:

- `CaseStudy/Code/ADNI_Y_Final.rds` (interval-censored outcomes)
- `CaseStudy/Code/ADNI_X_Final.rds` (SNP matrix)

`plot_adni_manhattan.R` additionally needs a PLINK `.bim`/`.map` file that maps
each SNP to its chromosome and base-pair position. By default it looks for

- `CaseStudy/Code/ADNI_1_GWAS_Plink/ADNI_cluster_01_forward_757LONI.bim`

If your file is elsewhere, edit `bim_file` in `plot_adni_manhattan.R`
(a `# PATH:` comment marks this line).

### SNP annotation inputs (Appendix C)

Before running the `snp_annotation` scripts, place these files in
`CaseStudy/Code/snp_annotation/`:

- `ADNI_top50_rsID.txt` - produced by `get_top50_genomewide.R`
- `results_from_Ensembl_Log_rank.txt` - Ensembl VEP (release 116) output for
  the Log-rank SIS top-50 SNPs. The ADD-SIS VEP output
  (`results_from_Ensembl_ADD.txt`) is not used.
- `VEP_Ensembl116_summary_by_SNP.csv` - VEP summary table

All paths in `config.py` are resolved relative to that file's own folder
(see the `# PATH:` comment); dbSNP and GWAS Catalog data are downloaded
automatically by the scripts.

## How to reproduce

Working directories matter. Simulation and plot scripts expect to run from the
**repository root**. The ADNI/Salary case-study scripts and
`get_top50_genomewide.R` expect the working directory `CaseStudy/Code`.

### 1. Simulations (Tables 1-7, Figure 1)

```r
source("Simulation/Example1-3/Code/LM-Log-rank.R")
source("Simulation/Example1-3/Code/GLM-Log-rank.R")
source("Simulation/Example1-3/Code/COX-Log-rank.R")
source("Simulation/Example1-3/Code/Proportion-Tables.R")   # needs the three runs above
source("Simulation/Outlier_Case/Code/Outlier_Case.R")
source("Plots/plot_screening_performance.R")               # Figure 1
```

### 2. Salary case study (Appendix B, Figure B.1)

1. Put the salary data in `CaseStudy/Data/`.
2. Knit `CaseStudy/Code/CaseStudyCode_Log-rank_Salary.Rmd` in RStudio.
3. From the repository root:

```r
source("Plots/plot_salary_cdf.R")
```

### 3. ADNI case study (Section 4, Tables 8-11, Figures 2-3, Appendix C)

1. Run `CaseStudy/Code/ADNI_Data_Preprocessing.R` (working directory
   `CaseStudy/Code`) to create `ADNI_Y_Final.rds` and `ADNI_X_Final.rds`.
2. Knit `CaseStudy/Code/CaseStudyCode_Log-rank_ADNI.Rmd` in RStudio.
3. From the repository root:

```bash
python Plots/plot_gof_sensitivity.py
```

4. From the repository root:

```r
source("Plots/plot_adni_manhattan.R")   # Figure 3 (runs the Log-rank screening automatically)
```

`get_top50_genomewide.R` can also be run on its own from `CaseStudy/Code/`
to regenerate the top-50 table; the plot script sources it automatically.

5. From `CaseStudy/Code/snp_annotation/`:

```bash
python fetch_dbsnp.py
python compare_and_summarize.py
python fetch_gwas_catalog.py download
python fetch_gwas_catalog.py traits
python fetch_gwas_catalog.py match
python enrich_and_nearest.py
```

## Output locations

- `Simulation/Example1-3/Intermediate/{LM,GLM,COX}_Log-rank/Full_res*.csv`
- `Simulation/Example1-3/Tables/Prop_Table_{LM,GLM,COX}.csv`
- `Simulation/Outlier_Case/Intermediate/Outlier_Case*.csv`
- `Screening_Performance_2x3_LargeFont.pdf` (Figure 1)
- `salary_analysis_plot.pdf` (Figure B.1)
- `GoF_Sensitivity_Analysis.pdf` and `GoF_Sensitivity_Analysis.png` (Figure 2)
- `Plots/ADNI_Manhattan_LogrankSIS.png` (Figure 3)
- `CaseStudy/Output/ADNI/adni_gof_by_d.csv` and `CaseStudy/Output/Salary/*.csv`
- `CaseStudy/Code/ADNI_top50_rsID.txt` and `ADNI_top50_genomewide_R.csv`
- `CaseStudy/Code/snp_annotation/*.csv` and `summary.txt` (Table C.1 inputs)

## Notes

- The repository contains only the scripts referenced by the manuscript. A
  standalone ADD-SIS genome-wide script used for an exploratory gene-level
  comparison is not included because that comparison does not appear in the
  paper.
- Syntax was verified for all `.R` files (`R parse`), all `.py` files
  (`py_compile`) and every R code chunk in both `.Rmd` files. Full numerical
  simulations and data analyses were not re-run during packaging.

## References

- Zhong, W., Qian, C., Liu, W., Zhu, L., & Li, R. (2023). Feature screening for
  interval-valued response with application to study association between posted
  salary and required skills. *JASA*, 118, 805-817.
- Fan, J., & Lv, J. (2008). Sure independence screening for ultrahigh dimensional
  feature space. *JRSS-B*, 70, 849-911.
- Sun, J. (1996). A non-parametric test for interval-censored failure time data
  with application to AIDS studies. *Statistics in Medicine*, 15, 1387-1395.
