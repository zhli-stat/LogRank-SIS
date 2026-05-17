# Log-rank SIS: Feature Screening for Ultrahigh-Dimensional Interval-Censored Data

This repository contains the R/Python code and data to reproduce all numerical
results in the paper:

> **A Computationally Efficient and Statistically Robust Feature Screening Procedure
> for Ultrahigh-Dimensional Interval-Censored Response**
> *Manuscript under review*, 2026.

---

## Overview

We propose a **Variance-Adjusted Generalized Log-rank Sure Independence Screening
(Log-rank SIS)** procedure for ultrahigh-dimensional Case II interval-censored data
with binary predictors. Compared with the existing ADD-SIS (Zhong et al., 2023),
the proposed method offers three key advantages:

- **Computationally efficient**: Requires only a single global NPMLE fit, reducing
  screening runtime from 3,767.89 s to 2.28 s at *p* = 8,000 (ADNI dataset).
- **Statistically robust**: Uses pooled failure and risk sets, avoiding instability
  from class-imbalanced conditional NPMLE estimation.
- **Model-free**: Avoids screening bias induced by model misspecification.

The core Log-rank SIS function has been implemented and made available in the
**MFSIS** R package on CRAN. The code in this repository reproduces all results
in the paper using the original source code, which provides full transparency
into the implementation details.

---

## R Package

The Log-rank SIS method is also available as a standalone function in the
**MFSIS** R package (Model-Free Sure Independence Screening):

```r
install.packages("MFSIS")
```

> **Note**: The simulation and case study scripts in this repository use the
> **original source code** rather than the MFSIS package directly. This ensures
> full reproducibility of all numerical results in the paper, including
> intermediate steps and benchmark comparisons. Users who wish to apply
> Log-rank SIS to their own data are encouraged to use the MFSIS package
> for convenience.

---

## Repository Structure

```
LogRank-SIS/
├── README.md
│
├── Plots/                                  # Code to reproduce all figures
│   ├── plot_salary_cdf.R                   # Figure 2: Conditional CDF plots (Salary)
│   ├── plot_adni_cdf.R                     # Figure 4: Conditional CDF plots (ADNI)
│   ├── plot_screening_performance.R        # Figure 1: MMS comparison (Simulation)
│   └── plot_gof_sensitivity.py             # Figure 3: GoF sensitivity across submodel sizes
│
├── CaseStudy/                              # Real data applications (Section 4)
│   ├── Code/
│   │   ├── ADNI_Data_Preprocessing.R           # Section 4.3: ADNI data preprocessing
│   │   ├── CaseStudyCode_Log-rank_Salary.Rmd   # Section 4.2: Job advertisement data
│   │   └── CaseStudyCode_Log-rank_ADNI.Rmd     # Section 4.3: ADNI genetic data
│   └── Data/                              # Salary dataset (see Data section below)
│       ├── REAL_DATA.csv                  # Original job posting text data
│       ├── Words.csv                      # 1,043 words extracted from postings
│       ├── Word_Filter.csv                # 677 words after low-frequency filtering
│       ├── WordMatrix_Filter.csv          # Binary design matrix (n x 677)
│       └── SelectedWord.csv              # Top selected words with English translation
│
└── Simulation/                            # Monte Carlo simulations (Section 3)
    ├── Example1-3/
    │   └── Code/
    │       ├── LM-Log-rank.R              # Example 1: Linear model (Table 1)
    │       ├── GLM-Log-rank.R             # Example 2: Poisson regression (Table 2)
    │       ├── COX-Log-rank.R             # Example 3: Cox model (Table 3)
    │       └── Proportion-Tables.R        # Examples 1–3: proportion tables from intermediate results
    └── Outlier_Case/
        └── Code/
            └── Outlier_Case.R             # Sensitivity analysis (Table 4)
```

---

## Data

###  Job Advertisement Data

This dataset was originally collected and processed by Zhong et al. (2023).
Please download it directly from their repository:
https://github.com/tsienchen/ADD-SIS/tree/main/CaseStudy

Place the the `Data` folder in the CaseStudy directory of Zhong et al. (2023) in this folder before running the case study code.

| File | Description |
|---|---|
| `REAL_DATA.csv` | Original job posting text data (*n* = 497 postings) |
| `Words.csv` | 1,043 words generated from the original text via word segmentation |
| `Word_Filter.csv` | 677 words after filtering out those appearing in fewer than 5 postings |
| `WordMatrix_Filter.csv` | Binary design matrix (*n* x 677), where each entry indicates presence (1) or absence (0) of a keyword |
| `SelectedWord.csv` | Top-ranked words selected by ADD-SIS with English translations |

The interval-valued response is the posted monthly salary (RMB), recorded as
salary brackets (e.g., 50k–70k RMB/month).

### ADNI Dataset (Section 4.3)

The Alzheimer's Disease Neuroimaging Initiative (ADNI) data cannot be redistributed.
Data access requires registration at the official portal:
[adni.loni.usc.edu](https://adni.loni.usc.edu)

The study cohort used in this paper consists of:

- *n* = 278 subjects with mild cognitive impairment (MCI)
- 155 (55.8%) interval-censored conversion events; 123 (44.2%) right-censored
- *p* = 530,437 SNPs after MAF >= 0.01 quality control filter
- Benchmark subset: top *p* = 8,000 SNPs by marginal empirical Bernoulli variance

---

## Requirements

### R Packages

```r
install.packages(c("MFSIS", "icenReg", "survival", "doParallel", "foreach",
                   "ggplot2", "dplyr"))
```

### Python (for Figure 3 only)

```python
pip install matplotlib numpy pandas
```

---

## How to Reproduce

### Simulation Results (Tables 1–3, Figure 1)

```r
# Example 1: Linear model -> reproduces Table 1
source("Simulation/Example1-3/Code/LM-Log-rank.R")

# Example 2: Poisson regression -> reproduces Table 2
source("Simulation/Example1-3/Code/GLM-Log-rank.R")

# Example 3: Cox proportional hazards -> reproduces Table 3
source("Simulation/Example1-3/Code/COX-Log-rank.R")

# Proportion tables for Examples 1–3
# Requires intermediate outputs from the three simulation scripts above
source("Simulation/Example1-3/Code/Proportion-Tables.R")
# Saves: Simulation/Example1-3/Tables/Prop_Table_{LM,GLM,COX}.csv
```

### Sensitivity Analysis (Table 4)

```r
source("Simulation/Outlier_Case/Code/Outlier_Case.R")
```

### Case Studies (Section 4)

Before running the ADNI case study, preprocess the raw ADNI data by running:

```r
source("CaseStudy/Code/ADNI_Data_Preprocessing.R")
```

This script processes the raw genotyping data downloaded from the ADNI portal
and prepares the input required by `CaseStudyCode_Log-rank_ADNI.Rmd`.

Open the following `.Rmd` files in **RStudio** and click **Knit**:

| File | Reproduces |
|---|---|
| `CaseStudy/Code/CaseStudyCode_Log-rank_Salary.Rmd` | Tables 5–8, Figure 2 |
| `CaseStudy/Code/CaseStudyCode_Log-rank_ADNI.Rmd` | Tables 9–12, Figures 3–4 |

> The ADNI case study requires the ADNI dataset. See the **Data** section above
> for access instructions.

### Figures

```r
# Figure 1: Screening performance comparison (MMS plot)
source("Plots/plot_screening_performance.R")

# Figure 2: Conditional CDF plots for salary keywords
source("Plots/plot_salary_cdf.R")

# Figure 4: Conditional CDF plots for ADNI
source("Plots/plot_adni_cdf.R")
```

```python
# Figure 3: GoF sensitivity analysis across submodel sizes
python Plots/plot_gof_sensitivity.py
```

---

## References

- Zhong, W., Qian, C., Liu, W., Zhu, L., & Li, R. (2023). Feature screening for
  interval-valued response with application to study association between posted
  salary and required skills. *JASA*, 118, 805–817.
- Fan, J., & Lv, J. (2008). Sure independence screening for ultrahigh dimensional
  feature space. *JRSS-B*, 70, 849–911.
- Sun, J. (1996). A non-parametric test for interval-censored failure time data
  with application to AIDS studies. *Statistics in Medicine*, 15, 1387–1395.