# Experimental design in R

R code from my graduate coursework in the design and analysis of experiments
(STAT 502, University of Washington). The scripts cover hypothesis tests, power
and sample size, ANOVA, contrasts and multiple comparisons, and block designs.

## Scripts

| Script | Topics |
|---|---|
| `R/01_two_sample_tests.R` | Pooled and Welch t-tests, randomization test, confidence interval for a difference in means |
| `R/02_power_and_sample_size.R` | Power of the paired t-test and the ANOVA F-test, sample size for a target power |
| `R/03_one_way_anova.R` | One-way ANOVA by hand and with `aov()`, residual checks, permutation test for F |
| `R/04_linear_model_matrices.R` | Design matrix under a sum constraint, hat matrix, least squares with matrix algebra |
| `R/05_transformations_and_contrasts.R` | Box-Cox transformation, planned contrasts with Bonferroni and Holm corrections |
| `R/06_blocking_and_paired_designs.R` | Randomized block design vs paired t-test (F = t^2), balance of an incomplete block design |

## Running the code

Open R in the repository folder and run a script, or run everything with

```r
source("run_all.R")
```

Packages: `ggplot2`, `MASS`, `multcomp`.

## Data

The small data sets used in the course are in `data/`: species diversity by zinc
level (`zinc.RDS`), platelet aggregation before and after smoking (`smoke.RDS`),
an incomplete block experiment (`ibd.RDS`), and lentil yields (`lentil.dat`).
