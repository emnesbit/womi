# R Pipeline Documentation

## Overview

The R pipeline has two scripts, which are run in this order:

1. `accuracy_analysis.R`: weight estimation accuracy (Results, Sections 1 and 4)
2. `mep_analysis.R`: MEP analyses (Results, Sections 2–4)

Both scripts follow the order of the Results. All statistics are printed to the console, and figures are saved as PNG files (300 dpi). Tables ST1–ST4 are printed in the format of the Supplementary Material.

**Packages** (see `requirements.R`): tidyverse, readxl, writexl, afex, emmeans, effectsize, broom, patchwork, RColorBrewer.

**Configuration**: set `bidsDir` at the top of each script to your BIDS directory.

---

## Shared conventions

- **Models**: all repeated-measures ANOVAs and ANCOVAs are fitted with `afex::aov_car` (Type III sums of squares). A covariate given in the formula is automatically crossed with all within-subject factors. `factorize = FALSE` keeps covariates numeric.
- **Sphericity**: the helper `rmTable()` prints Mauchly's test and the Greenhouse–Geisser ε, and applies the correction only to effects for which Mauchly's test is significant (p < 0.05).
- **Variables**: `id` is the participant, `accuracy_c` the mean-centred weight estimation accuracy (mean of AO and AOMI per participant).
- **Session order**: even participant IDs started with AO, odd IDs with AOMI.

---

## 1. `accuracy_analysis.R`

### Input

- `sub-XX/ses-{AO,AOMI}/beh/sub-XX_ses-{AO,AOMI}_womi_beh.tsv` (BIDS)
- `derivatives/vmiq_womi.xlsx` (VMIQ-2 scores, from OSF)

### Data preparation

For each session, one row per trial is created with the weight shown, the weight response and the TMS intensity. Control trials, the first (dummy) trial and trials without a response are removed. A trial is correct if the response matches the weight shown. Accuracy is the proportion of correct trials per participant and condition, and is saved to `accuracy_by_condition.xlsx`.

The helper `pct()` prints proportions as percentages with two decimals.

### Section 1: Weight estimation accuracy

| Analysis | Method |
|----------|--------|
| Accuracy vs. chance (1/3) | one-sided one-sample t-test per condition; Cohen's d against chance |
| AO vs. AOMI | paired t-test, Cohen's d; Pearson correlation between conditions |
| Condition × weight response | `accuracy ~ condition * weight + Error(id / (condition * weight))` |
| Weight responses | Bonferroni-corrected post-hoc comparisons, pooled over conditions |
| Response bias | chi-square goodness-of-fit test on response counts |

**Figure 2**: accuracy by weight response in AO and AOMI. Each participant has one colour in both panels; grey lines connect a participant's values; black diamonds show the group mean ± SD; brackets show the post-hoc tests and the condition effect.

### Section 4: Stability of weight estimation accuracy

| Analysis | Method |
|----------|--------|
| Session 1 vs. 2 | paired t-test (Figure S1) |
| AO first vs. AOMI first | Welch t-test on mean accuracy across conditions |
| TMS intensity | `accuracy ~ intensity + Error(id / intensity)` (Figure S2) |
| Imagery vividness | Pearson correlations of AO and AOMI accuracy with the VMIQ-2 mean score and subscales, Holm-corrected over eight tests |

### Output (`derivatives/behavioural_analysis/`)

- `accuracy_by_condition.xlsx`
- `Fig2_accuracy_by_weight.png`
- `FigS1_accuracy_by_session.png`
- `FigS2_accuracy_by_tms_intensity.png`

---

## 2. `mep_analysis.R`

### Input

- `derivatives/behavioural_analysis/accuracy_by_condition.xlsx` (from `accuracy_analysis.R`)
- `derivatives/mep_ztransform/mean_z_scores_by_intensity_and_weight.csv` (from `z_transformation.m`)
- `derivatives/mep_ztransform/mean_z_scores_control_vs_main.csv` (from `z_transformation.m`)
- `derivatives/vmiq_womi.xlsx` (from OSF)
- `derivatives/ADM_FDI_raw_MEPs_for_ANOVA.csv` (from OSF)

### Data preparation

- `dataLong`: one z-scored MEP per participant × condition × weight response × intensity (expected: 18 × 2 × 3 × 5 = 540 rows).
- `curves`: recruitment-curve slope per participant × condition × weight response, the least-squares slope of zMEP on intensity in z per 10% RMT.
- `gain`: slope averaged over weight responses, one value per participant × condition.

### Section 2: Effects of condition and weight (Hypotheses 1 and 2)

| Analysis | Model |
|----------|-------|
| Manipulation check | `zMEP ~ intensity + Error(id / intensity)`; Bonferroni-corrected tests of each intensity against the next |
| Hypotheses 1 and 2 | `slope ~ condition * weight + Error(id / (condition * weight))` |
| AOMI vs. AO | number of participants with a steeper slope in AOMI; paired Cohen's dz |

**Figure 3**: z-scored MEPs across intensities for each weight response, in AO (A) and AOMI (B); individual participants in grey, group mean ± SD in colour.

### Section 3: Exploratory, weight estimation accuracy

| Analysis | Model |
|----------|-------|
| Normality of the covariate | Shapiro–Wilk test |
| Slope ANCOVA (Table ST1) | `slope ~ accuracy_c + condition * weight + Error(id / (condition * weight))` |
| Simple slopes | `slope ~ accuracy` per condition, on `gain` |
| zMEP ANCOVA (Table ST2) | `zMEP ~ accuracy_c + condition * weight * intensity + Error(id / (condition * weight * intensity))` |
| Per intensity (Table ST3) | `zMEP ~ accuracy_c + condition * weight + Error(id / (condition * weight))` at each intensity; simple slopes `zMEP ~ accuracy_c` per condition on participant means across weights; Bonferroni correction over five intensities |

**Figure 4**: (A) recruitment-curve slope against accuracy in AO and AOMI, with linear fits; (B) simple slopes of accuracy on zMEP per intensity, with Bonferroni-adjusted 95% CIs.

### Section 4: Control analyses and robustness

| Analysis | Model |
|----------|-------|
| Imagery ability | slope ANCOVA with the centred VMIQ-2 mean score, or kinaesthetic score, as a second covariate |
| Session order (Table ST4) | `slope ~ order + accuracy_c + condition * weight + Error(id / (condition * weight))`; slope means and SDs by order × condition and order × weight |
| Control trials | paired t-test of z-scored control-trial MEPs, AO vs. AOMI session |
| FDI vs. ADM | `mep ~ muscle * condition * intensity + Error(id / (muscle * condition * intensity))` on raw MEPs (Figure S3, mean ± SD) |
| Residuals | Shapiro–Wilk test and Q–Q plot of the slope ANCOVA residuals (Figure S4) |
| Permutation test | correlation between accuracy and the AO − AOMI slope difference (equivalent to the accuracy × condition interaction), 10,000 permutations, two-sided, `set.seed(1)` |
| Leave-one-participant-out | slope ANCOVA re-run 18 times, each without one participant, with accuracy re-centred on the remaining 17 (Figure S5) |

### Output (`derivatives/ancova_analysis/`)

- `Fig3_zMEP_by_intensity.png`
- `Fig4_accuracy_condition.png`
- `FigS3_mep_by_muscle.png`
- `FigS4_residual_qq_slope.png`
- `FigS5_leave_one_out_slope.png`

---

## Pipeline dependencies

| Script | Needs |
|--------|-------|
| `accuracy_analysis.R` | BIDS behavioural files; `vmiq_womi.xlsx` |
| `mep_analysis.R` | `accuracy_by_condition.xlsx`; `z_transformation.m` outputs; `vmiq_womi.xlsx`; `ADM_FDI_raw_MEPs_for_ANOVA.csv` |

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | April 2026 | Initial BIDS-compatible version |
| 2.0 | November 2026 | Rewritten for the revised manuscript (see README) |

---

## Contact

Emma Nesbit: nesbit@cbs.mpg.de · GitHub Issues for bug reports
