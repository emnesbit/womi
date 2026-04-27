# R Pipeline Documentation

## Overview

The R pipeline processes TMS-EMG and behavioural data through two main stages:
1. **Behavioural Accuracy Analysis** (`accuracy.R`)
2. **MEP ANCOVA Analysis** (`main_3way_ANCOVA_with_accuracy.R`)

---

## Stage 1: Behavioural Accuracy Analysis (`accuracy.R`)

### Purpose
Extracts and analyses behavioural accuracy data from BIDS-formatted files, tests for session order and condition effects on weight discrimination performance, and generates publication-ready visualisations.

### Required Packages

```r
pacman::p_load(
  tidyverse,      # Data manipulation and plotting
  haven,          # Data import
  psych,          # Descriptive statistics
  writexl,        # Excel output
  afex,           # Repeated measures ANOVA
  emmeans,        # Estimated marginal means
  patchwork,      # Multi-panel figures
  effectsize,     # Cohen's d and eta squared
  RColorBrewer    # Colour palettes
)
```

### Input Files

The script reads directly from BIDS-formatted behavioural data:

```
/data/path/to/bids_directory/
├── sub-01/
│   ├── ses-AO/
│   │   └── beh/
│   │       └── sub-01_ses-AO_womi_beh.tsv
│   └── ses-AOMI/
│       └── beh/
│           └── sub-01_ses-AOMI_womi_beh.tsv
├── sub-02/
└── ...
```

**Behavioural File Format** (`*_beh.tsv`):
- `trial`: Trial number
- `event_type`: Event category (e.g., `weight_in_video`, `weight_estimation_response`, `Control_TMS_pulse`)
- `code`: Event value

### Configuration Parameters

```r
bidsDirectory  <- "/data/path/to/bids_directory/"
outputDirectory <- file.path(bidsDirectory, "derivatives", "behavioural_analysis")

participants   <- 1:18
sessions       <- c("AO", "AOMI")

# Session order: odd participants start with AOMI, even with AO
sessionOrder   <- ifelse(participants %% 2 == 0, "AO", "AOMI")
```

### Processing Steps

#### Part 1: Data Extraction
```r
aoData   <- extractBehaviouralData(p, "AO")
aomiData <- extractBehaviouralData(p, "AOMI")
```

For each participant and session, `extractBehaviouralData()`:
- Reads the BIDS `.tsv` behavioural file
- Identifies control trials via the presence of `Control_TMS_pulse` event type
- Extracts `weight_in_video` (weight shown) and `weight_estimation_response` (participant response) for main trials
- Excludes trial 1 (dummy), control trials, and rows with missing weight or response values
- Calculates trial-level accuracy (`correct = weightShown == weightResponse`)
- Returns overall accuracy, trial counts, and per-weight-response accuracy

**Trial Exclusion Criteria**:
| Criterion | Reason |
|-----------|--------|
| `isControl == TRUE` | Control trials have no weight judgment |
| `trial == 1` | First trial excluded as dummy |
| `is.na(weightShown)` | Missing stimulus information |
| `is.na(weightResponse)` | Missing participant response |

#### Part 2: Summary Tables

Per-participant accuracy is collated into three summary tables, each saved as `.xlsx`:

| Output File | Contents |
|-------------|----------|
| `accuracy_by_condition.xlsx` | Per-participant accuracy for AO and AOMI sessions |
| `accuracy_by_session_order.xlsx` | Accuracy by first vs. second session (session order) |
| `accuracy_by_weight.xlsx` | Accuracy broken down by weight response (1–3) for each session |

#### Part 3: Session Order Analysis (Session 1 vs. Session 2)

Tests whether accuracy differs between the first and second experimental session, regardless of which condition was completed first:

```r
# Normality check on difference scores
shapiro.test(Session1 - Session2)

# Paired t-test
t.test(Session1, Session2, paired = TRUE, conf.level = 0.95)

# Effect size
cohens_d(Session1, Session2, paired = TRUE)
```

**Output**: `session_order_plot.png`  
Paired line plot showing individual participant trajectories from Session 1 to Session 2, with points coloured by condition (AO/AOMI) and group mean ± 95% CI overlaid.

#### Part 4: Condition Analysis (AO vs. AOMI)

Tests whether accuracy differs between the Action Observation (AO) and Action Observation with Motor Imagery (AOMI) conditions:

```r
# Normality check on difference scores
shapiro.test(AO - AOMI)

# Paired t-test
t.test(AO, AOMI, paired = TRUE, conf.level = 0.95)

# Effect size
cohens_d(AO, AOMI, paired = TRUE)
```

**Output**: `condition_comparison_plot.png`, `condition_comparison_for_panel.png`  
Paired line plots with individual trajectories and group mean ± 95% CI. The panel version omits the title.

#### Part 5: Weight Response Analysis

Tests the effect of weight response category (Light/Moderate/Heavy) and condition on accuracy using repeated measures ANOVA:

**Two-Way Repeated Measures ANOVA** (condition × weight):
```r
aov_car(
  accuracy ~ condition * weight + 
    Error(participantID / (condition * weight)),
  data = weightData,
  type = 3
)
```

- Mauchly's test of sphericity is run for the `weight` factor and the `condition × weight` interaction (condition has only 2 levels and is exempt)
- Greenhouse–Geisser correction applied where sphericity is violated (p < .05)
- Partial η² reported as effect size

**One-Way Repeated Measures ANOVA** (weight only):
```r
aov_car(
  accuracy ~ weight + Error(participantID / weight),
  data = weightData,
  type = 3,
  fun_aggregate = mean
)
```

**Post-Hoc Comparisons** (Bonferroni-corrected pairwise):
```r
emmeans(anova1Way, pairwise ~ weight, adjust = "bonferroni")
```

**Chi-Square Test** (uniform response distribution):
```r
chisq.test(observed, p = c(1/3, 1/3, 1/3))
```

**Outputs**: `weight_by_condition_plot.png`, `weight_responses_plot.png`

#### Part 6: Combined Manuscript Figure

A two-panel figure combining the session order and weight response plots:

```r
behaviouralManuscriptFigure <- orderPlotA | weightPlotB
```

**Output**: `behavioural_manuscript_figure.png` (14 × 6 inches, 300 dpi)

| Panel | Content |
|-------|---------|
| A (left) | Session order effects (S1 vs. S2) |
| B (right) | Accuracy by weight response (Light/Moderate/Heavy) |

### Output Files

Saved to: `/derivatives/behavioural_analysis/`

```
accuracy_by_condition.xlsx
accuracy_by_session_order.xlsx
accuracy_by_weight.xlsx
session_order_plot.png
condition_comparison_plot.png
condition_comparison_for_panel.png
weight_by_condition_plot.png
weight_responses_plot.png
behavioural_manuscript_figure.png
```

### Colour Conventions

| Condition / Weight | Hex Code |
|--------------------|----------|
| AO | `#F8766D` (red) |
| AOMI | `#00BFC4` (blue) |
| Light weight | `#FDD0A2` |
| Moderate weight | `#FDAE6B` |
| Heavy weight | `#E6550D` |

### Expected Console Output

```
=== PART 1: EXTRACTING BEHAVIOURAL DATA ===

Processing participant 01...
Processing participant 02...
...
=== Data extraction complete ===
Total participants: 18
Total trials: XXXX

=== PART 2: CREATING SUMMARY TABLES ===
Mean accuracies:
  AO:   X.XXX (SD: X.XXX)
  AOMI: X.XXX (SD: X.XXX)

=== PART 3: SESSION ORDER ANALYSIS (S1 vs S2) ===
[Descriptives, Shapiro-Wilk, t-test, Cohen's d printed]

=== PART 4: CONDITION ANALYSIS (AO vs AOMI) ===
[Descriptives, Shapiro-Wilk, t-test, Cohen's d printed]

=== PART 5: WEIGHT RESPONSE ANALYSIS ===
[Descriptives, Mauchly's test, ANOVA tables, post-hoc comparisons printed]

=== PART 6: CREATING COMBINED MANUSCRIPT FIGURE ===
Two-panel figure saved as: behavioural_manuscript_figure.png
  - Panel A (left): Session order effects
  - Panel B (right): Accuracy by weight response
```

---

## Stage 2: MEP ANCOVA Analysis (`main_3way_ANCOVA_with_accuracy.R`)

### Purpose
Performs a three-way repeated measures ANCOVA (condition × weight response × TMS intensity) on z-transformed MEP amplitudes, with behavioural accuracy as a between-subject covariate. Conducts assumption checks, sphericity corrections, post-hoc analyses, and intensity-specific follow-up models.

### Required Packages

```r
pacman::p_load(
  tidyverse,    # Data manipulation and plotting
  readxl,       # Excel import
  afex,         # Repeated measures ANCOVA
  emmeans,      # Estimated marginal means and simple slopes
  effectsize,   # Partial eta squared
  car,          # Multivariate tests and VIF
  lme4,         # Linear mixed models
  lmerTest,     # p-values for mixed models
  patchwork,    # Multi-panel figures
  RColorBrewer  # Colour palettes
)
```

### Input Files

This script requires outputs from both Stage 1 (`accuracy.R`) and the MATLAB z-transformation pipeline (`z-transformation.m`):

```
/derivatives/
├── behavioural_analysis/
│   └── accuracy_by_condition.xlsx          ← From accuracy.R
└── mep_ztransform/
    └── mean_z_scores_by_intensity_and_weight.csv  ← From z-transformation.m
```

**Z-Scores File Format** (`mean_z_scores_by_intensity_and_weight.csv`):

Columns follow the pattern `{Session}_{TrialType}_{Weight}_{Intensity}`:
- Control conditions: `AO_Control_1` through `AOMI_Control_5` (2 sessions × 5 intensities)
- Main trial conditions: `AO_WtResp_1_1` through `AOMI_WtResp_3_5` (2 sessions × 3 weights × 5 intensities)

**Accuracy File Format** (`accuracy_by_condition.xlsx`):

| Column | Description |
|--------|-------------|
| participant | Participant number (1–18) |
| AO | Proportion correct in AO session |
| AOMI | Proportion correct in AOMI session |

### Configuration Parameters

```r
bidsDirectory    <- "/data/path/to/bids_directory/"
behaviouralDir   <- file.path(bidsDirectory, "derivatives", "behavioural_analysis")
zScoresDir       <- file.path(bidsDirectory, "derivatives", "mep_ztransform")
outputDirectory  <- file.path(bidsDirectory, "derivatives", "ancova_analysis")
```

### Processing Steps

#### Part 1: Load and Merge Data

```r
# Load accuracy (averaged across sessions)
dfAccuracy <- read_excel(accuracyFile) %>%
  mutate(accuracy = (AO + AOMI) / 2) %>%
  select(participant, accuracy)

# Load z-transformed MEP means
dfZScores <- read_csv(zScoresFile)

# Merge by participant ID
mergedData <- dfZScores %>%
  left_join(dfAccuracy, by = c("participant_num" = "participant"))
```

**Output**: `merged_data.csv`

#### Part 2: Reshape to Long Format

```r
dataLong <- mergedData %>%
  pivot_longer(
    cols = matches("^(AO|AOMI)_WtResp_\\d+_\\d+$"),
    names_to = c("condition", "weight", "intensity"),
    names_pattern = "(AO|AOMI)_WtResp_(\\d+)_(\\d+)",
    values_to = "zMEP"
  ) %>%
  mutate(
    condition  = factor(condition, levels = c("AO", "AOMI")),
    weight     = factor(weight, levels = c("1","2","3"), 
                        labels = c("Light","Moderate","Heavy")),
    intensity  = factor(intensity, levels = c("1","2","3","4","5"),
                        labels = c("90%","100%","110%","120%","130%")),
    accuracy_c = scale(accuracy, center = TRUE, scale = FALSE)[,1]
  )
```

Expected observations: 540 (18 participants × 2 conditions × 3 weights × 5 intensities). Missing cells (due to no trials in a combination) are dropped via `filter(!is.na(zMEP))`.

The covariate `accuracy_c` is mean-centred to reduce multicollinearity and improve interpretability of main effects.

#### Part 3: Assumption Checks

Four assumptions are tested and their outputs saved as plots:

**1. Normality of the Covariate**
```r
shapiro.test(dfAccuracy$accuracy)
```
Output: `assumption_accuracy_distribution.png`

**2. Multicollinearity (VIF)**
```r
modelVIF <- lm(zMEP ~ accuracy + condition + weight + intensity, data = dataLong)
vif(modelVIF)
```
VIF values > 10 trigger a warning. All factors are expected to be below this threshold.

**3. Linearity (Covariate–DV Relationship)**

Scatter plot of accuracy against z-MEP, with separate regression lines per condition:  
Output: `assumption_linearity.png`

**4. Homogeneity of Regression Slopes**
```r
modelSlopes <- lm(
  zMEP ~ accuracy + condition + weight + intensity +
    accuracy:condition + accuracy:weight + accuracy:intensity,
  data = dataLong
)
anova(modelSlopes)
```

If any interaction involving `accuracy` is significant (p < .05), those interaction terms are retained in the ANCOVA model (`includeAccuracyInteractions <- TRUE`). This determination is made automatically and governs the model formula used in Part 4.

#### Part 4: Repeated Measures ANCOVA

The model formula depends on the outcome of the homogeneity check:

**If slopes are homogeneous** (no significant interactions with accuracy):
```r
aov_car(
  zMEP ~ accuracy_c + condition * weight * intensity +
    Error(Participant_ID / (condition * weight * intensity)),
  data = dataLong, type = 3
)
```

**If slopes are heterogeneous** (significant accuracy interactions detected):
```r
aov_car(
  zMEP ~ accuracy_c + accuracy_c:condition + accuracy_c:intensity +
    condition * weight * intensity +
    Error(Participant_ID / (condition * weight * intensity)),
  data = dataLong, type = 3
)
```

**Sphericity**:

Mauchly's test is computed via a multivariate linear model (`car::Anova`) across all 30 within-subject cells. Greenhouse–Geisser correction is applied to any effect for which sphericity is violated:

```r
print(anova(anovaModel, correction = "GG", es = "pes"))
```

Effect sizes are reported as partial η²:
```r
eta_squared(anovaModel, ci = 0.95)
```

**Output**: `ancova_results.txt` (full results including assumption checks, Mauchly's tests, ANCOVA table, and effect sizes)

#### Part 5: Post-Hoc Analyses

**Pairwise comparisons for TMS intensity** (Bonferroni-corrected):
```r
emmeans(anovaModel, pairwise ~ intensity, adjust = "bonferroni")
```
Output: `pairwise_intensity.csv`

**Estimated marginal means for TMS intensity**:
```r
emmeans(anovaModel, ~ intensity)
```

#### Part 6a: Follow-Up Analysis at Highest Intensity (130% RMT)

A linear mixed model tests the accuracy × condition interaction at 130% RMT:

```r
lmer(
  zMEP ~ accuracy_c * condition + (1 | Participant_ID),
  data = dataHighest, REML = FALSE
)
```

Simple slopes are extracted per condition:
```r
emtrends(modelHighest, ~ condition, var = "accuracy_c")
```

**Output**: `simple_slopes_highest_intensity.csv`, `interaction_highest_intensity.png`

#### Part 6b: Accuracy × Condition Interaction Across All Intensities

The accuracy × condition interaction is tested at each TMS intensity level separately using linear mixed models:

```r
for (intLevel in c("90%", "100%", "110%", "120%", "130%")) {
  lmer(zMEP ~ accuracy_c * condition + (1 | Participant_ID), 
       data = filter(dataLong, intensity == intLevel), REML = FALSE)
}
```

For each intensity, the following are extracted and stored:

| Metric | Description |
|--------|-------------|
| `F_value` | F-statistic for accuracy × condition interaction |
| `df_num`, `df_denom` | Numerator and denominator degrees of freedom |
| `p_value` | p-value for the interaction |
| `effect_size` | Partial η² (SS_effect / (SS_effect + SS_residual)) |
| `AO_slope`, `AO_slope_se`, `AO_slope_p` | Simple slope for AO condition |
| `AOMI_slope`, `AOMI_slope_se`, `AOMI_slope_p` | Simple slope for AOMI condition |
| `slope_difference` | AO slope − AOMI slope |

**Outputs**:

| File | Description |
|------|-------------|
| `interaction_by_intensity.csv` | Summary table of interaction statistics per intensity |
| `interaction_all_intensities_faceted.png` | Faceted scatter plot across all intensities |
| `interaction_all_intensities_grid.png` | Grid of individual intensity plots |
| `effect_sizes_by_intensity.png` | Bar plot of partial η² with p-values annotated |
| `slopes_by_intensity.png` | Line plot of simple slopes (AO vs. AOMI) across intensities |

#### Part 7: Visualisations

All plots use a common publication theme (`pubTheme`) with 14pt base font size. The following are generated:

| Plot | File | Description |
|------|------|-------------|
| Main effect of condition | `main_effect_condition.png` | Violin + boxplot by AO/AOMI |
| Main effect of TMS intensity | `main_effect_intensity.png` | Violin + boxplot across intensity levels |
| Interaction at 130% | `interaction_highest_intensity.png` | Scatter + regression lines by condition |
| Combined poster plot | `combined_poster_plot.png` | Condition + intensity plots side by side |
| Marginal means: intensity | `marginal_means_intensity.png` | Adjusted means with 95% CI across intensities |
| Manuscript figure | `manuscript_figure_combined.png` | Three-panel figure for publication |

**Manuscript Figure Panel Layout** (`manuscript_figure_combined.png`, 16 × 11 inches):

| Panel | Content | Source |
|-------|---------|--------|
| A (top left) | Behavioural accuracy: AO vs. AOMI | `accuracy_by_condition.xlsx` |
| B (top right) | Z-MEP by TMS intensity (violin plot) | `dataLong` |
| C (bottom, full width) | Accuracy × Condition interaction, faceted by intensity | `dataLong` |

```r
manuscriptFigure <- (conditionPlotA | intensityPlotB) /
  wrap_elements(full = interactionPlotC) +
  plot_layout(heights = c(1, 1.3))
```

### Output Files

Saved to: `/derivatives/ancova_analysis/`

```
merged_data.csv
ancova_results.txt
pairwise_intensity.csv
simple_slopes_highest_intensity.csv
interaction_by_intensity.csv
assumption_accuracy_distribution.png
assumption_linearity.png
main_effect_condition.png
main_effect_intensity.png
interaction_highest_intensity.png
interaction_all_intensities_faceted.png
interaction_all_intensities_grid.png
effect_sizes_by_intensity.png
slopes_by_intensity.png
combined_poster_plot.png
marginal_means_intensity.png
manuscript_figure_combined.png
```

### Expected Console Output

```
=== WOMI MEP ANCOVA Analysis ===
Output directory: .../ancova_analysis

=== PART 1: LOADING DATA ===
Loaded accuracy data: 18 participants
Loaded z-scores data: 18 participants
Merged data saved.

=== PART 2: RESHAPING DATA ===
Reshaped data: XXX observations
Expected: 540 (18 × 2 × 3 × 5)
Data structure verification:
  Conditions: AO, AOMI
  Weights: Light, Moderate, Heavy
  Intensities: 90%, 100%, 110%, 120%, 130%

=== PART 3: ASSUMPTION CHECKS ===
1. Normality of accuracy (covariate): [Shapiro-Wilk printed]
2. Multicollinearity check (VIF): [VIF values printed]
   ✓ No multicollinearity issues (all VIF < 10)
3. Linearity check (accuracy × zMEP by condition): [plot saved]
4. Homogeneity of regression slopes: [ANOVA table printed]
   ✓ Homogeneity of slopes satisfied for all factors (all p ≥ .05)

=== PART 4: REPEATED MEASURES ANCOVA ===
=== Mauchly's Test of Sphericity ===
[Formatted W and p-values printed per effect]
=== ANCOVA Results (Greenhouse-Geisser corrected) ===
[ANCOVA table printed]
=== Effect Sizes (Partial η²) ===
[Effect size table printed]

=== PART 5: POST-HOC ANALYSES ===
1. Pairwise comparisons for TMS intensity (Bonferroni-corrected):
[Contrast table printed]
4. Estimated marginal means for TMS intensity:
[EMM table printed]

=== PART 6A: FOLLOW-UP ANALYSIS AT HIGHEST INTENSITY ===
Observations at 130% intensity: XXX
=== Mixed Model Results (130% intensity) ===
[Model summary and ANOVA table printed]
Simple slopes at 130% intensity:
[Slopes table with CIs and p-values printed]

=== PART 6B: ACCURACY × CONDITION INTERACTION ACROSS ALL INTENSITIES ===
--- Analyzing 90% RMT ---
--- Analyzing 100% RMT ---
...
--- Analyzing 130% RMT ---
=== Summary Table: Accuracy × Condition Interaction by Intensity ===
[Results table printed]
=== Intensity-Specific Analysis Complete ===

=== PART 7: GENERATING VISUALISATIONS ===
[Plots saved]

=== ANALYSIS COMPLETE ===
Results saved to: .../ancova_analysis
```

---

## Function Reference

### accuracy.R Functions

#### `extractBehaviouralData(participant, session)`
Reads a BIDS behavioural `.tsv` file and returns accuracy metrics for the specified participant and session.

**Arguments**:
- `participant`: Integer participant number (e.g., `1`)
- `session`: Character session label (`"AO"` or `"AOMI"`)

**Returns**: A named list with elements:
- `accuracy`: Overall proportion correct
- `totalTrials`: Number of valid main trials
- `correctTrials`: Number of correct trials
- `responseAccuracy`: Data frame of accuracy by weight response category (1–3)
- `trialData`: Trial-level data frame with `weightShown`, `weightResponse`, `correct`

#### `calcDescriptives(data, varName)`
Calculates mean, SD, and N for a given variable using tidy evaluation.

**Returns**: Single-row data frame with `mean`, `sd`, `n`.

---

### main_3way_ANCOVA_with_accuracy.R — Key Operations

This script does not define named helper functions but relies on the following pipeline operations:

#### Data Loading
```r
read_excel(accuracyFile)   # Loads behavioural accuracy from Stage 1
read_csv(zScoresFile)      # Loads z-transformed MEP means from MATLAB pipeline
```

#### Model Fitting
```r
aov_car(formula, data, type = 3)   # Repeated measures ANCOVA (afex)
lmer(formula, data, REML = FALSE)  # Linear mixed model (lme4/lmerTest)
```

#### Post-Hoc and Slopes
```r
emmeans(model, pairwise ~ intensity, adjust = "bonferroni")
emtrends(model, ~ condition, var = "accuracy_c")
```

#### Sphericity Testing
```r
car::Anova(mlm_model, idata = idata, idesign = ~condition * weight * intensity, type = 3)
summary(model_anova_sph, multivariate = FALSE)$sphericity.tests
```

---

## Code Modification Guidelines

### Changing Quality Thresholds for Accuracy Exclusions

```r
# Exclude participants below chance performance
accuracyData <- dfAccuracy %>% filter(accuracy > 0.33)
```

### Changing the Significance Threshold for Sphericity Correction

```r
# More conservative threshold
use_GG_effects <- rownames(sphericity_sph)[sphericity_sph[, "p-value"] < 0.01]
# Original: p < 0.05
```

### Modifying Plot Aesthetics

The publication theme is defined inline in Part 7 of `main_3way_ANCOVA_with_accuracy.R`:

```r
pubTheme <- theme_minimal(base_size = 14) +
  theme(
    plot.title  = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title  = element_text(size = 13),
    axis.text   = element_text(size = 12),
    legend.text = element_text(size = 11)
  )
```

To apply consistent styling across both scripts, extract this theme into a shared `theme_setup.R` source file.

---

## Pipeline Dependencies

| Stage | Script | Depends On |
|-------|--------|------------|
| 1 | `accuracy.R` | BIDS behavioural files |
| 2 | `main_3way_ANCOVA_with_accuracy.R` | `accuracy.R` outputs, MATLAB `z-transformation.m` outputs |

The ANCOVA script will stop with an informative error if either upstream file is missing:

```r
if (!file.exists(accuracyFile)) {
  stop("Accuracy file not found. Run behavioural analysis script first.")
}
if (!file.exists(zScoresFile)) {
  stop("Z-scores file not found. Run z-transformation script first.")
}
```

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | April 2026 | Initial BIDS-compatible version |

---

## Contact & Support

For questions about the R pipeline:
- Emma Nesbit: nesbit@cbs.mpg.de
- GitHub Issues for bug reports 