# ============================================================================
# WOMI MEP ANCOVA Analysis - BIDS Pipeline
# ============================================================================
# This script:
# 1. Performs repeated measures ANCOVA with behavioural accuracy as covariate
# 2. Tests condition × weight × intensity effects on z-transformed MEPs
# 3. Conducts assumption checks and post-hoc analyses
# 4. Generates publication-ready visualisations
# ============================================================================

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, readxl, afex, emmeans, effectsize, car, 
  lme4, lmerTest, patchwork, RColorBrewer
)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Set directories
bidsDirectory <- "/path/to/your/bids_directory/"
behaviouralDir <- file.path(bidsDirectory, "derivatives", "behavioural_analysis")
zScoresDir <- file.path(bidsDirectory, "derivatives", "mep_ztransform")
outputDirectory <- file.path(bidsDirectory, "derivatives", "ancova_analysis")

# Create output directory
if (!dir.exists(outputDirectory)) {
  dir.create(outputDirectory, recursive = TRUE)
}

cat("=== WOMI MEP ANCOVA Analysis ===\n")
cat(sprintf("Output directory: %s\n", outputDirectory))

# ============================================================================
# PART 1: LOAD AND MERGE DATA
# ============================================================================

cat("\n=== PART 1: LOADING DATA ===\n")

# Load behavioural accuracy data
accuracyFile <- file.path(behaviouralDir, "accuracy_by_condition.xlsx")
if (!file.exists(accuracyFile)) {
  stop("Accuracy file not found. Run behavioural analysis script first.")
}
dfAccuracy <- read_excel(accuracyFile) %>%
  mutate(accuracy = (AO + AOMI) / 2) %>%
  select(participant, accuracy)

cat(sprintf("Loaded accuracy data: %d participants\n", nrow(dfAccuracy)))

# Load z-transformed MEP data
zScoresFile <- file.path(zScoresDir, "mean_z_scores_by_intensity_and_weight.csv")
if (!file.exists(zScoresFile)) {
  stop("Z-scores file not found. Run z-transformation script first.")
}
dfZScores <- read_csv(zScoresFile, show_col_types = FALSE)

# Extract numeric participant ID from Participant column
dfZScores <- dfZScores %>%
  mutate(
    participant_num = as.numeric(str_extract(Participant, "\\d+"))
  )

cat(sprintf("Loaded z-scores data: %d participants\n", nrow(dfZScores)))

# Merge datasets
mergedData <- dfZScores %>%
  left_join(dfAccuracy, by = c("participant_num" = "participant")) %>%
  select(-participant_num) %>%
  rename(Participant_ID = Participant)

# Save merged data
write_csv(mergedData, file.path(outputDirectory, "merged_data.csv"))
cat("Merged data saved.\n")

# ============================================================================
# PART 2: RESHAPE TO LONG FORMAT
# ============================================================================

cat("\n=== PART 2: RESHAPING DATA ===\n")

# Reshape to long format for repeated measures ANCOVA
dataLong <- mergedData %>%
  pivot_longer(
    cols = matches("^(AO|AOMI)_WtResp_\\d+_\\d+$"),
    names_to = c("condition", "weight", "intensity"),
    names_pattern = "(AO|AOMI)_WtResp_(\\d+)_(\\d+)",
    values_to = "zMEP"
  ) %>%
  mutate(
    condition = factor(condition, levels = c("AO", "AOMI")),
    weight = factor(weight, levels = c("1", "2", "3"), 
                    labels = c("Light", "Moderate", "Heavy")),
    intensity = factor(intensity, levels = c("1", "2", "3", "4", "5"),
                       labels = c("90%", "100%", "110%", "120%", "130%")),
    Participant_ID = factor(Participant_ID),
    accuracy_c = scale(accuracy, center = TRUE, scale = FALSE)[,1]  # Centered accuracy
  ) %>%
  filter(!is.na(zMEP))

cat(sprintf("Reshaped data: %d observations\n", nrow(dataLong)))
cat(sprintf("Expected: %d (18 × 2 × 3 × 5)\n", 18 * 2 * 3 * 5))

# Verify data structure
cat("\nData structure verification:\n")
cat(sprintf("  Conditions: %s\n", paste(levels(dataLong$condition), collapse = ", ")))
cat(sprintf("  Weights: %s\n", paste(levels(dataLong$weight), collapse = ", ")))
cat(sprintf("  Intensities: %s\n", paste(levels(dataLong$intensity), collapse = ", ")))

# ============================================================================
# PART 3: ASSUMPTION CHECKS
# ============================================================================

cat("\n=== PART 3: ASSUMPTION CHECKS ===\n")

# 1. Normality of covariate
cat("1. Normality of accuracy (covariate):\n")
shapiroAccuracy <- shapiro.test(dfAccuracy$accuracy)
print(shapiroAccuracy)

accuracyHist <- ggplot(dfAccuracy, aes(x = accuracy)) +
  geom_histogram(bins = 10, fill = "skyblue", colour = "black", alpha = 0.7) +
  geom_density(aes(y = after_stat(count)), alpha = 0.3, fill = "blue") +
  labs(title = "Distribution of Accuracy", x = "Accuracy", y = "Frequency") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12)
  )

ggsave(file.path(outputDirectory, "assumption_accuracy_distribution.png"), 
       accuracyHist, width = 6, height = 5, dpi = 300)

# 2. Multicollinearity (VIF)
cat("\n2. Multicollinearity check (VIF):\n")
modelVIF <- lm(zMEP ~ accuracy + condition + weight + intensity, data = dataLong)
vifValues <- vif(modelVIF)
print(vifValues)

if (any(vifValues > 10)) {
  cat("\n⚠ WARNING: High multicollinearity detected (VIF > 10)\n\n")
} else {
  cat("\n✓ No multicollinearity issues (all VIF < 10)\n\n")
}

# 3. Linearity: Covariate-DV relationship
cat("3. Linearity check (accuracy × zMEP by condition):\n")

linearityPlot <- ggplot(dataLong, aes(x = accuracy, y = zMEP, colour = condition)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  labs(
    title = "Linearity: Accuracy × Condition on zMEP",
    x = "Accuracy (raw)", 
    y = "Z-transformed MEP",
    colour = "Condition"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12)
  )

ggsave(file.path(outputDirectory, "assumption_linearity.png"), 
       linearityPlot, width = 8, height = 6, dpi = 300)
print(linearityPlot)

# 4. Homogeneity of regression slopes
cat("\n4. Homogeneity of regression slopes:\n")
cat("Testing accuracy interactions with within-subject factors:\n\n")

modelSlopes <- lm(
  zMEP ~ accuracy + condition + weight + intensity +
    accuracy:condition + accuracy:weight + accuracy:intensity,
  data = dataLong
)
anovaSlopes <- anova(modelSlopes)
print(anovaSlopes)

sigInteractions <- c()
if (anovaSlopes["accuracy:condition", "Pr(>F)"] < 0.05) {
  sigInteractions <- c(sigInteractions, "accuracy × condition")
}
if (anovaSlopes["accuracy:weight", "Pr(>F)"] < 0.05) {
  sigInteractions <- c(sigInteractions, "accuracy × weight")
}
if (anovaSlopes["accuracy:intensity", "Pr(>F)"] < 0.05) {
  sigInteractions <- c(sigInteractions, "accuracy × intensity")
}

if (length(sigInteractions) > 0) {
  cat("\n⚠ WARNING: Homogeneity of slopes violated for:\n")
  cat(paste("  -", sigInteractions, collapse = "\n"))
  cat("\n   → These interactions will be retained in the model\n\n")
  includeAccuracyInteractions <- TRUE
} else {
  cat("\n✓ Homogeneity of slopes satisfied for all factors (all p ≥ .05)\n\n")
  includeAccuracyInteractions <- FALSE
}

# ============================================================================
# PART 4: REPEATED MEASURES ANCOVA
# ============================================================================

cat("\n=== PART 4: REPEATED MEASURES ANCOVA ===\n")

if (includeAccuracyInteractions) {
  cat("Running model with significant accuracy interactions retained\n\n")
  
  anovaModel <- aov_car(
    zMEP ~ accuracy_c + accuracy_c:condition + accuracy_c:intensity + 
      condition * weight * intensity + 
      Error(Participant_ID / (condition * weight * intensity)),
    data = dataLong,
    type = 3,
    factorize = FALSE
  )
} else {
  cat("Running model with accuracy as covariate only (no interactions)\n\n")
  
  anovaModel <- aov_car(
    zMEP ~ accuracy_c + condition * weight * intensity + 
      Error(Participant_ID / (condition * weight * intensity)),
    data = dataLong,
    type = 3,
    factorize = FALSE
  )
}

# Mauchly's tests for sphericity
cat("\n=== Mauchly's Test of Sphericity ===\n")

dataWide <- dataLong %>%
  mutate(
    intensity_clean = case_when(
      intensity == "90%" ~ "I90",
      intensity == "100%" ~ "I100",
      intensity == "110%" ~ "I110",
      intensity == "120%" ~ "I120",
      intensity == "130%" ~ "I130"
    )
  ) %>%
  unite("cell", condition, weight, intensity_clean, sep = "_") %>%
  select(Participant_ID, cell, zMEP) %>%
  pivot_wider(names_from = cell, values_from = zMEP)

allCombos <- expand.grid(
  condition = c("AO", "AOMI"),
  weight = c("Light", "Moderate", "Heavy"),
  intensity = c("I90", "I100", "I110", "I120", "I130")
) %>%
  mutate(colname = paste(condition, weight, intensity, sep = "_"))

mlm_cols <- allCombos$colname[allCombos$colname %in% names(dataWide)]
mlm_formula <- as.formula(paste("cbind(", paste("`", mlm_cols, "`", sep = "", collapse = ", "), ") ~ 1"))
mlm_model <- lm(mlm_formula, data = dataWide)

idata <- expand.grid(
  condition = factor(c("AO", "AOMI")),
  weight = factor(c("Light", "Moderate", "Heavy")),
  intensity = factor(c("I90", "I100", "I110", "I120", "I130"))
)

model_anova_sph <- car::Anova(mlm_model, idata = idata, 
                              idesign = ~condition * weight * intensity, type = 3)
mauchly_summary_sph <- summary(model_anova_sph, multivariate = FALSE)

if (!is.null(mauchly_summary_sph$sphericity.tests)) {
  sphericity_sph <- mauchly_summary_sph$sphericity.tests
  print(sphericity_sph)
  
  # Extract and format Mauchly's test results for reporting
  cat("\n=== Formatted Mauchly's Test Results ===\n")
  for (i in 1:nrow(sphericity_sph)) {
    effect_name <- rownames(sphericity_sph)[i]
    W_stat <- sphericity_sph[i, "Test statistic"]
    p_val <- sphericity_sph[i, "p-value"]
    
    cat(sprintf("%s: W = %.3f, p = %.3f\n", effect_name, W_stat, p_val))
  }
  
  use_GG_effects <- rownames(sphericity_sph)[sphericity_sph[, "p-value"] < 0.05]
  
  if (length(use_GG_effects) > 0) {
    cat("\nGreenhouse-Geisser correction will be applied for:\n")
    cat(paste("  -", use_GG_effects, collapse = "\n"))
    cat("\n")
  }
}

# Print ANCOVA results
cat("\n=== ANCOVA Results (Greenhouse-Geisser corrected) ===\n")
print(anova(anovaModel, correction = "GG", es = "pes"))

cat("\n=== Effect Sizes (Partial η²) ===\n")
etaSquared <- eta_squared(anovaModel, ci = 0.95)
print(etaSquared)

# Save results
capture.output(
  cat("=== Assumption Checks ===\n"),
  cat("\nNormality of accuracy:\n"),
  print(shapiroAccuracy),
  cat("\nMulticollinearity (VIF):\n"),
  print(vifValues),
  cat("\nHomogeneity of slopes:\n"),
  print(anovaSlopes),
  cat("\n\n=== Mauchly's Tests ===\n"),
  print(mauchly_summary_sph$sphericity.tests),
  cat("\n=== Formatted Mauchly's Results for Reporting ===\n"),
  if (!is.null(mauchly_summary_sph$sphericity.tests)) {
    for (i in 1:nrow(sphericity_sph)) {
      cat(sprintf("%s: W = %.3f, p = %.3f\n", 
                  rownames(sphericity_sph)[i],
                  sphericity_sph[i, "Test statistic"],
                  sphericity_sph[i, "p-value"]))
    }
  },
  cat("\n\n=== ANCOVA Results ===\n"),
  print(anova(anovaModel, correction = "GG", es = "pes")),
  cat("\n\n=== Effect Sizes ===\n"),
  print(etaSquared),
  file = file.path(outputDirectory, "ancova_results.txt")
)

# ============================================================================
# PART 5: POST-HOC ANALYSES
# ============================================================================

cat("\n=== PART 5: POST-HOC ANALYSES ===\n")

# 1. Pairwise comparisons for TMS intensity
cat("\n1. Pairwise comparisons for TMS intensity (Bonferroni-corrected):\n")
pairwiseIntensity <- emmeans(anovaModel, pairwise ~ intensity, adjust = "bonferroni")
print(pairwiseIntensity$contrasts)

write_csv(as.data.frame(pairwiseIntensity$contrasts), 
          file.path(outputDirectory, "pairwise_intensity.csv"))

# 2. Marginal means for TMS intensity
cat("\n4. Estimated marginal means for TMS intensity:\n")
emmIntensity <- emmeans(anovaModel, ~ intensity)
print(emmIntensity)

# ============================================================================
# PART 6A: FOLLOW-UP ANALYSIS AT HIGHEST INTENSITY
# ============================================================================

cat("\n=== PART 6A: FOLLOW-UP ANALYSIS AT HIGHEST INTENSITY ===\n")

dataHighest <- dataLong %>%
  filter(intensity == "130%")

cat(sprintf("Observations at 130%% intensity: %d\n", nrow(dataHighest)))

# Linear mixed model
modelHighest <- lmer(
  zMEP ~ accuracy_c * condition + (1 | Participant_ID),
  data = dataHighest,
  REML = FALSE
)

cat("\n=== Mixed Model Results (130% intensity) ===\n")
print(summary(modelHighest))
print(anova(modelHighest))

# Simple slopes at highest intensity
cat("\nSimple slopes at 130% intensity:\n")
slopesHighest <- emtrends(modelHighest, ~ condition, var = "accuracy_c")
print(slopesHighest)
print(confint(slopesHighest))

test(slopesHighest)

write_csv(as.data.frame(summary(slopesHighest)), 
          file.path(outputDirectory, "simple_slopes_highest_intensity.csv"))

# Plot for highest intensity
highestIntensityPlot <- ggplot(dataHighest, aes(x = accuracy, y = zMEP, colour = condition)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.5) +
  scale_colour_manual(values = c("AO" = "#F8766D", "AOMI" = "#00BFC4")) +
  labs(
    title = "Accuracy × Condition Interaction at 130% Intensity",
    x = "Accuracy",
    y = "Z-transformed MEP",
    colour = "Condition"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12)
  )

ggsave(file.path(outputDirectory, "interaction_highest_intensity.png"), 
       highestIntensityPlot, width = 8, height = 6, dpi = 300)


# ============================================================================
# PART 6B: ACCURACY × CONDITION INTERACTION ACROSS ALL INTENSITIES
# ============================================================================

cat("\n=== PART 6B: ACCURACY × CONDITION INTERACTION ACROSS ALL INTENSITIES ===\n")

# Initialize results storage
intensityLevels <- c("90%", "100%", "110%", "120%", "130%")
interactionResults <- list()

# List to store individual plots
intensityPlotsList <- list()

# Loop through each intensity level
for (intLevel in intensityLevels) {
  cat(sprintf("\n--- Analyzing %s RMT ---\n", intLevel))
  
  # Filter data for this intensity
  dataIntensity <- dataLong %>%
    filter(intensity == intLevel)
  
  n_trials <- nrow(dataIntensity)
  cat(sprintf("Number of observations: %d\n", n_trials))
  
  # Fit linear mixed model
  modelIntensity <- lmer(
    zMEP ~ accuracy_c * condition + (1 | Participant_ID),
    data = dataIntensity,
    REML = FALSE
  )
  
  # Get ANOVA results for interaction
  anovaIntensity <- anova(modelIntensity)
  
  # Extract interaction statistics
  interaction_F <- as.numeric(anovaIntensity["accuracy_c:condition", "F value"])
  interaction_df_num <- as.numeric(anovaIntensity["accuracy_c:condition", "NumDF"])
  interaction_df_denom <- as.numeric(anovaIntensity["accuracy_c:condition", "DenDF"])
  interaction_p <- as.numeric(anovaIntensity["accuracy_c:condition", "Pr(>F)"])
  
  # Calculate effect size (partial eta squared)
  SS_effect <- as.numeric(anovaIntensity["accuracy_c:condition", "Sum Sq"])
  SS_error <- sum(residuals(modelIntensity)^2)
  partial_eta_sq <- SS_effect / (SS_effect + SS_error)
  
  cat(sprintf("Interaction F(%d, %.1f) = %.2f, p = %.4f, ηp² = %.3f\n",
              interaction_df_num, interaction_df_denom, interaction_F, 
              interaction_p, partial_eta_sq))
  
  # Simple slopes analysis - use test() to get p-values directly
  slopesIntensity <- emtrends(modelIntensity, ~ condition, var = "accuracy_c")
  slopesSummary <- as.data.frame(summary(slopesIntensity))
  slopesTest <- as.data.frame(test(slopesIntensity))
  
  # Extract slope information - explicit numeric conversion
  AO_row <- slopesSummary[slopesSummary$condition == "AO", ]
  AOMI_row <- slopesSummary[slopesSummary$condition == "AOMI", ]
  
  AO_test <- slopesTest[slopesTest$condition == "AO", ]
  AOMI_test <- slopesTest[slopesTest$condition == "AOMI", ]
  
  # Check if rows exist
  if (nrow(AO_row) == 0 || nrow(AOMI_row) == 0) {
    cat("  WARNING: Could not extract slopes properly\n")
    next
  }
  
  # Extract AO values with explicit numeric conversion
  AO_slope <- as.numeric(AO_row$accuracy_c.trend)
  AO_slope_se <- as.numeric(AO_row$SE)
  AO_slope_t <- as.numeric(AO_test$t.ratio)
  AO_slope_df <- as.numeric(AO_test$df)
  AO_slope_p <- as.numeric(AO_test$p.value)
  
  # Extract AOMI values with explicit numeric conversion
  AOMI_slope <- as.numeric(AOMI_row$accuracy_c.trend)
  AOMI_slope_se <- as.numeric(AOMI_row$SE)
  AOMI_slope_t <- as.numeric(AOMI_test$t.ratio)
  AOMI_slope_df <- as.numeric(AOMI_test$df)
  AOMI_slope_p <- as.numeric(AOMI_test$p.value)
  
  # Slope difference
  slope_diff <- AO_slope - AOMI_slope
  
  cat(sprintf("  AO slope: β = %.2f, SE = %.2f, p = %.4f\n", 
              AO_slope, AO_slope_se, AO_slope_p))
  cat(sprintf("  AOMI slope: β = %.2f, SE = %.2f, p = %.4f\n", 
              AOMI_slope, AOMI_slope_se, AOMI_slope_p))
  cat(sprintf("  Slope difference: %.2f\n", slope_diff))
  
  # Store results in list
  interactionResults[[intLevel]] <- data.frame(
    intensity = intLevel,
    n_trials = n_trials,
    F_value = interaction_F,
    df_num = interaction_df_num,
    df_denom = interaction_df_denom,
    p_value = interaction_p,
    effect_size = partial_eta_sq,
    AO_slope = AO_slope,
    AO_slope_se = AO_slope_se,
    AO_slope_p = AO_slope_p,
    AOMI_slope = AOMI_slope,
    AOMI_slope_se = AOMI_slope_se,
    AOMI_slope_p = AOMI_slope_p,
    slope_difference = slope_diff,
    stringsAsFactors = FALSE
  )
  
  # Create plot for this intensity
  intensityPlot <- ggplot(dataIntensity, aes(x = accuracy, y = zMEP, colour = condition)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
    scale_colour_manual(values = c("AO" = "#F8766D", "AOMI" = "#00BFC4")) +
    labs(
      title = sprintf("%s RMT", intLevel),
      x = "Accuracy",
      y = "Z-transformed MEP",
      colour = "Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      legend.position = "right"
    )
  
  intensityPlotsList[[intLevel]] <- intensityPlot
}

# Combine results into a single dataframe
interactionResults <- bind_rows(interactionResults)

# Save results table
write_csv(interactionResults, 
          file.path(outputDirectory, "interaction_by_intensity.csv"))

cat("\n=== Summary Table: Accuracy × Condition Interaction by Intensity ===\n")
print(interactionResults)

# Create combined faceted plot with SAME y-axis for all panels
combinedIntensityPlot <- ggplot(dataLong, aes(x = accuracy, y = zMEP, colour = condition)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  facet_wrap(~ intensity, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("AO" = "#F8766D", "AOMI" = "#00BFC4")) +
  labs(
    title = "Accuracy × Condition Interaction Across TMS Intensities",
    x = "Accuracy",
    y = "Z-transformed MEP",
    colour = "Condition"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11)
  )

# Save combined plot
ggsave(file.path(outputDirectory, "interaction_all_intensities_faceted.png"), 
       combinedIntensityPlot, width = 16, height = 5, dpi = 300)

# Create individual plots arranged in grid (only if plots exist)
if (length(intensityPlotsList) > 0) {
  gridPlot <- wrap_plots(intensityPlotsList, ncol = 5) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  ggsave(file.path(outputDirectory, "interaction_all_intensities_grid.png"), 
         gridPlot, width = 18, height = 4, dpi = 300)
}

# Print combined plot
print(combinedIntensityPlot)

# Create summary visualisation of effect sizes and p-values
if (nrow(interactionResults) > 0) {
  # Reorder intensity as factor for proper plotting
  interactionResults$intensity <- factor(interactionResults$intensity, 
                                         levels = c("90%", "100%", "110%", "120%", "130%"))
  
  effectSizePlot <- ggplot(interactionResults, 
                           aes(x = intensity, y = effect_size)) +
    geom_col(fill = "#4575b4", alpha = 0.7) +
    geom_text(aes(label = sprintf("p = %.3f", p_value)), 
              vjust = -0.5, size = 3.5) +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_hline(yintercept = 0.06, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_hline(yintercept = 0.14, linetype = "dashed", color = "gray50", alpha = 0.5) +
    annotate("text", x = 0.7, y = 0.01, label = "Small", size = 3, hjust = 0) +
    annotate("text", x = 0.7, y = 0.06, label = "Medium", size = 3, hjust = 0) +
    annotate("text", x = 0.7, y = 0.14, label = "Large", size = 3, hjust = 0) +
    labs(
      title = "Effect Size (ηp²) of Accuracy × Condition Interaction by Intensity",
      x = "TMS Intensity (% RMT)",
      y = "Partial η²"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12)
    )
  
  ggsave(file.path(outputDirectory, "effect_sizes_by_intensity.png"), 
         effectSizePlot, width = 10, height = 6, dpi = 300)
  
  print(effectSizePlot)
  
  # Create slopes comparison plot
  slopesData <- interactionResults %>%
    select(intensity, AO_slope, AO_slope_se, AOMI_slope, AOMI_slope_se) %>%
    pivot_longer(
      cols = c(AO_slope, AOMI_slope),
      names_to = "condition_type",
      values_to = "slope"
    ) %>%
    mutate(
      condition = factor(ifelse(condition_type == "AO_slope", "AO", "AOMI"), 
                         levels = c("AO", "AOMI")),
      se = ifelse(condition_type == "AO_slope", AO_slope_se, AOMI_slope_se)
    ) %>%
    select(-condition_type)
  
  slopesPlot <- ggplot(slopesData, aes(x = intensity, y = slope, 
                                       color = condition, group = condition)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = slope - se, ymax = slope + se), 
                  width = 0.2, linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_color_manual(
      name = "Condition",
      values = c("AO" = "#F8766D", "AOMI" = "#00BFC4")
    ) +
    labs(
      title = "Simple Slopes: Accuracy Effect on zMEP by Condition and Intensity",
      x = "TMS Intensity (% RMT)",
      y = "Slope (β)",
      color = "Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12),
      legend.position = "bottom"
    )
  
  ggsave(file.path(outputDirectory, "slopes_by_intensity.png"), 
         slopesPlot, width = 10, height = 6, dpi = 300)
  
  print(slopesPlot)
}

cat("\n=== Intensity-Specific Analysis Complete ===\n")
cat("Generated files:\n")
cat("  - interaction_by_intensity.csv\n")
cat("  - interaction_all_intensities_faceted.png\n")
cat("  - interaction_all_intensities_grid.png\n")
cat("  - effect_sizes_by_intensity.png\n")
cat("  - slopes_by_intensity.png\n")

# ============================================================================
# PART 7: VISUALISATIONS
# ============================================================================

cat("\n=== PART 7: GENERATING VISUALISATIONS ===\n")

# Theme for publication
pubTheme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(size = 13, face = "plain"),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12)
  )

# 1. Main effect of condition
conditionPlot <- ggplot(dataLong, aes(x = condition, y = zMEP, fill = condition)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, colour = "black") +
  scale_fill_manual(values = c("AO" = "#F8766D", "AOMI" = "#00BFC4")) +
  labs(
    title = "Main Effect of Condition",
    x = "Condition",
    y = "Z-transformed MEP"
  ) +
  pubTheme +
  theme(legend.position = "none")

ggsave(file.path(outputDirectory, "main_effect_condition.png"), 
       conditionPlot, width = 6, height = 5, dpi = 300)

# 2. Main effect of TMS intensity
intensityColours <- brewer.pal(5, "Greens")

intensityPlot <- ggplot(dataLong, aes(x = intensity, y = zMEP, fill = intensity)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, colour = "black") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
               width = 0.2, colour = "black") +
  scale_fill_manual(
    name = "TMS Intensity",
    values = setNames(intensityColours, levels(dataLong$intensity))
  ) +
  labs(
    title = "Main Effect of TMS Intensity",
    x = "TMS Intensity (% RMT)",
    y = "Z-transformed MEP"
  ) +
  pubTheme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12))

ggsave(file.path(outputDirectory, "main_effect_intensity.png"), 
       intensityPlot, width = 8, height = 6, dpi = 300)

# 3. Accuracy × Condition interaction at highest intensity
interactionHighestPlot <- ggplot(dataHighest, aes(x = accuracy, y = zMEP, colour = condition)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.5) +
  scale_colour_manual(values = c("AO" = "#F8766D", "AOMI" = "#00BFC4")) +
  labs(
    title = "Accuracy × Condition Interaction at 130% Intensity",
    x = "Accuracy",
    y = "Z-transformed MEP",
    colour = "Condition"
  ) +
  pubTheme

ggsave(file.path(outputDirectory, "interaction_highest_intensity.png"), 
       interactionHighestPlot, width = 8, height = 6, dpi = 300)

# 4. Combined poster plot
posterPlot <- conditionPlot + intensityPlot + plot_layout(ncol = 2)

ggsave(file.path(outputDirectory, "combined_poster_plot.png"), 
       posterPlot, width = 14, height = 6, dpi = 300)

# 5. Marginal means for intensity
emmDF <- as.data.frame(emmIntensity)

marginalMeansPlot <- ggplot(emmDF, aes(x = intensity, y = emmean, group = 1)) +
  geom_line(colour = "blue", linewidth = 1.2) +
  geom_point(size = 3, colour = "blue") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, colour = "blue") +
  labs(
    title = "Estimated Marginal Means: TMS Intensity",
    x = "TMS Intensity (% RMT)",
    y = "Adjusted Z-transformed MEP"
  ) +
  pubTheme

ggsave(file.path(outputDirectory, "marginal_means_intensity.png"), 
       marginalMeansPlot, width = 8, height = 6, dpi = 300)

# 6. Combined three-panel figure for manuscript
# Load behavioural condition comparison data for Panel A
behaviouralConditionFile <- file.path(behaviouralDir, "accuracy_by_condition.xlsx")
if (!file.exists(behaviouralConditionFile)) {
  stop("Behavioural condition file not found. Run behavioural analysis script first.")
}
dfBehaviouralCondition <- read_excel(behaviouralConditionFile) %>%
  select(participant, AO, AOMI) %>%
  pivot_longer(cols = c(AO, AOMI), names_to = "condition", values_to = "accuracy") %>%
  mutate(condition = factor(condition, levels = c("AO", "AOMI")))

cat(sprintf("Loaded behavioural condition data for Panel A\n"))

# Panel A: Behavioural condition comparison (from behavioural script)
conditionPlotA <- ggplot(dfBehaviouralCondition, 
                         aes(x = condition, y = accuracy, group = participant)) +
  geom_line(colour = "grey", alpha = 0.7) +
  geom_point(aes(colour = condition), size = 3) +
  scale_colour_manual(
    name = "Condition",
    values = c("AO" = "#F8766D", "AOMI" = "#00BFC4")
  ) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               shape = 18, size = 4, colour = "black") +
  stat_summary(aes(group = 1), fun.data = mean_cl_normal, 
               geom = "errorbar", width = 0.1, colour = "black") +
  labs(
    title = "A",
    x = "Condition",
    y = "Proportion Correct"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 18, face = "bold", margin = margin(l = 0)),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )

# Panel B: TMS intensity plot (top right)
intensityPlotB <- ggplot(dataLong, aes(x = intensity, y = zMEP, fill = intensity)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, colour = "black") +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
               width = 0.2, colour = "black") +
  scale_fill_manual(
    name = "TMS Intensity",
    values = setNames(intensityColours, levels(dataLong$intensity))
  ) +
  labs(
    title = "B",
    x = "TMS Intensity (% RMT)",
    y = "Z-transformed MEP"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 18, face = "bold", margin = margin(l = 0)),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )

# Panel C: Accuracy × Condition interaction across intensities (bottom, full width)
interactionPlotC <- ggplot(dataLong, aes(x = accuracy, y = zMEP, colour = condition)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2) +
  facet_wrap(~ intensity, nrow = 1, scales = "free_x") +
  scale_colour_manual(
    name = "Condition",
    values = c("AO" = "#F8766D", "AOMI" = "#00BFC4")
  ) +
  labs(
    title = "C",
    x = "Proportion Correct",
    y = "Z-transformed MEP"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0, size = 18, face = "bold", margin = margin(l = 0)),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "top",
    legend.justification = "left",
    strip.text = element_text(face = "bold", size = 11)
  )

# Combine panels using patchwork
# Top row: A and B side by side
# Bottom row: C spanning full width (taller)
manuscriptFigure <- (conditionPlotA | intensityPlotB) / 
  wrap_elements(full = interactionPlotC) +
  plot_layout(heights = c(1, 1.3))

# Save the combined figure
ggsave(file.path(outputDirectory, "manuscript_figure_combined.png"), 
       manuscriptFigure, width = 16, height = 11, dpi = 300)

# Print key plots
print(conditionPlotA)
print(intensityPlotB)
print(interactionPlotC)
print(manuscriptFigure)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Results saved to: %s\n", outputDirectory))
cat("\nGenerated files:\n")
cat("  1. merged_data.csv\n")
cat("  2. ancova_results.txt\n")
cat("  3. pairwise_intensity.csv\n")
cat("  4. simple_slopes_highest_intensity.csv\n")
cat("  5. main_effect_condition.png\n")
cat("  6. main_effect_intensity.png\n")
cat("  7. interaction_highest_intensity.png\n")
cat("  8. combined_poster_plot.png\n")
cat("  9. marginal_means_intensity.png\n")
cat("  10. Assumption check plots (linearity, distribution)\n")

