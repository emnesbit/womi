# ============================================================================
# WOMI Behavioural Accuracy Analysis - BIDS Pipeline
# ============================================================================
# This script:
# 1. Reads BIDS behavioural data and calculates accuracy
# 2. Performs paired t-tests for session order and condition effects
# 3. Analyses accuracy by subjective weight responses
# 4. Generates publication-ready visualisations
# ============================================================================

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, haven, psych, writexl, afex, emmeans, 
  patchwork, effectsize, RColorBrewer
)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Set directories
bidsDirectory <- "/path/to/your/bids_directory/"
outputDirectory <- file.path(bidsDirectory, "derivatives", "behavioural_analysis")

# Create output directory
if (!dir.exists(outputDirectory)) {
  dir.create(outputDirectory, recursive = TRUE)
}

# Define participants and sessions
participants <- 1:18
sessions <- c("AO", "AOMI")

# Session order (odd participants start with AOMI, even with AO)
sessionOrder <- ifelse(participants %% 2 == 0, "AO", "AOMI")

# ============================================================================
# FUNCTIONS
# ============================================================================

# Extract behavioural data from BIDS format
extractBehaviouralData <- function(participant, session) {
  subID <- sprintf("sub-%02d", participant)
  
  # Construct BIDS behavioural file path
  behFile <- file.path(
    bidsDirectory, subID, paste0("ses-", session), "beh",
    sprintf("%s_ses-%s_womi_beh.tsv", subID, session)
  )
  
  if (!file.exists(behFile)) {
    warning(sprintf("Behavioural file not found: %s", behFile))
    return(NULL)
  }
  
  # Read BIDS behavioural TSV
  behData <- read_tsv(behFile, show_col_types = FALSE)
  uniqueTrials <- unique(behData$trial)
  
  # Initialise vectors
  weightShown <- rep(NA, length(uniqueTrials))
  weightResponse <- rep(NA, length(uniqueTrials))
  isControl <- rep(FALSE, length(uniqueTrials))
  
  # Process each trial
  for (i in seq_along(uniqueTrials)) {
    trialNum <- uniqueTrials[i]
    trialData <- behData %>% filter(trial == trialNum)
    
    # Check if control trial
    controlRow <- trialData %>% filter(event_type == "Control_TMS_pulse")
    isControl[i] <- nrow(controlRow) > 0
    
    if (!isControl[i]) {
      # Extract weight shown
      weightRow <- trialData %>% filter(event_type == "weight_in_video")
      if (nrow(weightRow) > 0) {
        weightShown[i] <- as.numeric(weightRow$code[1])
      }
      
      # Extract weight response
      responseRow <- trialData %>% filter(event_type == "weight_estimation_response")
      if (nrow(responseRow) > 0) {
        weightResponse[i] <- as.numeric(responseRow$code[1])
      }
    }
  }
  
  # Create trial dataframe
  trialDF <- data.frame(
    trial = uniqueTrials,
    weightShown = weightShown,
    weightResponse = weightResponse,
    isControl = isControl
  )
  
  # Filter to main experimental trials (exclude trial 1 and controls)
  trialDF <- trialDF %>%
    filter(!isControl, trial != 1, !is.na(weightShown), !is.na(weightResponse))
  
  if (nrow(trialDF) == 0) {
    warning(sprintf("No valid trials for %s session %s", subID, session))
    return(NULL)
  }
  
  # Calculate accuracy
  trialDF <- trialDF %>%
    mutate(correct = weightShown == weightResponse)
  
  # Overall accuracy
  overallAccuracy <- mean(trialDF$correct)
  totalTrials <- nrow(trialDF)
  correctTrials <- sum(trialDF$correct)
  
  # Accuracy by weight response
  responseAccuracy <- trialDF %>%
    group_by(weightResponse) %>%
    summarise(accuracy = mean(correct), .groups = "drop")
  
  return(list(
    accuracy = overallAccuracy,
    totalTrials = totalTrials,
    correctTrials = correctTrials,
    responseAccuracy = responseAccuracy,
    trialData = trialDF
  ))
}

# ============================================================================
# PART 1: DATA EXTRACTION
# ============================================================================

cat("\n=== PART 1: EXTRACTING BEHAVIOURAL DATA ===\n")

allResultsList <- list()
allTrialData <- list()

for (p in participants) {
  cat(sprintf("\nProcessing participant %02d...\n", p))
  
  # Process both sessions
  aoData <- extractBehaviouralData(p, "AO")
  aomiData <- extractBehaviouralData(p, "AOMI")
  
  # Create participant row
  participantRow <- data.frame(participant = p)
  
  # Add AO data
  if (!is.null(aoData)) {
    participantRow$AO_accuracy <- aoData$accuracy
    participantRow$AO_total_trials <- aoData$totalTrials
    participantRow$AO_correct_trials <- aoData$correctTrials
    
    for (resp in 1:3) {
      respAcc <- aoData$responseAccuracy %>% 
        filter(weightResponse == resp) %>% 
        pull(accuracy)
      participantRow[[paste0("AO_response", resp)]] <- 
        ifelse(length(respAcc) > 0, respAcc, NA)
    }
    
    # Store trial-level data
    allTrialData[[paste0("P", p, "_AO")]] <- aoData$trialData %>%
      mutate(participant = p, session = "AO")
  } else {
    participantRow[c("AO_accuracy", "AO_total_trials", "AO_correct_trials",
                     "AO_response1", "AO_response2", "AO_response3")] <- NA
  }
  
  # Add AOMI data
  if (!is.null(aomiData)) {
    participantRow$AOMI_accuracy <- aomiData$accuracy
    participantRow$AOMI_total_trials <- aomiData$totalTrials
    participantRow$AOMI_correct_trials <- aomiData$correctTrials
    
    for (resp in 1:3) {
      respAcc <- aomiData$responseAccuracy %>% 
        filter(weightResponse == resp) %>% 
        pull(accuracy)
      participantRow[[paste0("AOMI_response", resp)]] <- 
        ifelse(length(respAcc) > 0, respAcc, NA)
    }
    
    # Store trial-level data
    allTrialData[[paste0("P", p, "_AOMI")]] <- aomiData$trialData %>%
      mutate(participant = p, session = "AOMI")
  } else {
    participantRow[c("AOMI_accuracy", "AOMI_total_trials", "AOMI_correct_trials",
                     "AOMI_response1", "AOMI_response2", "AOMI_response3")] <- NA
  }
  
  allResultsList[[p]] <- participantRow
}

# Combine results
allResults <- bind_rows(allResultsList)
allTrials <- bind_rows(allTrialData)

cat("\n=== Data extraction complete ===\n")
cat(sprintf("Total participants: %d\n", nrow(allResults)))
cat(sprintf("Total trials: %d\n", nrow(allTrials)))

# ============================================================================
# PART 2: SUMMARY TABLES
# ============================================================================

cat("\n=== PART 2: CREATING SUMMARY TABLES ===\n")

# Accuracy by condition
accuracyByCondition <- allResults %>%
  select(participant, starts_with("AO"), starts_with("AOMI")) %>%
  select(participant, AO_accuracy, AOMI_accuracy, 
         AO_total_trials, AOMI_total_trials,
         AO_correct_trials, AOMI_correct_trials) %>%
  rename(AO = AO_accuracy, AOMI = AOMI_accuracy)

# Session order accuracy
accuracyBySessionOrder <- allResults %>%
  mutate(
    firstSession = sessionOrder[participant],
    session1Acc = ifelse(firstSession == "AO", AO_accuracy, AOMI_accuracy),
    session2Acc = ifelse(firstSession == "AO", AOMI_accuracy, AO_accuracy)
  ) %>%
  select(participant, session1Acc, session2Acc) %>%
  rename(Session1 = session1Acc, Session2 = session2Acc)

# Accuracy by weight response
accuracyByWeight <- allResults %>%
  select(participant, starts_with("AO_response"), starts_with("AOMI_response"))

# Save tables
write_xlsx(accuracyByCondition, 
           file.path(outputDirectory, "accuracy_by_condition.xlsx"))
write_xlsx(accuracyBySessionOrder, 
           file.path(outputDirectory, "accuracy_by_session_order.xlsx"))
write_xlsx(accuracyByWeight, 
           file.path(outputDirectory, "accuracy_by_weight.xlsx"))

# Print summary statistics
cat("\nMean accuracies:\n")
cat(sprintf("  AO:   %.3f (SD: %.3f)\n", 
            mean(accuracyByCondition$AO, na.rm = TRUE),
            sd(accuracyByCondition$AO, na.rm = TRUE)))
cat(sprintf("  AOMI: %.3f (SD: %.3f)\n", 
            mean(accuracyByCondition$AOMI, na.rm = TRUE),
            sd(accuracyByCondition$AOMI, na.rm = TRUE)))

# ============================================================================
# PART 3: SESSION ORDER ANALYSIS
# ============================================================================

cat("\n=== PART 3: SESSION ORDER ANALYSIS (S1 vs S2) ===\n")

# Descriptive statistics
descriptivesOrder <- accuracyBySessionOrder %>%
  summarise(
    mean_S1 = mean(Session1, na.rm = TRUE),
    sd_S1 = sd(Session1, na.rm = TRUE),
    mean_S2 = mean(Session2, na.rm = TRUE),
    sd_S2 = sd(Session2, na.rm = TRUE)
  )
cat("\nDescriptive statistics:\n")
print(descriptivesOrder)

# Check normality
diffOrder <- accuracyBySessionOrder$Session1 - accuracyBySessionOrder$Session2
shapiroOrder <- shapiro.test(diffOrder)
cat("\nShapiro-Wilk test (S1 - S2):\n")
print(shapiroOrder)

# Paired t-test
tTestOrder <- t.test(accuracyBySessionOrder$Session1, 
                     accuracyBySessionOrder$Session2, 
                     paired = TRUE, conf.level = 0.95)
cat("\nPaired t-test:\n")
print(tTestOrder)

# Effect size
cohenDOrder <- cohens_d(accuracyBySessionOrder$Session1, 
                        accuracyBySessionOrder$Session2, 
                        paired = TRUE)
cat("\nEffect size (Cohen's d):\n")
print(cohenDOrder)

# Visualisation
orderData <- accuracyBySessionOrder %>%
  mutate(
    firstTask = sessionOrder[participant],
    secondTask = ifelse(firstTask == "AO", "AOMI", "AO")
  ) %>%
  pivot_longer(
    cols = c(Session1, Session2),
    names_to = "session",
    values_to = "accuracy"
  ) %>%
  mutate(
    condition = factor(
      ifelse(session == "Session1", firstTask, secondTask),
      levels = c("AO", "AOMI")
    )
  )

orderPlot <- ggplot(orderData, aes(x = session, y = accuracy, group = participant)) +
  geom_line(colour = "grey70", alpha = 0.7) +
  geom_point(aes(colour = condition), size = 3) +
  scale_colour_manual(
    name = "Condition",
    values = c(
      "AO"   = "#F8766D",  # red
      "AOMI" = "#00BFC4"   # blue
    )) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               shape = 18, size = 4, colour = "black") +
  stat_summary(aes(group = 1), fun.data = mean_cl_normal, 
               geom = "errorbar", width = 0.1, colour = "black") +
  labs(
    title = "Accuracy and Session Order Effects",
    x = "Session",
    y = "Proportion Correct"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )

ggsave(file.path(outputDirectory, "session_order_plot.png"), 
       orderPlot, width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 4: CONDITION ANALYSIS
# ============================================================================

cat("\n=== PART 4: CONDITION ANALYSIS (AO vs AOMI) ===\n")

# Descriptive statistics
descriptivesCondition <- accuracyByCondition %>%
  summarise(
    mean_AO = mean(AO, na.rm = TRUE),
    sd_AO = sd(AO, na.rm = TRUE),
    mean_AOMI = mean(AOMI, na.rm = TRUE),
    sd_AOMI = sd(AOMI, na.rm = TRUE)
  )
cat("\nDescriptive statistics:\n")
print(descriptivesCondition)

# Check normality
diffCondition <- accuracyByCondition$AO - accuracyByCondition$AOMI
shapiroCondition <- shapiro.test(diffCondition)
cat("\nShapiro-Wilk test (AO - AOMI):\n")
print(shapiroCondition)

# Paired t-test
tTestCondition <- t.test(accuracyByCondition$AO, 
                         accuracyByCondition$AOMI, 
                         paired = TRUE, conf.level = 0.95)
cat("\nPaired t-test:\n")
print(tTestCondition)

# Effect size
cohenDCondition <- cohens_d(accuracyByCondition$AO, 
                            accuracyByCondition$AOMI, 
                            paired = TRUE)
cat("\nEffect size (Cohen's d):\n")
print(cohenDCondition)

# Visualisation - standard version
conditionData <- accuracyByCondition %>%
  select(participant, AO, AOMI) %>%
  pivot_longer(cols = c(AO, AOMI), names_to = "condition", values_to = "accuracy")

conditionPlot <- ggplot(conditionData, 
                        aes(x = condition, y = accuracy, group = participant)) +
  geom_line(colour = "grey", alpha = 0.7) +
  geom_point(aes(colour = condition), size = 3) +
  scale_colour_manual(values = c("AO" = "#F8766D", "AOMI" = "#00BFC4")) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               shape = 18, size = 4, colour = "black") +
  stat_summary(aes(group = 1), fun.data = mean_cl_normal, 
               geom = "errorbar", width = 0.1, colour = "black") +
  labs(
    title = "Paired Comparison of AO and AOMI",
    x = "Condition",
    y = "Proportion Correct"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )

ggsave(file.path(outputDirectory, "condition_comparison_plot.png"), 
       conditionPlot, width = 8, height = 6, dpi = 300)

# Visualisation - panel style version (for manuscript integration)
conditionPlotForPanel <- ggplot(conditionData, 
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
    x = "Condition",
    y = "Proportion Correct"
  ) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )

ggsave(file.path(outputDirectory, "condition_comparison_for_panel.png"), 
       conditionPlotForPanel, width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 5: WEIGHT RESPONSE ANALYSIS
# ============================================================================

cat("\n=== PART 5: WEIGHT RESPONSE ANALYSIS ===\n")

# Prepare data for ANOVA
weightData <- accuracyByWeight %>%
  pivot_longer(
    cols = -participant,
    names_to = "variable",
    values_to = "accuracy"
  ) %>%
  mutate(
    condition = if_else(str_detect(variable, "^AO_"), "AO", "AOMI"),
    weight = str_extract(variable, "\\d+$"),
    condition = factor(condition, levels = c("AO", "AOMI")),
    weight = factor(weight, levels = c("1", "2", "3"), labels = c("Light", "Moderate", "Heavy")),
    participantID = factor(participant)
  ) %>%
  select(participantID, condition, weight, accuracy) %>%
  filter(!is.na(accuracy))

# Descriptive statistics
descriptivesWeight <- weightData %>%
  group_by(condition, weight) %>%
  summarise(
    mean = mean(accuracy, na.rm = TRUE),
    sd = sd(accuracy, na.rm = TRUE),
    n = sum(!is.na(accuracy)),
    .groups = "drop"
  )
cat("\nDescriptive statistics by condition and weight:\n")
print(descriptivesWeight)

# Two-way repeated measures ANOVA
cat("\n=== Two-way repeated measures ANOVA ===\n")

anova2Way <- aov_car(
  accuracy ~ condition * weight + Error(participantID / (condition * weight)),
  data = weightData,
  type = 3,
  factorize = FALSE
)

# Mauchly's test
weightDataWide <- weightData %>%
  unite("cell", condition, weight, sep = "_") %>%
  pivot_wider(names_from = cell, values_from = accuracy, id_cols = participantID)

mlm_model <- lm(cbind(AO_Light, AO_Moderate, AO_Heavy, 
                      AOMI_Light, AOMI_Moderate, AOMI_Heavy) ~ 1,
                data = weightDataWide)

idata <- data.frame(
  condition = factor(rep(c("AO", "AOMI"), each = 3)),
  weight = factor(rep(c("Light", "Moderate", "Heavy"), 2))
)

model_anova <- car::Anova(mlm_model, idata = idata, 
                          idesign = ~condition * weight, type = 3)
mauchly_summary <- summary(model_anova, multivariate = FALSE)

cat("\nMauchly's Test of Sphericity:\n")
cat("Condition: Not applicable (only 2 levels)\n")

use_GG_weight <- FALSE
use_GG_interaction <- FALSE

if (!is.null(mauchly_summary$sphericity.tests)) {
  sphericity <- mauchly_summary$sphericity.tests
  
  if ("weight" %in% rownames(sphericity)) {
    cat(sprintf("Weight: Mauchly's W = %.3f, p = %.3f\n", 
                sphericity["weight", "Test statistic"],
                sphericity["weight", "p-value"]))
    use_GG_weight <- sphericity["weight", "p-value"] < 0.05
  }
  
  if ("condition:weight" %in% rownames(sphericity)) {
    cat(sprintf("Condition × Weight interaction: Mauchly's W = %.3f, p = %.3f\n", 
                sphericity["condition:weight", "Test statistic"],
                sphericity["condition:weight", "p-value"]))
    use_GG_interaction <- sphericity["condition:weight", "p-value"] < 0.05
  }
}

# ANOVA results
cat("\n2-Way ANOVA results:\n")
if (use_GG_interaction) {
  print(anova(anova2Way, correction = "GG", es = "pes"))
  cat("\n(Note: Greenhouse-Geisser correction applied for interaction)\n")
} else {
  print(anova(anova2Way, correction = "none", es = "pes"))
}

cat("\nEffect sizes (partial eta squared):\n")
print(eta_squared(anova2Way, ci = 0.95))

# Overall weight response counts
cat("\n=== Overall weight response counts ===\n")

weightResponseCounts <- allTrials %>%
  count(weightResponse) %>%
  arrange(weightResponse)

cat("\nWeight response distribution:\n")
print(weightResponseCounts)

observed <- weightResponseCounts$n
totalResponses <- sum(observed)

chiSquareTest <- chisq.test(observed, p = c(1/3, 1/3, 1/3))
cat(sprintf("\nχ²(%d) = %.3f, p < .001\n", 
            chiSquareTest$parameter, 
            chiSquareTest$statistic))
cat(sprintf("Total weight responses: %d\n", totalResponses))

# One-way ANOVA (weight only)
cat("\n=== One-way repeated measures ANOVA (weight only) ===\n")

anova1Way <- aov_car(
  accuracy ~ weight + Error(participantID / weight),
  data = weightData,
  type = 3,
  factorize = FALSE,
  fun_aggregate = mean
)

# Mauchly's test for one-way
weightData1Way <- weightData %>%
  group_by(participantID, weight) %>%
  summarise(accuracy = mean(accuracy, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = weight, values_from = accuracy)

mlm_model_1way <- lm(cbind(Light, Moderate, Heavy) ~ 1,
                     data = weightData1Way)

idata_1way <- data.frame(
  weight = factor(c("Light", "Moderate", "Heavy"))
)

model_anova_1way <- car::Anova(mlm_model_1way, idata = idata_1way, 
                               idesign = ~weight, type = 3)
mauchly_summary_1way <- summary(model_anova_1way, multivariate = FALSE)

cat("\nMauchly's Test of Sphericity:\n")

use_GG_1way <- FALSE

if (!is.null(mauchly_summary_1way$sphericity.tests)) {
  sphericity_1way <- mauchly_summary_1way$sphericity.tests
  
  if ("weight" %in% rownames(sphericity_1way)) {
    cat(sprintf("Weight: Mauchly's W = %.3f, p = %.3f\n", 
                sphericity_1way["weight", "Test statistic"],
                sphericity_1way["weight", "p-value"]))
    use_GG_1way <- sphericity_1way["weight", "p-value"] < 0.05
  }
}

cat("\nOne-way ANOVA results:\n")
if (use_GG_1way) {
  print(anova(anova1Way, correction = "GG", es = "pes"))
} else {
  print(anova(anova1Way, correction = "none", es = "pes"))
}

cat("\nEffect size (partial eta squared):\n")
print(eta_squared(anova1Way, ci = 0.95))

# Post-hoc comparisons
cat("\n=== Post-hoc comparisons (Bonferroni corrected) ===\n")
posthoc <- emmeans(anova1Way, pairwise ~ weight, adjust = "bonferroni")
print(posthoc$contrasts)

# Visualisations
weightColours <- c(
  "Light"    = "#FDD0A2",
  "Moderate" = "#FDAE6B",
  "Heavy"    = "#E6550D"
)

weightByConditionPlot <- ggplot(weightData, 
                                aes(x = weight, y = accuracy, group = participantID)) +
  geom_line(colour = "grey", alpha = 0.7) +
  geom_point(aes(colour = weight), size = 3) +
  scale_colour_manual(values = weightColours) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               shape = 18, size = 4, colour = "black") +
  stat_summary(aes(group = 1), fun.data = mean_cl_normal, 
               geom = "errorbar", width = 0.1, colour = "black") +
  facet_wrap(~ condition) +
  labs(
    title = "Accuracy by Condition and Weight Response",
    x = "Weight Estimation",
    y = "Proportion Correct"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "none"
  )

ggsave(file.path(outputDirectory, "weight_by_condition_plot.png"), 
       weightByConditionPlot, width = 10, height = 5, dpi = 300)

weightPlot <- ggplot(weightData,
                     aes(x = weight, y = accuracy,
                         group = interaction(participantID, condition))) +
  geom_line(colour = "grey80", alpha = 0.6) +
  geom_point(aes(colour = weight), size = 3) +
  scale_colour_manual(
    name = "Weight Response",
    values = weightColours,
    labels = c("Light", "Moderate", "Heavy")
  ) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               shape = 18, size = 4, colour = "black") +
  stat_summary(aes(group = 1), fun.data = mean_cl_normal, 
               geom = "errorbar", width = 0.1, colour = "black") +
  labs(
    title = "Accuracy by Weight Response",
    x = "Weight Response",
    y = "Proportion Correct"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12),
    legend.position = "right"
  )

ggsave(file.path(outputDirectory, "weight_responses_plot.png"), 
       weightPlot, width = 8, height = 6, dpi = 300)

# ============================================================================
# PART 6: COMBINED MANUSCRIPT FIGURE
# ============================================================================

cat("\n=== PART 6: CREATING COMBINED MANUSCRIPT FIGURE ===\n")

# Panel A: Session order plot (top left)
orderPlotA <- ggplot(orderData, aes(x = session, y = accuracy, group = participant)) +
  geom_line(colour = "grey70", alpha = 0.7) +
  geom_point(aes(colour = condition), size = 3) +
  scale_colour_manual(
    name = "Condition",
    values = c(
      "AO"   = "#F8766D",
      "AOMI" = "#00BFC4"
    )) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               shape = 18, size = 4, colour = "black") +
  stat_summary(aes(group = 1), fun.data = mean_cl_normal, 
               geom = "errorbar", width = 0.1, colour = "black") +
  labs(
    title = "A",
    x = "Session",
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

# Panel B: Weight response plot (top right)
weightPlotB <- ggplot(weightData,
                      aes(x = weight, y = accuracy,
                          group = interaction(participantID, condition))) +
  geom_line(colour = "grey80", alpha = 0.6) +
  geom_point(aes(colour = weight), size = 3) +
  scale_colour_manual(
    name = "Weight Response",
    values = weightColours,
    labels = c("Light", "Moderate", "Heavy")
  ) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", 
               shape = 18, size = 4, colour = "black") +
  stat_summary(aes(group = 1), fun.data = mean_cl_normal, 
               geom = "errorbar", width = 0.1, colour = "black") +
  labs(
    title = "B",
    x = "Weight Response",
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

# Combine panels using patchwork
behaviouralManuscriptFigure <- orderPlotA | weightPlotB

# Save the combined figure
ggsave(file.path(outputDirectory, "behavioural_manuscript_figure.png"), 
       behaviouralManuscriptFigure, width = 14, height = 6, dpi = 300)

cat("\n=== Behavioural Manuscript Figure Created ===\n")
cat("Two-panel figure saved as: behavioural_manuscript_figure.png\n")
cat("  - Panel A (left): Session order effects\n")
cat("  - Panel B (right): Accuracy by weight response\n")

# Print the combined figure
print(behaviouralManuscriptFigure)

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Results saved to: %s\n", outputDirectory))
cat("\nGenerated files:\n")
cat("  1. accuracy_by_condition.xlsx\n")
cat("  2. accuracy_by_session_order.xlsx\n")
cat("  3. accuracy_by_weight.xlsx\n")
cat("  4. session_order_plot.png\n")
cat("  5. condition_comparison_plot.png\n")
cat("  6. condition_comparison_for_panel.png (ready for manuscript panel A)\n")
cat("  7. weight_by_condition_plot.png\n")
cat("  8. weight_responses_plot.png\n")
cat("  9. behavioural_manuscript_figure.png\n")

