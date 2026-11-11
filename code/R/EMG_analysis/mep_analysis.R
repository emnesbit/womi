# ============================================================================
# WOMI - MEP analyses
# Follows the Results: 2 (Fig. 3), 3 (ST1-ST3, Fig. 4), 4 (ST4, Figs. S3-S5)
# Needs accuracy_by_condition.xlsx from the accuracy script.
# ============================================================================

pacman::p_load(tidyverse, readxl, afex, emmeans, effectsize, patchwork, broom)

bidsDir  <- "/data/path/to/bids_directory/"
derivDir <- file.path(bidsDir, "derivatives")
outDir   <- file.path(derivDir, "ancova_analysis")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

condColours       <- c(AO = "#F8766D", AOMI = "#00BFC4")
aoWeightColours   <- c(Light = "#FBC4BF", Moderate = "#F8766D", Heavy = "#C2412F")
aomiWeightColours <- c(Light = "#A6E6E8", Moderate = "#00BFC4", Heavy = "#007F83")
intensityLabels   <- c("90%", "100%", "110%", "120%", "130%")

# RM ANOVA table, Greenhouse-Geisser corrected only where Mauchly p < .05
rmTable <- function(m) {
  tab <- anova(m, correction = "none", es = "pes")
  s   <- summary(m$Anova, multivariate = FALSE)
  if (length(s$sphericity.tests) > 0) {
    sph <- s$sphericity.tests
    print(cbind(sph, GG_eps = s$pval.adjustments[rownames(sph), "GG eps"]))
    bad <- rownames(sph)[sph[, "p-value"] < .05]
    tab[bad, ] <- anova(m, correction = "GG", es = "pes")[bad, ]
  }
  tab
}

# Format as in the supplement: F(df1, df2) = F, p = p, ηp² = pes
fmtTable <- function(tab) {
  df    <- function(x) ifelse(abs(x - round(x)) < 1e-6, sprintf("%d", round(x)), sprintf("%.2f", x))
  small <- function(x) ifelse(x < .005, "< 0.01", sprintf("= %.2f", x))
  tibble(effect = rownames(tab),
         result = sprintf("F(%s, %s) %s, p %s, ηp² %s",
                          df(tab$`num Df`), df(tab$`den Df`), small(tab$F),
                          ifelse(tab$`Pr(>F)` < .001, "< 0.001", sprintf("= %.3f", tab$`Pr(>F)`)),
                          small(tab$pes)))
}

# Simple slope of accuracy on y, per condition (one value per participant)
simpleSlopes <- function(d, y, x = "accuracy") {
  d %>%
    group_by(condition) %>%
    group_modify(~ tidy(lm(reformulate(x, y), data = .x)) %>% filter(term == x)) %>%
    print()
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
accFile <- file.path(derivDir, "behavioural_analysis", "accuracy_by_condition.xlsx")
acc <- read_excel(accFile) %>% transmute(participant, accuracy = (AO + AOMI) / 2)

orderInfo <- tibble(participant = 1:18,
                    order = factor(if_else(participant %% 2 == 0, "AO first", "AOMI first")))

dataLong <- read_csv(file.path(derivDir, "mep_ztransform", "mean_z_scores_by_intensity_and_weight.csv"),
                     show_col_types = FALSE) %>%
  mutate(participant = as.numeric(str_extract(Participant, "\\d+"))) %>%
  pivot_longer(matches("^(AO|AOMI)_WtResp_\\d+_\\d+$"),
               names_to = c("condition", "weight", "intensity"),
               names_pattern = "(AO|AOMI)_WtResp_(\\d+)_(\\d+)",
               values_to = "zMEP") %>%
  filter(!is.na(zMEP)) %>%
  left_join(acc, by = "participant") %>%
  left_join(orderInfo, by = "participant") %>%
  mutate(id         = factor(Participant),
         condition  = factor(condition, levels = c("AO", "AOMI")),
         weight     = factor(weight, levels = 1:3, labels = c("Light", "Moderate", "Heavy")),
         intensity  = factor(intensity, levels = 1:5, labels = intensityLabels),
         accuracy_c = accuracy - mean(acc$accuracy))

cat(sprintf("Observations: %d (expected %d)\n", nrow(dataLong), 18 * 2 * 3 * 5))

# Recruitment-curve slope per participant x condition x weight (z per 10% RMT)
curves <- dataLong %>%
  mutate(intensity_num = as.numeric(sub("%", "", intensity))) %>%
  group_by(id, participant, accuracy, accuracy_c, order, condition, weight) %>%
  filter(n() == 5) %>%
  summarise(slope = coef(lm(zMEP ~ intensity_num))[2] * 10, .groups = "drop")

# Slope averaged over weights: one value per participant x condition
gain <- curves %>%
  group_by(id, accuracy, order, condition) %>%
  summarise(slope = mean(slope), .groups = "drop")

# ============================================================================
# 2. EFFECTS OF CONDITION AND WEIGHT ON CSE (Hypotheses 1 and 2)
# ============================================================================

# Manipulation check: zMEP by intensity
intensityMeans <- dataLong %>%
  group_by(id, intensity) %>%
  summarise(zMEP = mean(zMEP), .groups = "drop")
print(intensityMeans %>% group_by(intensity) %>% summarise(M = mean(zMEP), SD = sd(zMEP)))

intensityModel <- aov_car(zMEP ~ intensity + Error(id / intensity), data = intensityMeans)
print(rmTable(intensityModel))
print(emmeans(intensityModel, consec ~ intensity, adjust = "bonferroni")$contrasts)

# Slope ANOVA: condition x weight
print(curves %>% group_by(condition, weight) %>%
        summarise(M = mean(slope), SD = sd(slope), .groups = "drop"))
slopeModel <- aov_car(slope ~ condition * weight + Error(id / (condition * weight)), data = curves)
print(rmTable(slopeModel))

# AOMI vs AO on the slope
gainWide <- gain %>% select(id, condition, slope) %>%
  pivot_wider(names_from = condition, values_from = slope)
cat(sprintf("Steeper slope in AOMI than AO: %d of %d participants\n",
            sum(gainWide$AOMI > gainWide$AO), nrow(gainWide)))
print(cohens_d(gainWide$AOMI, gainWide$AO, paired = TRUE))

# Figure 3: zMEP across intensities by weight response, AO and AOMI
plotCurves <- function(data, colours) {
  d <- data %>%
    group_by(weight, id, intensity) %>%
    summarise(zMEP = mean(zMEP), .groups = "drop")
  ggplot(d, aes(intensity, zMEP)) +
    geom_line(aes(group = id), colour = "grey88", linewidth = 0.4) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", colour = "black", linewidth = 2) +
    stat_summary(aes(group = 1, colour = weight), fun = mean, geom = "line", linewidth = 1.3) +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "errorbar", width = 0.15, linewidth = 0.5) +
    stat_summary(aes(fill = weight), fun = mean, geom = "point", shape = 21, size = 3.5, stroke = 0.6) +
    facet_wrap(~ weight, nrow = 1) +
    scale_colour_manual(values = colours) +
    scale_fill_manual(values = colours) +
    labs(x = "TMS intensity (% RMT)", y = "Z-transformed MEP") +
    theme_minimal() +
    theme(legend.position = "none",
          strip.text = element_text(face = "bold", size = 12),
          axis.title = element_text(size = 13), axis.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

fig3 <- (plotCurves(filter(dataLong, condition == "AO"), aoWeightColours) + labs(title = "AO") |
           plotCurves(filter(dataLong, condition == "AOMI"), aomiWeightColours) + labs(title = "AOMI", y = NULL)) +
  plot_annotation(tag_levels = "A") &
  coord_cartesian(ylim = range(dataLong$zMEP)) &
  theme(plot.tag = element_text(size = 18, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
print(fig3)
ggsave(file.path(outDir, "Fig3_zMEP_by_intensity.png"), fig3, width = 18, height = 4.5, dpi = 300)

# ============================================================================
# 3. EXPLORATORY: WEIGHT ESTIMATION ACCURACY
# ============================================================================

print(shapiro.test(acc$accuracy))

# Table ST1: slope ANCOVA
slopeAccModel <- aov_car(slope ~ accuracy_c + condition * weight + Error(id / (condition * weight)),
                         data = curves, factorize = FALSE)
print(fmtTable(rmTable(slopeAccModel)), n = Inf)

# Simple slopes of accuracy on the slope, per condition (Fig. 4A)
simpleSlopes(gain, "slope")

# Table ST2: zMEP ANCOVA
zModel <- aov_car(zMEP ~ accuracy_c + condition * weight * intensity +
                    Error(id / (condition * weight * intensity)),
                  data = dataLong, factorize = FALSE)
print(fmtTable(rmTable(zModel)), n = Inf)
print(eta_squared(zModel, ci = .95) %>% filter(Parameter == "accuracy_c:condition:intensity"))

# Table ST3: accuracy x condition at each intensity, Bonferroni over 5 intensities
slopeFor <- function(d, cond) {
  d <- d %>% filter(condition == cond) %>%
    group_by(id, accuracy_c) %>% summarise(z = mean(zMEP), .groups = "drop")
  s <- summary(lm(z ~ accuracy_c, data = d))$coefficients["accuracy_c", ]
  c(beta = s[["Estimate"]], se = s[["Std. Error"]], p = s[["Pr(>|t|)"]], df = nrow(d) - 2)
}

st3 <- map_dfr(intensityLabels, function(lvl) {
  d   <- filter(dataLong, intensity == lvl)
  m   <- aov_car(zMEP ~ accuracy_c + condition * weight + Error(id / (condition * weight)),
                 data = d, factorize = FALSE)
  tab <- anova(m, correction = "none", es = "pes")["accuracy_c:condition", ]
  ao  <- slopeFor(d, "AO"); aomi <- slopeFor(d, "AOMI")
  tibble(intensity = lvl, F = tab$F, df1 = tab$`num Df`, df2 = tab$`den Df`,
         p = tab$`Pr(>F)`, pes = tab$pes,
         AO_beta = ao[["beta"]], AO_se = ao[["se"]], AO_p = ao[["p"]],
         AOMI_beta = aomi[["beta"]], AOMI_se = aomi[["se"]], AOMI_p = aomi[["p"]],
         slope_df = ao[["df"]])
}) %>%
  mutate(intensity   = factor(intensity, levels = intensityLabels),
         p_bonf      = p.adjust(p, "bonferroni"),
         AO_p_bonf   = p.adjust(AO_p, "bonferroni"),
         AOMI_p_bonf = p.adjust(AOMI_p, "bonferroni"))
print(st3 %>% mutate(across(where(is.numeric), ~ round(.x, 3))), width = Inf)

# Figure 4: A = slope vs accuracy; B = simple slopes per intensity (Bonferroni 95% CI)
fig4A <- ggplot(gain, aes(accuracy, slope, colour = condition, fill = condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(aes(group = id), colour = "grey80", linewidth = 0.4) +
  geom_smooth(method = "lm", formula = y ~ x, alpha = 0.15, linewidth = 1.2) +
  geom_point(shape = 21, size = 3, colour = "black", stroke = 0.5) +
  scale_colour_manual(name = "Condition", values = condColours) +
  scale_fill_manual(name = "Condition", values = condColours) +
  labs(x = "Weight estimation accuracy (proportion correct)",
       y = "Recruitment curve gain (z per 10% RMT)") +
  theme_minimal()

tCrit <- qt(1 - .05 / (2 * nrow(st3)), df = st3$slope_df[1])
betas <- bind_rows(
  st3 %>% transmute(intensity, condition = "AO",   beta = AO_beta,   se = AO_se),
  st3 %>% transmute(intensity, condition = "AOMI", beta = AOMI_beta, se = AOMI_se)) %>%
  mutate(condition = factor(condition, levels = c("AO", "AOMI")),
         lower = beta - tCrit * se, upper = beta + tCrit * se)
starData <- st3 %>% filter(p_bonf < .05) %>%
  mutate(y = max(betas$upper) + 0.06 * diff(range(c(betas$lower, betas$upper))))
pd <- position_dodge(width = 0.35)

fig4B <- ggplot(betas, aes(intensity, beta, colour = condition, group = condition)) +
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf, fill = "grey92") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth = 0.7, position = pd) +
  geom_line(linewidth = 1.1, position = pd) +
  geom_point(aes(fill = condition), shape = 21, size = 3.5, colour = "black", stroke = 0.6, position = pd) +
  geom_text(data = starData, aes(x = intensity, y = y), label = "*", inherit.aes = FALSE, size = 9) +
  scale_colour_manual(name = "Condition", values = condColours) +
  scale_fill_manual(name = "Condition", values = condColours) +
  labs(x = "TMS intensity (% RMT)", y = "Effect of accuracy on zMEP (β)") +
  theme_minimal()

fig4 <- (fig4A | fig4B) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12), legend.text = element_text(size = 12),
        axis.title = element_text(size = 13), axis.text = element_text(size = 12),
        plot.tag = element_text(size = 18, face = "bold"))
print(fig4)
ggsave(file.path(outDir, "Fig4_accuracy_condition.png"), fig4, width = 12, height = 5.5, dpi = 300)

# ============================================================================
# 4. CONTROL ANALYSES & ROBUSTNESS
# (stability of accuracy and VMIQ-2 correlations: see the accuracy script)
# ============================================================================

# Imagery ability as a second covariate
vmiq <- read_excel(file.path(derivDir, "vmiq_womi.xlsx")) %>%
  transmute(participant = Participant,
            vmiq_mean   = (External + Internal + Kinaesthetic) / 3,
            vmiq_kin    = Kinaesthetic) %>%
  mutate(across(starts_with("vmiq"), ~ .x - mean(.x)))
curvesVmiq <- left_join(curves, vmiq, by = "participant")

for (cov in c("vmiq_mean", "vmiq_kin")) {
  cat("\n--- Slope ANCOVA with accuracy +", cov, "---\n")
  m <- aov_car(as.formula(paste("slope ~ accuracy_c +", cov,
                                "+ condition * weight + Error(id / (condition * weight))")),
               data = curvesVmiq, factorize = FALSE)
  print(fmtTable(rmTable(m)), n = Inf)
}

# Table ST4: session order
slopeOrderModel <- aov_car(slope ~ order + accuracy_c + condition * weight +
                             Error(id / (condition * weight)),
                           data = curves, factorize = FALSE)
print(fmtTable(rmTable(slopeOrderModel)), n = Inf)

# Slope by order x condition and order x weight (M, SD across participants)
print(gain %>% group_by(order, condition) %>%
        summarise(M = mean(slope), SD = sd(slope), n = n(), .groups = "drop"))
print(curves %>% group_by(order, weight) %>%
        summarise(M = mean(slope), SD = sd(slope), n = n(), .groups = "drop"))

# No-video control trials: zMEP in AO vs AOMI session
ctrl <- read_csv(file.path(derivDir, "mep_ztransform", "mean_z_scores_control_vs_main.csv"),
                 show_col_types = FALSE)
print(t.test(ctrl$AO_Control, ctrl$AOMI_Control, paired = TRUE))

# FDI vs ADM: muscle x condition x intensity on raw MEPs
muscleData <- read_csv(file.path(derivDir, "ADM_FDI_raw_MEPs_for_ANOVA.csv"),
                       show_col_types = FALSE) %>%
  mutate(id        = factor(Participant),
         muscle    = factor(Muscle, levels = c("FDI", "ADM")),
         condition = factor(Condition, levels = c("AO", "AOMI")),
         intensity = factor(Intensity, levels = 1:5, labels = intensityLabels),
         mep       = MEP_Raw)
print(count(muscleData, muscle))
print(count(muscleData, muscle, intensity))

print(muscleData %>% group_by(id, muscle) %>% summarise(mep = mean(mep), .groups = "drop") %>%
        group_by(muscle) %>% summarise(M = mean(mep), SD = sd(mep)))

muscleModel <- aov_car(mep ~ muscle * condition * intensity +
                         Error(id / (muscle * condition * intensity)),
                       data = muscleData, fun_aggregate = mean)
print(fmtTable(rmTable(muscleModel)), n = Inf)
print(eta_squared(muscleModel, ci = .95) %>% filter(Parameter == "muscle"))

# Figure S3
muscleMeans <- muscleData %>%
  group_by(id, muscle, intensity) %>%
  summarise(mep = mean(mep), .groups = "drop")
pdS3 <- position_dodge(0.25)
figS3 <- ggplot(muscleMeans, aes(intensity, mep, colour = muscle, group = muscle)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2, position = pdS3) +
  stat_summary(fun = mean, geom = "point", size = 3, position = pdS3) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar",
               width = 0.15, position = pdS3) +
  scale_colour_manual(values = c(FDI = "#542788", ADM = "#B2ABD2")) +
  labs(x = "TMS intensity (% RMT)", y = "MEP amplitude (mV)", colour = "Muscle") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12))
print(figS3)
ggsave(file.path(outDir, "FigS3_mep_by_muscle.png"), figS3, width = 8, height = 6, dpi = 300)

# Statistical robustness -----------------------------------------------------

# Residuals of the slope ANCOVA (Figure S4)
slopeResiduals <- as.numeric(residuals(slopeAccModel))
print(shapiro.test(slopeResiduals))

figS4 <- ggplot(data.frame(residual = slopeResiduals), aes(sample = residual)) +
  stat_qq(alpha = 0.5, size = 1.8, colour = "grey30") +
  stat_qq_line(colour = "#E6550D", linewidth = 1) +
  labs(x = "Theoretical quantiles (normal distribution)",
       y = "Sample quantiles (slope ANCOVA residuals)") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12))
print(figS4)
ggsave(file.path(outDir, "FigS4_residual_qq_slope.png"), figS4, width = 6, height = 5, dpi = 300)

# Accuracy x condition = correlation of accuracy with the AO - AOMI slope difference
gainDiff <- gainWide %>%
  left_join(distinct(gain, id, accuracy), by = "id") %>%
  mutate(diff = AO - AOMI)
print(cor.test(gainDiff$accuracy, gainDiff$diff))

# Permutation test
set.seed(1)
obsR  <- cor(gainDiff$accuracy, gainDiff$diff)
permR <- replicate(10000, cor(gainDiff$accuracy, sample(gainDiff$diff)))
cat(sprintf("Permutation test: r = %.2f, p = %.4f (10,000 permutations, two-sided)\n",
            obsR, mean(abs(permR) >= abs(obsR))))

# Leave-one-participant-out (Figure S5)
loo <- map_dfr(levels(curves$id), function(pid) {
  d <- curves %>% filter(id != pid) %>% droplevels()
  d$accuracy_c <- d$accuracy - mean(distinct(d, id, accuracy)$accuracy)
  m <- suppressMessages(aov_car(slope ~ accuracy_c + condition * weight +
                                  Error(id / (condition * weight)),
                                data = d, factorize = FALSE))
  tab <- anova(m, correction = "none", es = "pes")["accuracy_c:condition", ]
  tibble(removed = pid, F = tab$F, df1 = tab$`num Df`, df2 = tab$`den Df`,
         p = tab$`Pr(>F)`, pes = tab$pes)
})
print(loo %>% arrange(desc(p)), n = Inf)
cat(sprintf("F = %.2f-%.2f, p = %.3f-%.3f, ηp² = %.2f-%.2f, significant in %d of %d runs\n",
            min(loo$F), max(loo$F), min(loo$p), max(loo$p),
            min(loo$pes), max(loo$pes), sum(loo$p < .05), nrow(loo)))

figS5 <- ggplot(loo, aes(fct_reorder(removed, p), p)) +
  geom_hline(yintercept = .05, linetype = "dashed", colour = "#E6550D") +
  geom_point(size = 3, colour = "grey20") +
  scale_y_log10() +
  labs(x = "Participant removed", y = "p, accuracy × condition on slope\n(log scale)") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1))
print(figS5)
ggsave(file.path(outDir, "FigS5_leave_one_out_slope.png"), figS5, width = 9, height = 5, dpi = 300)