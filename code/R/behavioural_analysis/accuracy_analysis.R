# ============================================================================
# WOMI - Weight estimation accuracy
# Follows the Results: 1 (Fig. 2), 4 "Stability of weight estimation
# accuracy" (Figs. S1, S2, VMIQ-2)
# Writes accuracy_by_condition.xlsx, which the MEP script reads.
# ============================================================================

pacman::p_load(tidyverse, readxl, afex, emmeans, effectsize, patchwork, RColorBrewer, writexl)

bidsDir <- "/data/path/to/bids_directory/"
outDir  <- file.path(bidsDir, "derivatives", "behavioural_analysis")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

condColours     <- c(AO = "#F8766D", AOMI = "#00BFC4")
intensityLabels <- c("90%", "100%", "110%", "120%", "130%")

# Proportion as percentage with two decimals
pct <- function(x) sprintf("%.2f%%", 100 * x)

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

# ---------------------------------------------------------------------------
# Data: one row per valid video trial (no control trials, no trial 1)
# ---------------------------------------------------------------------------
readSession <- function(participant, condition) {
  sub <- sprintf("sub-%02d", participant)
  f   <- file.path(bidsDir, sub, paste0("ses-", condition), "beh",
                   sprintf("%s_ses-%s_womi_beh.tsv", sub, condition))
  read_tsv(f, show_col_types = FALSE) %>%
    group_by(trial) %>%
    summarise(control   = any(event_type == "Control_TMS_pulse"),
              shown     = as.numeric(code[event_type == "weight_in_video"][1]),
              response  = as.numeric(code[event_type == "weight_estimation_response"][1]),
              intensity = as.numeric(code[event_type == "TMS_pulse_intensity"][1]),
              .groups = "drop") %>%
    filter(!control, trial != 1, !is.na(shown), !is.na(response)) %>%
    mutate(participant = participant, condition = condition, correct = shown == response)
}

trials <- expand_grid(participant = 1:18, condition = c("AO", "AOMI")) %>%
  pmap_dfr(readSession) %>%
  mutate(id        = factor(participant),
         condition = factor(condition, levels = c("AO", "AOMI")),
         order     = factor(if_else(participant %% 2 == 0, "AO first", "AOMI first")),
         weight    = factor(response, levels = 1:3, labels = c("Light", "Moderate", "Heavy")))

acc <- trials %>%
  group_by(participant, order, condition) %>%
  summarise(accuracy = mean(correct), .groups = "drop")

accWide <- acc %>% pivot_wider(names_from = condition, values_from = accuracy)
print(accWide, n = Inf)
write_xlsx(select(accWide, participant, AO, AOMI), file.path(outDir, "accuracy_by_condition.xlsx"))

# ============================================================================
# 1. WEIGHT ESTIMATION ACCURACY
# ============================================================================

# Accuracy vs chance (1/3), one-sided
chance <- map_dfr(c("AO", "AOMI"), function(cond) {
  x  <- accWide[[cond]]
  tt <- t.test(x, mu = 1/3, alternative = "greater")
  tibble(condition = cond, M = pct(mean(x)), SD = pct(sd(x)),
         t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
         d = (mean(x) - 1/3) / sd(x), n_above = sum(x > 1/3), min = pct(min(x)))
})
print(chance, width = Inf)

# AO vs AOMI
ttCond <- t.test(accWide$AO, accWide$AOMI, paired = TRUE)
print(ttCond)
cat(sprintf("Mean difference = %s, CI [%s, %s]\n",
            pct(ttCond$estimate), pct(ttCond$conf.int[1]), pct(ttCond$conf.int[2])))
print(cohens_d(accWide$AO, accWide$AOMI, paired = TRUE))
print(cor.test(accWide$AO, accWide$AOMI))

# Condition x weight response
weightAcc <- trials %>%
  group_by(id, condition, weight) %>%
  summarise(accuracy = mean(correct), .groups = "drop")

weightModel <- aov_car(accuracy ~ condition * weight + Error(id / (condition * weight)),
                       data = weightAcc)
print(rmTable(weightModel))
print(eta_squared(weightModel, ci = .95))

# Post-hoc: weight response, pooled over conditions (Bonferroni)
posthoc <- emmeans(weightModel, pairwise ~ weight, adjust = "bonferroni")
as.data.frame(posthoc$emmeans) %>% transmute(weight, M = pct(emmean)) %>% print()
as.data.frame(posthoc$contrasts) %>%
  transmute(contrast, MD = pct(abs(estimate)), SE = pct(SE),
            t = round(t.ratio, 2), p = round(p.value, 3)) %>%
  print()

# Response counts
counts <- count(trials, weight)
print(counts)
print(chisq.test(counts$n))

# Figure 2: accuracy by weight response, one colour per participant (same in AO and AOMI)
stars        <- function(p) case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ "n.s.")
weightLevels <- levels(weightAcc$weight)
idColours    <- setNames(scales::hue_pal(l = 65, c = 120)(nlevels(weightAcc$id)), levels(weightAcc$id))

brackets <- as.data.frame(summary(posthoc$contrasts)) %>%
  separate(contrast, into = c("from", "to"), sep = " - ") %>%
  mutate(x1 = match(from, weightLevels), x2 = match(to, weightLevels), label = stars(p.value)) %>%
  arrange(x2 - x1) %>%
  mutate(y = 1.02 + 0.07 * (row_number() - 1))

weightPanel <- function(cond) {
  ggplot(filter(weightAcc, condition == cond), aes(weight, accuracy, group = id)) +
    geom_hline(yintercept = 1/3, linetype = "dashed", colour = "grey40") +
    geom_line(colour = "grey80", alpha = 0.6) +
    geom_point(aes(colour = id), size = 3) +
    annotate("text", x = 0.5, y = 1/3, label = "chance (33%)", vjust = -0.5, hjust = 0,
             size = 3.5, colour = "grey40") +
    stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 18, size = 5, colour = "black") +
    stat_summary(aes(group = 1), fun.data = mean_sdl, fun.args = list(mult = 1),
                 geom = "errorbar", width = 0.1, colour = "black") +
    geom_segment(data = brackets, inherit.aes = FALSE, aes(x = x1, xend = x2, y = y, yend = y)) +
    geom_segment(data = brackets, inherit.aes = FALSE, aes(x = x1, xend = x1, y = y, yend = y - 0.015)) +
    geom_segment(data = brackets, inherit.aes = FALSE, aes(x = x2, xend = x2, y = y, yend = y - 0.015)) +
    geom_text(data = brackets, inherit.aes = FALSE,
              aes(x = (x1 + x2) / 2, y = y + 0.01, label = label), vjust = 0, size = 4.5) +
    scale_colour_manual(name = "Participant", values = idColours) +
    scale_y_continuous(limits = c(min(weightAcc$accuracy, 1/3) - 0.03, max(brackets$y) + 0.06),
                       breaks = seq(0, 1, 0.2)) +
    labs(title = cond, x = "Weight Response", y = "Proportion Correct") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 13), axis.text = element_text(size = 12))
}

conditionBracket <- ggplot() +
  annotate("segment", x = 0.25, xend = 0.75, y = 0.5, yend = 0.5) +
  annotate("segment", x = c(0.25, 0.75), xend = c(0.25, 0.75), y = 0.5, yend = 0.1) +
  annotate("text", x = 0.5, y = 0.65, vjust = 0, size = 4.5,
           label = stars(anova(weightModel)["condition", "Pr(>F)"])) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.2), expand = c(0, 0)) +
  theme_void()

fig2 <- conditionBracket / (weightPanel("AO") | weightPanel("AOMI")) +
  plot_layout(heights = c(0.08, 1), guides = "collect") &
  theme(legend.position = "right")
print(fig2)
ggsave(file.path(outDir, "Fig2_accuracy_by_weight.png"), fig2, width = 12, height = 6, dpi = 300)

# ============================================================================
# 4. CONTROL ANALYSES - Stability of weight estimation accuracy
# ============================================================================

# First vs second session
sessionAcc <- acc %>%
  mutate(session = if_else((order == "AO first") == (condition == "AO"), "Session 1", "Session 2"))
sessionWide <- sessionAcc %>%
  select(participant, session, accuracy) %>%
  pivot_wider(names_from = session, values_from = accuracy)
ttSession <- t.test(sessionWide$`Session 1`, sessionWide$`Session 2`, paired = TRUE)
print(ttSession)
cat(sprintf("Mean difference (session 1 - 2) = %s\n", pct(ttSession$estimate)))

# Figure S1
figS1 <- ggplot(sessionAcc, aes(session, accuracy, group = participant)) +
  geom_line(colour = "grey70", alpha = 0.7) +
  geom_point(aes(colour = condition), size = 3) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 18, size = 4) +
  stat_summary(aes(group = 1), fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.1) +
  scale_colour_manual(name = "Condition", values = condColours) +
  labs(x = "Session", y = "Proportion Correct") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12))
print(figS1)
ggsave(file.path(outDir, "FigS1_accuracy_by_session.png"), figS1, width = 8, height = 6, dpi = 300)

# AO first vs AOMI first (mean accuracy across conditions, Welch)
composite <- accWide %>% mutate(accuracy = (AO + AOMI) / 2)
print(composite %>% group_by(order) %>% summarise(M = pct(mean(accuracy)), SD = pct(sd(accuracy))))
print(t.test(accuracy ~ order, data = composite))

# TMS intensity
tmsAcc <- trials %>%
  filter(!is.na(intensity)) %>%
  mutate(intensity = factor(intensity, levels = 1:5, labels = intensityLabels)) %>%
  group_by(id, intensity) %>%
  summarise(accuracy = mean(correct), .groups = "drop")

tmsModel <- aov_car(accuracy ~ intensity + Error(id / intensity), data = tmsAcc)
print(rmTable(tmsModel))

# Figure S2
tmsColours <- setNames(brewer.pal(6, "Greens")[2:6], intensityLabels)
figS2 <- ggplot(tmsAcc, aes(intensity, accuracy, group = id)) +
  geom_line(colour = "#C7E9C0", alpha = 0.7) +
  geom_point(aes(colour = intensity), size = 3, show.legend = FALSE) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 18, size = 5, colour = "#00441B") +
  stat_summary(aes(group = 1), fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.15, colour = "#00441B") +
  scale_colour_manual(values = tmsColours) +
  labs(x = "TMS intensity (% RMT)", y = "Proportion Correct") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 12))
print(figS2)
ggsave(file.path(outDir, "FigS2_accuracy_by_tms_intensity.png"), figS2, width = 8, height = 6, dpi = 300)

# Accuracy vs imagery vividness (VMIQ-2), Holm-corrected over 8 tests
vmiq <- read_excel(file.path(bidsDir, "derivatives", "vmiq_womi.xlsx")) %>%
  transmute(participant  = Participant,
            external     = External,
            internal     = Internal,
            kinaesthetic = Kinaesthetic,
            mean         = (external + internal + kinaesthetic) / 3)

accWide %>%
  left_join(vmiq, by = "participant") %>%
  pivot_longer(c(AO, AOMI), names_to = "condition", values_to = "accuracy") %>%
  pivot_longer(c(mean, external, internal, kinaesthetic), names_to = "scale", values_to = "score") %>%
  group_by(condition, scale) %>%
  summarise(r = cor(accuracy, score), df = n() - 2,
            p = cor.test(accuracy, score)$p.value, .groups = "drop") %>%
  mutate(p_holm = p.adjust(p, "holm")) %>%
  print()