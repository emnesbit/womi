# ============================================================================
# WOMI - R package requirements
# Run once before the analysis scripts: source("requirements.R")
# ============================================================================

if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse,     # data handling and plotting (incl. ggplot2, readr, tidyr, dplyr, purrr)
  readxl,        # read .xlsx input files
  writexl,       # write .xlsx output files
  afex,          # repeated-measures ANOVA / ANCOVA (aov_car)
  emmeans,       # post-hoc comparisons
  effectsize,    # Cohen's d, partial eta squared with CIs
  broom,         # tidy regression output (simple slopes)
  patchwork,     # multi-panel figures
  RColorBrewer   # colour palettes
)
