# ============================================================================
# WOMI R Package Requirements
# ============================================================================
# This script installs all required R packages for the WOMI analysis pipeline
# Run this before starting the analysis
# ============================================================================

cat("=== WOMI2 R Package Installation ===\n\n")

# Check R version
r_version <- as.numeric(paste(R.version$major, R.version$minor, sep = "."))
if (r_version < 4.0) {
  warning("R version 4.0 or later is recommended. Current version: ", R.Version()$version.string)
  cat("Some packages may not install correctly.\n")
  cat("Consider updating R from: https://cran.r-project.org/\n\n")
} else {
  cat("R version check: OK (", R.Version()$version.string, ")\n\n")
}

# Install pacman if not available
if (!require("pacman", quietly = TRUE)) {
  cat("Installing pacman for package management...\n")
  install.packages("pacman", repos = "https://cloud.r-project.org")
}

library(pacman)

# Define required packages with descriptions
required_packages <- list(
  # Data manipulation and import/export
  tidyverse = "Data manipulation, visualization, and workflow tools",
  readxl = "Read Excel files (.xlsx, .xls)",
  writexl = "Write Excel files (.xlsx)",
  haven = "Import/export for SPSS, Stata, SAS files",
  
  # Statistical analysis
  afex = "Analysis of Factorial Experiments (ANOVA/ANCOVA)",
  emmeans = "Estimated Marginal Means (post-hoc tests)",
  effectsize = "Effect size calculations (Cohen's d, eta squared)",
  car = "Companion to Applied Regression (Levene's test, VIF)",
  lme4 = "Linear Mixed-Effects Models",
  lmerTest = "Tests for linear mixed models (p-values)",
  
  # Descriptive statistics
  psych = "Psychological research tools (descriptive statistics)",
  
  # Visualisation
  patchwork = "Combine multiple plots into one figure",
  RColorBrewer = "Color palettes for plots",
  gghalves = "Half-violin plots (optional, for specific visualizations)"
)

# Install all packages
cat("Installing required packages...\n")
cat("This may take several minutes on first run.\n\n")

# Use pacman to install/load all packages
packages_to_install <- names(required_packages)

# Try to load packages, install if missing
p_load(char = packages_to_install, install = TRUE, character.only = TRUE)

cat("\n=== Verifying Installation ===\n\n")

# Verify each package loaded successfully
all_loaded <- TRUE
for (pkg in packages_to_install) {
  if (pkg %in% loadedNamespaces()) {
    cat(sprintf("✓ %s: OK\n", pkg))
  } else {
    cat(sprintf("✗ %s: FAILED\n", pkg))
    all_loaded <- FALSE
  }
}

if (all_loaded) {
  cat("\n=== All packages installed successfully! ===\n\n")
  cat("You can now run the analysis scripts:\n")
  cat("  1. code/R/behavioural/accuracy.R\n")
  cat("  2. code/R/statistical_analysis/main_3way_ANCOVA_with_accuracy.R\n\n")
} else {
  cat("\n=== Some packages failed to install ===\n\n")
  cat("Troubleshooting steps:\n")
  cat("  1. Update R to the latest version\n")
  cat("  2. Try installing failed packages individually:\n")
  cat("     install.packages('package_name')\n")
  cat("  3. Check for system dependencies (Linux users)\n")
  cat("  4. Restart R and try again\n\n")
}

# Print package versions for reproducibility
cat("=== Installed Package Versions ===\n\n")
for (pkg in packages_to_install) {
  if (pkg %in% loadedNamespaces()) {
    version <- packageVersion(pkg)
    cat(sprintf("%s: %s\n", pkg, version))
  }
}

cat("\n=== System Information ===\n\n")
cat("R version:", R.Version()$version.string, "\n")
cat("Platform:", R.Version()$platform, "\n")
cat("OS:", Sys.info()["sysname"], Sys.info()["release"], "\n")
cat("RStudio:", ifelse(Sys.getenv("RSTUDIO") == "1", "Yes", "No"), "\n")

cat("\n=== Installation Complete ===\n")
cat("Save this output for troubleshooting if needed.\n")

# Optional: Create a session info file
session_info_file <- "session_info.txt"
sink(session_info_file)
cat("=== WOMI2 R Session Information ===\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
print(sessionInfo())
sink()
cat("\nSession info saved to:", session_info_file, "\n")