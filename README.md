# WOMI: Weight Observation and Motor Imagery Analysis Pipeline

## Overview

This repository contains the complete data processing and statistical analysis pipeline for the WOMI (Weight Observation and Motor Imagery) experiment. The study investigates motor cortex excitability changes during action observation (AO) and action observation combined with motor imagery (AOMI) using transcranial magnetic stimulation (TMS) and electromyography (EMG).

**Key Research Question**: Does combining motor imagery with action observation enhance motor cortex facilitation compared to action observation alone?

---

## Repository Structure

```
womi_code/
├── code/
│   ├── matlab/
│   │   ├── preprocessing/
│   │   │   └── main.m                            # MEP extraction and quality control
│   │   └── analysis/
│   │       └── z-transformation.m                # Within-participant standardisation
│   └── R/
│       ├── behavioural_analysis/
│       │   └── accuracy.R                        # Weight estimation accuracy analysis
│       └── EMG_analysis/
│           └── main_3way_ANCOVA_with_accuracy.R  # Main statistical model
├── docs/
│   ├── MATLAB_DOCUMENTATION.md                   # Complete MATLAB reference
│   └── R_DOCUMENTATION.md                        # Complete R reference
├── requirements.R                                # R package dependencies
├── LICENSE.md                                    # MIT license file
└── README.md                                     # This file
```

---

## Quick Start

### For First-Time Users

1. **Read the documentation**:
   - [`MATLAB_DOCUMENTATION.md`](docs/MATLAB_DOCUMENTATION.md) for MATLAB details
   - [`R_DOCUMENTATION.md`](docs/R_DOCUMENTATION.md) for R analysis details

2. **Install R dependencies**:
   ```r
   source("requirements.R")
   ```

3. **Set paths**: Update `bidsDirectory` in each script to point to your local data location (see Configuration below).

4. **Run in order** (see Pipeline section below).

### For Experienced Users

```matlab
% MATLAB: Extract and z-transform MEPs
cd('code/matlab/preprocessing/')
main()
cd('../analysis/')
z-transformation()
```

```r
# R: Behavioural accuracy then ANCOVA
source("code/R/behavioural_analysis/accuracy.R")
source("code/R/EMG_analysis/main_3way_ANCOVA_with_accuracy.R")
```

**Total runtime**: ~10 minutes for 18 participants.

---

## Experimental Design

### Design Structure

**Within-Subjects Design**: 2 (Condition) × 3 (Weight) × 5 (TMS Intensity)

| Factor | Levels | Description |
|--------|--------|-------------|
| **Condition** | AO, AOMI | Action Observation vs. AO + Motor Imagery |
| **Weight** | Light, Moderate, Heavy | Subjective weight estimation (1–3) |
| **TMS Intensity** | 90%, 100%, 110%, 120%, 130% | Percentage of resting motor threshold |

### Sample

- **N = 18** participants
- **Sessions**: 2 per participant (AO and AOMI)
- **Trials**: 122 per session (before exclusion)
- **Total trials**: 4,392 initial → 4,220 final (3.9% excluded)

### Counterbalancing

- Even participant IDs (2, 4, 6 ...): AO first, then AOMI
- Odd participant IDs (1, 3, 5 ...): AOMI first, then AO

Note: Only main trials (weight judgement trials) are included in the ANCOVA. Control trials are analysed separately at the aggregation stage.

---

## Data Requirements

This pipeline expects data in [BIDS](https://bids.neuroimaging.io/) format. Data are documented separately in the data repository README on OSF.

**Expected directory structure**:

```
/data/path/to/bids_directory/
├── sub-01/
│   ├── ses-AO/
│   │   ├── emg/
│   │   │   └── sub-01_ses-AO_womi_emg.tsv.gz
│   │   └── beh/
│   │       └── sub-01_ses-AO_womi_beh.tsv
│   └── ses-AOMI/
│       ├── emg/
│       │   └── sub-01_ses-AOMI_womi_emg.tsv.gz
│       └── beh/
│           └── sub-01_ses-AOMI_womi_beh.tsv
├── sub-02/
└── ...
```

**Key file contents**:

- `*_emg.tsv.gz`: `time_relative_to_TMS_pulse_in_ms`, `FDI_in_mV`
- `*_beh.tsv`: `trial`, `event_type`, `code`

---

## Analysis Pipeline

```
Raw BIDS Data (EMG + Behavioural)
         ↓
┌────────────────────────────────┐
│  MATLAB Stage 1: MEP Extract   │
│  - Peak-to-peak calculation    │
│  - Quality control filtering   │
│  - Baseline RMS check          │
└────────────────────────────────┘
         ↓
Cleaned MEP Data (per participant/session)
         ↓
┌────────────────────────────────┐
│  MATLAB Stage 2: Z-Transform   │
│  - Within-participant std      │
│  - Aggregate by conditions     │
└────────────────────────────────┘
         ↓
Standardised MEP Data + Summary CSVs
         ↓
┌────────────────────────────────┐
│  R Stage 1: Behavioural Acc    │
│  - Accuracy by condition       │
│  - Session order effects       │
│  - Weight response patterns    │
└────────────────────────────────┘
         ↓
Accuracy Tables + Visualisations
         ↓
┌────────────────────────────────┐
│  R Stage 2: Statistical Model  │
│  - Repeated measures ANCOVA    │
│  - Assumption checks           │
│  - Post-hoc analyses           │
└────────────────────────────────┘
         ↓
ANCOVA Results + Figures
```

---

## Configuration

Update the path at the top of each script before running:

| Script | Variable | Default |
|--------|----------|---------|
| `main.m` | `config.bidsDirectory` | `/data/path/to/bids_directory/` |
| `z-transformation.m` | `config.bidsDirectory` | `/data/path/to/bids_directory/` |
| `accuracy.R` | `bidsDirectory` | `/data/path/to/bids_directory/` |
| `main_3way_ANCOVA_with_accuracy.R` | `bidsDirectory` | `/data/path/to/bids_directory/` |

**Key analysis parameters** (MATLAB `main.m`):

| Parameter | Value | Description |
|-----------|-------|-------------|
| MEP window | 15–50 ms post-TMS | Peak-to-peak amplitude window |
| Baseline window | −40 to −10 ms pre-TMS | Baseline RMS window |
| Pre-contraction threshold | 0.15 mV | Max peak-to-peak before TMS |
| Baseline RMS threshold | 0.05 mV | Max RMS in baseline window |

---

## Output Files

All outputs are saved to BIDS-compatible derivatives directories.

### MATLAB Stage 1 → `/derivatives/mep_analysis/`

- `sub-XX_ses-{AO,AOMI}_mep.tsv` — trial-level MEP data (BIDS)
- `sub-XX_ses-{AO,AOMI}_mep.mat` — MATLAB format
- 72 files total (36 TSV + 36 MAT)

### MATLAB Stage 2 → `/derivatives/mep_ztransform/`

- `sub-XX_ses-{AO,AOMI}_zmep.tsv` — z-transformed MEP data (BIDS)
- `sub-XX_ses-{AO,AOMI}_zmep.mat` — MATLAB format
- `mean_z_scores_control_vs_main.csv`
- `mean_z_scores_by_intensity.csv`
- `mean_z_scores_by_intensity_and_weight.csv` ← input to R Stage 2

### R Stage 1 → `/derivatives/behavioural_analysis/`

Excel tables:
- `accuracy_by_condition.xlsx` ← input to R Stage 2
- `accuracy_by_session_order.xlsx`
- `accuracy_by_weight.xlsx`

Figures (PNG, 300 DPI):
- `session_order_plot.png`
- `condition_comparison_plot.png`
- `condition_comparison_for_panel.png`
- `weight_by_condition_plot.png`
- `weight_responses_plot.png`
- `behavioural_manuscript_figure.png`

### R Stage 2 → `/derivatives/ancova_analysis/`

Results:
- `merged_data.csv`
- `ancova_results.txt`
- `pairwise_intensity.csv`
- `simple_slopes_highest_intensity.csv`
- `interaction_by_intensity.csv`

Assumption check plots:
- `assumption_accuracy_distribution.png`
- `assumption_linearity.png`

Figures (PNG, 300 DPI):
- `main_effect_condition.png`
- `main_effect_intensity.png`
- `interaction_highest_intensity.png`
- `interaction_all_intensities_faceted.png`
- `interaction_all_intensities_grid.png`
- `effect_sizes_by_intensity.png`
- `slopes_by_intensity.png`
- `marginal_means_intensity.png`
- `combined_poster_plot.png`
- `manuscript_figure_combined.png`

---

## Statistical Methods

### Within-Participant Z-Transformation

MEPs are standardised within each participant across both sessions:

```
z = (MEP − participant_mean) / participant_SD
```

This removes individual baseline differences (e.g., skull thickness, coil distance) while preserving within-participant experimental effects.

### Repeated Measures ANCOVA

```
zMEP ~ accuracy_c + condition × weight × intensity
       + Error(participant / (condition × weight × intensity))
```

- **Within-subjects factors**: condition (2), weight (3), intensity (5)
- **Covariate**: behavioural accuracy (mean-centred)
- **Sphericity correction**: Greenhouse–Geisser where Mauchly's test is significant
- **Effect sizes**: partial η²
- **Post-hoc**: Bonferroni-corrected pairwise comparisons for intensity; simple slopes if accuracy interactions are significant

---

## Prerequisites

### Software

- **MATLAB** R2019b or later (Statistics and Machine Learning Toolbox required)
- **R** 4.0 or later (RStudio recommended)

### R Packages

```r
source("requirements.R")
```

Or manually:

```r
install.packages("pacman")
pacman::p_load(
  tidyverse, haven, psych, writexl, afex, emmeans,
  patchwork, effectsize, RColorBrewer, readxl,
  car, lme4, lmerTest
)
```

---

## Citation

If you use this pipeline, please cite:

```bibtex
@article{womi_2026,
  title   = {Individual Differences in Corticospinal Excitability During Combined Action Observation and Motor Imagery Depend on Task Approach},
  author  = {Lange, Yassamin and Nesbit, Emma and Hermsdoerfer, Joachim and Stadler, Waltraud},
  journal = {[Journal Name]},
  year    = {2026},
  doi     = {[DOI when available]}
}
```

---

## License

MIT — see `LICENSE.md`.

---

## Contributors

- **Study Design**: Yassamin Lange, Emma Nesbit, Waltraud Stadler
- **Data Collection**: Yassamin Lange, Emma Nesbit
- **Analysis Pipeline**: Emma Nesbit
- **Code Review**: Yassamin Lange, Waltraud Stadler

---

## Contact

**Emma Nesbit** — nesbit@cbs.mpg.de  
Max Planck Institute for Human Cognitive and Brain Sciences  
GitHub Issues for bug reports· GitHub Discussions for questions

---

## Acknowledgments

We thank the Chair of Human Movement Science at TUM for providing facilities and all participants who volunteered their time.

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | April 2026 | Initial release with complete BIDS pipeline |

---

**Tested with**: MATLAB R2023a · R 4.3.1