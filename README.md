# WOMI: Weight Observation and Motor Imagery Analysis Pipeline

## Overview

This repository contains the data processing and statistical analysis pipeline for the WOMI (Weight Observation and Motor Imagery) experiment, accompanying the article *Individual Task Approach Modulates Corticospinal Excitability During Combined Action Observation and Motor Imagery*.

The study compares corticospinal excitability (CSE), measured as motor evoked potentials (MEPs) elicited by transcranial magnetic stimulation (TMS), during action observation (AO) and action observation combined with motor imagery (AOMI) in a weight estimation task (WET).

**Hypotheses**:
1. AOMI leads to higher CSE than AO.
2. Observing or imagining the lifting of a heavier weight leads to higher CSE.

The main outcome is the slope (gain) of the MEP recruitment curve across five TMS intensities. Weight estimation accuracy is analysed as an exploratory covariate.

---

## Repository Structure

```
womi_code/
├── code/
│   ├── matlab/
│   │   ├── preprocessing/
│   │   │   └── main.m                    # MEP extraction and quality control
│   │   └── analysis/
│   │       └── z_transformation.m        # Within-participant z-transformation and aggregation
│   └── R/
│       ├── behavioural_analysis/
│       │   └── accuracy_analysis.R       # Weight estimation accuracy (Results 1, 4)
│       └── EMG_analysis/
│           └── mep_analysis.R            # MEP analyses (Results 2, 3, 4)
├── docs/
│   ├── MATLAB_DOCUMENTATION.md
│   └── R_DOCUMENTATION.md
├── requirements.R                        # R package dependencies
├── LICENSE.md                            # MIT licence
└── README.md
```

---

## Quick Start

1. **Get the data** from OSF (see Data Requirements) and place it in a BIDS directory.
2. **Install R packages**:
   ```r
   source("requirements.R")
   ```
3. **Set the path** to your BIDS directory at the top of each script (see Configuration).
4. **Run the scripts in this order**:

```matlab
% MATLAB
cd('code/matlab/preprocessing/');  main
cd('../analysis/');                z_transformation
```

```r
# R
source("code/R/behavioural_analysis/accuracy_analysis.R")   # writes accuracy_by_condition.xlsx
source("code/R/EMG_analysis/mep_analysis.R")                # needs accuracy_by_condition.xlsx
```

All statistics and figures are printed to the console and figures are saved as PNG files (300 dpi).

---

## Experimental Design

**Within-subject design**: 2 (condition) × 3 (weight response) × 5 (TMS intensity)

| Factor | Levels | Description |
|--------|--------|-------------|
| Condition | AO, AOMI | Action observation vs. action observation with motor imagery |
| Weight response | Light, Moderate, Heavy | Weight reported by the participant on each trial |
| TMS intensity | 90, 100, 110, 120, 130% RMT | Percentage of resting motor threshold |

- **N = 18** participants, two sessions each (AO and AOMI, on separate days)
- **Trials**: 122 per session; 4392 in total, 4220 after exclusions (3.9% excluded)
- **Session order**: even participant IDs started with AO, odd IDs with AOMI (n = 9 each)
- Each session also contained no-video control trials, which are analysed separately.

---

## Data Requirements

Data are available on OSF and follow the [BIDS](https://bids.neuroimaging.io/) format. The data repository has its own README.

```
/data/path/to/bids_directory/
├── sub-01/
│   ├── ses-AO/
│   │   ├── emg/sub-01_ses-AO_womi_emg.tsv.gz
│   │   └── beh/sub-01_ses-AO_womi_beh.tsv
│   └── ses-AOMI/
│       ├── emg/sub-01_ses-AOMI_womi_emg.tsv.gz
│       └── beh/sub-01_ses-AOMI_womi_beh.tsv
├── sub-02/
├── ...
└── derivatives/
    ├── vmiq_womi.xlsx                    # VMIQ-2 scores (from OSF)
    └── ADM_FDI_raw_MEPs_for_ANOVA.csv    # raw FDI and ADM MEPs (from OSF)
```

- `*_emg.tsv.gz`: `time_relative_to_TMS_pulse_in_ms`, `FDI_in_mV`, `ADM_in_mV`, `trial`
- `*_beh.tsv`: `trial`, `event_type`, `code`
- `vmiq_womi.xlsx`: `Participant`, `External`, `Internal`, `Kinaesthetic` (VMIQ-2 subscale scores)
- `ADM_FDI_raw_MEPs_for_ANOVA.csv`: `Participant`, `Condition`, `Muscle` (FDI/ADM), `Intensity` (1–5), `MEP_Raw` (mV); trial-level raw MEPs of both muscles, preprocessed with an RMS threshold of 0.10 mV for comparable trial retention

---

## Analysis Pipeline

```
BIDS data (EMG + behaviour)
        ↓
MATLAB 1  main.m               MEP extraction, quality control
        ↓
MATLAB 2  z_transformation.m   within-participant z-scores, aggregation
        ↓
R 1       accuracy_analysis.R  weight estimation accuracy
        ↓
R 2       mep_analysis.R       recruitment-curve slopes, exploratory ANCOVAs,
                               control and robustness analyses
```

---

## Configuration

Set the BIDS directory at the top of each script:

| Script | Variable |
|--------|----------|
| `main.m` | `config.bidsDirectory` |
| `z_transformation.m` | `config.bidsDirectory` |
| `accuracy_analysis.R` | `bidsDir` |
| `mep_analysis.R` | `bidsDir` |

**MEP preprocessing parameters** (`main.m`):

| Parameter | Value |
|-----------|-------|
| MEP window (peak-to-peak) | 15–50 ms after the TMS pulse |
| Pre-contraction window / threshold | −50 to 0 ms, 0.15 mV peak-to-peak |
| Baseline RMS window / threshold | −40 to −10 ms, 0.05 mV |
| Dummy trials | first video trial and first control trial of each session |

---

## Output Files

### MATLAB → `derivatives/mep_analysis/` and `derivatives/mep_ztransform/`

- `sub-XX_ses-{AO,AOMI}_mep.tsv/.mat`: trial-level MEPs after exclusions
- `sub-XX_ses-{AO,AOMI}_zmep.tsv/.mat`: trial-level z-scores
- `mean_z_scores_control_vs_main.csv`: control vs. video trials per session (control-trial analysis)
- `mean_z_scores_by_intensity.csv`
- `mean_z_scores_by_intensity_and_weight.csv`: input to `mep_analysis.R`

### `accuracy_analysis.R` → `derivatives/behavioural_analysis/`

- `accuracy_by_condition.xlsx`: accuracy per participant in AO and AOMI (input to `mep_analysis.R`)
- `Fig2_accuracy_by_weight.png`
- `FigS1_accuracy_by_session.png`
- `FigS2_accuracy_by_tms_intensity.png`

### `mep_analysis.R` → `derivatives/ancova_analysis/`

- `Fig3_zMEP_by_intensity.png`
- `Fig4_accuracy_condition.png`
- `FigS3_mep_by_muscle.png`
- `FigS4_residual_qq_slope.png`
- `FigS5_leave_one_out_slope.png`

Tables ST1–ST4 are printed to the console in the format of the Supplementary Material.

---

## Statistical Methods

All models are fitted with `afex::aov_car` (Type III sums of squares). Greenhouse–Geisser correction is applied only to effects for which Mauchly's test is significant (p < 0.05). Effect sizes are partial η² and Cohen's d. In the formulas, `id` is the participant and `accuracy_c` is the mean-centred weight estimation accuracy (mean of AO and AOMI). In `afex`, a covariate is automatically crossed with all within-subject factors, so `accuracy_c + condition * weight` also tests the interactions of accuracy with condition and weight.

**Z-transformation**: one mean and SD per participant, computed from all valid trials pooled across both sessions, all intensities and both trial types (video and control).

**Recruitment curve**: for each participant, condition and weight response, the least-squares slope of the z-scored MEPs on TMS intensity (z per 10% RMT).

| Analysis | Model |
|----------|-------|
| Hypotheses 1 and 2 | `slope ~ condition * weight + Error(id / (condition * weight))` |
| Accuracy (Table ST1) | `slope ~ accuracy_c + condition * weight + Error(id / (condition * weight))` |
| Full zMEP model (Table ST2) | `zMEP ~ accuracy_c + condition * weight * intensity + Error(id / (condition * weight * intensity))` |
| Per intensity (Table ST3) | `zMEP ~ accuracy_c + condition * weight + Error(id / (condition * weight))`, Bonferroni-corrected over five intensities |
| Session order (Table ST4) | `slope ~ order + accuracy_c + condition * weight + Error(id / (condition * weight))` |

Control and robustness analyses: imagery ability (VMIQ-2) as a second covariate, no-video control trials (paired t-test), FDI vs. ADM (muscle × condition × intensity on raw MEPs), Shapiro–Wilk test of the residuals, permutation test (10,000 permutations) and leave-one-participant-out analysis of the accuracy × condition interaction on the slope.

---

## Prerequisites

- **MATLAB** R2019b or later (Statistics and Machine Learning Toolbox)
- **R** 4.0 or later; packages listed in `requirements.R`

---

## Use of AI Tools

The MATLAB and R code was developed with the assistance of an AI tool (Claude, Anthropic), which was used for code generation, restructuring and debugging. All analyses were specified by the authors, and all code and results were reviewed and verified by the authors.

---

## Citation

```bibtex
@article{womi_2026,
  title   = {Individual Task Approach Modulates Corticospinal Excitability During Combined Action Observation and Motor Imagery},
  author  = {Lange, Yassamin and Nesbit, Emma and Hermsdoerfer, Joachim and Stadler, Waltraud},
  journal = {[Journal Name]},
  year    = {2026},
  doi     = {[DOI when available]}
}
```

---

## License

MIT, see `LICENSE.md`.

---

## Contributors

- **Study design**: Yassamin Lange, Emma Nesbit, Waltraud Stadler
- **Data collection**: Yassamin Lange, Emma Nesbit
- **Analysis pipeline**: Emma Nesbit
- **Code review**: Yassamin Lange, Waltraud Stadler

---

## Contact

**Emma Nesbit**: nesbit@cbs.mpg.de
Max Planck Institute for Human Cognitive and Brain Sciences
Use GitHub Issues for bug reports and GitHub Discussions for questions.

---

## Acknowledgments

We thank the Chair of Human Movement Science at TUM for providing facilities and all participants who volunteered their time.

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | April 2026 | Initial release with complete BIDS pipeline |
| 2.0.0 | November 2026 | Revised analyses for resubmission: recruitment-curve slope as main outcome; per-intensity mixed models replaced by Bonferroni-corrected ANCOVAs; Greenhouse–Geisser only where Mauchly's test is significant; added control and robustness analyses; scripts reorganised to follow the Results; figures renamed to match the manuscript |

**Tested with**: MATLAB R2023a · R 4.3.3
