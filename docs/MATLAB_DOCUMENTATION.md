# MATLAB Pipeline Documentation

## Overview

The MATLAB pipeline processes the raw TMS-EMG data in BIDS format in two stages:

1. **MEP extraction and quality control** (`preprocessing/main.m`)
2. **Z-transformation and aggregation** (`analysis/z_transformation.m`)

Run them in this order. The outputs of Stage 2 are the inputs to `mep_analysis.R` (see `R_DOCUMENTATION.md`).

The pipeline processes the FDI muscle only. The EMG files also contain the abductor digiti minimi (ADM) in the column `ADM_in_mV`. To create the file used for the muscle comparison in the Supplementary Material (`ADM_FDI_raw_MEPs_for_ANOVA.csv`), adapt `main.m` to extract the MEPs from both `FDI_in_mV` and `ADM_in_mV` and keep the raw (not z-transformed) MEPs. That file used a baseline RMS threshold of 0.10 mV for both muscles (`config.baselineRMSThreshold = 0.10`), so that trial retention was comparable between the muscles.

---

## Stage 1: MEP Extraction (`preprocessing/main.m`)

### Purpose

Extracts peak-to-peak MEP amplitudes from the raw EMG data, applies the quality control criteria and saves the remaining trials.

### Input Files

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
└── ...
```

**EMG file** (`*_emg.tsv.gz`):
- `time_relative_to_TMS_pulse_in_ms`: time in ms
- `FDI_in_mV`: first dorsal interosseous EMG in mV
- `ADM_in_mV`: abductor digiti minimi EMG in mV (not used by `main.m`)
- `trial`: trial number

**Behavioural file** (`*_beh.tsv`):
- `trial`: trial number
- `event_type`: event category (TMS_pulse_intensity, weight_in_video, etc.)
- `code`: event value

### Configuration

```matlab
config.bidsDirectory = '/data/path/to/bids_directory/';
config.subjects = 1:18;              % participant IDs
config.sessions = {'AO', 'AOMI'};    % sessions

% Analysis windows (s)
config.baselineWindow = [-0.040, -0.010];   % baseline RMS
config.mepWindow      = [0.015, 0.050];     % MEP amplitude
config.preconWindow   = [-0.050, 0];        % pre-contraction check

% Quality control thresholds
config.preconThreshold      = 0.15;         % max pre-contraction (mV, peak-to-peak)
config.baselineRMSThreshold = 0.05;         % max baseline RMS (mV)
config.dummyTrials          = [1, 12];      % first video trial and first control trial
```

### Exclusion Criteria

A trial is excluded if it meets any of these criteria:

1. **Dummy trials**: trials 1 and 12 of each session, i.e. the first video trial and the first control trial
2. **Noisy trials**: trials with artefacts identified by visual inspection (listed in `loadNoisyTrialsList()`)
3. **Pre-contraction**: peak-to-peak EMG > 0.15 mV from 50 ms before the TMS pulse to the pulse
4. **Baseline activity**: RMS EMG > 0.05 mV from 40 to 10 ms before the TMS pulse

A trial can meet more than one criterion, so the counts per criterion add up to more than the total number of excluded trials.

### Processing Steps

#### 1. Data loading
```matlab
[emgData, behData] = loadBIDSData(config, subID, sessionName);
```
- Loads the compressed EMG file and converts time from ms to s
- Checks that the required columns are present and removes NaN values
- Loads the behavioural events
- Removes the final row of the behavioural file, which would otherwise be counted as a 123rd trial

#### 2. MEP calculation
```matlab
results = processMEPs(emgData, behData, config);
```
- **MEP amplitude**: `max(segment) - min(segment)` from 15 to 50 ms after the TMS pulse
- **Baseline RMS**: `rms(segment)` from −40 to −10 ms
- **Pre-contraction**: `max(segment) - min(segment)` from −50 to 0 ms

#### 3. Behavioural data
```matlab
[intensity, weightShown, weightResponse, isControl] = extractBehaviouralData(behData, nTrials);
```
Per trial:
- **TMS intensity**: 1–5 (90, 100, 110, 120, 130% RMT)
- **Weight shown**: weight in the video (1 = light, 2 = moderate, 3 = heavy)
- **Weight response**: the participant's weight estimate (1–3)
- **Is control**: no-video control trial, identified by the event type `Control_TMS_pulse`. All other trials are video (weight estimation) trials.

#### 4. Exclusions
```matlab
[results, exclusionInfo] = applyExclusions(results, config, noisyTrials, subIdx, sessionName);
```
```matlab
excludeMask = dummyMask | noisyMask | preconMask | baselineMask;
results(excludeMask, :) = [];
```

### Output (`derivatives/mep_analysis/`)

- `sub-XX_ses-{AO,AOMI}_mep.tsv` (BIDS derivative)
- `sub-XX_ses-{AO,AOMI}_mep.mat`

| Column | Description | Unit |
|--------|-------------|------|
| MEP_Amplitude | peak-to-peak MEP | mV |
| TMS_Intensity | intensity level | 1–5 |
| Weight_Shown | weight in the video | 1–3 |
| Weight_Response | weight estimate | 1–3 |
| Is_Control | control trial | 0/1 |
| Precontraction | pre-TMS activity | mV |
| Baseline_RMS | baseline activity | mV |

### Noisy Trials List

Hard-coded in `loadNoisyTrialsList()` as `{participant_id, session, [trial numbers]}`, for example:

```matlab
noisyTrials = {
    'womi2_01', 'AO',   [];
    'womi2_01', 'AOMI', [67, 56, 57];
    'womi2_02', 'AO',   [56, 55];
    ...
};
```

### Console Output

```
=== Processing sub-01 ses-AO ===
Trials: 122 initial -> 120 final (2 excluded)
...
=== PROCESSING SUMMARY ===
Initial trials: 4392
Final trials: 4220
Total excluded: 172 (3.9%)

Exclusion breakdown:
  Dummy trials: 72
  Noisy trials: 41
  Precontraction: 19
  High baseline RMS: 61
Analysis complete.
```

---

## Stage 2: Z-Transformation (`analysis/z_transformation.m`)

### Purpose

Standardises the MEPs within each participant, which removes differences in overall MEP size between participants (e.g. due to scalp-to-cortex distance) and keeps the within-participant effects.

### Input

The Stage 1 output in `derivatives/mep_analysis/` (`sub-XX_ses-{AO,AOMI}_mep.tsv`).

### Z-Transformation

For each participant, **one mean and one SD are computed from all valid trials pooled across both sessions (AO and AOMI), all five intensities and both trial types (video and control trials)**:

```matlab
mergedData = [dataAO; dataAOMI];
participantMean = mean(mergedData.MEP_Amplitude);
participantStd  = std(mergedData.MEP_Amplitude);
mergedData.Z_MEP = (mergedData.MEP_Amplitude - participantMean) / participantStd;
```

The merged data are then split back into the two sessions and saved per session.

### Aggregation

| Part | Output | Content |
|------|--------|---------|
| 2 | `mean_z_scores_control_vs_main.csv` | mean zMEP per participant for control and video trials in each session (columns `AO_Control`, `AO_Main`, `AOMI_Control`, `AOMI_Main`); input to the control-trial analysis in `mep_analysis.R` |
| 3 | `mean_z_scores_by_intensity.csv` | mean zMEP per intensity for control and video trials in each session (`AO_Control_1` … `AOMI_Main_5`) |
| 4 | `mean_z_scores_by_intensity_and_weight.csv` | mean zMEP per session × weight response × intensity, plus control trials per intensity; **input to `mep_analysis.R`** |
| 5 | console | paired t-test of control-trial zMEPs, AO vs. AOMI session |

Column names in Part 4: `{Session}_{TrialType}_{Intensity}`, where TrialType is `Control` or `WtResp_{1,2,3}`. Example: `AO_WtResp_1_3` = AO session, light weight response, 100% RMT.

### Output (`derivatives/mep_ztransform/`)

- `sub-XX_ses-{AO,AOMI}_zmep.tsv` and `.mat`: trial-level data with the added `Z_MEP` column
- the three CSV files above

### Console Output

```
=== PART 1: Merging and z-transforming data ===
Participant sub-01: Mean = 0.712, SD = 0.987
...
=== PART 2: Mean z-scores for control vs main trials ===
=== PART 3: Mean z-scores by intensity ===
=== PART 4: Mean z-scores by intensity and subjective weight ===
=== PART 5: Paired t-test - Control trials AO vs AOMI ===
t(17) = X.XXX, p = X.XXX
Mean difference = X.XXX, 95% CI [X.XXX, X.XXX]
=== ALL ANALYSES COMPLETE ===
```

---

## Function Reference

### `main.m`

| Function | Description |
|----------|-------------|
| `initialiseStats()` / `updateStats(stats, info)` | track exclusion counts across participants |
| `loadBIDSData(config, subID, sessionName)` | load and check EMG (`.time`, `.fdi`) and behavioural data |
| `processMEPs(emgData, behData, config)` | MEP amplitude, baseline RMS and pre-contraction per trial |
| `extractBehaviouralData(behData, nTrials)` | trial-level intensity, weight shown, weight response and control flag |
| `parseCode(code)` | convert event codes to numbers (`str2double`) |
| `applyExclusions(results, config, noisyTrials, subIdx, sessionName)` | apply the exclusion criteria; returns the filtered table and the counts |
| `loadNoisyTrialsList()` | list of noisy trials |
| `saveResults(outputDir, subID, sessionName, results)` | save as TSV and MAT |
| `printSummary(stats)` | print the exclusion summary |

### `z_transformation.m`

| Function | Description |
|----------|-------------|
| `loadMEPData(inputDir, subID, sessionName)` | load Stage 1 output |
| `loadZTransformedData(outputDir, subID, sessionName)` | load z-transformed data |
| `saveZTransformedData(outputDir, subID, sessionName, data)` | save as TSV and MAT |

---

## Changing Parameters

Time windows and thresholds are set in `config` at the top of `main.m`, for example:

```matlab
config.mepWindow = [0.015, 0.060];      % default [0.015, 0.050]
config.baselineRMSThreshold = 0.10;     % default 0.05
```

Time values must match the timestamps in the EMG files.

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | April 2026 | Initial BIDS-compatible version |
| 2.0 | November 2026 | Documentation updated for the revised manuscript: dummy trials described correctly; control trials stated as part of the z-transformation; `z_transformation.m` file name; link to the R pipeline and OSF inputs |

---

## Contact

Emma Nesbit: nesbit@cbs.mpg.de · GitHub Issues for bug reports
