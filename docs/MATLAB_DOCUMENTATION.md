# MATLAB Pipeline Documentation

## Overview

The MATLAB pipeline processes raw TMS-EMG data from BIDS format through two main stages:
1. **MEP Extraction & Quality Control** (`main.m`)
2. **Z-Score Transformation & Aggregation** (`z-transformation.m`)

---

## Stage 1: MEP Extraction (`preprocessing/main.m`)

### Purpose
Extracts Motor Evoked Potentials (MEPs) from raw EMG data, applies quality control filters, and saves clean trial-level data.

### Input Files

The script expects BIDS-formatted data:

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

**EMG File Format** (`*_emg.tsv.gz`):
- `time_relative_to_TMS_pulse_in_ms`: Time in milliseconds
- `FDI_in_mV`: First Dorsal Interosseous muscle activity in millivolts

**Behavioural File Format** (`*_beh.tsv`):
- `trial`: Trial number
- `event_type`: Event category (TMS_pulse_intensity, weight_in_video, etc.)
- `code`: Event value

### Configuration Parameters

```matlab
config.bidsDirectory = '/data/path/to/bids_directory/';
config.subjects = 1:18;              % Participant IDs
config.sessions = {'AO', 'AOMI'};    % Session types

% Analysis windows (seconds)
config.baselineWindow = [-0.040, -0.010];   % Pre-TMS baseline
config.mepWindow = [0.015, 0.050];          % MEP measurement window
config.preconWindow = [-0.050, 0];          % Pre-contraction check

% Quality control thresholds
config.preconThreshold = 0.15;              % Max pre-contraction (mV)
config.baselineRMSThreshold = 0.05;         % Max baseline RMS (mV)
config.dummyTrials = [1, 12];               % Trials to exclude
```

### Quality Control Criteria

Trials are excluded if they meet ANY of these criteria:

1. **Dummy Trials**: Trials 1 and 12 (practice/buffer trials)
2. **Noisy Trials**: Manually identified trials with artifacts (see `loadNoisyTrialsList()`)
3. **Pre-contraction**: Peak-to-peak EMG > 0.15 mV in the 50ms before TMS
4. **High Baseline Activity**: RMS EMG > 0.05 mV in the 40-10ms before TMS

### Processing Steps

#### 1. Data Loading
```matlab
[emgData, behData] = loadBIDSData(config, subID, sessionName);
```
- Loads compressed EMG files (`.tsv.gz`)
- Converts time from milliseconds to seconds
- Validates required columns
- Removes NaN values
- Loads behavioural event data
- Removes final row (as safety fix, since it is counted as a 123rd trial otherwise)

#### 2. MEP Calculation
```matlab
results = processMEPs(emgData, behData, config);
```

**MEP Amplitude** (peak-to-peak):
```matlab
for i = 1:nTrials
    segment = emgData.fdi(startIdx:endIdx);
    mepAmplitudes(i) = max(segment) - min(segment);
end
```
- Window: 15-50 ms post-TMS pulse
- Captures both positive and negative deflections

**Baseline RMS** (quality metric):
```matlab
baselineRMS(i) = rms(baselineSegment);
```
- Window: -40 to -10 ms before TMS
- Used for quality control

**Pre-contraction** (muscle activity check):
```matlab
precontraction(i) = max(segment) - min(segment);
```
- Window: -50 to 0 ms before TMS
- Peak-to-peak amplitude

#### 3. Behavioural Data Extraction
```matlab
[intensity, weightShown, weightResponse, isControl] = extractBehaviouralData(behData, nTrials);
```

Extracts per trial:
- **TMS Intensity**: Stimulation level (1-5 scale)
- **Weight Shown**: Object weight in video (1=Light, 2=Moderate, 3=Heavy)
- **Weight Response**: Participant's weight judgment (1-3 scale)
- **Is Control**: Whether trial is a control condition (no weight judgment), detected via presence of event type Control_TMS_pulse in behavioural log. Trials without this event type are treated as main trials (weight judgement trials).

#### 4. Exclusion Application
```matlab
[results, exclusionInfo] = applyExclusions(results, config, noisyTrials, subIdx, sessionName);
```

Creates combined exclusion mask:
```matlab
excludeMask = dummyMask | noisyMask | preconMask | baselineMask;
results(excludeMask, :) = [];
```

### Output Files

Saved to: `/derivatives/mep_analysis/`

**TSV Format** (BIDS derivatives):
```
sub-01_ses-AO_mep.tsv
```

**MAT Format** (MATLAB native):
```
sub-01_ses-AO_mep.mat
```

**Output Table Structure**:
| Column | Description | Units |
|--------|-------------|-------|
| MEP_Amplitude | Peak-to-peak MEP | mV |
| TMS_Intensity | Stimulation level | 1-5 |
| Weight_Shown | Stimulus weight | 1-3 |
| Weight_Response | Participant judgment | 1-3 |
| Is_Control | Control trial flag | boolean |
| Precontraction | Pre-TMS activity | mV |
| Baseline_RMS | Baseline muscle activity | mV |

### Noisy Trials List

Hard-coded in `loadNoisyTrialsList()` function. Example entries:

```matlab
noisyTrials = {
    'womi2_01', 'AO', [];
    'womi2_01', 'AOMI', [67, 56, 57];
    'womi2_02', 'AO', [56, 55];
    ...
};
```

Format: `{participant_id, session, [trial_numbers]}`

### Expected Console Output

```
=== Processing sub-01 ses-AO === 

Trials: 122 initial -> 120 final (2 excluded)

=== Processing sub-01 ses-AOMI ===

Trials: 122 initial -> 117 final (5 excluded)


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

## Stage 2: Z-Score Transformation (`analysis/z-transformation.m`)

### Purpose
Performs within-participant standardisation to remove baseline MEP amplitude differences while preserving experimental effects.

### Input Files

From Stage 1 output:
```
/derivatives/mep_analysis/
├── sub-01_ses-AO_mep.tsv
├── sub-01_ses-AOMI_mep.tsv
├── sub-02_ses-AO_mep.tsv
└── ...
```

### Z-Score Transformation Method

**Within-Participant Standardisation**:

1. Merge both sessions for each participant
2. Calculate participant-specific parameters:
   ```matlab
   participantMean = mean(all_MEPs);
   participantStd = std(all_MEPs);
   ```
3. Transform all MEPs:
   ```matlab
   zMEP = (MEP - participantMean) / participantStd;
   ```

**Why Within-Participant?**
- Removes individual baseline differences (e.g., skull thickness, coil distance)
- Preserves within-participant experimental effects
- Makes participants directly comparable
- Standard approach in TMS-EMG research

### Processing Steps

#### Part 1: Merge and Transform
```matlab
% Load both sessions
dataAO = loadMEPData(config.inputDirectory, subID, 'AO');
dataAOMI = loadMEPData(config.inputDirectory, subID, 'AOMI');

% Merge
mergedData = [dataAO; dataAOMI];

% Calculate z-scores
mepValues = mergedData.MEP_Amplitude;
participantMean = mean(mepValues);
participantStd = std(mepValues);
zMEP = (mepValues - participantMean) / participantStd;

% Add to data
mergedData.Z_MEP = zMEP;

% Split back into sessions
nAO = height(dataAO);
zDataAO = mergedData(1:nAO, :);
zDataAOMI = mergedData(nAO+1:end, :);
```
- Both sessions are loaded, merged, z-transformed using combined parameters, then split back into session-level files. Individual z-transformed files are saved to /derivatives/mep_ztransform/.

#### Part 2: Control vs Main Trial Means
```matlab
resultTable1.AO_Control(k) = mean(dataAO.Z_MEP(dataAO.Is_Control == 1));
resultTable1.AO_Main(k) = mean(dataAO.Z_MEP(dataAO.Is_Control == 0));
```
- Calculates per-participant mean z-scores separately for control and main trials in each session. 

Output: `mean_z_scores_control_vs_main.csv`

| Participant | AO_Control | AO_Main | AOMI_Control | AOMI_Main |
|-------------|------------|---------|--------------|-----------|
| sub-01 | -0.234 | 0.156 | -0.189 | 0.201 |
| ... | ... | ... | ... | ... |

#### Part 3: Mean Z-Scores by Intensity
```matlab
for intensity = 1:5
    controlMask = data.Is_Control == 1 & data.TMS_Intensity == intensity;
    mainMask = data.Is_Control == 0 & data.TMS_Intensity == intensity;
    ...
end
```
- Calculates per-participant mean z-scores for each TMS intensity level (1-5) separately for control and main trials in each session.

Output: `mean_z_scores_by_intensity.csv`

Columns: 
- `AO_Control_1` through `AO_Control_5`
- `AO_Main_1` through `AO_Main_5`
- `AOMI_Control_1` through `AOMI_Control_5`
- `AOMI_Main_1` through `AOMI_Main_5`

#### Part 4: Full Factorial Aggregation
```matlab
for intensity = 1:5
    for weightResp = 1:3
        wtMask = data.Is_Control == 0 & 
                 data.Weight_Response == weightResp & 
                 data.TMS_Intensity == intensity;
        ...
    end
end
```
- Calculates mean z-scores for all combinations of session, trial type, weight response and TMS intensity.

Output: `mean_z_scores_by_intensity_and_weight.csv`

Columns:
- Control conditions (10 columns: 2 sessions × 5 intensities)
- Weight response conditions (30 columns: 2 sessions × 3 weights × 5 intensities)

Format: `{Session}_{TrialType}_{Weight}_{Intensity}`
- Example: `AO_WtResp_1_3` = AO session, Light weight response, 100% intensity

#### Part 5: Paired Samples t-test (AO vs. AOMI control trials)
```matlab
[~, p, ci, stats] = ttest(aoControlMeans(validIdx), aomiControlMeans(validIdx));
meanDiff = mean(aoControlMeans(validIdx) - aomiControlMeans(validIdx));
```
- Tests whether mean z-scored MEP amplitudes for control trials differ significantly between AO and AOMI sessions, using a withing participant paired t-test. NaN values are excluded from the analysis.

Console Output: 
```matlab
t(df) = X.XXX, p = X.XXX
Mean difference = X.XXX, 95% CI [X.XXX, X.XXX]
```

### Output Files

Saved to: `/derivatives/mep_ztransform/`

**Individual Z-Transformed Data**:
```
sub-01_ses-AO_zmep.tsv
sub-01_ses-AO_zmep.mat
sub-01_ses-AOMI_zmep.tsv
sub-01_ses-AOMI_zmep.mat
```

**Aggregated CSV Files**:
```
mean_z_scores_control_vs_main.csv          # Part 2 output
mean_z_scores_by_intensity.csv             # Part 3 output
mean_z_scores_by_intensity_and_weight.csv  # Part 4 output (→ R input)
```

### Expected Console Output

```
=== PART 1: Merging and z-transforming data ===
Participant sub-01: Mean = 0.712, SD = 0.987
Participant sub-02: Mean = 1.101, SD = 1.439
...
Within-participant z-transformation completed.

=== PART 2: Mean z-scores for control vs main trials ===
[Table displayed]
Result table 1 saved to .../mean_z_scores_control_vs_main.csv

=== PART 3: Mean z-scores by intensity ===
[Table displayed]
Result table 2 saved to .../mean_z_scores_by_intensity.csv

=== PART 4: Mean z-scores by intensity and subjective weight ===
[Table displayed]
Result table 3 saved to .../mean_z_scores_by_intensity_and_weight.csv

=== PART 5: Paired t-test - Control trials AO vs AOMI ===
t(17) = X.XXX, p = X.XXX
Mean difference = X.XXX, 95% CI [X.XXX, X.XXX]

=== ALL ANALYSES COMPLETE ===
```

---

## Function Reference

### Main Script Functions

#### `initialiseStats()`
Returns empty statistics structure for tracking exclusions.

#### `updateStats(stats, info)`
Accumulates exclusion counts across participants.

#### `loadBIDSData(config, subID, sessionName)`
Loads and validates EMG and behavioural data files.

**Returns**:
- `emgData`: Structure with `.time` and `.fdi` fields
- `behData`: Table with trial-level behavioural events

#### `processMEPs(emgData, behData, config)`
Calculates MEP amplitudes and quality metrics.

**Returns**: Table with columns:
- MEP_Amplitude, TMS_Intensity, Weight_Shown, Weight_Response
- Is_Control, Precontraction, Baseline_RMS

#### `extractBehaviouralData(behData, nTrials)`
Parses behavioural events into trial-level vectors.

#### `parseCode(code)`
Helper function within extractBehaviouralData(). Converts event code values to numeric format, handling both cell array and already-numeric inputs via str2double.

#### `applyExclusions(results, config, noisyTrials, subIdx, sessionName)`
Applies all quality control filters and removes bad trials.

**Returns**:
- `results`: Filtered data table
- `info`: Structure with exclusion counts

#### `loadNoisyTrialsList()`
Returns hard-coded cell array of manually identified noisy trials.

#### `saveResults(outputDir, subID, sessionName, results)`
Saves results in both TSV (BIDS) and MAT formats.

#### `printSummary(stats)`
Displays final exclusion statistics.

### Z-Transform Script Functions

#### `loadMEPData(inputDir, subID, sessionName)`
Loads processed MEP data from Stage 1 output.

#### `loadZTransformedData(outputDir, subID, sessionName)`
Loads z-transformed MEP data.

#### `saveZTransformedData(outputDir, subID, sessionName, data)`
Saves z-scores in TSV and MAT formats.

---

## Code Modification Guidelines

### Changing Time Windows

To modify analysis windows:

```matlab
% Example: Extend MEP window to 60ms
config.mepWindow = [0.015, 0.060];  % Original: [0.015, 0.050]
```

**Important**: Time values must exactly match EMG file timestamps.

### Adjusting Quality Thresholds

```matlab
% More lenient (keep more trials)
config.preconThreshold = 0.20;        % Original: 0.15
config.baselineRMSThreshold = 0.08;   % Original: 0.05

% More strict (exclude more trials)
config.preconThreshold = 0.10;
config.baselineRMSThreshold = 0.03;
```
---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | April 2026 | Initial BIDS-compatible version |

---

## Contact & Support

For questions about the MATLAB pipeline:
- Emma Nesbit: nesbit@cbs.mpg.de
- GitHub Issues for bug reports
