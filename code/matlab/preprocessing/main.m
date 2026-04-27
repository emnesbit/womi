%% MEP Analysis Pipeline for BIDS-Formatted TMS-EMG Data
% Processes EMG and behavioural data to extract Motor Evoked Potentials (MEPs)
% Author: Emma Nesbit
% Date: 22.04.2026

function analyse_mep_data_bids()
    %% Configuration
    config = struct();
    config.bidsDirectory = '/data/path/to/bids_directory/';
    config.outputDirectory = fullfile(config.bidsDirectory, 'derivatives', 'mep_analysis');
    config.subjects = 1:18;
    config.sessions = {'AO', 'AOMI'};
    
    % Analysis parameters
    config.baselineWindow = [-0.040, -0.010];  % seconds
    config.mepWindow = [0.015, 0.050];         % seconds
    config.preconWindow = [-0.050, 0];         % seconds
    config.preconThreshold = 0.15;             % mV
    config.baselineRMSThreshold = 0.05;        % mV 
    config.dummyTrials = [1, 12];              % trials to exclude
    
    % Create output directory
    if ~exist(config.outputDirectory, 'dir')
        mkdir(config.outputDirectory);
    end
    
    % Initialise tracking
    stats = initialiseStats();
    noisyTrials = loadNoisyTrialsList();
    
    %% Main processing loop
    for subIdx = config.subjects
        subID = sprintf('sub-%02d', subIdx);
        
        for sesIdx = 1:length(config.sessions)
            sessionName = config.sessions{sesIdx};
            
            fprintf('\n=== Processing %s ses-%s ===\n', subID, sessionName);
            
            try
                % Load data
                [emgData, behData] = loadBIDSData(config, subID, sessionName);
                
                % Process MEPs
                results = processMEPs(emgData, behData, config);
                
                % Apply exclusion criteria
                [results, exclusionInfo] = applyExclusions(results, config, noisyTrials, subIdx, sessionName);
                
                % Update statistics
                stats = updateStats(stats, exclusionInfo);
                
                % Save results
                saveResults(config.outputDirectory, subID, sessionName, results);
                
            catch ME
                fprintf('ERROR processing %s ses-%s: %s\n', subID, sessionName, ME.message);
                fprintf('Error in function: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
                continue;
            end
        end
    end
    
    % Print summary
    printSummary(stats);
end

%% Helper Functions

function stats = initialiseStats()
    stats.initialTrials = 0;
    stats.finalTrials = 0;
    stats.dummyExcluded = 0;
    stats.noisyExcluded = 0;
    stats.preconExcluded = 0;
    stats.baselineExcluded = 0;
end

function stats = updateStats(stats, info)
    stats.initialTrials = stats.initialTrials + info.initial;
    stats.finalTrials = stats.finalTrials + info.final;
    stats.dummyExcluded = stats.dummyExcluded + info.dummy;
    stats.noisyExcluded = stats.noisyExcluded + info.noisy;
    stats.preconExcluded = stats.preconExcluded + info.precon;
    stats.baselineExcluded = stats.baselineExcluded + info.baseline;
end

function [emgData, behData] = loadBIDSData(config, subID, sessionName)
    %% Load EMG data
    emgFile = fullfile(config.bidsDirectory, subID, ['ses-' sessionName], 'emg', ...
        sprintf('%s_ses-%s_womi_emg.tsv.gz', subID, sessionName));
    
    if ~exist(emgFile, 'file')
        error('EMG file not found: %s', emgFile);
    end
    
    emgTable = readtable(emgFile, 'FileType', 'text', 'Delimiter', '\t');
    
    % Validate columns
    requiredCols = {'time_relative_to_TMS_pulse_in_ms', 'FDI_in_mV'};
    if ~all(ismember(requiredCols, emgTable.Properties.VariableNames))
        error('EMG file missing required columns');
    end
    
    % Extract and convert data
    timeMS = emgTable.time_relative_to_TMS_pulse_in_ms;
    fdiMV = emgTable.FDI_in_mV;
    
    % Handle non-numeric data
    if ~isnumeric(timeMS)
        timeMS = str2double(timeMS);
    end
    if ~isnumeric(fdiMV)
        fdiMV = str2double(fdiMV);
    end
    
    % Remove NaN values
    validRows = ~isnan(timeMS) & ~isnan(fdiMV);
    if sum(~validRows) > 0
        warning('%d rows with NaN values removed', sum(~validRows));
    end
    
    emgData.time = timeMS(validRows) / 1000;
    emgData.fdi = fdiMV(validRows);
    
    %% Load behavioural data
    behFile = fullfile(config.bidsDirectory, subID, ['ses-' sessionName], 'beh', ...
        sprintf('%s_ses-%s_womi_beh.tsv', subID, sessionName));
    
    if ~exist(behFile, 'file')
        error('Behavioural file not found: %s', behFile);
    end
    
    behData = readtable(behFile, 'FileType', 'text', 'Delimiter', '\t');
    
    % Remove last row (safety fix)
    if height(behData) > 0
        behData(end, :) = [];
    end
end

function results = processMEPs(emgData, behData, config)
    % Find trial markers
    mepStartIdx = find(emgData.time == config.mepWindow(1));
    mepEndIdx = find(emgData.time == config.mepWindow(2));
    baselineStartIdx = find(emgData.time == config.baselineWindow(1));
    baselineEndIdx = find(emgData.time == config.baselineWindow(2));
    preconStartIdx = find(emgData.time == config.preconWindow(1));
    preconEndIdx = find(emgData.time == config.preconWindow(2));
    
    if isempty(mepStartIdx)
        error('No trial markers found at expected timepoints');
    end
    
    nTrials = length(mepStartIdx);
    
    % Calculate MEP amplitudes (peak-to-peak: 15-50 ms post-TMS)
    mepAmplitudes = zeros(nTrials, 1);
    for i = 1:nTrials
        startIdx = mepStartIdx(i);
        endIdx = mepEndIdx(i);
        segment = emgData.fdi(startIdx:endIdx);
        mepAmplitudes(i) = max(segment) - min(segment);
    end
    
    % Calculate baseline RMS for quality control
    baselineRMS = zeros(nTrials, 1);
    for i = 1:nTrials
        baselineSegment = emgData.fdi(baselineStartIdx(i):baselineEndIdx(i));
        baselineRMS(i) = rms(baselineSegment);
    end
    
    % Calculate precontraction
    precontraction = zeros(nTrials, 1);
    for i = 1:nTrials
        startIdx = preconStartIdx(i);
        endIdx = preconEndIdx(i);
        segment = emgData.fdi(startIdx:endIdx);
        precontraction(i) = max(segment) - min(segment);
    end
    
    % Extract behavioural information
    [intensity, weightShown, weightResponse, isControl] = extractBehaviouralData(behData, nTrials);
    
    % Combine results
    results = table(mepAmplitudes, intensity, weightShown, weightResponse, ...
        isControl, precontraction, baselineRMS, ...
        'VariableNames', {'MEP_Amplitude', 'TMS_Intensity', 'Weight_Shown', ...
        'Weight_Response', 'Is_Control', 'Precontraction', 'Baseline_RMS'});
end

function [intensity, weightShown, weightResponse, isControl] = extractBehaviouralData(behData, nTrials)
    
    % Initialise with defaults
    intensity = zeros(nTrials, 1);
    weightShown = nan(nTrials, 1);
    weightResponse = nan(nTrials, 1);
    isControl = false(nTrials, 1);
    
    uniqueTrials = unique(behData.trial);
    
    for i = 1:length(uniqueTrials)
        trialNum = uniqueTrials(i);
        
        if trialNum > nTrials
            continue;
        end
        
        trialData = behData(behData.trial == trialNum, :);
        
        % Check if control trial
        isControl(trialNum) = any(strcmp(trialData.event_type, 'Control_TMS_pulse'));
        
        % Extract TMS intensity
        intensityRow = trialData(strcmp(trialData.event_type, 'TMS_pulse_intensity'), :);
        if height(intensityRow) > 0
            intensity(trialNum) = parseCode(intensityRow.code(1));
        end
        
        % Extract weight information (only for non-control trials)
        if ~isControl(trialNum)
            weightRow = trialData(strcmp(trialData.event_type, 'weight_in_video'), :);
            if height(weightRow) > 0
                weightShown(trialNum) = parseCode(weightRow.code(1));
            end
            
            responseRow = trialData(strcmp(trialData.event_type, 'weight_estimation_response'), :);
            if height(responseRow) > 0
                weightResponse(trialNum) = parseCode(responseRow.code(1));
            end
        end
    end
end

function value = parseCode(code)
    if iscell(code)
        code = code{1};
    end
    if isnumeric(code)
        value = code;
    else
        value = str2double(code);
    end
end

function [results, info] = applyExclusions(results, config, noisyTrials, subIdx, sessionName)
    info.initial = height(results);
    
    % Dummy trials
    dummyMask = false(height(results), 1);
    dummyMask(config.dummyTrials) = true;
    info.dummy = sum(dummyMask);
    
    % Noisy trials (from manual inspection)
    noisyMask = false(height(results), 1);
    participantName = sprintf('womi2_%02d', subIdx);
    idx = find(strcmp(noisyTrials(:, 1), participantName) & strcmp(noisyTrials(:, 2), sessionName));
    if ~isempty(idx)
        noisyTrialList = noisyTrials{idx, 3};
        if ~isempty(noisyTrialList)
            validTrials = noisyTrialList(noisyTrialList <= height(results));
            noisyMask(validTrials) = true;
        end
    end
    info.noisy = sum(noisyMask);
    
    % Precontraction threshold mask
    preconMask = results.Precontraction > config.preconThreshold;
    info.precon = sum(preconMask);
    
    % Baseline RMS threshold
    baselineMask = results.Baseline_RMS > config.baselineRMSThreshold;
    info.baseline = sum(baselineMask);
    
    % Combine exclusions
    excludeMask = dummyMask | noisyMask | preconMask | baselineMask;
    
    % Apply exclusions
    results(excludeMask, :) = [];
    info.final = height(results);
    
    fprintf('  Trials: %d initial -> %d final (%d excluded)\n', ...
        info.initial, info.final, info.initial - info.final);
end

function noisyTrials = loadNoisyTrialsList()
    noisyTrials = {
        'womi2_01', 'AO', [];
        'womi2_01', 'AOMI', [67, 56, 57];
        'womi2_02', 'AO', [56, 55];
        'womi2_02', 'AOMI', [119, 120];
        'womi2_03', 'AO', [58, 86, 113, 114];
        'womi2_03', 'AOMI', [];
        'womi2_04', 'AO', [];
        'womi2_04', 'AOMI', [];
        'womi2_05', 'AO', [];
        'womi2_05', 'AOMI', [];
        'womi2_06', 'AO', [];
        'womi2_06', 'AOMI', [];
        'womi2_07', 'AO', [29, 31, 42];
        'womi2_07', 'AOMI', [];
        'womi2_08', 'AO', [];
        'womi2_08', 'AOMI', [];
        'womi2_09', 'AO', [59];
        'womi2_09', 'AOMI', [];
        'womi2_10', 'AO', [58, 78, 81];
        'womi2_10', 'AOMI', [];
        'womi2_11', 'AO', [];
        'womi2_11', 'AOMI', [113];
        'womi2_12', 'AO', [];
        'womi2_12', 'AOMI', [];
        'womi2_13', 'AO', [90];
        'womi2_13', 'AOMI', [];
        'womi2_14', 'AO', [];
        'womi2_14', 'AOMI', [112];
        'womi2_15', 'AO', [29];
        'womi2_15', 'AOMI', [57];
        'womi2_16', 'AO', [];
        'womi2_16', 'AOMI', [];
        'womi2_17', 'AO', [14, 44, 45, 66, 67];
        'womi2_17', 'AOMI', [44, 67, 68, 69, 106, 112, 113, 114];
        'womi2_18', 'AO', [15, 92, 107, 110];
        'womi2_18', 'AOMI', [79]
    };
end

function saveResults(outputDir, subID, sessionName, results)
    % Save as TSV (BIDS derivatives format)
    filename = sprintf('%s_ses-%s_mep.tsv', subID, sessionName);
    writetable(results, fullfile(outputDir, filename), 'FileType', 'text', 'Delimiter', '\t');
    
    % Also save as .mat for MATLAB users
    filename = sprintf('%s_ses-%s_mep.mat', subID, sessionName);
    save(fullfile(outputDir, filename), 'results');
end

function printSummary(stats)
    fprintf('\n=== PROCESSING SUMMARY ===\n');
    fprintf('Initial trials: %d\n', stats.initialTrials);
    fprintf('Final trials: %d\n', stats.finalTrials);
    totalExcluded = stats.initialTrials - stats.finalTrials;
    percentExcluded = 100 * totalExcluded / stats.initialTrials;
    fprintf('Total excluded: %d (%.1f%%)\n\n', totalExcluded, percentExcluded);
    
    fprintf('Exclusion breakdown:\n');
    fprintf('  Dummy trials: %d\n', stats.dummyExcluded);
    fprintf('  Noisy trials: %d\n', stats.noisyExcluded);
    fprintf('  Precontraction: %d\n', stats.preconExcluded);
    fprintf('  High baseline RMS: %d\n', stats.baselineExcluded);
    fprintf('\nAnalysis complete.\n');
end