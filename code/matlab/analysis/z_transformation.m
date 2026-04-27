%% Z-Transformation Analysis for BIDS-Formatted MEP Data
% Performs within-participant z-transformation and creates summary tables
% Author: Emma Nesbit
% Date: 22.04.2026

function analyse_mep_ztransform_bids()
    %% Configuration
    config = struct();
    config.bidsDirectory = '/data/path/to/bids_directory/';
    config.inputDirectory = fullfile(config.bidsDirectory, 'derivatives', 'mep_analysis');
    config.outputDirectory = fullfile(config.bidsDirectory, 'derivatives', 'mep_ztransform');
    config.subjects = 1:18;
    config.sessions = {'AO', 'AOMI'};
    
    % Create output directory
    if ~exist(config.outputDirectory, 'dir')
        mkdir(config.outputDirectory);
    end
    
    participants = cell(length(config.subjects), 1);
    for i = 1:length(config.subjects)
        participants{i} = sprintf('sub-%02d', config.subjects(i));
    end
    
    %% Part 1: Merge sessions and apply within-participant z-transformation
    fprintf('=== PART 1: Merging and z-transforming data ===\n');
    
    for k = 1:length(participants)
        subID = participants{k};
        
        try
            % Load both sessions
            dataAO = loadMEPData(config.inputDirectory, subID, 'AO');
            dataAOMI = loadMEPData(config.inputDirectory, subID, 'AOMI');
            
            % Merge sessions
            mergedData = [dataAO; dataAOMI];
            
            % Calculate within-participant z-scores
            mepValues = mergedData.MEP_Amplitude;
            participantMean = mean(mepValues);
            participantStd = std(mepValues);
            zMEP = (mepValues - participantMean) / participantStd;
            
            % Add z-scores to data
            mergedData.Z_MEP = zMEP;
            
            % Split back into sessions
            nAO = height(dataAO);
            zDataAO = mergedData(1:nAO, :);
            zDataAOMI = mergedData(nAO+1:end, :);
            
            % Save z-transformed data
            saveZTransformedData(config.outputDirectory, subID, 'AO', zDataAO);
            saveZTransformedData(config.outputDirectory, subID, 'AOMI', zDataAOMI);
            
            fprintf('Participant %s: Mean = %.3f, SD = %.3f\n', subID, participantMean, participantStd);
            
        catch ME
            fprintf('ERROR processing %s: %s\n', subID, ME.message);
            continue;
        end
    end
    
    fprintf('Within-participant z-transformation completed.\n\n');
    
    %% Part 2: Calculate mean z-scores for control vs main trials
    fprintf('=== PART 2: Mean z-scores for control vs main trials ===\n');
    
    resultTable1 = table();
    resultTable1.Participant = participants;
    resultTable1.AO_Control = zeros(length(participants), 1);
    resultTable1.AO_Main = zeros(length(participants), 1);
    resultTable1.AOMI_Control = zeros(length(participants), 1);
    resultTable1.AOMI_Main = zeros(length(participants), 1);
    
    for k = 1:length(participants)
        subID = participants{k};
        
        % Process AO session
        try
            dataAO = loadZTransformedData(config.outputDirectory, subID, 'AO');
            resultTable1.AO_Control(k) = mean(dataAO.Z_MEP(dataAO.Is_Control == 1));
            resultTable1.AO_Main(k) = mean(dataAO.Z_MEP(dataAO.Is_Control == 0));
        catch
            resultTable1.AO_Control(k) = NaN;
            resultTable1.AO_Main(k) = NaN;
        end
        
        % Process AOMI session
        try
            dataAOMI = loadZTransformedData(config.outputDirectory, subID, 'AOMI');
            resultTable1.AOMI_Control(k) = mean(dataAOMI.Z_MEP(dataAOMI.Is_Control == 1));
            resultTable1.AOMI_Main(k) = mean(dataAOMI.Z_MEP(dataAOMI.Is_Control == 0));
        catch
            resultTable1.AOMI_Control(k) = NaN;
            resultTable1.AOMI_Main(k) = NaN;
        end
    end
    
    disp(resultTable1);
    csvFile1 = fullfile(config.outputDirectory, 'mean_z_scores_control_vs_main.csv');
    writetable(resultTable1, csvFile1);
    fprintf('Result table 1 saved to %s\n\n', csvFile1);
    
    %% Part 3: Calculate mean z-scores by intensity
    fprintf('=== PART 3: Mean z-scores by intensity ===\n');
    
    % Create column names
    colNames = {'Participant'};
    for session = {'AO', 'AOMI'}
        for trialType = {'Control', 'Main'}
            for intensity = 1:5
                colNames{end+1} = sprintf('%s_%s_%d', session{1}, trialType{1}, intensity);
            end
        end
    end
    
    resultTable2 = array2table(zeros(length(participants), length(colNames)-1), ...
        'VariableNames', colNames(2:end));
    resultTable2.Participant = participants;
    resultTable2 = resultTable2(:, [end, 1:end-1]);
    
    for k = 1:length(participants)
        subID = participants{k};
        
        % Process AO session
        try
            dataAO = loadZTransformedData(config.outputDirectory, subID, 'AO');
            for intensity = 1:5
                controlMask = dataAO.Is_Control == 1 & dataAO.TMS_Intensity == intensity;
                mainMask = dataAO.Is_Control == 0 & dataAO.TMS_Intensity == intensity;
                
                if sum(controlMask) > 0
                    resultTable2.(sprintf('AO_Control_%d', intensity))(k) = mean(dataAO.Z_MEP(controlMask));
                else
                    resultTable2.(sprintf('AO_Control_%d', intensity))(k) = NaN;
                end
                
                if sum(mainMask) > 0
                    resultTable2.(sprintf('AO_Main_%d', intensity))(k) = mean(dataAO.Z_MEP(mainMask));
                else
                    resultTable2.(sprintf('AO_Main_%d', intensity))(k) = NaN;
                end
            end
        catch
            for intensity = 1:5
                resultTable2.(sprintf('AO_Control_%d', intensity))(k) = NaN;
                resultTable2.(sprintf('AO_Main_%d', intensity))(k) = NaN;
            end
        end
        
        % Process AOMI session
        try
            dataAOMI = loadZTransformedData(config.outputDirectory, subID, 'AOMI');
            for intensity = 1:5
                controlMask = dataAOMI.Is_Control == 1 & dataAOMI.TMS_Intensity == intensity;
                mainMask = dataAOMI.Is_Control == 0 & dataAOMI.TMS_Intensity == intensity;
                
                if sum(controlMask) > 0
                    resultTable2.(sprintf('AOMI_Control_%d', intensity))(k) = mean(dataAOMI.Z_MEP(controlMask));
                else
                    resultTable2.(sprintf('AOMI_Control_%d', intensity))(k) = NaN;
                end
                
                if sum(mainMask) > 0
                    resultTable2.(sprintf('AOMI_Main_%d', intensity))(k) = mean(dataAOMI.Z_MEP(mainMask));
                else
                    resultTable2.(sprintf('AOMI_Main_%d', intensity))(k) = NaN;
                end
            end
        catch
            for intensity = 1:5
                resultTable2.(sprintf('AOMI_Control_%d', intensity))(k) = NaN;
                resultTable2.(sprintf('AOMI_Main_%d', intensity))(k) = NaN;
            end
        end
    end
    
    disp(resultTable2);
    csvFile2 = fullfile(config.outputDirectory, 'mean_z_scores_by_intensity.csv');
    writetable(resultTable2, csvFile2);
    fprintf('Result table 2 saved to %s\n\n', csvFile2);
    
    %% Part 4: Calculate mean z-scores by intensity and subjective weight
    fprintf('=== PART 4: Mean z-scores by intensity and subjective weight ===\n');
    
    % Create column names
    colNames = {'Participant'};
    for session = {'AO', 'AOMI'}
        for trialType = {'Control', 'WtResp_1', 'WtResp_2', 'WtResp_3'}
            for intensity = 1:5
                colNames{end+1} = sprintf('%s_%s_%d', session{1}, trialType{1}, intensity);
            end
        end
    end
    
    resultTable3 = array2table(zeros(length(participants), length(colNames)-1), ...
        'VariableNames', colNames(2:end));
    resultTable3.Participant = participants;
    resultTable3 = resultTable3(:, [end, 1:end-1]);
    
    for k = 1:length(participants)
        subID = participants{k};
        
        % Process AO session
        try
            dataAO = loadZTransformedData(config.outputDirectory, subID, 'AO');
            
            for intensity = 1:5
                % Control trials
                controlMask = dataAO.Is_Control == 1 & dataAO.TMS_Intensity == intensity;
                if sum(controlMask) > 0
                    resultTable3.(sprintf('AO_Control_%d', intensity))(k) = mean(dataAO.Z_MEP(controlMask));
                else
                    resultTable3.(sprintf('AO_Control_%d', intensity))(k) = NaN;
                end
                
                % Weight response trials
                for weightResp = 1:3
                    wtMask = dataAO.Is_Control == 0 & dataAO.Weight_Response == weightResp & dataAO.TMS_Intensity == intensity;
                    if sum(wtMask) > 0
                        resultTable3.(sprintf('AO_WtResp_%d_%d', weightResp, intensity))(k) = mean(dataAO.Z_MEP(wtMask));
                    else
                        resultTable3.(sprintf('AO_WtResp_%d_%d', weightResp, intensity))(k) = NaN;
                    end
                end
            end
        catch
            for intensity = 1:5
                resultTable3.(sprintf('AO_Control_%d', intensity))(k) = NaN;
                for weightResp = 1:3
                    resultTable3.(sprintf('AO_WtResp_%d_%d', weightResp, intensity))(k) = NaN;
                end
            end
        end
        
        % Process AOMI session
        try
            dataAOMI = loadZTransformedData(config.outputDirectory, subID, 'AOMI');
            
            for intensity = 1:5
                % Control trials
                controlMask = dataAOMI.Is_Control == 1 & dataAOMI.TMS_Intensity == intensity;
                if sum(controlMask) > 0
                    resultTable3.(sprintf('AOMI_Control_%d', intensity))(k) = mean(dataAOMI.Z_MEP(controlMask));
                else
                    resultTable3.(sprintf('AOMI_Control_%d', intensity))(k) = NaN;
                end
                
                % Weight response trials
                for weightResp = 1:3
                    wtMask = dataAOMI.Is_Control == 0 & dataAOMI.Weight_Response == weightResp & dataAOMI.TMS_Intensity == intensity;
                    if sum(wtMask) > 0
                        resultTable3.(sprintf('AOMI_WtResp_%d_%d', weightResp, intensity))(k) = mean(dataAOMI.Z_MEP(wtMask));
                    else
                        resultTable3.(sprintf('AOMI_WtResp_%d_%d', weightResp, intensity))(k) = NaN;
                    end
                end
            end
        catch
            for intensity = 1:5
                resultTable3.(sprintf('AOMI_Control_%d', intensity))(k) = NaN;
                for weightResp = 1:3
                    resultTable3.(sprintf('AOMI_WtResp_%d_%d', weightResp, intensity))(k) = NaN;
                end
            end
        end
    end
    
    disp(resultTable3);
    csvFile3 = fullfile(config.outputDirectory, 'mean_z_scores_by_intensity_and_weight.csv');
    writetable(resultTable3, csvFile3);
    fprintf('Result table 3 saved to %s\n\n', csvFile3);
    
    %% Part 5: Paired samples t-test for control trials (AO vs AOMI)
    fprintf('=== PART 5: Paired t-test - Control trials AO vs AOMI ===\n');
    
    % Extract and filter valid data
    aoControlMeans = resultTable1.AO_Control;
    aomiControlMeans = resultTable1.AOMI_Control;
    validIdx = ~isnan(aoControlMeans) & ~isnan(aomiControlMeans);
    
    % Perform paired t-test
    [~, p, ci, stats] = ttest(aoControlMeans(validIdx), aomiControlMeans(validIdx));
    meanDiff = mean(aoControlMeans(validIdx) - aomiControlMeans(validIdx));
    
    % Display results
    fprintf('t(%d) = %.3f, p = %.3f\n', stats.df, stats.tstat, p);
    fprintf('Mean difference = %.3f, 95%% CI [%.3f, %.3f]\n\n', meanDiff, ci(1), ci(2));
    
    %% Summary
    fprintf('=== ALL ANALYSES COMPLETE ===\n');
    fprintf('Generated files:\n');
    fprintf('  1. %s\n', csvFile1);
    fprintf('  2. %s\n', csvFile2);
    fprintf('  3. %s\n', csvFile3);
end

%% Helper Functions

function data = loadMEPData(inputDir, subID, sessionName)
    % Load MEP data from TSV file
    filename = sprintf('%s_ses-%s_mep.tsv', subID, sessionName);
    filepath = fullfile(inputDir, filename);
    
    if ~exist(filepath, 'file')
        error('File not found: %s', filepath);
    end
    
    data = readtable(filepath, 'FileType', 'text', 'Delimiter', '\t');
end

function data = loadZTransformedData(outputDir, subID, sessionName)
    % Load z-transformed data
    filename = sprintf('%s_ses-%s_zmep.tsv', subID, sessionName);
    filepath = fullfile(outputDir, filename);
    
    if ~exist(filepath, 'file')
        error('File not found: %s', filepath);
    end
    
    data = readtable(filepath, 'FileType', 'text', 'Delimiter', '\t');
end

function saveZTransformedData(outputDir, subID, sessionName, data)
    % Save as TSV
    filename = sprintf('%s_ses-%s_zmep.tsv', subID, sessionName);
    filepath = fullfile(outputDir, filename);
    writetable(data, filepath, 'FileType', 'text', 'Delimiter', '\t');
    
    % Also save as MAT
    filename = sprintf('%s_ses-%s_zmep.mat', subID, sessionName);
    filepath = fullfile(outputDir, filename);
    save(filepath, 'data');
end