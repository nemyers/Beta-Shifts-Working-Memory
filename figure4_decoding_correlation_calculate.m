% figure4_decoding_correlation_calculate.m - Decode memory representations and correlate with beta power
%
% This script performs inverted encoding model (IEM) decoding to reconstruct
% memory representations, then correlates decoding fidelity with beta power.
%
% KEY ANALYSES:
%   1. Decode active and latent memory representations from ERPs
%   2. Calculate decoding fidelity (reconstruction quality)
%   3. Extract beta power (15-25Hz) at central channels
%   4. Correlate beta power with decoding fidelity across trials
%
% KEY FINDING:
%   Lower beta power predicts stronger decoding of the newly prioritized item.
%   Beta desynchronization tracks the quality of the active memory representation.
%
% INPUTS:
%   - config.m (path configuration)
%   - Subject EEG data (*_EEG_final.mat)
%   - Time-frequency data (*_HanningTF_Laplace.mat)
%   - Behavioral data (Study12_TFPow_Beh_Decoding.mat)
%
% OUTPUTS:
%   - results/Figure4_Decoding_BetaCorrelation_Results_[timestamp].mat
%     Contains:
%       output.decoding_fidelity - IEM decoding quality by trial
%       output.beta_power - Beta power by trial
%       output.correlations - Trial-wise correlations
%
% Author: Nicholas E. Myers
% Cleaned for GitHub: 2026
% For: Myers, Stokes, & Muhle-Karbe (2026) J. Neuroscience

clearvars;
close all;
clc;

%% 1. SETUP PATHS AND PARAMETERS

% Load configuration
paths = config();

% Subject information
subject_list = [2 4 5 7 11 14 15 16 17 18 19 21 22 23 27 29 31 32 35 36 2:31];
experiment_list = [ones(1,20)*1 ones(1,30)*2];
num_subjects = length(subject_list);

% Subject pairs
subject_pairs = [05 10 11 12 13 14 15; 
                 28 23 22 24 34 42 41];

fprintf('\n=== Figure 4: Decoding-Beta Correlation Analysis ===\n\n');

% Load behavioral data
behavioral_data_file = fullfile(paths.results, 'Study12_TFPow_Beh_Decoding.mat');
load(behavioral_data_file);  % loads 'tfdata'

%% 2. DECODING AND BETA ANALYSIS PARAMETERS

% Decoding parameters
decoding_channels = 45:61;  % Posterior channels for decoding
time_windows = 0.0:0.05:1.4;  % Time points for decoding (seconds)
window_length = [-0.450, 0];  % Temporal window for each decoding timepoint
temporal_step = 3;  % Downsampling step within window

% Beta power parameters
beta_freq_band = [15, 25];  % Hz
beta_time_window = [0.4, 0.8];  % seconds

% Baseline parameters
baseline_window = [-0.75, -0.25];  % seconds
erp_baseline_window = [-0.200, -0.050];  % seconds for ERP

%% 3. SUBJECT LOOP

for i_sub = 1:num_subjects
    
    fprintf('\nProcessing Subject %d/%d (ID: S%02d)...\n', ...
        i_sub, num_subjects, subject_list(i_sub));
    
    % Determine experiment
    experiment = experiment_list(i_sub);
    data_path = fullfile(paths.data, sprintf('Switch_EEG%d', experiment));
    tf_path = fullfile(paths.derivatives, 'TF', sprintf('Switch_EEG%d', experiment));
    
    % Subject-specific paths
    subject_id = sprintf('S%02d', subject_list(i_sub));
    subject_path = fullfile(data_path, subject_id);
    
    %% 3.1 LOAD DATA
    
    % Load EEG and behavioral data
    filename = fullfile(subject_path, sprintf('%s_EEG_final.mat', subject_id));
    if ~exist(filename, 'file')
        error('File does not exist: %s', filename);
    end
    load(filename);  % loads 'data' and 'behav'
    
    % Load time-frequency data
    tf_filename = fullfile(tf_path, subject_id, sprintf('%s_HanningTF_Laplace.mat', subject_id));
    load(tf_filename);  % loads 'datatfr'
    
    % Load artifact rejection
    artifact_file = fullfile(subject_path, sprintf('%s_semiautomaticAR.mat', subject_id));
    if exist(artifact_file, 'file')
        load(artifact_file, 'rejsemiautomatic');
    else
        rejsemiautomatic = [];
    end
    
    %% 3.2 IDENTIFY USABLE TRIALS
    
    num_trials = length(data.trial);
    reaction_time = behav.time;
    missed_trials = isnan(reaction_time);
    prev_missed = [false; missed_trials(1:end-1)];
    
    is_artifact = zeros(num_trials, 1);
    is_artifact(rejsemiautomatic) = 1;
    is_artifact(missed_trials) = 1;
    is_artifact(prev_missed) = 1;
    
    usable_trials = ~is_artifact & ~isinf(behav.time) & ~isinf([0; behav.time(1:end-1)]);
    
    fprintf('  Usable trials: %d (%.1f%%)\n', ...
        nnz(usable_trials), 100*nnz(usable_trials)/length(usable_trials));
    
    %% 3.3 EXTRACT BEHAVIORAL VARIABLES
    
    % Memory item orientations (to be decoded)
    % Column 1: active (cued) item
    % Column 2: latent (uncued) item
    item_angles = behav.relrule(usable_trials, :) * 2*pi/180;  % Convert to radians
    
    % Trial type
    trial_type = behav.tnum(usable_trials);
    trial_type(trial_type > 3) = 3;
    
    block = behav.block(usable_trials);
    block_trial = behav.blocktrial(usable_trials);
    
    %% 3.4 PREPROCESS TIME-FREQUENCY DATA FOR BETA POWER
    
    tf_power = datatfr.powspctrm(usable_trials, :, :, :);
    time_tf = datatfr.time;
    frequencies = datatfr.freq;
    
    % Log-transform
    tf_power = 10 * log10(tf_power);
    
    % Baseline correction
    baseline_idx = time_tf >= baseline_window(1) & time_tf <= baseline_window(2);
    baseline_power = mean(tf_power(:, :, :, baseline_idx), 4);
    tf_power = bsxfun(@minus, tf_power, baseline_power);
    
    % Extract beta power at central channels
    central_channels = ismember(datatfr.label, ...
        {'C1', 'Cz', 'C2', 'CP1', 'CPz', 'CP2', 'P1', 'Pz', 'P2'});
    
    beta_freq_idx = frequencies >= beta_freq_band(1) & frequencies <= beta_freq_band(2);
    beta_time_idx = time_tf >= beta_time_window(1) & time_tf <= beta_time_window(2);
    
    % Average beta power for each trial
    beta_power_by_trial = squeeze(mean(mean(mean(...
        tf_power(:, central_channels, beta_freq_idx, beta_time_idx), 2), 3), 4));
    
    %% 3.5 PREPROCESS ERP DATA FOR DECODING
    
    % Concatenate all trials
    erp_data = permute(cat(3, data.trial{usable_trials}), [3, 1, 2]);
    time_erp = data.time{1};
    sampling_rate = data.fsample;
    
    % Apply Gaussian smoothing
    smoothing_duration = 0.024;  % seconds
    smoothing_size = smoothing_duration * sampling_rate;
    erp_data = filtfast(erp_data, 3, [], 'gaussian', smoothing_size);
    
    % Re-epoch and baseline correct
    epoch_window = [-0.600, +2.200];  % seconds
    baseline_window_erp = [-0.200, -0.050];  % seconds
    
    epoch_samples = fix(epoch_window * sampling_rate);
    baseline_samples = fix(baseline_window_erp * sampling_rate);
    
    time_epoched = epoch_samples(1):epoch_samples(2);
    time_epoched = time_epoched / sampling_rate;
    num_timepoints = length(time_epoched);
    
    num_usable_trials = nnz(usable_trials);
    erp_epoched = nan(num_usable_trials, length(decoding_channels), num_timepoints);
    
    for i_trial = 1:num_usable_trials
        % Find cue onset
        cue_onset = find(data.time{i_trial} <= 0, 1, 'last');
        
        % Extract epoch
        time_indices = cue_onset + epoch_samples(1:end);
        time_indices = min(time_indices, length(data.time{i_trial}));
        
        erp_epoched(i_trial, :, :) = erp_data(i_trial, decoding_channels, time_indices);
        
        % Baseline correction
        baseline_indices = cue_onset + baseline_samples(1:end);
        baseline_mean = mean(erp_data(i_trial, decoding_channels, baseline_indices), 3);
        erp_epoched(i_trial, :, :) = squeeze(erp_epoched(i_trial, :, :)) - baseline_mean';
    end
    
    %% 3.6 PERFORM INVERTED ENCODING MODEL (IEM) DECODING
    
    fprintf('  Running IEM decoding...\n');
    
    % Define angle basis functions
    unique_angles = unique(round(item_angles(:,1) * 1000) / 1000);
    num_basis_functions = length(unique_angles) - 1;
    
    % Temporal window parameters
    window_samples = fix(window_length * sampling_rate);
    window_times = (window_samples(1):temporal_step:window_samples(2)) / sampling_rate;
    num_steps = length(window_times) - 1;
    
    % Cross-validation folds (based on blocks)
    fold_labels = ceil(block / 10);
    
    % Preallocate outputs
    num_time_windows = length(time_windows);
    decoding_fidelity = nan(num_usable_trials, num_time_windows, 2, 2);  % trials x time x train_ang x test_ang
    
    % Decode active and latent items
    for i_train_angle = 1:2  % 1=active, 2=latent
        for i_test_angle = 1:2
            
            fprintf('    Training angle %d, Testing angle %d...\n', ...
                i_train_angle, i_test_angle);
            
            % Select trials (exclude first trial in block)
            valid_trials = block_trial > 1;
            
            train_angles = item_angles(valid_trials, i_train_angle);
            test_angles = item_angles(valid_trials, i_test_angle);
            
            % Loop over time windows
            for i_window = 1:num_time_windows
                
                % Find time window onset
                [~, onset_idx] = min(abs(time_epoched - time_windows(i_window)));
                
                % Extract data in temporal window
                time_idx = onset_idx + window_samples(1:temporal_step:end);
                window_data = nan(nnz(valid_trials), length(decoding_channels), num_steps);
                
                for i_step = 1:num_steps
                    step_data = erp_epoched(valid_trials, :, time_idx(i_step):(time_idx(i_step) + temporal_step - 1));
                    window_data(:, :, i_step) = mean(step_data, 3);
                end
                
                % Demean across time
                window_data = bsxfun(@minus, window_data, mean(window_data, 3));
                
                % Reshape for decoding
                window_data = reshape(window_data, [nnz(valid_trials), ...
                    length(decoding_channels) * num_steps]);
                
                % Run IEM decoding with cross-validation
                [fidelity, ~] = mahal_iem_trialwise(window_data, train_angles, ...
                    struct('thetabasis', unique_angles', ...
                          'foldvar', fold_labels(valid_trials), ...
                          'testangles', test_angles));
                
                % Store fidelity (centered to remove mean)
                fidelity = fidelity - mean(fidelity);
                decoding_fidelity(valid_trials, i_window, i_train_angle, i_test_angle) = fidelity;
            end
        end
    end
    
    %% 3.7 CALCULATE TRIAL-WISE CORRELATIONS
    
    fprintf('  Calculating trial-wise correlations...\n');
    
    % Focus on switch and first repeat trials
    trial_types_to_analyze = [1, 2];  % Switch and Repeat1
    
    for i_trial_type = 1:length(trial_types_to_analyze)
        trials = trial_type == trial_types_to_analyze(i_trial_type) & block_trial > 1;
        
        % Average decoding fidelity across time windows
        mean_decoding(i_sub, :, :, i_trial_type) = squeeze(mean(...
            decoding_fidelity(trials, :, :, :), 1));
    end
    
    % Store outputs
    if i_sub == 1
        output = struct();
        output.time_windows = time_windows;
        output.beta_freq_band = beta_freq_band;
        output.beta_time_window = beta_time_window;
    end
    
    output.decoding(i_sub, :, :, :, :) = decoding_fidelity;
    output.beta_power(i_sub, :) = beta_power_by_trial;
    output.trial_type(i_sub, :) = trial_type;
    output.block_trial(i_sub, :) = block_trial;
    
end

%% 4. CALCULATE GROUP-LEVEL CORRELATIONS

fprintf('\n=== Calculating group-level correlations ===\n');

% For each subject, correlate beta power with decoding fidelity
for i_sub = 1:num_subjects
    for i_time = 1:num_time_windows
        for i_train = 1:2
            for i_test = 1:2
                % Get trials for this subject
                trials = output.trial_type(i_sub, :) == 1 & output.block_trial(i_sub, :) > 1;
                
                % Extract beta and decoding for these trials
                beta = output.beta_power(i_sub, trials);
                decoding = squeeze(output.decoding(i_sub, trials, i_time, i_train, i_test));
                
                % Calculate correlation
                valid = ~isnan(beta) & ~isnan(decoding);
                if nnz(valid) > 10
                    output.correlation(i_sub, i_time, i_train, i_test) = ...
                        corr(beta(valid)', decoding(valid));
                else
                    output.correlation(i_sub, i_time, i_train, i_test) = NaN;
                end
            end
        end
    end
end

%% 5. SAVE RESULTS

fprintf('\n=== Saving Results ===\n');

results_dir = paths.results;
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

timestamp = datestr(now, 'yyyymmdd-HHMM');
output_filename = fullfile(results_dir, ...
    sprintf('Figure4_Decoding_BetaCorrelation_Results_%s.mat', timestamp));

save(output_filename, 'output', '-v7.3');
fprintf('Results saved to: %s\n', output_filename);

fprintf('\n=== Analysis Complete! ===\n');
fprintf('Next step: Run figure4_decoding_correlation_plot.m to generate figures.\n');
