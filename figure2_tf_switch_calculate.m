% figure2_tf_switch_calculate.m - Calculate TF power differences: Switch vs. Repeat
%
% This script analyzes time-frequency power differences between switch and
% repeat trials, revealing oscillatory correlates of priority switching.
%
% KEY ANALYSES:
%   1. Load time-frequency decomposed data for each subject
%   2. Apply baseline correction (-750 to -250 ms)
%   3. Average power by trial type (Switch, Repeat 1, Repeat 2+)
%   4. Organize by response repetition (same vs. different response)
%
% KEY FINDING:
%   Switch trials show decreased beta (15-25Hz) and increased theta (4-8Hz)
%   power in the cue-probe interval.
%
% INPUTS:
%   - config.m (path configuration)
%   - Subject EEG data (*_EEG_final.mat)
%   - Time-frequency data (*_HanningTF_Laplace.mat)
%   - Artifact rejection files (*_semiautomaticAR.mat)
%
% OUTPUTS:
%   - results/Figure2_TF_SwitchVsRepeat_Results_[timestamp].mat
%     Contains:
%       output.tfmu - TF power by condition
%       output.timeTF, output.freq - Time and frequency axes
%       output.frontalchans, output.centralchans, output.postchans - Channel groups
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

% Subject pairs (participated in both experiments)
subject_pairs = [05 10 11 12 13 14 15; 
                 28 23 22 24 34 42 41];

fprintf('\n=== Figure 2: Time-Frequency Analysis (Switch vs. Repeat) ===\n\n');

%% 2. SUBJECT LOOP - PROCESS EACH SUBJECT

for i_sub = 1:num_subjects
    
    fprintf('Processing Subject %d/%d (ID: S%02d)...\n', i_sub, num_subjects, subject_list(i_sub));
    
    % Determine experiment
    experiment = experiment_list(i_sub);
    data_path = fullfile(paths.data, sprintf('Switch_EEG%d', experiment));
    tf_path = fullfile(paths.derivatives, 'TF', sprintf('Switch_EEG%d', experiment));
    
    % Subject-specific paths
    subject_id = sprintf('S%02d', subject_list(i_sub));
    subject_path = fullfile(data_path, subject_id);
    
    %% 2.1 LOAD DATA
    
    % Load EEG and behavioral data
    filename = fullfile(subject_path, sprintf('%s_EEG_final.mat', subject_id));
    if ~exist(filename, 'file')
        error('File does not exist: %s', filename);
    end
    load(filename);  % loads 'data' and 'behav'
    
    num_trials = length(data.trial);
    
    % Load time-frequency data
    tf_filename = fullfile(tf_path, subject_id, sprintf('%s_HanningTF_Laplace.mat', subject_id));
    load(tf_filename);  % loads 'datatfr'
    
    % Load artifact rejection
    artifact_file = fullfile(subject_path, sprintf('%s_semiautomaticAR.mat', subject_id));
    if exist(artifact_file, 'file')
        load(artifact_file, 'rejsemiautomatic');
    else
        warning('Artifact rejection file not found: %s', artifact_file);
        rejsemiautomatic = [];
    end
    
    %% 2.2 IDENTIFY USABLE TRIALS
    
    reaction_time = behav.time;
    missed_trials = isnan(reaction_time);
    prev_missed = [false; missed_trials(1:end-1)];
    
    % Mark artifacts
    is_artifact = zeros(num_trials, 1);
    is_artifact(rejsemiautomatic) = 1;
    is_artifact(missed_trials) = 1;
    is_artifact(prev_missed) = 1;
    
    usable_trials = ~is_artifact;
    usable_trials = usable_trials & ~isinf(behav.time) & ~isinf([0; behav.time(1:end-1)]);
    
    fprintf('  Usable trials: %d (%.1f%%)\n', ...
        nnz(usable_trials), 100*nnz(usable_trials)/length(usable_trials));
    
    num_usable = nnz(usable_trials);
    
    %% 2.3 EXTRACT BEHAVIORAL VARIABLES
    
    % Response and accuracy
    response = behav.resp;
    prev_response = [0; response(1:end-1)];
    response_repeat = response == prev_response;
    response_repeat = response_repeat(usable_trials);
    
    accuracy = behav.corr(usable_trials);
    block = behav.block(usable_trials);
    block_trial = behav.blocktrial(usable_trials);
    
    % Trial type (1=switch, 2=first repeat, 3+=subsequent repeats)
    trial_type = behav.tnum(usable_trials);
    trial_type(trial_type > 3) = 3;
    
    % Template distance (angular distance between items)
    rule_angles = behav.relrule(usable_trials, :);
    template_distance = circ_dist(rule_angles(:,1)*2*pi/180, rule_angles(:,2)*2*pi/180) * 180/pi/2;
    template_distance = abs(round(template_distance * 1000) / 1000);
    
    %% 2.4 PREPROCESS TIME-FREQUENCY DATA
    
    % Extract TF power for usable trials
    tf_power = datatfr.powspctrm(usable_trials, :, :, :);
    
    time_tf = datatfr.time;
    frequencies = datatfr.freq;
    num_frequencies = length(frequencies);
    num_timepoints = length(time_tf);
    
    % Log-transform power (convert to dB)
    tf_power = 10 * log10(tf_power);
    
    % Baseline correction
    baseline_window = [-0.75, -0.25];  % seconds
    baseline_idx = time_tf >= baseline_window(1) & time_tf <= baseline_window(2);
    baseline_power = mean(tf_power(:, :, :, baseline_idx), 4);
    tf_power = bsxfun(@minus, tf_power, baseline_power);
    
    % Channel information
    channel_labels = datatfr.label;
    
    % Define channel groups (on first subject)
    if i_sub == 1
        output = struct();
        output.timeTF = time_tf;
        output.freq = frequencies;
        output.subpairs = subject_pairs;
        
        output.frontalchans = ismember(channel_labels, ...
            {'AF3', 'AFz', 'AF4', 'F1', 'Fz', 'F2'});
        output.centralchans = ismember(channel_labels, ...
            {'C1', 'Cz', 'C2', 'CP1', 'CPz', 'CP2', 'P1', 'Pz', 'P2'});
        output.postchans = ismember(channel_labels, ...
            {'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'O1', 'Oz', 'O2'});
    end
    
    %% 2.5 AVERAGE POWER BY CONDITION
    
    % Conditions: [trial_type, response_repeat]
    % trial_type: 1=switch, 2=repeat1, 3=repeat2+
    % response_repeat: 0=different response, 1=same response
    conditions = [1 2 3 1 2 3;
                  0 0 0 1 1 1];
    
    num_conditions = size(conditions, 2);
    
    for i_cond = 1:num_conditions
        trials = (trial_type == conditions(1, i_cond)) & ...
                 (response_repeat == conditions(2, i_cond)) & ...
                 (block_trial > 1);  % exclude first trial in block
        
        if nnz(trials) > 0
            tf_mean(i_sub, :, :, :, i_cond) = mean(tf_power(trials, :, :, :), 1);
        else
            tf_mean(i_sub, :, :, :, i_cond) = NaN(1, size(tf_power, 2), num_frequencies, num_timepoints);
        end
    end
    
    % Store results on last subject
    if i_sub == num_subjects
        output.tfmu = tf_mean;
        output.condition_labels = {'Switch, Diff Resp', 'Repeat1, Diff Resp', 'Repeat2+, Diff Resp', ...
                                   'Switch, Same Resp', 'Repeat1, Same Resp', 'Repeat2+, Same Resp'};
    end
    
end

%% 3. SAVE RESULTS

fprintf('\n=== Saving Results ===\n');

results_dir = paths.results;
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

timestamp = datestr(now, 'yyyymmdd-HHMM');
output_filename = fullfile(results_dir, ...
    sprintf('Figure2_TF_SwitchVsRepeat_Results_%s.mat', timestamp));

save(output_filename, 'output', '-v7.3');
fprintf('Results saved to: %s\n', output_filename);

fprintf('\n=== Analysis Complete! ===\n');
fprintf('Next step: Run figure2_tf_switch_plot.m to generate figures.\n');

%% HELPER FUNCTION

function d = circ_dist(alpha, beta)
    % CIRC_DIST - Circular distance between angles
    d = angle(exp(1i * alpha) .* exp(-1i * beta));
end
