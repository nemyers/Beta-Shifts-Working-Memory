% figure1_behavior_calculate.m - Extract and organize behavioral data for analysis
%
% This script loads behavioral data from all subjects' EEG files and
% organizes it into a structure for subsequent analysis and plotting.
% The behavioral data includes response times, accuracy, decision variables,
% template distances, and trial types (switch vs. repeat).
%
% KEY ANALYSES:
%   1. Load behavioral data from preprocessed EEG files
%   2. Organize by subject, experiment, and trial type
%   3. Calculate derived variables (template distance, congruency, etc.)
%
% INPUTS:
%   - config.m (path configuration)
%   - Subject EEG data files (*_EEG_final.mat)
%     Each contains:
%       - data: EEG time series
%       - behav: behavioral data structure
%
% OUTPUTS:
%   - results/Study12_TFPow_Beh_Decoding.mat
%     Contains:
%       tfdata(isub).behaviour - behavioral data structure
%       tfdata(isub).experiment - which experiment (1 or 2)
%       tfdata(isub).subind - subject ID
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
% Experiments 1 and 2 are combined in the analysis
sublist = [2 4 5 7 11 14 15 16 17 18 19 21 22 23 27 29 31 32 35 36 2:31];
explist = [ones(1,20)*1 ones(1,30)*2];  % which experiment each subject belongs to
nsubs   = length(sublist);

% Subject pairs (subjects who participated in both experiments)
% Row 1: Experiment 1 ID, Row 2: Experiment 2 ID
subpairs = [05 10 11 12 13 14 15; 
            28 23 22 24 34 42 41];

fprintf('\n=== Extracting Behavioral Data for Figure 1 ===\n');
fprintf('Processing %d subjects...\n\n', nsubs);

%% 2. SUBJECT LOOP - EXTRACT BEHAVIORAL DATA

% Initialize output structure
tfdata = struct();

for isub = 1:nsubs
    
    fprintf('Processing Subject %d/%d (ID: S%02d)... ', isub, nsubs, sublist(isub));
    
    % Determine which experiment this subject participated in
    experiment = explist(isub);
    datapath   = fullfile(paths.data, sprintf('Switch_EEG%d', experiment));
    
    % Subject-specific paths
    substrg = sprintf('S%02d', sublist(isub));
    subpath = fullfile(datapath, substrg);
    
    %% 2.1 LOAD DATA
    
    % Load preprocessed EEG and behavioral data
    filename = fullfile(subpath, sprintf('%s_EEG_final.mat', substrg));
    if ~exist(filename, 'file')
        error('File does not exist: %s', filename);
    end
    load(filename);  % loads 'data' and 'behav' structures
    
    ntrials = length(data.trial);
    
    %% 2.2 EXTRACT BEHAVIORAL VARIABLES
    
    % Core behavioral measures
    % Note: These come directly from the behav structure in the EEG file
    
    behaviour = struct();
    
    % Response (1 or 2, corresponding to button press)
    behaviour.resp = behav.resp;
    
    % Reaction time (seconds)
    behaviour.rt = behav.time;
    
    % Accuracy (1 = correct, 0 = incorrect)
    behaviour.acc = behav.corr;
    
    % Decision variable: angular offset of probe from cued and uncued items
    % Column 1: offset from cued item
    % Column 2: offset from uncued item
    behaviour.dv = behav.dv;
    
    % Repetition number (1=switch, 2=first repeat, 3=second repeat, etc.)
    behaviour.repnum = behav.tnum;
    
    % Block number (each block contains 16 trials)
    behaviour.block = behav.block;
    
    % Trial number within block (1-16)
    behaviour.blocktrial = behav.blocktrial;
    
    % Template angles (orientations of the two memory items)
    % Column 1: cued (active) template
    % Column 2: uncued (latent) template
    behaviour.template_angles = behav.relrule;
    
    % Stimulus angles (actual probe orientations)
    if isfield(behav, 'angs')
        behaviour.stimulus_angles = behav.angs;
    end
    
    % Cue type (which item was cued)
    if isfield(behav, 'cue')
        behaviour.cue = behav.cue;
    end
    
    % Presentation side (experiment 2 only - lateral probe presentation)
    if isfield(behav, 'side')
        behaviour.side = behav.side;
    end
    
    %% 2.3 CALCULATE DERIVED VARIABLES
    
    % Template distance: angular distance between two memory items
    % This represents how similar the two items are (0-90 degrees)
    template_dist = circ_dist(behaviour.template_angles(:,1) * 2*pi/180, ...
                              behaviour.template_angles(:,2) * 2*pi/180);
    behaviour.template_distance = abs(template_dist * 180/pi/2);
    
    % Congruency: whether the correct response is the same for both items
    % If probe is between the two items, responses would differ (incongruent)
    % If probe is outside both items, responses are the same (congruent)
    behaviour.congruent = diff(sign(behaviour.dv), [], 2) == 0;
    
    %% 2.4 STORE IN OUTPUT STRUCTURE
    
    tfdata(isub).behaviour  = behaviour;
    tfdata(isub).experiment = experiment;
    tfdata(isub).subind     = sublist(isub);
    
    fprintf('done. (%d trials)\n', ntrials);
    
end

fprintf('\n=== Data extraction complete! ===\n');

%% 3. SAVE RESULTS

fprintf('\nSaving behavioral data...\n');

resultsdir = paths.results;
if ~exist(resultsdir, 'dir')
    mkdir(resultsdir);
end

% Save with standard filename (no timestamp for this file, as it's referenced by other scripts)
outfname = fullfile(resultsdir, 'Study12_TFPow_Beh_Decoding.mat');

save(outfname, 'tfdata', '-v7.3');
fprintf('Results saved to: %s\n', outfname);

%% 4. DISPLAY SUMMARY STATISTICS

fprintf('\n=== Summary Statistics ===\n');

% Calculate overall statistics
all_acc = [];
all_rt = [];
all_repnum = [];

for isub = 1:nsubs
    beh = tfdata(isub).behaviour;
    
    % Exclude missing responses
    valid_trials = ~isnan(beh.rt);
    
    all_acc = [all_acc; beh.acc(valid_trials)];
    all_rt = [all_rt; beh.rt(valid_trials & beh.acc == 1)];  % RT only for correct trials
    all_repnum = [all_repnum; beh.repnum(valid_trials)];
end

fprintf('\nOverall Performance:\n');
fprintf('  Mean accuracy: %.1f%%\n', 100 * mean(all_acc));
fprintf('  Median RT (correct): %.3f s\n', median(all_rt));

fprintf('\nBy Trial Type:\n');
for irep = 1:3
    trials_this_type = all_repnum == irep;
    acc_this_type = all_acc(trials_this_type);
    rt_this_type = all_rt(ismember(find(all_repnum == irep), find(all_acc == 1)));
    
    if irep == 1
        type_name = 'Switch';
    elseif irep == 2
        type_name = 'Repeat 1';
    else
        type_name = 'Repeat 2+';
    end
    
    fprintf('  %s: %.1f%% accuracy, %.3f s RT\n', ...
            type_name, 100 * mean(acc_this_type), median(rt_this_type));
end

% Subject count by experiment
n_exp1 = sum(explist == 1);
n_exp2 = sum(explist == 2);
n_both = size(subpairs, 2);

fprintf('\nSubjects:\n');
fprintf('  Experiment 1 only: %d\n', n_exp1 - n_both);
fprintf('  Experiment 2 only: %d\n', n_exp2 - n_both);
fprintf('  Both experiments: %d\n', n_both);
fprintf('  Total unique subjects: %d\n', n_exp1 + n_exp2 - n_both);

fprintf('\n=== Analysis Complete! ===\n');
fprintf('Next step: Run figure1_behavior_plot.m to generate figures.\n');

%% HELPER FUNCTIONS

function d = circ_dist(alpha, beta)
    % CIRC_DIST - Circular distance between angles
    %
    % Computes the signed angular distance between two angles.
    % Result is in the range [-pi, pi].
    %
    % INPUTS:
    %   alpha - first angle(s) in radians
    %   beta  - second angle(s) in radians
    %
    % OUTPUT:
    %   d - angular distance in radians
    %
    % Note: This is a simplified version. For full functionality,
    % use the CircStat toolbox: https://github.com/circstat/circstat-matlab
    
    d = angle(exp(1i * alpha) .* exp(-1i * beta));
end
