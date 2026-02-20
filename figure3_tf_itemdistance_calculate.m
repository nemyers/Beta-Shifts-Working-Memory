% figure3_tf_itemdistance_calculate.m - Calculate power modulation by item distance
%
% This script analyzes how oscillatory power (particularly beta-band) scales
% with the angular distance between the two memory items ("item distance").
% Larger item distances require larger-magnitude priority shifts.
%
% KEY ANALYSES:
%   1. GLM regression: template distance predicts TF power
%   2. Binning: average power for each item distance level
%
% KEY FINDING:
%   Beta power (15-25Hz) scales with item distance on Switch trials:
%   larger switches show greater beta desynchronization
%
% INPUTS:
%   - config.m (path configuration)
%   - Subject EEG data (*_EEG_final.mat)
%   - Time-frequency data (*_HanningTF_Laplace.mat)
%   - Artifact rejection files (*_semiautomaticAR.mat)
%
% OUTPUTS:
%   - results/Study12_Figure3_TF_ItemDistance_Results_[timestamp].mat
%     Contains:
%       output.betamat_templatedist - GLM betas for template distance regressor
%       output.tfmutempdist         - TF power binned by template distance
%       output.timeTF, output.freq  - Time and frequency axes
%       output.regressors           - Names of GLM regressors used
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

%% 2. SUBJECT LOOP - PROCESS EACH SUBJECT
for isub = 1:nsubs
    
    fprintf('\n=== Processing Subject %d/%d (ID: S%02d) ===\n', isub, nsubs, sublist(isub));
    
    % Determine which experiment this subject participated in
    experiment = explist(isub);
    datapath   = fullfile(paths.data, sprintf('Switch_EEG%d', experiment));
    tfpath     = fullfile(paths.derivatives, 'TF', sprintf('Switch_EEG%d', experiment));
    
    % Subject-specific paths
    substrg = sprintf('S%02d', sublist(isub));
    subpath = fullfile(datapath, substrg);
    
    %% 2.1 LOAD DATA
    
    % Load preprocessed EEG and behavioral data
    filename = fullfile(subpath, sprintf('%s_EEG_final.mat', substrg));
    if ~exist(filename, 'file')
        error('File does not exist: %s', filename);
    end
    fprintf('Loading EEG data...\n');
    load(filename);  % loads 'data' and 'behav' structures
    
    % Clean up data structure
    if isfield(data, 'nTrials')
        data = rmfield(data, 'nTrials');
    end
    
    ntrials = length(data.trial);
    time    = data.time{1};
    
    % Load time-frequency data
    fname = fullfile(tfpath, substrg, sprintf('%s_HanningTF_Laplace.mat', substrg));
    fprintf('Loading time-frequency data...\n');
    load(fname);  % loads 'datatfr' structure
    
    % Load artifact rejection info
    rejfile = fullfile(subpath, sprintf('%s_semiautomaticAR.mat', substrg));
    if exist(rejfile, 'file')
        load(rejfile, 'rejsemiautomatic');
    else
        warning('Artifact rejection file not found: %s', rejfile);
        rejsemiautomatic = [];
    end
    
    %% 2.2 IDENTIFY USABLE TRIALS
    
    % Trials to exclude:
    % - Marked as artifacts
    % - No response (RT is NaN)
    % - Previous trial had no response (to avoid contamination)
    % - Infinite RT values
    
    rt   = behav.time;
    miss = isnan(rt);
    missprev = [false; miss(1:end-1)];
    
    artefact = zeros(ntrials, 1);
    artefact(rejsemiautomatic) = 1;
    artefact(miss) = 1;
    artefact(missprev) = 1;
    
    usable = ~artefact;
    usable = usable & ~isinf(behav.time) & ~isinf([0; behav.time(1:end-1)]);
    
    fprintf('S%02d: %d usable trials (%.1f%%)\n', ...
        sublist(isub), nnz(usable), 100*nnz(usable)/length(usable));
    
    ntrials_usable = nnz(usable);
    
    %% 2.3 EXTRACT BEHAVIORAL VARIABLES
    
    % Trial type: 1=switch, 2=first repeat, 3=second+ repeat
    tnum = behav.tnum(usable);
    tnum(tnum > 3) = 3;  % collapse repeats beyond 2nd
    
    % Previous trial's repetition number (for GLM regressor)
    tnum_prev  = [0; behav.tnum(1:end-1)];
    tnum_prev2 = [0; 0; behav.tnum(1:end-2)];
    tnum_prev3 = [0; 0; 0; behav.tnum(1:end-3)];
    tnum_prev  = tnum_prev(usable);
    tnum_prev2 = tnum_prev2(usable);
    tnum_prev3 = tnum_prev3(usable);
    
    % For repeat trials, use the tnum from the corresponding switch trial
    tnum_prev(tnum == 2) = tnum_prev2(tnum == 2);
    tnum_prev(tnum == 3) = tnum_prev3(tnum == 3);
    
    % Other behavioral variables
    block     = behav.block(usable);
    acc       = behav.corr(usable);
    rt        = behav.time(usable);
    blktrl    = behav.blocktrial(usable);
    
    % Previous trial variables
    rtprev = [0; behav.time(1:end-1)];
    rtprev = rtprev(usable);
    acprev = [0; behav.corr(1:end-1)];
    acprev = acprev(usable);
    
    % Spatial side (experiment 2 only; lateral probe presentation)
    if experiment == 2
        side = behav.side(usable);
    else
        side = ones(nnz(usable), 1);
    end
    
    % Rule orientations (active and latent templates)
    rls = behav.relrule(usable, :);
    
    % TEMPLATE DISTANCE: Angular distance between the two memory items
    % This is the key variable for this analysis
    % It represents how different the two items are (range: 0-90 degrees)
    cong = circ_dist(rls(:,1)*2*pi/180, rls(:,2)*2*pi/180) * 180/pi/2;
    templatedist = abs(round(cong*1000)/1000);  % absolute angular distance
    
    utemplatedist = unique(templatedist);
    ntemplatedist = length(utemplatedist);
    
    fprintf('Template distance range: %.2f to %.2f degrees\n', ...
        min(utemplatedist), max(utemplatedist));
    
    %% 2.4 PREPROCESS TIME-FREQUENCY DATA
    
    % Extract TF power for usable trials
    dat = datatfr.powspctrm(usable, :, :, :);
    
    timeTF = datatfr.time;
    ntime  = length(timeTF);
    freq   = datatfr.freq;
    nfreq  = length(freq);
    nchans = size(dat, 2);
    
    % Log-transform power (convert to dB)
    dat = 10 * log10(dat);
    
    % ORTHOGONALIZE WITH RESPECT TO BASELINE
    % This removes variance explained by baseline power fluctuations
    % Important for regression analyses
    fprintf('Orthogonalizing TF power with respect to baseline...\n');
    
    baseline_win = [-0.75 -0.25];  % baseline window in seconds
    ib = timeTF >= baseline_win(1) & timeTF <= baseline_win(2);
    baseline_power = mean(dat(:, :, :, ib), 4);
    
    % Demean baseline power across trials
    baseline_power = bsxfun(@minus, baseline_power, mean(baseline_power, 1));
    
    % For each frequency and channel, regress out baseline power
    dat_ortho = nan(size(dat));
    for ifreq = 2:nfreq  % skip DC component
        for ichan = 1:nchans
            % Design matrix: intercept + baseline power
            X = [ones(ntrials_usable, 1), baseline_power(:, ichan, ifreq)];
            pX = pinv(X);
            
            % Data: power at all timepoints for this channel/frequency
            Y = squeeze(dat(:, ichan, ifreq, :));
            
            % Regress and get residuals
            beta = pX * Y;
            Yhat = X * beta;
            dat_ortho(:, ichan, ifreq, :) = Y - Yhat;
        end
    end
    dat = dat_ortho;
    clear dat_ortho;
    
    %% 2.5 DEFINE CHANNEL GROUPS
    
    chanlabels = datatfr.label;
    
    % Store channel info on first subject
    if isub == 1
        output = struct();
        output.timeTF = timeTF;
        output.freq   = freq;
        output.subpairs = subpairs;
        
        % Define anatomical channel groups
        output.frontalchans = ismember(chanlabels, ...
            {'AF3', 'AFz', 'AF4', 'F1', 'Fz', 'F2'});
        output.centralchans = ismember(chanlabels, ...
            {'C1', 'Cz', 'C2', 'CP1', 'CPz', 'CP2', 'P1', 'Pz', 'P2'});
        output.postchans = ismember(chanlabels, ...
            {'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'O1', 'Oz', 'O2'});
    end
    
    %% 2.6 GLM ANALYSIS: TEMPLATE DISTANCE EFFECT
    
    % Run separate GLM for each trial type (switch, repeat1, repeat2+)
    fprintf('Running GLM analyses...\n');
    
    for icond = 1:3
        % Select trials for this condition
        % Exclude first trial in block (no valid previous trial)
        itrl = (blktrl > 1) & (acc >= 0) & (tnum == icond);
        
        % Design matrix (all regressors z-scored)
        X = [templatedist, tnum_prev, rt, acc, rtprev, acprev, side, blktrl, block];
        X = X(itrl, :);
        y = dat(itrl, :, :, :);
        
        % Z-score for comparable regression coefficients
        y = zscore(y, [], 1);
        X = zscore(X, [], 1);
        
        % Mass univariate GLM
        % Output: betas has shape [nregressors x nchans x nfreqs x ntime]
        betas = massGLM(X, y, false);
        
        % Store betas: [nsubs x nchans x nfreqs x ntime x nregressors x nconds]
        betamat_templatedist(isub, :, :, :, :, icond) = permute(betas, [2 3 4 1]);
    end
    
    % Store regressor names on final subject
    if isub == nsubs
        output.betamat_templatedist = betamat_templatedist;
        output.regressors = {'templatedist', 'tnum_prev', 'rt', 'acc', ...
                             'rtprev', 'acprev', 'side', 'blktrl', 'block'};
    end
    
    %% 2.7 BIN POWER BY TEMPLATE DISTANCE
    
    % For visualization: bin TF power by template distance and trial type
    fprintf('Binning power by template distance...\n');
    
    for icond = 1:3  % switch, repeat1, repeat2+
        for itemp_dist = 1:ntemplatedist
            % Select trials
            itrl = (tnum == icond) & ...
                   (templatedist == utemplatedist(itemp_dist)) & ...
                   (blktrl > 1);
            
            % Average power across trials
            tfmutempdist(isub, :, :, :, itemp_dist, icond) = ...
                mean(dat(itrl, :, :, :), 1);
        end
    end
    
    % Store on final subject
    if isub == nsubs
        output.tfmutempdist = tfmutempdist;
        output.utemplatedist = utemplatedist;
    end
    
end  % end subject loop

%% 3. SAVE RESULTS

fprintf('\n=== Saving Results ===\n');

resultsdir = paths.results;
if ~exist(resultsdir, 'dir')
    mkdir(resultsdir);
end

timestamp = datestr(now, 'yyyymmdd-HHMM');
outfname = fullfile(resultsdir, ...
    sprintf('Study12_Figure3_TF_ItemDistance_Results_%s.mat', timestamp));

save(outfname, 'output', '-v7.3');
fprintf('Results saved to: %s\n', outfname);

fprintf('\n=== Analysis Complete! ===\n');
