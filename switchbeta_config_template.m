function paths = config()
% CONFIG - Configuration file for SwitchBeta analysis
%
% Returns a structure containing all path information needed for the analysis.
% 
% INSTRUCTIONS FOR SETUP:
%   1. Copy this file to config.m (which is git-ignored)
%   2. Edit the paths below to match your local system
%   3. Run config() to verify paths are correct
%
% OUTPUTS:
%   paths - Structure with fields:
%     .study       - Main project directory
%     .data        - Raw EEG data location
%     .derivatives - Preprocessed data (time-frequency, etc.)
%     .results     - Where to save analysis outputs
%     .scripts     - Location of analysis scripts
%     .fieldtrip   - FieldTrip toolbox path
%     .toolbox     - Custom analysis toolbox path
%
% Author: Nicholas Myers
% Created: 2026
% For: Myers, Stokes, & Muhle-Karbe (2026) JNeurosci

%% EDIT THESE PATHS FOR YOUR SYSTEM
% Example paths shown - replace with your actual paths

% Main project directory
paths.study = 'E:\Myers\OneDrive - The University of Nottingham\Projects\El-SwitchBeta';

% Data directories
paths.data        = fullfile(paths.study, 'data');
paths.derivatives = fullfile(paths.study, 'data', 'derivatives');
paths.results     = fullfile(paths.study, 'results');
paths.scripts     = fullfile(paths.study, 'scripts');

% Toolbox paths
paths.fieldtrip = '/Users/yourname/matlab/fieldtrip';
paths.toolbox   = '/Users/yourname/matlab/toolbox';

%% AUTOMATIC SETUP
% Add paths to MATLAB search path
addpath(paths.fieldtrip);
addpath(paths.toolbox);
addpath(paths.scripts);

% Initialize FieldTrip
ft_defaults;

%% VERIFY PATHS EXIST
% Check that all critical directories exist
critical_paths = {paths.study, paths.fieldtrip, paths.toolbox};
for i = 1:length(critical_paths)
    if ~exist(critical_paths{i}, 'dir')
        warning('Path does not exist: %s', critical_paths{i});
    end
end

% Create results directory if it doesn't exist
if ~exist(paths.results, 'dir')
    fprintf('Creating results directory: %s\n', paths.results);
    mkdir(paths.results);
end

fprintf('\nConfiguration loaded successfully!\n');
fprintf('Study path: %s\n', paths.study);
fprintf('FieldTrip:  %s\n', paths.fieldtrip);
fprintf('Results:    %s\n\n', paths.results);

end
