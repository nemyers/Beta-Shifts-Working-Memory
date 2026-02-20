# Code Cleaning Guide

This document explains how to clean up the remaining MATLAB scripts for GitHub release.

## What I've Done So Far

I've created cleaned versions of:
1. **README.md** - Complete project documentation
2. **config_template.m** - Path configuration template  
3. **figure3_tf_itemdistance_calculate.m** - Example cleaned calculation script
4. **.gitignore** - Standard gitignore for MATLAB projects

## Cleaning Checklist for Remaining Scripts

For each of the remaining scripts, follow this process:

### 1. Add Header Documentation

Replace:
```matlab
clearvars; close('all'); clc
```

With:
```matlab
% SCRIPT_NAME - Brief one-line description
%
% Detailed description of what this script does and why.
%
% KEY ANALYSES:
%   1. First main analysis
%   2. Second main analysis  
%
% KEY FINDING:
%   Brief summary of main result from the paper
%
% INPUTS:
%   - config.m (path configuration)
%   - List specific data files needed
%
% OUTPUTS:
%   - results/Output_File_Name_[timestamp].mat
%     Description of what's in the output
%
% Author: Nicholas E. Myers
% Cleaned for GitHub: 2026
% For: Myers, Stokes, & Muhle-Karbe (2026) J. Neuroscience

clearvars;
close all;
clc;
```

### 2. Replace Hardcoded Paths

**Old code:**
```matlab
workstation = 'nick_pc';
switch workstation
    case 'nick_pc'
        onedrivepath = 'E:\Myers\OneDrive - The University of Nottingham';
        ...
end
```

**New code:**
```matlab
%% 1. SETUP PATHS AND PARAMETERS
% Load configuration
paths = config();

% Subject information
sublist = [2 4 5 7 11 14 15 16 17 18 19 21 22 23 27 29 31 32 35 36 2:31];
explist = [ones(1,20)*1 ones(1,30)*2];
nsubs   = length(sublist);
```

### 3. Add Section Headers

Organize code into clear sections with headers:

```matlab
%% 1. SETUP PATHS AND PARAMETERS

%% 2. SUBJECT LOOP - PROCESS EACH SUBJECT
for isub = 1:nsubs
    
    %% 2.1 LOAD DATA
    
    %% 2.2 IDENTIFY USABLE TRIALS
    
    %% 2.3 EXTRACT BEHAVIORAL VARIABLES
    
    %% 2.4 RUN ANALYSIS
    
end

%% 3. GROUP-LEVEL STATISTICS

%% 4. SAVE RESULTS
```

### 4. Add Inline Comments

Add comments explaining:
- What each variable represents
- Why certain choices were made
- Links to methods in the paper

**Old code:**
```matlab
tnum = behav.tnum(usable);
tnum(tnum>3) = 3;
```

**New code:**
```matlab
% Trial type: 1=switch, 2=first repeat, 3=second+ repeat
tnum = behav.tnum(usable);
tnum(tnum > 3) = 3;  % collapse repeats beyond 2nd
```

### 5. Improve Variable Names (Where Possible)

Some improvements:
- `itrl` → `trial_idx` or keep if very common in your codebase
- `isub` → keep (standard in loops)
- `dv` → add comment explaining it's "decision variable" or similar
- `td` → add comment that it's "template distance"

### 6. Fix File Paths

**Old:**
```matlab
filename = sprintf('%s%s%s_EEG_final.mat',subpath,fs,substrg);
```

**New:**  
```matlab
filename = fullfile(subpath, sprintf('%s_EEG_final.mat', substrg));
```

### 7. Add Progress Indicators

**Old:**
```matlab
for isub = 1:nsubs
    load(filename);
```

**New:**
```matlab
for isub = 1:nsubs
    fprintf('\n=== Processing Subject %d/%d (ID: S%02d) ===\n', ...
        isub, nsubs, sublist(isub));
    
    fprintf('Loading EEG data...\n');
    load(filename);
```

### 8. Consistent Results Saving

**Old:**
```matlab
save_output = true;
if save_output
    resultsdir = sprintf('%s/results',studypath);
    outfname = sprintf('%s/Study12_...',resultsdir);
    outfname = sprintf('%s_%s.mat',outfname,datestr(now,'yyyymmdd-HHMM'));
    save(outfname,'output','-v7.3');
end
```

**New:**
```matlab
%% 3. SAVE RESULTS
fprintf('\n=== Saving Results ===\n');

resultsdir = paths.results;
if ~exist(resultsdir, 'dir')
    mkdir(resultsdir);
end

timestamp = datestr(now, 'yyyymmdd-HHMM');
outfname = fullfile(resultsdir, ...
    sprintf('Study12_Figure2_Results_%s.mat', timestamp));

save(outfname, 'output', '-v7.3');
fprintf('Results saved to: %s\n', outfname);

fprintf('\n=== Analysis Complete! ===\n');
```

## Script-Specific Notes

### Figure 1 (Behavior)
- Complex psychometric modeling with swap errors
- Consider breaking into subfunctions for model fitting
- Document the mixture model approach

### Figure 2 (TF Switch/Repeat)  
- Main oscillatory power analysis
- Document baseline correction approach
- Explain cluster-based permutation testing

### Figure 4 (Decoding Correlation)
- Decoding analysis setup
- Explain the Mahalanobis distance approach
- Document cross-validation scheme

## Testing Your Cleaned Scripts

Before committing:
1. **Test each script** runs without errors (with dummy/subset data if needed)
2. **Check all paths** resolve correctly with config.m
3. **Verify output files** are created in the right location
4. **Compare results** with original scripts on a test subject

## Naming Conventions

Stick to these file naming patterns:
- Calculation scripts: `figureN_description_calculate.m`
- Plotting scripts: `figureN_description_plot.m`
- Utility functions: `util_description.m`
- All lowercase with underscores

## Functions to Consider Creating

If you find repeated code blocks, consider extracting to functions:
- `load_subject_data.m` - Load EEG, behavior, TF data for one subject
- `identify_usable_trials.m` - Apply artifact rejection and trial exclusion
- `merge_experiments.m` - Merge data from subjects in both experiments

## Final Checks

Before publishing:
- [ ] All scripts have header documentation
- [ ] No hardcoded paths remain
- [ ] config_template.m is well documented  
- [ ] README.md is complete and accurate
- [ ] .gitignore prevents committing data files
- [ ] LICENSE file is included (MIT or CC-BY recommended)
- [ ] Test on a fresh clone of the repository

## Questions?

If unclear about any of these guidelines:
1. Look at the cleaned figure3_tf_itemdistance_calculate.m as a reference
2. Check MATLAB style guides (MathWorks has good ones)
3. Prioritize readability - future you will thank you!
