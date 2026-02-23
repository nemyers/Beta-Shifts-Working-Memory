# Quick Start Guide

New to this codebase? Start here!

## First-Time Setup (5 minutes)

### 1. Install Prerequisites

```matlab
% Check MATLAB version
ver

% Required: MATLAB R2019b or later
% Required: FieldTrip toolbox
```

### 2. Configure Paths

```bash
# Copy the template
cp switchbeta_config_template.m switchbeta_config.m

# Edit config.m with your paths
nano switchbeta_config.m  # or use your editor
```

In `switchbeta_config.m`, update:
```matlab
paths.study = '/path/to/your/project';      % Where you cloned this repo
paths.fieldtrip = '/path/to/fieldtrip';     % Where FieldTrip is installed
paths.toolbox = '/path/to/custom/toolbox';  % Your analysis toolbox
```

### 3. Verify Setup

```matlab
% Test configuration
paths = config();

% Should print your paths without errors
```

## Understanding the Code Structure

```
switchbeta-analysis/
│
├── README.md                    ← Start here for overview
├── switchbeta_config_template.m            ← Copy this to switchbeta_config.m
├── CLEANING_GUIDE.md           ← How to clean remaining scripts
├── DEPENDENCIES.md             ← What external functions are needed
│
├── figure1_behavior_calculate.m     ← Run behavioral analysis
├── figure1_behavior_plot.m          ← Plot Figure 1
│
├── figure2_tf_switch_calculate.m    ← Calculate TF power differences
├── figure2_tf_switch_plot.m         ← Plot Figure 2
│
├── figure3_tf_itemdistance_calculate.m  ← Calculate power by item distance
├── figure3_tf_itemdistance_plot.m       ← Plot Figure 3
│
├── figure4_decoding_correlation_calculate.m  ← Calculate decoding correlations
└── figure4_decoding_correlation_plot.m       ← Plot Figure 4
```

## Running Your First Analysis

### Test on a Single Subject

```matlab
% In any calculate script, change:
% for isub = 1:nsubs
% to:
for isub = 1:1  % Just first subject
    % ...
end
```

This lets you test the code quickly!

### Full Pipeline

```matlab
% 1. Calculate results (slow - may take hours)
figure1_behavior_calculate;
figure2_tf_switch_calculate;
figure3_tf_itemdistance_calculate;
figure4_decoding_correlation_calculate;

% 2. Generate plots (fast)
figure1_behavior_plot;
figure2_tf_switch_plot;
figure3_tf_itemdistance_plot;
figure4_decoding_correlation_plot;
```