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
cp config_template.m config.m

# Edit config.m with your paths
nano config.m  # or use your editor
```

In `config.m`, update:
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
├── config_template.m            ← Copy this to config.m
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

## What Each Figure Shows

**Figure 1** - Behavioral switch costs
- Do switches hurt performance?
- Does item similarity matter?

**Figure 2** - Oscillatory power changes during switches
- Beta decreases, theta increases on switches
- Spatial topography of effects

**Figure 3** - Beta power scales with update magnitude
- Bigger switches → more beta desynchronization
- Theta doesn't scale the same way

**Figure 4** - Beta power predicts decoding strength
- Lower beta → better decoding of new item
- Specific to newly prioritized item

## Common Issues

### "Function not found: massGLM"
→ See DEPENDENCIES.md for custom functions needed

### "File not found: S02_EEG_final.mat"  
→ Check your paths.data directory structure

### "Out of memory"
→ Process fewer subjects at once, or use higher-memory machine

### Results look different from paper
→ Did you merge subjects who did both experiments?
→ Are you using the same trial exclusion criteria?

## Tips for Success

1. **Start small** - Test on 1 subject before running all
2. **Check intermediate outputs** - Don't wait until the end
3. **Read the paper** - Understanding helps debugging
4. **Comment as you go** - Future you will thank you
5. **Version control** - Commit working code before experimenting

## Getting Help

1. Check the [README](README.md) for detailed documentation
2. Review [CLEANING_GUIDE.md](CLEANING_GUIDE.md) for code patterns
3. Look at [DEPENDENCIES.md](DEPENDENCIES.md) for missing functions
4. Email: nicholas.myers@nottingham.ac.uk

## Next Steps

Once basic setup works:
1. Run test analysis on subset of data
2. Compare with published results
3. Clean remaining scripts (see CLEANING_GUIDE.md)
4. Add your own analyses!

## Citation

Remember to cite the paper if you use this code:

> Myers, N.E., Stokes, M.G., & Muhle-Karbe, P.S. (2026). Human Beta 
> Oscillations Reflect Magnitude and Fidelity of Priority Shifts in 
> Working Memory. *Journal of Neuroscience*, DOI: 10.1523/JNEUROSCI.1548-25.2026
