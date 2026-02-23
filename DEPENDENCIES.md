# Toolbox Dependencies

This document lists the custom functions used in the analysis scripts that come from external toolboxes.

## Required Custom Functions

Some functions are called in the scripts but not part of standard MATLAB or FieldTrip. 
They are included in the utils folder.

### Core Analysis Functions

**`massGLM.m`**
- Used in: All calculate scripts
- Purpose: Mass univariate General Linear Model analysis
- Runs GLM across multiple channels/frequencies/timepoints

### Experiment-Specific Functions

**`Switch_MergeExperiments.m`**
- Used in: Plotting scripts
- Purpose: Average data from subjects who participated in both experiments