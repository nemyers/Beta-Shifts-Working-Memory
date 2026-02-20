# Toolbox Dependencies

This document lists the custom functions used in the analysis scripts that come from external toolboxes.

## Required Custom Functions

The following functions are called in the scripts but not part of standard MATLAB or FieldTrip. You'll need to either:
1. Include them from your custom toolbox
2. Replace them with standard MATLAB equivalents
3. Provide them in a `utils/` folder in the repository

### Core Analysis Functions

**`massGLM.m`**
- Used in: All calculate scripts
- Purpose: Mass univariate General Linear Model analysis
- Runs GLM across multiple channels/frequencies/timepoints
- Alternative: Could use MATLAB's `fitlm` in a loop, but would be slower

**`SetupToolboxPath.m`**  
- Used in: All scripts (via config.m)
- Purpose: Add subdirectories of toolbox to path
- Alternative: Replace with standard `addpath(genpath(...))` but be careful about conflicts

### Statistics Functions

**`circ_dist.m`** (from CircStat toolbox)
- Used in: Behavioral analysis scripts
- Purpose: Compute circular distance between angles
- Note: CircStat is open source: https://github.com/circstat/circstat-matlab
- You should include it as a git submodule or dependency

### Experiment-Specific Functions

**`Switch_MergeExperiments.m`**
- Used in: Plotting scripts
- Purpose: Average data from subjects who participated in both experiments
- You'll need to include this function or document how to handle duplicate subjects

### Optional/Plotting Functions

Functions that may be used in plotting scripts:
- Custom plotting utilities (colors, formatting, etc.)
- May want to replace with standard MATLAB plotting code

## Recommended Approach

### Option 1: Include a utils/ folder
```
switchbeta-analysis/
├── utils/
│   ├── massGLM.m
│   ├── Switch_MergeExperiments.m
│   └── ...
```

Then in config.m:
```matlab
addpath(fullfile(paths.study, 'utils'));
```

### Option 2: Document as prerequisites

In README.md, add:
```markdown
## Prerequisites

This code requires the following custom functions:
- `massGLM.m` - Available at: [link]
- CircStat toolbox - Available at: https://github.com/circstat/circstat-matlab

Install these before running the analysis.
```

### Option 3: Provide simplified versions

For functions that do simple operations, you could provide simplified alternatives:

**Example: Simplified circ_dist**
```matlab
function d = circ_dist(alpha, beta)
    % CIRC_DIST - Circular distance between angles
    % Simplified version for this analysis
    % For full functionality, use CircStat toolbox
    d = angle(exp(1i*alpha) .* exp(-1i*beta));
end
```

## CircStat Toolbox Functions Used

If using CircStat, these specific functions are called:
- `circ_dist.m` - Circular distance

You can either:
1. Add CircStat as a git submodule
2. Copy just the needed functions (check license allows this)
3. Direct users to install CircStat separately

## Checking for Missing Functions

To find all undefined functions in your scripts, run:

```matlab
% Check for undefined functions
all_scripts = dir('*.m');
for i = 1:length(all_scripts)
    fprintf('Checking %s...\n', all_scripts(i).name);
    [fList, pList] = matlab.codetools.requiredFilesAndProducts(all_scripts(i).name);
    
    % Find functions not in MATLAB path
    for j = 1:length(fList)
        if ~exist(fList{j}, 'file')
            warning('Missing function: %s', fList{j});
        end
    end
end
```

## Decision Matrix

| Function | Include? | Why |
|----------|----------|-----|
| massGLM | Yes | Core to analysis, likely custom |
| SetupToolboxPath | Replace | Simple path management |
| circ_dist | Submodule | Part of CircStat (open source) |
| Switch_MergeExperiments | Yes | Experiment-specific |
| Plotting helpers | Maybe | Depends on complexity |

## Next Steps

1. **Identify all custom functions** in your scripts
2. **Check licenses** - Can you redistribute them?
3. **Document requirements** in README.md
4. **Test on clean install** - Does everything run?

## Contact About Toolbox

If users have trouble with the custom toolbox:
- Point them to this documentation
- Consider releasing minimal required functions
- Provide contact for questions about specific functions
