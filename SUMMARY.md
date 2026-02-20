# Code Cleaning Summary

## What I've Created For You

I've prepared a clean, well-documented code repository structure for your GitHub release:

### Core Documentation
1. **README.md** - Comprehensive project overview
   - Installation instructions
   - Data structure requirements  
   - Analysis pipeline description
   - Citation information

2. **QUICKSTART.md** - Get-started-in-5-minutes guide
   - Setup instructions
   - Test workflows
   - Common troubleshooting

3. **CLEANING_GUIDE.md** - Step-by-step instructions for cleaning remaining scripts
   - Before/after code examples
   - Specific patterns to change
   - Quality checklist

4. **DEPENDENCIES.md** - Documentation of external function requirements
   - What custom functions are needed
   - How to handle them (include vs. reference)
   - CircStat and other toolboxes

### Configuration
5. **config_template.m** - Path configuration template
   - Users copy this to config.m
   - No hardcoded paths in any analysis scripts
   - Automatic path verification

6. **.gitignore** - Prevents committing data/results
   - config.m is git-ignored (contains local paths)
   - Data directories excluded
   - Standard MATLAB ignores

### Example Cleaned Scripts

7. **figure3_tf_itemdistance_calculate.m** - FULLY CLEANED calculation script
   - Complete header documentation
   - Section headers throughout
   - Inline comments explaining logic
   - Uses config.m for paths
   - Clear variable names
   - Progress indicators
   - Proper error handling
   - Timestamped output saving

8. **figure1_behavior_plot_TEMPLATE.m** - Template for plotting scripts
   - Shows structure for plot scripts
   - Loads from timestamped results
   - Helper functions included
   - Comments on what each panel shows

## What You Need To Do

### Immediate Next Steps

1. **Download the cleaned code**
   - All files are in the folder I've prepared
   - Review the structure and documentation

2. **Read the CLEANING_GUIDE.md**
   - Shows before/after examples
   - Explains each type of change needed

3. **Clean remaining scripts** using the guide:
   - Figure 1 calculate & plot
   - Figure 2 calculate & plot
   - Figure 4 calculate & plot
   - Any supplementary figure scripts

4. **Handle custom functions** (see DEPENDENCIES.md):
   - Decide: include, reference, or rewrite
   - massGLM.m is critical - include or document
   - CircStat - consider git submodule
   - Switch_MergeExperiments - probably include

5. **Test everything**:
   ```matlab
   % Create config.m from template
   paths = config();
   
   % Run cleaned scripts on test data
   figure3_tf_itemdistance_calculate;
   % etc.
   ```

6. **Create GitHub repository**:
   ```bash
   git init
   git add .
   git commit -m "Initial commit: cleaned analysis code"
   git remote add origin [your-repo-url]
   git push -u origin main
   ```

### Key Improvements Made

**Before (original code):**
- Hardcoded paths for specific computers
- Minimal documentation
- Unclear variable names
- No section organization
- Inconsistent file handling
- Mixed analysis and setup code

**After (cleaned code):**
- ✅ Portable configuration system
- ✅ Comprehensive documentation
- ✅ Inline comments explaining logic  
- ✅ Clear section headers
- ✅ Consistent use of fullfile()
- ✅ Separated concerns
- ✅ Progress indicators
- ✅ Error handling
- ✅ Timestamped outputs

## Example of Changes

### Path Configuration
**Old:**
```matlab
workstation = 'nick_pc';
switch workstation
    case 'nick_pc'
        onedrivepath = 'E:\Myers\OneDrive - The University of Nottingham';
        studypath = [onedrivepath '\Projects\El-SwitchBeta'];
end
```

**New:**
```matlab
%% 1. SETUP PATHS AND PARAMETERS
paths = config();  % Loads from config.m
```

### Documentation
**Old:**
```matlab
clearvars; close('all'); clc
for isub = 1:nsubs
    load(filename);
```

**New:**
```matlab
% figure3_tf_itemdistance_calculate.m - Calculate power modulation by item distance
%
% This script analyzes how oscillatory power scales with angular distance...
%
% KEY FINDING:
%   Beta power scales with item distance on Switch trials
%
% INPUTS/OUTPUTS: [detailed list]

%% 1. SETUP PATHS AND PARAMETERS
paths = config();

%% 2. SUBJECT LOOP
for isub = 1:nsubs
    fprintf('\n=== Processing Subject %d/%d ===\n', isub, nsubs);
    
    %% 2.1 LOAD DATA
    fprintf('Loading EEG data...\n');
```

### File Paths
**Old:**
```matlab
filename = sprintf('%s%s%s_EEG_final.mat',subpath,fs,substrg);
```

**New:**
```matlab
filename = fullfile(subpath, sprintf('%s_EEG_final.mat', substrg));
```

## Checklist Before Publishing

- [ ] All scripts cleaned following CLEANING_GUIDE.md
- [ ] config_template.m documented and tested
- [ ] README.md reviewed and accurate
- [ ] DEPENDENCIES.md lists all custom functions
- [ ] Custom functions either included or documented
- [ ] .gitignore prevents committing data
- [ ] Test on fresh checkout (different computer if possible)
- [ ] Add LICENSE file (I recommend MIT)
- [ ] Scripts run without errors (at least on test subject)
- [ ] Output files created in correct locations
- [ ] Consider adding example output structure

## Optional Enhancements

Once basics are working, consider:

1. **Add example data subset** (if ethically possible)
   - One subject's data (anonymized)
   - Lets users test code immediately

2. **Automated tests**
   - Simple scripts that verify outputs match expected structure

3. **Docker container** or **MATLAB project**
   - Makes environment reproducible

4. **Jupyter notebooks** or **live scripts**
   - Interactive demonstrations of key analyses

5. **Zenodo archive**
   - Get DOI for code
   - Cite alongside paper

## Timeline Suggestion

- **Day 1**: Review documentation, set up config.m
- **Day 2-3**: Clean remaining calculate scripts
- **Day 4**: Clean plotting scripts
- **Day 5**: Handle custom functions/dependencies
- **Day 6**: Test full pipeline
- **Day 7**: Create GitHub repo and publish!

## Questions?

If anything is unclear:
1. Check the example cleaned script (figure3_...)
2. Review CLEANING_GUIDE.md for specific patterns
3. Look at QUICKSTART.md for workflow
4. Feel free to ask me for clarification!

## What Makes This Code GitHub-Ready?

Your code will be:
- **Portable**: Runs on any computer with proper configuration
- **Documented**: Others can understand what it does
- **Organized**: Clear structure and file naming
- **Testable**: Can verify it works correctly
- **Citable**: Properly attributed to authors
- **Findable**: Good README helps discoverability
- **Professional**: Reflects well on your research group

Good luck with your GitHub release! 🚀
