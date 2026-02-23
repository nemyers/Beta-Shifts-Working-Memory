# Beta Oscillations and Working Memory Priority Shifts

Code accompanying: Myers, N.E., Stokes, M.G., & Muhle-Karbe, P.S. (2026). Human Beta Oscillations Reflect Magnitude and Fidelity of Priority Shifts in Working Memory. *Journal of Neuroscience*.

## Overview

This repository contains MATLAB code to reproduce the analyses and figures from our study.

## Requirements

### Software
- MATLAB (tested on R2019b or later)
- [FieldTrip toolbox](http://www.fieldtriptoolbox.org/) (tested with version 20200607)
- Circular statistics functions (CircStat toolbox), included as submodule

### Custom Toolboxes
The code depends on several custom analysis functions, which are included in utils. If a function is not found in your MATLAB path, it's likely from the custom toolbox referenced in the scripts. The most commonly used functions include:
- `massGLM.m` - Mass univariate GLM analysis
- `SetupToolboxPath.m` - Path setup utility


## Setup

1. **Clone this repository:**
   ```bash
   git clone https://github.com/yourusername/switchbeta-analysis
   cd switchbeta-analysis
   ```

2. **Configure paths:**
   - Copy `switchbeta_config_template.m` to `switchbeta_config.m`
   - Edit `switchbeta_config.m` to point to your local paths:
     - `paths.data` - Location of raw EEG data
     - `paths.derivatives` - Location of preprocessed time-frequency data
     - `paths.results` - Where to save analysis outputs
     - `paths.fieldtrip` - FieldTrip toolbox location
     - `paths.toolbox` - Custom analysis toolbox location

3. **Add FieldTrip to your path (done in switchbeta_config):**
   ```matlab
   addpath('/path/to/fieldtrip');
   ft_defaults;
   ```

## Data Structure

The analysis expects data organized as follows:
```
data/
├── Switch_EEG1/        # Experiment 1 data
│   ├── S02/            # Subject 02
│   │   ├── S02_EEG_final.mat # broadband EEG data and behavioural data
│   │   └── S02_semiautomaticAR.mat # vector of artifactual trials, to be rejected before final analyses
│   └── ...
├── Switch_EEG2/        # Experiment 2 data
│   └── ...
└── derivatives/
    └── TF/
        ├── Switch_EEG1/
        │   ├── S02/
        │   │   └── S02_HanningTF_Laplace.mat
        │   └── ...
        └── Switch_EEG2/
            └── ...
```

### Data Files
- `*_EEG_final.mat` - Preprocessed EEG data with behavioral information
- `*_semiautomaticAR.mat` - Artifact rejection indices
- `*_HanningTF_Laplace.mat` - Time-frequency decomposition (Hanning tapers, Laplacian filtered)

## Analysis Pipeline

The analysis is organized into separate calculation and plotting scripts for each figure:

### Figure 1: Behavioral Effects
- **Calculate:** `figure1_behavior_calculate.m`
- **Plot:** `figure1_behavior_plot.m`
- **Outputs:** Behavioral switch costs, psychometric fits, item distance effects

### Figure 2: Switch-Related Oscillatory Power Changes  
- **Calculate:** `figure2_tf_switch_calculate.m`
- **Plot:** `figure2_tf_switch_plot.m`
- **Outputs:** Time-frequency power differences (Switch vs. Repeat trials)
- **Key finding:** Beta power (15-25Hz) decreases on switch trials

### Figure 3: Power Scaling with Update Magnitude
- **Calculate:** `figure3_tf_itemdistance_calculate.m`  
- **Plot:** `figure3_tf_itemdistance_plot.m`
- **Outputs:** Correlation between beta power and item distance (update magnitude)
- **Key finding:** Larger priority switches → greater beta desynchronization

### Figure 4: Power Predicting Decoding Fidelity
- **Calculate:** `figure4_decoding_correlation_calculate.m`
- **Plot:** `figure4_decoding_correlation_plot.m`
- **Outputs:** Relationship between beta power and WM decoding strength
- **Key finding:** Lower beta power → stronger decoding of newly prioritized item

## Running the Analyses

### Quick Start
To reproduce all figures (assuming data are available):

```matlab
% 1. Configure paths
switchbeta_config;

% 2. Run calculations (may take several hours)
figure1_behavior_calculate;
figure2_tf_switch_calculate;
figure3_tf_itemdistance_calculate;
figure4_decoding_correlation_calculate;

% 3. Generate plots
figure1_behavior_plot;
figure2_tf_switch_plot;
figure3_tf_itemdistance_plot;
figure4_decoding_correlation_plot;
```

### Individual Analyses
Each calculation script can be run independently and will save results to the configured results directory with a timestamp.

## Key Analysis Parameters

### Subject Information
- **Experiment 1:** 20 participants (IDs: 2, 4, 5, 7, 11, 14, 15, 16, 17, 18, 19, 21, 22, 23, 27, 29, 31, 32, 35, 36)
- **Experiment 2:** 30 participants (IDs: 2-31)
- **Combined:** 43 unique participants (7 participated in both experiments)

### Channel Groups
- **Central channels:** C1, Cz, C2, CP1, CPz, CP2, P1, Pz, P2 (beta effects)
- **Frontal channels:** AF3, AFz, AF4, F1, Fz, F2 (theta effects)
- **Posterior channels:** PO7, PO3, POz, PO4, PO8, O1, Oz, O2 (alpha/decoding)

## Citation

If you use this code, please cite:

```bibtex
@article{myers2026beta,
  title={Human Beta Oscillations Reflect Magnitude and Fidelity of Priority Shifts in Working Memory},
  author={Myers, Nicholas E. and Stokes, Mark G. and Muhle-Karbe, Paul S.},
  journal={Journal of Neuroscience},
  year={2026},
  doi={10.1523/JNEUROSCI.1548-25.2026}
}
```

## License

This code is released under the MIT License. See LICENSE file for details.

## Contact

For questions about the code:
- Nicholas Myers: nicholas.myers@nottingham.ac.uk