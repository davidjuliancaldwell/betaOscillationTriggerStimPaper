### Code to analyze evoked potentials from beta oscillation triggered direct electrical stimulation in humans.

This repository contains MATLAB and R code to analyze the data from the beta oscillation triggered stimulation experiments.

The main script is ***master_script_betaStim.m***, which calls other sub scripts and analysis functions.

---

### Finding stimuli timing

The ***find_stims*** folder has code to figure out when the stimuli were delivered and build tables to help further data extraction.

---

### Phase R_analysis_scripts

The function ***B_phaseCalc_allChans_processed.m*** is the key function for the nonlinear sinusoid fits to estimate the phase of delivery across all the channels.

---

### Peak Extraction

The ***peak_extraction*** folder has scripts illustrating sweeping through the subjects and extracting peak to peak voltages.


---

### R analysis

The ***R_analysis_scripts*** folder contains the R scripts required to fit linear mixed models and generate statistical plots after the data structure generated from ***master_script_betaStim*** has been run

- ***betaStim_R_script.R*** — Main mixed effects analysis. Contains four model specifications with progressively improved random effects structures, effect size reporting via `emmeans::eff_size()`, and model comparison. See `statistical_audit.md` for rationale.
- ***R_compare_control_cond.R*** — Closed-loop vs playback control comparison (subject 7). Includes Cohen's d effect sizes and permutation tests (10,000 permutations, two-sided) for the dose-response interaction.
- ***R_compare_subject_6_random.R*** — Random vs phase-targeted comparison (subject 6).

---

### Filter analysis

Code to analyze the performance of the real time filter on the TDT is in the ***test_real_time_filter*** folder

---

David J Caldwell, BSD-3 License
