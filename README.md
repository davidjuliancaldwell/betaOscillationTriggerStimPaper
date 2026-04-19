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

The ***R_analysis_scripts*** folder contains the R scripts required to fit linear mixed models and generate statistical plots after the data structure generated from ***master_script_betaStim*** has been run.

- ***betaStim_R_script.R*** — Main mixed effects analysis.
  - **Primary models: 5a and 5a-gf2** (continuous circular phase via `sin(phaseDeg) + cos(phaseDeg)`, ordinal dose × sin/cos, random dose slope, baseline ANCOVA). 5a applies a channel-level phase-quality filter (`phaseVecLength ≥ 0.3`); 5a-gf2 adds a good-fit burst restriction (≥1 conditioning stim with R² > 0.7 and 12-20 Hz) and relaxes the phase-quality threshold to 0.2. Reported together as the primary pair — 5a is the conservative broad test, 5a-gf2 adds the "beta actually present" quality check.
  - **Supporting models: 3a-3e** (binary phaseClass 90/270). Retained as sensitivity checks; no longer primary because binary binning discards circular information and conflates distinct measured phases at multi-phase channels.
  - Effect-size reporting via `emmeans::eff_size()` and `effectsize::eta_squared()`; residual diagnostics and model-comparison tables exported to `output_plots/betaStim_statistical_tables.docx`. See `statistical_audit.md` (Finding 23) and `CLAUDE.md` for the rationale for the 5a/5a-gf2 switch.
- ***R_compare_control_cond.R*** — Within-subject control analyses.
  - **0b5a2e CL vs Playback**: matched-pair sign-flip permutation (exact for n ≤ 20, MC 10k otherwise) with Holm correction across 3 dose bins; paired-probe forest plot by channel × phase condition.
  - **ecb43e targeted vs random**: combined mixed model (`magnitude ~ numStims_ord * (sin_phase + cos_phase) + condType * (sin_phase + cos_phase) + baselineMag_c + (1|channel)`); Wald tests for targeted-specific phase effects beyond random.
  - **Null-burst vs Baseline (0b5a2e)**: validity check for null-burst control; 0/8 channels significant, aggregated perm p=0.226.
  - **Conditioned vs Baseline per-cell permutation** (all 7 subjects): two-sample label-shuffle permutation (10k MC, median diff) with bootstrap 95% CIs, tested on each `(subject × channel × phaseDeg_round × dose)` cell. BH FDR within (subject × dose) respects channel-in-subject nesting. Companion to 5a-gf2: LMM gives the population-level average; per-cell permutation describes the distribution across channels. Outputs include a phase-sorted forest plot with beta-trigger channels highlighted in pink y-axis labels.
  - **Percent modulation from baseline**: per subject × dose descriptive summary.
- ***R_compare_subject_6_random.R*** — Random vs phase-targeted comparison (subject 6).
- ***compare_three_models.R*** — Standalone script comparing Models 3a, 3b, 3c side-by-side with emmip and forest plots (legacy, superseded by the 5a/5a-gf2 framing).

---

### Filter analysis

Code to analyze the performance of the real time filter on the TDT is in the ***test_real_time_filter*** folder

---

David J Caldwell, BSD-3 License
