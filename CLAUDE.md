# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Analysis pipeline for the paper "Dose Dependent Enhancement of Cortically Evoked Potentials During Beta-Oscillation Phase Triggered Direct Cortical Stimulation of Human Cortex" (David J. Caldwell). Examines how beta oscillation phase at stimulation time modulates evoked potential magnitude across 8 human subjects with intracranial ECoG recordings.

## Running the Pipeline

**MATLAB** (primary): Open MATLAB in the repo root, then run `master_script_betaStim.m`. Set `generateIntermediateData = 1` for the full pipeline (stim table building, peak extraction, phase calculation), or `0` to skip data generation and only produce plots/tables.

**R** (statistical analysis): Run `R_analysis_scripts/betaStim_R_script.R` after MATLAB has generated the output CSV in `data/output_table/`. Requires packages: lme4, lmerTest, afex, emmeans, sjPlot, multcomp, ggplot2, plyr, dplyr, Hmisc, here.

## Pipeline Architecture

The pipeline runs in lettered stages (A, B, C) that must execute in order:

1. **A: Stim table building** (`find_stims/A_BuildStimTablesFirst6.m`, `A_BuildStimTablesSubj7andPlayback.m`) -- Identifies stimulation times from raw TDT trigger data, outputs to `data/stim_timing_data/`
2. **B: Neural data extraction** -- Two parallel tracks:
   - `peak_extraction/B_ExtractNeuralData_PP_reref.m` -- Extracts peak-to-peak CEP magnitudes with median CAR re-referencing, outputs to `data/EP_data/`
   - `B_phaseCalc_allChans_processed.m` -- Nonlinear sinusoid fitting (`analysis_functions/sinfit.m`) to estimate beta phase at stimulation across all channels, outputs to `data/phase_data/`. This step is computationally expensive.
3. **C: Aggregation & visualization** -- `multipleSubj_GLMM_script_PP.m` combines EP and phase data into CSV tables in `data/output_table/`, then plotting scripts generate figures
4. **R: Statistical modeling** -- Linear mixed effects models (LME4) with subject/channel random effects on the output CSV

## Key Configuration

- `setup_environment.m` / `Z_Constants.m`: Define subject IDs (SIDS), folder paths for all data directories. Both files set the same variables; `setup_environment.m` is called by `master_script_betaStim.m`.
- Subject IDs: `{'d5cd55', 'c91479', '7dbdec', '9ab7ab', '702d24', 'ecb43e', '0b5a2e', '0b5a2ePlayback'}`
- Subjects `702d24` and `0b5a2ePlayback` are excluded in the R analysis.

## Data Directories (under `data/`)

- `ECoG_data/` -- Raw ECoG recordings (not in repo, loaded at runtime)
- `stim_timing_data/` -- Stimulation timing tables (stage A output)
- `phase_data/` -- Phase calculation results per subject/channel (stage B output)
- `EP_data/` -- Peak-to-peak evoked potentials (stage B output)
- `output_table/` -- CSV tables for R analysis (stage C output). Multiple threshold variants (30, 50, 100 uV minimum)
- `coordinates/` -- Subject-specific MRI electrode coordinates (at repo root level)

## Key Analysis Parameters

- Beta band: ~10-30 Hz, extracted via nonlinear sinusoid fitting (not bandpass)
- CEP magnitude thresholds: 25-1500 uV range; current preferred minimum is 100 uV
- Phase binning: 45-degree bins (8 bins per cycle)
- Dose levels: Base, [1,2], [3,4], [5,inf) conditioning stimuli
- Phase classes: depolarizing (90) vs hyperpolarizing (270)

## Code Organization

| Directory | Purpose |
|---|---|
| `analysis_functions/` | Signal processing: sinfit, CEP amplitude extraction, filtering, re-referencing, permutation tests |
| `brain_plotting_functions/` | Cortical surface visualization with electrode overlays |
| `find_stims/` | Stage A: stimulation timing detection from TDT data |
| `peak_extraction/` | Stage B: peak-to-peak voltage extraction |
| `phase_visualizations/` | Phase distribution plots (circular histograms, cortex overlays) |
| `plotting_functions/` | General visualization utilities |
| `manuscript_generate_scripts/` | Publication figure generation |
| `test_real_time_filter/` | Validation of the real-time beta filter on TDT hardware |
| `R_analysis_scripts/` | Linear mixed effects models and statistical plots |

## Core Functions

- `analysis_functions/sinfit.m` -- Nonlinear least-squares sinusoidal fitting for phase/amplitude/period estimation. Central to the phase analysis.
- `analysis_functions/CEPamp.m` -- Extracts peak amplitudes and latencies of cortically evoked potentials
- `analysis_functions/rereference_CAR_median.m` -- Common average re-referencing with median subtraction
- `analysis_functions/phase_circstats_calc.m` -- Circular statistics for phase distributions

## Statistical Models

`R_analysis_scripts/betaStim_R_script.R` contains four models with progressively improved random effects. The recommended models are **Model 3** (summary-level, no singularity) and **Model 4** (trial-level, correct DF for all predictors):

```r
# Model 1: Original — random intercepts only. numStims DF inflated (~37K).
fit.intercepts.only = lmer(absDiff ~ numStims * phaseClass + betaLabels +
  numStims:betaLabels + (1|sid/channel), data=dataNoBaseline)

# Model 2: Adds random dose slopes per subject. Fixes numStims DF (~5),
# but phaseClass DF still inflated (~4K).
fit.trial.level = lmer(absDiff ~ numStims * phaseClass +
  (1|sid) + (0+numStims|sid) + (1|channel), data=dataNoBaseline)

# Model 3: Summary-level (one median per cell). All DF correct. No singularity.
fit.summary.level = lmer(magnitude ~ numStims * phaseClass +
  (1|sid) + (0+numStims|sid) + (1|channel), data=summaryNoBaseline)

# Model 4: Trial-level with condition nesting. Fixes phaseClass DF (~31)
# by adding the level at which phaseClass operates. Singular fit expected
# (single-phase subjects have one condition per channel).
fit.nested.condition = lmer(absDiff ~ numStims * phaseClass +
  (0+numStims|sid) + (1|channel) + (1|channel:setToDeliverPhase), data=dataNoBaseline)
```

Key data structure notes:
- `numStims` (dose) varies trial-to-trial within a channel — real trial-level predictor
- `phaseClass` is a channel-level constant — the circular mean of phase-at-delivery, binned to 90/270, replicated across all trials (`multipleSubj_GLMM_script_PP.m:190-191`)
- Channel IDs are unique per subject (subjectNum*100 + raw channel), so `(1|channel)` implicitly nests within subject
- 6 subjects, 31 channels, 49 channel x condition cells, ~37K trials after exclusions

`R_analysis_scripts/R_compare_control_cond.R` compares closed-loop (0b5a2e) vs playback control (0b5a2ePlayBack) with Cohen's d effect sizes and permutation tests. See `statistical_audit.md` for full findings.
