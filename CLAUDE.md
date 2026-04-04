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

`R_analysis_scripts/betaStim_R_script.R` contains four models with progressively improved random effects. **Model 3** (summary-level) is the primary reported model:

```r
# Model 1: Original — random intercepts only. numStims DF inflated (~37K).
fit.intercepts.only = lmer(absDiff ~ numStims * phaseClass + betaLabels +
  numStims:betaLabels + (1|sid/channel), data=dataNoBaseline)

# Model 2: Adds random dose slopes per subject. Fixes numStims DF (~5),
# but phaseClass DF still inflated (~4K).
fit.trial.level = lmer(absDiff ~ numStims * phaseClass +
  (1|sid) + (0+numStims|sid) + (1|channel), data=dataNoBaseline)

# Model 3 (primary): Summary-level (one median per cell). No singularity.
# Random intercepts only — dose slopes removed because 6 subjects cannot
# support a 3x3 covariance matrix (correlations hit 1.0). setToDeliverPhase
# is not used as a random effect (it is a fixed experimental condition).
fit.summary.level = lmer(magnitude ~ numStims * phaseClass +
  (1|sid) + (1|channel), data=summaryNoBaseline)

# Model 4: Trial-level with condition nesting. Kept for reference.
# Singular fit due to (1|channel:setToDeliverPhase) redundancy and
# near-saturated dose slope covariance.
fit.nested.condition = lmer(absDiff ~ numStims * phaseClass +
  (0+numStims|sid) + (1|channel) + (1|channel:setToDeliverPhase), data=dataNoBaseline)
```

Key data structure notes:
- `numStims` (dose) varies trial-to-trial within a channel — real trial-level predictor
- `phaseClass` is a channel-level constant — the circular mean of phase-at-delivery, binned to 90/270, replicated across all trials (`multipleSubj_GLMM_script_PP.m:190-191`)
- `setToDeliverPhase` is a fixed experimental condition (not a random grouping variable)
- Channel IDs are unique per subject (subjectNum*100 + raw channel), so `(1|channel)` implicitly nests within subject
- 6 subjects, 31 channels, 120 summary observations (median per cell), ~37K trials before aggregation
- phaseClass DF limitation: Satterthwaite assigns ~84 DF for phaseClass instead of the ideal ~30 (between-channel). Adding `(1|channel:setToDeliverPhase)` would correct this but causes singularity — redundant with `(1|channel)` for 22/31 single-phase channels. Does not affect conclusions (phaseClass p=0.42 at DF=84)

Model 3 key results (easystats reporting added to R script):
- numStims: F(2,84) = 9.60, **p = 0.0002**, partial eta² = 0.19 (large)
- phaseClass: F(1,84) = 0.64, p = 0.42, partial eta² = 0.008
- Interaction: F(2,84) = 1.85, p = 0.16, partial eta² = 0.04
- Performance: conditional R² = 0.998, marginal R² = 0.0005, ICC = 0.998
- Emmeans: dose effect at phase 270 ([5,inf) vs [1,2] = +11.3 uV, p=0.0001); no dose effect at phase 90; phase contrast at [5,inf) = 5.8 uV, p=0.063

`R_analysis_scripts/R_compare_control_cond.R` compares closed-loop (0b5a2e) vs playback control (0b5a2ePlayBack) with Cohen's d effect sizes and permutation tests. See `statistical_audit.md` for full findings.
