# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Analysis pipeline for the paper "Dose Dependent Enhancement of Cortically Evoked Potentials During Beta-Oscillation Phase Triggered Direct Cortical Stimulation of Human Cortex" (David J. Caldwell). Examines how beta oscillation phase at stimulation time modulates evoked potential magnitude across 8 human subjects with intracranial ECoG recordings.

## Running the Pipeline

**MATLAB** (primary): Open MATLAB in the repo root, then run `master_script_betaStim.m`. Set `generateIntermediateData = 1` for the full pipeline (stim table building, peak extraction, phase calculation), or `0` to skip data generation and only produce plots/tables.

**R** (statistical analysis): Run `R_analysis_scripts/betaStim_R_script.R` after MATLAB has generated the output CSV in `data/output_table/`. Requires packages: lme4, lmerTest, afex, emmeans, sjPlot, multcomp, ggplot2, plyr, dplyr, Hmisc, here, effectsize, performance, report, officer, flextable.

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
- Subject `0b5a2ePlayback` is excluded from the main R analysis (used only in the CL vs playback comparison). `702d24` is now included (previously excluded, 1 channel only — no convergence issues).

## Data Directories (under `data/`)

- `ECoG_data/` -- Raw ECoG recordings (not in repo, loaded at runtime)
- `stim_timing_data/` -- Stimulation timing tables (stage A output)
- `phase_data/` -- Phase calculation results per subject/channel (stage B output)
- `EP_data/` -- Peak-to-peak evoked potentials (stage B output)
- `output_table/` -- CSV tables for R analysis (stage C output). Multiple threshold variants (30, 50, 100 uV minimum)
- `coordinates/` -- Subject-specific MRI electrode coordinates (at repo root level)

## Key Analysis Parameters

- Beta band: ~10-30 Hz, extracted via nonlinear sinusoid fitting (not bandpass)
- CEP magnitude filtering (two stages):
  - **Channel-level** (MATLAB, `multipleSubj_GLMM_script_PP.m:185`): exclude channels where mean baseline EP < 100 uV (`epThresholdMag`). Screens out channels with weak/absent evoked potentials.
  - **Trial-level** (R, `betaStim_R_script.R:36-37`): exclude individual trials with magnitude < 25 uV or > 1500 uV. Removes artifacts and non-responses.
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

`R_analysis_scripts/betaStim_R_script.R` contains five summary-level models (3a-3e), plus trial-level models (1, 2, 4) kept for reference. All summary models collapse to one **median** per (subject x channel x phaseClass x numStims) cell to eliminate pseudoreplication. Median is used consistently throughout: cell summaries, baseline computation, and permutation test statistics.

```r
# Model 3a: absDiff (baseline median pre-subtracted per channel).
# Random intercepts only. Most powerful for dose (p=0.0002) but no random slopes.
fit.absDiff = lmer(absDiff ~ numStims * phaseClass +
  (1|sid) + (1|channel), data=summaryNB)  # 120 obs

# Model 3b: Raw magnitude with baseline as a fourth dose level.
# Random linear dose slope per subject (doseNum=0,1,2,3).
# Baseline-by-phase confound weakens interaction (p=0.34).
fit.modelD = lmer(magnitude ~ numStims * phaseClass +
  (1+doseNum|sid) + (1|channel), data=summaryAll)  # 151 obs

# Model 3c (ANCOVA): Raw magnitude, baseline as fixed covariate (centered,
# beta~1.03). Absorbs between-channel variance (channel SD 153→10 uV).
fit.ancova = lmer(magnitude ~ numStims * phaseClass + baselineMag_c +
  (1+doseNum|sid) + (1|channel), data=summaryNB_ancova)  # 120 obs

# Model 3d: Ordinal dose (polynomial .L/.Q contrasts). Reparameterization
# of 3c — identical fit. afex::mixed with expand_re tests .L and .Q
# separately; .Q random variance ~0, confirming linear slope sufficient.
fit.ordinal_afex = afex::mixed(magnitude ~ numStims_ord * phaseClass +
  baselineMag_c + (numStims_ord || sid) + (1|channel),
  data=summaryNB_ancova, expand_re=TRUE, per_parameter="numStims_ord")

# Model 3e: Fully numeric dose (fixed + random). Most parsimonious.
# Nested within 3c/3d — LRT confirms quadratic unnecessary (p=0.76).
fit.numeric = lmer(magnitude ~ doseNum * phaseClass + baselineMag_c +
  (1+doseNum|sid) + (1|channel), data=summaryNB_ancova)  # 120 obs
```

Key data structure notes:
- All cell summaries use **median** (magnitude, absDiff, percentDiff, baseline)
- Permutation tests also use **median** as the test statistic
- `numStims` (dose) varies trial-to-trial within a channel — real trial-level predictor
- `phaseClass` is a channel-level constant — the circular mean of phase-at-delivery, binned to 90/270
- Channel IDs are unique per subject (subjectNum*100 + raw channel), so `(1|channel)` implicitly nests within subject
- 6 subjects, 31 channels, 120 summary observations (no baseline) or 151 (with baseline)
- `doseNum` is a numeric encoding of the dose factor; `dose_linpoly` = `contr.poly(3)[,1]` is an equivalent linear rescaling used with ordinal models
- None of the five primary models are singular

Model comparison (all on median, after phaseClass fix):

| Model | AIC | Dose p | Interaction p | Phase contrast |
|-------|-----|--------|---------------|----------------|
| 3a (absDiff, no slopes) | 910 | **~0.0002** | ~0.86 | ns |
| 3c (ANCOVA, categorical) | 912 | 0.188 | **0.859** | ns |
| 3d (ordinal, .L only) | 914 | 0.188 | **0.616** (.L coeff) | ns |
| 3e (numeric, linear) | 917 | 0.106 | **0.615** | ns |

**Critical bug fix (2026-04-06)**: `multipleSubj_GLMM_script_PP.m` had a phase label swap for multi-phase subjects — both `phaseClass` and `setToDeliverPhase` were inverted for c91479, 0b5a2e, 0b5a2ePlayBack. The previously reported phase × dose interaction (p = 0.045-0.110) was an artifact. After fix, the interaction is non-significant (p = 0.62-0.86). The dose main effect remains significant in the no-random-slopes sensitivity analysis (Model 3a).

## Phase Label Verification

`verify_phase_consistency.m` validates that measured phase-at-delivery (circular mean from sinfit) is consistent with intended target phase for every subject's beta reference channel. Run from repo root in MATLAB.

### Beta reference channels (from `valueSet` in `multipleSubj_GLMM_script_PP.m`)

| Subject | betaChan | Type | desiredF | Hardware convention |
|---------|----------|------|----------|-------------------|
| d5cd55 | 53 | s | 180 | Single condition |
| c91479 | 64 | m | [0, 180] | stims(8)==1→pos(0°), stims(8)==0→neg(180°) |
| 7dbdec | 4 | s | 180 | Single condition |
| 9ab7ab | 51 | s | 270 | Single condition |
| 702d24 | 5 | m | [90, 270] | stims(8)==1→pos(90°), stims(8)==0→neg(270°) |
| ecb43e | 55 | t | [270, 90, rand, rand] | **Inverted**: stims(8)==0→pos(270°), stims(8)==1→neg(90°) |
| 0b5a2e | 31 | m | [90, 270] | stims(8)==1→pos(90°), stims(8)==0→neg(270°) |
| 0b5a2ePlayBack | 31 | m | [90, 270] | stims(8)==1→pos(90°), stims(8)==0→neg(270°) |

**Note**: 0b5a2e and 0b5a2ePlayBack share the same patient and betaChan=31. `B_ExtractNeuralData_PP_reref.m` previously had a stale value of betaChan=23 for 0b5a2e (fixed 2026-04-08); this variable was unused in the extraction pipeline so it had no effect on output data.

### Verification results (2026-04-08)

All beta reference channels show measured circular mean within 40° of target, **except** 0b5a2ePlayBack neg condition (ch31 circMean=0.1°, target=270°, 90° off). This is the playback condition where stimulation timing was replayed asynchronously from the live beta oscillation, so degraded phase targeting is expected.

Non-reference EP channels show spatial phase offsets of up to ~180° from the beta reference channel due to cortical beta phase propagation gradients. This is expected physiology, not a labeling error.

c91479 0/180 targets verified: all EP channels show pos condition closer to 0° and neg condition closer to 180° — no swap.

## Residual Diagnostics

`betaStim_R_script.R` includes residual diagnostics for all five summary-level models (3a-3e):
- **Shapiro-Wilk** normality test
- **Skewness** and **excess kurtosis** (moment-based, computed from raw residuals)
- **ggplot2 QQ plots** with skewness/kurtosis annotated in subtitles, saved to `output_plots/betaStim_qq_3{a-e}_*.png`
- **Residuals-vs-fitted plots** with loess smoother, saved to `output_plots/betaStim_resid_vs_fitted_3{a-e}_*.png`

Current findings: skewness is acceptable (|skew| < 1 across all models), but excess kurtosis is elevated (5-7), indicating heavy tails from a few channels with extreme CEP values. LME is relatively robust to leptokurtic residuals at these sample sizes.

## Manuscript .docx Export

The R script generates `output_plots/betaStim_statistical_tables.docx` using `officer` + `flextable`, containing:
- Residual diagnostics table (Shapiro-Wilk, skewness, kurtosis for all 5 models)
- Model comparison table (AIC, BIC, singularity status)
- Random effects, Type III ANOVA (Satterthwaite df), and fixed effects tables for models 3a, 3c, 3e
- EMM pairwise dose and phase contrasts with 95% CIs (Tukey-adjusted)
- Cohen's d effect size tables for dose contrasts

`R_analysis_scripts/R_compare_control_cond.R` compares closed-loop (0b5a2e) vs playback control (0b5a2ePlayBack) with Cohen's d effect sizes and permutation tests (using median as test statistic). See `statistical_audit.md` for full findings.
