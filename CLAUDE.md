# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Analysis pipeline for the paper "Dose Dependent Enhancement of Cortically Evoked Potentials During Beta-Oscillation Phase Triggered Direct Cortical Stimulation of Human Cortex" (David J. Caldwell). Examines how beta oscillation phase at stimulation time modulates evoked potential magnitude across 8 human subjects with intracranial ECoG recordings.

## Local Machine Setup (one-time, macOS)

The `coordinates/` folder contains only `surf/obj/` subdirectories — binary `.mat` files (electrode locations, cortex meshes) are gitignored and must be copied from OneDrive.

**OneDrive source**: `OneDrive-UCSF/UWGoogleDriveBackup/djcald_uw_backup_v2/Data-pistachio/Subjects/`

For each subject in `{d5cd55, c91479, 7dbdec, 9ab7ab, 702d24, ecb43e, 0b5a2e}`, copy:
- `trodes.mat` → `coordinates/<sid>/trodes.mat`
- `tail_trodes.mat` (7dbdec only) → `coordinates/<sid>/tail_trodes.mat`
- `other/` → `coordinates/<sid>/other/`
- `surf/*.mat` → `coordinates/<sid>/surf/`

**Required MATLAB Add-On**: Install **"Cyclic color map"** by Chad Greene from MATLAB Add-Ons Manager. Provides `phasemap`, `phasebar`, `phasewrap` used by `plotPhase_distributions_cortex.m` and `plot_phase_cortex.m`.

## Running the Pipeline

**MATLAB** (primary): Open MATLAB in the repo root, then run `master_script_betaStim.m`. Set `generateIntermediateData = 1` for the full pipeline (stim table building, peak extraction, phase calculation), or `0` to skip data generation and only produce plots/tables.

**R** (statistical analysis): Run `R_analysis_scripts/betaStim_R_script.R` after MATLAB has generated the output CSV in `data/output_table/`. Requires packages: lme4, lmerTest, afex, emmeans, sjPlot, multcomp, ggplot2, plyr, dplyr, Hmisc, here, effectsize, performance, report, car, splines, officer, flextable.

## Pipeline Architecture

The pipeline runs in lettered stages (A, B, C) that must execute in order:

1. **A: Stim table building** (`find_stims/A_BuildStimTablesFirst6.m`, `A_BuildStimTablesSubj7andPlayback.m`) -- Identifies stimulation times from raw TDT trigger data, outputs to `data/stim_timing_data/`
2. **B: Neural data extraction** -- Two parallel tracks:
   - `peak_extraction/B_ExtractNeuralData_PP_reref.m` -- Extracts peak-to-peak CEP magnitudes with median CAR re-referencing, outputs to `data/EP_data/`
   - `B_phaseCalc_allChans_processed.m` -- Nonlinear sinusoid fitting (`analysis_functions/sinfit.m`) to estimate beta phase at stimulation across all channels, outputs to `data/phase_data/`. Computationally expensive.
3. **C: Aggregation & visualization** -- `multipleSubj_GLMM_script_PP.m` combines EP and phase data into CSV tables in `data/output_table/`, then plotting scripts generate figures
4. **R: Statistical modeling** -- Linear mixed effects models (LME4) with subject/channel random effects on the output CSV

## Key Configuration

- `setup_environment.m` / `Z_Constants.m`: Define subject IDs (SIDS), folder paths for all data directories. Both files set the same variables; `setup_environment.m` is called by `master_script_betaStim.m`.
- Subject IDs: `{'d5cd55', 'c91479', '7dbdec', '9ab7ab', '702d24', 'ecb43e', '0b5a2e', '0b5a2ePlayback'}`
- `0b5a2ePlayback` excluded from main R analysis (CL vs playback comparison only). `702d24` included (1 channel, no convergence issues).
- `external_deps/` contains bundled `CircStat2012a/` and `savitzkyGolay.m`/`sgolayfilt_complete.m`. `setup_environment.m` uses `addpath(genpath(locationsDir))` only — no external paths required.

## Data Directories (under `data/`)

- `ECoG_data/` -- Raw ECoG recordings (gitignored, symlinked from OneDrive: `OneDrive-UCSF/Research/UW_research/betastim/Data-from-xps/ECoG_data`). Contains `{sid}_ECoG.mat` and `{sid}_forBetaPhase.mat` per subject.
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
- Phase classes: depolarizing (90°) vs hyperpolarizing (270°)
- All cell summaries use **median** (magnitude, absDiff, percentDiff, baseline). Median is robust to right-skewed EP distributions (overall skew=1.41). Permutation tests also use median as test statistic.

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

`R_analysis_scripts/betaStim_R_script.R` contains five summary-level models (3a-3e), plus trial-level models (1, 2, 4) kept for reference. All summary models collapse to one **median** per (subject × channel × phaseClass × numStims) cell to eliminate pseudoreplication.

```r
# Model 3a: absDiff (baseline median pre-subtracted per channel).
# Random intercepts only. Most powerful for dose (p=0.0002) but no random slopes.
fit.absDiff = lmer(absDiff ~ numStims * phaseClass +
  (1|sid) + (1|channel), data=summaryNB)  # 120 obs

# Model 3b: Raw magnitude with baseline as a fourth dose level.
# Random linear dose slope per subject.
fit.modelD = lmer(magnitude ~ numStims * phaseClass +
  (1+doseNum|sid) + (1|channel), data=summaryAll)  # 151 obs

# Model 3c (ANCOVA): Raw magnitude, baseline as fixed covariate (centered).
# Absorbs between-channel variance (channel SD 153→10 uV).
fit.ancova = lmer(magnitude ~ numStims * phaseClass + baselineMag_c +
  (1+doseNum|sid) + (1|channel), data=summaryNB_ancova)  # 120 obs

# Model 3d: Ordinal dose (polynomial .L/.Q contrasts). Reparameterization of 3c.
fit.ordinal_afex = afex::mixed(magnitude ~ numStims_ord * phaseClass +
  baselineMag_c + (numStims_ord || sid) + (1|channel),
  data=summaryNB_ancova, expand_re=TRUE, per_parameter="numStims_ord")

# Model 3e: Fully numeric dose (fixed + random). Most parsimonious.
fit.numeric = lmer(magnitude ~ doseNum * phaseClass + baselineMag_c +
  (1+doseNum|sid) + (1|channel), data=summaryNB_ancova)  # 120 obs
```

Key data structure notes:
- `numStims` (dose) varies trial-to-trial within a channel — real trial-level predictor
- `phaseClass` is a channel-level constant — the circular mean of phase-at-delivery, binned to 90/270
- Channel IDs are unique per subject (subjectNum*100 + raw channel), so `(1|channel)` implicitly nests within subject
- 6 subjects, 31 channels, 120 summary observations (no baseline) or 151 (with baseline)
- `doseNum` is a numeric encoding of the dose factor; `dose_linpoly` = `contr.poly(3)[,1]` is an equivalent linear rescaling used with ordinal models
- None of the five primary models are singular

Model comparison (all on median, after phaseClass fix):

| Model | AIC | Dose p | Interaction p |
|-------|-----|--------|---------------|
| 3a (absDiff, no slopes) | 910 | **~0.0002** | ~0.86 |
| 3c (ANCOVA, categorical) | 912 | 0.188 | 0.859 |
| 3d (ordinal, .L only) | 914 | 0.188 | 0.616 |
| 3e (numeric, linear) | 917 | 0.106 | 0.615 |

### Model 5: Continuous circular phase (sin/cos decomposition)

Replaces binary `phaseClass` (90/270) with `sin(phase)` and `cos(phase)`. The output CSV includes `phaseDeg` (channel-level circular mean in degrees, from `circ_mean` of stims with R² > 0.7), plus quality metrics `phaseVecLength`, `phaseCircStd`, `phaseOmnibusP`. Groups by `(sid, phaseDeg_round, numStims, channel)` — keeps conditions with different measured phases separate. Before phase-quality filters: 141 obs. With default filters (`minPhaseVecLength_5a=0.3`): 78 obs.

```r
# Model 5a (primary): Ordinal dose × sin/cos + betaLabels, ANCOVA.
# Channel-level phase. phaseVecLength >= 0.3. Non-singular.
fit.sincos.ordinal = lmer(magnitude ~ numStims_ord * (sin_phase + cos_phase) +
  betaLabels + baselineMag_c +
  (1+dose_linpoly|sid) + (1|channel), data=summaryNB_m5)  # 78 obs (r>=0.3)

# Model 5a-gf (secondary): Good beta fit restriction, per-burst phase.
# Groups by (sid, numStims, channel) — collapses across delivered phases.
# Intercepts only — random dose slope singular at 96 obs.
# Dose p is anti-conservative without the random slope.
fit.sincos.ordinal.gf = lmer(magnitude ~ numStims_ord * (sin_phase + cos_phase) +
  betaLabels + baselineMag_c +
  (1|sid) + (1|channel), data=summaryNB_gf)  # 96 obs

# Model 5a-gf2 (preferred secondary): Good beta fit, channel-level phase.
# phaseVecLength >= 0.2. Random dose slope non-singular.
fit.sincos.ordinal.gf2 = lmer(magnitude ~ numStims_ord * (sin_phase + cos_phase) +
  betaLabels + baselineMag_c +
  (1+dose_linpoly|sid) + (1|channel), data=summaryNB_gf2)  # 102 obs (r>=0.2)
```

Default filter thresholds: `minPhaseVecLength_5a = 0.3`, `minPhaseVecLength_gf2 = 0.2`, `minGoodBetaPerBurst = 1`, `minBurstVecLength_gf = 0`.

Model 5 results:
- **5a** (78 obs, 19 ch): Dose.L p=**0.068** (trend), cos_phase p=**0.090** (trend). Effect size 6.4 µV.
- **5a-gf** (96 obs, intercepts only): Dose.L p=0.047 — **anti-conservative**. Sensitivity check only.
- **5a-gf2** (102 obs, 25 ch): Dose.L p=**0.049**, sin_phase p=0.19, cos_phase p=0.13. Effect size 8.6 µV.

Sensitivity analysis (`output_plots/betaStim_phase_quality_sensitivity.csv`): Models 5a and 5a-gf2 fit at r ∈ {0, 0.1, 0.2, 0.3, 0.4}. Dose effect strengthens monotonically with r.

ANCOVA rationale: `baselineMag_c` (grand-mean-centered baseline median per channel) absorbs ~10-fold between-channel variance. Baseline trials excluded to avoid baseline-by-phase confound. `betaLabels` (1=beta reference channel) included as additive fixed effect — consistently non-significant (p=0.16-0.55).

Effect sizes: `d_total = estimate / sqrt(Var_int + E[x²]*Var_slope + Var_channel + Var_resid)`, where E[x²]=1/3 for `contr.poly(3)`.

### Model 6: Continuous dose spline (EXPLORATORY)

Uses `nCondStims` with natural splines `ns(nCondStims, df=3)`. Without random slopes (singular at 7 subjects), p-values are anti-conservative. Retained for descriptive visualization of the continuous dose-response shape.

### Within-subject analyses (`R_compare_control_cond.R`)

**0b5a2e CL vs Playback**: Probes matched by sequential position (CL probe N ↔ PB probe N). Combined model: `magnitude ~ numStims_ord * (sin_phase + cos_phase) + condition * (sin_phase + cos_phase) + baselineMag_c + (1|channel_raw)`. `condition` ref = PB.

Three phase quantities in play (NOT interchangeable):
1. **`setToDeliverPhase`** — hardware intended target (90° or 270°), channel-independent
2. **`phaseDeg`/`phaseDeg_round`** — channel-level circular mean of sinfit phases. EP channels show spatial offsets up to ~180° from beta reference (expected physiology). Some channels show nearly identical phases for both target conditions — they cannot be treated as providing two distinct phase conditions.
3. **`burstCircMean`/`cl_burst_phase`** — per-burst measured phase from `compute_burst_phase_precision.m` → `{sid}_burst_phase_precision.csv`

For 0b5a2e CL vs PB: matched-pair permutation (sign-flip) is the primary inferential tool; mixed model is used for effect sizes. After `minPhaseVecLength_clpb = 0.2` filter (13 of 16 channel × condition cells kept), none of the aggregated sign-flip permutations clear Holm correction across 3 dose bins.

**ecb43e targeted vs random**: 36 obs, severely underpowered. No phase modulation detected.

**Percent modulation from baseline**: Grand mean [1,2]=1.3%, [3,4]=3.1%, [5,inf)=4.9%. Strongest responder: c91479 (14-19%).

**Null-burst validity test (0b5a2e)**: 0/8 channels p<0.05 uncorrected; aggregated perm p=0.226. Null-burst EPs are statistically indistinguishable from baseline, validating the null-burst control.

## Phase Label Verification

`verify_phase_consistency.m` validates measured phase-at-delivery against intended target for each subject's beta reference channel.

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

All beta reference channels show measured circular mean within 40° of target, except 0b5a2ePlayBack neg condition (expected — asynchronous playback). Non-reference EP channels show spatial offsets up to ~180° due to cortical beta phase propagation gradients (expected physiology, not a labeling error).

**Phase label audit (2026-04-12)**: All pipeline stages verified consistent — `stims(8)` → burst type → phase variable → CSV column → R analysis. No phase-label swap bugs remain. The critical fix was `correctIdx = 2 - bt` in `multipleSubj_GLMM_script_PP.m` (2026-04-06), backported to `phase_vs_peak.m` (2026-04-12). Type 's' subjects (single condition) were trivially unaffected; type 't' (ecb43e, inverted hardware convention) was coincidentally unaffected.

## Residual Diagnostics

Models 3a-3e, 5a-5c: Shapiro-Wilk normality test, skewness/kurtosis, QQ plots (`output_plots/betaStim_qq_*.png`), residuals-vs-fitted plots (`output_plots/betaStim_resid_vs_fitted_*.png`). Current findings: |skew| < 1 across all models; excess kurtosis 5-7 (heavy tails from a few extreme channels — acceptable at these sample sizes).

## Manuscript .docx Export

`output_plots/betaStim_statistical_tables.docx` (from `betaStim_R_script.R`): residual diagnostics table, model comparison (AIC/BIC), ANOVA tables, EMM pairwise contrasts, Cohen's d effect sizes for models 3a, 3c, 3e, 5a.

`output_plots/betaStim_within_subject_tables.docx` (from `R_compare_control_cond.R`): CL vs PB ANOVA, ecb43e ANOVA, sign-flip permutation tables (exact + MC), null-burst validity tables, percent modulation tables.

CSV outputs in `output_plots/`: `betaStim_clpb_perm_chan_aggregate.csv`, `betaStim_clpb_perm_perchan_bycell.csv`, `betaStim_null_vs_base_perchan.csv`, `betaStim_null_vs_base_aggregate.csv`, percent modulation CSVs.

## Permutation test conventions

**Sign-flip (paired/within-subject)**: Each unit contributes one paired difference; signs are independently exchangeable under H₀. For n ≤ 20 units: exact enumeration (2^n configurations, fast matrix-multiply path for mean statistic). For n > 20: Monte Carlo 10,000 draws. Used for CL vs PB probe pairs and channel-level Null vs Base comparisons.

**Two-sample label shuffling**: Pool trials from both groups, shuffle labels without replacement, recompute difference of group medians. Used for unmatched within-channel tests (e.g., Null vs Base trial-level). Always Monte Carlo.

Holm step-down correction applied for clpb permutation tests across 3 dose bins.

| Test | Unit | n | Method |
|------|------|---|--------|
| Null vs Base aggregated (0b5a2e) | channel | 8 | Exact (2^8 = 256) |
| clpb channel × condition aggregated | channel × phase cell | 13–16 per dose | Exact (2^13–2^16) |
| Null vs Base per channel | trial | ~410 | MC 10k |
| clpb per-cell sign-flip | probe pair | ~30–100 | MC 10k |
