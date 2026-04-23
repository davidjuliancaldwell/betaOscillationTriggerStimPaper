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

### Re-reference / ECO-loading bugs: fixed 2026-04-21

All ECO-struct loading loops that used the narrow `achan ∈ {1,2}` or wider `achan ∈ {1,2,4,6}` "heuristic reload gate" have been replaced with a direct `grp ≠ prev_grp` group-change detector (initialized `prev_grp = -1` to force the first-iteration load). Affected files:

- `peak_extraction/B_ExtractNeuralData_PP_reref.m` (both loops, lines 158 and 228)
- `phase_visualizations/B_phaseCalc_allChans_processed.m` (loop at line 90; also removed debug `chans = 64;` override at line 59, restored `idxVec = [1:8]`, fixed save filename `12samps → 51samps`)
- `manuscript_generate_scripts/plot_EP_goodfit_by_phase.m` (reref loop, line 120)
- `manuscript_generate_scripts/BETA_ExtractNeuralDataCEPscreen.m` (both loops)
- `manuscript_generate_scripts/plot_example_dose_dependent_time_series.m` (reref loop)

Also restored `B_ExtractNeuralData_PP_reref.m` line 19 `for idx = 1:1` (debug override) → `for idx = 1:8` (all subjects).

**Verification**: EP pipeline regenerated for c91479; CSV float-identical (≤1e-12 µV) to pre-fix across all analyzed channels, confirming the bug was latent for the channels that actually entered the analysis. Primary R models (5a, 5a-gf2) reproduce CLAUDE.md numbers exactly. The `phase_data/*.mat` files were NOT regenerated — the Dec 2018 hyak-cluster files remain canonical. With the debug hardcodes removed, regenerating phase data is now safe (just slow — hours for the nonlinear sinusoid fits).

See `phase_data_origin_and_bugs.md` for the full diagnosis, safe-subject matrix, and verification details. Full pre-fix data backup at `data/_backup_20260421_161313/`.

## Key Configuration

- `setup_environment.m` / `Z_Constants.m`: Define subject IDs (SIDS), folder paths for all data directories. Both files set the same variables; `setup_environment.m` is called by `master_script_betaStim.m`.
- Subject IDs: `{'d5cd55', 'c91479', '7dbdec', '9ab7ab', '702d24', 'ecb43e', '0b5a2e', '0b5a2ePlayback'}`
- `0b5a2ePlayback` excluded from main R analysis (CL vs playback comparison only). `702d24` included (1 channel, no convergence issues).
- `external_deps/` contains bundled `CircStat2012a/` and `savitzkyGolay.m`/`sgolayfilt_complete.m`. `setup_environment.m` uses `addpath(genpath(locationsDir))` only — no external paths required.

## Data Directories (under `data/`)

- `ECoG_data/` -- Raw ECoG recordings (gitignored, symlinked from OneDrive: `OneDrive-UCSF/Research/UW_research/betastim/Data-from-xps/ECoG_data`). Contains `{sid}_ECoG.mat` and `{sid}_forBetaPhase.mat` per subject.
- `stim_timing_data/` -- Stimulation timing tables (stage A output)
- `phase_data/` -- Phase calculation results per subject/channel (stage B output). **Dec 2018 hyak-generated; do not regenerate without reading `phase_data_origin_and_bugs.md`.**
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
- Phase classes: **hyperpolarizing (90°) vs depolarizing (270°)** — per Zanos et al. convention (Curr Biol 2018, PIIS0960982218309084.pdf). *Earlier docs had the labels reversed; this is the correct convention for the manuscript.*
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

### Contrast coding (global)

Both `R_analysis_scripts/betaStim_R_script.R` and `R_analysis_scripts/R_compare_control_cond.R` set `options(contrasts = c("contr.sum", "contr.poly"))` at the top. This applies **sum-to-zero** coding to unordered factors (`phaseClass`, `betaLabels`, `numStims` when unordered) and polynomial coding to ordered factors (`numStims_ord`). Required for `anova(fit)` / `anova(fit, type = 3)` Type III tests to yield **marginal** main effects in models with interactions, rather than effects conditional on a reference level.

**Empirical impact** (verified by running the full script once with each coding, same data):
- **Models 5a / 5a-gf2**: *no change at all*. Coefficients, F-values, and p-values are byte-identical between contr.treatment and contr.sum. `numStims_ord` already uses `contr.poly`, sin/cos are continuous, and `betaLabels` is a simple main effect (not in an interaction) — none of the conditions that make Type III coding-dependent are present.
- **Models 3a–3e, 5c**: interaction p-values shift meaningfully (3c/3d/3e: ~0.86 → ~0.98) because under the old default (contr.treatment) the Type III "interaction" was conditional on the reference level of `phaseClass`. Dose main effects essentially unchanged.

contr.sum is retained even though it doesn't move 5a/5a-gf2 — it's the correct default for Type III interpretation and future-proofs any model that gains an unordered-factor interaction.

### Primary models (for manuscript)

**The primary inferential models are 5a and 5a-gf2** (continuous circular phase). They supersede Models 3a-3e (binary 90/270 phaseClass) for these reasons:

1. **Phase is circular, not binary.** Models 3a-3e bin the channel-level circular mean into 90 or 270, which treats phase as a two-level factor and ignores the actual delivered angle. 5a/5a-gf2 decompose phase into `sin(phaseDeg) + cos(phaseDeg)`, using the full circular information (Fisher 1993).
2. **Grouping by measured phase keeps distinct conditions separate.** Models 3a-3e collapse multiple measured phases into one bin per channel. 5a groups by `(sid, phaseDeg_round, numStims, channel)`, preserving both measured phases for multi-phase channels.
3. **5a-gf2 adds phase-fit quality control.** Restricts to bursts where beta was actually present (≥1 conditioning stim with R² > 0.7 and 12-20 Hz), a scientifically motivated filter consistent with the "phase-triggered stimulation requires an ongoing oscillation" hypothesis.
4. **Model 3 series retained as supporting/sensitivity analyses** — they produced the initial result but are now framed as robustness checks, not primary.

Per-cell permutation tests (`R_compare_control_cond.R`, see "Within-subject analyses" below) are reported **alongside** the LMM: the LMM estimates the population-level average effect; the permutation tests describe how that effect is distributed across individual channels.

### Primary: 5a and 5a-gf2 (continuous circular phase)

See "Model 5" subsection below for full specification. Key points (results are the same under contr.treatment or contr.sum — verified):
- **5a**: channel-level phase, phaseVecLength ≥ 0.3, 78 obs, 19 channels. numStims_ord.L p = 0.068 (trend, summary t-test); Type III omnibus dose F p = 0.123; cos_phase p = 0.090.
- **5a-gf2**: good-fit burst restriction + channel-level phase, phaseVecLength ≥ 0.2, 102 obs, 25 channels. numStims_ord.L p = 0.053 (trend, summary t-test); Type III omnibus dose F p = 0.061; cos_phase p = 0.098. Effect size ~8.1 µV.

### Supporting: Models 3a-3e (binary phaseClass, summary-level)

Retained as sensitivity/robustness checks. All summary models collapse to one **median** per (subject × channel × phaseClass × numStims) cell to eliminate pseudoreplication.

```r
# Model 3a (sensitivity): absDiff (baseline median pre-subtracted per channel).
# Random intercepts only. Strong dose effect (p=0.0002) but no random slopes.
# Retained as a robustness check; primary inference is Models 5a / 5a-gf2.
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
- Models 3a-3e: 7 subjects, 32 channels, 126 summary observations (no baseline) or ~157 (with baseline). Grouping by (sid, channel, phaseClass, numStims).
- Models 5a/5a-gf2: 7 subjects, 32 channels, 141 summary obs pre-filter (grouping by phaseDeg_round allows multiple measured phases per channel)
- `doseNum` is a numeric encoding of the dose factor; `dose_linpoly` = `contr.poly(3)[,1]` is an equivalent linear rescaling used with ordinal models
- None of the five primary models are singular

Model comparison (all on median, under contr.sum/contr.poly; Type III ANOVA):

| Model | AIC | Dose p | Interaction p |
|-------|-----|--------|---------------|
| 3a (absDiff, no slopes) | 961 | **0.00037** | 0.966 |
| 3c (ANCOVA, categorical) | 961 | 0.214 | 0.980 |
| 3d (ordinal, .L only) | 959 | 0.214 | 0.980 |
| 3e (numeric, linear) | 961 | 0.109 | 0.981 |

Note: Under the prior contr.treatment default, interaction p-values for 3c/3d/3e were ~0.6–0.86 because the Type III "interaction" was conditional on the reference level of `phaseClass`. After contr.sum, the interactions are essentially null (~0.97–0.98), consistent with the consistent "no dose×phase interaction" conclusion. Dose main effects are essentially unchanged.

### Model 5: Continuous circular phase (sin/cos decomposition)

Replaces binary `phaseClass` (90/270) with `sin(phase)` and `cos(phase)`. The output CSV includes `phaseDeg` (channel-level circular mean in degrees, from `circ_mean` of stims with R² > 0.7), plus quality metrics `phaseVecLength`, `phaseCircStd`, `phaseOmnibusP`. Groups by `(sid, phaseDeg_round, numStims, channel)` — keeps conditions with different measured phases separate. Before phase-quality filters: 141 obs. With default filters (`minPhaseVecLength_5a=0.3`): 78 obs.

```r
# Model 5a (PRIMARY): Ordinal dose × sin/cos + betaLabels, ANCOVA.
# Channel-level phase. phaseVecLength >= 0.3. Non-singular.
fit.sincos.ordinal = lmer(magnitude ~ numStims_ord * (sin_phase + cos_phase) +
  betaLabels + baselineMag_c +
  (1+dose_linpoly|sid) + (1|channel), data=summaryNB_m5)  # 78 obs (r>=0.3)

# Model 5a-gf (diagnostic): Good beta fit restriction, per-burst phase.
# Groups by (sid, numStims, channel) — collapses across delivered phases.
# Intercepts only — random dose slope singular at 96 obs.
# Dose p is anti-conservative without the random slope. Sensitivity check only.
fit.sincos.ordinal.gf = lmer(magnitude ~ numStims_ord * (sin_phase + cos_phase) +
  betaLabels + baselineMag_c +
  (1|sid) + (1|channel), data=summaryNB_gf)  # 96 obs

# Model 5a-gf2 (PRIMARY-complement): Good beta fit + channel-level phase.
# phaseVecLength >= 0.2. Random dose slope non-singular.
# Companion to 5a: adds a scientifically motivated quality filter
# (beta actually present during conditioning) without losing the random-slope
# design. Reported alongside 5a as the primary pair.
fit.sincos.ordinal.gf2 = lmer(magnitude ~ numStims_ord * (sin_phase + cos_phase) +
  betaLabels + baselineMag_c +
  (1+dose_linpoly|sid) + (1|channel), data=summaryNB_gf2)  # 102 obs (r>=0.2)
```

Default filter thresholds: `minPhaseVecLength_5a = 0.3`, `minPhaseVecLength_gf2 = 0.2`, `minGoodBetaPerBurst = 1`, `minBurstVecLength_gf = 0`.

Model 5 results (under contr.sum/contr.poly):
- **5a** (78 obs, 19 ch): numStims_ord.L p=**0.068** (trend, summary t), Type III dose F p=0.123; sin_phase p=0.29, cos_phase p=**0.090** (trend). Effect size ~6.4 µV.
- **5a-gf** (96 obs, intercepts only): numStims_ord.L p≈0.05 — **anti-conservative**. Sensitivity check only.
- **5a-gf2** (102 obs, 25 ch): numStims_ord.L p=**0.053** (marginal, summary t), Type III dose F p=0.061; sin_phase p=0.19, cos_phase p=**0.098** (trend). Effect size ~8.1 µV.

Sensitivity analysis (`output_plots/betaStim_phase_quality_sensitivity.csv`): Models 5a and 5a-gf2 fit at r ∈ {0, 0.1, 0.2, 0.3, 0.4}. Dose effect strengthens monotonically with r.

ANCOVA rationale: `baselineMag_c` (grand-mean-centered baseline median per channel) absorbs ~10-fold between-channel variance. Baseline trials excluded to avoid baseline-by-phase confound. `betaLabels` (1=beta reference channel) included as additive fixed effect — consistently non-significant (p≈0.5–0.7 under contr.sum).

Effect sizes: `d_total = estimate / sqrt(Var_int + E[x²]*Var_slope + Var_channel + Var_resid)`, where E[x²]=1/3 for `contr.poly(3)`.

### Average Marginal Effects (AME) complement — added 2026-04-22

`marginaleffects::avg_comparisons()` added as a complement to the at-90°/at-270° emmeans contrasts for Models 5a and 5a-gf2. Answers a different question: the at-90°/at-270° EMMs evaluate the dose effect at two specific phases (sin=±1, cos=0); the AME averages pairwise dose contrasts over the **observed joint distribution** of sin_phase, cos_phase, baselineMag_c, and betaLabels — the phase-weighted population-level dose effect.

```r
ame_5a_raw <- avg_comparisons(fit.sincos.ordinal,
  variables = list(numStims_ord = "pairwise"))
ame_5a_adj <- hypotheses(ame_5a_raw, multcomp = "single-step")   # max-t FWER
```

Same pattern for `fit.sincos.ordinal.gf2` → `ame_gf2_dose`. Multiplicity: `multcomp = "single-step"` (multcomp::glht max-t) — the closest available analog to the Tukey-Kramer adjustment used by neighboring emmeans `pairs()` output. marginaleffects 0.29 does not expose `"tukey"` directly. Inference uses fixed-effect covariance only (standard behavior for `lmerMod` under `re.form = NA`), consistent with `emmeans` defaults for these models.

AME results (single-step adjusted, under default phase-quality filters):

| Contrast | 5a estimate (µV) | 5a p (adj) | 5a-gf2 estimate (µV) | 5a-gf2 p (adj) |
|----------|------------------|------------|----------------------|----------------|
| [3,4] − [1,2] | 5.18 | 0.221 | 2.28 | 0.801 |
| [5,inf) − [1,2] | **9.09** | **0.052** | **11.01** | **0.048** |
| [5,inf) − [3,4] | 3.91 | 0.422 | **8.73** | **0.041** |

These converge on the at-90°/at-270° estimates and reinforce the dose effect without relying on two specific phase anchors.

**Dose-effect-vs-phase curves** (`marginaleffects::avg_comparisons()` inside a `lapply` over phase angles, paired sin/cos): new figures `output_plots/betaStim_model5a_phase_curve_r30_dose_effect_vs_phase.{png,eps}` and `betaStim_model5a_gf2_phase_curve_r20_dose_effect_vs_phase.{png,eps}`. Y-axis is the pairwise dose contrast (µV); x-axis is delivered phase (0–330° in 30° steps); one line per contrast with pointwise 95% CI ribbon. Complementary to the emmeans-based phase-response curves which show predicted **magnitude** at each phase × dose.

**emmeans plot marginalization fix (2026-04-22)**: The primary single-panel phase-response curves (`p_5a`, `p_5a_gf2`) previously pinned `betaLabels = "0"` (non-beta channels only). They now use `weights = "proportional"` so emmeans marginalizes over `betaLabels` using observed proportions — matches the AME/G-computation semantics and makes the single-panel plots population-level rather than non-beta-specific. Faceted supplementary plots (`p_5a_beta`, `p_5a_gf2_beta`) still show both `betaLabels` levels side-by-side. `baselineMag_c = 0` retained in the `at=...` spec because the covariate is grand-mean centered and enters linearly with no interactions → pinning at 0 is mathematically identical to averaging over the observed distribution (E[baselineMag_c] = 0 by construction). Subtitles reworded from "baselineMag_c = 0" to "baseline at grand mean" for clarity.

### Model 6: Continuous dose spline (EXPLORATORY)

Uses `nCondStims` with natural splines `ns(nCondStims, df=3)`. Without random slopes (singular at 7 subjects), p-values are anti-conservative. Retained for descriptive visualization of the continuous dose-response shape.

### Within-subject analyses (`R_compare_control_cond.R`)

**0b5a2e CL vs Playback**: Probes matched by sequential position (CL probe N ↔ PB probe N). Combined model: `magnitude ~ numStims_ord * (sin_phase + cos_phase) + condition * (sin_phase + cos_phase) + baselineMag_c + (1|channel_raw)`. `condition` ref = PB.

Three phase quantities in play (NOT interchangeable):
1. **`setToDeliverPhase`** — hardware intended target (90° or 270°), channel-independent
2. **`phaseDeg`/`phaseDeg_round`** — channel-level circular mean of sinfit phases. EP channels show spatial offsets up to ~180° from beta reference (expected physiology). Some channels show nearly identical phases for both target conditions — they cannot be treated as providing two distinct phase conditions.
3. **`burstCircMean`/`cl_burst_phase`** — per-burst measured phase from `compute_burst_phase_precision.m` → `{sid}_burst_phase_precision.csv`

For 0b5a2e CL vs PB: matched-pair permutation (sign-flip) is the primary inferential tool; mixed model is used for effect sizes. After `minPhaseVecLength_clpb = 0.2` filter (13 of 16 channel × condition cells kept), none of the aggregated sign-flip permutations clear Holm correction across 3 dose bins ([1,2] raw p=0.205 / Holm 0.614; [3,4] p=0.965 / 1.0; [5,inf) p=0.595 / 1.0). Combined CL+PB ANCOVA model: dose trend p=0.061; condition (CL vs PB) p=0.52; phase at PB p=0.77/0.42 (null as expected, asynchronous); CL-specific phase beyond PB p=0.73/0.76.

**ecb43e targeted vs random**: 36 obs, 4 channels. Combined model shows significant dose effect (numStims_ord p=0.0125) but no phase modulation (sin_phase p=0.86, cos_phase p=0.53), no condType (targeted vs random) main effect (p=0.22), and no phase × condType interactions (p>0.45). Underpowered for phase comparison.

**Percent modulation from baseline**: Grand mean [1,2]=1.3%, [3,4]=3.1%, [5,inf)=4.9%. Strongest responder: c91479 (14-19%).

**Conditioned vs Baseline per-cell permutation** (added 2026-04-14): Companion to Model 5a-gf2. For each `(sid × channel × phaseDeg_round × dose)` cell, a two-sample label-shuffle permutation (median diff, 10k MC) tests whether conditioned probe magnitudes differ from channel-level baseline magnitudes. Bootstrap 95% CIs (2k resamples) provided for the forest plot. **Filters match 5a-gf2 exactly**: channel-level `phaseVecLength ≥ 0.2` AND good-fit burst restriction (`nGoodBeta ≥ 1`, i.e., each conditioned trial's preceding burst had ≥1 stim with R² > 0.7 and frequency 12-20 Hz). Per-cell sample-size minima: ≥10 baseline probes, ≥5 conditioned probes. Baselines are exempt from the good-fit filter (no preceding burst). Primary correction: **BH FDR within (subject × dose)** — channels are nested within subjects, so the correction family is one subject's cells at one dose. Pooled FDR reported as a diagnostic CSV column only (not plotted).

Purpose: the LMM estimates the **average** dose effect; per-cell permutation describes how that average is **distributed** across channels. These answer different questions:
- 5a/5a-gf2 (LMM): is there a population-level dose × phase effect?
- Per-cell perm: in how many individual channels is the effect individually detectable?

Results (101 cells, 7 subjects, with good-fit filter): good-fit filter retains ~46% of conditioned trials. Per-dose condensed summary (within-subject FDR, collapsed across subjects): [1,2] 5/34 uncorr, 4/34 FDR (11.8%); [3,4] 7/33 uncorr, 7/33 FDR (21.2%); [5,inf) 12/34 uncorr, 7/34 FDR (20.6%). Median effect grows with dose: 7.6 → 8.6 → 13.6 µV. Subject-level presence: 2 / 3 / 2 of 7 subjects with ≥1 FDR-sig cell at [1,2] / [3,4] / [5,inf). c91479 dominates (3-4/4 FDR-sig cells at every dose, median effect ~48 µV); 9ab7ab shows a dose gradient (1/4 → 2/4 → 3/4 cells); ecb43e has 1/5 at [3,4]. 0b5a2e's 13-cell FDR family is harsh — the good-fit filter makes its [5,inf) effect stronger (median 24 µV, 4 uncorrected-sig cells) but within-subject FDR still rejects at 13 tests. Phase results (from LMM): cos_phase is trend-level in both 5a (p=0.090) and 5a-gf2 (p=0.098); sin_phase non-significant; no dose × phase interactions (all p>0.23). Outputs in `output_plots/`: `betaStim_cond_vs_base_perchan.csv` (master, includes within-subject `perm_q` and pooled `perm_q_pooled` columns), `betaStim_cond_vs_base_per_subject.csv`, `betaStim_cond_vs_base_subject_presence.csv`, `betaStim_cond_vs_base_summary.csv` (per-dose condensed counts), `betaStim_cond_vs_base_pooled_summary.csv`, `betaStim_cond_vs_base_forest.png/.eps`. Forest plot sorts rows by measured phase (0° at top); beta-trigger channels highlighted with pink y-axis labels; dots colored by significance category (ns / p<0.05 uncorr / FDR q<0.05). All red cells are also uncorrected-significant by construction (BH q ≥ p always).

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

`output_plots/betaStim_statistical_tables.docx` (from `betaStim_R_script.R`): residual diagnostics table, model comparison (AIC/BIC), ANOVA tables, EMM pairwise contrasts, Cohen's d effect sizes for models 3a, 3c, 3e, 5a, and — added 2026-04-22 — AME pairwise dose contrasts for Models 5a and 5a-gf2 (single-step max-t corrected, phase-marginalized population-level complement to the at-90°/at-270° EMM contrasts).

`output_plots/betaStim_within_subject_tables.docx` (from `R_compare_control_cond.R`): CL vs PB ANOVA, ecb43e ANOVA, sign-flip permutation tables (exact + MC), null-burst validity tables, percent modulation tables, conditioned-vs-baseline per-cell permutation tables (across-subject summary, per-subject breakdown, pooled-FDR sensitivity).

CSV outputs in `output_plots/`: `betaStim_clpb_perm_chan_aggregate.csv`, `betaStim_clpb_perm_perchan_bycell.csv`, `betaStim_null_vs_base_perchan.csv`, `betaStim_null_vs_base_aggregate.csv`, `betaStim_cond_vs_base_perchan.csv`, `betaStim_cond_vs_base_per_subject.csv`, `betaStim_cond_vs_base_subject_presence.csv`, `betaStim_cond_vs_base_summary.csv`, `betaStim_cond_vs_base_pooled_summary.csv`, percent modulation CSVs.

## Permutation test conventions

**Sign-flip (paired/within-subject)**: Each unit contributes one paired difference; signs are independently exchangeable under H₀. For n ≤ 20 units: exact enumeration (2^n configurations, fast matrix-multiply path for mean statistic). For n > 20: Monte Carlo 10,000 draws. Used for CL vs PB probe pairs and channel-level Null vs Base comparisons.

**Two-sample label shuffling**: Pool trials from both groups, shuffle labels without replacement, recompute difference of group medians. Used for unmatched within-channel tests (Null vs Base trial-level, Conditioned vs Baseline per-cell). Always Monte Carlo.

**Bootstrap 95% CIs**: Nonparametric resampling within each group (BCa not used; percentile intervals from 2,000 bootstrap resamples). Used to visualize per-cell effect uncertainty on forest plots. Complements the permutation p-value, which uses a different null hypothesis (exchangeability) than the bootstrap (observed distribution).

**FDR correction strategies**:
- **Holm step-down**: clpb aggregated tests across 3 dose bins.
- **BH within (subject × dose)**: Conditioned vs Baseline per-cell. Family = one subject's cells at one dose. Respects the nested design (channels within subjects).
- **BH pooled within dose**: Conditioned vs Baseline, reported as sensitivity alongside the within-subject FDR.

| Test | Unit | n | Method | Correction |
|------|------|---|--------|------------|
| Null vs Base aggregated (0b5a2e) | channel | 8 | Exact sign-flip (2^8 = 256) | none |
| clpb channel × condition aggregated | channel × phase cell | 13–16 per dose | Exact sign-flip (2^13–2^16) | Holm across 3 doses |
| Null vs Base per channel | trial | ~410 | Two-sample MC 10k | descriptive |
| clpb per-cell sign-flip | probe pair | ~30–100 | Sign-flip MC 10k | BH within dose (3 families of ~13–16 cells) |
| Conditioned vs Baseline per-cell | trial | ~50–200 cond vs ~40–60 base | Two-sample MC 10k | BH within (subj × dose); pooled within dose as sensitivity |
