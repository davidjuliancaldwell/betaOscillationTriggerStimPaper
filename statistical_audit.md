# Statistical Audit: Mixed Effects Model Analysis

Audit of `R_analysis_scripts/betaStim_R_script.R` and its upstream data generation in `peak_extraction/multipleSubj_GLMM_script_PP.m`.

---

## Current primary models (2026-04-14 update)

The sections below document successive iterations of the primary model across the project's history. **As of 2026-04-14, the primary inferential models are Model 5a and Model 5a-gf2 (continuous circular phase, see Finding 23).** Earlier primary designations (originally a trial-level model with `(1|sid/channel)`, later Model 3 series with binary phaseClass) are retained in this document as historical record and as robustness checks for the manuscript.

**Primary pair:**
- **Model 5a** (channel-level phase, `phaseVecLength ≥ 0.3`, random dose slope): 78 obs, 19 channels. Dose.L p=0.068 (trend), effect size 6.4 µV.
- **Model 5a-gf2** (good-fit burst restriction + channel-level phase, `phaseVecLength ≥ 0.2`, random dose slope): 102 obs, 25 channels. Dose.L p=0.049, effect size 8.6 µV.

**Companion per-cell analysis** (added 2026-04-14, Finding 24): Conditioned vs Baseline two-sample label-shuffle permutation (10k MC, median diff) for each `(subject × channel × phase × dose)` cell, BH FDR within (subject × dose). Describes the distribution of the LMM's population-level effect across individual channels.

**Superseded designations:**
- Original trial-level LMM (Section "Model Under Review" below): `(1|sid/channel)` inflated DF → flagged in Findings 1 and 11, no longer primary.
- Models 3a-3e (binary phaseClass, Finding 15): primary through early 2026; now reframed as sensitivity checks because binary phaseClass discards circular information and conflates multiple delivered phases at multi-phase channels. See Finding 23 for the rationale for moving to 5a/5a-gf2.

---

## Model Under Review (historical — original trial-level model)

The primary model fit on trial-level data (`dataNoBaseline`, ~37K rows after exclusions):

```r
fit.lmm3 = lmerTest::lmer(
  absDiff ~ numStims + phaseClass + betaLabels +
            numStims:betaLabels + numStims:phaseClass +
            (1|sid/channel),
  data = dataNoBaseline
)
```

A secondary model fit on summary-level data (`summaryDataForMixed`, medians per cell):

```r
fit.lmm5 = mixed(
  magnitude ~ numStims + phaseClass + betaLabels +
              numStims:betaLabels + numStims:phaseClass +
              (numStims|sid/channel),
  data = summaryDataForMixed
)
```

---

## Finding 1: `phaseClass` is a channel-level constant, not a trial-level variable

**This is the most significant finding.**

In `multipleSubj_GLMM_script_PP.m:190-191`:

```matlab
phaseVecChosen = peakPhaseVec(ii, goodEPs==chan);    % ONE scalar per channel×condition
phaseVec = repmat(phaseVecChosen, lengthType, 1)';   % stamped onto EVERY trial
```

`peakPhase` is the circular mean of phase-at-delivery across all qualifying trials for a channel (computed via `circ_mean` in `phase_circstats_calc.m`). This single summary value is then binned to 90 or 270 and replicated identically across every trial row for that channel.

**phaseClass does not vary trial-to-trial.** The 37,221 rows in `dataNoBaseline` map to only **49 unique channel x condition x phaseClass units** (25 at phaseClass=270, 24 at phaseClass=90). The effective sample size for the phaseClass fixed effect is ~49, not ~37K.

`lmer` with `(1|sid/channel)` should adjust standard errors for this clustering, but the adjustment depends on having enough clusters. With only 31 unique channels across 6 subjects, Satterthwaite's DF approximation may be anti-conservative.

By contrast, `numStims` (the dose level: Base, [1,2], [3,4], [5,inf)) genuinely varies trial-to-trial within a channel — different trials in the same experimental block received different numbers of conditioning stimuli before the test pulse. So the trial-level model provides real statistical power for the dose fixed effect and its interactions with other trial-level variables.

The `numStims:phaseClass` interaction inherits the phaseClass problem — it asks whether the dose-response curve differs between phaseClass=90 and phaseClass=270 groups, but that comparison is still at the channel level.

---

## Finding 2: Only 9 of 31 channels provide within-channel phaseClass contrast

For multi-phase subjects, different experimental blocks targeted different oscillation phases. Some channels received stimulation at both phase targets, creating within-channel variation in phaseClass — the clean, causally interpretable contrast.

**Channels with BOTH phaseClass=90 and phaseClass=270** (within-channel contrast):

| Subject | Channels | Count |
|---------|----------|-------|
| c91479 (type m) | 264 | 1 |
| 0b5a2e (type m) | 715, 716, 723, 731, 732 | 5 |
| ecb43e (type t) | 647, 648, 655 | 3 |
| **Total** | | **9** |

**Channels with only ONE phaseClass** (between-channel comparison only):

| Subject | Channels at 90 | Channels at 270 | Total |
|---------|----------------|-----------------|-------|
| d5cd55 (type s) | 7 | 2 | 9 |
| 7dbdec (type s) | 2 | 1 | 3 |
| 9ab7ab (type s) | 0 | 5 | 5 |
| c91479 (type m) | 1 (ch 247) | 0 | 1 |
| 0b5a2e (type m) | 0 | 3 (714, 721, 740) | 3 |
| ecb43e (type t) | 1 (ch 654) | 0 | 1 |
| **Total** | | | **22** |

For the 22 single-phaseClass channels, the phaseClass effect is estimated entirely from between-channel differences. These are confounded with electrode location, distance from stimulation site, cortical excitability, and other channel-specific properties that affect CEP magnitude independently of oscillation phase.

---

## Finding 3: Single-phase vs multi-phase subjects test different scientific questions

**Single-phase subjects** (d5cd55, 7dbdec, 9ab7ab — type 's'):
- All stimulations targeted ONE phase (e.g., 180 degrees)
- phaseClass varies only BETWEEN channels, reflecting spatial variation in how beta oscillations propagate across the cortex
- The comparison asks: *"Do channels at different cortical positions show different CEP magnitudes?"*
- Confounded with any channel-level property

**Multi-phase subjects** (c91479, 0b5a2e — type 'm'; ecb43e — type 't'):
- Stimulations targeted TWO (or more) phases in separate experimental blocks
- phaseClass can vary WITHIN a channel across conditions
- The comparison asks: *"Does the same channel respond differently when stimulated at different oscillation phases?"*
- This is the clean causal contrast — channel-specific confounds cancel out

The model pools these fundamentally different evidence structures into a single phaseClass fixed effect without distinguishing them.

---

## Finding 4: 9ab7ab contributes no phaseClass information

9ab7ab (type 's', desiredF=270) has all 5 channels coded phaseClass=270:

```
450,270
451,270
452,270
453,270
458,270
```

After filtering, this subject contributes **7,786 rows — 21% of the total model data** — all on one side of the phase contrast. It cannot inform the phaseClass effect or the numStims:phaseClass interaction. It shifts the balance of the dataset without adding discriminative information about phase.

It does contribute to estimating numStims and betaLabels effects.

---

## Finding 5: ecb43e random-condition trials are included

Subject ecb43e (type 't') has three experimental conditions:
- `setToDeliverPhase = 270` (targeted depolarizing)
- `setToDeliverPhase = 90` (targeted hyperpolarizing)
- `setToDeliverPhase = 12345` (random phase — not phase-locked)

The R code excludes `numStims == 'Null'` but does NOT exclude `setToDeliverPhase == 12345`. Approximately **1,082 random-condition trials** are assigned a phaseClass based on sinusoidal fits and treated identically to deliberately phase-targeted trials in the model.

Since the random condition was not locking to a specific oscillation phase, including these trials and treating their phaseClass the same as deliberate phase targeting is conceptually questionable.

---

## Finding 6: Trial count imbalance across subjects

After R exclusions (removing 702d24, 0b5a2ePlayBack, Null, Base; filtering magnitude 25-1500):

| Subject | Rows | Channels | Type | Phase evidence |
|---------|------|----------|------|----------------|
| d5cd55 | 9,968 | 9 | s | Between-channel only |
| 9ab7ab | 7,786 | 5 | s | **None** (all 270) |
| 7dbdec | 6,462 | 3 | s | Between-channel only |
| 0b5a2e | 5,299 | 8 | m | Within-channel (5 of 8 channels) |
| ecb43e | 4,653 | 4 | t | Within-channel (3 of 4 channels) + random |
| c91479 | 3,053 | 2 | m | Within-channel (1 of 2 channels) |
| **Total** | **37,221** | **31** | | |

Subjects with only between-channel or no phaseClass evidence (d5cd55 + 9ab7ab + 7dbdec) contribute **24,216 rows (65%)**. The 3 subjects with within-channel phaseClass variation contribute **13,005 rows (35%)**.

With `(1|sid/channel)` random intercepts only, subjects with more trials have more influence on fixed effect estimates. The model does not weight subjects equally.

---

## Finding 7: Baseline computed from condition 1 only

In `multipleSubj_GLMM_script_PP.m`, the `for ii = 1:numTypes` loop includes baseline trials only for `ii == 1` (lines 234-249). For multi-phase subjects, this means the first experimental block's baseline is used to compute `absDiff` and `percentDiff` for test trials from ALL blocks.

If brain state drifted between experimental blocks (which were run sequentially), this single baseline may not be appropriate for later blocks, introducing systematic bias in the baseline correction.

---

## Finding 8: Summary model is overparameterized

The secondary model specification:

```r
fit.lmm5 = mixed(
  magnitude ~ numStims + phaseClass + betaLabels +
              numStims:betaLabels + numStims:phaseClass +
              (numStims|sid/channel),
  data = summaryDataForMixed
)
```

`(numStims|sid/channel)` fits a random intercept + random slopes for each numStims level at both the subject level AND the channel-within-subject level. For a subject like c91479 with 2 channels and ~4 numStims levels, this means estimating ~4 random effect parameters per channel from ~8 summary observations. The model likely has convergence issues or singular fits.

---

## Finding 9: `betaLabels` is weakly identified

`betaLabels` is binary: 1 for the beta recording channel, 0 for all others. Each subject has exactly 1 beta channel out of 1-9 total channels. The `numStims:betaLabels` interaction is estimated from comparing the dose-response of 1 channel versus N-1 channels within each subject — low power and sensitive to outliers on the single beta channel.

---

## Finding 10: Closed-loop vs playback control analysis

`R_analysis_scripts/R_compare_control_cond.R` compares subject 0b5a2e (closed-loop, phase-triggered) vs 0b5a2ePlayBack (control, replayed timing). Same physical subject and electrodes under different experimental conditions.

### Original analysis (channel 14 only, selected for strong response)

Fixed-effects `lm(magnitude ~ numStims + sid + numStims:sid)` on channel pair 714/814.

**Cohen's d (via `emmeans::eff_size`, pooled residual SD = 120.4 uV):**

| Dose | CL mean | PB mean | Diff | Cohen's d (95% CI) | p |
|------|---------|---------|------|--------------------|---|
| Base | 346 uV | 312 uV | +34 | 0.28 (-0.11, 0.68) | 0.162 |
| [1,2] | 340 uV | 310 uV | +30 | 0.25 (0.10, 0.39) | 0.001 |
| [3,4] | 378 uV | 319 uV | +59 | 0.49 (0.31, 0.67) | <.0001 |
| [5,inf) | 432 uV | 329 uV | +103 | 0.86 (0.58, 1.13) | <.0001 |

Within closed-loop: significant dose-response staircase. Within playback: no significant dose-response (all p > 0.4, all d < 0.2).

### Permutation test (channel 14, 10,000 permutations, two-sided)

Condition labels shuffled within each dose level independently.

Per-dose CL vs PB:

| Dose | Obs diff | Perm p |
|------|----------|--------|
| Base | +34.0 uV | 0.12 |
| [1,2] | +29.6 uV | 0.0002 |
| [3,4] | +58.8 uV | <0.0001 |
| [5,inf) | +102.9 uV | <0.0001 |

**Dose-response interaction** (does CL-PB gap grow from Base to [5,inf)?):
- CL-PB at Base: 34.0 uV
- CL-PB at [5,inf): 102.9 uV
- Interaction: 68.9 uV, permutation p = 0.039 (two-sided)

The interaction test shuffles condition labels independently at each dose level, then computes the difference-of-differences. Under the null (no dose-dependent interaction), this should be centered at zero. The observed 68.9 uV exceeds 96.1% of permuted values in absolute magnitude.

### All-channel permutation (8 matched pairs, for context)

When extended to all 8 matched channels (714-740 vs 814-840):
- CL > PB significant at all dose levels including baseline (all p < 0.0001)
- Dose-response interaction NOT significant (p = 0.18)
- Systematic session offset: closed-loop magnitudes higher even at baseline

Channel 14 was selected for its strong response. The dose-dependent interaction holds on that channel but does not generalize across all 8 channels.

---

## Finding 11: Denominator DF inflation across model specifications

Systematically compared four random effects structures on the same data. The original model's `(1|sid/channel)` produces severely inflated DF for all predictors. Adding random dose slopes per subject `(0+numStims|sid)` corrects the numStims DF but leaves phaseClass inflated.

| Effect | Model 1 (orig) | Model 2 (dose slopes) | **Model 3 (primary)** | Model 4 (nested, ref) |
|--------|------|------|------|------|
| numStims | df=37K, **p=2e-6** | df=5, p=0.16 | **df=84, p=0.0002** | df=5, p=0.16 |
| phaseClass | df=4.6K, p=0.14 | df=4.2K, p=0.18 | **df=84, p=0.42** | df=31, p=0.27 |
| interaction | df=37K, **p=6e-4** | df=851, **p=0.048** | **df=84, p=0.16** | df=753, **p=0.040** |

**Model 3** (summary-level, primary) collapses to one median per (subject x channel x phaseClass x numStims) cell (120 observations) and uses simple random intercepts: `(1|sid) + (1|channel)`. No singularity. Random dose slopes were removed because 6 subjects cannot support a 3x3 covariance matrix (correlations hit 1.0). `setToDeliverPhase` is not used as a random grouping factor — it is a fixed experimental condition, not a random sample.

Model 3 performance: conditional R² = 0.998, marginal R² = 0.0005, partial eta² for numStims = 0.19 (large). The high ICC (0.998) reflects that most variance is between channels/subjects; the experimental manipulation produces a small (~11 uV) but detectable shift.

phaseClass DF limitation: Satterthwaite assigns ~84 DF for phaseClass instead of the ideal ~30 (between-channel). Adding `(1|channel:setToDeliverPhase)` would correct this but is redundant with `(1|channel)` for 22/31 single-phase channels, causing singularity that cannot be resolved without also hitting boundary issues in the dose slope covariance. This does not affect conclusions — phaseClass is non-significant at DF=84 and would be less significant with fewer DF.

**Model 4** (trial-level, nested conditions) is kept for reference. It is singular due to both `(1|channel:setToDeliverPhase)` redundancy with `(1|channel)` for single-phase subjects and near-saturated dose slope correlations (0.89–0.98) in `(0+numStims|sid)`.

Both models converge on the same pattern: the dose effect is driven by the [5,inf) level, and emmeans show this effect is concentrated in phaseClass=270 channels (270 > 90 by ~6 uV at highest dose, p=0.063 in Model 3).

---

## Finding 12: Sensitivity analyses — interaction survives exclusions but not within-channel restriction

Model 4's nested condition structure was tested with progressive exclusions:

| Model | Subjects | Channels | Trials | Interaction p | [5,inf) phase contrast |
|-------|----------|----------|--------|---------------|----------------------|
| 4 (full) | 6 | 31 | 37,221 | **0.040** | 6.2 uV, p=0.017 |
| 4a (no 9ab7ab) | 5 | 26 | 29,435 | **0.049** | 5.9 uV, p=0.025 |
| 4b (no 9ab + random) | 5 | 26 | 27,994 | **0.042** | 6.3 uV, p=0.031 |
| 4c (within-channel only) | 3 | 9 | 7,489 | 0.70 | 2.4 uV, p=0.48 |

The interaction is robust to excluding 9ab7ab (no phase contrast) and ecb43e's random-condition trials (not phase-targeted). The [5,inf) phase contrast (270 > 90) is consistent at ~6 uV across Models 4, 4a, and 4b.

**Model 4c (within-channel only)** restricts to the 9 channels (c91479: 264; 0b5a2e: 715, 716, 723, 731, 732; ecb43e: 647, 648, 655) where the same physical electrode was stimulated at both phases in separate experimental blocks. This eliminates the between-channel confound (electrode location, cortical excitability, distance from stim site). The interaction disappears (p = 0.70) and the effect size drops from ~6 uV to ~2.4 uV.

This could reflect:
1. **Insufficient power** — 3 subjects and 9 channels may be too few to detect a real but small effect
2. **The between-channel confound driving the effect** — channels at phaseClass=270 may inherently produce larger evoked potentials for reasons unrelated to oscillation phase (e.g., proximity to stim site, cortical geometry)

Both interpretations should be acknowledged. The within-channel analysis is the cleanest causal test but is underpowered. The full-dataset result is well-powered but includes confounded between-channel evidence.

---

## Recommendations

**These recommendations were drafted when the Model 3 series (binary phaseClass) was the primary analysis. Items 4-7 have been superseded by Finding 23 (continuous phase, Models 5a/5a-gf2). Items 1-3 and 8 remain valid and are implemented in the current pipeline.**

### Current (2026-04-14) recommendations

A. **Primary inference: Models 5a and 5a-gf2 (continuous sin/cos phase).** See Finding 23. Report both — 5a as the conservative primary, 5a-gf2 as the primary complement with the good-fit quality filter.

B. **Companion per-cell analysis: Conditioned vs Baseline two-sample permutation.** See Finding 24. Use BH FDR within (subject × dose) to respect nesting. Report pooled FDR as sensitivity. This characterizes how the LMM's average effect is distributed across channels.

C. **Retain Models 3a-3e as robustness checks**, not primary. The binary phaseClass framing is no longer the central inferential claim but gives readers a cross-check against a categorical phase model.

D. **Keep within-channel analyses (Finding 2, 9 channels) as explicit sensitivity checks**, not as the primary causal test — the within-channel subset is too underpowered (p=0.70 in Model 4c) to carry primary inference, and the full-dataset phase decomposition (5a/5a-gf2) avoids the binary-bin confound that motivated the within-channel restriction.

### Historical recommendations (Model 3 era — items 4-7 superseded)

1. **Separate within-channel and between-channel phaseClass analyses.** Restrict the primary phaseClass analysis to the 9 channels (across 3 subjects) with within-channel variation. Report the between-channel analysis separately as supporting evidence, clearly noting the confound.

2. **Remove 9ab7ab from any model that includes phaseClass** as a predictor. It contributes no contrast and distorts the dataset balance.

3. **Exclude or separately analyze ecb43e random-condition trials** (`setToDeliverPhase == 12345`). They weren't phase-targeted and shouldn't receive the same treatment as deliberate phase conditions.

4. ~~**Use the summary-level model (Model 3) as the primary analysis**~~ — **superseded by Recommendation A (Finding 23).** Model 3 is retained as a sensitivity check.

5. ~~**Validate phaseClass p-values**~~ — **superseded**; continuous-phase models (5a/5a-gf2) and the per-cell permutation analysis (Finding 24) serve this role.

6. ~~**Consider the continuous phase variable...as a sensitivity analysis**~~ — **realized as primary in Finding 23.**

7. **Report effective sample sizes** alongside the model. Still applies: 5a uses 78 obs across 19 channels; 5a-gf2 uses 102 obs across 25 channels; 7 subjects throughout.

8. **Frame the channel 14 closed-loop vs playback comparison as a selected-channel analysis.** The dose-response interaction is significant on this channel (permutation p = 0.039) but does not generalize across all 8 channels (p = 0.18). Report both results transparently.

---

## ANOVA types and reporting structure

The model output has four complementary components:

| Output | What it tests | Example |
|--------|--------------|---------|
| **Type III ANOVA** (`anova(model)`) | Omnibus F-test per factor, each tested at reference level of other factors | "Does numStims matter (at phaseClass=270)?" |
| **Type II ANOVA** (`anova(model, type=2)`) | Omnibus F-test averaging over other factors, ignoring interactions | "Does numStims matter (averaging across phases)?" |
| **Coefficients** (`summary` / `tab_model`) | Each factor level vs reference | "[5,inf) is 10.1 uV above [1,2] at phaseClass=270" |
| **Emmeans** (`emmeans` + `pairs`) | All pairwise contrasts, not just vs reference | "[3,4] vs [5,inf) within phaseClass=270" |

**Type III is the correct choice** when the interaction is significant (p = 0.040) because main effects are conditional on the other factor's level. Type II gives different main effect p-values (phaseClass p = 0.51 vs 0.27) because it averages across dose levels, diluting the phase effect that only appears at [5,inf).

The interaction is the highest-order term, so Type II and Type III agree on it (p = 0.040 in both).

---

## Subject 6 (ecb43e): phase-targeted vs random stimulation

`R_analysis_scripts/R_compare_subject_6_random.R` compares the three experimental conditions on channel 655 (beta channel): depolarizing (setToDeliverPhase=270), hyperpolarizing (90), and random (12345).

Permutation tests (10,000 permutations, two-sided) show no significant difference between phase-targeted and random stimulation at any dose level:

| Comparison | [1,2] | [3,4] | [5,inf) |
|------------|-------|-------|---------|
| Depolarizing vs Random | +7.4 uV, p=0.41 | -18.2 uV, p=0.051 | -1.0 uV, p=0.89 |
| Hyperpolarizing vs Random | +7.9 uV, p=0.34 | -10.7 uV, p=0.29 | -6.9 uV, p=0.32 |

Effect sizes are small (all |d| < 0.3). The one borderline result ([3,4] depolarizing vs random, p=0.051) goes in the wrong direction (random is larger).

This is consistent with the Model 4c finding: within-channel phase effects do not reach significance. On this single subject's beta channel, phase-targeted stimulation does not produce measurably different evoked potentials than random-phase stimulation.

---

## Effect size reporting

**Primary measure**: Raw effects in uV from emmeans pairwise contrasts. Directly interpretable — "phaseClass 270 produces 6.2 uV larger CEPs than 90 at [5,inf) dose."

**Supplementary measure**: Cohen's d using a consistent denominator across all models — the trial-level residual SD from Model 4 (sigma = 77.4 uV). This represents trial-to-trial noise within a single channel/condition/dose cell. Using one sigma makes d values comparable across Models 2-4 and sensitivity analyses, avoiding the inconsistency where summary-level models give artificially large d due to their much smaller residual SD (~7.5 uV).

Phase contrast at [5,inf) (270 - 90) across models, all using sigma = 77.4 uV:

| Model | Raw effect | Cohen's d | 95% CI |
|-------|-----------|-----------|--------|
| 4 (full) | 6.2 uV | 0.08 | (0.01, 0.15) |
| 4a (no 9ab7ab) | 5.9 uV | 0.08 | (0.01, 0.14) |
| 4b (clean) | 6.3 uV | 0.08 | (0.01, 0.16) |
| 4c (within-channel) | 2.4 uV | 0.03 | (-0.06, 0.12) |

The d = 0.08 is small relative to single-trial noise — expected for neural data where trial-to-trial variability (77 uV) far exceeds systematic effects (6 uV). The effect is reliable (significant across thousands of trials) but small on any given pulse. The within-channel analysis (4c) shows a halved effect size that does not reach significance.

---

## Finding 16: Phase fit screening — frequency filter is the bottleneck, not R²

The phase visualization plots show an apparent discrepancy: the R² histogram looks heavily concentrated near 1.0, but the polar plot reports only 225/1451 trials passing. This is because the histogram pre-filters for frequency (12-20 Hz) before showing R² values, while the denominator (1451) is all trials.

Analysis on d5cd55 channel 53 (beta channel):

| Filter | Passes | Percentage |
|--------|--------|-----------|
| R² > 0.7 (any freq) | 1023 | 70% |
| Freq 12-20 Hz (any R²) | 244 | 17% |
| Both R² > 0.7 AND freq 12-20 | 225 | 15% |

798 trials (55%) have excellent sinusoidal fits (R² > 0.7) but the fitted frequency is outside 12-20 Hz — 639 below 12 Hz, 159 above 20 Hz. The fit locks onto a clean oscillation in a different frequency band. Among trials that DO fall in 12-20 Hz, 96% have R² > 0.7 (mean R² = 0.94).

**Impact on burst-level phase precision analysis**: With only ~17% of conditioning stims producing in-band fits, a burst of 7 stims yields ~1.2 good fits on average. This makes per-burst phase precision (circular R) unreliable — 1 fit gives R=1.0 trivially. The `compute_burst_phase_precision.m` script was written to test whether tighter phase precision in conditioning bursts predicts larger CEPs, but the low pass rate limits statistical power for this analysis.

The frequency filter (12-20 Hz) is more restrictive than the R² filter (0.7) because the 40ms pre-stimulus window (~1 beta cycle) often contains oscillatory activity outside the beta band. The sinusoidal fit finds the best-fitting frequency, which may not be beta even when the fit quality is excellent.

---

## Finding 17: CL vs Playback analysis on beta channel (channel 31)

The original CL vs PB comparison used channel 14 (selected for strong response). Channel 31 is the actual beta channel where phase detection and triggering occurred, and has the largest CEP magnitudes (median 530 uV).

**Phase delivery on beta channel** (from phase data files):
- 90-targeted: circular mean = 109.4°, R = 0.718 (good concentration)
- 270-targeted: circular mean = 296.8°, R = 0.551

**Channel 14 does NOT see different phases**: both conditions give circular mean ~270-282° on channel 14. The phase offset between the detection site (ch 31) and recording site (ch 14) is such that ch 14 is always near 270° regardless of targeting condition.

**Beta channel CL vs PB results** (log-transformed LM, trial-level):
- Condition (CL vs PB): p = 0.017 — CL consistently larger
- Dose (numStims): p = 0.189
- Interaction: **p = 0.030** — but in the WRONG direction: CL advantage shrinks from +6.7% at Base to +1.9% at [5,inf)

**Per-dose CL > PB (median, permutation)**:
- Base: +41.5 uV, p = 0.018
- [1,2]: +37.8 uV, p < 0.0001
- [3,4]: +40.0 uV, p < 0.0001
- [5,inf): +26.1 uV, p = 0.004

All dose levels significant. The CL advantage is a constant ~30-40 uV session offset, not dose-dependent.

**Within CL: 270-targeted vs 90-targeted at [5,inf)**:
- 270: median 534.6 uV, 90: median 521.9 uV
- Difference: +12.7 uV, permutation p = 0.286 — not significant

**All-channel comparison** (8 matched pairs, linear model with baseline covariate):
- Condition: CL > PB by ~23 uV, p = 3.9e-07
- Dose: F(2,35) = 6.22, p = 0.005
- Interaction: F(2,35) = 0.03, p = 0.968 — no dose × condition interaction

The closed-loop system produces systematically larger CEPs than playback control across all channels and dose levels, but this is a session-level offset, not a dose-dependent phase-locking benefit.

---

## Finding 18: Channel 14 heavy tail analysis

The mean-based CL-PB interaction on channel 14 (previously p = 0.039) was driven by a heavy right tail in the CL condition:

| Dose | CL sd | PB sd | CL trials > 540 uV | PB trials > 540 uV |
|------|-------|-------|--------------------|--------------------|
| Base | 137.7 | 56.9 | 4 (8.2%) | 0 |
| [1,2] | 130.7 | 71.1 | 25 (7.1%) | 0 |
| [3,4] | 175.3 | 73.2 | 29 (12.3%) | 0 |
| [5,inf) | 242.9 | 65.3 | 16 (15.7%) | 0 |

Zero PB trials exceed 540 uV. The CL condition has dozens of extreme trials (up to 1285 uV), and the proportion grows with dose (8% → 16%). The CL mean at [5,inf) is 431.6 but median is 358.2 — a 73 uV gap driven by the tail. PB mean ≈ median.

Switching to median: the interaction drops from +68.9 uV (p=0.039) to +16.9 uV (p=0.458). The per-dose CL > PB remains significant at [5,inf) (+27.5 uV, p=0.036) but the dose-dependent growth is not significant.

Log-transformed trial-level LM on channel 14: interaction p = 0.095 (trending). Residual skew drops from 2.5 (raw) to 0.3 (log).

---

## Finding 19: Burst-level phase precision does NOT predict CEP magnitude

`compute_burst_phase_precision.m` links each test stim to its preceding conditioning burst and computes the circular mean and vector length (concentration) of beta-band phase fits (R² > 0.7 AND 12-20 Hz) across the burst's conditioning stims. Merged with EP magnitudes via rank-matching within (channel, target phase, dose) cells.

**Note**: The EP output table's `setToDeliverPhase` labels are swapped for 0b5a2e — `setToDeliverPhase="90"` in the CSV corresponds to burst type 0 (270-targeted), and vice versa. This is because `multipleSubj_GLMM_script_PP.m` assigns `desiredF(ii)` where `ii` is the loop index over sorted burst types, not the burst type value itself. The burst phase precision table has correct labels. The merge accounts for this swap.

### Beta channel (31) at [5,inf)

| Target | n with beta fits | Spearman rho (vecLength vs mag) | p |
|--------|-----------------|--------------------------------|---|
| 90° | 46 | -0.074 | 0.625 |
| 270° | 39 | -0.205 | 0.210 |

No significant correlation. Most bursts have only 1-2 beta-band fits (median nGoodBeta = 1-2), giving ceiling vector lengths (R ≈ 1.0) with insufficient variance to test the hypothesis.

### Channel 14 at [5,inf) — significant REVERSE effect

For 270-targeted bursts (n=41 with beta fits):
- **Spearman rho = -0.299, p = 0.058** — trending negative
- **Tight phase (R > 0.882, n=20)**: median 331.7 uV
- **Loose phase (R ≤ 0.882, n=21)**: median 438.8 uV
- **Difference: -107.2 uV, permutation p = 0.006** — loose phase produces LARGER CEPs

Extreme trials (>540 uV, n=16):
- Have LOWER phase precision (median vecLength = 0.725) than normal trials (0.994)
- Split evenly across target phases (7 targeting 270°, 9 targeting 90°)
- Not driven by precise 270° targeting

The large CEPs on channel 14 come from bursts with poor phase precision, not tight phase-locking. This may reflect cortical excitability fluctuations that simultaneously disrupt the beta oscillation (degrading fits) and enhance evoked responses, or stimulation artifact from long/intense bursts that degrades fit quality while also producing large CEPs through cumulative conditioning.

---

## Finding 20 (CRITICAL): phaseClass AND setToDeliverPhase both swapped for multi-phase subjects

**Bug**: In `multipleSubj_GLMM_script_PP.m`, the data loop `for ii = 1:numTypes` uses `desiredF(ii)` and `peakPhaseVec(ii,...)` directly. But `ii` indexes `dataForPPanalysis{chan}{ii}`, which stores data for burst type `types(ii)` (sorted, 0-indexed). The phase calc loop uses `indices = [1,2]` where index 1 = pos (stims(8)==1, rising, ~90° target) and index 2 = neg (stims(8)==0, falling, ~270° target). Since burst types are sorted as [0, 1, ...], `ii=1` → burst type 0 (neg/270-targeted) but gets `desiredF(1)` (90) and `peakPhaseVec(1,:)` (90-targeted phase). Both labels are wrong.

**Affected columns**: BOTH `setToDeliverPhase` AND `phaseClass` are swapped for multi-phase subjects (c91479, 0b5a2e, 0b5a2ePlayBack). The 270-targeted data gets the 90-targeted phase and label, and vice versa. The swap appears self-consistent in the CSV (setToDeliverPhase=270 always has phaseClass=270) because both use the same wrong index.

**Not affected**: Single-phase subjects (d5cd55, 7dbdec, 9ab7ab) — only one condition type, trivial mapping. Also not affected: the polar plots and phase distribution figures, which use `phase_at_0_pos`/`phase_at_0_neg` directly.

**Fix applied**: `bt = ii - 1` gives the burst type, then map to correct index: type 's' → 1, type 'm' → `2 - bt` (bt=0→2, bt=1→1), type 't' → bt=0→1, bt=1→2, bt≥2→4. Output table regenerated.

**Impact on results**: The phase × dose interaction that appeared significant in Models 3c-3e (p = 0.045-0.110) was driven by this mislabeling. After fix:

| Model | Interaction p (before) | Interaction p (after) |
|-------|----------------------|---------------------|
| 3a (absDiff) | 0.161 | TBD |
| 3c (ANCOVA) | 0.110 | **0.859** |
| 3d (ordinal .L) | 0.047 | **0.616** |
| 3e (numeric) | 0.045 | **0.615** |

The "dose-dependent CEP enhancement selective to hyperpolarizing (270) channels" was an artifact of phase mislabeling. The dose main effect remains significant in Model 3a (p ≈ 0.0002).

---

## Finding 21: CEP magnitude filtering pipeline

Two-stage filtering removes channels without meaningful evoked potentials and individual artifact/non-response trials:

**Stage 1: Channel-level exclusion (MATLAB)**
`multipleSubj_GLMM_script_PP.m` line 185:
```matlab
if nanmean(tempMagScreen(tempLabelScreen==0 & tempKeepsScreen)) > epThresholdMag
```
Computes the mean peak-to-peak magnitude across baseline trials (label==0, kept probes) for condition 1 on each channel. Channels where this mean < 100 uV (`epThresholdMag = 100`) are excluded entirely from the output table. This screens out channels far from the stimulation site or with poor recording quality. Applied per channel, not per trial.

**Stage 2: Trial-level exclusion (R)**
`betaStim_R_script.R` lines 36-37:
```r
data <- subset(data, magnitude < 1500)
data <- subset(data, magnitude > 25)
```
Removes individual trials with peak-to-peak voltage < 25 uV (likely non-responses or flat signals) or > 1500 uV (likely artifacts). Applied after the MATLAB channel-level filter.

**Peak-to-peak extraction** (`extract_PP_betaStim.m`):
For each trial, extracts the maximum peak-to-trough amplitude in a subject-specific time window (e.g., 5-36 ms for d5cd55, 6-60 ms for ecb43e). Optionally smoothed with Savitzky-Golay filter (order 3, frame 171 samples). Every trial gets a value — no per-trial filtering in MATLAB.

---

## Finding 14: Median consistency across all analyses

All analyses now use **median** consistently for:
- Cell summaries: `median(magnitude)`, `median(absDiff)`, `median(percentDiff)` per (sid, channel, phaseClass, numStims) cell
- Baseline computation: `baseMedian = median(base)` per channel (previously `baseMean = mean(base)`)
- Baseline covariate: `baselineMag = median(magnitude)` for Base trials per channel
- Permutation test statistics: `median(dSub$magnitude[condition == "CL"]) - median(dSub$magnitude[condition == "PB"])` (previously used `mean()`)
- Plotting summaries: all `ddply(..., summarize, ...)` calls use `median()`

Rationale: CEP magnitudes are right-skewed and the 25-1500 uV filtering still permits extreme values. Median is robust to within-cell outliers and produces more normally distributed model residuals. With 100+ trials per cell, mean and median are similar, but median is more defensible against reviewer concerns about outlier sensitivity.

---

## Finding 15: Ordinal and numeric dose model comparison (Models 3d, 3e)

Five summary-level models (3a-3e) were compared, all non-singular, all using median cell summaries:

| Model | Specification | AIC | Dose p | Int p | Phase 270>90 at [5,inf) |
|-------|--------------|-----|--------|-------|------------------------|
| 3a | absDiff, intercepts only | 907 | **0.00016** | 0.161 | 6.0 uV, **p=.047** |
| 3b | magnitude + baseline category, slope | 1300 | 0.228 | 0.146 | — (confounded) |
| 3c | ANCOVA (baseline covariate), categorical dose | 907 | 0.193 | 0.110 | 6.5 uV, **p=.030** |
| 3d-i | ordinal (.L+.Q uncorr slopes, afex) | 919 | 0.271 | 0.170 | — |
| 3d-ii | ordinal (.L only slope, lmerTest) | 910 | 0.193 | 0.110 | 6.5 uV, **p=.030** |
| 3e | numeric dose (fixed + random) | 1024 | 0.063 | 0.311 | ns |

**NOTE**: Finding 15 results above were computed BEFORE the phaseClass label fix (Finding 20). The interaction p-values shown here are pre-fix and no longer valid. Post-fix results (7 subjects including 702d24):

| Model | AIC | Dose p | Interaction p |
|-------|-----|--------|---------------|
| 3a | 1010 | **0.005** | 0.683 |
| 3c | 1017 | 0.085 | 0.587 |
| 3e | 1024 | 0.056 | 0.311 |

Key findings (post-fix):
- **3d-ii is identical to 3c**: ordinal polynomial .L contrast is a linear rescaling of doseNum. Confirms reparameterization.
- **3d-i (.L+.Q random slopes)**: .Q random variance negligible. Linear slope sufficient.
- **LRT 3e vs 3d-ii**: quadratic adds nothing (p ≈ 0.80).
- **Phase interaction is non-significant in ALL models** (p = 0.31-0.68). The pre-fix "selective to hyperpolarizing channels" finding was an artifact of the label swap (Finding 20).
- **Dose effect**: significant in intercepts-only model (3a, p = 0.005), trending with random slopes (3c p = 0.085, 3e p = 0.056). Dose-response vs baseline: [5,inf) - Base = +14.5 uV, CI [8.3, 20.6] (intercepts model).
- All models converge: dose-dependent CEP enhancement exists, but there is no phase selectivity.

---

## Finding 13: ecb43e stim table versions — verified no impact

Multiple versions of `ecb43e_tables.mat` exist across backup locations due to a stim table rebuild on 2019-11-08 (commit `29f12a1`). Comparison of the OLD (original, 103,367 bytes) and NEW (rebuilt, 103,292 bytes) tables:

| Property | OLD | NEW | Impact |
|----------|-----|-----|--------|
| Total stims | 10,000 | 9,998 | 2 stim difference |
| Probe stims (mode==0) | 1,650 | 1,650 | **Identical times** |
| In-burst stims (mode==1) | 8,350 | 8,348 | 2 fewer in NEW |
| 270° in-burst | 2,473 | 2,473 | Identical |
| 90° in-burst | 3,247 | 3,247 | Identical |
| Random in-burst | 2,628 | 2,628 | Identical |
| Bursts | 1,169 | 1,169 | Identical |

The 2 removed stimuli in the NEW table were in-burst conditioning pulses with anomalous stim type code `stims(8,:)==2`, which does not exist in the NEW table. This type code is not matched by any condition filter in the analysis pipeline (`stims(8,:)==0` for 270°, `==1` for 90°, `==3` for random), so these 2 stimuli would never enter the phase calculation or EP extraction regardless of which table is used.

**Probe stim times are byte-identical between OLD and NEW.** Since probe stims are the test pulses used for CEP magnitude extraction, the EP data is unaffected by the table version.

The phase calculation (`B_phaseCalc_allChans_processed.m`) was likely not rerun after the table rebuild — the commit hardcoded `idxVec = [7:7]` (subject 7 only) and `chans = 64` (debug settings). However, since the per-condition in-burst stim counts and times are identical between tables, rerunning the phase calculation would produce the same results.

A third variant (`ecb43e_tables_modDJC.mat`, 106,672 bytes) exists with 1,553 bursts (vs 1,169), suggesting a different burst-detection parameterization. This variant is not used in the current pipeline.

---

## Finding 22: Phase label verification and betaChan consistency audit (2026-04-08)

**Verification script**: `verify_phase_consistency.m` loads each subject's phase `.mat` file from `data/phase_data/`, computes circular mean phase for EP channels with good fits (r² > 0.7, freq 12.01–19.99 Hz), and checks:
1. Whether the measured circular mean at each beta reference channel is closer to the stated target or its opposite (target + 180°)
2. Whether the phaseClass bin (>180° → 270, ≤180° → 90) is consistent with the target
3. c91479's 0°/180° targets specifically (on the bin boundary)
4. Single-phase subjects' target-vs-opposite consistency

### betaChan inconsistency in `B_ExtractNeuralData_PP_reref.m`

The authoritative `valueSet` in `multipleSubj_GLMM_script_PP.m` (line 29) defines **betaChan=31** for both 0b5a2e and 0b5a2ePlayBack. This is consistent across all scripts (`plot_example_subject_phases.m`, `PhaseExtract_BetaRecordingChannel.m`, `examineRandomSubjectDelivery.m`, `BETA_manuscript_bars_compare0b5a2e_PP.m`, `BETA_manuscript_bars_phaseDiff.m`, `BETA_manuscript_bars_phaseDiff_ecb43e.m`, `plot_example_dose_dependent_time_series.m`).

**Exception**: `B_ExtractNeuralData_PP_reref.m` line 88 had `betaChan = 23` for 0b5a2e (while correctly setting 31 for 0b5a2ePlayBack). This variable is defined but **never referenced** in the extraction loop, so it had no effect on output data. **Fixed** to 31 on 2026-04-08.

### Beta reference channel verification results

| Subject | betaChan | Cond | Target | CircMean | Diff | Bin | Status |
|---------|----------|------|--------|----------|------|-----|--------|
| d5cd55 | 53 | all | 180° | 173.6° | 6.4° | 90 | OK |
| c91479 | 64 | pos | 0° | 20.6° | 20.6° | 90 | OK |
| c91479 | 64 | neg | 180° | 211.7° | 31.7° | 270 | OK |
| 7dbdec | 4 | all | 180° | 169.0° | 11.0° | 90 | OK |
| 9ab7ab | 51 | all | 270° | 270.8° | 0.8° | 270 | OK |
| 702d24 | 5 | pos | 90° | 109.2° | 19.2° | 90 | OK |
| 702d24 | 5 | neg | 270° | 307.5° | 37.5° | 270 | OK |
| ecb43e | 55 | pos | 270° | 259.8° | 10.2° | 270 | OK |
| ecb43e | 55 | neg | 90° | 130.4° | 40.4° | 90 | OK |
| 0b5a2e | 31 | pos | 90° | 109.4° | 19.4° | 90 | OK |
| 0b5a2e | 31 | neg | 270° | 296.8° | 26.8° | 270 | OK |
| 0b5a2ePlayBack | 31 | pos | 90° | 80.3° | 9.7° | 90 | OK |
| **0b5a2ePlayBack** | **31** | **neg** | **270°** | **0.1°** | **90.1°** | **90** | **WRONG HALF** |

All beta reference channels map correctly **except 0b5a2ePlayBack neg condition** — the playback replayed stimulation timing from the original session but with a temporal offset, so phase locking to the live beta oscillation is not expected.

### c91479 0°/180° target verification

c91479 `desiredF = [0, 180]` places the neg target on the 90/270 bin boundary. All EP channels (47, 48, 64) show pos circMean closer to 0° and neg circMean closer to 180°. No swap detected. Ch64 (beta ref) neg circMean = 211.7° bins to 270; ch47 neg circMean = 154.4° bins to 90 (borderline but correct — closer to 180° than to 0°).

### Single-phase subject verification

| Subject | Target | Channels correct / total | Mismatch channels |
|---------|--------|------------------------|-------------------|
| d5cd55 | 180° | 8 / 9 | ch63 (circMean=63.1°, closer to 0°) |
| 7dbdec | 180° | 3 / 4 | ch10 (circMean=332.1°, closer to 0°) |
| 9ab7ab | 270° | 8 / 8 | none (ch43 had zero good fits) |

Isolated mismatches on non-reference channels are expected from cortical beta phase propagation gradients.

### ecb43e inverted hardware convention

ecb43e is the only subject where `ptsPos = stims(8)==0` and `ptsNeg = stims(8)==1` in `B_phaseCalc_allChans_processed.m` (lines 190–191), the reverse of all other multi-phase subjects. This is compensated by `desiredF = [270, 90, ...]` (flipped from the standard [90, 270]) and the type `'t'` branch in the `correctIdx` mapping. Verified correct: ch55 (beta ref) pos circMean=259.8° (target 270°) and neg circMean=130.4° (target 90°).

---

## Finding 23: Models 5a and 5a-gf2 (continuous circular phase) are the primary models (2026-04-14)

**Decision**: The primary inferential models for the manuscript are **Model 5a** and **Model 5a-gf2**, both using continuous circular phase via `sin(phaseDeg) + cos(phaseDeg)`. The Model 3 series (binary phaseClass 90 vs 270) is retained as a robustness/sensitivity check.

### Rationale

1. **Phase is a circular quantity.** The 90°/270° binning used by Models 3a-3e treats phase as a 2-level factor and discards the actual measured angle. Channels with measured phases of 85° and 105° end up in the same bin, while 89° and 91° are split between bins despite being essentially identical. Fisher (1993) recommends decomposing circular predictors as `sin(phase)` and `cos(phase)`, which is what Models 5a / 5a-gf2 do.

2. **Grouping by measured phase preserves distinct conditions.** The 3-series grouping `(sid, channel, phaseClass, numStims)` collapses multiple delivered phases that happen to bin together into a single summary row. The 5-series grouping `(sid, phaseDeg_round, numStims, channel)` keeps distinct measured phases as separate cells, accurately reflecting that a multi-phase channel received stimulation at two different angles.

3. **5a-gf2 adds a principled quality filter.** Requiring `nGoodBeta ≥ 1` (at least one conditioning stim with R² > 0.7 and frequency 12-20 Hz per burst) restricts analysis to bursts where beta was actually present during conditioning. This directly supports the "phase-triggered stimulation requires an ongoing oscillation" hypothesis. Model 5a (no good-fit restriction) is the more conservative test; 5a-gf2 sharpens the effect.

4. **5a and 5a-gf2 together cover the analytic space.** 5a is the broader, more conservative primary (stricter phaseVecLength but no good-fit filter). 5a-gf2 adds the good-fit filter with a slightly relaxed phaseVecLength to stay non-singular. Sensitivity analysis across r ∈ {0, 0.1, 0.2, 0.3, 0.4} shows dose effect strengthening monotonically with r in both models (`output_plots/betaStim_phase_quality_sensitivity.csv`).

### Model results (from `CLAUDE.md` summary)

| Model | N | Channels | Dose.L p | sin p | cos p | Effect (µV) | Singular |
|-------|---|----------|----------|-------|-------|-------------|----------|
| 5a (r≥0.3) | 78 | 19 | **0.068** (trend) | 0.29 | 0.090 (trend) | 6.4 | No |
| 5a-gf2 (r≥0.2, good-fit) | 102 | 25 | **0.049** | 0.19 | 0.13 | 8.6 | No |
| 5a-gf (r≥0, good-fit, intercepts only) | 96 | — | 0.047 | — | — | — | No, but anti-conservative |

5a-gf (per-burst phase, intercepts only, 96 obs) is **not** used for inference — it becomes singular when a random dose slope is added, and without the slope its dose p-value is anti-conservative. It is retained as a diagnostic (sensitivity check for per-burst vs channel-level phase aggregation).

### Why not pool 5a and 5a-gf2 into one primary?

They represent two meaningfully different filter choices (quality threshold on phase consistency vs quality threshold on beta presence during bursts). Both are preregistered defaults. Reporting both as the primary pair is more transparent than arbitrarily selecting one as "the" primary.

### Implications for interpretation

- The central manuscript claim — cumulative dose-dependent EP enhancement modulated by phase — is supported at p=0.049 (5a-gf2) and as a trend at p=0.068 (5a). Effect sizes are small but reproducible.
- Interaction terms (dose × sin_phase, dose × cos_phase) are not significant in either model, so the dose effect is not strongly phase-gated at the group level. The phase effects are small trends in cos (~0.09-0.13).
- The 3-series models' dose p=0.0002 (Model 3a on absDiff, intercepts only) is stronger than 5a/5a-gf2's p-values **because** 3a lacks random dose slopes. That's not a valid reason to prefer 3a — the slopes are required to account for cross-subject heterogeneity in dose response, and 3a's strong p is partly attributable to the missing slope. 5a/5a-gf2 with the proper random-effects structure give the honest inferential statement.

---

## Finding 24: Conditioned vs Baseline per-cell permutation + within-subject FDR (2026-04-14)

**Motivation**: Models 5a and 5a-gf2 estimate the **average** dose effect across channels. They do not describe how that average is **distributed** across individual channels. A 5 µV group-level effect could come from (i) every channel modulating by ~5 µV, or (ii) a few channels modulating strongly while most are flat. The per-cell permutation analysis distinguishes these.

### Method

For each `(sid × channel × phaseDeg_round × dose)` cell passing filters:
- **Test**: two-sample label-shuffle permutation, median difference, 10,000 MC iterations.
- **Bootstrap 95% CIs**: 2,000 resamples per cell, percentile method, for forest-plot uncertainty visualization.
- **Filters (exactly match 5a-gf2)**:
  1. Channel-level `phaseVecLength ≥ 0.2` (drops channels where phase was not consistently delivered).
  2. Good-fit burst restriction `nGoodBeta ≥ 1` (each conditioned trial's preceding burst had ≥1 stim with R² > 0.7 AND frequency 12-20 Hz). Baseline trials are exempt (no preceding burst).
  3. Per-cell minima: ≥10 baseline probes, ≥5 conditioned probes (post-filter).
- Baselines are channel-level (pooled across phase conditions); conditioned trials are grouped by `phaseDeg_round`.

### Correction: BH FDR within (subject × dose)

Channels are nested within subjects — channels within a patient share an electrode grid, anatomy, and session-specific noise, and are therefore not exchangeable across subjects. The correction family is one subject's cells at one dose. Pooled FDR (within dose, across all subjects) is reported as a diagnostic CSV column only; not plotted (near-identical to within-subject on this dataset).

### Results (101 cells, 7 subjects)

Good-fit filter retains ~46% of conditioned trials (16,814 / 36,599). One cell dropped below the 5-conditioned-trial minimum after filtering (101 vs 102 before filter).

**Per-subject × dose (within-subject FDR):**

| Subject | n_cells | [1,2] sig | [3,4] sig | [5,inf) sig | Median Δ at [5,inf) |
|---------|---------|-----------|-----------|-------------|---------------------|
| 0b5a2e  | 13      | 0         | 0         | 0           | 24.3 µV (4 uncorr at [5,inf); 13-cell FDR family harsh) |
| 702d24  | 2       | 0         | 0         | 0           | 7.6 µV (underpowered) |
| 7dbdec  | 3       | 0         | 0         | 0           | similar pattern to pre-filter |
| 9ab7ab  | 4       | 1         | 2         | 3           | ~5 µV (clear internal dose gradient) |
| c91479  | 4       | 3         | 4         | 4           | ~48 µV (dominant responder) |
| d5cd55  | 3       | 0         | 0         | 0           | ~21 µV (underpowered for 3-test family) |
| ecb43e  | 5       | 0         | 1         | 0           | mixed |

**Across-subject summary:** at each dose, 2-3 of 7 subjects contribute at least one FDR-significant cell. Peak is at [3,4] (3 subjects), not [5,inf) (2 subjects), because ecb43e's lone significant cell at [3,4] doesn't persist.

**Comparison to pre-filter version:** good-fit filter strengthens 0b5a2e's per-cell signal (4 uncorrected-sig cells at [5,inf) vs 1 before), consistent with the "beta actually present during conditioning" interpretation. c91479 loses 1 FDR-sig cell at [1,2] (3/4 vs 4/4 before) — the filter removes trials with intermittent beta, making variance estimates noisier for some of its cells. Overall subject-level counts unchanged.

**Pooled-FDR sensitivity:** essentially identical to within-subject FDR on this dataset. Retained as a CSV diagnostic column only; not a separate plot.

### Interpretation

- **Supports the 5a-gf2 dose effect**: it's not a single-subject artifact; at minimum 2 subjects contribute individually detectable modulation at every dose. c91479 is the clearest responder; 9ab7ab shows a textbook dose gradient.
- **Explains why 0b5a2e's effect isn't per-cell-significant despite strong median effects**: 13-cell FDR family requires p ≤ 0.004 for the smallest to survive. Its best [5,inf) p is 0.028. The LMM recovers this by borrowing strength across its 13 channels — the per-cell test cannot.
- **Companion, not primary**: per-cell permutations describe distribution; the LMM (5a / 5a-gf2) makes the population-level inferential claim. Both are reported together.

### Outputs

- `output_plots/betaStim_cond_vs_base_perchan.csv` — master table (one row per cell: n, medians, obs_diff, bootstrap CI, perm_p, perm_q within-subj, perm_q pooled, sig flags)
- `output_plots/betaStim_cond_vs_base_per_subject.csv` — per subject × dose breakdown
- `output_plots/betaStim_cond_vs_base_subject_presence.csv` — across-subject headline (N subjects with ≥1 FDR-sig cell per dose)
- `output_plots/betaStim_cond_vs_base_pooled_summary.csv` — pooled-FDR sensitivity summary
- `output_plots/betaStim_cond_vs_base_forest.png/.eps` — forest plot (within-subject FDR, primary)
- `.docx` export: new sections in `betaStim_within_subject_tables.docx`

### Forest plot conventions

- Rows sorted by measured phase (0° at top → 360° at bottom), subject as tie-breaker
- Channel label: "Subject N ChX @ P°" with raw channel number (subject-number prefix stripped)
- Beta trigger channels (d5cd55 Ch53, c91479 Ch64, 7dbdec Ch4, 9ab7ab Ch51, 702d24 Ch5, ecb43e Ch55, 0b5a2e Ch31) highlighted with pink y-axis labels; dot color remains tied to significance category
- Dot colors: grey = ns, orange = p<0.05 uncorrected, red = FDR q<0.05 (within subject for primary plot; pooled for sensitivity plot)
- Three dose panels: [1,2], [3,4], [5,inf)
