# Statistical Audit: Mixed Effects Model Analysis

Audit of `R_analysis_scripts/betaStim_R_script.R` and its upstream data generation in `peak_extraction/multipleSubj_GLMM_script_PP.m`.

## Model Under Review

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

1. **Separate within-channel and between-channel phaseClass analyses.** Restrict the primary phaseClass analysis to the 9 channels (across 3 subjects) with within-channel variation. Report the between-channel analysis separately as supporting evidence, clearly noting the confound.

2. **Remove 9ab7ab from any model that includes phaseClass** as a predictor. It contributes no contrast and distorts the dataset balance.

3. **Exclude or separately analyze ecb43e random-condition trials** (`setToDeliverPhase == 12345`). They weren't phase-targeted and shouldn't receive the same treatment as deliberate phase conditions.

4. **Use the summary-level model (Model 3) as the primary analysis** to eliminate pseudoreplication. Collapse to one median per (sid, phaseClass, numStims, channel) cell (120 observations). Use `magnitude ~ numStims * phaseClass + (1|sid) + (1|channel)` — random intercepts only. Random dose slopes cause singularity with 6 subjects; `setToDeliverPhase` is a fixed experimental condition and should not be a random grouping variable.

5. **Validate phaseClass p-values** by comparing trial-level model standard errors against a channel-level analysis (e.g., permutation test or bootstrap at the channel level).

6. **Consider the continuous phase variable** (`phaseDeliveryBinned45`, 8 bins of 45 degrees) instead of the binary 90/270 classification as a sensitivity analysis, at least for within-channel contrasts where more phase resolution is available.

7. **Report effective sample sizes** alongside the model. The reader should know that the phaseClass effect is estimated from ~49 channel-condition units (or 9 within-channel contrasts), not 37K trials.

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
