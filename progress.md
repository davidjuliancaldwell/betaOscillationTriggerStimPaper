# Progress

## Completed

### Statistical audit of mixed effects model
- Audited `R_analysis_scripts/betaStim_R_script.R` and upstream MATLAB table generation (`peak_extraction/multipleSubj_GLMM_script_PP.m`)
- Identified 12 findings documented in `statistical_audit.md`
- Key finding: `phaseClass` is a channel-level constant (not trial-level), effective N = 49 units not 37K trials
- Only 9 of 31 channels have within-channel phaseClass variation
- 9ab7ab contributes no phaseClass contrast (all channels = 270)

### Closed-loop vs playback control analysis (`R_compare_control_cond.R`)
- Added Cohen's d effect size calculations using `emmeans::eff_size()` with 95% CIs
- Added permutation test (10,000 permutations, two-sided) for channel 14:
  - Per-dose CL vs PB: significant at [1,2], [3,4], [5,inf) (p < 0.001); Base not significant (p = 0.12)
  - Dose-response interaction (Base vs [5,inf) gap): 68.9 uV, p = 0.039
- Added plots to `output_plots/`: null distribution histogram, per-dose lollipop, Cohen's d bar plot
- Ran all-channel permutation (8 matched pairs): interaction not significant across all channels (p = 0.18)

### Mixed model improvements (`betaStim_R_script.R`)
Seven models with progressively improved random effects and sensitivity analyses:

| Model | Description | Random effects |
|-------|-------------|---------------|
| `fit.intercepts.only` | Original baseline | `(1\|sid/channel)` |
| `fit.trial.level` | + dose slopes per subject | `(0+numStims\|sid) + (1\|channel)` |
| `fit.summary.level` | Summary-level, no pseudoreplication | Same, on medians per cell |
| `fit.nested.condition` | + condition nesting (Model 4) | `+ (1\|channel:setToDeliverPhase)` |
| `fit.no.9ab7ab` | Model 4, exclude 9ab7ab (no phase contrast) | Same |
| `fit.clean` | Model 4, exclude 9ab7ab + ecb43e random trials | Same |
| `fit.within.channel` | Model 4, within-channel phaseClass only (9 channels, 3 subjects) | Same |

DF comparison:

| Effect | Model 1 | Model 4 | Model 4c |
|--------|---------|---------|----------|
| numStims | df=37K, p=2e-6 | df=5, p=0.16 | df=2, p=0.25 |
| phaseClass | df=4.6K, p=0.14 | df=31, p=0.27 | df=11, p=0.77 |
| interaction | df=37K, p=6e-4 | df=753, p=0.040 | df=2778, p=0.70 |

Key results:
- Interaction significant in Models 4/4a/4b (p = 0.04-0.05), driven by [5,inf) dose: phaseClass 270 > 90 by ~6 uV
- Interaction NOT significant in Model 4c (within-channel only, p = 0.70) — effect drops to 2.4 uV
- Could be power (3 subjects, 9 channels) or between-channel confound driving the full-dataset result

### Documentation
- `CLAUDE.md` — Pipeline architecture, all model specifications, data structure notes
- `statistical_audit.md` — 12 findings and 8 recommendations
- `progress.md` — This file
- `README.md` — R script descriptions

---

### Subject 6 (ecb43e) random vs phase-targeted analysis (`R_compare_subject_6_random.R`)
- Fixed: Windows path, CSV file (now uses 100_thresh), y-axis label, `summaryDataChan` column name bug, savePlot default
- Added: linear model (`lm`), Cohen's d via `emmeans::eff_size()`, permutation tests (targeted vs random at each dose)
- Result: no significant difference between phase-targeted and random stimulation on channel 655 (all perm p > 0.05, all |d| < 0.3)

### ANOVA type documentation
- Type III ANOVA is correct for these models (interaction is significant, main effects are conditional)
- Type II and Type III agree on the interaction (p = 0.040); differ on main effects (phaseClass: Type III p=0.27, Type II p=0.51)
- Documented the four complementary outputs (Type III ANOVA, coefficients/tab_model, emmeans, effect sizes)

---

## Not yet started

- [ ] Add permutation test validation for phaseClass p-values in main model
- [ ] Consider expanding `R_compare_control_cond.R` to all 8 matched channels with `lmer(magnitude ~ numStims * sid + (1|channel))`
