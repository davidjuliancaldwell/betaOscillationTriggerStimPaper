# Progress

## Completed

### Phase-label backport, Model 5 sin/cos, d5cd55 fix, bundled deps (2026-04-12)

#### Phase-label swap fix backported to `phase_vs_peak.m`
- `phase_vs_peak.m` had the same burst-type ↔ phase-variable swap as the 2026-04-06 fix in `multipleSubj_GLMM_script_PP.m`. For type 'm' subjects (c91479, 702d24, 0b5a2e, 0b5a2ePlayBack), `index=1` was paired with `phase_at_0_pos` (burst type 1 phases) but `dataForPPanalysis{chan}{1}` (burst type 0 EPs) — inverted.
- Fix: replaced combined if-elseif with explicit per-type-per-index mapping via intermediate variables (`rsq_use`, `phase_use`, `f_use`, `target`). Added missing `type=='t' && index==4` branch for ecb43e random condition.
- Verified empirically: c91479 ch64 dot positions swapped correctly after fix.
- Output file naming now encodes active filters (e.g., `_r30`, `_r20_gf1`). Config block uses `exist(..., 'var')` checks so a wrapper script can pre-set filter values.
- New output files: `phase_vs_peak_all_subj_median_r30.{png,eps}`, `phase_vs_peak_all_subj_median_r20_gf1.{png,eps}`, per-subject variants.

#### Full pipeline phase-label consistency audit
- Systematic trace of `stims(8)` → burst type → phase variable → CSV column → R analysis for all subjects and pipeline stages. All stages confirmed consistent. No remaining phase-label swap bugs.

#### Model 5 (sin/cos continuous phase) added to `betaStim_R_script.R`
- Three variants: 5a (channel-level phase, phaseVecLength ≥ 0.3, 78 obs, random dose slope), 5a-gf (good-fit trials, per-burst phase, 96 obs, intercepts only — random slope singular), 5a-gf2 (good-fit trials, channel-level phase, phaseVecLength ≥ 0.2, 102 obs, random dose slope).
- Results: 5a Dose.L p=0.068 (trend), 5a-gf2 Dose.L p=0.049 (nominally significant). Phase main effects weakly estimated (cos_phase p~0.09-0.13). 5a-gf dose p=0.047 is anti-conservative (intercepts only).
- Output CSV includes `phaseDeg`, `phaseVecLength`, `phaseCircStd`, `phaseOmnibusP` from MATLAB.
- Sensitivity analysis at r ∈ {0, 0.1, 0.2, 0.3, 0.4} → `output_plots/betaStim_phase_quality_sensitivity.csv`.
- Effect sizes use corrected `d_total` formula: `Var_int + E[x²]*Var_slope + Var_channel + Var_resid` where E[x²]=1/3 for `contr.poly(3)`.

#### Critical d5cd55 probeSample alignment bug fixed
- `multipleSubj_GLMM_script_PP.m` was reconstructing `probeStims` without applying d5cd55's time filter (`stims(2,:) > 36536266`), giving 1982 probes instead of 1563. Magnitudes/dose labels were correct; only `probeSample` was wrong (0% match with precision CSV before fix, 100% after).
- Fix: reproduce extraction's pts selector in un-shifted coordinates (`stims(2,:) > 36536266 - delayDelivery`). Mirrored in `phase_vs_peak.m` for the good-fit filter path.
- Impact: d5cd55 Model 5a-gf/5a-gf2 results changed (bad precision CSV merge); Models 3a/3c/3e/5a unaffected.

#### Critical 702d24 extraction bug fixed
- Three interrelated bugs caused 702d24 ch5 to produce nearly all-NaN extractions (11% valid → 72% valid after fix):
  1. `framelen=171` too wide for 702d24's 21 ms window → `ppFramelen=91` via new 8th arg in `extract_PP_betaStim.m`
  2. `t_min=3.8ms` not accounting for `delayDelivery=14` shift → `t_min=0.00323` for 702d24
  3. `plot_EP_goodfit_by_phase.m` (new file) now applies `delayDelivery` shift to match pipeline
- Regenerated `betaStim_outputTable_50_new_100_thresh.csv`.

#### External dependencies bundled
- Copied `CircStat2012a/` and `sgolayfilt_complete.m`/`savitzkyGolay.m` into `external_deps/`. `setup_environment.m` now uses `addpath(genpath(locationsDir))` only — no external paths needed.

#### CL vs PB analysis expanded (`R_compare_control_cond.R`)
- Burst-quality filter config block: `minGoodBetaPerBurst_clpb`, `minBurstVecLength_clpb`, `minPhaseVecLength_clpb = 0.2`. Filter removes 26% trials + 3 of 16 channel × condition cells.
- Exact sign-flip permutations (replacing Monte Carlo) for small-n aggregated tests: n=8 (Null vs Base) and n=13–16 (clpb channel × condition). Deterministic, reproducible.
- Null-burst vs baseline validity test: 0/8 channels significant; aggregated perm p=0.226. Null EPs indistinguishable from baseline.
- New visualizations: 8×3 grid `betaStim_clpb_grid_ch_x_dose.{png,eps}`, forest plot `betaStim_clpb_forest_paired_effects.{png,eps}`, 8-phase dose plots (45° bins).
- CSV outputs: `betaStim_clpb_perm_chan_aggregate.csv`, `betaStim_clpb_perm_perchan_bycell.csv`, `betaStim_null_vs_base_perchan.csv`, `betaStim_null_vs_base_aggregate.csv`, percent modulation CSVs.
- Added to `.docx`: 4 new sign-flip permutation tables.

#### Other bug fixes and cleanup
- `compute_burst_phase_precision.m`: fixed `load()` overwriting `sid` argument; fixed playback phase file path (was loading CL's file).
- `plotting_functions/SaveFig.m`: fixed Windows-only path logic to recognize Unix absolute paths.
- `plotting_functions/plot_phase_cortex.m`: fixed legend icon sizes via `findobj` after `drawnow`.
- `helper_functions/getSubjDir.m`: fixed hardcoded Windows backslashes → `fullfile()`.
- `peak_extraction/C_PlotBrains_PP.m` + `phase_vs_peak.m`: removed spurious ch63 from ecb43e goodEPs (was silently filtered anyway).
- `BETA_manuscript_bars_compare0b5a2e_PP.m`: trimmed dead-weight 6 unused SIDS/valueSet entries.
- `showTabModel = FALSE` now consistently gated in `R_compare_control_cond.R` and `R_compare_subject_6_random.R`.
- `R_burst_phase_analysis.R`: minor fix.
- CLAUDE.md: trimmed from 525 to ~300 lines (removed development log/bug-fix histories, kept current-state facts).

---

### Ordinal/numeric dose models + median consistency (2026-04-05)
- Switched ALL analyses to use median: cell summaries, baseline computation (baseMean→baseMedian), permutation test statistics, plotting summaries
- Added Model 3d (ordinal dose with polynomial contrasts):
  - 3d-i: `afex::mixed(expand_re=TRUE, per_parameter="numStims_ord")` with uncorrelated .L+.Q random slopes — .Q variance ~0, confirming linear slope sufficient
  - 3d-ii: `lmerTest::lmer` with manual `dose_linpoly` (.L only) — identical to Model 3c (reparameterization)
- Added Model 3e (fully numeric dose): `doseNum * phaseClass + baselineMag_c` — most parsimonious (interaction was p=0.045 before phase label fix, now p=0.31)
- LRT 3e vs 3d: p=0.76, quadratic component unnecessary
- All five models (3a-3e) non-singular, all converge on same pattern: dose-dependent enhancement selective to phase 270
- With median, all four ANCOVA-family models now show significant phase contrast at [5,inf): 3a p=.047, 3c p=.030, 3d-ii p=.030, 3e p=.038
- Files changed: betaStim_R_script.R, R_compare_control_cond.R, R_compare_subject_6_random.R, compare_three_models.R

### Model restructuring and effect size audit (2026-04-03)
- Audited all effect size calculations (partial eta², Cohen's d, Hedge's g) across R and MATLAB scripts
- Restructured summary-level models into three complementary specifications (3a, 3b, 3c)
- Model 3c (ANCOVA with baseline covariate + random linear dose slope) is now primary
- Key finding: baseline covariate (beta~1.03) absorbs channel variance far better than random intercepts
- Interaction coefficient significant in ANCOVA: `numStims[5,inf):phaseClass90 = -7.34, p = 0.047`
- Discovered phase 90 channels have higher baseline magnitudes (~349 vs 261 uV) — between-channel confound affects interpretation when baseline is modeled as a category vs covariate
- Added emmip contrast plots and forest plots for all three models
- Added Cohen's d for Models 3a and 3c using consistent trial-level sigma
- Commented out SVG saves (require svglite package not installed)
- Added `compare_three_models.R` standalone comparison script

### Final results with 702d24 included + burst phase analysis (2026-04-06)
- Included 702d24 (1 channel) in main analysis — no convergence issues. 7 subjects, 32 channels.
- Added `probeSample` column to output table for robust trial-level merge with burst precision data
- Created `R_burst_phase_analysis.R`: tests phase direction (270° bin vs others) and phase error across all 7 subjects
- **No phase effect**: 270° bin vs others p = 0.806, interaction p = 0.838 in mixed model. Per-subject phase error correlations disappear in mixed model (between-channel confounds).
- **Good-beta-fit trials show NO dose effect** (p = 0.85) and actually have SMALLER magnitudes than non-fit trials. Beta presence predicts smaller, not larger, CEPs.
- **Dose effect confirmed**: [5,inf) vs Base = +14.5 uV, CI [8.3, 20.6] (intercepts model, p = 1.5e-6). With random slopes: same effect size but underpowered (p = 0.085).
- Documented as Finding 21 (filtering pipeline) in statistical_audit.md

### **CRITICAL BUG FIX**: phaseClass AND setToDeliverPhase swap (2026-04-06)
- **Bug**: `multipleSubj_GLMM_script_PP.m` used loop index `ii` to index `desiredF(ii)` and `peakPhaseVec(ii,...)`, but `ii` iterates over sorted burst types (0,1,2) while `desiredF` and `peakPhaseVec` are indexed by the phase calc loop (index 1=pos/stims(8)==1, index 2=neg/stims(8)==0). For type 'm' subjects: burst type 0 (neg/270-targeted) got index 1 (pos/90-targeted phase and label), and vice versa.
- **Affected**: BOTH `setToDeliverPhase` AND `phaseClass` columns for multi-phase subjects (c91479, 0b5a2e, 0b5a2ePlayBack). Single-phase subjects unaffected.
- **Fix**: map `bt = ii - 1` (burst type) to correct index. For type 'm': `correctIdx = 2 - bt`. Applied and regenerated output table.
- **Impact on results**: The phase × dose interaction that was significant in Models 3c-3e (p = 0.045-0.110) is now **non-significant (p = 0.62-0.86)**. The "dose-dependent enhancement selective to hyperpolarizing channels" was an artifact of misassigned phase labels.
- Dose main effect remains (Model 3a p = 0.0002 sensitivity analysis)

### Burst phase precision analysis (2026-04-06)
- Created `compute_burst_phase_precision.m`: links each test stim to its burst's beta-band phase fits
- Ran for 0b5a2e: 9864 rows (9 channels x 1096 probes), exported to CSV
- Beta channel (31) at [5,inf): 85/102 trials have beta fits, but most have only 1-2 fits per burst (median R ≈ 1.0)
- Channel 14 at [5,inf): significant REVERSE effect — loose phase precision produces LARGER CEPs (diff = -107 uV, perm p = 0.006)
- Documented as Findings 19-20 in statistical_audit.md

### CL vs Playback deep dive + phase fit screening analysis (2026-04-06)
- Analyzed beta channel (31) for CL vs PB: CL > PB by ~30-40 uV at all doses (session offset), no dose-dependent interaction
- Channel 14 phase delivery: both conditions give ~270° (doesn't track targeting condition). Only beta channel (31) shows correct phase delivery (109° for 90-target, 297° for 270-target)
- Heavy tail analysis on channel 14: CL condition has extreme trials (up to 1285 uV, 0 PB > 540 uV). Mean-based interaction (p=0.039) driven by tail; median-based interaction p=0.458
- Log LM on channel 14: interaction p=0.095, skew drops 2.5→0.3
- All-channel linear model (8 pairs): CL > PB by 23 uV (p=3.9e-07), no interaction (p=0.968)
- Phase fit screening bottleneck: frequency filter (12-20 Hz) rejects 83% of trials; R² > 0.7 alone passes 70%. Among in-band fits, 96% have R² > 0.7.
- Created `compute_burst_phase_precision.m` for per-burst phase precision analysis
- Symlinked data directories to OneDrive backup
- Documented as Findings 16-18 in statistical_audit.md

### ecb43e stim table verification (2026-04-03)
- Compared OLD (original, pre-Nov 2019) and NEW (rebuilt) `ecb43e_tables.mat` across OneDrive backups
- OLD: 10,000 stims; NEW: 9,998 stims — difference is exactly 2 in-burst conditioning pulses with anomalous type code (stims(8,:)==2)
- Probe stim times (mode==0, used for CEP extraction) are byte-identical between versions: 1,650 probes
- Per-condition in-burst counts identical: 2,473 (270°), 3,247 (90°), 2,628 (random)
- Phase calculation was not rerun after table rebuild (commit hardcoded to subject 7 only), but would produce identical results since the relevant stim times are unchanged
- Documented as Finding 13 in `statistical_audit.md`
- Also traced full phase fitting pipeline for ecb43e: sinfit.m → phase_calculation.m → phase_circstats_calc.m → binary binning at 180° threshold

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

#### Primary summary-level models (all non-singular, one median per cell)

| Model | Outcome | Baseline handling | Random effects | AIC |
|-------|---------|-------------------|---------------|-----|
| 3a `fit.absDiff` | absDiff (pre-subtracted) | Subtracted before model | `(1\|sid) + (1\|channel)` | 916 |
| 3b `fit.modelD` | magnitude | Category in numStims | `(1+doseNum\|sid) + (1\|channel)` | 1300 |
| **3c `fit.ancova`** | magnitude | **Fixed covariate** | `(1+doseNum\|sid) + (1\|channel)` | **907** |

Key design decisions:
- Random dose slopes use numeric `doseNum` (0,1,2,3) — categorical slopes cause singularity with 6 subjects
- Fixed effects keep `numStims` categorical to capture non-linear dose patterns
- ANCOVA baseline covariate (beta~1.03) absorbs channel-level variance far more effectively than random intercepts alone (channel SD: 153→10 uV)

Model 3c (ANCOVA, primary) results:
- Interaction coefficient: `numStims[5,inf):phaseClass90 = -7.34, p = 0.047`
- Phase contrast at [5,inf): 270 > 90 by 6.5 uV, CI [0.6, 12.3], p = 0.030
- Dose main effect: p = 0.19 (underpowered with random slope, 6 subjects)
- Sensitivity (3a, no random slopes): dose p = 0.0002, eta² = 0.19

All three models converge: dose-dependent CEP enhancement, no phase selectivity (phase interaction was artifact of label swap, corrected 2026-04-06).

#### Trial-level reference models (kept for comparison)

| Model | Description | Random effects |
|-------|-------------|---------------|
| `fit.intercepts.only` | Original, DF inflated | `(1\|sid/channel)` |
| `fit.trial.level` | + dose slopes, singular | `(0+numStims\|sid) + (1\|channel)` |
| `fit.nested.condition` | + condition nesting, singular | `+ (1\|channel:setToDeliverPhase)` |
| `fit.no.9ab7ab` | Model 4, exclude 9ab7ab | Same |
| `fit.clean` | Model 4, exclude 9ab7ab + random | Same |
| `fit.within.channel` | Within-channel only (9 ch, 3 subj) | Same |

Within-channel analysis (4c): interaction disappears (p = 0.70, effect drops to 2.4 uV). Could be power or between-channel confound.

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
