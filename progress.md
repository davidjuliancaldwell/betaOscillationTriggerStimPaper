# Progress

## Completed

### marginaleffects AME complement + dose-effect-vs-phase plots (2026-04-22)

#### Motivation

`emmeans::emmeans(fit, at = list(sin_phase = 1, cos_phase = 0, ...))` evaluates the dose effect at two specific phases (90°, 270°). A reviewer could reasonably ask "but what about all the other phases the study sampled?". `marginaleffects::avg_comparisons()` answers that by averaging pairwise dose contrasts over the **observed joint distribution** of sin_phase, cos_phase, baselineMag_c, and betaLabels — the phase-weighted population-level dose effect (G-computation / Average Marginal Effect).

#### Additions to `betaStim_R_script.R`

1. **Library**: added `library('marginaleffects')` (version 0.29.0 installed).

2. **AME numeric blocks** (after Model 5a effect sizes at line ~1285, after Model 5a-gf2 effect sizes at ~1867):
   ```r
   ame_5a_raw <- avg_comparisons(fit.sincos.ordinal,
     variables = list(numStims_ord = "pairwise"))
   ame_5a_adj <- as.data.frame(hypotheses(ame_5a_raw, multcomp = "single-step"))
   ame_5a_adj$contrast <- as.data.frame(ame_5a_raw)$contrast
   ame_5a_dose <- ame_5a_adj[, c("term", "contrast", ...)]
   ```
   Same pattern for `ame_gf2_dose`. `hypotheses(multcomp = "single-step")` applies multcomp::glht max-t adjustment (FWER control via the joint multivariate-t distribution) — the closest analog to Tukey-Kramer available in marginaleffects 0.29. `hypotheses()` drops the `contrast` label column, so it's reattached from the raw object.

3. **Dose-effect-vs-phase plots** (after `p_5a`/`p_5a_gf2` ggsave blocks):
   ```r
   cmp_5a_curve_df <- do.call(rbind, lapply(phase_vals_curve, function(ph) {
     nd <- summaryNB_m5
     nd$sin_phase <- sin(ph * pi / 180)
     nd$cos_phase <- cos(ph * pi / 180)
     df <- as.data.frame(avg_comparisons(fit.sincos.ordinal,
       variables = list(numStims_ord = "pairwise"),
       newdata = nd, re.form = NA))
     df$phase_deg <- ph
     df
   }))
   ```
   At each of 12 phase angles (0–330° in 30° steps), override sin/cos on the observed data and compute AME. Verified numerically identical to the canonical `datagrid(grid_type = "counterfactual") + by=` marginaleffects idiom (|est diff| = 0, |SE diff| = 5e-17 on a toy model) — the loop is needed only because `phase_deg` isn't a model predictor and the sin/cos pairing would otherwise cross-product into 144 off-unit-circle combinations. New output files:
   - `output_plots/betaStim_model5a_phase_curve_r30_dose_effect_vs_phase.{png,eps}`
   - `output_plots/betaStim_model5a_gf2_phase_curve_r20_dose_effect_vs_phase.{png,eps}`

4. **docx AME flextables**: added `render_ame_contrasts()` helper; new sections "AME Dose Contrasts (marginalized over observed phase)" inserted into `betaStim_statistical_tables.docx` after the existing EMM phase-contrast tables for each primary model.

#### Emmeans plot marginalization fix

The primary single-panel phase-response curves (`p_5a`, `p_5a_gf2`) previously pinned `betaLabels = "0"` in the `at = ...` spec — so they showed the non-beta-channel phase-response, not a population-level curve. Now use `weights = "proportional"` so emmeans marginalizes over `betaLabels` using observed proportions:

```r
em <- emmeans(fit.sincos.ordinal, ~ numStims_ord,
  at = list(sin_phase = ..., cos_phase = ..., baselineMag_c = 0),
  weights = "proportional")   # was: betaLabels = "0" in at=...
```

Matches the AME/G-computation semantics used by the new marginaleffects plots and the AME numeric tables in the docx. Supplementary faceted plots (`p_5a_beta`, `p_5a_gf2_beta`) unchanged — they still show both betaLabels levels side-by-side. `baselineMag_c = 0` retained because the covariate is grand-mean centered and enters linearly with no interactions → pinning at 0 is mathematically identical to averaging.

Subtitles reworded from `"baselineMag_c = 0, betaLabels = 0"` → `"baseline at grand mean; betaLabels marginalized over observed proportions"` for all primary and faceted plots.

#### Results

AME (single-step adjusted):

| Contrast | 5a estimate | 5a p | 5a-gf2 estimate | 5a-gf2 p |
|----------|-------------|------|-----------------|----------|
| [3,4] − [1,2] | 5.18 µV | 0.221 | 2.28 µV | 0.801 |
| [5,inf) − [1,2] | **9.09 µV** | **0.052** | **11.01 µV** | **0.048** |
| [5,inf) − [3,4] | 3.91 µV | 0.422 | **8.73 µV** | **0.041** |

Converge on the at-90°/at-270° EMM estimates — the phase-marginalized AME reinforces the dose-effect story without relying on two specific phase anchors. 5a-gf2 now shows two adjacent-dose contrasts clearing the 5% threshold; 5a has one trending contrast.

#### Verification

Full `betaStim_R_script.R` re-run, exit 0, 517 parsed expressions, 4 new plot files + regenerated docx. No regressions in existing outputs.

#### Why not also swap existing emmeans → marginaleffects elsewhere

Critical analysis before implementation concluded `emmeans` stays as the primary tool: `eff_size()` for Cohen's d with conditional/total variance denominators has no marginaleffects analog; log-scale back-transformation to percent change (`fit.lm.log` in `R_compare_control_cond.R`) is idiomatic in emmeans; spline and sin/cos phase-response curves work fine in emmeans with `at = ...`. AME is the one genuinely different question marginaleffects answers better — scoped narrowly to the 5a/5a-gf2 primary models + one complementary plot each.

---

### Reref/ECO-loading bug class fix + EP pipeline rerun verification (2026-04-21)

#### Bug class: narrow `achan`-based reload gate

TDT ECoG data is split across four 16-channel structs (`ECO1`–`ECO4`, ~1.7 GB each) inside the per-subject `_ECoG.mat` file. To avoid loading all four, the pipeline's channel loops cache the most-recently-loaded struct in `dataStruct` and reload only when the loop crosses into a new group. The reload gate was implemented as a heuristic — "reload if `achan ∈ {1, 2}`" (or sometimes `{1, 2, 4, 6}`) on the theory that any new group's first channel would have `achan=1`. This assumption breaks in two ways:

1. **Stale initial state** — the `plot_EP_goodfit_by_phase.m` reref loop starts with `dataStruct` already loaded for the **target** `chanInt`'s group (not `ECO1`). For `chanInt=47` (c91479 ch47, `grp=2 → ECO3`), the first 13 `rerefChans` iterations (`rc=4..16`) hit `achan ∈ {4..16}`, the narrow gate misses, and the loop reads ECO3 columns instead of the intended ECO1 columns. Result: 23% of reref channels silently replaced, ~5% CEP drift in the median-CAR downstream.

2. **Second-loop stale state in the main pipeline** — `B_ExtractNeuralData_PP_reref.m` has two loops (reref at line 158, peak-extract at line 228). The first uses the broader `{1,2,4,6}` gate and catches every group transition for current subjects; the second uses just `{1,2}` and enters stale from the end of the first loop. For c91479 (whose `chans` list starts at chan=4 / `achan=4`), the first ~12 iterations read from the wrong struct. Empirically those channels are excluded by the 100 µV baseline threshold so the R CSV never saw corrupted numbers — but that was luck, not design.

3. **B_phaseCalc_allChans_processed.m also has the narrow gate** — latent for all 8 subjects because every group transition in their bads/rerefChans lists happens to land at `achan ∈ {1, 2}`. Would silently corrupt any future subject whose bads block channels 1 and 2 of any group.

Documented extensively in `phase_data_origin_and_bugs.md` before fixing.

#### Fix: `prev_grp` group-change detector

Replaced all heuristic gates with a direct group-change detector:

```matlab
prev_grp = -1;        % forces load on first iteration regardless of caller state
for chan = <list>
  grp = floor((chan-1)/16);
  if grp ~= prev_grp
    load(..., sprintf('ECO%d', grp+1));
    dataStruct = eval(sprintf('ECO%d', grp+1));
    prev_grp = grp;
  end
  eco = 4 * dataStruct.data(:, chan - grp*16)';
  % ...
end
```

Applied to 8 loops across 5 files:

- `peak_extraction/B_ExtractNeuralData_PP_reref.m` (both loops, lines 158 and 228)
- `phase_visualizations/B_phaseCalc_allChans_processed.m` (channel loop, line 90)
- `manuscript_generate_scripts/plot_EP_goodfit_by_phase.m` (reref loop, line 120)
- `manuscript_generate_scripts/BETA_ExtractNeuralDataCEPscreen.m` (both loops)
- `manuscript_generate_scripts/plot_example_dose_dependent_time_series.m` (reref loop)

#### Collateral fixes (debug-leftover cleanup)

- `B_ExtractNeuralData_PP_reref.m:19` — `for idx = 1:1` (debug, d5cd55-only since March 2020) → `for idx = 1:8`
- `B_phaseCalc_allChans_processed.m:59` — removed `chans = 64;` hardcode (debug from Nov 2019 commit that would have processed only channel 64 for every subject)
- `B_phaseCalc_allChans_processed.m:2` — `idxVec = [1:7]` → `idxVec = [1:8]` (include playback)
- `B_phaseCalc_allChans_processed.m:395` — save filename `..._12samps_...` typo → `..._51samps_...` (matches downstream modifier and existing Dec 2018 hyak files)

#### Pipeline rerun (c91479 only) and verification

- **Why c91479 only**: per the safe-subject matrix in `phase_data_origin_and_bugs.md`, c91479 is the only subject whose `chans` list starts at `achan ∉ {1, 2}` — the only subject the chans-loop narrow gate could have bitten. Other subjects' EP data provably unchanged, so no rerun needed.
- **Why not phase**: Dec 2018 hyak-cluster `.mat` files are canonical; regenerating (now safe with the debug hardcode removed) would take hours for numerically-identical output.
- **Stages run**: `B_ExtractNeuralData_PP_reref` (c91479 only, patched `for idx = 2:2` with `onCleanup` restore; took 6 s), then `multipleSubj_GLMM_script_PP` across all 7 subjects (7 min). `compute_burst_phase_precision` not rerun — its input (phase data) unchanged, so outputs identical.
- **CSV diff vs backup**: all 10 numeric columns bit-identical for every subject and channel except c91479 ch47 (15 rows) and ch64 (137 rows), all differing by ≤1e-12 µV (floating-point accumulation noise). Confirms the pre-fix CSV magnitudes for analyzed channels were already correct — the bug never touched the filtered-in data.
- **R primary models reproduce CLAUDE.md exactly**:
  - Model 5a (78 obs, 19 ch, 7 sid, non-singular): `numStims_ord.L` p = 0.0676, `cos_phase` p = 0.0899, Type III dose F p = 0.123, effect 6.36 µV
  - Model 5a-gf2 (102 obs, 25 ch, non-singular): `numStims_ord.L` p = 0.0531, `cos_phase` p = 0.0983, effect 8.10 µV
- `R_compare_control_cond.R` also completed cleanly.

#### plot_EP_goodfit_by_phase.m inclusion-indicator sgtitle

Added a pre-loop channel-level check that matches the R-CSV inclusion criterion (`mean Base PP ≥ 100 µV`, per `multipleSubj_GLMM_script_PP.m:185`). Sgtitle now shows:

- **INCLUDED** (black): `[INCLUDED: mean Base PP X.X ≥ 100 µV]`
- **EXCLUDED** (red): `[EXCLUDED from R analysis: mean Base PP X.X < 100 µV threshold]`

Ch48 for c91479 correctly flagged EXCLUDED — pre-fix it had Base PP = 0 µV on the plot (artifact of the old reref corruption trashing single-trial peak extraction within the narrow 5–36 ms window); post-fix it shows Base PP = 36 µV (real but weak signal, still below threshold and correctly dropped from R analysis).

#### Plot regeneration with fixed reref

Regenerated 11 plots with corrected reref:
- c91479: ch47 (INCLUDED, mean 139.4 µV), ch48 (EXCLUDED, mean 46.7 µV), ch64 (INCLUDED, mean 251.7 µV)
- 0b5a2e: ch14 (340), ch15 (241), ch16 (277), ch21 (416), ch23 (433), ch31 (518, beta ref), ch32 (205), ch40 (146) — all INCLUDED

c91479 ch47/ch64 per-cell PP annotations shifted <5% as predicted (median CAR is robust even at 23% wrong reref channels; mean CAR would have shifted much more).

#### Backup

`data/_backup_20260421_161313/{output_table,phase_data,EP_data}/` holds full local copies of pre-fix state. Safe to remove once reviewers confirm the post-fix state.

---

### Primary-model documentation + conditioned-vs-baseline per-cell permutation (2026-04-14)

#### Documentation: 5a / 5a-gf2 established as primary
- Updated `CLAUDE.md`, `statistical_audit.md`, `README.md` to make clear that **Model 5a** (channel-level phase, `phaseVecLength ≥ 0.3`) and **Model 5a-gf2** (good-fit burst restriction + channel-level phase, `phaseVecLength ≥ 0.2`) are the primary inferential models. Model 3 series (binary phaseClass) is retained as a sensitivity/robustness check but is no longer primary.
- Rationale captured in statistical_audit.md Finding 23: phase is circular (binary binning loses information), grouping by measured phase keeps distinct conditions separate, 5a-gf2 adds a principled quality filter ("beta actually present during conditioning"), 5a/5a-gf2 reported together cover the analytic space.
- `CLAUDE.md` code comments updated: Model 3a → "(sensitivity)"; 5a → "(PRIMARY)"; 5a-gf2 → "(PRIMARY-complement)"; 5a-gf → "(diagnostic)".

#### New analysis: Conditioned vs Baseline per-cell permutation (companion to 5a-gf2)
- **Purpose**: 5a/5a-gf2 estimate the average dose effect across channels; the per-cell analysis describes how that average is distributed across individual channels — answering "broad shallow effect or a few strong responders?"
- **Method**: For each `(sid × channel × phaseDeg_round × dose)` cell, two-sample label-shuffle permutation (10k MC, median diff) of conditioned probes vs channel-level baseline probes. Bootstrap 95% CIs (2k resamples) for forest-plot uncertainty.
- **Filters (exactly match 5a-gf2)**: channel-level `phaseVecLength ≥ 0.2` AND good-fit burst restriction (`nGoodBeta ≥ 1`, applied to conditioned trials only — baselines exempt because they have no preceding burst). Per-cell minima: ≥10 baseline probes, ≥5 conditioned probes. Good-fit filter retains ~46% of conditioned trials (16,814 / 36,599). Final dataset: **101 cells across 7 subjects**.
- **Correction**: BH FDR within (subject × dose) — respects that channels are nested within subjects (electrode grid, anatomy, session noise shared within subject). Pooled FDR kept as a CSV diagnostic column only (near-identical to within-subject on this dataset).
- **Results**: c91479 has 3-4/4 cells FDR-sig at every dose (median effect ~48 µV at [5,inf), dominant responder). 9ab7ab shows a clear dose gradient (1/4 → 2/4 → 3/4 cells). ecb43e has 1/5 cells sig at [3,4]. 0b5a2e's [5,inf) effect strengthens post-filter (4 uncorrected-sig cells vs 1 pre-filter, median 24 µV) but its 13-test FDR family still rejects. Per-dose subject-presence: 2 / 3 / 2 of 7 subjects with ≥1 FDR-sig cell at [1,2] / [3,4] / [5,inf).
- **Forest plot**: rows sorted by measured phase (0° at top → 360° at bottom), subject as tie-breaker. Channel labels use raw channel numbers (subject-prefix stripped, e.g. 714→14) with "Subject N" format. Beta trigger channels (sin-fit reference channels per subject) highlighted via pink y-axis labels; dot color signals significance (grey = ns, orange = p<0.05 uncorrected, red = FDR q<0.05 within subject). Font bumped for manuscript legibility (base_size 17, plot dimensions 16×16 in).
- **Outputs**: `betaStim_cond_vs_base_perchan.csv` (master, includes `perm_q` within-subj and `perm_q_pooled` diagnostic columns), `betaStim_cond_vs_base_per_subject.csv`, `betaStim_cond_vs_base_subject_presence.csv`, `betaStim_cond_vs_base_summary.csv` (per-dose condensed counts, within-subject FDR), `betaStim_cond_vs_base_pooled_summary.csv`, `betaStim_cond_vs_base_forest.png/.eps` (single plot; pooled version dropped). New sections in `betaStim_within_subject_tables.docx`.
- Documented as Finding 24 in `statistical_audit.md`.

#### Per-dose condensed summary table added (2026-04-15)
- `betaStim_cond_vs_base_summary.csv`: one row per dose, collapses across subjects. Reports n_cells, n_subjects, n_channels, n_sig_uncorr, pct_sig_uncorr, n_sig_fdr (within-subject), pct_sig_fdr, median_effect, mean_effect.
- Key pattern: uncorrected rises monotonically (5→7→12 cells, 14.7%→21.2%→35.3%), FDR plateaus at ~21% for [3,4] and [5,inf) (0b5a2e's 4 uncorrected-sig cells at [5,inf) can't survive its 13-test family). Effect-size gradient clear: median 7.6→8.6→13.6 µV.

#### Phase results from primary models (documented 2026-04-15)
- cos_phase is trend-level in both 5a (p=0.090, β=4.92 µV) and 5a-gf2 (p=0.098, β=5.54 µV) — modulation on the 0°/180° axis.
- sin_phase non-significant in both (5a p=0.29, 5a-gf2 p=0.19; sign flips between models).
- No dose × phase interactions (all p>0.23). Phase effect is additive, not dose-gated.

#### 5a-gf2 effect sizes and .docx restructuring
- Computed total-variance Cohen's d for 5a-gf2 (mirroring 5a's existing computation). total_sd_gf2 = 22.3 µV, resid_sd_gf2 = 13.1 µV.
- `betaStim_statistical_tables.docx` restructured: PRIMARY MODELS (5a + 5a-gf2 with ANOVA, fixed effects, EMMs, effect sizes) at TOP; SUPPORTING / SENSITIVITY MODELS (3a, 3c, 3e) below with explicit demotion language.
- All stale "(primary)" references to Model 3c in R script comments updated to "(sensitivity)".

#### FDR correction strategy discussion
- Considered three correction families: pooled (all 34 cells per dose), within-subject (channels within a patient), within-channel (phases × doses within a channel).
- Went with **within-subject FDR** as the plotted primary: family sizes 2-13 give BH real correction power while respecting the dominant source of clustering. Within-channel families (1-6) reduce to near-uncorrected; pooled ignores clustering entirely. Pooled kept in CSV for reference; within-channel not implemented (near-uncorrected in practice).

---

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
- `CLAUDE.md` — Pipeline architecture, all model specifications, data structure notes, primary model rationale
- `statistical_audit.md` — 24 findings and updated recommendations (5a/5a-gf2 primary, 3-series sensitivity)
- `progress.md` — This file
- `README.md` — R script descriptions (updated to reflect 5a/5a-gf2 as primary)

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

- [ ] Manuscript figure integration: incorporate forest plot (Fig 8) and per-dose summary table into manuscript draft
- [ ] Sensitivity analysis: run cond-vs-base at alternative phaseVecLength thresholds (r ∈ {0, 0.1, 0.3}) to confirm stability
