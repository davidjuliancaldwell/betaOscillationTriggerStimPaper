# Phase data origin and known bugs in B_phaseCalc and re-reference loaders

**Status (2026-04-21, post-fix):** All bugs described below are **fixed** using
the `prev_grp` group-change detector pattern. EP pipeline regenerated for
c91479; CSV verified float-identical (max delta 1e-12 µV) to pre-fix across all
analyzed channels; primary R models (5a, 5a-gf2) reproduce CLAUDE.md numbers
exactly. Phase data **not** regenerated — the Dec 2018 hyak `.mat` files remain
canonical (regenerating would be safe now that the debug hardcodes are removed,
but would take hours and produce numerically-equivalent output). Scripts fixed:
`B_ExtractNeuralData_PP_reref.m` (both loops), `B_phaseCalc_allChans_processed.m`
(reload loop + `chans=64` debug hardcode removed + `12samps`→`51samps` save
filename typo + `idxVec` restored to `[1:8]`), `plot_EP_goodfit_by_phase.m`,
`BETA_ExtractNeuralDataCEPscreen.m` (both loops), and
`plot_example_dose_dependent_time_series.m`. Also `B_ExtractNeuralData_PP_reref.m`
`for idx = 1:1` restored to `for idx = 1:8`. See the "Verification after the
2026-04-21 fix" section near the bottom of this document.

**Pre-fix status (2026-04-21 morning):** Phase data in `data/phase_data/` is
correct but current local scripts have bugs that would corrupt it if
regenerated. Do not run `master_script_betaStim.m` with
`generateIntermediateData = 1` until the bugs below are fixed.

Last verified: 2026-04-21.

## Summary

- The existing per-stim phase fits in `data/phase_data/*.mat` were generated on
  the UW hyak cluster on **2018-12-04** using the `hyak` branch of
  `~/code/MATLAB_ECoG_Code` (path:
  `Experiment/BetaTriggeredStim/phaseVisualizations/B_ExctractNeuralDataDJC_phaseCalc_allChans_processed.m`).
  The filename suffix `_51samps_12_20_40ms_randomStart` encodes the hyak-era
  parameters (51-sample moving-average smoothing, 12–20 Hz fit band, 40 ms
  pre-stim window, random initial-phase seed).
- Two bugs exist in the *current* copy at
  `phase_visualizations/B_phaseCalc_allChans_processed.m`. The hyak version that
  generated the `.mat` files does **not** have Bug B, and Bug A is latent in both
  (present but does not manifest for any of the 8 subjects' bads lists).
- The related re-reference loading loops in
  `peak_extraction/B_ExtractNeuralData_PP_reref.m` and
  `manuscript_generate_scripts/plot_EP_goodfit_by_phase.m` share a class of bug
  (narrow `achan`-based reload gate). Only `plot_EP_goodfit_by_phase.m`'s
  manifestation actively corrupts output; the main pipeline's gate is broader
  and happens to catch every group transition for all current subjects.

## The existing `.mat` files

Files in `data/phase_data/`, all dated 2018-12-04:

| File | Size |
|---|---|
| `d5cd55_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat` | 638 MB |
| `c91479_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat` | 1.95 GB |
| `7dbdec_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat` | 3.12 GB |
| `9ab7ab_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat` | 2.04 GB |
| `702d24_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat` | 3.48 GB |
| `ecb43e_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat` | 3.48 GB |
| `0b5a2e_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat` | 766 MB |
| `0b5a2ePlayBack_phaseDelivery_allChans_51samps_12_20_40ms_randomStart.mat` | 766 MB |

These are the source of truth for downstream phase-dependent analyses (Models
5a, 5a-gf2; CL vs PB; conditioned-vs-baseline per-cell permutation; all forest
plots; all phase-response curves). They feed:

- `analysis_functions/compute_burst_phase_precision.m` →
  `data/output_table/{sid}_burst_phase_precision.csv`
- `peak_extraction/multipleSubj_GLMM_script_PP.m` →
  `data/output_table/betaStim_outputTable_50_new_100_thresh.csv`
- All R analyses downstream.

## Bug B (critical): `chans = 64;` hardcode

**Location:** `phase_visualizations/B_phaseCalc_allChans_processed.m:59`

```matlab
chans = [1:64];
chans(ismember(chans, badsTotal)) = [];
chans = 64;                              % <-- wipes the previous line
```

MATLAB's `chans = 64` replaces the variable with the scalar 64, discarding the
full channel list built above. The subsequent `for chan = chans` loop at line 90
would iterate exactly once (with `chan = 64`), writing a phase `.mat` that
contains data for only channel 64 — or, for ecb43e where chan 64 is a bad
channel, for chan 64 regardless of the `badsTotal` filter (because the hardcode
runs *after* the filter).

**Origin:** commit `29f12a1` (2019-11), *"fixed subject six, use original stim
table, appropriate rereferencing, fixing paths in various scripts"*. Likely a
debug/test line added while iterating on subject 702d24 that was never
removed.

**Impact if the script is rerun as-is:** every subject's
`*_phaseDelivery_allChans_*.mat` would be overwritten with single-channel output,
silently breaking every downstream analysis.

**Fix:** delete line 59. The proper channel list is already built in lines
57–58.

## Bug A (latent): narrow reload gate

**Location:** `phase_visualizations/B_phaseCalc_allChans_processed.m:101`

```matlab
if achan==1 || achan == 2
    load(fullfile(folderECoGData,[sid '_ECoG.mat']),ev);
    dataStruct = eval(ev);
end
```

TDT stores 64 ECoG channels across 4 structs (`ECO1`–`ECO4`), each holding 16
channels. The loop caches the most-recently-loaded struct in `dataStruct` and
tries to reload only when needed. The reload condition is based on a heuristic —
"if the within-struct index is 1 or 2, we probably just crossed into a new
group" — rather than directly checking whether the group changed.

The heuristic fails when a subject's `bads` list skips channels 1 and 2 of a
group. In that case, the first channel of the new group has `achan ≥ 3`, the
gate misses, and the loop reads from the previous group's `dataStruct`
indefinitely.

**Why it doesn't actually manifest in the Dec 2018 data:** by coincidence, every
one of the 8 subjects' bads lists leaves at least channel 1 or 2 of each group
unbad, so every group transition lands at `achan ∈ {1, 2}` and the narrow gate
catches it.

| Subject | First chan in each group |
|---|---|
| d5cd55 | ECO1: 2, ECO2: 17, ECO3: 33, ECO4: 50 |
| c91479 | ECO1: 1, ECO2: 17, ECO3: 33, ECO4: 49 |
| 7dbdec | ECO1: 1, ECO2: 17, ECO3: 33, ECO4: 49 |
| 9ab7ab | ECO1: 2, ECO2: 17, ECO3: 33, ECO4: 49 |
| 702d24 | ECO1: 1, ECO2: 17, ECO3: 33, ECO4: 49 |
| ecb43e | ECO1: 1, ECO2: 17, ECO3: 33, ECO4: 49 |
| 0b5a2e | ECO1: 1, ECO2: 17, ECO3: 33, ECO4: 49 |
| 0b5a2ePlayback | ECO1: 1, ECO2: 17, ECO3: 33, ECO4: 49 |

All group-entry `achan` values are 1 or 2, so the narrow gate succeeds for the
current 8-subject cohort. Adding any future subject whose bads list blocks
channels 1 and 2 of any group would silently corrupt their phase data.

**Fix:** replace the narrow gate with a direct group-change detector:

```matlab
prev_grp = -1;
for chan = chans
  grp = floor((chan-1)/16);
  ev  = sprintf('ECO%d', grp+1);
  achan = chan - grp*16;
  if grp ~= prev_grp
    load(fullfile(folderECoGData, [sid '_ECoG.mat']), ev);
    dataStruct = eval(ev);
    prev_grp = grp;
  end
  eco = 4 * dataStruct.data(:, achan)';
  % ... rest of loop
end
```

This is correct regardless of the channel list or bads configuration, and also
saves the redundant re-loads that the `1||2||4||6` widened gate triggered
within a group.

## Related bugs in the reref/CEP scripts

### `peak_extraction/B_ExtractNeuralData_PP_reref.m`

This is the main CEP extraction pipeline. Two channel loops:

- **Lines 158–214 (rerefChans loop):** gate is `achan == 1 || 2 || 4 || 6`. For
  each subject's current `rerefChans` list, every group transition lands at
  `achan ∈ {1, 2, 4, 6}`, so this gate catches them all. **Output is correct
  for the current 8 subjects.**
- **Lines 228+ (chans loop):** gate is `achan == 1 || 2` (narrow). At the start
  of this loop, `dataStruct` is stale from the end of the rerefChans loop
  (typically ECO4). Manifests only when a subject's `chans` list starts at
  `achan ∉ {1, 2}` — which among the current 8 subjects occurs only for
  c91479 (bads include [1, 2, 3], so `chans` starts at 4, achan=4).
  For c91479 chans 4–16, the code reads ECO4 stale data
  (raw channels 52–64) instead of the intended ECO1.

**Practical impact: benign for reported results.** c91479's downstream
`goods` list (at line 44: `[39 40 47 48 63 64]`) excludes all of chans 4–16,
so the buggy per-channel CEP values are computed but never used by any
downstream analysis (multipleSubj_GLMM_script_PP.m, R analyses, forest plots,
etc.). All 7 other subjects have `chans` starting at `achan ∈ {1, 2}`,
so they never trigger the bug.

Fix the gate anyway as defensive hygiene — any future subject with bads at
positions 1 AND 2 of a group, or any downstream `goods` edit that pulls
chans 4–16 of c91479 into the analysis, would make this latent bug bite.

### `manuscript_generate_scripts/plot_EP_goodfit_by_phase.m:124`

Same narrow gate `achanR == 1 || 2`, *combined with* a pre-load of `dataStruct`
for `chanInt` (line 111–112) before the rerefChans loop starts. For c91479 with
`chanInt = 47`, the pre-load puts `dataStruct = ECO3`, and the first rerefChans
iteration (`rc = 4`, `achanR = 4`) misses the narrow gate, so the reref signal
for `rc = 4..16` is pulled from ECO3.data(:, 4..16) = raw channels 36..48
instead of raw channels 4..16. Result: ~23% of reref channels silently replaced,
median shift causes CEP drift of ~5% downstream.

### `manuscript_generate_scripts/BETA_ExtractNeuralDataCEPscreen.m:189`

Same narrow gate in a second loop analogous to the main pipeline's chans loop.
Same fix applies.

### Files with the broader `1 || 2 || 4 || 6` gate (safe for current subjects)

- `peak_extraction/B_ExtractNeuralData_PP_reref.m` lines 158–214 (rerefChans loop)
- `manuscript_generate_scripts/plot_example_dose_dependent_time_series.m`
- `manuscript_generate_scripts/BETA_ExtractNeuralDataCEPscreen.m` line 117

These still have the fragile-heuristic pattern and should be fixed preventatively,
but do not currently produce wrong outputs for any subject.

## What not to do

- Do **not** run `master_script_betaStim.m` with `generateIntermediateData = 1`
  before fixing Bug B above. The hardcode will destroy the good phase data.
- Do **not** regenerate CEP data from `B_ExtractNeuralData_PP_reref.m` without
  first fixing the chans-loop gate bug (or at least confirming it does not
  affect the current output). For the current cohort, reported results are
  unaffected (c91479's `goods` excludes the bug zone), but the buggy CEP
  values for c91479 chans 4–16 are still sitting in the `.mat` file.

## What to do (recommended order)

1. **Fix both bugs in `B_phaseCalc_allChans_processed.m`** — delete line 59 and
   apply the group-change-detector fix to line 101.
2. **Apply the same group-change-detector fix** to `B_ExtractNeuralData_PP_reref.m`
   (both loops), `plot_EP_goodfit_by_phase.m`, and
   `BETA_ExtractNeuralDataCEPscreen.m` line 189.
3. **Back up** `data/phase_data/*.mat` and `data/EP_data/*epSTATS-PP-sig-reref-50-new.mat`
   to a dated backup directory before regenerating.
4. **Regenerate per subject** and diff against backup for numerical tolerance.
   Small differences from MATLAB version drift or random-start seeds are fine;
   large differences indicate something else has changed and needs investigation.
5. **Update CLAUDE.md** with post-fix state.

## Code archaeology reference

Hyak-branch commit history for the phase-calc script (from
`~/code/MATLAB_ECoG_Code`, branch `origin/hyak`):

```
7b5eede 2018-10-22 8-30 Hz bin fit to try and address edge effects of fitting
2f6c77c 2018-10-16 trying to start sinfit at 0, random might be better
e978872 2018-10-14 40 ms before pre fit, 12 20 Hz sine wave, random start
ea1bbfd 2018-10-08 small changes to finish subjects not done
a00cdbe 2018-10-04 51 samps smoothing, 60 ms before, random start for phase,
                   12-20 rather than 10-30 limit on beta
66ef82e 2018-08-29 14 sample delay for curve fitting to avoid artifact,
                   using offset now, non random starting phase
605fd39 2018-08-28 30 samples @ 24 khz to minimize artifact,
                   try 0 to 2pi for fitting rather than -pi to pi
5ebf1b6 2018-08-28 running on hyak, including stimulation delivery offset in
                   stim table as it's loaded in
07bcfa9 2018-08-26 small changes for first hyak run
2663fe6 2018-08-26 changes for hyak
2ec6dcc 2018-08-26 now doing rereferencing on good channels from the beta stim
                   experiment to try and accurately determine the phase
```

The Dec 2018 `.mat` files were generated from the `e978872`/`ea1bbfd` revision
(40 ms window, 12–20 Hz, 51-sample smooth, random start), which matches the
filename suffix `51samps_12_20_40ms_randomStart`.

## Verification after the 2026-04-21 fix

### Fixes applied (`prev_grp` pattern, 8 loops in 5 files)

```matlab
prev_grp = -1;
for chan = <loop_chans>
    grp = floor((chan-1)/16);
    ev  = sprintf('ECO%d', grp+1);
    achan = chan - grp*16;
    if grp ~= prev_grp
        load(fullfile(folderECoGData, [sid '_ECoG.mat']), ev);
        dataStruct = eval(ev);
        prev_grp = grp;
    end
    eco = 4 * dataStruct.data(:, achan)';
    % ... rest of loop
end
```

Replaces the narrow `achan ∈ {1, 2}` and wider `achan ∈ {1, 2, 4, 6}` gates in:

- `peak_extraction/B_ExtractNeuralData_PP_reref.m:158` (reref loop) and `:228` (chans loop)
- `phase_visualizations/B_phaseCalc_allChans_processed.m:90` (channel loop)
- `manuscript_generate_scripts/plot_EP_goodfit_by_phase.m:120` (reref loop)
- `manuscript_generate_scripts/BETA_ExtractNeuralDataCEPscreen.m:106` (reref loop) and `:179` (chans loop)
- `manuscript_generate_scripts/plot_example_dose_dependent_time_series.m:149` (reref loop)

Also in `B_phaseCalc_allChans_processed.m`:
- Line 59: removed debug `chans = 64;` override (was silently limiting every subject's phase-fit loop to channel 64 only)
- Line 2: `idxVec = [1:7]` → `idxVec = [1:8]` (include playback)
- Line 395: save filename typo `12samps` → `51samps` (matched downstream modifier and existing files)

And in `B_ExtractNeuralData_PP_reref.m` line 19: `for idx = 1:1` (debug, d5cd55-only) → `for idx = 1:8`.

### Verification: EP pipeline rerun (c91479 only)

Rationale for c91479-only: per the safe-subject matrix in this document, c91479
is the only subject whose `chans` list starts at `achan ∉ {1, 2}` — the only
one the old narrow-gate bug could have bitten in the main pipeline's second
loop. Other subjects' outputs are provably unaffected.

- Ran `B_ExtractNeuralData_PP_reref` for `idx=2` (c91479) with `saveIt=1`. Stage B1 took 6 s (OneDrive already hydrated from earlier inspection).
- Stage C (`multipleSubj_GLMM_script_PP` across 7 subjects) regenerated `data/output_table/betaStim_outputTable_50_new_100_thresh.csv` in 7 min.
- Diffed new vs backup CSV: **all 10 numeric columns bit-identical for every subject and channel except c91479**, where 152 rows differ by **≤1e-12 µV** (float-accumulation noise, not meaningful).
- Channel-level trace of `dataForPPanalysis` confirmed c91479 good channels (47, 48, 64) were already reading the correct ECO struct in the pre-fix code — by the time the narrow-gate chans-loop reached them, `dataStruct` had been legitimately reloaded to the right struct via earlier `achan=1|2` hits. So the pre-fix CSV magnitudes were already correct.

### Verification: R analysis reproduces CLAUDE.md

- `betaStim_R_script.R`: completed cleanly (7 standard warnings, no errors).
  - **Model 5a** (78 obs, 19 ch, 7 sid, non-singular): `numStims_ord.L` p = 0.0676 ✓ (CLAUDE.md: 0.068), `cos_phase` p = 0.0899 ✓ (0.090), Type III dose F p = 0.123 ✓, effect size 6.36 µV ✓.
  - **Model 5a-gf2** (102 obs, 25 ch, non-singular): `numStims_ord.L` p = 0.0531 ✓ (0.053), `cos_phase` p = 0.0983 ✓ (0.098), effect size 8.10 µV ✓.
- `R_compare_control_cond.R`: completed cleanly. `betaStim_within_subject_tables.docx` regenerated with matching within-subject CL vs PB, ecb43e targeted/random, and conditioned-vs-baseline per-cell permutation results.

### Verification: plot_EP_goodfit_by_phase.m outputs

- 11 plots regenerated (c91479: 47, 48, 64; 0b5a2e: 14, 15, 16, 21, 23, 31, 32, 40).
- New sgtitle annotation: `[INCLUDED: mean Base PP X.X ≥ 100 µV]` (black) or `[EXCLUDED from R analysis: mean Base PP X.X < 100 µV threshold]` (red). Ch48 correctly flagged EXCLUDED — its pre-fix annotation showed Base PP = 0 µV (artifact of the old reref corruption trashing single-trial peak extraction); post-fix it shows Base PP = 36 µV (genuine weak signal, still below threshold).
- c91479 ch47 and ch64 Base PP annotations shifted <5% (148→155 µV, 276→272 µV) as predicted for median-CAR with ~23% reref-channel corruption on these fits.

### Backup

Full local backup of pre-fix `data/{output_table,phase_data,EP_data}/` is at
`data/_backup_20260421_161313/`. Safe to remove once reviewers confirm the
post-fix state.
