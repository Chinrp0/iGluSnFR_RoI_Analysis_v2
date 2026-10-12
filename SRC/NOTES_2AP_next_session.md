# 2AP (paired-pulse / PPF) analysis — handoff notes

Dataset: `E:\...\iglu3fast_NGR\PPF\PPF_2026-06-30--12\E_raw_mean`
Pilot: iGluSnFR3fast, Doc2b WT vs R213W. Stat unit = **coverslip**.
WT = 2 coverslips (Cs1: c2/c4, Cs2: c2/c3), R213W = 1 coverslip (Cs1: c2/c5), 1 animal/genotype.
ISIs: 30, 60, 90, 120, 250, 500, 1000 ms. Mostly 1 sweep/condition (a few have -01/-02).
Latest figures: `...\PPF\E_raw_mean_figures_20260707_7\`

---

## KEY FINDING (context for everything below)
The apparent **PPR < 1 ("depression") at all ISIs was a selection / regression-to-the-mean
artifact**, NOT biology. PPR was a mean of per-ROI ratios computed only on pulse-1
responders (gating on A1 inflates the denominator). Proof: gating on A1 gives A2/A1 < 1
AND gating on A2 gives A1/A2 < 1 — impossible for a real effect. The unbiased
**ratio-of-means** (coverslip average trace) shows real **facilitation peaking ~60–120 ms**
(R213W/Cs1 pk ≈ 1.71 at 90 ms) decaying to ~0.85 by 1000 ms (indicator rundown). The 30 ms
point is unreliable because the two pulses fuse at 6-frame spacing.

## DONE this session
- **Schmitt alignment**: 2AP now uses the SAME dual-threshold validation as 1AP
  (`window_peak` counts run >= LOWER threshold; responder needs peak >= UPPER).
  Config: `lower_threshold_sigma` now used; `min_frames_above_upper` -> `min_frames_above_lower` (=2).
  Effect: pulse-1 responders ~doubled (recovered fast transients the old 2-frames-above-upper rule discarded).
- **AUC metric** added (`auc1`, `auc2`, `ppr_auc`, `config.twoAP.auc_window_ms=40`); in CSV + summary;
  replaced the unstable sync panel in the PPR-vs-ISI figure (sync still in CSV).
- **Montage**: wider window (`montage_pre_ms=250`, `montage_post_ms=300`), shaded peak-search
  windows, lower-threshold line, `max=/auc=` per-tile titles.
- **Coverslip grand-average QC figure** (`plot_twoap_cs_average` + `trace_ppr`): mean of ALL ROIs
  per coverslip x ISI; title pk/auc = ratio-of-means PPR. This is the honest view.

Files touched: `config/tracenorm_config.m`, `tracenorm/twoap_analyzer.m`, `twoap_batch_analysis.m`.

### DONE 2026-07-14 (ratio-of-means primary + amplitude figure)
- **Ratio-of-means is now the PRIMARY PPR estimator** (TODO #1, done). New `compute_rom_by_isi`
  runs `trace_ppr` on each `gd.cs_avg` trace -> per-coverslip peak/AUC PPR + absolute A1/A2/AUC
  at each ISI. `plot_twoap_ppr` rebuilt: panels 1-2 = ratio-of-means peak & AUC (bold "primary"),
  panel 3 = per-ROI responder-gated peak-max (italic "diagnostic", labelled selection-biased).
  `trace_ppr` now returns `[ppr_pk, ppr_au, A1, A2, au1, au2]`. Results struct gains
  `ppr_by_isi.<group>.ratio_of_means`; summary prints ROM peak/auc + A1/A2 with gated alongside.
  Validated on real data: R213W ROM peak = 1.71 @ 90 ms (matches prior finding); gated stuck 0.6-0.9.
- **NEW amplitude figure** `plot_twoap_amplitudes` (`plot_handles.amplitudes`): GROUPED BAR chart,
  one panel/group, three bars per ISI = Pulse 1 (A1) / Pulse 2 incremental (above pre-stim2
  residual) / Pulse 2 from-baseline (above pre-stim1 baseline). mean+/-SEM across coverslips,
  shared y-axis. `trace_ppr` now also returns `A2_abs` (7th out); `compute_rom_by_isi` collects
  `a2_abs`. Shows facilitation directly (A1 flat ~0.009-0.011, A2-incremental peaks ~90 ms) AND
  the residual pedestal: incremental << from-baseline at 30 ms (0.009 vs 0.017), converged by
  250 ms. NB incremental subtracts the mean of the 2 frames just before stim2 (NOT the trough).

Shared plot helpers added: `draw_ppr_series`, `finish_ppr_axis`.
- **Gating diagnostic UNIFIED with the primary ROM (2026-07-14)**: previously the gating figure's
  "no-gate" was a POOLED per-ROI ratio (sumA2/sumA1) which did NOT equal the Fig-1 avg-trace ROM
  (confusingly both were called "ratio of means"; e.g. R213W 90ms Fig1=1.71 vs old no-gate=1.20).
  Now `compute_gating_bias` builds per-coverslip AVERAGE TRACES restricted to each responder subset
  (`gd.cs_avg_gate`, accumulated in aggregate_group_2ap alongside cs_avg) and takes peak2/peak1 of
  each -> **no-gate is now IDENTICAL to the primary ROM curve** (verified exact match), gate_A1 <
  no-gate < gate_A2 fan preserved and wider. gate_A2 reconstructs P2 responders as
  `peak2_amp >= upper_threshold - prestim_base` (= upper_sigma*sd; no config needed). Helper
  `colsum` added. Fig ylabel now "PPR (A2/A1 of avg trace)", no-gate label "= primary ROM".
- **cs_avg now excludes noise-QC-failed ROIs** (`qc = isfinite(upper_threshold)` = analyzable &
  not noise-excluded; NOT response-gated, so no selection bias). Tightens SEM on the primary ROM +
  amplitude figures; mean barely moves (R213W 90ms 1.71->1.70), confirming noise was averaging out.
- **Figure saving refactored to shared io/ helpers**: `io/make_figures_dir.m` (dated auto-increment
  dir next to data) + `io/save_open_figures.m` (PNG+FIG of given/all open figs). `run_pipeline_fixed`
  now calls these (local `make_output_dir`/`save_all_figures` deleted). `twoap_batch_analysis` gained
  `options.saveFigures` (default FALSE so run_pipeline_fixed doesn't double-save) -> saves its own
  handles + sets `results.output_dir`. Normal run path (run_pipeline_fixed) already saved the new
  figures automatically; opt-in is for direct/standalone calls.
- **Responder QC tightened in `twoap_analyzer.m`** (noise ROIs like ROI189 with negative PPR):
  1) `window_peak` now uses **findpeaks** (genuine local max >= upper + MinPeakDistance), matching 1AP,
     instead of a plain window max -> rejects monotonic noise ramps. New config `min_peak_distance_ms=10`.
  2) **Baseline-excursion floor**: pulse-1 responder must have `peak1_amp > k * prctile(pre-stim
     excursion, pct)`. New config `baseline_excursion_pct=99`, `peak_over_baseline_k=1.0` (k=0 disables).
     Exported `res.base_exc_ceil` for QC.
  3) **Negative A2 clip**: if incremental A2 (or AUC2) < 0 (pulse 2 didn't respond / fused on P1 decay),
     clip to 0 and set is_responder2=false -> per-ROI PPR floored at 0, no more negative ratios.
  Effect on real data: WT responders 2367->2216 (-6.4%), R213W 1278->1213 (-5.1%); negative peak-PPR
  84->0, negative AUC-PPR 715->0. IMPORTANT: the PRIMARY ROM figures (ppr_vs_isi, amplitudes) use
  cs_avg = mean of ALL ROIs (ungated by design) so they are UNCHANGED; QC only cleans responder-based
  outputs (montage, counts, gated diagnostic, gating-bias, per-ROI CSV). Diagnostic scripts used to
  pick the metric are in scratchpad (baseline-event-rate FAILED to discriminate; std-SNR ~= MAD-SNR;
  peak-vs-baseline-excursion was the discriminating metric). NB many 30 ms auc=0.00 in montage = P2
  fused on P1 decay -> clipped; item #4 (decay subtraction) would recover those.

---

## TODO next session

### 1. ~~Ratio of means as the PRIMARY PPR estimator~~ — DONE (see above)

### 2. ~~Explore gating on A1 vs A2 vs either/both (bias exploration)~~ — DONE 2026-07-14
- `compute_gating_bias(gd, config)` + `plot_twoap_gating_bias` (`plot_handles.gating_bias`) +
  `print_gating_bias`. Population PPR = sum(A2)/sum(A1) pooled per group x ISI under 5 gates:
  no_gate / gate_A1 (is_responder1) / gate_A2 / gate_either / gate_both. One figure panel per group
  + a printed table in the run log. Result confirmed the artifact at EVERY ISI:
  **gate-A1 < no-gate < gate-A2** (e.g. R213W 90 ms: 0.78 < 1.20 < 1.60; 30 ms: 0.47 < 1.22 < 2.09).
- gate_A2 reconstructs the pulse-2 responder as `peak2_amp >= upper_threshold_sigma*prestim_sd`
  (NOT is_responder2, per the caveat below) — A2 is already above the local pre-stim2 residual so
  the k*sigma comparison is direct.
- Original caveat (now handled): `is_responder2` is unreliable at short ISI (pulse 2 sits on
  pulse-1's elevated decay, always above the LOWER threshold), so gate_A2 uses the reconstruction.

### 3. ~~Rework the AUC calculation~~ — DONE 2026-07-21 (total-charge, encompassing window)
- **Design decision** (with user): the AUC is no longer a per-pulse ratio. It is now a SINGLE
  TOTAL CHARGE over ONE window that ENCOMPASSES BOTH pulses: `q_pair = ∫ dF/F over
  [stim1, stim2 + auc_post_ms]`, above the SINGLE pre-stim1 baseline (no residual, no per-pulse
  split). New config `config.twoAP.auc_post_ms = 250` (matches the 1AP auc window; captures the
  full iGluSnFR3fast decay after pulse 2 and stays inside 600 frames at 1000 ms: stim2~=466 -> end
  516). Replaced `auc_window_ms`. This kills the old short-ISI pathologies outright: no per-pulse
  window shrink (`min(auc_frames, isi-1)` is gone), no negative-AUC2-above-residual, no NaN blowups.
- **Why total charge, not a ratio**: a single encompassing window gives one number = total glutamate
  charge of the pair; you can't recover A2/A1 from a blended integral. The window MUST scale with ISI
  (it has to reach stim2), so its length = ISI + auc_post_ms and is NOT the same across ISIs — that
  is unavoidable, not a flaw. User chose **absolute Q, NO normalization** (options considered:
  per-coverslip long-ISI pulse-1 charge, per-trace pulse-1, from the 1AP dataset — all deferred).
- **Code**: `twoap_analyzer.m` computes `res.q_pair` (per ROI, UNGATED — for all analyzable ROIs).
  `trace_ppr` returns `[ppr_pk, q_pair, A1, A2, A2_abs]` (dropped ppr_au/au1/au2). `compute_rom_by_isi`
  carries `q_pair`. Figures/CSV/summary updated: PPR-fig **panel 2 is now "Total charge Q" on its own
  y-axis** (absolute dF/F*ms, no PPR=1 line, NOT shared with the ratio panels 1&3). Montage +
  cs_average tile titles show `Q=` (was `AUCr=`). CSV col `q_pair` (dropped auc1/auc2/ppr_auc).
  Summary prints `Q` per ISI. Removed `ppr_by_isi.<grp>.auc` and `res.n_ppr_valid_auc` (unused).
- **Verified**: checkcode clean; synthetic trace q_pair == independent trapz exactly at every ISI;
  full real-data run OK (15465/15540 finite q_pair; 75 NaN = zero-noise skipped ROIs).
- **KEY CAVEAT (drift confound)**: on the real coverslip AVERAGE traces, R213W Q DECREASES with ISI
  (1.70 @30ms -> 0.62 @1000ms). This is partly METHODOLOGICAL, not biology: integrating above a FIXED
  pre-stim1 baseline over a window that LENGTHENS with ISI lets photobleaching/rundown baseline droop
  after the stim accumulate NEGATIVELY over the long-ISI windows. Consistent with the documented
  long-ISI rundown, but it means raw Q-vs-ISI is confounded. To make Q-vs-ISI a clean facilitation
  readout, revisit normalization (deferred above) AND/OR bleach-correct the baseline before integrating.

### 4. Decay subtraction at short ISI (30–90 ms)
- Subtract pulse-1's decay under pulse 2 before measuring A2, so the fused 30 ms response is
  measurable. Fit pulse-1 tail (single/double exp; `tracenorm/fit_decay_double.m` exists) or
  template-subtract the average pulse-1 waveform, then measure A2 above the extrapolation.
- This is the one change that makes the 30–90 ms points trustworthy for BOTH peak and AUC PPR.

### 5. ~~Consistent y-axis on the coverslip average-trace figure~~ — DONE 2026-07-14 (option a)
- `plot_twoap_cs_average` now links y per coverslip ROW (`linkaxes(rowAxes,'y')`), so pulse-2 vs
  pulse-1 heights are comparable across ISIs within a coverslip without a global scale squashing the
  smaller coverslip. Stim markers switched to `xline` so they stay full-height after linking.
- Optional future nicety (not done): option (c) NORMALIZED variant (divide each trace by its
  pulse-1 peak -> every A1=1.0, facilitation readable on a common 0-2 scale).

### 6. ~~Figure labeling / consistency pass~~ — DONE 2026-07-14
- Every 2AP figure now states its ROI basis + stat unit in the title: PPR primary panels &
  amplitudes = "all ROIs" (avg trace); PPR panel 3 = "P1 responders only"; gating = "pooled ROIs";
  cs_average = "mean of ALL ROIs"; montage = "example pulse-1 responders (one reference file)".
- PPR-vs-ISI figure: y-min pinned to 0 + shared y across its 3 panels. Gating figure: shared y
  across group panels, y-min 0 (matches PPR fig). Terminology unified: PPR = A2/A1, AUCr = AUC
  ratio (montage + cs_average tiles relabeled from the ambiguous max=/auc= that read as amplitude).
  Stim colors consistent everywhere (blue=stim1, red=stim2).

---

## Data caveats to keep in mind
- Pilot N (2 vs 1 coverslip, 1 animal/genotype): NO genotype stats. Use for method-locking +
  effect-size/variance to power the real study. Frame figures accordingly.
- Long-ISI PPR < 1 (~0.85) with well-separated peaks = likely indicator photobleaching/rundown,
  a real but separate effect from PPR. Consider bleach correction.
- iGluSnFR PPR magnitude (~1.2–1.7) is compressed vs electrophysiology (2–3) — indicator
  saturation + averaging over boutons of mixed release probability. Direction is what matters here.
