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
  dir next to data) + `io/save_open_figures.m` (PNG+FIG of given/all open figs). `run_pipeline`
  now calls these (local `make_output_dir`/`save_all_figures` deleted). `twoap_batch_analysis` gained
  `options.saveFigures` (default FALSE so run_pipeline doesn't double-save) -> saves its own
  handles + sets `results.output_dir`. Normal run path (run_pipeline) already saved the new
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

### DONE 2026-10-15 (new PRIMARY PPR: pooled gate-either, decay-subtracted)
- **Shared amplitude function** `tracenorm/twoap_pulse_amplitudes.m`: A1/A2 five ways (peak_max,
  win_mean [amp_window_ms 5-20], sync, p1lock, decay_sub). Used by twoap_analyzer (per ROI) AND
  trace_ppr (avg traces) -> same frames / local-max rule / baseline (trace_ppr now uses prestim_idx,
  i.e. drops the 100 ms guard). Effect on old ROM curve <= 0.01. Test: `tests/test_twoap_pulse_amplitudes.m`.
- **Decay subtraction**: pulse-1 template = group's longest-ISI (>= 250 ms) cs_avg trace, stim1..stim2,
  peak-normalised, zeroed once < `decay_template_min_frac` (0.05). Scaled to each trace's pre-stim2
  residual, subtracted in [stim2, stim2+peak_search_ms]; A2 = largest local max of the subtracted window.
  Synthetic: peak_max underestimates true PPR 18-23% at 30 ms; decay_sub exact.
- **PRIMARY PPR** (fig panel 1, summary, `ppr_by_isi.<grp>.primary`, `twoap_per_coverslip.csv`):
  per coverslip x ISI, SumA2_decay / SumA1 over **gate-either** ROIs (is_responder1 OR incremental
  peak2_amp >= upper_sigma*sd), mean across coverslips. Panel 3 = QC ref: ungated avg trace, decay-sub.
- **Diagnostic** (`plot_twoap_ppr_diagnostic`, `_distributions`, `_decay_examples`, `print_ppr_diagnostic`):
  mean of per-ROI ratios is inflated 1.2-11x (small A1) -> never use; pooled ~= avg trace; per-ROI
  p1lock biased low (latency picked on noise); sync unstable at short ISI.
- Per-ROI CSV gained is_responder2_incr, gate_either, a2_decay_sub, ppr_decay_sub.
- Datasets in scope: 3fast NGR 5 ms + SF4f 10 ms (skip the 5 ms SF4f file); pdgfr lacks Recordings.csv.

### FIXED 2026-10-15: pulse-1 gate collapsed at 30 ms ISI
- At 30 ms the pulse-1 window is cut to stim2-1 (10 ms frames: 0/10/20 ms). Pulse 1 peaks ~20 ms =
  the window's last frame, and findpeaks never counts endpoints -> is_responder1 rate 5-7% at 30 ms vs
  ~45% at other ISIs (SF4f). Gate-either then ~= gate-A2 -> PPR biased UP (SF4f 30 ms 2.1 vs 1.3-1.6
  no-gate). 3fast 30 ms milder (15-18% vs 20-34%). Fix candidate: when the window is truncated by
  stim2, allow the window-end max to pass the P1 gate (amplitude + Schmitt run), or drop 30 ms.
- **Fix**: when stim2 cuts the pulse-1 window short, its LAST frame counts as a peak candidate if it
  rose into it (twoap_analyzer window_peak `allow_end` for the gate, twoap_pulse_amplitudes local_peak
  for A1); same >= upper / Schmitt-run / baseline-floor checks. NOT caused by the 10.0075 ms frame time
  (30/10.0075 and 30/10 both round to 3 frames; Recordings.csv refined stims are 3 frames apart).
  Effect: only 30 ms rows change. SF4f P1 responders 5-7% -> 18-52%; primary PPR 2.09/2.05/2.12 ->
  1.67/1.44/1.53 (R213W/Doc2b-WT/WT; noGate 1.61/1.48/1.29). 3fast 30 ms -0.08 to -0.17. SF4f 30 ms is
  still low-confidence: pulses fuse in many ROIs (trace rising through stim2).

### FINDING 2026-10-15: what gate-either excludes (3fast NGR)
- Diagnostic figs (gating-bias, estimator comparison, per-ROI distributions) now behind
  `options.diagnosticPlots` (default false). Amplitude fig (Fig 2) "incremental" is still LEGACY
  (peak-max above flat residual, ungated avg trace), not decay-subtracted.
- Excluded = QC-passed ROIs failing both gates: ~55-70% of ROIs, ~half of total SumA1/SumA2. NOT noise:
  their avg trace responds to both pulses, but A1 is ~1/3 of included ROIs (0.006-0.008 vs 0.016-0.022 dF/F)
  and they FACILITATE MORE (avg-trace decay-sub PPR, excl vs incl: WT 1.75/1.42/1.09/1.25 vs
  1.39/1.08/1.07/0.97; R213W 2.14/1.62/1.69/1.83 vs 1.38/1.07/1.16/1.13 at 30/60/90/120 ms; ~equal or
  lower at 1000 ms). Consistent with low-Pr boutons. Their avg trace also decays slower after pulse 2.
- So no-gate PPR > gate-either PPR because the gate drops a weak, strongly facilitating subpopulation,
  not because of regression to the mean.
- Per-ROI pooled PPR UNDERSTATES the excluded set (1.11 vs 1.42 avg trace at WT 60 ms): peak-picking a
  sub-threshold trace adds a noise bump to both A1 and A2, compressing the ratio toward 1. For included
  ROIs pooled == avg trace (1.08 vs 1.08). => sub-threshold ROIs need average-trace estimates.
- **Reporting decision (user, 2026-10-15)**: gate-either 5 sigma stays PRIMARY; the excluded ROIs are
  reported as a separate **sub-threshold** population, PPR from their per-coverslip AVERAGE trace
  (decay-sub) -> dashed curve in PPR fig panel 1, `ppr_diag.sub_threshold`, CSV `ppr_sub_threshold`,
  summary `subThr`. Excluded ROIs are 1.6x noisier at baseline (median sd 0.0081 vs 0.0050) AND smaller
  (median ~3 sigma vs ~7). Fig 2 (amplitudes) now = mean A1 / decay-sub A2 per gate-either ROI per coverslip
  (ratio == primary PPR); legacy incremental / from-baseline bars removed.

### DONE 2026-10-15: every 2AP figure on the decay-subtracted basis
- New helper `trace_ppr_primary` (avg-trace A2/A1 with the group's decay template; peak-max only if
  the group has no template). Used by the gating figure and the coverslip-average tiles.
- **Gating figure** now decay-subtracted: its no-gate curve == PPR fig panel 3 EXACTLY (verified 0 diff,
  3fast + SF4f 10 ms). Labels say "= PPR fig panel 3" (was the stale "= primary ROM").
- **Coverslip-average tiles**: title = `PPR <primary, gate-either pooled, from cs_rows> (trace <this avg
  trace, decay-sub = panel 3>) Q <panel 2>`. Across-ROI SEM band REMOVED (ROIs are pseudoreplicates);
  grey region (c#) means show the spread. Shading now = blue/red peak-search windows (same as montage and
  decay tiles); Q-window end = grey dotted line. ISIs split at <=120 / >120 ms (none dropped any more).
- **Montage**: shows gate-either ROIs, tag [A1|A2|A1+A2]; title = A1 and decay-sub A2 (what the primary
  pools), NOT a per-ROI ratio (A2-only ROIs have tiny A1 -> ratios of 14, NaN). Orange diamond =
  decay-sub A2 (as in decay tiles); p1-lock marker changed d -> ^ so diamond means one thing.
- Coverslip-average + decay tiles use `make_time_axis` (same x unit as the montage). PPR fig white bg.
- NB MATLAB TeX has no `\diamond` (it errors); use the word.
- **Q window single-sourced**: `tracenorm/twoap_charge_window.m` returns the integrated frames
  [stim1, stim2 + round(auc_post_ms / frame)]; used by twoap_analyzer (per-ROI q_pair), trace_ppr
  (avg-trace Q) AND the grey band in the coverslip-average + montage figures, so the shading IS the
  integrated window (verified: band ends at stim2 + 500.9 ms = 100 frames at 5.009 ms). Q values unchanged.
  Dropped the analyzer's hidden 250 ms default for auc_post_ms (now required from config).

### DONE 2026-10-15: sync / async event binning (3fast only)
- `tracenorm/twoap_sync_async_events.m`: every event in a window = local max >= upper (5 sigma),
  Schmitt run, RISE >= `events.min_prominence_sigma` (5) sd above the trough BEFORE it (left trough
  only; findpeaks' two-sided prominence killed real pulse-1 peaks at 30 ms ISI, where the window is
  cut by stim2 -> P(sync) P1 0.04-0.09 instead of ~0.2). `allow_end` for the cut pulse-1 window.
- Windows: P1 = [stim1, min(stim2-1, stim1+100 ms)], P2 = [stim2, stim2+100 ms]; spont = pre-stim1.
  Config `config.twoAP.events` (sync_cutoff_ms, async_window_ms, min_prominence_sigma, max_frame_ms=6:
  10 ms SF4f is skipped).
- **Cutoff on PEAK latency, 30 ms** (user, 2026-10-15). Data: 3fast first-event peaks spread 5-30 ms
  (mode 15 ms, rise + frame integration), then a flat async tail ~0.001-0.002 ev/ROI/frame to 100 ms.
  A 10 ms peak cutoff would call ~60% of sync responses async. Detector robust: 3/4/5 sigma prominence
  give 0.33/0.31/0.30 events/ROI. Spont rate ~0.001-0.002 ev/s/ROI (negligible correction).
- Batch: `compute_sync_async` (all QC-passed ROIs; the detector is the gate) -> per coverslip x ISI
  P(sync) P1/P2, sync PPF, async events/ROI and /s (spont-subtracted), async PPF (equal windows only),
  async share of P2 events, mean event rise. Figure `plot_handles.sync_async`, summary block, CSV
  `twoap_sync_async_per_coverslip.csv`, per-ROI CSV gains n_sync1/n_async1/n_sync2/n_async2/n_spont.
  Test: `tests/test_twoap_sync_async_events.m`.
- First look (3fast NGR, pilot N): sync PPF 1.8 @30 ms -> ~1 @120 -> 0.5-0.7 @>=250 (same shape as
  amplitude PPR); async after P2 rises steeply at short ISI (WT 3.1/s @30 vs ~0.2/s @>=250; async share
  of P2 events 0.37 -> 0.05-0.09). WT and R213W overlap.
- Caveats: at short ISI the P2 async bin also holds late pulse-1 async events and possibly pulse-2
  sync peaks delayed by riding on pulse-1's decay. async PPF is a ratio of a few events per coverslip
  (noisy, 0.2-8). P(sync) is per ROI, so ROI segmentation (boutons per ROI) affects it.
- Example-trace figure `plot_handles.sync_async_ex` (2026-10-16): per group, seeded random ROIs with
  async events + sync-only ROIs (`visualization.montage.twoAP_sync_async_rois = [8 4]`), markers =
  the analyzer's own events, gap between pulse windows cut out at long ISI.
- Short-ISI P2 async (2026-10-16): the largest P2 event per ROI lands past 30 ms in 13-34% of ROIs at
  ISI <= 120 ms vs 2-9% at >= 250 ms. **Resolved by the avg-trace check**: decay-subtracted P2 avg
  trace (all QC ROIs) peaks at 14-20 ms at EVERY ISI (= P1, 16-17 ms; 50% rise 3-7 ms), but its
  centroid moves 29-35 -> 39-48 ms and the area share after 30 ms 0.40-0.52 -> 0.56-0.69 at ISI 30-90.
  So the sync peak does not shift; short ISI adds LATE charge = real extra late release (residual Ca),
  not a delayed sync peak.

### CHANGED 2026-10-16: sync / async binned on event ONSET <= 10 ms (user)
- 3fast kinetics from these data (pulse-1 avg trace, ISI >= 250, per coverslip): time to peak 15-20 ms,
  10-90% rise 12-14 ms, half-decay 16-20 ms, single-exp tau 18-22 ms, FWHM 28-32 ms (upper bounds:
  5 ms frame integration + release jitter). Peak-aligned isolated single events are NOT usable (aligning
  thresholded noisy events on their max makes a 1-frame spike; only 11 isolated spont events).
- Onset = first frame of the final climb above `events.onset_frac` (0.2) of trough -> peak, trough =
  lowest frame in `events.onset_trough_ms` (25) before the peak, not past the previous event. The first
  try (nearest local min before the peak) put the onset at the last noise dip of a jagged rise
  (often = peak frame) -> biased async. `sync_cutoff_ms` = 10 (2 frames). Event lists carry `onset_ms`;
  histograms are onset latency; example figure marks the onset frame (open square + faint line).
- Effect (3fast NGR): P(sync) P1 ~0.15-0.27 (was 0.2-0.4 with peak <= 30); async share of P2 events
  0.18-0.28 at ISI >= 250, 0.37-0.61 at ISI <= 120; sync PPF shape unchanged (1.8-1.9 @30 -> 0.6-0.8
  @>=250). Onset histogram: high at 0-5 ms, a 10-25 ms shoulder, then a ~10x drop near 30 ms to a flat
  tail -> most "async" events are the 10-25 ms shoulder. At 5 ms frames onset is +-1 frame, so events
  rising ~10 ms land on either side of the cutoff.

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
