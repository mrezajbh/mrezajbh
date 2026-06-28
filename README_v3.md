# OMIE price-maker stress test — v3

`omie_price_maker_model_v3.py` is a backward-compatible upgrade of the v2 model used in
the manuscript *"When Cheap Electricity Is No Longer Cheap"*. It targets the four
limitations stated in the paper while leaving the published results intact.

## Governing principle: v3 is a strict superset of v2

Run with default arguments, v3 reproduces v2's output **byte-for-byte**. This was
verified two ways:

1. The price-responsive clearing routine, fed a vertical (inelastic) bid, returns the
   *exact* same price and unserved quantity as v2's `price_anchored_net_supply` across
   25 randomised curves × 7 order sizes (`tests_v3.py::test_inelastic_bid_equals_anchored`).
2. A full v2-vs-v3 pipeline run on identical inputs produces 11 identical core CSVs
   (`impact_results_long`, `impact_results_wide`, `depth_metrics`, `schedule_summary`,
   `schedule_detail`, `model_comparison_quantiles`, `rank_reversal_summary`,
   `critical_mass_summary`, `validation_summary`, and both diagnostics files).

The bootstrap also reproduces the paper's headline statistics on the real pilot results:
median impact at 1000 MWh = **3.775 €/MWh**, p95 = **17.524 €/MWh**.

Every addition uses **only your own OMIE data** plus standard mathematical/statistical
methods. No elasticity, coefficient, or empirical magnitude is hard-coded into the model.

## What each limitation gets, and what stays out of scope

### 1. Not EUPHEMIA (anchoring)
v3 does **not** re-implement EUPHEMIA — full market-coupling re-clear with complex
orders, welfare maximisation, and network constraints remains a separate, much larger
project, exactly as the paper states. What v3 adds is making the anchoring assumption
*auditable* rather than asserted:
- `curve_validation.csv` — checks the constructed net-supply curve `N(p)` for
  monotonicity on its own grid and reports the local slope `dN/dp` at the clearing price
  (the "depth" in MWh per €/MWh).
- `anchoring_sensitivity.csv` — recomputes each impact estimate with the anchor shifted
  by **± that period's own raw-crossing residual** `|raw − official|`, and reports how
  far the estimate moves. The perturbation size comes from the data, not from outside.
- `impact_bracket.csv` / `impact_bracket_summary.csv` — reorganises the three mechanisms
  into an explicit lower / central / upper bound, in the spirit of the supply-curve-bounds
  approach in the robust self-scheduling literature already cited in the paper.

**Honest limit:** these quantify and bound the anchoring error; they do not eliminate it.

### 2. Sample size
- **Fault-tolerant batch loading** (`load_days_robust`): a missing file, an unparseable
  curve, or a missing official price is logged in `data_coverage.csv` and skipped, so a
  full-year or multi-year run completes instead of aborting on a single bad day.
- **Both Iberian zones in one run**: `--zones ES,PT` writes per-zone subfolders.
- **Bootstrap confidence intervals** (`confidence_intervals.csv`): every headline number
  (median / mean / p95 impact, share ≥ 5 €/MWh, rank-reversal share, naive understatement,
  impact-aware savings) is reported with a percentile CI obtained by resampling **whole
  days** with replacement. Because periods within a day are serially correlated, the day
  is the natural resampling unit (a cluster bootstrap). The point estimate is unchanged;
  only the uncertainty is now stated explicitly.

**Honest limit:** the code cannot fabricate data you don't have. On the one-week pilot the
CIs are wide (e.g. median impact at 1000 MWh: 3.775 with 95% CI ≈ [2.81, 5.16]) — which is
precisely the proof-of-concept caveat made quantitative. Point the loader at a full year of
files and the same machinery delivers the econometric-grade run.

### 3. Behavioural assumption (inelastic demand)
The added order can now be a **price-responsive demand bid** instead of a vertical line.
`PeriodCurve.price_for_demand_bid` solves `P = min{p ≥ p0 : N(p) − N(p0) ≥ b(p)}` for any
non-increasing bid `b(p)`. Supported shapes (all opt-in via `--demand-bid`):
- `inelastic` (default) — vertical bid; reproduces v2.
- `piecewise` — a user-supplied `(price, quantity)` curve via `--demand-bid-file`.
- `elastic` — constant-elasticity `q(p) = q_ref·(p/p_ref)^(−ε)`, with `ε` and `p_ref`
  **supplied by you** (`--demand-elasticity`, `--demand-ref-price`).

Results go to `impact_price_responsive.csv` (cleared vs withdrawn MWh) and never touch the
core inelastic tables.

**Honest limit:** the bid *shape* is an input grounded in your data or in the literature;
the model implements the mechanism but does not invent an elasticity. Block / exclusive-group
orders are not yet implemented — their acceptance economics require an explicit surplus
criterion and are the natural next extension.

### 4. Hourly vs 15-minute resolution
Period counts are inferred from the data (never assumed to be 24), so a 96-period day flows
through unchanged. `--period-minutes 15` documents the market time unit and is recorded in
`resolution_report.csv`. The `HxQx` quarter-hour period labels are already parsed.

**Honest limit:** the exact 15-minute OMIE file schema (column layout / sentinels) should be
confirmed against a real downloaded sample before production use — none was available to test
against here, so the reader is written to be tolerant rather than hard-coded to a guessed
format.

## Usage

```bash
# Reproduces v2 exactly (default inelastic, single zone):
python omie_price_maker_model_v3.py \
  --start 2024-01-01 --end 2024-01-07 --country ES \
  --daily-energy 1000 --max-power 250 --step 25 \
  --cache-dir raw_cache --no-download --out out_v3

# Full year, both Iberian zones, fault-tolerant, 2000 bootstrap reps:
python omie_price_maker_model_v3.py \
  --start 2024-01-01 --end 2024-12-31 --zones ES,PT \
  --bootstrap-reps 2000 --ci-level 0.95 \
  --cache-dir raw_cache --no-download --out out_v3_year

# Price-responsive buyer (elasticity supplied by you):
python omie_price_maker_model_v3.py \
  --start 2024-01-01 --end 2024-01-07 --country ES \
  --demand-bid elastic --demand-elasticity 0.3 --demand-ref-price 60 \
  --cache-dir raw_cache --no-download --out out_v3_elastic

# Quarter-hour resolution (when 15-minute inputs exist):
python omie_price_maker_model_v3.py \
  --start 2025-06-01 --end 2025-06-07 --country ES --period-minutes 15 \
  --cache-dir raw_cache_15min --no-download --out out_v3_qh
```

New CLI flags: `--zones`, `--period-minutes`, `--bootstrap-reps`, `--ci-level`, `--seed`,
`--demand-bid`, `--demand-elasticity`, `--demand-ref-price`, `--demand-bid-file`.

## New output files (additions only)

| File | Limitation | Contents |
|---|---|---|
| `data_coverage.csv` | 2 | per-day load status / skip reason |
| `confidence_intervals.csv` | 2 | day-cluster bootstrap CIs for every headline statistic |
| `resolution_report.csv` | 4 | inferred periods-per-day and configured MTU |
| `curve_validation.csv` | 1 | `N(p)` monotonicity + local slope at `p0` |
| `anchoring_sensitivity.csv` | 1 | impact under ± own-residual anchor shifts |
| `impact_bracket.csv`, `impact_bracket_summary.csv` | 1 | lower/central/upper mechanism bounds |
| `impact_price_responsive.csv` | 3 | cleared vs withdrawn MWh under an elastic/piecewise bid (only if `--demand-bid` ≠ inelastic) |

## Tests

```bash
python tests_v3.py
```
Seven checks: inelastic ≡ anchored (exact), elastic monotonicity, piecewise validation,
bootstrap reproducing the paper numbers on the real CSV, full synthetic pipeline incl. a
skipped missing day, quarter-hour handling, and the price-responsive output path.

## Files

- `omie_price_maker_model_v3.py` — the upgraded model.
- `omie_price_maker_model_v2_reference.py` — the original v2, unchanged, for diffing.
- `tests_v3.py` — the test suite.
- `README_v3.md` — this file.

## Note on the archive format

This bundle is delivered as a **.zip**. A true RAR archive needs the proprietary `rar`
encoder, which was not available in the build environment. A `.zip` opens identically in
WinRAR, 7-Zip, and the OS archiver; rename or re-compress to `.rar` on your machine if a
RAR container is specifically required.
