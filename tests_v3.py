#!/usr/bin/env python3
"""Tests for omie_price_maker_model_v3.

These verify three things:
  1. Backward compatibility: the price-responsive generalisation reduces EXACTLY to
     the v2 anchored inelastic estimator, so the paper's numbers are preserved.
  2. The new statistics reproduce known quantities on the REAL pilot result CSV.
  3. The full pipeline (parse -> curves -> impact -> schedule -> bootstrap ->
     diagnostics -> zip) runs end to end on synthetic OMIE-format files, including a
     missing-day skip and a quarter-hour (96-period) day.

No synthetic number here is part of any reported result; these fixtures only exercise
code paths. Run:  python tests_v3.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

import omie_price_maker_model_v3 as m


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

def make_synthetic_curve(p0: float = 50.0, seed: int = 0, bid_price: float = 3000.0) -> m.PeriodCurve:
    """Build a PeriodCurve directly from a controlled offer ladder (no file IO)."""
    rng = np.random.default_rng(seed)
    sell_prices = np.array([10, 25, 40, 50, 65, 80, 100, 130, 170, 220], dtype=float)
    sell_mwh = rng.uniform(80, 300, size=sell_prices.size)
    buy_prices = np.array([260, 180, 120, 90, 65, 45, 25], dtype=float)
    buy_mwh = rng.uniform(80, 300, size=buy_prices.size)

    rows = []
    for pr, q in zip(sell_prices, sell_mwh):
        rows.append({"date": pd.Timestamp("2024-01-01"), "period": 1, "zone": "ES", "side": "V", "price": float(pr), "mwh": float(q), "flag": "O"})
    for pr, q in zip(buy_prices, buy_mwh):
        rows.append({"date": pd.Timestamp("2024-01-01"), "period": 1, "zone": "ES", "side": "C", "price": float(pr), "mwh": float(q), "flag": "O"})
    offered = pd.DataFrame(rows)
    # Matched: sells at/below p0 and buys at/above p0 (only needed so residual curve exists).
    matched = offered.copy()
    matched = matched[((matched.side == "V") & (matched.price <= p0)) | ((matched.side == "C") & (matched.price >= p0))].copy()
    matched["flag"] = "C"
    return m.PeriodCurve(date=pd.Timestamp("2024-01-01"), period=1, zone="ES", p0=float(p0), offered=offered, matched=matched, bid_price=bid_price)


def write_synthetic_omie_day(cache: Path, d: str, periods: int = 24, p0: float = 50.0) -> None:
    """Write curva_pbc and marginalpdbc files in the formats the readers expect."""
    ds = d.replace("-", "")
    yyyy, mm, dd = d.split("-")
    # --- curve file: 2 ignored lines, then header, comma-decimal, dd/mm/YYYY ---
    header = "Hora;Fecha;Pais;Unidad;Tipo Oferta;Energia Compra/Venta;Precio Compra/Venta;Ofertada (O)/Casada (C)"
    lines = ["OMIE - synthetic curve fixture", "(second ignored line)", header]
    date_es = f"{int(dd):02d}/{int(mm):02d}/{yyyy}"
    sell_prices = [10, 25, 40, 50, 65, 80, 100, 130, 170, 220]
    buy_prices = [260, 180, 120, 90, 65, 45, 25]
    for per in range(1, periods + 1):
        for pr in sell_prices:
            lines.append(f"{per};{date_es};ES;U1;V;{150.0:.1f};{float(pr):.1f};O".replace(".", ","))
            if pr <= p0:
                lines.append(f"{per};{date_es};ES;U1;V;{150.0:.1f};{float(pr):.1f};C".replace(".", ","))
        for pr in buy_prices:
            lines.append(f"{per};{date_es};ES;U1;C;{150.0:.1f};{float(pr):.1f};O".replace(".", ","))
            if pr >= p0:
                lines.append(f"{per};{date_es};ES;U1;C;{150.0:.1f};{float(pr):.1f};C".replace(".", ","))
    (cache / f"curva_pbc_{ds}.1").write_text("\n".join(lines) + "\n", encoding="latin1")

    # --- price file: 1 ignored header line, dot-decimal, year;month;day;period;pt;es; ---
    plines = ["MARGINALPDBC;synthetic"]
    for per in range(1, periods + 1):
        plines.append(f"{yyyy};{int(mm)};{int(dd)};{per};{p0 - 2.0:.2f};{p0:.2f};")
    plines.append("*")
    (cache / f"marginalpdbc_{ds}.1").write_text("\n".join(plines) + "\n", encoding="latin1")


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

def test_inelastic_bid_equals_anchored():
    """CRITICAL: price-responsive clearing of a vertical bid == v2 anchored estimator."""
    fails = 0
    for seed in range(25):
        p0 = float(40 + seed)
        pc = make_synthetic_curve(p0=p0, seed=seed)
        for q in (10, 50, 100, 250, 500, 1000, 5000):
            p_old, u_old = pc.price_anchored_net_supply(q)
            p_new, u_new = pc.price_for_demand_bid(m.InelasticBid(cap=pc.bid_price, q=float(q)))
            if abs(p_old - p_new) > 1e-9 or abs(u_old - u_new) > 1e-6:
                fails += 1
                print(f"  MISMATCH seed={seed} q={q}: old=({p_old},{u_old}) new=({p_new},{u_new})")
    assert fails == 0, f"{fails} inelastic-vs-anchored mismatches"
    print("PASS test_inelastic_bid_equals_anchored (vertical bid reduces to v2 exactly)")


def test_elastic_bid_withdraws():
    """A price-responsive buyer clears no more than nominal and withdraws as price rises."""
    pc = make_synthetic_curve(p0=50.0, seed=3)
    q_ref = 1000.0
    p_inel, _ = pc.price_for_demand_bid(m.InelasticBid(cap=pc.bid_price, q=q_ref))
    last_price = None
    for e in (0.0, 0.2, 0.5, 1.0):
        bid = m.ConstantElasticityBid(cap=pc.bid_price, q_ref=q_ref, p_ref=pc.p0, elasticity=e)
        p_after, unserved = pc.price_for_demand_bid(bid)
        cleared = max(0.0, bid.quantity_at(p_after) - unserved)
        assert cleared <= q_ref + 1e-6, "cleared exceeds nominal reference quantity"
        assert p_after >= pc.p0 - 1e-9, "price below p0"
        if last_price is not None:
            assert p_after <= last_price + 1e-9, "more elastic buyer should clear at <= price"
        last_price = p_after
    # zero elasticity must coincide with the inelastic clearing price.
    bid0 = m.ConstantElasticityBid(cap=pc.bid_price, q_ref=q_ref, p_ref=pc.p0, elasticity=0.0)
    p0e, _ = pc.price_for_demand_bid(bid0)
    assert abs(p0e - p_inel) < 1e-9, "elasticity=0 should match inelastic"
    print("PASS test_elastic_bid_withdraws (price-responsive demand behaves monotonically)")


def test_piecewise_validation_and_interp():
    ok = m.PiecewiseLinearBid(cap=3000.0, points=((10.0, 1000.0), (50.0, 600.0), (100.0, 200.0)))
    assert abs(ok.quantity_at(30.0) - 800.0) < 1e-9  # midpoint of first segment
    assert ok.quantity_at(5.0) == 1000.0 and ok.quantity_at(200.0) == 200.0  # held outside range
    for bad in (((10.0, 100.0), (50.0, 200.0)),  # increasing quantity
                ((50.0, 100.0), (10.0, 200.0))):  # unsorted prices
        try:
            m.PiecewiseLinearBid(cap=3000.0, points=bad)
            assert False, "should have rejected invalid piecewise bid"
        except ValueError:
            pass
    print("PASS test_piecewise_validation_and_interp")


def test_bootstrap_on_real_results():
    """Bootstrap reproduces the paper's point estimates and brackets them."""
    real = Path("/home/claude/bundle/price_maker_paper_final/results/impact_results_long.csv")
    if not real.exists():
        print("SKIP test_bootstrap_on_real_results (real CSV not present)")
        return
    d = pd.read_csv(real)
    sub = d[(d.mechanism == "anchored_net_supply") & (d.order_mwh == 1000.0)]
    med = float(sub["dp_eur_mwh"].median())
    p95 = float(sub["dp_eur_mwh"].quantile(0.95))
    # Known paper values (compared at the 3-decimal precision the paper reports).
    assert round(med, 3) == 3.775, f"median {med} != 3.775"
    assert round(p95, 3) == 17.524, f"p95 {p95} (rounds to {round(p95,3)}) != 17.524"
    rng = np.random.default_rng(0)
    point, lo, hi, n = m.day_cluster_bootstrap(sub, lambda x: x["dp_eur_mwh"].median(), B=500, level=0.95, rng=rng)
    assert abs(point - med) < 1e-9, "bootstrap point estimate must equal the plain estimate"
    assert lo <= point + 1e-9 <= hi + 1e-9 or lo <= point <= hi, "point not inside CI"
    assert lo <= hi, "CI low above high"
    assert n == 7, f"expected 7 days, got {n}"
    print(f"PASS test_bootstrap_on_real_results (median={point:.3f}, 95% CI=[{lo:.3f},{hi:.3f}], n_days={n})")


def test_full_pipeline_synthetic():
    """End-to-end run on synthetic files, plus missing-day skip."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        cache = tmp / "cache"
        cache.mkdir()
        write_synthetic_omie_day(cache, "2024-01-01", periods=24, p0=50.0)
        write_synthetic_omie_day(cache, "2024-01-02", periods=24, p0=55.0)
        # Day 2024-01-03 intentionally absent -> must be skipped, not crash.
        out = tmp / "out"
        rc = m.main([
            "--start", "2024-01-01", "--end", "2024-01-03",
            "--country", "ES",
            "--order-sizes", "50,250,1000",
            "--daily-energy", "1000", "--max-power", "250", "--step", "250",
            "--bootstrap-reps", "50",
            "--cache-dir", str(cache), "--no-download",
            "--out", str(out),
        ])
        assert rc == 0
        core = ["impact_results_long.csv", "schedule_summary.csv", "depth_metrics.csv"]
        added = ["confidence_intervals.csv", "data_coverage.csv", "curve_validation.csv",
                 "anchoring_sensitivity.csv", "impact_bracket.csv", "resolution_report.csv"]
        for f in core + added:
            assert (out / f).exists(), f"missing output {f}"
        cov = pd.read_csv(out / "data_coverage.csv")
        assert (cov.status == "loaded").sum() == 2, "expected 2 loaded days"
        assert (cov.status != "loaded").sum() == 1, "expected 1 skipped day"
        brk = pd.read_csv(out / "impact_bracket.csv")
        assert (brk["dp_bracket_low"] <= brk["dp_central"] + 1e-9).all(), "bracket low > central"
        assert (brk["dp_central"] <= brk["dp_bracket_high"] + 1e-9).all(), "central > bracket high"
        print("PASS test_full_pipeline_synthetic (ran end to end; 1 missing day skipped; bracket ordered)")


def test_quarter_hour_resolution():
    """A 96-period synthetic day flows through; resolution is inferred, not assumed."""
    assert m.period_to_int("H3Q2") == 10  # 15-min label parsing
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        cache = tmp / "cache"
        cache.mkdir()
        write_synthetic_omie_day(cache, "2025-06-01", periods=96, p0=60.0)
        out = tmp / "out_qh"
        rc = m.main([
            "--start", "2025-06-01", "--end", "2025-06-01",
            "--country", "ES", "--order-sizes", "100,1000",
            "--daily-energy", "1000", "--max-power", "500", "--step", "500",
            "--period-minutes", "15", "--bootstrap-reps", "20",
            "--cache-dir", str(cache), "--no-download", "--out", str(out),
        ])
        assert rc == 0
        res = pd.read_csv(out / "resolution_report.csv")
        assert int(res["periods_in_day"].iloc[0]) == 96, "did not infer 96 periods"
        assert int(res["configured_period_minutes"].iloc[0]) == 15
        print("PASS test_quarter_hour_resolution (96-period day handled; period_to_int OK)")


def test_price_responsive_output():
    """Opt-in elastic analysis writes its own file and clears <= nominal."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        cache = tmp / "cache"
        cache.mkdir()
        write_synthetic_omie_day(cache, "2024-01-01", periods=24, p0=50.0)
        out = tmp / "out_el"
        rc = m.main([
            "--start", "2024-01-01", "--end", "2024-01-01",
            "--country", "ES", "--order-sizes", "250,1000",
            "--daily-energy", "1000", "--max-power", "500", "--step", "500",
            "--demand-bid", "elastic", "--demand-elasticity", "0.3", "--demand-ref-price", "50",
            "--bootstrap-reps", "20",
            "--cache-dir", str(cache), "--no-download", "--out", str(out),
        ])
        assert rc == 0
        pr = pd.read_csv(out / "impact_price_responsive.csv")
        assert (pr["cleared_mwh"] <= pr["nominal_order_mwh"] + 1e-6).all()
        assert (pr["dp_eur_mwh"] >= -1e-9).all()
        # The core inelastic table must still exist and be untouched by the elastic run.
        assert (out / "impact_results_long.csv").exists()
        print("PASS test_price_responsive_output (elastic file written; core table intact)")


def main() -> int:
    tests = [
        test_inelastic_bid_equals_anchored,
        test_elastic_bid_withdraws,
        test_piecewise_validation_and_interp,
        test_bootstrap_on_real_results,
        test_full_pipeline_synthetic,
        test_quarter_hour_resolution,
        test_price_responsive_output,
    ]
    failed = 0
    for t in tests:
        try:
            t()
        except AssertionError as exc:
            failed += 1
            print(f"FAIL {t.__name__}: {exc}")
        except Exception as exc:  # noqa
            failed += 1
            import traceback
            print(f"ERROR {t.__name__}: {exc}")
            traceback.print_exc()
    print("\n" + ("ALL TESTS PASSED" if failed == 0 else f"{failed} TEST(S) FAILED"))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
