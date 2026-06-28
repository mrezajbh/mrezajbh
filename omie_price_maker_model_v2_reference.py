#!/usr/bin/env python3
"""
OMIE price-maker stress test, v2
================================

This script is a research model for one question:

    When an optimization paper says "buy electricity in the cheapest hours",
    what happens if the order itself is large enough to move the day-ahead
    auction curve?

Compared with v1, this version adds a more economic counterfactual. It does not
only look at residual sell offers above the official price. It also uses the
published aggregate buy curve. For a new inelastic purchase q, it asks how far the
published aggregate net-supply curve must move from the observed clearing price.

Core empirical object
---------------------
For period t and price p, build from OMIE offered curves:

    S_t(p) = total sell volume offered at prices <= p
    D_t(p) = total buy volume offered at prices >= p
    N_t(p) = S_t(p) - D_t(p)          # aggregate net supply

Let p0_t be the official OMIE marginal day-ahead price. Because the full official
European clearing algorithm includes market coupling, complex orders, network
constraints and local rules, directly re-clearing simple aggregate curves will
not always reproduce p0_t. So the default estimator is deliberately anchored:

    P_t(q) = min{p >= p0_t : N_t(p) - N_t(p0_t) >= q}.

That says: to clear q additional MWh of price-inelastic demand, the market must
move up the observed aggregate net-supply curve until q more MWh of net supply is
available. This uses both extra sellers entering and pre-existing buyers leaving.

The script also reports two robustness mechanisms:

    residual_supply  : v1-style unmatched sell volume above p0_t.
    raw_reclear_delta: raw supply-demand reclearing delta, then anchored to p0_t.

Caveat
------
This is not a reimplementation of EUPHEMIA. It is a transparent local
counterfactual using public OMIE curves. It should be described as a price-impact
stress test or partial-equilibrium auction-curve perturbation.

Data files
----------
- curva_pbc_YYYYMMDD.1        aggregate day-ahead supply/demand curves
- marginalpdbc_YYYYMMDD.1     official day-ahead marginal prices

Example
-------
python omie_price_maker_model_v2.py \
  --start 2024-01-01 --end 2024-01-07 \
  --country ES \
  --daily-energy 1000 --max-power 250 --step 25 \
  --cache-dir raw_cache --no-download \
  --out out_v2

Dependencies
------------
pip install pandas numpy matplotlib requests
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import zipfile
from dataclasses import dataclass
from datetime import date, datetime, timedelta
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
except Exception:  # pragma: no cover
    plt = None

try:
    import requests
except Exception:  # pragma: no cover
    requests = None


OMIE_DOWNLOAD = "https://www.omie.es/en/file-download?filename={filename}&parents={parent}"
MECHANISMS = ("anchored_net_supply", "residual_supply", "raw_reclear_delta")


@dataclass(frozen=True)
class Config:
    start: date
    end: date
    country: str = "ES"
    order_sizes: Tuple[float, ...] = (10, 50, 100, 250, 500, 1000)
    daily_energy: float = 1000.0
    max_power: float = 250.0
    step: float = 25.0
    bid_price: float = 3000.0
    mechanisms: Tuple[str, ...] = ("anchored_net_supply", "residual_supply")
    out: Path = Path("out_omie_v2")
    cache_dir: Optional[Path] = None
    no_download: bool = False
    overwrite: bool = False


@dataclass
class PeriodCurve:
    date: pd.Timestamp
    period: int
    zone: str
    p0: float
    offered: pd.DataFrame
    matched: pd.DataFrame
    bid_price: float = 3000.0

    def __post_init__(self) -> None:
        off = self.offered
        self._sell = off[off.side == "V"].groupby("price")["mwh"].sum().sort_index()
        self._buy = off[off.side == "C"].groupby("price")["mwh"].sum().sort_index()
        self._sell_prices = self._sell.index.to_numpy(dtype=float)
        self._sell_cum = self._sell.cumsum().to_numpy(dtype=float)
        self._buy_prices = self._buy.index.to_numpy(dtype=float)
        self._buy_vol = self._buy.to_numpy(dtype=float)
        self._buy_cum_less = np.concatenate(([0.0], np.cumsum(self._buy_vol)))
        self._buy_total = float(self._buy_vol.sum())
        # Grid on which prices can move. Include p0 and the buyer bid price even
        # if they are not existing offer-price ticks.
        grid = set(float(x) for x in off["price"].dropna().to_numpy(dtype=float))
        grid.add(float(self.p0))
        grid.add(float(self.bid_price))
        self._grid = np.array(sorted(grid), dtype=float)
        self._net0 = self.net_supply(self.p0)
        self._raw0 = self.raw_clearing_price(0.0)
        self._residual_curve = self._build_residual_sell_curve()

    def supply(self, price: float) -> float:
        idx = int(np.searchsorted(self._sell_prices, price, side="right"))
        if idx <= 0:
            return 0.0
        return float(self._sell_cum[idx - 1])

    def demand(self, price: float) -> float:
        # Sum buy offers with bid price >= price.
        idx_less = int(np.searchsorted(self._buy_prices, price, side="left"))
        return float(self._buy_total - self._buy_cum_less[idx_less])

    def net_supply(self, price: float) -> float:
        return self.supply(price) - self.demand(price)

    def incremental_net_supply(self, price: float) -> float:
        return self.net_supply(price) - self._net0

    def _build_residual_sell_curve(self) -> pd.Series:
        offered_sell = self.offered[self.offered.side == "V"].groupby("price")["mwh"].sum()
        matched_sell = self.matched[self.matched.side == "V"].groupby("price")["mwh"].sum()
        residual = (offered_sell - matched_sell.reindex(offered_sell.index, fill_value=0.0)).clip(lower=0.0)
        residual = residual[residual.index >= self.p0 - 1e-9].sort_index()
        return residual[residual > 1e-9]

    def price_anchored_net_supply(self, q_mwh: float) -> Tuple[float, float]:
        """Default model: min p >= p0 such that N(p)-N(p0) >= q.

        Returns (price, unserved_mwh). unserved_mwh is nonzero only if the buyer's
        bid-price cap is reached before enough net supply appears.
        """
        q = float(q_mwh)
        if q <= 1e-12:
            return float(self.p0), 0.0
        candidates = self._grid[(self._grid >= self.p0 - 1e-9) & (self._grid <= self.bid_price + 1e-9)]
        if candidates.size == 0:
            return float(self.p0), q
        target = self._net0 + q
        # Loop is fine: per period grids are a few thousand points.
        for p in candidates:
            if self.net_supply(float(p)) + 1e-9 >= target:
                return float(max(p, self.p0)), 0.0
        max_depth = max(0.0, self.net_supply(float(candidates[-1])) - self._net0)
        return float(candidates[-1]), max(0.0, q - max_depth)

    def price_residual_supply(self, q_mwh: float) -> Tuple[float, float]:
        """v1 robustness model: unmatched offered sell curve above p0."""
        q = float(q_mwh)
        if q <= 1e-12:
            return float(self.p0), 0.0
        residual = self._residual_curve
        if residual.empty:
            return float(self.p0), q
        cum = residual.cumsum().to_numpy(dtype=float)
        prices = residual.index.to_numpy(dtype=float)
        # Respect bid cap.
        valid = prices <= self.bid_price + 1e-9
        if not valid.any():
            return float(self.p0), q
        cum_valid = cum[valid]
        prices_valid = prices[valid]
        idx = int(np.searchsorted(cum_valid, q, side="left"))
        if idx >= len(prices_valid):
            return float(prices_valid[-1]), max(0.0, q - float(cum_valid[-1]))
        return float(max(prices_valid[idx], self.p0)), 0.0

    def raw_clearing_price(self, q_mwh: float) -> float:
        """Raw simple-curve reclearing: min p with S(p) >= D(p)+q.

        This usually does not reproduce official p0 because the official clearing
        algorithm contains details not visible in aggregate simple curves. We use
        it only as a robustness delta.
        """
        q = float(q_mwh)
        candidates = self._grid[self._grid <= self.bid_price + 1e-9]
        if candidates.size == 0:
            return float(self.p0)
        for p in candidates:
            if self.supply(float(p)) + 1e-9 >= self.demand(float(p)) + q:
                return float(p)
        return float(candidates[-1])

    def price_raw_reclear_delta(self, q_mwh: float) -> Tuple[float, float]:
        q = float(q_mwh)
        if q <= 1e-12:
            return float(self.p0), 0.0
        raw_q = self.raw_clearing_price(q)
        p = max(float(self.p0), float(self.p0 + raw_q - self._raw0))
        if p > self.bid_price:
            return float(self.bid_price), np.nan
        return float(p), 0.0

    def price(self, q_mwh: float, mechanism: str) -> Tuple[float, float]:
        if mechanism == "anchored_net_supply":
            return self.price_anchored_net_supply(q_mwh)
        if mechanism == "residual_supply":
            return self.price_residual_supply(q_mwh)
        if mechanism == "raw_reclear_delta":
            return self.price_raw_reclear_delta(q_mwh)
        raise KeyError(f"Unknown mechanism: {mechanism}")

    def depth_to_move(self, delta_eur_mwh: float, mechanism: str) -> float:
        """How many MWh can be added before price reaches p0+delta?

        Larger means deeper/safer. This is a diagnostic, not an optimization input.
        """
        threshold = min(self.p0 + float(delta_eur_mwh), self.bid_price)
        if mechanism == "anchored_net_supply":
            return max(0.0, self.net_supply(threshold) - self._net0)
        if mechanism == "residual_supply":
            residual = self._residual_curve
            if residual.empty:
                return 0.0
            return float(residual[residual.index <= threshold + 1e-9].sum())
        if mechanism == "raw_reclear_delta":
            # Approximate using the anchored net-supply diagnostic; raw-reclear is
            # for robustness, not the main depth measure.
            return max(0.0, self.net_supply(threshold) - self._net0)
        raise KeyError(mechanism)


# ---------------------------------------------------------------------------
# IO utilities
# ---------------------------------------------------------------------------

def daterange(start: date, end: date) -> Iterable[date]:
    d = start
    while d <= end:
        yield d
        d += timedelta(days=1)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def download_omie_file(parent: str, filename: str, dest: Path, overwrite: bool = False) -> Path:
    if dest.exists() and dest.stat().st_size > 0 and not overwrite:
        return dest
    if requests is None:
        raise RuntimeError("requests is not installed. Install it or run with --no-download and --cache-dir.")
    url = OMIE_DOWNLOAD.format(filename=filename, parent=parent)
    r = requests.get(url, timeout=90)
    r.raise_for_status()
    if len(r.content) < 50:
        raise RuntimeError(f"Downloaded file is suspiciously small: {url}")
    dest.write_bytes(r.content)
    return dest


def fetch_day_files(d: date, cache_dir: Path, no_download: bool, overwrite: bool) -> Tuple[Path, Path]:
    ds = d.strftime("%Y%m%d")
    curve_name = f"curva_pbc_{ds}.1"
    price_name = f"marginalpdbc_{ds}.1"
    curve_path = cache_dir / curve_name
    price_path = cache_dir / price_name
    if no_download:
        missing = [str(p) for p in (curve_path, price_path) if not p.exists()]
        if missing:
            raise FileNotFoundError("Missing OMIE files in cache: " + ", ".join(missing))
    else:
        download_omie_file("curva_pbc", curve_name, curve_path, overwrite=overwrite)
        download_omie_file("marginalpdbc", price_name, price_path, overwrite=overwrite)
    return curve_path, price_path


def period_to_int(x) -> int:
    if isinstance(x, str):
        xs = x.strip()
        m = re.match(r"H(\d+)Q(\d+)", xs)
        if m:
            h, q = int(m.group(1)), int(m.group(2))
            return (h - 1) * 4 + q
    return int(float(x))


def read_curve_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=";",
        encoding="latin1",
        skiprows=2,
        decimal=",",
        thousands=".",
        engine="python",
    )
    df.columns = [str(c).strip() for c in df.columns]
    df = df.loc[:, ~df.columns.str.contains(r"^Unnamed")]

    rename_map = {}
    for c in df.columns:
        cs = str(c).strip()
        if cs in ("Hora", "Periodo"):
            rename_map[c] = "period"
        elif cs == "Fecha":
            rename_map[c] = "date"
        elif cs in ("Pais", "País"):
            rename_map[c] = "zone"
        elif cs == "Tipo Oferta":
            rename_map[c] = "side"
        elif cs in ("Energía Compra/Venta", "Energia Compra/Venta", "Potencia Compra/Venta"):
            rename_map[c] = "mwh"
        elif cs == "Precio Compra/Venta":
            rename_map[c] = "price"
        elif cs.startswith("Ofertada"):
            rename_map[c] = "flag"
    df = df.rename(columns=rename_map)

    required = ["period", "date", "zone", "side", "mwh", "price", "flag"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns {missing} in {path}. Found: {df.columns.tolist()}")

    df = df.dropna(subset=required).copy()
    df["period"] = df["period"].map(period_to_int).astype(int)
    df["date"] = pd.to_datetime(df["date"], format="%d/%m/%Y", errors="coerce")
    df["zone"] = df["zone"].astype(str).str.strip()
    df["side"] = df["side"].astype(str).str.strip()  # C=Compra/buy, V=Venta/sell
    df["flag"] = df["flag"].astype(str).str.strip()  # O=offered, C=matched
    df["mwh"] = pd.to_numeric(df["mwh"], errors="coerce")
    df["price"] = pd.to_numeric(df["price"], errors="coerce")
    df = df.dropna(subset=["date", "mwh", "price"])
    df = df[df["mwh"] > 0].copy()
    return df[required]


def read_price_file(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(
        path,
        sep=";",
        skiprows=1,
        header=None,
        names=["year", "month", "day", "period", "price_portugal", "price_spain", "_"],
        engine="python",
    )
    raw = raw.dropna(subset=["period"]).copy()
    raw = raw[raw["year"].astype(str).str.match(r"^\d{4}$", na=False)]
    raw["date"] = pd.to_datetime(
        dict(year=raw["year"].astype(int), month=raw["month"].astype(int), day=raw["day"].astype(int))
    )
    raw["period"] = raw["period"].map(period_to_int).astype(int)
    raw["price_spain"] = pd.to_numeric(raw["price_spain"], errors="coerce")
    raw["price_portugal"] = pd.to_numeric(raw["price_portugal"], errors="coerce")
    return raw[["date", "period", "price_spain", "price_portugal"]]


def choose_zone_for_country(g: pd.DataFrame, country: str = "ES") -> str:
    country = country.upper()
    zones = set(g["zone"].unique())
    if country in zones:
        return country
    if "MI" in zones:
        return "MI"
    if len(zones) == 1:
        return next(iter(zones))
    return str(g["zone"].value_counts().index[0])


def official_price_col(country: str) -> str:
    return "price_spain" if country.upper() == "ES" else "price_portugal"


# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------

def build_period_curves(curves: pd.DataFrame, prices: pd.DataFrame, cfg: Config) -> Dict[Tuple[pd.Timestamp, int], PeriodCurve]:
    price_col = official_price_col(cfg.country)
    price_map = prices.set_index(["date", "period"])[price_col]
    out: Dict[Tuple[pd.Timestamp, int], PeriodCurve] = {}
    for (dt, period), gp in curves.groupby(["date", "period"], sort=True):
        if (dt, period) not in price_map.index:
            continue
        zone = choose_zone_for_country(gp, cfg.country)
        gz = gp[gp.zone == zone]
        p0 = float(price_map.loc[(dt, period)])
        offered = gz[gz.flag == "O"].copy()
        matched = gz[gz.flag == "C"].copy()
        if offered.empty:
            continue
        out[(pd.Timestamp(dt), int(period))] = PeriodCurve(
            date=pd.Timestamp(dt), period=int(period), zone=zone, p0=p0,
            offered=offered, matched=matched, bid_price=cfg.bid_price,
        )
    return out


def build_impact_tables(period_curves: Mapping[Tuple[pd.Timestamp, int], PeriodCurve], cfg: Config) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    wide_rows = []
    for key in sorted(period_curves):
        pc = period_curves[key]
        base = {
            "date": pc.date,
            "period": pc.period,
            "zone": pc.zone,
            "p0": pc.p0,
            "offered_buy_mwh": float(pc.offered[pc.offered.side == "C"].mwh.sum()),
            "offered_sell_mwh": float(pc.offered[pc.offered.side == "V"].mwh.sum()),
            "matched_buy_mwh": float(pc.matched[pc.matched.side == "C"].mwh.sum()),
            "matched_sell_mwh": float(pc.matched[pc.matched.side == "V"].mwh.sum()),
            "net_supply_at_p0": float(pc.net_supply(pc.p0)),
            "raw_simple_curve_price_q0": float(pc._raw0),
            "raw_minus_official_eur_mwh": float(pc._raw0 - pc.p0),
        }
        wide = dict(base)
        for mech in cfg.mechanisms:
            for q in cfg.order_sizes:
                p_after, unserved = pc.price(q, mech)
                row = dict(base)
                row.update(
                    {
                        "mechanism": mech,
                        "order_mwh": float(q),
                        "p_after": float(p_after),
                        "dp_eur_mwh": float(p_after - pc.p0),
                        "unserved_mwh_at_bid_cap": float(unserved) if pd.notna(unserved) else np.nan,
                    }
                )
                rows.append(row)
                wide[f"{mech}_p_{int(q)}"] = float(p_after)
                wide[f"{mech}_dp_{int(q)}"] = float(p_after - pc.p0)
                wide[f"{mech}_unserved_{int(q)}"] = float(unserved) if pd.notna(unserved) else np.nan
        wide_rows.append(wide)
    return pd.DataFrame(rows), pd.DataFrame(wide_rows)


def build_depth_table(period_curves: Mapping[Tuple[pd.Timestamp, int], PeriodCurve], cfg: Config) -> pd.DataFrame:
    rows = []
    for key in sorted(period_curves):
        pc = period_curves[key]
        for mech in cfg.mechanisms:
            row = {"date": pc.date, "period": pc.period, "zone": pc.zone, "p0": pc.p0, "mechanism": mech}
            for d in (1, 2, 5, 10, 20):
                row[f"depth_before_plus_{d}_eur_mwh"] = pc.depth_to_move(d, mech)
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Scheduling models
# ---------------------------------------------------------------------------

def price_taker_schedule(day_prices: pd.DataFrame, energy_mwh: float, max_power_mwh: float) -> Dict[int, float]:
    remaining = float(energy_mwh)
    sched = {int(p): 0.0 for p in day_prices["period"]}
    for _, row in day_prices.sort_values(["p0", "period"]).iterrows():
        take = min(float(max_power_mwh), remaining)
        if take <= 1e-9:
            break
        sched[int(row["period"])] = take
        remaining -= take
    if remaining > 1e-6:
        raise ValueError("Not enough periods/max_power to allocate daily energy.")
    return sched


def impact_aware_dp(
    period_curves_for_day: Mapping[int, PeriodCurve],
    periods: Sequence[int],
    mechanism: str,
    energy_mwh: float,
    max_power_mwh: float,
    step_mwh: float,
) -> Dict[int, float]:
    """Discrete DP minimizing invoice cost sum_t q_t * P_t(q_t)."""
    K_float = energy_mwh / step_mwh
    M_float = max_power_mwh / step_mwh
    K = int(round(K_float))
    M = int(round(M_float))
    if not math.isclose(K_float, K, abs_tol=1e-9):
        raise ValueError("daily-energy must be an integer multiple of --step.")
    if not math.isclose(M_float, M, abs_tol=1e-9):
        raise ValueError("max-power must be an integer multiple of --step.")

    n = len(periods)
    dp = np.full((n + 1, K + 1), np.inf)
    prev = np.full((n + 1, K + 1), -1, dtype=int)
    dp[0, 0] = 0.0
    for j, period in enumerate(periods, start=1):
        pc = period_curves_for_day[int(period)]
        costs = []
        for m in range(M + 1):
            q = m * step_mwh
            p_after, unserved = pc.price(q, mechanism)
            if unserved and unserved > 1e-7:
                costs.append(np.inf)
            else:
                costs.append(q * p_after)
        for k in range(K + 1):
            best = np.inf
            best_m = -1
            max_m = min(M, k)
            for m in range(max_m + 1):
                val = dp[j - 1, k - m] + costs[m]
                if val < best:
                    best = val
                    best_m = m
            dp[j, k] = best
            prev[j, k] = best_m

    if not np.isfinite(dp[n, K]):
        raise ValueError(f"No feasible impact-aware schedule under mechanism={mechanism}")

    sched = {int(p): 0.0 for p in periods}
    k = K
    for j in range(n, 0, -1):
        m = int(prev[j, k])
        p = int(periods[j - 1])
        sched[p] = m * step_mwh
        k -= m
    return sched


def evaluate_schedule(
    period_curves_for_day: Mapping[int, PeriodCurve],
    sched: Mapping[int, float],
    mechanism: str,
    policy: str,
) -> pd.DataFrame:
    rows = []
    for period, q in sorted(sched.items()):
        pc = period_curves_for_day[int(period)]
        p_after, unserved = pc.price(float(q), mechanism)
        rows.append(
            {
                "date": pc.date,
                "period": int(period),
                "zone": pc.zone,
                "mechanism": mechanism,
                "policy": policy,
                "q_mwh": float(q),
                "p0": pc.p0,
                "p_after": float(p_after),
                "dp_eur_mwh": float(p_after - pc.p0),
                "unserved_mwh_at_bid_cap": float(unserved) if pd.notna(unserved) else np.nan,
                "naive_cost_eur": float(q) * pc.p0,
                "actual_invoice_cost_eur": float(q) * float(p_after),
            }
        )
    return pd.DataFrame(rows)


def run_scheduling(period_curves: Mapping[Tuple[pd.Timestamp, int], PeriodCurve], cfg: Config) -> Tuple[pd.DataFrame, pd.DataFrame]:
    by_day: Dict[pd.Timestamp, Dict[int, PeriodCurve]] = {}
    for (dt, period), pc in period_curves.items():
        by_day.setdefault(dt, {})[period] = pc

    summary_rows = []
    detail_frames = []
    for dt in sorted(by_day):
        day = by_day[dt]
        periods = sorted(day)
        day_prices = pd.DataFrame({"period": periods, "p0": [day[p].p0 for p in periods]})
        pt_sched = price_taker_schedule(day_prices, cfg.daily_energy, cfg.max_power)

        for mech in cfg.mechanisms:
            pt_eval = evaluate_schedule(day, pt_sched, mech, "price_taking")
            ia_sched = impact_aware_dp(day, periods, mech, cfg.daily_energy, cfg.max_power, cfg.step)
            ia_eval = evaluate_schedule(day, ia_sched, mech, "impact_aware")
            detail_frames.extend([pt_eval, ia_eval])

            naive = float(pt_eval["naive_cost_eur"].sum())
            pt_actual = float(pt_eval["actual_invoice_cost_eur"].sum())
            ia_naive = float(ia_eval["naive_cost_eur"].sum())
            ia_actual = float(ia_eval["actual_invoice_cost_eur"].sum())
            summary_rows.append(
                {
                    "date": dt,
                    "mechanism": mech,
                    "daily_energy_mwh": cfg.daily_energy,
                    "max_power_mwh_per_period": cfg.max_power,
                    "step_mwh": cfg.step,
                    "price_taking_naive_cost_eur": naive,
                    "price_taking_actual_invoice_cost_eur": pt_actual,
                    "impact_aware_naive_cost_eur": ia_naive,
                    "impact_aware_actual_invoice_cost_eur": ia_actual,
                    "naive_understatement_pct": (pt_actual / naive - 1.0) * 100 if naive else np.nan,
                    "impact_aware_savings_vs_price_taking_pct": (1.0 - ia_actual / pt_actual) * 100 if pt_actual else np.nan,
                    "price_taking_periods": ",".join(str(k) for k, v in pt_sched.items() if v > 0),
                    "impact_aware_periods": ",".join(str(k) for k, v in ia_sched.items() if v > 0),
                }
            )
    detail = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame()
    return pd.DataFrame(summary_rows), detail


# ---------------------------------------------------------------------------
# Figures and reporting
# ---------------------------------------------------------------------------

def mechanism_label(x: str) -> str:
    return {
        "anchored_net_supply": "anchored net supply",
        "residual_supply": "residual supply",
        "raw_reclear_delta": "raw reclear delta",
    }.get(x, x)


def plot_outputs(impact_long: pd.DataFrame, depth: pd.DataFrame, sched_summary: pd.DataFrame, sched_detail: pd.DataFrame, cfg: Config) -> None:
    if plt is None:
        print("matplotlib is not installed; skipping figures", file=sys.stderr)
        return
    fig_dir = cfg.out / "figures"
    ensure_dir(fig_dir)

    # Figure 1: price-impact quantiles by mechanism.
    qdf = (
        impact_long.groupby(["mechanism", "order_mwh"])["dp_eur_mwh"]
        .quantile([0.50, 0.75, 0.95])
        .rename("dp")
        .reset_index()
        .pivot_table(index=["mechanism", "order_mwh"], columns="level_2", values="dp")
        .reset_index()
    )
    for mech in cfg.mechanisms:
        sub = qdf[qdf.mechanism == mech].sort_values("order_mwh")
        if sub.empty:
            continue
        plt.figure(figsize=(7.4, 4.9))
        plt.plot(sub["order_mwh"], sub[0.50], marker="o", label="median")
        plt.plot(sub["order_mwh"], sub[0.75], marker="o", label="75th percentile")
        plt.plot(sub["order_mwh"], sub[0.95], marker="o", label="95th percentile")
        plt.xlabel("Added inelastic demand per period (MWh)")
        plt.ylabel("Price impact (€/MWh)")
        plt.title(f"Price impact: {mechanism_label(mech)}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(fig_dir / f"fig1_price_impact_quantiles_{mech}.png", dpi=220)
        plt.close()

    # Figure 2: side-by-side total schedule cost.
    for mech in cfg.mechanisms:
        sub = sched_summary[sched_summary.mechanism == mech]
        if sub.empty:
            continue
        totals = pd.Series(
            {
                "price-taking\nnaive": sub["price_taking_naive_cost_eur"].sum(),
                "price-taking\nactual": sub["price_taking_actual_invoice_cost_eur"].sum(),
                "impact-aware\nactual": sub["impact_aware_actual_invoice_cost_eur"].sum(),
            }
        )
        plt.figure(figsize=(6.6, 4.9))
        totals.plot(kind="bar")
        plt.ylabel("Total buyer invoice cost (€)")
        plt.title(f"Flexible-load scheduling under {mechanism_label(mech)}")
        plt.xticks(rotation=0)
        plt.tight_layout()
        plt.savefig(fig_dir / f"fig2_schedule_costs_{mech}.png", dpi=220)
        plt.close()

    # Figure 3: one dramatic counterfactual curve for the default mechanism.
    default_mech = cfg.mechanisms[0]
    wide = impact_long[impact_long.mechanism == default_mech]
    if not wide.empty:
        idx = wide[wide.order_mwh == max(cfg.order_sizes)]["dp_eur_mwh"].idxmax()
        chosen = wide.loc[idx]
        sub = wide[(wide.date == chosen.date) & (wide.period == chosen.period)].sort_values("order_mwh")
        plt.figure(figsize=(7.2, 4.9))
        plt.plot(sub["order_mwh"], sub["p_after"], marker="o")
        plt.axhline(float(chosen["p0"]), linestyle="--", linewidth=1)
        plt.xlabel("Added demand in the period (MWh)")
        plt.ylabel("Counterfactual clearing price (€/MWh)")
        plt.title(f"A cheap hour has finite depth: {pd.Timestamp(chosen.date).date()}, period {int(chosen.period)}")
        plt.tight_layout()
        plt.savefig(fig_dir / f"fig3_selected_counterfactual_curve_{default_mech}.png", dpi=220)
        plt.close()

    # Figure 4: daily understatement by mechanism.
    if not sched_summary.empty:
        plt.figure(figsize=(8.0, 4.9))
        for mech in cfg.mechanisms:
            sub = sched_summary[sched_summary.mechanism == mech].sort_values("date")
            if not sub.empty:
                plt.plot(pd.to_datetime(sub["date"]), sub["naive_understatement_pct"], marker="o", label=mechanism_label(mech))
        plt.ylabel("Price-taking cost understatement (%)")
        plt.xlabel("Date")
        plt.title("How much the cheapest-hour model underprices its own order")
        plt.legend()
        plt.tight_layout()
        plt.savefig(fig_dir / "fig4_daily_understatement.png", dpi=220)
        plt.close()

    # Figure 5: depth distribution before a +5 €/MWh move.
    if not depth.empty:
        plt.figure(figsize=(7.2, 4.9))
        data = []
        labels = []
        for mech in cfg.mechanisms:
            sub = depth[depth.mechanism == mech]["depth_before_plus_5_eur_mwh"].dropna()
            if not sub.empty:
                data.append(sub.to_numpy())
                labels.append(mechanism_label(mech))
        if data:
            
            try:
                plt.boxplot(data, tick_labels=labels, showfliers=False)
            except TypeError:  # older matplotlib
                plt.boxplot(data, labels=labels, showfliers=False)
            plt.ylabel("MWh before price rises by 5 €/MWh")
            plt.title("Market depth is not constant across periods")
            plt.xticks(rotation=15)
            plt.tight_layout()
            plt.savefig(fig_dir / "fig5_depth_before_plus_5.png", dpi=220)
            plt.close()


def write_model_note(cfg: Config, impact_long: pd.DataFrame, depth: pd.DataFrame, sched_summary: pd.DataFrame) -> Path:
    note = cfg.out / "MODEL_NOTE.md"
    lines: List[str] = []
    lines.append("# OMIE price-maker model v2\n")
    lines.append("## Critical question\n")
    lines.append("> If the optimizer selects the cheapest electricity hours, are those hours still cheap after the optimizer submits a large MWh-scale order?\n")
    lines.append("## Main model\n")
    lines.append("For each period, the script constructs offered supply S(p), offered demand D(p), and net supply N(p)=S(p)-D(p). ")
    lines.append("The default counterfactual is anchored to the official OMIE price p0 and sets P(q)=min{p>=p0: N(p)-N(p0)>=q}. ")
    lines.append("This is intentionally not a full EUPHEMIA rerun; it is a transparent local perturbation of the public aggregate auction curves.\n")
    lines.append("## Why this is better than v1\n")
    lines.append("v1 used only unmatched sell volume above p0. v2 uses both sides of the book: as price rises, extra sellers enter and some historical buyers leave. ")
    lines.append("That makes the price-impact estimate more disciplined and usually less alarmist, while still testing whether the price-taking assumption survives scale.\n")
    if not impact_long.empty:
        lines.append("## Price-impact quantiles\n")
        q = (
            impact_long.groupby(["mechanism", "order_mwh"])["dp_eur_mwh"]
            .quantile([0.5, 0.75, 0.95])
            .unstack()
            .reset_index()
        )
        lines.append(q.to_markdown(index=False, floatfmt=".3f"))
        lines.append("\n")
    if not sched_summary.empty:
        lines.append("## Scheduling summary\n")
        cols = [
            "mechanism", "price_taking_naive_cost_eur", "price_taking_actual_invoice_cost_eur",
            "impact_aware_actual_invoice_cost_eur", "naive_understatement_pct",
            "impact_aware_savings_vs_price_taking_pct",
        ]
        totals = sched_summary.groupby("mechanism")[cols[1:]].sum(numeric_only=True).reset_index()
        # Percent columns need weighted recomputation.
        recomputed = []
        for mech, g in sched_summary.groupby("mechanism"):
            naive = g["price_taking_naive_cost_eur"].sum()
            pt = g["price_taking_actual_invoice_cost_eur"].sum()
            ia = g["impact_aware_actual_invoice_cost_eur"].sum()
            recomputed.append({
                "mechanism": mech,
                "price_taking_naive_cost_eur": naive,
                "price_taking_actual_invoice_cost_eur": pt,
                "impact_aware_actual_invoice_cost_eur": ia,
                "naive_understatement_pct": (pt/naive - 1)*100 if naive else np.nan,
                "impact_aware_savings_vs_price_taking_pct": (1 - ia/pt)*100 if pt else np.nan,
            })
        total_df = pd.DataFrame(recomputed)
        lines.append(total_df.to_markdown(index=False, floatfmt=".3f"))
        lines.append("\n")
    if not depth.empty:
        lines.append("## Depth diagnostic\n")
        lines.append("The depth table asks: how many extra MWh can be added before the period price rises by 1, 2, 5, 10, or 20 €/MWh? A shallow 'cheap' period is dangerous for price-taking optimization.\n")
    note.write_text("\n".join(lines), encoding="utf-8")
    return note




def fmt_num(x: float) -> str:
    return str(int(round(float(x)))) if abs(float(x) - round(float(x))) < 1e-9 else str(float(x)).replace(".", "p")

def write_additional_diagnostics(impact_long: pd.DataFrame, impact_wide: pd.DataFrame, cfg: Config) -> None:
    """Write validation, rank-reversal, critical-mass, and comparison diagnostics."""
    # 1) Model-comparison quantiles.
    quant = (
        impact_long.groupby(["mechanism", "order_mwh"])["dp_eur_mwh"]
        .agg(
            mean_dp="mean",
            median_dp="median",
            p75_dp=lambda x: x.quantile(0.75),
            p95_dp=lambda x: x.quantile(0.95),
            max_dp="max",
            share_positive_dp=lambda x: (x > 0).mean(),
            share_ge_5eur=lambda x: (x >= 5).mean(),
            share_ge_10eur=lambda x: (x >= 10).mean(),
        )
        .reset_index()
    )
    quant.to_csv(cfg.out / "model_comparison_quantiles.csv", index=False)

    # 2) Validation of raw simple curve crossing. This justifies anchoring.
    if "raw_minus_official_eur_mwh" in impact_wide:
        err = impact_wide["raw_minus_official_eur_mwh"].dropna()
        validation = pd.DataFrame([
            {
                "n_periods": int(len(err)),
                "raw_simple_crossing_mae_eur_mwh": float(err.abs().mean()) if len(err) else np.nan,
                "raw_simple_crossing_rmse_eur_mwh": float(np.sqrt((err ** 2).mean())) if len(err) else np.nan,
                "raw_simple_crossing_median_abs_error_eur_mwh": float(err.abs().median()) if len(err) else np.nan,
                "raw_simple_crossing_p90_abs_error_eur_mwh": float(err.abs().quantile(0.90)) if len(err) else np.nan,
                "raw_simple_crossing_bias_eur_mwh": float(err.mean()) if len(err) else np.nan,
                "share_abs_error_le_1eur": float((err.abs() <= 1).mean()) if len(err) else np.nan,
                "share_abs_error_le_5eur": float((err.abs() <= 5).mean()) if len(err) else np.nan,
            }
        ])
    else:
        validation = pd.DataFrame()
    validation.to_csv(cfg.out / "validation_summary.csv", index=False)

    # 3) Rank reversal: does the published-price cheapest period remain cheapest after q?
    rank_rows = []
    for (mech, q), gmq in impact_long.groupby(["mechanism", "order_mwh"], sort=True):
        for dt, g in gmq.groupby("date", sort=True):
            base = g.sort_values(["p0", "period"]).iloc[0]
            best = g.sort_values(["p_after", "period"]).iloc[0]
            rank_rows.append(
                {
                    "date": dt,
                    "mechanism": mech,
                    "order_mwh": float(q),
                    "published_cheapest_period": int(base["period"]),
                    "published_cheapest_p0": float(base["p0"]),
                    "published_cheapest_after_price": float(base["p_after"]),
                    "impact_adjusted_cheapest_period": int(best["period"]),
                    "impact_adjusted_cheapest_price": float(best["p_after"]),
                    "cheapest_period_changes": int(int(base["period"]) != int(best["period"])),
                    "avoidable_price_error_eur_mwh": max(0.0, float(base["p_after"]) - float(best["p_after"])),
                }
            )
    rank = pd.DataFrame(rank_rows)
    rank.to_csv(cfg.out / "rank_reversal_diagnostics.csv", index=False)
    if not rank.empty:
        rank_summary = (
            rank.groupby(["mechanism", "order_mwh"])
            .agg(
                reversal_share_pct=("cheapest_period_changes", lambda x: 100.0 * x.mean()),
                median_avoidable_price_error=("avoidable_price_error_eur_mwh", "median"),
                p75_avoidable_price_error=("avoidable_price_error_eur_mwh", lambda x: x.quantile(0.75)),
                max_avoidable_price_error=("avoidable_price_error_eur_mwh", "max"),
            )
            .reset_index()
        )
    else:
        rank_summary = pd.DataFrame()
    rank_summary.to_csv(cfg.out / "rank_reversal_summary.csv", index=False)

    # 4) Critical mass: smallest tested q that creates a threshold price move.
    thresholds = [1.0, 5.0, 10.0, 20.0]
    crit_rows = []
    for (dt, period, mech), g in impact_long.groupby(["date", "period", "mechanism"], sort=True):
        row = {"date": dt, "period": int(period), "mechanism": mech, "p0": float(g["p0"].iloc[0])}
        for th in thresholds:
            hit = g[g["dp_eur_mwh"] >= th].sort_values("order_mwh")
            row[f"critical_mwh_for_plus_{fmt_num(th)}eur"] = float(hit["order_mwh"].iloc[0]) if not hit.empty else np.nan
        crit_rows.append(row)
    crit = pd.DataFrame(crit_rows)
    crit.to_csv(cfg.out / "critical_mass_diagnostics.csv", index=False)

    crit_summary_rows = []
    for mech, g in crit.groupby("mechanism") if not crit.empty else []:
        for th in thresholds:
            col = f"critical_mwh_for_plus_{fmt_num(th)}eur"
            values = g[col].dropna()
            crit_summary_rows.append(
                {
                    "mechanism": mech,
                    "threshold_eur_per_mwh": float(th),
                    "share_reached_within_tested_orders_pct": float(100.0 * g[col].notna().mean()),
                    "median_critical_mwh_conditional": float(values.median()) if len(values) else np.nan,
                    "p25_critical_mwh_conditional": float(values.quantile(0.25)) if len(values) else np.nan,
                    "p75_critical_mwh_conditional": float(values.quantile(0.75)) if len(values) else np.nan,
                }
            )
    crit_summary = pd.DataFrame(crit_summary_rows)
    crit_summary.to_csv(cfg.out / "critical_mass_summary.csv", index=False)

    # Extra figures.
    if plt is not None:
        fig_dir = cfg.out / "figures"
        ensure_dir(fig_dir)
        if not quant.empty:
            plt.figure(figsize=(7.4, 4.9))
            for mech, g in quant.groupby("mechanism"):
                g = g.sort_values("order_mwh")
                plt.plot(g["order_mwh"], g["median_dp"], marker="o", label=mechanism_label(mech))
            plt.xlabel("Added demand per period (MWh)")
            plt.ylabel("Median price impact (€/MWh)")
            plt.title("Median impact across counterfactual mechanisms")
            plt.legend()
            plt.tight_layout()
            plt.savefig(fig_dir / "fig6_median_impact_model_comparison.png", dpi=220)
            plt.close()

        if not rank_summary.empty:
            plt.figure(figsize=(7.4, 4.9))
            for mech, g in rank_summary.groupby("mechanism"):
                g = g.sort_values("order_mwh")
                plt.plot(g["order_mwh"], g["reversal_share_pct"], marker="o", label=mechanism_label(mech))
            plt.xlabel("Added demand in the published-cheapest period (MWh)")
            plt.ylabel("Days where cheapest period changes (%)")
            plt.title("Does the published-price cheapest period remain cheapest?")
            plt.legend()
            plt.tight_layout()
            plt.savefig(fig_dir / "fig7_rank_reversal_share.png", dpi=220)
            plt.close()

        main = "anchored_net_supply"
        cm = crit_summary[crit_summary["mechanism"] == main].sort_values("threshold_eur_per_mwh") if not crit_summary.empty else pd.DataFrame()
        if not cm.empty:
            plt.figure(figsize=(7.4, 4.9))
            plt.plot(cm["threshold_eur_per_mwh"], cm["share_reached_within_tested_orders_pct"], marker="o")
            plt.xlabel("Price-impact threshold (€/MWh)")
            plt.ylabel("Periods reaching threshold within tested orders (%)")
            plt.title("Critical mass: how often scale moves the price")
            plt.tight_layout()
            plt.savefig(fig_dir / "fig8_critical_mass_reach_anchored_net_supply.png", dpi=220)
            plt.close()

    # Append compact diagnostics to the model note.
    note = cfg.out / "MODEL_NOTE.md"
    if note.exists():
        text = note.read_text(encoding="utf-8")
        if not validation.empty:
            text += "\n\n## Validation of the raw simple crossing\n\n"
            text += (
                "The raw simple supply-demand crossing is reported only as a diagnostic. "
                f"For this run, its MAE against the official OMIE price is "
                f"{validation['raw_simple_crossing_mae_eur_mwh'].iloc[0]:.2f} €/MWh, "
                "which justifies anchoring the v2 perturbation at the official price rather than pretending to re-run EUPHEMIA.\n"
            )
        if not rank_summary.empty:
            text += "\n## Rank-reversal diagnostic\n\n"
            text += "The clean optimization question is whether the published-price cheapest period remains cheapest after the order. Summary:\n\n"
            text += rank_summary.to_markdown(index=False, floatfmt=".3f")
            text += "\n"
        if not crit_summary.empty:
            text += "\n## Critical-mass diagnostic\n\n"
            text += "Critical mass is the smallest tested order size that creates at least a given price movement. Summary:\n\n"
            text += crit_summary.to_markdown(index=False, floatfmt=".3f")
            text += "\n"
        note.write_text(text, encoding="utf-8")

def print_console_summary(impact_long: pd.DataFrame, sched_summary: pd.DataFrame, cfg: Config) -> None:
    print("\nPrice-impact quantiles by mechanism (€/MWh):")
    for mech in cfg.mechanisms:
        sub = impact_long[impact_long.mechanism == mech]
        if sub.empty:
            continue
        table = sub.pivot_table(index="order_mwh", values="dp_eur_mwh", aggfunc=["median", lambda x: x.quantile(0.75), lambda x: x.quantile(0.95), "max"])
        table.columns = ["median", "p75", "p95", "max"]
        print(f"\n[{mech}]")
        print(table.to_string(float_format=lambda x: f"{x:,.3f}"))

    if not sched_summary.empty:
        print("\nScheduling totals:")
        rows = []
        for mech, g in sched_summary.groupby("mechanism"):
            naive = g["price_taking_naive_cost_eur"].sum()
            pt = g["price_taking_actual_invoice_cost_eur"].sum()
            ia = g["impact_aware_actual_invoice_cost_eur"].sum()
            rows.append({
                "mechanism": mech,
                "price_taking_naive": naive,
                "price_taking_actual": pt,
                "impact_aware_actual": ia,
                "understatement_%": (pt/naive - 1)*100 if naive else np.nan,
                "savings_%": (1 - ia/pt)*100 if pt else np.nan,
            })
        print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:,.3f}"))


def make_zip(out_dir: Path, zip_path: Optional[Path] = None) -> Path:
    zip_path = zip_path or out_dir.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in sorted(out_dir.rglob("*")):
            if p.is_file():
                zf.write(p, arcname=p.relative_to(out_dir.parent))
    return zip_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[List[str]] = None) -> Config:
    parser = argparse.ArgumentParser(description="OMIE price-maker electricity optimization stress test, v2")
    parser.add_argument("--start", required=True, help="start date, YYYY-MM-DD")
    parser.add_argument("--end", required=True, help="end date, YYYY-MM-DD")
    parser.add_argument("--country", default="ES", choices=["ES", "PT"], help="buyer price zone")
    parser.add_argument("--order-sizes", default="10,50,100,250,500,1000", help="comma-separated MWh shocks")
    parser.add_argument("--daily-energy", type=float, default=1000.0, help="MWh/day in flexible-load example")
    parser.add_argument("--max-power", type=float, default=250.0, help="max MWh per period")
    parser.add_argument("--step", type=float, default=25.0, help="MWh discretization step for dynamic program")
    parser.add_argument("--bid-price", type=float, default=3000.0, help="buyer bid price cap in €/MWh")
    parser.add_argument(
        "--mechanisms",
        default="anchored_net_supply,residual_supply",
        help="comma-separated mechanisms: anchored_net_supply,residual_supply,raw_reclear_delta,all",
    )
    parser.add_argument("--out", default="out_omie_v2", help="output directory")
    parser.add_argument("--cache-dir", default=None, help="directory containing/downloading OMIE files; default OUT/cache")
    parser.add_argument("--no-download", action="store_true", help="use files already present in cache-dir")
    parser.add_argument("--overwrite", action="store_true", help="re-download files if present")
    args = parser.parse_args(argv)

    mechs_raw = [x.strip() for x in args.mechanisms.split(",") if x.strip()]
    if "all" in mechs_raw:
        mechanisms = MECHANISMS
    else:
        unknown = [m for m in mechs_raw if m not in MECHANISMS]
        if unknown:
            raise ValueError(f"Unknown mechanisms {unknown}. Allowed: {MECHANISMS} or all")
        mechanisms = tuple(mechs_raw)
    order_sizes = tuple(float(x.strip()) for x in args.order_sizes.split(",") if x.strip())
    out = Path(args.out)
    cache_dir = Path(args.cache_dir) if args.cache_dir else out / "cache"
    return Config(
        start=datetime.strptime(args.start, "%Y-%m-%d").date(),
        end=datetime.strptime(args.end, "%Y-%m-%d").date(),
        country=args.country.upper(),
        order_sizes=order_sizes,
        daily_energy=args.daily_energy,
        max_power=args.max_power,
        step=args.step,
        bid_price=args.bid_price,
        mechanisms=mechanisms,
        out=out,
        cache_dir=cache_dir,
        no_download=args.no_download,
        overwrite=args.overwrite,
    )


def main(argv: Optional[List[str]] = None) -> int:
    cfg = parse_args(argv)
    ensure_dir(cfg.out)
    assert cfg.cache_dir is not None
    ensure_dir(cfg.cache_dir)

    curve_frames = []
    price_frames = []
    for d in daterange(cfg.start, cfg.end):
        curve_path, price_path = fetch_day_files(d, cfg.cache_dir, cfg.no_download, cfg.overwrite)
        curve_frames.append(read_curve_file(curve_path))
        price_frames.append(read_price_file(price_path))

    curves = pd.concat(curve_frames, ignore_index=True)
    prices = pd.concat(price_frames, ignore_index=True)
    period_curves = build_period_curves(curves, prices, cfg)
    if not period_curves:
        raise RuntimeError("No period curves built. Check country, dates and files.")

    impact_long, impact_wide = build_impact_tables(period_curves, cfg)
    depth = build_depth_table(period_curves, cfg)
    sched_summary, sched_detail = run_scheduling(period_curves, cfg)

    impact_long.to_csv(cfg.out / "impact_results_long.csv", index=False)
    impact_wide.to_csv(cfg.out / "impact_results_wide.csv", index=False)
    depth.to_csv(cfg.out / "depth_metrics.csv", index=False)
    sched_summary.to_csv(cfg.out / "schedule_summary.csv", index=False)
    sched_detail.to_csv(cfg.out / "schedule_detail.csv", index=False)

    plot_outputs(impact_long, depth, sched_summary, sched_detail, cfg)
    write_model_note(cfg, impact_long, depth, sched_summary)
    write_additional_diagnostics(impact_long, impact_wide, cfg)
    zip_path = make_zip(cfg.out)

    print_console_summary(impact_long, sched_summary, cfg)
    print(f"\nSaved outputs in: {cfg.out.resolve()}")
    print(f"Saved bundle: {zip_path.resolve()}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
