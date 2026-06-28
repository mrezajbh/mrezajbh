#!/usr/bin/env python3
"""
OMIE price-maker stress test, v3
================================

This script is a research model for one question:

    When an optimization paper says "buy electricity in the cheapest hours",
    what happens if the order itself is large enough to move the day-ahead
    auction curve?

Compared with v1, this version adds a more economic counterfactual. It does not
only look at residual sell offers above the official price. It also uses the
published aggregate buy curve. For a new inelastic purchase q, it asks how far the
published aggregate net-supply curve must move from the observed clearing price.

What is new in v3
-----------------
v3 is a strict superset of v2. With default arguments it reproduces the v2 core
results exactly (same impact tables, same scheduling numbers). All additions are
either opt-in flags or extra diagnostic files that do not touch the core
estimates. Every addition uses only the user's own OMIE data plus standard
mathematical or statistical methods; no external coefficient, elasticity, or
empirical magnitude is hard-coded into the model.

The additions target the four documented limitations of the v2 paper:

1. Anchoring vs full EUPHEMIA. v3 does not re-implement EUPHEMIA (that remains a
   separate, much larger project). Instead it makes the anchoring assumption
   auditable: it validates the constructed net-supply curve for monotonicity,
   reports the local curve slope at the clearing price, and quantifies how much
   each impact estimate moves when the anchor price is perturbed by the period's
   own raw-crossing residual (anchoring-band sensitivity). The three mechanisms
   are reorganised into an explicit lower/central/upper bracket, in the spirit of
   the supply-curve-bounds approach used in the robust self-scheduling literature
   that the paper already cites.

2. Sample size. The data loader is now fault-tolerant: days with missing or
   unparseable files are logged and skipped rather than aborting the whole run,
   so a full-year or multi-year batch completes. A data-coverage report is
   written. Both Iberian price zones can be run in one invocation. Crucially, all
   headline statistics are reported with bootstrap confidence intervals computed
   by resampling whole days (a cluster bootstrap that respects within-day serial
   correlation), so however long the sample is, the uncertainty is stated
   honestly rather than implied.

3. Behavioural assumption. The added order can now be a price-responsive demand
   bid instead of a vertical inelastic line. The clearing price of a monotone
   non-increasing demand bid against the incremental net-supply curve is solved
   exactly; the inelastic case is recovered as the special case of a vertical bid
   at the bid-price cap. The bid shape (piecewise-linear breakpoints, or a
   constant-elasticity form whose elasticity the user supplies) is an input, not
   an assumption baked into the code.

4. Market time unit. Period counts are inferred from the data rather than fixed
   at 24, and an explicit --period-minutes flag documents the resolution, so the
   same diagnostic runs at quarter-hour resolution when 15-minute inputs are
   provided. (The exact 15-minute OMIE file schema should be confirmed against a
   real sample before production use.)

Core empirical object

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
# Reproduces v2 exactly (inelastic demand, single zone):
python omie_price_maker_model_v3.py \
  --start 2024-01-01 --end 2024-01-07 \
  --country ES \
  --daily-energy 1000 --max-power 250 --step 25 \
  --cache-dir raw_cache --no-download \
  --out out_v3

# Full-year, both Iberian zones, fault-tolerant, with 2000 bootstrap reps:
python omie_price_maker_model_v3.py \
  --start 2024-01-01 --end 2024-12-31 \
  --zones ES,PT \
  --bootstrap-reps 2000 --ci-level 0.95 \
  --cache-dir raw_cache --no-download \
  --out out_v3_year

# Price-responsive buyer (constant-elasticity demand bid; elasticity supplied by user):
python omie_price_maker_model_v3.py \
  --start 2024-01-01 --end 2024-01-07 --country ES \
  --demand-bid elastic --demand-elasticity 0.3 --demand-ref-price 60 \
  --cache-dir raw_cache --no-download --out out_v3_elastic

# Quarter-hour resolution (when 15-minute inputs are available):
python omie_price_maker_model_v3.py \
  --start 2025-06-01 --end 2025-06-07 --country ES \
  --period-minutes 15 \
  --cache-dir raw_cache_15min --no-download --out out_v3_qh

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
from dataclasses import dataclass, replace
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
    # --- v3 additions (defaults reproduce v2 behaviour exactly) ---
    zones: Tuple[str, ...] = ()       # if multiple, run each Iberian zone and combine
    period_minutes: int = 60          # market time unit; 15 for quarter-hour inputs
    bootstrap_reps: int = 1000        # day-cluster bootstrap replications
    ci_level: float = 0.95            # confidence level for reported intervals
    seed: int = 0                     # RNG seed for reproducible bootstrap
    demand_bid: str = "inelastic"     # inelastic | elastic | piecewise
    demand_elasticity: float = 0.0    # used only when demand_bid == elastic (user supplied)
    demand_ref_price: float = 0.0     # reference price for the elasticity form (user supplied)
    demand_bid_file: Optional[Path] = None  # CSV of (price,quantity) points for piecewise


@dataclass(frozen=True)
class DemandBid:
    """Base class for a buyer's willingness-to-pay schedule.

    A demand bid maps a price to the quantity the buyer still wants at that price.
    It must be non-increasing in price. The default model uses the inelastic
    special case; richer shapes are supplied by the user, never invented here.
    """

    cap: float = math.inf  # buyer bid-price cap in EUR/MWh

    def quantity_at(self, price: float) -> float:  # pragma: no cover - abstract
        raise NotImplementedError

    def label(self) -> str:  # pragma: no cover - abstract
        raise NotImplementedError


@dataclass(frozen=True)
class InelasticBid(DemandBid):
    """Vertical demand line: q MWh wanted at any price up to the cap.

    This reproduces the v2 inelastic perturbation exactly.
    """

    q: float = 0.0

    def quantity_at(self, price: float) -> float:
        return self.q if price <= self.cap + 1e-9 else 0.0

    def label(self) -> str:
        return f"inelastic({self.q:g} MWh)"


@dataclass(frozen=True)
class PiecewiseLinearBid(DemandBid):
    """User-supplied piecewise-linear demand curve.

    ``points`` is a sequence of (price, quantity) breakpoints with non-increasing
    quantity. Linear interpolation is used between breakpoints; outside the range
    the nearest endpoint quantity is held. Every number here comes from the user.
    """

    points: Tuple[Tuple[float, float], ...] = ()

    def __post_init__(self) -> None:
        if len(self.points) < 2:
            raise ValueError("PiecewiseLinearBid needs at least two (price, quantity) points.")
        ps = [float(p) for p, _ in self.points]
        if any(b < a for a, b in zip(ps, ps[1:])):
            raise ValueError("PiecewiseLinearBid prices must be sorted ascending.")
        qs = [float(q) for _, q in self.points]
        if any(b > a + 1e-9 for a, b in zip(qs, qs[1:])):
            raise ValueError("PiecewiseLinearBid quantity must be non-increasing in price.")

    def quantity_at(self, price: float) -> float:
        xs = [float(p) for p, _ in self.points]
        ys = [float(q) for _, q in self.points]
        if price <= xs[0]:
            return ys[0]
        if price >= xs[-1]:
            return ys[-1]
        return float(np.interp(price, xs, ys))

    def label(self) -> str:
        return f"piecewise({len(self.points)} pts)"


@dataclass(frozen=True)
class ConstantElasticityBid(DemandBid):
    """Constant-elasticity demand bid q(p) = q_ref * (p / p_ref) ** (-elasticity).

    The elasticity and reference point are supplied by the user. The model does not
    assume any particular elasticity; it only implements the functional form so the
    user can test a hypothesis grounded in their own data or in the literature.
    """

    q_ref: float = 0.0
    p_ref: float = 1.0
    elasticity: float = 0.0  # >= 0; 0 reproduces inelastic demand at q_ref

    def __post_init__(self) -> None:
        if self.p_ref <= 0:
            raise ValueError("ConstantElasticityBid requires p_ref > 0.")
        if self.elasticity < 0:
            raise ValueError("ConstantElasticityBid elasticity must be non-negative.")

    def quantity_at(self, price: float) -> float:
        if price <= self.cap + 1e-9:
            p = max(float(price), 1e-9)
            return float(self.q_ref * (p / self.p_ref) ** (-self.elasticity))
        return 0.0

    def label(self) -> str:
        return f"elastic(e={self.elasticity:g}, q_ref={self.q_ref:g}@{self.p_ref:g})"


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

    def price_for_demand_bid(self, bid: "DemandBid") -> Tuple[float, float]:
        """Clearing price for a price-responsive buyer (limitation 3).

        Generalises ``price_anchored_net_supply``. The buyer submits a demand bid
        b(p): a non-increasing quantity-as-a-function-of-price schedule. The market
        must move up the observed incremental net-supply curve until the extra net
        supply available equals the quantity the buyer still demands at that price:

            P = min{ p >= p0 : N(p) - N(p0) >= b(p) }.

        Because N(p) - N(p0) is non-decreasing in p and b(p) is non-increasing in p,
        the gap [N(p) - N(p0)] - b(p) crosses zero at most once on the price grid, so
        the solution is well defined. For a vertical bid b(p) = q (for p below the
        buyer's cap) this returns exactly the same price as the inelastic anchored
        model, which is asserted in the test suite. No elasticity or bid shape is
        assumed here; the shape is carried entirely by the supplied ``bid`` object.

        Returns (price, unserved_mwh), where unserved_mwh is the residual demand that
        could not be cleared below the buyer's bid-price cap.
        """
        cap = float(min(self.bid_price, bid.cap))
        q0 = float(bid.quantity_at(self.p0))
        if q0 <= 1e-12:
            return float(self.p0), 0.0
        candidates = self._grid[(self._grid >= self.p0 - 1e-9) & (self._grid <= cap + 1e-9)]
        if candidates.size == 0:
            return float(self.p0), q0
        for p in candidates:
            avail = self.net_supply(float(p)) - self._net0
            need = float(bid.quantity_at(float(p)))
            if avail + 1e-9 >= need:
                return float(max(p, self.p0)), 0.0
        p_last = float(candidates[-1])
        unserved = max(0.0, float(bid.quantity_at(p_last)) - (self.net_supply(p_last) - self._net0))
        return p_last, unserved

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
# v3: demand bids (limitation 3)
# ---------------------------------------------------------------------------

def build_demand_bid(cfg: Config, order_mwh: float) -> DemandBid:
    """Construct the demand bid for a nominal order size, per the run configuration.

    The inelastic default reproduces v2. The elastic and piecewise forms use only
    parameters the user supplied on the command line or in a file.
    """
    kind = cfg.demand_bid.lower()
    if kind == "inelastic":
        return InelasticBid(cap=cfg.bid_price, q=float(order_mwh))
    if kind == "elastic":
        p_ref = cfg.demand_ref_price if cfg.demand_ref_price > 0 else None
        if p_ref is None:
            raise ValueError("--demand-bid elastic requires --demand-ref-price > 0 (supplied by the user).")
        return ConstantElasticityBid(
            cap=cfg.bid_price, q_ref=float(order_mwh), p_ref=float(p_ref), elasticity=float(cfg.demand_elasticity)
        )
    if kind == "piecewise":
        if cfg.demand_bid_file is None or not Path(cfg.demand_bid_file).exists():
            raise ValueError("--demand-bid piecewise requires --demand-bid-file pointing to a (price,quantity) CSV.")
        pts_df = pd.read_csv(cfg.demand_bid_file)
        cols = {c.lower(): c for c in pts_df.columns}
        if "price" not in cols or "quantity" not in cols:
            raise ValueError("Piecewise bid file must have 'price' and 'quantity' columns.")
        pts = tuple(
            (float(r[cols["price"]]), float(r[cols["quantity"]]))
            for _, r in pts_df.sort_values(cols["price"]).iterrows()
        )
        return PiecewiseLinearBid(cap=cfg.bid_price, points=pts)
    raise ValueError(f"Unknown demand bid kind: {cfg.demand_bid}")


def write_price_responsive_analysis(
    period_curves: Mapping[Tuple[pd.Timestamp, int], "PeriodCurve"], cfg: Config
) -> Optional[Path]:
    """Opt-in: clear a price-responsive demand bid against each period's net supply.

    Only runs when --demand-bid is not 'inelastic'. Writes its own CSV and never
    touches the core impact tables that feed the paper. The nominal order size sets
    the reference quantity of the bid; the realised cleared quantity can be smaller
    because a price-responsive buyer withdraws as the price rises.
    """
    if cfg.demand_bid.lower() == "inelastic":
        return None
    rows = []
    for key in sorted(period_curves):
        pc = period_curves[key]
        for q in cfg.order_sizes:
            bid = build_demand_bid(cfg, q)
            p_after, unserved = pc.price_for_demand_bid(bid)
            cleared = max(0.0, float(bid.quantity_at(p_after)) - unserved)
            rows.append(
                {
                    "date": pc.date,
                    "period": pc.period,
                    "zone": pc.zone,
                    "p0": pc.p0,
                    "bid_kind": cfg.demand_bid.lower(),
                    "bid_label": bid.label(),
                    "nominal_order_mwh": float(q),
                    "p_after": float(p_after),
                    "dp_eur_mwh": float(p_after - pc.p0),
                    "cleared_mwh": float(cleared),
                    "withdrawn_mwh": float(max(0.0, q - cleared)),
                    "unserved_mwh_at_bid_cap": float(unserved),
                }
            )
    df = pd.DataFrame(rows)
    path = cfg.out / "impact_price_responsive.csv"
    df.to_csv(path, index=False)
    return path


# ---------------------------------------------------------------------------
# v3: fault-tolerant data loading + coverage (limitation 2)
# ---------------------------------------------------------------------------

def load_days_robust(cfg: Config) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load every requested day, skipping (and logging) days that fail.

    A single missing file, an unparseable curve, or a missing official price no
    longer aborts the whole run. This is what makes a full-year or multi-year batch
    feasible. Returns (curves, prices, coverage). Behaviour on a fully clean range is
    identical to the v2 loader.
    """
    assert cfg.cache_dir is not None
    curve_frames: List[pd.DataFrame] = []
    price_frames: List[pd.DataFrame] = []
    coverage_rows: List[dict] = []
    for d in daterange(cfg.start, cfg.end):
        status, reason = "loaded", ""
        try:
            curve_path, price_path = fetch_day_files(d, cfg.cache_dir, cfg.no_download, cfg.overwrite)
            cdf = read_curve_file(curve_path)
            pdf = read_price_file(price_path)
            if cdf.empty or pdf.empty:
                raise ValueError("empty curve or price frame after parsing")
            curve_frames.append(cdf)
            price_frames.append(pdf)
        except FileNotFoundError as exc:
            status, reason = "skipped_missing_file", str(exc)
        except ValueError as exc:
            status, reason = "skipped_parse_error", str(exc)
        except Exception as exc:  # network or unexpected
            status, reason = "skipped_other", f"{type(exc).__name__}: {exc}"
        coverage_rows.append({"date": pd.Timestamp(d), "status": status, "reason": reason[:300]})

    coverage = pd.DataFrame(coverage_rows)
    curves = pd.concat(curve_frames, ignore_index=True) if curve_frames else pd.DataFrame()
    prices = pd.concat(price_frames, ignore_index=True) if price_frames else pd.DataFrame()
    return curves, prices, coverage


def summarise_resolution(period_curves: Mapping[Tuple[pd.Timestamp, int], "PeriodCurve"], cfg: Config) -> pd.DataFrame:
    """Report the market time unit actually seen in the data (limitation 4).

    Period counts are inferred from the data, not assumed to be 24, so the same code
    runs at quarter-hour resolution when 15-minute inputs are supplied.
    """
    by_day: Dict[pd.Timestamp, int] = {}
    for (dt, _period) in period_curves:
        by_day[dt] = by_day.get(dt, 0) + 1
    rows = [
        {
            "date": dt,
            "periods_in_day": n,
            "configured_period_minutes": cfg.period_minutes,
            "implied_minutes_per_period_if_24h": round(24 * 60 / n, 3) if n else np.nan,
        }
        for dt, n in sorted(by_day.items())
    ]
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# v3: statistical inference by day-cluster bootstrap (limitation 2)
# ---------------------------------------------------------------------------

def _percentile_ci(samples: np.ndarray, level: float) -> Tuple[float, float]:
    samples = samples[np.isfinite(samples)]
    if samples.size == 0:
        return (np.nan, np.nan)
    alpha = (1.0 - level) / 2.0
    lo = float(np.percentile(samples, 100.0 * alpha))
    hi = float(np.percentile(samples, 100.0 * (1.0 - alpha)))
    return (lo, hi)


def day_cluster_bootstrap(
    df: pd.DataFrame, value_fn, B: int, level: float, rng: np.random.Generator
) -> Tuple[float, float, float, int]:
    """Percentile bootstrap CI obtained by resampling whole days with replacement.

    Periods within a day are serially correlated, so the natural resampling unit is
    the day, not the individual period (a cluster bootstrap). ``value_fn`` maps a
    pooled resample to a scalar statistic. Returns (point_estimate, ci_low, ci_high,
    n_days). This adds honest uncertainty to every headline number without changing
    the point estimate, and uses only the user's own observations.
    """
    if df.empty:
        return (np.nan, np.nan, np.nan, 0)
    point = float(value_fn(df))
    groups = [g for _, g in df.groupby("date", sort=True)]
    n = len(groups)
    if n == 0:
        return (point, np.nan, np.nan, 0)
    stats = np.empty(B, dtype=float)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        sample = pd.concat([groups[i] for i in idx], ignore_index=True)
        try:
            stats[b] = float(value_fn(sample))
        except Exception:
            stats[b] = np.nan
    lo, hi = _percentile_ci(stats, level)
    return (point, lo, hi, n)


def write_statistical_inference(
    impact_long: pd.DataFrame, sched_summary: pd.DataFrame, rank_summary_src: pd.DataFrame, cfg: Config
) -> Path:
    """Write confidence intervals for the headline statistics (limitation 2)."""
    rng = np.random.default_rng(cfg.seed)
    B, lvl = cfg.bootstrap_reps, cfg.ci_level
    rows = []

    # 1) Price-impact location statistics per (mechanism, order size).
    for (mech, q), g in impact_long.groupby(["mechanism", "order_mwh"], sort=True):
        for stat_name, fn in (
            ("median_dp_eur_mwh", lambda x: x["dp_eur_mwh"].median()),
            ("mean_dp_eur_mwh", lambda x: x["dp_eur_mwh"].mean()),
            ("p95_dp_eur_mwh", lambda x: x["dp_eur_mwh"].quantile(0.95)),
            ("share_dp_ge_5eur", lambda x: (x["dp_eur_mwh"] >= 5).mean()),
        ):
            point, lo, hi, n = day_cluster_bootstrap(g, fn, B, lvl, rng)
            rows.append(
                {
                    "family": "price_impact",
                    "mechanism": mech,
                    "order_mwh": float(q),
                    "statistic": stat_name,
                    "point_estimate": point,
                    "ci_low": lo,
                    "ci_high": hi,
                    "ci_level": lvl,
                    "bootstrap_reps": B,
                    "n_days": n,
                }
            )

    # 2) Rank-reversal share per (mechanism, order size), bootstrapped over days.
    if rank_summary_src is not None and not rank_summary_src.empty:
        for (mech, q), g in rank_summary_src.groupby(["mechanism", "order_mwh"], sort=True):
            point, lo, hi, n = day_cluster_bootstrap(
                g, lambda x: 100.0 * x["cheapest_period_changes"].mean(), B, lvl, rng
            )
            rows.append(
                {
                    "family": "rank_reversal",
                    "mechanism": mech,
                    "order_mwh": float(q),
                    "statistic": "reversal_share_pct",
                    "point_estimate": point,
                    "ci_low": lo,
                    "ci_high": hi,
                    "ci_level": lvl,
                    "bootstrap_reps": B,
                    "n_days": n,
                }
            )

    # 3) Scheduling understatement and savings per mechanism, bootstrapped over days.
    if sched_summary is not None and not sched_summary.empty:
        for mech, g in sched_summary.groupby("mechanism", sort=True):
            def understatement(x):
                naive = x["price_taking_naive_cost_eur"].sum()
                actual = x["price_taking_actual_invoice_cost_eur"].sum()
                return (actual / naive - 1.0) * 100.0 if naive else np.nan

            def savings(x):
                pt = x["price_taking_actual_invoice_cost_eur"].sum()
                ia = x["impact_aware_actual_invoice_cost_eur"].sum()
                return (1.0 - ia / pt) * 100.0 if pt else np.nan

            for stat_name, fn in (
                ("naive_understatement_pct", understatement),
                ("impact_aware_savings_pct", savings),
            ):
                point, lo, hi, n = day_cluster_bootstrap(g, fn, B, lvl, rng)
                rows.append(
                    {
                        "family": "scheduling",
                        "mechanism": mech,
                        "order_mwh": np.nan,
                        "statistic": stat_name,
                        "point_estimate": point,
                        "ci_low": lo,
                        "ci_high": hi,
                        "ci_level": lvl,
                        "bootstrap_reps": B,
                        "n_days": n,
                    }
                )

    ci = pd.DataFrame(rows)
    path = cfg.out / "confidence_intervals.csv"
    ci.to_csv(path, index=False)
    return path


# ---------------------------------------------------------------------------
# v3: curve validation, anchoring sensitivity, impact bracket (limitation 1)
# ---------------------------------------------------------------------------

def _impact_at_anchor(pc: "PeriodCurve", anchor: float, q: float) -> float:
    """Anchored impact if the curve were anchored at ``anchor`` instead of p0.

    Used only to measure sensitivity to the anchoring choice. The perturbation size
    is the period's own raw-crossing residual, so no external magnitude is introduced.
    """
    if q <= 1e-12:
        return 0.0
    net_anchor = pc.net_supply(float(anchor))
    target = net_anchor + q
    cand = pc._grid[(pc._grid >= anchor - 1e-9) & (pc._grid <= pc.bid_price + 1e-9)]
    if cand.size == 0:
        return 0.0
    for p in cand:
        if pc.net_supply(float(p)) + 1e-9 >= target:
            return float(max(p, anchor) - anchor)
    return float(cand[-1] - anchor)


def write_curve_diagnostics(
    period_curves: Mapping[Tuple[pd.Timestamp, int], "PeriodCurve"], cfg: Config
) -> Tuple[Path, Path, Path]:
    """Make the anchoring assumption auditable (limitation 1).

    Writes three files:
      curve_validation.csv     - monotonicity of N(p) and the local slope dN/dp at p0.
      anchoring_sensitivity.csv- how much the anchored impact moves when the anchor is
                                 shifted by +/- the period's own raw-crossing residual.
      impact_bracket.csv       - the three mechanisms reorganised as lower/central/upper
                                 bounds, in the spirit of supply-curve-bounds robust
                                 self-scheduling.
    """
    grid_probe_eur = 1.0  # probe +1 EUR/MWh to estimate the local slope; not an output magnitude

    val_rows = []
    sens_rows = []
    brk_rows = []
    for key in sorted(period_curves):
        pc = period_curves[key]
        # --- monotonicity of the constructed net-supply curve on its own grid ---
        gp = pc._grid[(pc._grid >= pc.p0 - 1e-9) & (pc._grid <= pc.bid_price + 1e-9)]
        ns = np.array([pc.net_supply(float(p)) for p in gp], dtype=float)
        diffs = np.diff(ns) if ns.size > 1 else np.array([0.0])
        non_monotone = float((diffs < -1e-6).sum())
        max_decrease = float(-diffs.min()) if diffs.size and diffs.min() < 0 else 0.0
        slope = (pc.net_supply(pc.p0 + grid_probe_eur) - pc._net0) / grid_probe_eur  # MWh per EUR/MWh near p0
        val_rows.append(
            {
                "date": pc.date,
                "period": pc.period,
                "zone": pc.zone,
                "p0": pc.p0,
                "net_supply_at_p0": float(pc._net0),
                "grid_points_at_or_above_p0": int(gp.size),
                "non_monotone_steps": non_monotone,
                "max_net_supply_decrease_mwh": max_decrease,
                "local_depth_mwh_per_eur": float(slope),
                "raw_minus_official_eur_mwh": float(pc._raw0 - pc.p0),
            }
        )

        # --- anchoring-band sensitivity and three-mechanism bracket ---
        residual = abs(float(pc._raw0 - pc.p0))
        for q in cfg.order_sizes:
            dp_central = _impact_at_anchor(pc, pc.p0, q)
            dp_low_anchor = _impact_at_anchor(pc, pc.p0 - residual, q)
            dp_high_anchor = _impact_at_anchor(pc, pc.p0 + residual, q)
            sens_rows.append(
                {
                    "date": pc.date,
                    "period": pc.period,
                    "order_mwh": float(q),
                    "anchor_residual_eur_mwh": residual,
                    "dp_central_eur_mwh": dp_central,
                    "dp_anchor_minus_eur_mwh": dp_low_anchor,
                    "dp_anchor_plus_eur_mwh": dp_high_anchor,
                    "anchor_band_width_eur_mwh": abs(dp_high_anchor - dp_low_anchor),
                }
            )
            dps = {m: pc.price(q, m)[0] - pc.p0 for m in MECHANISMS}
            vals = [v for v in dps.values() if pd.notna(v)]
            brk_rows.append(
                {
                    "date": pc.date,
                    "period": pc.period,
                    "order_mwh": float(q),
                    "dp_anchored_net_supply": dps["anchored_net_supply"],
                    "dp_residual_supply": dps["residual_supply"],
                    "dp_raw_reclear_delta": dps["raw_reclear_delta"],
                    "dp_bracket_low": float(min(vals)) if vals else np.nan,
                    "dp_central": dps["anchored_net_supply"],
                    "dp_bracket_high": float(max(vals)) if vals else np.nan,
                }
            )

    val = pd.DataFrame(val_rows)
    sens = pd.DataFrame(sens_rows)
    brk = pd.DataFrame(brk_rows)
    p_val = cfg.out / "curve_validation.csv"
    p_sens = cfg.out / "anchoring_sensitivity.csv"
    p_brk = cfg.out / "impact_bracket.csv"
    val.to_csv(p_val, index=False)
    sens.to_csv(p_sens, index=False)
    brk.to_csv(p_brk, index=False)

    # Compact bracket summary by order size (central with min/max across mechanisms).
    if not brk.empty:
        bracket_summary = (
            brk.groupby("order_mwh")
            .agg(
                median_dp_central=("dp_central", "median"),
                median_dp_bracket_low=("dp_bracket_low", "median"),
                median_dp_bracket_high=("dp_bracket_high", "median"),
                p95_dp_central=("dp_central", lambda x: x.quantile(0.95)),
                p95_dp_bracket_high=("dp_bracket_high", lambda x: x.quantile(0.95)),
            )
            .reset_index()
        )
        bracket_summary.to_csv(cfg.out / "impact_bracket_summary.csv", index=False)

    return p_val, p_sens, p_brk


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[List[str]] = None) -> Config:
    parser = argparse.ArgumentParser(description="OMIE price-maker electricity optimization stress test, v3")
    parser.add_argument("--start", required=True, help="start date, YYYY-MM-DD")
    parser.add_argument("--end", required=True, help="end date, YYYY-MM-DD")
    parser.add_argument("--country", default="ES", choices=["ES", "PT"], help="buyer price zone (single-zone runs)")
    parser.add_argument("--zones", default=None, help="comma-separated zones to run and combine, e.g. ES,PT (overrides --country)")
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
    parser.add_argument("--out", default="out_omie_v3", help="output directory")
    parser.add_argument("--cache-dir", default=None, help="directory containing/downloading OMIE files; default OUT/cache")
    parser.add_argument("--no-download", action="store_true", help="use files already present in cache-dir")
    parser.add_argument("--overwrite", action="store_true", help="re-download files if present")
    # v3 options
    parser.add_argument("--period-minutes", type=int, default=60, help="market time unit in minutes (15 for quarter-hour inputs)")
    parser.add_argument("--bootstrap-reps", type=int, default=1000, help="day-cluster bootstrap replications for confidence intervals")
    parser.add_argument("--ci-level", type=float, default=0.95, help="confidence level for reported intervals")
    parser.add_argument("--seed", type=int, default=0, help="RNG seed for reproducible bootstrap")
    parser.add_argument("--demand-bid", default="inelastic", choices=["inelastic", "elastic", "piecewise"], help="buyer demand-bid shape")
    parser.add_argument("--demand-elasticity", type=float, default=0.0, help="elasticity for --demand-bid elastic (user supplied)")
    parser.add_argument("--demand-ref-price", type=float, default=0.0, help="reference price for --demand-bid elastic (user supplied)")
    parser.add_argument("--demand-bid-file", default=None, help="CSV with price,quantity columns for --demand-bid piecewise")
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

    if args.zones:
        zones = tuple(z.strip().upper() for z in args.zones.split(",") if z.strip())
        for z in zones:
            if z not in ("ES", "PT"):
                raise ValueError(f"Unknown zone {z}. Allowed: ES, PT")
    else:
        zones = (args.country.upper(),)

    if args.demand_bid == "elastic" and args.demand_ref_price <= 0:
        raise ValueError("--demand-bid elastic requires --demand-ref-price > 0 (a value you supply).")
    if args.demand_bid == "piecewise" and not args.demand_bid_file:
        raise ValueError("--demand-bid piecewise requires --demand-bid-file.")

    return Config(
        start=datetime.strptime(args.start, "%Y-%m-%d").date(),
        end=datetime.strptime(args.end, "%Y-%m-%d").date(),
        country=zones[0],
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
        zones=zones,
        period_minutes=args.period_minutes,
        bootstrap_reps=args.bootstrap_reps,
        ci_level=args.ci_level,
        seed=args.seed,
        demand_bid=args.demand_bid,
        demand_elasticity=args.demand_elasticity,
        demand_ref_price=args.demand_ref_price,
        demand_bid_file=Path(args.demand_bid_file) if args.demand_bid_file else None,
    )


def run_one_zone(curves: pd.DataFrame, prices: pd.DataFrame, cfg: Config) -> Optional[pd.DataFrame]:
    """Run the full analysis for a single zone into cfg.out. Returns impact_long.

    The first five CSVs are byte-for-byte what v2 produced for the same inputs. The
    remaining files are v3 additions and never alter those core outputs.
    """
    ensure_dir(cfg.out)
    period_curves = build_period_curves(curves, prices, cfg)
    if not period_curves:
        print(f"[{cfg.country}] No period curves built; skipping zone.", file=sys.stderr)
        return None

    # --- v2-identical core outputs ---
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

    # --- v3 additions ---
    summarise_resolution(period_curves, cfg).to_csv(cfg.out / "resolution_report.csv", index=False)
    # rank diagnostics were written by write_additional_diagnostics; reuse them for CIs.
    rank_path = cfg.out / "rank_reversal_diagnostics.csv"
    rank_src = pd.read_csv(rank_path) if rank_path.exists() else pd.DataFrame()
    write_statistical_inference(impact_long, sched_summary, rank_src, cfg)
    write_curve_diagnostics(period_curves, cfg)
    write_price_responsive_analysis(period_curves, cfg)
    print_console_summary(impact_long, sched_summary, cfg)
    return impact_long


def main(argv: Optional[List[str]] = None) -> int:
    cfg = parse_args(argv)
    ensure_dir(cfg.out)
    assert cfg.cache_dir is not None
    ensure_dir(cfg.cache_dir)

    # Fault-tolerant load once (curve files contain all zones); reuse across zones.
    curves, prices, coverage = load_days_robust(cfg)
    coverage.to_csv(cfg.out / "data_coverage.csv", index=False)
    n_loaded = int((coverage["status"] == "loaded").sum())
    n_total = int(len(coverage))
    if n_loaded == 0:
        raise RuntimeError("No days could be loaded. Check dates, cache-dir and files (see data_coverage.csv).")
    if n_loaded < n_total:
        print(f"Loaded {n_loaded}/{n_total} days; {n_total - n_loaded} skipped (see data_coverage.csv).", file=sys.stderr)

    zones = cfg.zones if cfg.zones else (cfg.country,)
    multi = len(zones) > 1
    for zone in zones:
        zcfg = replace(cfg, country=zone, zones=(zone,), out=(cfg.out / zone if multi else cfg.out),
                       cache_dir=cfg.cache_dir)
        if multi:
            ensure_dir(zcfg.out)
        run_one_zone(curves, prices, zcfg)

    zip_path = make_zip(cfg.out)
    print(f"\nSaved outputs in: {cfg.out.resolve()}")
    print(f"Saved bundle: {zip_path.resolve()}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
