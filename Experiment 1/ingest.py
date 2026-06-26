"""Load cached chunks into a normalized container the analysis stages share.

We deliberately stay schema-agnostic: CDAWeb variable names differ per product,
so we pick variables heuristically by dimensionality and unit hints rather than
hard-coding names. Override `pick_*` if your dataset uses unusual names.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
import numpy as np

from .config import Product, KIND_TIMESERIES, KIND_CHANNELS


@dataclass
class Series:
    """Normalized product data. Times are float seconds since 1970 (UTC)."""
    product: Product
    times: np.ndarray                          # (N,) epoch seconds
    values: np.ndarray                         # (N,) scalar  OR (N,F) channels
    freqs: Optional[np.ndarray] = None         # (F,) for KIND_CHANNELS
    fs: Optional[float] = None                 # sampling rate (Hz) for timeseries
    meta: Dict = field(default_factory=dict)

    @property
    def empty(self) -> bool:
        return self.times.size == 0


def _epoch_seconds(time_da) -> np.ndarray:
    t = np.asarray(time_da.values)
    if np.issubdtype(t.dtype, np.datetime64):
        return t.astype("datetime64[ns]").astype("int64") / 1e9
    return t.astype(float)


def _load_xr(outdir: Path):
    import xarray as xr
    files = sorted(outdir.glob("*.nc"))
    if not files:
        return None
    parts = []
    for f in files:
        try:
            parts.append(xr.open_dataset(f))
        except Exception:
            continue
    if not parts:
        return None
    # concat along the primary record dim (first dim of first data var)
    dim = list(parts[0].dims)[0]
    try:
        return xr.concat(parts, dim=dim)
    except Exception:
        return parts[0]


def _pick_time(ds):
    for name in ds.coords:
        if np.issubdtype(np.asarray(ds[name].values).dtype, np.datetime64):
            return ds[name]
    for cand in ("Epoch", "epoch", "time", "Time"):
        if cand in ds:
            return ds[cand]
    return ds[list(ds.dims)[0]]


def _pick_scalar(ds, unit_hints):
    """Best 1-D numeric data variable, preferring matching unit hints."""
    best, best_score = None, -1
    for name, da in ds.data_vars.items():
        if da.ndim != 1:
            continue
        units = str(da.attrs.get("UNITS", da.attrs.get("units", ""))).lower()
        score = sum(h in units for h in unit_hints) + (0.1 if "mag" in name.lower() else 0)
        if score > best_score:
            best, best_score = da, score
    return best


def _pick_channels(ds):
    """A 2-D (time, freq) variable + its frequency axis, for spectrum analyzers."""
    chan = None
    for name, da in ds.data_vars.items():
        if da.ndim == 2:
            chan = da
            break
    if chan is None:
        return None, None
    freqs = None
    for cand in ("frequency", "Frequency", "freq", "FREQUENCY"):
        if cand in ds.coords or cand in ds:
            freqs = np.asarray(ds[cand].values, dtype=float).ravel()
            break
    if freqs is None:  # last resort: nominal 16-channel PWS-SA grid (Hz)
        freqs = 10.0 * 10 ** (np.arange(16) / 4.0)
    return chan, freqs


def load_series(product: Product, raw_root: Path) -> Series:
    outdir = raw_root / product.name
    ds = _load_xr(outdir)
    if ds is None:
        return Series(product, np.array([]), np.array([]))

    time_da = _pick_time(ds)
    times = _epoch_seconds(time_da)

    if product.kind == KIND_CHANNELS:
        chan, freqs = _pick_channels(ds)
        if chan is None:
            return Series(product, np.array([]), np.array([]))
        vals = np.asarray(chan.values, dtype=float)
        if vals.shape[0] != times.shape[0] and vals.shape[-1] == times.shape[0]:
            vals = vals.T
        return Series(product, times, vals, freqs=freqs,
                      meta={"source": str(outdir)})

    unit_hints = ["nt"] if product.quantity == "B" else ["v/m", "v", "volt"]
    da = _pick_scalar(ds, unit_hints)
    if da is None:
        return Series(product, np.array([]), np.array([]))
    vals = np.asarray(da.values, dtype=float)
    fs = None
    if times.size > 4:
        dt = np.median(np.diff(times))
        if dt > 0:
            fs = 1.0 / dt
    return Series(product, times, vals, fs=fs, meta={"source": str(outdir)})
