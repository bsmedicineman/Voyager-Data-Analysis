"""Anomaly detection.

Two complementary detectors:
  * spectral events  -> for any time series with usable AC content. Flags
    spectrogram frames whose band power exceeds a robust rolling threshold,
    groups consecutive flagged frames into events, and reports each event's
    center frequency, bandwidth (FWHM), start/stop and duration.
  * DC baseline shifts -> for the magnetometer's near-DC layer, where the
    interesting features (shocks, heliopause) are step changes / variance
    spikes in the field magnitude rather than oscillations at a frequency.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from typing import List
import numpy as np

from .config import Product, Band
from .ingest import Series
from . import spectra


@dataclass
class Anomaly:
    product: str
    spacecraft: str
    instrument: str
    quantity: str             # "E" | "B"
    t_start: float            # epoch s
    t_stop: float
    duration_s: float
    center_freq_hz: float     # 0.0 for DC baseline shifts
    bandwidth_hz: float
    amplitude: float          # band power (E, relative) or |dB| step (B, nT)
    kind: str                 # "spectral" | "dc_shift"


def _robust_z(x: np.ndarray, window: int) -> np.ndarray:
    """Rolling median/MAD z-score; robust to the very spikes we want to find."""
    n = len(x)
    if n == 0:
        return x
    w = max(8, min(window, n))
    z = np.zeros(n)
    half = w // 2
    for i in range(n):
        a, b = max(0, i - half), min(n, i + half)
        seg = x[a:b]
        med = np.median(seg)
        mad = np.median(np.abs(seg - med)) + 1e-30
        z[i] = 0.6745 * (x[i] - med) / mad
    return z


def _group(flags: np.ndarray) -> List[tuple]:
    events, i, n = [], 0, len(flags)
    while i < n:
        if flags[i]:
            j = i
            while j < n and flags[j]:
                j += 1
            events.append((i, j - 1))
            i = j
        else:
            i += 1
    return events


def _fwhm(f: np.ndarray, spec: np.ndarray) -> tuple:
    base = np.median(spec)
    spec = spec - base
    k = int(np.argmax(spec))
    peak = spec[k]
    if peak <= 0:
        return float(f[k]), 0.0
    half = peak / 2.0
    lo = k
    while lo > 0 and spec[lo] > half:
        lo -= 1
    hi = k
    while hi < len(spec) - 1 and spec[hi] > half:
        hi += 1
    return float(f[k]), float(f[hi] - f[lo])


def detect_spectral(series: Series, z_thresh: float, window: int) -> List[Anomaly]:
    p = series.product
    if series.empty or series.fs is None:
        return []
    x = series.values.astype(float)
    x = x[np.isfinite(x)]
    if x.size < 64:
        return []
    fs = series.fs
    f, t, Sxx = spectra.spectrogram(x, fs)
    if Sxx.size == 0:
        return []
    bp = spectra.band_power(f, Sxx, p.band)
    z = _robust_z(np.log10(bp + 1e-30), window)
    out: List[Anomaly] = []
    t0 = series.times[0]
    frame_dt = (t[1] - t[0]) if len(t) > 1 else (len(x) / fs)
    for a, b in _group(z > z_thresh):
        spec = Sxx[:, a:b + 1].mean(axis=1)
        cf, bw = _fwhm(f, spec)
        ts, te = t0 + t[a] - frame_dt / 2, t0 + t[b] + frame_dt / 2
        out.append(Anomaly(
            product=p.name, spacecraft=p.spacecraft, instrument=p.instrument,
            quantity=p.quantity, t_start=ts, t_stop=te, duration_s=te - ts,
            center_freq_hz=cf, bandwidth_hz=bw,
            amplitude=float(bp[a:b + 1].max()), kind="spectral"))
    return out


def detect_dc_shifts(series: Series, z_thresh: float, window: int) -> List[Anomaly]:
    """Step changes / variance spikes in a near-DC magnetometer magnitude."""
    p = series.product
    if series.empty or p.quantity != "B":
        return []
    x = series.values.astype(float)
    m = np.isfinite(x)
    x, t = x[m], series.times[m]
    if x.size < 16:
        return []
    dx = np.abs(np.diff(x, prepend=x[0]))      # level changes
    z = _robust_z(dx, window)
    out: List[Anomaly] = []
    for a, b in _group(z > z_thresh):
        out.append(Anomaly(
            product=p.name, spacecraft=p.spacecraft, instrument=p.instrument,
            quantity="B", t_start=float(t[a]), t_stop=float(t[b]),
            duration_s=float(t[b] - t[a]), center_freq_hz=0.0, bandwidth_hz=0.0,
            amplitude=float(dx[a:b + 1].max()), kind="dc_shift"))
    return out


def detect(series: Series, z_thresh: float, window: int) -> List[Anomaly]:
    out = detect_spectral(series, z_thresh, window)
    if series.product.band.lo <= 0:               # DC-inclusive product
        out += detect_dc_shifts(series, z_thresh, window)
    return out


def to_rows(anoms: List[Anomaly]) -> List[dict]:
    return [asdict(a) for a in anoms]
