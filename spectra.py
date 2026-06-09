"""Broad-spectrum signal processing.

For a true time series (MAG, PWS waveform) we compute a spectrogram across the
*entire* valid instrument band (no narrow pre-filter) so anomalies are found
wherever they occur. A bandpass helper is provided for optional isolation; note
a band whose low edge is <= 0 becomes a low-pass (the correct way to "keep DC").
"""
from __future__ import annotations
from typing import Tuple
import numpy as np
from scipy import signal

from .config import Band

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


def make_bandpass(fs: float, band: Band, order: int = 4):
    nyq = fs / 2.0
    lo = max(band.lo, 0.0)
    hi = min(band.hi, nyq * 0.999)
    if lo <= 0 and hi < nyq * 0.999:
        return signal.butter(order, hi / nyq, btype="low", output="sos")
    if lo > 0 and hi >= nyq * 0.999:
        return signal.butter(order, lo / nyq, btype="high", output="sos")
    if lo > 0:
        return signal.butter(order, [lo / nyq, hi / nyq], btype="band", output="sos")
    return None  # band spans the whole spectrum -> no filtering needed


def bandpass(x: np.ndarray, fs: float, band: Band) -> np.ndarray:
    sos = make_bandpass(fs, band)
    if sos is None:
        return x
    return signal.sosfiltfilt(sos, x)


def welch_psd(x: np.ndarray, fs: float, nperseg: int = 1024
              ) -> Tuple[np.ndarray, np.ndarray]:
    nperseg = int(min(nperseg, len(x)))
    if nperseg < 8:
        return np.array([]), np.array([])
    f, pxx = signal.welch(x, fs=fs, nperseg=nperseg, detrend="constant")
    return f, pxx


def spectrogram(x: np.ndarray, fs: float, nperseg: int = 1024, overlap: float = 0.5
                ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (freqs, frame_times[s], Sxx[freq,time]) power spectral density."""
    nperseg = int(min(nperseg, len(x)))
    if nperseg < 8:
        return np.array([]), np.array([]), np.zeros((0, 0))
    f, t, Sxx = signal.spectrogram(
        x, fs=fs, nperseg=nperseg, noverlap=int(nperseg * overlap),
        detrend="constant", mode="psd")
    return f, t, Sxx


def band_power(f: np.ndarray, Sxx: np.ndarray, band: Band) -> np.ndarray:
    """Integrate PSD over [band.lo, band.hi] -> one value per spectrogram frame."""
    lo = max(band.lo, f[0]) if f.size else band.lo
    sel = (f >= lo) & (f <= band.hi)
    if not sel.any():
        return np.zeros(Sxx.shape[1])
    return _trapz(Sxx[sel, :], f[sel], axis=0)
