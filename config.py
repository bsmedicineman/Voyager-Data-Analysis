"""Configuration: data products, their real instrument bands, and paths.

The band limits here are hardware facts, not preferences:
  * MAG (fluxgate) responds from DC (0 Hz). The low-field sensor's passband is
    0-8.3 Hz; archived averaged products sample far slower (1.92 s -> ~0.26 Hz
    Nyquist, 48 s -> ~0.01 Hz, hourly -> ~1.4e-4 Hz).
  * PWS waveform receiver is AC-coupled through a 40 Hz - 12 kHz analog bandpass
    sampled at 28800 S/s. There is NO data below 40 Hz; it was never digitized.
  * PWS spectrum analyzer = 16 fixed log-spaced channels, 10 Hz - 56.2 kHz,
    one spectrum / 4 s. These are intensities, not a waveform you can filter.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional


@dataclass
class Band:
    """Analysis band in Hz. lo <= 0 means 'keep DC' (low-pass only)."""
    lo: float
    hi: float

    def describe(self) -> str:
        lo = "DC" if self.lo <= 0 else f"{self.lo:g} Hz"
        return f"{lo} - {self.hi:g} Hz"


# How a product's samples are organized, which drives the analysis path.
KIND_TIMESERIES = "timeseries"   # uniformly-sampled scalar/vector (MAG, PWS waveform)
KIND_CHANNELS = "channels"       # pre-binned spectral channels (PWS SA)


@dataclass
class Product:
    name: str                    # short id used in filenames/tables
    spacecraft: str              # "voyager1" | "voyager2"
    instrument: str              # "MAG" | "PWS"
    kind: str                    # KIND_*
    band: Band                   # full instrument band we analyze across
    quantity: str                # "B" (nT, calibrated) | "E" (relative, uncalibrated)
    # Keywords used to auto-resolve the CDAWeb dataset id from the catalog.
    catalog_keywords: List[str] = field(default_factory=list)
    # Optional explicit override if you already know the exact CDAWeb id.
    dataset_id: Optional[str] = None
    notes: str = ""


@dataclass
class Paths:
    root: Path = Path("./voyager_em_out")

    @property
    def raw(self) -> Path: return self.root / "raw"          # cached CDF/netcdf chunks
    @property
    def kernels(self) -> Path: return self.root / "kernels"  # SPICE kernels
    @property
    def tables(self) -> Path: return self.root / "tables"    # anomaly + merged tables
    @property
    def figures(self) -> Path: return self.root / "figures"

    def ensure(self) -> None:
        for p in (self.raw, self.kernels, self.tables, self.figures):
            p.mkdir(parents=True, exist_ok=True)


def default_products() -> List[Product]:
    """One entry per (spacecraft, product). dataset_id left None -> auto-resolve."""
    prods: List[Product] = []
    for sc, scnum in (("voyager1", 1), ("voyager2", 2)):
        prods += [
            Product(
                name=f"{sc}_mag",
                spacecraft=sc, instrument="MAG", kind=KIND_TIMESERIES,
                band=Band(0.0, 8.3), quantity="B",
                catalog_keywords=[f"voyager{scnum}", "mag", "magnetic"],
                notes="Calibrated field (nT). Keep DC; this is where the "
                      "termination shock / heliopause signatures live.",
            ),
            Product(
                name=f"{sc}_pws_sa",
                spacecraft=sc, instrument="PWS", kind=KIND_CHANNELS,
                band=Band(10.0, 56200.0), quantity="E",
                catalog_keywords=[f"voyager{scnum}", "pws", "spectrum", "analyzer"],
                notes="16 log-spaced channels, 4 s cadence. Discrete intensities.",
            ),
            Product(
                name=f"{sc}_pws_wf",
                spacecraft=sc, instrument="PWS", kind=KIND_TIMESERIES,
                band=Band(40.0, 12000.0), quantity="E",
                catalog_keywords=[f"voyager{scnum}", "pws", "waveform"],
                notes="28800 S/s, RELATIVE amplitude (AGC, uncalibrated). "
                      "Sparse: ~weekly 48 s captures.",
            ),
        ]
    return prods


@dataclass
class Config:
    paths: Paths = field(default_factory=Paths)
    products: List[Product] = field(default_factory=default_products)
    # Full-mission default window; downloader chunks this internally.
    t_start: str = "1977-09-05T00:00:00Z"
    t_stop: str = "2025-01-01T00:00:00Z"
    download_chunk_days: int = 91
    # Anomaly detector knobs (robust z on rolling baseline).
    anomaly_z: float = 6.0
    anomaly_window: int = 256
