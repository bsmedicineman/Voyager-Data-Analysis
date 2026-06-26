"""Command-line entry point.

  python -m voyager_em selftest      # offline end-to-end on synthetic data
  python -m voyager_em kernels       # fetch SPICE kernels (needs network)
  python -m voyager_em discover      # resolve CDAWeb dataset ids
  python -m voyager_em download      # pull all products (chunked, resumable)
  python -m voyager_em analyze       # spectra + anomaly detection -> table
  python -m voyager_em map           # attach positions + render figures
  python -m voyager_em all           # download -> analyze -> map
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np

from .config import Config, Product, Band, KIND_TIMESERIES
from .ingest import Series, load_series
from . import anomaly, mapping


def _analyze(cfg: Config):
    rows = []
    for p in cfg.products:
        s = load_series(p, cfg.paths.raw)
        if s.empty:
            print(f"[{p.name}] no cached data"); continue
        a = anomaly.detect(s, cfg.anomaly_z, cfg.anomaly_window)
        print(f"[{p.name}] {len(a)} anomalies ({p.band.describe()})")
        rows += anomaly.to_rows(a)
    return rows


def _map(cfg: Config, rows):
    from .position import attach_distance
    try:
        rows = attach_distance(rows, cfg.paths.kernels)
    except Exception as e:
        print(f"[position] skipped ({e}); run 'kernels' first for distances")
    df = mapping.build_table(rows, cfg.paths.tables)
    figs = mapping.plot_map(df, cfg.paths.figures)
    print(f"table: {cfg.paths.tables/'anomalies.csv'}")
    for f in figs:
        print(f"figure: {f}")
    return df


def _selftest(cfg: Config):
    """Synthetic E/B data + fake trajectory -> exercises the whole pipeline."""
    cfg.paths.ensure()
    rng = np.random.default_rng(0)
    rows = []

    # Synthetic PWS waveform: broadband noise + frequency-drifting tone bursts.
    fs = 28800.0
    n = int(fs * 40)
    t = np.arange(n) / fs
    x = rng.normal(0, 1, n)
    for (start, dur, f0) in [(5, 1.5, 80), (15, 1.0, 140), (28, 2.0, 60)]:
        m = (t >= start) & (t < start + dur)
        x[m] += 12 * np.sin(2 * np.pi * f0 * t[m])
    wf = Series(Product("voyager1_pws_wf", "voyager1", "PWS", KIND_TIMESERIES,
                        Band(40, 12000), "E"),
                times=1.0e9 + t, values=x, fs=fs)

    # Synthetic MAG magnitude: slow 1/r falloff + ULF wiggle + a shock step.
    fsB = 1 / 1.92
    nB = int(fsB * 3600 * 24 * 10)
    tB = np.arange(nB) / fsB
    B = 0.5 + 2.0 / (1 + tB / 5e5) + 0.02 * np.sin(2 * np.pi * 0.03 * tB)
    B += rng.normal(0, 0.01, nB)
    B[nB // 2:] += 0.4                         # termination-shock-like step
    mag = Series(Product("voyager1_mag", "voyager1", "MAG", KIND_TIMESERIES,
                         Band(0, 8.3), "B"),
                 times=1.0e9 + tB, values=B, fs=fsB)

    for s in (wf, mag):
        a = anomaly.detect(s, cfg.anomaly_z, cfg.anomaly_window)
        print(f"[{s.product.name}] {len(a)} anomalies")
        rows += anomaly.to_rows(a)

    # Fake heliocentric distances (AU) so the map renders without SPICE.
    for r in rows:
        r["distance_au"] = float(20 + (r["t_start"] - 1.0e9) / 5e6 + rng.normal(0, 0.3))

    df = mapping.build_table(rows, cfg.paths.tables)
    figs = mapping.plot_map(df, cfg.paths.figures)
    print(f"\n{len(df)} total anomalies -> {cfg.paths.tables/'anomalies.csv'}")
    for f in figs:
        print(f"figure: {f}")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="voyager_em")
    ap.add_argument("cmd", choices=["selftest", "kernels", "discover",
                                    "download", "analyze", "map", "all"])
    ap.add_argument("--out", default="./voyager_em_out")
    ap.add_argument("--only", nargs="*", help="limit to product names")
    args = ap.parse_args(argv)

    cfg = Config()
    cfg.paths.root = Path(args.out)
    cfg.paths.ensure()

    if args.cmd == "selftest":
        _selftest(cfg)
    elif args.cmd == "kernels":
        from .position import fetch_kernels
        fetch_kernels(cfg.paths.kernels)
    elif args.cmd == "discover":
        from cdasws import CdasWs
        from .catalog import resolve_dataset_ids
        ids = resolve_dataset_ids(cfg, CdasWs(), refresh=True)
        for k, v in ids.items():
            print(f"  {k:18s} -> {v}")
    elif args.cmd == "download":
        from .download import download
        download(cfg, only=args.only)
    elif args.cmd == "analyze":
        _map  # noqa
        rows = _analyze(cfg)
        mapping.build_table(rows, cfg.paths.tables)
    elif args.cmd == "map":
        rows = _analyze(cfg)
        _map(cfg, rows)
    elif args.cmd == "all":
        from .download import download
        download(cfg, only=args.only)
        _map(cfg, _analyze(cfg))


if __name__ == "__main__":
    main()
