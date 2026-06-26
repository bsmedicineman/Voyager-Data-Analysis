"""Recompiled output: one table of every anomaly + the EM map figures.

Figure 1 (the chart you described): x = heliocentric distance (AU),
y = fluctuation amplitude, color = anomaly center frequency (log Hz), with
DC magnetometer shifts drawn distinctly along the bottom.
Figure 2: distance vs frequency, color = amplitude -- the full broad-spectrum
picture along each spacecraft's route.
"""
from __future__ import annotations
from pathlib import Path
from typing import List
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def build_table(rows: List[dict], tables_dir: Path) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    tables_dir.mkdir(parents=True, exist_ok=True)
    if not df.empty:
        df = df.sort_values(["spacecraft", "t_start"]).reset_index(drop=True)
        for col, base in (("t_start", "iso_start"), ("t_stop", "iso_stop")):
            df[base] = pd.to_datetime(df[col], unit="s", utc=True)
    df.to_csv(tables_dir / "anomalies.csv", index=False)
    try:
        df.to_parquet(tables_dir / "anomalies.parquet", index=False)
    except Exception:
        pass
    return df


def plot_map(df: pd.DataFrame, figures_dir: Path) -> List[Path]:
    figures_dir.mkdir(parents=True, exist_ok=True)
    out: List[Path] = []
    if df.empty or "distance_au" not in df:
        return out
    d = df.dropna(subset=["distance_au"]).copy()
    if d.empty:
        return out

    spectral = d[d["kind"] == "spectral"]
    dc = d[d["kind"] == "dc_shift"]
    markers = {"voyager1": "o", "voyager2": "^"}

    # ---- Figure 1: distance vs amplitude, colored by frequency ----
    fig, ax = plt.subplots(figsize=(11, 6))
    if not spectral.empty:
        cf = np.log10(spectral["center_freq_hz"].clip(lower=1e-3))
        for sc, mk in markers.items():
            s = spectral[spectral["spacecraft"] == sc]
            if s.empty:
                continue
            sc_cf = np.log10(s["center_freq_hz"].clip(lower=1e-3))
            sca = ax.scatter(s["distance_au"], s["amplitude"], c=sc_cf,
                             cmap="viridis", marker=mk, s=28, alpha=0.8,
                             vmin=cf.min(), vmax=cf.max(), label=sc)
        cb = fig.colorbar(sca, ax=ax)
        cb.set_label("anomaly center frequency  log10(Hz)")
    if not dc.empty:
        for sc, mk in markers.items():
            s = dc[dc["spacecraft"] == sc]
            if not s.empty:
                ax.scatter(s["distance_au"], s["amplitude"], marker=mk, s=34,
                           facecolors="none", edgecolors="crimson",
                           label=f"{sc} DC shift")
    ax.set_xlabel("heliocentric distance (AU)")
    ax.set_ylabel("fluctuation amplitude\n(PWS: relative power | MAG DC: |ΔB| nT)")
    ax.set_yscale("log")
    ax.set_title("Voyager broad-spectrum EM anomalies along the heliospheric route")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    f1 = figures_dir / "distance_amplitude_frequency.png"
    fig.tight_layout(); fig.savefig(f1, dpi=140); plt.close(fig)
    out.append(f1)

    # ---- Figure 2: distance vs frequency, colored by amplitude ----
    if not spectral.empty:
        fig, ax = plt.subplots(figsize=(11, 6))
        amp = np.log10(spectral["amplitude"].clip(lower=1e-30))
        for sc, mk in markers.items():
            s = spectral[spectral["spacecraft"] == sc]
            if s.empty:
                continue
            sca = ax.scatter(s["distance_au"], s["center_freq_hz"],
                             c=np.log10(s["amplitude"].clip(lower=1e-30)),
                             cmap="magma", marker=mk, s=28, alpha=0.8,
                             vmin=amp.min(), vmax=amp.max(), label=sc)
        ax.set_yscale("log")
        cb = fig.colorbar(sca, ax=ax); cb.set_label("amplitude  log10")
        ax.set_xlabel("heliocentric distance (AU)")
        ax.set_ylabel("center frequency (Hz)")
        ax.set_title("Frequency of EM anomalies vs distance from the Sun")
        ax.legend(fontsize=8); ax.grid(alpha=0.25)
        f2 = figures_dir / "distance_frequency.png"
        fig.tight_layout(); fig.savefig(f2, dpi=140); plt.close(fig)
        out.append(f2)
    return out
