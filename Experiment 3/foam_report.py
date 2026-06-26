#!/usr/bin/env python3
r"""
Foam Score Report Generator
===========================

Reads the large scored-output file produced by the Voyager foam analyzer
(`*__streaming_scores.csv.gz`) and produces a human-readable report WITHOUT
ever loading the whole file into memory.

Why this is needed:
    The analyzer appended the file one row at a time with gzip, which creates
    millions of tiny stacked gzip "members". That bloats the file on disk
    (a few-hundred-MB dataset can balloon past 2 GB) and makes it impossible
    to open in a spreadsheet. But the file is still a normal multi-member
    gzip stream, so pandas can read every row by streaming in chunks. This
    script does exactly that with a fixed, small memory footprint:
        - foam_score distribution is tracked with a 1000-bin histogram
          (percentiles are estimated from it; no need to store every value)
        - the highest-scoring rows are kept in a fixed-size top-N heap
        - everything else is running sums / counts

Outputs (written next to wherever you run it, in ./foam_analysis_report/):
    1. foam_report.md          - the full readable report
    2. foam_top_rows.csv        - the N highest foam_score rows (opens in Excel)
    3. foam_score_hist.png      - histogram of foam_score (if matplotlib present)

Usage:
    python3 foam_report.py --csv "anomalies_full__streaming_scores.csv.gz"
    python3 foam_report.py            # will prompt for the path

Works on either the .csv.gz (gzipped) or a plain .csv - it sniffs the file.
"""

from pathlib import Path
import argparse
import heapq
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
OUT_DIR = Path("./foam_analysis_report")
CHUNK_ROWS = 200_000          # rows per streaming chunk
SCORE_BINS = 1000             # resolution of the foam_score histogram (0..1)
TOP_N = 500                   # how many top-scoring rows to keep + export
PROGRESS_EVERY = 2_000_000    # print a progress line every N rows

# Score thresholds to report exact counts for
SCORE_THRESHOLDS = [0.01, 0.05, 0.10, 0.20, 0.30, 0.40,
                    0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99]

# Frequency bands (Hz) for a coarse breakdown
FREQ_BAND_EDGES = [-np.inf, 0, 8, 20, 100, 1_000, 10_000, np.inf]
FREQ_BAND_LABELS = ["< 0 Hz", "0-8 Hz", "8-20 Hz", "20-100 Hz",
                    "100 Hz-1 kHz", "1-10 kHz", "> 10 kHz"]

EXPECTED_COLS = ["row_number", "timestamp", "center_freq_hz", "amplitude",
                 "distance_au", "foam_score", "spacecraft", "source"]


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def sniff_compression(path: Path) -> str:
    """Return 'gzip' if the file starts with the gzip magic bytes, else 'infer'."""
    try:
        with open(path, "rb") as f:
            magic = f.read(2)
        return "gzip" if magic == b"\x1f\x8b" else "infer"
    except Exception:
        return "infer"


def percentile_from_hist(edges: np.ndarray, counts: np.ndarray, q: float) -> float:
    """Estimate the q-th percentile (0-100) from a histogram via linear interpolation."""
    total = counts.sum()
    if total == 0:
        return float("nan")
    target = (q / 100.0) * total
    cum = np.cumsum(counts)
    idx = int(np.searchsorted(cum, target, side="left"))
    idx = min(idx, len(counts) - 1)
    lower_cum = cum[idx - 1] if idx > 0 else 0.0
    bin_count = counts[idx] if counts[idx] > 0 else 1
    frac = (target - lower_cum) / bin_count
    frac = min(max(frac, 0.0), 1.0)
    return float(edges[idx] + frac * (edges[idx + 1] - edges[idx]))


def fmt_int(n) -> str:
    return f"{int(n):,}"


# ----------------------------------------------------------------------
# Aggregator (constant memory)
# ----------------------------------------------------------------------
class Aggregator:
    def __init__(self):
        self.n_rows = 0

        # foam_score
        self.score_hist_edges = np.linspace(0.0, 1.0, SCORE_BINS + 1)
        self.score_hist = np.zeros(SCORE_BINS, dtype=np.int64)
        self.score_n = 0
        self.score_sum = 0.0
        self.score_sumsq = 0.0
        self.score_min = np.inf
        self.score_max = -np.inf
        self.thresh_counts = {t: 0 for t in SCORE_THRESHOLDS}

        # top-N rows by foam_score (min-heap of (score, tiebreak, row_tuple))
        self.top_heap = []
        self._tie = 0

        # center_freq_hz
        self.f_n = 0
        self.f_sum = 0.0
        self.f_min = np.inf
        self.f_max = -np.inf
        self.f_band_counts = np.zeros(len(FREQ_BAND_LABELS), dtype=np.int64)
        self.f_near8 = 0      # within +/-10% of 8 Hz
        self.f_near20 = 0     # within +/-10% of 20 Hz

        # amplitude (orders of magnitude -> track log10)
        self.a_n = 0
        self.a_logsum = 0.0
        self.a_min = np.inf
        self.a_max = -np.inf

        # distance_au
        self.d_finite = 0
        self.d_nan = 0
        self.d_sum = 0.0
        self.d_min = np.inf
        self.d_max = -np.inf

        # categoricals
        self.spacecraft = {}
        self.source = {}

        # timestamps (fixed-width strings sort chronologically) + row_number
        self.ts_min = None
        self.ts_max = None
        self.ts_epoch1970 = 0
        self.rn_min = np.inf
        self.rn_max = -np.inf

    # --- per-chunk update ---
    def update(self, df: pd.DataFrame):
        n = len(df)
        self.n_rows += n

        # ---- foam_score ----
        s = pd.to_numeric(df.get("foam_score"), errors="coerce").to_numpy(dtype=float)
        s = s[np.isfinite(s)]
        if s.size:
            sc = np.clip(s, 0.0, 1.0)
            self.score_hist += np.histogram(sc, bins=self.score_hist_edges)[0]
            self.score_n += s.size
            self.score_sum += float(s.sum())
            self.score_sumsq += float(np.square(s).sum())
            self.score_min = min(self.score_min, float(s.min()))
            self.score_max = max(self.score_max, float(s.max()))
            for t in SCORE_THRESHOLDS:
                self.thresh_counts[t] += int((s >= t).sum())

        # ---- top-N rows ----
        sc_full = pd.to_numeric(df.get("foam_score"), errors="coerce").to_numpy(dtype=float)
        order = np.argsort(sc_full)  # ascending; tail = highest
        take = order[-min(TOP_N, len(order)):]
        for i in take:
            val = sc_full[i]
            if not np.isfinite(val):
                continue
            if len(self.top_heap) < TOP_N:
                heapq.heappush(self.top_heap, (val, self._tie, self._row_tuple(df, i)))
                self._tie += 1
            elif val > self.top_heap[0][0]:
                heapq.heapreplace(self.top_heap, (val, self._tie, self._row_tuple(df, i)))
                self._tie += 1

        # ---- center_freq_hz ----
        f = pd.to_numeric(df.get("center_freq_hz"), errors="coerce").to_numpy(dtype=float)
        f = f[np.isfinite(f)]
        if f.size:
            self.f_n += f.size
            self.f_sum += float(f.sum())
            self.f_min = min(self.f_min, float(f.min()))
            self.f_max = max(self.f_max, float(f.max()))
            self.f_band_counts += np.histogram(f, bins=FREQ_BAND_EDGES)[0]
            self.f_near8 += int((np.abs(f - 8.0) <= 0.8).sum())
            self.f_near20 += int((np.abs(f - 20.0) <= 2.0).sum())

        # ---- amplitude ----
        a = pd.to_numeric(df.get("amplitude"), errors="coerce").to_numpy(dtype=float)
        a = a[np.isfinite(a) & (a > 0)]
        if a.size:
            self.a_n += a.size
            self.a_logsum += float(np.log10(a).sum())
            self.a_min = min(self.a_min, float(a.min()))
            self.a_max = max(self.a_max, float(a.max()))

        # ---- distance_au ----
        d = pd.to_numeric(df.get("distance_au"), errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(d)
        self.d_nan += int((~finite).sum())
        df_ = d[finite]
        if df_.size:
            self.d_finite += df_.size
            self.d_sum += float(df_.sum())
            self.d_min = min(self.d_min, float(df_.min()))
            self.d_max = max(self.d_max, float(df_.max()))

        # ---- categoricals ----
        if "spacecraft" in df:
            for k, v in df["spacecraft"].astype(str).value_counts().items():
                self.spacecraft[k] = self.spacecraft.get(k, 0) + int(v)
        if "source" in df:
            for k, v in df["source"].astype(str).value_counts().items():
                self.source[k] = self.source.get(k, 0) + int(v)

        # ---- timestamps + row_number ----
        if "timestamp" in df:
            ts = df["timestamp"].astype(str)
            cmn, cmx = ts.min(), ts.max()
            self.ts_min = cmn if self.ts_min is None else min(self.ts_min, cmn)
            self.ts_max = cmx if self.ts_max is None else max(self.ts_max, cmx)
            self.ts_epoch1970 += int(ts.str.startswith("1970-01-01").sum())
        if "row_number" in df:
            rn = pd.to_numeric(df["row_number"], errors="coerce").to_numpy(dtype=float)
            rn = rn[np.isfinite(rn)]
            if rn.size:
                self.rn_min = min(self.rn_min, float(rn.min()))
                self.rn_max = max(self.rn_max, float(rn.max()))

    def _row_tuple(self, df, i):
        return tuple(df.iloc[i].get(c, "") for c in EXPECTED_COLS)

    # --- finalize ---
    def top_rows_df(self) -> pd.DataFrame:
        rows = sorted(self.top_heap, key=lambda x: x[0], reverse=True)
        return pd.DataFrame([r[2] for r in rows], columns=EXPECTED_COLS)


# ----------------------------------------------------------------------
# Report writer
# ----------------------------------------------------------------------
def build_report(agg: Aggregator, csv_path: Path, elapsed: float) -> str:
    L = []
    A = L.append

    score_mean = agg.score_sum / agg.score_n if agg.score_n else float("nan")
    score_var = (agg.score_sumsq / agg.score_n - score_mean ** 2) if agg.score_n else float("nan")
    score_std = np.sqrt(max(score_var, 0.0)) if agg.score_n else float("nan")
    pct = {q: percentile_from_hist(agg.score_hist_edges, agg.score_hist, q)
           for q in [1, 5, 25, 50, 75, 90, 95, 99, 99.9]}

    f_mean = agg.f_sum / agg.f_n if agg.f_n else float("nan")
    a_geo = 10 ** (agg.a_logsum / agg.a_n) if agg.a_n else float("nan")
    d_mean = agg.d_sum / agg.d_finite if agg.d_finite else float("nan")

    A("# Foam Score Report\n")
    A(f"**Source file:** `{csv_path.name}`  ")
    A(f"**File size on disk:** {csv_path.stat().st_size / 1e9:.2f} GB  ")
    A(f"**Rows scored:** {fmt_int(agg.n_rows)}  ")
    A(f"**Processed in:** {elapsed:.1f}s\n")

    # coverage
    A("## Coverage\n")
    A(f"- Input row numbers span **{fmt_int(agg.rn_min)} -> {fmt_int(agg.rn_max)}** "
      f"(these are the original line numbers from `anomalies_full.csv`).")
    A(f"- Rows with a usable distance value: **{fmt_int(agg.d_finite)}**; "
      f"rows with missing distance (the PARTIAL rows): **{fmt_int(agg.d_nan)}**.")
    if agg.ts_min is not None:
        A(f"- Timestamp range (as stored): `{agg.ts_min}` -> `{agg.ts_max}`.")
        if agg.ts_epoch1970:
            A(f"  - Note: {fmt_int(agg.ts_epoch1970)} rows carry a `1970-01-01` (Unix-epoch) "
              f"timestamp, i.e. the original time field did not parse into a real date for those rows.")
    A("")

    # the headline: foam_score distribution
    A("## Foam score distribution\n")
    A("| Statistic | Value |")
    A("|---|---|")
    A(f"| count | {fmt_int(agg.score_n)} |")
    A(f"| min | {agg.score_min:.6f} |")
    A(f"| max | {agg.score_max:.6f} |")
    A(f"| mean | {score_mean:.6f} |")
    A(f"| std dev | {score_std:.6f} |")
    A(f"| 1st pct | {pct[1]:.6f} |")
    A(f"| 5th pct | {pct[5]:.6f} |")
    A(f"| 25th pct | {pct[25]:.6f} |")
    A(f"| median (50th) | {pct[50]:.6f} |")
    A(f"| 75th pct | {pct[75]:.6f} |")
    A(f"| 90th pct | {pct[90]:.6f} |")
    A(f"| 95th pct | {pct[95]:.6f} |")
    A(f"| 99th pct | {pct[99]:.6f} |")
    A(f"| 99.9th pct | {pct[99.9]:.6f} |")
    A("")

    A("### Rows at or above each score threshold\n")
    A("| Threshold | Rows >= | % of scored |")
    A("|---|---|---|")
    for t in SCORE_THRESHOLDS:
        c = agg.thresh_counts[t]
        pctg = 100 * c / agg.score_n if agg.score_n else 0
        A(f"| {t:.2f} | {fmt_int(c)} | {pctg:.4f}% |")
    A("")

    # frequency
    A("## Frequency (center_freq_hz)\n")
    A(f"- Range: **{agg.f_min:.3f} Hz -> {fmt_int(agg.f_max)} Hz**, mean **{f_mean:.2f} Hz** "
      f"(over {fmt_int(agg.f_n)} rows with a finite frequency).")
    A(f"- Within +/-10% of 8 Hz: **{fmt_int(agg.f_near8)}** rows; "
      f"within +/-10% of 20 Hz: **{fmt_int(agg.f_near20)}** rows.")
    A("")
    A("| Band | Rows |")
    A("|---|---|")
    for lab, c in zip(FREQ_BAND_LABELS, agg.f_band_counts):
        A(f"| {lab} | {fmt_int(c)} |")
    A("")

    # amplitude
    A("## Amplitude\n")
    A(f"- Range: **{agg.a_min:.3e} -> {agg.a_max:.3e}** "
      f"(over {fmt_int(agg.a_n)} positive, finite values).")
    A(f"- Geometric mean (typical order of magnitude): **{a_geo:.3e}**.")
    A("")

    # distance
    A("## Distance (distance_au)\n")
    if agg.d_finite:
        A(f"- Over the {fmt_int(agg.d_finite)} rows that have a value: "
          f"**{agg.d_min:.2f} -> {agg.d_max:.2f} AU**, mean **{d_mean:.2f} AU**.")
    else:
        A("- No rows had a finite distance value.")
    A("")

    # spacecraft
    A("## Rows by spacecraft\n")
    A("| Spacecraft | Rows | % |")
    A("|---|---|---|")
    for k, v in sorted(agg.spacecraft.items(), key=lambda x: -x[1]):
        A(f"| {k} | {fmt_int(v)} | {100*v/agg.n_rows:.2f}% |")
    A("")

    # sources
    A("## Top 25 sources\n")
    A("| Source | Rows |")
    A("|---|---|")
    for k, v in sorted(agg.source.items(), key=lambda x: -x[1])[:25]:
        A(f"| {k} | {fmt_int(v)} |")
    n_src = len(agg.source)
    if n_src > 25:
        A(f"\n_({fmt_int(n_src)} distinct sources in total; showing the 25 largest.)_")
    A("")

    # top rows
    A(f"## Top {min(TOP_N, len(agg.top_heap))} highest-scoring rows\n")
    A(f"Exported in full to `foam_top_rows.csv`. The 15 highest:\n")
    top = agg.top_rows_df().head(15)
    A("| foam_score | freq_hz | amplitude | distance_au | spacecraft | source | input_row |")
    A("|---|---|---|---|---|---|---|")
    for _, r in top.iterrows():
        A(f"| {float(r['foam_score']):.6f} | {r['center_freq_hz']} | {r['amplitude']} | "
          f"{r['distance_au']} | {r['spacecraft']} | {r['source']} | {r['row_number']} |")
    A("")

    return "\n".join(L)


def maybe_plot(agg: Aggregator, out_png: Path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False
    centers = 0.5 * (agg.score_hist_edges[:-1] + agg.score_hist_edges[1:])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2))
    ax1.bar(centers, agg.score_hist, width=1.0 / SCORE_BINS, color="#3b6ea5")
    ax1.set_title("foam_score distribution")
    ax1.set_xlabel("foam_score"); ax1.set_ylabel("rows")
    ax2.bar(centers, np.maximum(agg.score_hist, 0), width=1.0 / SCORE_BINS, color="#a5453b")
    ax2.set_yscale("log")
    ax2.set_title("foam_score distribution (log y - shows the tail)")
    ax2.set_xlabel("foam_score"); ax2.set_ylabel("rows (log)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    return True


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Generate a report from the foam scores file.")
    p.add_argument("--csv", "-c", help="Path to *__streaming_scores.csv.gz (or .csv)")
    p.add_argument("--outdir", "-o", default=str(OUT_DIR))
    args = p.parse_args()

    csv_path = args.csv or input("Enter path to the scores file (.csv.gz): ").strip()
    csv_path = csv_path.strip().strip('"').strip("'")
    csv_path = Path(os.path.normpath(os.path.expanduser(csv_path)))
    if not csv_path.exists():
        print(f"ERROR: file not found: {csv_path}")
        sys.exit(2)

    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    comp = sniff_compression(csv_path)
    print(f"[REPORT] Reading: {csv_path}")
    print(f"[REPORT] Detected compression: {comp}")
    print(f"[REPORT] Size on disk: {csv_path.stat().st_size / 1e9:.2f} GB")
    print(f"[REPORT] Streaming in chunks of {CHUNK_ROWS:,} rows (constant memory)...\n")

    agg = Aggregator()
    t0 = time.time()
    try:
        reader = pd.read_csv(csv_path, compression=comp, chunksize=CHUNK_ROWS,
                             low_memory=False)
        last_print = 0
        for chunk in reader:
            agg.update(chunk)
            if agg.n_rows - last_print >= PROGRESS_EVERY:
                last_print = agg.n_rows
                rate = agg.n_rows / max(time.time() - t0, 1e-9)
                print(f"[REPORT] rows={agg.n_rows:,}  rate={rate:,.0f}/s  "
                      f"elapsed={time.time()-t0:.0f}s")
    except Exception as e:
        print(f"[REPORT] Stopped after {agg.n_rows:,} rows due to: {e}")
        print("[REPORT] Writing a partial report from what was read so far.")

    elapsed = time.time() - t0
    if agg.n_rows == 0:
        print("[REPORT] No rows were read - nothing to report.")
        sys.exit(1)

    # write outputs
    report_md = build_report(agg, csv_path, elapsed)
    report_path = out_dir / "foam_report.md"
    report_path.write_text(report_md, encoding="utf-8")

    top_path = out_dir / "foam_top_rows.csv"
    agg.top_rows_df().to_csv(top_path, index=False)

    png_path = out_dir / "foam_score_hist.png"
    drew = maybe_plot(agg, png_path)

    print(f"\n[REPORT] DONE. {agg.n_rows:,} rows in {elapsed:.1f}s")
    print(f"[REPORT] Report : {report_path}")
    print(f"[REPORT] Top rows: {top_path}")
    if drew:
        print(f"[REPORT] Histogram: {png_path}")
    print("\n----- quick console summary -----")
    print(f"  rows scored        : {agg.n_rows:,}")
    print(f"  foam_score min/max : {agg.score_min:.4f} / {agg.score_max:.4f}")
    sm = agg.score_sum / agg.score_n if agg.score_n else float('nan')
    print(f"  foam_score mean    : {sm:.4f}")
    print(f"  rows >= 0.50       : {agg.thresh_counts[0.50]:,}")
    print(f"  rows >= 0.90       : {agg.thresh_counts[0.90]:,}")
    print(f"  highest single row : {agg.score_max:.4f}")


if __name__ == "__main__":
    main()
