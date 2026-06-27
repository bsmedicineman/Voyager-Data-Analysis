#!/usr/bin/env python3
"""
analyze_em_foam.py
==================
Stream CDF field data window-by-window and write one row of EM / spectral
features per window straight to an output CSV. Looking for narrowband,
persistent, drifting spectral lines as candidate spacetime-foam signatures in
magnetometer / electric-field / search-coil data.

Streams input -> output: never holds more than one analysis window (plus a
small read block) in memory, so RAM stays flat no matter how big the archive is.

Per the spec, EVERY row carries:  timestamp (start/end)  |  EM features  |  location.

Usage
-----
  python analyze_em_foam.py run --input  /data/cdf \\
                                --output /data/foam_features.csv

  python analyze_em_foam.py run -i file.cdf -o out.csv \\
        --window 65536 --overlap 0.5 --fmin 1e-3 --fmax 0.5 \\
        --p-fa 1e-3 --q-min 20 --known-lines 0.00417,0.0083

  python analyze_em_foam.py selftest          # verify the detector math

Resume: re-running the same command skips input files already recorded in
<output>.done and appends; a crash loses at most the current window.

IMPORTANT (read the README): a flagged window is a *candidate anomaly*, i.e. a
narrowband feature that survives a colored-noise null test. Calling it "foam"
requires ruling out spacecraft interference, instrument lines, and natural
plasma waves (ion-cyclotron, Langmuir, whistler...). This tool helps you find
and characterize candidates; it cannot by itself attribute them to new physics.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time

import numpy as np

import em_foam_common as C

LOG = C.get_logger("analyze")

CSV_COLUMNS = [
    "file", "dataset", "variable", "component",
    "t_start_utc", "t_end_utc",
    "pos_r", "pos_x", "pos_y", "pos_z",
    "n_samples", "fs_hz",
    "mean_abs", "std", "psd_slope", "spectral_flatness",
    "peak_freq_hz", "peak_ratio", "peak_pvalue", "peak_q", "coherence_s",
    "freq_drift_hz_per_win", "n_lines",
    "candidate", "candidate_reason",
]


# --------------------------------------------------------------------------- #
#  Per-file streaming analysis
# --------------------------------------------------------------------------- #
def analyze_file(path: str, writer: csv.writer, args) -> int:
    """Stream one CDF, emit feature rows. Returns number of rows written."""
    import cdflib
    try:
        cdf = cdflib.CDF(path)
    except Exception as e:
        LOG.warning("skip (open failed) %s: %s", os.path.basename(path), e)
        return 0

    em_vars = C.find_em_variables(cdf)
    if not em_vars:
        LOG.info("no EM variables in %s", os.path.basename(path))
        return 0

    pos_var = C.find_position_variable(cdf)
    dataset = os.path.basename(path).split("_20")[0]   # crude dataset tag
    rows = 0

    for var in em_vars:
        tvar = C.find_time_variable(cdf, near=var)
        if not tvar:
            continue
        nrec = C._record_count(cdf, var)
        window = _fit_window(args.window, nrec)
        if window < args.min_window:
            continue
        step = max(1, int(window * (1.0 - args.overlap)))
        fill = C.fillval(cdf, var)
        reader = C.make_reader(cdf, tvar, var, fill)
        pos_reader = (C.make_reader(cdf, tvar, pos_var, C.fillval(cdf, pos_var))
                      if pos_var else None)

        recent_peak = []        # tiny rolling buffer for drift estimate (O(1) mem)

        for t_ns, data in C.iter_windows(reader, nrec, window, step):
            if np.isnan(data).mean() > args.max_nan:
                continue
            fs = C.median_sample_rate_hz(t_ns)
            if not np.isfinite(fs) or fs <= 0:
                continue

            mag = _magnitude(data)               # rotation-invariant series
            mag = _fill_gaps(mag)
            if mag is None:
                continue

            fmax = min(args.fmax, fs / 2.0 * 0.95)
            if fmax <= args.fmin:
                continue

            nperseg = max(256, window // 8)
            f, pxx, kseg = C.welch_psd(mag, fs, nperseg)
            peaks = C.narrowband_scan(f, pxx, kseg, args.p_fa,
                                      args.fmin, fmax)
            slope = C.powerlaw_slope(f, pxx, args.fmin, fmax)
            flat = C.spectral_flatness(pxx, f, args.fmin, fmax)

            top = peaks[0] if peaks else None
            drift = _drift(recent_peak, top.freq if top else None)

            pos = _sample_position(pos_reader, window) if pos_reader else None
            cand, reason = _candidate(top, args, fmax)

            t0 = np.datetime64(int(t_ns[0]), "ns")
            t1 = np.datetime64(int(t_ns[-1]), "ns")
            writer.writerow([
                os.path.basename(path), dataset, var, "|B|",
                str(t0), str(t1),
                _fmt(pos[0] if pos else None),
                _fmt(pos[1] if pos else None),
                _fmt(pos[2] if pos else None),
                _fmt(pos[3] if pos else None),
                len(mag), round(fs, 6),
                _fmt(float(np.nanmean(np.abs(mag)))),
                _fmt(float(np.nanstd(mag))),
                _fmt(slope), _fmt(flat),
                _fmt(top.freq if top else None),
                _fmt(top.ratio if top else None),
                _fmt(top.p_value if top else None),
                _fmt(top.q if top else None),
                _fmt(top.coherence_s if top else None),
                _fmt(drift), len(peaks),
                int(cand), reason,
            ])
            rows += 1

    return rows


# --------------------------------------------------------------------------- #
#  small helpers
# --------------------------------------------------------------------------- #
def _fit_window(req: int, nrec: int) -> int:
    """Largest power-of-two window <= requested and <= record count."""
    w = min(req, nrec)
    if w < 2:
        return 0
    return 1 << int(np.floor(np.log2(w)))


def _magnitude(data: np.ndarray) -> np.ndarray:
    if data.shape[1] >= 3:
        return np.sqrt(np.nansum(data[:, :3] ** 2, axis=1))
    return data[:, 0]


def _fill_gaps(x: np.ndarray) -> np.ndarray | None:
    """Linear-interpolate isolated NaNs; bail if mostly gaps."""
    good = np.isfinite(x)
    if good.sum() < max(16, len(x) // 2):
        return None
    if good.all():
        return x
    idx = np.arange(len(x))
    return np.interp(idx, idx[good], x[good])


def _drift(buf: list, freq) -> float | None:
    """Drift of dominant line vs previous window (Hz per window)."""
    if freq is None:
        return None
    d = None
    if buf:
        d = freq - buf[-1]
    buf.append(freq)
    if len(buf) > 4:
        buf.pop(0)
    return d


def _sample_position(pos_reader, window):
    """Read position only at the window midpoint -> r, x, y, z (memory-light)."""
    try:
        _, p = pos_reader(0, min(window, 2))
        v = p[0]
        if v.shape[0] >= 3:
            x, y, z = float(v[0]), float(v[1]), float(v[2])
            return (float(np.sqrt(x * x + y * y + z * z)), x, y, z)
        return (float(v[0]), None, None, None)
    except Exception:
        return None


def _candidate(top, args, fmax):
    """Apply the foam-candidate gate and explain the verdict."""
    if top is None:
        return False, "no_line"
    if top.p_value >= args.p_fa:
        return False, "not_significant"
    if not np.isfinite(top.q) or top.q < args.q_min:
        return False, f"low_q({top.q:.1f})"
    for fl in args.known_lines:
        if abs(top.freq - fl) <= args.line_tol * max(fl, 1e-9):
            return False, f"known_line~{fl:g}Hz"
    return True, "candidate"


def _fmt(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return v


# --------------------------------------------------------------------------- #
#  driver
# --------------------------------------------------------------------------- #
def collect_inputs(path: str) -> list[str]:
    if os.path.isfile(path):
        return [path]
    out = []
    for root, _, files in os.walk(path):
        for fn in sorted(files):
            if fn.lower().endswith(".cdf"):
                out.append(os.path.join(root, fn))
    return out


def run(args) -> int:
    inputs = collect_inputs(args.input)
    if not inputs:
        LOG.error("no .cdf files under %s", args.input)
        return 2

    done_log = args.output + ".done"
    done = set()
    if args.resume and os.path.exists(done_log):
        with open(done_log) as fh:
            done = {ln.strip() for ln in fh if ln.strip()}

    new_file = not os.path.exists(args.output) or os.path.getsize(args.output) == 0
    out = open(args.output, "a", newline="")
    writer = csv.writer(out)
    if new_file:
        writer.writerow(CSV_COLUMNS)
        out.flush()

    total_rows = 0
    t0 = time.time()
    try:
        for i, path in enumerate(inputs, 1):
            if path in done:
                continue
            LOG.info("[%d/%d] %s", i, len(inputs), os.path.basename(path))
            n = analyze_file(path, writer, args)
            out.flush()
            os.fsync(out.fileno())
            with open(done_log, "a") as dh:        # mark file complete (resume)
                dh.write(path + "\n")
            total_rows += n
            LOG.info("    -> %d windows (running total %d)", n, total_rows)
    finally:
        out.close()

    dt = time.time() - t0
    cand = _count_candidates(args.output)
    LOG.info("done: %d files, %d feature rows, %d candidate windows in %.1fs",
             len(inputs), total_rows, cand, dt)
    LOG.info("output: %s", args.output)
    if cand:
        LOG.info("NOTE: candidates are anomalies vs a noise null -- vet against "
                 "spacecraft/instrument/plasma lines before any foam claim.")
    return 0


def _count_candidates(csv_path: str) -> int:
    n = 0
    try:
        with open(csv_path) as fh:
            r = csv.DictReader(fh)
            for row in r:
                if row.get("candidate") == "1":
                    n += 1
    except Exception:
        pass
    return n


# --------------------------------------------------------------------------- #
#  self-test: inject a known line into noise and confirm we recover it
# --------------------------------------------------------------------------- #
def selftest(_args) -> int:
    LOG.info("self-test: synthetic narrowband line in colored noise")
    rng = np.random.default_rng(0)
    fs = 1.0                          # 1 Hz sampling
    n = 1 << 16
    t = np.arange(n) / fs
    f0 = 0.1                          # injected line at 0.1 Hz
    noise = np.cumsum(rng.standard_normal(n))        # red-ish noise
    noise -= noise.mean()
    noise /= noise.std()
    signal = 0.6 * np.sin(2 * np.pi * f0 * t)
    x = noise + signal

    f, pxx, k = C.welch_psd(x, fs, nperseg=4096)
    peaks = C.narrowband_scan(f, pxx, k, p_fa=1e-3, fmin=1e-3, fmax=0.49)

    ok = bool(peaks) and abs(peaks[0].freq - f0) < 0.005
    if ok:
        p = peaks[0]
        LOG.info("PASS: recovered f=%.4f Hz (truth %.3f), ratio=%.1f, "
                 "p=%.2e, Q=%.1f, %d total lines",
                 p.freq, f0, p.ratio, p.p_value, p.q, len(peaks))
    else:
        LOG.error("FAIL: did not recover the injected line (got %s)",
                  [round(p.freq, 4) for p in peaks])

    # negative control: pure noise should rarely flag anything at p_fa=1e-3
    fn, pn, kn = C.welch_psd(noise, fs, nperseg=4096)
    false = C.narrowband_scan(fn, pn, kn, p_fa=1e-3, fmin=1e-3, fmax=0.49)
    LOG.info("negative control (pure noise): %d false lines at p_fa=1e-3", len(false))
    return 0 if ok else 1


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Stream CDF EM data and flag narrowband foam-candidate lines.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="analyze a file or directory of CDFs")
    r.add_argument("-i", "--input", required=True,
                   help="input .cdf file OR directory to walk")
    r.add_argument("-o", "--output", required=True, help="output CSV path")
    r.add_argument("--window", type=int, default=65536,
                   help="samples per analysis window (rounded down to power of 2)")
    r.add_argument("--overlap", type=float, default=0.5,
                   help="window overlap fraction [0,1)")
    r.add_argument("--fmin", type=float, default=1e-3, help="min analysis freq (Hz)")
    r.add_argument("--fmax", type=float, default=0.5,
                   help="max analysis freq (Hz); clipped to 0.95*Nyquist")
    r.add_argument("--p-fa", type=float, default=1e-3,
                   help="per-window false-alarm probability for line detection")
    r.add_argument("--q-min", type=float, default=20.0,
                   help="minimum quality factor for a foam candidate")
    r.add_argument("--known-lines", type=_floatlist, default=[],
                   help="comma-separated Hz to veto (spin harmonics, reaction "
                        "wheels, heaters...)")
    r.add_argument("--line-tol", type=float, default=0.02,
                   help="fractional tolerance when matching known lines")
    r.add_argument("--max-nan", type=float, default=0.2,
                   help="skip a window if more than this fraction is fill/NaN")
    r.add_argument("--min-window", type=int, default=256,
                   help="skip variables/files shorter than this many samples")
    r.add_argument("--no-resume", dest="resume", action="store_false",
                   help="reprocess everything instead of skipping <output>.done")
    r.set_defaults(resume=True, func=run)

    s = sub.add_parser("selftest", help="verify the detector on a synthetic signal")
    s.set_defaults(func=selftest)
    return p


def _floatlist(s: str) -> list[float]:
    return [float(x) for x in s.split(",") if x.strip()]


def main() -> int:
    args = build_parser().parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
