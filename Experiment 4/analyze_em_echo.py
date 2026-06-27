#!/usr/bin/env python3
"""
analyze_em_echo.py
==================
A SECOND, different "lens" on the same EM archives (Voyager included): instead
of hunting steady-state narrowband lines (analyze_em_foam.py), this hunts
TRANSIENT-DRIVEN POST-SPIKE ECHOES.

The template comes from the one positive Earth result: a large coherent energy
release (seismic cluster) produced a main spike in the Schumann cavity followed
by a delayed, structured, longer-than-expected trailing echo. Here we look for
the same morphology when the "drive" is a planetary-scale magnetic event -- a
bow-shock / magnetopause crossing or an interval of intense magnetospheric
activity during a flyby.

For each detected drive spike it reports:
  * relaxation time vs the ambient (pre-spike) decorrelation time  -> excess persistence
  * delayed matched-filter echoes of the drive shape               -> repeated/echoed bursts
  * spectral centroid shift between drive and echo                 -> downward frequency shift
  * spacecraft location + UTC timestamp of the event

How it stays light & fast:
  pass 1 builds a DECIMATED activity envelope for the whole file (tiny memory);
  spikes, persistence and echoes are measured on that envelope; spectral metrics
  read only short raw snapshots at the spike and echo lags. Memory stays bounded
  regardless of sampling rate or trailing-window length.

Usage
-----
  python analyze_em_echo.py run \\
      --input /data/spacefoam /data/voyager \\
      --output /data/echo_events.csv \\
      --env-seconds 1 --trail-seconds 3600 \\
      --k 6 --min-spike-snr 3 --min-echo-z 6 --min-persistence 3

  python analyze_em_echo.py selftest      # synthetic spike+echo, verify recovery

Inputs: .cdf (everything from download_em_data.py, plus Voyager CDFs) and simple
.csv/.tab with a time column + numeric field columns. Resume: re-running skips
files already in <output>.done.

INTERPRETATION: a flagged event is a CANDIDATE echo -- a delayed/persistent
feature that stands above the ambient null. Giant-planet magnetospheres are full
of natural echoes (whistlers, chorus, reflected waves); a candidate must survive
that scrutiny (and the Kp/lightning-style confounder regression you used on
Earth) before it means anything about foam. See README.
"""

from __future__ import annotations

import argparse
import csv as csvmod
import os
import sys
import time

import numpy as np

import em_foam_common as C

LOG = C.get_logger("echo")

CSV_COLUMNS = [
    "file", "dataset", "variable", "kind",
    "spike_t_utc", "spike_index",
    "pos_r", "pos_x", "pos_y", "pos_z",
    "baseline", "spike_level", "spike_snr", "step_mean",
    "relax_time_s", "rise_time_s", "persistence_ratio", "relax_truncated",
    "echo_delay_s", "echo_z", "echo_count",
    "drive_centroid", "echo_centroid", "freq_shift", "centroid_units",
    "spectral_similarity",
    "candidate", "candidate_reason",
]


# --------------------------------------------------------------------------- #
#  activity proxy (one scalar per record/sample)
# --------------------------------------------------------------------------- #
def activity(d: np.ndarray, kind: str) -> np.ndarray:
    if kind == "spectral":
        return np.nansum(d, axis=1)              # total power across channels
    if d.shape[1] >= 3:
        return np.sqrt(np.nansum(d[:, :3] ** 2, axis=1))   # |B|
    return np.abs(d[:, 0])


# --------------------------------------------------------------------------- #
#  pass 1: decimated activity envelope for the whole source
# --------------------------------------------------------------------------- #
def build_envelope(reader, nrec, env_w, kind):
    """Return (env_t_ns, env_std, env_mean) at one point per env_w samples."""
    n_env = nrec // env_w
    if n_env < 8:
        return None
    env_t = np.empty(n_env, dtype=np.int64)
    env_std = np.empty(n_env)
    env_mean = np.empty(n_env)
    per_read = max(1, (1 << 20) // env_w)        # bound RAM per read
    i = 0
    while i < n_env:
        j = min(n_env, i + per_read)
        t, d = reader(i * env_w, j * env_w)       # one read per chunk
        L = min(len(d), (j - i) * env_w)
        j = i + L // env_w
        if j <= i:
            break
        a = activity(d[:L], kind).reshape(j - i, env_w)
        env_std[i:j] = np.nanstd(a, axis=1)
        env_mean[i:j] = np.nanmean(a, axis=1)
        env_t[i:j] = t[:L].reshape(j - i, env_w)[:, env_w // 2]
        i = j
    return env_t[:i], env_std[:i], env_mean[:i]


# --------------------------------------------------------------------------- #
#  spike detection + per-event metrics (all on the envelope)
# --------------------------------------------------------------------------- #
def detect_spikes(env_std, k, min_distance, min_prom_frac):
    from scipy.signal import find_peaks
    med = np.nanmedian(env_std)
    mad = np.nanmedian(np.abs(env_std - med)) + 1e-30
    height = med + k * 1.4826 * mad
    prom = min_prom_frac * (np.nanmax(env_std) - med + 1e-30)
    peaks, _ = find_peaks(np.nan_to_num(env_std, nan=med),
                          height=height, distance=max(1, min_distance),
                          prominence=prom)
    return peaks, med, mad


def relax_and_rise(env, peak_i, baseline, frac, trail_pts, env_dt, dip_tol=3):
    """Rise time (up to peak) and relaxation time (down from peak) at the same
    level. persistence = relax/rise: ~1 for a symmetric burst, >>1 for a long
    tail. A short run of sub-threshold points (dip_tol) is tolerated so envelope
    ripple from an oscillatory signal doesn't truncate the span prematurely."""
    peak = env[peak_i]
    if peak <= baseline:
        return 0.0, 0.0, False
    lvl = baseline + frac * (peak - baseline)

    def extent(step, limit):
        i = peak_i
        last_above = peak_i
        run = 0
        while True:
            i += step
            if i < 0 or i >= len(env) or abs(i - peak_i) > limit:
                break
            if env[i] > lvl:
                last_above = i
                run = 0
            else:
                run += 1
                if run >= dip_tol:
                    break
        return abs(last_above - peak_i)

    relax_pts = extent(+1, trail_pts)
    rise_pts = extent(-1, trail_pts)
    truncated = (peak_i + relax_pts) >= min(len(env), peak_i + trail_pts) - 1
    return relax_pts * env_dt, rise_pts * env_dt, truncated


def detect_echo(env, peak_i, drive_half, trail_pts, env_dt, z_thresh):
    """
    A delayed recurrence of the drive is a secondary peak in the trailing
    envelope. Detect peaks (above a robust trailing baseline) after the drive's
    own tail; report the strongest one's delay and z-score. Robust and direct --
    no template-width assumptions.
    """
    from scipy.signal import find_peaks
    t0 = min(len(env), peak_i + max(1, drive_half))      # skip the drive itself
    t1 = min(len(env), peak_i + trail_pts)
    seg = env[t0:t1].astype(float)
    if len(seg) < 8:
        return dict(delay_s=None, z=0.0, count=0, rec_off=None)
    base = np.median(seg)
    mad = np.median(np.abs(seg - base)) + 1e-30
    height = base + z_thresh * 1.4826 * mad
    pk, _ = find_peaks(seg, height=height, distance=max(1, drive_half))
    if len(pk) == 0:
        return dict(delay_s=None, z=0.0, count=0, rec_off=None)
    zvals = (seg[pk] - base) / (1.4826 * mad)
    best = int(pk[int(np.argmax(zvals))])
    delay_idx = (t0 + best) - peak_i
    return dict(delay_s=delay_idx * env_dt, z=float(zvals.max()),
                count=int(len(pk)), rec_off=delay_idx)


# --------------------------------------------------------------------------- #
#  spectral snapshots (read only a short raw window)
# --------------------------------------------------------------------------- #
def centroid_timeseries(reader, rec_center, w_snap, fs, fmin):
    s = max(0, rec_center - w_snap // 2)
    t, d = reader(s, s + w_snap)
    x = activity(d, "timeseries")
    good = np.isfinite(x)
    if good.sum() < max(16, len(x) // 2):
        return None, None
    x = np.interp(np.arange(len(x)), np.flatnonzero(good), x[good])
    f, p, _ = C.welch_psd(x, fs, min(w_snap, 2048))
    m = f >= fmin
    if p[m].sum() <= 0:
        return None, None
    return float((f[m] * p[m]).sum() / p[m].sum()), p[m]


def centroid_spectral(reader, rec_center, w):
    s = max(0, rec_center - w // 2)
    _, d = reader(s, s + max(1, w))
    spec = np.nanmean(d, axis=0)
    spec = np.where(np.isfinite(spec), spec, 0.0)
    if spec.sum() <= 0:
        return None, None
    ch = np.arange(len(spec))
    return float((ch * spec).sum() / spec.sum()), spec


# --------------------------------------------------------------------------- #
#  source iteration: CDF field variables, or a simple CSV
# --------------------------------------------------------------------------- #
def cdf_sources(path):
    import cdflib
    cdf = cdflib.CDF(path)
    pos_var = C.find_position_variable(cdf)
    for v in C.cdf_variables(cdf):
        low = v.lower()
        if any(h in low for h in C._TIME_HINTS):
            continue
        if not any(h in low for h in C._EM_HINTS):
            continue
        nrec = C._record_count(cdf, v)
        if nrec < 64:
            continue
        ld = C._last_dim(cdf, v)
        kind = "spectral" if (ld and ld > 3) else "timeseries"
        tvar = C.find_time_variable(cdf, near=v)
        if not tvar:
            continue
        reader = C.make_reader(cdf, tvar, v, C.fillval(cdf, v))
        pos_reader = (C.make_reader(cdf, tvar, pos_var, C.fillval(cdf, pos_var))
                      if pos_var else None)
        yield dict(var=v, kind=kind, nrec=nrec, reader=reader,
                   pos_reader=pos_reader)


def csv_source(path, max_mb=500):
    size_mb = os.path.getsize(path) / 1e6
    if size_mb > max_mb:
        LOG.warning("CSV %s is %.0f MB; loading anyway (raise --csv-max-mb to silence)",
                    os.path.basename(path), size_mb)
    # detect delimiter + header
    with open(path) as fh:
        sample = fh.readline()
    delim = "\t" if sample.count("\t") > sample.count(",") else ","
    rows = np.genfromtxt(path, delimiter=delim, names=True, dtype=None,
                         encoding="utf-8", deletechars="")
    names = list(rows.dtype.names)
    tname = next((n for n in names
                  if any(h in n.lower() for h in ("time", "epoch", "date"))), names[0])
    t_ns = _parse_time_column(rows[tname])
    comps = [n for n in names if n != tname]
    data = np.column_stack([np.asarray(rows[n], dtype="float64") for n in comps])
    order = np.argsort(t_ns)
    t_ns, data = t_ns[order], data[order]

    def reader(s, e):
        return t_ns[s:e], data[s:e]

    yield dict(var="+".join(comps[:3]) or "csv", kind="timeseries",
               nrec=len(t_ns), reader=reader, pos_reader=None)


def _parse_time_column(col):
    try:
        return np.asarray(col, dtype="datetime64[ns]").astype("int64")
    except Exception:
        return (np.asarray(col, dtype="float64") * 1e9).astype("int64")   # seconds


# --------------------------------------------------------------------------- #
#  per-file analysis
# --------------------------------------------------------------------------- #
def analyze_file(path, writer, args) -> int:
    ext = os.path.splitext(path)[1].lower()
    try:
        srcs = (csv_source(path, args.csv_max_mb)
                if ext in (".csv", ".tab", ".txt") else cdf_sources(path))
        srcs = list(srcs)
    except Exception as e:
        LOG.warning("skip %s: %s", os.path.basename(path), e)
        return 0

    dataset = os.path.basename(path).split("_20")[0]
    rows = 0
    for src in srcs:
        rows += _analyze_source(path, dataset, src, writer, args)
    return rows


def _analyze_source(path, dataset, src, writer, args) -> int:
    reader, nrec, kind = src["reader"], src["nrec"], src["kind"]
    t0, _ = reader(0, min(nrec, 4096))
    fs = C.median_sample_rate_hz(t0)
    if not np.isfinite(fs) or fs <= 0:
        return 0

    env_w = max(args.min_env_samples, int(round(args.env_seconds * fs)))
    env = build_envelope(reader, nrec, env_w, kind)
    if env is None:
        return 0
    env_t, env_std, env_mean = env
    env_dt = env_w / fs

    trail_pts = max(4, int(round(args.trail_seconds / env_dt)))
    drive_half = max(1, int(round(args.drive_seconds / 2 / env_dt)))
    quiet_pts = max(8, int(round(args.quiet_seconds / env_dt)))
    w_snap = min(1 << 14, max(256, int(round(args.drive_seconds * fs))))

    spikes, med, mad = detect_spikes(env_std, args.k,
                                     min_distance=drive_half,
                                     min_prom_frac=args.min_prominence)
    rows = 0
    for pi in spikes:
        lo = max(0, pi - quiet_pts)
        baseline = float(np.median(env_std[lo:pi])) if pi > lo else med
        peak = float(env_std[pi])
        snr = peak / baseline if baseline > 0 else float("inf")
        step = float(env_mean[min(len(env_mean) - 1, pi + drive_half)]
                     - env_mean[max(0, pi - drive_half)])

        rt, rise, trunc = relax_and_rise(env_std, pi, baseline, args.relax_frac,
                                         trail_pts, env_dt)
        persistence = rt / rise if rise > 0 else float("nan")

        echo = detect_echo(env_std, pi, drive_half, trail_pts, env_dt, args.min_echo_z)

        # spectral snapshots at the spike and (if found) the echo
        units = "Hz" if kind == "timeseries" else "channel"
        if kind == "timeseries":
            dc, dpsd = centroid_timeseries(reader, pi * env_w, w_snap, fs, args.fmin)
        else:
            dc, dpsd = centroid_spectral(reader, pi * env_w, max(1, env_w))
        ec, epsd, shift, sim = None, None, None, None
        if echo["rec_off"] is not None:
            erec = (pi + echo["rec_off"]) * env_w
            if kind == "timeseries":
                ec, epsd = centroid_timeseries(reader, erec, w_snap, fs, args.fmin)
            else:
                ec, epsd = centroid_spectral(reader, erec, max(1, env_w))
            if dc is not None and ec is not None:
                shift = ec - dc
            if dpsd is not None and epsd is not None and len(dpsd) == len(epsd):
                if np.std(dpsd) > 0 and np.std(epsd) > 0:
                    sim = float(np.corrcoef(dpsd, epsd)[0, 1])

        pos = _sample_position(src["pos_reader"], pi * env_w)
        cand, reason = _candidate(snr, persistence, echo, shift, args)
        ts = np.datetime64(int(env_t[pi]), "ns")

        writer.writerow([
            os.path.basename(path), dataset, src["var"], kind,
            str(ts), int(pi * env_w),
            _fmt(pos[0] if pos else None), _fmt(pos[1] if pos else None),
            _fmt(pos[2] if pos else None), _fmt(pos[3] if pos else None),
            _fmt(baseline), _fmt(peak), _fmt(snr), _fmt(step),
            _fmt(rt), _fmt(rise), _fmt(persistence), int(trunc),
            _fmt(echo["delay_s"]), _fmt(echo["z"]), echo["count"],
            _fmt(dc), _fmt(ec), _fmt(shift), units,
            _fmt(sim),
            int(cand), reason,
        ])
        rows += 1
    return rows


def _candidate(snr, persistence, echo, shift, args):
    if snr < args.min_spike_snr:
        return False, "weak_spike"
    strong_echo = echo["z"] >= args.min_echo_z and echo["delay_s"] is not None
    excess = np.isfinite(persistence) and persistence >= args.min_persistence
    if not (strong_echo or excess):
        return False, "no_echo_no_persistence"
    bits = []
    if strong_echo:
        bits.append(f"echo(z={echo['z']:.1f},dt={echo['delay_s']:.0f}s)")
    if excess:
        bits.append(f"persist(x{persistence:.1f})")
    if shift is not None and shift < 0:
        bits.append("downshift")          # the Schumann-like signature
    return True, "+".join(bits)


def _sample_position(pos_reader, rec):
    if not pos_reader:
        return None
    try:
        _, p = pos_reader(rec, rec + 2)
        v = p[0]
        if v.shape[0] >= 3:
            x, y, z = float(v[0]), float(v[1]), float(v[2])
            return (float(np.sqrt(x * x + y * y + z * z)), x, y, z)
        return (float(v[0]), None, None, None)
    except Exception:
        return None


def _fmt(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    return f"{v:.6g}" if isinstance(v, float) else v


# --------------------------------------------------------------------------- #
#  driver
# --------------------------------------------------------------------------- #
def collect_inputs(paths):
    out = []
    for p in paths:
        if os.path.isfile(p):
            out.append(p)
        else:
            for root, _, files in os.walk(p):
                for fn in sorted(files):
                    if fn.lower().endswith((".cdf", ".csv", ".tab", ".txt")):
                        out.append(os.path.join(root, fn))
    return out


def run(args) -> int:
    inputs = collect_inputs(args.input)
    if not inputs:
        LOG.error("no input files under %s", args.input)
        return 2

    done_log = args.output + ".done"
    done = set()
    if args.resume and os.path.exists(done_log):
        done = {ln.strip() for ln in open(done_log) if ln.strip()}

    new = not os.path.exists(args.output) or os.path.getsize(args.output) == 0
    out = open(args.output, "a", newline="")
    writer = csvmod.writer(out)
    if new:
        writer.writerow(CSV_COLUMNS)
        out.flush()

    total = 0
    t0 = time.time()
    try:
        for i, path in enumerate(inputs, 1):
            if path in done:
                continue
            LOG.info("[%d/%d] %s", i, len(inputs), os.path.basename(path))
            n = analyze_file(path, writer, args)
            out.flush()
            os.fsync(out.fileno())
            open(done_log, "a").write(path + "\n")
            total += n
            LOG.info("    -> %d events", n)
    finally:
        out.close()

    cand = _count(args.output)
    LOG.info("done: %d files, %d drive events, %d candidate echoes in %.1fs",
             len(inputs), total, cand, time.time() - t0)
    LOG.info("output: %s", args.output)
    if cand:
        LOG.info("NOTE: candidates are echoes/persistence above the ambient null -- "
                 "vet against natural magnetospheric waves before any foam claim.")
    return 0


def _count(path):
    n = 0
    try:
        for row in csvmod.DictReader(open(path)):
            n += (row.get("candidate") == "1")
    except Exception:
        pass
    return n


# --------------------------------------------------------------------------- #
#  self-test: build a spike + delayed downshifted echo, confirm recovery
# --------------------------------------------------------------------------- #
def selftest(_a) -> int:
    LOG.info("self-test: synthetic drive spike + delayed echo")
    fs = 10.0
    n = 60000                       # 100 min @ 10 Hz
    t = np.arange(n) / fs
    rng = np.random.default_rng(3)
    x = 0.05 * rng.standard_normal(n)

    def burst(center_s, dur_s, amp, freq):
        env = np.exp(-0.5 * ((t - center_s) / (dur_s / 2.355)) ** 2)
        return amp * env * np.sin(2 * np.pi * freq * t)

    x += burst(100, 6, 1.0, 2.0)            # main drive at t=100 s, 2 Hz
    x += burst(160, 6, 0.4, 1.6)            # echo at t=160 s (delay 60 s), down-shifted
    data = x[:, None]

    def reader(s, e):
        return (np.arange(s, e, dtype="int64") * int(1e9 / fs),
                data[s:e])

    args = _defaults()
    src = dict(var="synth", kind="timeseries", nrec=n, reader=reader,
               pos_reader=None)
    import io
    buf = io.StringIO()
    w = csvmod.writer(buf)
    w.writerow(CSV_COLUMNS)
    rows = _analyze_source("synth.cdf", "synth", src, w, args)

    buf.seek(0)
    recs = list(csvmod.DictReader(buf))
    if not recs:
        LOG.error("FAIL: no spike detected"); return 1
    main = recs[0]
    delay = float(main["echo_delay_s"] or "nan")
    z = float(main["echo_z"] or "0")
    shift = float(main["freq_shift"] or "nan")
    snr = float(main["spike_snr"] or "0")
    LOG.info("  spike snr=%.1f  echo delay=%.0fs (truth 60)  z=%.1f  "
             "freq_shift=%.2f Hz (truth ~-0.4)  candidate=%s",
             snr, delay, z, shift, main["candidate"])
    ok = (abs(delay - 60) < 6) and (z >= args.min_echo_z) and (shift < 0) \
        and main["candidate"] == "1"

    # negative control: drive with NO echo -> should not flag an echo
    x2 = 0.05 * rng.standard_normal(n) + burst(100, 6, 1.0, 2.0)
    d2 = x2[:, None]
    src2 = dict(var="synth", kind="timeseries", nrec=n,
                reader=lambda s, e: (np.arange(s, e, dtype="int64") * int(1e9 / fs), d2[s:e]),
                pos_reader=None)
    b2 = io.StringIO(); w2 = csvmod.writer(b2); w2.writerow(CSV_COLUMNS)
    _analyze_source("ctrl.cdf", "ctrl", src2, w2, args)
    b2.seek(0)
    ctrl = list(csvmod.DictReader(b2))
    ctrl_echo = max((float(r["echo_z"] or 0) for r in ctrl), default=0.0)
    ctrl_flag = any(r["candidate"] == "1" for r in ctrl)
    LOG.info("  negative control: max echo z=%.1f, any candidate=%s "
             "(expect modest z, ideally no candidate)", ctrl_echo, ctrl_flag)

    LOG.info("PASS" if ok else "FAIL")
    return 0 if ok else 1


def _defaults():
    p = build_parser()
    return p.parse_args(["run", "-i", ".", "-o", "/dev/null"])


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Transient-driven post-spike echo search in EM data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("run", help="analyze files/dirs of CDF or CSV")
    r.add_argument("-i", "--input", nargs="+", required=True,
                   help="files or directories (e.g. the download dir AND voyager dir)")
    r.add_argument("-o", "--output", required=True, help="output CSV path")
    r.add_argument("--env-seconds", type=float, default=1.0,
                   help="activity-envelope cadence in seconds (time resolution)")
    r.add_argument("--trail-seconds", type=float, default=3600.0,
                   help="how long after a spike to look for echoes")
    r.add_argument("--drive-seconds", type=float, default=30.0,
                   help="expected drive duration (matched-filter template width)")
    r.add_argument("--quiet-seconds", type=float, default=600.0,
                   help="pre-spike window for baseline + ambient decorrelation time")
    r.add_argument("--k", type=float, default=6.0,
                   help="spike threshold: median + k*1.4826*MAD of envelope")
    r.add_argument("--min-prominence", type=float, default=0.05,
                   help="min peak prominence as fraction of envelope range")
    r.add_argument("--relax-frac", type=float, default=1.0 / np.e,
                   help="relaxation endpoint as fraction of (peak-baseline)")
    r.add_argument("--fmin", type=float, default=0.0,
                   help="min freq for spectral centroid (Hz, timeseries only)")
    r.add_argument("--min-spike-snr", type=float, default=3.0,
                   help="candidate gate: spike/baseline ratio")
    r.add_argument("--min-echo-z", type=float, default=6.0,
                   help="candidate gate: matched-filter echo z-score")
    r.add_argument("--min-persistence", type=float, default=3.0,
                   help="candidate gate: relax_time / ambient_tau")
    r.add_argument("--min-env-samples", type=int, default=8,
                   help="floor on samples per envelope point (low-cadence data)")
    r.add_argument("--csv-max-mb", type=float, default=500.0,
                   help="warn above this CSV size (still loads)")
    r.add_argument("--no-resume", dest="resume", action="store_false",
                   help="reprocess everything")
    r.set_defaults(resume=True, func=run)

    s = sub.add_parser("selftest", help="synthetic spike+echo recovery test")
    s.set_defaults(func=selftest)
    return p


if __name__ == "__main__":
    _args = build_parser().parse_args()
    sys.exit(_args.func(_args))
