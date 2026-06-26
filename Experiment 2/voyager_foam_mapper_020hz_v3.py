#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
voyager_foam_mapper_020hz_v3.py  —  Experiment 2, rewritten
============================================================
Maps the 0-20 Hz "foam-oscillation" score over Voyager plasma-wave / magnetometer
telemetry. This is a rewrite of voyager_foam_mapper_020hz_version2.py with three
deliberate changes:

  1. INPUT and OUTPUT are now command-line options (-i / --input, -o / --output).
  2. NASA CDF files are read with `cdflib` (the previous version used xarray, which
     cannot open NASA CDF, so every real Voyager .cdf was silently skipped).
  3. There is NO silent synthetic fallback. If you point it at real data and nothing
     usable is found, it stops with an error. Synthetic test data is generated ONLY
     when you explicitly pass --selftest, and it is clearly labelled as synthetic in
     every output.

The frequency used for the resonance term is taken from a real frequency column in
the data (e.g. the PWS spectrum-analyzer channel frequencies). The previous version
set frequency_hz equal to the band-pass-filtered amplitude, which was not a frequency
at all; that has been corrected.

Memory: input is streamed in two passes (scan for global amplitude bounds, then score),
so arbitrarily large inputs use roughly constant memory.

USAGE
-----
  # Real data (a folder of CDF/CSV files, or a single file):
  python3 voyager_foam_mapper_020hz_v3.py --input  /path/to/raw_data \
                                           --output /path/to/analysis_out

  # Explicit synthetic self-test (clearly labelled, no real data needed):
  python3 voyager_foam_mapper_020hz_v3.py --selftest --output ./selftest_out

OPTIONS
-------
  -i, --input        Path to a data file OR a directory of data files.
  -o, --output       Output directory (created if missing). Default: ./voyager_output
  -v, --voyager      Spacecraft number, 1 or 2 (sets launch date). Default: 1
      --fs           Sampling rate (Hz) used by the 0-20 Hz band-pass. Default: 1000
      --min-freq     Band-pass low edge (Hz). Default: 0.0
      --max-freq     Band-pass high edge (Hz). Default: 20.0
      --chunk-size   Rows per streaming chunk for CSV/TAB. Default: 500000
      --max-rows     Optional hard cap on total rows processed (0 = no cap). Default: 0
      --selftest     Generate labelled synthetic data instead of reading input.
"""

import argparse
import csv as csvlib
import glob
import os
import shutil
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

warnings.filterwarnings("ignore")

# Optional specialized readers
try:
    import cdflib
except Exception:
    cdflib = None
try:
    import xarray as xr
except Exception:
    xr = None

# Column-name patterns (lowercased substring match)
TIME_PATTERNS = ["epoch", "timestamp", "time", "t_start", "date", "datetime", "ut", "scet"]
AMP_PATTERNS  = ["amplitude", "mag_amplitude", "pws_amplitude", "e_field", "b_field",
                 "field", "efield", "bfield", "intensity", "power", "spectral",
                 "mag_total", "pws_total", "value", "amp", "mag", "pws"]
FREQ_PATTERNS = ["center_freq", "frequency_hz", "frequency", "freq_hz", "freq", "channel_freq"]
DIST_PATTERNS = ["distance_au", "distance", "radial", "range_au", "helio", "r_au"]

LAUNCH = {1: pd.Timestamp("1977-09-05"), 2: pd.Timestamp("1977-08-20")}
AU_KM = 149_597_870.7


# ----------------------------------------------------------------------------
class FoamMapper:
    def __init__(self, fs=1000.0, min_freq=0.0, max_freq=20.0, voyager=1):
        self.fs = float(fs)
        self.min_freq = float(min_freq)
        self.max_freq = float(max_freq)
        self.voyager = int(voyager)
        self.launch = LAUNCH.get(self.voyager, LAUNCH[1])

    # ---- band-pass (faithful to the original) ----
    def bandpass(self, data):
        nyq = 0.5 * self.fs
        if nyq <= 0:
            return data
        low = self.min_freq / nyq
        high = self.max_freq / nyq
        try:
            if low <= 0 and 0 < high < 1.0:
                b, a = butter(4, high, btype="low")
            elif high >= 1.0 and 0 < low < 1.0:
                b, a = butter(4, low, btype="high")
            else:
                low = max(low, 1e-6)
                high = min(max(high, low + 1e-6), 0.999999)
                b, a = butter(4, [low, high], btype="band")
            return filtfilt(b, a, data)
        except Exception:
            return data  # graceful: leave unfiltered

    @staticmethod
    def _find_col(columns, patterns):
        cols = [str(c) for c in columns]
        low = [c.lower() for c in cols]
        for pat in patterns:
            for i, c in enumerate(low):
                if pat in c:
                    return cols[i]
        return None

    # ---- loaders -----------------------------------------------------------
    def _detect_header_row(self, path, delimiter):
        """Skip a preamble / acknowledgement block: find the first line whose field
        count matches the rows below it and whose tokens are not all numeric."""
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                lines = [next(f) for _ in range(300)]
        except StopIteration:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.readlines()
        except Exception:
            return 0
        def nfields(s):
            return len(s.rstrip("\n").split(delimiter))
        def mostly_numeric(s):
            toks = [t.strip() for t in s.rstrip("\n").split(delimiter) if t.strip() != ""]
            if not toks:
                return True
            num = 0
            for t in toks:
                try:
                    float(t); num += 1
                except ValueError:
                    pass
            return num / len(toks) > 0.6
        for i in range(len(lines) - 2):
            if lines[i].lstrip().startswith("#"):
                continue
            nf = nfields(lines[i])
            if nf < 2:
                continue
            # header candidate: same field count as next two lines, not mostly numeric,
            # and the lines below ARE mostly numeric (i.e. data)
            if (nfields(lines[i + 1]) == nf and nfields(lines[i + 2]) == nf
                    and not mostly_numeric(lines[i])
                    and (mostly_numeric(lines[i + 1]) or mostly_numeric(lines[i + 2]))):
                return i
        return 0

    def load_csv(self, path):
        # sniff delimiter
        delimiter = ","
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                sample = f.read(8192)
            delimiter = csvlib.Sniffer().sniff(sample, delimiters=",;\t|").delimiter
        except Exception:
            for d in [",", "\t", ";", "|"]:
                if d in sample:
                    delimiter = d
                    break
        skip = self._detect_header_row(path, delimiter)
        frames = []
        try:
            for chunk in pd.read_csv(path, sep=delimiter, skiprows=skip, comment="#",
                                     chunksize=self.chunk_size, low_memory=False,
                                     on_bad_lines="skip", engine="c"):
                frames.append(chunk)
        except Exception:
            try:
                df = pd.read_csv(path, sep=delimiter, skiprows=skip, engine="python",
                                 on_bad_lines="skip")
                frames = [df]
            except Exception as e:
                print(f"   [WARN] could not parse {os.path.basename(path)}: {e}")
                return None
        if not frames:
            return None
        df = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
        return df if not df.empty else None

    def load_cdf(self, path):
        if cdflib is None:
            print(f"   [WARN] cdflib not installed — cannot read CDF {os.path.basename(path)} "
                  f"(pip install cdflib)")
            return None
        try:
            cdf = cdflib.CDF(path)
            info = cdf.cdf_info()
            zvars = list(getattr(info, "zVariables", []) or []) + list(getattr(info, "rVariables", []) or [])
            if not zvars:
                return None
            # time variable
            tname = self._find_col(zvars, TIME_PATTERNS)
            if tname is None:
                return None
            epoch = cdf.varget(tname)
            # Convert by the variable's DECLARED CDF epoch type, not by array dtype.
            # cdflib.to_datetime auto-detects from dtype (float->CDF_EPOCH, int->TT2000),
            # which can misfire on real files; force the dtype to match the declared type.
            try:
                dtcode = getattr(cdf.varinq(tname), "Data_Type", None)
            except Exception:
                dtcode = None
            try:
                if dtcode == 31:        # CDF_EPOCH (float ms since year 0)
                    times = pd.to_datetime(cdflib.cdfepoch.to_datetime(np.asarray(epoch, dtype=np.float64)))
                elif dtcode == 33:      # CDF_TIME_TT2000 (int64 ns since J2000)
                    times = pd.to_datetime(cdflib.cdfepoch.to_datetime(np.asarray(epoch, dtype=np.int64)))
                else:                   # CDF_EPOCH16 or unknown — let cdflib decide
                    times = pd.to_datetime(cdflib.cdfepoch.to_datetime(epoch))
            except Exception:
                times = pd.to_datetime(np.asarray(epoch), errors="coerce")
            ntime = len(times)
            # frequency variable (channel centres)
            fname = self._find_col([v for v in zvars if v != tname], FREQ_PATTERNS)
            freqs = None
            if fname is not None:
                try:
                    freqs = np.asarray(cdf.varget(fname), dtype=float).ravel()
                except Exception:
                    freqs = None
            # primary data variable = largest numeric var that is not time/freq
            data_name, data_arr = None, None
            for v in zvars:
                if v in (tname, fname):
                    continue
                try:
                    arr = np.asarray(cdf.varget(v))
                except Exception:
                    continue
                if arr is None or arr.dtype.kind not in "fiu":
                    continue
                if data_arr is None or arr.size > data_arr.size:
                    data_name, data_arr = v, arr
            if data_arr is None:
                return None
            # distance variable (optional)
            dname = self._find_col([v for v in zvars if v not in (tname, fname, data_name)], DIST_PATTERNS)
            dist = None
            if dname is not None:
                try:
                    dist = np.asarray(cdf.varget(dname), dtype=float).ravel()
                except Exception:
                    dist = None
            # reshape to long format
            if data_arr.ndim == 2 and freqs is not None and data_arr.shape[0] == ntime \
                    and data_arr.shape[1] == len(freqs):
                # [time x channel] -> one row per (time, channel)
                nt, nf = data_arr.shape
                df = pd.DataFrame({
                    "timestamp": np.repeat(times.values, nf),
                    "amplitude": data_arr.reshape(-1).astype(float),
                    "frequency_hz": np.tile(freqs, nt),
                })
                if dist is not None and len(dist) == nt:
                    df["distance_au"] = np.repeat(dist, nf)
            else:
                vals = data_arr.reshape(ntime, -1)[:, 0] if data_arr.ndim > 1 else data_arr.ravel()
                n = min(ntime, len(vals))
                df = pd.DataFrame({
                    "timestamp": times.values[:n],
                    "amplitude": np.asarray(vals[:n], dtype=float),
                    "frequency_hz": (freqs[0] if (freqs is not None and freqs.size == 1) else np.nan),
                })
                if dist is not None and len(dist) >= n:
                    df["distance_au"] = dist[:n]
            df["source_file"] = os.path.basename(path)
            print(f"   [INFO] CDF {os.path.basename(path)}: vars time='{tname}', "
                  f"data='{data_name}', freq='{fname}' -> {len(df):,} rows")
            return df if not df.empty else None
        except Exception as e:
            print(f"   [WARN] CDF load failed {os.path.basename(path)}: {e}")
            return None

    def load_nc(self, path):
        if xr is None:
            return None
        try:
            ds = xr.open_dataset(path)
            df = ds.to_dataframe().reset_index()
            return df if not df.empty else None
        except Exception:
            return None

    # ---- normalize an arbitrary frame to the columns we need ----
    def normalize(self, df, source):
        df.columns = [str(c).strip().lower().replace(" ", "_").replace(".", "_") for c in df.columns]
        out = pd.DataFrame()

        # amplitude (required)
        amp_col = self._find_col(df.columns, AMP_PATTERNS)
        if amp_col is None:
            num = df.select_dtypes(include=[np.number]).columns
            amp_col = num[0] if len(num) else None
        if amp_col is None:
            return None
        amp = pd.to_numeric(df[amp_col], errors="coerce")
        valid = np.isfinite(amp.values)
        if valid.sum() < 10:
            return None
        out["amplitude"] = amp.values

        # band-pass filtered copy (the 0-20 Hz focus)
        af = np.full(len(amp), np.nan)
        try:
            af[valid] = self.bandpass(amp.values[valid])
        except Exception:
            af[valid] = amp.values[valid]
        out["amplitude_filtered"] = af

        # frequency (real column only; NaN if absent — never faked)
        if "frequency_hz" in df.columns:
            out["frequency_hz"] = pd.to_numeric(df["frequency_hz"], errors="coerce").values
        else:
            fcol = self._find_col(df.columns, FREQ_PATTERNS)
            out["frequency_hz"] = pd.to_numeric(df[fcol], errors="coerce").values if fcol else np.nan

        # timestamp
        tcol = self._find_col(df.columns, TIME_PATTERNS)
        if tcol is not None:
            out["timestamp"] = pd.to_datetime(df[tcol], errors="coerce", utc=True).dt.tz_localize(None)
        else:
            out["timestamp"] = pd.NaT

        # distance
        dcol = self._find_col(df.columns, DIST_PATTERNS)
        out["distance_au"] = pd.to_numeric(df[dcol], errors="coerce").values if dcol else np.nan

        out["source_file"] = source
        return out

    # ---- foam score (faithful formula; global amplitude bounds) ----
    def foam_score(self, df, amp_min, amp_max):
        n = len(df)
        amp = pd.to_numeric(df["amplitude"], errors="coerce").fillna(0.0).astype(float).values
        if amp_max > amp_min:
            amp_norm = np.clip((amp - amp_min) / (amp_max - amp_min), 0.0, 1.0)
        else:
            amp_norm = np.ones(n)
        freq = pd.to_numeric(df["frequency_hz"], errors="coerce").values
        rw = np.zeros(n)
        finite_f = np.isfinite(freq)
        for fp in (8.0, 20.0):
            bw = max(0.5, fp * 0.1)
            d = np.where(finite_f, np.abs(freq - fp), np.inf)
            rw += np.where(finite_f, np.exp(-(np.clip(d, 0, 1e6) ** 2) / (2 * bw ** 2)), 0.0)
        rw = np.minimum(rw / 2.0, 1.0)

        days = ((df["timestamp"] - self.launch).dt.total_seconds() / 86400.0).values
        days = np.where(np.isfinite(days), days, 0.0)
        solar = 0.5 + 0.5 * np.sin(2 * np.pi * days / (365.25 * 11.0))

        dist = pd.to_numeric(df["distance_au"], errors="coerce").values
        if np.isfinite(dist).any():
            dmax = np.nanmax(dist) if np.nanmax(dist) > 0 else 1.0
            dnorm = dist / dmax
            deff = np.where(np.isfinite(dist), 0.7 + 0.3 * np.exp(-((dnorm - 0.5) ** 2) / 0.1), 1.0)
        else:
            deff = np.ones(n)

        doy = df["timestamp"].dt.dayofyear.values.astype(float)
        doy = np.where(np.isfinite(doy), doy, 1.0)
        annual = 0.8 + 0.2 * np.sin(2 * np.pi * doy / 365.25)

        coh = np.clip(solar * deff * annual, 0.1, 0.9)
        mismatch = 1.0 / (rw + 0.01)
        score = np.clip(coh * amp_norm / (mismatch + 0.1), 0.0, 1.0)
        return score

    @staticmethod
    def classify_bands(freq_series):
        f = pd.to_numeric(freq_series, errors="coerce").fillna(-9999).values
        bands = np.full(f.shape, "Other", dtype=object)
        bands[(f >= 0) & (f < 1)] = "ULF (0-1 Hz)"
        bands[(f >= 1) & (f < 3)] = "VLF (1-3 Hz)"
        bands[(f >= 3) & (f < 8)] = "LF (3-8 Hz)"
        bands[(f >= 8.5) & (f < 15)] = "Mid (8.5-15 Hz)"
        bands[(f >= 15) & (f <= 20)] = "High (15-20 Hz)"
        bands[(f >= 7.5) & (f < 8.5)] = "8 Hz Resonance"
        bands[(f >= 19.5) & (f < 20.5)] = "20 Hz Resonance"
        return bands

    # ---- synthetic self-test (explicit only, clearly labelled) ----
    def synthetic(self, n=1_000_000):
        print("[SELFTEST] Generating labelled SYNTHETIC data (not Voyager telemetry).")
        t = pd.date_range(start=self.launch, periods=n, freq="s")
        tt = np.arange(n) / max(self.fs, 1.0)
        mag = np.sin(2 * np.pi * 8.0 * tt) + 0.3 * np.sin(2 * np.pi * 20 * tt) + 0.1 * np.random.normal(0, 1, n)
        df = pd.DataFrame({
            "timestamp": t,
            "amplitude": mag,
            "amplitude_filtered": mag,
            "frequency_hz": np.tile([8.0, 20.0, 17.8, 31.1], int(np.ceil(n / 4)))[:n],
            "distance_au": 1.0 + 3.42 * (np.arange(n) / float(n)),
            "source_file": "SYNTHETIC_SELFTEST",
        })
        return df


# ----------------------------------------------------------------------------
def collect_files(input_path):
    p = Path(input_path)
    if p.is_file():
        return [str(p)]
    if p.is_dir():
        exts = ["*.csv", "*.CSV", "*.txt", "*.TXT", "*.dat", "*.DAT", "*.tab", "*.TAB",
                "*.tsv", "*.TSV", "*.cdf", "*.CDF", "*.nc", "*.NC"]
        files = []
        for e in exts:
            files.extend(glob.glob(str(p / e)))
        return sorted(set(files))
    return []


def load_one(mapper, path):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".cdf",):
        raw = mapper.load_cdf(path)
    elif ext in (".nc",):
        raw = mapper.load_nc(path)
    else:
        raw = mapper.load_csv(path)
    if raw is None:
        return None
    if "source_file" in raw.columns and ext == ".cdf":
        # CDF loader already produced normalized columns
        norm = mapper.normalize(raw, os.path.basename(path)) if "amplitude" not in raw.columns else raw
        if "amplitude_filtered" not in norm.columns:
            norm = mapper.normalize(norm, os.path.basename(path))
        return norm
    return mapper.normalize(raw, os.path.basename(path))


def run(args):
    mapper = FoamMapper(fs=args.fs, min_freq=args.min_freq, max_freq=args.max_freq, voyager=args.voyager)
    mapper.chunk_size = int(args.chunk_size)
    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / "_tmp_chunks"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    cap = int(args.max_rows) if args.max_rows else None

    # -------- acquire normalized frames (real or synthetic) --------
    sources, n_rows = [], 0
    amp_min, amp_max = np.inf, -np.inf
    t_min, t_max = None, None
    d_min, d_max = np.inf, -np.inf
    chunk_paths = []

    def stash(norm, idx):
        nonlocal n_rows, amp_min, amp_max, t_min, t_max, d_min, d_max
        a = pd.to_numeric(norm["amplitude"], errors="coerce").values
        a = a[np.isfinite(a)]
        if a.size:
            amp_min = min(amp_min, float(a.min())); amp_max = max(amp_max, float(a.max()))
        ts = pd.to_datetime(norm["timestamp"], errors="coerce")
        if ts.notna().any():
            lo, hi = ts.min(), ts.max()
            t_min = lo if t_min is None else min(t_min, lo)
            t_max = hi if t_max is None else max(t_max, hi)
        dd = pd.to_numeric(norm["distance_au"], errors="coerce").values
        dd = dd[np.isfinite(dd)]
        if dd.size:
            d_min = min(d_min, float(dd.min())); d_max = max(d_max, float(dd.max()))
        path = tmp_dir / f"chunk_{idx:05d}.pkl"
        norm.to_pickle(path)
        chunk_paths.append(path)
        n_rows += len(norm)

    if args.selftest:
        norm = mapper.synthetic(cap if cap else 1_000_000)
        sources = ["SYNTHETIC_SELFTEST"]
        stash(norm, 0)
    else:
        files = collect_files(args.input)
        if not files:
            shutil.rmtree(tmp_dir, ignore_errors=True)
            sys.exit(f"[ERROR] No data files found at: {args.input}\n"
                     f"        Give --input a CSV/TAB/CDF file or a folder of them, "
                     f"or use --selftest for synthetic data.")
        print(f"[LOAD] {len(files)} file(s) at {args.input}")
        idx = 0
        for fp in files:
            print(f"  [{idx+1}/{len(files)}] {os.path.basename(fp)}")
            norm = load_one(mapper, fp)
            if norm is None or norm.empty:
                print("        (no usable rows — skipped)")
                continue
            if cap and n_rows + len(norm) > cap:
                norm = norm.iloc[: max(0, cap - n_rows)]
            sources.append(os.path.basename(fp))
            stash(norm, idx); idx += 1
            if cap and n_rows >= cap:
                print(f"  [INFO] reached --max-rows cap ({cap:,}); stopping.")
                break

        if n_rows == 0:
            shutil.rmtree(tmp_dir, ignore_errors=True)
            sys.exit("[ERROR] Files were found but none contained usable amplitude data.\n"
                     "        Nothing was written. (No synthetic substitution — use --selftest "
                     "if you actually want synthetic data.)")

    # -------- score + write (stream the stashed chunks) --------
    scored_path = out_dir / f"voyager_{args.voyager}_foam_scored.csv"
    if scored_path.exists():
        scored_path.unlink()
    score_max, score_sum, score_cnt = -np.inf, 0.0, 0
    ge = {t: 0 for t in (0.01, 0.05, 0.1, 0.2)}
    wrote_header = False
    for path in chunk_paths:
        norm = pd.read_pickle(path)
        norm["days_since_launch"] = (norm["timestamp"] - mapper.launch).dt.total_seconds() / 86400.0
        norm["distance_km"] = pd.to_numeric(norm["distance_au"], errors="coerce") * AU_KM
        norm["year"] = pd.to_datetime(norm["timestamp"], errors="coerce").dt.year
        norm["foam_oscillation_score"] = mapper.foam_score(norm, amp_min, amp_max)
        norm["frequency_band"] = mapper.classify_bands(norm["frequency_hz"])
        s = norm["foam_oscillation_score"].values
        score_max = max(score_max, float(np.nanmax(s))) if len(s) else score_max
        score_sum += float(np.nansum(s)); score_cnt += int(np.isfinite(s).sum())
        for t in ge:
            ge[t] += int((s >= t).sum())
        cols = ["timestamp", "source_file", "frequency_hz", "frequency_band", "amplitude",
                "amplitude_filtered", "distance_au", "distance_km", "days_since_launch",
                "year", "foam_oscillation_score"]
        norm[cols].to_csv(scored_path, mode="a", header=not wrote_header, index=False)
        wrote_header = True
    shutil.rmtree(tmp_dir, ignore_errors=True)

    # -------- honest summary --------
    is_syn = bool(args.selftest)
    summary_path = out_dir / f"voyager_{args.voyager}_data_summary.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(f"Voyager {args.voyager} data summary\n")
        f.write("=" * 60 + "\n")
        if is_syn:
            f.write("DATA SOURCE   : *** SYNTHETIC SELF-TEST (NOT Voyager telemetry) ***\n")
        else:
            f.write("DATA SOURCE   : REAL input files\n")
            f.write(f"Files read    : {len(sources)}\n")
            for s in sources[:25]:
                f.write(f"                - {s}\n")
            if len(sources) > 25:
                f.write(f"                ... (+{len(sources)-25} more)\n")
        f.write(f"Rows          : {n_rows:,}\n")
        if t_min is not None:
            f.write(f"Time range    : {t_min} -- {t_max}\n")
        else:
            f.write("Time range    : (no parseable timestamps)\n")
        if np.isfinite(d_min):
            f.write(f"Distance      : {d_min:.3f} -- {d_max:.3f} AU\n")
        else:
            f.write("Distance      : (no distance column in input)\n")
        f.write(f"Amplitude     : {amp_min:.3e} -- {amp_max:.3e} (global, for normalization)\n")
        f.write("\n--- foam_oscillation_score ---\n")
        mean = score_sum / score_cnt if score_cnt else float("nan")
        f.write(f"max           : {score_max:.6f}\n")
        f.write(f"mean          : {mean:.6f}\n")
        for t in (0.01, 0.05, 0.1, 0.2):
            pct = 100 * ge[t] / score_cnt if score_cnt else 0
            f.write(f"rows >= {t:<4} : {ge[t]:,} ({pct:.4f}%)\n")
        f.write(f"\nScored output : {scored_path.name}\n")
        if is_syn:
            f.write("\nNOTE: This run used SYNTHETIC self-test data. Do not interpret these\n"
                    "      numbers as Voyager measurements.\n")
    print(f"\n[DONE] rows={n_rows:,}  source={'SYNTHETIC' if is_syn else 'REAL'}")
    print(f"[DONE] scored CSV : {scored_path}")
    print(f"[DONE] summary    : {summary_path}")


def build_cli():
    p = argparse.ArgumentParser(
        description="Voyager 0-20 Hz foam-oscillation mapper (Experiment 2, rewritten).",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", help="Data file OR directory of data files.")
    p.add_argument("-o", "--output", default="./voyager_output", help="Output directory.")
    p.add_argument("-v", "--voyager", type=int, default=1, choices=[1, 2], help="Spacecraft (1 or 2).")
    p.add_argument("--fs", type=float, default=1000.0, help="Sampling rate (Hz) for the band-pass.")
    p.add_argument("--min-freq", type=float, default=0.0, help="Band-pass low edge (Hz).")
    p.add_argument("--max-freq", type=float, default=20.0, help="Band-pass high edge (Hz).")
    p.add_argument("--chunk-size", type=int, default=500_000, help="Rows per CSV streaming chunk.")
    p.add_argument("--max-rows", type=int, default=0, help="Cap total rows (0 = no cap).")
    p.add_argument("--selftest", action="store_true",
                   help="Generate labelled SYNTHETIC data instead of reading input.")
    return p


def main():
    args = build_cli().parse_args()
    if not args.selftest and not args.input:
        build_cli().error("provide --input PATH (a file or folder), or --selftest for synthetic data.")
    print(f"[INIT] Voyager {args.voyager}  |  band {args.min_freq}-{args.max_freq} Hz  |  fs={args.fs}")
    print(f"[INIT] input  = {args.input if args.input else '(synthetic self-test)'}")
    print(f"[INIT] output = {args.output}\n")
    run(args)


if __name__ == "__main__":
    main()
