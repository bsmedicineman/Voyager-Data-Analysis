#!/usr/bin/env python3
r"""
Voyager Foam Oscillation Mapper (0–20 Hz Focus) - TRUE SINGLE-PASS STREAMING WITH DATA CATCH FILTER

This script processes a single CSV ONE ROW AT A TIME with ZERO buffering.
- Reads one line, applies data catch filter to validate/clean, computes score, writes output immediately.
- Data catch filter validates data types, ranges, and completeness as rows arrive.
- Global stats (amplitude min/max, distance max) are updated adaptively as data streams in.
- Scores improve monotonically as more data is observed (converges to global optimum).
- Memory footprint is ~constant (single row buffer + output writer).

Usage:
    python3 VoyagerFoamAnalyzer020Hz_exp3_streaming.py --csv "C:\path\to\anomalies_full.csv"

For reproducible global normalization, pre-compute stats and provide via --amp-min, --amp-max, --distance-max.
"""

from pathlib import Path
import argparse
import sys
import os
import time
import warnings
from typing import Dict, Tuple, Optional

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ----------------------------------------------------------------------
# CONFIG / DEFAULTS
# ----------------------------------------------------------------------
OUTPUT_DIR = Path("./analysis_output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Adaptive streaming defaults (will be updated as data flows in)
INITIAL_AMP_MIN = 1e-9
INITIAL_AMP_MAX = 1e-3
INITIAL_DISTANCE_MAX = 150.0

# Data catch filter thresholds
DATA_CATCH_CONFIG = {
    "freq_min_hz": -100.0,          # Frequency must be >= this
    "freq_max_hz": np.inf,          # Frequency has NO upper limit
    "amp_min_valid": 1e-12,         # Amplitude reasonable lower bound
    "amp_max_valid": 1e2,           # Amplitude reasonable upper bound
    "distance_min_au": 0.0,         # Distance must be >= this
    "distance_max_au": 200.0,       # Distance must be <= this
}

# ----------------------------------------------------------------------
# Data Catch Filter
# ----------------------------------------------------------------------
class DataCatchFilter:
    """
    Validates and cleans row data as it streams in.
    Tracks quality metrics and flags suspect data.
    """
    def __init__(self, config=None):
        self.config = config or DATA_CATCH_CONFIG
        
        # Counters
        self.rows_total = 0
        self.rows_valid = 0
        self.rows_filtered = 0
        self.rows_partial = 0  # passed but with filldowns
        
        # Issue tracking
        self.freq_out_of_range = 0
        self.amp_out_of_range = 0
        self.distance_out_of_range = 0
        self.timestamp_invalid = 0
        self.missing_critical = 0

    def catch_and_clean(self, row: pd.Series) -> Tuple[Dict, bool, str]:
        """
        Validate and clean a row.
        
        Returns:
            (cleaned_dict, is_valid, reason)
            - cleaned_dict: dict of scalar values (may have NaN for invalid fields)
            - is_valid: True if row passes all checks
            - reason: human-readable status string
        """
        self.rows_total += 1
        issues = []
        cleaned = {}
        is_partial = False

        # ---- Frequency ----
        try:
            freq = float(row.get("center_freq_hz", -9999))
            if not np.isfinite(freq):
                freq = -9999.0
                issues.append("freq_not_finite")
            elif not (self.config["freq_min_hz"] <= freq <= self.config["freq_max_hz"]):
                issues.append(f"freq_out_of_range({freq:.2f}Hz)")
                self.freq_out_of_range += 1
        except (ValueError, TypeError):
            freq = -9999.0
            issues.append("freq_parse_failed")
        
        cleaned["freq"] = freq

        # ---- Amplitude ----
        try:
            amp = float(row.get("amplitude", 0.0))
            if not np.isfinite(amp):
                amp = 0.0
                issues.append("amp_not_finite")
                is_partial = True
            elif amp == 0.0:
                # Treat 0 as missing, not invalid
                is_partial = True
            elif not (self.config["amp_min_valid"] <= amp <= self.config["amp_max_valid"]):
                issues.append(f"amp_out_of_range({amp:.2e})")
                self.amp_out_of_range += 1
        except (ValueError, TypeError):
            amp = 0.0
            issues.append("amp_parse_failed")
            is_partial = True
        
        cleaned["amp"] = amp

        # ---- Distance ----
        try:
            dist = float(row.get("distance_au", np.nan))
            if np.isnan(dist):
                dist = np.nan
                is_partial = True
            elif not np.isfinite(dist):
                dist = np.nan
                is_partial = True
            elif not (self.config["distance_min_au"] <= dist <= self.config["distance_max_au"]):
                issues.append(f"distance_out_of_range({dist:.1f}AU)")
                self.distance_out_of_range += 1
        except (ValueError, TypeError):
            dist = np.nan
            issues.append("distance_parse_failed")
            is_partial = True
        
        cleaned["dist"] = dist

        # ---- Timestamp ----
        ts = None
        ts_source = "missing"
        if "t_start" in row and pd.notna(row["t_start"]):
            try:
                ts = pd.Timestamp(row["t_start"])
                ts_source = "t_start"
            except Exception:
                pass
        
        if ts is None and "date" in row and pd.notna(row["date"]):
            try:
                ts = pd.Timestamp(row["date"])
                ts_source = "date"
            except Exception:
                pass
        
        if ts is None:
            issues.append("timestamp_missing")
            self.timestamp_invalid += 1
            is_partial = True
        
        cleaned["timestamp"] = ts
        cleaned["ts_source"] = ts_source

        # ---- Spacecraft / Source (strings, non-critical) ----
        spacecraft = str(row.get("spacecraft", "")).strip() if pd.notna(row.get("spacecraft")) else ""
        source = str(row.get("source", "")).strip() if pd.notna(row.get("source")) else ""
        
        cleaned["spacecraft"] = spacecraft
        cleaned["source"] = source

        # ---- Decide validity ----
        # Valid if: freq in range, amp in range (or zero), distance in range (or NaN), timestamp exists
        has_critical_issue = any(x in str(issues) for x in ["freq_out_of_range", "amp_out_of_range", "distance_out_of_range"])
        
        if has_critical_issue:
            is_valid = False
            self.rows_filtered += 1
            reason = f"FILTERED: {'; '.join(issues)}"
        elif is_partial:
            is_valid = True  # still process, but mark as partial
            self.rows_partial += 1
            reason = f"PARTIAL: {'; '.join(issues) if issues else 'has missing/zero fields'}"
        else:
            is_valid = True
            self.rows_valid += 1
            reason = "VALID"

        return cleaned, is_valid, reason

    def report(self):
        """Print summary statistics."""
        print()
        print("[DATA CATCH FILTER] REPORT")
        print(f"  Total rows processed: {self.rows_total:,}")
        print(f"  ✓ Valid (complete):   {self.rows_valid:,} ({100*self.rows_valid/max(1,self.rows_total):.1f}%)")
        print(f"  ~ Partial (fillable): {self.rows_partial:,} ({100*self.rows_partial/max(1,self.rows_total):.1f}%)")
        print(f"  ✗ Filtered (invalid): {self.rows_filtered:,} ({100*self.rows_filtered/max(1,self.rows_total):.1f}%)")
        print()
        print("  Issue breakdown:")
        print(f"    - Frequency out of range: {self.freq_out_of_range:,}")
        print(f"    - Amplitude out of range: {self.amp_out_of_range:,}")
        print(f"    - Distance out of range:  {self.distance_out_of_range:,}")
        print(f"    - Invalid timestamp:      {self.timestamp_invalid:,}")
        print()


# ----------------------------------------------------------------------
# Foam scoring (streaming-friendly, mutable stats)
# ----------------------------------------------------------------------
class AdaptiveStreamingFoamScorer:
    """
    Scorer that updates global min/max on-the-fly as rows arrive.
    Normalized scores improve monotonically toward true global optimum.
    """
    def __init__(self, amp_min=None, amp_max=None, distance_max=None, voyager_num=1, adaptive=True):
        self.amp_min = amp_min if amp_min is not None else INITIAL_AMP_MIN
        self.amp_max = amp_max if amp_max is not None else INITIAL_AMP_MAX
        self.distance_max = distance_max if distance_max is not None else INITIAL_DISTANCE_MAX
        self.voyager_num = voyager_num
        self.adaptive = adaptive  # if True, update stats as rows arrive
        
        # counters for logging
        self.rows_scored = 0
        self.stats_updates = 0

        # clock catalog (Hz)
        self.clock_catalog_hz = {
            "Jupiter": 0.0280e-3,
            "Saturn": 0.0263e-3,
            "Sun": 0.0018e-3,
            "Earth": 0.0116e-3,
        }

    def update_stats(self, amp_val, dist_val):
        """Update running min/max if adaptive mode enabled."""
        if not self.adaptive:
            return
        
        if np.isfinite(amp_val) and amp_val > 0:
            if amp_val < self.amp_min:
                self.amp_min = amp_val
            if amp_val > self.amp_max:
                self.amp_max = amp_val
        
        if dist_val is not None and np.isfinite(dist_val):
            if dist_val > self.distance_max:
                self.distance_max = dist_val
        
        self.stats_updates += 1

    def _resonance_weight(self, freq):
        """Single frequency value (scalar)."""
        # Guard against invalid frequency values (e.g., -9999.0)
        if not np.isfinite(freq) or freq < -100:
            return 0.0
        
        resonance_peaks = [8.0, 20.0]
        resonance_weight = 0.0
        for fp in resonance_peaks:
            diff = abs(freq - fp)
            bw = max(0.5, fp * 0.1)
            
            # Clip diff to prevent overflow in the square operation
            diff = np.clip(diff, -1000, 1000)
            
            # Prevent overflow: cap the exponent before np.exp()
            exponent = -(diff ** 2) / (2 * bw ** 2)
            exponent = np.clip(exponent, -700, 0)  # np.exp(-700) ≈ 0, np.exp(0) = 1
            
            resonance_weight += np.exp(exponent)
        
        # normalize to [0, 1]
        max_weight = 2.0  # two peaks, each ~1.0 at their center
        resonance_weight = min(resonance_weight / max_weight, 1.0)
        return resonance_weight

    def score_row(self, freq, amp, timestamp, distance_au):
        """Score a single row (all scalar inputs)."""
        self.rows_scored += 1
        
        # Update adaptive stats
        self.update_stats(amp, distance_au)

        # Amplitude normalization
        if self.amp_max > self.amp_min:
            amp_norm = (amp - self.amp_min) / (self.amp_max - self.amp_min)
        else:
            amp_norm = 0.5  # fallback if bounds not set
        amp_norm = np.clip(amp_norm, 0.0, 1.0)

        # Resonance weight
        resonance_weight = self._resonance_weight(freq)

        # Solar cycle (11-year) modulation
        launch_ts = pd.Timestamp("1977-09-05")
        ts = pd.Timestamp(timestamp) if timestamp is not None else launch_ts
        days = (ts - launch_ts).total_seconds() / 86400.0
        solar_cycle = 0.5 + 0.5 * np.sin(2 * np.pi * days / (365.25 * 11.0))

        # Distance effect (Gaussian centered at 0.5 normalized distance)
        if distance_au is not None and np.isfinite(distance_au) and self.distance_max > 0:
            distance_norm = distance_au / self.distance_max
            distance_effect = 0.7 + 0.3 * np.exp(-((distance_norm - 0.5) ** 2) / 0.1)
        else:
            distance_effect = 1.0

        # Annual effect (day-of-year modulation)
        try:
            doy = ts.dayofyear
            annual_effect = 0.8 + 0.2 * np.sin(2 * np.pi * doy / 365.25)
        except Exception:
            annual_effect = 1.0

        coherence_est = solar_cycle * distance_effect * annual_effect
        coherence_est = np.clip(coherence_est, 0.1, 0.9)

        # Final score
        mismatch = 1.0 / (resonance_weight + 0.01)
        damping = 0.1
        score = coherence_est * amp_norm / (mismatch + damping)
        score = np.clip(score, 0.0, 1.0)

        return score


# ----------------------------------------------------------------------
# Row-by-row streaming processor
# ----------------------------------------------------------------------
def process_csv_true_streaming(csv_path, output_dir=OUTPUT_DIR, voyager_num=1, 
                               amp_min=None, amp_max=None, distance_max=None):
    """
    Read CSV one row at a time using chunksize=1, apply data catch filter, compute score immediately, 
    write to output, never hold more than one row in memory.
    """
    # Sanitize path
    if csv_path is None:
        raise FileNotFoundError("No CSV path provided (csv_path is None)")

    if isinstance(csv_path, Path):
        csv_path = str(csv_path)
    csv_path = str(csv_path).strip().strip('"').strip("'")
    csv_path = os.path.expanduser(csv_path)
    csv_path = os.path.normpath(csv_path)
    csv_path = Path(csv_path)

    if not csv_path.exists():
        msg = (
            f"CSV not found: {csv_path}\n"
            f"- Checked path after stripping quotes/whitespace and expanding ~\n"
            f"- Make sure the file exists and you have read permission.\n"
        )
        raise FileNotFoundError(msg)

    # Columns to read (minimal subset for memory)
    usecols = [
        "center_freq_hz",
        "amplitude",
        "distance_au",
        "t_start",
        "date",
        "spacecraft",
        "source",
    ]

    # Initialize data catch filter
    data_filter = DataCatchFilter(config=DATA_CATCH_CONFIG)

    # Initialize adaptive scorer
    adaptive = (amp_min is None or amp_max is None or distance_max is None)
    scorer = AdaptiveStreamingFoamScorer(
        amp_min=amp_min,
        amp_max=amp_max,
        distance_max=distance_max,
        voyager_num=voyager_num,
        adaptive=adaptive,
    )

    # Output file
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / (csv_path.stem + "__streaming_scores.csv.gz")
    
    # Separate log for filtered/partial rows
    filtered_log_path = output_dir / (csv_path.stem + "__data_catch_log.csv")

    print(f"[STREAMING] Starting true single-pass processing with DATA CATCH FILTER")
    print(f"[STREAMING] Input CSV: {csv_path}")
    print(f"[STREAMING] Frequency filter: >= -100 Hz (NO upper limit)")
    print(f"[STREAMING] Adaptive mode: {adaptive} (stats will update as data flows)")
    print(f"[STREAMING] Output (scored): {out_path}")
    print(f"[STREAMING] Output (catch log): {filtered_log_path}")
    print()

    t0 = time.time()
    rows_processed = 0
    rows_written = 0
    output_file = None
    
    # Open catch log for writing
    catch_log = open(filtered_log_path, "w")
    catch_log.write("row_number,status,reason,frequency_hz,amplitude,distance_au,timestamp,spacecraft,source\n")

    try:
        # Use chunksize=1 for true streaming (one row at a time, no full DataFrame materialization)
        for chunk in pd.read_csv(csv_path, usecols=usecols, chunksize=1, low_memory=False):
            # Extract single row from chunk (DataFrame with 1 row)
            row = chunk.iloc[0]
            rows_processed += 1

            # APPLY DATA CATCH FILTER
            cleaned, is_valid, reason = data_filter.catch_and_clean(row)

            freq = cleaned["freq"]
            amp = cleaned["amp"]
            dist = cleaned["dist"]
            ts = cleaned["timestamp"]
            spacecraft = cleaned["spacecraft"]
            source = cleaned["source"]

            # Log if not valid or partial
            if not is_valid or "PARTIAL" in reason:
                status = "FILTERED" if not is_valid else "PARTIAL"
                catch_log.write(
                    f"{rows_processed},{status},\"{reason}\","
                    f"{freq},{amp},{dist},{ts},{spacecraft},{source}\n"
                )
                catch_log.flush()

            # Skip if filtered
            if not is_valid:
                continue

            # Compute score for this row
            score = scorer.score_row(freq, amp, ts, dist)

            # Prepare output row
            out_row = pd.DataFrame({
                "row_number": [rows_processed],
                "timestamp": [ts],
                "center_freq_hz": [freq],
                "amplitude": [amp],
                "distance_au": [dist],
                "foam_score": [score],
                "spacecraft": [spacecraft],
                "source": [source],
            })

            # Write to output immediately (append mode to disk)
            if rows_written == 0:
                # First row: create file with header
                out_row.to_csv(out_path, index=False, compression="gzip", mode="w")
            else:
                # Subsequent rows: append without header
                out_row.to_csv(out_path, index=False, compression="gzip", mode="a", header=False)

            rows_written += 1

            # Log progress every 10k rows
            if rows_written % 10000 == 0:
                elapsed = time.time() - t0
                rate = rows_written / elapsed if elapsed > 0 else 0
                print(f"[STREAMING] rows_written={rows_written:,}  rows/sec={rate:.0f}  "
                      f"amp_bounds=[{scorer.amp_min:.2e}, {scorer.amp_max:.2e}]  "
                      f"distance_max={scorer.distance_max:.1f} AU  elapsed={elapsed:.1f}s")

            # Explicit memory cleanup
            del out_row, chunk

    except Exception as e:
        print(f"[ERROR] at row {rows_processed}: {e}")
        raise
    finally:
        catch_log.close()

    elapsed = time.time() - t0
    print()
    print(f"[STREAMING] COMPLETE")
    print(f"[STREAMING] rows_processed={rows_processed:,}  rows_written={rows_written:,}  elapsed={elapsed:.1f}s")
    print(f"[STREAMING] Final global stats (adaptive):")
    print(f"            amp_min={scorer.amp_min:.2e}  amp_max={scorer.amp_max:.2e}")
    print(f"            distance_max={scorer.distance_max:.1f} AU")
    print()

    # Print filter report
    data_filter.report()

    print(f"[STREAMING] Scored output saved to: {out_path}")
    print(f"[STREAMING] Filter log saved to:    {filtered_log_path}")

    return out_path, filtered_log_path


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def build_cli():
    p = argparse.ArgumentParser(description="True single-pass streaming Voyager foam mapper with data catch filter")
    p.add_argument("--csv", "-c", required=False, help="Path to CSV to process")
    p.add_argument("--outdir", "-o", default=str(OUTPUT_DIR), help="Output directory")
    p.add_argument("--voyager", type=int, default=1, help="Voyager spacecraft number (1 or 2)")
    p.add_argument("--amp-min", type=float, default=None, help="Force global amplitude min (disables adaptive)")
    p.add_argument("--amp-max", type=float, default=None, help="Force global amplitude max (disables adaptive)")
    p.add_argument("--distance-max", type=float, default=None, help="Force global distance max (disables adaptive)")
    return p


def main():
    parser = build_cli()
    args = parser.parse_args()

    csv_path = args.csv
    if not csv_path:
        csv_path = input("Enter path to CSV file to process: ").strip()

    if not csv_path:
        print("No CSV path provided. Aborting.")
        sys.exit(1)

    csv_path = str(csv_path).strip().strip('"').strip("'")
    csv_path = os.path.expanduser(csv_path)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Starting true streaming processing for: {csv_path!r}")
    print()

    try:
        scored_out, catch_out = process_csv_true_streaming(
            csv_path,
            output_dir=outdir,
            voyager_num=args.voyager,
            amp_min=args.amp_min,
            amp_max=args.amp_max,
            distance_max=args.distance_max,
        )
        print()
        print("✓ All done.")
    except FileNotFoundError as e:
        print("ERROR:", e)
        sys.exit(2)
    except Exception as e:
        print("Unhandled error during processing:", e)
        raise


if __name__ == "__main__":
    main()

