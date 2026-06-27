# Spacefoam EM Pipeline

Programs for testing whether spacetime foam leaves detectable signatures in
spacecraft electromagnetic data:

1. **`download_em_data.py`** — concurrently downloads and sorts EM (magnetometer /
   E-field / search-coil) CDF data from the inner-heliosphere & Mercury missions.
2. **`analyze_em_foam.py`** — *steady-state lens*. Streams each file and flags
   narrowband, persistent, drifting spectral lines as candidate foam signatures.
3. **`analyze_em_echo.py`** — *transient lens*. Finds large "drive" spikes (bow
   shocks, magnetopause crossings, intense activity) and looks in the trailing
   interval for delayed structured echoes, excess persistence, and downward
   frequency shifts — the Earth Schumann template applied to flyby spikes.

Shared code lives in `em_foam_common.py`. Keep all files in the same folder.

Covered missions (per your source list): **Parker Solar Probe, Solar Orbiter,
MESSENGER, Helios 1/2** download automatically via NASA CDAWeb. **BepiColombo**
lives in ESA's PSA (a different protocol) and is flagged for a manual/extension
path rather than silently skipped.

---

## Install

```bash
pip install -r requirements.txt
```

`cdasws` (NASA's CDAWeb client) is only needed for the downloader's discovery
step; the analyzer needs only `cdflib numpy scipy`.

---

## Workflow

### 1. (optional) confirm current dataset IDs
CDAWeb IDs drift over time, so check what exists right now:

```bash
python download_em_data.py list-datasets --missions psp solo messenger helios
```

### 2. download
```bash
python download_em_data.py fetch \
    --missions psp solo messenger helios \
    --start 2022-01-01 --end 2022-02-01 \
    --output-dir /data/spacefoam \
    --workers 8
```

Files land under `/<output-dir>/<MISSION>/<DATASET_ID>/<file>.cdf`.
Preview first with `--dry-run`. Already have URLs? Use `--url-list urls.txt`.

### 3. analyze
```bash
python analyze_em_foam.py run \
    --input  /data/spacefoam \
    --output /data/foam_features.csv \
    --window 65536 --overlap 0.5 \
    --fmin 1e-3 --fmax 0.5 \
    --p-fa 1e-3 --q-min 20 \
    --known-lines 0.00417,0.0083     # spin harmonics, reaction wheels, heaters...
```

Output CSV carries, for every window: **timestamps** (start/end UTC), **location**
(`pos_r` + x/y/z when the file has ephemeris), and **EM features** (mean/std,
PSD slope, spectral flatness, dominant line frequency/strength/significance,
Q-factor, coherence time, frequency drift, and a `candidate` flag with reason).

---

## Why it meets the spec

| Requirement | How |
|---|---|
| **Negligible memory** | Downloads stream in 64 KB chunks; analysis keeps only one window (+ a few-MB read block) in RAM. A 2 M-record file analyzes at ~constant ~200 MB RSS — almost all of which is the Python/scipy import baseline, not the data. |
| **Fast** | Parallel download threads saturate bandwidth; analysis reads in large blocks and slides windows in-memory to minimize per-call overhead. |
| **Self-recovering** | Downloader keeps a SQLite manifest + atomic `.part` files with HTTP-Range resume; just re-run the same command to continue. Deleted/missing completed files are detected and re-queued. Analyzer logs completed inputs to `<output>.csv.done` and skips them on rerun; a crash loses at most the current window. |
| **CLI input/output paths** | `download ... -o/--output-dir`, `analyze ... -i/--input -o/--output`. |
| **EM only: data + timestamps + location** | The analyzer extracts field variables, epoch, and position; every CSV row has all three. |
| **Stream loaded file → output** | Each window's features are written and flushed straight to the CSV; nothing accumulates. |

Graceful pause/resume: hit **Ctrl-C** during a download — in-flight chunks finish,
the manifest is saved, and the next run picks up exactly where it stopped.

---

## What the "foam signature" detector actually does

For each window it computes a Welch PSD of `|B|`, fits a robust running-median
**continuum** (so the turbulent background doesn't hide lines), and searches for
**narrowband excess** above that continuum. The null model is colored noise:
with *K* averaged Welch segments, PSD/background follows a Gamma(*K*, 1/*K*)
distribution, so the detection threshold is set for a chosen per-window
false-alarm probability (`--p-fa`, Bonferroni-corrected across bins). For the
strongest line it estimates a **quality factor** *Q = f₀/Δf* (persistence /
coherence time) and tracks the line's **frequency drift** across windows.

A window is flagged `candidate=1` only if its dominant line is (a) statistically
significant, (b) sufficiently narrow/persistent (`--q-min`), and (c) not within
tolerance of a user-supplied **known instrument line**.

### Honest interpretation
A `candidate` is an **anomaly that survives a colored-noise null test** — not a
foam detection. Real narrowband features in this data are usually mundane:
spacecraft interference (reaction wheels, heaters, spin harmonics), instrument
lines, or natural plasma waves (ion-cyclotron, Langmuir, whistler). The tool
helps you *find and characterize* candidates and veto known lines, but
attributing one to new physics requires ruling those sources out — ideally with
a pre-registered test, multi-instrument coincidence, and independent passes.
That epistemic burden is on the analysis you build around these candidates, not
on the flag itself.

---

## The transient lens (`analyze_em_echo.py`)

A different lens on the same archives — and the one motivated by your single
positive Earth result. Instead of steady-state lines it hunts **transient-driven
post-spike echoes**: a large coherent drive (bow shock, magnetopause crossing,
intense magnetospheric interval) followed by a delayed, structured, longer-than-
expected trailing response.

```bash
python analyze_em_echo.py run \
    --input /data/spacefoam /data/voyager \    # multiple roots OK
    --output /data/echo_events.csv \
    --env-seconds 1 --trail-seconds 3600 \
    --drive-seconds 30 \
    --k 6 --min-spike-snr 3 --min-echo-z 6 --min-persistence 3
```

One CSV row per detected drive spike, each carrying timestamp, location, and:

- `spike_snr`, `step_mean` — spike strength and the mean-|B| jump (shock signature)
- `relax_time_s`, `rise_time_s`, `persistence_ratio` — relaxation vs the drive's
  own rise time; `persistence_ratio ≫ 1` means a long tail beyond the drive's
  intrinsic timescale (excess persistence)
- `echo_delay_s`, `echo_z`, `echo_count` — delayed recurrences of the drive in
  the trailing envelope (secondary peaks above a robust trailing baseline)
- `drive_centroid`, `echo_centroid`, `freq_shift` — spectral centroid of drive vs
  echo; **negative `freq_shift` is the Schumann-like downward shift** you noted
- `spectral_similarity` — correlation of drive vs echo spectra (morphology match)
- `candidate`, `candidate_reason` — flagged if the spike is real **and** shows a
  strong delayed echo or excess persistence; the reason string spells out which

**How it stays light:** pass 1 builds a *decimated activity envelope* for the
whole file (one point per `--env-seconds`), and spikes, persistence, and echoes
are all measured on that envelope. Spectral metrics read only short raw snapshots
at the spike and echo lags. Memory stays bounded (~140 MB measured on a 400k-
record file) regardless of sampling rate or how long the trailing window is — so
high-rate PWS waveform data over hour-long trails is fine.

**Inputs:** `.cdf` (everything from the downloader, plus Voyager CDFs from
PDS/CDAWeb) and simple `.csv`/`.tab` (a time column + numeric field columns).
For 2-D spectral products (e.g. PWS spectrum-analyzer channels), the activity
proxy is total power across channels and the centroid/shift are reported in
**channel-index** units (labeled `centroid_units`) since per-channel Hz needs the
instrument frequency table.

**Tuning for Voyager flybys:** point `--input` at your Voyager files, set
`--env-seconds` to a few seconds (or larger for low-cadence cruise data — it
auto-floors at `--min-env-samples`), and `--drive-seconds` to the expected event
duration. Start with one or two well-studied events (e.g. V1 Jupiter inbound bow
shock) before scaling up, exactly as your reframing note suggests.

### Honest interpretation (same standard as the steady-state lens)
A `candidate` is a delayed/persistent feature that stands above the ambient null
— **not** a foam detection. Giant-planet magnetospheres are full of natural
echoes: whistlers, chorus, reflected and mode-converted waves, and the rich
post-shock turbulence itself. A candidate must survive that scrutiny and the same
confounder regression (Kp/lightning-style controls) you applied on Earth before
it says anything about new physics. The tool finds and characterizes candidates;
the inferential burden stays with the follow-up analysis.

---

## Self-tests (run anytime, no network)

```bash
python analyze_em_foam.py selftest      # injects a known line, confirms recovery + clean null
python analyze_em_echo.py selftest      # injects a spike+echo, confirms delay/shift recovery + clean control
python download_em_data.py selftest     # simulates a crash, confirms correct resume
```

---

## Notes & extension points

- **Position / ephemeris:** many L2 field CDFs don't embed spacecraft position.
  The analyzer captures it when present and otherwise leaves `pos_*` blank; for
  rigorous heliocentric coordinates, merge a separate ephemeris/SPICE product by
  time. (Hook: extend `_sample_position` / `find_position_variable`.)
- **BepiColombo / ESA PSA:** add a discovery function returning CDF URLs from the
  PSA TAP service and feed them through `--url-list`; the download/manifest path
  is identical.
- **Search-coil / spectral-matrix products** (already-PSD data) are skipped by
  the time-series analyzer; they need a separate reader.








Run these one at a time. Each is a single line — use the copy button on each block so nothing wraps or picks up a stray character.
1. See how big the light survey is (downloads nothing):
python download_em_data.py fetch-all -o spacefoam_data --light --dry-run
2. Download the light survey and leave it (small, fast, all four missions):
python download_em_data.py fetch-all -o spacefoam_data --light
3. When you want the full-resolution data — check the size first:
python download_em_data.py fetch-all -o spacefoam_data --dry-run
4. Then download the full set and leave it (many GB, resumable):
python download_em_data.py fetch-all -o spacefoam_data
5. Analyze — steady-state lens (narrowband foam lines):
python analyze_em_foam.py run -i spacefoam_data -o foam_features.csv
6. Analyze — transient lens (post-spike echoes):
python analyze_em_echo.py run -i spacefoam_data -o echo_events.csv
That's the whole workflow. Pick 2 or 4 for the download (not both — 4 is a superset). Re-running any download command resumes where it left off, so if it stops you just paste the same line again.







