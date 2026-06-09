# voyager_em — broad-spectrum EM mapping of Voyager 1 & 2

Downloads Voyager plasma-wave (PWS) and magnetometer (MAG) data, analyzes each
across its full instrument band, detects anomalies, locates each one in space,
and renders a map of EM fluctuation vs distance from the Sun.

## Install

```bash
pip install -r requirements.txt
```

## Run

```bash
# 1. Prove the pipeline works offline (synthetic data, no network):
python -m voyager_em selftest

# 2. On a networked machine, fetch positions + resolve datasets + pull data:
python -m voyager_em kernels        # SPICE kernels (leap-sec, de440, Voyager SPK)
python -m voyager_em discover       # resolve live CDAWeb dataset ids
python -m voyager_em download       # chunked, resumable, multi-decade
python -m voyager_em map            # analyze + locate + render figures

# or the whole chain:
python -m voyager_em all
```

Outputs land in `./voyager_em_out/`: `raw/` (cached chunks), `tables/anomalies.csv`
(every anomaly with frequency, bandwidth, start/stop, duration, distance), and
`figures/` (the distance-vs-amplitude-vs-frequency maps).

## Data sources (all public, free)

* CDAWeb / SPDF (NASA Goddard) — CDF + REST API, accessed via `cdasws`.
* PDS-PPI node (UCLA) — authoritative archive.
* U. Iowa Plasma Wave Group — docs + PWS calibration tables.
* NAIF SPICE kernels — spacecraft trajectory for heliocentric distance.

`discover` resolves the exact CDAWeb dataset ids at runtime (they drift), so
nothing is hard-coded; override `Product.dataset_id` in `config.py` if needed.

## The physics this design respects (read before trusting a result)

* **MAG is your DC-to-~8 Hz layer.** A fluxgate magnetometer measures the field
  from 0 Hz up; the heliospheric structure (termination shock, heliopause,
  sector boundaries) shows up as near-DC steps and variance spikes, not as a
  tone at some frequency. The detector handles these as `dc_shift` events.
  Archived MAG products are heavily averaged (1.92 s → ~0.26 Hz Nyquist), so
  true 1–8 Hz content needs the full-resolution 60 ms data.

* **PWS has no data below 40 Hz — at all.** The waveform receiver is AC-coupled
  through a 40 Hz–12 kHz *analog* bandpass; sub-40 Hz was never digitized and
  cannot be recovered by any software. The spectrum analyzer's lowest channel is
  10 Hz (16 discrete log-spaced channels, 4 s cadence — intensities, not a
  waveform). So PWS contributes 40 Hz upward only.

* **PWS amplitude is relative, not calibrated.** Automatic gain control with no
  gain in telemetry → use it to find/locate features, not for absolute V/m.
  MAG is properly calibrated in nT.

* **E and B are different instruments and units** (V/m relative vs nT). They are
  analyzed and plotted as parallel layers, never merged into one filtered band.

## Module map

| stage | file | does |
|-------|------|------|
| config    | `config.py`   | products, real instrument bands, paths |
| discover  | `catalog.py`  | resolve CDAWeb dataset ids from keywords |
| download  | `download.py` | chunked resumable CDF pull |
| ingest    | `ingest.py`   | cached chunks → normalized `Series` |
| spectra   | `spectra.py`  | bandpass / Welch PSD / spectrogram |
| anomaly   | `anomaly.py`  | spectral events + DC baseline shifts |
| position  | `position.py` | SPICE heliocentric distance |
| mapping   | `mapping.py`  | merged table + the EM map figures |
| CLI       | `__main__.py` | orchestration + offline self-test |

## Tuning

`Config.anomaly_z` (default 6.0) sets detection strictness; `anomaly_window` is
the rolling-baseline length. The synthetic self-test is the fastest way to see
the effect of changes before committing to a multi-decade run.
