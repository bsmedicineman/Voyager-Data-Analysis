# Voyager Data Analysis

Analysis of electromagnetic (EM) data from NASA's Voyager 1 and Voyager 2 spacecraft.

This repository contains tools and experiments for processing plasma wave (PWS) and magnetometer (MAG) data from the Voyager missions, focused on detecting transients and potential spacetime foam signatures in the heliosphere.

## Repository Structure

```
Voyager-Data-Analysis/
├── Archive/                  # Archived and deprecated files
├── Experiment 1/             # Early foam-oscillation mapping experiments
├── Experiment 2/             # Streaming / line-by-line analysis experiments
├── Experiment 3/             # Newer and actively updated experiment work
├── __main__.py               # CLI entrypoint and orchestration
├── anomaly.py                # Spectral event and DC shift detection
├── catalog.py                # CDAWeb dataset discovery
├── config.py                 # Instrument configuration and frequency bands
├── download.py               # Chunked, resumable data downloads
├── ingest.py                 # Data normalization and ingestion
├── mapping.py                # Results merging and visualization
├── position.py               # SPICE heliocentric distance calculations
├── spectra.py                # Spectrogram and spectral analysis
├── requirements.txt
├── LICENSE
└── README.md
```

## Core Pipeline

The root-level modules provide a modular pipeline for Voyager EM data analysis:

| Module        | Purpose                                           |
|---------------|---------------------------------------------------|
| `config.py`   | Defines instrument bands and data paths           |
| `catalog.py`  | Resolves current CDAWeb dataset IDs               |
| `download.py` | Downloads large CDF files in resumable chunks     |
| `ingest.py`   | Converts cached data into clean time series       |
| `spectra.py`  | Generates spectrograms and computes band power    |
| `anomaly.py`  | Detects spectral anomalies and DC baseline shifts |
| `position.py` | Calculates spacecraft heliocentric distance       |
| `mapping.py`  | Merges data and generates visualizations          |
| `__main__.py` | Command-line interface and workflow orchestration |

### Quick Start

```bash
pip install -r requirements.txt

# Run offline self-test (no internet required)
python -m voyager_em selftest

# Full pipeline (download + analyze)
python -m voyager_em all
```

## Experiment Folders

| Folder        | Description                                      | Status     |
|---------------|--------------------------------------------------|------------|
| `Experiment 1/` | Early “classic” mapper with foam-oscillation scoring | Legacy     |
| `Experiment 2/` | Streaming / memory-efficient analysis approach     | Legacy     |
| `Experiment 3/` | Newer and actively developed experiment work       | Active     |
| `Archive/`      | Deprecated files and early drafts                  | Archived   |

## Important Physics Notes

- **MAG** data is primarily sensitive to near-DC fields and step changes. Many heliospheric features appear as baseline shifts or variance increases.
- **PWS** waveform receiver has a hard ~40 Hz analog high-pass cutoff. There is no recoverable data below this frequency.
- PWS amplitudes are relative (due to automatic gain control) and should be used for detection and localization rather than absolute field strength.
- Electric (PWS) and magnetic (MAG) data are analyzed as separate layers.

## Current Status

**Work in progress** — June 2026

The repository is being actively reorganized. The modern pipeline at the root is the primary development focus. Legacy experiment code is preserved in the `Experiment X/` folders for reference.

## License

This project is released under the **Knight Industries Proprietary Research License 2026**.

See the `LICENSE` file for full terms. Commercial use requires prior written permission.
```
