# Voyager Data Analysis

Analysis of electromagnetic data from NASA's Voyager 1 and Voyager 2 spacecraft.

This repository contains tools and experiments for processing and analyzing plasma wave (PWS) and magnetometer (MAG) data from the Voyager missions, with the goal of identifying transient events and potential signatures related to spacetime foam research.

## Repository Structure
Voyager-Data-Analysis/
├── Archive/                  # Archived / deprecated files and early experiments
├── Experiment 1/             # Early foam-oscillation mapping work
├── Experiment 2/             # Streaming / line-by-line analysis experiments
├── Experiment 3/             # Newer / updated experiment work
├── main.py               # CLI orchestration
├── anomaly.py                # Anomaly and transient detection
├── catalog.py                # CDAWeb dataset discovery
├── config.py                 # Configuration and instrument bands
├── download.py               # Chunked data downloading
├── ingest.py                 # Data ingestion and normalization
├── mapping.py                # Visualization and mapping
├── position.py               # SPICE heliocentric positioning
├── spectra.py                # Spectral analysis tools
├── requirements.txt
├── LICENSE
└── README.md
text## Core Pipeline

The root-level modules form a modular pipeline for working with Voyager EM data:

- **`config.py`** — Defines frequency bands and data paths for PWS and MAG instruments
- **`catalog.py`** — Resolves current dataset IDs from CDAWeb
- **`download.py`** — Handles large, resumable downloads from CDAWeb
- **`ingest.py`** — Converts raw data into clean time series
- **`spectra.py`** — Spectrogram generation and spectral processing
- **`anomaly.py`** — Detects spectral events and DC-level shifts
- **`position.py`** — Computes spacecraft distance from the Sun using SPICE
- **`mapping.py`** — Combines results and generates visualizations
- **`__main__.py`** — Command-line interface

### Basic Usage

```bash
pip install -r requirements.txt

# Test the pipeline without downloading data
python -m voyager_em selftest

# Full run (download + analyze)
python -m voyager_em all
Experiment Folders
This repository contains multiple iterations of analysis:






























FolderDescriptionStatusExperiment 1/Early "classic" mapper with foam-oscillation scoringLegacyExperiment 2/Streaming / memory-efficient analysis approachLegacyExperiment 3/Newer / updated experiment workActiveArchive/Deprecated files and early draftsArchived
The Experiment 3/ folder is the most recent area of active work.
Important Instrument Notes

MAG primarily captures near-DC and low-frequency variations. Many heliospheric structures appear as step changes or variance increases rather than narrowband tones.
PWS waveform data has a hard analog high-pass cutoff near 40 Hz. There is no usable data below this frequency in the waveform receiver products.
PWS amplitudes are relative due to automatic gain control and are best used for detection rather than absolute calibration.

Current Status
Work in progress — June 2026
The repository is being actively reorganized. The modern pipeline at the root level is the main focus going forward. Legacy experiment code remains available in the Experiment X/ folders for reference.
Results and methods from this work support broader research into resonant spacetime phenomena.

License
Released under the Knight Industries Proprietary Research License 2026.
See the LICENSE file for terms. Commercial use or redistribution requires prior written permission from Knight Industries.
