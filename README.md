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


Voyager‑Data Experiments 

The Voyager‑Data repository contains (currently 06.13.2026) three self‑contained Python experiments that work directly on the raw Voyager plasma‑wave (PWS) and magnetic‑field (MAG) telemetry that NASA released in the Planetary Data System (PDS). All three experiments produce event catalogs (times, frequencies, scores, candidate foamy signatures) that can be fed into the Foam‑Modeling repo (baseline_foam_map.py, solar_foam_animation.py, craft_trajectory_plotter.py).  


📂 Repository layout (Voyager‑Data)

          
            
            
          
          voyager-data/

├─ Voyager data download program

│
│

├─  # Experiment 1 – “classic” mapper    

├─  # Experiment 2 – true single‑pass streaming

├─  # Experiment 3 – chunked loading + diagnostics  ***NOT CURRENTLY UPLOADED TO GITHUB

│
│
│
│

├─ data/                                      # Raw Voyager EM/MAG files (CSV, TAB, NC, …)  ***NOT UPLOADED TO GITHUB

│

└─ docs/                                      # Optional notebooks / plots   ***NOT CURRENTLY UPLOADED TO GITHUB, IN PROGRESS
      
1️⃣ Experiment 1 – Classic 0‑20 Hz Foam‑Oscillation Mapper
File: voyager_foam_mapper_020hz.txt



What it does
How it works
Main outputs



Loads every raw Voyager EM/MAG file (CSV, TAB, NC, LBLX, BSP, etc.) using a robust multi‑delimiter parser.
• Detects the most likely amplitude column (e.g. mag_amplitude, pws_amplitude). <br>• Applies a 0‑20 Hz Butterworth band‑pass filter (order 4). <br>• Generates a synthetic “frequency” column from the filtered amplitude (the same trick used in the original research).
voyager_1_foam_events_020hz.csv – a table of foam‑oscillation events (start/end time, duration, peak score, mean frequency, distance from Sun, frequency‑band label). <br>voyager_1_data_summary.txt – quick sanity check (rows, columns, time range, distance range).


Computes the original “foam‑oscillation score” (Section 5 of the Master‑Equation reference).
Normalises amplitude, adds resonance weighting around 8 Hz and 20 Hz, modulates by a solar‑cycle term, distance‑effect Gaussian, and an annual sinusoid.
The score column (foam_oscillation_score) is stored in the DataFrame and later used for thresholding.


Detects bursts by sliding a 5‑year window, taking the 95 th percentile of the score, and grouping contiguous high‑score points (≥ 1 day gap).
Event‑processing routine (_process_event) computes duration, distance travelled, frequency statistics, and a handful of “stiffness/flow” proxy metrics.
Event‑level CSV (see above) plus a set of frequency‑band statistics (voyager_1_frequency_band_stats.csv).


Why it matters for the my other Foam‑Modeling repo

The CSV of foam events can be loaded as a source of “real” foam disturbances and plotted alongside the synthetic foam field generated by solar_foam_animation.py.  
The event timestamps and frequencies provide ground‑truth anchors for the Resonance Engine (the driven‑Mathieu model) – you can compare the predicted mode amplitudes at the Voyager location with the observed peaks.  
The distance‑versus‑time information lets the Foam‑Detector (core of the detection pipeline) evaluate whether a candidate matches the expected gravitational‑foam scaling (Eq. 9.1).


2️⃣ Experiment 2 – True Single‑Pass Streaming with Data‑Catch Filter
File: voyager_foam_analyzer020hz_exp3.txt



What it does
How it works
Main outputs



Processes the raw CSV line‑by‑line (chunksize=1) – constant O(1) memory usage.
• DataCatchFilter validates each row (frequency range, amplitude limits, distance limits, timestamp sanity). <br>• Rows that fail are logged to *_data_catch_log.csv. <br>• Valid rows are passed to AdaptiveStreamingFoamScorer.
*_streaming_scores.csv.gz – one row per valid observation, containing row_number, timestamp, center_freq_hz, amplitude, distance_au, and the online foam score (foam_score).


Updates global amplitude/min‑max and distance‑max adaptively as data streams in, so the normalisation converges toward the true global optimum.
The scorer keeps running minima/maxima (amp_min, amp_max, distance_max) and recomputes the normalised amplitude and resonance weight for every new sample.
Final adaptive statistics printed at the end of the run (e.g. amp_min=2.3e‑12, amp_max=1.4e‑2).


Provides a “live” score that monotonically improves as more data are read, ideal for on‑board or near‑real‑time monitoring.
No post‑hoc batch processing; the score is emitted instantly after each row.
The compressed CSV can be visualised as a time‑series in any plotting tool (Matplotlib, Plotly) without loading the whole dataset into RAM.


Why it matters for the Foam‑Modeling repo

The streaming scores are a real‑time analogue of the static foam‑oscillation score used by the detection pipeline. They can be fed directly into the Signal Processor (matched‑filter) as a time‑ordered series, bypassing the need for a full DataFrame.  
The Data‑Catch Log gives a clean “quality‑control” mask that can be applied when you later compute the foam‑gradient force for Voyager’s position (the log tells you which rows have reliable distance_au and center_freq_hz).  
Because the script works with a single CSV (the Voyager PWS “anomalies_full.csv” that many researchers already use), you can run it on the same file that you later load into the Foam‑Detector for a side‑by‑side comparison of online vs offline scores.


3️⃣ Experiment 3 – Version 2: Chunked Loading, Diagnostics & Synthetic Fallback
File: voyager_foam_mapper_020hz_version_2.txt



What it does
How it works
Main outputs



Scans the raw‑data folder, loads every supported file type (CSV/TAB, NetCDF, ASCII, BSP) in chunks of 1 000 000 rows to keep memory usage bounded.
• Uses the same multi‑delimiter logic as Experiment 1. <br>• Detects the amplitude column, applies the 0‑20 Hz band‑pass filter, writes a temporary “amplitude_filtered” column. <br>• Stores each chunk to a temporary pickle/CSV file, then concatenates all chunks (down‑sampling if total rows exceed max_rows_in_memory).
A single consolidated DataFrame (self.data) that contains the full Voyager data set (or a synthetic fallback if no raw files are found).


Runs a diagnostics routine (run_diagnostics()) that reports file‑type distribution, sample sizes, and attempts a manual parse of the first three files.
Helpful when the raw folder contains legacy PDS/VICAR formats that the automatic loader cannot handle.
Diagnostic log printed to STDOUT; any failures fall back to synthetic data (a 1‑year, 1 Hz synthetic mission) so the rest of the pipeline never crashes.


Performs the same foam‑oscillation detection as Experiment 1, but on the concatenated data set.
Uses the same scoring function (_calculate_foam_score_020hz) and the same 5‑year sliding‑window event detection.
Identical event CSV/summary files as Experiment 1 (voyager_1_foam_events_020hz.csv, summary_*.txt/.csv/.png).


Adds a new “spacetime‑foam / GW candidate” detector (detect_spacetime_foam).
• Computes Welch PSDs for MAG and PWS in a moving 7‑day window. <br>• Applies five criteria from the Master‑Equation (peak SNR, cross‑layer coherence, clock‑frequency exclusion, scaling test, clock‑catalog distance). <br>• Produces a final candidate list (spacetime_foam_candidates_020hz.csv).
A list of potential gravitational‑wave / foam signatures that satisfy the full 5‑criterion pipeline.


Produces visualisations (2‑D path map, frequency‑band analysis, 5‑year summaries).
Uses Matplotlib to generate PNGs that can be inspected locally or uploaded to the Foam‑Modeling repo for comparison.
voyager_1_2d_foam_map_020hz.png, voyager_1_frequency_bands.png, summary_*.png.


Why it matters for the Foam‑Modeling repo

Unified data set – The concatenated DataFrame (self.data) is the canonical Voyager PWS/MAG time‑series that the Foam‑Detector expects as input. You can simply import the CSV produced by this experiment into the detection pipeline (telemetry_loader.load_telemetry).  
Spacetime‑foam candidates – The candidate CSV (spacetime_foam_candidates_020hz.csv) is a pre‑filtered list of frequencies that already passed the first four criteria (SNR, coherence, clock‑exclusion, scaling). The Resonance Engine can now focus only on those narrow bands when building the matched‑filter bank, dramatically speeding up the Signal Processor.  
Diagnostics & synthetic fallback – Even if you later add new Voyager data (e.g., a newly released PDS volume), the diagnostics will automatically tell you whether the loader can handle the format. If not, the script creates a synthetic Voyager mission that still lets the foam‑modeling visualisations run, ensuring the repository remains self‑contained.  
Visualization assets – The PNGs generated here are perfect for the Foam‑Modeling repo’s README or for a Jupyter notebook that compares the observed Voyager foam events with the predicted foam field from solar_foam_animation.py.


🔗 How the three experiments feed the Foam‑Modeling pipeline



Foam‑Modeling component
Input it needs
Which Voyager experiment supplies it



Baseline foam grid (baseline_foam_map.py)
None – creates a static 3‑D foam matrix.
Independent (shared by all experiments).


Time‑varying foam field (solar_foam_animation.py)
Positions & masses of Solar‑System bodies (already hard‑coded).
Independent – provides the reference foam field that we compare Voyager events against.


Spacecraft trajectory visualiser (craft_trajectory_plotter.py)
CSV/JSON with spacecraft state (time, position, velocity, attitude).
Experiment 1 and Experiment 3 already output Voyager event times & distances; you can convert those into a trajectory file (e.g., voyager1_trajectory.csv) and feed it to the plotter so the Voyager path is shown together with the synthetic foam.


Telemetry for detection (telemetry_loader in the detection pipeline)
Raw PWS/MAG time series (timestamp, frequency, phase, range, SNR).
Experiment 2 (streaming scores) – gives a compact, line‑by‑line CSV that can be read directly. <br>Experiment 3 – its consolidated DataFrame (self.data) can be exported to voyager_full_telemetry.csv and used as the input for the detection stage.


Foam‑event catalog (used as “ground truth” for validation)
List of timestamps, frequencies, scores, distances.
Experiment 1 (voyager_1_foam_events_020hz.csv) and Experiment 3 (voyager_1_foam_events_020hz.csv, identical output).


Spacetime‑foam candidate list (narrow‑band seeds for the Resonance Engine)
Frequency, PSD, coherence, distance for each candidate.
Experiment 3 (spacetime_foam_candidates_020hz.csv).


Visualization assets (PNG/summary)
Static images for documentation or notebooks.
Experiment 3 (2‑D map, frequency‑band plots, 5‑year summaries).


In practice the workflow looks like this:

Run Experiment 3 → obtain a single CSV that contains every cleaned Voyager observation (voyager_full_telemetry.csv) and the event catalog (voyager_1_foam_events_020hz.csv).  
Feed the telemetry CSV to the Detection pipeline (run_foam_detection.py – the conceptual driver described in the Foam‑Modeling repo).  
The pipeline uses the candidate frequencies (spacetime_foam_candidates_020hz.csv) to initialise the Resonance Engine, then runs the matched‑filter on the telemetry.  
The posterior foam map (foam_post_XXXXX.vts) is written and can be visualised together with the Voyager trajectory (generated from the event catalog or from the original PWS ephemeris) using craft_trajectory_plotter.py.  
The final Detection Report (PDF/HTML) can embed the PNGs produced by Experiment 3, giving a side‑by‑side view of observed foam events vs. modelled foam distortions.














