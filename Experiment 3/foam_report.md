# Foam Score Report

**Source file:** `anomalies_full__streaming_scores.csv.gz`  
**File size on disk:** 2.52 GB  
**Rows scored:** 14,429,424  
**Processed in:** 240.5s

## Coverage

- Input row numbers span **1 -> 14,429,424** (these are the original line numbers from `anomalies_full.csv`).
- Rows with a usable distance value: **14,417,750**; rows with missing distance (the PARTIAL rows): **11,674**.
- Timestamp range (as stored): `1970-01-01 00:00:00.240940722` -> `1970-01-01 00:00:00.287476`.
  - Note: 14,429,424 rows carry a `1970-01-01` (Unix-epoch) timestamp, i.e. the original time field did not parse into a real date for those rows.

## Foam score distribution

| Statistic | Value |
|---|---|
| count | 14,429,424 |
| min | 0.000000 |
| max | 0.215273 |
| mean | 0.001220 |
| std dev | 0.001201 |
| 1st pct | 0.000048 |
| 5th pct | 0.000241 |
| 25th pct | 0.001061 |
| median (50th) | 0.001419 |
| 75th pct | 0.001777 |
| 90th pct | 0.001992 |
| 95th pct | 0.003700 |
| 99th pct | 0.005532 |
| 99.9th pct | 0.005959 |

### Rows at or above each score threshold

| Threshold | Rows >= | % of scored |
|---|---|---|
| 0.01 | 1,453 | 0.0101% |
| 0.05 | 327 | 0.0023% |
| 0.10 | 86 | 0.0006% |
| 0.20 | 2 | 0.0000% |
| 0.30 | 0 | 0.0000% |
| 0.40 | 0 | 0.0000% |
| 0.50 | 0 | 0.0000% |
| 0.60 | 0 | 0.0000% |
| 0.70 | 0 | 0.0000% |
| 0.80 | 0 | 0.0000% |
| 0.90 | 0 | 0.0000% |
| 0.95 | 0 | 0.0000% |
| 0.99 | 0 | 0.0000% |

## Frequency (center_freq_hz)

- Range: **-9999.000 Hz -> 177,827,941,003,892,282,464,666,329,044,544,350,433,461,841,993,555,425,302,819,213,811,781,960,363,881,908,293,009,773,095,947,472,054,642,624,823,151,102,215,353,416,768,185,912,439,251,732,020,762,732,641,835,392,185,020,240,736,560,040,281,759,666,003,866,472,845,087,673,665,468,854,956,796,323,669,433,933,121,607,641,259,702,987,255,869,961,538,025,987,457,967,480,249,004,057,700,357,996,251,019,476,992 Hz**, mean **inf Hz** (over 14,429,424 rows with a finite frequency).
- Within +/-10% of 8 Hz: **0** rows; within +/-10% of 20 Hz: **0** rows.

| Band | Rows |
|---|---|
| < 0 Hz | 2,689,920 |
| 0-8 Hz | 1 |
| 8-20 Hz | 52,792 |
| 20-100 Hz | 29,787 |
| 100 Hz-1 kHz | 772,744 |
| 1-10 kHz | 1,298,815 |
| > 10 kHz | 9,585,365 |

## Amplitude

- Range: **5.941e-07 -> 1.500e+01** (over 14,429,424 positive, finite values).
- Geometric mean (typical order of magnitude): **1.769e-01**.

## Distance (distance_au)

- Over the 14,417,750 rows that have a value: **1.01 -> 5.13 AU**, mean **4.54 AU**.

## Rows by spacecraft

| Spacecraft | Rows | % |
|---|---|---|
| voyager1 | 13,844,765 | 95.95% |
| voyager2 | 584,658 | 4.05% |
| unknown | 1 | 0.00% |

## Top 25 sources

| Source | Rows |
|---|---|
| vg1_pws_wf_1978-10-10T10_v1.0.cdf | 733,190 |
| vg1_pws_wf_1978-12-13T21_v1.0.cdf | 511,033 |
| vg1_pws_wf_1979-01-18T21_v1.0.cdf | 264,155 |
| vg1_pws_wf_1979-01-21T20_v1.0.cdf | 239,780 |
| vg1_pws_wf_1979-01-07T23_v1.0.cdf | 212,779 |
| vg1_pws_wf_1979-01-14T22_v1.0.cdf | 192,784 |
| vg1_pws_wf_1979-01-15T21_v1.0.cdf | 167,128 |
| vg1_pws_wf_1979-01-25T20_v1.0.cdf | 159,990 |
| vg1_pws_wf_1979-01-26T20_v1.0.cdf | 156,618 |
| vg1_pws_wf_1979-01-16T21_v1.0.cdf | 155,273 |
| vg1_pws_wf_1979-01-22T20_v1.0.cdf | 142,950 |
| vg1_pws_wf_1979-02-09T10_v1.0.cdf | 134,693 |
| vg1_pws_wf_1979-01-20T21_v1.0.cdf | 134,365 |
| vg1_pws_wf_1979-01-12T22_v1.0.cdf | 133,463 |
| vg1_pws_wf_1979-02-06T23_v1.0.cdf | 133,167 |
| vg1_pws_wf_1978-09-26T12_v1.0.cdf | 132,368 |
| vg1_pws_wf_1979-01-10T22_v1.0.cdf | 131,741 |
| vg1_pws_wf_1978-09-26T19_v1.0.cdf | 131,732 |
| vg1_pws_wf_1979-01-28T19_v1.0.cdf | 130,108 |
| vg1_pws_wf_1979-02-07T01_v1.0.cdf | 129,418 |
| vg1_pws_wf_1979-01-11T22_v1.0.cdf | 128,998 |
| vg1_pws_wf_1979-01-24T20_v1.0.cdf | 128,125 |
| vg1_pws_wf_1979-01-17T21_v1.0.cdf | 123,960 |
| vg1_pws_wf_1979-01-08T23_v1.0.cdf | 119,256 |
| vg1_pws_wf_1979-01-19T21_v1.0.cdf | 117,820 |

_(1,165 distinct sources in total; showing the 25 largest.)_

## Top 500 highest-scoring rows

Exported in full to `foam_top_rows.csv`. The 15 highest:

| foam_score | freq_hz | amplitude | distance_au | spacecraft | source | input_row |
|---|---|---|---|---|---|---|
| 0.215273 | 17.799999237060547 | 0.003182461252436 | nan | voyager2 | vg2pws_lr_19770820_v5.30.cdf | 175 |
| 0.215273 | 17.799999237060547 | 0.0806223452091217 | nan | voyager1 | vg1pws_lr_19770905_v5.30.cdf | 7302 |
| 0.157269 | 17.799999237060547 | 0.1414427161216735 | 3.3098116782462057 | voyager2 | vg2pws_lr_19780619_v5.30.cdf | 2053925 |
| 0.157267 | 17.799999237060547 | 0.1414427161216735 | 3.3059267115787723 | voyager2 | vg2pws_lr_19780618_v5.30.cdf | 2052239 |
| 0.157266 | 17.799999237060547 | 0.1414427161216735 | 3.302696455192963 | voyager2 | vg2pws_lr_19780618_v5.30.cdf | 2050034 |
| 0.157255 | 17.799999237060547 | 0.1414427161216735 | 3.2768821655289737 | voyager2 | vg2pws_lr_19780614_v5.30.cdf | 2030990 |
| 0.157254 | 17.799999237060547 | 0.1414427161216735 | 3.2739757940028853 | voyager2 | vg2pws_lr_19780614_v5.30.cdf | 2029660 |
| 0.157252 | 17.799999237060547 | 0.1414427161216735 | 3.2699586294638765 | voyager2 | vg2pws_lr_19780613_v5.30.cdf | 2026499 |
| 0.157251 | 17.799999237060547 | 0.1414427161216735 | 3.2669613718700137 | voyager2 | vg2pws_lr_19780613_v5.30.cdf | 2024506 |
| 0.157249 | 17.799999237060547 | 0.1414427161216735 | 3.261901333636008 | voyager2 | vg2pws_lr_19780612_v5.30.cdf | 2017841 |
| 0.157239 | 17.799999237060547 | 0.1414427161216735 | 3.2379661929499224 | voyager2 | vg2pws_lr_19780609_v5.30.cdf | 2003306 |
| 0.157237 | 17.799999237060547 | 0.1414427161216735 | 3.232948283395139 | voyager2 | vg2pws_lr_19780608_v5.30.cdf | 2001074 |
| 0.157225 | 17.799999237060547 | 0.1414427161216735 | 3.203828154618242 | voyager2 | vg2pws_lr_19780604_v5.30.cdf | 1981541 |
| 0.157203 | 17.799999237060547 | 0.1414427161216735 | 3.150309856187428 | voyager2 | vg2pws_lr_19780528_v5.30.cdf | 1939432 |
| 0.157147 | 17.799999237060547 | 0.1414427161216735 | 3.017433150756819 | voyager2 | vg2pws_lr_19780510_v5.30.cdf | 1801617 |
