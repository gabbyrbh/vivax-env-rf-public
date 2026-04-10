# vivax-env-rf-public

**Associations between weather and *Plasmodium vivax* malaria in an Amazonian elimination setting: a distributed lag analysis from 2016–2024**

Benjamin-Chung Lab | First author: Gabriella Barratt Heitmann | PI: Jade Benjamin-Chung

---

## Overview

This repository contains the public analysis pipeline for a study examining lagged associations between weekly weather exposures (minimum temperature, maximum temperature, total precipitation) and *P. vivax* malaria incidence across 136 communities in Loreto, Peru. The primary method is distributed lag nonlinear modeling (DLNM) with Poisson/log link, 2–16 week lookback, followed by generalized additive models (GAMs) at the primary lag with bootstrap confidence intervals.

**Study period**: 2016–2024 (2020–2021 excluded due to COVID-19 disruptions)
**Software**: R v4.2.1, Python 3, DLNM v2.4.7, mgcv

> **Note on reproducibility**: Community centroid coordinates (lat/long) are not published for participant confidentiality reasons. The ERA5 weather extraction script (`0-data-cleaning/2-era5-extraction.py`) is provided for reference only and cannot be run from this repository. The public analysis pipeline begins at step `1-data-processing/` using weather data already merged into the public dataset at `6-public-data/output/vivax-env-erf-public.csv`.

---

## Repository Structure

```
vivax-env-rf-public/
├── 0-config.R                        # Shared paths and packages — edit before running
├── 0-base-functions.R                # Shared utility functions
│
├── 0-data-cleaning/                  # Reference only — cannot be run publicly
│   └── 2-era5-extraction.py          # ERA5-Land extraction via GEE (requires coordinates)
│
├── 1-data-processing/                # Requires original Box data — for reference
│   ├── 1-aggregate-era5.R            # Kelvin/meters → °C/mm, daily → weekly
│   └── 2-merge-weather.R             # Merge incidence + weather → analysis dataset
│
├── 2-analysis/                       # Reproducible with public dataset
│   └── 1-fit-dlnm-gam.R             # Primary DLNM-GAMs (produces DLNM slices + IR estimates)
│
├── 3-sherlock/                       # Bootstrap GAMs — run on HPC or locally (slow)
│   ├── gam-temp-min.R / run-gam-temp-min.sh
│   ├── gam-temp-max.R / run-gam-temp-max.sh
│   ├── gam-precip-total.R / run-gam-precip-total.sh
│   ├── gam-{exposure}-enso.R / .sh       (×3)
│   └── gam-{exposure}-comm-type.R / .sh  (×3)
│
├── 4-figure-scripts/                 # Manuscript figures
│   ├── fig2-map-study-area.R         # Requires geographic data — cannot be reproduced publicly
│   ├── fig3-time-series.R
│   ├── fig4-plot-main-gam.R          # Requires bootstrap results from 3-sherlock/
│   └── fig5-plot-sub-groups.R        # Requires bootstrap results from 3-sherlock/
│
├── 5-table-scripts/                  # Manuscript tables
│   └── sub-group-table.R             # Requires bootstrap results from 3-sherlock/
│
├── 6-public-data/                    # Public dataset
│   ├── create-public-dataset.R       # Script used to generate the public dataset
│   └── output/
│       └── vivax-env-erf-public.csv  # Public analysis dataset (start here)
│
├── data/                             # Small intermediate objects
├── results/                          # Analysis outputs (DLNM slices, IR estimates)
├── figures/                          # Generated figures
└── tables/                           # Generated tables
```

---

## Setup

### 1. Configure paths

Open `0-config.R` and set `data_dir` to the directory where you have stored the data files:

```r
data_dir <- "/path/to/your/data/"
```

If you are using only the public dataset (recommended), no external data directory is needed — data are read from `6-public-data/output/` automatically.

### 2. Install R packages

```r
install.packages(c(
  "dplyr", "magrittr", "purrr", "tidyr", "stringr", "ggplot2", "here",
  "glue", "data.table", "lubridate", "forestplot", "grid", "gridExtra",
  "cowplot", "lemon", "foreach", "parallel", "doParallel",
  "mgcv", "lme4", "dlnm", "splines", "glue", "broom", "boot",
  "scales", "cowplot", "gKRLS", "kableExtra"
))
```

---

## Public Dataset

The primary analysis dataset is at **`6-public-data/output/vivax-env-erf-public.csv`**. Community IDs have been randomly permuted and geographic identifiers removed.

### Variables

| Variable | Description |
|----------|-------------|
| `comm_id` | Anonymized community identifier (scrambled) |
| `year` | Calendar year |
| `week` | Epidemiological week within year |
| `monday_date` | Date of Monday starting each study week |
| `n_cases` | Weekly *P. vivax* case count |
| `population` | Community population estimate |
| `years_since_itn_delivery` | Years since insecticide-treated net distribution |
| `num_weeks` | Total weeks of observation for community |
| `cum_weeks` | Cumulative study week number |
| `study_week` | Sequential study week |
| `temp_wk_min` | Weekly minimum temperature (°C, ERA5-Land) |
| `temp_wk_max` | Weekly maximum temperature (°C, ERA5-Land) |
| `precip_wk_total` | Weekly total precipitation (mm, ERA5-Land) |
| `oni_index` | Oceanic Niño Index (ONI) |
| `comm_type` | Community type (`highway` or `riverine`) |

### Reading the public dataset into analysis scripts

Replace the Box data-loading block at the top of each analysis script with:

```r
nonlagged_data <- read.csv(paste0(here::here(), "/6-public-data/output/vivax-env-erf-public.csv")) %>%
  rename(time = "week") %>%
  mutate(population = as.numeric(population),
         n_cases    = as.numeric(n_cases),
         comm_id    = as.factor(comm_id),
         oni_index  = as.numeric(oni_index)) %>%
  filter(!is.na(population)) %>%
  mutate(year = factor(as.character(year),
                       levels = c("2016","2017","2018","2019","2020","2021","2022","2023","2024"))) %>%
  filter(year != "2020") %>%
  filter(year != "2021")
```

---

## Running the Pipeline

### Scripts reproducible with the public dataset

| Script | Output | Notes |
|--------|--------|-------|
| `2-analysis/1-fit-dlnm-gam.R` | DLNM slices, IR estimates | Run locally |
| `3-sherlock/gam-*.R` | GAM base results + bootstrap CIs | See HPC instructions below |
| `4-figure-scripts/fig3-time-series.R` | Figure 3 | Run locally |
| `4-figure-scripts/fig4-plot-main-gam.R` | Figure 4 | Requires bootstrap results |
| `4-figure-scripts/fig5-plot-sub-groups.R` | Figure 5 | Requires bootstrap results |
| `5-table-scripts/sub-group-table.R` | Tables 1 & 2 | Requires bootstrap results |

> **Note:** Because community IDs are scrambled, community-level fixed effect estimates will differ from those in the manuscript. All exposure-response curves, IRs, and subgroup analyses are unaffected.

> **Note:** Figure 2 (study area map) requires community centroid coordinates and cannot be reproduced from the public dataset.

---

## Bootstrap GAM Scripts (3-sherlock/)

The 9 scripts in `3-sherlock/` fit GAMs with 1,000 bootstrap iterations and are computationally intensive (~2–4 hours each with parallelization). They can be run locally or on an HPC cluster.

### Running locally

Each script can be run directly in R. Edit the number of workers at the top of each script to match your machine:

```r
r_workers <- 4  # set to number of available CPU cores - 1
```

Then run:
```r
source("3-sherlock/gam-temp-min.R")
```

Or from the terminal:
```bash
Rscript 3-sherlock/gam-temp-min.R
```

### Running on an HPC (SLURM)

Each `.R` script has a corresponding `.sh` SLURM submission script. To run on a SLURM-based HPC:

**1. Transfer scripts to your HPC:**
```bash
scp 3-sherlock/gam-*.R 3-sherlock/run-gam-*.sh \
  <user>@<hpc-login-node>:<your-project-dir>/
```

**2. Edit the `.sh` files** to match your cluster's account, partition, and module names. Key lines to update:
```bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
module load R/4.2.1   # adjust to your cluster's R module
```

**3. Transfer the data** your cluster needs — at minimum the public dataset:
```bash
scp 6-public-data/output/vivax-env-erf-public.csv \
  <user>@<hpc-login-node>:<your-project-dir>/6-public-data/output/
```

**4. Submit all 9 jobs:**
```bash
sbatch run-gam-temp-min.sh
sbatch run-gam-temp-max.sh
sbatch run-gam-precip-total.sh
sbatch run-gam-temp-min-enso.sh
sbatch run-gam-temp-max-enso.sh
sbatch run-gam-precip-total-enso.sh
sbatch run-gam-temp-min-comm-type.sh
sbatch run-gam-temp-max-comm-type.sh
sbatch run-gam-precip-total-comm-type.sh
```

**5. Retrieve outputs** — each script saves two files to `results/base results/` and `bs/`:
- `base-{predictor}.RDS` — marginal incidence curve from the fitted model
- `bs-{predictor}.RDS` — bootstrap distribution (1,000 iterations × predictor sequence)

Copy these back before running figure and table scripts:
```bash
scp <user>@<hpc-login-node>:<your-project-dir>/results/base\ results/*.RDS \
  results/base\ results/
scp <user>@<hpc-login-node>:<your-project-dir>/bs/*.RDS bs/
```

---

## Manuscript Outputs

| Figure/Table | Script |
|-------------|--------|
| Figure 2 — Study area map | `4-figure-scripts/fig2-map-study-area.R` *(requires geographic data)* |
| Figure 3 — Weather & incidence time series | `4-figure-scripts/fig3-time-series.R` |
| Figure 4 — Main GAM results with bootstrap CIs | `4-figure-scripts/fig4-plot-main-gam.R` |
| Figure 5 — Subgroup results (ENSO, community type) | `4-figure-scripts/fig5-plot-sub-groups.R` |
| Table 1 — Incidence ratios by community type | `5-table-scripts/sub-group-table.R` |
| Table 2 — Incidence ratios by ENSO period | `5-table-scripts/sub-group-table.R` |

---

## Citation

*Manuscript in preparation.* Please contact the authors for citation information.
