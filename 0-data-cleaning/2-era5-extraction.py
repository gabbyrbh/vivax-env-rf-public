# ------------------------------------------------------------------------------
# @Organization: Benjamin-Chung Lab
# @Project     : Environmental Risk Factors for Vivax Malaria in Peru
# @File        : 2-era5-extraction.py
#
# *** NOTE FOR PUBLIC REPOSITORY ***
# This script is provided for reference only and CANNOT be run from the
# public repository. It requires:
#   (1) Community centroid coordinates (lat/long) which are not published
#       for participant confidentiality reasons.
#   (2) A Google Earth Engine (GEE) service account key with access to a
#       GCP project authorized for the Earth Engine API.
#
# If you wish to adapt this script for your own data, you will need to:
#   - Supply your own community centroid CSV with comm_id, latitude, longitude
#   - Create a GCP project, enable the Earth Engine API, and generate a
#     service account key (see https://developers.google.com/earth-engine)
#   - Update PROJECT_ID and SERVICE_ACCOUNT_FILENAME in the CONFIGURATION
#     section below
#
# The public analysis pipeline begins at step 1-data-processing/ using the
# pre-extracted weather data already merged into the public dataset at:
#   6-public-data/output/vivax-env-erf-public.csv
# ***********************************
#
# @Description : Extract ERA5-Land hourly weather for study community centroids,
#                aggregate to daily summaries, and save CSVs by 6-month period.
#
#                Outputs to:
#                  raw-era5-data/
#                  comm_era5_daily_{year}_1.csv  (January–June)
#                  comm_era5_daily_{year}_2.csv  (July–December)
#
#                Columns: comm_id, date, temperature_2m_mean, temperature_2m_max,
#                         temperature_2m_min, total_precipitation_sum
#                Units:   temperatures in Kelvin; precipitation in meters
#                         (unit conversion to °C and mm is done in
#                          1-data-processing/1-aggregate-era5.R)
# ------------------------------------------------------------------------------

import ee
import os
import json
from pathlib import Path
from datetime import datetime
import pandas as pd


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PROJECT_ID               = "YOUR_GCP_PROJECT_ID"        # replace with your GCP project ID
SERVICE_ACCOUNT_FILENAME = "YOUR_SERVICE_ACCOUNT.json"   # replace with your service account key filename
RESOLUTION               = 11132   # meters (~0.1 degree, ERA5-Land native)
ERA5_BANDS               = ["temperature_2m", "total_precipitation"]
START_YEAR               = 2016    # 2016 needed for lagged analysis of 2017 data
END_YEAR                 = 2024    # last complete year of outcome data


def get_box_root():
    """Find Box directory across different systems."""
    possible_paths = [
        Path.home() / "Library/CloudStorage/Box-Box",  # macOS
        Path.home() / "Box",                            # Windows/Linux
        Path.home() / "Box Sync",                       # Older Box
    ]
    for path in possible_paths:
        if path.exists():
            return str(path)
    raise RuntimeError(
        "Could not find Box directory. Expected one of:\n"
        "  ~/Library/CloudStorage/Box-Box\n"
        "  ~/Box\n"
        "  ~/Box Sync"
    )


BOX_ROOT             = get_box_root()
BOX_FLAME_ERF        = os.path.join(BOX_ROOT, "FLAME-Env-Risk-factors")
CREDENTIALS_DIR      = os.path.join(BOX_FLAME_ERF, "credentials")
SERVICE_ACCOUNT_KEY  = os.path.join(CREDENTIALS_DIR, SERVICE_ACCOUNT_FILENAME)
CENTROIDS_CSV        = os.path.join(BOX_FLAME_ERF, "centroids_final.csv")
OUTPUT_DIR           = os.path.join(BOX_FLAME_ERF, "raw-era5-data")


# ==============================================================================
# EARTH ENGINE INITIALIZATION
# ==============================================================================

def initialize_earth_engine():
    """Initialize Earth Engine with service account authentication."""

    if not os.path.exists(BOX_FLAME_ERF):
        raise FileNotFoundError(
            f"FLAME-Env-Risk-factors folder not found: {BOX_FLAME_ERF}\n\n"
            f"Ensure Box is synced and the folder exists under: {BOX_ROOT}"
        )

    if not os.path.exists(SERVICE_ACCOUNT_KEY):
        raise FileNotFoundError(
            f"Service account key not found: {SERVICE_ACCOUNT_KEY}\n\n"
            f"Contact lab administrator for the key file '{SERVICE_ACCOUNT_FILENAME}'\n"
            f"and save it to: {CREDENTIALS_DIR}"
        )

    with open(SERVICE_ACCOUNT_KEY, "r") as f:
        service_account_info = json.load(f)

    credentials = ee.ServiceAccountCredentials(
        service_account_info["client_email"],
        SERVICE_ACCOUNT_KEY
    )
    ee.Initialize(credentials, project=PROJECT_ID)

    print(f"Earth Engine initialized")
    print(f"  Service Account : {service_account_info['client_email']}")
    print(f"  Project         : {PROJECT_ID}\n")


# ==============================================================================
# DATA LOADING
# ==============================================================================

def load_community_centroids():
    """
    Load community centroid coordinates from Box.

    Reads centroids_final.csv generated by 0-data-cleaning/1-all-district-centroids.R.
    Expected columns: comm_id, long, lat

    Returns:
        DataFrame with columns: comm_id (zero-padded 3-digit string), long (float), lat (float)
    """
    if not os.path.exists(CENTROIDS_CSV):
        raise FileNotFoundError(
            f"Community centroids file not found: {CENTROIDS_CSV}\n\n"
            f"Run 0-data-cleaning/1-all-district-centroids.R first to generate it."
        )

    df = pd.read_csv(CENTROIDS_CSV)

    required = {"comm_id", "long", "lat"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Centroids CSV is missing expected columns: {missing}\n"
            f"Found columns: {list(df.columns)}"
        )

    df = df[["comm_id", "long", "lat"]].dropna().reset_index(drop=True)
    df["comm_id"] = df["comm_id"].apply(lambda x: f"{int(x):03d}")

    print(f"Loaded {len(df)} community centroids\n")
    return df


# ==============================================================================
# ERA5 EXTRACTION
# ==============================================================================

def aggregate_hourly_to_daily(start_date, end_date):
    """
    Aggregate ERA5-Land hourly data to daily summaries server-side.

    For each day computes:
      - temperature_2m_mean, temperature_2m_max, temperature_2m_min  (Kelvin)
      - total_precipitation_sum                                        (meters)

    Args:
        start_date: Start date string (YYYY-MM-DD), inclusive
        end_date:   End date string (YYYY-MM-DD), inclusive

    Returns:
        ee.ImageCollection of daily aggregated images, one image per day
    """
    era5_hourly = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY").select(ERA5_BANDS)

    start   = ee.Date(start_date)
    end_excl = ee.Date(end_date).advance(1, "day")
    n_days  = end_excl.difference(start, "day")

    def make_daily_image(i):
        day      = start.advance(i, "day")
        next_day = day.advance(1, "day")
        daily    = era5_hourly.filterDate(day, next_day)

        temp       = daily.select("temperature_2m")
        temp_mean  = temp.mean().rename("temperature_2m_mean")
        temp_max   = temp.max().rename("temperature_2m_max")
        temp_min   = temp.min().rename("temperature_2m_min")

        precip_sum = daily.select("total_precipitation").sum().rename("total_precipitation_sum")

        return (temp_mean
                .addBands([temp_max, temp_min, precip_sum])
                .set("date", day.format("YYYY-MM-dd")))

    days = ee.List.sequence(0, n_days.subtract(1))
    return ee.ImageCollection(days.map(make_daily_image))


def extract_values_for_points(df_centroids, daily_images, print_every=10):
    """
    Extract daily ERA5 values at each community centroid.

    Args:
        df_centroids:  DataFrame with comm_id, long, lat
        daily_images:  ee.ImageCollection of daily aggregated images
        print_every:   Print progress every N days

    Returns:
        DataFrame with columns: comm_id, date, temperature_2m_mean,
        temperature_2m_max, temperature_2m_min, total_precipitation_sum
        Units are unchanged from ERA5 (Kelvin, meters).
    """
    features = [
        ee.Feature(
            ee.Geometry.Point(row["long"], row["lat"]),
            {"comm_id": row["comm_id"]}
        )
        for _, row in df_centroids.iterrows()
    ]
    fc = ee.FeatureCollection(features)

    daily_list = daily_images.toList(daily_images.size())
    n_days     = int(daily_images.size().getInfo())

    rows = []

    for i in range(n_days):
        if i % print_every == 0 or i == n_days - 1:
            print(f"    Day {i + 1}/{n_days}...")

        try:
            img  = ee.Image(daily_list.get(i))
            date = img.get("date").getInfo()

            result       = img.reduceRegions(fc, reducer=ee.Reducer.first(), scale=RESOLUTION)
            features_out = result.getInfo().get("features", [])

            for feature in features_out:
                props = feature.get("properties", {})
                rows.append({
                    "comm_id":                  props.get("comm_id"),
                    "date":                     date,
                    "temperature_2m_mean":      props.get("temperature_2m_mean"),
                    "temperature_2m_max":       props.get("temperature_2m_max"),
                    "temperature_2m_min":       props.get("temperature_2m_min"),
                    "total_precipitation_sum":  props.get("total_precipitation_sum"),
                })

        except Exception as e:
            print(f"    Warning: failed on day {i + 1} ({date if 'date' in dir() else '?'}): {e}")
            continue

    return pd.DataFrame(rows)


# ==============================================================================
# PERIOD GENERATION
# ==============================================================================

def generate_6month_periods():
    """
    Generate 6-month periods from START_YEAR to today.

    Returns:
        List of tuples: (start_date, end_date, part_number, year)
        part_number = 1 for Jan-Jun, 2 for Jul-Dec
        End date is capped at today for the current in-progress period.
    """
    periods = []
    for year in range(START_YEAR, END_YEAR + 1):
        # H1: January–June
        periods.append((f"{year}-01-01", f"{year}-06-30", 1, year))

        # H2: July–December
        periods.append((f"{year}-07-01", f"{year}-12-31", 2, year))

    return periods


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("ERA5 EXTRACTION — Vivax Malaria Environmental Risk Factors")
    print(f"  Box folder  : {BOX_FLAME_ERF}")
    print(f"  Centroids   : {CENTROIDS_CSV}")
    print(f"  Output dir  : {OUTPUT_DIR}")
    print(f"  Years       : {START_YEAR} – {END_YEAR}\n")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    initialize_earth_engine()
    df_centroids = load_community_centroids()

    periods = generate_6month_periods()
    print(f"Processing {len(periods)} half-year periods\n")

    for start_date, end_date, part, year in periods:
        out_file = os.path.join(OUTPUT_DIR, f"comm_era5_daily_{year}_{part}.csv")

        if os.path.exists(out_file):
            print(f"Skipping {year} part {part} — file already exists: {out_file}")
            continue

        half = "Jan–Jun" if part == 1 else "Jul–Dec"
        print(f"Processing {year} part {part} ({half}): {start_date} → {end_date}")

        daily_data = aggregate_hourly_to_daily(start_date, end_date)
        df_period  = extract_values_for_points(df_centroids, daily_data)

        df_period.to_csv(out_file, index=False)
        print(f"  Saved {len(df_period):,} rows → {out_file}\n")

    print("EXTRACTION COMPLETE")


if __name__ == "__main__":
    main()
