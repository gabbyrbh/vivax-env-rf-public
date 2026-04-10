# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : 0-config.R
# @Description  : Shared configuration for all scripts: loads packages, sets
#                 paths, and defines analysis-wide constants.
#
# SETUP: Replace "YOUR_FILE_PATH" below with the path to your local copy of
#        the data directory before running any scripts.
# ------------------------------------------------------------------------------

# Load libraries ----------------------------------------------------------
library(dplyr)
library(magrittr)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(here)
library(glue)
library(data.table)
library(lubridate)
library(forestplot)
library(grid)
library(gridExtra)
library(cowplot)
library(lemon)
library(foreach)
library(parallel)
library(doParallel)

# File paths ------------------------------------------------------------------
# Set this to the directory where you have stored the data files.
# This can be a local folder or a cloud-synced directory (e.g., Box, Dropbox).
# Example: data_dir <- "/Users/yourname/path/to/data/"
data_dir <- "YOUR_FILE_PATH"

# Path to the main project data folder (non-lagged analysis dataset,
# bootstrap results, centroid files, etc.)
box_path_flame_erf <- paste0(data_dir, "FLAME-Env-Risk-factors/")

# Path to shared spatial/raster files — required only for Figure 2
# (land use raster). Not needed if using the public dataset.
box_path_ucsf_shared <- paste0(data_dir, "YOUR_SHARED_FOLDER/")

# To use the public dataset (recommended), data are already in:
public_data_path <- paste0(here::here(), "/6-public-data/output/")

# Local output paths (no changes needed below this line) ----------------------
data_path    <- paste0(here::here(), "/data/")
results_path <- paste0(here::here(), "/results/")
tables_path  <- paste0(here::here(), "/tables/")
figure_path  <- paste0(here::here(), "/figures/")
