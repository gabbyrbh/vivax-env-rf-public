# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : 1-aggregate-era5.R
# @Description  : Convert ERA5-Land daily data (Kelvin/meters) to Celsius/mm,
#                 aggregate to weekly summaries, and save era_{year}_weekly.RDS
#                 to Box. Also builds the ONI index and saves mean thresholds.
# ------------------------------------------------------------------------------
rm(list=ls())
source(paste0(here::here(), "/0-config.R"))
library(mgcv)
library(lme4)
library(lubridate)
library(zoo)
select <- dplyr::select
summarize <- dplyr::summarize

# combine 2021, 2022, and 2023 -- had to be split bc of processing timeouts in jupyter notebook
# Load and combine raw daily ERA5 data ------------------------------------
combine_era_parts = function(year) {
  part_1 = read.csv(paste0(box_path_flame_erf, "raw-era5-data/comm_era5_daily_", year, "_1.csv"))
  part_2 = read.csv(paste0(box_path_flame_erf, "raw-era5-data/comm_era5_daily_", year, "_2.csv"))
  full = rbind(part_1, part_2)
  full = full %>% 
    mutate(data_year = year(date)) %>%
    rename("year" = "data_year")
  return(full)
}

years <- c(2016:2024)
df_all_years <- map_dfr(years, ~combine_era_parts(.x))

for (this_year in years) {
  this_year_df <- df_all_years %>%
    filter(year == this_year)
  write.csv(this_year_df, paste0(box_path_flame_erf, "raw-era5-data/comm_era5_daily_", this_year, ".csv"))
}

# Build ONI index ---------------------------------------------------------
# manually create El Niño dataset -- from https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php
oni_df <- data.frame(
  year = 2010:2024,
  Jan = c(1.5, -1.4, -0.9, -0.4, -0.4, 0.5, 2.5, -0.3, -0.9, 0.7, 0.5, -1.0, -1.0, -0.7, 1.8),
  Feb = c(1.2, -1.2, -0.7, -0.4, -0.5, 0.5, 2.1, -0.2, -0.9, 0.7, 0.5, -0.9, -0.9, -0.4, 1.5),
  Mar = c(0.8, -0.9, -0.6, -0.3, -0.3, 0.5, 1.6, 0.1, -0.7, 0.7, 0.4, -0.8, -1.0, -0.1, 1.1),
  Apr = c(0.4, -0.7, -0.5, -0.3, 0.0, 0.7, 0.9, 0.2, -0.5, 0.7, 0.2, -0.7, -1.1, 0.2, 0.7),
  May = c(-0.2, -0.6, -0.3, -0.4, 0.2, 0.9, 0.4, 0.3, -0.2, 0.5, -0.1, -0.5, -1.0, 0.5, 0.4),
  Jun = c(-0.7, -0.4, 0.0, -0.4, 0.2, 1.2, -0.1, 0.3, 0.0, 0.5, -0.3, -0.4, -0.9, 0.8, 0.2),
  Jul = c(-1.0, -0.5, 0.2, -0.4, 0.0, 1.5, -0.4, 0.1, 0.1, 0.3, -0.4, -0.4, -0.8, 1.1, 0.0),
  Aug = c(-1.3, -0.6, 0.4, -0.3, 0.1, 1.9, -0.5, -0.1, 0.2, 0.1, -0.6, -0.5, -0.9, 1.3, -0.1),
  Sep = c(-1.6, -0.8, 0.4, -0.3, 0.2, 2.2, -0.6, -0.4, 0.5, 0.2, -0.9, -0.7, -1.0, 1.6, -0.2),
  Oct = c(-1.6, -1.0, 0.3, -0.2, 0.5, 2.4, -0.7, -0.7, 0.8, 0.3, -1.2, -0.8, -1.0, 1.8, -0.3),
  Nov = c(-1.6, -1.1, 0.1, -0.2, 0.6, 2.6, -0.7, -0.8, 0.9, 0.5, -1.3, -1.0, -0.9, 1.9, -0.4),
  Dec = c(-1.6, -1.0, -0.2, -0.3, 0.7, 2.6, -0.6, -1.0, 0.8, 0.5, -1.2, -1.0, -0.8, 2.0, -0.5)
)

oni_df_long <- oni_df %>%
  pivot_longer(cols = -year, 
               names_to = "month_abb", 
               values_to = "oni_index") %>%
  mutate(month = match(month_abb, month.abb)) %>%
  select(year, month, oni_index) %>%
  arrange(year, month)


# Aggregate to weekly summaries and save ---------------------------------
process_era5 = function(year) {
  df = read.csv(paste0(box_path_flame_erf, "raw-era5-data/comm_era5_daily_", year, ".csv"))
  
  df = df %>% mutate(total_precipitation_sum_mm = total_precipitation_sum * 1000,
                     across(c(temperature_2m_max, temperature_2m_min, temperature_2m_mean), ~ .x - 273.15))
  df = df %>% mutate_at(vars(everything()), ~ replace(., . < 0, NA)) # replace neg values as missing
  
  # calculate daily temp range
  df = df %>% mutate(temp_day_range = temperature_2m_max - temperature_2m_min)
  
  # calculate percentile daily precipitation for binary heavy rain vars
  precip_90 = quantile(df$total_precipitation_sum_mm, probs = 0.9, na.rm = T)
  precip_10 = quantile(df$total_precipitation_sum_mm, probs = 0.1, na.rm = T)
  
  heat_90 = quantile(df$temperature_2m_mean, probs = 0.90, na.rm = T)

    # aggregate weekly
  df = df %>% mutate(week = lubridate::isoweek(date),
                     month = lubridate::month(date))
  
  # generate weekly measures
  df_proc = df %>% 
    group_by(comm_id, year, month, week) %>% 
    mutate(comm_id = sprintf("%03d", comm_id)) %>%
    summarize(temp_wk_mean = mean(temperature_2m_mean),
              temp_wk_min = min(temperature_2m_min),
              temp_wk_max = max(temperature_2m_max),
              temp_wk_mean_range = mean(temp_day_range),
              heat_wk_high = ifelse(any(temperature_2m_mean > heat_90), 1, 0),
              precip_wk_max = max(total_precipitation_sum_mm),
              precip_wk_total = sum(total_precipitation_sum_mm),
              precip_wk_mean = mean(total_precipitation_sum_mm),
              precip_wk_heavy = ifelse(any(total_precipitation_sum_mm > precip_90), 1, 0))
  
  # add in El Nino
  df_proc <- df_proc %>%
    left_join(oni_df_long, by = c("year", "month"))
  
  saveRDS(df_proc, paste0(box_path_flame_erf, "era5-weekly-agg/era_", year, "_weekly.RDS"))
  
  precip_df = data.frame(year = year, type = "precip", threshold = precip_90)
  temp_df = data.frame(year = year, type = "temp", threshold = heat_90)
  threshold_df = rbind(precip_df, temp_df)
  
  return(threshold_df)
}

  
thresholds <- map_dfr(.x = years, ~ process_era5(.x))

mean_thresholds <- thresholds %>% group_by(type) %>% summarize(mean_threshold = mean(threshold))
saveRDS(mean_thresholds, paste0(data_path, "mean-thresholds.RDS"))

# Check for negative values (proxy for missing/corrupt data) -------------
check_missing <- function(year) {
  df = read.csv(paste0(box_path_flame_erf, "raw-era5-data/comm_era5_daily_", year, ".csv"))
  # process temp & precip into C and mm, respectively
  df = df %>% mutate(total_precipitation_sum_mm = total_precipitation_sum * 1000,
                     across(c(temperature_2m_max, temperature_2m_min, temperature_2m_mean), ~ .x - 273.15))
  missing_df <- df %>% 
    filter(if_any(everything(), ~ . < 0))
  
  return(missing_df)
}

missing_data <- map_dfr(.x = years, ~ check_missing(.x))
        
