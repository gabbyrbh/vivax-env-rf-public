# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : 2-merge-weather.R
# @Description  : Merge weekly ERA5 weather summaries with community-level
#                 P. vivax incidence data; create continuous study-week index;
#                 save non-lagged-analysis-data_ext.RDS to Box.
# ------------------------------------------------------------------------------

# Clear environment and load configuration
rm(list = ls())
source(paste0(here::here(), "/0-config.R"))

# Load required libraries
library(readxl)    # For reading Excel files
library(lubridate) # For date manipulation
library(ISOweek)   # For ISO week calculations

# Avoid namespace conflicts with select function
select <- dplyr::select

# Define functions -------------------------------------------------------

#' Process Weather Data
#' 
#' Reads and processes ERA5 weekly weather data for a specific year
#' 
#' @param file_year Year of data to process (numeric)
#' @param box_path_flame_erf File path to data directory
#' @return Processed weather dataframe with standardized year and week columns
process_weather_data <- function(file_year, box_path_flame_erf) {
  # Read the ERA5 weekly aggregated weather data
  weather_data <- readRDS(paste0(box_path_flame_erf, "era5-weekly-agg/era_", file_year, "_weekly.RDS"))
  
  # Standardize data types for merging
  processed_data <- weather_data %>% 
    mutate(year = as.character(year),    # Convert year to character for consistency
           week = as.numeric(week))      # Ensure week is numeric
  
  return(processed_data)
}

# Load and aggregate weather data ----------------------------------------

# Define study period years
years <- c(2016:2024)

# Process and combine weather data for all study years
df_all_years <- map_dfr(years, ~process_weather_data(file_year = .x, box_path_flame_erf))

# Aggregate weather variables by community, year, and week
# Take mean values where multiple observations exist
df_all_years <- df_all_years %>% 
  group_by(comm_id, year, week) %>%
  summarize(across(c("temp_wk_min", "temp_wk_max", "precip_wk_total", "oni_index"), mean)) %>%
  mutate(week = as.numeric(week)) %>% 
  ungroup() %>%
  # Calculate number of weeks per year for study week indexing
  group_by(year) %>%
  mutate(num_weeks = length(unique(week))) %>%
  ungroup()

# Create study week index ------------------------------------------------

# Calculate cumulative weeks to create continuous study week variable
# This accounts for varying numbers of weeks per year
study_week_df <- df_all_years %>%
  ungroup() %>%
  group_by(year) %>%
  mutate(num_weeks = length(unique(week))) %>% 
  ungroup() %>%
  select(year, num_weeks) %>%
  filter(year != "2016") %>%  # Exclude 2016 as reference year
  distinct() %>%
  mutate(cum_weeks = cumsum(num_weeks)) %>%
  mutate(cum_weeks = lag(cum_weeks)) %>%  # Lag to get previous year cumulative weeks
  select(year, num_weeks, cum_weeks)

# Add study week variable to weather data
# Study week provides continuous numbering across all study years
weather_df <- df_all_years %>%
  left_join(study_week_df) %>%
  mutate(study_week = case_when(
    year == "2016" ~ week - num_weeks,  # 2016 gets negative values (baseline)
    year == "2017" ~ week,              # 2017 starts at week 1
    TRUE ~ week + cum_weeks             # Subsequent years continue sequence
  ))

# Load incidence data ----------------------------------------------------

# Load and process weekly incidence data
weekly_inc_agg = readRDS(paste0(box_path_flame_erf, "all_districts_weekly_incidence.RDS")) %>%
  mutate(week = as.numeric(week),
         year = as.character(year)) %>%
  # Apply same study week calculation as weather data
  left_join(study_week_df) %>%
  mutate(study_week = case_when(
    year == "2016" ~ week - num_weeks,
    year == "2017" ~ week,
    TRUE ~ week + cum_weeks
  )) %>%
  ungroup() %>%
  # Remove geographic identifier columns not needed for analysis
  select(-PROVINCIA, -DISTRITO, -localidades_inei)

# Get list of communities with incidence data
inc_comms <- weekly_inc_agg %>% 
  pull(comm_id) %>% 
  unique()

# Merge weather and incidence data ---------------------------------------

# Identify communities with weather data available
weather_comms <- weather_df %>%
  pull(comm_id) %>%
  unique()

# Flag communities with weather data in incidence dataset
weekly_inc_agg <- weekly_inc_agg %>%
  mutate(weather_comm = ifelse(comm_id %in% weather_comms, 1, 0))

# Identify communities with missing weather data for quality assessment
missing_data <- anti_join(weekly_inc_agg, weather_df, by = c("comm_id", "year", "week", "study_week"))

# Visualize temporal pattern of missing data
ggplot(missing_data, aes(x = monday_date, y = n_cases)) + geom_point()

# Perform left join to combine incidence and weather data
# This preserves all incidence observations, adding weather data where available
combined_study_week_df <- left_join(weekly_inc_agg, 
                                    weather_df,
                                    by = join_by(comm_id, year, week, num_weeks, cum_weeks, study_week))

# Save combined dataset for analysis
saveRDS(combined_study_week_df, paste0(box_path_flame_erf, "non-lagged-analysis-data_ext.RDS"))

# Data quality checks ----------------------------------------------------

# Analyze incidence patterns in low-population communities
# Small populations may have unstable incidence rates due to small denominators
low_pop_inc <- combined_study_week_df %>%
  filter(population < 15) %>%
  mutate(inc = n_cases/population) %>%  # Calculate incidence rate
  group_by(year) %>%
  summarize(mean_inc_10k = mean(inc)*10000,           # Mean incidence per 10,000
            num_communities = length(unique(comm_id))) # Number of communities

# Analyze incidence patterns in high-population communities for comparison
high_pop_inc <- combined_study_week_df %>%
  filter(population > 15) %>%
  mutate(inc = n_cases/population) %>%
  group_by(year) %>%
  summarize(mean_inc_10k = mean(inc)*10000,
            num_communities = length(unique(comm_id)))