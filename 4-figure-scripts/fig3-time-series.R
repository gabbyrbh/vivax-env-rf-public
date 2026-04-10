# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : fig3-time-series.R
# @Description  : Generate Figure 3: weekly weather variables (min/max temp,
#                 total precipitation) and P. vivax incidence time series.
# ------------------------------------------------------------------------------
rm(list=ls())

source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/0-base-functions.R"))
library(mgcv)
library(lme4)
library(lmtest)
library(sandwich)
library(gKRLS)
library(kableExtra)
library(ISOweek)
library(tidyverse)
library(cowplot)

select <- dplyr::select
summarize <- dplyr::summarize

# Load data ---------------------------------------------------------------
nonlagged_data <- read.csv(paste0(public_data_path, "vivax-env-erf-public.csv")) %>%
  rename(time = week) %>%
  mutate(
    population = as.numeric(population),
    n_cases = as.numeric(n_cases),
    comm_id = as.factor(comm_id),
    oni_index = as.numeric(oni_index),
    year = factor(as.character(year), levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))
  ) %>%
  filter(population > 15)

# Set median pop community as reference
median_pop_comm <- nonlagged_data %>%
  filter(year == 2019) %>%
  group_by(comm_id) %>%
  summarize(avg_population = mean(population), .groups = "drop") %>%
  mutate(diff_from_median = abs(avg_population - median(avg_population))) %>%
  slice_min(diff_from_median) %>%
  pull(comm_id) %>% 
  first() %>% 
  as.character()

nonlagged_data <- nonlagged_data %>% 
  mutate(comm_id = relevel(comm_id, ref = median_pop_comm))

# Read in comm_type csv
# comm_type is included in the public dataset
nonlagged_data <- nonlagged_data %>%
  mutate(comm_type = ifelse(is.na(comm_type) | comm_type == "", "riverine", comm_type)) %>%
  filter(!is.na(comm_type)) %>%
  mutate(
    comm_type = factor(comm_type, levels = c("highway", "riverine"), labels = c("Highway", "Riverine")),
    highway  = if_else(comm_type == "Highway", 1, 0),
    riverine = if_else(comm_type == "Riverine", 1, 0)
  ) %>%
  mutate(
    enso_period = case_when(
      oni_index >= 0.5 ~ "El Niño periods",
      oni_index <= -0.5 ~ "La Niña periods",
      TRUE ~ "Neutral ENSO periods"
    ),
    enso_period = factor(enso_period, levels = c("Neutral ENSO periods", "El Niño periods", "La Niña periods"))
  )

# Create ENSO shading data
enso_shading <- nonlagged_data %>%
  select(monday_date, oni_index) %>%
  distinct() %>%
  arrange(monday_date) %>%
  mutate(
    enso_type = case_when(
      oni_index >= 0.5 ~ "El Niño",
      oni_index <= -0.5 ~ "La Niña",
      TRUE ~ "Neutral"
    )
  ) %>%
  filter(enso_type != "Neutral") %>%
  # Create groups for consecutive periods
  mutate(
    group = cumsum(enso_type != lag(enso_type, default = first(enso_type)))
  ) %>%
  group_by(group, enso_type) %>%
  summarize(
    xmin = min(monday_date),
    xmax = max(monday_date),
    .groups = "drop"
  ) %>%
  mutate(
    fill_color = case_when(
      enso_type == "El Niño" ~ "orange", 
      enso_type == "La Niña" ~ "cornflowerblue"
    )
  )

df_long <- nonlagged_data %>% 
  select(monday_date, comm_id, n_cases, temp_wk_min, temp_wk_max, precip_wk_total, enso_period) %>% 
  group_by(monday_date) %>%
  summarize(across(c(temp_wk_min, temp_wk_max, precip_wk_total), mean), .groups = "drop") %>% 
  pivot_longer(
    cols = c(temp_wk_min, temp_wk_max, precip_wk_total), 
    names_to = "measure", 
    values_to = "value"
  ) %>%
  mutate(plot_type = if_else(measure == "precip_wk_total", "precip", "temp")) %>%
  mutate(plot_type = factor(plot_type, levels = c("temp", "precip")))

# Create plot labels
plot_labs <- c(
  "Minimum and maximum daily temperature (C)",
  "Total weekly precipitation (mm)"
)

plot_vars <- c("temp", "precip")

label_function <- function(var_list) {
  letters <- LETTERS[seq_along(var_list)]
  paste0(letters, ") ", var_list) %>%
    set_names(plot_vars)
}

labs <- label_function(plot_labs)

plot_data <- df_long

rf_plot_weekly <- ggplot(plot_data, aes(x = monday_date, y = value, group = measure, color = measure)) +
  # Add ENSO shading first 
  geom_rect(
    data = enso_shading,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = enso_type),
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  geom_line(na.rm = FALSE) +
  scale_x_date(
    date_breaks = "month", 
    date_labels = "%b %Y", 
    limits = c(as.Date("2017-01-04"), as.Date("2024-12-25")), 
    expand = c(0, 5)
  ) + 
  facet_wrap(~ plot_type, scales = "free_y", ncol = 1, labeller = labeller(plot_type = labs)) +
  scale_color_manual(values = c("temp_wk_min" = "#fca903", "temp_wk_max" = "#bf3b32", "precip_wk_total" = "#58a5e0")) +
  scale_fill_manual(
    values = c("El Niño" = "orange", "La Niña" = "cornflowerblue"),
    guide = "none"
  ) +
  theme_light() + 
  labs(x = NULL, y = NULL) +
  theme(
    plot.margin = margin(t = 0.1, r = 0.1, b = 0, l = 0.1, unit = "in"),
    strip.background = element_rect(fill = "white", color = "gray50"),
    strip.text.x = element_text(family = "serif", size = 10, color = "gray40", face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(family = "serif", size = 8),
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  expand_limits(x = c(as.Date("2017-01-04"), as.Date("2024-12-25")))

n_plot_dat <- nonlagged_data %>%
  group_by(monday_date) %>%
  summarize(
    n_cases = (sum(n_cases, na.rm = TRUE) / sum(population, na.rm = TRUE)) * 10000) %>%
  filter(monday_date >= as.Date("2017-01-01"), monday_date <= as.Date("2024-12-31")) %>%
  mutate(n_color = "black") %>%
  add_row(monday_date = as.Date("2017-01-04"), n_cases = 0, n_color = "black") %>%
  add_row(monday_date = as.Date("2024-12-25"), n_cases = 0, n_color = "black") %>%
  mutate(year = year(monday_date))

n_plot_weekly <- ggplot(n_plot_dat) + 
  # Add ENSO shading
  geom_rect(
    data = enso_shading,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = enso_type),
    alpha = 0.2,
    inherit.aes = FALSE) +
  geom_col(aes(x = monday_date, y = n_cases), fill = "black", na.rm = TRUE, alpha = 0.8, width = 6) + 
  scale_x_date(
    date_breaks = "year", 
    date_labels = "%Y", 
    limits = c(as.Date("2017-01-04"), as.Date("2024-12-25")), 
    expand = c(0, 5)
  ) + 
  labs(y = NULL, x = NULL, fill = "El Niño Southern Oscillation Cycle") + 
  scale_fill_manual(
    values = c("black" = "black", "El Niño" = "orange", "La Niña" = "cornflowerblue", "Neutral" = NA),
  ) +
  facet_wrap(~n_color, ncol = 1, labeller = as_labeller(c("black" = "C) Malaria Incidence per 10,000 person-weeks"))) +
  theme_light() + 
  theme(
    plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "in"),
    strip.background = element_rect(fill = "white", color = "gray50"),
    strip.text.x = element_text(family = "serif", size = 10, color = "gray40", face = "bold"),
    axis.text = element_text(family = "serif", size = 8),
    axis.text.x = element_text(family = "serif", angle = 90, size = 8),
    legend.position = "bottom",
    legend.text = element_text(family = "serif", size = 8),
    legend.title = element_text(family = "serif", size = 8),
    panel.grid.minor = element_blank()
  )

plot <- plot_grid(
  rf_plot_weekly, 
  n_plot_weekly, 
  ncol = 1, 
  rel_heights = c(2, 1.3), 
  axis = "lr", 
  align = "hv", 
  labels = NULL, 
  label_size = 1
)

ggsave(plot, filename = paste0(figure_path, "figure-3.pdf"),
       height = 8, width = 6, units = "in", dpi = 300)
