# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : fig4-plot-main-gam.R
# @Description  : Generate Figure 4: DLNM exposure-response and lag-response
#                 slices for the primary analyses (min temp 7-wk, max temp 9-wk,
#                 total precipitation 11-wk) with bootstrap 95% CIs.
# ------------------------------------------------------------------------------
rm(list = ls())
source(paste0(here::here(), "/0-base-functions.R"))
source(paste0(here::here(), "/0-config.R"))
library(mgcv)
library(lme4)
library(glue)
library(broom)
library(gridExtra)  # For grid.arrange()
library(tidyr)      # For pivot_longer()
library(stringr)
select <- dplyr::select
summarize <- dplyr::summarize

# Load data ---------------------------------------------------------------
nonlagged_data <- read.csv(paste0(public_data_path, "vivax-env-erf-public.csv")) %>%
  rename(time = "week") %>%
  mutate(population = as.numeric(population),
         n_cases = as.numeric(n_cases),
         comm_id = as.factor(comm_id),
         oni_index = as.numeric(oni_index)) %>%
  filter(!is.na(population)) %>%
  mutate(year = factor(as.character(year), levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))) %>%
  filter(year != "2020") %>%
  filter(year != "2021") 

# set median pop community as reference
median_pop_comm <- nonlagged_data %>%
  filter(year == 2019) %>%
  group_by(comm_id) %>%
  summarize(avg_population = mean(population)) %>%
  ungroup() %>%
  mutate(diff_from_median = abs(avg_population - median(avg_population))) %>%
  slice_min(diff_from_median) %>%
  pull(comm_id) %>% unique() %>% first() %>% as.character()

nonlagged_data <- nonlagged_data %>% mutate(comm_id = relevel(comm_id, ref = median_pop_comm))

# comm_type is included in the public dataset
nonlagged_data <- nonlagged_data %>%
  mutate(comm_type = ifelse(is.na(comm_type) | comm_type == "", "riverine", comm_type))

# Calculate incidence ----------------------------------------------------------
incidence = sum(nonlagged_data$n_cases)/sum(nonlagged_data$population)*10000

# Read in DLNM results ---------------------------------------------------------
plots_temp_min <- readRDS(paste0(results_path, "dlnm/temp_wk_min-slices.RDS"))

plots_temp_max <- readRDS(paste0(results_path, "dlnm/temp_wk_max-slices.RDS"))

plots_precip_total <- readRDS(paste0(results_path, "dlnm/precip_wk_total-slices.RDS"))

ir_results <- read.csv(here("results", "dlnm-ir-estimates.csv"))

# Read in GAM results ----------------------------------------------------------
gam_temp_min <- readRDS(paste0(here::here(), "/results/base results/base-temp-min.RDS")) %>% 
  mutate(incidence_10k = incidence * 10000) %>%
  mutate(mean_incidence_10k = mean(incidence) * 10000,
         incidence_difference_10k = ((max(incidence) - min(incidence))*10000),
         incidence_ratio = max(incidence)/min(incidence))

gam_temp_max <- readRDS(paste0(here::here(), "/results/base results/base-temp-max.RDS")) %>% 
  mutate(incidence_10k = incidence * 10000)

gam_precip_total <- readRDS(paste0(here::here(), "/results/base results/base-precip-total.RDS")) %>% 
  mutate(incidence_10k = incidence * 10000)

# Load in bootstrap results from Sherlock --------------------------------------
# Load bootstrap results from Sherlock
bs_temp_wk_min <- readRDS(paste0(bs_path, "bs-temp-min.RDS"))

bs_temp_wk_max <- readRDS(paste0(bs_path, "bs-temp-max.RDS"))

bs_precip_wk_total <- readRDS(paste0(bs_path, "bs-precip-total.RDS"))

# Define plotting function ------------------------------------------------
plot_gam <- function(base_result, boot_results, predictor, predictor_name, lag_value, color_code, plot_tags, axis_breaks, axis_limits) {
  # Compute confidence intervals
  ci_data <- boot_results %>%
    group_by(!!sym(glue("{predictor}_{lag_value}_weeklag"))) %>%
    summarize(
      lower_ci = quantile(incidence, 0.025),
      upper_ci = quantile(incidence, 0.975),
      .groups = "drop"
    )
  
  # Combine base result with confidence intervals
  final_result <- base_result %>%
    left_join(ci_data, by = glue("{predictor}_{lag_value}_weeklag"))
  
  hist_data <- nonlagged_data %>%
    filter(!!sym(predictor) >= axis_limits[1] & 
             !!sym(predictor) <= axis_limits[2]) %>%
    filter(!is.na(!!sym(predictor)))
  
  title <- paste0("Association at ", lag_value, " Weeks")
  
  p <- ggplot(final_result, aes(x = !!sym(glue::glue("{predictor}_{lag_value}_weeklag")),
                                  y = incidence*10000)) +
    geom_line(color = color_code) +
    geom_ribbon(aes(ymin = lower_ci*10000, ymax = upper_ci*10000),
                alpha = 0.25, fill = color_code) +
    geom_histogram(data = hist_data, 
                   binwidth = ifelse(predictor == "precip_wk_total", 100, 0.2),
                   aes(x = !!sym(predictor), y = after_stat(count)/max(after_stat(count)) * 5),
                   color = NA, fill = "black", alpha = 0.3, inherit.aes = FALSE) +
    labs(title = title,
         x = predictor_name,
         y = "Marginal Incidence",
         tags = plot_tags) +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(family = "serif", size = 10, face = "bold", hjust = 0.5, vjust = 3.0),
          plot.margin = margin(10, 10, 5, 10),
          plot.tag = element_text(family = "serif", size = 10, face = "bold", hjust = 0.5, vjust = -0.5),
          axis.line = element_line(color = "black", linewidth = 0.2),
          axis.title = element_text(family = "serif", size = 10),
          axis.text = element_text(family = "serif", size = 10))
  
  return(p)
}

# Generate and save figures -----------------------------------------------
gam_plot_temp_min <- plot_gam(base_result = gam_temp_min,
                              boot_results = bs_temp_wk_min,
                              predictor = "temp_wk_min",
                              lag_value = 7,
                              predictor_name = "Minimum Temperature (C)",
                              color_code = "#FF9F4D",
                              plot_tags = "B)",
                              axis_limits = c(20, 24),
                              axis_breaks = c(20, 21, 22, 23, 24))

gam_plot_temp_max <- plot_gam(base_result = gam_temp_max, 
                              boot_results = bs_temp_wk_max,
                              predictor = "temp_wk_max",
                              lag_value = 9,
                              predictor_name = "Maximum Temperature (C)",
                              color_code = "#E07A5F",
                              plot_tags = "E)",
                              axis_limits = c(29.5, 34.5),
                              axis_breaks = c(30, 31, 32, 33, 34))

gam_plot_precip_total <- plot_gam(base_result = gam_precip_total,
                                  boot_results = bs_precip_wk_total,
                                  predictor = "precip_wk_total",
                                  lag_value = 11,
                                  predictor_name = "Total Precipitation (mm)",
                                  color_code = "#5BA3C7",
                                  plot_tags = "H)",
                                  axis_limits = c(0, 1200),
                                  axis_breaks = c(0, 200, 400, 600, 800, 1000, 1200))

# Create combined dlnm-gam plots -----------------------------------------------
temp_min <- arrangeGrob(
  plots_temp_min[[1]], gam_plot_temp_min, plots_temp_min[[3]],
  widths = c(1, 1, 1.2), ncol = 3,
  top = textGrob("Minimum Temperature", 
                 gp = gpar(fontsize = 12, fontface = "italic", fontfamily = "serif"),
                 just = "left", x = 0.02, y = 0))
ggsave(paste0(figure_path, "figure-4-temp-min.pdf"),
       temp_min, bg = "white",
       width = 20, height = 6, units = "cm", dpi = 300)

temp_max <- arrangeGrob(
  plots_temp_max[[1]], gam_plot_temp_max, plots_temp_max[[3]],
  widths = c(1, 1, 1.2), ncol = 3,
  top = textGrob("Maximum Temperature", 
                 gp = gpar(fontsize = 12, fontface = "italic", fontfamily = "serif"),
                 just = "left", x = 0.02, y = 0))
ggsave(paste0(figure_path, "figure-4-temp-max.pdf"),
       temp_max, bg = "white",
       width = 20, height = 6, units = "cm", dpi = 300)

precip_total <- arrangeGrob(
  plots_precip_total[[1]], gam_plot_precip_total, plots_precip_total[[3]],
  widths = c(1, 1, 1.2), ncol = 3,
  top = textGrob("Total Precipitation", 
                 gp = gpar(fontsize = 12, fontface = "italic", fontfamily = "serif"),
                 just = "left", x = 0.02, y = 0))
ggsave(paste0(figure_path, "figure-4-precip-total.pdf"),
       precip_total, bg = "white",
       width = 20, height = 6, units = "cm", dpi = 300)

# Combine rows
figure_4 <- arrangeGrob(temp_min, temp_max, precip_total, nrow = 3)

ggsave(paste0(figure_path, "figure-4.pdf"),
       figure_4, bg = "white",
       width = 20, height = 15, units = "cm", dpi = 300)

