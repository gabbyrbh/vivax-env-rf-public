# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : fig5-plot-sub-groups.R
# @Description  : Generate Figure 5: DLNM subgroup results for ENSO period
#                 (El Nino / La Nina / Neutral) and community type
#                 (highway / riverine) with bootstrap 95% CIs.
# ------------------------------------------------------------------------------
rm(list = ls())
source(paste0(here::here(), "/0-base-functions.R"))
source(paste0(here::here(), "/0-config.R"))
library(mgcv)
library(lme4)
library(lmtest)
library(sandwich)
library(gKRLS)
library(dlnm)
library(splines)
library(kableExtra)
library(glue)
library(data.table)
library(broom)
library(boot)
library(scales) 
library(ggplot2)  
library(dplyr)    
library(cowplot)  
select <- dplyr::select
summarize <- dplyr::summarize

nonlagged_data <- read.csv(paste0(public_data_path, "vivax-env-erf-public.csv")) %>%
  rename(time = "week") %>%
  mutate(population = as.numeric(population),
         n_cases = as.numeric(n_cases),
         comm_id = as.factor(comm_id),
         oni_index = as.numeric(oni_index)) %>%
  filter(!is.na(population)) %>%
  mutate(year = factor(as.character(year), levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))) %>%
  filter(year != "2020") %>%
  filter(year != "2021") %>%
  filter(population > 15)

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

# Load Community Type Data ----------------------------------------------------
base_precip_total_comm <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-precip-wk-total-highway.RDS")) %>% 
    mutate(comm_type = "highway"),
  readRDS(file = paste0(here::here(), "/results/base results/base-precip-wk-total-riverine.RDS")) %>% 
    mutate(comm_type = "riverine")
)

bs_precip_total_comm <- rbind(
  readRDS(paste0(bs_path, "bs-precip-wk-total-highway.RDS")) %>% mutate(comm_type = "highway"),
  readRDS(paste0(bs_path, "bs-precip-wk-total-riverine.RDS")) %>% mutate(comm_type = "riverine")
)

base_temp_max_comm <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-wk-max-highway.RDS")) %>% 
    mutate(comm_type = "highway"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-wk-max-riverine.RDS")) %>% 
    mutate(comm_type = "riverine"))

bs_temp_max_comm <- rbind(
  readRDS(paste0(bs_path, "bs-temp-wk-max-highway.RDS")) %>% mutate(comm_type = "highway"),
  readRDS(paste0(bs_path, "bs-temp-wk-max-riverine.RDS")) %>% mutate(comm_type = "riverine")
)

base_temp_min_comm <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-wk-min-highway.RDS")) %>% 
    mutate(comm_type = "highway"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-wk-min-riverine.RDS")) %>% 
    mutate(comm_type = "riverine")
)

bs_temp_min_comm <- rbind(
  readRDS(paste0(bs_path, "bs-temp-wk-min-highway.RDS")) %>% mutate(comm_type = "highway"),
  readRDS(paste0(bs_path, "bs-temp-wk-min-riverine.RDS")) %>% mutate(comm_type = "riverine")
)

# Load ENSO Data --------------------------------------------------------------
temp_min_base_enso <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-min-enso-normal.RDS")) %>% mutate(enso_period = "Neutral"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-min-enso-el-nino.RDS")) %>% mutate(enso_period = "El Niño"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-min-enso-la-nina.RDS")) %>% mutate(enso_period = "La Niña")
) %>% mutate(incidence_10k = incidence*10000)

temp_max_base_enso <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-max-enso-normal.RDS")) %>% mutate(enso_period = "Neutral"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-max-enso-el-nino.RDS")) %>% mutate(enso_period = "El Niño"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-max-enso-la-nina.RDS")) %>% mutate(enso_period = "La Niña")
) %>% mutate(incidence_10k = incidence*10000)

precip_total_base_enso <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-precip-total-enso-normal.RDS")) %>% mutate(enso_period = "Neutral"),
  readRDS(file = paste0(here::here(), "/results/base results/base-precip-total-enso-el-nino.RDS")) %>% mutate(enso_period = "El Niño"),
  readRDS(file = paste0(here::here(), "/results/base results/base-precip-total-enso-la-nina.RDS")) %>% mutate(enso_period = "La Niña")
) %>% mutate(incidence_10k = incidence*10000)

temp_min_bs_enso <- rbind(
  readRDS(file = paste0(bs_path, "bs-temp-min-enso-normal.RDS")) %>% 
    mutate(enso_period = "Neutral"),
  readRDS(file = paste0(bs_path, "bs-temp-min-enso-el-nino.RDS")) %>% 
    mutate(enso_period = "El Niño"),
  readRDS(file = paste0(bs_path, "bs-temp-min-enso-la-nina.RDS")) %>% 
    mutate(enso_period = "La Niña")
)

temp_max_bs_enso <- rbind(
  readRDS(file = paste0(bs_path, "bs-temp-max-enso-normal.RDS")) %>% mutate(enso_period = "Neutral"),
  readRDS(file = paste0(bs_path, "bs-temp-max-enso-el-nino.RDS")) %>% mutate(enso_period = "El Niño"),
  readRDS(file = paste0(bs_path, "bs-temp-max-enso-la-nina.RDS")) %>% mutate(enso_period = "La Niña")
)

precip_total_bs_enso <- rbind(
  readRDS(file = paste0(bs_path, "bs-precip-total-enso-normal.RDS")) %>% mutate(enso_period = "Neutral"),
  readRDS(file = paste0(bs_path, "bs-precip-total-enso-el-nino.RDS")) %>% mutate(enso_period = "El Niño"),
  readRDS(file = paste0(bs_path, "bs-precip-total-enso-la-nina.RDS")) %>% mutate(enso_period = "La Niña")
)

# Function for Community Type plots
plot_gam_comm_type <- function(base_result, 
                               boot_results,
                               predictor, predictor_name, lag_value,
                               axis_breaks, axis_limits, plot_tag, show_legend = FALSE) {
  
  ci_data <- boot_results %>%
    group_by(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")), comm_type) %>%
    summarize(
      lower_ci = quantile(incidence, 0.025),
      upper_ci = quantile(incidence, 0.975),
      .groups = "drop"
    )
  
  final_result <- base_result %>%
    left_join(ci_data, by = c(glue::glue("{predictor}_{lag_value}_weeklag"), "comm_type"))
  
  title <- paste0(predictor_name, ":\nAssociation at ", lag_value, "-week lag")
  
  p <- ggplot(final_result, aes(x = !!sym(glue("{predictor}_{lag_value}_weeklag")), 
                                y = incidence*10000)) +
    geom_line(aes(color = comm_type), linewidth = 0.5) +
    geom_ribbon(aes(ymin = lower_ci*10000, ymax = upper_ci*10000, fill = comm_type),
                alpha = 0.2) +
    scale_color_manual(values = c("highway" = "peru", 
                                  "riverine" = "cyan4"),
                       labels = c("highway" = "Highway",
                                  "riverine" = "Riverine"),
                       name = "Community Type") +
    scale_fill_manual(values = c("highway" = "peru", 
                                 "riverine" = "cyan4"),
                      labels = c("highway" = "Highway",
                                 "riverine" = "Riverine"),
                      name = "Community Type") +
    scale_y_continuous(limits = c(0, 40)) +
    labs(title = title,
         x = predictor_name,
         y = "Marginal Incidence",
         tag = plot_tag) +
    theme_minimal() +
    theme(plot.title = element_text(family = "serif", size = 12, hjust = 0.5, vjust = 3.0, face = "bold"),
          plot.margin = margin(5, 15, 5, 0),
          plot.tag = element_text(family = "serif", size = 12, face = "bold", hjust = 0.2, vjust = -0.5),
          axis.line = element_line(color = "black", linewidth = 0.2),
          axis.title.y = element_text(family = "serif", size = 12, margin = margin(r = 5)),
          axis.title.x = element_text(family = "serif", size = 12, margin = margin(b = 5)),
          axis.text = element_text(family = "serif", size = 12),
          legend.position = ifelse(show_legend, "bottom", "none"),
          legend.margin = margin(1, 1, 5, 5),
          legend.title = element_text(family = "serif", size = 12),
          legend.text = element_text(family = "serif", size = 12),
          legend.key.height = unit(0.3, "in"),
          legend.key.width = unit(0.3, "in"),
          legend.key.spacing = unit(0.05, "in"),
          panel.grid.major = element_line(linewidth = 1),
          panel.grid.minor = element_line(linewidth = 0.5))
  
  return(p)
}

# Function for ENSO plots
plot_gam_enso <- function(base_result, 
                          boot_results,
                          predictor, predictor_name, lag_value,
                          axis_breaks, axis_limits, plot_tag, show_legend = FALSE) {
  
  # Check if this is a precipitation dataset
  is_precip <- predictor == "precip_wk_total"
  
  if (is_precip) {
    predictor_col <- glue::glue("{predictor}_{lag_value}_weeklag")
    
    base_unique_vals <- unique(base_result[[predictor_col]])
    bs_unique_vals <- unique(boot_results[[predictor_col]])
    
    base_index <- data.frame(
      idx = 1:length(base_unique_vals),
      pred_val = base_unique_vals
    )
    
    bs_index <- data.frame(
      idx = 1:length(bs_unique_vals),
      pred_val = bs_unique_vals
    )
    
    ci_data <- boot_results %>%
      rename(pred_val = !!predictor_col) %>%
      left_join(bs_index, by = "pred_val") %>%
      group_by(idx, enso_period) %>%
      summarize(
        lower_ci = quantile(incidence, 0.025),
        upper_ci = quantile(incidence, 0.975),
        .groups = "drop"
      ) 
    
    final_result <- base_result %>%
      rename(pred_val = !!predictor_col) %>%
      left_join(base_index, by = "pred_val") %>%
      left_join(ci_data, by = c("idx", "enso_period")) %>%
      rename(!!predictor_col := pred_val)
    
  } else {
    # For temperature: apply sprintf formatting
    ci_data <- boot_results %>%
      mutate(pred_val = sprintf("%0.1f", !!sym(glue::glue("{predictor}_{lag_value}_weeklag")))) %>%
      group_by(pred_val, enso_period) %>%
      summarize(
        lower_ci = quantile(incidence, 0.025),
        upper_ci = quantile(incidence, 0.975),
        .groups = "drop"
      ) 
    
    final_result <- base_result %>%
      mutate(pred_val = sprintf("%0.1f", !!sym(glue::glue("{predictor}_{lag_value}_weeklag")))) %>%
      left_join(ci_data, by = c("pred_val", "enso_period"))
  }
  
  title <- paste0(predictor_name, ":\nAssociation at ", lag_value, "-week lag")
  
  p <- ggplot(final_result, aes(x = !!sym(glue("{predictor}_{lag_value}_weeklag")), 
                                y = incidence*10000)) +
    geom_line(aes(color = enso_period), linewidth = 0.5) +
    geom_ribbon(aes(ymin = lower_ci*10000, ymax = upper_ci*10000, fill = enso_period),
                alpha = 0.2) + 
    scale_color_manual(values = c("Neutral" = "darkgreen", 
                                  "El Niño" = "orange",
                                  "La Niña" = "cornflowerblue"),
                       name = "ENSO Period") +
    scale_fill_manual(values = c("Neutral" = "darkgreen", 
                                 "El Niño" = "orange",
                                 "La Niña" = "cornflowerblue"),
                      name = "ENSO Period") +
    #scale_x_continuous(limits = axis_limits, breaks = axis_breaks) +
    scale_y_continuous(limits = c(0, 40)) +
    labs(title = title,
         x = predictor_name,
         y = "Marginal Incidence",
         tag = plot_tag) +
    theme_minimal() +
    theme(plot.title = element_text(family = "serif", size = 12, hjust = 0.5, vjust = 3.0, face = "bold"),
          plot.margin = margin(5, 10, 5, 5),
          plot.tag = element_text(family = "serif", size = 12, face = "bold", hjust = 0.2, vjust = -0.5),
          axis.line = element_line(color = "black", linewidth = 0.2),
          axis.title.y = element_text(family = "serif", size = 12, margin = margin(r = 5)),
          axis.title.x = element_text(family = "serif", size = 12, margin = margin(b = 5)),
          axis.text = element_text(family = "serif", size = 12),
          legend.position = ifelse(show_legend, "bottom", "none"),
          legend.margin = margin(1, 1, 5, 5),
          legend.title = element_text(family = "serif", size = 12),
          legend.text = element_text(family = "serif", size = 12),
          legend.key.height = unit(0.3, "in"),
          legend.key.width = unit(0.3, "in"),
          legend.key.spacing = unit(0.05, "in"),
          panel.grid.major = element_line(linewidth = 1),
          panel.grid.minor = element_line(linewidth = 0.5)
    )
  
  return(p)
}

# Generate Community Type plots (top row)
plot_temp_min_comm <- plot_gam_comm_type(base_result = base_temp_min_comm,
                                         boot_results = bs_temp_min_comm,
                                         predictor = "temp_wk_min",
                                         predictor_name = "Minimum Temperature (C)",
                                         lag_value = 7,
                                         axis_limits = c(20, 25),
                                         axis_breaks = c(20, 21, 22, 23, 24),
                                         plot_tag = "A)",
                                         show_legend = FALSE)

plot_temp_max_comm <- plot_gam_comm_type(base_result = base_temp_max_comm,
                                         boot_results = bs_temp_max_comm,
                                         predictor = "temp_wk_max",
                                         predictor_name = "Maximum Temperature (C)",
                                         lag_value = 9,
                                         axis_limits = c(30, 34),
                                         axis_breaks = c(30, 31, 32, 33, 34),
                                         plot_tag = "B)",
                                         show_legend = TRUE)

plot_precip_total_comm <- plot_gam_comm_type(base_result = base_precip_total_comm,
                                             boot_results = bs_precip_total_comm,
                                             predictor = "precip_wk_total",
                                             lag_value = 11,
                                             predictor_name = "Total Precipitation (mm)",
                                             axis_limits = c(0, 1200),
                                             axis_breaks = c(0, 200, 400, 600, 800, 1000, 1200),
                                             plot_tag = "C)",
                                             show_legend = FALSE)

# Generate ENSO plots (bottom row)
plot_temp_min_enso <- plot_gam_enso(base_result = temp_min_base_enso,
                                    boot_results = temp_min_bs_enso,
                                    predictor = "temp_wk_min",
                                    predictor_name = "Minimum Temperature (C)",
                                    lag_value = 7,
                                    axis_limits = c(20, 25),
                                    axis_breaks = c(20, 21, 22, 23, 24),
                                    plot_tag = "D)",
                                    show_legend = FALSE)

plot_temp_max_enso <- plot_gam_enso(base_result = temp_max_base_enso,
                                    boot_results = temp_max_bs_enso,
                                    predictor = "temp_wk_max",
                                    predictor_name = "Maximum Temperature (C)",
                                    lag_value = 9,
                                    axis_limits = c(30, 34),
                                    axis_breaks = c(30, 31, 32, 33, 34),
                                    plot_tag = "E)",
                                    show_legend = TRUE)

plot_precip_total_enso <- plot_gam_enso(base_result = precip_total_base_enso,
                                        boot_results = precip_total_bs_enso,
                                        predictor = "precip_wk_total",
                                        predictor_name = "Total Precipitation (mm)",
                                        lag_value = 11,
                                        axis_limits = c(0, 1200),
                                        axis_breaks = c(0, 200, 400, 600, 800, 1000, 1200),
                                        plot_tag = "F)",
                                        show_legend = FALSE)

# Extract legends separately
legend_comm <- cowplot::get_legend(
  plot_gam_comm_type(base_result = base_temp_min_comm,
                     boot_results = bs_temp_min_comm,
                     predictor = "temp_wk_min",
                     predictor_name = "Minimum Temperature (C)",
                     lag_value = 7,
                     axis_limits = c(20, 24),
                     axis_breaks = c(20, 21, 22, 23, 24),
                     plot_tag = "",
                     show_legend = TRUE)
)

legend_enso <- cowplot::get_legend(
  plot_gam_enso(base_result = temp_min_base_enso,
                boot_results = temp_min_bs_enso,
                predictor = "temp_wk_min",
                predictor_name = "Minimum Temperature (C)",
                lag_value = 7,
                axis_limits = c(20, 24),
                axis_breaks = c(20, 21, 22, 23, 24),
                plot_tag = "",
                show_legend = TRUE)
)

# Combine plots in 1x3 grid
comm_plots <- cowplot::plot_grid(
  plot_temp_min_comm, plot_temp_max_comm, plot_precip_total_comm,
  nrow = 1, ncol = 3, align = "hv"
)
ggsave(comm_plots, file = paste0(figure_path, "figure-5-comm-type.pdf"),
       width = 12, height = 5, units = "in", dpi = 300)

# Combine plots in 1x3 grid
enso_plots <- cowplot::plot_grid(
  plot_temp_min_enso, plot_temp_max_enso, plot_precip_total_enso,
  nrow = 1, ncol = 3, align = "hv"
)
ggsave(enso_plots, file = paste0(figure_path, "figure-5-enso.pdf"),
       width = 12, height = 5, units = "in", dpi = 300)

# Combine plots in 2x3 grid
plots_grid <- cowplot::plot_grid(
  plot_temp_min_comm, plot_temp_max_comm, plot_precip_total_comm,
  plot_temp_min_enso, plot_temp_max_enso, plot_precip_total_enso,
  nrow = 2, ncol = 3, align = "hv"
)

ggsave(plots_grid, file = paste0(figure_path, "figure-5.pdf"),
       width = 12, height = 8, dpi = 300)
