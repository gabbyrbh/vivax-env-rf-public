# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : 1-fit-dlnm-gam.R
# @Description  : Fit primary DLNM-GAMs for minimum temperature, maximum
#                 temperature, and total precipitation (Poisson/log link,
#                 natural splines 2df, 1 internal lag knot, 2-16 week lookback).
#                 Saves DLNM slice objects to results/dlnm/.
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
  mutate(year = factor(as.character(year), 
                       levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))) %>%
  filter(!year %in% c("2020", "2021")) %>%
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

# Define DLNM fitting function --------------------------------------------
fit_dlnm <- function(data, predictor, predictor_type, var_value, lag_value, covariate, predictor_name, color_code, plot_tags, outcome = "n_cases", offset_var = "population") {
  
  # Trim to 5th–95th percentile of predictor to avoid unstable estimates at extremes
  analysis_data <- data %>%
    mutate(!!sym(predictor) := as.numeric(!!sym(predictor)),
           !!sym(covariate) := as.numeric(!!sym(covariate))) %>%
    filter(!!sym(predictor) > quantile(data[[predictor]], p = 0.05) &
             !!sym(predictor) < quantile(data[[predictor]], p = 0.95))

  # Place knots at equal quantiles of the observed predictor distribution.
  # df = 2 → 1 internal knot in the exposure dimension (2 df natural spline).
  predictor_varknots <- equalknots(analysis_data[[predictor]], fun = "ns", df = 2)
  covariate_varknots <- equalknots(analysis_data[[covariate]], fun = "ns", df = 2)
  # 1 internal knot over the 2–16 week lag range (nk = 1).
  lagknots <- equalknots(x = c(2, 16), nk = 1)

  # Build cross-basis matrix for the primary predictor.
  # lag = c(2, 16): 2–16 week lookback; lags 0–1 excluded to reduce reverse-
  #   causality bias (observed cases can't cause prior-week weather).
  # bylag = 1: evaluate the lag-response at weekly resolution.
  # argvar: natural splines (fun = "ns") constrain the exposure-response to be
  #   linear beyond the boundary knots, reducing extrapolation artefacts.
  # arglag: same natural spline basis applied across the lag dimension.
  cb_predictor <- crossbasis(analysis_data[[predictor]], lag = c(2, 16), bylag = 1,
                             argvar = list(fun = "ns", knots = predictor_varknots),
                             arglag = list(knots = lagknots))

  # Mutual-adjustment cross-basis for the covariate (same specification).
  cb_covariate <- crossbasis(analysis_data[[covariate]], lag = c(2, 16), bylag = 1,
                             argvar = list(fun = "ns", knots = covariate_varknots),
                             arglag = list(knots = lagknots))

  # Poisson GAM with log link.
  # cb_predictor + cb_covariate: joint cross-bases allow simultaneous adjustment
  #   for lagged effects of both exposure and covariate.
  # year: fixed effect for secular trends (COVID years already excluded upstream).
  # comm_id: community fixed effects absorb time-invariant community confounding.
  # s(oni_index): smooth term for ENSO, captures non-linear climate cycle effects.
  # years_since_itn_delivery: categorical covariate for ITN distribution timing.
  # offset(log(population)): converts expected counts to incidence rates.
  formula <- as.formula(paste(outcome, "~ cb_predictor + cb_covariate + year + comm_id + s(oni_index) + years_since_itn_delivery + offset(log(", offset_var, "))"))

  model <- gam(data = analysis_data, formula = formula, family = "poisson")
  
  # save cov results
  coef_summary <- summary(model)$p.table
  coef_summary <- as.data.frame(coef_summary)
  coef_summary$term <- rownames(coef_summary)
  colnames(coef_summary) <- c("estimate", "std_error", "z_value", "p_value", "term")
  
  cov_results <- coef_summary %>%
    mutate(
      IR = exp(estimate),
      IR_Lower = exp(estimate - 1.96 * std_error),
      IR_Upper = exp(estimate + 1.96 * std_error),
      IR_Estimate = sprintf("%.2f", IR),
      IR_CI_Lower = sprintf("%.2f", IR_Lower),
      IR_CI_Upper = sprintf("%.2f", IR_Upper),
      IR_95CI = paste0(IR_Estimate, " (", IR_CI_Lower, " - ", IR_CI_Upper, ")")
    ) %>%
    select(Variable = term, IR_Estimate, IR_CI_Lower, IR_CI_Upper, IR_95CI) %>%
    filter(!str_detect(Variable, "^cb")) %>%
    filter(Variable != "(Intercept)")
  
  write.csv(cov_results, here("results", glue("covariate_ir_estimates_{predictor}.csv")), row.names = FALSE)
  
  # crosspred: generate predictions from the fitted cross-basis over the
  # observed exposure range.
  # cen: the reference value against which all IRs are computed.
  #   Temperature: minimum observed value (IRs express risk relative to coolest
  #   observed week; appropriate when relationship is monotone-increasing).
  #   Precipitation: 0 mm (IRs express risk relative to a dry week).
  # by: prediction grid resolution (0.1°C for temp; 50 mm for precip).
  # cumul = TRUE: also compute cumulative (overall) effects summed across all lags.
  if (predictor_type == "temp") {
    if (predictor == "temp_wk_max") {
      pred <- crosspred(cb_predictor, model,
                        cen = min(analysis_data[[predictor]], na.rm = TRUE),
                        by = 0.1, cumul = T)
    } else {
      pred <- crosspred(cb_predictor, model,
                        cen = min(analysis_data[[predictor]], na.rm = TRUE),
                        by = 0.1, cumul = T)
    }
  } else {
    pred <- crosspred(cb_predictor, model, cen = 0, by = 50, cumul = T)
  }
  
  plot_data <- data.frame(pred["predvar"], 
                          pred[["matRRfit"]]) %>% 
    pivot_longer(cols = starts_with("lag"),
                 names_to = "lag",
                 names_prefix = "lag",
                 values_to = "fit")
  
  lower_ci_data <- data.frame(pred["predvar"], 
                              pred["matRRlow"]) %>% 
    pivot_longer(cols = starts_with("matRRlow."),
                 names_to = "lag",
                 names_prefix = "matRRlow.lag",
                 values_to = "lower_ci")
  plot_data <- left_join(plot_data, lower_ci_data, by = c("predvar", "lag"))
  
  upper_ci_data <- data.frame(pred["predvar"], 
                              pred["matRRhigh"]) %>% 
    pivot_longer(cols = starts_with("matRRhigh."),
                 names_to = "lag",
                 names_prefix = "matRRhigh.lag",
                 values_to = "upper_ci")
  if (predictor == "temp_wk_max") {
    plot_data <- left_join(plot_data, upper_ci_data, by = c("predvar", "lag")) %>% 
      mutate(lag = as.numeric(lag)) %>% 
      mutate(predvar = round(predvar, 1))
  } else {
    plot_data <- left_join(plot_data, upper_ci_data, by = c("predvar", "lag")) %>% 
      mutate(lag = as.numeric(lag)) %>% 
      mutate(predvar = round(predvar, 1)) %>% 
      filter(predvar > quantile(data[[predictor]], p = 0.05) &
               predvar < quantile(data[[predictor]], p = 0.95)) # truncate limits to data range!
  }
  
  point_title = paste0("Association at Median")
  
  plot_list <- list()
  result_data <- data.frame()
  this_plot_data <- plot_data %>% filter(predvar == var_value)
  
  if (predictor_name == "Total Precipitation (mm)") {
      point_title = "90th percentile vs. 0"
  } else {
    point_title = "90th percentile vs. min"
  }
  
  plot <- ggplot(this_plot_data, aes(x = lag, y = fit)) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), 
                  width = 0.2, color = color_code) +
    geom_point(color = color_code) +
    geom_hline(yintercept = 1, linetype = "longdash", color = "black") +
    labs(title = point_title,
         tag = plot_tags[1],
         x = "Lag (weeks)", 
         y = "Incidence Ratio") +
    scale_y_continuous(trans = "log", 
                       breaks = scales::breaks_extended()) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(family = "serif", size = 10, face = "bold", hjust = 0.5, vjust = 3.0),
          plot.margin = margin(10, 10, 5, 10),
          plot.tag = element_text(family = "serif", size = 10, face = "bold", hjust = 0.2, vjust = -1.0),
          axis.line = element_line(color = "black", linewidth = 0.2),
          axis.title = element_text(family = "serif", size = 10), 
          axis.text = element_text(family = "serif", size = 10))
  result_data <- rbind(result_data, this_plot_data)
  plot_list[[1]] <- plot
  
  this_plot_data <- plot_data %>% filter(lag == lag_value)
  lag_plot <- ggplot(this_plot_data, mapping = aes(x = predvar, y = fit)) +
    geom_ribbon(mapping = aes(ymin = lower_ci, ymax = upper_ci), 
                alpha = 0.2, color = NA, fill = color_code) +
    geom_line(color = color_code) +
    geom_hline(yintercept = 1, linetype = "longdash", color = "black") +
    labs(title = paste0("Association at ", lag_value, "-week lag"),
         tag = plot_tags[2],
         x = predictor_name, 
         y = "Incidence Ratio") +
    scale_y_continuous(trans = "log", 
                       breaks = scales::breaks_extended()) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(family = "serif", size = 10, face = "bold", hjust = 0.5, vjust = 3.0),
          plot.margin = margin(10, 10, 5, 10),
          plot.tag = element_text(family = "serif", size = 10, face = "bold", hjust = 0.2, vjust = -1.0),
          axis.line = element_line(color = "black", linewidth = 0.2),
          axis.title = element_text(family = "serif", size = 10),
          axis.text = element_text(family = "serif", size = 10))
  result_data <- rbind(result_data, this_plot_data)
  plot_list[[2]] <- lag_plot
  
  # create contour plot
  contour_plot <- ggplot(plot_data, aes(x = predvar, y = lag, z = fit)) +
    geom_raster(aes(fill = fit)) +
    geom_contour(color = "white", alpha = 0.5) +
    scale_fill_gradient2(
      low = "#4DAF4A",        
      mid = "white",       
      high = "#C51B7D",        
      midpoint = 1,
      limits = c(0.9, 1.2),      # Fixed limits for all plots
      breaks = c(0.9, 1.0, 1.1, 1.2),  # Fixed breaks for all plots
      labels = c("0.9", "1.0", "1.1", "1.2"),  # Fixed labels
      oob = scales::squish  # Handle values outside limits
    ) +
    labs(title = "Predictor-Lag Surface",
         tag = plot_tags[3],
         x = predictor_name,
         y = "Lag (weeks)",
         fill = "Incidence\nRatio"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(family = "serif", size = 10, face = "bold", hjust = 0.5, vjust = 3.0),
          plot.margin = margin(10, 10, 5, 10),
          plot.tag = element_text(family = "serif", size = 10, face = "bold", hjust = 0.2, vjust = -1.0),
          axis.title = element_text(family = "serif", size = 10),
          axis.text = element_text(family = "serif", size = 10),
          legend.text = element_text(family = "serif", size = 8),
          legend.title = element_text(family = "serif", size = 8),
          legend.key.size = unit(0.6, "cm"),
          legend.key.width = unit(0.3, "cm"),
          panel.grid = element_blank())
  plot_list[[3]] <- contour_plot
  
  result_data <- result_data %>% mutate(Measure = predictor_name,
                                        Covariate = covariate,
                                        `Predictor Value` = predvar,
                                        `Lag` = lag,
                                        `Incidence Ratio (95% CI)` = paste0(sprintf("%.02f", fit), " (", sprintf("%.02f", lower_ci), " - ", sprintf("%.02f", upper_ci), ")")) %>% select(Measure, Covariate, Lag, `Predictor Value`, `Incidence Ratio (95% CI)`)
  
  plot_list[[4]] <- result_data
  # print(summary(model))
  return(plot_list)
}

# Fit models and save outputs ---------------------------------------------
plots_temp_min <- fit_dlnm(data = nonlagged_data,
                           predictor = "temp_wk_min",
                           predictor_type = "temp",
                           var_value = 23.5, lag_value = 7,
                           covariate = "precip_wk_total",
                           predictor_name = "Minimum Temperature (C)",
                           color_code = "#FF9F4D",
                           plot_tags = c("A)", "B)", "C)"))
saveRDS(plots_temp_min, paste0(results_path, "dlnm/temp_wk_min-slices.RDS"))

plots_temp_max <- fit_dlnm(data = nonlagged_data, 
                           predictor = "temp_wk_max",
                           predictor_type = "temp",
                           var_value = 33.7, 
                           lag_value = 9,
                           covariate = "precip_wk_total",
                           predictor_name = "Maximum Temperature (C)",
                           color_code = "#E07A5F",
                           plot_tags = c("D)", "E)", "F)"))
saveRDS(plots_temp_max, paste0(results_path, "dlnm/temp_wk_max-slices.RDS"))

plots_precip_total <- fit_dlnm(data = nonlagged_data, predictor = "precip_wk_total",
                               predictor_type = "precip",
                               var_value = 1000, lag_value = 11,
                               covariate = "temp_wk_min",
                               predictor_name = "Total Precipitation (mm)",
                               color_code = "#5BA3C7",
                               plot_tags = c("G)", "H)", "I)"))

saveRDS(plots_precip_total, paste0(results_path, "dlnm/precip_wk_total-slices.RDS"))

# save point estimates with CIs
ir_results <- rbind(plots_temp_min[[4]], plots_temp_max[[4]], plots_precip_total[[4]])
write.csv(ir_results, here("results", "dlnm-ir-estimates.csv"))

