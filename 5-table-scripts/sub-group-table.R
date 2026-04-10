# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : sub-group-table.R
# @Description  : Generate Table 1: incidence ratios by community type
#                 (highway / riverine) from GLM with interaction term.
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
select <- dplyr::select
summarize <- dplyr::summarize

# Load data ---------------------------------------------------------------
nonlagged_data <- readRDS(paste0(box_path_flame_erf,
                                 "non-lagged-analysis-data_ext.RDS")) %>%
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

pop_ref <- nonlagged_data %>% 
  select(year, comm_id, population) %>%
  distinct()

# Set median pop community as reference
median_pop_comm <- nonlagged_data %>%
  filter(year == 2019) %>%
  group_by(comm_id) %>%
  summarize(avg_population = mean(population)) %>%
  ungroup() %>%
  mutate(diff_from_median = abs(avg_population - median(avg_population))) %>%
  slice_min(diff_from_median) %>%
  pull(comm_id) %>% unique() %>% first() %>% as.character()

nonlagged_data <- nonlagged_data %>% mutate(comm_id = relevel(comm_id, ref = median_pop_comm))

# Read in comm_type csv
comm_type <- read.csv(paste0(box_path_flame_erf, "all_district_centroids_comm_type.csv")) %>%
  mutate(comm_id = sprintf("%03d", comm_id)) %>%
  mutate(comm_type = ifelse(comm_type == "", "riverine", comm_type)) %>%
  filter(!is.na(comm_type)) %>%
  mutate(comm_id = as.factor(comm_id)) %>%
  mutate(comm_id = relevel(comm_id, ref = median_pop_comm)) %>%
  mutate(highway = ifelse(comm_type == "highway", 1, 0),
         riverine = ifelse(comm_type == "riverine", 1, 0))

nonlagged_data <- nonlagged_data %>% 
  left_join(comm_type, by = "comm_id") %>%
  filter(!is.na(comm_type)) %>%
  mutate(comm_type = factor(comm_type, levels = c("highway", "riverine"), labels = c("Highway", "Riverine")))

# Add in ENSO
nonlagged_data <- nonlagged_data %>%
  mutate(
    enso_period = case_when(
      oni_index >= 0.5 ~ "El Niño periods",
      oni_index <= -0.5 ~ "La Niña periods",
      TRUE ~ "Normal ENSO periods")) %>%
  mutate(enso_period = factor(enso_period, 
                              levels = c("Normal ENSO periods", "El Niño periods", "La Niña periods"))) %>% 
  mutate(el_nino = case_when(
    enso_period == "El Niño periods" ~ "El Niño periods",
    enso_period == "Normal ENSO periods" ~ "Normal ENSO periods",
    TRUE ~ NA_character_),
    la_nina = case_when(
      enso_period == "La Niña periods" ~ "La Niña periods",
      enso_period == "Normal ENSO periods" ~ "Normal ENSO periods",
      TRUE ~ NA_character_
    )) %>%
  mutate(
    el_nino = factor(el_nino, levels = c("Normal ENSO periods", "El Niño periods")),
    la_nina = factor(la_nina, levels = c("Normal ENSO periods", "La Niña periods"))
  )

enso_count <- nonlagged_data %>%
  select(year, study_week, la_nina, el_nino) %>%
  distinct()

# Calculate incidence
comm_type_inc <- nonlagged_data %>%
  group_by(comm_type) %>%
  summarize(incidence_10k = sum(n_cases)/sum(population)*10000)

enso_inc <- nonlagged_data %>%
  group_by(enso_period) %>%
  summarize(incidence_10k = sum(n_cases)/sum(population)*10000)

# Function to lag data
lag_data <- function(data, predictor, covariate, lag_value) {
  lagged_data <- data %>%
    group_by(comm_id) %>%
    mutate(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) := lag(!!sym(predictor), n = lag_value),
           !!sym(glue::glue("{covariate}_{lag_value}_weeklag")) := lag(!!sym(covariate), n = lag_value)) %>%
    ungroup() %>%
    mutate(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) := 
             as.numeric(!!sym(glue::glue("{predictor}_{lag_value}_weeklag"))),
           !!sym(glue::glue("{covariate}_{lag_value}_weeklag")) := 
             as.numeric(!!sym(glue::glue("{covariate}_{lag_value}_weeklag")))) %>%
    rename("longitude" = "long")
  
  return(lagged_data)
}

fit_intx_glm <- function(data, predictor, lag_value, covariate, intx_var, intx_var_name) {
  
  lagged_data <- lag_data(data, predictor, covariate, lag_value)
  
  analysis_data <- lagged_data %>%
    filter(!is.na(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")))) %>%
    filter(!is.na(!!sym(intx_var))) %>%
    mutate(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) := 
             as.numeric(!!sym(glue::glue("{predictor}_{lag_value}_weeklag"))),
           !!sym(glue::glue("{covariate}_{lag_value}_weeklag")) := 
             as.numeric(!!sym(glue::glue("{covariate}_{lag_value}_weeklag"))))
  
  analysis_data <- analysis_data %>%
    filter(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) > 
             quantile(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]], p = 0.05) &
             !!sym(glue::glue("{predictor}_{lag_value}_weeklag")) < 
             quantile(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]], p = 0.95))
  
  model_intx <- glm(data = analysis_data,
                    as.formula(glue::glue(
                      "n_cases ~
                   {predictor}_{lag_value}_weeklag +
                   {covariate}_{lag_value}_weeklag +
                   {predictor}_{lag_value}_weeklag*{intx_var} +
                   {intx_var} +
                   oni_index +
                   comm_id +
                   years_since_itn_delivery +
                   year")),
                    offset = log(population),
                    family = "poisson")
  
  intx_df <- tidy(model_intx) %>%
    filter(term %in% c(glue::glue("{predictor}_{lag_value}_weeklag"),
                       glue::glue("{predictor}_{lag_value}_weeklag:{intx_var}{intx_var_name}")))
  
  beta_main <- intx_df %>%
    filter(term == glue::glue("{predictor}_{lag_value}_weeklag")) %>%
    pull(estimate)
  
  se_main <- intx_df %>%
    filter(term == glue::glue("{predictor}_{lag_value}_weeklag")) %>%
    pull(std.error)
  
  beta_intx <- intx_df %>%
    filter(term == glue::glue("{predictor}_{lag_value}_weeklag:{intx_var}{intx_var_name}")) %>%
    pull(estimate)
  
  se_intx <- intx_df %>%
    filter(term == glue::glue("{predictor}_{lag_value}_weeklag:{intx_var}{intx_var_name}")) %>%
    pull(std.error)
  
  p_intx <- intx_df %>%
    filter(term == glue::glue("{predictor}_{lag_value}_weeklag:{intx_var}{intx_var_name}")) %>%
    pull(p.value)
  
  # Calculate estimates and SEs first
  base_estimate <- ifelse(grepl("precip", predictor), beta_main * 500, beta_main)
  base_se <- ifelse(grepl("precip", predictor), se_main * 500, se_main)
  base_IR <- exp(base_estimate)
  base_CI_lower <- exp(base_estimate - 1.96 * base_se)
  base_CI_upper <- exp(base_estimate + 1.96 * base_se)
  
  # For interaction term - proper SE calculation using delta method
  intx_estimate <- ifelse(grepl("precip", predictor), 
                          (beta_main + beta_intx) * 500, 
                          beta_main + beta_intx)
  
  # Get variance-covariance matrix for proper SE calculation
  vcov_matrix <- vcov(model_intx)
  main_term <- glue::glue("{predictor}_{lag_value}_weeklag")
  intx_term <- glue::glue("{predictor}_{lag_value}_weeklag:{intx_var}{intx_var_name}")
  
  # Calculate SE for sum of coefficients using delta method
  var_sum <- vcov_matrix[main_term, main_term] + 
    vcov_matrix[intx_term, intx_term] + 
    2 * vcov_matrix[main_term, intx_term]
  
  intx_se <- ifelse(grepl("precip", predictor), 
                    sqrt(var_sum) * 500, 
                    sqrt(var_sum))
  
  intx_IR <- exp(intx_estimate)
  intx_CI_lower <- exp(intx_estimate - 1.96 * intx_se)
  intx_CI_upper <- exp(intx_estimate + 1.96 * intx_se)
  
  # Create data frames
  base_IR_df <- data.frame(
    intx_var = intx_var,
    predictor = predictor,
    intx_var_name = ifelse(intx_var == "comm_type", "Highway", "Normal ENSO periods"),
    IR_CI = paste0(sprintf("%.2f", base_IR), " (", 
                   sprintf("%.2f", base_CI_lower), " - ", 
                   sprintf("%.2f", base_CI_upper), ")"),
    p_value = " ",
    stringsAsFactors = FALSE
  )
  
  intx_IR_df <- data.frame(
    intx_var = intx_var, 
    predictor = predictor,
    intx_var_name = intx_var_name,
    IR_CI = paste0(sprintf("%.2f", intx_IR), " (", 
                   sprintf("%.2f", intx_CI_lower), " - ", 
                   sprintf("%.2f", intx_CI_upper), ")"),
    p_value = sprintf("%.3f", p_intx),
    stringsAsFactors = FALSE
  )
  
  if (intx_var_name == "La Niña periods") {
    df <- intx_IR_df
  } else {
    df <- rbind(base_IR_df, intx_IR_df)
  }
  
  df <- df %>%
    mutate(p_value = ifelse(p_value == "0.000", "<0.001", p_value))
  
  return(df)
}

# Fit GLM models
glm_temp_max_comm_type <- fit_intx_glm(data = nonlagged_data,
                                       predictor = "temp_wk_max",
                                       lag_value = 9,
                                       covariate = "precip_wk_total",
                                       intx_var = "comm_type",
                                       intx_var_name = "Riverine")

glm_temp_min_comm_type <- fit_intx_glm(data = nonlagged_data,
                                       predictor = "temp_wk_min",
                                       lag_value = 7,
                                       covariate = "precip_wk_total",
                                       intx_var = "comm_type",
                                       intx_var_name = "Riverine")

glm_precip_comm_type <- fit_intx_glm(data = nonlagged_data,
                                     predictor = "precip_wk_total",
                                     lag_value = 11,
                                     covariate = "temp_wk_min",
                                     intx_var = "comm_type",
                                     intx_var_name = "Riverine")

glm_temp_max_el_nino <- fit_intx_glm(data = nonlagged_data,
                                     predictor = "temp_wk_max",
                                     lag_value = 9,
                                     covariate = "precip_wk_total",
                                     intx_var = "enso_period",
                                     intx_var_name = "El Niño periods")

glm_temp_min_el_nino <- fit_intx_glm(data = nonlagged_data,
                                     predictor = "temp_wk_min",
                                     lag_value = 7,
                                     covariate = "precip_wk_total",
                                     intx_var = "enso_period",
                                     intx_var_name = "El Niño periods")

glm_precip_el_nino <- fit_intx_glm(data = nonlagged_data,
                                   predictor = "precip_wk_total",
                                   lag_value = 11,
                                   covariate = "temp_wk_min",
                                   intx_var = "enso_period",
                                   intx_var_name = "El Niño periods")

glm_temp_max_la_nina <- fit_intx_glm(data = nonlagged_data,
                                     predictor = "temp_wk_max",
                                     lag_value = 9,
                                     covariate = "precip_wk_total",
                                     intx_var = "enso_period",
                                     intx_var_name = "La Niña periods")
glm_temp_max_enso <- rbind(glm_temp_max_el_nino,
                           glm_temp_max_la_nina)

glm_temp_min_la_nina <- fit_intx_glm(data = nonlagged_data,
                                     predictor = "temp_wk_min",
                                     lag_value = 7,
                                     covariate = "precip_wk_total",
                                     intx_var = "enso_period",
                                     intx_var_name = "La Niña periods")
glm_temp_min_enso <- rbind(glm_temp_min_el_nino,
                           glm_temp_min_la_nina)

glm_precip_la_nina <- fit_intx_glm(data = nonlagged_data,
                                   predictor = "precip_wk_total",
                                   lag_value = 11,
                                   covariate = "temp_wk_min",
                                   intx_var = "enso_period",
                                   intx_var_name = "La Niña periods")
glm_precip_enso <- rbind(glm_precip_el_nino,
                           glm_precip_la_nina)

# Load base results for predictor and incidence ranges ------------------------
## Community Type ----------------------------------------------------
base_precip_total_comm <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-precip-wk-total-highway.RDS")) %>% 
    mutate(comm_type = "Highway"),
  readRDS(file = paste0(here::here(), "/results/base results/base-precip-wk-total-riverine.RDS")) %>% 
    mutate(comm_type = "Riverine")) %>%
  mutate(incidence_10k = incidence*10000) %>%
  rename("intx_var_name" = "comm_type")

base_temp_max_comm <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-wk-max-highway.RDS")) %>% 
    mutate(comm_type = "Highway"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-wk-max-riverine.RDS")) %>% 
    mutate(comm_type = "Riverine")) %>%
  mutate(incidence_10k = incidence*10000) %>%
  rename("intx_var_name" = "comm_type")

base_temp_min_comm <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-wk-min-highway.RDS")) %>% 
    mutate(comm_type = "Highway"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-wk-min-riverine.RDS")) %>% 
    mutate(comm_type = "Riverine")) %>%
  mutate(incidence_10k = incidence*10000)  %>%
  rename("intx_var_name" = "comm_type")

## ENSO Data --------------------------------------------------------------
temp_min_base_enso <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-min-enso-normal.RDS")) %>%
    mutate(enso_period = "Normal ENSO periods"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-min-enso-el-nino.RDS")) %>%
    mutate(enso_period = "El Niño periods"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-min-enso-la-nina.RDS")) %>%
    mutate(enso_period = "La Niña periods")) %>%
  mutate(incidence_10k = incidence*10000) %>%
  rename("intx_var_name" = "enso_period")

temp_max_base_enso <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-max-enso-normal.RDS")) %>%
    mutate(enso_period = "Normal ENSO periods"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-max-enso-el-nino.RDS")) %>%
    mutate(enso_period = "El Niño periods"),
  readRDS(file = paste0(here::here(), "/results/base results/base-temp-max-enso-la-nina.RDS")) %>%
    mutate(enso_period = "La Niña periods")) %>%
  mutate(incidence_10k = incidence*10000) %>%
  rename("intx_var_name" = "enso_period")

precip_total_base_enso <- rbind(
  readRDS(file = paste0(here::here(), "/results/base results/base-precip-total-enso-normal.RDS")) %>%
    mutate(enso_period = "Normal ENSO periods"),
  readRDS(file = paste0(here::here(), "/results/base results/base-precip-total-enso-el-nino.RDS")) %>%
    mutate(enso_period = "El Niño periods"),
  readRDS(file = paste0(here::here(), "/results/base results/base-precip-total-enso-la-nina.RDS")) %>%
    mutate(enso_period = "La Niña periods")) %>%
  mutate(incidence_10k = incidence*10000) %>%
  rename("intx_var_name" = "enso_period")

# Table 1: Community type results
comm_type_table <- rbind(
  glm_temp_min_comm_type,
  glm_temp_max_comm_type,
  glm_precip_comm_type
) %>%
  left_join(
    rbind(
      base_temp_min_comm %>% 
        group_by(intx_var_name) %>%
        summarize(predictor = "temp_wk_min",
                  intx_var = "comm_type",
                  min_inc = min(incidence_10k),
                  max_inc = max(incidence_10k),
                  min_pred = sprintf("%0.1f", min(temp_wk_min_7_weeklag)),
                  max_pred = sprintf("%0.1f", max(temp_wk_min_7_weeklag))),
      base_temp_max_comm %>% 
        group_by(intx_var_name) %>%
        summarize(predictor = "temp_wk_max",
                  intx_var = "comm_type",
                  min_inc = min(incidence_10k),
                  max_inc = max(incidence_10k),
                  min_pred = sprintf("%0.1f", min(temp_wk_max_9_weeklag)),
                  max_pred = sprintf("%0.1f", max(temp_wk_max_9_weeklag))),
      base_precip_total_comm %>%
        group_by(intx_var_name) %>%
        summarize(predictor = "precip_wk_total",
                  intx_var = "comm_type",
                  min_inc = min(incidence_10k),
                  max_inc = max(incidence_10k),
                  min_pred = round(min(precip_wk_total_11_weeklag)),
                  max_pred = round(max(precip_wk_total_11_weeklag)))
    ),
    by = c("predictor", "intx_var", "intx_var_name")
  ) %>%
  mutate(
    predictor = case_when(
      predictor == "temp_wk_min" ~ "Minimum temperature",
      predictor == "temp_wk_max" ~ "Maximum temperature",
      predictor == "precip_wk_total" ~ "Precipitation"
    ),
    pred_range = paste0(min_pred, " - ", max_pred),
    inc_range = paste0(sprintf("%.1f", min_inc), " - ", sprintf("%.1f", max_inc))
  ) %>%
  select(predictor, intx_var_name, pred_range, inc_range, IR_CI, p_value) %>%
  rename(
    "Community type" = intx_var_name,
    "Predictor range" = pred_range,
    "Incidence range" = inc_range,
    "Incidence Ratio (95% CI)" = IR_CI,
    "Interaction p-value" = p_value
  )

# Table 2: ENSO results
enso_table <- rbind(
  glm_temp_min_enso,
  glm_temp_max_enso,
  glm_precip_enso
) %>%
  left_join(
    rbind(
      temp_min_base_enso %>% 
        group_by(intx_var_name) %>%
        summarize(predictor = "temp_wk_min",
                  intx_var = "enso_period",
                  min_inc = min(incidence_10k),
                  max_inc = max(incidence_10k),
                  min_pred = sprintf("%0.1f", min(temp_wk_min_7_weeklag)),
                  max_pred = sprintf("%0.1f", max(temp_wk_min_7_weeklag))),
      temp_max_base_enso %>% 
        group_by(intx_var_name) %>%
        summarize(predictor = "temp_wk_max",
                  intx_var = "enso_period",
                  min_inc = min(incidence_10k),
                  max_inc = max(incidence_10k),
                  min_pred = sprintf("%0.1f", min(temp_wk_max_9_weeklag)),
                  max_pred = sprintf("%0.1f", max(temp_wk_max_9_weeklag))),
      precip_total_base_enso %>% 
        group_by(intx_var_name) %>%
        summarize(predictor = "precip_wk_total",
                  intx_var = "enso_period",
                  min_inc = min(incidence_10k),
                  max_inc = max(incidence_10k),
                  min_pred = round(min(precip_wk_total_11_weeklag)),
                  max_pred = round(max(precip_wk_total_11_weeklag)))
    ),
    by = c("predictor", "intx_var", "intx_var_name")
  ) %>%
  mutate(
    predictor = case_when(
      predictor == "temp_wk_min" ~ "Minimum temperature",
      predictor == "temp_wk_max" ~ "Maximum temperature",
      predictor == "precip_wk_total" ~ "Precipitation"
    ),
    intx_var_name = case_when(
      intx_var_name == "Normal ENSO periods" ~ "Normal",
      intx_var_name == "El Niño periods" ~ "El Nino",
      intx_var_name == "La Niña periods" ~ "La Nina"
    ),
    pred_range = paste0(min_pred, " - ", max_pred),
    inc_range = paste0(sprintf("%.1f", min_inc), " - ", sprintf("%.1f", max_inc))
  ) %>%
  select(predictor, intx_var_name, pred_range, inc_range, IR_CI, p_value) %>%
  rename(
    "ENSO period" = intx_var_name,
    "Predictor range" = pred_range,
    "Incidence range" = inc_range,
    "Incidence Ratio (95% CI)" = IR_CI,
    "Interaction p-value" = p_value
  )

# Write tables to CSV
write.csv(comm_type_table, paste0(tables_path, "table-1.csv"), row.names = FALSE)
write.csv(enso_table, paste0(tables_path, "table-2.csv"), row.names = FALSE)
