# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : gam-temp-max-enso.R
# @Description  : DLNM-GAM for maximum temperature stratified by ENSO period
#                 with 1,000 bootstrap CIs. Runs on Stanford Sherlock HPC.
#                 Outputs bootstrap results to results/bs/.
# ------------------------------------------------------------------------------
rm(list = ls())
library(mgcv)
library(lme4)
library(lmtest)
library(sandwich)
library(gKRLS)
library(splines)
library(glue)
library(data.table)
library(broom)
library(boot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(parallel)
library(doParallel)
library(foreach)
library(furrr)
library(zoo)

# Control parallelization properly
slurm_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = detectCores()))
r_workers <- max(1, min(8, slurm_cpus %/% 4))  # Cap at 8 workers, use 1/4 of cores
cores_per_worker <- slurm_cpus %/% r_workers

cat("Slurm allocated CPUs:", slurm_cpus, "\n")
cat("Using R workers:", r_workers, "\n")
cat("Cores per worker:", cores_per_worker, "\n")

# Set threading environment variables
Sys.setenv(OMP_NUM_THREADS = as.character(cores_per_worker))
Sys.setenv(MKL_NUM_THREADS = as.character(cores_per_worker))
Sys.setenv(OPENBLAS_NUM_THREADS = as.character(cores_per_worker))

# Paths (defined relative to project root — edit in 0-config.R if needed)
public_data_path <- paste0(here::here(), "/6-public-data/output/")
results_path     <- paste0(here::here(), "/results/")
bs_path          <- paste0(here::here(), "/bs/")
dir.create(bs_path, showWarnings = FALSE, recursive = TRUE)

# Load data ---------------------------------------------------------------
nonlagged_data <- read.csv(paste0(public_data_path, "vivax-env-erf-public.csv")) %>%
  rename(time = "study_week") %>%
  mutate(population = as.numeric(population),
         n_cases = as.numeric(n_cases),
         comm_id = as.factor(comm_id)) %>%
  filter(!is.na(population)) %>%
  mutate(year = factor(as.character(year), levels = c("2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))) %>%
  filter(year != "2020") %>%
  filter(year != "2021") %>%
  filter(population > 15) %>%
  arrange(time) %>%
  mutate(
    oni_above_threshold = abs(oni_index) >= 0.5,
    rolling_sum = rollapply(oni_above_threshold, width = 5, FUN = sum, align = "right", fill = NA),
    abnormal_oni = ifelse(rolling_sum >= 5, 1, 0),
    el_nino_threshold = oni_index >= 0.5,
    el_nino_rolling_sum = rollapply(el_nino_threshold, width = 5, FUN = sum, align = "right", fill = NA),
    el_nino = ifelse(el_nino_rolling_sum >= 5, 1, 0),
    la_nina_threshold = oni_index <= -0.5,
    la_nina_rolling_sum = rollapply(la_nina_threshold, width = 5, FUN = sum, align = "right", fill = NA),
    la_nina = ifelse(la_nina_rolling_sum >= 5, 1, 0)
  ) %>%
  mutate(oni_period = case_when(
    abnormal_oni == 0 ~ "Normal",
    el_nino == 1 ~ "El Niño",
    la_nina == 1 ~ "La Niña"
  )) %>%
  ungroup()

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

normal_data <- nonlagged_data %>% filter(oni_period == "Normal")
el_nino_data <- nonlagged_data %>% filter(oni_period == "El Niño")
la_nina_data <- nonlagged_data %>% filter(oni_period == "La Niña")

# Function to lag data for a specific variable and lag
lag_data <- function(data, predictor, covariate, lag_value) {
  lagged_data <- data %>%
    group_by(comm_id) %>%
    mutate(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) := lag(!!sym(predictor), n = lag_value),
           !!sym(glue::glue("{covariate}_{lag_value}_weeklag")) := lag(!!sym(covariate), n = lag_value)) %>%
    ungroup()
  
  return(lagged_data)
}

# Function to compute marginal incidence
compute_marginal_incidence <- function(model, data, pred_seq, predictor, lag_value) {
  pred_col <- glue("{predictor}_{lag_value}_weeklag")

  results <- lapply(pred_seq, function(pval) {
    newdata <- data
    newdata[[pred_col]] <- pval
    preds <- predict(model, newdata = newdata, type = "response")
    data.frame(
      pred_col_val = pval,
      incidence = mean(preds, na.rm = TRUE)
    )
  })

  plot_data <- do.call(rbind, results)
  colnames(plot_data)[1] <- pred_col
  return(plot_data)
}

# Fixed Bootstrap function with all variables properly exported
bootstrap_incidence_parallel <- function(data, pred_seq, predictor, lag_value, covariate, cores_per_worker, n_boot = 1000, seed = 123) {

  set.seed(seed)

  boot_results <- foreach(i = 1:n_boot,
                          .combine = rbind,
                          .packages = c("mgcv", "dplyr", "glue"),
                          .export = c("compute_marginal_incidence", "predictor", "lag_value", "covariate", "cores_per_worker", "pred_seq")) %dopar% {

                            Sys.setenv(OMP_NUM_THREADS = as.character(cores_per_worker))
                            Sys.setenv(MKL_NUM_THREADS = as.character(cores_per_worker))
                            Sys.setenv(OPENBLAS_NUM_THREADS = as.character(cores_per_worker))

                            boot_indices <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
                            boot_data <- data[boot_indices, ]

                            boot_model <- gam(data = boot_data,
                                              as.formula(glue::glue(
                                                "n_cases ~
                    s({predictor}_{lag_value}_weeklag, k = 3) +
                    s({covariate}_{lag_value}_weeklag, k = 3) +
                    comm_id +
                    years_since_itn_delivery +
                    year")),
                                              offset = log(population),
                                              family = "poisson")

                            result <- compute_marginal_incidence(boot_model, boot_data, pred_seq, predictor, lag_value)
                            result$boot_id <- i
                            return(result)
                          }

  return(boot_results)
}

# Set up parallel backend once
cl <- makeCluster(r_workers)
registerDoParallel(cl)
cat("Parallel backend ready with", getDoParWorkers(), "workers\n")

# Analysis for NORMAL data
data = normal_data
predictor = "temp_wk_max"
lag_value = 9
covariate = "precip_wk_total"

cat("Starting analysis for normal data...\n")

lagged_data <- lag_data(data, predictor, covariate, lag_value)

analysis_data <- lagged_data %>% 
  filter(!is.na(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")))) %>%
  mutate(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) := as.numeric(!!sym(glue::glue("{predictor}_{lag_value}_weeklag"))),
         !!sym(glue::glue("{covariate}_{lag_value}_weeklag")) := as.numeric(!!sym(glue::glue("{covariate}_{lag_value}_weeklag"))))

analysis_data <- analysis_data %>%
  filter(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) > quantile(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]], p = 0.05) &
           !!sym(glue::glue("{predictor}_{lag_value}_weeklag")) < quantile(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]], p = 0.95))

cat("Fitting GAM model for normal data...\n")
model <- gam(data = analysis_data,
             as.formula(glue::glue(
               "n_cases ~
                   s({predictor}_{lag_value}_weeklag, k = 3) + 
                   s({covariate}_{lag_value}_weeklag, k = 3) + 
                   comm_id + 
                   years_since_itn_delivery +
                   year")),
             offset = log(population),
             family = "poisson")

# Create a sequence of values for the predictor
pred_seq <-
  if(predictor == "precip_wk_total") {
    seq(
      min(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]),
      max(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]),
      length = 8
    )
  } else {
    seq(
      min(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]), 
      max(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]), 
      by = 0.2
    )
  }
base_result <- compute_marginal_incidence(model, analysis_data, pred_seq, predictor, lag_value)
saveRDS(base_result, file = paste0(results_path, "base-temp-max-enso-normal.RDS"))

# Perform bootstrap for normal data
cat("Starting bootstrap for normal data...\n")
boot_results <- bootstrap_incidence_parallel(analysis_data, pred_seq, predictor, lag_value, covariate, cores_per_worker, n_boot = 1000)
saveRDS(boot_results, file = paste0(bs_path, "bs-temp-max-enso-normal.RDS"))
cat("Normal data bootstrap completed.\n")

# Analysis for El Niño data
data = el_nino_data
predictor = "temp_wk_max"
lag_value = 9
covariate = "precip_wk_total"

cat("Starting analysis for El Niño data...\n")

lagged_data <- lag_data(data, predictor, covariate, lag_value)

analysis_data <- lagged_data %>% 
  filter(!is.na(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")))) %>%
  mutate(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) := as.numeric(!!sym(glue::glue("{predictor}_{lag_value}_weeklag"))),
         !!sym(glue::glue("{covariate}_{lag_value}_weeklag")) := as.numeric(!!sym(glue::glue("{covariate}_{lag_value}_weeklag"))))

analysis_data <- analysis_data %>%
  filter(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) > quantile(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]], p = 0.05) &
           !!sym(glue::glue("{predictor}_{lag_value}_weeklag")) < quantile(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]], p = 0.95))

cat("Fitting GAM model for abnormal data...\n")
model <- gam(data = analysis_data,
             as.formula(glue::glue(
               "n_cases ~
                   s({predictor}_{lag_value}_weeklag, k = 3) + 
                   s({covariate}_{lag_value}_weeklag, k = 3) + 
                   comm_id + 
                   years_since_itn_delivery +
                   year")),
             offset = log(population),
             family = "poisson")

# Create a sequence of values for the predictor
pred_seq <-
  if(predictor == "precip_wk_total") {
    seq(
      min(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]),
      max(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]),
      length = 8
    )
  } else {
    seq(
      min(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]), 
      max(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]), 
      by = 0.2
    )
  }
base_result <- compute_marginal_incidence(model, analysis_data, pred_seq, predictor, lag_value)
saveRDS(base_result, file = paste0(results_path, "base-temp-max-enso-el-nino.RDS"))

# Perform bootstrap for El Niño data
cat("Starting bootstrap for El Niño data...\n")
boot_results <- bootstrap_incidence_parallel(analysis_data, pred_seq, predictor, lag_value, covariate, cores_per_worker, n_boot = 1000)
saveRDS(boot_results, file = paste0(bs_path, "bs-temp-max-enso-el-nino.RDS"))
cat("El Niño data bootstrap completed.\n")

# Analysis for La Niña data
data = la_nina_data
predictor = "temp_wk_max"
lag_value = 9
covariate = "precip_wk_total"

cat("Starting analysis for La Niña data...\n")

lagged_data <- lag_data(data, predictor, covariate, lag_value)

analysis_data <- lagged_data %>% 
  filter(!is.na(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")))) %>%
  mutate(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) := as.numeric(!!sym(glue::glue("{predictor}_{lag_value}_weeklag"))),
         !!sym(glue::glue("{covariate}_{lag_value}_weeklag")) := as.numeric(!!sym(glue::glue("{covariate}_{lag_value}_weeklag"))))

analysis_data <- analysis_data %>%
  filter(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) > quantile(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]], p = 0.05) &
           !!sym(glue::glue("{predictor}_{lag_value}_weeklag")) < quantile(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]], p = 0.95))

cat("Fitting GAM model for abnormal data...\n")
model <- gam(data = analysis_data,
             as.formula(glue::glue(
               "n_cases ~
                   s({predictor}_{lag_value}_weeklag, k = 3) + 
                   s({covariate}_{lag_value}_weeklag, k = 3) + 
                   comm_id + 
                   years_since_itn_delivery +
                   year")),
             offset = log(population),
             family = "poisson")

# Create a sequence of values for the predictor
pred_seq <-
  if(predictor == "precip_wk_total") {
    seq(
      min(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]),
      max(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]),
      length = 8
    )
  } else {
    seq(
      min(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]), 
      max(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]]), 
      by = 0.2
    )
  }
base_result <- compute_marginal_incidence(model, analysis_data, pred_seq, predictor, lag_value)
saveRDS(base_result, file = paste0(results_path, "base-temp-max-enso-la-nina.RDS"))

# Perform bootstrap for La Niña data
cat("Starting bootstrap for La Niña data...\n")
boot_results <- bootstrap_incidence_parallel(analysis_data, pred_seq, predictor, lag_value, covariate, cores_per_worker, n_boot = 1000)
saveRDS(boot_results, file = paste0(bs_path, "bs-temp-max-enso-la-nina.RDS"))
cat("La Niña data bootstrap completed.\n")