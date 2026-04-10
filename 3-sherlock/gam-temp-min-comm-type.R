# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : gam-temp-min-comm-type.R
# @Description  : DLNM-GAM for minimum temperature stratified by community type
#                 (highway / riverine) with 1,000 bootstrap CIs.
#                 Runs on Stanford Sherlock HPC. Outputs to results/bs/.
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

data_path <- "/home/groups/jadebc/vivax-env-rf/data/"
figure_path <- "/home/groups/jadebc/vivax-env-rf/figures/"
results_path <- "/home/groups/jadebc/vivax-env-rf/results/"

# Load data ---------------------------------------------------------------
nonlagged_data <- readRDS(paste0(data_path,
                                 "non-lagged-analysis-data_ext.RDS")) %>%
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

# Read in comm_type csv
comm_type <- read.csv(paste0(data_path, "all_district_centroids_comm_type.csv")) %>%
  mutate(comm_type = ifelse(comm_type == "", "riverine", comm_type)) %>%
  mutate(comm_id = sprintf("%03d", comm_id)) %>%
  mutate(comm_id = as.factor(comm_id)) %>%
  mutate(comm_id = relevel(comm_id, ref = median_pop_comm))

nonlagged_data <- nonlagged_data %>% left_join(comm_type, by = "comm_id")

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
                    s(oni_index, k = 3) +
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

# Main analysis function
run_analysis <- function(data, predictor, lag_value, covariate, comm_type_value) {
  
  cat(glue("Starting analysis for predictor: {predictor}, comm_type: {comm_type_value}, lag: {lag_value}\n"))
  
  # Filter data by community type
  filtered_data <- data %>% filter(comm_type == comm_type_value)
  
  # Create lagged data
  lagged_data <- lag_data(filtered_data, predictor, covariate, lag_value)
  
  # Prepare analysis data
  analysis_data <- lagged_data %>% 
    filter(!is.na(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")))) %>%
    mutate(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) := as.numeric(!!sym(glue::glue("{predictor}_{lag_value}_weeklag"))),
           !!sym(glue::glue("{covariate}_{lag_value}_weeklag")) := as.numeric(!!sym(glue::glue("{covariate}_{lag_value}_weeklag"))))
  
  # Filter extreme values
  analysis_data <- analysis_data %>%
    filter(!!sym(glue::glue("{predictor}_{lag_value}_weeklag")) > quantile(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]], p = 0.05) &
             !!sym(glue::glue("{predictor}_{lag_value}_weeklag")) < quantile(analysis_data[[glue("{predictor}_{lag_value}_weeklag")]], p = 0.95))
  
  cat("Fitting GAM model...\n")
  
  # Fit the main model
  model <- gam(data = analysis_data,
               as.formula(glue::glue(
                 "n_cases ~
                     s({predictor}_{lag_value}_weeklag, k = 3) + 
                     s({covariate}_{lag_value}_weeklag, k = 3) + 
                     s(oni_index, k = 3) +
                     comm_id + 
                     years_since_itn_delivery +
                     year")),
               offset = log(population),
               family = "poisson")
  
  # Create prediction sequences
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
  
  # Compute base result
  cat("Computing base result...\n")
  base_result <- compute_marginal_incidence(model, analysis_data, pred_seq, predictor, lag_value)

  # Create filename components
  predictor_clean <- gsub("_", "-", predictor)
  comm_type_clean <- gsub("_", "-", comm_type_value)

  # Save base result
  base_filename <- glue("base-{predictor_clean}-{comm_type_clean}.RDS")
  saveRDS(base_result, file = paste0(results_path, base_filename))
  cat(glue("Base result saved as: {base_filename}\n"))

  # Perform bootstrap
  cat("Starting bootstrap...\n")
  boot_results <- bootstrap_incidence_parallel(analysis_data, pred_seq, predictor, lag_value, covariate, cores_per_worker, n_boot = 1000)
  
  # Save bootstrap result
  bootstrap_filename <- glue("bs-{predictor_clean}-{comm_type_clean}.RDS")
  saveRDS(boot_results, file = paste0(results_path, bootstrap_filename))
  cat(glue("Bootstrap result saved as: {bootstrap_filename}\n"))
  
  cat(glue("Analysis completed for predictor: {predictor}, comm_type: {comm_type_value}\n\n"))
  
  return(list(base_result = base_result, 
              boot_results = boot_results,
              base_filename = base_filename,
              bootstrap_filename = bootstrap_filename))
}

# Set up parallel backend once
cl <- makeCluster(r_workers)
registerDoParallel(cl)
cat("Parallel backend ready with", getDoParWorkers(), "workers\n")

# Run analysis for min temp with riverine communities
results_temp_max_riverine <- run_analysis(
  data = nonlagged_data,
  predictor = "temp_wk_min",
  lag_value = 7,
  covariate = "precip_wk_total",
  comm_type_value = "riverine"
)

# Run analysis for min temp with highway communities
results_temp_max_road <- run_analysis(
  data = nonlagged_data,
  predictor = "temp_wk_min", 
  lag_value = 7,
  covariate = "precip_wk_total",
  comm_type_value = "highway"
)

# Clean up
stopCluster(cl)