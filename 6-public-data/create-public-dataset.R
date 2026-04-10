# ------------------------------------------------------------------------------
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : create-public-dataset.R
# @Description  : Scramble community IDs using a random key and drop geographic
#                 identifiers (long, lat, weather_comm) to produce a public
#                 dataset safe for publication. Saves the anonymized dataset and
#                 the ID mapping key separately.
#
# Outputs:
#   - public-data/vivax-env-erf-public.RDS   : anonymized dataset
#   - public-data/comm-id-key.csv            : mapping of original -> scrambled IDs
#                                              (keep private or publish separately)
# ------------------------------------------------------------------------------
rm(list = ls())
source(paste0(here::here(), "/0-config.R"))
library(dplyr)

set.seed(42)  # Set seed for reproducibility of the scramble key

# Load data -------------------------------------------------------------------
data <- readRDS(paste0(box_path_flame_erf, "non-lagged-analysis-data_ext.RDS"))

# Join comm_type before scrambling so original IDs match ---------------------
comm_type <- read.csv(paste0(box_path_flame_erf, "all_district_centroids_comm_type.csv")) %>%
  mutate(comm_id   = sprintf("%03d", comm_id),
         comm_type = ifelse(comm_type == "", "riverine", comm_type)) %>%
  select(comm_id, comm_type)

data <- data %>% left_join(comm_type, by = "comm_id")

# Create scrambled ID key -----------------------------------------------------
original_ids <- sort(unique(data$comm_id))

# Randomly shuffle the original IDs to create anonymized IDs
# Result: original comm "001" maps to some other existing ID, preventing
# geographic inference without the key
scrambled_ids <- sample(original_ids, size = length(original_ids), replace = FALSE)

id_key <- data.frame(
  original_comm_id = original_ids,
  public_comm_id   = scrambled_ids,
  stringsAsFactors = FALSE
)

# Apply key to dataset --------------------------------------------------------
public_data <- data %>%
  left_join(id_key, by = c("comm_id" = "original_comm_id")) %>%
  mutate(comm_id = public_comm_id) %>%
  select(-public_comm_id) %>%
  # Drop geographic identifiers
  select(-long, -lat, -weather_comm)

# Verify no original IDs leaked through and dimensions are preserved ----------
stopifnot(nrow(public_data) == nrow(data))
stopifnot(ncol(public_data) == ncol(data) - 3)  # dropped long, lat, weather_comm (comm_type already joined)
stopifnot(all(sort(unique(public_data$comm_id)) == sort(scrambled_ids)))
cat("Rows:    ", nrow(public_data), "\n")
cat("Columns: ", ncol(public_data), "\n")
cat("Columns in public dataset:\n")
print(names(public_data))
cat("Unique communities (should be", length(original_ids), "):",
    length(unique(public_data$comm_id)), "\n")

# Save outputs ----------------------------------------------------------------
output_dir <- paste0(here::here(), "/6-public-data/output/")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(public_data, file = paste0(output_dir, "vivax-env-erf-public.csv"), row.names = FALSE)
write.csv(id_key,      file = paste0(output_dir, "comm-id-key.csv"),          row.names = FALSE)

cat("\nSaved public dataset to:  ", paste0(output_dir, "vivax-env-erf-public.csv"), "\n")
cat("Saved ID key to:          ", paste0(output_dir, "comm-id-key.csv"), "\n")
cat("\nRemember: keep comm-id-key.csv private or publish separately from the dataset.\n")
