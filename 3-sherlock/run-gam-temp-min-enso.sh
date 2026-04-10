#!/bin/bash
# ==============================================================================
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : run-gam-temp-min-enso.sh
# @Description  : SLURM submission script for gam-temp-min-enso.R (DLNM-GAM
#                 for minimum temperature stratified by ENSO period) on
#                 Sherlock HPC.
# ==============================================================================

#SBATCH --job-name=run-gam-temp-min-enso
#SBATCH --begin=now
#SBATCH --dependency=singleton
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --mem=64G
#SBATCH --output=00-run-gam-temp-max-enso.out
#SBATCH --time=30:00:00

cd /home/groups/jadebc/vivax-env-rf/
# Export the CPU count for R to use
export SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
  
module purge 

module load math
module load devel
module load system
module load curl
module load physics
module load gcc
module load R/4.2.0

R CMD BATCH --no-save gam-temp-min-enso.R gam-temp-min-enso.out