#!/bin/bash
# ==============================================================================
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : run-gam-temp-max-enso.sh
# @Description  : SLURM submission script for gam-temp-max-enso.R (DLNM-GAM
#                 for maximum temperature stratified by ENSO period) on
#                 Sherlock HPC.
# ==============================================================================

#SBATCH --job-name=run-gam-temp-max-enso
#SBATCH --begin=now
#SBATCH --dependency=singleton
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --mem=64G
#SBATCH --output=00-run-gam-temp-max-enso.out
#SBATCH --time=30:00:00

# Export the CPU count for R to use
export SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

cd /home/groups/jadebc/vivax-env-rf/
  
module purge 

module load math
module load devel
module load system
module load curl
module load physics
module load gcc
module load R/4.2.0

R CMD BATCH --no-save gam-temp-max-enso.R gam-temp-max-enso.out