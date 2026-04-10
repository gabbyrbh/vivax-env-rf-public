#!/bin/bash
# ==============================================================================
# @Organization : Benjamin-Chung Lab
# @Project      : Environmental Risk Factors for Vivax Malaria in Peru
# @File         : run-gam-temp-min-comm-type.sh
# @Description  : SLURM submission script for gam-temp-min-comm-type.R (DLNM-GAM
#                 for minimum temperature stratified by community type) on
#                 Sherlock HPC.
# ==============================================================================

#SBATCH --job-name=run-gam-temp-min-comm-type
#SBATCH --begin=now
#SBATCH --dependency=singleton
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=32
#SBATCH --mem=40G
#SBATCH --mem=40G
#SBATCH --output=00-run-gam-temp-min-comm-type.out
#SBATCH --time=20:00:00

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

R CMD BATCH --no-save gam-temp-min-comm-type.R gam-temp-min-comm-type.out