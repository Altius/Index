#!/bin/bash
#SBATCH --job-name=ML_overlap.%a
#SBATCH -n 1                       # Number of cores
#SBATCH --error=error/ML_overlap.%a.txt
#SBATCH --output=out/ML_overlap.%a.txt

##########################################################################################################
#                                                                                                        #
# Method for generating a master list / Index of DNaseI hypersensitivity sites.                          #
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.                         #
#                                                                                                        #
# Version: WM20180313                                                                                    #
#                                                                                                        #
##########################################################################################################

mkdir -p error out Rout
module load R/3.3.3
R CMD BATCH --no-save --no-restore "--args chunknum=${SLURM_ARRAY_TASK_ID}" \
  code_overlap.R Rout/output_overlap_chunk_${SLURM_ARRAY_TASK_ID}.Rout


