#!/bin/bash
#SBATCH --job-name=ML_matrix.%a
#SBATCH -n 1                       # Number of cores
#SBATCH --error=error/ML_matrix.%a.txt
#SBATCH --output=out/ML_matrix.%a.txt

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
R CMD BATCH --no-save --no-restore "--args chunknum=${SLURM_ARRAY_TASK_ID} type='all'" \
  code_matrix.R Rout/output_matrix_all_chunk_${SLURM_ARRAY_TASK_ID}.Rout
R CMD BATCH --no-save --no-restore "--args chunknum=${SLURM_ARRAY_TASK_ID} type='nonovl_core'" \
  code_matrix.R Rout/output_matrix_nonovl_core_chunk_${SLURM_ARRAY_TASK_ID}.Rout
R CMD BATCH --no-save --no-restore "--args chunknum=${SLURM_ARRAY_TASK_ID} type='nonovl_any'" \
  code_matrix.R Rout/output_matrix_nonovl_any_chunk_${SLURM_ARRAY_TASK_ID}.Rout


