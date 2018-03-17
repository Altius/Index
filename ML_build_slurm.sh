#!/bin/bash
#SBATCH --job-name=ML_build.%a
#SBATCH -n 1                       # Number of cores
#SBATCH --array=1-5000
#SBATCH --error=error/ML_build.%a.txt
#SBATCH --output=out/ML_build.%a.txt
#SBATCH --partition=queue1

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
FILEPATH="/net/seq/data/projects/ENCODE3_publications/erynes/ENCODEpaper2017/SubsetOf665/varWidthPeaks/chunks_10kb_min500elementsPerChunk/";
R CMD BATCH --no-save --no-restore "--args chunknum=${SLURM_ARRAY_TASK_ID} filepath='${FILEPATH}'" \
  code_build.R Rout/output_build_chunk_${SLURM_ARRAY_TASK_ID}.Rout


