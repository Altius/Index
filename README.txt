# Method for generating a master list / Index of DNaseI hypersensitivity sites.
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.
#####################################################################################################

### How to run:

## Generates chunked versions of the master list, for each of ~5000 chunks (run in parallel)
sbatch --array=1-5000 ML_build_slurm.sh 
# wait for jobs to finish

## Resolves overlaps for alternative master list versions where no overlap is required
sbatch --array=1-5000 ML_overlap_slurm.sh 
# wait for jobs to finish

## Generates the final concatenated versions, as well as browser tracks
./code_gen_masterlist.sh <ID/DATE>
# ID/DATE will be the name of the masterlist -- the latest I generated is called 'WM20180313'

## Optional: generate presence/absence and signal matrices for the three versions of the master list chunks
sbatch --array=1-5000 ML_matrix_slurm.sh 
# wait for jobs to finish
## Optional: Generate the final concatenated versions.
./code_construct_matrix.sh <ID/DATE>

#####################################################################################################

### Files of interest:

README.txt
code_build.R
code_compare_masterlists.sh
code_construct_matrix.sh
code_gen_masterlist.sh
code_matrix.R
code_ML.R
code_overlap.R
code_stats.R
ML_overlap_slurm.sh
ML_build_slurm.sh
ML_matrix_slurm.sh

##########################################################################################################

### Misc notes:

For naming DHSs, I adapted this from Eric Haugen:
 /home/ehaugen/work/encode3/naming/src/run_name_master_list.sh input_file.txt

## Check to see if number of accounted peaks is the same (it is!):

[meuleman@sched0 WM20180220_detailed_ML_analysis]$ cat peaks_all/chunk*| wc -l
68698924
[meuleman@sched0 WM20180220_detailed_ML_analysis]$ cat peaks_nonovl_any/chunk*| wc -l
68698924



