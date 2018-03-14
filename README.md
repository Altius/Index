Method for generating a master list / Index of DNaseI hypersensitivity sites.
All code, implementation details and design by Wouter Meuleman and Eric Rynes.

# How to run:

## Generates chunked versions of the master list, for each of ~5000 chunks (run in parallel)
`sbatch --array=1-5000 ML_build_slurm.sh`
wait for jobs to finish

## Resolves overlaps for alternative master list versions where no overlap is required
`sbatch --array=1-5000 ML_overlap_slurm.sh`
wait for jobs to finish

## Generates the final concatenated versions, as well as browser tracks
`./code_gen_masterlist.sh <ID/DATE>`
ID/DATE will be the name of the masterlist -- the latest I generated is called 'WM20180313'

## Optional: generate presence/absence and signal matrices for the three versions of the master list chunks
`sbatch --array=1-5000 ML_matrix_slurm.sh`
wait for jobs to finish
## Optional: Generate the final concatenated versions.
`./code_construct_matrix.sh <ID/DATE>`


# Files of interest:

```
code_ML.R

ML_build_slurm.sh
code_build.R

ML_overlap_slurm.sh
code_overlap.R
code_gen_masterlist.sh

ML_matrix_slurm.sh
code_matrix.R
code_construct_matrix.sh

code_compare_masterlists.sh
```


