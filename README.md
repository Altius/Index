Method for generating a master list / Index of DNaseI hypersensitivity sites.
All code, implementation details and design by Wouter Meuleman and Eric Rynes.

## How to run within Docker:
This approach is not recommended, as it will execute command serially, which will take an enormous amount of time for any real-life dataset.

1. Run the Docker script

`./run_from_scratch.sh`

## How to run with SLURM:

1. Generates chunked versions of the master list, for each of ~5000 chunks (run in parallel)

`sbatch ML_build_slurm.sh`
wait for jobs to finish

2. Resolves overlaps for alternative master list versions where no overlap is required

`sbatch ML_overlap_slurm.sh`
wait for jobs to finish

3. Generates the final concatenated versions, as well as browser tracks

`./code_gen_masterlist.sh <ID/DATE> <numchunks>`

`ID/DATE` will be the name of the masterlist -- the latest I generated is called 'WM20180313'.

`numchunks` is the number of genomic chunks that was used to process the R code in parallel.

## Prerequisites:

This code depends on the availability of individual "chunk" files, each containing hotspot2 peak calls for a subset of the genome, across all samples of interest.
Starting from hotspot2, the following steps need to be taken to obtain these:

TODO Eric Rynes

## Files of interest:

| File | Purpose |
| --- | --- |
| `code_ML.R` | common routines |
| `code_build.R` | code used for converting a genomic chunk of peak calls into tentative DHSs |
| `ML_build_slurm.sh` | SLURM submission script for `code_build.R` |
| `code_overlap.R` | code used to detect and resolve overlapping elements, if so desired |
| `ML_overlap_slurm.sh` | SLURM submission script for `code_overlap.R` |
| `code_gen_masterlist.sh` | code used to concatenate the output of all chunks and generate browser tracks |



