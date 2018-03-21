Method for generating a master list / Index of DNaseI hypersensitivity sites.

All code, implementation details and design by Wouter Meuleman and Eric Rynes. Docker implementation by Jemma Nelson.

## How to run within Docker:
This approach is not recommended, as it will execute command serially, which will take an enormous amount of time for any real-life dataset.

1. Build the docker image from this repo, with `docker build . --tag=masterlist`
2. Change to the directory that contains your input files
3. Run `docker run --mount type=bind,source="$(pwd)",target=/data masterlist run_sequential.sh my_output_dir chrom_sizes.bed listOfFiles.txt`

`my_output_dir` will be created inside of the working directory. `chrom_sizes.bed` is a standard bed file containing the sizes of the chromosomes. `listOfFiles.txt` is a list of peak files (which are `starch`-formatted), containing one peak file per line. A relative path should be used to a file inside the working directory - absolute paths are not supported, due to the way the files are mounted in the Docker container.

## How to run with SLURM:

From any directory, execute

`[path to]ML_build_slurm.sh workdir outputID chrom_sizes.bed listOfFiles.txt partitionName memSize`

where

* messages and intermediate files will be written into `workdir` (it will be created if it doesn't already exist)
* `outputID` is an identifier to be used in the output filenames (e.g., containing the date, number of samples, organism, etc.)
* `chrom_sizes.bed` is a 0-based 3-column BED file containing the lengths of the relevant chromosomes
* `listOfFiles.txt` is a plain text file containing the paths to variable-width peak files, one per biological sample, one path per line
* `partitionName` is the name of the cluster on which the SLURM jobs will run (`--partition=partitionName`)
* `memSize` is the amount of memory to require for each step, best if tailored to the needs of the collation of peak files (`--mem=memSize`)

Output files will be written to the current directory.  Messages and intermediate files will be written into `workdir`.

## Files of interest:

| File | Purpose |
| --- | --- |
| `chunk_bed.awk` | script for partitioning the peaks into genomic islands or "chunks" that can be processed in parallel |
| `code_ML.R` | common routines |
| `code_build.R` | code used for converting a genomic chunk of peak calls into tentative DHSs |
| `code_overlap.R` | code used to detect and resolve overlapping elements, if so desired |
| `code_gen_masterlist.sh` | code used to concatenate the output of all chunks and generate browser tracks |
| `ML_build_slurm.sh` | SLURM submission script which executes each of the above in sequence |




