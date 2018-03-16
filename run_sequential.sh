#!/bin/bash

set -e -o pipefail

dir=$1
chrom_file=$2
outname=${3:-masterlist}

outdir=Rout
mkdir -p "$outdir"
numchunks=$(find "$dir" -maxdepth 1 -name 'chunk*.bed' | wc -l)

# Do something
i=1
for chunkfile in "$dir"/chunk*.bed ; do
    R CMD BATCH \
        --no-save --no-restore \
        "--args chunknum=$i $chunkfile" \
        ./code_build.R \
        "$outdir/output_build_chunk_$i.Rout"
    ((i++))
done

# Do something else
for i in $(seq 1 "$numchunks") ; do
    R CMD BATCH --no-save --no-restore "--args chunknum=${i}" \
        ./code_overlap.R "$outdir/output_overlap_chunk_${i}.Rout"
done

# Generate masterlist
./code_gen_masterlist.sh "$(basename "$outname")" "$chrom_file"

# Copy output files out of the container
cp -r \
    DHS* \
    Rout \
    masterlist_DHSs* \
    peaks_* \
    "$dir"
# Make a matrix
# TODO: This comes later
# ./code_construct_matrix.sh "$(basename "$outname")" "$numchunks"
