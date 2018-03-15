#!/bin/bash

dir=$1
outname=${2:-masterlist}

outdir=Rout
mkdir -p "$outdir"
numchunks=$(find "$dir" -maxdepth 1 -name 'chunk*.bed' | wc -l)

# Do something
i=1
for chunkfile in "$dir"/chunk*.bed ; do
    R CMD BATCH \
        --no-save --no-restore \
        "--args chunknum=$i filepath=$chunkfile" \
        ./code_build.R \
        "$outdir/output_build_chunk_$i.Rout"
    (i++)
done

# Do something else
for i in $(seq 1 "$numchunks") ; do
    R CMD BATCH --no-save --no-restore "--args chunknum=${i}" \
        code_overlap.R "Rout/output_overlap_chunk_${i}.Rout"
done

# Generate masterlist
./code_gen_masterlist.sh "$outname"

# Make a matrix
./code_construct_matrix.sh "$outname" "$numchunks"
