#!/bin/bash

set -e -o pipefail

dir=$1
chrom_file=$(readlink -f "$2")
in_file_list=$(readlink -f "$3")
outname=${4:-masterlist}

starting_dir=$(pwd)
here=$(dirname "$(readlink -f "$0")")

if [[ $# -lt 3 ]] ; then
    echo "Usage: $0 output_dir chrom_file in_file_list [outname]"
    exit 2
fi

mkdir -p "$dir"
dir=$(readlink -f "$dir")

rout_dir="$dir"/Rout
mkdir -p "$rout_dir"

# Collate peaks together
OUTFILE=allPeaks.starch
{
    # Merge peaks
    echo -e "Collating peaks..."
    declare -a infiles
    readarray -t infiles < "$in_file_list"
    {
        for f in "${infiles[@]}" ; do
            unstarch "$f"
        done
    } \
        | sort-bed - \
        | starch - \
        > "${dir}/$OUTFILE"
    if ! mkdir -p "${dir}" ; then
        echo -e "Failed to create directory ${dir}."
        exit 2
    fi

    # Chunk peaks up
    echo -e '"Chunking" collated peaks for eventual use in master list creation...'
    cd "${dir}" || exit 2
    unstarch "$OUTFILE" | awk \
        -v minChunkSep=10000 \
        -v minChunkSep2=4000 \
        -v minPerChunk=500 \
        -v maxPerChunk=100000 \
        -f "$here"/chunk_bed.awk
}

cd "$starting_dir"

# Do something
i=1
numchunks=$(find "$dir" -maxdepth 1 -name 'chunk*.bed' | wc -l)
cd "$here"
for chunkfile in "$dir"/chunk*.bed ; do
    Rscript "$here"/code_build.R \
        chunknum=$i \
        filepath="'$dir'" \
        workdir="'$dir'" \
        sourcedir="'$here'" \
        > "$rout_dir/output_build_chunk_$i.Rout" \
        2>&1
    ((i++))
done

# Do something else
for i in $(seq 1 "$numchunks") ; do
    Rscript "$here"/code_overlap.R \
        chunknum="$i" \
        workdir="'$dir'" \
        sourcedir="'$here'" \
        > "$rout_dir/output_overlap_chunk_$i.Rout" \
        2>&1
done
cd "$dir"

# Generate masterlist
"$here"/code_gen_masterlist.sh "$(basename "$outname")" "$chrom_file" .

# Make a matrix
# TODO: This comes later
# ./code_construct_matrix.sh "$(basename "$outname")" "$numchunks"
