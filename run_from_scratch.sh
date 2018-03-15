#! /bin/bash

set -e -u -o pipefail

if [[ $# != 4 && $# != 5 ]]; then
  echo -e "Usage:  $0 fileOfBAMfiles outdir chrom_sizes.bed mappable.starch [min_varWidth_peak_width]"
  exit 2
fi

FILE_OF_BAM_FILES=$1
OUTDIR=$2
CHROM_FILTER=$3
MAPPABLE_REGIONS=$4

if [ $# = 5 ]; then
  MIN_VARWIDTH_PEAK_WIDTH=$5
else
  MIN_VARWIDTH_PEAK_WIDTH="20"
fi

# Sanity check on input filenames.
linenum=1
while read -r line
do
  fname=$line
  if [ ! -s "$fname" ]; then
    echo -e "Error:  on line $linenum of $FILE_OF_BAM_FILES, file not found, or the file is empty ($fname)."
    exit 2
  fi
  ((linenum++))
done < "$FILE_OF_BAM_FILES"

mkdir -p "$OUTDIR"

# Input files and parameter values
#CENTER_SITES=/home/erynes/topics/ENCODEpaper2017/CenterSitesFiles/centerSites_halfWin100bp_K36mapOnlyMinusBlacklist_chrs1-22XY.hg38.starch
#MAPPABLE_REGIONS=/home/erynes/topics/ENCODEpaper2017/GRCh38_36merMappableOnly_minusBlacklist_chrs1-22XY.starch
#CHROM_FILTER=/home/erynes/topics/ENCODEpaper2017/chromLengths_GRCh38_only1-22XY.bed3
CALL_THRESHOLD=1.0 # write all sites to disk, including sites with FDR = 1
HOTSPOT_FDR_THRESHOLD=0.0010 # 0.1% (not 1%)
CENTER_SITES="$OUTDIR/centersites.starch"

function exe_check(){
  for exe in "$@"; do
    if [ ! -x "$exe" ]; then
      echo -e "Error:  $exe was not found, or execution privileges for it are not set."
      exit 2
    fi
  done
}

# Executables and scripts
HOTSPOT2=$(which hotspot2.sh)
EXTRACT=$(which extractCenterSites.sh)

exe_check "$HOTSPOT2" "$EXTRACT"

# Create center sites
echo "Creating center sites"
"$EXTRACT" -c "$CHROM_FILTER" -o "$CENTER_SITES"

ls -l "$CENTER_SITES"


# Run hotspot2.sh
while read -r BAM
do
  bname=$(basename "$BAM" .bam)
  varWidthSpec="varWidth_${MIN_VARWIDTH_PEAK_WIDTH}_${bname}"
  {
    unset TMPDIR
    "$HOTSPOT2" -c "$CHROM_FILTER" -C "$CENTER_SITES" -M "$MAPPABLE_REGIONS" -p "$varWidthSpec" -f "$HOTSPOT_FDR_THRESHOLD" -F "$CALL_THRESHOLD" "$BAM" "$OUTDIR"
  }
done < "$FILE_OF_BAM_FILES"


OUTFILE=allPeaks.starch
{
    echo -e "Collating peaks..."
    bedops -u "${OUTDIR}"/*.peaks.starch \
       | starch - \
       > "${OUTDIR}/$OUTFILE"
    if ! mkdir -p "${OUTDIR}/subdir" ; then
       echo -e "Failed to create directory ${OUTDIR}/subdir."
       exit 2
     fi
    echo -e '"Chunking" collated peaks for eventual use in master list creation...'
    cd "${OUTDIR}/subdir" || exit 2
    unstarch "../${OUTFILE}" \
       | awk -v minChunkSep=10000 -v minChunkSep2=4000 -v minPerChunk=500 -v maxPerChunk=100000 \
          'BEGIN{n=1; nWritten=0; outf="chunk0001.bed"}{
             if(NR > 1){
                if(\$1 != prevChr || \$2-prev3 > minChunkSep){
                   if(nWritten >= minPerChunk){
                      outf = sprintf("chunk%04d.bed", ++n);
                      nWritten = 0;
                    }
                  }
                else{
                   if(\$2-prev3 > minChunkSep2 && nWritten > maxPerChunk){
                      outf = sprintf("chunk%04d.bed",++n);
                      nWritten=0;
                    }
                  }
                }
             print \$0 > outf;
             nWritten++; prevChr = \$1; prev3 = \$3;
           }'
    echo -e "Done!"
}

# Run the masterlist stuff
./run_sequential.sh "$OUTDIR" "$OUTDIR"
