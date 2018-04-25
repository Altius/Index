#!/bin/bash

##########################################################################################################
#                                                                                                        #
# Method for generating a master list / Index of DNaseI hypersensitivity sites.                          #
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.                         #
#                                                                                                        #
# Version: WM20180313                                                                                    #
#                                                                                                        #
##########################################################################################################

usage () {
    echo -e "Usage:  $0 output_dir outputID chromSizes.bed fileOfPeakfileNames.txt partitionName memSize"
    echo -e "where output_dir is where all messages and intermediate output will be written,"
    echo -e "outputID is a \"version ID\" to be used in the output filenames (e.g. containing the date and number of samples),"
    echo -e "chromSizes.bed is a 3-column 0-based BED file containing the lengths of the relevant chromosomes,"
    echo -e "and fileOfPeakfileNames.txt contains paths to variable-width peak files, one file / biological sample"
    echo -e "per line, from which the Index will be constructed."
    echo -e "output_dir will be created if it does not already exist."
    echo -e "partitionName is the name of the cluster on which the SLURM jobs will run (--partition=partitionName),"
    echo -e "and memSize is the amount of memory to require for collating the files in fileOfPeakfileNames.txt (--mem=memSize)."
    echo -e "Output files will be written to the directory from which the script is invoked, i.e. the current directory."
    exit 2
}

if [[ $# != 6 ]]; then
    usage
fi

SOURCE_DIR=`realpath $0 | awk '{len=split($1,x,"/");printf("/%s",x[2]);for(i=3;i<len;i++){printf("/%s",x[i])}printf("\n")}'`
export SOURCE_DIR
OUTDIR=$1
OUTPUT_VERSION_ID="$2"
export OUTPUT_VERSION_ID
CHROM_SIZES=$3
export CHROM_SIZES
PEAKFILES=$4
export PEAKFILES
PARTITION=$5
export PARTITION
memSize=$6
export memSize

if [ ! -s $CHROM_SIZES ]; then
    echo -e "Error:  File $CHROM_SIZES was not found, or it is empty."
    usage
fi
if [ ! -s $PEAKFILES ]; then
    echo -e "Error:  File $PEAKFILES was not found, or it is empty."
    usage
fi

R_BUILD_SCRIPT=${SOURCE_DIR}/code_build.R
export R_BUILD_SCRIPT
if [ ! -s "$R_BUILD_SCRIPT" ]; then
    echo -e "Error:  Required file $R_BUILD_SCRIPT was not found, or it is empty."
    exit 2
fi
R_OVERLAP_SCRIPT=${SOURCE_DIR}/code_overlap.R
export R_OVERLAP_SCRIPT
if [ ! -s "$R_OVERLAP_SCRIPT" ]; then
    echo -e "Error:  Required file $R_OVERLAP_SCRIPT was not found, or it is empty."
    exit 2
fi
GEN_ML_SCRIPT=${SOURCE_DIR}/code_gen_masterlist.sh
export GEN_ML_SCRIPT
if [ ! -s "$GEN_ML_SCRIPT" ]; then
    echo -e "Error:  Required file $GEN_ML_SCRIPT was not found, or it is empty."
    exit 2
fi
CHUNK_SCRIPT=${SOURCE_DIR}/chunk_bed.awk
export CHUNK_SCRIPT
if [ ! -s "$CHUNK_SCRIPT" ]; then
    echo -e "Error:  Required file $CHUNK_SCRIPT was not found, or it is empty."
    exit 2
fi

mkdir -p $OUTDIR
if [ $? != "0" ]; then
    echo -e "Error:  Failed to create output directory $OUTDIR."
    usage
fi
OUTDIR=`realpath $OUTDIR`
export OUTDIR
STDERR_DIR=${OUTDIR}/errorMessages
export STDERR_DIR
STDOUT_DIR=${OUTDIR}/outputMessages
export STDOUT_DIR
ROUT_DIR=${OUTDIR}/R_Messages
export ROUT_DIR
mkdir -p $STDERR_DIR
if [ $? != "0" ]; then
    echo -e "Error:  Failed to create directory $STDERR_DIR for general messages sent to stderr."
    exit 2
fi
mkdir -p $STDOUT_DIR
if [ $? != "0" ]; then
    echo -e "Error:  Failed to create directory $STDOUT_DIR for general messages sent to stdout."
    exit 2
fi
mkdir -p $ROUT_DIR
if [ $? != "0" ]; then
    echo -e "Error:  Failed to create directory $ROUT_DIR for messages sent from the R code."
    exit 2
fi

# Collate the peak files and "chunk" them into genomic "islands" that can be processed in parallel.
COLLATED_INPUT=`pwd`
COLLATED_INPUT="${COLLATED_INPUT}/allPeaks.bed"
export COLLATED_INPUT
jobName=CollateAndChunk
jobID1=$(sbatch --parsable --partition=$PARTITION --export=ALL --job-name=$jobName --output=${STDOUT_DIR}/${jobName}.o%j --error=${STDERR_DIR}/${jobName}.e%j --mem=$memSize <<EOF1
#! /bin/bash
   module load bedops
   if ! which bedops >/dev/null 2>/dev/null ; then
      echo -e "Error:  bedops must be available on computing cluster $PARTITION; it was not found."
      exit 2
   fi
   echo -e "Collating peaks..."
   for f in \`cat $PEAKFILES\` ; do unstarch \$f ; done \
      | sort-bed - \
      > $COLLATED_INPUT
   if [ \$? != "0" ]; then
      echo -e "An error occurred while collating the peaks in $PEAKFILE."
      exit 2
   fi
   cd $OUTDIR
   echo -e "\"Chunking\" collated peaks into smaller files..."
   awk -v minChunkSep=10000 -v minChunkSep2=4000 -v minPerChunk=500 -v maxPerChunk=100000 -f $CHUNK_SCRIPT $COLLATED_INPUT
   if [ \$? != "0" ]; then
      echo -e "An error occurred while executing $CHUNK_SCRIPT."
      exit 2
   else
      rm -f $COLLATED_INPUT
   fi
EOF1
     )

dependency2="afterok:${jobID1}"

jobName="ProcessChunks"
jobID2=$(sbatch --parsable --dependency=$dependency2 --partition=$PARTITION --export=ALL --job-name=$jobName --output=${STDOUT_DIR}/${jobName}.o%j --error=${STDERR_DIR}/${jobName}.e%j --mem=$memSize <<"EOF2"
#! /bin/bash
   numJobs=`ls -1 ${OUTDIR}/chunk*bed | wc -l`
   jobName3a="ML_build.%a"
   jobID3a=$(sbatch --parsable --partition=$PARTITION --export=ALL --job-name=${jobName3a} -n 1 --array=1-$numJobs --output=${STDOUT_DIR}/${jobName3a}.o%j --error=${STDERR_DIR}/${jobName3a}.e%j --mem=$memSize <<"EOF3a"
#! /bin/bash
      module load R/3.3.3
      if ! which R >/dev/null 2>/dev/null ; then
         echo -e "Error:  R must be available on computing cluster $PARTITION; it was not found."
         exit 2
      fi
      echo -e "Building the initial Index..."
      R CMD BATCH --no-save --no-restore "--args chunknum=${SLURM_ARRAY_TASK_ID} filepath='${OUTDIR}' workdir='${OUTDIR}' sourcedir='${SOURCE_DIR}'" \
         $R_BUILD_SCRIPT $ROUT_DIR/output_build_chunk_${SLURM_ARRAY_TASK_ID}.Rout
EOF3a
      )
   dependency3b="afterok:${jobID3a}"
   jobName3b="ML_overlap.%a"
   jobID3b=$(sbatch --parsable --partition=$PARTITION --dependency=$dependency3b --export=ALL --job-name=${jobName3b} -n 1 --array=1-$numJobs --output=${STDOUT_DIR}/${jobName3b}.o%j --error=${STDERR_DIR}/${jobName3b}.e%j --mem=$memSize <<"EOF3b"
#! /bin/bash
      module load R/3.3.3
      if ! which R >/dev/null 2>/dev/null ; then
         echo -e "Error:  R must be available on computing cluster $PARTITION; it was not found."
         exit 2
      fi
      echo -e "Resolving overlaps in the initial Index..."
      R CMD BATCH --no-save --no-restore "--args chunknum=${SLURM_ARRAY_TASK_ID} workdir='${OUTDIR}' sourcedir='${SOURCE_DIR}'" \
         $R_OVERLAP_SCRIPT $ROUT_DIR/output_overlap_chunk_${SLURM_ARRAY_TASK_ID}.Rout
EOF3b
      )
   dependency3c="afterok:${jobID3b}"
   jobName3c=Finishing
   jobID3c=$(sbatch --parsable --partition=$PARTITION --dependency=$dependency3c --export=ALL --job-name=${jobName3c} --output=${STDOUT_DIR}/${jobName3c}.o%j --error=${STDERR_DIR}/${jobName3c}.e%j --mem=$memSize <<"EOF3c"
#! /bin/bash
      module load bedops
      if ! which bedops >/dev/null 2>/dev/null ; then
         echo -e "Error:  bedops must be available on computing cluster $PARTITION; it was not found."
         exit 2
      fi
      module load kentutil
      if ! which bedToBigBed >/dev/null 2>/dev/null ; then
         echo -e "Error:  bedToBigBed must be available on computing cluster $PARTITION; it was not found."
         exit 2
      fi
      echo -e "Creating the final Index files..."
      $GEN_ML_SCRIPT $OUTPUT_VERSION_ID $CHROM_SIZES $OUTDIR
      if [ $? != "0" ]; then
         echo -e "An error occurred while executing $GEN_ML_SCRIPT."
         exit 2
      fi
EOF3c
      )
EOF2
      )

exit 0
