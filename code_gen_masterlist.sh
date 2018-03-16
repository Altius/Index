##########################################################################################################
#                                                                                                        #
# Method for generating a master list / Index of DNaseI hypersensitivity sites.                          #
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.                         #
#                                                                                                        #
# Version: WM20180313                                                                                    #
#                                                                                                        #
##########################################################################################################

if [[ $# -eq 0 ]] ; then
    echo 'Provide the name/version of the file, e.g. WM20180301'
    exit 1
fi

NAME=$1;

TYPES="all nonovl_any nonovl_core";
for TYPE in ${TYPES}; do
  echo $TYPE;

  FILE_CHUNKIDS="masterlist_DHSs_${NAME}_${TYPE}_chunkIDs.txt";
  FILE_INDEXIDS="masterlist_DHSs_${NAME}_${TYPE}_indexIDs.txt";
  FILE_BED="masterlist_DHSs_${NAME}_${TYPE}.bed";
  FILE_BIGBED="masterlist_DHSs_${NAME}_${TYPE}.bb";

  ### Create final master list
  if ! [[ -f "$FILE_CHUNKIDS" ]] ; then
    echo "Concatenating DHS chunks"
    cat DHSs_${TYPE}/* | sort-bed - > ${FILE_CHUNKIDS}
  fi
  
  #### Generate label mapping
  #if ! [[ -f "masterlist_DHSs_${NAME}_all_chunkIDs2indexIDs.txt" && ${TYPE}=="all" ]] ; then
  #  echo "Generating unique DHS identifiers"
  #  ./run_name_master_list.sh ${FILE_CHUNKIDS};
  #fi
  #
  #### Apply mapping
  #if ! [[ -f "$FILE_INDEXIDS" ]] ; then
  #  echo "Mapping chunk identifiers to final DHS identifiers"
  #  awk 'BEGIN{ FS = OFS = "\t" }
  #    FNR == NR { split($0, f, /\t/); map[f[2]] = f[1]; next } 
  #    { if ($4 in map) { $4 = map[$4] } } 
  #    { print }' masterlist_DHSs_${NAME}_all_chunkIDs2indexIDs.txt ${FILE_CHUNKIDS} > ${FILE_INDEXIDS}
  #fi
  FILE_INDEXIDS=${FILE_CHUNKIDS}

  ### Create browser loadable BED12 files
  if ! [[ -f "$FILE_BED" ]] ; then
    echo "Constructing BED12 file"
    #echo "browser position chr6:26020208-26022677" > ${FILE_BED}
    #echo "track name='Master list DHSs ${NAME} ${TYPE}' description='Master list DHSs ${NAME} ${TYPE}' visibility=2 useScore=1" >> ${FILE_BED}
    cat ${FILE_INDEXIDS} | awk '
    function round(x) { if(x=="NA") { return 0 } else { return int(x + 0.5) } }
    BEGIN { OFS="\t"; } 
    {
      if ($10 == "NA" || $11 == "NA") {
        thickStart=$2
        thickEnd=$3
        blockCount=1
        blockSizes=$3-$2
        blockStarts=0
      } else {
        $9=($9 <= $2 ? $2+1 : $9)
        $9=($9 >= $3 ? $3-1 : $9)
        $10=($10 <= $2 ? $2+1 : $10)
        $11=($11 >= $3 ? $3-1 : $11)
    
        thickStart= $9-1
        thickEnd  = $9+1
    
        blockCount=1
        blockSizes=1
        blockStarts=0
    
        blockSize2=$9-$10
        blockStart2=$10-$2
        if (blockSize2 != 0 && blockStart2 < $3-$2-1) { 
          blockCount+=1; 
          blockSizes=blockSizes","blockSize2; 
          blockStarts=blockStarts","blockStart2 
        }
    
        blockSize3=$11-$9
        blockStart3=$9-$2
        if (blockSize3 != 0 && blockStart3 < $3-$2-1) { 
          blockCount+=1; 
          blockSizes=blockSizes","blockSize3; 
          blockStarts=blockStarts","blockStart3 
        }
    
        blockSize4=1
        blockStart4=$3-$2-1
        blockCount+=1;
        blockSizes=blockSizes","blockSize4; 
        blockStarts=blockStarts","blockStart4
      }
  
      score=round(log($5+1)/log(10)*500)
      score=(score > 1000 ? 1000 : score)
      print $1, $2, $3, $4, score, ".", thickStart, thickEnd, "0,0,0", blockCount, blockSizes, blockStarts
    }' > ${FILE_BED}
  fi
  
  if ! [[ -f "$FILE_BIGBED" ]] ; then
    echo "Converting BED file to BIGBED file"
    bedToBigBed -type=bed12 ${FILE_BED} <(tail -n +2 hg38.genome) ${FILE_BIGBED}
  fi

  echo "Load this track in the UCSC browser using the following:"
  echo "track type=bigBed name=master_list_${NAME}_${TYPE} useScore=1 visibility=2 itemRgb='On' bigDataUrl=https://encode:collabor8@resources.altius.org/~meuleman/ML_tracks/${FILE_BIGBED}"

done;


