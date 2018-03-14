##########################################################################################################
#                                                                                                        #
# Method for generating a master list / Index of DNaseI hypersensitivity sites.                          #
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.                         #
#                                                                                                        #
# Version: WM20180313                                                                                    #
#                                                                                                        #
##########################################################################################################

# Compare two lists:

if [[ $# -lt 2 ]] ; then
    echo 'Provide the name/version of each of the two files you like to compare, e.g. WM20180301'
    exit 1
fi

NAME1=$1;
NAME2=$2;

LIST1="masterlist_DHSs_${NAME1}_chunkIDs.txt"
LIST2="masterlist_DHSs_${NAME2}_chunkIDs.txt"

diff -y --suppress-common-lines \
  <(cut -f 1-3 ${LIST1} | awk 'function round(x){return int(x + 0.5)}{print $1, round($2/100)*100, round($3/100)*100}') \
  <(cut -f 1-3 ${LIST2} | awk 'function round(x){return int(x + 0.5)}{print $1, round($2/100)*100, round($3/100)*100}') | \
  awk 'BEGIN{OFS="\t"}{if (($0 != "") && ($1 != ">")) print $1, $2-500, $3+500}' | bedops -m - > \
  diffs_${NAME1}_vs_${NAME2}.bed




