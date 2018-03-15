#!/bin/bash
##########################################################################################################
#                                                                                                        #
# Method for generating a master list / Index of DNaseI hypersensitivity sites.                          #
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.                         #
#                                                                                                        #
# Version: WM20180313                                                                                    #
#                                                                                                        #
##########################################################################################################


if [[ $# -le 2 ]] ; then
    echo 'Provide the name/version of the file, e.g. WM20180301, and the number of chunks'
    exit 1
fi

NAME=$1;
numchunks=$2

rm -f "matrix_bin_all_${NAME}.txt.gz"
for i in $(seq -f "%04g" 1 "$numchunks"); do
  echo "$i"
  if [[ -f "DHSs_all/chunk${i}.bed" ]] ; then
    #paste <(cut -f 1-3 DHSs_all/chunk${i}.bed) matrix_bin_all/chunk${i}.bed | gzip -c - >> matrix_bin_all_${NAME}.txt.gz
    paste <(cut -f 4 "DHSs_all/chunk${i}.bed") "matrix_bin_all/chunk${i}.bed" | gzip -c - >> "matrix_bin_all_${NAME}.txt.gz"
  fi
done

rm -f "matrix_bin_nonovl_core_${NAME}.txt.gz"
for i in $(seq -f "%04g" 1 "$numchunks"); do
  echo "$i"
  if [[ -f "DHSs_nonovl_core/chunk${i}.bed" ]] ; then
    #paste <(cut -f 1-3 DHSs_nonovl_core/chunk${i}.bed) matrix_bin_nonovl_core/chunk${i}.bed | gzip -c - >> matrix_bin_nonovl_core_${NAME}.txt.gz
    paste <(cut -f 4 "DHSs_nonovl_core/chunk${i}.bed") "matrix_bin_nonovl_core/chunk${i}.bed" | gzip -c - >> "matrix_bin_nonovl_core_${NAME}.txt.gz"
  fi
done

rm -f "matrix_bin_nonovl_any_${NAME}.txt.gz"
for i in $(seq -f "%04g" 1 "$numchunks"); do
  echo "$i"
  if [[ -f "DHSs_nonovl_any/chunk${i}.bed" ]] ; then
    #paste <(cut -f 1-3 DHSs_nonovl_any/chunk${i}.bed) matrix_bin_nonovl_any/chunk${i}.bed | gzip -c - >> matrix_bin_nonovl_any_${NAME}.txt.gz
    paste <(cut -f 4 "DHSs_nonovl_any/chunk${i}.bed") "matrix_bin_nonovl_any/chunk${i}.bed" | gzip -c - >> "matrix_bin_nonovl_any_${NAME}.txt.gz"
  fi
done
