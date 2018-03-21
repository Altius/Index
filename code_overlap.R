##########################################################################################################
#                                                                                                        #
# Method for generating a master list / Index of DNaseI hypersensitivity sites.                          #
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.                         #
#                                                                                                        #
# Version: WM20180313                                                                                    #
#                                                                                                        #
##########################################################################################################

library(caTools)
#source("code_ML.R")

######################################################################################################################################

args=(commandArgs(TRUE))
if (length(args)==0) {
  stop("No arguments supplied.")
} else {
  eval(parse(text=args[[1]])) # parse first argument: chunknum
  eval(parse(text=args[[2]])) # parse second argument: workdir
  eval(parse(text=args[[3]])) # parse third argument: sourcedir
}

setwd(workdir)

source(paste(sourcedir,"code_ML.R",sep="/"))

chunk <- paste("chunk", sprintf("%04d", chunknum), ".bed", sep="");
print(chunk)

options(scipen=20)
#options(warn=0)
options(warn=2)

dir.create("DHSs_nonovl_core", showWarnings=FALSE, recursive=TRUE)
dir.create("peaks_nonovl_core", showWarnings=FALSE, recursive=TRUE)

dir.create("DHSs_nonovl_any", showWarnings=FALSE, recursive=TRUE)
dir.create("peaks_nonovl_any", showWarnings=FALSE, recursive=TRUE)

dir.create("DHSs_nonovl_stats", showWarnings=FALSE, recursive=TRUE)

######################################################################################################################################

## Example line:
# chr11   52983767        52983861        id-708944       0.0622013       52983811        LN4866  52983810
#chunks <- dir("peaks", pattern="^chunk.*.bed");

#for (chunk in chunks) {
  # Load in data per chunk, each separated by at least 10kb
  peaks <- read.delim(paste("peaks_all", chunk, sep="/"), header=FALSE, as.is=T)
  #colnames(peaks) <- c("seqname", "start", "end", "ID", "score", "density_summit", "sampleID", "wavelet_summit")
  colnames(peaks) <- c("seqname", "start", "end", "sampleID", "score", "density_summit", "wavelet_summit", "ID")

  DHSs <- read.delim(paste("DHSs_all", chunk, sep="/"), header=FALSE, as.is=T)
  colnames(DHSs) <- c("seqname", "start", "end", "ID", "score", "numsamples", "numpeaks", "width", "summit", "disp_low", "disp_high");

  # Resolve overlaps across DHS summit core regions
  DHSs_nonovl_core <- DHSs; peaks_nonovl_core <- peaks;
  num <- 1;
  while(num > 0) {
    ovl_rm <- merge_overlap(DHSs_nonovl_core, peaks_nonovl_core, type="core", force=TRUE)
    DHSs_nonovl_core <- ovl_rm$DHSs
    peaks_nonovl_core <- ovl_rm$peaks
    num <- ovl_rm$num
    print(num)
  }
  write.table(DHSs_nonovl_core, file=paste("DHSs_nonovl_core", chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(peaks_nonovl_core, file=paste("peaks_nonovl_core", chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

  # Resolve overlaps across DHS overall
  DHSs_nonovl_any <- DHSs; peaks_nonovl_any <- peaks;
  num <- 1;
  while(num > 0) {
    ovl_rm <- merge_overlap(DHSs_nonovl_any, peaks_nonovl_any, type="any", force=TRUE)
    DHSs_nonovl_any <- ovl_rm$DHSs
    peaks_nonovl_any <- ovl_rm$peaks
    num <- ovl_rm$num
    print(num)
  }
  write.table(DHSs_nonovl_any, file=paste("DHSs_nonovl_any", chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(peaks_nonovl_any, file=paste("peaks_nonovl_any", chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

  # Register level of overlap to original list
  peak_map <- unique(cbind(peaks$ID, peaks_nonovl_core$ID))
  peak_map_list <- peak_map[,2]; names(peak_map_list) <- peak_map[,1]
  DHSs$ovl_core <- peak_map_list[DHSs$ID];

  peak_map <- unique(cbind(peaks$ID, peaks_nonovl_any$ID))
  peak_map_list <- peak_map[,2]; names(peak_map_list) <- peak_map[,1]
  DHSs$ovl_any <- peak_map_list[DHSs$ID];

  write.table(DHSs, file=paste("DHSs_nonovl_stats", chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
#}

######################################################################################################################################

