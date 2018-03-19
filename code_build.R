#!/usr/bin/env Rscript
##########################################################################################################
#                                                                                                        #
# Method for generating a master list / Index of DNaseI hypersensitivity sites.                          #
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.                         #
#                                                                                                        #
# Version: WM20180313                                                                                    #
#                                                                                                        #
##########################################################################################################

library(caTools)
source("code_ML.R")

######################################################################################################################################

args=(commandArgs(TRUE))
if (length(args)==0) {
  stop("No arguments supplied.")
} else {
  eval(parse(text=args[[1]])) # parse first argument: chunknum
  eval(parse(text=args[[2]])) # parse second argument: filepath
  eval(parse(text=args[[3]])) # parse second argument: workdir
}

setwd(workdir)

chunk <- paste("chunk", sprintf("%04d", chunknum), ".bed", sep="");
print(chunk)

options(scipen=20)
#options(warn=0)
options(warn=2)

######################################################################################################################################

## Example line:
# chr11   52983767        52983861        id-708944       0.0622013       52983811        LN4866  52983810

dir.create("DHSs_all", showWarnings=FALSE, recursive=TRUE)
dir.create("peaks_all", showWarnings=FALSE, recursive=TRUE)

#for (chunk in chunks) {
  # Load in data per chunk, each separated by at least 10kb
  peaks <- read.delim(filepath, header=FALSE, as.is=T)
  #colnames(peaks) <- c("seqname", "start", "end", "ID", "score", "density_summit", "sampleID", "wavelet_summit")
  colnames(peaks) <- c("seqname", "start", "end", "sampleID", "score", "density_summit", "wavelet_summit")
  peaks <- peaks[order(peaks$wavelet_summit),] # Order peaks by wavelet summit first

  DHSs <- NULL; # Final list of delineated DHSs

  # Overwrite peak IDs, chunked based on a 20bp+ separation of summits
  loc_ID_seps <- c(0, which(diff(peaks$wavelet_summit) > 20), nrow(peaks))
  loc_IDs <- rep(1:(length(loc_ID_seps)-1), times=diff(loc_ID_seps))
  peaks$ID <- loc_IDs;

  # Iterate over peak clumps
  for (loc_ID in unique(peaks$ID)) {
    # if (loc_ID == "id-1431") stop() # debugging

    # localize data, in particular summit coordinates
    sel <- which(peaks$ID == loc_ID);
    loc_peaks <- peaks[sel,]
    xlim <- range(loc_peaks$wavelet_summit)
    loc_peaks$wavelet_summit_loc <- loc_peaks$wavelet_summit - xlim[1] + 1
    n <- diff(xlim)+1

    # Obtain cut points, to be processed independently
    cuts <- get_cut_points(loc_peaks$wavelet_summit_loc);

    # Select localized data for FWHM determination 
    loc_DHSs <- NULL
    for (i in 1:(length(cuts)-1)) {
      idx <- which(loc_peaks$wavelet_summit_loc > cuts[i] & loc_peaks$wavelet_summit_loc <= cuts[i+1])
      if (length(idx) == 0) next;
      new_DHS <- get_FWHM(loc_peaks[idx,])
      new_DHS$ID <- paste(loc_ID, i, sep="_")
      loc_peaks$ID[idx] <- paste(loc_ID, i, sep="_")
      loc_DHSs <- rbind(loc_DHSs, new_DHS) # Save resultant (candidate) DHS
    }
    
    if (nrow(loc_DHSs) > 1) {
      # Resolve summit overlaps (summit of one element in FWHM of another)
      num <- 1;
      while(num > 0) {
        ovl_rm <- merge_overlap(loc_DHSs, loc_peaks, type="summit")
        loc_DHSs <- ovl_rm$DHSs
        loc_peaks <- ovl_rm$peaks
        num <- ovl_rm$num
      }
    }

    DHSs <- rbind(DHSs, loc_DHSs)
    peaks$ID[sel] <- loc_peaks$ID
  }
  
  # Final pass across all data in this chunk
  if (nrow(DHSs) > 1) {
    # Resolve summit overlaps across originally defined "peak clumps"
    num <- 1;
    while(num > 0) {
      ovl_rm <- merge_overlap(DHSs, peaks, type="summit")
      DHSs <- ovl_rm$DHSs
      peaks <- ovl_rm$peaks
      num <- ovl_rm$num
      print(num)
    }
  }
  DHSs <- DHSs[order(DHSs$seqname, DHSs$start, DHSs$end),]

  peaks$ID <- paste(gsub(".bed", "", chunk), peaks$ID, sep="_")
  DHSs$ID <- paste(gsub(".bed", "", chunk), DHSs$ID, sep="_")

  write.table(DHSs, file=paste("DHSs_all", chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  write.table(peaks, file=paste("peaks_all", chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
#}

######################################################################################################################################

