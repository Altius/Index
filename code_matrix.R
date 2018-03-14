##########################################################################################################
#                                                                                                        #
# Method for generating a master list / Index of DNaseI hypersensitivity sites.                          #
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.                         #
#                                                                                                        #
# Version: WM20180313                                                                                    #
#                                                                                                        #
##########################################################################################################

source("../../code/general.R")
library(caTools)
library(Matrix)
source("code_ML.R")

######################################################################################################################################

args=(commandArgs(TRUE))
if (length(args)==0) {
  stop("No arguments supplied.")
} else {
  eval(parse(text=args[[1]])) # parse first argument: chunknum
  eval(parse(text=args[[2]])) # parse second argument: type
}

chunk <- paste("chunk", sprintf("%04d", chunknum), ".bed", sep="");
print(chunk)

options(scipen=20)
#options(warn=0)
options(warn=2)

dir.create(paste("matrix", type, sep="_"), showWarnings=FALSE, recursive=TRUE)
dir.create(paste("matrix_bin", type, sep="_"), showWarnings=FALSE, recursive=TRUE)

######################################################################################################################################

sample_info <- read.table("/net/seq/data/projects/ENCODE3_publications/erynes/ENCODEpaper2017/SubsetOf665/nameMapping665_LN_DS_sample.txt", as.is=TRUE, header=FALSE)
library_nums <- sample_info[,1];
sample_nams <- paste(sample_info[,3], sample_info[,2], sep=".")

#chunks <- dir("peaks", pattern="^chunk.*.bed");
#for (chunk in chunks) {
  # Load in data per chunk, each separated by at least 10kb
  peaks <- read.delim(paste(paste("peaks", type, sep="_"), chunk, sep="/"), header=FALSE, as.is=T)
  colnames(peaks) <- c("seqname", "start", "end", "sampleID", "score", "density_summit", "wavelet_summit", "ID")
  peaks$sampleID <- factor(peaks$sampleID, levels=library_nums)

  DHSs <- read.delim(paste(paste("DHSs", type, sep="_"), chunk, sep="/"), header=FALSE, as.is=T)
  colnames(DHSs) <- c("seqname", "start", "end", "ID", "score", "numsamples", "numpeaks", "width", "summit", "disp_low", "disp_high");

  ### Construct the signal matrix
  dat <- matrix(0, nrow=nrow(DHSs), ncol=length(sample_nams), dimnames=list(DHSs$ID, library_nums))
  # Iterate over DHSs
  for (DHS_ID in DHSs$ID) {
    sel <- which(peaks$ID == DHS_ID);
    dat[DHS_ID,] <- tapply(peaks$score[sel], peaks$sampleID[sel], max)
  }
  dat[is.na(dat)] <- 0
  write.table(dat, file=paste(paste("matrix", type, sep="_"), chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

  # Now the binary matrix
  write.table((dat > 0)+0, file=paste(paste("matrix_bin", type, sep="_"), chunk, sep="/"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  #indexes <- which(dat > 1, arr.ind=TRUE)
  #coords_paste <- paste(DHSs$seqname, ":", DHSs$start, "-", DHSs$end, sep="")
  #dat_bin <- sparseMatrix(i=indexes[,1], j=indexes[,2], dims=dim(dat), dimnames=list(coords_paste, sample_nams));
  #save(dat_bin, file=paste("matrix_bin/", gsub(".bed", "", chunk), ".RData", sep=""));

#}

######################################################################################################################################

