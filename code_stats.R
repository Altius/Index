##########################################################################################################
#                                                                                                        #
# Method for generating a master list / Index of DNaseI hypersensitivity sites.                          #
# All code, implementation details and design by Wouter Meuleman and Eric Rynes.                         #
#                                                                                                        #
# Version: WM20180313                                                                                    #
#                                                                                                        #
##########################################################################################################

source("../../code/general.R")

#annot <- read.delim("master_list_stats_WM20180130.txt", as.is=TRUE);
annot <- read.delim("masterlist_DHSs_WM20180302_overlap_annot.txt", as.is=T, header=FALSE)
colnames(annot) <- c("seqname", "start", "end", "sumND", "numsamples", "width", "summit", "nonoverlapping", "ID", "nonoverlapping")
stats <- read.delim("stats_per_sumND.txt", as.is=TRUE) # Obtained from Eric Rynes

for (i in c(0,1,2,3,4,5,10,15,20)) { 
  print(summary(annot$DHS_width[annot$total_signal > i])) 
}

plotfile("sumND_vs_numelem", type="pdf")
par(mar=c(5,7,3,1))
plot(stats$sumND, stats$numelem/1e6, log="y", type="b", pch=16, cex=2, lwd=3,
     xlab="signal", ylab="", cex.lab=2, cex.axis=2, yaxt="n")
axis(2, las=2, cex.axis=2)
mtext("# DHSs (millions)", side=2, cex=2, line=4)
dev.off()

plotfile("sumND_vs_coverage_perc", type="pdf")
par(mar=c(5,7,3,1))
plot(stats$sumND, stats$coverage_perc, type="b", pch=16, cex=2, lwd=3,
     xlab="signal", ylab="", cex.lab=2, cex.axis=2, yaxt="n")
axis(2, las=2, cex.axis=2)
mtext("genome coverage (%)", side=2, cex=2, line=4)
dev.off()

plotfile("DHS_widths", type="pdf")
par(mar=c(5,7,3,1))
plot(density(annot$DHS_width), xlab="DHS widths", ylab="relative frequency", 
     lwd=3, cex.lab=2, cex.axis=2, main="")
abline(v=median(annot$DHS_width), lwd=3, lty=2, col="grey")
dev.off()

plotfile("DHS_widths_jitter", type="pdf")
par(mar=c(5,7,3,1))
plot(density(annot$DHS_width+sample(-10:10, nrow(annot), replace=TRUE)), main="",
     xlab="DHS widths", ylab="relative frequency", lwd=3, cex.lab=2, cex.axis=2)
abline(v=median(annot$DHS_width), lwd=3, lty=2, col="grey")
dev.off()


