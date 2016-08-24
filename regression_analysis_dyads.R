#!/usr/bin/env Rscript
library("optparse")

options = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Dyads alignment file [e.g. -f Avg_Dyads_Plus1.ExpLabel.RData]")
) 

opt_parser = OptionParser(option_list=options)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("The alignment file name has not been supplied!", call.=FALSE)
}


##################
# Initialization #
##################
# Results files
res.file = opt$file

# Label
label = paste(strsplit(res.file, ".", fixed = TRUE)[[1]][2],
              strsplit(res.file, ".", fixed = TRUE)[[1]][3],
              sep =".")
sample.label = strsplit(res.file, ".", fixed = TRUE)[[1]][2]
lengths.label = gsub("_", "-", strsplit(res.file, ".", fixed = TRUE)[[1]][3])

# Type of alignments, i.e. reference points
selectedReference = strsplit(strsplit(res.file, ".", fixed = TRUE)[[1]][1], "_")[[1]][3]

# Load the necessary R packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(caTools)
  library(colorRamps)
})


load(res.file)
Smoothed_Avg_Dyads = runmean(Avg_Dyads, 21, alg="C", endrule="mean")


########################
## Regression analysis #
########################
Position = -beforeRef:afterRef
Signal = as.numeric(Smoothed_Avg_Dyads)
suppressPackageStartupMessages(library(pracma))
x = findpeaks(Signal, sortstr=TRUE, minpeakdistance=147)

pdf(paste("Regression_", selectedReference, ".", label, ".pdf", sep=""), w=12, h=5)
par(mfrow=c(1,2), mar=c(5,5,1.5,2), oma = c(0, 0, 2, 0)) 
plot(Position, Signal, axes=FALSE, typ='l', ann=FALSE, col = "blue",
     xlim=c(-beforeRef, afterRef), xaxs="i", yaxs="i", ylim = c(0, 1.1*max(Signal)),
     panel.first = c(abline(v = seq(-beforeRef+500, afterRef-500, 500), col = "lightgray", lty = "dashed")))
# x axis
par(tcl= -0.2)
axis(1, at=seq(-beforeRef, afterRef, by=100), labels=FALSE, lwd=1, lwd.ticks=1)
par(tcl= -0.5)
axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
par(tcl= -0.5)
axis(2, at=seq(0,5,0.2), cex.axis=1.25)
switch(selectedReference, 
       TSS={
         title(main="", xlab="Position relative to TSS (bp)", ylab="Average dyad density", cex.main=1.4, cex.lab=1.4)
       },
       Plus1={
         title(main="", xlab="Position relative to +1 nuc. (bp)", ylab="Average dyad density", cex.main=1.4, cex.lab=1.4)
       })
downstreamPeaks = (Position[x[, 2]] >= -50)
upstreamPeaks = (Position[x[, 2]] < -50)
points(Position[x[, 2]][downstreamPeaks], Signal[x[, 2]][downstreamPeaks], pch=20, cex = 2, col="blue")

nucRank = 1:5
peakLoc = sort(Position[x[, 2]][downstreamPeaks])
peakLoc = peakLoc[1:5]
plot(nucRank, peakLoc, axes=FALSE, ann=FALSE, pch=20, cex = 2, col = "blue",
     xaxs="i", yaxs="i", xlim=c(0.5, 5.5), ylim = c(-100, 1.05*max(peakLoc)),
     panel.first = c(abline(h = seq(0, 1000, 200), col = "lightgray", lty = "dashed")))
fit = lm(peakLoc ~ nucRank)
abline(fit, col="red")
par(tcl= -0.5)
axis(1, at=seq(0,6,1), xlim=c(0.5, 5.5), labels=c("0","+1","+2","+3","+4","+5","+6"), cex.axis=1.25)
par(tcl= -0.2)
axis(2, at=seq(-1000, 1000, by=20), labels=F, lwd=1, lwd.ticks=1)
par(tcl= -0.5)
axis(2, at=seq(-1000, 1000, by=100), lwd=0, lwd.ticks=1, cex.axis=1.25)
switch(selectedReference, 
       TSS={
         title(main="", ylab="Position relative to TSS (bp)", xlab="Nucleosome", cex.main=1.4, cex.lab=1.4)
       },
       Plus1={
         title(main="", ylab="Position relative to +1 nuc. (bp)", xlab="Nucleosome", cex.main=1.4, cex.lab=1.4)
       })

mtext(paste("Sample: ", sample.label, "; Fragment lengths: ", lengths.label, " bp", sep=""), outer = TRUE, cex = 1.5)

cf = round(coef(fit), 2)
eq = paste0("y = ", cf[2], " x ", ifelse(sign(cf[1])==1, " + ", " - "), abs(cf[1]))
text(5.5, 100, labels = eq, pos=2, cex=1.25, col="red")
garbage = dev.off()
## End of the regression analysis 