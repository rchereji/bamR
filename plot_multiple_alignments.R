#!/usr/bin/env Rscript
library("optparse")

options = list(
  make_option(c("-f", "--files"), type="character", default=NULL, 
              help="Alignment file names, separated by commas only [e.g. -f Avg_Occ_1.RData,Avg_Occ_2.Rdata]"),
  make_option(c("-l", "--labels"), type="character", default=NULL, 
              help="Labels for the plot, separated by commas only [e.g. -l Exp_1,Exp_2]"),
  make_option(c("-t", "--type"), type="character", default="occ", 
              help="Type of average plot to be generated [options: occ, dyads; default = %default]"),
  make_option(c("-r", "--reference"), type="character", default="TSS", 
              help="Reference points to align [options: TSS (default), TTS, Plus1]")
) 

opt_parser = OptionParser(option_list=options)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least the alignment file names must be supplied.", call.=FALSE)
}

##################
# Initialization #
##################
# Results files
res.files = opt$files
res.files = strsplit(res.files, ',')[[1]]
noFiles = length(res.files)

# Labels
if (! is.null(opt$labels)){
  line.labels = strsplit(opt$labels, ',')[[1]]
} else {
  line.labels = paste("Data", seq(1,noFiles), sep="_")
}

# Type of alignments, i.e. reference points
selectedReference = opt$reference

plot.type = toupper(opt$type)

# Load the necessary R packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(caTools)
  library(colorRamps)
})


if ("OCC" %in% plot.type){
  
  Profile = list()
  for (f in 1:noFiles){
    load(res.files[f])
    Profile[[f]] = Avg_Occ
    }
  
  
  # Plot the figure
  switch(selectedReference, 
         TSS={
           pdf(paste("Avg_Occ_TSS.", paste(line.labels, collapse = "__"), ".", Lmin, "_", Lmax, ".pdf", sep=""), w=5, h=4)
           
           # 1D Occ
           par(mar=c(5,5,2,2))
           
           maxValue = max(c(Profile, recursive=TRUE))
           plot(0, 0, axes=FALSE, type="n", ann=FALSE, col = "blue",
                xlim=c(-beforeRef, afterRef), xaxs="i", yaxs="i", ylim = c(0, 1.1*maxValue),
                panel.first = c(abline(v = seq(-beforeRef+500, afterRef-500, 500), col = "lightgray", lty = "dashed")))
           
           all.colors = rainbow(noFiles)
           
           for (i in 1:noFiles){
             lines(-beforeRef:afterRef, Profile[[i]], lty = 1, col = all.colors[i], type = 'l')
           }
           
           box()
           # x axis
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=F, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           # y axis
           axis(2, at=seq(0,5,0.2), cex.axis=1.25)
           title(main=paste(line.labels, collapse = "__"), font.main = 1, xlab="Position relative to TSS (bp)", ylab="Average occupancy", cex.lab=1.4)
           legend("bottomright", inset=c(0.01, 0.01), legend=line.labels, col=all.colors, lty=1, cex=1, bg = "white")
           garbage = dev.off()
         },
         TTS={
           pdf(paste("Avg_Occ_TTS.", paste(line.labels, collapse = "__"), ".", Lmin, "_", Lmax, ".pdf", sep=""), w=5, h=4)
           
           # 1D Occ
           par(mar=c(5,5,2,2))
           
           maxValue = max(c(Profile, recursive=TRUE))
           plot(0, 0, axes=FALSE, type="n", ann=FALSE, col = "blue",
                xlim=c(-beforeRef, afterRef), xaxs="i", yaxs="i", ylim = c(0, 1.1*maxValue),
                panel.first = c(abline(v = seq(-beforeRef+500, afterRef-500, 500), col = "lightgray", lty = "dashed")))
           
           all.colors = rainbow(noFiles)
           
           for (i in 1:noFiles){
             lines(-beforeRef:afterRef, Profile[[i]], lty = 1, col = all.colors[i], type = 'l')
           }
           
           box()
           # x axis
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=F, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           # y axis
           axis(2, at=seq(0,5,0.2), cex.axis=1.25)
           title(main=paste(line.labels, collapse = "__"), font.main = 1, xlab="Position relative to TTS (bp)", ylab="Average occupancy", cex.lab=1.4)
           legend("bottomright", inset=c(0.01, 0.01), legend=line.labels, col=all.colors, lty=1, cex=1, bg = "white")
           garbage = dev.off()
         },
         Plus1={
           pdf(paste("Avg_Occ_Plus1.", paste(line.labels, collapse = "__"), ".", Lmin, "_", Lmax, ".pdf", sep=""), w=5, h=4)
           
           # 1D Occ
           par(mar=c(5,5,2,2))
           
           maxValue = max(c(Profile, recursive=TRUE))
           plot(0, 0, axes=FALSE, type="n", ann=FALSE, col = "blue",
                xlim=c(-beforeRef, afterRef), xaxs="i", yaxs="i", ylim = c(0, 1.1*maxValue),
                panel.first = c(abline(v = seq(-beforeRef+500, afterRef-500, 500), col = "lightgray", lty = "dashed")))
           
           all.colors = rainbow(noFiles)
           
           for (i in 1:noFiles){
             lines(-beforeRef:afterRef, Profile[[i]], lty = 1, col = all.colors[i], type = 'l')
           }
           
           box()
           # x axis
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=F, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           # y axis
           axis(2, at=seq(0,5,0.2), cex.axis=1.25)
           title(main=paste(line.labels, collapse = "__"), font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Average occupancy", cex.lab=1.4)
           legend("bottomright", inset=c(0.01, 0.01), legend=line.labels, col=all.colors, lty=1, cex=1, bg = "white")
           garbage = dev.off()
         }
  )
}

if ("DYADS" %in% plot.type){
  
  Profile = list()
  for (f in 1:noFiles){
    load(res.files[f])
    Profile[[f]] = runmean(Avg_Dyads, 21, alg="C", endrule="mean")
  }
  
  # Plot the figure
  switch(selectedReference, 
         TSS={
           pdf(paste("Avg_Dyads_TSS.", paste(line.labels, collapse = "__"), ".", Lmin, "_", Lmax, ".pdf", sep=""), w=5, h=4)
           
           # 1D Occ
           par(mar=c(5,5,2,2))
           maxValue = max(c(Profile, recursive=TRUE))
           plot(0, 0, axes=FALSE, type="n", ann=FALSE, col = "blue",
                xlim=c(-beforeRef, afterRef), xaxs="i", yaxs="i", ylim = c(0, 1.1*maxValue),
                panel.first = c(abline(v = seq(-beforeRef+500, afterRef-500, 500), col = "lightgray", lty = "dashed")))
           
           all.colors = rainbow(noFiles)
           
           for (i in 1:noFiles){
             lines(-beforeRef:afterRef, Profile[[i]], lty = 1, col = all.colors[i], type = 'l')
           }
           
           box()
           # x axis
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=F, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           # y axis
           axis(2, at=seq(0,5,0.2), cex.axis=1.25)
           title(main=paste(line.labels, collapse = "__"), font.main = 1, xlab="Position relative to TSS (bp)", ylab="Average dyads density", cex.lab=1.4)
           legend("topleft", inset=c(0.01, 0.01), legend=line.labels, col=all.colors, lty=1, cex=1, bg = "white")
           garbage = dev.off()
         },
         TTS={
           pdf(paste("Avg_Dyads_TTS.", paste(line.labels, collapse = "__"), ".", Lmin, "_", Lmax, ".pdf", sep=""), w=5, h=4)
           
           # 1D Occ
           par(mar=c(5,5,2,2))
           
           maxValue = max(c(Profile, recursive=TRUE))
           plot(0, 0, axes=FALSE, type="n", ann=FALSE, col = "blue",
                xlim=c(-beforeRef, afterRef), xaxs="i", yaxs="i", ylim = c(0, 1.1*maxValue),
                panel.first = c(abline(v = seq(-beforeRef+500, afterRef-500, 500), col = "lightgray", lty = "dashed")))
           
           all.colors = rainbow(noFiles)
           
           for (i in 1:noFiles){
             lines(-beforeRef:afterRef, Profile[[i]], lty = 1, col = all.colors[i], type = 'l')
           }
           
           box()
           # x axis
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=F, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           # y axis
           axis(2, at=seq(0,5,0.2), cex.axis=1.25)
           title(main=paste(line.labels, collapse = "__"), font.main = 1, xlab="Position relative to TTS (bp)", ylab="Average dyads density", cex.lab=1.4)
           legend("topleft", inset=c(0.01, 0.01), legend=line.labels, col=all.colors, lty=1, cex=1, bg = "white")
           garbage = dev.off()
         },
         Plus1={
           pdf(paste("Avg_Dyads_Plus1.", paste(line.labels, collapse = "__"), ".", Lmin, "_", Lmax, ".pdf", sep=""), w=5, h=4)
           
           # 1D Occ
           par(mar=c(5,5,2,2))
           
           maxValue = max(c(Profile, recursive=TRUE))
           plot(0, 0, axes=FALSE, type="n", ann=FALSE, col = "blue",
                xlim=c(-beforeRef, afterRef), xaxs="i", yaxs="i", ylim = c(0, 1.1*maxValue),
                panel.first = c(abline(v = seq(-beforeRef+500, afterRef-500, 500), col = "lightgray", lty = "dashed")))
           
           all.colors = rainbow(noFiles)
           
           for (i in 1:noFiles){
             lines(-beforeRef:afterRef, Profile[[i]], lty = 1, col = all.colors[i], type = 'l')
           }
           
           box()
           # x axis
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=F, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           # y axis
           axis(2, at=seq(0,5,0.2), cex.axis=1.25)
           title(main=paste(line.labels, collapse = "__"), font.main = 1, xlab="Position relative to +1 nuc. (bp)", ylab="Average dyads density", cex.lab=1.4)
           legend("topleft", inset=c(0.01, 0.01), legend=line.labels, col=all.colors, lty=1, cex=1, bg = "white")
           garbage = dev.off()
         }
  )
}