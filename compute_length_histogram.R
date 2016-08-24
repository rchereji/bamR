#!/usr/bin/env Rscript
library("optparse")

options = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Dataset file name (BAM format)"),
  make_option(c("-o", "--outputs"), type="character", default="pdf,csv,RData", 
              help="Types of outputs to be generated, separated by commas only [options: pdf, csv, RData; default = %default]")
) 

opt_parser = OptionParser(option_list=options)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least the dataset file name must be supplied.", call.=FALSE)
}

##################
# Initialization #
##################
# Load the necessary R packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(caTools)
  library(colorRamps)
  library(Rsamtools)
})

#########################################
# Import the paired-end sequencing data #
#########################################
# Data file name
inputFilename = opt$file
sample.name = sub(".bam", "", inputFilename)
all_fields = c("rname", "pos", "isize")
param = ScanBamParam(what = all_fields, 
                     flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, 
                                        isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                                        isMinusStrand = FALSE, isMateMinusStrand = TRUE,
                                        isNotPassingQualityControls = FALSE))
bam = scanBam(inputFilename, param=param)

# Keep only the proper reads, with the length > 0
posStrandReads = (bam[[1]]$isize > 0)

reads = GRanges(seqnames=Rle(bam[[1]]$rname[posStrandReads]),
                ranges = IRanges(start=bam[[1]]$pos[posStrandReads], width=bam[[1]]$isize[posStrandReads]),
                strand = "*")
rm(bam)
readLength = width(reads)
TotalNoReads = length(reads)

#########################################
# Compute the fragment length histogram #
#########################################
output.files = toupper(opt$outputs)
filesToGenerate = strsplit(output.files, ',')[[1]]

# Compute the histogram
h = hist(readLength, breaks=seq(from = 0.5, to = 1000.5, by = 1), plot=FALSE)

# Create folder
dir.create("Length_histograms", showWarnings = FALSE)

if ("PDF" %in% filesToGenerate){
  # Plot the histogram
  pdf(paste("Length_histograms/Length_histogram.", sample.name, ".pdf", sep=""), w=5, h=4)
  plot(h$mids, 100*h$density, axes=FALSE, typ='l', ann=FALSE, col = "blue",
       xlim=c(0,500), xaxs="i", yaxs="i",
       panel.first = c(abline(v = seq(100, 400, 100), col = "lightgray", lty = "dashed")))
  par(tcl= -0.2)
  axis(1, at=seq(0, 500, by=10), labels=F, lwd=1, lwd.ticks=1)
  par(tcl= -0.5)
  axis(1, at=seq(0, 500, by=100), lwd=0, lwd.ticks=1, cex.axis=1.25)
  par(tcl= -0.5)
  axis(2, cex.axis=1.25)
  title(main=sample.name, xlab="Fragment length (bp)", ylab="Percentage (%)", cex.main=1.4, cex.lab=1.4)
  garbage = dev.off()
}

if ("CSV" %in% filesToGenerate){
  # Save the histogram in a CSV format
  write.csv(data.frame(Length=1:1000, Percentage=100*h$density), 
            file=paste("Length_histograms/Length_histogram.", sample.name, ".csv", sep=""), 
            row.names=FALSE)
}

if ("RDATA" %in% filesToGenerate){
  fragmentLength = 1:1000
  Percentage = 100*h$density
  save(fragmentLength, Percentage, TotalNoReads, file=paste("Length_histograms/Length_histogram.", sample.name, ".RData", sep=""))
}
