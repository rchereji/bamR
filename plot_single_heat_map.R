#!/usr/bin/env Rscript
library("optparse")

options = list(
  make_option(c("-f", "--files"), type="character", default=NULL,
              help="Data file names; replicates separated by comma [e.g.: -c data_1.bam,data_2.bam,data_3.bam]"),
  make_option("--sampleLabel", type="character", default="Control",
              help="Sample label [default = %default]"),
  make_option(c("-t", "--type"), type="character", default="occ", 
              help="Type of heat map; multiple types separated by commas only [options: occ, dyads; default = %default]"),
  make_option(c("-r", "--reference"), type="character", default="TSS",
              help="Reference points to align [options: TSS (default), TTS, Plus1]"),
  make_option(c("-p", "--presortedList"), type="character", default=NULL,
              help="Presorted list of genes (csv file)"),
  make_option(c("-l", "--minLength"), type="integer", default=50,
              help="The smallest DNA fragment to be considered [default = %default]"),
  make_option(c("-L", "--maxLength"), type="integer", default=300,
              help="The largest DNA fragment to be considered [default = %default]"),
  make_option(c("-u", "--upstreamPlotWindow"), type="integer", default=1000,
              help="Length of the upstream region that will be plotted [default = %default]"),
  make_option(c("-d", "--downstreamPlotWindow"), type="integer", default=1000,
              help="Length of the downstream region that will be plotted [default = %default]"),
  make_option("--colorScale", type="double", default=2,
              help="Maximum color scale in single condition heat maps [default = %default]")
)

opt_parser = OptionParser(option_list=options)
opt = parse_args(opt_parser)

if (is.null(opt$files)){
  print_help(opt_parser)
  stop("At least the dataset file names must be supplied.", call.=FALSE)
}

###############
# Test parameters
# opt = list(control="AGH0220_1.bam,AGH0220_2.bam,AGH0220_3.bam",
#            treatment="AGH0220_4.bam,AGH0220_5.bam,AGH0220_6.bam",
#            reference="Plus1",
#            presortedList=NULL,
#            minLength=50,
#            maxLength=200,
#            upstreamPlotWindow=1000,
#            downstreamPlotWindow=1000,
#            upstreamSortWindow=500,
#            downstreamSortWindow=100,
#            colorScaleDiff=0.5)
###############

##################
# Initialization #
##################
# Label
label = opt$sampleLabel

# Type of alignments, i.e. reference points
selectedReference = opt$reference

# Size selection parameters: specify the interval of lengths to be analyzed
Lmin = opt$minLength  # the smallest DNA fragment to be considered
Lmax = opt$maxLength  # the largest DNA fragment to be considered

# Window selection parameters
beforeRef = opt$upstreamPlotWindow    # length of the upstream region that will be plotted
afterRef  = opt$downstreamPlotWindow  # length of the downstream region that will be plotted

plot.type = toupper(opt$type)
plot.type = strsplit(plot.type, ',')[[1]]

# Load the necessary R packages
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(caTools)
  library(colorRamps)
  library(Rsamtools)
})


######################################
# Function to align multiple regions #
######################################
# Input arguments:
#   Profile           - The profile that needs to be aligned, e.g. coverage/occupancy profile
#   ReferenceGRanges  - GRanges with the windows centered on the reference points
#
# Output:
#   AlignedProfiles   - a matrix with each row corresponding to an aligned locus

AlignRegions = function(Profile, ReferenceGRanges)
{
  # Create Views with all the ReferenceGRanges
  chrName = unique(as.character(seqnames(ReferenceGRanges)))
  myViews = Views(Profile[chrName], as(ReferenceGRanges, "RangesList")[chrName])
  AlignedProfilesList = lapply(myViews, function(gr) t(viewApply(gr, as.vector)))
  AlignedProfiles = do.call("rbind", AlignedProfilesList)
  
  ## Get the index of ReferenceGRanges, which were reorganized by as(ReferenceGRanges, "RangesList")
  listInd = split(1:length(ReferenceGRanges), as.factor(seqnames(ReferenceGRanges)))
  idx = do.call("c", listInd)
  
  rownames(AlignedProfiles) = idx
  AlignedProfiles = AlignedProfiles[order(idx),]
  
  ## Flip regions from the Crick strand
  CrickInd = which(as.character(strand(ReferenceGRanges)) == "-")
  AlignedProfiles[CrickInd,] = AlignedProfiles[CrickInd, ncol(AlignedProfiles):1]
  return(AlignedProfiles)
}

##############################################
# Function for drawing the image color scale #
##############################################
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}


# Load annotations
sacCer3transcripts = read.csv("sacCer3_annotations.csv", header=TRUE, stringsAsFactors=FALSE)
rownames(sacCer3transcripts) = sacCer3transcripts$ORF

if (! is.null(opt$presortedList)){
  # Get the presorted genes
  presorted_genes = read.csv(opt$presortedList, header = FALSE, stringsAsFactors = FALSE)
  presorted_genes = presorted_genes$V1
  Indices = match(presorted_genes, rownames(sacCer3transcripts))
  sacCer3transcripts = sacCer3transcripts[Indices,]
}

switch(selectedReference, 
       TSS={
         ReferencePos = sacCer3transcripts$TSS
         ReferenceChr = sacCer3transcripts$Chr
         RefStrand = sacCer3transcripts$Strand
         ORF = sacCer3transcripts$ORF
         
         Watson = (RefStrand == 1)
         leftEdge = ReferencePos-beforeRef
         rightEdge = ReferencePos+afterRef
         
         leftEdgeCrick = ReferencePos-afterRef
         rightEdgeCrick = ReferencePos+beforeRef
         
         leftEdge[!Watson] = leftEdgeCrick[!Watson]
         rightEdge[!Watson] = rightEdgeCrick[!Watson]
         
         ReferenceGRanges = GRanges(seqnames=ReferenceChr,
                                    IRanges(start=leftEdge,
                                            end=rightEdge),
                                    strand=RefStrand)
         
         # Construct GRanges for the entire chromosomes...
         chrLen = c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
         names(chrLen) = c('chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI')
         wholeChr = GRanges(seqnames=names(chrLen),
                            IRanges(start=rep(1, length(chrLen)),
                                    end=chrLen))
         # ... and remove the windows that fall outside of the chromosome edges
         Hits = as.matrix(findOverlaps(ReferenceGRanges, wholeChr, type="within", ignore.strand=TRUE))
         idxToKeep = Hits[,"queryHits"]
         ReferenceGRanges = ReferenceGRanges[idxToKeep]
         ORF = ORF[idxToKeep]
       },
       TTS={
         ReferencePos = sacCer3transcripts$TTS
         ReferenceChr = sacCer3transcripts$Chr
         RefStrand = sacCer3transcripts$Strand
         ORF = sacCer3transcripts$ORF
         
         Watson = (RefStrand == 1)
         leftEdge = ReferencePos-beforeRef
         rightEdge = ReferencePos+afterRef
         
         leftEdgeCrick = ReferencePos-afterRef
         rightEdgeCrick = ReferencePos+beforeRef
         
         leftEdge[!Watson] = leftEdgeCrick[!Watson]
         rightEdge[!Watson] = rightEdgeCrick[!Watson]
         
         ReferenceGRanges = GRanges(seqnames=ReferenceChr,
                                    IRanges(start=leftEdge,
                                            end=rightEdge),
                                    strand=RefStrand)
         
         # Construct GRanges for the entire chromosomes...
         chrLen = c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
         names(chrLen) = c('chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI')
         wholeChr = GRanges(seqnames=names(chrLen),
                            IRanges(start=rep(1, length(chrLen)),
                                    end=chrLen))
         # ... and remove the windows that fall outside of the chromosome edges
         Hits = as.matrix(findOverlaps(ReferenceGRanges, wholeChr, type="within", ignore.strand=TRUE))
         idxToKeep = Hits[,"queryHits"]
         ReferenceGRanges = ReferenceGRanges[idxToKeep]
         ORF = ORF[idxToKeep]
       },
       Plus1={
         ReferencePos = sacCer3transcripts$Plus1
         ReferenceChr = sacCer3transcripts$Chr
         RefStrand = sacCer3transcripts$Strand
         ORF = sacCer3transcripts$ORF
         
         validValues = ! is.nan(ReferencePos)
         ReferencePos = ReferencePos[validValues]
         ReferenceChr = ReferenceChr[validValues]
         RefStrand = RefStrand[validValues]
         ORF = ORF[validValues]
         
         Watson = (RefStrand == 1)
         leftEdge = ReferencePos-beforeRef
         rightEdge = ReferencePos+afterRef
         
         leftEdgeCrick = ReferencePos-afterRef
         rightEdgeCrick = ReferencePos+beforeRef
         
         leftEdge[!Watson] = leftEdgeCrick[!Watson]
         rightEdge[!Watson] = rightEdgeCrick[!Watson]
         
         ReferenceGRanges = GRanges(seqnames=ReferenceChr,
                                    IRanges(start=leftEdge,
                                            end=rightEdge),
                                    strand=RefStrand)
         
         # Construct GRanges for the entire chromosomes...
         chrLen = c(230218,813184,316620,1531933,576874,270161,1090940,562643,439888,745751,666816,1078177,924431,784333,1091291,948066)
         names(chrLen) = c('chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII','chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI')
         wholeChr = GRanges(seqnames=names(chrLen),
                            IRanges(start=rep(1, length(chrLen)),
                                    end=chrLen))
         # ... and remove the windows that fall outside of the chromosome edges
         Hits = as.matrix(findOverlaps(ReferenceGRanges, wholeChr, type="within", ignore.strand=TRUE))
         idxToKeep = Hits[,"queryHits"]
         ReferenceGRanges = ReferenceGRanges[idxToKeep]
         ORF = ORF[idxToKeep]
       }
)


#########################################
# Import the paired-end sequencing data #
#########################################
# Read data
allDataFiles = opt$files
allDataFiles = strsplit(allDataFiles, ",")[[1]]
noFiles = length(allDataFiles)
reads = GRanges()

for (f in 1:noFiles){
  # Data file name
  inputFilename = allDataFiles[f]
  all_fields = c("rname", "pos", "isize")
  param = ScanBamParam(what = all_fields, 
                       flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, 
                                          isUnmappedQuery = FALSE, hasUnmappedMate = FALSE, 
                                          isMinusStrand = FALSE, isMateMinusStrand = TRUE,
                                          isNotPassingQualityControls = FALSE))
  bam = scanBam(inputFilename, param=param)
  
  # Keep only the proper reads, with the length > 0
  posStrandReads = (bam[[1]]$isize > 0)
  
  reads0 = GRanges(seqnames=Rle(bam[[1]]$rname[posStrandReads]),
                   ranges = IRanges(start=bam[[1]]$pos[posStrandReads], width=bam[[1]]$isize[posStrandReads]),
                   strand = "*")
  rm(bam)
  reads = c(reads, reads0)
}

# Discard the reads from the rDNA region, chrXII:451000-469000
rDNAregion = GRanges(seqnames = "chrXII",
                     ranges = IRanges(start=451000, end=469000),
                     strand = "*")
rDNAInd = overlapsAny(reads, rDNAregion, ignore.strand=TRUE)
reads = reads[!rDNAInd]

# Discard the reads from the yeast mitochondrial DNA
reads = reads[seqnames(reads) != 'chrM']


##################
# Size selection #
##################
# Eliminate the reads that are shorter than Lmin or longer than Lmax
readLength = width(reads)
goodReadsInd = ((readLength >= Lmin) & (readLength <= Lmax))
reads = reads[goodReadsInd]
readLength = width(reads)
TotalNoReads = length(reads)


if ("OCC" %in% plot.type){
  #################
  # Normalization #
  #################
  # Compute rescaling coefficients, such that after rescaling the average occupancy is 1, for each chromosome 
  Occ = coverage(reads)
  chrLabel = seqlevels(Occ)
  noChr = length(chrLabel)
  coverageWeight = list()
  for(chr in chrLabel)
  {
    coverageWeight[[chr]] = 1/mean(Occ[[chr]])
  }
  
  Occ = coverage(reads, weight=coverageWeight)
  
  
  ############################
  ## Align occupancy profile #
  ############################
  dir.create("Heatmap_Occ", showWarnings = FALSE, recursive = TRUE)
  
  # Align occupancy profile
  AlignedProfile = AlignRegions(Occ, ReferenceGRanges)
  AlignedProfile[AlignedProfile < 0] = 0
  rownames(AlignedProfile) = ORF
  colnames(AlignedProfile) = seq(from=-beforeRef, to=afterRef, by=1)
  
  SmoothedAlignedProfile = runmean(AlignedProfile, 21, alg="C", endrule="mean")
  SmoothedAlignedProfile = t(runmean(t(SmoothedAlignedProfile), 21, alg="C", endrule="mean"))
  SmoothedAlignedProfile[SmoothedAlignedProfile < 0] = 0
  
  sample.name.title = label
  if (! is.null(opt$presortedList)){
    sample.name = paste(label, ".presorted", sep="")
  } else {
    sample.name = label
  }
  
  maxScale = opt$colorScale
  
  # Save the aligned occupancy
  switch(selectedReference, 
         TSS={
           save(AlignedProfile, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, file=paste("Heatmap_Occ/Aligned_Occ_TSS.", sample.name, ".", Lmin, "_", Lmax, ".RData", sep=""))
         },
         TTS={
           save(AlignedProfile, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, file=paste("Heatmap_Occ/Aligned_Occ_TTS.", sample.name, ".", Lmin, "_", Lmax, ".RData", sep=""))
         },
         Plus1={
           save(AlignedProfile, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, file=paste("Heatmap_Occ/Aligned_Occ_Plus1.", sample.name, ".", Lmin, "_", Lmax, ".RData", sep=""))
         }
  )
  
  #####################
  # Plot the heat map #
  #####################
  
  # Plot the heat map
  switch(selectedReference, 
         TSS={
           pdf(paste("Heatmap_Occ/Heatmap_Occ_TSS.", sample.name, ".", Lmin, "_", Lmax, ".pdf", sep=""), w=6, h=8)
           layout(matrix(c(1,2), ncol=2), widths=c(5,1))
           par(mar=c(5,2,2,2))
           image(-beforeRef:afterRef, 1:nrow(SmoothedAlignedProfile), t(SmoothedAlignedProfile), col=matlab.like(101), ylim=c(nrow(SmoothedAlignedProfile)+0.5,0.5), breaks = c(seq(0, maxScale, length.out = 101), max(maxScale + 0.001, max(SmoothedAlignedProfile))), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=FALSE, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           title(main=sample.name.title, xlab="Position relative to TSS (bp)", ylab="", cex.main=1.4, cex.lab=1.4)
           abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 3)
           box()
           
           #Add scale
           par(mar=c(25,0,2,4))
           image.scale(col = matlab.like(100), breaks = seq(0, maxScale, length.out = 101), horiz=FALSE, xlab="", yaxt="n")
           if (maxScale <= 3) {
             axis(4, at=seq(0, maxScale, 0.5), las=2)
           } else {
             axis(4, at=seq(0, maxScale, 1), las=2)
           }
           box()
           garbage = dev.off()
         },
         TTS={
           pdf(paste("Heatmap_Occ/Heatmap_Occ_TTS.", sample.name, ".", Lmin, "_", Lmax, ".pdf", sep=""), w=6, h=8)
           layout(matrix(c(1,2), ncol=2), widths=c(5,1))
           par(mar=c(5,2,2,2))
           image(-beforeRef:afterRef, 1:nrow(SmoothedAlignedProfile), t(SmoothedAlignedProfile), col=matlab.like(101), ylim=c(nrow(SmoothedAlignedProfile)+0.5,0.5), breaks = c(seq(0, maxScale, length.out = 101), max(maxScale + 0.001, max(SmoothedAlignedProfile))), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=FALSE, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           title(main=sample.name.title, xlab="Position relative to TTS (bp)", ylab="", cex.main=1.4, cex.lab=1.4)
           abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 3)
           box()
           
           #Add scale
           par(mar=c(25,0,2,4))
           image.scale(col = matlab.like(100), breaks = seq(0, maxScale, length.out = 101), horiz=FALSE, xlab="", yaxt="n")
           if (maxScale <= 3) {
             axis(4, at=seq(0, maxScale, 0.5), las=2)
           } else {
             axis(4, at=seq(0, maxScale, 1), las=2)
           }
           box()
           garbage = dev.off()
         },
         Plus1={
           pdf(paste("Heatmap_Occ/Heatmap_Occ_Plus1.", sample.name, ".", Lmin, "_", Lmax, ".pdf", sep=""), w=6, h=8)
           layout(matrix(c(1,2), ncol=2), widths=c(5,1))
           par(mar=c(5,2,2,2))
           image(-beforeRef:afterRef, 1:nrow(SmoothedAlignedProfile), t(SmoothedAlignedProfile), col=matlab.like(101), ylim=c(nrow(SmoothedAlignedProfile)+0.5,0.5), breaks = c(seq(0, maxScale, length.out = 101), max(maxScale + 0.001, max(SmoothedAlignedProfile))), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=FALSE, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           title(main=sample.name.title, xlab="Position relative to +1 nuc. (bp)", ylab="", cex.main=1.4, cex.lab=1.4)
           abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 3)
           box()
           
           #Add scale
           par(mar=c(25,0,2,4))
           image.scale(col = matlab.like(100), breaks = seq(0, maxScale, length.out = 101), horiz=FALSE, xlab="", yaxt="n")
           if (maxScale <= 3) {
             axis(4, at=seq(0, maxScale, 0.5), las=2)
           } else {
             axis(4, at=seq(0, maxScale, 1), las=2)
           }
           box()
           garbage = dev.off()
         }
  )
}

if ("DYADS" %in% plot.type){
  #################
  # Normalization #
  #################
  # Compute rescaling coefficients, such that after rescaling the average occupancy is 1, for each chromosome
  reads = resize(reads, width = 1, fix = "center")
  Dyads = coverage(reads)
  chrLabel = seqlevels(Dyads)
  noChr = length(chrLabel)
  coverageWeight = list()
  for(chr in chrLabel)
  {
    coverageWeight[[chr]] = 1/mean(Dyads[[chr]])
  }
  
  Dyads = coverage(reads, weight=coverageWeight)
  
  ########################
  ## Align dyads profile #
  ########################
  dir.create("Heatmap_Dyads", showWarnings = FALSE, recursive = TRUE)
  
  # Align dyads profile
  AlignedProfile = AlignRegions(Dyads, ReferenceGRanges)
  AlignedProfile[AlignedProfile < 0] = 0
  rownames(AlignedProfile) = ORF
  colnames(AlignedProfile) = seq(from=-beforeRef, to=afterRef, by=1)
  
  SmoothedAlignedProfile = runmean(AlignedProfile, 21, alg="C", endrule="mean")
  SmoothedAlignedProfile = t(runmean(t(SmoothedAlignedProfile), 21, alg="C", endrule="mean"))
  SmoothedAlignedProfile[SmoothedAlignedProfile < 0] = 0;
  
  sample.name.title = label
  if (! is.null(opt$presortedList)){
    sample.name = paste(label, ".presorted", sep="")
  } else {
    sample.name = label
  }
  
  maxScale = opt$colorScale
  
  # Save the aligned dyads
  switch(selectedReference, 
         TSS={
           save(AlignedProfile, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, file=paste("Heatmap_Dyads/Aligned_Dyads_TSS.", sample.name, ".", Lmin, "_", Lmax, ".RData", sep=""))
         },
         TTS={
           save(AlignedProfile, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, file=paste("Heatmap_Dyads/Aligned_Dyads_TTS.", sample.name, ".", Lmin, "_", Lmax, ".RData", sep=""))
         },
         Plus1={
           save(AlignedProfile, Lmin, Lmax, beforeRef, afterRef, TotalNoReads, file=paste("Heatmap_Dyads/Aligned_Dyads_Plus1.", sample.name, ".", Lmin, "_", Lmax, ".RData", sep=""))
         }
  )
  
  #####################
  # Plot the heat map #
  #####################
  
  # Plot the heat map
  switch(selectedReference, 
         TSS={
           pdf(paste("Heatmap_Dyads/Heatmap_Dyads_TSS.", sample.name, ".", Lmin, "_", Lmax, ".pdf", sep=""), w=6, h=8)
           layout(matrix(c(1,2), ncol=2), widths=c(5,1))
           par(mar=c(5,2,2,2))
           image(-beforeRef:afterRef, 1:nrow(SmoothedAlignedProfile), t(SmoothedAlignedProfile), col=matlab.like(101), ylim=c(nrow(SmoothedAlignedProfile)+0.5,0.5), breaks = c(seq(0, maxScale, length.out = 101), max(maxScale + 0.001, max(SmoothedAlignedProfile))), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=FALSE, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           title(main=sample.name.title, xlab="Position relative to TSS (bp)", ylab="", cex.main=1.4, cex.lab=1.4)
           abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 3)
           box()
           
           #Add scale
           par(mar=c(25,0,2,4))
           image.scale(col = matlab.like(100), breaks = seq(0, maxScale, length.out = 101), horiz=FALSE, xlab="", yaxt="n")
           if (maxScale <= 3) {
             axis(4, at=seq(0, maxScale, 0.5), las=2)
           } else {
             axis(4, at=seq(0, maxScale, 1), las=2)
           }
           box()
           garbage = dev.off()
         },
         TTS={
           pdf(paste("Heatmap_Dyads/Heatmap_Dyads_TTS.", sample.name, ".", Lmin, "_", Lmax, ".pdf", sep=""), w=6, h=8)
           layout(matrix(c(1,2), ncol=2), widths=c(5,1))
           par(mar=c(5,2,2,2))
           image(-beforeRef:afterRef, 1:nrow(SmoothedAlignedProfile), t(SmoothedAlignedProfile), col=matlab.like(101), ylim=c(nrow(SmoothedAlignedProfile)+0.5,0.5), breaks = c(seq(0, maxScale, length.out = 101), max(maxScale + 0.001, max(SmoothedAlignedProfile))), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=FALSE, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           title(main=sample.name.title, xlab="Position relative to TTS (bp)", ylab="", cex.main=1.4, cex.lab=1.4)
           abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 3)
           box()
           
           #Add scale
           par(mar=c(25,0,2,4))
           image.scale(col = matlab.like(100), breaks = seq(0, maxScale, length.out = 101), horiz=FALSE, xlab="", yaxt="n")
           if (maxScale <= 3) {
             axis(4, at=seq(0, maxScale, 0.5), las=2)
           } else {
             axis(4, at=seq(0, maxScale, 1), las=2)
           }
           box()
           garbage = dev.off()
         },
         Plus1={
           pdf(paste("Heatmap_Dyads/Heatmap_Dyads_Plus1.", sample.name, ".", Lmin, "_", Lmax, ".pdf", sep=""), w=6, h=8)
           layout(matrix(c(1,2), ncol=2), widths=c(5,1))
           par(mar=c(5,2,2,2))
           image(-beforeRef:afterRef, 1:nrow(SmoothedAlignedProfile), t(SmoothedAlignedProfile), col=matlab.like(101), ylim=c(nrow(SmoothedAlignedProfile)+0.5,0.5), breaks = c(seq(0, maxScale, length.out = 101), max(maxScale + 0.001, max(SmoothedAlignedProfile))), axes=FALSE, xlab="", ylab="", useRaster=TRUE)
           par(tcl= -0.2)
           axis(1, at=seq(-beforeRef, afterRef, by=100), labels=FALSE, lwd=1, lwd.ticks=1)
           par(tcl= -0.5)
           axis(1, at=seq(-beforeRef, afterRef, by=500), labels=TRUE, lwd=0, lwd.ticks=1, cex.axis=1.25)
           title(main=sample.name.title, xlab="Position relative to +1 nuc. (bp)", ylab="", cex.main=1.4, cex.lab=1.4)
           abline(v = 0, untf = FALSE, col = "white", lty = "longdash", lwd = 3)
           box()
           
           #Add scale
           par(mar=c(25,0,2,4))
           image.scale(col = matlab.like(100), breaks = seq(0, maxScale, length.out = 101), horiz=FALSE, xlab="", yaxt="n")
           if (maxScale <= 3) {
             axis(4, at=seq(0, maxScale, 0.5), las=2)
           } else {
             axis(4, at=seq(0, maxScale, 1), las=2)
           }
           box()
           garbage = dev.off()
         }
  )
}