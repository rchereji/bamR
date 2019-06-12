#!/usr/bin/env Rscript
library("optparse")

options = list(
  make_option(c("-c", "--controlFiles"), type="character", default=NULL,
              help="Control file names, separated by comma [e.g.: -c control_1.bam,control_2.bam,control_3.bam]"),
  make_option(c("-t", "--treatmentFiles"), type="character", default=NULL,
              help="Treatment file names, separated by comma [e.g.: -c treatment_1.bam,treatment_2.bam,treatment_3.bam]"),
  make_option("--controlLabel", type="character", default="Control",
              help="Control label [default = %default]"),
  make_option("--treatmentLabel", type="character", default="Treatment",
              help="Treatment label [default = %default]"),
  make_option(c("-l", "--minLength"), type="integer", default=50, 
              help="The smallest DNA fragment to be considered [default = %default]"),
  make_option(c("-L", "--maxLength"), type="integer", default=300, 
              help="The largest DNA fragment to be considered [default = %default]"),
  make_option(c("-p", "--presortedList"), type="character", default=NULL,
              help="Presorted list of genes (csv file)"),
  make_option("--colorScaleDiff", type="double", default=1,
              help="Maximum color scale in difference heat map [default = %default]"),
  make_option("--annotationsFile", type="character", default="sacCer3_annotations.csv",
              help="Filename (csv format) containing the annotations [default = %default]")
) 

opt_parser = OptionParser(option_list=options)
opt = parse_args(opt_parser)

if (is.null(opt$controlFiles) | is.null(opt$treatmentFiles)){
  print_help(opt_parser)
  stop("At least the dataset file names must be supplied.", call.=FALSE)
}

##################
# Initialization #
##################
# Labels
controlLabel = opt$controlLabel
treatmentLabel = opt$treatmentLabel

# Size selection parameters: specify the interval of lengths to be analyzed
Lmin = opt$minLength   # the smallest DNA fragment to be considered
Lmax = opt$maxLength  # the largest DNA fragment to be considered

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
# Read control data
allDataFiles = opt$controlFiles
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

Occ_control = coverage(reads, weight=coverageWeight)


# Read treatment data
allDataFiles = opt$treatmentFiles
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

Occ_treatment = coverage(reads, weight=coverageWeight)


#####################################################
# Function to compute averages for multiple regions #
#####################################################
# Input arguments:
#   Profile           - The profile that needs to be aligned, e.g. coverage/occupancy profile
#   ReferenceGRanges  - GRanges with the windows centered on the reference points
#
# Output:
#   AvgOcc            - a vector with the average occupancy for each gene body

Average_Occ_over_promoter = function(Profile, ReferenceGRanges)
{
  # Create Views with all the ReferenceGRanges
  chrName = unique(as.character(seqnames(ReferenceGRanges)))
  myViews = Views(Profile[chrName], as(ReferenceGRanges, "IntegerRangesList")[chrName])
  AvgProfileList = viewMeans(myViews)
  AvgProfile = unlist(AvgProfileList)
  
  ## Get the index of ReferenceGRanges, which were reorganized by as(ReferenceGRanges, "IntegerRangesList")
  listInd = split(1:length(ReferenceGRanges), as.factor(seqnames(ReferenceGRanges)))
  idx = do.call("c", listInd)
  
  names(AvgProfile) = idx
  AvgProfile = AvgProfile[order(idx)]
  names(AvgProfile) = ReferenceGRanges$ORF
  return(AvgProfile)
}


# Load annotations
sacCer3transcripts = read.csv(opt$annotationsFile, header=TRUE, stringsAsFactors=FALSE)
rownames(sacCer3transcripts) = sacCer3transcripts$ORF

if (! is.null(opt$presortedList)){
  # Get the presorted genes
  presorted_genes = read.csv(opt$presortedList, header = FALSE, stringsAsFactors = FALSE)
  presorted_genes = presorted_genes$V1
  Indices = match(presorted_genes, rownames(sacCer3transcripts))
  sacCer3transcripts = sacCer3transcripts[Indices,]
}

Minus1 = sacCer3transcripts$Minus1
Plus1 = sacCer3transcripts$Plus1
ReferenceChr = sacCer3transcripts$Chr
RefStrand = sacCer3transcripts$Strand
ORF = sacCer3transcripts$ORF

Watson = (RefStrand == 1)
leftEdge = Minus1 - 73
rightEdge = Plus1 + 73

leftEdgeCrick = Plus1 - 73
rightEdgeCrick = Minus1 + 73

leftEdge[!Watson] = leftEdgeCrick[!Watson]
rightEdge[!Watson] = rightEdgeCrick[!Watson]

# Remove genes without +1/-1 annotations
available_annot = !is.nan(leftEdge)
ReferenceChr = ReferenceChr[available_annot]
leftEdge = leftEdge[available_annot]
rightEdge = rightEdge[available_annot]
RefStrand = RefStrand[available_annot]
ORF = ORF[available_annot]

ReferenceGRanges = GRanges(seqnames=ReferenceChr,
                           IRanges(start=leftEdge,
                                   end=rightEdge),
                           strand=RefStrand,
                           ORF=ORF)

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


##############################
## Compute average occupancy #
##############################
dir.create("Avg_Occ_over_promoters", showWarnings = FALSE, recursive = TRUE)

# Compute average occupancy
Avg_Occ_control = Average_Occ_over_promoter(Occ_control, ReferenceGRanges)
Avg_Occ_control[Avg_Occ_control < 0] = 0 # eliminate rounding errors for very small numbers

Avg_Occ_treatment = Average_Occ_over_promoter(Occ_treatment, ReferenceGRanges)
Avg_Occ_treatment[Avg_Occ_treatment < 0] = 0 # eliminate rounding errors for very small numbers

Avg_Occ_diff = Avg_Occ_treatment - Avg_Occ_control
Avg_Occ = Avg_Occ_diff


sample.name.title = paste(treatmentLabel, "-", controlLabel, sep="")
if (! is.null(opt$presortedList)){
  sample.name = paste(treatmentLabel, "-", controlLabel, ".presorted", sep="")
} else {
  sample.name = paste(treatmentLabel, "-", controlLabel, sep="")
}


# Save the average occupancy
save(Avg_Occ_diff, Lmin, Lmax, file=paste("Avg_Occ_over_promoters/Avg_Occ_diff_over_promoters.", sample.name, ".", Lmin, "_", Lmax, ".RData", sep=""))
write.csv(data.frame(Gene=names(Avg_Occ_diff), Average_Occ_diff=Avg_Occ_diff), 
          file=paste("Avg_Occ_over_promoters/Avg_Occ_diff_over_promoters.", sample.name, ".", Lmin, "_", Lmax, ".csv", sep=""), row.names=FALSE)


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


Smoothed_Avg_Occ = cbind(runmean(Avg_Occ, 21, alg="C", endrule="mean"), runmean(Avg_Occ, 21, alg="C", endrule="mean"))
maxScaleDiff = opt$colorScaleDiff
pdf(paste("Avg_Occ_over_promoters/Avg_Occ_diff_over_promoters.", sample.name, ".", Lmin, "_", Lmax, ".pdf", sep=""), w=3, h=8)
layout(matrix(c(1,2), ncol=2), widths=c(3,2))
par(mar=c(5,5,2,2))
image(c(1,2), 1:nrow(Smoothed_Avg_Occ), t(Smoothed_Avg_Occ), col=matlab.like(102), 
      ylim=c(nrow(Smoothed_Avg_Occ)+0.5,0.5), 
      breaks = c(min(-0.001-maxScaleDiff, min(Smoothed_Avg_Occ)), seq(-maxScaleDiff, maxScaleDiff, length.out = 101), max(0.001+maxScaleDiff,max(Smoothed_Avg_Occ))), 
      axes=FALSE, xlab="", ylab="Average occupancy difference over promoters", useRaster=TRUE, cex.lab=1.4)
title(main=sample.name.title, cex.main=1.4)
box()

#Add scale
par(mar=c(25,0,2,4))
image.scale(col = matlab.like(100), breaks = seq(-maxScaleDiff, maxScaleDiff, length.out = 101), horiz=FALSE, xlab="", yaxt="n")
axis(4, at=seq(-maxScaleDiff, maxScaleDiff, length.out = 5), las=2)
box()
garbage = dev.off()
