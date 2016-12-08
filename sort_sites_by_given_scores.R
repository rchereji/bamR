#!/usr/bin/env Rscript
library("optparse")

options = list(
  make_option(c("-a", "--annotFilename"), type="character", default=NULL, 
              help="Annotations of the sites to be aligned (csv format)"),
  make_option(c("-s", "--sortByFilename"), type="character", default=NULL, 
              help="Scores for the sites to be aligned (csv format)"),
  make_option(c("-o", "--order"), type="character", default="decreasing", 
              help="sorting order [options: increasing, decreasing; default = %default]")
) 

opt_parser = OptionParser(option_list=options)
opt = parse_args(opt_parser)

if (is.null(opt$annotFilename)){
  print_help(opt_parser)
  stop("Annotations file not provided!", call.=FALSE)
}

if (is.null(opt$sortByFilename)){
  print_help(opt_parser)
  stop("Scores file not provided!", call.=FALSE)
}

annotFilename = opt$annotFilename
sortByFilename = opt$sortByFilename
sortOrder = opt$order

# annotFilename = 'All_Gcn4_peaks.csv'
# sortByFilename = 'Avg_Occ_per_site.Gcn4_sites.Gcn4.50_300.csv'
# sortOrder = 'decreasing'

data = read.csv(sortByFilename, stringsAsFactors = FALSE)
if (sortOrder == 'increasing'){
  data = data[order(data$Average, decreasing = FALSE),]
} else {
  data = data[order(data$Average, decreasing = TRUE),]
}

sortedSites = data$siteID

annot = read.csv(annotFilename, stringsAsFactors = FALSE)
annot = annot[match(sortedSites, annot$siteID), ]

library(tools)
write.csv(annot, file=paste(file_path_sans_ext(annotFilename), sortOrder, file_ext(annotFilename), sep="."), row.names = FALSE)
