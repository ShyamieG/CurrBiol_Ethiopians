# AUTHOR: SSG
# This script filters distributions of RoH and IBD segments that overlap certain regions that are known to produce false positives
# It tests for a sufficient level of overlap between a given IBD segment and a region that must be excluded
# The required arguments are: Input file - the full name of the file (.match or .ibd file)
#                           : A table of regions to exclude, one per line, the first three columns being the chromosome, start, and end position of each region (ex. output of the get_IBD_outlier_regions.R script)
# The optional arguments are: Percent overlap - the proportion of overlap necessary between a segment and an excluded region to flag it as a possible spurious IBD segment (default is 0.95)
# USAGE:
# Rscript filter_RoH-IBD_segs.R CHABU_RoH_results_25_2.match exclude_regions_hg19_forIBDNe.txt

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Defaults
perc.overlap <- 0.95

# Error handling and updating default values
if (length(args)<2) {
  stop("At least two arguments must be supplied", call.=FALSE)
} else if (length(args)==3) {
  perc.overlap <- as.numeric(args[3])
} else if (length(args)>3) {
  stop("Too many arguments", call.=FALSE)
}

options(stringsAsFactors=F)
fileNameSplit <- strsplit(args[1], split=".", fixed=T)[[1]]
ext.ind <- which(fileNameSplit %in% c("match", "hom", "ibd"))
stem <- paste(fileNameSplit[-ext.ind], collapse=".")
ext <- fileNameSplit[ext.ind]

if (ext=="hom") {
  seg.file <- read.table(args[1], header=T)
  colnames(seg.file)[7:8] <- c("start", "end")
} else if (ext=="match") {
  seg.file <- read.table(args[1])
  if (ncol(seg.file)==13) {
    colnames(seg.file)[3:5] <- c("CHR", "start", "end")
  } else if (ncol(seg.file)==15) {
    colnames(seg.file)[5:7] <- c("CHR", "start", "end")
  }
} else if (ext=="ibd") {
  seg.file <- read.table(args[1])
  colnames(seg.file)[5:7] <- c("CHR", "start", "end")
}

excl.regs <- read.table(args[2])
  colnames(excl.regs) <- c("CHR", "start", "end", "type")

suppressMessages(library("GenomicRanges"))
Excl.Regs <- makeGRangesFromDataFrame(excl.regs[,c("CHR", "start", "end")], seqnames.field="CHR")
Seg.File <- makeGRangesFromDataFrame(seg.file[ ,c("CHR", "start", "end")], seqnames.field="CHR")

Hits <- findOverlaps(Seg.File, Excl.Regs) # find all hits
Overlaps <- pintersect(Seg.File[queryHits(Hits)], Excl.Regs[subjectHits(Hits)])
percentOverlap <- width(Overlaps)/width(Seg.File[queryHits(Hits)])
HighOverlapHits <- Hits[percentOverlap > perc.overlap] # Which have high overlap? Some will be removed

new.seg.file <- seg.file[-queryHits(HighOverlapHits),]

if (ext=="hom") {
  COLNAMES <- c(colnames(new.seg.file)[1:6], "POS1", "POS2", colnames(new.seg.file)[9:13])
  write.table(new.seg.file, "", row.names=F, col.names=COLNAMES, quote=F, sep="\t")
} else {
  write.table(new.seg.file, "", row.names=F, col.names=F, quote=F, sep="\t")
}

