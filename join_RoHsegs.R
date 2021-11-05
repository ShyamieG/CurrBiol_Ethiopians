# AUTHOR: SSG
# This script joins runs of homozygosity that span stretches of low SNP density
# It requires a .bim file from which to estimate SNP density and a PLINK file or GARLIC containing RoH segments to potentially 'repair'
# The required arguments are: Bim file
#                           : RoH file (the full name the file)
#                           : The window size (1 Mb recommended)
# The optional arguments are: Percentage threshold - what cutoff to use to define low density regions? Default is the top 5% lowest density regions on the chromsome or 1/50000 * window size (whichever is higher)
#                           : Percentage overlap - what is the minimum overlap required between low density region and a RoH gap in order to join the flanking segments? Default is 0.8
# USAGE:
# Rscript join_RoH_segs.R Ethiopian_MEGA_aut_geno-freqFiltered.bim Ethiopian_MEGA_RoH_miss2_het1_138SNPs.hom 1e6 0.05 0.8

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Defaults
perc.thresh <- 0.05
perc.overlap <- 0.8

# Error handling and updating default values
if (length(args)<3) {
  stop("At least three arguments must be supplied", call.=FALSE)
} else if (length(args)==4) {
  perc.thresh <- as.numeric(args[4])
} else if (length(args)==5) {
  perc.thresh <- as.numeric(args[4])
  perc.overlap <- as.numeric(args[5])
} else if (length(args)>5) {
  stop("Too many arguments", call.=FALSE)
}

options(stringsAsFactors=F)
Bim <- read.table(args[1])
WinSize <- as.numeric(args[3])

# Determine type of input segment file
cat("Determining input file type... ")
fileNameSplit <- strsplit(args[2], split=".", fixed=T)[[1]]
ext.ind <- match("hom", fileNameSplit)

if (length(fileNameSplit) > 1 & !is.na(ext.ind)) {
  stem <- paste(fileNameSplit[-ext.ind], collapse=".")
  ext <- "hom"
} else {
  stem <- args[2]
  ext="garlic"
}

if (ext=="hom") {
  Seg <- read.table(args[2], header=T)
  if (ncol(Seg)==13) {
    fileType <- "PLINK"
    colnames(Seg)[9] <- "LENGTH"
  }
  cat(fileType, "file type detected\n")
} else {
  Seg <- read.table(args[2], skip=1)
  Seg <- cbind(Seg, rep("placeholder_ID", nrow(Seg)))
  if (ncol(Seg)==10) {
    fileType <- "GARLIC"
    cat("WARNING: GARLIC file type detected - header line dropped\n")
    colnames(Seg)[c(1,2,3,4,5,10)] <- c("CHR","POS1","POS2","CLASS","LENGTH","IID")
    for (i in 1:22) {
      Seg[which(Seg$CHR==paste("chr", i, sep="")), "CHR"] <- i
    }
    Seg$CHR <- as.numeric(Seg$CHR)
    class.ranges <- list(); length(class.ranges) <- 3; names(class.ranges) <- c("A","B","C")
    for (class in c("A","B","C")) {
      class.ranges[[class]] <- range(as.numeric(Seg[which(Seg$CLASS==class), "LENGTH"]))
    }
    class.ranges[["A"]][2] <- mean(class.ranges[["A"]][2], class.ranges[["B"]][1])
    class.ranges[["B"]][1] <- mean(class.ranges[["A"]][2], class.ranges[["B"]][1]) + 0.5
    class.ranges[["B"]][2] <- mean(class.ranges[["B"]][2], class.ranges[["C"]][1]) 
    class.ranges[["C"]][1] <- mean(class.ranges[["B"]][2], class.ranges[["C"]][1]) + 0.5
  }
}

# Create copy of input file - eventual output file
Seg.new <- Seg

genome.SNP.count <- c()
chr.hists <- list()
  length(chr.hists) <- length(unique(Bim$V1))
  names(chr.hists) <- unique(Bim$V1)

for (chr in unique(Bim$V1)) {
  chr.data <- Bim[Bim$V1==chr,]
  start <- chr.data[1, "V4"]
  end <- chr.data[nrow(chr.data), "V4"]
  chr.hists[[as.character(chr)]] <- hist(chr.data$V4, breaks=seq(start, end+WinSize, WinSize), plot=F)
  genome.SNP.count <- c(genome.SNP.count, chr.hists[[as.character(chr)]]$counts)
}
min.SNPs <- max(sort(genome.SNP.count)[length(genome.SNP.count)*perc.thresh], 1/5e4*WinSize)

suppressMessages(library("GenomicRanges"))
for (chr in sort(unique(Seg$CHR))) {
  cat("Starting chr", chr, "\n")
  # Calculate SNP density across the chromosome in windows
  Hist <- chr.hists[[as.character(chr)]]
  Counts <- Hist$counts
  Breaks <- Hist$breaks
  # Store count information per region
  if (fileType=="PLINK") {
    chrWinInfo <- as.data.frame(matrix(ncol=4, nrow=length(Counts)))
      colnames(chrWinInfo) <- c("chr", "start", "end", "nSNPs")
  } else if (fileType=="GARLIC") {
    chrWinInfo <- as.data.frame(matrix(ncol=3, nrow=length(Counts)))
      colnames(chrWinInfo) <- c("chr", "start", "end")
  }
  chrWinInfo$chr <- chr
  chrWinInfo$start <- Breaks[1:length(Counts)]
  chrWinInfo$end <- Breaks[2:length(Breaks)]
  chrWinInfo$nSNPs <- Counts
  # Find low density regions
  low.dens.wins <- chrWinInfo[which(chrWinInfo$nSNPs < min.SNPs), 1:3]
  # If low density regions are found join continguous windows into regions
  if (nrow(low.dens.wins) > 0) {
    Low.Density.Regions <- low.dens.wins[1, ]
    if (nrow(low.dens.wins) > 1) {
      for (x in 2:nrow(low.dens.wins)) {
        if (low.dens.wins[x, "start"] == Low.Density.Regions[nrow(Low.Density.Regions), "end"]) {
          Low.Density.Regions[nrow(Low.Density.Regions), "end"] <- low.dens.wins[x, "end"]
        } else {
          Low.Density.Regions[nrow(Low.Density.Regions)+1, ] <- low.dens.wins[x, ]
        }
      }
    }
    Low.Density.Ranges <- makeGRangesFromDataFrame(Low.Density.Regions)
    # Join segments in input file
    for (Ind in unique(Seg$IID)) {
      chrSegInd <- Seg[Seg$CHR==chr&Seg$IID==Ind,]
        chrSegInd <- chrSegInd[order(chrSegInd$POS1),]
      if (nrow(chrSegInd)>1) {
        RoH.gaps <- as.data.frame(matrix(ncol=3,nrow=0))
          colnames(RoH.gaps) <- c("chr","start","end")
        for (y in 2:max(2, nrow(chrSegInd))) {
          RoH.gaps[y-1, ] <- c(chr, chrSegInd[y-1, "POS2"], chrSegInd[y, "POS1"])
        }
        RoH.gaps <- RoH.gaps[which(RoH.gaps$end > RoH.gaps$start),]
        if (nrow(RoH.gaps) > 0) {
          RoH.gap.Ranges <- makeGRangesFromDataFrame(RoH.gaps)
          # Do any gaps overlap a low density region?
          Hits <- findOverlaps(RoH.gap.Ranges, Low.Density.Ranges) # find all hits
          Overlaps <- pintersect(Low.Density.Ranges[subjectHits(Hits)], RoH.gap.Ranges[queryHits(Hits)])
          percentOverlap <- width(Overlaps) / width(RoH.gap.Ranges[queryHits(Hits)]) # calculate overlap
          Hits <- Hits[percentOverlap > perc.overlap] # sufficient overlap?
          joinOverGaps <- RoH.gaps[queryHits(Hits),]
          if (nrow(joinOverGaps)>0) {
            for (z in 1:nrow(joinOverGaps)) {
              seg1 <- rownames(chrSegInd[which(chrSegInd$POS2 == joinOverGaps[z, "start"]),])[1]
              seg2 <- rownames(chrSegInd[which(chrSegInd$POS1 == joinOverGaps[z, "end"]),])[1]
              joinedSeg <- chrSegInd[as.character(seg1),]
              joinedSeg$POS2 <- chrSegInd[as.character(seg2),"POS2"]
              joinedSeg$LENGTH <- chrSegInd[as.character(seg1),"LENGTH"] + chrSegInd[as.character(seg2),"LENGTH"]
              if (fileType=="PLINK") {
                joinedSeg$SNP2 <- chrSegInd[as.character(seg2),"SNP2"]
                joinedSeg$NSNP <- chrSegInd[as.character(seg1),"NSNP"] + chrSegInd[as.character(seg2),"NSNP"]
              } else if (fileType=="GARLIC") {
                if (joinedSeg$LENGTH < class.ranges[["A"]][1]) {
                  joinedSeg$CLASS <- "A"
                } else if (joinedSeg$LENGTH >= class.ranges[["A"]][1] & joinedSeg$LENGTH <= class.ranges[["A"]][2]) {
                  joinedSeg$CLASS <- "A"
                } else if (joinedSeg$LENGTH >= class.ranges[["B"]][1] & joinedSeg$LENGTH <= class.ranges[["B"]][2]) {
                  joinedSeg$CLASS <- "B"
                } else if (joinedSeg$LENGTH >= class.ranges[["C"]][1] & joinedSeg$LENGTH <= class.ranges[["C"]][2]) {
                  joinedSeg$CLASS <- "C"
                } else if (joinedSeg$LENGTH > class.ranges[["C"]][2]) {
                  joinedSeg$CLASS <- "C"
                }
              }
              Seg.new <- rbind(Seg.new, joinedSeg)
              # Remove the two split segments from the Seg.joined object
              Seg.new <- Seg.new[-match(as.character(seg1), rownames(Seg.new)),]
              Seg.new <- Seg.new[-match(as.character(seg2), rownames(Seg.new)),]
            }
          }
        }
      }
    }
  }
}

if (fileType=="PLINK") {
  colnames(Seg.new)[9] <- "KB"
  write.table(Seg.new, paste(stem, "_segsJoined.hom", sep=""), quote=F, row.names=F, col.names=T)
} else if (fileType=="GARLIC") {
  Seg.new <- Seg.new[,-ncol(Seg.new)]
  Seg.new$CHR <- paste("chr", Seg.new$CHR, sep="")
  write.table(Seg.new, paste(stem, "_segsJoined.garlic", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
}
