
# Get functions
library(monocle3)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
# modStatus <- loadMonoclePackages()
library(tidyr)
library(ggplot2)
# Note: Need bedtools loaded for this! bedtools/2.29.2
library(cicero)

library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-c", "--ciceroFile"), type="character", 
  			# default="Cicero_Defaults_cds_p_W144.heart.apex.s1_ciceroConnections.RDS", 
        default="Cicero_Defaults_cds_p_allHeartATACFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_ciceroConnections.RDS", 
              help="Name of cds to read in and process", metavar="character"),

  make_option(c("-p", "--ciceroPath"), type="character", 
  			default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_24_ciceroWork_nobackup/rdsOutputs/", 
              help="Path to cds to process", metavar="character"),

  make_option(c("-f", "--ciceroCDS"), type="character", 
        default="Cicero_Defaults_cds_p_allHeartATACFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_ciceroCDS.RDS", 
              help="Path to cds to process", metavar="character"),  

  make_option(c("-i", "--ciceroInputCDS"), type="character", 
        # default="Cicero_Defaults_cds_p_W144.heart.apex.s1_ciceroInput_CDS.RDS", 
        default="Cicero_Defaults_cds_p_allHeartATACFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_ciceroInput_CDS.RDS", 
              help="CDS before binning", metavar="character"),  
 # "_ciceroInput_CDS.RDS" 

  make_option(c("-g", "--genomeFile"), type="character", 
  			# default="/net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished", 
        # The default is a soft link to "/net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished"
  			default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/data/genomeData/bbiHg38.fasta",
              help="Path to genome fasta to use", metavar="character"),

  make_option(c("-t", "--geneAnnotationsForTSS"), type="character", 
        default="/net/bbi/vol1/data/genomes_stage/human/human_atac/gene_bodies.bed.gz",
              help="Path to genome fasta to use", metavar="character"),

  make_option(c("-u", "--promoterUpstream"), type="numeric", 
        default=1000,
              help="Bases upstream of TSS to use for input features", metavar="numeric"),
  make_option(c("-d", "--promoterDownstream"), type="numeric", 
        default=200,
              help="Bases downstream of TSS to use for input features", metavar="numeric"),

  make_option(c("-k", "--coaccessCutoff"), type="numeric", 
        default=.05,
              help="Cutoff for keeping coaccessible peaks", metavar="numeric"),

  make_option(c("-n", "--maxNdistalSites"), type="numeric", 
        default=5,
              help="Max number of sites to link to a gene's promoter", metavar="numeric"),

  make_option(c("-q", "--peakSize"), type="numeric", 
        default=600,
              help="-1 -> Keep as is, or enter a positive integer to make all peaks the same size", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



getPromoterDF <- function(gzTSSfile, opt){
  promotersDF = as.data.frame(read.table(gzfile(gzTSSfile)))

  # Get the upstream and downstream regions
  promotersDF$start = ifelse(promotersDF$V6 == "+", 
                      (promotersDF$V2 - opt$promoterUpstream), 
                      (promotersDF$V3 - opt$promoterDownstream))
  promotersDF$end = ifelse(promotersDF$V6 == "+", 
                      (promotersDF$V2 + opt$promoterDownstream), 
                      (promotersDF$V3 + opt$promoterUpstream))

  # Format into a DF that can be written as a bed file
  tssBed = data.frame("chr" = promotersDF$V1, "start" = promotersDF$start, "end" = promotersDF$end, "geneName" = promotersDF$V4,
                    "V5"=promotersDF$V5, "orientation"=promotersDF$V6)

  return(tssBed)
}


getCiceroPeaksAsBed <- function(ciceroDF){
  ciceroBaseSites = data.frame("Peaks" = levels(as.factor(ciceroDF$Peak1)))

  # Get peaks into bed format
  ciceroPeaksAsBed = data.frame(separate(data=ciceroBaseSites, col = "Peaks", into = c("chr", "start", "end"), sep="_", remove=TRUE))
  ciceroPeaksAsBed = separate(data=ciceroPeaksAsBed, col = "chr", into = c("Null", "chr"), sep="chr", remove=TRUE)
  ciceroPeaksAsBed$name = "Peak"
  ciceroPeaksAsBed = ciceroPeaksAsBed[c("chr", "start", "end", "name")]

  return(ciceroPeaksAsBed)
}

getCiceroPromoterIntersections <- function(opt){
  # Sort the promoter and cicero peak bed files
  sortedATAC = paste0("./fileOutputs/ATAC_Peaks_", opt$ciceroFile, "_sorted.bed")
  system(paste0("bedtools sort -i ", outFileATAC_sites, " > ", sortedATAC))
  sortedTSS = paste0("./fileOutputs/promoterRegions.bed")
  system(paste0("bedtools sort -i ", outFileTSS, " > ", sortedTSS))

  # Intersect
  intersectBed = paste0("./fileOutputs/TSS_Intersected_Peaks_", opt$ciceroFile, ".bed")
  system(paste0("bedtools intersect -b ", sortedTSS, " -a ", sortedATAC, " -wa -wb > ", intersectBed))

  # Get the intersected sites
  intersectedDF = as.data.frame(read.table(intersectBed))
  intersectedDF$formattedPeak = paste0("chr", intersectedDF$V5, "_", intersectedDF$V2, "_", intersectedDF$V3)

  intersectedDF$promoterCoord = paste0("chr", intersectedDF$V9, "_", intersectedDF$V6, "_", intersectedDF$V7)

  # Rename the GeneID
  colnames(intersectedDF) = c("PeakChr", "PeakStart", "PeakEnd", "Type", "PromoterChr", "PromoterStart", 
                            "PromoterEnd", "GeneID", "V9", "orientation", "formattedPeak", "promoterCoord")

  return(intersectedDF)
}


# Read the human gene_bodies.bed file into a dataframe
promoterDF = getPromoterDF(opt$geneAnnotationsForTSS, opt)

# Write this output
outFileTSS = paste0("./fileOutputs/Promoters_Up", as.character(opt$promoterUpstream),
              "Down", as.character(opt$promoterDownstream), ".bed")
write.table(promoterDF, file=outFileTSS, 
    quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

# Get the cicero connection file to use here
ciceroDF = readRDS(paste0(opt$ciceroPath, opt$ciceroFile))
ciceroPeaksAsBed = getCiceroPeaksAsBed(ciceroDF)

# Output this
outFileATAC_sites = paste0("./fileOutputs/ATAC_Peaks_", opt$ciceroFile, ".bed")
write.table(ciceroPeaksAsBed, file=outFileATAC_sites, 
    quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

# Get the intersection
ciceroPromoterIntersection =  getCiceroPromoterIntersections(opt)




# Get the subset of cicero links going to/from a TSS
ciceroTSS_DF = ciceroDF[(ciceroDF$Peak1 %in% ciceroPromoterIntersection$formattedPeak) ,]
ciceroTSS_DF$Peak2 = as.character(ciceroTSS_DF$Peak2)
# Get links above a cutoff
ciceroTSS_DF = ciceroTSS_DF[ciceroTSS_DF$coaccess > opt$coaccessCutoff,]  

# Set up a dataframe to hold the output coordinates of promoters + linked sites
promoterDistalDF = data.frame("GeneID" = ciceroPromoterIntersection$GeneID,
                              "chr" = ciceroPromoterIntersection$PromoterChr,
                              "PromoterPos" = ciceroPromoterIntersection$promoterCoord,
                              "orientation" = ciceroPromoterIntersection$orientation,
                              "linkedSites" = 0)
# Get rid of any redundant entries
promoterDistalDF = promoterDistalDF[!duplicated(promoterDistalDF),]

# Add columns for as many distal sites as can be linked, at most.
for (eachCol in 1:opt$maxNdistalSites){
  promoterDistalDF[paste0("Distal_", as.character(eachCol))] = "None"
}

# Get rid of NA's
ciceroTSS_DF = na.omit(ciceroTSS_DF)

# Now we need to go through and, for every TSS, get the subset of 
for (rowNum in 1:nrow(promoterDistalDF)){ # V8 gives the ensemble ID for a gene

  # Get the subset df
  thisGene = promoterDistalDF[rowNum, "GeneID"]
  peaksInProm = ciceroPromoterIntersection[ciceroPromoterIntersection$GeneID == thisGene,]
  peaksInProm = peaksInProm$formattedPeak
  subsetCiceroDF = ciceroTSS_DF[ciceroTSS_DF$Peak1 %in% peaksInProm, ]

  # If there are linked sites, order and get 
  if (nrow(subsetCiceroDF) != 0){
    # Sort
    subsetCiceroDF = subsetCiceroDF[order(-subsetCiceroDF$coaccess),]
    # Get up to 5
    sitesToGet = min(opt$maxNdistalSites, nrow(subsetCiceroDF))
    for (eachCoaccessInd in 1:sitesToGet){
      promoterDistalDF[rowNum, paste0("Distal_", as.character(eachCoaccessInd))] = (
            subsetCiceroDF$Peak2[eachCoaccessInd])
    }
    # Note how many this is
    promoterDistalDF[rowNum, "linkedSites"] = sitesToGet
  }

}

# Plot a histogram of how many sites per promoter
png(paste0("./plots/promoterRegionLinkedSitesHist", opt$maxNdistalSites, "maxSites", 
          opt$promoterUpstream, "upstream", opt$promoterDownstream, "down", opt$coaccessCutoff, "coacCut", ".png" ),
        width=1000, height=1000, res =200)
myPlot = ggplot(promoterDistalDF, aes(x=linkedSites)) + geom_histogram()
print(myPlot)
dev.off()

dfName = paste0("Gene_Prom_Plus_Distal_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff, ".RDS")

# Save this DF
saveRDS(promoterDistalDF, file=paste0("./rdsOutputs/", dfName))



updatedPromDistalDF = promoterDistalDF

# Check if need to resize distal peaks before getting sequence
if (opt$peakSize != -1){
  # Get the columns of distal coordiantes
  colToProcess = c()
  for (eachCol in 1:opt$maxNdistalSites){
    colToProcess = c(colToProcess, paste0("Distal_", as.character(eachCol)))
  }

  # Loop
  for (eachCol in colToProcess){
    # For each item
    promCount = 1
    for (eachEntry in strsplit(promoterDistalDF[[eachCol]], "_")){
      # browser()
      # Make sure not entry
      if (eachEntry[1] == "None"){
        # Return as is
        updatedPromDistalDF[promCount, eachCol] = "None"
      } else {
        # Get the midpoint
        midPoint = (as.numeric(eachEntry[2]) + as.numeric(eachEntry[3])) / 2.0
        updatedPromDistalDF[promCount, eachCol] = paste0(
              eachEntry[1], "_", as.character(as.integer(midPoint - opt$peakSize/2.0)),
                            "_", as.character(as.integer(midPoint + opt$peakSize/2.0)))
      }

      promCount = promCount + 1
    }
  }
}
promoterDistalDF = updatedPromDistalDF


writeBedFromCol <- function(promoterDistalDF, tempFile, colToWrite){
  coordDF = promoterDistalDF[,c(colToWrite, "GeneID", "orientation")]
  # Drop NAs
  coordDF = na.omit(coordDF)
  # Drop the entries with "None" written in
  coordDF = coordDF[!(coordDF[[colToWrite]] == "None"),]

  # Format
  coordDF = separate(coordDF, colToWrite, into=c("chr", "start", "end"),
                  sep="_", remove=TRUE)
  coordDF$chr = gsub("chr", "", coordDF$chr)
  coordDF$placeHold = 1

  coordDF = coordDF[,c("chr", "start", "end", "GeneID", "placeHold", "orientation")]
  # Make this into a bed file
  write.table(coordDF, file=tempFile, 
    quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
  # Done
}

library(dplyr)

getPromoterDistalSeq <- function(promoterDistalDF, opt){
  # Process each column
  colToProcess = c("PromoterPos")
  for (eachCol in 1:opt$maxNdistalSites){
    colToProcess = c(colToProcess, paste0("Distal_", as.character(eachCol)))
  }

  returnDF = promoterDistalDF[,!(names(promoterDistalDF) %in% colToProcess)]
  # Loop. For each of these, 
  #    1) write a bed file out of the coordiantes
  #    2) Get the fasta from this using bedtools
  #    3) Read it back in and set up a new DF
  tempFile = "./fileOutputs/tempBed.bed"
  seqFile = "./fileOutputs/tempSeq.fa"
  for (eachCol in colToProcess){
    print(paste0("Working on ", eachCol))
    writeBedFromCol(promoterDistalDF, tempFile, eachCol)
    # Get sequence
    system(paste0("bedtools getfasta -fi ", opt$genomeFile, " -bed ", tempFile, 
        " -tab > ", seqFile ))

    # Now read in the table
    peaksWithSeq = read.table(seqFile, sep='\t', header=FALSE,  col.names=c(eachCol, "seq"))
    # Fix formatting to match the original DF/peak format
    peaksWithSeq[[eachCol]] = paste0("chr", peaksWithSeq[[eachCol]])
    peaksWithSeq[[eachCol]] = gsub("-", "_", peaksWithSeq[[eachCol]])
    peaksWithSeq[[eachCol]] = gsub(":", "_", peaksWithSeq[[eachCol]])

    # Intersect in
    seqMergedDF = left_join(promoterDistalDF, peaksWithSeq[!(duplicated(peaksWithSeq[[eachCol]])),], by=eachCol)

    # Replace coordinates with the sequence
    returnDF[[eachCol]] = seqMergedDF$seq
  }
  
  return(returnDF)
}


# Part 2: Get sequence for all of these sites
promoterDistalDFWithSeqs = getPromoterDistalSeq(promoterDistalDF, opt)


# Save the output
dfName = paste0("Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize)

# Save this DF
saveRDS(promoterDistalDFWithSeqs, file=paste0("./rdsOutputs/", dfName, ".RDS"))








makeFastaFromSites <- function(promoterDistalDFWithSeqs, opt, dfName){
  # Loop through the promoter and all distal sites
  colToProcess = c("PromoterPos")
  for (eachCol in 1:opt$maxNdistalSites){
    colToProcess = c(colToProcess, paste0("Distal_", as.character(eachCol)))
  }

  # Loop through all
  sink(paste0("./fileOutputs/", dfName, ".fa"))
  for (eachCol in colToProcess){
    # For each item
    promCount = 1
    for (eachEntry in promoterDistalDFWithSeqs[[eachCol]] ){

      # If not NA, write to output
      if (!is.na(eachEntry)){
        thisGene = promoterDistalDFWithSeqs[promCount, "GeneID"]
        # Write the name of this sequence
        cat(paste0(">", thisGene, ":", eachCol, "\n"))
        # Write the sequence
        cat(eachEntry)
        cat("\n")
      }
      promCount = promCount + 1
    }
  }
  sink()
}



# Finally, go ahead and make a FASTA file for this sequence to feed into FIMO

makeFastaFromSites(promoterDistalDFWithSeqs, opt, dfName)

























