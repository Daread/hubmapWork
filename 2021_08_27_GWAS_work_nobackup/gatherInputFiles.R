
# Get functions
library(monocle3)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
# modStatus <- loadMonoclePackages()
library(tidyr)
library(ggplot2)
# Note: Need bedtools loaded for this! bedtools/2.29.2
library(cicero)

# Add this to stop R from formatting into scientific notation, which causes errors in bedtools (which expects explciit genomic coordinates)
options(scipen=999)

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
        # default="/net/bbi/vol1/data/genomes_stage/human/human_atac/gene_bodies.bed.gz", # This one from the ATAC dir isn't matching many RNA entries
        default="/net/bbi/vol1/data/genomes_stage/human/human_rna/latest.genes.bed",
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

opt$variableParams = paste0(opt$promoterUpstream, "_", opt$promoterDownstream, "_",
                           opt$coaccessCutoff, "_", opt$maxNdistalSites, "_", opt$peakSize )


getPromoterDF <- function(gzTSSfile, opt){
  promotersDF = as.data.frame(read.table(gzfile(gzTSSfile)))

  # # Get the upstream and downstream regions
  # promotersDF$start = ifelse(promotersDF$V6 == "+", 
  #                     (promotersDF$V2 - opt$promoterUpstream), 
  #                     (promotersDF$V3 - opt$promoterDownstream))
  # promotersDF$end = ifelse(promotersDF$V6 == "+", 
  #                     (promotersDF$V2 + opt$promoterDownstream), 
  #                     (promotersDF$V3 + opt$promoterUpstream))

  # Format into a DF that can be written as a bed file
  tssBed = data.frame("chr" = promotersDF$V1, "start" = promotersDF$V2, "end" = promotersDF$V3, "geneName" = promotersDF$V4,
                    "V5"=promotersDF$V5, "orientation"=promotersDF$V6)
  return(tssBed)
}

getPromoterDFnotGZ <- function(inputTSSfile, opt){
  promotersDF = as.data.frame(read.table((inputTSSfile), sep='\t', header=FALSE))

  # Get the upstream and downstream regions
  # promotersDF$start = ifelse(promotersDF$V6 == "+", 
  #                     ((promotersDF$V2 - opt$promoterUpstream)), 
  #                     ((promotersDF$V3 - opt$promoterDownstream)))
  # promotersDF$end = ifelse(promotersDF$V6 == "+", 
  #                     (promotersDF$V2 + opt$promoterDownstream), 
  #                     (promotersDF$V3 + opt$promoterUpstream))

  # Format into a DF that can be written as a bed file
  tssBed = data.frame("chr" = promotersDF$V1, "start" = promotersDF$V2, "end" = promotersDF$V3, "geneName" = promotersDF$V4,
                    "V5"=promotersDF$V5, "orientation"=promotersDF$V6)
  chromToKeep = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
                "17", "18", "19", "20", "21", "22", "X", "Y")
  tssBed = tssBed[tssBed$chr %in% chromToKeep,]

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

getCiceroPromoterIntersections <- function(opt, outFileTSS, outFileATAC_sites){
  # Sort the promoter and cicero peak bed files
  sortedATAC = paste0("./fileOutputs/ATAC_Peaks_", opt$variableParams, opt$ciceroFile, "_sorted.bed")
  if (!(file.exists(sortedATAC))){
    system(paste0("bedtools sort -i ", outFileATAC_sites, " > ", sortedATAC))
  }
  sortedTSS = paste0("./fileOutputs/promoterRegions", opt$variableParams, ".bed")
  if (!(file.exists(sortedTSS))){
    system(paste0("bedtools sort -i ", outFileTSS, " > ", sortedTSS))
  }
  

  # Intersect
  intersectBed = paste0("./fileOutputs/TSS", opt$variableParams, "_Intersected_Peaks_", opt$ciceroFile, ".bed")
  if (!(file.exists(intersectBed))){
    system(paste0("bedtools intersect -b ", sortedTSS, " -a ", sortedATAC, " -wa -wb > ", intersectBed))
  }

  # Get the intersected sites
  intersectedDF = as.data.frame(read.table(intersectBed))
  intersectedDF$formattedPeak = paste0("chr", intersectedDF$V5, "_", intersectedDF$V2, "_", intersectedDF$V3)
  # browser()

  intersectedDF$promoterCoord = paste0("chr", intersectedDF$V5, "_", intersectedDF$V6, "_", intersectedDF$V7)

  # Rename the GeneID
  colnames(intersectedDF) = c("PeakChr", "PeakStart", "PeakEnd", "Type", "PromoterChr", "PromoterStart", 
                            "PromoterEnd", "GeneID", "V9", "orientation", "formattedPeak", "promoterCoord")

  return(intersectedDF)
}

# Read the human gene_bodies.bed file into a dataframe
# promoterDF = getPromoterDF(opt$geneAnnotationsForTSS, opt) # for the atac one, in gz form
promoterDF = getPromoterDFnotGZ(opt$geneAnnotationsForTSS, opt) # for RNA-based gene annotations


# Write this output
outFileTSS = paste0("./fileOutputs/GeneBodies.bed")
if (!(file.exists(outFileTSS))){
  write.table(promoterDF, file=outFileTSS, 
    quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
}

# Get the cicero connection file to use here
ciceroDF = readRDS(paste0(opt$ciceroPath, opt$ciceroFile))
ciceroPeaksAsBed = getCiceroPeaksAsBed(ciceroDF)

# Output this
outFileATAC_sites = paste0("./fileOutputs/ATAC_Peaks_", opt$ciceroFile, opt$variableParams, ".bed")
if (!(file.exists(outFileATAC_sites))){
  write.table(ciceroPeaksAsBed, file=outFileATAC_sites, 
    quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
}


# Get the intersection.
# This holds peaks that intersected with a promoter
ciceroPromoterIntersection =  getCiceroPromoterIntersections(opt, outFileTSS, outFileATAC_sites)




# Get the subset of cicero links going to/from a TSS
ciceroTSS_DF = ciceroDF[(ciceroDF$Peak1 %in% ciceroPromoterIntersection$formattedPeak) ,]
ciceroTSS_DF$Peak2 = as.character(ciceroTSS_DF$Peak2)
# Get links above a cutoff
ciceroTSS_DF = ciceroTSS_DF[ciceroTSS_DF$coaccess > opt$coaccessCutoff,]  

# Set up a dataframe to hold the output coordinates of promoters + linked sites
# promoterDistalDF = data.frame("GeneID" = ciceroPromoterIntersection$GeneID,
#                               "chr" = ciceroPromoterIntersection$PromoterChr,
#                               "PromoterPos" = ciceroPromoterIntersection$promoterCoord,
#                               "orientation" = ciceroPromoterIntersection$orientation,
#                               "linkedSites" = 0)

promoterDistalDF = data.frame("GeneID" = promoterDF$geneName,
                              "chr" = promoterDF$chr,
                              "PromoterPos" = paste0("chr", promoterDF$chr, "_", as.character(promoterDF$start), 
                                                    "_", as.character(promoterDF$end)) ,
                              "orientation" = promoterDF$orientation,
                              "linkedSites" = 0)
# Note if a gene's promoter even had a peak assigned
promoterDistalDF$AnyPeakInProm = ifelse(promoterDistalDF$GeneID %in% ciceroPromoterIntersection$GeneID,
                                       1, 0)

# Get rid of any redundant entries
# promoterDistalDF = promoterDistalDF[!duplicated(promoterDistalDF),]

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
  # Make sure to drop links between peaks that are both intersecting with the promoter
  subsetCiceroDF = subsetCiceroDF[!(subsetCiceroDF$Peak2 %in% peaksInProm),]

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
          opt$promoterUpstream, "upstream", opt$promoterDownstream, "down", opt$coaccessCutoff, "coacCut", opt$peakSize, "peaks", ".png" ),
        width=1000, height=1000, res =200)
myPlot = ggplot(promoterDistalDF, aes(x=linkedSites)) + geom_histogram()
print(myPlot)
dev.off()


dfName = paste0("Gene_Prom_Plus_Distal_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff, opt$peakSize, "peaks",  ".RDS")



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





# Save this DF
saveRDS(promoterDistalDF, file=paste0("./rdsOutputs/", dfName))


# First I need to get a list of genes that represent highly expressed genes in each of the cell types I want to test for enrichment

getRNAdf <- function(opt, cellTypes){
  rnaDir = paste0("../2021_07_26_setup_ATAC_to_expr_data_nobackup/rdsOutputs/RNA_Data/" )
  rnaDF = readRDS(paste0(rnaDir, "Pseudobulked_RNA_countsPerMillionNotDeduplicated.RDS"))

  # # print(paste0(cellTypes, opt$predictionTask))
  # returnDF = rnaDF[,c("id", "gene_short_name", 
  #   paste0(cellTypes, opt$predictionTask))]

  return(rnaDF)
}

cellTypes = c("Adipocytes", "B_Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic_Endothelium", "Macrophage", "Mast_Cell", "Neuronal", "T_Cell",
                "Vascular_Endothelium", "VSM_and_Pericyte")

rnaData = getRNAdf(opt, cellTypes)

log2CPMcutoff = 2.0
log2RatioVsMeanCutoff = .1

library(reshape2)

writeBedFromDF <- function(promoterDistalDF, outputFile, opt){
  colToProcess = c("PromoterPos")
  for (eachCol in 1:opt$maxNdistalSites){
    colToProcess = c(colToProcess, paste0("Distal_", as.character(eachCol)))
  }

  coordDF = promoterDistalDF[,c(colToProcess, "GeneID", "orientation")]
  # Melt into format with one entry per site, instead of one row per gene+all its sites
  coordDF = melt(coordDF, id=c("GeneID", "orientation"))

  # Drop NAs
  coordDF = na.omit(coordDF)
  # Drop the entries with "None" written in
  coordDF = coordDF[!(coordDF[["value"]] == "None"),]

  # Format
  coordDF = separate(coordDF, value, into=c("chr", "start", "end"),
                  sep="_", remove=TRUE)
  coordDF$chr = gsub("chr", "", coordDF$chr)
  coordDF$placeHold = 1

  coordDF = coordDF[,c("chr", "start", "end", "GeneID", "placeHold", "orientation")]
  # Make this into a bed file
  write.table(coordDF, file=outputFile, 
    quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
}

library(dplyr)


outBedDir = paste0("./fileOutputs/bedFilesForGWAS_", opt$variableParams, "/")
dir.create(outBedDir)

# I need to get:
# 1: A single file holding all 
outputMasterFile = paste0(outBedDir, opt$variableParams, "_masterGeneAndPeakFile.bed")
writeBedFromDF(promoterDistalDF, outputMasterFile, opt)
# 2: A bed file for each of the 
for (eachCelltype in cellTypes){
  print(paste0("Working on ", eachCelltype))
  # Get a dataframe corresponding to this data
  miniDF = rnaData[,c("id", paste0(eachCelltype, c("_log2_CPM", "_log2_ratio_Vs_AllTypeMean")))]
  # Subset down
  miniDF = miniDF[miniDF[[paste0(eachCelltype, "_log2_CPM")]] > log2CPMcutoff,]
  miniDF = miniDF[miniDF[[paste0(eachCelltype, "_log2_ratio_Vs_AllTypeMean")]] > log2RatioVsMeanCutoff,]
  print(str(miniDF))
  # if (nrow(miniDF) == 0){
  #   next
  # }
  # Write the subset of the main DF that matches these genes
  thisOutputFile = paste0(outBedDir, opt$variableParams, "_", eachCelltype, "GeneAndPeakFile.bed")
  writeBedFromDF(promoterDistalDF[promoterDistalDF$GeneID %in% miniDF$id,], thisOutputFile, opt)
}









