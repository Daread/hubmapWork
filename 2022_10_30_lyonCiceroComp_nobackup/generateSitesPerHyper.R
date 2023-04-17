
.libPaths('1')


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
        # default="Cicero_Defaults_cds_p_allHeartATACFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_ciceroConnections.RDS",
        default="Cicero_Defaults_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p.rds_ciceroConnections.RDS", 
              help="Name of cds to read in and process", metavar="character"),

  make_option(c("-p", "--ciceroPath"), type="character", 
  			default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_04_redoCiceroWork_nobackup/rdsOutputs/", 
              help="Path to cds to process", metavar="character"),

  make_option(c("-l", "--linkSet"), type="character", 
        default="Cicero", 
              help="Path to cds to process", metavar="character"),

  make_option(c("-g", "--genomeFile"), type="character", 
  			# default="/net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished", 
        # The default is a soft link to "/net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished"
  			# default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/data/genomeData/bbiHg38.fasta",
        default = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/data/genomeData/backupBBIHumanGenome_nobackup/Homo_sapiens.GRCh38.dna.toplevel.fa.finished",
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


getCiceroPathData <- function(opt){
  # Set opt$peakLinkPath and opt$peakLinkFile based on defined inputs
  if (opt$linkSet == "Cicero"){
    opt$peakLinkPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_04_redoCiceroWork_nobackup/rdsOutputs/"
    opt$peakLinkFile = "Cicero_Defaults_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p.rds_ciceroConnections.RDS"
  } else if (opt$linkSet == "LyonV1"){
    opt$peakLinkPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/lyonFiles/"
    opt$peakLinkFile = "lyon_heart_cons_merged.Rds"
  }

  return(opt)
}


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

getPromoterDFnotGZ <- function(inputTSSfile, opt){
  promotersDF = as.data.frame(read.table((inputTSSfile), sep='\t', header=FALSE))

  # Get the upstream and downstream regions
  promotersDF$start = ifelse(promotersDF$V6 == "+", 
                      ((promotersDF$V2 - opt$promoterUpstream)), 
                      ((promotersDF$V3 - opt$promoterDownstream)))
  promotersDF$end = ifelse(promotersDF$V6 == "+", 
                      (promotersDF$V2 + opt$promoterDownstream), 
                      (promotersDF$V3 + opt$promoterUpstream))

  # Format into a DF that can be written as a bed file
  tssBed = data.frame("chr" = promotersDF$V1, "start" = promotersDF$start, "end" = promotersDF$end, "geneName" = promotersDF$V4,
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
  sortedATAC = paste0("./fileOutputs/ATAC_Peaks_", opt$variableParams, opt$peakLinkFile, "_sorted.bed")
  if (!(file.exists(sortedATAC))){
    system(paste0("bedtools sort -i ", outFileATAC_sites, " > ", sortedATAC))
  }
  sortedTSS = paste0("./fileOutputs/promoterRegions", opt$variableParams, ".bed")
  if (!(file.exists(sortedTSS))){
    system(paste0("bedtools sort -i ", outFileTSS, " > ", sortedTSS))
  }
  

  # Intersect
  intersectBed = paste0("./fileOutputs/TSS", opt$variableParams, "_Intersected_Peaks_", opt$peakLinkFile, ".bed")
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
# Only keep promoters that are 

# dir.create("./fileOutputs/")

# Write this output
outFileTSS = paste0("./fileOutputs/Promoters_Up", as.character(opt$promoterUpstream),
              "Down", as.character(opt$promoterDownstream), "_", opt$variableParams, ".bed")
if (!(file.exists(outFileTSS))){
  write.table(promoterDF, file=outFileTSS, 
    quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
}

# Get the cicero connection file to use here

opt = getCiceroPathData(opt)

ciceroDF = readRDS(paste0(opt$peakLinkPath, opt$peakLinkFile))

ciceroPeaksAsBed = getCiceroPeaksAsBed(ciceroDF)

# Output this
outFileATAC_sites = paste0("./fileOutputs/ATAC_Peaks_", opt$peakLinkFile, opt$variableParams, ".bed")
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
    subsetCiceroDF = subsetCiceroDF[!(duplicated(subsetCiceroDF$Peak2)),]
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
