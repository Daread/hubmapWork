
# 4-7-22 add

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
        default="LyonV1",   #"Cicero", 
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
        default=2000,
              help="Bases upstream of TSS to use for input features", metavar="numeric"),
  make_option(c("-d", "--promoterDownstream"), type="numeric", 
        default=1000,
              help="Bases downstream of TSS to use for input features", metavar="numeric"),

  make_option(c("-k", "--coaccessCutoff"), type="numeric", 
        default=.015,
              help="Cutoff for keeping coaccessible peaks", metavar="numeric"),

  make_option(c("-n", "--maxNdistalSites"), type="numeric", 
        default=5,
              help="Max number of sites to link to a gene's promoter", metavar="numeric"),

  make_option(c("-i", "--seqInfOnly"), type="logical", 
        default=FALSE,
              help="Only Output Info On linked sites", metavar="logical"),

  make_option(c("-q", "--peakSize"), type="numeric", 
        default=600,
              help="-1 -> Keep as is, or enter a positive integer to make all peaks the same size", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(opt)

opt$variableParams = paste0(opt$promoterUpstream, "_", opt$promoterDownstream, "_",
                           opt$coaccessCutoff, "_", opt$maxNdistalSites, "_", opt$peakSize )


getProteinCodingGenes <- function(rnaData, opt){
  # Get the gtf file
  gtfFile = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/fileOutputs/protCodingOnlyHumanGTF.gtf"
  gtfData = read.table((gtfFile), sep="\t")
  splitData = separate(gtfData, "V9", c("geneTxt", "idTxt", "ensemblName"))

  # Get unique names
  ensemblIDsProtein = unique(splitData[["ensemblName"]])
  protRNA = rnaData[rnaData$GeneID %in% ensemblIDsProtein,]
  return(protRNA)
}

# testRowForSig <- function(x, distSiteCount) {
# 	allDistNames = paste0("Distal_", as.character(1:distSiteCount))

# 	for (eachName in allDistNames){

# 	}
#   if (character == FALSE) {
#     x ^ 2
#   } else {
#     as.character(x ^2)
#   }
# }

assignAgeRelatedGenes = function(inputDF, distSiteCount){
	# Get DE Genes
	ageDEfile = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/plots/DE_Summaries/Combined_DE_by_Age_q_0.1.csv"
	ageDE = read.csv(ageDEfile)
	# Get ensembl names
	heartCDS = readRDS("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/allCells_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes.rds")
	shortNameVec = rownames(heartCDS)
	names(shortNameVec) = rowData(heartCDS)$gene_short_name

	ageDE$ensemblName = shortNameVec[ageDE$gene]

	ageGenes = unique(ageDE$ensemblName)

	inputDF$IsAgeDE = ifelse(inputDF$GeneID %in% ageGenes, "Yes", "No")
	

	# Get scores of distal sites' age-covariance
	agePeakFile = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/rdsOutputs/RNA_Data/dr_age_specific.Rds"
	agePeakData = readRDS(agePeakFile)

	ageScores = agePeakData$effect_by_error
	names(ageScores) = agePeakData$peak
	ageScores = setNames(c(ageScores, NA), c(names(ageScores), "None"))

	dir.create("./plots/lyonPlots/")
	png("./plots/lyonPlots/ageScoreHist.png")
	myPlot = ggplot(agePeakData[(agePeakData$effect_by_error < 5) & (agePeakData$effect_by_error > -5),], aes(x=effect_by_error)) + 
		geom_histogram()
	print(myPlot)
	dev.off()


	# Loop for each distal site
	for (eachSite in 1:distSiteCount){
		siteName = paste0("Distal_", as.character(eachSite))
		print(siteName)
		inputDF[[paste0(siteName, "_Score")]] = ageScores[inputDF[[siteName]]]
	}
	

	# Get count with 
	ageDEonly = inputDF[inputDF$IsAgeDE == "Yes",]
	allDistNames = paste0("Distal_", as.character(1:distSiteCount), "_Score")
	ageDEonly$ScoredSites = rowSums(!is.na(ageDEonly[allDistNames]))
	ageDEonly$HasScoredSites = ageDEonly$ScoredSites > 0

	print("Linked to any:")
	print(sum(as.numeric(ageDEonly$HasScoredSites)))

	replaceNAdf = ageDEonly
	replaceNAdf[is.na(replaceNAdf)] = 0
	replaceNAdf$AgeSpecPeak = rowSums(abs(ageDEonly[allDistNames]) > 1.96)

	

	ageSpecPeakCount = replaceNAdf$AgeSpecPeak[!is.na(replaceNAdf$AgeSpecPeak)]
	print("Linked to sig")
	print(sum(ageSpecPeakCount > 0))

	siteVec = character()
	allDistSites = paste0("Distal_", as.character(1:distSiteCount))
	for (eachCol in allDistSites){
		siteVec = c(siteVec, ageDEonly[[eachCol]])
	}

	browser()

	uniqueSites = unique(siteVec)
	siteFile = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_10_30_lyonCiceroComp_nobackup/fileInputs/ageDEsites.rds"
	saveRDS(uniqueSites, file =siteFile)

	# Get sites with >1.96 test statistic
	statSigSites = agePeakData$peak[abs(agePeakData$effect_by_error) > 1.96]
	uniqueSigSites = uniqueSites[uniqueSites %in% statSigSites]
	siteFile = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_10_30_lyonCiceroComp_nobackup/fileInputs/sigAgeDEsites.rds"
	saveRDS(uniqueSigSites, file =siteFile)

	# ageDEonly = ageDEonly[c("")]

	print("Hey")
	print("Hey")
	print("Hey")
	print("Hey")
	print("Hey")

	return(inputDF)

}

outputSeqContent = function(inputDF, opt, distSiteCount){
  distalCols = paste0("Distal_", as.character(1:distSiteCount))

  # Only count protein coding genes
  print("Getting Prot Coding Genes")
  inputDF = getProteinCodingGenes(inputDF)
  # Get empty entry counts
  miniDF = inputDF[distalCols]

  emptyEntries = sum(miniDF == "None")
  linkedSites = (distSiteCount * nrow(miniDF)) - emptyEntries

  linkedSeq = linkedSites * as.numeric(opt$peakSize)

  outputDF = data.frame("Upstream" = opt$promoterUpstream, "Downstream" = opt$promoterDownstream,
                  "maxSites" = opt$maxNdistalSites, "Cutoff" = opt$coaccessCutoff,
                  "Peak"= opt$peakSize, "LinkedCount" = linkedSites, "LinkedBP" = linkedSeq)

  # browser()

  #Output the result
  outFile = paste0("./plots/promoterRegionLinkedSiteCounts_", opt$linkSet, "_", opt$maxNdistalSites, "maxSites", 
          opt$promoterUpstream, "upstream", opt$promoterDownstream, "down", opt$coaccessCutoff, "coacCut", opt$peakSize, "peaks", ".csv" )
  write.csv(outputDF, file=outFile)

  # 11-1-22: Look at age-related content
  inputDF = assignAgeRelatedGenes(inputDF, distSiteCount)

  # Plot the histogram distribution


}

outputSeqContent(promoterDistalDF, opt, opt$maxNdistalSites)



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
print("Getting link data")
opt = getCiceroPathData(opt)

ciceroDF = readRDS(paste0(opt$peakLinkPath, opt$peakLinkFile))

ciceroPeaksAsBed = getCiceroPeaksAsBed(ciceroDF)

# Output this
outFileATAC_sites = paste0("./fileOutputs/ATAC_Peaks_", opt$peakLinkFile, opt$variableParams, ".bed")
if (!(file.exists(outFileATAC_sites))){
  write.table(ciceroPeaksAsBed, file=outFileATAC_sites, 
    quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)
}


print("Finding intersection")
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

print("Finding linked sites per TSS")

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

# dir.create("./plots/")
# Plot a histogram of how many sites per promoter
png(paste0("./plots/promoterRegionLinkedSitesHist", opt$maxNdistalSites, "maxSites", 
          opt$promoterUpstream, "upstream", opt$promoterDownstream, "down", opt$coaccessCutoff, "coacCut", opt$peakSize, "peaks", ".png" ),
        width=1000, height=1000, res =200)
myPlot = ggplot(promoterDistalDF, aes(x=linkedSites)) + geom_histogram()
print(myPlot)
dev.off()

if (opt$seqInfOnly){
  print("Getting Seq Info")
  outputSeqContent(promoterDistalDF, opt, opt$maxNdistalSites)
  quit(save="no")
} else{
  print("Getting and outputting sequence for linked sites")
}






# dir.create("./rdsOutputs/")

dfName = paste0("Gene_Prom_Plus_Distal_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff, opt$peakSize, "peaks",  ".RDS")

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
  tempFile = paste0( "./fileOutputs/tempBed", opt$variableParams, ".bed")
  seqFile = paste0("./fileOutputs/tempSeq", opt$variableParams, ".fa")
  for (eachCol in colToProcess){
    # See if this is an empty column. If so, skip
    # if (levels(as.factor(promoterDistalDF[[eachCol]])) == "None"){
    if ( length(unique(promoterDistalDF[[eachCol]])) == 1 ){
      print(paste0("No links for ", eachCol))
      returnDF[[eachCol]] = "None"
    } else {
      print(paste0("Working on ", eachCol))
      writeBedFromCol(promoterDistalDF, tempFile, eachCol)
      # Get sequence
      system(paste0("bedtools getfasta -s -fi ", opt$genomeFile, " -bed ", tempFile, 
          " -tab > ", seqFile ))
      # Now read in the table
      peaksWithSeq = read.table(seqFile, sep='\t', header=FALSE,  col.names=c(eachCol, "seq"))

      # Fix formatting to match the original DF/peak format
      peaksWithSeq[[eachCol]] = paste0("chr", peaksWithSeq[[eachCol]])
      peaksWithSeq[[eachCol]] = gsub("-", "_", peaksWithSeq[[eachCol]])
      peaksWithSeq[[eachCol]] = gsub(":", "_", peaksWithSeq[[eachCol]])
      # Strandedness option adds a (+) or (-) to the end
      peaksWithSeq[[eachCol]] = gsub("\\(\\+\\)", "", peaksWithSeq[[eachCol]])
      peaksWithSeq[[eachCol]] = gsub("\\(_\\)", "", peaksWithSeq[[eachCol]])

      # Intersect in
      seqMergedDF = left_join(promoterDistalDF, peaksWithSeq[!(duplicated(peaksWithSeq[[eachCol]])),], by=eachCol)

      # Replace coordinates with the sequence
      returnDF[[eachCol]] = seqMergedDF$seq
    }

  }  
  return(returnDF)
}


# Part 2: Get sequence for all of these sites
promoterDistalDFWithSeqs = getPromoterDistalSeq(promoterDistalDF, opt)

addTrainTestDesignation <- function(promoterDistalDFWithSeqs){
  promoterDistalDFWithSeqs$Model_Set = "Undefined"
  promoterDistalDFWithSeqs$Model_Set = ifelse(promoterDistalDFWithSeqs$chr %in% c("6", "16"), # 2334 genes
              "Test", promoterDistalDFWithSeqs$Model_Set)
  promoterDistalDFWithSeqs$Model_Set = ifelse(promoterDistalDFWithSeqs$chr %in% c("5", "15"), # 1958 genes
              "Validation", promoterDistalDFWithSeqs$Model_Set)
  promoterDistalDFWithSeqs$Model_Set = ifelse(promoterDistalDFWithSeqs$chr %in% 
                c("1", "2", "3", "4", "7", "8", "9", "10", "11", "12", "13", "14", "17", "18", "19", "20", "21", "22", "X"),
              "Train", promoterDistalDFWithSeqs$Model_Set)

  return(promoterDistalDFWithSeqs)
}

# Add notes on if train/validation/test
promoterDistalDFWithSeqs = addTrainTestDesignation(promoterDistalDFWithSeqs)

# Save the output
dfName = paste0("Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize, "_Links_", opt$linkSet)

# Save this DF
saveRDS(promoterDistalDFWithSeqs, file=paste0("./rdsOutputs/", dfName, ".RDS"))

# Save as a csv for easier swapping over into python for CNN work
write.csv(promoterDistalDFWithSeqs, paste0("./rdsOutputs/", dfName, ".csv"))






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


























