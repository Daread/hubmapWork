
readHeartCDS <- function(pathToBBIdirs){
	# List of samples to read in. Note that W137 failed in this prep, so I'm reading in all 14 others
	samplesToRead <- c("W134.Apex", "W135.Left.Vent", "W136.Apex", "W136.Left.Vent",
                   "W139.Apex", "W139.Left.Vent", "W139.Right.Vent",
                   "W139.Septum", "W142.Left.Vent", "W144.Apex", "W145.Apex",
                   "W145.Left.Vent", "W146.Apex", "W146.Left.Vent")
	# Get the data
	return(makeCDS_fromBBIoutputDirs(samplesToRead, pathToSample=pathToBBIdirs))
}

getQCParams = function(){
	# Set up parameters that need to be specified
	qcParams = list()
	qcParams["useScrublet"] = TRUE
	qcParams["scrubCutoff"] = .2
	qcParams["writeMatrixForScrub"] = TRUE
	qcParams["umiCellMin"] = 100
	qcParams["bgUMImax"] = 15
	qcParams["useMitoCutoff"] = TRUE
	qcParams["mitoCutoff"] = 10.0
	# Background Correction Options
	qcParams["usePackerBGcorrection"] = FALSE
	qcParams["useSouopX"] = FALSE
	qcParams["soupX_assumption"] = .2
	qcParams["useMNN"] = TRUE
	qcParams["mnnCategory"] = "sampleName"
	return(qcParams)
}

discardLowUMI_highMito = function(inputCDS, inputQClist){
	# Simplest part is to only keep cells above a particular UMI
	inputCDS = inputCDS[,colData(inputCDS)$n.umi >= inputQClist[["umiCellMin"]]]
	# Toss mitochondrial
	if (inputQClist[["useMitoCutoff"]]){
		inputCDS = inputCDS[,colData(inputCDS)$perc_mitochondrial_umis < inputQClist[["mitoCutoff"]]]
	}
	return(inputCDS)
}

getCorrectedCDS = function(inputCellCDS, inputFullCDS, qcParamList){
	# Set up a background CDS
	backgroundCDS = inputFullCDS[,colData(inputFullCDS)$n.umi < qcParamList[["bgUMImax"]]]

	# Disard cells with high scrublet scores from the cellCDS
	inputCellCDS = inputCellCDS[,colData(inputCellCDS)$new_scrub_score < qcParamList[["scrubCutoff"]]]

	# Now run the Packer BG correction method.
	if (qcParamList[["usePackerBGcorrection"]]){
		bgCorrection_List = findBackground_Packer_Method(inputCellCDS, backgroundCDS, returnOnlyMatrix=FALSE)
		inputCellCDS = (bgCorrection_List[["CDS"]])
	}
	# Return the updated CDS
	return(inputCellCDS)
}

addPostScrubletProcNote <- function(processingNote, inputQClist){
	# Scrublet?
	if (inputQClist[["useScrublet"]]){
		processingNote = paste0(processingNote, "Scrub=", as.character(inputQClist[["scrubCutoff"]]))
	} else{
		processingNote = paste0(processingNote, "NoScrub")
	}
	# Packer BG?
	if (inputQClist[["usePackerBGcorrection"]]){
		processingNote = paste0(processingNote, "packerBG")
	} else{
		processingNote = paste0(processingNote, "noPacker")
	}
	return(processingNote)
}

###########################################################################################
# Load relevant packages
source("../../../sharedProjectCode/utility/singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()
print(modStatus)
source("../../../sharedProjectCode/scRNA_Seq_Background_Code/backgroundMethods.R")
processingNote = "HM10"

# Set processing/QC parameters and read in the initial CDS
qcParams = getQCParams()
fullCDS <- readHeartCDS("../../../fibrosis/results/2021_03_29_HM10_Human15_ResetPlates_nobackup/")
colData(fullCDS)$log10_umi = log10(colData(fullCDS)$n.umi)
set.seed(7)

# 9-29-21: Make a matrix to give to Matt at Hubmap
mattCDS = fullCDS[,fullCDS$sampleName == "W144.Apex"]
mattCDS = mattCDS[,mattCDS$n.umi > 99]

# Output
mattDir = "./mattOutputs/"
dir.create(mattDir)
writeMM(exprs(mattCDS), paste0(mattDir, "W144_Apex_countMatrix.MM"))
write.csv(as.data.frame(colData(mattCDS)), paste0(mattDir, "W144_cellMetadata.csv") )
write.csv(as.data.frame(rowData(mattCDS)), paste0(mattDir, "W144_geneMetadata.csv") )

# Make some core plots of QC metrics before any correction/filtering
# makeQCplots(fullCDS[,colData(fullCDS)$n.umi > 99], processingNote="UMI>=100", makeUMAP=TRUE,
# 	genesToShow = c("GSN", "LDB2", "KCNAB1", "RBPJ", "TTN", "THEMIS"), geneSetName="MiscHeartMarkers")
# makeQCplots(fullCDS[,colData(fullCDS)$n.umi > 99], processingNote="UMI>=100", makeUMAP=TRUE)

# The easiest way to run Scrublet is, annoyingly, via Python. This can handle all parts of QC apart from the Scrublet call
cellCDS = discardLowUMI_highMito(fullCDS, qcParams)
processingNote = paste0(processingNote, "UMI=", as.character(qcParams[["umiCellMin"]]))
if (qcParams[["useMitoCutoff"]]) {
	processingNote = paste0(processingNote, "_mito=", as.character(qcParams[["mitoCutoff"]]))
}

# Write MM File for scrublet
writeMM(t(exprs(cellCDS)), paste0("./scrubletIntermediates/", processingNote, "filteredCellExpr.MM"))
# Run scrublet
scrubletCommand = paste0("python ../../../sharedProjectCode/utility/runScrublet.py ", processingNote, " ./scrubletIntermediates/")
print("source activate scrublet_anaconda3")
scrubletCommand

# In parallel, run the scrublet helper function command listed above, then proceed:
newScrubResults = read.csv(paste0("./scrubletIntermediates/", processingNote, "_scrubletScoresPreds.csv"))
colData(cellCDS)$new_scrub_score = newScrubResults$New_Scrub_Score

# Now, do background correction based on discarding high scrublet scores, SoupX, and Packer methods
# backgroundCDS = fullCDS[,colData(fullCDS)$n.umi < qcParams[["bgUMImax"]]]
cellCDS = getCorrectedCDS(cellCDS, fullCDS, qcParams)
processingNote = addPostScrubletProcNote(processingNote, qcParams)

# First, need to re-find size factors after discarding so many low-quality "nuclei"/"cells"
set.seed(7)
cellCDS <- monocle3::estimate_size_factors(cellCDS)
cellCDS <- preprocess_cds(cellCDS)

# Use MNN?
set.seed(7)
if (qcParams[["useMNN"]] ){
	cellCDS <- align_cds(cellCDS, alignment_group = qcParams[["mnnCategory"]])
	processingNote = paste0(processingNote, "MNN=", qcParams[["mnnCategory"]])
} else{
	processingNote = paste0(processingNote, "MNN=none" )
}
cellCDS <- reduce_dimension(cellCDS)

# QC of the remaining cells
makeQCplots(cellCDS[,colData(cellCDS)$n.umi > 99], processingNote=processingNote, makeUMAP=TRUE,
	genesToShow = c("GSN", "LDB2", "KCNAB1", "RBPJ", "TTN", "THEMIS"), geneSetName="MiscHeartMarkers")

plotUMAP_Monocle(cellCDS, processingNote, "sampleName", 
	outputPath=paste0("./plots/QC_", processingNote, "/"), show_labels=FALSE)

##################################################################################################################################################################

handleClustering = function(inputCDS, processingNote, inputKval){
	inputCDS = cluster_cells(inputCDS, k = inputKval)
	colData(inputCDS)$partition_label = as.character(partitions(inputCDS))
	colData(inputCDS)$cluster_label = as.character(clusters(inputCDS))
	# Automatically generate plots for the clusters and partitions
	dir.create(file.path("./plots", paste0(processingNote, "K=", as.character(inputKval))), showWarnings=FALSE)
	outputPath = paste0("./plots/", paste0(processingNote, "K=", as.character(inputKval)), "/")
	# Make the plots
	plotUMAP_Monocle(inputCDS, paste0("k=", as.character(inputKval), "_", processingNote), "cluster_label", show_labels=TRUE, outputPath=outputPath)
	plotUMAP_Monocle(inputCDS, paste0("k=", as.character(inputKval), "_", processingNote), "partition_label", show_labels=TRUE, outputPath=outputPath)
	# Return
	return(inputCDS)
}

# Cluster cells
print("Clustering Cells Now")
set.seed(7)
kVal = 40
cellCDS = handleClustering(cellCDS, processingNote, kVal)



# Try working with some high-level cell types
library(garnett)

garnettMarkerPath= "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/data/garnettModels/heartModels/"
markerFileToUse =  "heartBroadMarkV2" #"heartBroadMarkers"

garnettModelPath = makeGarnettModelHuman(cellCDS, garnettMarkerPath, markerFileToUse,
                 processingNote, returnPath=TRUE)

# garnettModelPath = paste0(garnettMarkerPath, markerFileToUse, "_trainedOn_",
# 			processingNote, ".rds")

# Cluster, so I can use cluster_extend
colData(cellCDS)$garnett_cluster = colData(cellCDS)$partition_label
cellCDS = applyGarnettModelHuman(cellCDS, garnettModelPath)

# Take a look at the results
outputPath = paste0("./plots/", paste0(processingNote, "K=", as.character(kVal)), "/")
plotUMAP_Monocle(cellCDS, paste0("k=", as.character(kVal), "_", processingNote, markerFileToUse), "cluster_ext_type", show_labels=FALSE, outputPath=outputPath)

# Takee a look at the distribution of the cell_type calls vs. clusters
propDF = plotGroupedProportions(cellCDS, paste0(processingNote, markerFileToUse), "partition_label", "cluster_ext_type", 
        pathToPlot=paste0("./plots/", processingNote, "K=", as.character(kVal), "/"), 
        widthToUse=1800)
# propTotDF = (propDF %>% group_by(partition_label) %>% summarise(propSum=sum(Freq)))

















# Write an output of the cells as they've been processed. Make a full CDS and one
#    with just the endothelial + vascular smooth muscle cells
processingNoteWithK = paste0(processingNote, "K=", as.character(kVal))
saveRDS(cellCDS, paste0("./rdsOutput/allCells_", processingNoteWithK, ".rds"))

# Now get the endothelial CDS
endothOrVSMclusts = c("4", "7", "8", "6", "13")
endothCDS = cellCDS[,colData(cellCDS)$partition_label %in% endothOrVSMclusts]
saveRDS(endothCDS, paste0("./rdsOutput/endothelium_", processingNoteWithK, ".rds"))
























