

###########################################################################################
# Load relevant packages
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()
print(modStatus)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/scRNA_Seq_Background_Code/backgroundMethods.R")
source("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_12_scalableSubtyping_nobackup/subtypeUtilFunctions.R")


# Get the passed parameters
library("optparse")

option_list = list(
  make_option(c("-n", "--name"), type="character",  default="Macrophage", # default="Endothelium", # default="TestRun", 
              help="Name of this subset analysis", metavar="character"),
    make_option(c("-c", "--clusters"), type="character", default="3,9,10,15", #default="7,12,16,22,13",#  default="1,2", 
              help="comma-separated clusters to analyze", metavar="character"),
    make_option(c("-m", "--markerFile"), type="character", default="macrophageSubtypesV2", #default="heartBroadMarkV2", 
              help="Marker file name", metavar="character"),
    make_option(c("-k", "--kval"), type="numeric", default=20, 
              help="Name of this subset analysis", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

processingNote = paste0("Heart_Subset_", opt$name)

# With arguments, read in data
rdsPath = "../2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/"
oldProcNote = "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40"
allCellCDS = readRDS(paste0(rdsPath, "allCells_", oldProcNote, ".rds"))

# Based off of this original CDS name, make a new directory for outputs
oldCDS_Path = paste0("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_12_scalableSubtyping_nobackup/plots/", oldProcNote, "/")
dir.create(file.path(oldCDS_Path), showWarnings=FALSE)
outputPath = paste0(oldCDS_Path, processingNote, "/")
dir.create(file.path(outputPath))

# Get the subset of clusters I want to look at here
clustToKeep = as.numeric(strsplit(opt$clusters, ",")[[1]])
myCDS = allCellCDS[,colData(allCellCDS)$cluster_label %in% clustToKeep]

# Re-do the pre-processing
myCDS = redoProcessing(myCDS, useMNN=TRUE, mnnCategory="sample")

# 10-14-21: Hard-assign parameters to plot outputs for 
source("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_12_scalableSubtyping_nobackup/subtypeUtilFunctions.R")
if (opt$name %in% c("Perivascular", "Endothelium", "Macrophage")){
	plotHardAssignedTypes(myCDS, opt, outputPath)
}


# Run the clustering at different levels
kToTry = c(10, 20, 50, 100)
kToKeep = as.numeric(opt$kval)
myCDS = runAllClustering(myCDS, kToTry, kToKeep, processingNote, outputPath)

processingNote = paste0(processingNote, "k=", as.character(kToKeep))

# Look at Garnett markers

garnettMarkerPath= "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/data/garnettModels/heartModels/"
markerFileToUse =  as.character(opt$markerFile) #"heartBroadMarkers"

garnettModelPath = makeGarnettModelHuman(myCDS, garnettMarkerPath, markerFileToUse,
                 processingNote, returnPath=TRUE)

# Get the markers to work with separately
garnettMarkers = getMarkerListFromGarnettFile(myCDS, garnettMarkerPath, markerFileToUse)
plotGarnettMarkersOnUMAP(myCDS, garnettMarkers, processingNote, outputPath)

# Now try running Garnett on the data and plot the results. Make sure to give it cluster data first
myCDS$garnett_cluster = colData(myCDS)$cluster_label
myCDS = applyGarnettModelHuman(myCDS, garnettModelPath)

# Plot
plotUMAP_Monocle(myCDS, paste0(processingNote, "_", markerFileToUse), "cluster_ext_type",
					show_labels=FALSE, outputPath=outputPath)
propDF = plotGroupedProportions(myCDS, paste0(processingNote, "_", markerFileToUse), 
					"cluster_label", "cluster_ext_type", pathToPlot=outputPath,
					widthToUse=2000)

# Also go ahead and find some markers, since the ones from literature may be uninformative

myTestResClust = runDEtestingToID_markers(myCDS, processingNote, "cluster_label",
									howManyGenesToTest = 50, outputPath=outputPath)


myTestRes = runDEtestingToID_markers(myCDS, processingNote, "partition_label",
									howManyGenesToTest = 50, outputPath=outputPath)




# 10-15-21: David working here

# clustRes = myTestResClust$marker_test_res
# clustRes = clustRes[clustRes$cell_group=="5",]
# clustRes = clustRes[order(clustRes$marker_score, decreasing=TRUE),]

# plotUMAP_Monocle_genes(myCDS, processingNote, c("DIAPH3", "TUBB", "TUBA1B"), paste0("Cytoskel_Markers"), 
#           outputPath = outputPath)


# plotUMAP_Monocle_genes(myCDS, processingNote, c("CLEC10A", "CD209A"), paste0("Dendritic_Markers"), 
#           outputPath = outputPath)


# plotUMAP_Monocle_genes(myCDS, processingNote, c("RBPJ"), paste0("General_Mac"), 
#           outputPath = outputPath)




