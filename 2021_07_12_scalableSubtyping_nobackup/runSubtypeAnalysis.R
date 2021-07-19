

###########################################################################################
# Load relevant packages
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()
print(modStatus)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/scRNA_Seq_Background_Code/backgroundMethods.R")

redoProcessing <- function(inputCDS, useMNN=TRUE, mnnCategory="sample"){
	set.seed(7)
	inputCDS <- monocle3::estimate_size_factors(inputCDS)
	inputCDS <- preprocess_cds(inputCDS)

	if (useMNN){
		inputCDS = align_cds(inputCDS, alignment_group = mnnCategory)
	}
	inputCDS = reduce_dimension(inputCDS)

	return(inputCDS)
}

runAllClustering <- function(inputCDS, kToTry, kToKeep, processingNote, outputPath){
	# Loop through all k that aren't the kToKeep. Just make plots
	for (eachK in kToTry){
		if (eachK != kToKeep){
			print(paste0("Clustering with k = ", as.character(eachK)))
			inputCDS = cluster_cells(inputCDS, k=eachK)
			colData(inputCDS)$cluster_label = as.character(clusters(inputCDS))
			# Plot
			plotUMAP_Monocle(inputCDS, paste0(processingNote, "k=", as.character(eachK)),
						 "cluster_label", outputPath=outputPath)
		}
	}

	# Now do it for the correct k
	print(paste0("Clustering with k = ", as.character(kToKeep)))
	inputCDS = cluster_cells(inputCDS, k=kToKeep)
	colData(inputCDS)$cluster_label = as.character(clusters(inputCDS))
	colData(inputCDS)$partition_label = as.character(partitions(inputCDS))
	# Plot
	plotUMAP_Monocle(inputCDS, paste0(processingNote, "k=", as.character(kToKeep)),
				 "cluster_label", outputPath=outputPath)
	plotUMAP_Monocle(inputCDS, paste0(processingNote, "k=", as.character(kToKeep)),
				 "partition_label", outputPath=outputPath)
	# Return
	return(inputCDS)
}

getMarkerListFromGarnettFile <- function(inputCDS, inputMarkerFilePath, inputMarkerFileName){

    # Load Garnett    
    load_all("~/bin/garnett")
    library(org.Hs.eg.db)

    # Initial test of marker utility
    marker_file_path = paste0(inputMarkerFilePath, inputMarkerFileName, ".txt")
    marker_check <- check_markers(inputCDS, marker_file_path,
                              db=org.Hs.eg.db,
                              # cds_gene_id_type = "SYMBOL",
                              cds_gene_id_type = "ENSEMBL",
                              marker_file_gene_id_type = "SYMBOL")
    return(marker_check)
}


plotGarnettMarkersOnUMAP <- function(inputCDS, garnettMarkers, processingNote, outputPath){
	
	for (eachCelltype in as.character(levels(as.factor(garnettMarkers$cell_type)))){
		# Get the genes
		subsetDF = garnettMarkers[garnettMarkers$cell_type == eachCelltype,]
		theseGenes = as.character(subsetDF$marker_gene )
		# Plot
		plotUMAP_Monocle_genes(inputCDS, processingNote, theseGenes, paste0(eachCelltype,"_Markers"), 
					outputPath = outputPath)
	}
}



# Get the passed parameters
library("optparse")

option_list = list(
  make_option(c("-n", "--name"), type="character", default="TestRun", 
              help="Name of this subset analysis", metavar="character"),
    make_option(c("-c", "--clusters"), type="character", default="1,2", 
              help="comma-separated clusters to analyze", metavar="character"),
    make_option(c("-m", "--markerFile"), type="character", default="heartBroadMarkV2", 
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

myTestRes = runDEtestingToID_markers(myCDS, processingNote, "cluster_label",
									howManyGenesToTest = 50, outputPath=outputPath)

myTestRes = runDEtestingToID_markers(myCDS, processingNote, "partition_label",
									howManyGenesToTest = 50, outputPath=outputPath)
















