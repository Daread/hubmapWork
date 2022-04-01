


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
			set.seed(7)
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




 hardAssignSubclusters <- function(inputCDS, cellTypeName, kvalToUse){
 	# Endothelial assignments:
 	if ((cellTypeName == "Endothelium") & (kvalToUse == 10)){
 		colData(inputCDS)$Subtype = "Unassigned"
 		colData(inputCDS)$Subtype = ifelse(colData(inputCDS)$cluster_label %in% c(2,6,14,12),
 						"Venous", colData(inputCDS)$Subtype)
 		colData(inputCDS)$Subtype = ifelse(colData(inputCDS)$cluster_label %in% c(7,11),
 						"Arterial", colData(inputCDS)$Subtype)
 		colData(inputCDS)$Subtype = ifelse(colData(inputCDS)$cluster_label %in% c(5,4,1,13,15),
 						"Capillary", colData(inputCDS)$Subtype)
 		colData(inputCDS)$Subtype = ifelse(colData(inputCDS)$cluster_label %in% c(10,9),
 						"Lymphatic", colData(inputCDS)$Subtype)
 		colData(inputCDS)$Subtype = ifelse(colData(inputCDS)$cluster_label %in% c(8,3),
 						"Endocardium", colData(inputCDS)$Subtype)
 	}
 	if ((cellTypeName == "Perivascular") & (kvalToUse == 10)){
 		colData(inputCDS)$Subtype = "Unassigned"
 		colData(inputCDS)$Subtype = ifelse(colData(inputCDS)$cluster_label %in% c(1,2,4),
 						"Pericyte", colData(inputCDS)$Subtype)
 		colData(inputCDS)$Subtype = ifelse(colData(inputCDS)$cluster_label %in% c(3),
 						"Vascular Smooth Muscle", colData(inputCDS)$Subtype)
 	}
 	 	if ((cellTypeName == "Macrophage") & (kvalToUse == 20)){
 		colData(inputCDS)$Subtype = "Undetermined"
 		colData(inputCDS)$Subtype = ifelse(colData(inputCDS)$cluster_label %in% c(1,2,3,11,9),
 						"M-S1: COLEC12+", colData(inputCDS)$Subtype)
 		colData(inputCDS)$Subtype = ifelse(colData(inputCDS)$cluster_label %in% c(4,5),
 						"M-S2: DOCK10+", colData(inputCDS)$Subtype)
 		colData(inputCDS)$Subtype = ifelse(colData(inputCDS)$cluster_label %in% c(6),
 						"DIAPH3+", colData(inputCDS)$Subtype)
 	}



 	return(inputCDS)
 }


plotHardAssignedTypes <- function(inputCDS, opt, outputPath){
	set.seed(7)
	# Cluster with the right k value
	if (opt$name %in%c("Perivascular", "Endothelium")){
		kvalToUse = 10
		inputCDS = cluster_cells(inputCDS, k=kvalToUse)
		colData(inputCDS)$cluster_label = as.character(clusters(inputCDS))
	}
	if (opt$name %in%c("Macrophage")){
		kvalToUse = 20
		inputCDS = cluster_cells(inputCDS, k=kvalToUse)
		colData(inputCDS)$cluster_label = as.character(clusters(inputCDS))
	}
	plotUMAP_Monocle(inputCDS, paste0(opt$name, "_k=", as.character(kvalToUse)),
				 "cluster_label", outputPath=outputPath)

	# Hard assign the cluster labels
	inputCDS = hardAssignSubclusters(inputCDS, opt$name, kvalToUse)

	# Now plot these
	png(paste0(outputPath, opt$name, "_SubtypesLabeled.png"), res=300, height = 1500, width=1500)
	myPlot <- (plot_cells(inputCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by="Subtype", label_cell_groups=FALSE,
          cell_stroke=.1 # , group_label_size=textSize        
                ))
    myPlot = (myPlot + theme(text=element_text(size=12)))
    print(myPlot)
	dev.off()
}

