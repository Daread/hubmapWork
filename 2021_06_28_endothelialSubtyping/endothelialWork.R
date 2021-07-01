
###########################################################################################
# Load relevant packages
source("../../../sharedProjectCode/utility/singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()
print(modStatus)
source("../../../sharedProjectCode/scRNA_Seq_Background_Code/backgroundMethods.R")
processingNote = "HM10_Endoth_and_SMC"

# Load the saved file in
rdsPath = "../2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/"
oldProcNote = "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40"

# Read in the endothelium
endothCDS = readRDS(paste0(rdsPath, "endothelium_", oldProcNote, ".rds"))

# Get the markers
garnettMarkerPath = "../../data/garnettModels/heartModels/"
markerFileToUse = "endothSubtypesLitvV1"

# Try running garnett
set.seed(7)
garnettModelPath = makeGarnettModelHuman(endothCDS, garnettMarkerPath, markerFileToUse,
                 processingNote, returnPath=TRUE)

# Try this out
colData(endothCDS)$garnett_cluster = colData(endothCDS)$partition_label
endothCDS = applyGarnettModelHuman(endothCDS, garnettModelPath)

# Take a look at the results
outputPath = paste0("./plots/", processingNote, "/")
dir.create(file.path(outputPath), showWarnings=FALSE)
plotUMAP_Monocle(endothCDS, paste0(processingNote, markerFileToUse), "cluster_ext_type", show_labels=FALSE, outputPath=outputPath)

# Takee a look at the distribution of the cell_type calls vs. clusters
propDF = plotGroupedProportions(endothCDS, paste0(processingNote, markerFileToUse), "partition_label", "cluster_ext_type", 
        pathToPlot=paste0("./plots/", processingNote, "/"), 
        widthToUse=1800)

plotUMAP_Monocle(endothCDS, paste0(processingNote, markerFileToUse), "partition_label", show_labels=FALSE, outputPath=outputPath)
plotUMAP_Monocle(endothCDS, paste0(processingNote, markerFileToUse), "cluster_label", show_labels=FALSE, outputPath=outputPath)

propDF = plotGroupedProportions(endothCDS, paste0(processingNote, markerFileToUse), "cluster_label", "cluster_ext_type", 
        pathToPlot=paste0("./plots/", processingNote, "/"), 
        widthToUse=1800)

hardAssignTentativeClustLabels = function(inputCDS, processingNote, oldProcNote){
	# Make sure using a processing note that I've hard-set
	if ((oldProcNote == "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40") & (processingNote == "HM10_Endoth_and_SMC")){
		# Assign
		colData(inputCDS)$hardCodedCellType = "Unknown"
		colData(inputCDS)$hardCodedCellType = ifelse(colData(inputCDS)$cluster_label == "11", 
						"Pericytes + Vascular Smooth Muscle", colData(inputCDS)$hardCodedCellType)
		colData(inputCDS)$hardCodedCellType = ifelse(colData(inputCDS)$cluster_label == "12", 
						"Venous Endothelium", colData(inputCDS)$hardCodedCellType)
		colData(inputCDS)$hardCodedCellType = ifelse(colData(inputCDS)$cluster_label == "13", 
						"Endocardium", colData(inputCDS)$hardCodedCellType)
		colData(inputCDS)$hardCodedCellType = ifelse(colData(inputCDS)$cluster_label == "16", 
						"Lymphatic + Venous", colData(inputCDS)$hardCodedCellType)
		colData(inputCDS)$hardCodedCellType = ifelse(colData(inputCDS)$cluster_label == "22", 
						"Lymphatic Endothelium", colData(inputCDS)$hardCodedCellType)
		colData(inputCDS)$hardCodedCellType = ifelse(colData(inputCDS)$cluster_label == "7", 
						"Arterial + Capillary", colData(inputCDS)$hardCodedCellType)
	}

	return(inputCDS)
}

endothCDS = hardAssignTentativeClustLabels(endothCDS, processingNote, oldProcNote)
plotUMAP_Monocle(endothCDS, paste0(processingNote, markerFileToUse), "hardCodedCellType", show_labels=TRUE, outputPath=outputPath, textSize=5)


#--------------------------------------------------------------------------------------------------------------------------------------
# Testing for differences by donor and anatomical site
#--------------------------------------------------------------------------------------------------------------------------------------


hardAssignAnatomicalSites <- function(inputCDS){
	colData(inputCDS)$Anatomical_Site = "Not_Specified"
	colData(inputCDS)$Anatomical_Site = ifelse(grepl("Apex", colData(inputCDS)$sample),
										"Apex", colData(inputCDS)$Anatomical_Site)
	colData(inputCDS)$Anatomical_Site = ifelse(grepl("Septum", colData(inputCDS)$sample),
										"Septum", colData(inputCDS)$Anatomical_Site)
	colData(inputCDS)$Anatomical_Site = ifelse(grepl("Left.Vent", colData(inputCDS)$sample),
										"Left_Vent", colData(inputCDS)$Anatomical_Site)
	colData(inputCDS)$Anatomical_Site = ifelse(grepl("Right.Vent", colData(inputCDS)$sample),
										"Right_Vent", colData(inputCDS)$Anatomical_Site)
	return(inputCDS)
}
hardAssignAnatomicalSitesDF <- function(inputDF){

	inputDF$Anatomical_Site = "Not_Specified"
	inputDF$Anatomical_Site = ifelse(grepl("Apex", inputDF$sample),
										"Apex", inputDF$Anatomical_Site)
	inputDF$Anatomical_Site = ifelse(grepl("Septum", inputDF$sample),
										"Septum", inputDF$Anatomical_Site)
	inputDF$Anatomical_Site = ifelse(grepl("Left.Vent", inputDF$sample),
										"Left_Vent", inputDF$Anatomical_Site)
	inputDF$Anatomical_Site = ifelse(grepl("Right.Vent", inputDF$sample),
										"Right_Vent", inputDF$Anatomical_Site)
	return(inputDF)
}
hardAssignDonors <- function(inputCDS){
	colData(inputCDS)$Donor = "Not_Specified"
	colData(inputCDS)$Donor = ifelse(grepl("W134", colData(inputCDS)$sample),
										"W134", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W135", colData(inputCDS)$sample),
										"W135", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W136", colData(inputCDS)$sample),
										"W136", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W139", colData(inputCDS)$sample),
										"W139", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W142", colData(inputCDS)$sample),
										"W142", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W144", colData(inputCDS)$sample),
										"W144", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W145", colData(inputCDS)$sample),
										"W145", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W146", colData(inputCDS)$sample),
										"W146", colData(inputCDS)$Donor)
	return(inputCDS)
}

colData(endothCDS)$sample = colData(endothCDS)$sampleName
endothCDS = hardAssignAnatomicalSites(endothCDS)
endothCDS = hardAssignDonors(endothCDS)


# Plot proportions of overall endothelium and sub-types. Split up to see variation by donor and subt-type

makeFacetedProportionPlot_withFill <- function(inputCDS, processingNote, colForGroup, colForProps, 
			fillCol = colForGroup, colorCol = colForGroup,
			subsetPropsToShow=levels(as.factor(colData(inputCDS)[[colForProps]]))){
	# browser()
	# Get the proportions
	countTable = with(colData(inputCDS), table(get(colForGroup), get(colForProps)))
	proportionTable = prop.table(countTable, margin=1)
	propDF = as.data.frame(proportionTable)

	# Get the alignment of colForGroup and fillCol
	countTable = with(colData(inputCDS), table(get(colForGroup), get(fillCol)))
	proportionTable = prop.table(countTable, margin=1)
	fillDF = as.data.frame(proportionTable)
	# Keep the ones that are synonymous
	fillDF = fillDF[fillDF$Freq == 1,]
	# Get the name vec
	groupToFillVec = as.character(fillDF$Var2)
	names(groupToFillVec) = as.character(fillDF$Var1)

########################
	# Get the alignment of colForGroup and colorCol
	countTable = with(colData(inputCDS), table(get(colForGroup), get(colorCol)))
	proportionTable = prop.table(countTable, margin=1)
	fillDF = as.data.frame(proportionTable)
	# Keep the ones that are synonymous
	fillDF = fillDF[fillDF$Freq == 1,]
	# Get the name vec
	groupToColorVec = as.character(fillDF$Var2)
	names(groupToColorVec) = as.character(fillDF$Var1)
##########################
#
	# Subset for only the subsets to show (by default this is all) and re-add the fill info
	colnames(propDF) = c(colForGroup, colForProps, "Proportion")
	propDF[[colForProps]] = as.character(propDF[[colForProps]])
	propDF[[colForGroup]] = as.character(propDF[[colForGroup]])
	propDF = propDF[propDF[[colForProps]] %in% subsetPropsToShow,]
	propDF[[fillCol]] = groupToFillVec[propDF[[colForGroup]]]
	propDF[[colorCol]] = groupToColorVec[propDF[[colForGroup]]]
	# browser()

	# Make the plot
	png(paste0("./plots/", processingNote, colForGroup, "propsBy", colForProps,
					"groupedBy", fillCol, "colorBy", colorCol, ".png"),
					width=1600, height=1200, res=200 )
	myPlot <- (ggplot(propDF, 
        aes_string(x=fillCol, y="Proportion", color=colorCol,
        						 fill=colForProps)) +
		geom_point() +
         facet_wrap(vars(get(colForProps)), scales="free") +
        # facet_wrap(vars(get(colForProps)), scales="free") +
      # geom_boxplot() + 
        # geom_jitter() + 
        # facet_wrap(vars(get(colForProps)), scales="free") +
        ylab("Proportion of Total") + #xlab("") +
        theme(text = element_text(size = 12))    +
         theme(axis.text.x = element_text(angle = 90, hjust = 1))  
         #  +
         # theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
         )
	print(myPlot)
	dev.off()

	return(myPlot)
}


hardAssignHighLevelCellTypes <- function(inputCDS){
	# Assign
	colData(inputCDS)$highLevelCellType = "Unassigned"
	colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label == "1", 
					"Fibroblast", colData(inputCDS)$highLevelCellType)
	colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("4", "7", "8", "13"), 
					"Endothelium", colData(inputCDS)$highLevelCellType)
	colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("6"), 
					"VSM_and_Pericyte", colData(inputCDS)$highLevelCellType)
	colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("2"), 
					"Macrophage", colData(inputCDS)$highLevelCellType)
	colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("3"), 
					"Cardiomyocyte", colData(inputCDS)$highLevelCellType)
	colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("5"), 
					"T_Cell", colData(inputCDS)$highLevelCellType)
	colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("11"), 
					"Adipocytes", colData(inputCDS)$highLevelCellType)
	colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("10"), 
					"Neuronal", colData(inputCDS)$highLevelCellType)
	

	return(inputCDS)
}

# Need a CDS with all the cells of this dataset, so I can estimate the fraction of all cells that are endothelium, etc.
allCellCDS = readRDS(paste0(rdsPath, "allCells_", oldProcNote, ".rds"))
allCellCDS = hardAssignHighLevelCellTypes(allCellCDS)

colData(allCellCDS)$sample = colData(allCellCDS)$sampleName
allCellCDS = hardAssignAnatomicalSites(allCellCDS)
allCellCDS = hardAssignDonors(allCellCDS)

allCellCDS = hardAssignTentativeClustLabels(allCellCDS, processingNote, oldProcNote)

# Overall endothelium, by donor
thisPlot = makeFacetedProportionPlot_withFill(allCellCDS, processingNote, "sample", "highLevelCellType", 
							fillCol = "Donor", subsetPropsToShow=c("Endothelium"), colorCol = "Anatomical_Site" )

endothCDS = hardAssignHighLevelCellTypes(endothCDS)

pureEndothCDS = endothCDS[,colData(endothCDS)$highLevelCellType=="Endothelium"]
thisPlot = makeFacetedProportionPlot_withFill(pureEndothCDS, processingNote, "sample", "hardCodedCellType", 
							fillCol = "Donor", colorCol = "Anatomical_Site" #, subsetPropsToShow=c("Endothelium")
							)











# -----------------------------------------------------------------------------------------
# DE testing by endothelial sub-types
# ---------------------------------------------------------------------------------------

# Within the endothelial group, I want to try running DE testing comparing each sub-type against the others.
#    Since I'm trying to find what differentiates each individual subtype from the rest, I need to split up and
#    Run a fit for each subtype individually, rather than a single fit using dummy variables.
#

#   That's a bit different than what I've tended to do before 
genesToTest = 1000
endotheliumMarkerTestList = runDEtestingToID_markers(pureEndothCDS, processingNote, "hardCodedCellType", howManyGenesToTest = genesToTest)

markerTestRes = endotheliumMarkerTestList[["marker_test_res"]]

# Save the results
write.csv(markerTestRes, paste0("./rdsOutput/", processingNote, "_endothelialSubtyping", as.character(genesToTest), "_genesPerType", ".csv"))



outputPath = paste0("./plots/", "OnlyEndothelium_", processingNote, "/")
dir.create(file.path(outputPath), showWarnings=FALSE)

arterialMarks = c("SEMA3G", "EFNB2", "DLL4")
venousMarks = c("NR2F2", "ACKR1")
endocardMarks = c("SMOC1", "NPR3")
lymphatMarks = c("PROX1", "TBX1", "PDPN")

# Also, just plot some marker gene sets that are known
plotUMAP_Monocle_genes(pureEndothCDS, processingNote, arterialMarks,
                     "Litv_Arterial_Markers", outputPath = outputPath, plotSetTotals=FALSE)
plotUMAP_Monocle_genes(pureEndothCDS, processingNote, venousMarks,
                     "Litv_Venous_Markers", outputPath = outputPath, plotSetTotals=FALSE)
plotUMAP_Monocle_genes(pureEndothCDS, processingNote, endocardMarks,
                     "Litv_Endocardial_Markers", outputPath = outputPath, plotSetTotals=FALSE)
plotUMAP_Monocle_genes(pureEndothCDS, processingNote, lymphatMarks,
                     "Litv_Lymphatic_Markers", outputPath = outputPath, plotSetTotals=FALSE)


allMarkersLitv = c(arterialMarks, venousMarks, endocardMarks, lymphatMarks)
onlyMarkerEndoCDS = pureEndothCDS[rowData(pureEndothCDS)$gene_short_name %in% allMarkersLitv,]

endotheliumMarkerTestList = runDEtestingToID_markers(onlyMarkerEndoCDS,
								 "litvKnownSubtypeMarkersOnly_", "hardCodedCellType")

markerTestRes = endotheliumMarkerTestList[["marker_test_res"]]

# Save the results
write.csv(markerTestRes, paste0("./rdsOutput/",  "litvKnownSubtypeMarkersOnly_", ".csv"))



