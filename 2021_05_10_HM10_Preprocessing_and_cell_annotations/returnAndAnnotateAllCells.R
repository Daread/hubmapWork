

###########################################################################################
# Load relevant packages
source("../../../sharedProjectCode/utility/singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()
print(modStatus)
source("../../../sharedProjectCode/scRNA_Seq_Background_Code/backgroundMethods.R")


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

hardAssignHighLevelCellTypes <- function(inputCDS, origProcessing){
	if (origProcessing == "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40"){
		# Assign
		colData(inputCDS)$highLevelCellType = "Unassigned"
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label == "1", 
						"Fibroblast", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("4"), 
						"Vascular_Endothelium", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("7"), 
						"Endocardium", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("8", "13"), 
						"Lymphatic_Endothelium", colData(inputCDS)$highLevelCellType)
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
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("12", "14"), 
						"B_Cell", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("9"), 
						"Mast_Cell", colData(inputCDS)$highLevelCellType)
	}
	
	return(inputCDS)
}





# 2021-07-16 Work: Returning to sub-type the last few clusters
# hard set processing note, and read in the CDS by a specific name...?


# Load the saved file in
rdsPath = "../2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/"
oldProcNote = "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40"
allCellCDS = readRDS(paste0(rdsPath, "allCells_", oldProcNote, ".rds"))

processingNote = paste0(oldProcNote, "addAllTypes")
outputPath = paste0("./plots/", oldProcNote, "/")

myMarkerRes = runDEtestingToID_markers(allCellCDS, processingNote, "partition_label",
									howManyGenesToTest = 25, outputPath=outputPath)

# From this, it looks like 14 and 12 are both B cell subsets, while 9 might be mast cells
colData(allCellCDS)$sample = colData(allCellCDS)$sampleName
allCellCDS = hardAssignHighLevelCellTypes(allCellCDS, oldProcNote)
allCellCDS = hardAssignDonors(allCellCDS)
allCellCDS = hardAssignAnatomicalSites(allCellCDS)


# Save
saveRDS(allCellCDS, paste0("./rdsOutput/allCells_", processingNote, ".rds"))


# Make a plot of the cell types
plotOutput = paste0("./plots/", oldProcNote, "/")
plotUMAP_Monocle(allCellCDS, processingNote, "highLevelCellType", outputPath=plotOutput, textSize=4)


# Make a plot of the cell types
plotOutput = paste0("./plots/", oldProcNote, "/")
plotUMAP_Monocle(allCellCDS, processingNote, "partition_label", outputPath=plotOutput, textSize=4)









