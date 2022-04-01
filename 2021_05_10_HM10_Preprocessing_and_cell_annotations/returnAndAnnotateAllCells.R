

###########################################################################################
# Load relevant packages
# source("../../../sharedProjectCode/utility/singleCellUtilFuncs.R")
# modStatus <- loadMonoclePackages()
# print(modStatus)
# source("../../../sharedProjectCode/scRNA_Seq_Background_Code/backgroundMethods.R")


# load requirements
suppressPackageStartupMessages({
  # library(ArchR)
  library(monocle3)
  library(dplyr)
  library(ggrastr)
  library(ggplot2)
  # library(Seurat)
})
source("../../../sharedProjectCode/utility/singleCellUtilFuncs.R")
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


# 10-4-21 mod: Change the clusters field by hard-coding in partitions so that plots look correct


hardAssignDonorAges <- function(inputCDS){
  colData(inputCDS)$Age = 0
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W134", 43, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W135", 60, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W136", 43, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W137", 49, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W139", 45, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W142", 55, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W144", 53, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W145", 51, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W146", 25, colData(inputCDS)$Age)

  colData(inputCDS)$Log10Age = log10(colData(inputCDS)$Age)

  return(inputCDS)
}


hardAssignDonorSexes <- function(inputCDS){
  colData(inputCDS)$Sex = "Not_Set"
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W134", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W135", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W136", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W137", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W139", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W142", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W144", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W145", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W146", "F", colData(inputCDS)$Sex)

  return(inputCDS)
}



hardAssignHighLevelCellTypes <- function(inputCDS, origProcessing){
	if (origProcessing == "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40"){
		# Assign
		colData(inputCDS)$highLevelCellType = "Unassigned"
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label == "1", 
						"Fibroblast", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("4"), 
						"Vascular Endothelium", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("7"), 
						"Endocardium", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("8", "13"), 
						"Lymphatic Endothelium", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("6"), 
						"Perivascular Cell", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("2"), 
						"Macrophage", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("3"), 
						"Cardiomyocyte", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("5"), 
						"T Cell", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("11"), 
						"Adipocyte", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("10"), 
						"Neuronal", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("12", "14"), 
						"B Cell", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("9"), 
						"Mast Cell", colData(inputCDS)$highLevelCellType)
	}
	
	return(inputCDS)
}


allCellCDS = hardAssignHighLevelCellTypes(allCellCDS, oldProcNote)
# monocle3::clusters(allCellCDS) = partitions(allCellCDS)


plotUMAP_Monocle_modded <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, textSize=10, outputPath = "./plots/",
                    returnPlotObj=FALSE, inUMAPscale=0.5){ #, xVal, yVal){
    png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored.png"),
             width=1400, height=1000, res=200)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels, group_cells_by="partition",
          cell_stroke=.1 , group_label_size=(textSize * inUMAPscale)        
                ))
    myPlot = (myPlot + theme(text=element_text(size=textSize)))
    print(myPlot)
    dev.off()   

    if (returnPlotObj){
        return(myPlot)
    }

}



plotUMAP_Monocle_formatted <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, textSize=10, outputPath = "./plots/",
                    returnPlotObj=FALSE){ #, xVal, yVal){

	colData(dataCDS)[[catToColor]] = gsub("[.]", " ", colData(dataCDS)[[catToColor]])
	# Shuffle data to not overplot with the last subset
	dataCDS = dataCDS[,sample(ncol(dataCDS))]

    png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored_formatted.png"),
             width=1400, height=1000, res=200)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels,
          cell_stroke=.1 , group_label_size=textSize        
                ))
    myPlot = (myPlot + theme(text=element_text(size=textSize)))
    print(myPlot)
    dev.off()   

    if (returnPlotObj){
        return(myPlot)
    }

}




# Make a plot of the cell types
plotOutput = paste0("./plots/", oldProcNote, "/")
plotUMAP_Monocle_modded(allCellCDS, paste0(processingNote, "_LabByPartition"), "highLevelCellType", outputPath=plotOutput, textSize=12)



# Make a plot of the cell types
plotOutput = paste0("./plots/", oldProcNote, "/")
plotUMAP_Monocle(allCellCDS, processingNote, "partition_label", outputPath=plotOutput, textSize=4)






# 2-24-22: Formatting for publication versions
colData(allCellCDS)$Sample = colData(allCellCDS)$sampleName
plotUMAP_Monocle_formatted(allCellCDS, processingNote, "Sample", outputPath=plotOutput, textSize=24, show_labels=FALSE)

# Make a plot of the cell types
plotOutput = paste0("./plots/", oldProcNote, "/")
plotUMAP_Monocle_modded(allCellCDS, paste0(processingNote, "_LabByPartition_Formattted"), "highLevelCellType", outputPath=plotOutput, textSize=24, inUMAPscale=.3)





# 10-15-21: Return and make plots showing representation of different sexes and donors for each cell type
# Shuffle the CDS
allCellCDS = allCellCDS[,sample(ncol(allCellCDS))]

allCellCDS = hardAssignDonorSexes(allCellCDS)

# Now plot by sample and sex
plotUMAP_Monocle(allCellCDS, processingNote, "Sex", outputPath=plotOutput, textSize=12, show_labels=FALSE)

colData(allCellCDS)$Sample = colData(allCellCDS)$sampleName
plotUMAP_Monocle(allCellCDS, processingNote, "Sample", outputPath=plotOutput, textSize=12, show_labels=FALSE)


sexVec = colData(allCellCDS)$Sex 
sum(sexVec == "F")
sum(sexVec == "M")




















