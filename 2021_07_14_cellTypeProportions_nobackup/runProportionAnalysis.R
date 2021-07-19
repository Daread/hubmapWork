

###########################################################################################
# Load relevant packages
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()
print(modStatus)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/scRNA_Seq_Background_Code/backgroundMethods.R")



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
	}
	
	

	return(inputCDS)
}


# Get the passed parameters
library("optparse")

option_list = list(
  make_option(c("-n", "--name"), type="character", default="TestRun", 
              help="Name of this subset analysis", metavar="character"),
    make_option(c("-c", "--colToUse"), type="character", default="highLevelCellType", 
              help="Column for which to run proportion tests", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

processingNote = paste0(opt$name, "_by_", opt$colToUse)

# With arguments, read in data
rdsPath = "../2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/"
oldProcNote = "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40"
allCellCDS = readRDS(paste0(rdsPath, "allCells_", oldProcNote, ".rds"))

# Assign the subtypes
colData(allCellCDS)$sample = colData(allCellCDS)$sampleName
allCellCDS = hardAssignHighLevelCellTypes(allCellCDS, oldProcNote)
allCellCDS = hardAssignDonors(allCellCDS)
allCellCDS = hardAssignAnatomicalSites(allCellCDS)

# Based off of this original CDS name, make a new directory for outputs
oldCDS_Path = paste0("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_14_cellTypeProportions_nobackup/plots/", oldProcNote, "/")
dir.create(file.path(oldCDS_Path), showWarnings=FALSE)
outputPath = paste0(oldCDS_Path, processingNote, "/")
dir.create(file.path(outputPath))


set.seed(7)
# Loop over each unique value in the colToUse
for (eachGroup in as.character(levels(as.factor(colData(allCellCDS)[[opt$colToUse]])))){
	# Plot proportions of this
	makeFacetedProportionPlot_withFill(allCellCDS, paste0(processingNote, "_", eachGroup), "sample", opt$colToUse, 
            fillCol = "Donor", subsetPropsToShow=c(eachGroup), colorCol="Anatomical_Site",
            outputPath=outputPath )
	
}

library(reshape2)
library("dplyr")
# Next, look and see if there are obvious trends in cell composition
# cellPropDF = (as.data.frame(colData(allCellCDS)) %>% )
countTable = with(as.data.frame(colData(allCellCDS)), table(get("sample"), get(opt$colToUse)))
propTable = prop.table(countTable, margin=1)
propDF = as.data.frame(propTable)
unmeltedPropDF = dcast(propDF, Var2 ~ Var1)
colnames(unmeltedPropDF)[1] = opt$colToUse
rownames(unmeltedPropDF) = unmeltedPropDF[[opt$colToUse]]
unmeltedPropDF = unmeltedPropDF[, -which(names(unmeltedPropDF) %in% c(opt$colToUse))]


# From https://r-coder.com/correlation-plot-r/

# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    Cor <- (cor(x, y)) # Remove abs function if desired
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
    if(missing(cex.cor)) {
        cex.cor <- 0.4 / strwidth(txt)
    }
    text(0.5, 0.5, txt,
         cex = 1 + cex.cor * abs(Cor)) # Resize the text by level of correlation
}
######


png(paste0(outputPath, "scatterPlotMatrix_", opt$colToUse, ".png"),
				width=2000, height=2000, res=200)
thisPlot = pairs(t(unmeltedPropDF),
				upper.panel = panel.cor)
print(thisPlot)
dev.off()


# Just make a plot to show the UMAP
plotUMAP_Monocle(allCellCDS, processingNote, "highLevelCellType", textSize=4,
				show_labels=TRUE, outputPath=outputPath)















