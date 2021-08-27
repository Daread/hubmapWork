

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
# oldProcNote = "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40"
oldProcNote = "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes"
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




# 8-17-21: Get a dataframe with proportions of each cell type for each sample, labeling donor and site.
#.         Then use this to look at inter-site vs inter-donor variation, as well as 

library(data.table)

cellTypePropDF = transpose(unmeltedPropDF)
rownames(cellTypePropDF) = colnames(unmeltedPropDF)
colnames(cellTypePropDF) = rownames(unmeltedPropDF)

cellTypePropDF$sampleName = rownames(cellTypePropDF)

meltedPropDF = reshape::melt(cellTypePropDF, id="sampleName")
colnames(meltedPropDF) = c("sampleName", "Cell_Type", "Proportion")


library(tidyr)
meltedPropDF = tidyr::separate(meltedPropDF, col="sampleName", into=c("Donor", "Anatomical_Site"), 
				sep="[.]", remove=FALSE, extra="merge")


# Now, time to get 
varianceComparisonDF = data.frame("Mean_Value" = double(),
								"Comparison" = character(),
								  "Cell_Type" = character())

cellTypes = c("Adipocytes", "B_Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic_Endothelium", "Macrophage", "Mast_Cell", "Neuronal", "T_Cell",
                "Vascular_Endothelium", "VSM_and_Pericyte")



getInterDonorSetSite <- function(subsetDF, siteToUse){
	thisDF = subsetDF[subsetDF$Anatomical_Site == siteToUse,]
	# Loop over all combinations
	comboMeanVec = numeric()

	for (eachInd in 1:(nrow(thisDF) - 1)){
		for (eachComp in (eachInd + 1):nrow(thisDF)){
			thisAbsDiff = abs(thisDF[eachInd, "Proportion"] - thisDF[eachComp, "Proportion"] )

			comboMeanVec = c(comboMeanVec, thisAbsDiff)
		}
	}
	return ( mean(comboMeanVec))
}


getIntraDonorDiff <- function(subsetDF, siteOne, siteTwo){
	thisDF = subsetDF[subsetDF$Anatomical_Site %in% c(siteOne, siteTwo),]

	interSiteVec = numeric()
	# Get all the intra donor combos
	for (eachDonor in levels(as.factor(subsetDF$Donor))){
		miniDF = thisDF[thisDF$Donor == eachDonor,]
		# If only one, skip
		if (nrow(miniDF) == 1){
			next 
		}
		# Otherwise, get the comparison
		interSiteDiff = abs(miniDF[1, "Proportion"] - miniDF[2, "Proportion"])
		interSiteVec = c(interSiteVec, interSiteDiff)
	}

	return(mean(interSiteVec))
}



for (eachCellType in cellTypes){
	subsetDF = meltedPropDF[meltedPropDF$Cell_Type == eachCellType,]
	interLV = getInterDonorSetSite(subsetDF, "Left.Vent")
	interApex = getInterDonorSetSite(subsetDF, "Apex")

	interSiteIntraDonor = getIntraDonorDiff(subsetDF, "Left.Vent", "Apex")

	cellTypeDF = data.frame("Mean_Value" = c(interLV, interApex, interSiteIntraDonor))
	cellTypeDF$Comparison = c("Inter_Donor_LV", "Inter_Donor_Apex", "Inter_Site_Within_Donor")
	cellTypeDF$Cell_Type = eachCellType
	varianceComparisonDF = rbind(varianceComparisonDF, cellTypeDF)

}

# Make a plot


png(paste0(outputPath, "Compare_InterSite_vs_InterDonor_Differences", ".png"),
				width=2000, height=2000, res=200)
thisPlot = ggplot(varianceComparisonDF, aes_string(x="Cell_Type", y="Mean_Value", color="Comparison")) +
				geom_point()+ 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(thisPlot)
dev.off()



# Calculate the t test for LV vs apex proportions
tTestDF = data.frame("Cell_Type" = cellTypes)
rownames(tTestDF) = tTestDF$Cell_Type

for (eachCellType in cellTypes){
	subsetDF = meltedPropDF[meltedPropDF$Cell_Type == eachCellType,]
	subsetDF = subsetDF[subsetDF$Anatomical_Site %in% c("Left.Vent", "Apex"),]

	# Get the t test
	xVec = subsetDF[subsetDF$Anatomical_Site == "Left.Vent",]$Proportion
	yVec = subsetDF[subsetDF$Anatomical_Site == "Apex",]$Proportion
	tRes = t.test(xVec, yVec)

	tTestDF[eachCellType, "P_Value_LV_vs_Apex"] = tRes$p.value

}

# Save this

write.csv(tTestDF, paste0(outputPath, "tTestsForProportions.csv"))

























