
library(ggplot2)
library(dplyr)
library("optparse")
library(stringr)
library(tidyr)
library(tidyverse)

print("Libraries loaded, starting now")


getCombinedCovariateCSV <- function(opt, inputCelltypes, inputCovariateList){
	# Loop and get a combined CSV for each covariate
	for (eachCovariate in names(inputCovariateList)){
		# Read the csv, combine into a larger on
		myDF = data.frame()
		for (eachType in inputCelltypes){
			# inputDEfile = paste0("./plots/", eachType, opt$modelNotes, "/", eachType, "_", eachCovariate, "_AllTestsRun_Table.csv")
			inputDEfile = paste0("./plots/", eachType, "/", eachType, "_", eachCovariate, "_AllTestsRun_Table.csv")
			thisCSV = read.csv(inputDEfile)
			# Keep only the pathway name, enrichment, and p val. 
			# browser()
			thisCSV = thisCSV[c("gene", "coefficientValue", "q_val")]
			# Add a column for the cell type, for tidy format plotting later
			thisCSV$cellType = eachType
			myDF = rbind(myDF, thisCSV)
		}

		# Add columns for absolute value of effect size and direction of effect (for plotting)
		myDF$Effect_Magnitude = abs(myDF$coefficientValue)
		myDF$Effect_Direction = ifelse(myDF$coefficientValue > 0.0, "Positive", "Negative")

		myDF = myDF[myDF$q_val < opt$padjCutoff,]

		inputCovariateList[[eachCovariate]] = myDF
	}

	# Return the whole csv
	return(inputCovariateList)
}


getCombinedCelltypeCSV <- function(opt, inputCelltypes, inputPath = "./plots/cellTypeSpec/Cell_Type_Specificity_for"){
	# Read the csv, combine into a larger on
	myDF = data.frame()
	for (eachType in inputCelltypes){
		# inputDEfile = paste0(inputPath, eachType, opt$modelNotes, "/", eachType, "_is_", eachType, "_AllTestsRun_Table.csv")
		inputDEfile = paste0(inputPath, eachType, "/", eachType, "_is_", eachType, "_AllTestsRun_Table.csv")
		thisCSV = read.csv(inputDEfile)
		# Keep only the pathway name, enrichment, and p val. 
		# browser()
		thisCSV = thisCSV[c("gene", "coefficientValue", "q_val")]
		# Add a column for the cell type, for tidy format plotting later
		thisCSV$cellType = eachType
		myDF = rbind(myDF, thisCSV)
	}

	# Add columns for absolute value of effect size and direction of effect (for plotting)
	myDF$Effect_Magnitude = abs(myDF$coefficientValue)
	myDF$Effect_Direction = ifelse(myDF$coefficientValue > 0.0, "Positive", "Negative")

	myDF = myDF[myDF$q_val < opt$padjCutoff,]

	# Return the whole csv
	return(myDF)
}




getDEcounts <- function(inputCSVlist, inputCovariates){
	countList = vector(mode='list', length = length(inputCovariates))
	names(countList) = inputCovariates

	# Group each by cell type
	for (eachCovariate in inputCovariates){
		countList[[eachCovariate]] = (inputCSVlist[[eachCovariate]] %>% count(cellType))
		countList[[eachCovariate]]$covariate = eachCovariate
	}

	combinedDF = bind_rows(countList)
}




plotCellMarkers <- function(inputCSV, plotNote, markersToPlot, cellTypesToPlot, opt, onlyPosCoefs = TRUE){
	# If only plotting enrichments, throw out all negative coefficient entries
	if (onlyPosCoefs){
		inputCSV = inputCSV[inputCSV$coefficientValue > 0,]
	}

	# Get the subset for entries corresponding to marker genes
	inputCSV = inputCSV[inputCSV$gene %in% markersToPlot,]

	# Set order
	inputCSV$cellType = factor(inputCSV$cellType, levels=cellTypesToPlot)
	inputCSV$gene = factor(inputCSV$gene, levels=markersToPlot)

	# Plot these markers
	outputDir = "./plots/atacMarkers/"
	dir.create(outputDir)
	outputFile = paste0(outputDir, plotNote, opt$modelNotes, ".png")
	png(outputFile, res=200, width = 1400, height=1200)
	myPlot = ggplot(inputCSV, aes_string(x="cellType", y="gene", size="coefficientValue", color="coefficientValue	")) + 
					geom_point() + 
					theme(axis.text.x = element_text(angle = 45, hjust=1))+
            theme(text = element_text(size = 20)) + xlab("Cell Type") + 
            	ylab("TF") + 
            	 guides(size=guide_legend(title="Enrichment"))+ 
            	 guides(color=FALSE)

	print(myPlot)
	dev.off()

}

formatCellType <- function(inputColumn){

	inputColumn = ifelse(inputColumn == "T_Cell", "T Cell", inputColumn)
	inputColumn = ifelse(inputColumn == "VSM_and_Pericyte", "Perviascular Cell", inputColumn)
	inputColumn = ifelse(inputColumn == "Vascular_Endothelium", "Vascular Endothelium", inputColumn)
	inputColumn = ifelse(inputColumn == "Lymphatic_Endothelium", "Lymphatic Endothelium", inputColumn)
	inputColumn = ifelse(inputColumn == "Mast_Cell", "Mast Cell", inputColumn)

	return(inputColumn)
}




# Get the passed parameters
option_list = list(
  make_option(c("-m", "--modelNotes"), type="character", 
  			default="_fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult",
              help="Processing note from model fitting", metavar="character"),
  # make_option(c("-b", "--covariate"), type="character", 
  # 			default="SexM",   #"Age", # "SexM"
  #             help="Covariate to plot GSEA summary plot", metavar="character"),
  make_option(c("-p", "--padjCutoff"), type="numeric", 
  			default=0.1,
              help="Max padj value to plot as significant", metavar="numeric"),
  make_option(c("-c", "--cellType"), type="character", 
  			default="Vascular_Endothelium",
              help="Cell type for which the model was fit", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cellTypes = c("Vascular_Endothelium", "Cardiomyocyte", "Macrophage", "T_Cell", "VSM_and_Pericyte", "Fibroblast", "Endocardium", "Mast_Cell", "Adipocytes", "Lymphatic_Endothelium")

covariatesToPlot = c("SexM", "Age")
covariateCSVs = vector(mode='list', length=length(covariatesToPlot))
names(covariateCSVs) = covariatesToPlot

combinedCSVlist = getCombinedCovariateCSV(opt, cellTypes, covariateCSVs)


combinedCellTypeCSV = getCombinedCelltypeCSV(opt, cellTypes)
combinedCellTypeCSV$cellType = formatCellType(combinedCellTypeCSV$cellType)


markerSet = c("Sox17", "FOSL1::JUN(var.2)", "MEF2A", "MEF2B", "SPI1", "SPIC", "RUNX3", "RXRA::VDR",  "ETV2", "CEBPA", "EBF1", "SOX9", "MITF", "Cebpa", "Ddit3::Cebpa", "E2F7", "E2F8" )
plotCellMarkers(combinedCellTypeCSV, "canonicalMarkers", markerSet, formatCellType(cellTypes), opt)








# Counts of the DE genes
deCounts = getDEcounts(combinedCSVlist, covariatesToPlot)



# Simplest version: only plot those that 
outDir = "./plots/DE_Summaries/"
dir.create(outDir)



# Output:
outfile = paste0("DE_hits_qVal_", as.character(opt$padjCutoff), "_by_", paste0(covariatesToPlot, collapse="_"), ".png" )
png(paste0(outDir, outfile), res=200, width=1200, height=1000)
myPlot = ggplot(deCounts, aes_string(x="cellType", y="n", fill="covariate")) +
			geom_bar(position="dodge", stat="identity") +
			 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(myPlot)
dev.off()


# Output a CSV holding these hits
for (eachCovariate in covariatesToPlot){
	fullDEfile = paste0(outDir, "Combined_DE_by_", eachCovariate, "_q_", as.character(opt$padjCutoff), ".csv")
	write.csv(combinedCSVlist[[eachCovariate]], file=fullDEfile)

} 


