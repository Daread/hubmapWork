
library(tidyr)
library(ggplot2)
library(glmnet)
library(ROCR)
library(data.table)
library(monocle3)
library(gridExtra)

library(RColorBrewer)

# Get shared functions between this and runRegressionModel.R
source("./utilityFunctions.R")

source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

library("optparse")
# Get the passed parameters
option_list = list(
   make_option(c("-c", "--cdsRNA"), type="character", 
        default="allCells_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes.rds", 
              help="Name of RNA cds to use", metavar="character"),
  make_option(c("-r", "--rnaPath"), type="character", 
        default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_19_largeScaleMixedModeling_nobackup/plots/", 
              help="Path to RNA cds", metavar="character"),

  make_option(c("-p", "--pvalCutoff"), type="numeric", 
        default=0.1, 
              help="Path to RNA cds", metavar="numeric"),

  make_option(c("-n", "--atacNotes"), type="character", 
        default="_fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult", 
              help="Name of RNA cds to use", metavar="character"),
  make_option(c("-a", "--atacPath"), type="character", 
        default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_05_redoPeakMotifRegression_nobackup/plots/", 
              help="Path to RNA cds", metavar="character")

)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

getATACcoefs <- function(opt, cellTypes, covariatesToGet=c("SexM", "Age")){
	# Set up the list to store the DFs
	covariateList = vector(mode="list", length = length(covariatesToGet))
	names(covariateList) = covariatesToGet
	# Loop and retrieve each covariate's DF for each cell type desired
	for (eachCovar in covariatesToGet){
		coefDF = data.frame()
		# Read in the csv for each 
		for (eachType in cellTypes){
			# thisFile = paste0(opt$atacPath, eachType, opt$atacNotes, "/", eachType, "_", eachCovar, "_AllTestsRun_Table.csv")
			thisFile = paste0(opt$atacPath, eachType, "/", eachType, "_", eachCovar, "_AllTestsRun_Table.csv")
			thisDF = read.csv(thisFile)
			thisDF$CellType = eachType
			coefDF = rbind(coefDF, thisDF)
			# Set order of cell types as a factor
			# coefDF$CellType = factor(coefDF$CellType, levels = cellTypes)

		}
		# Save this combined df into the list slot for this coefficient
		covariateList[[eachCovar]] = coefDF
	}
	
	return(covariateList)
}


checkForZeroHitCelltypes <- function(plotDF){
	# Check by each type, see if sum of coefficients is zero. If so, remove
	possibleCellTypes = levels(as.factor(plotDF$CellType))
	typesToKeep = character()
	for (eachType in possibleCellTypes){
		coefSum = sum(abs(plotDF[plotDF$CellType == eachType,][["coefficientValue"]]))
		if (coefSum > 0){
			typesToKeep = c(typesToKeep, eachType)
		}
	}
	# browser()

	plotDF = plotDF[plotDF$CellType %in% typesToKeep,]
	return(plotDF)
}


plotCoefBars <- function(coefList, cellTypes, opt, plotCovar="Age", geneSet = c("MEF2C"), 
				setLabel = "JustMEF2C", cutoffByP = TRUE){
	# Get the df for the proper coef
	plotDF = coefList[[plotCovar]]
	# Go down to the genes to be plotted
	plotDF = plotDF[plotDF$gene %in% geneSet,]

	plotDF$CellType = factor(plotDF$CellType, levels=formatCellType(cellTypes))

	# Set the custom pallete
	cellColors = setNames(brewer.pal(length(cellTypes), "Set1"), levels(plotDF$CellType))

	# Trim by p value?
	if (cutoffByP){
		# Changed on 2-22-22. Just set coefficients to zero if not significant. Makes plotting behave better
		# plotDF = plotDF[plotDF$q_val < opt$pvalCutoff,]
		plotDF$coefficientValue = ifelse(plotDF$q_val < opt$pvalCutoff, plotDF$coefficientValue, 0)
		# If one cell type doesn't have any significant results, remove it entirely from the plotting DF
		plotDF = checkForZeroHitCelltypes(plotDF)
		plotTitle = paste0(plotCovar, " Coefficients, q < ", as.character(opt$pvalCutoff))
	} else {
		plotTitle = paste0(plotCovar, " Coefficients")
	}

	# covarName = ifelse(plotCovar == "SexM", "Sex", plotCovar)
	yLabToUse = ifelse(plotCovar == "SexM", "Log Motif Enrichment, Male/Female",
										"Enrichment with Age")

	# Change TFs into a factor so that order is specified
	plotDF$gene = factor(plotDF$gene, levels=geneSet)

	# Plot the coefficients
	outDir = paste0("./plots/MM_Coef_Plots/")
	dir.create(outDir)
	png(paste0(outDir, setLabel, "_Coefficients.png"), 
				res=200, height=1200, width = 1400)
	myPlot = ggplot(plotDF, aes_string(x="gene", y="coefficientValue", fill="CellType")) + 
			# geom_bar(position="dodge", stat="identity") + 
			geom_col(width=0.5, position=position_dodge(0.5)) +
			#ggtitle(plotTitle) + 
			 theme(text = element_text(size = 20))  + #ylab(paste0("Coefficient for ", covarName)) +
			 ylab(yLabToUse) +
			 xlab("TF") + 
			 guides(fill=guide_legend(title="Cell Type")) +
			 theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  				# scale_fill_brewer(palette = "Dark2")
  				scale_fill_manual(values=cellColors, drop = TRUE, limits=force)
	print(myPlot)
	dev.off()
}

formatCellType <- function(inputColumn){
	inputColumn = ifelse(inputColumn == "T_Cell", "T Cell", inputColumn)
	inputColumn = ifelse(inputColumn == "VSM_and_Pericyte", "Perviascular Cell", inputColumn)
	inputColumn = ifelse(inputColumn == "Vascular_Endothelium", "Vascular Endothelium", inputColumn)

	return(inputColumn)
}

# cellTypes= c("Macrophage", "Cardiomyocyte", "Endocardium", "Fibroblast", "T_Cell", "Vascular_Endothelium")
cellTypes = c("Vascular_Endothelium", "Cardiomyocyte", "Macrophage", "Fibroblast", "T_Cell", "VSM_and_Pericyte",  "Endocardium", "Mast_Cell") #, "Adipocytes", "Lymphatic_Endothelium", "B_Cell", "Neuronal")


# rnaResults = getRNAcoefs(opt)
atacResults = getATACcoefs(opt, cellTypes)

# Format the celltype names
for (eachCovar in names(atacResults)){
	atacResults[[eachCovar]][["CellType"]] = formatCellType(atacResults[[eachCovar]][["CellType"]])
}

# Plot desired coefficients
plotCoefBars(atacResults, cellTypes, opt, plotCovar="Age", geneSet = c("IRF7", "ATF7", "STAT1"), 
				setLabel = "AgeTestSet")


# plotCoefBars(atacResults, cellTypes, opt, plotCovar="SexM", geneSet = c("Smad4", "JDP2", "SNAI2", "HIF1A", "ARNT::HIF1A", "RORA"), 
# 				setLabel = "Sex_TGFB_and_Metab_Regs")

plotCoefBars(atacResults, cellTypes, opt, plotCovar="SexM", geneSet = c("SMAD3", "JDP2",  "HIF1A", "MYC", "RORA"), 
				setLabel = "Sex_TGFB_and_Metab_Regs")

# plotCoefBars(atacResults, cellTypes, opt, plotCovar="Age", geneSet = c("IRF1", "IRF7", "NFKB1", "STAT1", "STAT2", "CUX1"), 
# 				setLabel = "Age_Immunity_and_Senesc")

# c("SMAD3", "JDP2",  "HIF1A", "ARNT::HIF1A", "RORA"), 



# plotCoefBars(atacResults, cellTypes, opt, plotCovar="Age", geneSet = c("IRF1", "IRF7", "NFKB1", "CUX1", "JUNB", "JUND"), 
# 				setLabel = "Age_Immunity_and_Senesc")


plotCoefBars(atacResults, cellTypes, opt, plotCovar="Age", geneSet = c("IRF1", "IRF7", "NFKB2", "FOS::JUN(var.2)", "CUX1", "JUND"), 
				setLabel = "Age_Immunity_and_Senesc")












