
library(ggplot2)
library(dplyr)
library("optparse")
library(stringr)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

print("Libraries loaded, starting now")


monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}



formatCellType <- function(inputColumn){

	inputColumn = ifelse(inputColumn == "T_Cell", "T Cell", inputColumn)
	inputColumn = ifelse(inputColumn == "VSM_and_Pericyte", "Perivascular Cell", inputColumn)
	inputColumn = ifelse(inputColumn == "Vascular_Endothelium", "Vascular Endothelium", inputColumn)
	inputColumn = ifelse(inputColumn == "Lymphatic_Endothelium", "Lymphatic Endothelium", inputColumn)
	inputColumn = ifelse(inputColumn == "Mast_Cell", "Mast Cell", inputColumn)

	inputColumn = gsub("_", " ", inputColumn)

	return(inputColumn)
}


getCombinedCSV <- function(opt, inputCelltypes, inputCovariateList, filterbyQval = TRUE){

	# Loop and get a combined CSV for each covariate
	for (eachCovariate in names(inputCovariateList)){
		# Read the csv, combine into a larger on
		myDF = data.frame()
		for (eachType in inputCelltypes){
			inputDEfile = paste0("./plots/", eachType, opt$modelNotes, "/", eachType, "_", eachCovariate, "_AllTestsRun_Table.csv")
			thisCSV = read.csv(inputDEfile)
			# Keep only the pathway name, enrichment, and p val. 
			# browser()
			thisCSV = thisCSV[c("gene", "coefficientValue", "q_val", "pval")]
			# Add a column for the cell type, for tidy format plotting later
			thisCSV$cellType = eachType
			myDF = rbind(myDF, thisCSV)
		}

		# Add columns for absolute value of effect size and direction of effect (for plotting)
		myDF$Effect_Magnitude = abs(myDF$coefficientValue)
		myDF$Effect_Direction = ifelse(myDF$coefficientValue > 0.0, "Positive", "Negative")

		if (filterbyQval){
			myDF = myDF[myDF$q_val < opt$padjCutoff,]
		}
		
		inputCovariateList[[eachCovariate]] = myDF
	}

	
	# Return the whole csv
	return(inputCovariateList)
}


outputFullCSV = function(inputList){
	outDir = "./plots/DE_Summaries/"
	covarToOutput = c("Age", "SexM")
	for (eachCovar in covarToOutput){
		write.csv(inputList[[eachCovar]], file=paste0(outDir, eachCovar, "_Combined_DE_All.csv"))
	}
}


# Get the passed parameters
option_list = list(
  make_option(c("-m", "--modelNotes"), type="character", 
  			default="_fix_Anatomical_Site,Age,Sex,DataSource,log10_umi_rand_Donor_NucleiOnlySharedGenesCDS_noAtr_",
              help="Processing note from model fitting", metavar="character"),
  # make_option(c("-b", "--covariate"), type="character", 
  # 			default="SexM",   #"Age", # "SexM"
  #             help="Covariate to plot GSEA summary plot", metavar="character"),
  make_option(c("-p", "--padjCutoff"), type="numeric", 
  			default=0.1,
              help="Max padj value to plot as significant", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# cellTypes = c("Vascular_Endothelium", "Cardiomyocyte", "Macrophage", "T_Cell", "VSM_and_Pericyte", "Fibroblast", "Endocardium")
cellTypes = c("Endothelium", "Fibroblast", "Lymphocyte", "Myeloid", "Perivascular", "Ventricular_Cardiomyocytes", "Neuron", "Adipocytes")


covariatesToPlot = c("SexM", "Age")
# covariatesToPlot = c("SexM", "Age", "Anatomical_SiteApex", 
# 					"Anatomical_SiteRV", "Anatomical_SiteSeptum") 

covariateCSVs = vector(mode='list', length=length(covariatesToPlot))
names(covariateCSVs) = covariatesToPlot

combinedCSVlist = getCombinedCSV(opt, cellTypes, covariateCSVs)



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


# Counts of the DE genes
deCounts = getDEcounts(combinedCSVlist, covariatesToPlot)

deCounts$cellType = formatCellType(deCounts$cellType)
# Format the coviarate column
deCounts$Covariate = ifelse(deCounts$covariate == "SexM", "Sex", deCounts$covariate)
deCounts$Covariate = gsub("Anatomical_Site", "", deCounts$covariate) #ifelse(deCounts$covariate == "SexM", "Sex", deCounts$covariate)

# Simplest version: only plot those that 
outDir = "./plots/DE_Summaries/"
dir.create(outDir)


deCounts$cellType = formatCellType(deCounts$cellType)

# Output:
outfile = paste0("DE_hits_qVal_", as.character(opt$padjCutoff), "_by_", paste0(covariatesToPlot, collapse="_"), ".png" )
png(paste0(outDir, outfile), res=200, width=1200, height=1000)
# png(paste0(outDir, outfile), res=200, width=1600, height=1000)
myPlot = ggplot(deCounts, aes_string(x="cellType", y="n", fill="Covariate")) +
			geom_bar(position="dodge", stat="identity") +
			monocle_theme_opts() + 
			 # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			 theme(axis.text.x = element_text(angle = 45,  hjust=1)) + 
			 theme(text = element_text(size = 20))  + ylab("DE Genes") +
			 xlab("Cell Type") + 
			 scale_fill_brewer(palette = "Accent")

print(myPlot)
dev.off()


# Output a CSV holding these hits
for (eachCovariate in covariatesToPlot){
	fullDEfile = paste0(outDir, "Combined_DE_by_", eachCovariate, "_q_", as.character(opt$padjCutoff), ".csv")
	write.csv(combinedCSVlist[[eachCovariate]], file=fullDEfile)

} 







# Look at all results, not q value filtered:

combinedCSVlistFull = getCombinedCSV(opt, cellTypes, covariateCSVs, filterbyQval = FALSE)

sexResFull = combinedCSVlistFull$SexM

outputFullCSV(combinedCSVlistFull)

# Plot P value distribution

plotPvalDists = function(inputList){
	outDir = "./plots/DE_Summaries/p_val_dist/"
	dir.create(outDir)

	covarToOutput = c("Age", "SexM")
	cellTypes = c("Endothelium", "Fibroblast", "Lymphocyte", "Myeloid", "Perivascular", "Ventricular_Cardiomyocytes", "Neuron", "Adipocytes")
	for (eachCovar in covarToOutput){
		thisCovarCSV = inputList[[eachCovar]]

		for (eachType in cellTypes){
			print(paste0("Working on ", eachType))
			miniCSV = thisCovarCSV[thisCovarCSV$cellType == eachType,]
			# Plot a histogram of the p values here
			png(paste0(outDir, eachType, "_by_", eachCovar, "_pvals.png"), res=200, height=800, width=800)
			myPlot = ggplot(miniCSV, aes_string(x="pval")) + 
					geom_histogram() +
					ggtitle(paste0(eachType, " by ", eachCovar))
			print(myPlot)
			dev.off()
		}
	}


}




plotPvalDists(combinedCSVlistFull)










