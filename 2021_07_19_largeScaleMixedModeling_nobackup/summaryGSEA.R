
library(ggplot2)
library(dplyr)
library("optparse")
library(stringr)
library(msigdbr)
library(fgsea)
library(tidyr)
library(tidyverse)

print("Libraries loaded, starting now")

# Get the passed parameters
option_list = list(
  make_option(c("-m", "--modelNotes"), type="character", 
  			default="_fix_Anatomical_Site,Age,Sex_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult",
              help="Processing note from model fitting", metavar="character"),
  make_option(c("-b", "--covariate"), type="character", 
  			default="SexM",   #"Age", # "SexM"
              help="Covariate to plot GSEA summary plot", metavar="character"),
  make_option(c("-p", "--padjCutoff"), type="numeric", 
  			default=0.1,
              help="Max padj value to plot as significant", metavar="numeric"),
  make_option(c("-c", "--cellType"), type="character", 
  			default="Vascular_Endothelium",
              help="Cell type for which the model was fit", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cellTypes = c("Vascular_Endothelium", "Cardiomyocyte", "Macrophage", "T_Cell", "VSM_and_Pericyte", "Fibroblast", "Endocardium")

getCombinedCSV <- function(opt, inputCelltypes){
	# Read the csv, combine into a larger on
	myDF = data.frame()
	for (eachType in inputCelltypes){
		inputGSEAfile = paste0("./plots/", eachType, opt$modelNotes, "/GSEA_", eachType, "_", opt$covariate, ".csv")
		thisCSV = read.csv(inputGSEAfile)
		# Keep only the pathway name, enrichment, and p val. 
		# browser()
		thisCSV = thisCSV[c("pathway", "NES", "padj")]
		# Add a column for the cell type, for tidy format plotting later
		thisCSV$cellType = eachType
		myDF = rbind(myDF, thisCSV)
	}

	# Add columns for absolute value of effect size and direction of effect (for plotting)
	myDF$Effect_Magnitude = abs(myDF$NES)
	myDF$Effect_Direction = ifelse(myDF$NES > 0.0, "Positive", "Negative")

	# Rename any pathways to remove the "HALLMARK_" prefix
	myDF$pathway = sub("HALLMARK_", "", myDF$pathway)

	# Return the whole csv
	return(myDF)
}


combinedCSV = getCombinedCSV(opt, cellTypes)

# Simplest version: only plot those that 
outDir = "./plots/GSEA_Summaries/"
dir.create(outDir)

# Output:
outfile = paste0("OnlySigHits_p", as.character(opt$padjCutoff), "_by_", opt$covariate, ".png" )
png(paste0(outDir, outfile), res=200, width=1200, height=1000)
myPlot = ggplot(combinedCSV[combinedCSV$padj < opt$padjCutoff,], aes_string(x="pathway", y="cellType", col="Effect_Direction")) +
			geom_point(aes(size=Effect_Magnitude)) +
			 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(myPlot)
dev.off()




