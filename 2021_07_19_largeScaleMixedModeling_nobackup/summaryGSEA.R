
library(ggplot2)
library(dplyr)
library("optparse")
library(stringr)
library(msigdbr)
library(fgsea)
library(tidyr)
library(tidyverse)
library(RColorBrewer)

print("Libraries loaded, starting now")

# Get the passed parameters
option_list = list(
  make_option(c("-m", "--modelNotes"), type="character", 
  			default="_fix_Anatomical_Site,Age,Sex_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult",
              help="Processing note from model fitting", metavar="character"),
  make_option(c("-b", "--covariate"), type="character", 
  			default=     "Age",   #"Age", # "SexM"
              help="Covariate to plot GSEA summary plot", metavar="character"),
  make_option(c("-p", "--padjCutoff"), type="numeric", 
  			default=0.1,
              help="Max padj value to plot as significant", metavar="numeric")#,
  # make_option(c("-c", "--cellType"), type="character", 
  # 			default="Vascular_Endothelium",
  #             help="Cell type for which the model was fit", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cellTypes = c("Vascular_Endothelium", "Cardiomyocyte", "Macrophage", "T_Cell", "VSM_and_Pericyte", "Fibroblast", "Endocardium")


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
	myDF$Enrichment = abs(myDF$NES)
	myDF$Effect_Direction = ifelse(myDF$NES > 0.0, "Positive", "Negative")

	# Rename any pathways to remove the "HALLMARK_" prefix
	myDF$pathway = sub("HALLMARK_", "", myDF$pathway)

	# Return the whole csv
	return(myDF)
}


formatCellType <- function(inputColumn){

	inputColumn = ifelse(inputColumn == "T_Cell", "T Cell", inputColumn)
	inputColumn = ifelse(inputColumn == "VSM_and_Pericyte", "Perivascular Cell", inputColumn)
	inputColumn = ifelse(inputColumn == "Vascular_Endothelium", "Vascular Endothelium", inputColumn)

	return(inputColumn)
}

getFormattedEffectDirection <- function(combinedCSV, opt){

	if (opt$covariate == "SexM"){
		combinedCSV$Effect_Direction = ifelse(combinedCSV$Effect_Direction == "Positive", "Higher in Males", 
																				"Higher in Females")
	}

	if (opt$covariate == "Age"){
		combinedCSV$Effect_Direction = ifelse(combinedCSV$Effect_Direction == "Positive", "Higher in Age", 
																				"Lower in Age")
	}

	return(combinedCSV)
}


getQvalAsterisk <- function(combinedCSV, opt){
	# Add a column with one asterisk if < .1, two if < .05, three if < .01
	combinedCSV$Significance = ifelse(combinedCSV$padj >= .1, "",
								ifelse(combinedCSV$padj >= .05, "*", 
									ifelse(combinedCSV$padj >= .01, "**", "***")))
	return(combinedCSV)
}



combinedCSV = getCombinedCSV(opt, cellTypes)

combinedCSV$cellType = formatCellType(combinedCSV$cellType)

combinedCSV$pathway = gsub("_", " ", combinedCSV$pathway)

# Format the effect direction labels to indicate sex meant by a pos/negative result

combinedCSV = getFormattedEffectDirection(combinedCSV, opt)

# Simplest version: only plot those that 
outDir = "./plots/GSEA_Summaries/"
dir.create(outDir)

if (opt$covariate == "SexM"){
	figWidth = 1800
	colorPaletteToUse = "Dark2"
} else{
	figWidth = 2600
	colorPaletteToUse = "Set1"
}

# Add col with labels based on p value
combinedCSV = getQvalAsterisk(combinedCSV, opt)

# Output:
outfile = paste0("OnlySigHits_p", as.character(opt$padjCutoff), "_by_", opt$covariate, ".png" )
png(paste0(outDir, outfile), res=200, width=figWidth, height=1400)
#myPlot = ggplot(combinedCSV[combinedCSV$padj < opt$padjCutoff,], aes_string(x="pathway", y="cellType", col="Effect_Direction")) +
# myPlot = ggplot(combinedCSV[combinedCSV$padj < opt$padjCutoff,],
# 					 aes_string(x="pathway", y="cellType", col="Effect_Direction", label="Significance")) +
# 			geom_point(aes(size=Enrichment)) +
# 			geom_text(vjust=0, nudge_y=.1, colour="Black") +
# 			monocle_theme_opts() + 
# 			 theme(axis.text.x = element_text(angle = 45, hjust=1)) +              #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
# 			 theme(text = element_text(size = 18))  + ylab("Cell Type") +
# 			 xlab("Pathway") + 
# 			 guides(col=guide_legend(title="Effect Direction"))+
# 			 scale_color_brewer(palette=colorPaletteToUse)

myPlot = ggplot(combinedCSV[combinedCSV$padj < opt$padjCutoff,],
					 aes_string(x="pathway", y="cellType", col="Effect_Direction", label="Significance")) +
			geom_point(aes(size=Enrichment)) +
			monocle_theme_opts() + 
			geom_text(vjust=0, nudge_y=.1, colour="Black", size=6) +
			 theme(axis.text.x = element_text(angle = 45, hjust=1)) +              #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
			 theme(text = element_text(size = 18))  + 
			 ylab("Cell Type") +
			 xlab("Pathway") + 
			 guides(col=guide_legend(title="Effect Direction"))+
			 scale_color_brewer(palette=colorPaletteToUse)

print(myPlot)
dev.off()




