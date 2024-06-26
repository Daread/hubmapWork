
library(ggplot2)
library(dplyr)
library("optparse")
library(stringr)
library(ggrepel)
# options(ggrepel.max.overlaps = Inf)

print("Libraries loaded, starting now")


# Get the passed parameters
option_list = list(
  make_option(c("-m", "--modelNotes"), type="character", 
  			# default="_fix_Anatomical_Site_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult", 
  			# default="_fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult",
  			# default="_fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult",
  			default="_fix_Anatomical_Site,Age,Sex_rand_Donor_PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyLabels_regressProtocad_cds_p_MMresult",
              help="Processing note from model fitting", metavar="character"),
  # make_option(c("-c", "--cellType"), type="character", 
  # 			# default="Endocardium", 
  # 			default="Vascular_Endothelium",
  #             help="Cell type for which the model was fit", metavar="character"),
   make_option(c("-p", "--padjCutoff"), type="numeric", 
        default=0.05,
              help="Max padj value to plot as significant", metavar="numeric"),

   make_option(c("-f", "--fetalInput"), type="character", 
        default="./fileInputs/motifs_all_tissues.for_website.csv",
              help="Max padj value to plot as significant", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



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



getCombinedCSV <- function(opt, inputCelltypes, inputcellTypeList){

    # Read the per cell type csv, combine into a larger on
  myDF = data.frame()
  for (eachType in inputCelltypes){
    # inputDAfile = paste0("./plots/cellTypeSpec/Cell_Type_Specificity_for", eachType, opt$modelNotes, "/", eachType, "_is_", eachType, "_AllTestsRun_Table.csv")
    inputDAfile = paste0("./plots/cellTypeSpec/Cell_Type_Specificity_for", eachType, "/", eachType, "_is_", eachType, "_AllTestsRun_Table.csv")
    thisCSV = read.csv(inputDAfile)
    # Keep only the pathway name, enrichment, and p val. 
    # browser()
    thisCSV = thisCSV[c("gene", "coefficientValue", "q_val")]
    # Add a column for the cell type, for tidy format plotting later
    thisCSV$cellType = eachType
    myDF = rbind(myDF, thisCSV)
  }

  # browser()
  # Add columns for absolute value of effect size and direction of effect (for plotting)
  myDF$Effect_Magnitude = abs(myDF$coefficientValue)
  myDF$Effect_Direction = ifelse(myDF$coefficientValue > 0.0, "Positive", "Negative")

  # myDF = myDF[myDF$q_val < opt$padjCutoff,]
  
  # Return the whole csv
  return(myDF)
}




getFetalCSV <- function(opt, typesToKeep){
  # Read in the file as it originally was set up
  rawCSV = read.csv(opt$fetalInput)
  # Keep the subset of types that have matches in adult data
  # Myeloid cells = closest to macrophages, at the level of high-level cell types in the fetal atlas
  # Thymocytes = closest to mature T cells
  # Watch out for "Cardiomyocyte lineage/Endothelial cell". Maybe a better fit for endocardium?
  # Use "Epicardial fat cells" for adipocytes
  # "Lymphatic endothelial cells"
  # "Myeloid/Lymphoid cell"
  # "Purkinje neurons" for neurons

  # # Backup, 3-31
  # fetalTypesToUse = c("Cardiomyocytes", "Vascular endothelial cells", "Endocardial cells", "Myeloid cells", "Smooth muscle cells", "Stromal cells")
  # equivalentAdult = c( "Cardiomyocyte", "Vascular_Endothelium",  "Endocardium", "Macrophage",  "VSM_and_Pericyte", "Fibroblast") 
  # names(equivalentAdult) = fetalTypesToUse
  # "Endocardial cells"

  # Backup, 3-31
  fetalTypesToUse = c("Cardiomyocytes", "Vascular endothelial cells", "Endocardial cells", "Myeloid cells", "Smooth muscle cells", "Stromal cells",
  						"Epicardial fat cells", "Purkinje neurons", "Thymocytes")
  equivalentAdult = c( "Cardiomyocyte", "Vascular_Endothelium",  "Endocardium",                                  "Macrophage",  "VSM_and_Pericyte",      "Fibroblast",
  						"Adipocytes", "Neuronal", "T_Cell") 
  names(equivalentAdult) = fetalTypesToUse


  # Map fetal atlas cell type names to closest equivalent in adult tissue. Add this info, then return the desired types
  returnCSV = rawCSV[rawCSV$cell_type %in% fetalTypesToUse,]
  returnCSV$adultMatchType = as.vector(equivalentAdult[returnCSV$cell_type])

  # Remove hits below the desired q value cutoff
  # returnCSV = returnCSV[returnCSV$qvalue < opt$padjCutoff,]
  return(returnCSV)
}


defineGenesToPlot <- function(inputDF, cellTypeHere){

	inputDF$labelGene = ""

	# For genes that I want to show on the scatter plot, add the motif name to the labelGene column
	if (cellTypeHere == "Vascular_Endothelium"){
		# Show clusters of genes above and below the y=x line. Also show a group at the end of highest enrichments
		# High in adult, down in fetal
		inputDF$labelGene = ifelse((inputDF$Adult_Fold_Change > 1.08 &  inputDF$Fetal_Fold_Change < .975), inputDF$Motif, inputDF$labelGene)
		# High in fetal, down in adult
		# inputDF$labelGene = ifelse((inputDF$Adult_Fold_Change <.8 &  inputDF$Fetal_Fold_Change >1.025), inputDF$Motif, inputDF$labelGene)
		# # High in both
		# inputDF$labelGene = ifelse(	(inputDF$Fetal_Fold_Change > 1.1) | 
		# 						(inputDF$Adult_Fold_Change > 1.3 ), inputDF$Motif, inputDF$labelGene)
		# High in fetal, little chang in adult
		inputDF$labelGene = ifelse(	inputDF$Fetal_Fold_Change > 1.075 & inputDF$Adult_Fold_Change < 1.05, inputDF$Motif, inputDF$labelGene)
	}
	# Cardiomyocytes
	if (cellTypeHere == "Cardiomyocyte"){
		# Show clusters of genes above and below the y=x line. Also show a group at the end of highest enrichments
		# High in adult, down in fetal
		# inputDF$labelGene = ifelse((inputDF$Adult_Fold_Change > 1.08 &  inputDF$Fetal_Fold_Change < .975), inputDF$Motif, inputDF$labelGene)
		# # High in fetal, down in adult
		# inputDF$labelGene = ifelse((inputDF$Adult_Fold_Change <.8 &  inputDF$Fetal_Fold_Change >1.025), inputDF$Motif, inputDF$labelGene)
		# High in both
		inputDF$labelGene = ifelse(	inputDF$Fetal_Fold_Change > 1.3 & inputDF$Adult_Fold_Change > 1.275, inputDF$Motif, inputDF$labelGene)
	}
	# Macrophages
	if (cellTypeHere == "Macrophage"){
		# Show clusters of genes above and below the y=x line. Also show a group at the end of highest enrichments
		# High in adult, down in fetal
		# inputDF$labelGene = ifelse((inputDF$Adult_Fold_Change > 1.08 &  inputDF$Fetal_Fold_Change < .975), inputDF$Motif, inputDF$labelGene)
		# # High in fetal, down in adult
		inputDF$labelGene = ifelse((inputDF$Adult_Fold_Change >1.1 &  inputDF$Fetal_Fold_Change <.9), inputDF$Motif, inputDF$labelGene)
		# High in both
		inputDF$labelGene = ifelse(	inputDF$Fetal_Fold_Change > 1.12 | inputDF$Adult_Fold_Change > 1.25, inputDF$Motif, inputDF$labelGene)
	}
	# T cells
	# if (cellTypeHere == "T_Cell"){
	# 	# High in both
	# 	inputDF$labelGene = ifelse(	inputDF$Fetal_Fold_Change < .9 & inputDF$Adult_Fold_Change > 1.1, inputDF$Motif, inputDF$labelGene)
	# }
	# T cells
	if (cellTypeHere == "Fibroblast"){
		# High in both
		inputDF$labelGene = ifelse(	inputDF$Adult_Fold_Change > 1.2, inputDF$Motif, inputDF$labelGene)
	}
	if (cellTypeHere == "Neuronal"){
		# High in both
		inputDF$labelGene = ifelse(	inputDF$Fetal_Fold_Change > 1.3 | inputDF$Adult_Fold_Change > 1.3, inputDF$Motif, inputDF$labelGene)
	}

	return(inputDF)
}


findMotifOverlap <- function(adultCombinedCSV, fetalCombinedCSV, opt, cellTypes){
	# # Filter down to the right q value
	# adultCombinedCSV = adultCombinedCSV[adultCombinedCSV$q_val < opt$padjCutoff,]
	# fetalCombinedCSV = fetalCombinedCSV[fetalCombinedCSV$qvalue < opt$padjCutoff,]

	overlapDF = data.frame(Motif_Count = numeric(), Overlap_Type=character(), Cell_Type=character())
	overlapTypes = c("Adult_Only", "Fetal_Only", "Both", "Random")
	# Loop. For each cell type, find the motifs that overlap and are unique to each
	for (eachType in cellTypes){
		# motifsInAdult = adultCombinedCSV[adultCombinedCSV$cellType == eachType,][["gene"]]
		# motifsInFetal = fetalCombinedCSV[fetalCombinedCSV$adultMatchType == eachType,][["motif"]]

		testsInAdult = adultCombinedCSV[adultCombinedCSV$cellType == eachType,]
		testsInFetal = fetalCombinedCSV[fetalCombinedCSV$adultMatchType == eachType,]

		# Get the number of tests run
		adultTestCount = nrow(testsInAdult)
		fetalTestCount = nrow(testsInFetal)

		# browser()
		# Only use positive results
		testsInAdult = testsInAdult[testsInAdult$coefficientValue > 0.0,]
		testsInFetal = testsInFetal[testsInFetal$fold_change > 1.0,]

		# Filter by p value
		motifsInAdult = testsInAdult[testsInAdult$q_val < opt$padjCutoff,][["gene"]]
		motifsInFetal = testsInFetal[testsInFetal$qvalue < opt$padjCutoff,][["motif"]]

		# Get the expected random overlap number
		expectedOverlap = length(motifsInFetal) * (1.0 * length(motifsInAdult) / adultTestCount)

		# Get the overlaps
		adultOnly = length(motifsInAdult) - sum(as.numeric(motifsInAdult %in% motifsInFetal))
		fetalOnly = length(motifsInFetal) - sum(as.numeric(motifsInFetal %in% motifsInAdult))
		inBoth = sum(as.numeric(motifsInAdult %in% motifsInFetal))
		# browser()

		miniDF = data.frame(Motif_Count = c(adultOnly, fetalOnly, inBoth, expectedOverlap), 
							Overlap_Type = overlapTypes, Cell_Type = eachType)
		overlapDF = rbind(overlapDF, miniDF)
	}

	return(overlapDF)
}

plotOverlaps <- function(overlapDF, opt, outDir = "./plots/"){

	# Make a bar plot
	outputFile = paste0(outDir, "Fetal_and_Adult_Motif_Overlap_qval_", as.character(opt$padjCutoff), ".png")
	png(outputFile, res=200, width=1200, height=1000)
	myPlot = ggplot(overlapDF, aes_string(x="Cell_Type", y="Motif_Count", fill="Overlap_Type")) + 
			geom_bar(position="dodge", stat="identity") + 
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

	print(myPlot)
	dev.off()	


	overlapDF = overlapDF[overlapDF$Overlap_Type %in% c("Both", "Random"),]
	overlapDF$Overlap_Type = ifelse(overlapDF$Overlap_Type == "Both", "Observed Intersection", "Random Intersection")

	outputFile = paste0(outDir, "Fetal_and_Adult_Motif_Only_Overlap_qval_", as.character(opt$padjCutoff), ".png")
	png(outputFile, res=200, width=1200, height=1000)
	myPlot = ggplot(overlapDF, aes_string(x="Cell_Type", y="Motif_Count", fill="Overlap_Type")) + 
			geom_bar(position="dodge", stat="identity") + 
			theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
            theme(text = element_text(size = 20)) + xlab("Cell Type") + 
            	ylab("TF Count") + 
            	 guides(fill=guide_legend(title="Intersection"))

	print(myPlot)
	dev.off()

}


formatCellType <- function(inputColumn){

	inputColumn = ifelse(inputColumn == "T_Cell", "T Cell", inputColumn)
	inputColumn = ifelse(inputColumn == "VSM_and_Pericyte", "Perivascular Cell", inputColumn)
	inputColumn = ifelse(inputColumn == "Vascular_Endothelium", "Vascular Endothelium", inputColumn)

	return(inputColumn)
}



plotCSVcomparisons <- function(adultCSV, fetalCSV, opt, cellTypes, outDir="./plots/"){

	correlationDF = data.frame("Type"=cellTypes, "Correlation"=100, "pValue" = 100)
	rownames(correlationDF) = cellTypes

	# See how the fold changes correlate
	for (eachType in cellTypes){
		# Get and format the csv entries for this cell type
		miniFetal = fetalCSV[fetalCSV$adultMatchType == eachType,]
		miniFetal = miniFetal[c("motif", "fold_change", "qvalue")]
		colnames(miniFetal) = c("Motif", "Fetal_Fold_Change", "Fetal_Qval")

		miniAdult = adultCSV[adultCSV$cellType == eachType,]
		miniAdult = miniAdult[c("gene", "coefficientValue", "q_val")]
		miniAdult$coefficientValue = exp(miniAdult$coefficientValue)
		colnames(miniAdult) = c("Motif", "Adult_Fold_Change", "Adult_Qval")

		# Merge
		miniMerge = merge(miniAdult, miniFetal, by="Motif")

		# Get the -log(q_val)
		miniMerge$adultNegLogQ = -log10(miniMerge$Adult_Qval)
		miniMerge$fetalNegLogQ = -log10(miniMerge$Fetal_Qval)

		# Plot these
		qValCor = cor(miniMerge$adultNegLogQ, miniMerge$fetalNegLogQ)
		thisPlotFile = paste0(outDir, eachType, "q_", as.character(opt$padjCutoff), "_motif_qval_correlation.png")
		png(thisPlotFile, res=200, width=1200, height=1000)
		myPlot = ggplot(miniMerge, aes_string(x="adultNegLogQ", y="fetalNegLogQ")) +
					geom_point() + ggtitle(paste0(eachType, " -Log10(QVal) Correlation = ", as.character(qValCor)))
		print(myPlot)
		dev.off()

		# Keep those that are sig in at least one
		# miniMerge = miniMerge[(miniMerge$Adult_Qval < opt$padjCutoff),]
		miniMerge = miniMerge[(miniMerge$Adult_Qval < opt$padjCutoff | miniMerge$Fetal_Qval < opt$padjCutoff),]

		miniMerge = defineGenesToPlot(miniMerge, eachType)

		thisCorr = cor(miniMerge$Adult_Fold_Change, miniMerge$Fetal_Fold_Change)
		thisPVal = cor.test(miniMerge$Adult_Fold_Change, miniMerge$Fetal_Fold_Change)

		# browser()
		correlationDF[eachType,"Correlation"] = thisCorr
		correlationDF[eachType, "pValue"] = thisPVal$p.value

		# Debug: 2-28-2022
		# if (eachType == "Vascular_Endothelium"){
		# 	browser()
		# }

		# Plot the correlation for these
		thisPlotFile = paste0(outDir, eachType, "q_", as.character(opt$padjCutoff), "_motif_coef_correlation.png")
		png(thisPlotFile, res=200, width=1200, height=1000)
		myPlot = ggplot(miniMerge, aes_string(x="Fetal_Fold_Change", y="Adult_Fold_Change")) +
					geom_point(color="blue") + #ggtitle(paste0(eachType, " Fold Change Correlation = ", as.character(thisCorr)))
					monocle_theme_opts() + 
					theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
					xlab("Log Enrichment in Fetal") + 
					ylab("Log Enrichment in Adult")+
            		theme(text = element_text(size = 24)) +
            		geom_smooth(stat="smooth", method="lm", se=FALSE, color="black", linetype='dashed', alpha=0.5) + 
            		geom_text_repel(aes(label = labelGene), box.padding =1) 
		print(myPlot)
		dev.off()
	}

	# Save the corr df
	write.csv(correlationDF, file= paste0(outDir, "Correlations", as.character(opt$padjCutoff), "_motif_coef_correlation.csv"))

	return(correlationDF)
}



dir.create("./plots/fetalAtlasComparison/")
outputDir = paste0("./plots/fetalAtlasComparison/", opt$modelNotes, "/")
dir.create(outputDir)
# cellTypes = c("Vascular_Endothelium", "Cardiomyocyte", "Macrophage",  "VSM_and_Pericyte", "Fibroblast", "Endocardium") 

cellTypes = c( "Cardiomyocyte", "Vascular_Endothelium",  "Endocardium",                                  "Macrophage",  "VSM_and_Pericyte",      "Fibroblast",
  						"Adipocytes", "Neuronal", "T_Cell") 

# 1-28-22: Re-run "T_Cell", as naming conventions were off for that. Looks like I didn't properly re-run that after changing file name conventions for output
cellTypeCSVs = vector(mode='list', length=length(cellTypes))
names(cellTypeCSVs) = cellTypes

adultCombinedCSV = getCombinedCSV(opt, cellTypes, cellTypeCSVs)

fetalCombinedCSV = getFetalCSV(opt, cellTypes)


# 2024_05_28: Update to output data for figure plots
library(readr)
common_out <- "../combined_processed_data/"
dir.create(common_out)

write_csv(adultCombinedCSV,  paste0(common_out, "adult_motif_enrichments.csv"))
write_csv(fetalCombinedCSV,  paste0(common_out, "fetal_motif_enrichments.csv"))
####################


correlationDF = plotCSVcomparisons(adultCombinedCSV, fetalCombinedCSV, opt, cellTypes, outDir=outputDir)



# overlapDF = findMotifOverlap(adultCombinedCSV, fetalCombinedCSV, opt, cellTypes)

# overlapDF$Cell_Type = formatCellType(overlapDF$Cell_Type)






# plotOverlaps(overlapDF, opt, outDir = outputDir)







