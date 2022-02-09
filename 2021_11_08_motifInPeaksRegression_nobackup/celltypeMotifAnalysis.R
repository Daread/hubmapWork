
library(ggplot2)
library(dplyr)
library("optparse")
library(stringr)


print("Libraries loaded, starting now")


# Get the passed parameters
option_list = list(
  make_option(c("-m", "--modelNotes"), type="character", 
  			# default="_fix_Anatomical_Site_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult", 
  			default="_fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult",
  			# default="_fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult",
  			# default="_fix_Anatomical_Site,Age,Sex,log10umi_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult",
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




getCombinedCSV <- function(opt, inputCelltypes, inputcellTypeList){

    # Read the per cell type csv, combine into a larger on
  myDF = data.frame()
  for (eachType in inputCelltypes){
    inputDAfile = paste0("./plots/cellTypeSpec/Cell_Type_Specificity_for", eachType, opt$modelNotes, "/", eachType, "_is_", eachType, "_AllTestsRun_Table.csv")
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
  fetalTypesToUse = c("Cardiomyocytes", "Vascular endothelial cells", "Endocardial cells", "Myeloid cells", "Smooth muscle cells", "Stromal cells")
  equivalentAdult = c( "Cardiomyocyte", "Vascular_Endothelium",  "Endocardium", "Macrophage",  "VSM_and_Pericyte", "Fibroblast") 
  names(equivalentAdult) = fetalTypesToUse

  # Map fetal atlas cell type names to closest equivalent in adult tissue. Add this info, then return the desired types
  returnCSV = rawCSV[rawCSV$cell_type %in% fetalTypesToUse,]
  returnCSV$adultMatchType = as.vector(equivalentAdult[returnCSV$cell_type])

  # Remove hits below the desired q value cutoff
  # returnCSV = returnCSV[returnCSV$qvalue < opt$padjCutoff,]

  return(returnCSV)
}


plotCSVcomparisons <- function(adultCSV, fetalCSV, opt, cellTypes, outDir="./plots/"){
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

		# browser()


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
		miniMerge = miniMerge[(miniMerge$Adult_Qval < opt$padjCutoff),]
		# miniMerge = miniMerge[(miniMerge$Adult_Qval < opt$padjCutoff | miniMerge$Fetal_Qval < opt$padjCutoff),]

		thisCorr = cor(miniMerge$Adult_Fold_Change, miniMerge$Fetal_Fold_Change)

		# Plot the correlation for these
		thisPlotFile = paste0(outDir, eachType, "q_", as.character(opt$padjCutoff), "_motif_coef_correlation.png")
		png(thisPlotFile, res=200, width=1200, height=1000)
		myPlot = ggplot(miniMerge, aes_string(x="Fetal_Fold_Change", y="Adult_Fold_Change")) +
					geom_point() + ggtitle(paste0(eachType, " Fold Change Correlation = ", as.character(thisCorr)))
		print(myPlot)
		dev.off()
	}
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
}



outputDir = paste0("./plots/fetalAtlasComparison/", opt$modelNotes, "/")
dir.create(outputDir)
cellTypes = c("Vascular_Endothelium", "Cardiomyocyte", "Macrophage",  "VSM_and_Pericyte", "Fibroblast", "Endocardium") 
# 1-28-22: Re-run "T_Cell", as naming conventions were off for that. Looks like I didn't properly re-run that after changing file name conventions for output
cellTypeCSVs = vector(mode='list', length=length(cellTypes))
names(cellTypeCSVs) = cellTypes

adultCombinedCSV = getCombinedCSV(opt, cellTypes, cellTypeCSVs)

fetalCombinedCSV = getFetalCSV(opt, cellTypes)



plotCSVcomparisons(adultCombinedCSV, fetalCombinedCSV, opt, cellTypes, outDir=outputDir)



overlapDF = findMotifOverlap(adultCombinedCSV, fetalCombinedCSV, opt, cellTypes)

plotOverlaps(overlapDF, opt, outDir = outputDir)







