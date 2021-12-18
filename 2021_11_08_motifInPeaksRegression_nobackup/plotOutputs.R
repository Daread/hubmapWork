
library(ggplot2)
library(dplyr)
library("optparse")
library(stringr)


print("Libraries loaded, starting now")


# Get the passed parameters
option_list = list(
  make_option(c("-m", "--modelNotes"), type="character", 
  			# default="_fix_Anatomical_Site_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult", 
  			# default="_fix_Anatomical_Site,Age,Sex_rand_Donor_FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20_MMresult",
  			default="_fix_Anatomical_Site,Age,Sex_rand_Donor_Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_MMresult",
              help="Processing note from model fitting", metavar="character"),
  make_option(c("-c", "--cellType"), type="character", 
  			# default="Endocardium", 
  			default="Vascular_Endothelium",
              help="Cell type for which the model was fit", metavar="character"),

  make_option(c("-p", "--rdsPathMod"), type="character", 
  			# default="Endocardium", 
  			default="",  # "./rdsOutput/mixedModels/cellTypeSpec",
              help="Cell type for which the model was fit", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cellType = opt$cellType
modelFitDescription = opt$modelNotes

# Read in the fit
rdsPath = paste0("./rdsOutput/mixedModels/", opt$rdsPathMod)
fitResults = readRDS(paste0(rdsPath, cellType, modelFitDescription))

# Get the coefficients and p-values for all the fixed effects
fixedCoefsToGet = rownames(fitResults[1,"model_summary"][[1]]$coefficients)
fixedCoefsToGet = fixedCoefsToGet[!(fixedCoefsToGet %in% c("(Intercept)"))]

# Store the fits for different fixed effects as separate DFs, for now
fixedEffectDFlist = vector(mode="list", length=length(fixedCoefsToGet))
counter = 1
for (eachCoef in fixedCoefsToGet){
	# Set up DF
	#thisDF = data.frame()
	thisDF = NULL
	thisDF$gene = fitResults[["gene_short_name"]]
	thisDF = as.data.frame(thisDF)
	thisDF$coefficientName = eachCoef

	# Loop through and get the relevant coefficients and p values
	for (eachInd in 1:nrow(fitResults)){

		thisDF[eachInd,"coefficientValue"] = (
			fitResults[["model_summary"]][[eachInd]]$coefficients[eachCoef,"Estimate"])
		thisDF[eachInd,"pval"] = (
			fitResults[["model_summary"]][[eachInd]]$coefficients[eachCoef,"Pr(>|z|)"])
		thisDF[eachInd,"z_value"] = (
			fitResults[["model_summary"]][[eachInd]]$coefficients[eachCoef,"z value"])
		thisDF[eachInd,"coef_std_error"] = (
			fitResults[["model_summary"]][[eachInd]]$coefficients[eachCoef,"Std. Error"])
		thisDF[eachInd,"negLog10Pval"] = -log10(thisDF[eachInd,"pval"])
		thisDF[eachInd, "Donor_Stddev"] = attr(fitResults[["model_summary"]][[eachInd]]$varcor[["Donor"]], "stddev")

		# 8-9-21 addition
		thisDF[eachInd, "overDispersion"] = as.numeric(
			str_match_all(fitResults[["model_summary"]][[eachInd]]$family, "(?<=\\().+?(?=\\))")[[1]][1,1])
		thisDF[eachInd, "fracCellsExpr"] = fitResults[["fracCellsExpr"]][[eachInd]]
			
	}

	# Also get the adjust p values by benjamini hochberg
	thisDF$q_val = p.adjust(thisDF$pval, method = "BH")

	# Save
	fixedEffectDFlist[[counter]] = thisDF
	counter	= counter + 1
}
names(fixedEffectDFlist) = fixedCoefsToGet

dir.create(paste0("./plots/", opt$rdsPathMod))
outputDir = paste0("./plots/", opt$rdsPathMod, cellType, modelFitDescription, "/")
dir.create(outputDir)
print(outputDir)

qValCutoff = .1
# Make some plots
for (eachCoef in fixedCoefsToGet){
	dfHere = fixedEffectDFlist[[eachCoef]]

	# Not different levels of q value

	dfHere$q_value_range = ifelse(dfHere$q_val > .2, ">.2",
					ifelse(dfHere$q_val > .1, ".1_to_.2", 
					ifelse(dfHere$q_val > .05, ".05_to_.1", "<.05")))


	print(str(dfHere))
	# Volcano plot
	png(paste0(outputDir, cellType, "_volcano_", eachCoef, ".png"),
		height=1000, width=1000, res=200)
	myPlot = ggplot(dfHere, 
		aes_string(x="coefficientValue", y="negLog10Pval", color="q_value_range")) +
		xlab("Coefficient") + ylab("-log10(P value)") + 
		ggtitle(paste0("Volcano plot for ", eachCoef)) + geom_point()
	print(myPlot)
	dev.off()

	# Save a csv with q < .1 to take a look at later
	candidateHits = dfHere[dfHere$q_val < qValCutoff,]
	# Sort by decending qval
	candidateHits = candidateHits[order(abs(candidateHits$coefficientValue), decreasing=TRUE),]
	write.csv(candidateHits, paste0(outputDir, cellType, "_", eachCoef, "_Table.csv"))

	# Also write a csv with all tests run, not just those passing a qval threshold
	allHits = dfHere[order(abs(dfHere$coefficientValue), decreasing=TRUE),]
	write.csv(allHits, paste0(outputDir, cellType, "_", eachCoef, "_AllTestsRun_Table.csv"))

	# Don't plot points with huge std errors
	stderrorDF = dfHere[dfHere$coef_std_error < 50,]

	# 8-19-21 Add:
	# Also try plotting coefficeints vs. std error. Scuffed version of the value-suppressing color pallete idea
	png(paste0(outputDir, cellType, "_coefVsStderror_", eachCoef, ".png"),
		height=1000, width=1000, res=200)
	myPlot = ggplot(stderrorDF, 
		aes_string(x="coefficientValue", y="coef_std_error", color="q_value_range")) +
		xlab("Coefficient") + ylab("Std Error of Coefficient") + 
		ggtitle(paste0("Coef vs. Std Error for ", eachCoef)) + geom_point()
	print(myPlot)
	dev.off()
}

# Histogram of overdispersion estimates
# Clip to be in the 0-2 range
overdisperseMax = 2.0
dfHere$overDispersion = ifelse(dfHere$overDispersion > overdisperseMax, overdisperseMax, dfHere$overDispersion )

png(paste0(outputDir, cellType, "_NB_Overdispersion.png"))
myPlot = ggplot(dfHere, 
	aes_string("overDispersion")) + geom_histogram() + 
	xlab("Overdispersion from NB Fitting") +
	ggtitle(paste0("Gene overdispersion in ", cellType))
print(myPlot)
dev.off()


print("All Done")































