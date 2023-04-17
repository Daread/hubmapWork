
# Run to set up modules:

#########
# module unload R/4.0.0
# module unload pcre2/10.35
# module unload proj/4.9.3
# module unload gdal/2.4.1
# module unload hdf5/1.10.1

# module add sqlite/3.32.3 proj/7.1.0 gdal/3.5.2 hdf5/1.10.5 pcre2/10.39 R/4.1.2
##########


# library(nebula)
# #> Warning: package 'nebula' was built under R version 4.1.2
# data(sample_data)




setwd("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_09_08_meta_mixed_DE_nobackup")
###########################################################################################
# Load relevant packages
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
library(monocle3)
# print(modStatus)
library("optparse")

library("nebula") # Requres R4.1.2 

# library("lme4")
# Source the draft code for mixed model fitting
# source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/mixedModelFittingMonocle3Style.R")

filterByCellsExpressed <- function(inputCDS, cellPropMin){
	rowData(inputCDS)$fracCellsExpr = (rowSums(exprs(inputCDS) != 0) /
								nrow(colData(inputCDS)) )
	inputCDS = inputCDS[rowData(inputCDS)$fracCellsExpr > cellPropMin,]
	return(inputCDS)
}

getMixedModelFormula <- function(inputFixedEffects, inputRandomEffects){
	# Formatted Random Effects
	fixedFormulaPart = paste(inputFixedEffects, collapse = " + ")

	# For now, assume there are random effects (otherwise we'd just use the fit_models code in monocle3...)
	randomFormulaPart = paste(inputRandomEffects, collapse = ") + (1|")
	# Add the front and back of the random intercepts formula
	randomFormulaPart = paste0("(1|", randomFormulaPart, ")")

	fullFormula = paste0("~", fixedFormulaPart, " + ", randomFormulaPart)
	return(fullFormula)
}

getGeneSubsetCDS <- function(inputCDS, opt){
  nGenes = nrow(inputCDS)
  lowerGeneCount = as.integer((nGenes * ((opt$split -  1) / opt$nSplits)) + 1)
  upperGeneCount = as.integer(nGenes * ((opt$split) / opt$nSplits))
  return(inputCDS[lowerGeneCount:upperGeneCount,])
}

fit_nebula <- function(inputCDS, modelFormula, fixedEffectVec){
  # Requires cells to be ordered by individual
  inputCDS = inputCDS[,order(colData(inputCDS)$Donor)]
  # Get design matrix to input to nebula
  nebulaDF = model.matrix(as.formula(paste0("~", paste(fixedEffectVec, collapse = "+"))), data = as.data.frame(colData(inputCDS)))
  offsetVal = colData(inputCDS)$size_factor 

  print(str(nebulaDF))
  print(str(offsetVal))

  resList = vector(mode="list", length=5)

  browser()

  print("Default fitting:")
  startTime = Sys.time()
  nebResDefault = nebula(exprs(inputCDS), colData(inputCDS)$Donor, pred=nebulaDF, offset = colData(inputCDS)$umi, covariance = TRUE)
  resList[1] = nebResDefault
	endTime = Sys.time()
	print("Run time:")
	print(endTime - startTime)

  print("Fitting with no offset")
  startTime = Sys.time()
  nebResNoOffset = nebula(exprs(inputCDS), colData(inputCDS)$Donor, pred=nebulaDF, covariance = TRUE)
  resList[2] = nebResNoOffset
	endTime = Sys.time()
	print("Run time:")
	print(endTime - startTime)

	print("Fitting with log offset")
  startTime = Sys.time()
  nebResLogOffset = nebula(exprs(inputCDS), colData(inputCDS)$Donor, pred=nebulaDF, offset = log(colData(inputCDS)$umi), covariance = TRUE)
  resList[3] = nebResLogOffset
	endTime = Sys.time()
	print("Run time:")
	print(endTime - startTime)

	print("Fitting with size factor offset")
  startTime = Sys.time()
  nebResSFoffset = nebula(exprs(inputCDS), colData(inputCDS)$Donor, pred=nebulaDF, offset = offsetVal, covariance = TRUE)
  resList[4] = nebResSFoffset
	endTime = Sys.time()
	print("Run time:")
	print(endTime - startTime)

	print("Fitting with UMI offset, NBLMM")
	 startTime = Sys.time()
  nebResNBLMM = nebula(exprs(inputCDS), colData(inputCDS)$Donor, pred=nebulaDF,
  							 offset = colData(inputCDS)$umi, model="NBLMM", covariance = TRUE)
  resList[5] = nebResNBLMM
	endTime = Sys.time()
	print("Run time:")
	print(endTime - startTime)

  return(resList)
}

addShortNames <- function(mixedTestRes, inputCDS){
  shortNames = rowData(inputCDS)$gene_short_name
  names(shortNames) = rownames(rowData(inputCDS))
  mixedTestRes[[1]]["gene_short_name"] = shortNames[mixedTestRes[[1]][["gene"]]]
  return(mixedTestRes)
}

set.seed(7)
# Get the passed parameters
option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
  			default="NucleiOnlySharedGenesCDS", 
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-c", "--cellType"), type="character", 
               default="Endothelium",  #"Adipocytes", #"Ventricular_Cardiomyocytes", #  default="Adipocytes", #default="Endocardium", 
              help="Cell Type to subset for testing", metavar="character"),
    make_option(c("-r", "--randomEffects"), type="character", default="Donor",  # Syntax would be "Donor,Batch" if including multiple
              help="Comma-separated string of variables to model with random effects", metavar="character"),
    make_option(c("-f", "--fixedEffects"), type="character",
               # default="Anatomical_Site", # Syntax would be "Anatomical_Site,Age" if including multiple
              default="Anatomical_Site,Age,Sex,DataSource,log10_umi",
              help="Comma-separated string of variables to model with fixed effects", metavar="character"),
    make_option(c("-g", "--groupColumn"), type="character", default="Cell_Shared_Label", 
              help="column in colData that cellType uses to select", metavar="character"),
    make_option(c("-m", "--minGeneExprPercent"), type="numeric", default=.01, 
              help="Proportion of cells that must express a gene for it to be tested", metavar="numeric"),
    make_option(c("-a", "--noAtrium"), type="logical", default = TRUE,
              help="Exclude samples from atrium", metavar="logical"),
    make_option(c("-n", "--nSplits"), type="numeric", default = 10,
              help="Number of sub-jobs to split the DE testing task into", metavar="numeric"),
    make_option(c("-k", "--split"), type="numeric", default=1,
              help="Number of total splits to be run here", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# Parse random and fixed effects
fixedEffectVec = strsplit(opt$fixedEffects, ",")[[1]]
randomEffectVec = strsplit(opt$randomEffects, ",")[[1]]

print("Fixed effects are:")
print(fixedEffectVec)
print("Cell type:")
print(opt$cellType)

# Get the Model formula from these vectors
modelFormula = getMixedModelFormula(fixedEffectVec, randomEffectVec)

processingNote = paste0(opt$cellType, "_fix_", opt$fixedEffects, "_rand_", opt$randomEffects, "_", opt$processingNote)

# Get a CDS where I've already annotated the high level cell types
rdsPath = "../2022_08_22_addNewSamples_nobackup/formattedData/"
oldProcNote = opt$processingNote
allCellCDS = readRDS(paste0(rdsPath, oldProcNote, ".rds"))
print("Read in CDS")

colData(allCellCDS)$Age = as.numeric(colData(allCellCDS)$Age)

# Filter down to the cell type being tested here
testCDS = allCellCDS[,colData(allCellCDS)[[opt$groupColumn]] == opt$cellType]

if (opt$noAtrium){
  testCDS = testCDS[,!(colData(testCDS)$Anatomical_Site %in% c("LA", "RA"))]
  processingNote = paste0(processingNote, "_noAtr")
}

# Added 7-20-21: Try to reduce memory footprint. Hitting issues with batch submissions
allCellCDS = NULL

# Filter down to a minimum expression level. Try 1% of cells within this subtype
testCDS = filterByCellsExpressed(testCDS, opt$minGeneExprPercent)

# 9-16-21: 
# To fix errors in the fitting package, need to raise the min expr level for some of the rarer cell types
abundantCellTypes = c("Endothelium", "Fibroblast", "Lymphocyte", "Myeloid", "Perivascular", "Ventricular_Cardiomyocytes")
if (!(opt$cellType %in% abundantCellTypes)){
  # Update to require 5% expressing cells
  testCDS = filterByCellsExpressed(testCDS, .05)
}

# veryRareCellTypes = c("Adipocytes", "Neuron" )
# if ((opt$cellType %in% veryRareCellTypes)){
#   # Update to require 5% expressing cells
#   testCDS = filterByCellsExpressed(testCDS, .1)
# }

# Re-estimate size factors for this subset, in case there is some odd UMI distribution behavior specific to a cell type/subset
testCDS = estimate_size_factors(testCDS)
colData(testCDS)$Size_Factor = size_factors(testCDS)
# Re-set UMIs, as some genes may have been dropped when merging datasets
colData(testCDS)$umi = colSums(exprs(testCDS))
colData(testCDS)$log10_umi = log10(colData(testCDS)$umi)
names(colData(testCDS)$umi) = NULL
names(colData(testCDS)$log10_umi) = NULL
names(colData(testCDS)$Size_Factor) = NULL

# Relevel so LV (most common) is baseline
colData(testCDS)$Anatomical_Site = relevel(as.factor(colData(testCDS)$Anatomical_Site), ref="LV")

# # # # Make this a mini run to test
# # ############################################################## Testing only 7-19-21
processingNote = paste0(processingNote, "miniTest")
testCDS = testCDS[1:200,]
# # ###################################################################################

colData(testCDS)$Donor = paste0(colData(testCDS)$DataSource, "_", colData(testCDS)$Donor)

# Format so that only testing the desired subset of data
testCDS = getGeneSubsetCDS(testCDS, opt)

startTime = Sys.time()
# Now test for these genes
print("Fitting Models Now")
# mixedTestRes = fit_mixed_models(testCDS, modelFormula, #c(fixedEffectVec, randomEffectVec),
#                             expression_family="mixed-negbinomial", clean_model=T)
mixedTestRes = fit_nebula(testCDS, modelFormula, fixedEffectVec)
endTime = Sys.time()
print("Run time:")
print(endTime - startTime)

names(mixedTestRes) = c("Default", "No_Offset", "Log_Offset", "SF_Offset", "NBLMM")

library(ggplot2)
library(ggpubr)

# Make somoe plots testing consistency of results across methods
qcDir = paste0("./plots/NebulaTesting/")
dir.create(qcDir)
for (eachCol in c("logFC_(Intercept)", "logFC_log10_umi", "logFC_Age", "se_Age", "se_SexM", "logFC_SexM")){
	panelList = vector(mode="list", length=4)
	for (eachIndex in 2:5){
		# Set up mini DF with basic and comparison outputs
		miniDF = data.frame("Default"=mixedTestRes[[1]][[eachCol]], 
					"Compare" = mixedTestRes[[eachIndex]][[eachCol]])
		thisScatter = ggplot(miniDF, aes_string(x="Default", y="Compare")) + 
						geom_point() + ggtitle(paste0(eachCol, " vs ", names(mixedTestRes)[eachIndex]))
		panelList[[eachIndex - 1]] = thisScatter
	}
	print(paste0("Plotting ", eachCol))

	# Plot one comparison

	# Plot the multi-panel output
	png(paste0(qcDir, eachCol, "_compareToDefault.png"), res=200, height=1200, width=1200)
	myPlot = ggarrange(plotlist = panelList, ncol=2, nrow=2)
	print(myPlot)
	dev.off()
}







mixedTestRes = addShortNames(mixedTestRes, testCDS)


# # testResDF = as.data.frame(mixedTestRes)
# # testResDF = testResDF[, !(names(testResDF) %in% c("model"))]

# # Save the output
# dir.create("./rdsOutput/")
# dir.create(paste0("./rdsOutput/nebulaFits/"))

# # Update processing note with sub-job info
# processingNote = paste0(processingNote, "_", as.character(opt$split), "of", as.character(opt$nSplits))
# outPath = paste0("./rdsOutput/nebulaFits/", processingNote, "_MMresult.rds")



# saveRDS(mixedTestRes, outPath)

# print("All Done")
