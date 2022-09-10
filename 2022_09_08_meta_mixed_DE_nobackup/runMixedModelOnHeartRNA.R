


###########################################################################################
# Load relevant packages
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
library(monocle3)
# print(modStatus)
library("optparse")

library("nebula")

# library("lme4")
# Source the draft code for mixed model fitting
# source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/mixedModelFittingMonocle3Style.R")

filterByCellsExpressed = function(inputCDS, cellPropMin){
	rowData(inputCDS)$fracCellsExpr = (rowSums(exprs(inputCDS) != 0) /
								nrow(colData(inputCDS)) )
	inputCDS = inputCDS[rowData(inputCDS)$fracCellsExpr > cellPropMin,]
	return(inputCDS)
}

getMixedModelFormula = function(inputFixedEffects, inputRandomEffects){
	# Formatted Random Effects
	fixedFormulaPart = paste(inputFixedEffects, collapse = " + ")

	# For now, assume there are random effects (otherwise we'd just use the fit_models code in monocle3...)
	randomFormulaPart = paste(inputRandomEffects, collapse = ") + (1|")
	# Add the front and back of the random intercepts formula
	randomFormulaPart = paste0("(1|", randomFormulaPart, ")")

	fullFormula = paste0("~", fixedFormulaPart, " + ", randomFormulaPart)
	return(fullFormula)
}

fit_nebula <- function(inputCDS, modelFormula, fixedEffectVec){
  # Requires cells to be ordered by individual
  inputCDS = inputCDS[,order(colData(inputCDS)$Donor)]
  # Get design matrix to input to nebula
  nebulaDF = model.matrix(as.formula(paste0("~", paste(fixedEffectVec, collapse = "+"))), data = as.data.frame(colData(inputCDS)))
  offsetVal = colData(inputCDS)$size_factor 

  print(str(nebulaDF))
  print(str(offsetVal))

  nebRes = nebula(exprs(inputCDS), colData(inputCDS)$Donor, pred=nebulaDF, offset = (1.0/colData(inputCDS)$Size_Factor))
  return(nebRes)
}

# Get the passed parameters
option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
  			default="NucleiOnlySharedGenesCDS", 
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-c", "--cellType"), type="character", 
               default="Adipocytes", #"Ventricular_Cardiomyocytes", #  default="Adipocytes", #default="Endocardium", 
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
              help="Exclude samples from atrium", metavar="logical")
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

# Relevel so LV (most common) is baseline
colData(testCDS)$Anatomical_Site = relevel(as.factor(colData(testCDS)$Anatomical_Site), ref="LV")

# # # Make this a mini run to test
# # ############################################################## Testing only 7-19-21
processingNote = paste0(processingNote, "miniTest")
testCDS = testCDS[1:20,]
# # ###################################################################################

colData(testCDS)$Donor = paste0(colData(testCDS)$Donor, "_", colData(testCDS)$DataSource)

# source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/utils-MM.R")
# source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/expr_models-MM.R")
# source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/mixedModelFittingMonocle3Style.R")

startTime = Sys.time()
# Now test for these genes
print("Fitting Models Now")
# mixedTestRes = fit_mixed_models(testCDS, modelFormula, #c(fixedEffectVec, randomEffectVec),
#                             expression_family="mixed-negbinomial", clean_model=T)
mixedTestRes = fit_nebula(testCDS, modelFormula, fixedEffectVec)
endTime = Sys.time()

print(endTime - startTime)



testResDF = as.data.frame(mixedTestRes)
testResDF = testResDF[, !(names(testResDF) %in% c("model"))]

# Save the output
dir.create("./rdsOutput/")
dir.create(paste0("./rdsOutput/nebulaFits/"))
outPath = paste0("./rdsOutput/nebulaFits/", processingNote, "_MMresult")





saveRDS(testResDF, outPath)

print("All Done")




# for (eachRowInd in 1:nrow(testCDS)){
#   print(eachRowInd)
#   testFit = fit_mixed_models(testCDS[eachRowInd,], modelFormula, #c(fixedEffectVec, randomEffectVec),
#                             expression_family="mixed-negbinomial", clean_model=T)
# }





















