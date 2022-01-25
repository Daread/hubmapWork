
hardAssignDonorAges <- function(inputCDS){
  colData(inputCDS)$Age = 0
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W134", 43, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W135", 60, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W136", 43, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W137", 49, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W139", 45, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W142", 55, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W144", 53, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W145", 51, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W146", 25, colData(inputCDS)$Age)

  colData(inputCDS)$Log10Age = log10(colData(inputCDS)$Age)

  return(inputCDS)
}


hardAssignDonorSexes <- function(inputCDS){
  colData(inputCDS)$Sex = "Not_Set"
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W134", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W135", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W136", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W137", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W139", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W142", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W144", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W145", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W146", "F", colData(inputCDS)$Sex)

  return(inputCDS)
}



###########################################################################################
# Load relevant packages
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()
print(modStatus)
library("optparse")
library("lme4")
# Source the draft code for mixed model fitting
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/mixedModelFittingMonocle3Style.R")

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

# Get the passed parameters
option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
  			default="HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes", 
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-c", "--cellType"), type="character", default="Adipocytes", #default="Endocardium", 
              help="Cell Type to subset for testing", metavar="character"),
    make_option(c("-r", "--randomEffects"), type="character", default="Donor",  # Syntax would be "Donor,Batch" if including multiple
              help="Comma-separated string of variables to model with random effects", metavar="character"),
    make_option(c("-f", "--fixedEffects"), type="character",
               # default="Anatomical_Site", # Syntax would be "Anatomical_Site,Age" if including multiple
              default="Anatomical_Site,Age,Sex",
              help="Comma-separated string of variables to model with fixed effects", metavar="character"),
    make_option(c("-g", "--groupColumn"), type="character", default="highLevelCellType", 
              help="column in colData that cellType uses to select", metavar="character"),
    make_option(c("-m", "--minGeneExprPercent"), type="numeric", default=.01, 
              help="Proportion of cells that must express a gene for it to be tested", metavar="numeric")
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
rdsPath = "../2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/"
# oldProcNote = "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes"
oldProcNote = opt$processingNote
allCellCDS = readRDS(paste0(rdsPath, "allCells_", oldProcNote, ".rds"))
print("Read in CDS")

#Added 7-28-21
allCellCDS = hardAssignDonorAges(allCellCDS)

# Added 8-16-21
allCellCDS = hardAssignDonorSexes(allCellCDS)

# Filter down to the cell type being tested here

# Getting problems 
testCDS = allCellCDS[,colData(allCellCDS)[[opt$groupColumn]] == opt$cellType]


# Added 7-20-21: Try to reduce memory footprint. Hitting issues with batch submissions
allCellCDS = NULL

# Filter down to a minimum expression level. Try 1% of cells within this subtype
testCDS = filterByCellsExpressed(testCDS, opt$minGeneExprPercent)

# 9-16-21: 
# To fix errors in the fitting package, need to raise the min expr level for some of the rarer cell types
abundantCellTypes = c("Fibroblast", "Macrophage", "Cardiomyocyte", "T_Cell", "VSM_and_Pericyte", "Vascular_Endothelium")
if (!(opt$cellType %in% abundantCellTypes)){
  # Update to require 5% expressing cells
  testCDS = filterByCellsExpressed(testCDS, .05)
}
veryRareCellTypes = c("Adipocytes")
if ((opt$cellType %in% veryRareCellTypes)){
  # Update to require 5% expressing cells
  testCDS = filterByCellsExpressed(testCDS, .1)
}

# Re-estimate size factors for this subset, in case there is some odd UMI distribution behavior specific to a cell type/subset
testCDS = estimate_size_factors(testCDS)

# Relevel so Apex is always the lowest anatomical site (and can readily be compared to the Left_ventricle in a clear coefficient)
colData(testCDS)$Anatomical_Site = relevel(as.factor(colData(testCDS)$Anatomical_Site), ref="Left_Vent")
# 8-27-21: Changed from "Aepx" to "Left_Vent" as baseline. Makes for easier LV vs. RV interpretation

# # # # Make this a mini run to test
# # # ############################################################## Testing only 7-19-21
# processingNote = paste0(processingNote, "miniTest")
# testCDS = testCDS[1:20,]
# # # ###################################################################################

source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/utils-MM.R")
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/expr_models-MM.R")
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/mixedModelFittingMonocle3Style.R")

# Now test for these genes
print("Fitting Models Now")
mixedTestRes = fit_mixed_models(testCDS, modelFormula, #c(fixedEffectVec, randomEffectVec),
                            expression_family="mixed-negbinomial", clean_model=T)

testResDF = as.data.frame(mixedTestRes)
testResDF = testResDF[, !(names(testResDF) %in% c("model"))]

# Save the output
dir.create(paste0("./rdsOutput/mixedModels/"))
outPath = paste0("./rdsOutput/mixedModels/", processingNote, "_MMresult")





saveRDS(testResDF, outPath)

print("All Done")




# for (eachRowInd in 1:nrow(testCDS)){
#   print(eachRowInd)
#   testFit = fit_mixed_models(testCDS[eachRowInd,], modelFormula, #c(fixedEffectVec, randomEffectVec),
#                             expression_family="mixed-negbinomial", clean_model=T)
# }





















