
# load requirements
suppressPackageStartupMessages({
  # library(ArchR)
  library(monocle3)
  library(dplyr)
  library(ggrastr)
  # library(Seurat)
})
library(tidyr)
library(plyr)
library(ggplot2)
library(Matrix)
library("lme4")

assignDonorAndSite <- function(inputCDS){
	updatedColdata = separate(data=as.data.frame(colData(inputCDS)), col="sampleName",
							into=c("Donor", "Tissue", "Site"), sep="[.]", remove=FALSE)
	colData(inputCDS)$Donor = updatedColdata$Donor
	colData(inputCDS)$Anatomical_Site = updatedColdata$Site
	return(inputCDS)
}

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



source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

library("optparse")

# Get the passed parameters
option_list = list(

    make_option(c("-w", "--cdsToSave"), type="numeric", default=NULL,   # Peak_CDS
              help="Features of cds to save", metavar="character"),

    make_option(c("-k", "--kVal"), type="numeric", default=20, 
              help="K value for clustering", metavar="character"),
    make_option(c("-m", "--matrixPath"), type="character", default="/net/trapnell/vol1/HuBMAP/novaseq/210111_Riza_sciATAC3_split_sample/analyze_out/", 
              help="Path to peak x motif matrices", metavar="character"),
    make_option(c("-z", "--pcToUse"), type="numeric", default=50, 
              help="Principel components/LSI coords to use", metavar="numeric"),

     make_option(c("-n", "--motifSetNote"), type="character", default="",
          ##default="/net/trapnell/vol1/HuBMAP/novaseq/210111_Riza_sciATAC3_split_sample/analyze_out/", 
              help="Note about motif collection tested", metavar="character"),

  make_option(c("-b", "--pvalToUse"), type="character", default=NULL, 
              help="Principel components/LSI coords to use", metavar="character"),

    make_option(c("-c", "--cellType"), type="character", default="Cardiomyocyte", #default="Endocardium", 
              help="Cell Type to subset for testing", metavar="character"),

    make_option(c("-p", "--cdsPath"), type="character",# default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/archr/results/NB9/All_Cells/", # "peakMat" for greg/riza matrix, "bMat" for archr bins
                                                                              # "gMat" for archr activity scores
                                  default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/archr/results/NB11/All_Cells/",
              help="Path to where the ATAC cds is saved", metavar="character"),
    make_option(c("-g", "--groupColumn"), type="character", default="harmonyKNN_type",     # default="Assigned_Cell_Type", 

              help="column in colData that cellType uses to select", metavar="character"),

    make_option(c("-r", "--randomEffects"), type="character", default="Donor",  # Syntax would be "Donor,Batch" if including multiple
              help="Comma-separated string of variables to model with random effects", metavar="character"),
    make_option(c("-f", "--fixedEffects"), type="character",
               # default="Anatomical_Site", # Syntax would be "Anatomical_Site,Age" if including multiple
              default="Anatomical_Site,Age,Sex",
              help="Comma-separated string of variables to model with fixed effects", metavar="character"),

    make_option(c("-a", "--ATACprocNote"), type="character", #default="FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20", 
                    default="Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds",
              help="How ATAC Cells were filtered", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# Parse random and fixed effects
fixedEffectVec = strsplit(opt$fixedEffects, ",")[[1]]
randomEffectVec = strsplit(opt$randomEffects, ",")[[1]]

# Add in the indicator for the cell type being tested
fixedEffectVec = c(fixedEffectVec, paste0("is_", opt$cellType))


# Get the Model formula from these vectors
modelFormula = getMixedModelFormula(fixedEffectVec, randomEffectVec)

processingNote = paste0("Cell_Type_Specificity_for", opt$cellType, "_fix_", opt$fixedEffects, "_rand_", opt$randomEffects, "_", opt$ATACprocNote)


samplesATACnames = c(
  "W134.heart.apex.s1",
  "W135.heart.LV.s1", 
  "W136.heart.apex.s1", "W136.heart.LV.s1", 
  "W139.heart.apex.s1", "W139.heart.LV.s1", "W139.heart.RV.s1", "W139.heart.septum.s1", 
  "W142.heart.LV.s1", 
  "W144.heart.apex.s1", 
  "W145.heart.apex.s1", "W145.heart.LV.s1", 
  "W146.heart.apex.s1", "W146.heart.LV.s1")
names(samplesATACnames) = c("W134.Apex",
  "W135.Left.Vent", 
  "W136.Apex", "W136.Left.Vent", 
  "W139.Apex", "W139.Left.Vent", "W139.Right.Vent", "W139.Septum", 
  "W142.Left.Vent", 
  "W144.Apex", 
  "W145.Apex", "W145.Left.Vent", 
  "W146.Apex", "W146.Left.Vent")

getPeakMotifMat <- function(opt, sampleNames){

	# All samples have the same matrix. Get one and return it (will be a dgTmatrix)
	eachSample = sampleNames[1]
	thisMatrixPath = paste0(opt$matrixPath, eachSample, "/motif_matrices/", eachSample, "-peak_motif_matrix.mtx.gz")
	thisMatrix = readMM(gzfile(thisMatrixPath))

	matrixCol = read.table(paste0(opt$matrixPath, eachSample, "/motif_matrices/", eachSample, "-peak_motif_matrix.columns.txt"),
							header=FALSE)
	matrixRow = read.table(paste0(opt$matrixPath, eachSample, "/motif_matrices/", eachSample, "-peak_motif_matrix.rows.txt"),
							header=FALSE)

	colnames(thisMatrix) = paste0("chr", matrixCol$V1 )
	rownames(thisMatrix) = matrixRow$V1

	return(thisMatrix)
}




getPeakMotifMatCustomPval <- function(opt){

  # Get the matrix file name based on pvalue and the path
  matrixPath = paste0(opt$matrixPath, opt$motifSetNote, "peak_x_motif_matrix_pVal", as.character(opt$pvalToUse))
  thisMatrix = readMM(gzfile(paste0(matrixPath, ".mtx")))

  matrixCol = read.table(paste0(matrixPath, ".columns.txt"), header=FALSE)
  matrixRow = read.table(paste0(matrixPath, ".rows.txt"), header=FALSE)

  colnames(thisMatrix) = paste0("chr", matrixCol$V1 )
  rownames(thisMatrix) = matrixRow$V1

  return(thisMatrix)
}




# First, read in the peak x motif matrices and make sure we have a single matrix for the shared peaks
if (  is.null(opt$pvalToUse) ){
  print("Using default peak matrix")
  peakMotifMat = getPeakMotifMat(opt, samplesATACnames)
  } else {
    peakMotifMat = getPeakMotifMatCustomPval(opt)
    processingNote = paste0("Cell_Type_Specificity_for", opt$cellType, "_fix_", opt$fixedEffects, "_rand_", opt$randomEffects, "_", opt$ATACprocNote, '_p', as.character(opt$pvalToUse), opt$motifSetNote)
  }


# Now get the cds holding a binary peak matrix
peakCDS = readRDS(paste0(opt$cdsPath, opt$ATACprocNote, ".rds"))

# Now multiply the motif/peak matrix by the cell x peak matrix to get cell x motif counts
# Get the matching peaks
peakMotifMat = peakMotifMat[,colnames(peakMotifMat) %in% rownames(exprs(peakCDS))]
peakMotifMat = peakMotifMat[,rownames(exprs(peakCDS))]

motifCellMat = as(peakMotifMat, "dgCMatrix") %*% exprs(peakCDS)

# Make a new cds out of this
motifDF = data.frame("Motif" = rownames(motifCellMat))
rownames(motifDF) = rownames(motifCellMat)

motifCellCDS = new_cell_data_set(motifCellMat,
								cell_metadata = as.data.frame(colData(peakCDS)),
								gene_metadata = motifDF)

# Set up a column to track if a cell is the cell type of interest or not
colData(motifCellCDS)[[paste0("is_", opt$cellType)]] = ifelse(colData(motifCellCDS)[[opt$groupColumn]] == opt$cellType, 1, 0)

testCDS = motifCellCDS

# Filter down to a minimum expression level. Try 1% of cells within this subtype
testCDS = filterByCellsExpressed(testCDS, .01)

testCDS = assignDonorAndSite(testCDS)

#Added 7-28-21
testCDS = hardAssignDonorAges(testCDS)

# Added 8-16-21
testCDS = hardAssignDonorSexes(testCDS)

# Re-estimate size factors for this subset, in case there is some odd UMI distribution behavior specific to a cell type/subset
testCDS = estimate_size_factors(testCDS)

# Relevel so Apex is always the lowest anatomical site (and can readily be compared to the Left_ventricle in a clear coefficient)
colData(testCDS)$Anatomical_Site = relevel(as.factor(colData(testCDS)$Anatomical_Site), ref="LV")

rowData(testCDS)$gene_short_name = rowData(testCDS)$Motif

source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/utils-MM.R")
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/expr_models-MM.R")
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/mixedModelFittingMonocle3Style.R")


# # # # Make this a mini run to test
# # # ############################################################## Testing only 7-19-21
# processingNote = paste0(processingNote, "miniTest")
# testCDS = testCDS[1:20,]
# # # ###################################################################################


# Now test for these genes
print("Fitting Models Now")
mixedTestRes = fit_mixed_models(testCDS, modelFormula, #c(fixedEffectVec, randomEffectVec),
                            expression_family="mixed-negbinomial", clean_model=T)

testResDF = as.data.frame(mixedTestRes)

testResDF = testResDF[, !(names(testResDF) %in% c("model"))]



######### Debugging: 11-9-21
# source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/Hierarchical_DE_Methods.R")

# oldFormatDF = as.data.frame(colData(testCDS))
# oldFormatDF$Expr = exprs(testCDS)[1,]


# oldFormatRes = RunMixedModel_lme4(oldFormatDF,"Expr", fixedEffectVec, "Donor",
#         nAGQ=0, errorModel="NegBinom")



###############














cellTypePath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/rdsOutput/mixedModels/cellTypeSpec/"
dir.create(cellTypePath)


# Save the output
# dir.create(paste0("./rdsOutput/mixedModels/"))
outPath = paste0(cellTypePath, processingNote, "_MMresult")




saveRDS(testResDF, outPath)

print("All Done")












# # > sum(peakMotifMat["SPI1",])
# # [1] 2079

# as.vector(peakMotifMat["SPI1",])



# # # Debuggin/investigating odd results: 2-8-22

# for (eachType in levels(as.factor(colData(testCDS)[[opt$groupColumn]]))){

#   miniCDS = testCDS[,colData(testCDS)[[opt$groupColumn]] == eachType]
#   print(eachType)

  
#   totalReadsAv = sum(exprs(miniCDS)) / ncol(miniCDS)
#   print(paste0("Average motifs in peaks is ",  as.character(totalReadsAv)))
  

#   thisAv = sum(exprs(miniCDS[rowData(miniCDS)$Motif == "SPI1",])) / ncol(miniCDS)
#   print(paste0("SPI1 expression is ", thisAv))

#   exprNormd = thisAv *1000 / totalReadsAv
#   print(paste0("Normalized expr is ", exprNormd))


#   umiAv = sum(colData(miniCDS)$umi) / ncol(miniCDS)
#   print(paste0("UMI Average ", as.character(umiAv)))
#   print(" ")
# }






