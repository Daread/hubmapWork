
library(tidyr)
library(ggplot2)
library(glmnet)

library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-p", "--pValFIMOcutoff"), type="numeric", 
  			default=1e-6, 
              help="Max FIMO p value to retain match", metavar="character"),

  make_option(c("-f", "--featureSelection"), type="character", 
        default="Binary_Combined_Motif_Counts", # "Binary_Combined_Motif_Counts", "Binary_PromOnly_Motif_Counts"
              help="Option of which features to use for input", metavar="character"),

  make_option(c("-t", "--predictionTask"), type="character", 
        default="_log2_ratio_Vs_AllTypeMean", 
              help="Option of which RNA data to predict", metavar="character"),

  make_option(c("-u", "--promoterUpstream"), type="numeric", 
        default=1000,
              help="Bases upstream of TSS to use for input features", metavar="numeric"),
  make_option(c("-d", "--promoterDownstream"), type="numeric", 
        default=200,
              help="Bases downstream of TSS to use for input features", metavar="numeric"),

  make_option(c("-k", "--coaccessCutoff"), type="numeric", 
        default=.05,
              help="Cutoff for keeping coaccessible peaks", metavar="numeric"),

  make_option(c("-n", "--maxNdistalSites"), type="numeric", 
        default=5,
              help="Max number of sites to link to a gene's promoter", metavar="numeric"),

  make_option(c("-q", "--peakSize"), type="numeric", 
        default=600,
              help="-1 -> Keep as is, or enter a positive integer to make all peaks the same size", metavar="numeric")


)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


outDirName = paste0("./plots/",
            "Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize,"pVal", as.character(opt$pValFIMOcutoff), "/" )
dir.create(outDirName)

getInputDF <- function(opt){
  # Get the path
  dirName = paste0("../2021_07_26_setup_ATAC_to_expr_data_nobackup/fileOutputs/", 
            "Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize, "/" )
  # dirName = paste0("../2021_07_26_setup_ATAC_to_expr_data_nobackup/fileOutputs/backupAllOldIntermediatesBeforeFixRaceConditionBug/", 
  #           "Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
  #               "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
  #                   "peakSize", opt$peakSize, "/" )

  # Get previous output
  dfName = paste0("Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize, "pVal", as.character(opt$pValFIMOcutoff))

  # Read in according to the featureSelection specified
  featureDF = read.csv(paste0(dirName, dfName, "_", opt$featureSelection, ".csv"))

  return(featureDF)
}

getRNAdf <- function(opt, cellTypes){
  rnaDir = paste0("../2021_07_26_setup_ATAC_to_expr_data_nobackup/rdsOutputs/RNA_Data/" )
  rnaDF = readRDS(paste0(rnaDir, "Pseudobulked_RNA_countsPerMillionNotDeduplicated.RDS"))

  print(paste0(cellTypes, opt$predictionTask))
  returnDF = rnaDF[,c("id", "gene_short_name", 
    paste0(cellTypes, opt$predictionTask))]

  return(returnDF)
}


cellTypes = c("Adipocytes", "B_Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic_Endothelium", "Macrophage", "Mast_Cell", "Neuronal", "T_Cell",
                "Vascular_Endothelium", "VSM_and_Pericyte")

inputFeatures = getInputDF(opt)

rnaData = getRNAdf(opt, cellTypes)

names(rnaData)[names(rnaData) == "id"] <- "GeneID"


getCombinedDF <- function(inputFeatures, rnaData, opt, cellTypes){
  # Need to merge these. Go down to genes that are in both of these sets
  inputFeatures = inputFeatures[inputFeatures$GeneID %in% rnaData$GeneID,]
  # ...and vice-versa
  rnaData = rnaData[rnaData$GeneID %in% inputFeatures$GeneID,]

  # 
  comboDF = merge(inputFeatures, rnaData, by="GeneID")

  comboDF = comboDF[, !(colnames(comboDF) %in% c("X", "chr", "orientation", "linkedSites"))]

  rownames(comboDF) = comboDF$Gene_ID
  return(comboDF)
}

getSubsetOfCombined <- function(combinedDF, opt, setToUse){
  combinedSubset = combinedDF[combinedDF$Model_Set == setToUse,]
  combinedSubset = combinedSubset[, !(colnames(combinedSubset) %in% c("Model_Set", "GeneID", "gene_short_name"))]

  return(combinedSubset)
}

getAccuracy <-function(fitCV_Model, inputX, inputY){
  # Get the prediction
  predictedVal = predict(fitCV_Model, newx=as.matrix(inputX), s = "lambda.min")
  # Get the r_squared
  rSquared = (cor(predictedVal, inputY) ^ 2)
  return(rSquared)
}


combinedDF = getCombinedDF(inputFeatures, rnaData, opt, cellTypes)
combinedTrain = getSubsetOfCombined(combinedDF, opt, "Train")


trainX = combinedTrain[, !(colnames(combinedTrain) %in% c(paste0(cellTypes, opt$predictionTask)))]
trainY = combinedTrain[, (colnames(combinedTrain) %in% c(paste0(cellTypes, opt$predictionTask)))]

combinedVal = getSubsetOfCombined(combinedDF, opt, "Validation")
valX = combinedVal[, !(colnames(combinedVal) %in% c(paste0(cellTypes, opt$predictionTask)))]
valY = combinedVal[, (colnames(combinedVal) %in% c(paste0(cellTypes, opt$predictionTask)))]


# thisCellType = "Fibroblast"
thisAlpha = .5
alphaToUse = c(0.0, .25, .5, .75, 1.0)

fitAccuracyDF = data.frame("Cell_Type" = character(),
                            "Alpha" = double(),
                            "CV_R_squared" = double(),
                            "Train_R_Squared" = double(),
                            "Val_R_Squared" = double())

# rownames(fitAccuracyDF) = fitAccuracyDF$Cell_Type

# cellTypes = c("Fibroblast", "Cardiomyocyte")

set.seed(7)
for (thisCellType in cellTypes){
  print(paste0("Working on ", thisCellType))
  miniY = trainY[[paste0(thisCellType, opt$predictionTask)]]
  trainFit = glmnet(as.matrix(trainX), miniY)

  for (thisAlpha in alphaToUse){
    print(paste0("Working on alpha = ", thisAlpha))
    trainCVfit = cv.glmnet(as.matrix(trainX), miniY)
    minLambda = trainCVfit$lambda.min 

    thisDF = data.frame("dev.ratio" = trainCVfit$glmnet.fit$dev.ratio,
                        "lambda" = trainCVfit$glmnet.fit$lambda)
    thisDF = thisDF[thisDF$lambda == minLambda,]
    CVrSquared_bestLambda = thisDF[1, "dev.ratio"]

    valRsquared = getAccuracy(trainCVfit, valX, valY[[paste0(thisCellType, opt$predictionTask]])
    trainRsquared = getAccuracy(trainCVfit, trainX, miniY)
    newRow = data.frame(thisCellType, thisAlpha, CVrSquared_bestLambda, trainRsquared, valRsquared)
    names(newRow) = c("Cell_Type", "Alpha", "CV_R_squared", "Train_R_Squared", "Val_R_Squared")
    fitAccuracyDF = rbind(fitAccuracyDF, newRow)
    # fitAccuracyDF[thisCellType, "CV_R_squared"] = rSquared_bestLambda
    
  }

}



fitAccuracyDF$Alpha = as.character(fitAccuracyDF$Alpha)


png(paste0(outDirName, "FitAccuracy", "glmnet_varyAlpha_",
      opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
      height=1000, width=1000, res=200)

fitPlot = (ggplot(fitAccuracyDF, aes_string(x="Cell_Type", y="CV_R_squared", color="Alpha")) +
     # geom_point() + 
     ggtitle(paste0("Fit Accuracies by CV in glmnet")) + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
     geom_jitter())
print(fitPlot)
dev.off()


png(paste0(outDirName, "FitAccuracy_OnValidation_", "glmnet_varyAlpha_",
      opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
      height=1000, width=1000, res=200)
fitPlot = (ggplot(fitAccuracyDF, aes_string(x="Cell_Type", y="Val_R_Squared", color="Alpha")) +
     # geom_point() + 
     ggtitle(paste0("Fit Accuracies by CV in glmnet")) + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
     geom_jitter())
print(fitPlot)
dev.off()



png(paste0(outDirName, "FitAccuracy_OnTrainSet_", "glmnet_varyAlpha_",
      opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
      height=1000, width=1000, res=200)
fitPlot = (ggplot(fitAccuracyDF, aes_string(x="Cell_Type", y="Train_R_Squared", color="Alpha")) +
     # geom_point() + 
     ggtitle(paste0("Fit Accuracies by CV in glmnet")) + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
     geom_jitter())
print(fitPlot)
dev.off()























