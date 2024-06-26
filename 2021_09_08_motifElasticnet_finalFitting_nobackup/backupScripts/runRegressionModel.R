
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

  make_option(c("-a", "--alphaToUse"), type="numeric", 
        default=0.5, 
              help="Alpha value to use in glment", metavar="numeric"),

  make_option(c("-k", "--coaccessCutoff"), type="numeric", 
        default=.05,
              help="Cutoff for keeping coaccessible peaks", metavar="numeric"),

  make_option(c("-n", "--maxNdistalSites"), type="numeric", 
        default=5,
              help="Max number of sites to link to a gene's promoter", metavar="numeric"),

  make_option(c("-q", "--peakSize"), type="numeric", 
        default=600,
              help="-1 -> Keep as is, or enter a positive integer to make all peaks the same size", metavar="numeric"),

  make_option(c("-g", "--gtfFile"), type="character", 
        default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/fileOutputs/protCodingOnlyHumanGTF.gtf",
              help="Path and file for human gtf", metavar="character")


)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


outDirName = paste0("./plots/regression/",
            "Prot_Only_Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize,"pVal", as.character(opt$pValFIMOcutoff), "/" )
dir.create(outDirName)

cellTypes = c("Adipocytes", "B_Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic_Endothelium", "Macrophage", "Mast_Cell", "Neuronal", "T_Cell",
                "Vascular_Endothelium", "VSM_and_Pericyte")

inputFeatures = getInputDF(opt)
rnaData = getRNAdf(opt, cellTypes)
names(rnaData)[names(rnaData) == "id"] <- "GeneID"

# Added 9-1-21
rnaData = getProteinCodingGenes(rnaData, opt)

# Get data that combines rna expression levels and input motif data
combinedDF = getCombinedDF(inputFeatures, rnaData, opt, cellTypes)
combinedTrain = getSubsetOfCombined(combinedDF, opt, "Train")


trainX = combinedTrain[, !(colnames(combinedTrain) %in% c(paste0(cellTypes, opt$predictionTask)))]
trainY = combinedTrain[, (colnames(combinedTrain) %in% c(paste0(cellTypes, opt$predictionTask)))]

combinedVal = getSubsetOfCombined(combinedDF, opt, "Validation")
valX = combinedVal[, !(colnames(combinedVal) %in% c(paste0(cellTypes, opt$predictionTask)))]
valY = combinedVal[, (colnames(combinedVal) %in% c(paste0(cellTypes, opt$predictionTask)))]




fitAccuracyDF = data.frame("Cell_Type" = character(),
                            "Alpha" = double(),
                            "CV_R_squared" = double(),
                            "Train_R_Squared" = double(),
                            "Val_R_Squared" = double())

# rownames(fitAccuracyDF) = fitAccuracyDF$Cell_Type

# cellTypes = c("Fibroblast", "Cardiomyocyte")
exactPredictionsDF = data.frame("Placeholder" = rep(0, nrow(valY)))

set.seed(7)
for (thisCellType in cellTypes){
  print(paste0("Working on ", thisCellType))
  miniY = trainY[[paste0(thisCellType, opt$predictionTask)]]
  # trainFit = glmnet(as.matrix(trainX), miniY)

  for (thisAlpha in alphaToUse){
    print(paste0("Working on alpha = ", thisAlpha))
    trainCVfit = cv.glmnet(as.matrix(trainX), miniY, alpha=thisAlpha)
    minLambda = trainCVfit$lambda.min 

    thisDF = data.frame("dev.ratio" = trainCVfit$glmnet.fit$dev.ratio,
                        "lambda" = trainCVfit$glmnet.fit$lambda)
    thisDF = thisDF[thisDF$lambda == minLambda,]
    CVrSquared_bestLambda = thisDF[1, "dev.ratio"]

    valRsquared = getAccuracy(trainCVfit, valX, valY[[paste0(thisCellType, opt$predictionTask)]])
    trainRsquared = getAccuracy(trainCVfit, trainX, miniY)
    newRow = data.frame(thisCellType, thisAlpha, CVrSquared_bestLambda, trainRsquared, valRsquared)
    names(newRow) = c("Cell_Type", "Alpha", "CV_R_squared", "Train_R_Squared", "Val_R_Squared")
    fitAccuracyDF = rbind(fitAccuracyDF, newRow)
    # fitAccuracyDF[thisCellType, "CV_R_squared"] = rSquared_bestLambda

    # 8-31-21 edit: Add in plotting code for residuals
    if (abs(thisAlpha - .5) < .01 ){
      valAccDF = getValAccuracyDF(trainCVfit, valX,  valY[[paste0(thisCellType, opt$predictionTask)]])
      exactPredictionsDF[[paste0(thisCellType, "_Actual")]] = valAccDF[["Actual"]]
      exactPredictionsDF[[paste0(thisCellType, "_Predicted")]] = valAccDF[["Predicted"]]
    }
  }
}


# Plot the residuals
residualDir = paste0(outDirName, "residuals/")
dir.create(residualDir)

for (eachCelltype in cellTypes){
  # Plot prediction vs. actual
  png(paste0(residualDir, "Actual_vs_Pred", eachCelltype, "_",
      opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
      height=1000, width=1000, res=200)
  fitPlot = (ggplot(exactPredictionsDF, aes_string(x=paste0(thisCellType, "_Actual"),
                               y=paste0(thisCellType, "_Predicted"))) +
       ggtitle(paste0("Expression vs Prediction ", eachCelltype, "_", opt$predictionTask, "_", opt$featureSelection)) + 
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
       geom_point() )
  print(fitPlot)
  dev.off()

  exactPredictionsDF[["Residuals"]] = (exactPredictionsDF[[paste0(thisCellType, "_Actual")]] - exactPredictionsDF[[paste0(thisCellType, "_Predicted")]])
  # Plot prediction vs. actual
  png(paste0(residualDir, "Actual_vs_Residuals", eachCelltype, "_",
      opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
      height=1000, width=1000, res=200)
  fitPlot = (ggplot(exactPredictionsDF, aes_string(x=paste0(thisCellType, "_Actual"),
                               y="Residuals") )+
       ggtitle(paste0("Expression vs Prediction ", eachCelltype, "_", opt$predictionTask, "_", opt$featureSelection)) + 
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
       geom_point())
  print(fitPlot)
  dev.off()
}




fitAccuracyDF$Alpha = as.character(fitAccuracyDF$Alpha)

png(paste0(outDirName, "FitAccuracy", "glmnet_varyAlpha_",
      opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
      height=1000, width=1000, res=200)

fitPlot = (ggplot(fitAccuracyDF, aes_string(x="Cell_Type", y="CV_R_squared", color="Alpha")) +
     # geom_point() + 
     ggtitle(paste0("Fit Accuracies by CV in glmnet ", opt$predictionTask, "_", opt$featureSelection)) + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
     geom_jitter())
print(fitPlot)
dev.off()


png(paste0(outDirName, "FitAccuracy_OnValidation_", "glmnet_varyAlpha_",
      opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
      height=1000, width=1000, res=200)
fitPlot = (ggplot(fitAccuracyDF, aes_string(x="Cell_Type", y="Val_R_Squared", color="Alpha")) +
     # geom_point() + 
     ggtitle(paste0("Fit Accuracies by CV in glmnet ", opt$predictionTask, "_", opt$featureSelection)) + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
     geom_jitter())
print(fitPlot)
dev.off()



png(paste0(outDirName, "FitAccuracy_OnTrainSet_", "glmnet_varyAlpha_",
      opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
      height=1000, width=1000, res=200)
fitPlot = (ggplot(fitAccuracyDF, aes_string(x="Cell_Type", y="Train_R_Squared", color="Alpha")) +
     # geom_point() + 
     ggtitle(paste0("Fit Accuracies by CV in glmnet ", opt$predictionTask, "_", opt$featureSelection)) + 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
     geom_jitter())
print(fitPlot)
dev.off()





