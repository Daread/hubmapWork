
library(tidyr)
library(ggplot2)
library(glmnet)
library(ROCR)

# Get shared functions between this and runRegressionModel.R
source("./utilityFunctions.R")

library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-p", "--pValFIMOcutoff"), type="numeric", 
  			default=1e-5, 
              help="Max FIMO p value to retain match", metavar="character"),

  make_option(c("-f", "--featureSelection"), type="character", 
        default="Binary_Combined_Motif_Counts", # "Binary_Combined_Motif_Counts", "Binary_PromOnly_Motif_Counts"
              help="Option of which features to use for input", metavar="character"),

  make_option(c("-t", "--predictionTask"), type="character", 
        default="_log2_CPM",   # "_log2_ratio_Vs_AllTypeMean", #  "_log2_CPM"
              help="Option of which RNA data to predict", metavar="character"),

  make_option(c("-z", "--highCutoff"), type="numeric", 
        default=0, 
              help="Minimum value, above which label = 1", metavar="numeric"),

  make_option(c("-l", "--lowCutoff"), type="numeric", 
        default=0, 
              help="Max value, below which label = 0", metavar="numeric"),

  make_option(c("-a", "--alphaToUse"), type="numeric", 
        default=0.5, 
              help="Alpha value to use in glment", metavar="numeric"),

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
              help="-1 -> Keep as is, or enter a positive integer to make all peaks the same size", metavar="numeric"),

  make_option(c("-g", "--gtfFile"), type="character", 
        default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/fileOutputs/protCodingOnlyHumanGTF.gtf",
              help="Path and file for human gtf", metavar="character")


)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$variableParams = paste0(opt$promoterUpstream, "_", opt$promoterDownstream, "_",
                           opt$coaccessCutoff, "_", opt$maxNdistalSites, "_", opt$peakSize )

outDirName = paste0("./plots/classification/",
            "Prot_Only_Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize,"pVal", as.character(opt$pValFIMOcutoff), "/" )
dir.create(outDirName)



cellTypes = c("Adipocytes", "B_Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic_Endothelium", "Macrophage", "Mast_Cell", "Neuronal", "T_Cell",
                "Vascular_Endothelium", "VSM_and_Pericyte")
cellTypes = c( "Macrophage", 
                "Vascular_Endothelium")



# Get the RNA data that'll be used in all models
rnaData = getRNAdf(opt, cellTypes)
names(rnaData)[names(rnaData) == "id"] <- "GeneID"
# Only use genes that are annotated as protein-coding
rnaData = getProteinCodingGenes(rnaData, opt)
rnaData = binarizeRNAdata(rnaData, opt, cellTypes)



#
getModelFitsAndAccuracy <- function(featureSelection, setsForTraining, setsForTesting, opt, rnaData, cellTypes){
  opt$featureSelection = featureSelection
  # Get the correct input features
  inputFeatures = getInputDF(opt, specifyFeatures=TRUE, specifiedFeatures = featureSelection)  # (opt, specifyFeatures=FALSE, specifiedFeatures="NULL")
  combinedDF = getCombinedDF(inputFeatures, rnaData, opt, cellTypes)

  # Now get the subset that gets used for training.
  combinedTrain = getSubsetOfCombined(combinedDF, opt, setsForTraining)
  combinedEvaluate = getSubsetOfCombined(combinedDF, opt, setsForTesting)

  trainX = combinedTrain[, !(colnames(combinedTrain) %in% c(paste0(cellTypes, opt$predictionTask)))]
  trainY = combinedTrain[, (colnames(combinedTrain) %in% c(paste0(cellTypes, opt$predictionTask)))]
  # And for evaluation data
  valX = combinedEvaluate[, !(colnames(combinedEvaluate) %in% c(paste0(cellTypes, opt$predictionTask)))]
  valY = combinedEvaluate[, (colnames(combinedEvaluate) %in% c(paste0(cellTypes, opt$predictionTask)))]

  # Make a column that makes sense for how this model is being evaluated
  evalNames = c("AUC","R^2")
  names(evalNames) = c("Classification", "Regression")
  evalCol = paste0("Eval_Set_", evalNames[opt$predictionFraming]) 
  fitAccuracyDF = data.frame("Cell_Type" = character(),
                            "Alpha" = double(),
                            evalCol = double(),
                            "Best_Lambda" = double() ) #,
                            # "Feature_Set" = featureSelection)

  set.seed(7)
  for (thisCellType in cellTypes){
    print(paste0("Working on ", thisCellType))
    miniY = trainY[[paste0(thisCellType, opt$predictionTask)]]

    # Fit the model using cross-validation to pick a best lambda
    if (opt$predictionFraming == "Classification"){
      trainCVfit = cv.glmnet(as.matrix(trainX), miniY, family="binomial", alpha=opt$alphaToUse)
    } else if (opt$predictionFraming == "Regression"){
      trainCVfit = cv.glmnet(as.matrix(trainX), miniY, alpha=opt$alphaToUse)
    }
    
    minLambda = trainCVfit$lambda.min 

    # Classification or regression?
    if (opt$predictionFraming == "Classification"){
      evalValue = getAUC(trainCVfit, valX, valY[[paste0(thisCellType, opt$predictionTask)]])
    } else if (opt$predictionFraming == "Regression"){
      evalValue = getAccuracy(trainCVfit, valX, valY[[paste0(thisCellType, opt$predictionTask)]])
    }
    
    newRow = data.frame(thisCellType,  opt$alphaToUse, evalValue, minLambda)
    names(newRow) = c("Cell_Type", "Alpha", evalCol, "Best_Lambda")
    fitAccuracyDF = rbind(fitAccuracyDF, newRow)

  }

  fitAccuracyDF$Feature_Set = featureSelection
  return(fitAccuracyDF)

}



# Now run this
opt$predictionFraming = "Classification" # Or "Regression"
evalNames = c("AUC","R^2")
names(evalNames) = c("Classification", "Regression")
evalCol = paste0("Eval_Set_", evalNames[opt$predictionFraming]) 

trainAndVal_CombinedResults = getModelFitsAndAccuracy("Binary_Combined_Motif_Counts", c("Train", "Validation"), c("Test"),
                                                         opt, rnaData, cellTypes)

trainAndVal_PromoterResults = getModelFitsAndAccuracy("Binary_PromOnly_Motif_Counts", c("Train", "Validation"), c("Test"),
                                                         opt, rnaData, cellTypes)

fullFitDF = rbind(trainAndVal_CombinedResults, trainAndVal_PromoterResults)


# Now plot
png(paste0("./plots/classification/", opt$predictionFraming, "_", opt$variableParams, "_Test_Performance_Prom_Vs_Combined.png" ),
      height=1400, width=1600, res = 200)
myPlot = ggplot(fullFitDF, aes_string(x="Cell_Type", y=evalCol, color="Feature_Set")) +
      geom_point() + ggtitle("Promoter-Only Vs. Promoter+Distal Performance") 
print(myPlot)
dev.off()




# inputFeatures = getInputDF(opt)


# # Get data that combines rna expression levels and input motif data
# combinedDF = getCombinedDF(inputFeatures, rnaData, opt, cellTypes)
# combinedTrain = getSubsetOfCombined(combinedDF, opt, "Train")


# trainX = combinedTrain[, !(colnames(combinedTrain) %in% c(paste0(cellTypes, opt$predictionTask)))]
# trainY = combinedTrain[, (colnames(combinedTrain) %in% c(paste0(cellTypes, opt$predictionTask)))]

# combinedVal = getSubsetOfCombined(combinedDF, opt, "Validation")
# valX = combinedVal[, !(colnames(combinedVal) %in% c(paste0(cellTypes, opt$predictionTask)))]
# valY = combinedVal[, (colnames(combinedVal) %in% c(paste0(cellTypes, opt$predictionTask)))]


# # Define an alpha value. 


# fitAccuracyDF = data.frame("Cell_Type" = character(),
#                             "Alpha" = double(),
#                             "Train_AUROC" = double(),
#                             "Val_AUROC" = double())

# exactPredictionsDF = data.frame("Placeholder" = rep(0, nrow(valY)))



# ###########

# set.seed(7)
# # First, loop through and get accuracy with 

# for (thisCellType in cellTypes){
#   print(paste0("Working on ", thisCellType))
#   miniY = trainY[[paste0(thisCellType, opt$predictionTask)]]
#   # trainFit = glmnet(as.matrix(trainX), miniY, family="binomial")


#   print(paste0("Working on alpha = ", opt$alphaToUse))
#   trainCVfit = cv.glmnet(as.matrix(trainX), miniY, family="binomial", alpha=opt$alphaToUse)
#   minLambda = trainCVfit$lambda.min 

#   valAUC = getAUC(trainCVfit, valX, valY[[paste0(thisCellType, opt$predictionTask)]])
#   trainAUC = getAUC(trainCVfit, trainX, miniY)

#   newRow = data.frame(thisCellType, opt$alphaToUse, trainAUC, valAUC)
#   names(newRow) = c("Cell_Type", "Alpha", "Train_AUC", "Val_AUC")
#   fitAccuracyDF = rbind(fitAccuracyDF, newRow)
#   # fitAccuracyDF[thisCellType, "CV_R_squared"] = rSquared_bestLambda

#   # # 8-31-21 edit: Add in plotting code for residuals
#   if (abs(opt$alphaToUse - .5) < .01 ){
#     valAccDF = getValAccuracyDF(trainCVfit, valX,  valY[[paste0(thisCellType, opt$predictionTask)]])
#     exactPredictionsDF[[paste0(thisCellType, "_Actual")]] = valAccDF[["Actual"]]
#     exactPredictionsDF[[paste0(thisCellType, "_Predicted")]] = valAccDF[["Predicted"]]
#   }
# }






# # Plot the residuals
# residualDir = paste0(outDirName, "residuals/")
# dir.create(residualDir)

# for (eachCelltype in cellTypes){
#   # Plot prediction vs. actual
#   png(paste0(residualDir, "Actual_vs_Pred", eachCelltype, "_",
#       opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
#       height=1000, width=1000, res=200)
#   exactPredictionsDF[[paste0(thisCellType, "_Actual")]] = as.factor(exactPredictionsDF[[paste0(thisCellType, "_Actual")]])
#   fitPlot = (ggplot(exactPredictionsDF, aes_string(x=paste0(thisCellType, "_Actual"),
#                                y=paste0(thisCellType, "_Predicted"))) +
#        ggtitle(paste0("Expression vs Prediction ", eachCelltype, "_", opt$predictionTask, "_", opt$featureSelection)) + 
#        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#        geom_violin() )
#   print(fitPlot)
#   dev.off()

#   # exactPredictionsDF[["Residuals"]] = (exactPredictionsDF[[paste0(thisCellType, "_Actual")]] - exactPredictionsDF[[paste0(thisCellType, "_Predicted")]])
#   # # Plot prediction vs. actual
#   # png(paste0(residualDir, "Actual_vs_Residuals", eachCelltype, "_",
#   #     opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
#   #     height=1000, width=1000, res=200)
#   # fitPlot = (ggplot(exactPredictionsDF, aes_string(x=paste0(thisCellType, "_Actual"),
#   #                              y="Residuals") )+
#   #      ggtitle(paste0("Expression vs Prediction ", eachCelltype, "_", opt$predictionTask, "_", opt$featureSelection)) + 
#   #      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#   #      geom_point())
#   # print(fitPlot)
#   # dev.off()
# }




# fitAccuracyDF$Alpha = as.character(fitAccuracyDF$Alpha)


# png(paste0(outDirName, "FitAUC_OnValidation_", "glmnet_varyAlpha_",
#       opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
#       height=1000, width=1000, res=200)
# fitPlot = (ggplot(fitAccuracyDF, aes_string(x="Cell_Type", y="Val_AUC", color="Alpha")) +
#      # geom_point() + 
#      ggtitle(paste0("AUC on Validation ", opt$predictionTask, "_", opt$featureSelection)) + 
#      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#      geom_jitter())
# print(fitPlot)
# dev.off()



# png(paste0(outDirName, "FitAccuracy_OnTrainSet_", "glmnet_varyAlpha_",
#       opt$predictionTask, "_from_", opt$featureSelection, "fitting.png" ),
#       height=1000, width=1000, res=200)
# fitPlot = (ggplot(fitAccuracyDF, aes_string(x="Cell_Type", y="Train_R_Squared", color="Alpha")) +
#      # geom_point() + 
#      ggtitle(paste0("Fit Accuracies by CV in glmnet ", opt$predictionTask, "_", opt$featureSelection)) + 
#      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#      geom_jitter())
# print(fitPlot)
# dev.off()
