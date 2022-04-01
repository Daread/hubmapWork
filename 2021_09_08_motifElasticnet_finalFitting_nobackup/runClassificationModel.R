
library(tidyr)
library(ggplot2)
library(glmnet)
library(ROCR)
library(data.table)

# Get shared functions between this and runRegressionModel.R
source("./utilityFunctions.R")

library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-p", "--pValFIMOcutoff"), type="numeric", 
  			default=1e-4, 
              help="Max FIMO p value to retain match", metavar="character"),

  make_option(c("-f", "--featureSelection"), type="character", 
        default="Binary_Combined_Motif_Counts", # "Binary_Combined_Motif_Counts", "Binary_PromOnly_Motif_Counts"
              help="Option of which features to use for input", metavar="character"),

  make_option(c("-t", "--predictionTask"), type="character", 
        default="_log2_CPM",   # "_log2_ratio_Vs_AllTypeMean", #  "_log2_CPM"
              help="Option of which RNA data to predict", metavar="character"),

  make_option(c("-y", "--predictionFraming"), type="character", 
        default="Regression",# "Classification",   # "Classification" or "Regression"
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
        default=1500,
              help="Bases upstream of TSS to use for input features", metavar="numeric"),
  make_option(c("-d", "--promoterDownstream"), type="numeric", 
        default=500,
              help="Bases downstream of TSS to use for input features", metavar="numeric"),

  make_option(c("-k", "--coaccessCutoff"), type="numeric", 
        default=.015,
              help="Cutoff for keeping coaccessible peaks", metavar="numeric"),

  make_option(c("-n", "--maxNdistalSites"), type="numeric", 
        default=20,
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

outDirName = paste0("./plots/", opt$predictionFraming, "/",
            "Prot_Only_Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize,"pVal", as.character(opt$pValFIMOcutoff), "/" )
dir.create(outDirName)



cellTypes = c("Adipocytes", "B_Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic_Endothelium", "Macrophage", "Mast_Cell", "Neuronal", "T_Cell",
                "Vascular_Endothelium", "VSM_and_Pericyte")
# cellTypes = c( "Macrophage", 
#                 "Vascular_Endothelium")


# Get the RNA data that'll be used in all models
rnaData = getRNAdf(opt, cellTypes)
names(rnaData)[names(rnaData) == "id"] <- "GeneID"
# Only use genes that are annotated as protein-coding
rnaData = getProteinCodingGenes(rnaData, opt)

# If classification, turn all the RNA labels into high/low labels
if (opt$predictionFraming == "Classification"){
  rnaData = binarizeRNAdata(rnaData, opt, cellTypes)
}



# Make the full, empty dataframe to hold all coefficients
outDir = paste0("./fileOutputs/", opt$predictionFraming, "/")
dir.create(outDir)
outFile = paste0(outDir,  "With_and_without_promoter_", opt$predictionFraming, "_FitModel.csv")


evalNames = c("AUC","R2")
names(evalNames) = c("Classification", "Regression")
evalCol = paste0("Eval_Set_", evalNames[opt$predictionFraming]) 


# If I've already fit and saved the results, just load that:
myRun = tryCatch({
  fullFitDF = read.csv(outFile)
  colnames(fullFitDF) = c("X", "Cell_Type", "Alpha", evalCol, "Best_Lambda", "Feature_Set")
  # evalCol = "Eval_Set_Accuracy"
  print("File read in successfully")
  },
  error=function(cond){
    print(cond)
    print("Fit not already run. Running now.")

    trainAndVal_CombinedResults = getModelFitsAndAccuracy("Binary_Combined_Motif_Counts", c("Train", "Validation"), c("Test"),
                                                             opt, rnaData, cellTypes)
    trainAndVal_PromoterResults = getModelFitsAndAccuracy("Binary_PromOnly_Motif_Counts", c("Train", "Validation"), c("Test"),
                                                             opt, rnaData, cellTypes)
    fullFitDF = rbind(trainAndVal_CombinedResults, trainAndVal_PromoterResults)

    # Save in case needed later
    write.csv(fullFitDF, outFile)
  }
  )






# Now run this
# opt$predictionFraming = "Classification" # Or "Regression"
# evalNames = c("AUC","R2")
# names(evalNames) = c("Classification", "Regression")
# evalCol = paste0("Eval_Set_", evalNames[opt$predictionFraming]) 

# trainAndVal_CombinedResults = getModelFitsAndAccuracy("Binary_Combined_Motif_Counts", c("Train", "Validation"), c("Test"),
#                                                          opt, rnaData, cellTypes)
# trainAndVal_PromoterResults = getModelFitsAndAccuracy("Binary_PromOnly_Motif_Counts", c("Train", "Validation"), c("Test"),
#                                                          opt, rnaData, cellTypes)
# fullFitDF = rbind(trainAndVal_CombinedResults, trainAndVal_PromoterResults)

# # Save in case needed later
# write.csv(fullFitDF, outFile)



# colnames(fullFitDF) = c("Cell_Type", "Alpha", "Eval_Set_R2", "Best_Lambda", "Feature_Set")


# Now plot
png(paste0("./plots/", opt$predictionFraming, "/", opt$predictionFraming, "_", opt$variableParams, "_Test_Performance_Prom_Vs_Combined.png" ),
      height=1400, width=1600, res = 200)
myPlot = ggplot(fullFitDF, aes_string(x="Cell_Type", y=evalCol, color="Feature_Set")) +
      geom_point() + ggtitle("Promoter-Only Vs. Promoter+Distal Performance") + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(myPlot)
dev.off()




png(paste0("./plots/", opt$predictionFraming, "/", opt$predictionFraming, "_", opt$variableParams, "_Barplot_Test_Performance_Prom_Vs_Combined.png" ),
      height=1400, width=2000, res = 200)
myPlot = ggplot(fullFitDF, aes_string(x="Cell_Type", y=evalCol, fill="Feature_Set")) +
      geom_bar(position='dodge', stat='identity') + # ggtitle("Promoter-Only Vs. Promoter+Distal Performance") + 
      # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
       theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      scale_fill_discrete(name = "Sequence Used", labels = c("Promoter + Distal", "Promoter Alone"))+
      scale_x_discrete(breaks=cellTypes, labels=c("Adipocyte", "B Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic Endothelium", "Macrophage", "Mast Cell", "Neuronal", "T Cell",
                "Vascular Endothelium", "Perivascular Cell")) + ylab("R^2 on Test Data") +
                 xlab("Cell Type")+theme(text=element_text(size=21))
print(myPlot)
dev.off()



formatCellType <- function(inputColumn){

  inputColumn = ifelse(inputColumn == "T_Cell", "T Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "VSM_and_Pericyte", "Perviascular Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "Vascular_Endothelium", "Vascular Endothelium", inputColumn)
  inputColumn = ifelse(inputColumn == "Lymphatic_Endothelium", "Lymphatic Endothelium", inputColumn)
  inputColumn = ifelse(inputColumn == "Mast_Cell", "Mast Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "B_Cell", "B Cell", inputColumn)

  return(inputColumn)
}


makeFitResults_vs_proportion_plot <- function(fullFitDF, opt, cellPropCol="RNA"){
  fullFitDF$Cell_Type = as.character(fullFitDF$Cell_Type)
  fullFitDF$Cell_Type = formatCellType(fullFitDF$Cell_Type)
  # For each cell type, get plots of proportion vs. 1) R^2 promoter only
                              #                     2) R^2 prom+distal/prom ratio
  # browser()

  propAndAccuracyDF = data.frame("Cell_Type" = fullFitDF$Cell_Type)
  # Get accuracies of promoters only
  promoterOnlyDF = fullFitDF[fullFitDF$Feature_Set == 'Binary_PromOnly_Motif_Counts',]
  promAccuracies = promoterOnlyDF$Eval_Set_R2
  names(promAccuracies) = promoterOnlyDF$Cell_Type
  propAndAccuracyDF[["Promoter_Only_R2"]] = promAccuracies[propAndAccuracyDF$Cell_Type]
  #... and the ratio of when distal info is added
  combinedSeqDF = fullFitDF[fullFitDF$Feature_Set == "Binary_Combined_Motif_Counts",]
  combinedAccuracies = combinedSeqDF$Eval_Set_R2
  names(combinedAccuracies) = combinedSeqDF$Cell_Type
  propAndAccuracyDF[["Promoter_Plus_Distal_R2"]] = combinedAccuracies[propAndAccuracyDF$Cell_Type]

  # Now read in cell type proportions
  cellPropFile = "../2021_07_15_Greg_ATAC_Code_nobackup/archr/results/finalPlots/RNA_vs_ATAC_CellType_Proportions.csv"
  cellPropDF = read.csv(cellPropFile)
  cellProps = cellPropDF[[paste0(cellPropCol, "_Prop")]]
  names(cellProps) = cellPropDF$Cell_Type
  propAndAccuracyDF[["Cell_Type_Proportion"]] = cellProps[propAndAccuracyDF$Cell_Type]

  propAndAccuracyDF[["Combined_Over_Promoter_Ratio"]] = propAndAccuracyDF$Promoter_Plus_Distal_R2 / propAndAccuracyDF$Promoter_Only_R2

  # Now make plots from the two desired comparisons
  plotFile = paste0("./plots/", opt$predictionFraming, "/", opt$predictionFraming, "_", 
                    opt$variableParams, "_Scatter_Promoter_Accuracy_Vs_Celltype_Prop.png" )
  png(plotFile, res=200, height = 1000, width=1000)
  myPlot = ggplot(propAndAccuracyDF, aes_string(x="Promoter_Only_R2", y="Cell_Type_Proportion")) + 
            geom_point() + 
            xlab("R^2 Using Promoter Only") + 
            ylab(paste0("Cell type proportion in ", cellPropCol)) + 
            theme(text=element_text(size=18))
  print(myPlot)
  dev.off()

  # Gain vs. abundance
  plotFile = paste0("./plots/", opt$predictionFraming, "/", opt$predictionFraming, "_", 
                    opt$variableParams, "_Scatter_Combined_Gain_Accuracy_Vs_Celltype_Prop.png" )
  png(plotFile, res=200, height = 1000, width=1000)
  myPlot = ggplot(propAndAccuracyDF, aes_string(x="Combined_Over_Promoter_Ratio", y="Cell_Type_Proportion")) + 
            geom_point() + 
            xlab("Distal+Promoter / Promoter R^2") + 
            ylab(paste0("Cell type proportion in ", cellPropCol)) + 
            theme(text=element_text(size=18))
  print(myPlot)
  dev.off()

}

# Make a plot comparing accuracy gains vs. cell type abundances
makeFitResults_vs_proportion_plot(fullFitDF, opt)









#
getModelFitsAccuracyAndModel_noCV <- function(featureSelection, setsForTraining, setsForTesting, opt, rnaData, cellTypes,
							bestLambdaList, saveModel=TRUE){
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
                            "Alpha" = double() ) #,
                            # "Feature_Set" = featureSelection)
  fitModelList = vector(mode="list", length(cellTypes))
  typeCount = 1

  set.seed(7)
  for (thisCellType in cellTypes){
    print(paste0("Working on ", thisCellType))
    miniY = trainY[[paste0(thisCellType, opt$predictionTask)]]

    # Fit the model using cross-validation to pick a best lambda
    if (opt$predictionFraming == "Classification"){
      fitModel = glmnet(as.matrix(trainX), miniY, family="binomial",
      				 alpha=opt$alphaToUse, lambda = bestLambdaList[thisCellType])
    } else if (opt$predictionFraming == "Regression"){
      fitModel = glmnet(as.matrix(trainX), miniY,
      					  alpha=opt$alphaToUse, lambda = bestLambdaList[thisCellType])
    }

    # # Classification or regression?
    # if (opt$predictionFraming == "Classification"){
    #   evalValue = getAUC(trainCVfit, valX, valY[[paste0(thisCellType, opt$predictionTask)]])
    # } else if (opt$predictionFraming == "Regression"){
    #   evalValue = getAccuracy(trainCVfit, valX, valY[[paste0(thisCellType, opt$predictionTask)]])
    # }
    
    newRow = data.frame(thisCellType,  opt$alphaToUse)
    names(newRow) = c("Cell_Type", "Alpha")
    fitAccuracyDF = rbind(fitAccuracyDF, newRow)

    # Save the model
    fitModelList[[typeCount]] = fitModel
    typeCount = typeCount + 1
  }
  # browser()
  fitAccuracyDF$Fit_Models = fitModelList
  fitAccuracyDF$Feature_Set = featureSelection
  return(fitAccuracyDF)
}


# Feed in a named list of the optimal lambda values found in the train+validation set. These will be used for the final model fitting

bestLambdaList = trainAndVal_CombinedResults$Best_Lambda
names(bestLambdaList) = trainAndVal_CombinedResults$Cell_Type


combinedModelRes = getModelFitsAccuracyAndModel_noCV("Binary_Combined_Motif_Counts", c("Train", "Validation", "Test"), c("Test"),
                                                         opt, rnaData, cellTypes, bestLambdaList)

# Make the full, empty dataframe to hold all coefficients
outDir = paste0("./fileOutputs/", opt$predictionFraming, "/")
dir.create(outDir)


allCoefDF = data.frame("Motif" = character(),
						"EmptyFit" = numeric())
for (eachCelltype in cellTypes){
  subsetResults = combinedModelRes[combinedModelRes$Cell_Type == eachCelltype,]
  thisFit = subsetResults[1,"Fit_Models"]
  # Now get the nonzero coefficients
  # coefficientDF = as.data.frame(as.matrix(coef(thisFit[[1]], s= "lambda.min")))
  coefficientDF = as.data.frame(as.matrix(coef(thisFit[[1]])))
  colnames(coefficientDF) = c("Coefficient")
  coefficientDF["Motif"] = rownames(coefficientDF)
  coefficientDF["AbsVal_Motif"] = abs(coefficientDF[["Coefficient"]])
  # Drop intercept
  coefficientDF = coefficientDF[coefficientDF$Motif != "(Intercept)",]
  # Sort
  # browser()
  coefficientDF = coefficientDF[order(-coefficientDF$AbsVal_Motif),]

  # Write this as an output alone
  outFile = paste0(outDir, eachCelltype, "_", opt$predictionFraming, "_Final_FullFit_Coefficients.csv")
  write.csv(coefficientDF, outFile)

  # Add this into the dataframe holding all results
  coefficientDF[[eachCelltype]] = coefficientDF$Coefficient
  allCoefDF = merge(allCoefDF, coefficientDF[c("Motif", eachCelltype)], by="Motif", all=TRUE)

}


# From the allCoefDF, drop the placeholder
allCoefDF = allCoefDF[, -which(names(allCoefDF) %in% c("EmptyFit"))]
# Save this to a new CSV for convenience
outFile = paste0(outDir,  "All_Celltype_Results_", opt$predictionFraming, "_Final_FullFit_Coefficients.csv")
write.csv(allCoefDF, outFile)




# Make a plot showing motifs x cell types for some high-abs motifs
makeMotifByTypePlot = function(allCoefDF, cellTypes, definedMotifs=FALSE, motifsToUse=NULL, nPerCelltype=3, plotNote="",
								minZerosToPlot=0){
	# For each celltype, get n top hits to show
	if (definedMotifs == FALSE){
		print("Finding top motifs by cell type")
		hitMotifs = c()
		# Filter down to rows with enough zeros
		allCoefDF = allCoefDF[rowSums(allCoefDF == 0) >= minZerosToPlot,]

		# Get the top nPerCellType TF motifs by coefficient magnitude
		for (eachCelltype in cellTypes){
			miniDF = data.frame(theseCoefs = allCoefDF[[eachCelltype]],
								motifs = allCoefDF$Motif)
			miniDF$absCoef = abs(miniDF$theseCoefs)
			# Reorder
			miniDF = miniDF[order(miniDF$absCoef, decreasing=TRUE),]
			# Get the first few coefs
			sortedMotifs = miniDF$motifs 
			hitMotifs = c(hitMotifs, sortedMotifs[1:nPerCelltype])
		}
		# Now get the unique ones
		uniqueHitMotifs = unique(hitMotifs)
		fileName = paste0("./plots/", opt$predictionFraming, "/", plotNote, "_", opt$predictionFraming, "_", minZerosToPlot, "zeros_",
						 "_n", nPerCelltype, "_", opt$variableParams, "_MotifCoefDotplot.png")
	} else {
		fileName =paste0("./plots/", opt$predictionFraming, "/", plotNote, "_", opt$predictionFraming,  opt$variableParams, "_MotifCoefDotplot.png")
		uniqueHitMotifs = motifsToUse
	}
	

	# Get a formatted dataframe
	meltedDF = data.table::melt(allCoefDF, id.vars=c("Motif"), variable.name="Cell_Type", value.name="Coefficient")
	meltedDF$absCoefficient = abs(meltedDF$Coefficient)
	meltedDF$Coef_Sign = ifelse(meltedDF$Coefficient > 0, "Positive", "Negative")

	# Now make a plot of coefficients
	png(fileName,
			res=100, width=1500, height=1200)
	myPlot = ggplot(meltedDF[meltedDF$Motif %in% uniqueHitMotifs,],
						 aes(x=Cell_Type, y=Motif, 
						 	col=Coef_Sign, size=ifelse(absCoefficient==0, NA, absCoefficient))) +
						geom_point()+ 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	print(myPlot)
	dev.off()

}


makeMotifByTypePlot(allCoefDF, cellTypes, nPerCelltype=10, minZerosToPlot=2)
makeMotifByTypePlot(allCoefDF, cellTypes, nPerCelltype=10, minZerosToPlot=5)

makeMotifByTypePlot(allCoefDF, cellTypes)


makeMotifByTypePlot(allCoefDF, cellTypes, nPerCelltype=5)


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("Smad2..Smad3_AllSeq", "SMAD2..SMAD3..SMAD4_AllSeq",
																		    "SMAD3_AllSeq", "Smad4_AllSeq", "FOXH1_AllSeq",
																		    "FOS_AllSeq", "FOS..JUN_AllSeq", "FOS..JUN.var.2._AllSeq",
																		    "JUN_AllSeq", "SNAI1_AllSeq"), plotNote="TGFBeta_Set")

makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("GATA4_AllSeq", "GATA5_AllSeq", "GATA6_AllSeq", "MEF2_AllSeq",
																"FOXO_AllSeq", "NKX2.5_AllSeq", "YY1_AllSeq", "HEY2_AllSeq",
																"MITF_AllSeq"), plotNote="Hypertrophy_Set")



makeMotifByTypePlot(allCoefDF, cellTypes, nPerCelltype=30)
























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
