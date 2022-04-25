
getInputDF <- function(opt, specifyFeatures=FALSE, specifiedFeatures="NULL"){
  # Get the path
  dirName = paste0("../2022_04_07_redoATACtoExprData_nobackup/fileOutputs/", 
            "Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize, "/" )
  # dirName = paste0("../2022_04_07_redoATACtoExprData_nobackup/fileOutputs/backupAllOldIntermediatesBeforeFixRaceConditionBug/", 
  #           "Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
  #               "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
  #                   "peakSize", opt$peakSize, "/" )

  # Get previous output
  dfName = paste0("Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize, "pVal", as.character(opt$pValFIMOcutoff))

  # Read in according to the featureSelection specified
  if (specifyFeatures){
    featureDF = read.csv(paste0(dirName, dfName, "_", specifiedFeatures, ".csv"))
  } else {
    featureDF = read.csv(paste0(dirName, dfName, "_", opt$featureSelection, ".csv"))
  }
  

  return(featureDF)
}

getRNAdf <- function(opt, cellTypes){
  rnaDir = paste0("../2022_04_07_redoATACtoExprData_nobackup/rdsOutputs/RNA_Data/" )
  rnaDF = readRDS(paste0(rnaDir, "Pseudobulked_RNA_countsPerMillionNotDeduplicated.RDS"))

  print(paste0(cellTypes, opt$predictionTask))
  returnDF = rnaDF[,c("id", "gene_short_name", 
    paste0(cellTypes, opt$predictionTask))]

  return(returnDF)
}




getProteinCodingGenes <- function(rnaData, opt){
  # Get the gtf file
  gtfData = read.table((opt$gtfFile), sep="\t")
  splitData = separate(gtfData, "V9", c("geneTxt", "idTxt", "ensemblName"))

  # Get unique names
  ensemblIDsProtein = unique(splitData[["ensemblName"]])
  protRNA = rnaData[rnaData$GeneID %in% ensemblIDsProtein,]
  return(protRNA)
}


binarizeRNAdata <- function(rnaData, opt, cellTypes){
  # Loop through each of the entries
  for (eachCol in paste0(cellTypes, opt$predictionTask)){
    rnaData[[eachCol]] = ifelse(rnaData[[eachCol]] > opt$highCutoff, 1, 
                    ifelse(rnaData[[eachCol]] < opt$lowCutoff, 0, -1))
  }

  # Return the values
  return(rnaData)
}



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

getSubsetOfCombined <- function(combinedDF, opt, setsToUse){
  combinedSubset = combinedDF[combinedDF$Model_Set %in% setsToUse,]
  combinedSubset = combinedSubset[, !(colnames(combinedSubset) %in% c("Model_Set", "GeneID", "gene_short_name"))]

  return(combinedSubset)
}

getAccuracy <-function(fitCV_Model, inputX, inputY){
  # Get the prediction
  predictedVal = predict(fitCV_Model, newx=as.matrix(inputX), s = "lambda.min")
  # Get the r_squared
  rSquared = (cor(predictedVal[,1], inputY) ^ 2)
  return(rSquared)
}

getValAccuracyDF <-function(fitCV_Model, inputX, inputY){
  # Get the prediction
  predictedVal = predict(fitCV_Model, newx=as.matrix(inputX), s = "lambda.min")
  # Make a DF out of this
  predictDF = data.frame("Actual" = inputY, "Predicted" = predictedVal[,1])
  # browser()
  return(predictDF)
}


getAUC <- function(fitCV_Model, inputX, inputY){
  # Get the prediction
  predictedVal = predict(fitCV_Model, newx=as.matrix(inputX), s = "lambda.min")
  predROCR = prediction(predictedVal, inputY)
  aucValue = performance(predROCR, measure="auc")@y.values[[1]]

  return(aucValue)
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
  evalNames = c("AUC","R2")
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



   
monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}


plot_genes_violin_DFR_customPseudo <- function (cds_subset,
                               group_cells_by = NULL,
                               min_expr = 0,
                               nrow = NULL,
                               ncol = 1,
                               panel_order = NULL,
                               label_by_short_name = TRUE,
                               normalize = TRUE,
                               log_scale = TRUE,
                               pseudocount = 0, 
                               leaveFillBlank=FALSE) {

  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))

  if(!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% names(colData(cds_subset)),
                            msg = paste("group_cells_by must be a column in",
                                        "the colData table"))
  }

  assertthat::assert_that(assertthat::is.number(min_expr))

  if(!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }

  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(assertthat::is.number(pseudocount))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)),
                            msg = paste("When label_by_short_name = TRUE,",
                                        "rowData must have a column of gene",
                                        "names called gene_short_name."))
  }
  if(!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in%
                                    rowData(cds_subset)$gene_short_name))
    } else {
      assertthat::assert_that(all(panel_order %in%
                                    row.names(rowData(cds_subset))))
    }
  }

  assertthat::assert_that(is.logical(normalize))
  assertthat::assert_that(is.logical(log_scale))


  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100,
                          msg = paste("cds_subset has more than 100 genes -",
                                      "pass only the subset of the CDS to be",
                                      "plotted."))

  if (pseudocount > 0) {
    cds_exprs <- SingleCellExperiment::counts(cds_subset) + pseudocount
  } else {
    cds_exprs <- SingleCellExperiment::counts(cds_subset)
  }
  if (normalize) {
    cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  } else {
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  }

  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr

  cds_exprs <- merge(cds_exprs, rowData(cds_subset), by.x = "f_id",
                     by.y = "row.names")
  # cds_exprs <- merge(cds_exprs, colData(cds_subset), by.x = "Cell",  by.y = "row.names")
  # browser()
  # row.names(x)
  # cds_exprs$Cell = as.character(cds_exprs$Cell)
  cds_exprs <- merge(as.data.frame(cds_exprs), as.data.frame(colData(cds_subset)),
         by.x = "Cell",  by.y = "row.names")

  if (label_by_short_name) {
    if (!is.null(cds_exprs$gene_short_name)) {
      cds_exprs$feature_label <- cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    } else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  } else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }

  if (!is.null(panel_order)) {
    cds_exprs$feature_label = factor(cds_exprs$feature_label,
                                     levels = panel_order)
  }

  cds_exprs[,group_cells_by] <- as.factor(cds_exprs[,group_cells_by])

  q <- ggplot(aes_string(x = group_cells_by, y = "expression"),
              data = cds_exprs) +
    monocle_theme_opts()

  cds_exprs[,group_cells_by] <- as.factor(cds_exprs[,group_cells_by])
  if (!leaveFillBlank){
    q <- q + geom_violin(aes_string(fill = group_cells_by), scale="width") +
    guides(fill=FALSE)
    } else{
      q <- q + geom_violin(scale="width") +
           guides(fill=FALSE)
    }
  

  q <- q + stat_summary(fun=mean, geom="point", size=1, color="black")
  q <- q + facet_wrap(~feature_label, nrow = nrow,
                      ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }

  q <- q + ylab("Expression") + xlab(group_cells_by)

  if (log_scale){
    q <- q + scale_y_log10()
  }
  q
}




