
getInputDF <- function(opt, specifyFeatures=FALSE, specifiedFeatures="NULL"){
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
  if (specifyFeatures){
    featureDF = read.csv(paste0(dirName, dfName, "_", specifiedFeatures, ".csv"))
  } else {
    featureDF = read.csv(paste0(dirName, dfName, "_", opt$featureSelection, ".csv"))
  }
  

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

