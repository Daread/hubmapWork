
library(tidyr)
library(ggplot2)
library(glmnet)
library(ROCR)
library(data.table)
library(viridis)

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


# Get the motif information
getMotifMatrix <- function(featureSelection, setsToGet, opt, rnaData, cellTypes){
  opt$featureSelection = featureSelection
  # Get the correct input features
  inputFeatures = getInputDF(opt, specifyFeatures=TRUE, specifiedFeatures = featureSelection)  # (opt, specifyFeatures=FALSE, specifiedFeatures="NULL")
  combinedDF = getCombinedDF(inputFeatures, rnaData, opt, cellTypes)

  # Now get the subset that gets used for training.
  combinedData = getSubsetOfCombined(combinedDF, opt, setsToGet)

  motifData = combinedData[, !(colnames(combinedData) %in% c(paste0(cellTypes, opt$predictionTask)))]

  return(motifData)
}



combinedSeq_allSets = getMotifMatrix("Binary_Combined_Motif_Counts", c("Train", "Validation", "Test"),
                                                         opt, rnaData, cellTypes)

motifCorrelationMatrix = cor(combinedSeq_allSets, method="pearson")



# Read in from the csv where outputs were made at the end of final fitting
outDir = paste0("./fileOutputs/", opt$predictionFraming, "/")
inputFile = paste0(outDir,  "All_Celltype_Results_", opt$predictionFraming, "_Final_FullFit_Coefficients.csv")
allCoefDF = read.csv(inputFile)
# Drop row numbers
allCoefDF = allCoefDF[,-which(colnames(allCoefDF) == "X")]


corrPlotDir = "./plots/correlationMatModified/"
dir.create(corrPlotDir)

getCorrModifiedDF <- function(allCoefDF, motifCorrelationMatrix, cellTypes, normalizeCorr=FALSE,
                    cutoffVal=0.0){
  # Reduce values under a cutoff to zero
  motifCorrelationMatrix[motifCorrelationMatrix < cutoffVal] = 0

  # Normalize the correlation matrix?
  if (normalizeCorr){
    motifCorrelationMatrix = apply(motifCorrelationMatrix, 1, function(x)(x/sum(x)))
  }
  # browser()

  # Loop and multiply
  for (eachCelltype in cellTypes){
    allCoefDF[eachCelltype] = as.numeric(motifCorrelationMatrix %*% allCoefDF[[eachCelltype]])
  }

  return(allCoefDF)
}


corrMultipliedDF = getCorrModifiedDF(allCoefDF, motifCorrelationMatrix, cellTypes)
normdCorrMultDF = getCorrModifiedDF(allCoefDF, motifCorrelationMatrix, cellTypes, normalizeCorr = TRUE)


normdCorrMultDF_highCorrOnly = getCorrModifiedDF(allCoefDF, motifCorrelationMatrix, cellTypes, 
                                        normalizeCorr = TRUE, cutoffVal = .5)


# Show the correlation matrix itself
meltedCorr = melt(as.data.frame(motifCorrelationMatrix))
png(paste0(corrPlotDir, "Correlation_Histogram.png"), res=200, height=1000, width=1000)
myPlot = ggplot(meltedCorr[meltedCorr$value < .9999999,], aes(x=value)) +
               geom_histogram() +
             # ggtitle("Correlations of motifs") +
             theme(text=element_text(size=24)) + xlab("Correlation") +
             ylab("Motif Pair Count")
print(myPlot)
dev.off()

# Make a heatmap of the correlation matrix
meltedMatrix = melt(motifCorrelationMatrix)
# Get the order for clustering

set.seed(7)
hclustRes = hclust(dist(as.data.frame(motifCorrelationMatrix), method="euclidean"), method="ward.D")
motifOrder = hclustRes$order

# # Order factors
# meltedMatrix$Var1 = as.factor(sapply(strsplit(as.character(meltedMatrix$Var1), "_"), `[`, 1))
# meltedMatrix$Var2 = as.factor(sapply(strsplit(as.character(meltedMatrix$Var2), "_"), `[`, 1))
# meltedMatrix$Var1 <- factor(meltedMatrix$Var1, levels = as.character(meltedMatrix$Var1)[motifOrder])
# meltedMatrix$Var2 <- factor(meltedMatrix$Var2, levels = as.character(meltedMatrix$Var12)[motifOrder])


meltedMatrix$Var1 <- factor(meltedMatrix$Var1, levels = colnames(motifCorrelationMatrix)[motifOrder])
meltedMatrix$Var2 <- factor(meltedMatrix$Var2, levels = colnames(motifCorrelationMatrix)[motifOrder])

png(paste0(corrPlotDir, "Correlation_Heatmap.png"), res=200, height=14000, width=14000)
myPlot = ggplot(meltedMatrix, aes(x=Var1, y=Var2, fill=value)) + geom_tile() + ggtitle("Correlations of motifs") +
      scale_fill_viridis() + 
      theme(axis.text.x=element_text(angle=45, hjust=1))
print(myPlot)
dev.off()

png(paste0(corrPlotDir, "Zerod_Correlation_Heatmap.png"), res=200, height=14000, width=14000)
myPlot = ggplot(meltedMatrix, aes(x=Var1, y=Var2, fill=ifelse(value<.5, 0, value))) +
     geom_tile() + ggtitle("Correlations of motifs") +
      scale_fill_viridis() + 
      theme(axis.text.x=element_text(angle=45, hjust=1))
print(myPlot)
dev.off()


# Loop, and for each cell type make
# 1: A histogram of the coefficients
# 2: A histogram of the cofficients after multiplying by the correlation matrix
for (eachCelltype in cellTypes){

  # Coefficients as is
  png(paste0(corrPlotDir, eachCelltype, "_coefficients.png"), res=200, height=800, width=1000)
  myPlot = ggplot(allCoefDF, aes_string(x=eachCelltype)) + geom_histogram() +
          ggtitle(paste0(eachCelltype, " Coefficients"))
  print(myPlot)
  dev.off()

 # Coefficients multiplied by unscaled correlation
  png(paste0(corrPlotDir, eachCelltype, "_unscaleMult_coeffs.png"), res=200, height=800, width=1000)
  myPlot = ggplot(corrMultipliedDF, aes_string(x=eachCelltype)) + geom_histogram() +
          ggtitle(paste0(eachCelltype, " Coeffs mult by Unscaled Corr"))
  print(myPlot)
  dev.off()


 # Coefficients multiplied by scaled correlation
  png(paste0(corrPlotDir, eachCelltype, "_scaledMult_coeffs.png"), res=200, height=800, width=1000)
  myPlot = ggplot(normdCorrMultDF, aes_string(x=eachCelltype)) + geom_histogram() +
          ggtitle(paste0(eachCelltype, " Coeffs mult by Scaled Corr"))
  print(myPlot)
  dev.off()

   # Coefficients multiplied by scaled correlation and only using corr > .5
  png(paste0(corrPlotDir, eachCelltype, "_cutoff_0.5_scaledMult_coeffs.png"), res=200, height=800, width=1000)
  myPlot = ggplot(normdCorrMultDF_highCorrOnly, aes_string(x=eachCelltype)) + geom_histogram() +
          ggtitle(paste0(eachCelltype, " Coeffs mult by Scaled Corr"))
  print(myPlot)
  dev.off()



}

















matrixPartOnly = normdCorrMultDF_highCorrOnly[,-which(colnames(normdCorrMultDF_highCorrOnly) == "Motif")]

matrixPartOnly = as.matrix(matrixPartOnly)











