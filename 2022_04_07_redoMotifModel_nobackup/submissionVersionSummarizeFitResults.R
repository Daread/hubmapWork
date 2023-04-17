
library(tidyr)
library(ggplot2)
library(glmnet)
library(ROCR)

library("optparse")
# Get the passed parameters
option_list = list(

  make_option(c("-f", "--featureSelection"), type="character", 
        default="Binary_Combined_Motif_Counts", # "Binary_Combined_Motif_Counts", "Binary_PromOnly_Motif_Counts"
              help="Option of which features to use for input", metavar="character"),

  make_option(c("-t", "--predictionTask"), type="character", 
        default="_log2_CPM",   # "_log2_ratio_Vs_AllTypeMean", #  "_log2_CPM"
              help="Option of which RNA data to predict", metavar="character"),

  make_option(c("-y", "--predictionFraming"), type="character", 
        default= "regression", #"classification",   # "classification" or "regression"
              help="Option of which RNA data to predict", metavar="character"),

  make_option(c("-s", "--summarizationPolicy"), type="character", 
        default = "alpha0.5_average", #"alpha0.5_max", # "alpha0.5_average",   # "alpha0.5_max" "alpha0.5_average"
              help="How to quantify the accuracy of a whole run (which comprises many alpha vals and cell types)", metavar="character"),

  make_option(c("-p", "--paramFile"), type="character", 
        default="VaryPromotersAllRunsMade",   # "VaryPromotersAllRunsMade.csv",  "VaryDistalParams",  #"VaryPromoters",   #"VaryDistalParams",   # "VaryPromoters" or "VaryDistalParams"
              help="File of which parameter sets to summarize", metavar="character")

)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$variableParams = paste0(opt$promoterUpstream, "_", opt$promoterDownstream, "_",
                           opt$coaccessCutoff, "_", opt$maxNdistalSites, "_", opt$peakSize )

paramFileName = opt$paramFile
paramSetsToPlotFile = paste0("./fileOutputs/", paramFileName, ".csv")

outDirName = paste0("./plots/summaryPlotsResubmitCheck/")
dir.create(outDirName)

# Read in the sets of hyperparameters we want to summarize
parameterSets = read.csv(paramSetsToPlotFile)


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



getFitDF <- function(inputParams, opt){
  # Read the file
  subDir = paste0(  "Prot_Only_Gene_Prom_Plus_Distal_WithSequence_Sites_Max", inputParams[1,"maxNsites"], "_Upstream",inputParams[1,"Upstream"],
                "_Downstream", inputParams[1,"Downstream"], "_cicCuf", inputParams[1,"coaccessCutoff"],
                    "peakSize", inputParams[1,"peakSize"],"pVal", as.character(inputParams[1,"pVal"]), "/",
                    "fitVals", opt$predictionTask, "_from_", opt$featureSelection, ".csv" )
  thisResultFile = paste0("./plots/", opt$predictionFraming, "/", subDir )
  thisFitResult = read.csv(thisResultFile)
  return(thisFitResult)
}

# I want to condense down a whole lot of run information into a single number.
# Otherwise I'll end up with annoyances that some parameter sets are ever-so-slightly better in one cell type than another
# One policy is to simply grab the average accuracy value (AUC or R^2) across all fits, all fits at a given alpha, etc
getSummarizedDF <- function(inputDF, opt){
  evalNames = c("Val_AUC","Val_R_Squared")
  names(evalNames) = c("classification", "regression")

  # Run the policy to summarize
  if (opt$summarizationPolicy == "alpha0.5_average"){
    subsetDF = inputDF[inputDF$Alpha == 0.5,]
    # Average these runs
    # browser()
    thisSummary = mean(subsetDF[[evalNames[opt$predictionFraming]]])
  }

  # Run the policy to summarize
  if (opt$summarizationPolicy == "alpha0.5_max"){
    subsetDF = inputDF[inputDF$Alpha == 0.5,]
    # Average these runs
    # browser()
    thisSummary = max(subsetDF[[evalNames[opt$predictionFraming]]])
  }

  # Return
  return(thisSummary)
}




dfSummaries = data.frame(matrix(ncol=(length(colnames(parameterSets)) + 1),nrow=0))
colnames(dfSummaries) = c(opt$summarizationPolicy, colnames(parameterSets))

for (eachRowNum in 1:nrow(parameterSets)){
  thisParamDF = parameterSets[eachRowNum,]
  fitDF = getFitDF(thisParamDF, opt)
  # browser()
  # Now, summarize this fitting acrossGroups
  thisParamDF[[opt$summarizationPolicy]] = getSummarizedDF(fitDF, opt)
  # Add this to the summary data frame
  dfSummaries = rbind(dfSummaries, thisParamDF)
}





# Plot, based on the input file used
makePlotThisComparison <- function(paramFileName, dfSummaries, opt){

  if (paramFileName == "VaryPromoters") {
    # Show varied upstrea/downstream sizes and a point per alpha
    dfSummaries$Parameter_Set = paste0(dfSummaries$Upstream, "/", dfSummaries$Downstream)
    dfSummaries$pVal = as.character(dfSummaries$pVal)

    if (opt$summarizationPolicy == "alpha0.5_max"){
      yLabToUse = "Highest R^2 Obtained"
    } else {
      yLabToUse = "Average R^2 Obtained"
    }

    png(paste0("./plots/summaryPlots/", opt$predictionFraming,
           "_", paramFileName, "_", opt$summarizationPolicy, ".png"),
        height=1200, width=1400,res=200)
    myPlot = ggplot(dfSummaries, aes_string(x="Parameter_Set", y=opt$summarizationPolicy, color="pVal") ) + 
          geom_point() + 
          xlab("TSS Up/Downstream Lengths") + #+ ggtitle(paste0(paramFileName, "_", opt$summarizationPolicy)) +
            theme(text = element_text(size = 20)) + 
              ylab(yLabToUse) + 
               guides(color=guide_legend(title="Motif p Value"))+ 
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
    print(myPlot)
    dev.off()
  }

  if (paramFileName == "VaryPromotersAllRunsMade") {
    # Show varied upstrea/downstream sizes and a point per alpha
    dfSummaries$Parameter_Set = paste0(dfSummaries$Upstream, "/", dfSummaries$Downstream)
    dfSummaries$pVal = as.character(dfSummaries$pVal)

    if (opt$summarizationPolicy == "alpha0.5_max"){
      yLabToUse = "Highest R^2 Obtained"
    } else {
      yLabToUse = "Average R^2 Obtained"
    }

    png(paste0("./plots/summaryPlots/", opt$predictionFraming,
           "_", paramFileName, "_", opt$summarizationPolicy, ".png"),
        height=1200, width=1400,res=200)
    myPlot = ggplot(dfSummaries, aes_string(x="Parameter_Set", y=opt$summarizationPolicy, color="pVal") ) + 
          geom_point() + 
          xlab("TSS Up/Downstream Lengths") + #+ ggtitle(paste0(paramFileName, "_", opt$summarizationPolicy)) +
            theme(text = element_text(size = 20)) + 
              ylab(yLabToUse) + 
               guides(color=guide_legend(title="Motif p Value"))+ 
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
    print(myPlot)
    dev.off()
  }


  if (paramFileName %in% c("VaryDistalAllRunsMade", "VaryDistalParams")) {
    # Show varied upstrea/downstream sizes and a point per alpha
    dfSummaries$Parameter_Set = paste0(dfSummaries$coaccessCutoff, "/", dfSummaries$maxNsites, "/", dfSummaries$peakSize)
    dfSummaries$pVal = as.character(dfSummaries$pVal)

    png(paste0("./plots/summaryPlots/", opt$predictionFraming,
           "_", paramFileName, "_", opt$summarizationPolicy, ".png"),
        height=1000, width=1600,res=200)
    myPlot = ggplot(dfSummaries, aes_string(x="Parameter_Set", y=opt$summarizationPolicy, color="pVal") ) + 
          geom_point() + 
          xlab("Cicero Cutoff/Max Distal Sites Linked/Distal Site Size in bp") + ylab("R^2 Average") +
          #ggtitle(paste0(paramFileName, "_", opt$summarizationPolicy))+ 
          guides(color=guide_legend(title="Motif p Value"))+ 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    print(myPlot)
    dev.off()
  }
}


print("Plotting Comparison Now")
# Plot, based on the input file used
makePlotThisComparison(paramFileName, dfSummaries, opt)







dfFullData = data.frame(matrix(ncol=(length(colnames(parameterSets)) + ncol(fitDF)),nrow=0))
colnames(dfFullData) = c(colnames(fitDF), colnames(parameterSets))

for (eachRowNum in 1:nrow(parameterSets)){
  thisParamDF = parameterSets[eachRowNum,]
  fitDF = getFitDF(thisParamDF, opt)

  fitDFparams = data.frame(do.call('rbind',strsplit(as.character(fitDF$hyperParams), "_", fixed=TRUE)))
  colnames(fitDFparams) = c("Upstream", "Downstream", "coaccessCutoff", "maxNsites", "peakSize")
  fitDF = cbind(fitDF, fitDFparams)
  fitDF$pVal = thisParamDF$pVal

  # Now, summarize this fitting acrossGroups
  # thisParamDF[[opt$summarizationPolicy]] = getSummarizedDF(fitDF, opt)
  # Add this to the summary data frame
  dfFullData = rbind(dfFullData, fitDF)
}

dfFullData = dfFullData[dfFullData$Alpha == 0.5,]


# fitDFparams = data.frame(do.call('rbind',strsplit(as.character(fitDF$hyperParams), "_", fixed=TRUE)))
# colnames(fitDFparams) = c("Upstream", "Downstream", "coaccessCutoff", "maxNsites", "pVal")

# fitDF = cbind(fitDF, fitDFparams)
# # fitDFparams = data.frame(do.call('rbind',strsplit(as.character(fitDF$hyperParams), "_", fixed=TRUE)))





# Plot, based on the input file used
makePlotUnsummarizedThisComparison <- function(paramFileName, fitDF, opt){

  if (paramFileName == "VaryPromoters") {
    # Show varied upstrea/downstream sizes and a point per alpha
    fitDF$Parameter_Set = paste0(fitDF$Upstream, "/", fitDF$Downstream)
    fitDF$pVal = as.character(fitDF$pVal)

    yLabToUse = "Validation R^2"

    png(paste0("./plots/summaryPlots/ShowCellTypes_", opt$predictionFraming,
           "_", paramFileName,  ".png"),
        height=1200, width=1400,res=200)
    myPlot = ggplot(fitDF, aes_string(x="Parameter_Set", y="Val_R_Squared") ) + 
          geom_boxplot(aes(fill=pVal)) + 
          geom_point(position=position_dodge(width=0.75), aes(group=pVal)) +
          monocle_theme_opts() + 
          xlab("TSS Up/Downstream Lengths") + #+ ggtitle(paste0(paramFileName, "_", opt$summarizationPolicy)) +
            theme(text = element_text(size = 20)) + 
              ylab(yLabToUse) + 
               guides(fill=guide_legend(title="Motif p Value"))+ 
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
    print(myPlot)
    dev.off()
  }

  if (paramFileName == "VaryPromotersAllRunsMade") {
    # Show varied upstrea/downstream sizes and a point per alpha
    fitDF$Parameter_Set = paste0(fitDF$Upstream, "/", fitDF$Downstream)
    fitDF$pVal = as.character(fitDF$pVal)

    yLabToUse = "Validation R^2"

    png(paste0("./plots/summaryPlots/ShowCellTypes_", opt$predictionFraming,
           "_", paramFileName, "_", opt$featureSelection,  ".png"),
        height=1800, width=2100,res=300)
    myPlot = ggplot(fitDF, aes_string(x="Parameter_Set", y="Val_R_Squared") ) + 
          geom_boxplot(aes(fill=pVal)) + 
          geom_point(position=position_dodge(width=0.75), aes(group=pVal)) +
          monocle_theme_opts() + 
          xlab("TSS Up/Downstream Lengths") + #+ ggtitle(paste0(paramFileName, "_", opt$summarizationPolicy)) +
            theme(text = element_text(size = 20)) + 
              ylab(yLabToUse) + 
               guides(fill=guide_legend(title="Motif p Value"))+ 
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
          scale_fill_brewer(palette="Set1")
    print(myPlot)
    dev.off()
  }

  if (paramFileName == "VaryDistalParams") {
    # Show varied upstrea/downstream sizes and a point per alpha
    fitDF$Parameter_Set = paste0(fitDF$coaccessCutoff, "/", fitDF$maxNsites, "/", fitDF$peakSize)
    fitDF$pVal = as.character(fitDF$pVal)

    yLabToUse = "Validation R^2"

    png(paste0("./plots/summaryPlots/ShowCellTypes_", opt$predictionFraming,
           "_", paramFileName,  ".png"),
        height=1200, width=1800,res=200)
    myPlot = ggplot(fitDF, aes_string(x="Parameter_Set", y="Val_R_Squared") ) + 
          geom_boxplot(aes(fill=pVal)) + 
          geom_point(position=position_dodge(width=0.75), aes(group=pVal)) +
          monocle_theme_opts() + 
          xlab("Cicero Cutoff/Max Sites Linked/Distal Site Size in BP") + #+ ggtitle(paste0(paramFileName, "_", opt$summarizationPolicy)) +
            theme(text = element_text(size = 20)) + 
              ylab(yLabToUse) + 
               guides(fill=guide_legend(title="Motif p Value"))+ 
          theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
          scale_fill_brewer(palette="Set1")
    print(myPlot)
    dev.off()
  }
}


# Plot, based on the input file used
makePlotUnsummarizedThisComparison(paramFileName, dfFullData, opt)
