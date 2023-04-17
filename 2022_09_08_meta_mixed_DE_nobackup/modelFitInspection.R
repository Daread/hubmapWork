
library(ggplot2)
library(dplyr)
library("optparse")
library(stringr)
library(ggpubr)

print("Libraries loaded, starting now")

getFitResults <- function(opt){
	rdsPath = "./rdsOutput/nebulaFits/"

	# Loop through each of the splits
	resultDF = data.frame()
  covarDF = data.frame()
	for (eachSplit in 1:(opt$nSplits)){
		inFile = paste0(rdsPath, opt$cellType, opt$modelNotes, as.character(eachSplit), "of", as.character(opt$nSplits), "_MMresult.rds")
		thisRes = readRDS(inFile)
		# Get the first list entry, holding model coefficients/pvalues
		resultDF = rbind(resultDF, thisRes[[1]])
    # Also get the covariance matrix information
    thisCovarDF = thisRes[[5]]
    rownames(thisCovarDF) = thisRes[[1]]$gene
    covarDF = rbind(covarDF, thisCovarDF)
    # browser()
	}
	return(list(resultDF, covarDF))
}

formatCovarRow = function(x, fixedEffCount = 12){
  # From Nebula github: https://github.com/lhe17/nebula
  cov= matrix(NA,fixedEffCount,fixedEffCount)
  cov[lower.tri(cov,diag=T)] = as.numeric(x)
  cov[upper.tri(cov)] = t(cov)[upper.tri(cov)]

  # Return
  return(cov)
}

getCovarMatDF <- function(covarInput, fixedCount, formatCovarRow){
  # browser()
  #
  covarDF = apply(covarInput, 1, formatCovarRow, fixedEffCount = fixedCount)

  return(covarDF)
}

getCovarNames <-function(fitResults, fixedCovariateCount){
  covarNames = colnames(fitResults)[1:fixedCovariateCount]
  # Format to remove logFC_
  covarNames = sub("logFC_", "", covarNames)
  return(covarNames)
}

makeVarCovarHists <- function(inputCovarDF, covarNames){

  covarCount = length(covarNames)
  histList = vector(mode="list", length=covarCount*covarCount)
  dataIndex = 1
  for (yPos in 1:covarCount){
    for (xPos in 1:covarCount){
      # Label this
      if (yPos == xPos){
        thisLabel = paste0(covarNames[yPos], " Variance")
      } else{
        thisLabel = paste0("Cov ", covarNames[yPos], " w/ ", covarNames[xPos])
      }



      miniDF = data.frame("Data" = log10(as.numeric(inputCovarDF[dataIndex,])))

      # Get this hist
      thisHist = (ggplot(miniDF, aes(x=Data)) + geom_histogram() 
          + xlab(thisLabel) + xlim(-10,10))
      histList[[dataIndex]] = thisHist
      dataIndex = dataIndex + 1
    }
  }

  print("Plotting hist now")
  # browser()
  # Plot the multi-panel plot
  png(paste0("./plots/", "CovarHistSingleIntVar.png"), res=200, height = 2400, width=2400)
  print(histList[[1]])
  dev.off()

  scaleVar = 3.5
  # Plot the multi-panel plot
  png(paste0("./plots/", "CovarHistPanel.png"), res=200, height = 2400*scaleVar, width=2400*scaleVar)
  arrangedData = ggarrange(plotlist = histList, ncol=covarCount, nrow=covarCount)
  print(arrangedData)
  dev.off()
}


# Get the passed parameters
option_list = list(
  make_option(c("-m", "--modelNotes"), type="character", 
  			# default="_fix_Anatomical_Site_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult", 
  			default="_fix_Anatomical_Site,Age,Sex,DataSource,log10_umi_rand_Donor_NucleiOnlySharedGenesCDS_noAtr_",
              help="Processing note from model fitting", metavar="character"),
  make_option(c("-c", "--cellType"), type="character", 
  			# default="Endocardium", 
  			default="Endothelium",
              help="Cell type for which the model was fit", metavar="character"),
    make_option(c("-n", "--nSplits"), type="numeric", default = 10,
              help="Number of sub-jobs to split the DE testing task into", metavar="numeric"),
    make_option(c("-v", "--covariates"), type="character", 
  			# default="Endocardium", 
  			default="SexM,Age",
              help="Covariates to summarize", metavar="character")
    )
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

dir.create("./plots/")

cellType = opt$cellType
modelFitDescription = opt$modelNotes

# Read in the fit
fitResultsList = getFitResults(opt)
fitResults = fitResultsList[[1]]
covarResults = fitResultsList[[2]]

# Nebula outputs only the non-redundant parts of the matrix, by default. Format it
fixedCovariateCount = sum(str_detect(colnames(fitResults), "logFC"))
covarMatDF = getCovarMatDF(covarResults, fixedCovariateCount, formatCovarRow)


covarNames = getCovarNames(fitResults, fixedCovariateCount)

# Now plot the histogram of variance/covariance
makeVarCovarHists(covarMatDF, covarNames)























