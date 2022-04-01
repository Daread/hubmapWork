#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)

print("Starting Now")
library("SymSim")
source("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_29_simulations_and_power_nobackup/simFuncsDavidMod.R")
library(distr)
# source("/net/trapnell/vol1/home/readdf/trapLabDir/fibrosis/results/HeirarchalDE_Methods.R")
library("optparse")
library(ggplot2)

outDir = "./plots/downsampleTests/"
dir.create(outDir)

# # Get input params
# paramList = processParamArgs()
# print(paramList)


# Get the passed parameters
option_list = list(
   make_option(c("-c", "--cells"), type="numeric", 
  			default=300, 
              help="Number of cells to simulate", metavar="numeric"),
  make_option(c("-i", "--indiv"), type="numeric", 
  			default=3, 
              help="Number of individuals per group", metavar="numeric"),
  make_option(c("-g", "--ngenes"), type="numeric", 
  			default=2000, 
              help="Number of genes to simulate", metavar="numeric"),
  make_option(c("-s", "--indSD"), type="numeric", 
  			default=.1, 
              help="Stddev on log10 scale between individuals", metavar="numeric"),
  make_option(c("-x", "--cellSD"), type="numeric", 
  			default=.02, 
              help="Number of cells to simulate", metavar="numeric"),
  make_option(c("-b", "--bimod"), type="numeric", 
  			default=5.0, 
              help="Bimodalisty Constant", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# David Read
# Wrapper script that runs large numbers of parameterss (like fold-changes) and
#    runs power analysis for each.

# Set up some defaults
# foldChangesToTry = c(1.0, 1.5, 2.0, 5.0, 7.5, 10.0) #c(2.0, 5.0) #c(1.1, 1.25, 1.5, 2, 3, 4, 5, 7.5, 10.0)
downSampleToTry = c(.1, .05,  .01, .02, .03, .04, .005, .001 )

cellsPerIndividual = opt$cells #paramList[["cells"]]   # 300
indivPerGroup = opt$indiv # paramList[["indiv"]]
ngenes = opt$ngenes #paramList[["ngenes"]]
zeiselDen=getZeiselDensities()
propHGE_toUse=0
noiseMean=0.0
noiseSDindiv =  opt$indSD # paramList[["indSD"]]
noiseSDcell =  opt$cellSD #paramList[["cellSD"]]
bimodToUse = opt$bimod

# Store the results in a DF to save
# resultsDF = data.frame(foldChange=foldChangesToTry, posResult=rep(NA,length(foldChangesToTry)),
# 		 negResult=rep(NA,length(foldChangesToTry)), noConvergence=rep(NA,length(foldChangesToTry))) 
# rownames(resultsDF) = as.character(foldChangesToTry)
randSeed=7
set.seed(randSeed)
indivNoiseDist=Norm(mean=noiseMean, sd=noiseSDindiv)
cellNoiseDist=Norm(mean=noiseMean, sd=noiseSDcell)
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace=TRUE) #replace is false in the SymSim vignette but I want to be able to make arbitrarily large


getPercentileVsMeanPlot <- function(observedCounts, eachSample, outDir, percentilesToPlot,
					UMI_median, processingNote){
	# Get the means
	observedMeans = rowSums(observedCounts) / ncol(observedCounts)
	observedMeans = sort(observedMeans)

	# Get percentile positions
	percentilePositions = percentilesToPlot * length(observedMeans)

	plotDF = data.frame("Percentile" = as.character(percentilesToPlot),
						"Gene_Mean" = observedMeans[percentilePositions])
	# browser()
	# Plot
	png(paste0(outDir, "percentileVsMean_Downsample=", as.character(eachSample), processingNote, ".png"),
				width=1000, height=1000, res=200)
	myPlot = ggplot(plotDF, aes_string(x="Percentile", y="Gene_Mean")) + 
				ggtitle(paste0("Percentile vs mean for ", as.character(eachSample), " UMI = ", as.character(UMI_median))) + 
				geom_point() + ylab("Observed mean counts")
	print(myPlot)
	dev.off()
}


getPercentileViolinPlot <-function(observedCounts, eachSample, outDir, percentilesToPlot,
					UMI_median, processingNote){
	# Sort observed Counts
	observedCounts = observedCounts[order(rowSums(observedCounts)),]

	# Get percentile positions
	percentilePositions = percentilesToPlot * nrow(observedCounts)

	# Get this data for violin plotting
	plotDF = data.frame(observedCounts[percentilePositions,])
	plotDF$Percentile = as.character(percentilesToPlot)
	# Get into plotting form
	plotDF = melt(plotDF, id.vars = "Percentile")

	# browser()
	# Plot the violin
	png(paste0(outDir, "percentile_violins_Downsample=", as.character(eachSample), processingNote, ".png"),
				width=1000, height=1000, res=200)
	myPlot = ggplot(plotDF, aes_string(x="Percentile", y="value")) + 
				ggtitle(paste0("Percentile violins for ", as.character(eachSample), " UMI = ", as.character(UMI_median))) + 
				geom_violin() + ylab("Observed counts")
	print(myPlot)
	dev.off()


}


processingNote = paste0("bimod=", as.character(bimodToUse), "indSD=", as.character(noiseSDindiv))

percentilesToPlot =  c(.2, .3, .4, .5, .6, .7, .8, .9, .99, .995, .999, .9995)

###############################################################################################
# Run simulations:

for (eachSample in downSampleToTry){

	print(paste0("Running downsample rate: ", as.character(eachSample)))
	foldChangeBetweenGroups = 1.0
	drawnParams <- getCellParamsFromPops(foldChangeBetweenGroups, zeiselDen, indivNoiseDist,
                            cellNoiseDist, indivPerGroup, cellsPerIndividual,
                             nGenes=ngenes, randomSeed=randSeed, bimod=bimodToUse)
	print("Getting True counts")
	trueCounts <- getTrueCountsGivenParams(drawnParams,
                ngenes=ngenes, prop_hge=propHGE_toUse, randseed=randSeed)

	print("Getting Observed Counts now")
	observedCounts <- GetObsCounts_Given_Param_In(true_counts=trueCounts[[1]], 
			protocol="UMI", alpha_mean=eachSample, alpha_sd=0.00, # alpha_sd=0.05,
			gene_len=gene_len, 
				depth_mean=max(2e4,(ngenes*160)), # Depth mean-> seq depth. Assume saturation
				depth_sd=1e3)

	# Get the cell-wide means. This can be used to estimate UMIs per cell, if I simulate enough genes
	UMI_median = median( colSums(observedCounts[[1]]) )
	# Scale so acts like ~20k genes in a transcriptome
	UMI_median = (UMI_median * (20000 / ngenes))

	getPercentileVsMeanPlot(observedCounts[[1]], eachSample, outDir, percentilesToPlot, UMI_median,
				processingNote)

	getPercentileViolinPlot(observedCounts[[1]], eachSample, outDir, percentilesToPlot, UMI_median,
				processingNote)

	# getMeanVsVarPlot(observedCounts[[1]], eachSample, outDir)
}

# # browser()
# # Add in the lists of pvalues/coefficients
# resultsDF$pValsMM = pvalMMList
# resultsDF$coefMM  = coefMMList
# resultsDF$pValsFM = pvalFMList
# resultsDF$coefFM  = coefFMList


# ###############################################################################################
# # Save output
# paramString = paste(as.character(foldChangesToTry), collapse="_")
# paramString = paste0(paramString,"_folds_", as.character(cellsPerIndividual), "_cells_",
# 					as.character(ngenes), '_genes_',
# 					as.character(noiseSDindiv), "sd_indiv_",
# 					as.character(noiseSDcell), "sd_cell",
# 					as.character(indivPerGroup), "_indiv.rds")
# # Save the Df to pull up later
# saveRDS(resultsDF, 
#   paste0( "/net/trapnell/vol1/home/readdf/trapLabDir/fibrosis/results/2020_12_29_SymSim_Power_Work/rdsOutputs_nobackup/",
#   			 paramString))



# print("All Done")
# # browser()

