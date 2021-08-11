#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)

print("Starting Now")
library("SymSim")
source("/net/trapnell/vol1/home/readdf/trapLabDir/fibrosis/results/2020_12_29_SymSim_Power_Work/simFuncsDavidMod.R")
library(distr)
source("/net/trapnell/vol1/home/readdf/trapLabDir/fibrosis/results/HeirarchalDE_Methods.R")
library("optparse")

# Get input params
paramList = processParamArgs()
print(paramList)


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
              help="Number of cells to simulate", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# David Read
# Wrapper script that runs large numbers of parameterss (like fold-changes) and
#    runs power analysis for each.

# Set up some defaults
foldChangesToTry = c(1.0, 1.5, 2.0, 5.0, 7.5, 10.0) #c(2.0, 5.0) #c(1.1, 1.25, 1.5, 2, 3, 4, 5, 7.5, 10.0)


cellsPerIndividual = opt$cells #paramList[["cells"]]   # 300
indivPerGroup = opt$indiv # paramList[["indiv"]]
ngenes = opt$ngenes #paramList[["ngenes"]]
zeiselDen=getZeiselDensities()
propHGE_toUse=0
noiseMean=0.0
noiseSDindiv =  opt$indSD # paramList[["indSD"]]
noiseSDcell =  opt$cellSD #paramList[["cellSD"]]

# Store the results in a DF to save
resultsDF = data.frame(foldChange=foldChangesToTry, posResult=rep(NA,length(foldChangesToTry)),
		 negResult=rep(NA,length(foldChangesToTry)), noConvergence=rep(NA,length(foldChangesToTry))) 
rownames(resultsDF) = as.character(foldChangesToTry)
randSeed=7
set.seed(randSeed)
indivNoiseDist=Norm(mean=noiseMean, sd=noiseSDindiv)
cellNoiseDist=Norm(mean=noiseMean, sd=noiseSDcell)
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace=TRUE) #replace is false in the SymSim vignette but I want to be able to make arbitrarily large
pvalMMList = vector(mode="list", length=length(foldChangesToTry))
coefMMList = vector(mode="list", length=length(foldChangesToTry))
pvalFMList = vector(mode="list", length=length(foldChangesToTry))
coefFMList = vector(mode="list", length=length(foldChangesToTry))





###############################################################################################
# Run simulations:
foldCount = 1
for (eachFold in foldChangesToTry){
	print(paste0("Running fold: ", as.character(eachFold)))
	drawnParams <- getCellParamsFromPops(eachFold, zeiselDen, indivNoiseDist,
                            cellNoiseDist, indivPerGroup, cellsPerIndividual,
                             nGenes=ngenes, randomSeed=randSeed)
	trueCounts <- getTrueCountsGivenParams(drawnParams,
                ngenes=ngenes, prop_hge=propHGE_toUse, randseed=randSeed)
	print("Getting Observed Counts now")
	observedCounts <- GetObsCounts_Given_Param_In(true_counts=trueCounts[[1]], 
			protocol="UMI", alpha_mean=0.1, alpha_sd=0.05, gene_len=gene_len, 
				depth_mean=max(2e4,(ngenes*160)), # Depth mean-> seq depth. Assume saturation
				depth_sd=1e3)

	# Run fixed effects:
	deResultsFM = runDEtest(observedCounts, ngenes, cellsPerIndividual, indivPerGroup,
						testType="fixedModel", errorForm="NegBinom")
	# Get the outputs
	resultsDF[as.character(eachFold),"posResultFM"] = deResultsFM["posRate"]
	resultsDF[as.character(eachFold), "negResultFM"] = deResultsFM["negRate"]
	resultsDF[as.character(eachFold), "noConvergenceFM"] = deResultsFM["noConvRate"]
	pvalFMList[[foldCount]] = deResultsFM[["pvals"]]
	coefFMList[[foldCount]] = deResultsFM[["coefs"]]

	# Mixed model
	deResultsMM = runDEtest(observedCounts, ngenes, cellsPerIndividual, indivPerGroup,
						testType="mixedModel", errorForm="NegBinom")	
	# Get the outputs
	resultsDF[as.character(eachFold),"posResultMM"] = deResultsMM["posRate"]
	resultsDF[as.character(eachFold), "negResultMM"] = deResultsMM["negRate"]
	resultsDF[as.character(eachFold), "noConvergenceMM"] = deResultsMM["noConvRate"]
	pvalMMList[[foldCount]] = deResultsMM[["pvals"]]
	coefMMList[[foldCount]] = deResultsMM[["coefs"]]

	foldCount = foldCount + 1
}

# browser()
# Add in the lists of pvalues/coefficients
resultsDF$pValsMM = pvalMMList
resultsDF$coefMM  = coefMMList
resultsDF$pValsFM = pvalFMList
resultsDF$coefFM  = coefFMList


###############################################################################################
# Save output
paramString = paste(as.character(foldChangesToTry), collapse="_")
paramString = paste0(paramString,"_folds_", as.character(cellsPerIndividual), "_cells_",
					as.character(ngenes), '_genes_',
					as.character(noiseSDindiv), "sd_indiv_",
					as.character(noiseSDcell), "sd_cell",
					as.character(indivPerGroup), "_indiv.rds")
# Save the Df to pull up later
saveRDS(resultsDF, 
  paste0( "/net/trapnell/vol1/home/readdf/trapLabDir/fibrosis/results/2020_12_29_SymSim_Power_Work/rdsOutputs_nobackup/",
  			 paramString))



print("All Done")
# browser()

