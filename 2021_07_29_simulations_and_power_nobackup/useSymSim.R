
# Use SymSim to generate simlated datasets for use in DE testing and power analysis
print("Starting Now")

library("SymSim")

source("./simFuncsDavidMod.R")
library(distr)


# Testing out, per vignette
ngenes <- 100
true_counts_res <- SimulateTrueCounts(ncells_total=300, ngenes=ngenes,
				 evf_type="one.population", Sigma=0.4, randseed=0)

# Get observed counts
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], 
	meta_cell=true_counts_res[[3]], 
	protocol="nonUMI", alpha_mean=0.1, alpha_sd=0.05, gene_len=gene_len, depth_mean=1e5, depth_sd=3e3)



getTree <- function(treePath){
	thisPhyla = read.tree(treePath)
	return(thisPhyla)
}





source("./simFuncsDavidMod.R")

# Testing out my sim functions
zeiselDen = getZeiselDensities()
foldToDraw = 10.0
kineticDraw = drawKineticParams(foldToDraw, zeiselDen)

meanToUse = 0
sdToUse = .2
indivNoiseDist = Norm(mean=meanToUse, sd=sdToUse)

indivToDraw = 10
indivMeanSvals = drawFromDist(indivToDraw, indivNoiseDist)

######################################################################



# Test 1-5-20:
ngenes <- 50
# true_counts_res <- getTrueCountsGivenParams(ncells_total=300, ngenes=ngenes,
# 				 evf_type="one.population", Sigma=0.4, randseed=0)
zeiselDen = getZeiselDensities()
foldToDraw = 10.0

meanToUse = 0.0
sdToUse = .2
indivNoiseDist = Norm(mean=meanToUse, sd=sdToUse)
cellNoiseDist = Norm(mean=meanToUse, sd=(sdToUse/2.0))

indivPerGroup = 5
cellsPerSample = 100

randSeed=7
drawnParams <-  getCellParamsFromPops(foldToDraw, zeiselDen, indivNoiseDist,
                            cellNoiseDist, indivPerGroup, cellsPerSample,
                             nGenes=ngenes, randomSeed=randSeed)


source("./simFuncsDavidMod.R")
# Draw actual counts
randSeed=7
trueCountRes <- getTrueCountsGivenParams(#ncells_total,min_popsize,
                drawnParams,
                ngenes=ngenes, 
                prop_hge=0,
                               randseed=randSeed)

source("./simFuncsDavidMod.R")
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
observed_counts <- GetObsCounts_Given_Param_In(true_counts=trueCountRes[[1]], 
	protocol="UMI", alpha_mean=0.1, 
	alpha_sd=0.05, gene_len=gene_len, 
	depth_mean=2e4, depth_sd=1e3)




# Make a dataframe (analagous to coldata in a cds) holding individual info
groupLabels = rep(c("Control", "Case"), each=(indivPerGroup*cellsPerSample))
indivLabels = rep((1:(indivPerGroup*2)), each= cellsPerSample)

cellMetaDF = data.frame(Group=groupLabels,Individual=indivLabels)
cellMetaDF$log10UMI = colSums(trueCountRes[[1]])


# getCellParamsFromPops <- function(meanFoldChange, kineticDensities, indivNoiseDist,
#                             cellNoiseDist, indivsPerPop, cellsPerIndiv,
#                             popOneSetExpressionLevel = NULL, nGenes=ngenes)

# Try fitting a model
source("../HeirarchalDE_Methods.R")

genesToFit = 50
myFitVec = vector(mode="list", length=genesToFit)
pvalVec = rep(NA,genesToFit)
for (eachGeneInd in 1:genesToFit){
	print(eachGeneInd)
	thisDF = cellMetaDF
	thisDF["GeneCount"] = as.vector(trueCountRes[[1]][eachGeneInd,])
	thisFit= RunMixedModel(thisDF, "GeneCount",
		         c("log10UMI", "Group"), "Individual")
	myFitVec[eachGeneInd] = coef(summary(thisFit))

	thisPval = coef(summary(thisFit))[[1]][,4][["GroupControl"]]
	thisCoef = coefficients(thisFit)[["Group"]]

	pvalVec[eachGeneInd] = thisPval
}






# Testing the vignette-defined workflow
#############################################
myTreePath = "./5v5_two_group_tree.txt"
myTree = getTree(myTreePath)


# testTree = getTree("./testTree.txt")
cellsToUse = 1000
ngenes = 300
mySim_true_count_res = SimulateTrueCounts(ncells_total=cellsToUse, min_popsize=50,
									 i_minpop=2, ngenes=ngenes, nevf=10, 
									evf_type="discrete", n_de_evf=9, vary="s",
								 Sigma=0.4, phyla=myTree, randseed=0)

data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
mySim_obs_counts = True2ObservedCounts(true_counts=mySim_true_count_res[[1]],
					meta_cell=mySim_true_count_res[[3]], protocol="UMI", alpha_mean=0.05,
					 alpha_sd=0.02, gene_len=gene_len, depth_mean=1e4, depth_sd=2e3)

# Make a plot
tsne_UMI_counts <- PlotTsne(meta=mySim_obs_counts[[2]],
			 data=log2(mySim_obs_counts[[1]]+1),
			  evf_type="discrete", n_pc=20, label='pop', 
			  saving = TRUE, plotname="./obsCountstSNE.pdf")


# Group individuals by condition
grouopOneLabs = c(1,2,3,4,5)
groupTwoLabs = c(6,7,8,9,10)
mySim_true_count_res$cell_meta$origPop = mySim_true_count_res$cell_meta$pop
mySim_true_count_res$cell_meta$pop = ifelse(
			(mySim_true_count_res$cell_meta$origPop %in% grouopOneLabs), "cont", "disease" )

# Get some DE genes.
deRes <- getDEgenes(true_counts_res = mySim_true_count_res, popA="cont", popB="disease")



print("All Done")