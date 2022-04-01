
###########################################################################################
# Load relevant packages
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()
print(modStatus)
library("optparse")


# library("lme4")
# Source the draft code for mixed model fitting
# source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/DE_Code/mixedModelFittingMonocle3Style.R")

filterByCellsExpressed = function(inputCDS, cellPropMin){
	rowData(inputCDS)$fracCellsExpr = (rowSums(exprs(inputCDS) != 0) /
								nrow(colData(inputCDS)) )
	inputCDS = inputCDS[rowData(inputCDS)$fracCellsExpr > cellPropMin,]
	return(inputCDS)
}


# From https://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
# RowVar <- function(x, ...) {
#   rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
# }


# Get the passed parameters
option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
  			default="HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes", 
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-c", "--cellType"), type="character", default="Cardiomyocyte", 
              help="Cell Type to subset for testing", metavar="character"),
    make_option(c("-g", "--groupColumn"), type="character", default="highLevelCellType", 
              help="column in colData that cellType uses to select", metavar="character"),
    make_option(c("-m", "--minGeneExprPercent"), type="numeric", default=.01, 
              help="Proportion of cells that must express a gene for it to be tested", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
# Parse random and fixed effects
# fixedEffectVec = strsplit(opt$fixedEffects, ",")[[1]]
# randomEffectVec = strsplit(opt$randomEffects, ",")[[1]]

processingNote = paste0(opt$processingNote, "_", opt$cellType)

# Get a CDS where I've already annotated the high level cell types
rdsPath = "../2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/"
# oldProcNote = "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes"
oldProcNote = opt$processingNote
allCellCDS = readRDS(paste0(rdsPath, "allCells_", oldProcNote, ".rds"))
print("Read in CDS")

# Subset down to a cell type
myCDS = allCellCDS[,colData(allCellCDS)$highLevelCellType == opt$cellType]

# Make an output directory for work on this cell type
outputDir = paste0("./plots/", opt$cellType, "/")
dir.create(outputDir)




plotMeanVar <- function(inputCDS, processingNote, outputDir, outliersToDrop =4){

  # First, just plot the overall mean/variance relationship
  meanVarDF = data.frame("gene_means" = rowSums(exprs(inputCDS))/ncol(exprs(inputCDS)))
  meanVarDF$gene_var = apply(exprs(inputCDS), 1, var)

  # Skip a few outliers
  meanVarDF = meanVarDF[order(meanVarDF$gene_means),]
  outliersToDrop = outliersToDrop
  meanVarDF = meanVarDF[outliersToDrop:(nrow(meanVarDF)-outliersToDrop),]

  png(paste0(outputDir, processingNote, "_meanVsVarAllGenes.png"),
        width=800, height=800, res=200)
  myPlot = (ggplot(meanVarDF, aes(gene_means, gene_var)) + geom_point() 
    # + geom_spline()
    # + geom_line(aes(meanVarDF$gene_means, meanVarDF$gene_var) data=data.frame(spline(gene_means, gene_var))) 
    )
  print(myPlot)
  dev.off()

  return(meanVarDF)
}

# Plot for all genes
fullMeanVarRes = plotMeanVar(myCDS, processingNote, outputDir)
# Now for one sample at a atime
for (eachSample in as.character(levels(as.factor(myCDS$sampleName )))){
  # Subset and plot
  thisCDS = myCDS[,colData(myCDS)$sampleName == eachSample]
  meanVarRes = plotMeanVar(thisCDS, paste0(processingNote, "_", eachSample), outputDir)
}



myCDS = estimate_size_factors(myCDS)

# Take a look at some genes' distributions. Pick some examples to plot based on percentiles
percentilesToPlot = c(.2, .3, .4, .5, .6, .7, .8, .9, .99, .995, .999, .9995)

# percentilesToPlot = c(.25, .5, .75)
# Now, take a look at individual genes' variation. 
rowData(myCDS)$meanExpr = rowSums(exprs(myCDS)) / ncol(exprs(myCDS))
myCDS = myCDS[order(rowData(myCDS)$meanExpr),]

# Get the few genes at percentiles
percentilePos = as.integer(percentilesToPlot * nrow(exprs(myCDS)))
percentileCDS = myCDS[percentilePos,]

# Plot gene counts. As for mean/variance plots, do this for all cells and then split up by sample
plotIndivGenes <- function(inputCDS, processingNote, outputDir){
  png(paste0(outputDir, processingNote, "_PercentileGeneViolinsNormAndLog.png"),
    height=1200, width=1200, res = 200)
  myPlot = plot_genes_violin(inputCDS, group_cells_by="Donor",
                        min_expr=.001, nrow=4, ncol=3,
                         panel_order = as.character(rowData(inputCDS)$gene_short_name))
  print(myPlot)
  dev.off()

  # Now one that isn't normalized
  png(paste0(outputDir, processingNote, "_PercentileGeneViolinsLogNoNorm.png"),
    height=1200, width=1200, res = 200)
  myPlot = plot_genes_violin(inputCDS, group_cells_by="Donor",
                        min_expr=.001, nrow=4, ncol=3, normalize=FALSE,
                         panel_order = as.character(rowData(inputCDS)$gene_short_name))
  print(myPlot)
  dev.off()

  # Non-normalized plus linear scale:
  png(paste0(outputDir, processingNote, "_PercentileGeneViolinsNoLogNoNorm.png"),
    height=1200, width=1200, res = 200)
  myPlot = plot_genes_violin(inputCDS, group_cells_by="Donor",
                        min_expr=.001, nrow=4, ncol=3, normalize=FALSE, log_scale=FALSE,
                         panel_order = as.character(rowData(inputCDS)$gene_short_name))
  print(myPlot)
  dev.off()

  return(myPlot)
}

allCellsViolin = plotIndivGenes(percentileCDS, processingNote, outputDir)








# png(paste0(outputDir, processingNote, "_PercentileGeneViolinsNoGroup.png"))
# myPlot = monocle3::plot_genes_violin(percentileCDS, min_expr=.01)
# print(myPlot)
# dev.off()






# Take a look at the mean count for various percentiles
percentileMeans = data.frame("geneMeans" = rowData(percentileCDS)$meanExpr,
                            "percentile" = as.character(percentilesToPlot))

# Plot
png(paste0(outputDir, processingNote, "_Percentile_vs_Mean_Expr.png"), 
            height=1200, width=1200, res = 200)
myPlot = ggplot(percentileMeans, aes_string(x="percentile", y="geneMeans")) +  
        geom_point() + ggtitle("Percentile vs. mean expression")
print(myPlot)
dev.off()



















