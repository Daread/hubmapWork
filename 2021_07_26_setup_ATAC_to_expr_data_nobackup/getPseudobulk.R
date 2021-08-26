
library(monocle3)
library("optparse")
library(ggplot2)

outDir = "rdsOutputs/RNA_Data/"
dir.create(outDir)

# Get the passed parameters
option_list = list(
   make_option(c("-c", "--cdsRNA"), type="character", 
  			default="allCells_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes.rds", 
              help="Name of RNA cds to use", metavar="character"),
  make_option(c("-p", "--rnaPath"), type="character", 
  			default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/", 
              help="Path to RNA cds", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

processingNote = "countsPerMillionNotDeduplicated"


##
# load RNA cds
cellCDS = readRDS(paste0(opt$rnaPath, opt$cdsRNA))

# Set up a dataframe that'll hold the pseudobulked 
pseudobulkedDF = as.data.frame(rowData(cellCDS)[c("gene_short_name", "id", "chromosome", "bp1", "bp2", "gene_strand")])

pseudoCount = .125

# Loop and get the data
for (eachCelltype in levels(as.factor(colData(cellCDS)$highLevelCellType))){
	# Get a subset
	subsetCDS = cellCDS[,colData(cellCDS)$highLevelCellType == eachCelltype]
	# Add the pseudobulked data.
	countsPerMillion = (rowSums(exprs(subsetCDS)) / sum(exprs(subsetCDS))) * 1e6
	pseudobulkedDF[paste0(eachCelltype, "_CPM")] = countsPerMillion
  pseudobulkedDF[paste0(eachCelltype, "_log2_CPM")] = log2(countsPerMillion + pseudoCount)
}


# See a histogram of what these counts look like
for (eachCelltype in levels(as.factor(colData(cellCDS)$highLevelCellType))){
  png(paste0("./plots/", eachCelltype, "pseudoBulkHistoramLinear.png"), 
          width=1000, height=1000, res=200)
  myPlot = ggplot(pseudobulkedDF, aes_string(x=paste0(eachCelltype, "_CPM"))) + 
          geom_histogram()
  print(myPlot)        
  dev.off()

  # Log
  png(paste0("./plots/", eachCelltype, "pseudoBulkHistoramLog2.png"), 
          width=1000, height=1000, res=200)
  myPlot = ggplot(pseudobulkedDF, aes_string(x=paste0(eachCelltype, "_log2_CPM"))) + 
          geom_histogram()
  print(myPlot)        
  dev.off()

}



# # Save the output
# saveRDS(pseudobulkedDF, paste0(outDir, "Pseudobulked_RNA_", processingNote, ".RDS"))


cellCPMmatrix = as.matrix(pseudobulkedDF[,paste0(levels(as.factor(colData(cellCDS)$highLevelCellType)), "_CPM")])
cellCPMmatrix = cellCPMmatrix + pseudoCount

geneMeans = rowSums(cellCPMmatrix) / length(levels(as.factor(colData(cellCDS)$highLevelCellType)))

# Additionally, calculate a version where the values are in terms of ratio vs. the mean, taken across all cell types
for (eachCelltype in levels(as.factor(colData(cellCDS)$highLevelCellType))) {
  # Get pseudobulked data
  pseudobulkedDF[paste0(eachCelltype, "_log2_ratio_Vs_AllTypeMean")] = (
    log2(cellCPMmatrix[,paste0(eachCelltype, "_CPM")] / geneMeans)
    )
}



# Save the output
saveRDS(pseudobulkedDF, paste0(outDir, "Pseudobulked_RNA_", processingNote, ".RDS"))
write.csv(pseudobulkedDF, paste0(outDir, "Pseudobulked_RNA_", processingNote, ".csv"))

print("All Done")

# rowData(cellCDS)$geneSums = rowSums(exprs(cellCDS))
# rowData(cellCDS[duplicated(rowData(cellCDS)$gene_short_name, fromLast=TRUE),])
# rowData(cellCDS[rowData(cellCDS)$gene_short_name == "MATR3",])












