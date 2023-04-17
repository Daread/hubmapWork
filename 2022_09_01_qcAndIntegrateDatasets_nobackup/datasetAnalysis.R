
# library("devtools")
# library('spdep')
# load_all("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/misc/monocle3mmVersion/monocle3")

library("monocle3")
library(ggplot2)

# Load helper functions
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

formatCelltypes <- function(fullCDS){
  # Get rid of cells that weren't annotated across datasets
  typesToRemove = c("Epicardium", "Native_Cell", "Doublets", "Unassigned", "Atrial_Cardiomyocytes", "Mast_Cells")
  fullCDS = fullCDS[,!(colData(fullCDS)$Cell_Shared_Label %in% typesToRemove )]

  # Format names
  colData(fullCDS)$Cell_Shared_Label = gsub("_", " ", colData(fullCDS)$Cell_Shared_Label) 

  return(fullCDS)
}


set.seed(7)
library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-c", "--cdsRNA"), type="character", 
        # default="../2022_08_22_addNewSamples_nobackup/formattedData/allDatasetCDS.rds", 
        default="../2022_08_22_addNewSamples_nobackup/formattedData/NucleiOnlySharedGenesCDS.rds", 
              help="Name of RNA cds to use", metavar="character"),
  make_option(c("-d", "--downSamp"), type="numeric", 
        default=5,   
        help="Down sample rate of cells to speed processing", metavar="numeric"),
  make_option(c("-v", "--varToMNN"), type="character", 
        default="DataSource", 
              help="Name of RNA cds to use", metavar="character"),
  make_option(c("-n", "--nuclOnly"), type="logical", 
        default=TRUE, 
              help="Name of RNA cds to use", metavar="character"),
  make_option(c("-a", "--removeAtrium"), type = "logical",
        default=TRUE,
        help="Remove LA and RA samples?")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

processingNote = "SharedGenes"

if (!is.na(opt$varToMNN)){
  mnnVec = strsplit(opt$varToMNN, ",", fixed=TRUE)[[1]]
}

loadCDS = FALSE


# Read in data
cdsPath = opt$cdsRNA
fullCDS = readRDS(cdsPath)

# Shuffle order for plotting, and select the downsampled data
fullCDS = fullCDS[,sample(ncol(fullCDS))]
fullCDS = fullCDS[,1:as.integer(ncol(fullCDS) / opt$downSamp)]


if (opt$nuclOnly){
	processingNote = paste0("Nuclei_Only_", processingNote)
	fullCDS = fullCDS[,colData(fullCDS)$DataType == "Nuclei"]
}

if (opt$removeAtrium){
  processingNote = paste0("No_Atrium_", processingNote)
  fullCDS = fullCDS[,!(colData(fullCDS)$Anatomical_Site %in% c("LA", "RA")) ]
  fullCDS = fullCDS[,!(colData(fullCDS)$Cell_Shared_Label == "Atrial_Cardiomyocytes")]
}

# Format to remove unwanted cell types
fullCDS = formatCelltypes(fullCDS)



# Try the minimalist strategy of UMAPing all together, batch correct by source
fullCDS = estimate_size_factors(fullCDS)
fullCDS = preprocess_cds(fullCDS, num_dim=30)


if (is.na(opt$varToMNN)){
  processingNote = paste0("No_MNN_Downsamp_", as.character(opt$downSamp), processingNote)
} else if (length(mnnVec) == 1){
	
  processingNote = paste0("MNN_By_", paste(mnnVec, collapse = "_"), "_Downsamp_", as.character(opt$downSamp), processingNote)
  fullCDS = align_cds(fullCDS, alignment_group= mnnVec[1])
} else{
  processingNote = paste0("MNN_By_", paste(mnnVec, collapse = "_"), "_Downsamp_", as.character(opt$downSamp), processingNote)
	fullCDS = align_cds(fullCDS, alignment_group= mnnVec[1],
		residual_model_formula_str = paste0("~ ", paste(mnnVec[2:length(mnnVec)], collapse = " + ")))
  }


fullCDS = reduce_dimension(fullCDS, preprocess_method="Aligned")




dir.create("./rdsOutputs/")

# Save this cds if we want to come back later & not re-run the embedding steps
saveRDS(fullCDS, file=paste0("./rdsOutputs/", processingNote, "_CDS.rds"))

# Load RDS?
if (loadCDS){
  fullCDS = readRDS(paste0("./rdsOutputs/", processingNote, "_CDS.rds"))
}





# Plot outputs
dir.create('./plots/')

colData(fullCDS)$Publication = colData(fullCDS)$DataSource
colData(fullCDS)["Cell Type"] = colData(fullCDS)$Cell_Shared_Label

# For figure panels

plotUMAP_Monocle(fullCDS, processingNote, "Publication", 
  outputPath= "./plots/", show_labels=FALSE, textSize = 22)


plotUMAP_Monocle(fullCDS, processingNote, "Cell Type", 
  outputPath= "./plots/", show_labels=FALSE, textSize = 22)







plotUMAP_Monocle(fullCDS, processingNote, "DataSource", 
	outputPath= "./plots/", show_labels=FALSE)


plotUMAP_Monocle(fullCDS, processingNote, "Cell_Shared_Label", 
	outputPath= "./plots/", show_labels=FALSE)


plotUMAP_Monocle(fullCDS, processingNote, "Anatomical_Site", 
	outputPath= "./plots/", show_labels=FALSE)


plotUMAP_Monocle(fullCDS, processingNote, "log10_umi", 
	outputPath= "./plots/", show_labels=FALSE)

plotUMAP_Monocle(fullCDS, processingNote, "DataType", 
	outputPath= "./plots/", show_labels=FALSE)





plotUMAP_Monocle_faceted <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, textSize=10, outputPath = "./plots/",
                    returnPlotObj=FALSE){ #, xVal, yVal){
    png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored.png"),
             width=2000, height=2000, res=200)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels,
          cell_stroke=.1 , group_label_size=textSize        
                ))
    myPlot = (myPlot + theme(text=element_text(size=textSize)) + facet_wrap(vars(get("Cell_Shared_Label")), scales="free"))
    print(myPlot)
    dev.off()   

    if (returnPlotObj){
        return(myPlot)
    }

}


plotUMAP_Monocle_faceted(fullCDS, paste0("Faceted_", processingNote), "DataSource", 
	outputPath= "./plots/", show_labels=FALSE)



plotUMAP_Monocle_faceted(fullCDS, paste0("Faceted_", processingNote), "Cell_Shared_Label", 
  outputPath= "./plots/", show_labels=FALSE)




