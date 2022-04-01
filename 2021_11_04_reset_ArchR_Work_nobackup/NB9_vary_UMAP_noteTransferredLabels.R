
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
dir.create(paste0(basepath, "archr/results/NB9/"))
out_dir = paste0(basepath, "archr/results/NB9/")
setwd(out_dir)
set.seed(7)

# load requirements
suppressPackageStartupMessages({
  # library(ArchR)
  library(monocle3)
  library(dplyr)
  library(ggrastr)
  library(Seurat)
})
library(tidyr)
library(plyr)
library(ggplot2)


source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

library("optparse")

# Get the passed parameters
option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
  			default="HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes", 
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-s", "--sampleRNAname"), type="character",  default= "W144.Apex", #  "All_Cells",  # default="W142.Left.Vent", 
              help="Sample to integrate", metavar="character"),
    make_option(c("-u", "--useSubset"), type="character", default="Subset", 
              help="Subset or All", metavar="character"),

    make_option(c("-m", "--useMNN"), type="character", default="useMNN", 
              help="Subset or All", metavar="character"),

    make_option(c("-c", "--fripMin"), type="numeric", default=0.2, 
              help="Min FRIP value to permit", metavar="numeric"),

    # make_option(c("-o", "--lowDimOne"), type="numeric", default=1, 
    #           help="Coordinate to start in LSI", metavar="numeric"),

    # make_option(c("-t", "--lowDimTwo"), type="numeric", default=50, 
    #           help="Highest coord in LSI", metavar="numeric"),

    make_option(c("-z", "--pcToUse"), type="numeric", default=10, 
              help="Principel components/LSI coords to use", metavar="numeric"),

    make_option(c("-f", "--featureSet"), type="character", default="peakMat", # "peakMat" for greg/riza matrix, "bMat" for archr bins
                                                                              # "gMat" for archr activity scores
              help="peakMat, bMat, or gMat", metavar="character"),

    make_option(c("-a", "--ATACprocNote"), type="character", default="FRIP=0.2_FRIT=0.05UMI=1000DL=0.5", 
              help="How ATAC Cells were filtered", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

samplesATACnames = c(
  "W134.heart.apex.s1",
  "W135.heart.LV.s1", 
  "W136.heart.apex.s1", "W136.heart.LV.s1", 
  "W139.heart.apex.s1", "W139.heart.LV.s1", "W139.heart.RV.s1", "W139.heart.septum.s1", 
  "W142.heart.LV.s1", 
  "W144.heart.apex.s1", 
  "W145.heart.apex.s1", "W145.heart.LV.s1", 
  "W146.heart.apex.s1", "W146.heart.LV.s1")
names(samplesATACnames) = c("W134.Apex",
  "W135.Left.Vent", 
  "W136.Apex", "W136.Left.Vent", 
  "W139.Apex", "W139.Left.Vent", "W139.Right.Vent", "W139.Septum", 
  "W142.Left.Vent", 
  "W144.Apex", 
  "W145.Apex", "W145.Left.Vent", 
  "W146.Apex", "W146.Left.Vent")

PPIval = 300

# Get the ATAC cds
atacProcNote =  opt$ATACprocNote # "FRIP=0.2_FRIT=0.05UMI=1000"


# }
plotUMAP_MonocleModded <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, textSize=10, outputPath = "./plots/"){ #, xVal, yVal){
    png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored.png"),
             width=1400, height=1000, res=200)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels,
          cell_stroke=.1 , group_label_size=textSize        
                )) 
    myPlot = (myPlot + theme(text=element_text(size=textSize)) +facet_wrap(~predicted.id))
    print(myPlot)
    dev.off()   

}


# Get the saved CDS
rdsPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/rdsOutput/"
# fileName = paste0(opt$sampleRNAname, "_", opt$ATACprocNote, "gMatrixCDS_postTransfer.RDS")
fileName = paste0("All_Cells", "_", opt$ATACprocNote, "gMatrixCDS_postTransfer.RDS")

# Read in the data
cds.a.g = readRDS(paste0(rdsPath, fileName))
# bMatFile = paste0(opt$sampleRNAname, "_", opt$ATACprocNote, "bMatrixCDS_postTransfer.RDS")
bMatFile = paste0("All_Cells", "_", opt$ATACprocNote, "bMatrixCDS_postTransfer.RDS")
cds.a.b = readRDS(paste0(rdsPath, bMatFile))

# Get data from Greg's filtering process
# Added 7-29-21: Filter by checking names in a cds filtered by greg/riza's method
filterCDSName = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/cds_objects/"
processingNote = opt$ATACprocNote

oldProcNote = "FRIP=0.2_FRIT=0.05UMI=1000"
# cds_p holds data
load(paste0(filterCDSName, "cds_p_allHeartATAC", oldProcNote))


# Added 10-26-21: Get the CDS with gene activity scores calculated by cicero links
ciceroActivPath = paste0("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_24_ciceroWork_nobackup/fileOutputs/", 
                       "CiceroActivityCDS_", "Cicero_Defaults_cds_p_allHeartATAC", opt$ATACprocNote, "_ciceroConnections.RDS")
ciceroActivCDS = readRDS(ciceroActivPath)

##################################################################################### cicero activity cds formatting
# Reformat sampleName
colDataDF = separate(data=as.data.frame(colData(ciceroActivCDS)), col = "sampleName", 
                        into = c("sampleName", "Processing"), sep="FRIP", remove=TRUE)
# colData(ciceroActivCDS)$sampleName = colDataDF$sampleName
ciceroActivCDS$sampleName = colDataDF$sampleName

# Reformat the colnames of the activity matrix to match the bin/arrowActivity cds's
colDataDF$rowname = rownames(colDataDF)
colDataDF = separate(data=colDataDF, col = "rowname", 
                        into = c("Part1", "Part2", "Part3"), sep="_", remove=TRUE)
rownames(colDataDF) = paste0(colDataDF$Part1, "_", colDataDF$Part2, "_", colDataDF$Part3)


rownames(colData(ciceroActivCDS)) = paste0(colData(ciceroActivCDS)$sampleName, "#", rownames((colDataDF)))

colnames(ciceroActivCDS) = rownames(colData(ciceroActivCDS))

#####################################################################################

################################################################################# cds_p formatting
# Get the sampleName reformated
colDataDF = separate(data=as.data.frame(colData(cds_p)), col = "sampleName", 
                        into = c("sampleName", "Processing"), sep="FRIP", remove=TRUE)
cds_p$sampleName = colDataDF$sampleName
# Reformat the colnames of cds_p (original/peak matrix) to match the others
colDataDF$rowname = rownames(colDataDF)
colDataDF = separate(data=colDataDF, col = "rowname", 
                        into = c("Part1", "Part2", "Part3"), sep="_", remove=TRUE)
rownames(colDataDF) = paste0(colDataDF$Part1, "_", colDataDF$Part2, "_", colDataDF$Part3)


rownames(colData(cds_p)) = paste0(colData(cds_p)$sampleName, "#", rownames((colDataDF)))
colnames(cds_p) = rownames(colData(cds_p))
#####################################################################################





# Now get the portion of cds_p we want, in order
ciceroActivCDS = ciceroActivCDS[,colnames(cds.a.g)]

cds_p = cds_p[,colnames(cds.a.g)]
colData(cds.a.g)$sampleName = colData(cds_p)$sampleName
colData(cds.a.b)$sampleName = colData(cds_p)$sampleName
# Get data from 

cdsList = list(cds_p, cds.a.g, cds.a.b, ciceroActivCDS)
names(cdsList) = c("Peak_CDS", "Activity_CDS", "Bin_CDS", "CiceroActiv_CDS")

preprocessMethods = c("LSI", "PCA", "LSI", "PCA")
names(preprocessMethods) = c("Peak_CDS", "Activity_CDS", "Bin_CDS", "CiceroActiv_CDS")


# opt$sampleRNAname = "W144.heart.apex.s1"
# opt$useMNN = "NoMNN"

out_dir = paste0(basepath, "archr/results/NB9/", opt$sampleRNAname, "/")
dir.create(out_dir)
setwd(out_dir)

# Loop, make a UMAP for each of these and 
for (eachCDSname in names(cdsList)){
  thisCDS = cdsList[[eachCDSname]]
  colData(thisCDS)$Sample = colData(thisCDS)$sampleName
  colData(thisCDS)$FRIP = colData(cds_p)$FRIP 
  colData(thisCDS)$umi = colData(cds_p)$umi
  colData(thisCDS)$log10umi = log10(colData(cds_p)$umi)
  colData(thisCDS)$doublet_likelihood = colData(cds_p)$doublet_likelihood

  # Subset?
  if (!(opt$sampleRNAname == "All_Cells")){
    # thisCDS = thisCDS[,colData(thisCDS)$sampleName == opt$sampleRNAname]
    # miniCDS_g = cds.a.g[,colData(cds.a.g)$sampleName == opt$sampleRNAname]
    thisCDS = thisCDS[,colData(thisCDS)$sampleName == samplesATACnames[opt$sampleRNAname]]
    miniCDS_g = cds.a.g[,colData(cds.a.g)$sampleName == samplesATACnames[opt$sampleRNAname]]
  } else {
    miniCDS_g = cds.a.g
  }

  # Subset by FRIP
  miniCDS_g$FRIP = thisCDS$FRIP
  miniCDS_g = miniCDS_g[,colData(miniCDS_g)$FRIP > opt$fripMin]
  thisCDS = thisCDS[,colData(thisCDS)$FRIP > opt$fripMin]

  print(paste0("Working on ", eachCDSname))
  print(str(colData(thisCDS)))

  thisCDS = estimate_size_factors(thisCDS)
  thisCDS = preprocess_cds(thisCDS, method=preprocessMethods[eachCDSname], num_dim=opt$pcToUse)

  # MNN?
  if (opt$useMNN == "useMNN"){
    thisCDS = align_cds(thisCDS, preprocess_method = preprocessMethods[eachCDSname], alignment_group = "Sample", 
                                         residual_model_formula_str = ~log10umi + FRIP + doublet_likelihood)
    thisCDS = reduce_dimension(thisCDS, reduction_method = 'UMAP', preprocess_method = "Aligned")
  } else{
    thisCDS = reduce_dimension(thisCDS, reduction_method = 'UMAP', preprocess_method = preprocessMethods[eachCDSname])
  }
  # Reduce dimension

  # Give these over to the g matrix
  thisProcNote = paste0("UMAP_by_", eachCDSname, "_", opt$useMNN, "FRIP_", as.character(opt$fripMin), "PC", as.character(opt$pcToUse))

  reducedDims(miniCDS_g)$UMAP = reducedDims(thisCDS)$UMAP
  plotUMAP_MonocleModded(miniCDS_g, paste0(thisProcNote, "_faceted"), "sampleName",
                        outputPath = out_dir, show_labels=FALSE)

  plotUMAP_Monocle(miniCDS_g, thisProcNote, "sampleName",
                        outputPath = out_dir, show_labels=FALSE)

  # Now color by some marker genes
  plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("TTN", "RBPJ","GSN", "LDB2", "MYH11", "THEMIS" ),
                        "generalMarkers",
                        outputPath = out_dir)
  # 
  plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("PECAM1", "CD106", "CD62E", "SELE", "KDR", "ENG" ),
                        "coleEndothMark",
                        outputPath = out_dir)
  plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("MYH11", "SMTN", "CALD1", "CNN1", "CNN2" ),
                        "coleSmoothMuscMark",
                        outputPath = out_dir)
  plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("CD68", "CD14", "FCGR1", "MS4A7", "FCER1A", "CD163", "LY6C1", "FCN1", "MERTK" ),
                        "coleMonoMacMark",
                        outputPath = out_dir)
  plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("CD3E", "CD3G", 'CD3D", "CD4", "CD8A', "CD8B", "BCL11B" ),
                        "coleTCellMark",
                        outputPath = out_dir)
  plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("CD90", "THY1", "VIM", "DES" ),
                        "coleFibroblastMark",
                        outputPath = out_dir)
}




# rownames(colDataDF) = colDataDF$rowname

# rownames(colData(cds.a.g)) = paste0(rownames(colData(cds.a.g)), "_1")
# rownames(colData(cds.a.b)) = paste0(rownames(colData(cds.a.b)), "_1")

# 







thisCDS = cdsList[["Activity_CDS"]]
colData(thisCDS)$Sample = colData(thisCDS)$sampleName
colData(thisCDS)$FRIP = colData(cds_p)$FRIP 
colData(thisCDS)$umi = colData(cds_p)$umi
colData(thisCDS)$log10umi = log10(colData(cds_p)$umi)
colData(thisCDS)$doublet_likelihood = colData(cds_p)$doublet_likelihood

# Save this output
# saveRDS(thisCDS, paste0(out_dir, "cds.g_with_peakMetadata", opt$ATACprocNote, "_", opt$sampleRNAname, ".rds"))




thisCDS = cdsList[["CiceroActiv_CDS"]]
colData(thisCDS)$Sample = colData(thisCDS)$sampleName
colData(thisCDS)$FRIP = colData(cds_p)$FRIP 
colData(thisCDS)$umi = colData(cds_p)$umi
colData(thisCDS)$log10umi = log10(colData(cds_p)$umi)
colData(thisCDS)$doublet_likelihood = colData(cds_p)$doublet_likelihood

# Save this output
# saveRDS(thisCDS, paste0(out_dir, "ciceroActivityCDS_with_peakMetadata", opt$ATACprocNote, "_", opt$sampleRNAname, ".rds"))















# qcInOrigFiltering = c("FRIP", "FRIT", "umi", "doublet_score", "doublet", "doublet_likelihood")
# # 
# dfTransferredColdata = as.data.frame(colData(cds.a.g))
# dfOrigFilteredColdata = as.data.frame(colData(cds_p))

# if (opt$sampleRNAname != "All_Cells"){
#   dfOrigFilteredColdata = dfOrigFilteredColdata[c(qcInOrigFiltering, "cell", genesToPlotActivity)]
#   dfOrigFilteredColdata$cell = paste0(samplesATACnames[opt$sampleRNAname], "#", dfOrigFilteredColdata$cell)
# } else {
#   dfOrigFilteredColdata = separate(data=dfOrigFilteredColdata, col = "sampleName", 
#                         into = c("sampleName", "Processing"), sep="FRIP", remove=TRUE)
#   dfOrigFilteredColdata = dfOrigFilteredColdata[c(qcInOrigFiltering, "cell", "sampleName")]
#   dfOrigFilteredColdata$cell = paste0(dfOrigFilteredColdata$sampleName, "#", dfOrigFilteredColdata$cell)
# }


# dfTransferredColdata$cell = rownames(dfTransferredColdata)


# # mergedDF = merge(dfTransferredColdata, dfOrigFilteredColdata, by="cell")
# mergedDF = join(dfTransferredColdata, dfOrigFilteredColdata, by="cell")
# rownames(mergedDF) = rownames(dfTransferredColdata)
# # Put into a new cds
# annotatedCDS_gMat = new_cell_data_set(expression_data = exprs(cds.a.g), 
#                             cell_metadata = mergedDF, 
#                             gene_metadata = as.data.frame(rowData(cds.a.g)))

# # reducedDims(annotatedCDS)$UMAP = reducedDims(cds.a.g)$UMAP
# annotatedCDS_bMat = new_cell_data_set(expression_data = exprs(cds.a.b), 
#                             cell_metadata = mergedDF, 
#                             gene_metadata = as.data.frame(rowData(cds.a.b)))



# # Get the subset of the original cells that were in the annotated DF.
# dfTransferredColdata <- tidyr::separate(dfTransferredColdata, "cell", into=c("SamplePart", "OrigLabelPart"),
#                 sep="#", remove=FALSE)
# dfTransferredColdata$OrigLabelPart = paste0(dfTransferredColdata$OrigLabelPart)

# rownames(mergedDF) = paste0(dfTransferredColdata$OrigLabelPart, "_1")
# origCDS_withLabeledCells = cds_p[, colData(cds_p)$cell %in% dfTransferredColdata$OrigLabelPart]

# dfSubsettedOrigFilteredColdata = as.data.frame(colData(origCDS_withLabeledCells))
# dfSubsettedOrigFilteredColdata  = separate(data=dfSubsettedOrigFilteredColdata, col = "sampleName", 
#                         into = c("sampleName", "Processing"), sep="FRIP", remove=TRUE)
# dfSubsettedOrigFilteredColdata$cell = paste0(dfSubsettedOrigFilteredColdata$sampleName, "#", dfSubsettedOrigFilteredColdata$cell)

# mergedDF = join(dfSubsettedOrigFilteredColdata, dfTransferredColdata, by="cell")

# rownames(mergedDF) = rownames(dfSubsettedOrigFilteredColdata)
# origCDS_withLabeledCells = new_cell_data_set(expression_data = exprs(origCDS_withLabeledCells), 
#                             cell_metadata = mergedDF, 
#                             gene_metadata = as.data.frame(rowData(origCDS_withLabeledCells)))



# # Decide which CDS to use for downstream plotting
# if (opt$featureSet == "peakMat"){
#   cdsToPlot = origCDS_withLabeledCells
#   # Re-order to have the same order as the label-transferred data
#   colnames(cdsToPlot) = cdsToPlot$cell
#   orderedCDS = cdsToPlot[,annotatedCDS_gMat$cell]
#   cdsToPlot = orderedCDS
# } else if (opt$featureSet == "gMat"){
#   cdsToPlot = annotatedCDS_gMat
# } else if (opt$featureSet == "bMat"){
#   cdsToPlot = annotatedCDS_bMat
# }






# colData(cdsToPlot)$log10umi = log10(colData(cdsToPlot)$umi)

# processingNote = paste0(processingNote, "_", opt$useMNN, "_", opt$featureSet, "_LSI", opt$lowDimOne, "_", opt$lowDimTwo)

# out_dir = paste0(out_dir, processingNote, "/")
# dir.create(out_dir)
# setwd(out_dir)


# # Try making a UMAP embedding using LSI, dropping first component, then 
# set.seed(7)
# plotDir = out_dir
# cdsToPlot = estimate_size_factors(cdsToPlot)
# cdsToPlot = preprocess_cds(cdsToPlot, method = "LSI", num_dimensions=opt$lowDimTwo)
# reducedDim(cdsToPlot) <- reducedDim(cdsToPlot)[,opt$lowDimOne:opt$lowDimTwo] 


# print("MNN now")
# if (opt$useMNN == "useMNN"){
#   # processingNote = paste0(processingNote, "MNN_Used")
#   cdsToPlot = align_cds(cdsToPlot, preprocess_method = "LSI",
#                                           alignment_group = "Sample",
#                                          residual_model_formula_str = ~log10umi + FRIP + doublet_likelihood)
#   cdsToPlot = reduce_dimension(cdsToPlot, reduction_method = 'UMAP', preprocess_method = "Aligned")
# }


# plotUMAP_Monocle(cdsToPlot, paste0("No_FACET_Peak_based", opt$sampleRNAname,  "_", processingNote), 
#             "sampleName", outputPath="./", show_labels=FALSE)

# for (eachGene in genesToPlotActivity){
#   # colData(cdsToPlot)[[eachGene]] = log2(colData(cdsToPlot)[[eachGene]])
#   plotUMAP_Monocle(cdsToPlot, paste0("No_FACET_Peak_based", opt$sampleRNAname,  "_", processingNote), 
#             eachGene, outputPath="./", show_labels=FALSE)
# }


# # Plot reduced dim coordinates
# lsiDir = "./lsiDir/"
# dir.create(lsiDir)
# for (eachLSI in 1:(opt$lowDimTwo - opt$lowDimOne)){
#   colData(cdsToPlot)[paste0("LSI_Coord_", as.character(eachLSI))] = reducedDim(cdsToPlot)[,eachLSI]
#   plotUMAP_Monocle(cdsToPlot, paste0(opt$sampleRNAname,  "_", processingNote), 
#             paste0("LSI_Coord_", as.character(eachLSI)), outputPath=lsiDir, show_labels=FALSE)
# }


# colData(cdsToPlot)$monoFragProp = colData(cdsToPlot)$nMonoFrags / colData(cdsToPlot)$nFrags 
# colData(cdsToPlot)$diFragProp = colData(cdsToPlot)$nDiFrags / colData(cdsToPlot)$nFrags 


# for (eachCol in c("sampleName", qcInOrigFiltering, "log10umi", "NucleosomeRatio", "TSSEnrichment", "PromoterRatio", "BlacklistRatio", 
#             "monoFragProp", "diFragProp")){
#   plotUMAP_Monocle(cdsToPlot, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", processingNote),
#                   eachCol, outputPath="./", show_labels=FALSE)
# }


# # Try hard-plotting this using the monocle genes function
# umapAssignedGmat = annotatedCDS_gMat
# annotatedCDS_gMat <- estimate_size_factors(annotatedCDS_gMat)
# reducedDims(annotatedCDS_gMat)$UMAP = reducedDims(cdsToPlot)$UMAP

# plotUMAP_Monocle_genes(annotatedCDS_gMat, processingNote, genesToPlotActivity,
#         "VariousMarkers", outputPath=plotDir)



# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("TTN", "FHL2"), "CardiomyoctyeMarks", outputPath="./")
# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("GSN", "CRISPLD2", "NEGR1"), "FibroblastMarks", outputPath="./")
# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("RBPJ", "LYVE1"), "MacrophageMarks", outputPath="./")
# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("LDB2", "SYNE1"), "EndothelialMarks", outputPath="./")
# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("MYH11"), "VSM_Pericyte_Marks", plotSetTotals=TRUE,  outputPath="./")

# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("THEMIS", "CD34", "CD247", "SKAP1", "CD3D", "CD3E"), "T_Cell", plotSetTotals=TRUE,  outputPath="./")

# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("DCN", "PDGFRA"), "FibroblastMarksLitvin", outputPath="./")
# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("RGS5", "ABCC9", "KCNJ8", "TAGLN", "MYH11", "ACTA2"), plotSetTotals=TRUE,  "VSM_Peri_MarksLitvin", outputPath="./")

# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("SMOC1", "NPR3" ), "EndocardMarks", outputPath="./")
# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("SEMA3G", "EFNB2", "DLL4", "NR2F2", "ACKR1", "RGCC", "CA4" ), "Vasc_Endothelial_Marks",
#                       plotSetTotals=TRUE, outputPath="./")

# plotUMAP_Monocle_genes(annotatedCDS_gMat, paste0("gMat_Activity", processingNote),
#                     c("MYH11"), "VSM_Pericyte_Marks", plotSetTotals=TRUE,  outputPath="./")





# # reducedDims(annotatedCDS)$UMAP = reducedDims(cds.a.g)$UMAP



# # }
# plotUMAP_MonocleModded <- function(dataCDS, processingNote, catToColor,
#                     show_labels=TRUE, textSize=10, outputPath = "./plots/"){ #, xVal, yVal){
#     png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored.png"),
#              width=1400, height=1000, res=200)
#     myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
#         color_cells_by=catToColor, label_cell_groups=show_labels,
#           cell_stroke=.1 , group_label_size=textSize        
#                 )) 
#     myPlot = (myPlot + theme(text=element_text(size=textSize)) +facet_wrap(~predicted.id))
#     print(myPlot)
#     dev.off()   

# }

# plotUMAP_MonocleModded(cdsToPlot, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
#                   "predicted.id", outputPath="./", show_labels=FALSE)


# kVal=20
# cdsToPlot = cluster_cells(cdsToPlot, k=kVal)

# colData(cdsToPlot)$cluster_label = clusters(cdsToPlot)
# colData(cdsToPlot)$cluster_label = as.character(colData(cdsToPlot)$cluster_label)
# colData(cdsToPlot)$partition_label = as.character(partitions(cdsToPlot))



# hardAssignATAC_clusters <- function(inputCDS, processingNote, kValUsed){

#   # 8-6-21: Using k = 15 I made a prelim assignment
#   if ((processingNote == "FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_peakMat_LSI1_50") & (kValUsed == 20)){
#     colData(inputCDS)$Assigned_Cell_Type = "Ambiguous"

#     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(22), 
#                     "Endocardium", colData(inputCDS)$Assigned_Cell_Type)
#     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(6, 23), 
#                     "Macrophage", colData(inputCDS)$Assigned_Cell_Type)
#     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(17,25), 
#                     "VSM_and_Perictye", colData(inputCDS)$Assigned_Cell_Type)
#     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(16,2,11), 
#                     "Fibroblast", colData(inputCDS)$Assigned_Cell_Type)
#     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(24,14,10,8,9,19,12), 
#                     "Vascular_Endothelium", colData(inputCDS)$Assigned_Cell_Type)
#     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(1,15,21,4,13,20,7,3,5,18), 
#                     "Cardiomyocyte", colData(inputCDS)$Assigned_Cell_Type)
#   }

#   return(inputCDS)
# }



# plotUMAP_Monocle(cdsToPlot, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
#             "cluster_label", outputPath="./", show_labels=TRUE)
# plotUMAP_Monocle(cdsToPlot, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
#             "partition_label", outputPath="./", show_labels=TRUE)



# cdsToPlot = hardAssignATAC_clusters(cdsToPlot, processingNote, kVal)
# plotUMAP_Monocle(cdsToPlot, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
#             "Assigned_Cell_Type", outputPath="./", show_labels=FALSE)







# # thisOutput = paste0(out_dir, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),  "Doublet_Likelihood_Hist.png")
# # png(thisOutput, width=1000, height=1000, res = 200)
# # myPlot = ggplot(as.data.frame(colData(annotatedCDS)) , aes_string("doublet_likelihood")) + 
# #       geom_histogram()
# # print(myPlot)
# # dev.off()






# # # Plot doublet scores since I wonder if this is important in blurring cluster borders
# # # 8-6-21 addition of this plot
# # plotDir = out_dir
# # thisOutput = paste0(out_dir, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),  "Doublet_Likelihood_Hist.png")
# # png(thisOutput, width=1000, height=1000, res = 200)
# # myPlot = ggplot(as.data.frame(colData(annotatedCDS)) , aes_string("doublet_likelihood")) + 
# #       geom_histogram()
# # print(myPlot)
# # dev.off()


# # # # }
# # # plotUMAP_MonocleModded <- function(dataCDS, processingNote, catToColor,
# # #                     show_labels=TRUE, textSize=10, outputPath = "./plots/"){ #, xVal, yVal){
# # #     png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored.png"),
# # #              width=1400, height=1000, res=200)
# # #     myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
# # #         color_cells_by=catToColor, label_cell_groups=show_labels,
# # #           cell_stroke=.1 , group_label_size=textSize        
# # #                 )) 
# # #     myPlot = (myPlot + theme(text=element_text(size=textSize)) +facet_wrap(~predicted.id))
# # #     print(myPlot)
# # #     dev.off()   

# # # }

# # plotUMAP_MonocleModded(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                   "predicted.id", outputPath="./", show_labels=FALSE)


# # #shuffle
# # # annotatedCDS = annotatedCDS[sample(1:nrow(annotatedCDS)), sample(1:ncol(annotatedCDS))]


# # colData(annotatedCDS)$log10UMI = log10(colData(annotatedCDS)$umi)
# # colData(annotatedCDS)$monoFragProp = colData(annotatedCDS)$nMonoFrags / colData(annotatedCDS)$nFrags 
# # colData(annotatedCDS)$diFragProp = colData(annotatedCDS)$nDiFrags / colData(annotatedCDS)$nFrags 


# # for (eachCol in c("sampleName", qcInOrigFiltering, "log10UMI", "NucleosomeRatio", "TSSEnrichment", "PromoterRatio", "BlacklistRatio", 
# #             "monoFragProp", "diFragProp")){

# #   plotUMAP_MonocleModded(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                   eachCol, outputPath="./", show_labels=FALSE)
# # }


# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("TTN", "FHL2"), "CardiomyoctyeMarks", outputPath="./")
# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("GSN", "CRISPLD2", "NEGR1"), "FibroblastMarks", outputPath="./")
# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("RBPJ", "LYVE1"), "MacrophageMarks", outputPath="./")
# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("LDB2", "SYNE1"), "EndothelialMarks", outputPath="./")
# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("MYH11"), "VSM_Pericyte_Marks", plotSetTotals=TRUE,  outputPath="./")

# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("THEMIS", "CD34", "CD247", "SKAP1", "CD3D", "CD3E"), "T_Cell", plotSetTotals=TRUE,  outputPath="./")

# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("DCN", "PDGFRA"), "FibroblastMarksLitvin", outputPath="./")
# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("RGS5", "ABCC9", "KCNJ8", "TAGLN", "MYH11", "ACTA2"), plotSetTotals=TRUE,  "VSM_Peri_MarksLitvin", outputPath="./")

# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("SMOC1", "NPR3" ), "EndocardMarks", outputPath="./")
# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("SEMA3G", "EFNB2", "DLL4", "NR2F2", "ACKR1", "RGCC", "CA4" ), "Vasc_Endothelial_Marks",
# #                       plotSetTotals=TRUE, outputPath="./")

# # plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     c("MYH11"), "VSM_Pericyte_Marks", plotSetTotals=TRUE,  outputPath="./")


# # kVal=10
# # annotatedCDS = cluster_cells(annotatedCDS, k=kVal)

# # colData(annotatedCDS)$cluster_label = clusters(annotatedCDS)
# # colData(annotatedCDS)$cluster_label = as.character(colData(annotatedCDS)$cluster_label)

# # plotUMAP_Monocle(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
# #             "cluster_label", outputPath="./", show_labels=TRUE)


# # # Single big plot showing sample source
# # # plotUMAP_MonocleModded(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# # #                   eachCol, outputPath="./", show_labels=FALSE)


# # plotUMAP_Monocle(annotatedCDS, paste0("ATAC_With_Transfer_No_Facet_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
# #             "sampleName", outputPath="./", show_labels=FALSE)



# # hardAssignATAC_clusters <- function(inputCDS, processingNote, kValUsed){

# #   # 8-6-21: Using k = 15 I made a prelim assignment
# #   if ((processingNote == "ATAC_With_Transfer_All_Cells_FRIP=0.2_FRIT=0.05UMI=1000") & (kValUsed == 15)){
# #     colData(inputCDS)$Assigned_Cell_Type = "Ambiguous"

# #     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(13), 
# #                     "Endocardium", colData(inputCDS)$Assigned_Cell_Type)
# #     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(7, 30), 
# #                     "Macrophage", colData(inputCDS)$Assigned_Cell_Type)
# #     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(5), 
# #                     "VSM_and_Perictye", colData(inputCDS)$Assigned_Cell_Type)
# #     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(28,33,10,11), 
# #                     "Fibroblast", colData(inputCDS)$Assigned_Cell_Type)
# #     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(24,27,14,25,32,18), 
# #                     "Vascular_Endothelium", colData(inputCDS)$Assigned_Cell_Type)
# #     colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(1,20,21,3,31,2,26,9,15,19,4,34,23,6), 
# #                     "Cardiomyocyte", colData(inputCDS)$Assigned_Cell_Type)
# #   }

# #   return(inputCDS)
# # }


# # annotatedCDS = hardAssignATAC_clusters(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                     kVal)


# # plotUMAP_Monocle(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
# #             "Assigned_Cell_Type", outputPath="./", show_labels=FALSE)


# # # Given this  set of annotations, do we see representation of all samples within the different clusters?


# # makeFacetedProportionPlot_withFill(annotatedCDS, 
# #     paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
# #       "sampleName" , "Assigned_Cell_Type", 
# #       subsetPropsToShow=c("Cardiomyocyte", "Endocardium", "Fibroblast", "Macrophage", "Vascular_Endothelium", "VSM_and_Perictye"),
# #              outputPath="./" )


# # # Assign site/donor info to the dataframe
# # sampleInfoDF = tidyr::separate(data=as.data.frame(colData(annotatedCDS)), col = "sampleName", 
# #                         into = c("Donor", "Organ", "Anatomical_Site", "Replicate"), sep="\\.", remove=FALSE)


# # colData(annotatedCDS)$Donor = sampleInfoDF$Donor
# # colData(annotatedCDS)$Anatomical_Site = sampleInfoDF$Anatomical_Site



# # makeFacetedProportionPlot_withFill(annotatedCDS, 
# #     paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
# #       "sampleName" , "Assigned_Cell_Type", fillCol =  "Anatomical_Site", colorCol =  "Anatomical_Site",
# #       subsetPropsToShow=c("Cardiomyocyte", "Endocardium", "Fibroblast", "Macrophage", "Vascular_Endothelium", "VSM_and_Perictye"),
# #              outputPath="./" )





# # plotUMAP_Monocle(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote), 
# #             "doublet_likelihood", outputPath="./", show_labels=FALSE)





















# # plotUMAP_MonocleModded(cds.a.g, paste0("Not_annotated_ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
# #                   "predicted.id", outputPath="./", show_labels=FALSE)





