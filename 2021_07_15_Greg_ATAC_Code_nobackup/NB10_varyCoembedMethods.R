
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
# dir.create(paste0(basepath, "archr/results/NB6/"))
out_dir = paste0(basepath, "archr/results/NB6/")
setwd(out_dir)

# load requirements
suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(dplyr)
  library(ggrastr)
  library(Seurat)
})

library(uwot)

source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

library("optparse")
# Get the passed parameters

option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
        default="HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes", 
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-s", "--sampleRNAname"), type="character",  default="W144.Apex",  # default="W142.Left.Vent",  # "All_Cells"
              help="Sample to integrate", metavar="character"),
    make_option(c("-u", "--useSubset"), type="character", default="Subset", 
              help="Subset or All", metavar="character"),

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


##
# load cds.rna 
# cds.rna <- readRDS(paste0(basepath, "cds_objects/sciRNA_heart/basicCellTyped_HM10UMI=100BG=15MNN=sample_cds.RDS"))
rnaPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/"
inputRNA_CDSname = "allCells_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes.rds"
cds.rna = readRDS(paste0(rnaPath, inputRNA_CDSname))
#adjust rownames to match ATAC geneScore matrix
rd = rowData(cds.rna)
cm = counts(cds.rna)
cd = colData(cds.rna)
row.names(rd) <- rd$gene_short_name
row.names(cm) <- rd$gene_short_name
cds.rna = new_cell_data_set(expression_data = cm, 
                            cell_metadata = cd, 
                            gene_metadata = rd)
# Alias, referred to interchaneably in the original NB6
cds.r = cds.rna

# Get the saved CDS for g matrix
rdsPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/archr/results/NB9/"
fileName = paste0("cds.g_with_peakMetadata", atacProcNote, ".rds")

# Read in the data
cds.a.g = readRDS(paste0(rdsPath, fileName))

# bMatFile = paste0("All_Cells", "_", opt$ATACprocNote, "bMatrixCDS_postTransfer.RDS")
# cds.a.b = readRDS(paste0(rdsPath, bMatFile))


# Get the subsets for these ones only
if (!(opt$sampleRNAname == "All_Cells")){
  cds.a.g = cds.a.g[, colData(cds.a.g)$Sample == samplesATACnames[opt$sampleRNAname]]
  cds.rna = cds.rna[,colData(cds.rna)$sampleName == opt$sampleRNAname]
  # cds.a.b = cds.a.b[, colData(cds.a.b)$Sample == samplesATACnames[opt$sampleRNAname]]
} 



# Make an output directory, if it doesn't exist
out_dir = paste0(basepath, "archr/results/NB10/")
dir.create(out_dir)
out_dir = paste0(out_dir, opt$sampleRNAname, "/")
dir.create(out_dir)
setwd(out_dir)

################################################################################################################################################


# }
plotUMAP_MonocleModded <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, textSize=10, outputPath = "./plots/"){ #, xVal, yVal){
    png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored.png"),
             width=1400, height=1000, res=200)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels,
          cell_stroke=.1 , group_label_size=textSize        
                )) 
    myPlot = (myPlot + theme(text=element_text(size=textSize)) +facet_wrap(~tech))
    print(myPlot)
    dev.off()   

}


# Try a co-embedding one of two ways: First, try MNN
rnaCDS_genesInATAC = cds.rna[rowData(cds.rna)$gene_short_name %in% rowData(cds.a.g)$gene_short_name,]
colData(rnaCDS_genesInATAC)$tech = "RNA"
atacCDS_genesInRNA = cds.a.g[rowData(cds.a.g)$gene_short_name %in% rowData(cds.rna)$gene_short_name,]
colData(atacCDS_genesInRNA)$tech = "ATAC"
# Combine

###########################################################################################
# Plot expression distributions

getPercentileVsMeanPlot <- function(observedCounts, eachSample, outDir, percentilesToPlot,
           processingNote){
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
        ggtitle(paste0("Percentile vs mean for ", as.character(eachSample))) + 
        geom_point() + ylab("Observed mean counts")
  print(myPlot)
  dev.off()
}

getPercentileViolinPlot <-function(observedCounts, eachSample, outDir, percentilesToPlot,
          processingNote){
  # Sort observed Counts
  observedCounts = observedCounts[order(rowSums(observedCounts)),]

  # Get percentile positions
  percentilePositions = percentilesToPlot * nrow(observedCounts)

  # Get this data for violin plotting
  plotDF = data.frame(as.matrix(observedCounts[percentilePositions,]))
  plotDF$Percentile = as.character(percentilesToPlot)
  # Get into plotting form
  plotDF = melt(plotDF, id.vars = "Percentile")

  # browser()
  # Plot the violin
  png(paste0(outDir, "percentile_violins_Downsample=", as.character(eachSample), processingNote, ".png"),
        width=1000, height=1000, res=200)
  myPlot = ggplot(plotDF, aes_string(x="Percentile", y="value")) + 
        ggtitle(paste0("Percentile violins for ", as.character(eachSample))) + 
        geom_violin() + ylab("Observed counts")
  print(myPlot)
  dev.off()
}


percentilesToPlot =  c(.2, .3, .4, .5, .6, .7, .8, .9, .99, .995, .999, .9995)

getPercentileViolinPlot(exprs(rnaCDS_genesInATAC), "RNA Data", out_dir, percentilesToPlot, 
           paste0( "RNA_Based_", opt$ATACprocNote) )
getPercentileViolinPlot(exprs(atacCDS_genesInRNA), "ATAC Data", out_dir, percentilesToPlot, 
           paste0( "ATAC_Activity_Based_", opt$ATACprocNote) )

getPercentileVsMeanPlot(exprs(rnaCDS_genesInATAC), "RNA Data", out_dir, percentilesToPlot, 
           paste0( "RNA_Based_", opt$ATACprocNote) )
getPercentileVsMeanPlot(exprs(atacCDS_genesInRNA), "ATAC Data", out_dir, percentilesToPlot, 
           paste0( "ATAC_Activity_Based_", opt$ATACprocNote) )


#################################################################################################













# Combine the data and take a look at the UMAPs that result under different alignment/correction strategies
colData(rnaCDS_genesInATAC)$umi = colData(rnaCDS_genesInATAC)$n.umi
atacAndRNA_cds = combine_cds(list(rnaCDS_genesInATAC, atacCDS_genesInRNA))

colData(atacAndRNA_cds)$log10umi = log10(colData(atacAndRNA_cds)$umi)

# # Try the approach of aligning/covariate correcting 
atacAndRNA_cds = estimate_size_factors(atacAndRNA_cds)
atacAndRNA_cds = preprocess_cds(atacAndRNA_cds, method="PCA")


rnaCoPCA = atacAndRNA_cds[,colData(atacAndRNA_cds)$tech == "RNA"]
rnaCDS_genesInATAC_aligned = align_cds(rnaCoPCA,
                         preprocess_method="PCA", alignment_group = "sampleName", 
                        residual_model_formula_str = ~log10umi )

atacCoPCA = atacAndRNA_cds[,colData(atacAndRNA_cds)$tech == "ATAC"]
atacCDS_genesInRNA_aligned = align_cds(atacCoPCA,
                         preprocess_method="PCA", alignment_group = "sampleName", 
                        residual_model_formula_str =  ~log10umi + FRIP + doublet_likelihood )


reducedDims(rnaCDS_genesInATAC_aligned)$PCA = reducedDims(rnaCDS_genesInATAC_aligned)$Aligned
reducedDims(atacCDS_genesInRNA_aligned)$PCA = reducedDims(atacCDS_genesInRNA_aligned)$Aligned

combinedAfterCorr = combine_cds(list(rnaCDS_genesInATAC_aligned, atacCDS_genesInRNA_aligned))
# Hard set the reduced dim
combinedAfterCorr = estimate_size_factors(combinedAfterCorr)
reducedDims(combinedAfterCorr)$PCA = rbind(reducedDims(rnaCDS_genesInATAC_aligned)$Aligned, reducedDims(atacCDS_genesInRNA_aligned)$Aligned)

# Align by tech
combinedAfterCorr = align_cds(combinedAfterCorr, preprocess_method="PCA", alignment_group="tech")

# Reduce dim
combinedAfterCorr = reduce_dimension(combinedAfterCorr, preprocess_method="Aligned",
                            reduction_method="UMAP")

plotUMAP_MonocleModded(combinedAfterCorr, paste0(opt$sampleRNAname,"iterativeMNN_and_Covariate"),
                       "highLevelCellType", show_labels=FALSE,
                            outputPath = out_dir)


# Reduce dim
combinedAfterCorr = reduce_dimension(combinedAfterCorr, preprocess_method="PCA",
                            reduction_method="UMAP")

plotUMAP_MonocleModded(combinedAfterCorr, paste0(opt$sampleRNAname,"iterativeMNN_and_Covariate_butNotByTech"),
                     "highLevelCellType", show_labels=FALSE,
                            outputPath = out_dir)


# Try skipping the whole intermediate, separate correction set of steps
atacAndRNA_cds = align_cds(atacAndRNA_cds, preprocess_method="PCA", alignment_group="tech")
atacAndRNA_cds = reduce_dimension(atacAndRNA_cds, preprocess_method="Aligned",
                            reduction_method="UMAP")

plotUMAP_MonocleModded(atacAndRNA_cds, paste0(opt$sampleRNAname, "only_align_by_Tech"),
                             "highLevelCellType", show_labels=FALSE,
                            outputPath = out_dir)


####




















#####






# #########################################################################################################################
# # Try running parametric UMAP and embedding ATAC data based on an embedding defined for RNA data
# atacAndRNA_cds_genesIntersected = estimate_size_factors(atacAndRNA_cds_genesIntersected)
# atacAndRNA_cds_genesIntersected = preprocess_cds(atacAndRNA_cds_genesIntersected)

# # Get the PC coordinates
# rnaOnlyPostPCA_cds = atacAndRNA_cds_genesIntersected[,colData(atacAndRNA_cds_genesIntersected)$tech == "RNA"]
# rnaPCmatrix = reducedDims(rnaOnlyPostPCA_cds)$PCA
# rnaUMAPres = uwot::umap(rnaPCmatrix, ret_model=TRUE)

# # Apply to ATAC
# atacOnlyPostPCA_cds = atacAndRNA_cds_genesIntersected[,colData(atacAndRNA_cds_genesIntersected)$tech == "ATAC"]
# atacPCmatrix = reducedDims(atacOnlyPostPCA_cds)$PCA
# atacUMAPres = umap_transform(atacPCmatrix, rnaUMAPres)


# plotCombinedUMAP <- function(rnaUMAPres, atacUMAPres, processingNote, sampleName){
#   # Get ATAC into a dataframe
#   atacDF = as.data.frame(atacUMAPres)
#   colnames(atacDF) = c("UMAP1", "UMAP2")
#   atacDF$tech = "ATAC"
#   # Get from RNA
#   rnaDF = as.data.frame(rnaUMAPres[["embedding"]])
#   colnames(rnaDF) = c("UMAP1", "UMAP2")
#   rnaDF$tech = "RNA"
#   # Combine
#   comboDF = rbind(atacDF, rnaDF)

#   # Plot
#   png(paste0("./", processingNote, "_UMAP_Parametric_Coembed_", sampleName, ".png"),
#           width=1000, height=1000, res=200)
#   myPlot = ggplot(comboDF, aes_string(x="UMAP1", y="UMAP2", color="tech")) + 
#           ggtitle(paste0("Parametric embed atac from RNA, ", sampleName)) + geom_point()
#   print(myPlot)
#   dev.off()
#   return(comboDF )
# }

# # Combine a plot of these
# parametricCoembedRes = plotCombinedUMAP(rnaUMAPres, atacUMAPres, processingNote, opt$sampleRNAname)

# # 7-26-21: Added. 
# # Use the recorded model to 
# #.  1: Embed all cells (ATAC and RNA)
# #   2: Use "knn" from the 'class' package to transfer labels from the RNA data onto the ATAC data.
# #   3: Put this data back into the CDS and proceed onward to co-embedding based on MNN and simultaneous UMAP embedding
# library('class')

# allCellsEmbed = umap_transform(reducedDims(atacAndRNA_cds_genesIntersected)$PCA, rnaUMAPres)

# allCellsEmbed = as.data.frame(allCellsEmbed)
# colnames(allCellsEmbed) = c("UMAP1", "UMAP2")
# # Label
# allCellsEmbed$highLevelCellType = colData(atacAndRNA_cds_genesIntersected)$highLevelCellType
# allCellsEmbed$tech = colData(atacAndRNA_cds_genesIntersected)$tech

# # KNN transfer:
# trainCells = allCellsEmbed[allCellsEmbed$tech == "RNA",]
# testCells  = allCellsEmbed[allCellsEmbed$tech == "ATAC",]

# # KNN for ALL cells.
# knnRes = class::knn(trainCells[,c("UMAP1", "UMAP2")], allCellsEmbed[,c("UMAP1", "UMAP2")], #testCells[,c("UMAP1", "UMAP2")],
#                     trainCells$highLevelCellType,  k = 7)

# # Plot this on the UMAP
# reducedDims(atacAndRNA_cds_genesIntersected)$UMAP = as.matrix(allCellsEmbed[,c("UMAP1", "UMAP2")])
# colData(atacAndRNA_cds_genesIntersected)$knnCellTypeCall = knnRes 

# # Plot
# plotUMAP_Monocle(atacAndRNA_cds_genesIntersected, paste0(processingNote, "_ParametricCoembedded"), 
#                   "knnCellTypeCall", show_labels=FALSE,
#               outputPath="./")




# #########################################################################################################################








# # Code for MNN coembed combining ATAC+RNA
# ##############################################################################################################
# mnnCoembedATAC_and_RNA <- function(inputCDS, processingNote, inputSample, 
#                               mnnResidStr = NULL){
#   set.seed(7)
#   # Need to re-estimate size factors and so on
#   inputCDS = estimate_size_factors(inputCDS)
#   inputCDS = preprocess_cds(inputCDS)
#   inputCDS = align_cds(inputCDS, alignment_group="tech", residual_model_formula_str=mnnResidStr)
#   # UMAP it
#   inputCDS = reduce_dimension(inputCDS)
#   # Plot
#   plotUMAP_Monocle(inputCDS, paste0(processingNote, "_MNN_Coembed_", inputSample), "tech",
#                   show_labels =FALSE, outputPath = "./")
#   return(inputCDS)
# }

# # First, try MNN based co-embedding
# mnnRes = mnnCoembedATAC_and_RNA(atacAndRNA_cds_genesIntersected, processingNote, opt$sampleRNAname)
# # Shuffle the order
# mnnRes = mnnRes[sample(1:nrow(mnnRes)), sample(1:ncol(mnnRes))]

# # Make a few more plots off this
# generalHeartMarkers = c("GSN", "LDB2", "KCNAB1", "TTN", "RBPJ", "THEMIS", "MYH11")
# plotUMAP_Monocle_genes(mnnRes, paste0(processingNote, "_MNN_Coembed_", opt$sampleRNAname),
#            generalHeartMarkers, "GeneralMarkers",outputPath = "./")

# # Plot lots of QC plots
# columnsToColor = c("knnCellTypeCall", "tech", "log10_umi", "NucleosomeRatio",
#          "PromoterRatio", "ReadsInTSS", "TSSEnrichment", "highLevelCellType")
# for (eachCol in columnsToColor){
#   plotUMAP_Monocle(mnnRes, paste0(processingNote, "_MNN_Coembed_", opt$sampleRNAname),
#            eachCol, outputPath = "./", show_labels=FALSE)
# }
# ##################################################################################################################





# # }
# plotUMAP_MonocleModded <- function(dataCDS, processingNote, catToColor,
#                     show_labels=TRUE, textSize=10, outputPath = "./plots/"){ #, xVal, yVal){
#     png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored.png"),
#              width=1400, height=1000, res=200)
#     myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
#         color_cells_by=catToColor, label_cell_groups=show_labels,
#           cell_stroke=.1 , group_label_size=textSize        
#                 )) 
#     myPlot = (myPlot + theme(text=element_text(size=textSize)) +facet_wrap(~tech))
#     print(myPlot)
#     dev.off()   

# }



# # Plot lots of QC plots
# columnsToColor = c("knnCellTypeCall")
# for (eachCol in columnsToColor){
#   plotUMAP_MonocleModded(mnnRes, paste0(processingNote, "_MNN_Coembed_FACET", opt$sampleRNAname),
#            eachCol, outputPath = "./", show_labels=FALSE)
# }

# # Plot
# plotUMAP_MonocleModded(atacAndRNA_cds_genesIntersected, paste0(processingNote, "_ParametricCoembedded"), 
#                   "knnCellTypeCall", show_labels=FALSE,
#               outputPath="./")


































# ####################################################################################
# # Get the UMAP one of different ways for each
# for (eachName in names(seObjList)){
#   print(paste0("Working on ", eachName))
#   thisSeObj = seObjList[[eachName]]

#   print("Getting UMAP as Greg set up with Seurat pre-processing")
#   if (eachName == "ATAC"){
#     umapRes = getUMAP_GregNB6_ATAC_way(thisSeObj, processingNote, opt$sampleRNAname)
#   }
#   # If RNA, use Greg's method for RNA as done in NB 6
#   if (eachName == "RNA"){
#     umapRes = getUMAP_GregNB6_RNA_way(thisSeObj, processingNote, opt$sampleRNAname)
#   }
  
#   print("Getting UMAP with minimalist Seurat pre-processing")
#   if (eachName == "ATAC"){
#     DefaultAssay(thisSeObj) <- "ACTIVITY"
#   }
#   umapResult = getUMAP_minimalistSeurat(thisSeObj, processingNote, opt$sampleRNAname, eachName)

#   # Now Monocle Style
#   if (eachName == "RNA"){
#     thisCDS = cds.rna
#   } else {
#     thisCDS = cds.a.g
#     colData(thisCDS)$sampleName = opt$sampleRNAname
#   }
#   umapresult = getMonocleUMAPfromCDS(thisCDS, processingNote, opt$sampleRNAname, eachName)

# }


#   # png(paste0("Minimalist_Method_", eachName, "_UMAP_", processingNote, opt$sampleRNAname), 
#   #       width = 1000, height = 900, res=200)
#   # print(DimPlot(umapResult, reduction = "umap", group.by = "tech"))
#   # dev.off()




































