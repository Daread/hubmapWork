
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/"
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
library("harmony")
library('class')
# Get the passed parameters

option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
        default="HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes", 
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-s", "--sampleRNAname"), type="character",  default= "All_Cells", # "W144.Apex",  # default="W142.Left.Vent",  # "All_Cells"
              help="Sample to integrate", metavar="character"),
    make_option(c("-u", "--useSubset"), type="character", default="Subset", 
              help="Subset or All", metavar="character"),

    make_option(c("-g", "--geneActiv"), type="character", default="ArchR",  # "ArchR", "Cicero"
              help="ArchR or Cicero activity scores for genes", metavar="character"),

    make_option(c("-r", "--saveHarmonyCDS"), type="logical", default=TRUE,  # "ArchR", "Cicero"
              help="Save a cds of RNA+ATAC data after running harmony alignment", metavar="character"),

    make_option(c("-n", "--numPCs"), type="numeric", default=20,  # "ArchR", "Cicero"
              help="ArchR or Cicero activity scores for genes", metavar="character"),

    make_option(c("-a", "--ATACprocNote"), type="character", default="FRIP=0.1_FRIT=0.08UMI=1000DL=0.5", 
              help="How ATAC Cells were filtered", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# samplesATACnames = c(
#   "W134.heart.apex.s1",
#   "W135.heart.LV.s1", 
#   "W136.heart.apex.s1", "W136.heart.LV.s1", 
#   "W139.heart.apex.s1", "W139.heart.LV.s1", "W139.heart.RV.s1", "W139.heart.septum.s1", 
#   "W142.heart.LV.s1", 
#   "W144.heart.apex.s1", 
#   "W145.heart.apex.s1", "W145.heart.LV.s1", 
#   "W146.heart.apex.s1", "W146.heart.LV.s1")
# names(samplesATACnames) = c("W134.Apex",
#   "W135.Left.Vent", 
#   "W136.Apex", "W136.Left.Vent", 
#   "W139.Apex", "W139.Left.Vent", "W139.Right.Vent", "W139.Septum", 
#   "W142.Left.Vent", 
#   "W144.Apex", 
#   "W145.Apex", "W145.Left.Vent", 
#   "W146.Apex", "W146.Left.Vent")

PPIval = 300


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

# }
plotUMAP_MonocleModded_predicted <- function(dataCDS, processingNote, catToColor,
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

hardAssignAnatomicalSites <- function(inputCDS){
  colData(inputCDS)$Anatomical_Site = "Not_Specified"
  colData(inputCDS)$Anatomical_Site = ifelse(grepl("Apex", colData(inputCDS)$sample),
                    "Apex", colData(inputCDS)$Anatomical_Site)
  colData(inputCDS)$Anatomical_Site = ifelse(grepl("Septum", colData(inputCDS)$sample),
                    "Septum", colData(inputCDS)$Anatomical_Site)
  colData(inputCDS)$Anatomical_Site = ifelse(grepl("Left.Vent", colData(inputCDS)$sample),
                    "Left_Vent", colData(inputCDS)$Anatomical_Site)
  colData(inputCDS)$Anatomical_Site = ifelse(grepl("Right.Vent", colData(inputCDS)$sample),
                    "Right_Vent", colData(inputCDS)$Anatomical_Site)
  return(inputCDS)
}

hardAssignDonors <- function(inputCDS){
  colData(inputCDS)$Donor = "Not_Specified"
  colData(inputCDS)$Donor = ifelse(grepl("W134", colData(inputCDS)$sample),
                    "W134", colData(inputCDS)$Donor)
  colData(inputCDS)$Donor = ifelse(grepl("W135", colData(inputCDS)$sample),
                    "W135", colData(inputCDS)$Donor)
  colData(inputCDS)$Donor = ifelse(grepl("W136", colData(inputCDS)$sample),
                    "W136", colData(inputCDS)$Donor)
  colData(inputCDS)$Donor = ifelse(grepl("W139", colData(inputCDS)$sample),
                    "W139", colData(inputCDS)$Donor)
  colData(inputCDS)$Donor = ifelse(grepl("W142", colData(inputCDS)$sample),
                    "W142", colData(inputCDS)$Donor)
  colData(inputCDS)$Donor = ifelse(grepl("W144", colData(inputCDS)$sample),
                    "W144", colData(inputCDS)$Donor)
  colData(inputCDS)$Donor = ifelse(grepl("W145", colData(inputCDS)$sample),
                    "W145", colData(inputCDS)$Donor)
  colData(inputCDS)$Donor = ifelse(grepl("W146", colData(inputCDS)$sample),
                    "W146", colData(inputCDS)$Donor)
  return(inputCDS)
}


regressProtocadherinAndReumap <- function(inputCDS, opt){
  # Get some sums of protocadherins and UGT1A transcripts (which have also popped up as a nuisance gene set)
  pcdhaSet = c("PCDHA6", "PCDHA5", "PCDHA2")
  colData(inputCDS)$pcdhaSet = colSums(exprs(inputCDS[rowData(inputCDS)$gene_short_name %in% pcdhaSet,]))
  #
  ugt1aSet = c("UGT1A9", "UGT1A7", "UGT1A6")
  colData(inputCDS)$ugt1aSet = colSums(exprs(inputCDS[rowData(inputCDS)$gene_short_name %in% ugt1aSet,]))
  #
  pcdhgSet = c("PCDHGA1", "PCDHGA2", "PCDHGB1", "PCDHGB2", "PCDHGC3", "PCDHGC4")
  colData(inputCDS)$pcdhgSet = colSums(exprs(inputCDS[rowData(inputCDS)$gene_short_name %in% pcdhgSet,]))

  # Regress these out
  # inputCDS = align_cds(inputCDS, residual_model_formula_str = ~pcdhaSet + ugt1aSet + pcdhgSet + n.umi)
  inputCDS = align_cds(inputCDS, preprocess_method="PCA", residual_model_formula_str = ~pcdhaSet + ugt1aSet + pcdhgSet )
  # UMAP again
  inputCDS = reduce_dimension(inputCDS, preprocess_method="Aligned",
                            reduction_method="UMAP" )

  return (inputCDS)
}



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
if (opt$geneActiv == "ArchR"){
  rdsPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/cds_objects/"
  fileName = paste0("HM10_all_heart_fullFilter_", opt$ATACprocNote, "_ATAC_cds_g.RDS")
  opt$ATACprocNote = paste0(opt$ATACprocNote, "_ArchRActiv")
} else {
  print("Using Cicero")
  # rdsPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/archr/results/NB9/All_Cells/ciceroActivityCDS_with_peakMetadata"
  # fileName = paste0(opt$ATACprocNote, "_All_Cells.rds")
  # opt$ATACprocNote = paste0(opt$ATACprocNote, "_CiceroActiv")
}
# Read in the data
cds.a.g = readRDS(paste0(rdsPath, fileName))

# bMatFile = paste0("All_Cells", "_", opt$ATACprocNote, "bMatrixCDS_postTransfer.RDS")
# cds.a.b = readRDS(paste0(rdsPath, bMatFile))
# # Get the subsets for these ones only
# if (!(opt$sampleRNAname == "All_Cells")){
#   cds.a.g = cds.a.g[, colData(cds.a.g)$Sample == samplesATACnames[opt$sampleRNAname]]
#   cds.rna = cds.rna[,colData(cds.rna)$sampleName == opt$sampleRNAname]
#   # cds.a.b = cds.a.b[, colData(cds.a.b)$Sample == samplesATACnames[opt$sampleRNAname]]
# } 



# Make an output directory, if it doesn't exist
out_dir = paste0(basepath, "archr/results/NB10/")
dir.create(out_dir)
out_dir = paste0(out_dir, opt$sampleRNAname, "/")
dir.create(out_dir)
setwd(out_dir)



# Get the genes in common
rnaCDS_genesInATAC = cds.rna[rowData(cds.rna)$gene_short_name %in% rowData(cds.a.g)$gene_short_name,]
colData(rnaCDS_genesInATAC)$tech = "RNA"
atacCDS_genesInRNA = cds.a.g[rowData(cds.a.g)$gene_short_name %in% rowData(cds.rna)$gene_short_name,]
colData(atacCDS_genesInRNA)$tech = "ATAC"
# Combine
colData(cds.rna)$tech = "RNA"
colData(cds.a.g)$tech = "ATAC"


# # Get an idea of UMI distributions for each of these
# colData(atacCDS_genesInRNA)$activSum = colSums(exprs(atacCDS_genesInRNA))
# colData(atacCDS_genesInRNA)$log10_activSum = log10(colData(atacCDS_genesInRNA)$activSum)

# # Get the sums by gene
# rowData(rnaCDS_genesInATAC)$log10_geneSum = log10(rowSums(exprs(rnaCDS_genesInATAC)))
# rowData(atacCDS_genesInRNA)$log10_geneSum = log10(rowSums(exprs(atacCDS_genesInRNA)))



knnRNAtoATAC <- function(inputCDS){

  # Set up a df
  allCellsEmbed = as.data.frame(colData(inputCDS))
  allCellsEmbed$UMAP1 = reducedDims(inputCDS)$UMAP[,1]
  allCellsEmbed$UMAP2 = reducedDims(inputCDS)$UMAP[,2]

  rnaCells = allCellsEmbed[allCellsEmbed$tech == "RNA",]

  knnRes = class::knn(rnaCells[,c("UMAP1", "UMAP2")], allCellsEmbed[,c("UMAP1", "UMAP2")], 
                          rnaCells$highLevelCellType, k =7)
  colData(inputCDS)$harmonyKNN_type = knnRes

  return(inputCDS)
}


#################################################################################################

processingSetup = paste0("PC_", as.character(opt$numPCs), "_", opt$sampleRNAname, opt$ATACprocNote)

# Combine the data and take a look at the UMAP from harmony alignment
colData(rnaCDS_genesInATAC)$umi = colData(rnaCDS_genesInATAC)$n.umi
atacAndRNA_cds = combine_cds(list(rnaCDS_genesInATAC, atacCDS_genesInRNA))

colData(atacAndRNA_cds)$log10umi = log10(colData(atacAndRNA_cds)$umi)

set.seed(7)
# # Try the approach of aligning/covariate correcting 
atacAndRNA_cds = estimate_size_factors(atacAndRNA_cds)
atacAndRNA_cds = preprocess_cds(atacAndRNA_cds, method="PCA", num_dim=opt$numPCs)

# Harmony-based:
# Try using harmony between PC and UMAP levels of reduction and see if that works better?
harmonyPCs = harmony::HarmonyMatrix(reducedDims(atacAndRNA_cds)$PCA, 
                    as.data.frame(colData(atacAndRNA_cds)), c("tech", "sampleName"), do_pca = FALSE, verbose=TRUE)

harmonyCDS = atacAndRNA_cds
reducedDims(harmonyCDS)$PCA = harmonyPCs

harmonyCDS = reduce_dimension(harmonyCDS, preprocess_method="PCA", reduction_method="UMAP")

plotUMAP_MonocleModded(harmonyCDS, paste0(processingSetup,  "harmonyAlign"),
                             "highLevelCellType", show_labels=FALSE,
                            outputPath = out_dir)

# Save this?
if (opt$saveHarmonyCDS){

  # Assign by knn
  harmonyCDS = knnRNAtoATAC(harmonyCDS)
  testCDS = knnRNAtoATAC(harmonyCDS)

  plotUMAP_MonocleModded(harmonyCDS, paste0(processingSetup,  "harmonyAlign_labelKNN"),
                             "harmonyKNN_type", show_labels=FALSE,
                            outputPath = out_dir)

  # Set the path
  outFile = paste0(processingSetup,  "harmonyAligned_cds.rds")
  saveRDS(harmonyCDS, paste0(rdsPath, outFile))
}



# Transfer this onto cds p
atacOnlyPostHarmony = harmonyCDS[,colData(harmonyCDS)$tech == "ATAC"]

cds_p = readRDS(paste0(rdsPath, paste0("HM10_all_heart_fullFilter_", opt$ATACprocNote, "_ATAC_cds_p.RDS")))
cds_p = cds_p[,colnames(atacOnlyPostHarmony)]
colData(cds_p)$harmonyAlign_labelKNN = atacOnlyPostHarmony$harmonyAlign_labelKNN

# UMAP and color
set.seed(7)
cds_p =estimate_size_factors(cds_p)
cds_p = preprocess_cds(cds_p, method="LSI")
cds_p = align_cds(cds_p, preprocess_method = "LSI", alignment_group = "sampleName", 
                        residual_model_formula_str = ~log10umi + FRIP + doublet_likelihood)
# Plot the output
plotUMAP_Monocle(cds_p, "CDS_p_with_harmonyLabelTrans", "harmonyKNN_type", show_labels = FALSE,
                      outputPath = out_dir)


outFile = paste0(processingSetup,  "harmonyAligned_cds_p.rds")
saveRDS(cds_p, paste0(rdsPath, outFile))





regressProtocadherinAndReumap <- function(inputCDS, opt){
  # Get some sums of protocadherins and UGT1A transcripts (which have also popped up as a nuisance gene set)
  pcdhaSet = c("PCDHA6", "PCDHA5", "PCDHA2")
  colData(inputCDS)$pcdhaSet = colSums(exprs(inputCDS[rowData(inputCDS)$gene_short_name %in% pcdhaSet,]))
  #
  ugt1aSet = c("UGT1A9", "UGT1A7", "UGT1A6")
  colData(inputCDS)$ugt1aSet = colSums(exprs(inputCDS[rowData(inputCDS)$gene_short_name %in% ugt1aSet,]))
  #
  pcdhgSet = c("PCDHGA1", "PCDHGA2", "PCDHGB1", "PCDHGB2", "PCDHGC3", "PCDHGC4")
  colData(inputCDS)$pcdhgSet = colSums(exprs(inputCDS[rowData(inputCDS)$gene_short_name %in% pcdhgSet,]))

  # Regress these out
  # inputCDS = align_cds(inputCDS, residual_model_formula_str = ~pcdhaSet + ugt1aSet + pcdhgSet + n.umi)
  inputCDS = align_cds(inputCDS, preprocess_method="PCA", residual_model_formula_str = ~pcdhaSet + ugt1aSet + pcdhgSet )
  # UMAP again
  inputCDS = reduce_dimension(inputCDS, preprocess_method="Aligned",
                            reduction_method="UMAP" )

  return (inputCDS)
}



# Test out also regressing out repetitive gene signal
cadherinSubtractedCDS = regressProtocadherinAndReumap(harmonyCDS)

plotUMAP_MonocleModded(cadherinSubtractedCDS, paste0(processingSetup,  "harmonyAlign_labelKNN_protocadherinRegressed"),
                             "harmonyKNN_type", show_labels=FALSE,
                            outputPath = out_dir)
outFile = paste0(processingSetup,  "harmonyAligned_protocad_reg_cds.rds")
saveRDS(cadherinSubtractedCDS, paste0(rdsPath, outFile))







################################################################################################################################################

# }


# # If using the cicero CDS, need to add in gene_short_names
# if (opt$geneActiv == "Cicero"){
#   shortNameVec = rowData(cds.rna)$gene_short_name
#   names(shortNameVec) = rowData(cds.rna)$id
#   rowData(cds.a.g)$gene_short_name = shortNameVec[rowData(cds.a.g)$GeneName]
# }


###########################################################################################
# Plot expression distributions

# getPercentileVsMeanPlot <- function(observedCounts, eachSample, outDir, percentilesToPlot,
#            processingNote){
#   # Get the means
#   observedMeans = rowSums(observedCounts) / ncol(observedCounts)
#   observedMeans = sort(observedMeans)

#   # Get percentile positions
#   percentilePositions = percentilesToPlot * length(observedMeans)

#   plotDF = data.frame("Percentile" = as.character(percentilesToPlot),
#             "Gene_Mean" = observedMeans[percentilePositions])
#   # browser()
#   # Plot
#   png(paste0(outDir, "percentileVsMean_Downsample=", as.character(eachSample), processingNote, ".png"),
#         width=1000, height=1000, res=200)
#   myPlot = ggplot(plotDF, aes_string(x="Percentile", y="Gene_Mean")) + 
#         ggtitle(paste0("Percentile vs mean for ", as.character(eachSample))) + 
#         geom_point() + ylab("Observed mean counts")
#   print(myPlot)
#   dev.off()
# }

# getPercentileViolinPlot <-function(observedCounts, eachSample, outDir, percentilesToPlot,
#           processingNote){
#   # Sort observed Counts
#   observedCounts = observedCounts[order(rowSums(observedCounts)),]

#   # Get percentile positions
#   percentilePositions = percentilesToPlot * nrow(observedCounts)

#   # Get this data for violin plotting
#   plotDF = data.frame(as.matrix(observedCounts[percentilePositions,]))
#   plotDF$Percentile = as.character(percentilesToPlot)
#   # Get into plotting form
#   plotDF = melt(plotDF, id.vars = "Percentile")

#   # browser()
#   # Plot the violin
#   png(paste0(outDir, "percentile_violins_Downsample=", as.character(eachSample), processingNote, ".png"),
#         width=1000, height=1000, res=200)
#   myPlot = ggplot(plotDF, aes_string(x="Percentile", y="value")) + 
#         ggtitle(paste0("Percentile violins for ", as.character(eachSample))) + 
#         geom_violin() + ylab("Observed counts")
#   print(myPlot)
#   dev.off()
# }


# getCDS_hist <- function(metadataToPlot, eachSample, outDir, processingNote, colToPlot){
#   # Make a hist of the preferred colData field

#   png(paste0(outDir, "histogram_for_", colToPlot, "_", as.character(eachSample), processingNote, ".png"),
#         width=1000, height=1000, res=200)
#   myPlot = ggplot(as.data.frame(metadataToPlot), aes_string(x=colToPlot)) +
#               geom_histogram() + 
#                 ggtitle(paste0("Histogram for ", colToPlot, " ", eachSample))
#   print(myPlot)
#   dev.off()
# }




# Get an idea of UMI distributions for each of these
# colData(atacCDS_genesInRNA)$activSum = colSums(exprs(atacCDS_genesInRNA))
# colData(atacCDS_genesInRNA)$log10_activSum = log10(colData(atacCDS_genesInRNA)$activSum)

# # Get the sums by gene
# rowData(rnaCDS_genesInATAC)$log10_geneSum = log10(rowSums(exprs(rnaCDS_genesInATAC)))
# rowData(atacCDS_genesInRNA)$log10_geneSum = log10(rowSums(exprs(atacCDS_genesInRNA)))



# percentilesToPlot =  c(.2, .3, .4, .5, .6, .7, .8, .9, .99, .995, .999, .9995)

# getPercentileViolinPlot(exprs(rnaCDS_genesInATAC), "RNA Data", out_dir, percentilesToPlot, 
#            paste0( "RNA_Based_", opt$ATACprocNote) )
# getPercentileViolinPlot(exprs(atacCDS_genesInRNA), "ATAC Data", out_dir, percentilesToPlot, 
#            paste0( "ATAC_Activity_Based_", opt$ATACprocNote) )

# getPercentileVsMeanPlot(exprs(rnaCDS_genesInATAC), "RNA Data", out_dir, percentilesToPlot, 
#            paste0( "RNA_Based_", opt$ATACprocNote) )
# getPercentileVsMeanPlot(exprs(atacCDS_genesInRNA), "ATAC Data", out_dir, percentilesToPlot, 
#            paste0( "ATAC_Activity_Based_", opt$ATACprocNote) )

# getCDS_hist(colData(rnaCDS_genesInATAC), "RNA Data", out_dir, 
#            paste0( "RNA_Based_", opt$ATACprocNote), "log10_umi" )
# getCDS_hist(colData(atacCDS_genesInRNA), "Cicero Activitiy Data", out_dir, 
#            paste0( "CiceroActiv_Based_", opt$ATACprocNote), "log10_activSum" )

# getCDS_hist(rowData(rnaCDS_genesInATAC), "RNA Data", out_dir, 
#            paste0( "RNA_Based_", opt$ATACprocNote), "log10_geneSum" )
# getCDS_hist(rowData(atacCDS_genesInRNA), "Cicero Activitiy Data", out_dir, 
#            paste0( "CiceroActiv_Based_", opt$ATACprocNote), "log10_geneSum" )





#################################################################################################

# generateUMAP_and_qc <- function(inputCDS, outDirName, colsToColor, cdsNote, opt, pcToUse=50, kVal=20){
#   plotLabel = paste0(cdsNote, as.character(pcToUse), opt$sampleRNAname, opt$ATACprocNote, "ATAC_Alone")
#   inputCDS = estimate_size_factors(inputCDS)
#   inputCDS = preprocess_cds(inputCDS, method='PCA', num_dim=pcToUse)
#   inputCDS = align_cds(inputCDS, preprocess_method='PCA', alignment_group="sampleName",
#                                   residual_model_formula_str =  ~log10umi + FRIP + doublet_likelihood)
#   inputCDS = reduce_dimension(inputCDS, preprocess_method='Aligned',
#                             reduction_method="UMAP")

#   # Plot the outputs
#   for (eachCol in colsToColor){
#     plotUMAP_Monocle(inputCDS, plotLabel,
#                      eachCol, show_labels=FALSE,
#                             outputPath = outDirName)
#   }
#   plotUMAP_MonocleModded_predicted(inputCDS, paste0(plotLabel, "_faceted"),
#                      "sampleName", show_labels=FALSE,
#                             outputPath = outDirName)
#   plotUMAP_Monocle_genes(inputCDS, plotLabel,
#                       c("GSN", "LDB2", "KCNAB1", "RBPJ", "TTN", "THEMIS", "FHL2"), 
#                       "genMarkers", 
#                             outputPath = outDirName)

#   # 11-1-21 Add clustering and take a look
#   inputCDS = cluster_cells(inputCDS, k=kVal)
#   colData(inputCDS)$cluster_label = as.character(clusters(inputCDS))
#       # Plot
#     plotUMAP_Monocle(inputCDS, paste0(plotLabel, "k=", as.character(kVal)),
#              "cluster_label", outputPath=outDirName)
#   # Find markers within these clusters
#   myTestResClust = runDEtestingToID_markers(inputCDS, plotLabel, "cluster_label",
#                   howManyGenesToTest = 50, outputPath=outDirName)

#   return(myTestResClust)
# }


# seuratFromCDS <- function(inputCDS){
#   countData = assay(inputCDS)
#   metaDF = as.data.frame(colData(inputCDS))
#   seuratObj = CreateSeuratObject(counts=countData, project="RNA_Data", assay="RNA",
#                 meta.data = metaDF)
#   return(seuratObj)
# }



# getVariableGenes <- function(inputCDS){
#   inputSO = seuratFromCDS(inputCDS)
#   inputSO = NormalizeData(inputSO )
#   inputSO <- FindVariableFeatures(inputSO, selection.method = "vst", nfeatures = 2000)
#   return(VariableFeatures(object = inputSO))
# }


#################################################################################################

# Get some UMAPs out of just ATAC data alone, colored a few different ways

# atacAloneOutput = paste0(out_dir, "varyATAC_alone/")
# dir.create(atacAloneOutput)

# set.seed(7)

# colData(cds.a.g)$activSum = colSums(exprs(cds.a.g))
# colData(cds.a.g)$log10_activSum = log10(colData(cds.a.g)$activSum)

# colData(cds.a.g)$sample = colData(cds.a.g)$sampleName
# cds.a.g = hardAssignDonors(cds.a.g)

# pcsHere = opt$numPCs
# cdsAndDEtest = generateUMAP_and_qc(cds.a.g, atacAloneOutput, 
#             c("log10_activSum", "Donor", "sampleName", "predicted.id", "FRIP", "doublet_likelihood", "TSSEnrichment"),
#                       "fullActivity", opt, pcToUse=pcsHere)

# # Save the test result
# DEtestResFile = paste0(atacAloneOutput, "fullActivity", as.character(pcsHere), opt$sampleRNAname, opt$ATACprocNote, "DE_Testing.csv")
# write.csv(cdsAndDEtest$marker_test_res[order(cdsAndDEtest$marker_test_res$cell_group),], DEtestResFile)


# generateUMAP_and_qc(cds.a.g, atacAloneOutput, c("log10_activSum", "sampleName", "predicted.id", "FRIP", "doublet_likelihood", "TSSEnrichment"),
#                       "fullActivity", opt, pcToUse=50)
# generateUMAP_and_qc(cds.a.g, atacAloneOutput, c("log10_activSum", "sampleName", "predicted.id", "FRIP", "doublet_likelihood", "TSSEnrichment"),
#                       "fullActivity", opt, pcToUse=35)
# generateUMAP_and_qc(cds.a.g, atacAloneOutput, c("log10_activSum", "sampleName", "predicted.id", "FRIP", "doublet_likelihood", "TSSEnrichment"),
#                       "fullActivity", opt, pcToUse=20)


# generateUMAP_and_qc(atacCDS_genesInRNA, atacAloneOutput, c("log10_activSum", "sampleName", "predicted.id", "FRIP", "doublet_likelihood", "TSSEnrichment"),
#                       "genesAlsoInRNA", opt, pcToUse=50)
# generateUMAP_and_qc(atacCDS_genesInRNA, atacAloneOutput, c("log10_activSum", "sampleName", "predicted.id", "FRIP", "doublet_likelihood", "TSSEnrichment"),
#                       "genesAlsoInRNA", opt, pcToUse=10)

# Get the highly variable genes from the RNA or ATAC cds's, and try using only those as input to the PCA/UMAP process
# varFromRNA = getVariableGenes(cds.rna)
# varSubset = atacCDS_genesInRNA[rownames(atacCDS_genesInRNA) %in% varFromRNA,]


# generateUMAP_and_qc(varSubset, atacAloneOutput, c("log10_activSum", "sampleName", "predicted.id", "FRIP", "doublet_likelihood", "TSSEnrichment"),
#                       "variableInRNA", opt, pcToUse=50)
# generateUMAP_and_qc(varSubset, atacAloneOutput, c("log10_activSum", "sampleName", "predicted.id", "FRIP", "doublet_likelihood", "TSSEnrichment"),
#                       "variableInRNA", opt, pcToUse=10)

# varFromATAC = getVariableGenes(cds.a.g)
# varSubset = cds.a.g[rownames(cds.a.g) %in% varFromATAC,]


# generateUMAP_and_qc(varSubset, atacAloneOutput, c("log10_activSum", "sampleName", "predicted.id", "FRIP", "doublet_likelihood", "TSSEnrichment"),
#                       "variableInATAC", opt, pcToUse=50)
# generateUMAP_and_qc(varSubset, atacAloneOutput, c("log10_activSum", "sampleName", "predicted.id", "FRIP", "doublet_likelihood", "TSSEnrichment"),
#                       "variableInATAC", opt, pcToUse=10)

# # Try one using ATAC activity data only
# atacCDS_genesInRNA = estimate_size_factors(atacCDS_genesInRNA)
# atacCDS_genesInRNA = preprocess_cds(atacCDS_genesInRNA, method='PCA')
# atacCDS_genesInRNA = align_cds(atacCDS_genesInRNA, preprocess_method='PCA', alignment_group="sampleName",
#                                 residual_model_formula_str =  ~log10umi + FRIP + doublet_likelihood)
# atacCDS_genesInRNA = reduce_dimension(atacCDS_genesInRNA, preprocess_method='Aligned',
#                           reduction_method="UMAP")

# plotUMAP_MonocleModded(atacCDS_genesInRNA, paste0(opt$sampleRNAname, opt$ATACprocNote, "ATAC_Alone"),
#                      "log10_activSum", show_labels=FALSE,
#                             outputPath = out_dir)




# knnRNAtoATAC <- function(inputCDS){

#   # Set up a df
#   allCellsEmbed = as.data.frame(colData(inputCDS))
#   allCellsEmbed$UMAP1 = reducedDims(inputCDS)$UMAP[,1]
#   allCellsEmbed$UMAP2 = reducedDims(inputCDS)$UMAP[,2]

#   rnaCells = allCellsEmbed[allCellsEmbed$tech == "RNA",]

#   knnRes = class::knn(rnaCells[,c("UMAP1", "UMAP2")], allCellsEmbed[,c("UMAP1", "UMAP2")], 
#                           rnaCells$highLevelCellType, k =7)
#   colData(inputCDS)$harmonyKNN_type = knnRes

#   return(inputCDS)
# }



# ### #  opt$numPCs = 50


# #################################################################################################

# processingSetup = paste0("PC_", as.character(opt$numPCs), "_", opt$sampleRNAname, opt$ATACprocNote)

# # Combine the data and take a look at the UMAPs that result under different alignment/correction strategies
# colData(rnaCDS_genesInATAC)$umi = colData(rnaCDS_genesInATAC)$n.umi
# atacAndRNA_cds = combine_cds(list(rnaCDS_genesInATAC, atacCDS_genesInRNA))

# colData(atacAndRNA_cds)$log10umi = log10(colData(atacAndRNA_cds)$umi)

# set.seed(7)
# # # Try the approach of aligning/covariate correcting 
# atacAndRNA_cds = estimate_size_factors(atacAndRNA_cds)
# atacAndRNA_cds = preprocess_cds(atacAndRNA_cds, method="PCA", num_dim=opt$numPCs)

# # Harmony-based:
# # Try using harmony between PC and UMAP levels of reduction and see if that works better?
# harmonyPCs = harmony::HarmonyMatrix(reducedDims(atacAndRNA_cds)$PCA, 
#                     as.data.frame(colData(atacAndRNA_cds)), c("tech", "sampleName"), do_pca = FALSE, verbose=TRUE)

# harmonyCDS = atacAndRNA_cds
# reducedDims(harmonyCDS)$PCA = harmonyPCs

# harmonyCDS = reduce_dimension(harmonyCDS, preprocess_method="PCA", reduction_method="UMAP")

# plotUMAP_MonocleModded(harmonyCDS, paste0(processingSetup,  "harmonyAlign"),
#                              "highLevelCellType", show_labels=FALSE,
#                             outputPath = out_dir)

# # Save this?
# if (opt$saveHarmonyCDS){

#   # Assign by knn
#   harmonyCDS = knnRNAtoATAC(harmonyCDS)
#   testCDS = knnRNAtoATAC(harmonyCDS)

#   plotUMAP_MonocleModded(harmonyCDS, paste0(processingSetup,  "harmonyAlign_labelKNN"),
#                              "harmonyKNN_type", show_labels=FALSE,
#                             outputPath = out_dir)

#   # Set the path
#   outFile = paste0(processingSetup,  "harmonyAligned_cds.rds")
#   saveRDS(harmonyCDS, paste0(out_dir, outFile))
# }




# # Alternately, try using MNN & co to correct a co-embedding

# rnaCoPCA = atacAndRNA_cds[,colData(atacAndRNA_cds)$tech == "RNA"]
# rnaCDS_genesInATAC_aligned = align_cds(rnaCoPCA,
#                          preprocess_method="PCA", alignment_group = "sampleName", 
#                         residual_model_formula_str = ~log10umi )

# atacCoPCA = atacAndRNA_cds[,colData(atacAndRNA_cds)$tech == "ATAC"]
# atacCDS_genesInRNA_aligned = align_cds(atacCoPCA,
#                          preprocess_method="PCA", alignment_group = "sampleName", 
#                         residual_model_formula_str =  ~log10umi + FRIP + doublet_likelihood )


# reducedDims(rnaCDS_genesInATAC_aligned)$PCA = reducedDims(rnaCDS_genesInATAC_aligned)$Aligned
# reducedDims(atacCDS_genesInRNA_aligned)$PCA = reducedDims(atacCDS_genesInRNA_aligned)$Aligned

# combinedAfterCorr = combine_cds(list(rnaCDS_genesInATAC_aligned, atacCDS_genesInRNA_aligned))
# # Hard set the reduced dim
# combinedAfterCorr = estimate_size_factors(combinedAfterCorr)
# reducedDims(combinedAfterCorr)$PCA = rbind(reducedDims(rnaCDS_genesInATAC_aligned)$Aligned, reducedDims(atacCDS_genesInRNA_aligned)$Aligned)

# # Align by tech
# combinedAfterCorr = align_cds(combinedAfterCorr, preprocess_method="PCA", alignment_group="tech")

# # Reduce dim
# combinedAfterCorr = reduce_dimension(combinedAfterCorr, preprocess_method="Aligned",
#                             reduction_method="UMAP")

# plotUMAP_MonocleModded(combinedAfterCorr, paste0(processingSetup, "iterativeMNN_and_Covariate"),
#                        "highLevelCellType", show_labels=FALSE,
#                             outputPath = out_dir)


# # Reduce dim
# combinedAfterCorr = reduce_dimension(combinedAfterCorr, preprocess_method="PCA",
#                             reduction_method="UMAP")

# plotUMAP_MonocleModded(combinedAfterCorr, paste0(processingSetup, "iterativeMNN_and_Covariate_butNotByTech"),
#                      "highLevelCellType", show_labels=FALSE,
#                             outputPath = out_dir)


# # Try skipping the whole intermediate, separate correction set of steps
# atacAndRNA_cds = align_cds(atacAndRNA_cds, preprocess_method="PCA", alignment_group="tech")
# atacAndRNA_cds = reduce_dimension(atacAndRNA_cds, preprocess_method="Aligned",
#                             reduction_method="UMAP")

# plotUMAP_MonocleModded(atacAndRNA_cds, paste0(processingSetup,  "only_align_by_Tech"),
#                              "highLevelCellType", show_labels=FALSE,
#                             outputPath = out_dir)




# ####




# runDimRedRNA_seurat <- function(inputSO){
#   # inputSO = NormalizeData(inputSO, normalization.method = "LogNormalize", scale.factor = 10000 )
#   inputSO = NormalizeData(inputSO )
#   inputSO <- FindVariableFeatures(inputSO, selection.method = "vst", nfeatures = 2000)
#   all.genes <- rownames(inputSO)
#   inputSO <- ScaleData(inputSO, features = all.genes)
#   # reduce dimensions
#   inputSO <- RunPCA(inputSO, features = VariableFeatures(object = inputSO), ncomponents=10)
#   # Cluster cells
#   inputSO <- FindNeighbors(inputSO, dims = 1:10)
#   inputSO <- FindClusters(inputSO, resolution = 0.5)
#   inputSO <- RunUMAP(inputSO, dims = 1:10)

#   return(inputSO)
# }

# # runCoembed <- function(inputRNA, inputATAC){


# # }


# runSeuratStrategy <- function(inputRNAcds, inputATACcds, inputProcNote, out_dir,
#                      strategyLabel, impute="impute", activityType="ArchR"){
#   if (activityType == "Cicero"){
#     rownames(inputRNAcds) = rowData(inputRNAcds)$id
#   }
  
#   # Make a seurat object out of the cds
#   seuratRNA =  seuratFromCDS(inputRNAcds)
#   # Name the RNA rows after gene short names

#   seuratATAC = seuratFromCDS(inputATACcds)

#   # Reduce and process
#   seuratRNA = runDimRedRNA_seurat(seuratRNA)
#   seuratATAC = runDimRedRNA_seurat(seuratATAC)

#   ###########################
#   # display projections of ATAC and RNA datasets
#   p1 <- DimPlot(seuratATAC, group.by = "tech") + ggtitle(paste0("Gene activity ", inputProcNote))
#   p2 <- DimPlot(seuratRNA, group.by = "tech") + ggtitle("scRNA-seq")
#   png(paste0(out_dir, "seurat_UMAPs_RNA_ATAC", opt$sampleRNAname, "_", inputProcNote, strategyLabel, ".png"), 
#           width=1600, height=600, res=200)
#   myPlot = CombinePlots(plots = list(p1, p2))
#   print(myPlot)
#   dev.off()

#   ##########################
#   # transfer labels 
#   transfer.anchors <- FindTransferAnchors(
#     reference = seuratRNA, 
#     query = seuratATAC, 
#     features = VariableFeatures(object = seuratRNA), 
#     reference.assay = "RNA", 
#     query.assay = "RNA", 
#     reduction = "cca")

#   CellType_preds <- TransferData(
#     anchorset = transfer.anchors, 
#     refdata = seuratRNA$highLevelCellType,
#     weight.reduction = seuratATAC[["pca"]])

#   seuratATAC <- AddMetaData(seuratATAC, metadata = CellType_preds)

#   ###########################
#   # display projections of ATAC and RNA datasets (With transferred cell type annotations)
#   p1 <- DimPlot(seuratATAC, group.by = "predicted.id") + ggtitle(paste0("Gene activity ", inputProcNote))
#   p2 <- DimPlot(seuratRNA, group.by = "highLevelCellType") + ggtitle("scRNA-seq")
#   png(paste0(out_dir, "seurat_UMAPs_RNA_ATAC_CTAnnotations", opt$sampleRNAname, "_", inputProcNote,strategyLabel, ".png"), 
#           width=2000, height=600, res=200)
#   myPlot = CombinePlots(plots = list(p1, p2))
#   print(myPlot)
#   dev.off()

#   ##########################

#   genes.use <- VariableFeatures(seuratRNA)
#   refdata <- GetAssayData(seuratRNA, assay = "RNA", slot = "data")[genes.use, ]
#   imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = seuratATAC[["pca"]])

#   # this line adds the imputed data matrix to the pbmc.atac object
#   if (impute == "impute"){
#     seuratATAC[["RNA"]] <- imputation
#   }
#   coembed <- merge(x = seuratRNA, y = seuratATAC)

#   # Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
#   # datasets
#   # coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
#   # coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
#   # coembed <- RunUMAP(coembed, dims = 1:30)
#   # coembed$cellType <- ifelse(!is.na(coembed$highLevelCellType), coembed$highLevelCellType, coembed$predicted.id)

#   # Variation: Find with the same way I do UMAPs for the indidivudal RNA and ATAC datasets
#   coembed = NormalizeData(coembed )
#   coembed <- FindVariableFeatures(coembed, selection.method = "vst", nfeatures = 2000)
#   all.genes <- rownames(coembed)
#   coembed <- ScaleData(coembed, features = all.genes)
#   # reduce dimensions
#   coembed <- RunPCA(coembed, features = VariableFeatures(object = coembed))
#   coembed <- RunUMAP(coembed, dims = 1:10)
#   coembed$cellType <- ifelse(!is.na(coembed$highLevelCellType), coembed$highLevelCellType, coembed$predicted.id)

#   # plot
#   p1 <- DimPlot(coembed, group.by = "tech")
#   p2 <- DimPlot(coembed, group.by = "cellType")
#   png(paste0(out_dir, "coembed_ATAC_RNA_Annotated_", impute,  opt$sampleRNAname, "_", inputProcNote, strategyLabel, ".png"), 
#           width=2400, height=1000, res=200)
#   myPlots = CombinePlots(list(p1, p2))
#   print(myPlots)
#   dev.off()
#   # ggsave(filename = ,
#   #      width = 12, height = 5, units="in", dpi=PPIval)
# }



# runSeuratStrategy( cds.rna, cds.a.g, opt$ATACprocNote, out_dir,
#                 "fullGenes", impute="impute", activityType=opt$geneActiv)



# runSeuratStrategy( rnaCDS_genesInATAC, atacCDS_genesInRNA, opt$ATACprocNote,
#                   out_dir, "matchingGenes", impute="impute", activityType=opt$geneActiv)
# runSeuratStrategy( rnaCDS_genesInATAC, atacCDS_genesInRNA,opt$ATACprocNote, 
#                  out_dir, "matchingGenes", impute="no_impute", activityType=opt$geneActiv)

# runSeuratStrategy( cds.rna, cds.a.g, opt$ATACprocNote, out_dir,
#                 "fullGenes", impute="no_impute", activityType=opt$geneActiv)







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




































