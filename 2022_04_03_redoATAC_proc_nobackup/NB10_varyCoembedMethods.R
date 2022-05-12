
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
  unModProcNote = opt$ATACprocNote
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



knnRNAtoATAC <- function(inputCDS,kToUse=7){

  # Set up a df
  allCellsEmbed = as.data.frame(colData(inputCDS))
  allCellsEmbed$UMAP1 = reducedDims(inputCDS)$UMAP[,1]
  allCellsEmbed$UMAP2 = reducedDims(inputCDS)$UMAP[,2]

  rnaCells = allCellsEmbed[allCellsEmbed$tech == "RNA",]

  knnRes = class::knn(rnaCells[,c("UMAP1", "UMAP2")], allCellsEmbed[,c("UMAP1", "UMAP2")], 
                          rnaCells$highLevelCellType, k =kToUse)
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

  plotUMAP_MonocleModded(harmonyCDS, paste0(processingSetup,  "harmonyAlign_labelKNN"),
                             "harmonyKNN_type", show_labels=FALSE,
                            outputPath = out_dir)

  # Set the path
  outFile = paste0(processingSetup,  "harmonyAligned_cds.rds")
  saveRDS(harmonyCDS, paste0(rdsPath, outFile))
}

# Transfer this onto cds p
atacOnlyPostHarmony = harmonyCDS[,colData(harmonyCDS)$tech == "ATAC"]

# Fix the formatting on the colnames
colnames(atacOnlyPostHarmony) = sub('_[^_]*$', '', colnames(atacOnlyPostHarmony))

cds_p = readRDS(paste0(rdsPath, paste0("HM10_all_heart_fullFilter_", unModProcNote, "_ATAC_cds_p.RDS")))
cds_p = cds_p[,colnames(atacOnlyPostHarmony)]

colData(cds_p)$harmonyKNN_type = colData(atacOnlyPostHarmony)$harmonyKNN_type

# UMAP and color
set.seed(7)
cds_p =estimate_size_factors(cds_p)
cds_p = preprocess_cds(cds_p, method="LSI")
cds_p = align_cds(cds_p, preprocess_method = "LSI", alignment_group = "sampleName", 
                        residual_model_formula_str = ~log10umi + FRIP + doublet_likelihood)
cds_p = reduce_dimension(cds_p, , reduction_method = 'UMAP', preprocess_method = "Aligned")
# Plot the output

plotUMAP_Monocle(cds_p, "CDS_p_with_harmonyLabelTrans", "harmonyKNN_type", show_labels = FALSE,
                      outputPath = out_dir)


outFile = paste0(processingSetup,  "harmonyLabels_cds_p.rds")
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

# plotUMAP_MonocleModded(cadherinSubtractedCDS, paste0(processingSetup,  "harmonyAlign_labelKNN_protocadherinRegressed"),
#                              "harmonyKNN_type", show_labels=FALSE,
#                             outputPath = out_dir)

outFile = paste0(processingSetup,  "harmonyAligned_protocad_reg_cds.rds")

# saveRDS(cadherinSubtractedCDS, paste0(rdsPath, outFile))

# cadherinSubtractedCDS = readRDS(paste0(rdsPath, outFile))

cadherinSubtractedCDS = knnRNAtoATAC(cadherinSubtractedCDS, kToUse=11)

# Plot
plotUMAP_MonocleModded(cadherinSubtractedCDS, paste0(processingSetup,  "harmonyAlign_labelKNN_protocadherinRegressed"),
                             "harmonyKNN_type", show_labels=FALSE,
                            outputPath = out_dir)
plotUMAP_MonocleModded(cadherinSubtractedCDS, paste0(processingSetup,  "harmonyAlign_labelKNN_protocadherinRegressed"),
                             "highLevelCellType", show_labels=FALSE,
                            outputPath = out_dir)

# harmonyCDS = readRDS(paste0(rdsPath, processingSetup,  "harmonyAligned_cds.rds"))
saveRDS(cadherinSubtractedCDS, paste0(rdsPath, outFile))



# Save a cds_p with the knn labels after regressing out protocadherin



# Transfer this onto cds p
atacOnlyRegressProto = harmonyCDS[,colData(cadherinSubtractedCDS)$tech == "ATAC"]

# Fix the formatting on the colnames
colnames(atacOnlyRegressProto) = sub('_[^_]*$', '', colnames(atacOnlyRegressProto))

cds_p = readRDS(paste0(rdsPath, paste0("HM10_all_heart_fullFilter_", unModProcNote, "_ATAC_cds_p.RDS")))
cds_p = cds_p[,colnames(atacOnlyRegressProto)]

colData(cds_p)$harmonyKNN_type = colData(atacOnlyRegressProto)$harmonyKNN_type


outFile = paste0(processingSetup,  "harmonyLabels_regressProtocad_cds_p.rds")
saveRDS(cds_p, paste0(rdsPath, outFile))


# UMAP and color
set.seed(7)
cds_p =estimate_size_factors(cds_p)
cds_p = preprocess_cds(cds_p, method="LSI")
cds_p = align_cds(cds_p, preprocess_method = "LSI", alignment_group = "sampleName", 
                        residual_model_formula_str = ~log10umi + FRIP + doublet_likelihood)
cds_p = reduce_dimension(cds_p, , reduction_method = 'UMAP', preprocess_method = "Aligned")
# Plot the output

plotUMAP_Monocle(cds_p, "CDS_p_with_harmonyLabelTrans", "harmonyKNN_type", show_labels = FALSE,
                      outputPath = out_dir)













