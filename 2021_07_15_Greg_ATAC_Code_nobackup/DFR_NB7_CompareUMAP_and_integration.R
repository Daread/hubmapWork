
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
dir.create(paste0(basepath, "archr/results/NB6/"))
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

library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
  			default="HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes", 
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-s", "--sampleRNAname"), type="character", default="W144.Apex", 
              help="Sample to integrate", metavar="character")
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

# Get the ATAC cds
cds.a.g = readRDS("HM10_all_heart_apex_ATAC_cds_g.RDS")
cds.a.b = readRDS("HM10_all_heart_apex_ATAC_cds_b.RDS")

# Get the subsets for these ones only
cds.a.g = cds.a.g[, colData(cds.a.g)$Sample == samplesATACnames[opt$sampleRNAname]]
cds.rna = cds.rna[,colData(cds.rna)$sampleName == opt$sampleRNAname]
cds.a.b = cds.a.b[, colData(cds.a.b)$Sample == samplesATACnames[opt$sampleRNAname]]
cds.atac.g = cds.a.g
cds.atac.b = cds.a.b

# Make an output directory, if it doesn't exist
out_dir = paste0(basepath, "archr/results/NB7/")
dir.create(out_dir)
out_dir = paste0(out_dir, opt$sampleRNAname, "/")
dir.create(out_dir)
setwd(out_dir)

################################################################################################################################################


# Integrating outside of the ArchR framework 
# Going directly to Seurat
###########
# put peak matrix in ATAC assay 
# put activity matrix in "ACTIVITY" assay
peaks = counts(cds.a.b)
activity.matrix = counts(cds.a.g)

se.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = paste0(opt$sampleRNAname, "_ATAC"))
se.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

# add metaData to Seurat object
meta = as.data.frame(colData(cds.a.b))
se.atac <- AddMetaData(se.atac, metadata = meta)
se.atac$tech <- "atac"


##########################
# convert RNA cds to seurat object
cds.rna = cds.r

rna.matrix <- counts(cds.rna)
se.rna <- CreateSeuratObject(counts = rna.matrix, assay = "RNA", project = paste0(opt$sampleRNAname, "_RNA"))
#se.rna <- subset(se.rna, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
# add metaData to Seurat object
meta.rna = as.data.frame(colData(cds.rna))
se.rna <- AddMetaData(se.rna, metadata = meta.rna)
se.rna$tech <- "rna"


################################################################################################################################################

getUMAP_GregNB6_ATAC_way <- function(inputSeObj, processingNote, sampleName){
  # preprocess to find anchors (highly variable features) between ATAC and gene activity matrices
  DefaultAssay(inputSeObj) <- "ACTIVITY"
  inputSeObj <- FindVariableFeatures(inputSeObj)
  inputSeObj <- NormalizeData(inputSeObj)
  inputSeObj <- ScaleData(inputSeObj)

  # reduce dimensions
  DefaultAssay(inputSeObj) <- "ATAC"
  VariableFeatures(inputSeObj) <- names(which(Matrix::rowSums(inputSeObj) > 10))
  inputSeObj <- RunLSI(inputSeObj, n = 50, scale.max = NULL)
  inputSeObj <- RunUMAP(inputSeObj, reduction = "lsi", dims = 2:50)

  png(paste0("Greg_ATC_Method_UMAP_", processingNote, sampleName, ".png" ), width = 5, height = 4.5)
  print(DimPlot(inputSeObj, reduction = "umap", group.by = "tech"))
  dev.off()

  return(inputSeObj)
}

getUMAP_GregNB6_RNA_way <- function(inputSeObj, processingNote, sampleName){

  # preprocess the RNA counts matrix with seurat
  inputSeObj <- NormalizeData(inputSeObj, normalization.method = "LogNormalize", scale.factor = 10000)
  inputSeObj <- FindVariableFeatures(inputSeObj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(inputSeObj)
  inputSeObj <- ScaleData(inputSeObj, features = all.genes)
  # reduce dimensions
  inputSeObj <- RunPCA(inputSeObj, features = VariableFeatures(object = inputSeObj))
  # Cluster cells
  inputSeObj <- FindNeighbors(inputSeObj, dims = 1:10)
  inputSeObj <- FindClusters(inputSeObj, resolution = 0.5)
  inputSeObj <- RunUMAP(inputSeObj, dims = 1:10)

  png(paste0("Greg_RNA_Method_UMAP_", processingNote, sampleName, ".png" ), width = 5, height = 4.5)
  print(DimPlot(inputSeObj, reduction = "umap", group.by = "highLevelCellType"))
  dev.off()

}

getUMAP_minimalistSeurat <- function(inputSeObj, processingNote, sampleName, assayName){
    # preprocess the RNA counts matrix with seurat
  inputSeObj <- NormalizeData(inputSeObj) #, normalization.method = "LogNormalize", scale.factor = 10000)
  inputSeObj <- FindVariableFeatures(inputSeObj, selection.method = "vst", nfeatures = nrow(inputSeObj))
  all.genes <- rownames(inputSeObj)
  inputSeObj <- ScaleData(inputSeObj, features = all.genes)
  # inputSeObj <- SetAssayData(inputSeObj, assay=assay, slot="scale.data", new.data="normalized.data")
  # reduce dimensions
  inputSeObj <- RunPCA(inputSeObj) #, features = VariableFeatures(object = inputSeObj))
  # Cluster cells
  inputSeObj <- FindNeighbors(inputSeObj) #, dims = 1:10)
  inputSeObj <- FindClusters(inputSeObj) #, resolution = 0.5)
  inputSeObj <- RunUMAP(inputSeObj, dims = 1:50)

  png(paste0("Minimalist_Method_", assayName, "_UMAP_", processingNote, sampleName, ".png" ), width = 5, height = 4.5)
  print(DimPlot(inputSeObj, reduction = "umap", group.by = "tech"))
  dev.off()

}

# DFR: New code from here on, versus NB6.5
#    
DefaultAssay(se.atac) <- "ACTIVITY" # Originally did this with the code for pre-processing ATAC fo the solo UMAP
processingNote = "CompareUMAPs"
seObjList = list(se.atac, se.rna)
names(seObjList) = c("ATAC", "RNA")


# Get the UMAP one of different ways for each
for (eachName in names(seObjList)){
  print(paste0("Working on ", eachName))
  thisSeObj = seObjList[[eachName]]

  # print("Getting UMAP as Greg set up with Seurat pre-processing")
  # if (eachName == "ATAC"){
  #   umapRes = getUMAP_GregNB6_ATAC_way(thisSeObj, processingNote, opt$sampleRNAname)
  # }
  # # If RNA, use Greg's method for RNA as done in NB 6
  # if (eachName == "RNA"){
  #   umapRes = getUMAP_GregNB6_RNA_way(thisSeObj, processingNote, opt$sampleRNAname)
  # }
  
  print("Getting UMAP with minimalist Seurat pre-processing")
  if (eachName == "ATAC"){
    DefaultAssay(thisSeObj) <- "ACTIVITY"
  }
  umapRes = getUMAP_minimalistSeurat(thisSeObj, processingNote, opt$sampleRNAname, eachName)

}







































