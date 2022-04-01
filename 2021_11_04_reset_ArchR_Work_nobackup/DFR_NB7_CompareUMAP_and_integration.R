
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

library(uwot)

source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

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
cds.r = cds.rna

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
# cds.rna = cds.r

rna.matrix <- counts(cds.rna)
se.rna <- CreateSeuratObject(counts = rna.matrix, assay = "RNA", project = paste0(opt$sampleRNAname, "_RNA"))
#se.rna <- subset(se.rna, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
# add metaData to Seurat object
meta.rna = as.data.frame(colData(cds.rna))
se.rna <- AddMetaData(se.rna, metadata = meta.rna)
se.rna$tech <- "rna"


################################################################################################################################################

getUMAP_GregNB6_ATAC_way <- function(inputSeObj, processingNote, sampleName){
  set.seed(7)
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

  png(paste0("Greg_ATC_Method_UMAP_", processingNote, sampleName, ".png" ), 
          width = 1000, height = 900, res=200)
  print(DimPlot(inputSeObj, reduction = "umap", group.by = "tech"))
  dev.off()

  return(inputSeObj)
}

getUMAP_GregNB6_RNA_way <- function(inputSeObj, processingNote, sampleName){
  set.seed(7)

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

  png(paste0("Greg_RNA_Method_UMAP_", processingNote, sampleName, ".png" ), 
          width = 1000, height = 900, res=200)
  print(DimPlot(inputSeObj, reduction = "umap", group.by = "highLevelCellType"))
  dev.off()

}

getUMAP_minimalistSeurat <- function(inputSeObj, processingNote, sampleName, assayName){
  set.seed(7)
    # preprocess the RNA counts matrix with seurat
  inputSeObj <- NormalizeData(inputSeObj) #, normalization.method = "LogNormalize", scale.factor = 10000)
  inputSeObj <- FindVariableFeatures(inputSeObj, selection.method = "vst", nfeatures = nrow(inputSeObj))
  all.genes <- rownames(inputSeObj)
  inputSeObj <- ScaleData(inputSeObj, features = all.genes)
  # inputSeObj <- SetAssayData(inputSeObj, assay=assay, slot="scale.data", new.data="normalized.data")
  # reduce dimensions
  inputSeObj <- RunPCA(inputSeObj) #, features = VariableFeatures(object = inputSeObj))
  # Cluster cells
  inputSeObj <- FindNeighbors(inputSeObj, dims = 1:10)
  inputSeObj <- FindClusters(inputSeObj, resolution = 0.5)
  inputSeObj <- RunUMAP(inputSeObj, dims = 1:20)

  png(paste0("Minimalist_Method_", assayName, "_UMAP_", processingNote, sampleName, ".png" ), 
          width = 1000, height = 900, res=200)
  print(DimPlot(inputSeObj, reduction = "umap", group.by = "tech"))
  dev.off()

  return(inputSeObj)
}

getMonocleUMAPfromCDS<- function(inputCDS, processingNote, sampleName, assayName){
  set.seed(7)
  inputCDS = estimate_size_factors(inputCDS)
  inputCDS = preprocess_cds(inputCDS)
  inputCDS = reduce_dimension(inputCDS)
  # Now plot
  plotUMAP_Monocle(inputCDS, paste0("MonocleProcessUMAP_", processingNote, "_", assayName, "_", sampleName),
             "sampleName", outputPath = "./", show_labels=FALSE)
  generalHeartMarkers = c("GSN", "LDB2", "KCNAB1", "TTN", "RBPJ", "THEMIS", "MYH11")
  plotUMAP_Monocle_genes(inputCDS, paste0("MonocleProcessUMAP_", processingNote, "_", assayName, "_", sampleName),
           generalHeartMarkers, "GeneralMarkers",
          outputPath = "./")
  return(inputCDS)

}




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




# DFR: New code from here on, versus NB6.5
#    
DefaultAssay(se.atac) <- "ACTIVITY" # Originally did this with the code for pre-processing ATAC fo the solo UMAP
processingNote = "CompareUMAPs"

seObjList = list(se.atac, se.rna)
names(seObjList) = c("ATAC", "RNA")





# Try a co-embedding one of two ways: First, try MNN
rnaCDS_genesInATAC = cds.rna[rowData(cds.rna)$gene_short_name %in% rowData(cds.a.g)$gene_short_name,]
colData(rnaCDS_genesInATAC)$tech = "RNA"
atacCDS_genesInRNA = cds.a.g[rowData(cds.a.g)$gene_short_name %in% rowData(cds.rna)$gene_short_name,]
colData(atacCDS_genesInRNA)$tech = "ATAC"
# Combine
atacAndRNA_cds_genesIntersected = combine_cds(list(rnaCDS_genesInATAC, atacCDS_genesInRNA))





#########################################################################################################################
# Try running parametric UMAP and embedding ATAC data based on an embedding defined for RNA data
atacAndRNA_cds_genesIntersected = estimate_size_factors(atacAndRNA_cds_genesIntersected)
atacAndRNA_cds_genesIntersected = preprocess_cds(atacAndRNA_cds_genesIntersected)

# Get the PC coordinates
rnaOnlyPostPCA_cds = atacAndRNA_cds_genesIntersected[,colData(atacAndRNA_cds_genesIntersected)$tech == "RNA"]
rnaPCmatrix = reducedDims(rnaOnlyPostPCA_cds)$PCA
rnaUMAPres = uwot::umap(rnaPCmatrix, ret_model=TRUE)

# Apply to ATAC
atacOnlyPostPCA_cds = atacAndRNA_cds_genesIntersected[,colData(atacAndRNA_cds_genesIntersected)$tech == "ATAC"]
atacPCmatrix = reducedDims(atacOnlyPostPCA_cds)$PCA
atacUMAPres = umap_transform(atacPCmatrix, rnaUMAPres)


plotCombinedUMAP <- function(rnaUMAPres, atacUMAPres, processingNote, sampleName){
  # Get ATAC into a dataframe
  atacDF = as.data.frame(atacUMAPres)
  colnames(atacDF) = c("UMAP1", "UMAP2")
  atacDF$tech = "ATAC"
  # Get from RNA
  rnaDF = as.data.frame(rnaUMAPres[["embedding"]])
  colnames(rnaDF) = c("UMAP1", "UMAP2")
  rnaDF$tech = "RNA"
  # Combine
  comboDF = rbind(atacDF, rnaDF)

  # Plot
  png(paste0("./", processingNote, "_UMAP_Parametric_Coembed_", sampleName, ".png"),
          width=1000, height=1000, res=200)
  myPlot = ggplot(comboDF, aes_string(x="UMAP1", y="UMAP2", color="tech")) + 
          ggtitle(paste0("Parametric embed atac from RNA, ", sampleName)) + geom_point()
  print(myPlot)
  dev.off()
  return(comboDF )
}

# Combine a plot of these
parametricCoembedRes = plotCombinedUMAP(rnaUMAPres, atacUMAPres, processingNote, opt$sampleRNAname)

# 7-26-21: Added. 
# Use the recorded model to 
#.  1: Embed all cells (ATAC and RNA)
#   2: Use "knn" from the 'class' package to transfer labels from the RNA data onto the ATAC data.
#   3: Put this data back into the CDS and proceed onward to co-embedding based on MNN and simultaneous UMAP embedding
library('class')

allCellsEmbed = umap_transform(reducedDims(atacAndRNA_cds_genesIntersected)$PCA, rnaUMAPres)

allCellsEmbed = as.data.frame(allCellsEmbed)
colnames(allCellsEmbed) = c("UMAP1", "UMAP2")
# Label
allCellsEmbed$highLevelCellType = colData(atacAndRNA_cds_genesIntersected)$highLevelCellType
allCellsEmbed$tech = colData(atacAndRNA_cds_genesIntersected)$tech

# KNN transfer:
trainCells = allCellsEmbed[allCellsEmbed$tech == "RNA",]
testCells  = allCellsEmbed[allCellsEmbed$tech == "ATAC",]

# KNN for ALL cells.
knnRes = class::knn(trainCells[,c("UMAP1", "UMAP2")], allCellsEmbed[,c("UMAP1", "UMAP2")], #testCells[,c("UMAP1", "UMAP2")],
                    trainCells$highLevelCellType,  k = 7)

# Plot this on the UMAP
reducedDims(atacAndRNA_cds_genesIntersected)$UMAP = as.matrix(allCellsEmbed[,c("UMAP1", "UMAP2")])
colData(atacAndRNA_cds_genesIntersected)$knnCellTypeCall = knnRes 

# Plot
plotUMAP_Monocle(atacAndRNA_cds_genesIntersected, paste0(processingNote, "_ParametricCoembedded"), 
                  "knnCellTypeCall", show_labels=FALSE,
              outputPath="./")




#########################################################################################################################








# Code for MNN coembed combining ATAC+RNA
##############################################################################################################
mnnCoembedATAC_and_RNA <- function(inputCDS, processingNote, inputSample, 
                              mnnResidStr = NULL){
  set.seed(7)
  # Need to re-estimate size factors and so on
  inputCDS = estimate_size_factors(inputCDS)
  inputCDS = preprocess_cds(inputCDS)
  inputCDS = align_cds(inputCDS, alignment_group="tech", residual_model_formula_str=mnnResidStr)
  # UMAP it
  inputCDS = reduce_dimension(inputCDS)
  # Plot
  plotUMAP_Monocle(inputCDS, paste0(processingNote, "_MNN_Coembed_", inputSample), "tech",
                  show_labels =FALSE, outputPath = "./")
  return(inputCDS)
}

# First, try MNN based co-embedding
mnnRes = mnnCoembedATAC_and_RNA(atacAndRNA_cds_genesIntersected, processingNote, opt$sampleRNAname)
# Shuffle the order
mnnRes = mnnRes[sample(1:nrow(mnnRes)), sample(1:ncol(mnnRes))]

# Make a few more plots off this
generalHeartMarkers = c("GSN", "LDB2", "KCNAB1", "TTN", "RBPJ", "THEMIS", "MYH11")
plotUMAP_Monocle_genes(mnnRes, paste0(processingNote, "_MNN_Coembed_", opt$sampleRNAname),
           generalHeartMarkers, "GeneralMarkers",outputPath = "./")

# Plot lots of QC plots
columnsToColor = c("knnCellTypeCall", "tech", "log10_umi", "NucleosomeRatio",
         "PromoterRatio", "ReadsInTSS", "TSSEnrichment", "highLevelCellType")
for (eachCol in columnsToColor){
  plotUMAP_Monocle(mnnRes, paste0(processingNote, "_MNN_Coembed_", opt$sampleRNAname),
           eachCol, outputPath = "./", show_labels=FALSE)
}
##################################################################################################################





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



# Plot lots of QC plots
columnsToColor = c("knnCellTypeCall")
for (eachCol in columnsToColor){
  plotUMAP_MonocleModded(mnnRes, paste0(processingNote, "_MNN_Coembed_FACET", opt$sampleRNAname),
           eachCol, outputPath = "./", show_labels=FALSE)
}

# Plot
plotUMAP_MonocleModded(atacAndRNA_cds_genesIntersected, paste0(processingNote, "_ParametricCoembedded"), 
                  "knnCellTypeCall", show_labels=FALSE,
              outputPath="./")


































####################################################################################
# Get the UMAP one of different ways for each
for (eachName in names(seObjList)){
  print(paste0("Working on ", eachName))
  thisSeObj = seObjList[[eachName]]

  print("Getting UMAP as Greg set up with Seurat pre-processing")
  if (eachName == "ATAC"){
    umapRes = getUMAP_GregNB6_ATAC_way(thisSeObj, processingNote, opt$sampleRNAname)
  }
  # If RNA, use Greg's method for RNA as done in NB 6
  if (eachName == "RNA"){
    umapRes = getUMAP_GregNB6_RNA_way(thisSeObj, processingNote, opt$sampleRNAname)
  }
  
  print("Getting UMAP with minimalist Seurat pre-processing")
  if (eachName == "ATAC"){
    DefaultAssay(thisSeObj) <- "ACTIVITY"
  }
  umapResult = getUMAP_minimalistSeurat(thisSeObj, processingNote, opt$sampleRNAname, eachName)

  # Now Monocle Style
  if (eachName == "RNA"){
    thisCDS = cds.rna
  } else {
    thisCDS = cds.a.g
    colData(thisCDS)$sampleName = opt$sampleRNAname
  }
  umapresult = getMonocleUMAPfromCDS(thisCDS, processingNote, opt$sampleRNAname, eachName)

}


  # png(paste0("Minimalist_Method_", eachName, "_UMAP_", processingNote, opt$sampleRNAname), 
  #       width = 1000, height = 900, res=200)
  # print(DimPlot(umapResult, reduction = "umap", group.by = "tech"))
  # dev.off()




































