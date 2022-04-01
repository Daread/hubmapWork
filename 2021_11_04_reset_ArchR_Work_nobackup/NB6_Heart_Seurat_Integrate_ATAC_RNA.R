# Note: This only runs if you get the right resources upon qlogin/qsub. Run:
#  qlogin -q trapnell-login.q -l mfree=4G -pe serial 16
# to start up

# Trying:
# qlogin -l mfree=20G -q trapnell-login.q -pe serial 16 -l centos=7

# Trying this on 7-20-21. (Didn't work, going back to 20g)
# qlogin -l mfree=40G -q trapnell-login.q -pe serial 16 -l centos=7

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

# set genome
# addArchRGenome("hg19")
addArchRGenome("hg38")

# addArchRThreads(threads = 32) 
addArchRThreads(threads = 16) 

# load filtered ArchR project
prj = loadArchRProject(path = paste0(basepath, "archr/Heart_filtered"),
                       showLogo = FALSE)

########################
# Convert ATAC data to cds object
# both bin and gene score matrices


# Trying to run with a single thread right now
# addArchRThreads(threads = 1) 


# load PeakMatrix into memory
bMat = getMatrixFromProject(
  ArchRProj = prj, 
  useMatrix = "TileMatrix", 
  binarize = TRUE, 
  threads = getArchRThreads())

# format rowData 
rd = data.frame(rowData(bMat)) %>% 
  dplyr::mutate(stop = start + 499,
                bin = paste(seqnames, start, stop, sep = "_")) %>% 
  dplyr::select(bin)
row.names(rd) <- rd$bin
row.names(bMat) <- rd$bin

# Create CDS from tile Matrix (SummarizedExperiment)
cds_b = monocle3::new_cell_data_set(
  assays(bMat)$TileMatrix, 
  cell_metadata = colData(bMat),
  gene_metadata = rd)


# load gene score Matrix into memory
gMat = getMatrixFromProject(
  ArchRProj = prj, 
  useMatrix = "GeneScoreMatrix", 
  binarize = FALSE, 
  threads = getArchRThreads())

# format rowData 
rd = data.frame(rowData(gMat)) %>% 
  dplyr::mutate(bin = paste(seqnames, start, end, sep = "_"), 
                gene_short_name = name) %>% 
  dplyr::select(bin, name, gene_short_name)
row.names(rd) <- rd$name
row.names(gMat) <- rd$name

# Create CDS from gene score Matrix (SummarizedExperiment)
cds_g = monocle3::new_cell_data_set(
  assays(gMat)$GeneScoreMatrix, 
  cell_metadata = colData(gMat),
  gene_metadata = rd)


##
# load cds.rna 
# cds.rna <- readRDS(paste0(basepath, "cds_objects/sciRNA_heart/basicCellTyped_HM10UMI=100BG=15MNN=sample_cds.RDS"))
rnaPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/"
cds.rna = readRDS(paste0(rnaPath, "allCells_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes.rds"))
#adjust rownames to match ATAC geneScore matrix
rd = rowData(cds.rna)
cm = counts(cds.rna)
cd = colData(cds.rna)
row.names(rd) <- rd$gene_short_name
row.names(cm) <- rd$gene_short_name
cds.rna = new_cell_data_set(expression_data = cm, 
                            cell_metadata = cd, 
                            gene_metadata = rd)

# Added 7-29-21: Filter by checking names in a cds filtered by greg/riza's method
filterByGregCDS = TRUE
filterCDSName = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/cds_objects/"
processingNote = "FRIP=0.2_FRIT=0.05UMI=1000"
# processingNote = "FRIP=0.3_FRIT=0.1UMI=1000"
filterByDL = TRUE
DLcutoff = .5


if (filterByGregCDS){
  library(tidyr)
  # cds_p holds data
  load(paste0(filterCDSName, "cds_p_allHeartATAC", processingNote))
  sampleNameDF = separate(data=as.data.frame(colData(cds_p)), col = "sampleName", 
            into = c("sampleName", "procNote"), sep="FRIP", remove=TRUE)
  colData(cds_p)$sampleName = sampleNameDF$sampleName 
  sampleNameDF = NULL

  # Filter by doublet likelihood here?
  if (filterByDL){
    cds_p = cds_p[,colData(cds_p)$doublet_likelihood < DLcutoff]
    processingNote = paste0(processingNote, "DL=", as.character(DLcutoff))
  }
  namesToSave = paste0(colData(cds_p)$sampleName,"#", (colData(cds_p)$cell))

  cds_b = cds_b[,(rownames(colData(cds_b)) %in% namesToSave) ]
  cds_g = cds_g[,(rownames(colData(cds_g)) %in% namesToSave) ]
} else{
  processingNote = ""
}


############################
# Modification from Greg: Use all samples, not W144
############################
# get restricted cells and features (smaller file size)
# cds.atac.b <- cds_b[, colData(cds_b)$Sample == "W144.heart.apex.s1"]

#
processingNote = paste0("fixedHg38_", processingNote)


cds.atac.b = cds_b
cds.atac.b <- detect_genes(cds.atac.b)
ncells = ncol(cds.atac.b)*0.002
cds.atac.b = cds.atac.b[rowData(cds.atac.b)$num_cells_expressed > ncells,]
# cds.atac.g = cds_g[, colData(cds_g)$Sample == "W144.heart.apex.s1"]
# cds.rna = cds.rna[,colData(cds.rna)$sampleName == "W144.Apex"]
cds.atac.g = cds_g
cds.rna = cds.rna

# save these cds objects for future use (saves time)
saveRDS(cds.atac.b, file=paste0("HM10_all_heart", processingNote, "_ATAC_cds_b.RDS"))
saveRDS(cds.atac.g, file=paste0("HM10_all_heart", processingNote, "_ATAC_cds_g.RDS"))

saveRDS(cds.rna, file="HM10_all_heart_apex_RNA_cds.RDS")

########################################################
########################################################
# can start from here

# Load cds objects
cds.r = readRDS("HM10_all_heart_apex_RNA_cds.RDS")
cds.a.g = readRDS("HM10_all_heart_apex_ATAC_cds_g.RDS")
cds.a.b = readRDS("HM10_all_heart_apex_ATAC_cds_b.RDS")

# Integrating outside of the ArchR framework 
# Going directly to Seurat
###########
# put peak matrix in ATAC assay 
# put activity matrix in "ACTIVITY" assay
peaks = counts(cds.a.b)
activity.matrix = counts(cds.a.g)

se.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "AllHeart_ATAC")
se.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

# add metaData to Seurat object
meta = as.data.frame(colData(cds.a.b))
se.atac <- AddMetaData(se.atac, metadata = meta)
se.atac$tech <- "atac"

# preprocess to find anchors (highly variable features) between ATAC and gene activity matrices
DefaultAssay(se.atac) <- "ACTIVITY"
se.atac <- FindVariableFeatures(se.atac)
se.atac <- NormalizeData(se.atac)
se.atac <- ScaleData(se.atac)

# reduce dimensions
DefaultAssay(se.atac) <- "ATAC"
VariableFeatures(se.atac) <- names(which(Matrix::rowSums(se.atac) > 10))
se.atac <- RunLSI(se.atac, n = 50, scale.max = NULL)
se.atac <- RunUMAP(se.atac, reduction = "lsi", dims = 2:50)

pdf("DimPlot_ATAC.pdf", width = 5, height = 4.5)
DimPlot(se.atac, reduction = "umap", group.by = "tech")
dev.off()

##########################
# convert RNA cds to seurat object
cds.rna = cds.r

rna.matrix <- counts(cds.rna)
se.rna <- CreateSeuratObject(counts = rna.matrix, assay = "RNA", project = "AllHeart_RNA")
#se.rna <- subset(se.rna, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
# add metaData to Seurat object
meta.rna = as.data.frame(colData(cds.rna))
se.rna <- AddMetaData(se.rna, metadata = meta.rna)
se.rna$tech <- "rna"

# preprocess the RNA counts matrix with seurat
se.rna <- NormalizeData(se.rna, normalization.method = "LogNormalize", scale.factor = 10000)
se.rna <- FindVariableFeatures(se.rna, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(se.rna)
se.rna <- ScaleData(se.rna, features = all.genes)
# reduce dimensions
se.rna <- RunPCA(se.rna, features = VariableFeatures(object = se.rna))
# Cluster cells
se.rna <- FindNeighbors(se.rna, dims = 1:10)
se.rna <- FindClusters(se.rna, resolution = 0.5)
se.rna <- RunUMAP(se.rna, dims = 1:10)

##############################################
#only run on 



######################################



pdf("DimPlot_RNA.pdf", width = 5.5, height = 3.5)
DimPlot(se.rna, reduction = "umap", group.by = "highLevelCellType")
dev.off()


###########################
# display projections of ATAC and RNA datasets
p1 <- DimPlot(se.atac, group.by = "tech") + ggtitle("scATAC-seq")
p2 <- DimPlot(se.rna, group.by = "tech") + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))
ggsave(filename = "seurat_UMAPs_RNA_ATAC.pdf", width = 8, height = 3)


# DFR Add: CCA takes forever on full data. Save here.
saveRDS(se.rna, "HM10_All_Heart_RNA_preCCA.rds")
saveRDS(se.atac, "HM10_All_Heart_ATAC_preCCA.rds")


# Can resume here
# DFR
se.rna = readRDS("HM10_All_Heart_RNA_preCCA.rds")

se.atac = readRDS("HM10_All_Heart_ATAC_preCCA.rds")


##########################
# transfer labels 
transfer.anchors <- FindTransferAnchors(
  reference = se.rna, 
  query = se.atac, 
  features = VariableFeatures(object = se.rna), 
  reference.assay = "RNA", 
  query.assay = "ACTIVITY", 
  reduction = "cca")

CellType_preds <- TransferData(
  anchorset = transfer.anchors, 
  refdata = se.rna$highLevelCellType,
  weight.reduction = se.atac[["lsi"]])

se.atac <- AddMetaData(se.atac, metadata = CellType_preds)

###########################
# display projections of ATAC and RNA datasets (With transferred cell type annotations)
p1 <- DimPlot(se.atac, group.by = "predicted.id") + ggtitle("scATAC-seq")
p2 <- DimPlot(se.rna, group.by = "highLevelCellType") + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))
ggsave(filename = "seurat_UMAPs_RNA_ATAC_CTAnnotations.pdf", width = 10, height = 3)

##########################
# co-embedding data 
# note that we restrict the imputation to variable genes from scRNA-seq, 
# but could impute the full transcriptome if we wanted to
genes.use <- VariableFeatures(se.rna)
refdata <- GetAssayData(se.rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = se.atac[["lsi"]])

# this line adds the imputed data matrix to the pbmc.atac object
se.atac[["RNA"]] <- imputation
coembed <- merge(x = se.rna, y = se.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$cellType <- ifelse(!is.na(coembed$highLevelCellType), coembed$highLevelCellType, coembed$predicted.id)

# plot
p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "cellType")
CombinePlots(list(p1, p2))
ggsave(filename = "coembed_ATAC_RNA_Annotated.pdf", width = 12, height = 5)

###################
###################
# add cell-type annotations back into monocle cds object

atac.anno.preds = data.frame(se.atac@meta.data) %>% 
  tibble::rownames_to_column(var = "cell") %>% 
  select(cell, predicted.id)

cdat = data.frame(colData(cds.a.b)) %>%
  tibble::rownames_to_column(var = "cell") %>% 
  left_join(atac.anno.preds, by = "cell")

colData(cds.a.b)$predicted.id = cdat$predicted.id

# monocle3 UMAPs (W144.heart.apex only)
set.seed(2017)
cds_pl <- estimate_size_factors(cds.a.b)
cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
reducedDim(cds_pl) <- reducedDim(cds_pl)[,2:50]
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "LSI")
cds_pl = cluster_cells(cds_pl, reduction_method = 'UMAP')

# prepare relevant custom UMAPs 
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))

# remove low abundace cell types 
# lowAbundance = c("T Cells", "B Cells", "Neuronal") # DFR commented out 7-18-21

# by predicted ID
pdf("Monocle3_bMat_UMAP_All_Heart_Annotations.pdf", width = 3, height = 1.75)
filter(TCdat_pr, !(predicted.id %in% lowAbundance)) %>% 
ggplot(aes(x = UMAP_1, y = UMAP_2, color = predicted.id)) +
  geom_point_rast(size=0.4, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
dev.off()

### UMAP of RNA data from W144.heart.apex only####
set.seed(2017)
cds.r <- detect_genes(cds.r)
cds_pl2 <- estimate_size_factors(cds.r)
cds_pl2 = preprocess_cds(cds_pl2, method = "PCA", num_dim = 50)
cds_pl2 = reduce_dimension(cds_pl2, preprocess_method = "PCA")
cds_pl2 = cluster_cells(cds_pl2, reduction_method = 'UMAP')

# prepare relevant custom UMAPs 
colData(cds_pl2)$UMAP_1 <- reducedDims(cds_pl2)$UMAP[,1]
colData(cds_pl2)$UMAP_2 <- reducedDims(cds_pl2)$UMAP[,2]
colData(cds_pl2)$cluster <- cds_pl2@clusters$UMAP$clusters
TCdat_rna2 = data.frame(colData(cds_pl2))

# remove low abundace cell types 
# lowAbundance2 = c("20", "Adipocytes", "Neuronal") # DFR commented out 7-18-21

# by monocle cluster
pdf("Monocle3_sciRNA_UMAP_All_Heart_Annotations.pdf", width = 3, height = 1.75)
filter(TCdat_rna2, !(highLevelCellType %in% lowAbundance2)) %>% 
ggplot( aes(x = UMAP_1, y = UMAP_2, color = highLevelCellType)) +
  geom_point_rast(size=0.4, stroke = 0) + #rasterizes scatter plot while keeping axes etc. vectorized
  theme(legend.position = "right", text = element_text(size = 6),  
        legend.key.width = unit(0.2,"line"), legend.key.height = unit(0.25,"line")) + 
  xlab("Component 1") +
  ylab("Component 2") +
  guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  monocle3:::monocle_theme_opts() +
  theme(axis.text = element_blank(), axis.line = element_blank(), axis.line.x = element_blank(), 
        axis.line.y = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
dev.off()


