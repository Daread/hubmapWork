
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
dir.create(paste0(basepath, "archr/results/NB6/"))
out_dir = paste0(basepath, "archr/results/NB6/")
setwd(out_dir)
set.seed(7)

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

PPIval = 300

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


# Subset names

# Get the subsets for these ones only
cds.a.g = cds.a.g[, colData(cds.a.g)$Sample == samplesATACnames[opt$sampleRNAname]]
cds.rna = cds.rna[,colData(cds.rna)$sampleName == opt$sampleRNAname]
cds.a.b = cds.a.b[, colData(cds.a.b)$Sample == samplesATACnames[opt$sampleRNAname]]
cds.atac.g = cds.a.g
cds.atac.b = cds.a.b
cds.r = cds.rna

################################################################################################################################################


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

png(paste0("DimPlotATAC", opt$sampleRNAname, ".png" ), width = 5*PPIval, height = 4.5*PPIval, res=PPIval)
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


png(paste0("DimPlot_RNA", opt$sampleRNAname, ".png"), width = 5.5 * PPIval, height = 3.5 * PPIval, res=PPIval)
DimPlot(se.rna, reduction = "umap", group.by = "highLevelCellType")
dev.off()


###########################
# display projections of ATAC and RNA datasets
p1 <- DimPlot(se.atac, group.by = "tech") + ggtitle("scATAC-seq")
p2 <- DimPlot(se.rna, group.by = "tech") + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))
ggsave(filename = paste0("seurat_UMAPs_RNA_ATAC", opt$sampleRNAname, ".png"),
		 width = 8, height = 3, units="in", dpi=PPIval)




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
ggsave(filename = paste0("seurat_UMAPs_RNA_ATAC_CTAnnotations", opt$sampleRNAname, ".png"),
			 width = 10, height = 3, units="in", dpi=PPIval)

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
ggsave(filename = paste0("coembed_ATAC_RNA_Annotated", opt$sampleRNAname, ".png"),
		 width = 12, height = 5, units="in", dpi=PPIval)

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
lowAbundance = c("Unknown")

# by predicted ID
png(paste0("Monocle3_bMat_UMAP_All_Heart_Annotations", opt$sampleRNAname, ".png"),
			 width = 3*PPIval, height = 1.75*PPIval, res=PPIval)
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
lowAbundance2 = c("Unknown")

# by monocle cluster
png(paste0("Monocle3_sciRNA_UMAP_All_Heart_Annotations", opt$sampleRNAname, ".png"),
			 width = 3*PPIval, height = 1.75*PPIval, res=PPIval)
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



































