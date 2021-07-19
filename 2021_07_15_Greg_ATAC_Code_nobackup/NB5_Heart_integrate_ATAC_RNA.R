# Greg Booth 2021

# run with qlogin session on multiple cores: 
# qlogin -q trapnell-login.q -l mfree=25G -pe serial 4 

# This script integrates sciATAC and sciRNA from HEART using ArchR 
# The sciATAC data is already stored in an ArchR project 
# The sciRNA data is pulled in as a cds object and converted to seurat for integration

basepath = "~/projects/ATAC/210125_sciATAC_hubmap_NOVASEQ/nobackup/analysis/split_samples/"
dir.create(paste0(basepath, "archr/results/heart/NB5/"))
out_dir = paste0(basepath, "archr/results/heart/NB5/")
setwd(out_dir)

# load requirements
suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(dplyr)
  library(ggrastr)
})

# set genome
addArchRGenome("hg19")
addArchRThreads(threads = 32) 

# load filtered ArchR project
prj = loadArchRProject(path = paste0(basepath, "archr/Heart_filtered"),
                       showLogo = FALSE)

prj <- addIterativeLSI(
  ArchRProj = prj, 
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  varFeatures = 250000,
  LSIMethod = 3,
  dimsToUse = 2:30,
  iterations = 2,
  force = TRUE)


# quickly inspect ArchR UMAP
prj <- addUMAP(
  ArchRProj = prj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

p1 <- plotEmbedding(
  ArchRProj = prj, 
  colorBy = "cellColData", 
  name = "Sample", 
  embedding = "UMAP")

plotPDF(p1, 
        name = "ArchRPlot_LSI_custom_UMAP_ATACbySample.pdf", 
        addDOC = FALSE,
        width = 5,
        height = 5)


########################
# load cds.rna and convert to seurat object 
# need to change feature names to gene_short_names for 
# matching to ATAC activity matrix.
cds.rna <- readRDS(paste0(basepath, "cds_objects/sciRNA_heart/basicCellTyped_HM10UMI=100BG=15MNN=sample_cds.RDS"))
rd = rowData(cds.rna)
cm = counts(cds.rna)
cd = colData(cds.rna)
row.names(rd) <- rd$gene_short_name
row.names(cm) <- rd$gene_short_name
cds.rna = new_cell_data_set(expression_data = cm, 
                            cell_metadata = cd, 
                            gene_metadata = rd)

# Create Seurat object for RNA dataset
rna.matrix <- counts(cds.rna)
seRNA <- CreateSeuratObject(counts = rna.matrix, assay = "RNA", project = "SC2_RNA")
seRNA <- subset(seRNA, subset = nFeature_RNA > 200)

# add metaData to Seurat object
meta.rna = as.data.frame(colData(cds.rna))
seRNA <- AddMetaData(seRNA, metadata = meta.rna)
seRNA$tech <- "rna"

##########################
# unconstrained integration of scRNA with scATAC
prj <- addGeneIntegrationMatrix(
  ArchRProj = prj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "cluster_cell_type", # transfer these labels
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un", 
  force = TRUE
)

###################
# port to monocle for dim reduction and plotting 
# archR UMAPS don't work with these data 
# are transferred labels grouped?

# create monocle UMAP from bin x cell matrix
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

# Create CDS from peak Matrix (SummarizedExperiment)
cds_b = monocle3::new_cell_data_set(
  assays(bMat)$TileMatrix, 
  cell_metadata = colData(bMat),
  gene_metadata = rd)

#############################
# alternatively, make cds from the gene sore matrix 
gMat = getMatrixFromProject(
  ArchRProj = prj, 
  useMatrix = "GeneScoreMatrix", 
  binarize = FALSE, 
  threads = getArchRThreads())

rd_g = data.frame(rowData(gMat)) %>% 
  dplyr::mutate(bin = paste(seqnames, start, end, sep = "_")) %>% 
  dplyr::select(bin, gene_short_name = name)
row.names(rd_g) <- rd_g$bin
row.names(gMat) <- rd_g$bin

# Create CDS from peak Matrix (SummarizedExperiment)
cds_g = monocle3::new_cell_data_set(
  assays(gMat)$GeneScoreMatrix, 
  cell_metadata = colData(gMat),
  gene_metadata = rd_g)
##################################

# preprocess monocle3 cds

# reduce dimensions
set.seed(2017)
#cds_b <- detect_genes(cds_b)
#ncells = ncol(cds_b)*0.002
#cds_pl = cds_b[rowData(cds_b)$num_cells_expressed > ncells,]
cds_pl <- detect_genes(cds_g)
cds_pl = cds_pl[rowData(cds_pl)$num_cells_expressed > 0,]
cds_pl <- estimate_size_factors(cds_pl)
#cds_pl = preprocess_cds(cds_pl, method = "LSI", num_dimensions=50)
cds_pl = preprocess_cds(cds_pl, method = "PCA", num_dimensions=50)
### 
#alignment and MNN correction
#cds_pl = align_cds(cds_pl, preprocess_method = "LSI", residual_model_formula_str = "~log(nFrags)", alignment_group="Sample")
cds_pl = align_cds(cds_pl, preprocess_method = "PCA", residual_model_formula_str = "~log(nFrags)", alignment_group="Sample")

####
#reducedDim(cds_pl) <- reducedDim(cds_pl)[,2:50]
cds_pl = reduce_dimension(cds_pl, reduction_method = 'UMAP', preprocess_method = "Aligned")
cds_pl = cluster_cells(cds_pl, reduction_method = 'UMAP')

# prepare relevant custom UMAPs 
colData(cds_pl)$UMAP_1 <- reducedDims(cds_pl)$UMAP[,1]
colData(cds_pl)$UMAP_2 <- reducedDims(cds_pl)$UMAP[,2]
colData(cds_pl)$cluster <- cds_pl@clusters$UMAP$clusters
TCdat_pr = data.frame(colData(cds_pl))

# by monocle cluster
pdf("Monocle3_gMat_UMAP_AlignedKNN_samples.pdf", width = 4, height = 3.75)
ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = Sample)) +
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

###########
# quick plots to color cells by markers
mGenes = c(
  #cardiomyocyte_marks
  "MYL4","MYL7", "SLN",
  # endothelial_marks
  "PECAM1", "SELE" , "KDR", 
  # fibroblast makers
  "FN1") 

pdf("Monocle3_gMat_UMAP_markerGenes.pdf", width = 8, height = 8)
print(
  plot_cells(cds_pl, genes = mGenes,
             label_cell_groups = FALSE,
             label_groups_by_cluster = FALSE)
)
dev.off()

###########



### UMAP of RNA data from heart ####
cds.rna = cds.rna
# has alread been processed and dims reduced (by David)

# prepare relevant custom UMAPs 
colData(cds.rna)$UMAP_1 <- reducedDims(cds.rna)$UMAP[,1]
colData(cds.rna)$UMAP_2 <- reducedDims(cds.rna)$UMAP[,2]
TCdat_rna = data.frame(colData(cds.rna))

# by monocle cluster
pdf("Monocle3_sciRNA_UMAP_cellAnno.pdf", width = 4, height = 3.75)
ggplot(TCdat_rna, aes(x = UMAP_1, y = UMAP_2, color = cluster_cell_type)) +
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


############################
# try integration with just one sample (W144.heart.apex) from each data set:
############################

prj_W144_heart_apex = prj[prj$Sample == "W144.heart.apex.s1",]
cds.atac_W144_heart_apex = cds_b[, colData(cds_b)$Sample == "W144.heart.apex.s1"]
cds.rna_W144_heart_apex = cds.rna[,colData(cds.rna)$sampleName == "W144.Apex"]

prj_W144_heart_apex <- addIterativeLSI(
  ArchRProj = prj_W144_heart_apex, 
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  varFeatures = 100000,
  LSIMethod = 3,
  dimsToUse = 2:30,
  iterations = 2,
  force = TRUE)

# quick look at ArchR UMAP (W144.heart.apex only)
prj_W144_heart_apex <- addUMAP(
  ArchRProj = prj_W144_heart_apex, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

# cluster cells 
prj_W144_heart_apex <- addClusters(
  input = prj_W144_heart_apex,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

p2 <- plotEmbedding(
  ArchRProj = prj_W144_heart_apex, 
  colorBy = "cellColData", 
  name = "Clusters", 
  embedding = "UMAP")

plotPDF(p2, 
        name = "ArchRPlot_LSIcustom_UMAP_ATAC_W144apex_clusters.pdf", 
        addDOC = FALSE,
        width = 5,
        height = 5)

# add imputed gene scores (MAGIC-based)
prj_W144_heart_apex <- addImputeWeights(prj_W144_heart_apex)

#quick look at monocle3 UMAPs (W144.heart.apex only)
set.seed(2017)
cds.atac_W144_heart_apex <- detect_genes(cds.atac_W144_heart_apex)
ncells = ncol(cds.atac_W144_heart_apex)*0.002
cds_pl2 = cds.atac_W144_heart_apex[rowData(cds.atac_W144_heart_apex)$num_cells_expressed > ncells,]
cds_pl2 <- estimate_size_factors(cds_pl2)
cds_pl2 = preprocess_cds(cds_pl2, method = "LSI", num_dimensions=50)
reducedDim(cds_pl2) <- reducedDim(cds_pl2)[,2:50]
cds_pl2 = reduce_dimension(cds_pl2, reduction_method = 'UMAP', preprocess_method = "LSI")
cds_pl2 = cluster_cells(cds_pl2, reduction_method = 'UMAP')

# prepare relevant custom UMAPs 
colData(cds_pl2)$UMAP_1 <- reducedDims(cds_pl2)$UMAP[,1]
colData(cds_pl2)$UMAP_2 <- reducedDims(cds_pl2)$UMAP[,2]
colData(cds_pl2)$cluster <- cds_pl2@clusters$UMAP$clusters
TCdat_pr2 = data.frame(colData(cds_pl2))

# by monocle cluster
pdf("Monocle3_bMat_UMAP_W144Apex_clusters.pdf", width = 2, height = 1.75)
ggplot(TCdat_pr2, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
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
cds.rna_W144_heart_apex <- detect_genes(cds.rna_W144_heart_apex)
cds_pl2 <- estimate_size_factors(cds.rna_W144_heart_apex)
cds_pl2 = preprocess_cds(cds_pl2, method = "PCA", num_dim = 50)
cds_pl2 = reduce_dimension(cds_pl2, preprocess_method = "PCA")
cds_pl2 = cluster_cells(cds_pl2, reduction_method = 'UMAP')

# prepare relevant custom UMAPs 
colData(cds_pl2)$UMAP_1 <- reducedDims(cds_pl2)$UMAP[,1]
colData(cds_pl2)$UMAP_2 <- reducedDims(cds_pl2)$UMAP[,2]
colData(cds_pl2)$cluster <- cds_pl2@clusters$UMAP$clusters
TCdat_rna2 = data.frame(colData(cds_pl2))

# by monocle cluster
pdf("Monocle3_sciRNA_UMAP_W144apex_cellAnno.pdf", width = 2.5, height = 1.75)
ggplot(TCdat_rna2, aes(x = UMAP_1, y = UMAP_2, color = cluster_cell_type)) +
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

### Use ArchR to Transfer Cell type labels to ATAC cells from RNA cells #######

# convert RNA cds to seurat object
rd = rowData(cds.rna_W144_heart_apex)
cm = counts(cds.rna_W144_heart_apex)
cd = colData(cds.rna_W144_heart_apex)
row.names(rd) <- rd$gene_short_name
row.names(cm) <- rd$gene_short_name
cds.rna_W144_heart_apex_2 = new_cell_data_set(expression_data = cm, 
                            cell_metadata = cd, 
                            gene_metadata = rd)

# Create Seurat object for RNA dataset
rna.matrix <- counts(cds.rna_W144_heart_apex_2)
seRNA <- CreateSeuratObject(counts = rna.matrix, assay = "RNA", project = "SC2_RNA")
#seRNA <- subset(seRNA, subset = nFeature_RNA > 200)

# add metaData to Seurat object
meta.rna = as.data.frame(colData(cds.rna_W144_heart_apex_2))
seRNA <- AddMetaData(seRNA, metadata = meta.rna)
seRNA$tech <- "rna"

##########################
# unconstrained integration of scRNA with scATAC
prj_W144_heart_apex <- addGeneIntegrationMatrix(
  ArchRProj = prj_W144_heart_apex, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "cluster_cell_type", # transfer these labels
  nameCell = "predictedCell_Un",
  nameGroup = "predictedAnno_Un",
  nameScore = "predictedScore_Un", 
  force = TRUE
)




