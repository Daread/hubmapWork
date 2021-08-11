


source("../../../sharedProjectCode/utility/singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()
print(modStatus)
# sessionInfo()

source("../../../sharedProjectCode/scRNA_Seq_Background_Code/backgroundMethods.R")
processingNote = "DanaRulaTest_6-7-21"

useScrublet = FALSE
writeMatrixForScrub = TRUE
scrubCutoff = .2
umiCellMin = 50
bgUMImax = 15

print("Starting Now")

dataPath <- "./final-output/"
#myCDS <- readRDS(paste0(dataPath, "cds.RDS"))
# myCDS <- load.cds(paste0(dataPath, "UMI.count.matrix"), paste0(dataPath, "gene.annotations"),
#          paste0(dataPath, "cell.annotations"))
# print("CDS Loaded")

samplesToRead <- c("Spleen", "Pancreas", "Lung", "Liver")

# myCDS <- readRDS(paste0(dataPath, "cdsWithHash.RDS"))
myCDS <- makeCDS_fromBBIoutputDirs(samplesToRead)

colData(myCDS)$Cell = rownames(colData(myCDS))
colData(myCDS)$cell = colData(myCDS)$Cell


# Get the background and main cell subsets
bgCDS = myCDS[,colData(myCDS)$n.umi <= bgUMImax]
fullCDS = myCDS
myCDS = myCDS[,colData(myCDS)$n.umi >= umiCellMin]
processingNote = paste0(processingNote, "UMI=", as.character(umiCellMin), 
    "BG=", as.character(bgUMImax))
cellNames = colData(myCDS)$Cell

# Color by the PCR well
splitCellCol = (separate(as.data.frame(colData(myCDS)), cell, c("P7", "P5", "RT_Well",
         "Lig_Well", "SampleHolderNum"),
        sep="_"))
colData(myCDS)$P7_ind = splitCellCol$P7
colData(myCDS)$P5_ind = splitCellCol$P5
colData(myCDS)$RT_Well = splitCellCol$RT_Well
colData(myCDS)$Lig_Well = splitCellCol$Lig_Well

colData(myCDS)$sample = as.character(colData(myCDS)$sampleName)


print("CDS Loaded")

# Log10 scale is easier to plot
colData(myCDS)$log_10_umi <- log10(colData(myCDS)$n.umi)


################################################################################################################
# Barnyard plot and species labeling
print("Making species assignments and barnyard plot")
human.genes = rownames(rowData(myCDS))[grepl("^ENSG", rowData(myCDS)$id)]
mouse.genes = rownames(rowData(myCDS))[grepl("^ENSMUSG", rowData(myCDS)$id)]
print(str(human.genes))
colData(myCDS)$n_human_umi = Matrix::colSums(exprs(myCDS)[human.genes,])
colData(myCDS)$n_mouse_umi = Matrix::colSums(exprs(myCDS)[mouse.genes,])
print("Loaded species-specific UMIs")

# ggplot(colData(myCDS))  + 
#   geom_point(aes(x = n_human_umi, y = n_mouse_umi), size = .75, stroke = 0) +
#   scale_color_manual(values = c("TRUE" = "dimgrey", "FALSE" = "black")) +
#   xlim(0,40000) + ylim(0,40000) + 
#   monocle3:::monocle_theme_opts()
# ggsave("./plots/Lung_barnyard.png")

colData(myCDS)$Species <- ifelse(colData(myCDS)$n_human_umi/colData(myCDS)$n_mouse_umi > 10.0, "Human",
                            ifelse(colData(myCDS)$n_mouse_umi/colData(myCDS)$n_human_umi > 10.0, "Mouse", "Mixed"))

colData(myCDS)$sample <- as.character(colData(myCDS)$sample)

boxplot_stat_by_X(myCDS, processingNote, "log_10_umi", "sample")
# # See UMI distributions
# ggplot(as.data.frame(colData(myCDS)), aes(x=sample, y=log_10_umi)) + 
#     geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ggsave("./plots/UMIs_by_Sample.png")
# # ggsave("./plots/Scrambled_UMIs_by_Sample.png")

# boxplot_stat_by_X(myCDS, processingNote, "log_10_umi", "sample")
counts = (as.data.frame(colData(myCDS)) %>% dplyr::group_by(sample) %>%
             dplyr::summarise(n=n()))
png(paste0("./plots/sampleCellCounts", processingNote, ".png"))
ggplot(as.data.frame(counts), aes_string(x="sample", y="n")) +
    geom_bar(position="dodge", stat="identity")+ 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title="Cells per sample")
dev.off()

preUMAPCDS = myCDS
preUMAPNotes = processingNote

# tuckerHeartMarkers = c("RBPJ", "LDB2", "DLC1", "NRXN1", "NEGR1",
#                     "GSN", "FHL2", "TTN")

# Make UMAPs w/ color-coding for each sample, alone
for (eachSample in samplesToRead){
    subsetCDS = myCDS[,colData(myCDS)$sample == eachSample]
    # PCA, UMAP, show some heart markers
    print(paste0("Processing ", eachSample,  " alone"))
    if (nrow(colData(subsetCDS)) > 1){
        set.seed(7)
        subsetCDS <- preprocess_cds(subsetCDS)
        subsetCDS <- reduce_dimension(subsetCDS)
        plotUMAP_Monocle(subsetCDS, paste0(eachSample, "only",processingNote),
             "sample", show_labels=FALSE)
    }
}

# Get some summary plots
medianStats = (as.data.frame(colData(myCDS)) %>%  group_by(sample) 
				%>% summarize(the_medians = median(n.umi)))

# # Pre-process and plot
set.seed(7)
myCDS = preprocess_cds(myCDS) #, num_dim = 50)
print("CDS Pre-processed")
# png(paste0("./plots/", processingNote, "_varianceExplained.png"))
# plot_pc_variance_explained(myCDS)
# dev.off()

tryAlign = TRUE
set.seed(7)
if (tryAlign){
    print("Aligning data by sample")
    alignBy = "sample"
    myCDS = monocle3::align_cds(myCDS, num_dim=50, alignment_group = alignBy)
    processingNote <- paste0(processingNote, "MNN=", alignBy)
} else{
    processingNote <- paste0(processingNote, "_noMNN")
}

set.seed(7)
myCDS = reduce_dimension(myCDS, umap.fast_sgd=FALSE)

# See some plots
plotUMAP_Monocle(myCDS, processingNote, "sample", show_labels=FALSE)

plotUMAP_Monocle_genes(myCDS, processingNote, 
        tuckerHeartMarkers, "Markers_TuckerEtAl")




plotUMAP_Monocle_genes(myCDS, processingNote, 
        c("NEGR1", "ELN", "GSN", "NRXN1", 
            "NRXN3", "CADM2", "ATP1B3", "CRISPLD2"), "_Fibroblast_vs_neuronal_TuckerEtAl")
plotUMAP_Monocle_genes(myCDS, processingNote, 
        c("PROX1", "FLT4", "PDPN", "BMX", "NPR3"), "_endothelial_subtyping_TuckerEtAl")
plotUMAP_Monocle_genes(myCDS, processingNote, 
        c("PRKG1", "RBMS3", "PLA2G5", "GMDS", "RP11-767I20.1"), "_endothelial_2_aka_endocardium_TuckerEtAl")

plotUMAP_Monocle_genes(myCDS, processingNote, 
        c("SYNE1", "PTPRM", "LDB2", "LINC00486"), "_pandEndothelial_TuckerEtAl")

plotUMAP_Monocle_genes(myCDS, processingNote, 
        c("KCNAB1", "CARMN", "ITGA8", "MYH11"), "_vascularSmoothMus_TuckerEtAl")

plotUMAP_Monocle_genes(myCDS, processingNote, 
        c("ACACB", "GPAM", "PDE3B", "EBF1"), "_adipocytes_TuckerEtAl")




# Search for markers
deRes = runDEtestingToID_markers(myCDS, processingNote, "cluster_label")

# Code for hard-code labeling of heart clusters

# Cluster for partition-based coloring
print("Clustering Cells Now")
set.seed(7)
kVal = 40
myCDS = cluster_cells(myCDS, k = kVal)#, resolution=c(10^seq(-6,-1))) #, 
# myCDS = cluster_cells(myCDS), k = 3) #, resolution=c(10^seq(-6,-1)) )#, 
       #python_home="/net/gs/vol3/software/modules-sw-python/2.7.3/louvain")# /381b7db/Linux/RHEL6/x86_64/lib/python2.7/site-packages/python_louvain-0.13-py2.7.egg/")

        #community/")

# Get the partitions data
colData(myCDS)$partitionLabel <- partitions(myCDS) # alternately, clusters(myCDS)
colData(myCDS)$partitionLabel <- as.character(colData(myCDS)$partitionLabel)
colData(myCDS)$partitionGuess <- as.character(colData(myCDS)$partitionLabel)
colData(myCDS)$cluster_label = as.character(clusters(myCDS))
plotUMAP_Monocle(myCDS, paste0("k=", as.character(kVal), "_", processingNote), "cluster_label", show_labels=TRUE)
plotUMAP_Monocle(myCDS, paste0("k=", as.character(kVal), "_", processingNote), "partitionLabel", show_labels=TRUE)

# testTable = with(colData(myCDS), table(sample, partitionLabel))

myCDS <- hardAssignClustLabels(myCDS, kVal)
names(colData(myCDS)$cluster_cell_type) = rownames(colData(myCDS))
myCDS@clusters[["UMAP"]]$clusters <- as.factor(colData(myCDS)$cluster_cell_type)
plotUMAP_Monocle(myCDS, processingNote, "cluster_cell_type", show_labels=TRUE)



# Save an RDS
# saveRDS(myCDS, paste0("./rdsOutput/basicCellTyped_", processingNote, "_cds.RDS"))


# Note anatomical sites
myCDS = readRDS(paste0("./rdsOutput/basicCellTyped_", processingNote, "_cds.RDS"))
myCDS = hardAssignAnatomicalSites(myCDS)
# Plot the distributions

# testSplit = (as.data.frame(colData(myCDS)) %>% 
# 			group_by(sample, Anatomical_Site, cluster_cell_type) 
# 			 )
# str(testSplit)

countTable = with(colData(myCDS), table(sample, cluster_cell_type))
proportionTable <- prop.table(countTable, margin=1)
propDF = as.data.frame(proportionTable)

propDF

propDF <- hardAssignAnatomicalSitesDF(propDF)
# Remove the Unassigned clusters
propDF$cluster_cell_type = as.character(propDF$cluster_cell_type)
propDF = propDF[!(propDF$cluster_cell_type %in% c("20", "23")),]
propDF$cluster_cell_type = ifelse(propDF$cluster_cell_type == "Vascular Smooth Muscle",
				"Vascular Smooth\nMuscle", propDF$cluster_cell_type)
colnames(propDF) = c("sample", "cluster_cell_type", "Freq", "Anatomical_Site")

png(paste0("./plots/", processingNote, "propsByAnatomicalSite.png"))
myPlot <- (ggplot(propDF, 
		aes_string(x="cluster_cell_type", y="Freq", fill="Anatomical_Site",
					color="Anatomical_Site")) +
		geom_boxplot() + 
		ylab("Proportion of Total") + xlab("Cell Type") +
		guide_legend(title="Anatomical Site") +
		#scale_fill_discrete(name = "Anatomical Site")+
		theme(text = element_text(size = 12))    +
		 theme(axis.text.x = element_text(angle = 90, hjust = 1)))
print(myPlot)
dev.off()


plotUMAP_Monocle <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, textSize=14){ #, xVal, yVal){
    png(paste0("./plots/", processingNote, "_UMAP_", catToColor, "colored.png"),
             width=1400, height=1000, res=200)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels,
          cell_stroke=.1 , group_label_size=textSize        
                ))
    myPlot = (myPlot + theme(text=element_text(size=textSize)))
    print(myPlot)
    dev.off()   
}
plotUMAP_Monocle(myCDS, processingNote, "cluster_cell_type", show_labels=TRUE)



































# library(rlang)

# samVar = "sample"
# parVar = "partitionLabel"

# testTableTwo = with(colData(myCDS), table(get(samVar), get(parVar) ))

# testTable = with(colData(myCDS), table(sample, partitionLabel))


# countTable = (as.data.frame(colData(myCDS)) %>% 
#                 dplyr::groupy(c("sample", "partitionLabel")) %>%
#              dplyr::summarise(n=n()))



plotGroupedProportions(myCDS, processingNote, "sample", "partitionLabel")




if (setToShow == "Heart"){
    plotUMAP_Monocle(myCDS, processingNote, "donor", show_labels=FALSE)
    if (tryAlign == FALSE){
        myCDS = cluster_cells(myCDS, k = 6)
        colData(myCDS)$cluster_label = as.character(clusters(myCDS))
        colData(myCDS)$cluster_cell_type = colData(myCDS)$cluster_label


        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type == "1",
                                            "Cardiomyocyte", colData(myCDS)$cluster_cell_type)
        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type == "2",
                                            "Fibroblasts", colData(myCDS)$cluster_cell_type)
        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type == "3",
                                            "Macrophages", colData(myCDS)$cluster_cell_type) #RBPJ
        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type == "4",
                                            "Endothelium", colData(myCDS)$cluster_cell_type) #LDB@
        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type == "5",
                                            "Pericytes", colData(myCDS)$cluster_cell_type) #DLC1
    }
    if ((tryAlign == TRUE) & (alignBy == "donor")){
        myCDS = cluster_cells(myCDS, k = 4)
        colData(myCDS)$cluster_label = as.character(clusters(myCDS))
        colData(myCDS)$cluster_cell_type = colData(myCDS)$cluster_label

        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type %in% c("1", "6"),
                                            "Cardiomyocyte", colData(myCDS)$cluster_cell_type)
        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type %in% c("3"),
                                            "Fibroblast I & II", colData(myCDS)$cluster_cell_type)
        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type %in% c("9", "8", "4"),
                                            "Fibroblast III", colData(myCDS)$cluster_cell_type)
        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type %in% c("12", "10", "2"),
                                            "Macrophages", colData(myCDS)$cluster_cell_type)
        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type == "5",
                                            "Endothelium", colData(myCDS)$cluster_cell_type)
        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type == "7",
                                            "Pericytes", colData(myCDS)$cluster_cell_type)
        colData(myCDS)$cluster_cell_type = ifelse(colData(myCDS)$cluster_cell_type == "11",
                                            "Neuronal", colData(myCDS)$cluster_cell_type)

    }
    names(colData(myCDS)$cluster_cell_type) = rownames(colData(myCDS))
    myCDS@clusters[["UMAP"]]$clusters <- as.factor(colData(myCDS)$cluster_cell_type)
    plotUMAP_Monocle(myCDS, processingNote, "cluster_cell_type", show_labels=TRUE)
}

# Capitalized looks nicer
colData(myCDS)$Donor = colData(myCDS)$donor
myCDS@clusters[["UMAP"]]$clusters <- as.factor(colData(myCDS)$cluster_cell_type)

axisLim = 14
plotUMAP_Monocle_small_panel <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, labelNotByColor=FALSE,
                insteadOfColorLabel="sample", forceLegendShow=FALSE){ #, xVal, yVal){
    png(paste0("./plots/", processingNote, "_UMAP_", catToColor, "colored.png"),
        width=2.5, height=1.75, units="in", res=300)
    print(plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, 
        group_cells_by= "cluster",
        label_cell_groups=show_labels,
          cell_stroke=.1, group_label_size=2,
          labelNotByColor=labelNotByColor,
                insteadOfColorLabel=insteadOfColorLabel,
                legendOverride=forceLegendShow              ) + 
        
        theme(text=element_text(size=6)) + coord_fixed()
        # + xlim(-1*axisLim, axisLim) + ylim(-1*axisLim,axisLim) 
        + xlim(-10, axisLim) + ylim(-10,axisLim) 
        + theme(legend.key.size = unit(.02, "in")) 
         +theme(legend.box.margin=margin(-10,0,-10,-10), legend.margin=margin(0,0,0,0)
            ))
    dev.off()   
}

plotUMAP_Monocle_small_panel(myCDS, "For_IPR", "Donor", labelNotByColor=TRUE,
                insteadOfColorLabel="cluster_cell_type", forceLegendShow=TRUE)

monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}


plot_cells <- function(cds,
                       x=1,
                       y=2,
                       reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
                       color_cells_by="cluster",
                       group_cells_by=c("cluster", "partition"),
                       genes=NULL,
                       show_trajectory_graph=TRUE,
                       trajectory_graph_color="grey28",
                       trajectory_graph_segment_size=0.75,
                       norm_method = c("log", "size_only"),
                       label_cell_groups = TRUE,
                       label_groups_by_cluster=TRUE,

                       labelNotByColor=FALSE,
                       insteadOfColorLabel ="sample",

                       group_label_size=2,
                       labels_per_group=1,
                       label_branch_points=TRUE,
                       label_roots=TRUE,
                       label_leaves=TRUE,
                       graph_label_size=2,
                       cell_size=0.35,
                       cell_stroke= I(cell_size / 2),
                       alpha = 1,
                       min_expr=0.1,
                       rasterize=FALSE,
                       legendOverride=FALSE,
                       scale_to_range=FALSE) {
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
                          msg = paste("No dimensionality reduction for",
                                      reduction_method, "calculated.",
                                      "Please run reduce_dimension with",
                                      "reduction_method =", reduction_method,
                                      "before attempting to plot."))
  low_dim_coords <- reducedDims(cds)[[reduction_method]]
  assertthat::assert_that(ncol(low_dim_coords) >=max(x,y),
                          msg = paste("x and/or y is too large. x and y must",
                                      "be dimensions in reduced dimension",
                                      "space."))
  if(!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster", "partition",
                                                  "pseudotime") |
                              color_cells_by %in% names(colData(cds)),
                            msg = paste("color_cells_by must one of",
                                        "'cluster', 'partition', 'pseudotime,",
                                        "or a column in the colData table."))

    if(color_cells_by == "pseudotime") {
      tryCatch({pseudotime(cds, reduction_method = reduction_method)},
               error = function(x) {
                 stop(paste("No pseudotime for", reduction_method,
                            "calculated. Please run order_cells with",
                            "reduction_method =", reduction_method,
                            "before attempting to color by pseudotime."))})

    }
  }
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
                          msg = paste("Either color_cells_by or markers must",
                                      "be NULL, cannot color by both!"))

  norm_method = match.arg(norm_method)
  group_cells_by=match.arg(group_cells_by)
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
                          msg = paste("Either color_cells_by or genes must be",
                                      "NULL, cannot color by both!"))

  if (show_trajectory_graph &&
      is.null(principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }

  gene_short_name <- NA
  sample_name <- NA
  #sample_state <- colData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA
  if (rasterize){
    plotting_func <- ggrastr::geom_point_rast
  }else{
    plotting_func <- ggplot2::geom_point
  }

  S_matrix <- reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[,c(x,y)])

  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)

  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  if (group_cells_by == "cluster"){
    data_df$cell_group <-
      tryCatch({clusters(cds,
                         reduction_method = reduction_method)[
                           data_df$sample_name]},
               error = function(e) {NULL})
  } else if (group_cells_by == "partition") {
    data_df$cell_group <-
      tryCatch({partitions(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]},
               error = function(e) {NULL})
  } else{
    stop("Error: unrecognized way of grouping cells.")
  }

  if (color_cells_by == "cluster"){
    data_df$cell_color <-
      tryCatch({clusters(cds,
                         reduction_method = reduction_method)[
                           data_df$sample_name]},
               error = function(e) {NULL})
  } else if (color_cells_by == "partition") {
    data_df$cell_color <-
      tryCatch({partitions(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]},
               error = function(e) {NULL})
  } else if (color_cells_by == "pseudotime") {
    data_df$cell_color <-
      tryCatch({pseudotime(cds,
                           reduction_method = reduction_method)[
                             data_df$sample_name]}, error = function(e) {NULL})
  } else{
    data_df$cell_color <- colData(cds)[data_df$sample_name,color_cells_by]
  }

  ## Graph info
  if (show_trajectory_graph) {

    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
      as.data.frame() %>%
      dplyr::select_(prin_graph_dim_1 = x, prin_graph_dim_2 = y) %>%
      dplyr::mutate(sample_name = rownames(.),
                    sample_state = rownames(.))

    dp_mst <- cds@principal_graph[[reduction_method]]

    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select_(source = "from", target = "to") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(
                           source="sample_name",
                           source_prin_graph_dim_1="prin_graph_dim_1",
                           source_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "source") %>%
      dplyr::left_join(ica_space_df %>%
                         dplyr::select_(
                           target="sample_name",
                           target_prin_graph_dim_1="prin_graph_dim_1",
                           target_prin_graph_dim_2="prin_graph_dim_2"),
                       by = "target")
  }

  ## Marker genes
  markers_exprs <- NULL
  expression_legend_label <- NULL
  if (!is.null(genes)) {
    if (!is.null(dim(genes)) && dim(genes) >= 2){
      markers = unlist(genes[,1], use.names=FALSE)
    } else {
      markers = genes
    }
    markers_rowData <- rowData(cds)[(rowData(cds)$gene_short_name %in% markers) |
                                    (row.names(rowData(cds)) %in% markers),,drop=FALSE]
    markers_rowData <- as.data.frame(markers_rowData)
    if (nrow(markers_rowData) == 0) {
      stop("None of the provided genes were found in the cds")
    }
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), ,drop=FALSE]
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))

      if (!is.null(dim(genes)) && dim(genes) >= 2){
        genes = as.data.frame(genes)
        row.names(genes) = genes[,1]
        genes = genes[row.names(cds_exprs),]

        agg_mat = as.matrix(aggregate_gene_expression(cds, genes, norm_method=norm_method, scale_agg_values=FALSE))
        markers_exprs = agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        if (is.factor(genes[,2]))
          markers_exprs$feature_id = factor(markers_exprs$feature_id,
                                            levels=levels(genes[,2]))

        markers_exprs$feature_label <- markers_exprs$feature_id
        norm_method = "size_only"
        expression_legend_label = "Expression score"
      } else {
        cds_exprs@x = round(10000*cds_exprs@x)/10000
        markers_exprs = matrix(cds_exprs, nrow=nrow(markers_rowData))
        colnames(markers_exprs) = colnames(SingleCellExperiment::counts(cds))
        row.names(markers_exprs) = row.names(markers_rowData)
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c('feature_id','cell_id')
        markers_exprs <- merge(markers_exprs, markers_rowData,
                               by.x = "feature_id", by.y="row.names")
        if (is.null(markers_exprs$gene_short_name)) {
          markers_exprs$feature_label <-
            as.character(markers_exprs$feature_id)
        } else {
          markers_exprs$feature_label <-
            as.character(markers_exprs$gene_short_name)
        }

        markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) | !as.character(markers_exprs$feature_label) %in% markers,
                                              as.character(markers_exprs$feature_id),
                                              as.character(markers_exprs$feature_label))

        markers_exprs$feature_label <- factor(markers_exprs$feature_label,
                                              levels = markers)
        if (norm_method == "size_only")
          expression_legend_label = "Expression"
        else
          expression_legend_label = "log10(Expression)"
      }

      if (scale_to_range){
        markers_exprs = dplyr::group_by(markers_exprs, feature_label) %>%
          dplyr::mutate(max_val_for_feature = max(value),
                        min_val_for_feature = min(value)) %>%
          dplyr::mutate(value = 100 * (value - min_val_for_feature) / (max_val_for_feature - min_val_for_feature))
        expression_legend_label = "% Max"
      }
    }
  }

  if (label_cell_groups && is.null(color_cells_by) == FALSE){
    if (is.null(data_df$cell_color)){
      if (is.null(genes)){
        message(paste(color_cells_by, "not found in colData(cds), cells will",
                      "not be colored"))
      }
      text_df = NULL
      label_cell_groups = FALSE
    }else{
      if(is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
       #browser()

        # DFR Check: Override the "color" value with something else to define
        #   labels?
        if (labelNotByColor){
            backupColor = data_df$cell_color
            data_df$cell_color = data_df[[insteadOfColorLabel]]
        }
        #####################################

        if (label_groups_by_cluster && is.null(data_df$cell_group) == FALSE){
          text_df = data_df %>%
            dplyr::group_by(cell_group) %>%
            dplyr::mutate(cells_in_cluster= dplyr::n()) %>%
            dplyr::group_by(cell_color, add=TRUE) %>%
            dplyr::mutate(per=dplyr::n()/cells_in_cluster)
          median_coord_df = text_df %>%
            dplyr::summarize(fraction_of_group = dplyr::n(),
                             text_x = stats::median(x = data_dim_1),
                             text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df,
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_group) %>%
            dplyr::top_n(labels_per_group, per)
        } else {
          text_df = data_df %>% dplyr::group_by(cell_color) %>%
            dplyr::mutate(per=1)
          median_coord_df = text_df %>%
            dplyr::summarize(fraction_of_group = dplyr::n(),
                             text_x = stats::median(x = data_dim_1),
                             text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>%
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df,
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_color) %>%
            dplyr::top_n(labels_per_group, per)
        }

        text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
        # I feel like there's probably a good reason for the bit below, but I
        # hate it and I'm killing it for now.
        # text_df$label <- paste0(1:nrow(text_df))
        # text_df$process_label <- paste0(1:nrow(text_df), '_',
        # as.character(as.matrix(text_df[, 1])))
        # process_label <- text_df$process_label
        # names(process_label) <- as.character(as.matrix(text_df[, 1]))
        # data_df[, group_by] <-
        #  process_label[as.character(data_df[, group_by])]
        # text_df$label = process_label
      } else {
        message(paste("Cells aren't colored in a way that allows them to",
                      "be grouped."))
        text_df = NULL
        label_cell_groups = FALSE
      }
    }
  }

  # David: Add this to reset color
  data_df$cell_color = backupColor
  ###############################

  if (!is.null(markers_exprs) && nrow(markers_exprs) > 0){
    data_df <- merge(data_df, markers_exprs, by.x="sample_name",
                     by.y="cell_id")
    data_df$value <- with(data_df, ifelse(value >= min_expr, value, NA))
    ya_sub <- data_df[!is.na(data_df$value),]
    na_sub <- data_df[is.na(data_df$value),]
    if(norm_method == "size_only"){
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
        plotting_func(aes(data_dim_1, data_dim_2), size=I(cell_size),
                      stroke = I(cell_stroke), color = "grey80", alpha = alpha,
                      data = na_sub) +
        plotting_func(aes(color=value), size=I(cell_size),
                      stroke = I(cell_stroke),
                      data = ya_sub[order(ya_sub$value),]) +
        viridis::scale_color_viridis(option = "viridis",
                                     name = expression_legend_label,
                                     na.value = NA, end = 0.8,
                                     alpha = alpha) +
        guides(alpha = FALSE) + facet_wrap(~feature_label)
    } else {
      g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) +
        plotting_func(aes(data_dim_1, data_dim_2), size=I(cell_size),
                      stroke = I(cell_stroke), color = "grey80",
                      data = na_sub, alpha = alpha) +
        plotting_func(aes(color=log10(value+min_expr)),
                      size=I(cell_size), stroke = I(cell_stroke),
                      data = ya_sub[order(ya_sub$value),],
                      alpha = alpha) +
        viridis::scale_color_viridis(option = "viridis",
                                     name = expression_legend_label,
                                     na.value = NA, end = 0.8,
                                     alpha = alpha) +
        guides(alpha = FALSE) + facet_wrap(~feature_label)
    }
  } else {
    g <- ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2))

    # We don't want to force users to call order_cells before even being able
    # to look at the trajectory, so check whether it's null and if so, just
    # don't color the cells
    if(color_cells_by %in% c("cluster", "partition")){
      if (is.null(data_df$cell_color)){
        g <- g + geom_point(color=I("gray"), size=I(cell_size),
                            stroke = I(cell_stroke), na.rm = TRUE,
                            alpha = I(alpha))
        message(paste("cluster_cells() has not been called yet, can't",
                      "color cells by cluster"))
      } else{
        g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                            stroke = I(cell_stroke), na.rm = TRUE,
                            alpha = alpha)
      }
      g <- g + guides(color = guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    } else if (class(data_df$cell_color) == "numeric"){
      g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + viridis::scale_color_viridis(name = color_cells_by, option="C")
    } else {
      g <- g + geom_point(aes(color = cell_color), size=I(cell_size),
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + guides(color = guide_legend(title = color_cells_by,
                                           override.aes = list(size = 4)))
    }

  }
  if (show_trajectory_graph){
    g <- g + geom_segment(aes_string(x="source_prin_graph_dim_1",
                                     y="source_prin_graph_dim_2",
                                     xend="target_prin_graph_dim_1",
                                     yend="target_prin_graph_dim_2"),
                          size=trajectory_graph_segment_size,
                          color=I(trajectory_graph_color),
                          linetype="solid",
                          na.rm=TRUE,
                          data=edge_df)


    if (label_branch_points){
      mst_branch_nodes <- branch_nodes(cds, reduction_method)
      branch_point_df <- ica_space_df %>%
        dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>%
        dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))

      g <- g +
        geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                   shape = 21, stroke=I(trajectory_graph_segment_size),
                   color="white",
                   fill="black",
                   size=I(graph_label_size * 1.5),
                   na.rm=TRUE, branch_point_df) +
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                             label="branch_point_idx"),
                  size=I(graph_label_size), color="white", na.rm=TRUE,
                  branch_point_df)
    }

    if (label_leaves){
      mst_leaf_nodes <- leaf_nodes(cds, reduction_method)
      leaf_df <- ica_space_df %>%
        dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>%
        dplyr::mutate(leaf_idx = seq_len(dplyr::n()))

      g <- g +
        geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                   shape = 21, stroke=I(trajectory_graph_segment_size),
                   color="black",
                   fill="lightgray",
                   size=I(graph_label_size * 1.5),
                   na.rm=TRUE,
                   leaf_df) +
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                             label="leaf_idx"),
                  size=I(graph_label_size), color="black", na.rm=TRUE, leaf_df)
    }

    if (label_roots){
      mst_root_nodes <- root_nodes(cds, reduction_method)
      root_df <- ica_space_df %>%
        dplyr::slice(match(names(mst_root_nodes), sample_name)) %>%
        dplyr::mutate(root_idx = seq_len(dplyr::n()))

      g <- g +
        geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                   shape = 21, stroke=I(trajectory_graph_segment_size),
                   color="black",
                   fill="white",
                   size=I(graph_label_size * 1.5),
                   na.rm=TRUE,
                   root_df) +
        geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2",
                             label="root_idx"),
                  size=I(graph_label_size), color="black", na.rm=TRUE, root_df)
    }
  }

  if(label_cell_groups) {
    g <- g + ggrepel::geom_text_repel(data = text_df,
                                      mapping = aes_string(x = "text_x",
                                                           y = "text_y",
                                                           label = "label"),
                                      size=I(group_label_size))
    # If we're coloring by gene expression, don't hide the legend
    if (is.null(markers_exprs) & !legendOverride)
      g <- g + theme(legend.position="none")
  }

  g <- g +
    #scale_color_brewer(palette="Set1") +
    monocle_theme_opts() +
    xlab(paste(reduction_method, x)) +
    ylab(paste(reduction_method, y)) +
    #guides(color = guide_legend(label.position = "top")) +
    theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill='white'))
  g
}


















































# From HM6, for reference
# axisLim = 11
# plotUMAP_Monocle <- function(dataCDS, processingNote, catToColor,
#                     show_labels=TRUE){ #, xVal, yVal){
#     png(paste0("./plots/", processingNote, "_UMAP_", catToColor, "colored.png"),
#         width=2.5, height=1.75, units="in", res=300)
#     print(plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
#         color_cells_by=catToColor, label_cell_groups=show_labels,
#           cell_stroke=.1, group_label_size=4              ) + 
#         theme(text=element_text(size=6)) + coord_fixed()
#         + xlim(-1*axisLim, axisLim) + ylim(-1*axisLim,axisLim) 
#         + theme(legend.key.size = unit(.02, "in")) 
#          +theme(legend.box.margin=margin(-10,0,-10,-10), legend.margin=margin(0,0,0,0)
#             ))
#     dev.off()   
# }








png(paste0("./plots/new_LungUMAP_colorBySample", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="sample", label_cell_groups=FALSE)
plot_cells(myCDS, color_cells_by="sample_split_by_hash", label_cell_groups=FALSE)
dev.off()

# By RT
png(paste0("./plots/LungUMAP_With_HEK_colorByRT", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="sample", label_cell_groups=FALSE)
plot_cells(myCDS, color_cells_by="RT_Sample", label_cell_groups=FALSE)
dev.off()

# png("./plots/Lung_upper_airway.png", width=1400, height=1000, res=200)
png(paste0("./plots/Lung_upper_airway", processingNotes, ".png"), width=1000, height=1000, res=200)
plot_cells(myCDS, genes=c("Cdhr3", "Tff2", "Muc5ac", "Muc5b", "Lipf", "Krt13", "Krt4"),  label_cell_groups=FALSE)
dev.off()

# png("./plots/Lung_Epithelial_Cells.png", width=1400, height=1000, res=200)
png(paste0("./plots/Lung_Epithelial_Cells", processingNotes, ".png"), width=1000, height=1000, res=200)
plot_cells(myCDS, genes=c("Ager", "Hopx", "Sftpc", "Slc34a2", "Abca3"),  label_cell_groups=FALSE)
dev.off()

# png("./plots/Lung_Macrophages.png", width=1400, height=1000, res=200)
png(paste0("./plots/Lung_Macrophages", processingNotes, ".png"), width=1000, height=1000, res=200)
plot_cells(myCDS, genes=c("Cd68", "Fcgr1", "Tgfb1", "Trem1", "Csf2ra", "Adgre1", "Siglecf", "Lrp1"),  label_cell_groups=FALSE)
dev.off()

# From jonathan markers
# png("./plots/Lung_Type_1_Pneumocyte.png", width=1400, height=1000, res=200)
# plot_cells(myCDS, genes=c("Ager", "Cav1", "Emp2","Aqp5"),  label_cell_groups=FALSE)
# dev.off()

# png("./plots/Lung_Type_2_Pneumocyte.png", width=1400, height=1000, res=200)
# plot_cells(myCDS, genes=c("Sftpa1", "Sftpa2", "Sftpc", "Sftpb"),  label_cell_groups=FALSE)
# dev.off()

# png("./plots/Lung_Ciliated_Epith.png", width=1400, height=1000, res=200)
# plot_cells(myCDS, genes=c("Foxj1", "Rsph1", "Sntn", "Tppp3"),  label_cell_groups=FALSE)
# dev.off()

# png("./plots/Lung_Basal_Epith.png", width=1400, height=1000, res=200)
# plot_cells(myCDS, genes=c("Krt5", "Tp63", "Ngfr"),  label_cell_groups=FALSE)
# dev.off()

# png("./plots/Lung_Club_Epith.png", width=1400, height=1000, res=200)
# plot_cells(myCDS, genes=c("Scgb1a1", "Scgb3a1", "Muc5ac", "Muc5b"),  label_cell_groups=FALSE)
# dev.off()

# png("./plots/Lung_Fibroblast.png", width=1400, height=1000, res=200)
# plot_cells(myCDS, genes=c("Pdgfra", "Pdgfrb", "Thy1", "Lum", "Col1a1"),  label_cell_groups=FALSE)
# dev.off()

png("./plots/Lung_Jonathan_Mac.png", width=1400, height=1000, res=200)
plot_cells(myCDS, genes=c("Cd68", "Apoc1", "Ftl", "C1qa", "Cd74", "Mcemp1"),  label_cell_groups=FALSE)
dev.off()

# Olfactory markers
png("./plots/Olfactory_Epithelium.png", width=1400, height=1000, res=200)
plot_cells(myCDS, genes=c("MUC5AC","Muc8", "TP63", "TUBB"),  label_cell_groups=FALSE)
dev.off()

png("./plots/Lung_Human_Fibroblast.png", width=1400, height=1000, res=200)
plot_cells(myCDS, genes=c("PDGFRA", "PDGFRB", "THY1", "LUM", "COL1A1"),  label_cell_groups=FALSE)
dev.off()

png("./plots/Lung_Human_Ciliated_Epith.png", width=1400, height=1000, res=200)
plot_cells(myCDS, genes=c("FOXJ1", "RSPH1", "SNTN", "TPPP3"),  label_cell_groups=FALSE)
dev.off()

png("./plots/Lung_Human_upper_airway.png", width=1400, height=1000, res=200)
plot_cells(myCDS, genes=c("CDHR3", "TFF2", "MUC5AC", "MUC5B", "LIPF", "KRT13", "KRT4"),  label_cell_groups=FALSE)
dev.off()

png("./plots/Lung_Human_Selected_upper_airway.png", width=1400, height=1000, res=200)
plot_cells(myCDS, genes=c( "MUC5B", "KRT13", "KRT4"),  label_cell_groups=FALSE)
dev.off()

png("./plots/Lung_Human_Basal_Epith.png", width=1400, height=1000, res=200)
plot_cells(myCDS, genes=c("KRT5", "TP63", "NGFR" ),  label_cell_groups=FALSE)
dev.off()

# png("./plots/Debley_Markers_Cell_Types.png", width=1400, height=1000, res=200)
png(paste0("./plots/Debley_Markers_Cell_Types", processingNotes, ".png"), width=1000, height=1000, res=200)
plot_cells(myCDS, genes=c("TRPM5", "POU2F3", "ALOX5AP", "IL25", "FOXI1", "CFTR", "MUC5AC", "CDHR3", "NFIA"),  label_cell_groups=FALSE)
dev.off()

# png("./plots/Debley_Markers_General.png", width=1400, height=1000, res=200)
png(paste0("./plots/Debley_Markers_General", processingNotes, ".png"), width=1000, height=1000, res=200)
plot_cells(myCDS, genes=c("TGFB2", "TGFB1", "FSTL3", "INHBA", "TSLP", "IL25", "IFIH1"),  label_cell_groups=FALSE)
dev.off()

# png("./plots/TGFB_Polarity.png", width=1400, height=1000, res=200)
png(paste0("./plots/TGFB_Polarity", processingNotes, ".png"), width=1000, height=1000, res=200)
plot_cells(myCDS, genes=c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3"),  label_cell_groups=FALSE)
dev.off()

# png("./plots/EMT_Signature.png", width=1400, height=1000, res=200)
png(paste0("./plots/EMT_Signature", processingNotes, ".png"), width=1000, height=1000, res=200)
plot_cells(myCDS, genes=c("FN1", "VIM", "CDH1", "TJP1", "ACTA2", "COL1A1"),  label_cell_groups=FALSE)
dev.off()


png("./plots/TH2_Signature.png", width=1400, height=1000, res=200)
plot_cells(myCDS, genes=c("POSTN", "SERPINB2", "CLCA1", "TJP1", "CLDN18", "CLDN1", "CLDN7", "CLDN9"),  label_cell_groups=FALSE)
dev.off()

# png("./plots/Debley_Type_Markers.png", width=1400, height=1000, res=200)
png(paste0("./plots/Debley_Type_Markers", processingNotes, ".png"), width=1000, height=1000, res=200)
plot_cells(myCDS, genes=c("TP63", "KRT5", "KRT14", "KRT8", "TUBB4A","FOXJ1", "MUC5AC",  "MUC5B",  "IL25",
                        "TSLP"),  label_cell_groups=FALSE)
dev.off()

png("./plots/Cell_Cycle_Human.png", width=1400, height=1000, res=200)
plot_cells(myCDS, genes=c( "CCNE1", "CCNA1", "CCNB1"),  label_cell_groups=FALSE)
dev.off()


# Plot cells by patient-of-origin, showing some key markers that Jason is interested in
for (eachPatient in c("T178", "T175", "T390")){
    patientCDS = myCDS[,colData(myCDS)$patient_origin == eachPatient]
    str(colData(patientCDS))
    print(paste0("./plots/", eachPatient, "_R01_Genes", processingNotes, ".png"))

    png(paste0("./plots/", eachPatient, "_R01_Genes", processingNotes, ".png"), width=1000, height=1000, res=400)
    plot_cells(patientCDS, genes=c("TGFB1", "TGFB2", "TGFBR2", "FN1", "VIM", "CDH1", "ACTA2", "COL1A1", "TJP1", "POSTN"),
                  label_cell_groups=FALSE)
    dev.off()
}


# png("./plots/Lung_Human_Ciliated_Epith.png", width=1400, height=1000, res=200)
# plot_cells(myCDS, genes=c("Foxj1", "Rsph1", "Sntn", "Tppp3"),  label_cell_groups=FALSE)
# dev.off()

# Cluster for partition-based coloring
print("Clustering Cells Now")
myCDS = cluster_cells(myCDS, k = 5)#, resolution=c(10^seq(-6,-1))) #, 
# myCDS = cluster_cells(myCDS), k = 3) #, resolution=c(10^seq(-6,-1)) )#, 
       #python_home="/net/gs/vol3/software/modules-sw-python/2.7.3/louvain")# /381b7db/Linux/RHEL6/x86_64/lib/python2.7/site-packages/python_louvain-0.13-py2.7.egg/")

        #community/")

# Get the partitions data
colData(myCDS)$partitionLabel <- partitions(myCDS) # alternately, clusters(myCDS)
colData(myCDS)$partitionLabel <- as.character(colData(myCDS)$partitionLabel)
colData(myCDS)$partitionGuess <- as.character(colData(myCDS)$partitionLabel)

# colData(myCDS)$partitionGuess[which(colData(myCDS)$partitionGuess %in% c("1", "17","19","24","12","13"))] = "Ciliated"
# # colData(myCDS)$partitionGuess[which(colData(myCDS)$partitionGuess %in% c("2", "9","21","5","21","10","3","23"))] = "Basal"
# colData(myCDS)$partitionGuess[which(colData(myCDS)$partitionGuess %in% c("9","21","5","21","10","3","23"))] = "Basal"
# colData(myCDS)$partitionGuess[which(colData(myCDS)$partitionGuess %in% c("6", "8", "11", "22"))] = "Goblet" # Maybe 2?
# colData(myCDS)$partitionGuess[which(colData(myCDS)$partitionGuess %in% c("4", "20","7"))] = "Basal_to_Ciliated_Trans?"
# colData(myCDS)$partitionGuess[which(colData(myCDS)$partitionGuess %in% c("15", "14", "16", "18"))] = "Unkown"


# png("./plots/Lung_partitionPlot.png", width=1400, height=1000, res=200 )
png(paste0("./plots/Lung_partitionPlot", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="partition", group_cells_by="partition")
plot_cells(myCDS, color_cells_by="partitionLabel", group_cells_by="partition")
dev.off()

# Learn a graph within this
myCDS <- learn_graph(myCDS, use_partition = FALSE)
png(paste0("./plots/graph_LungUMAP_colorBySample", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="sample", label_cell_groups=FALSE)
plot_cells(myCDS, color_cells_by="sample_split_by_hash", label_cell_groups=FALSE, label_leaves=FALSE,
                label_branch_points=FALSE)
dev.off()

png(paste0("./plots/graph_labeled_branches_LungUMAP_colorBySample", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="sample", label_cell_groups=FALSE)
plot_cells(myCDS, color_cells_by="sample_split_by_hash", label_cell_groups=FALSE, label_leaves=TRUE,
                label_branch_points=TRUE)
dev.off()

# Order. 
# Looks like leaf 7 is the start I want to try
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds){  #, time_bin="130-170"){
  #cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  cell_ids <- which(colData(cds)[, "sample"] == "ALI_5")
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}



# Trying to figure a way to pick based on designating a node
################################################################################################
myCDS <- order_cells(myCDS, root_pr_nodes=c(names(leaf_nodes(myCDS)[7])[1] ) )#, root_pr_nodes=get_earliest_principal_node(myCDS))
# myCDS <- order_cells(myCDS, root_pr_nodes=c("Y_165") )#, root_pr_nodes=get_earliest_principal_node(myCDS))


png(paste0("./plots/graph_ordered_LungUMAP_colorBySample", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="sample", label_cell_groups=FALSE)
plot_cells(myCDS, color_cells_by="pseudotime", label_cell_groups=FALSE, label_leaves=FALSE,
                label_branch_points=FALSE)
dev.off()
################################################################################################

# Alternately, I can try using Hyeon Jin's method of picking cells based on a UMAP coordinate
root_cells_filter = myCDS@reducedDims$UMAP[,1] > 3
root_cells <- rownames(as.data.frame(myCDS@reducedDims$UMAP[root_cells_filter, ]))
myCDS = order_cells(myCDS, root_cells = root_cells)
myCDS = order_cells(myCDS, 
                root_pr_nodes = unique(myCDS@principal_graph_aux[["UMAP"]]$root_pr_nodes))

png(paste0("./plots/pseudotime_HJ_rootmethod_LungUMAP", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="sample", label_cell_groups=FALSE)
plot_cells(myCDS, color_cells_by="pseudotime", label_cell_groups=FALSE, label_leaves=FALSE,
                label_branch_points=FALSE)
dev.off()

################################################################################

# > leaf_nodes(myCDS)
#  Y_10  Y_46  Y_58  Y_78 Y_141 Y_143 Y_165 Y_196 
#    10    46    58    78   141   143   165   196 

myCDS <- order_cells(myCDS, root_pr_nodes=c("Y_58"))#, root_pr_nodes=get_earliest_principal_node(myCDS))

png(paste0("./plots/diff_start_graph_LungUMAP_colorBySample", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="sample", label_cell_groups=FALSE)
plot_cells(myCDS, color_cells_by="pseudotime", label_cell_groups=FALSE, label_leaves=FALSE,
                label_branch_points=FALSE)
dev.off()



# Marker testing
###############################################################################################################

# Graph Test

graph_test_res <- graph_test(myCDS, neighbor_graph="principal_graph", cores=4)
underQ_05 <- row.names(subset(graph_test_res, q_value < 0.05))





# Load Garnett
# library("~/bin/garnett")
# (base) [readdf@t003 bin]$ git clone --single-branch --branch monocle3 https://github.com/cole-trapnell-lab/garnett.git
load_all("~/bin/garnett")
# load_all("~/bin/org.Hs.eg.db")
library(org.Hs.eg.db)

# myDB = library(org.Hs.eg.db)

# marker_file_path <- system.file("extdata",
#          "/net/trapnell/vol1/home/readdf/trapLabDir/fibrosis/data/garnettModels/ALI_Custom_Markers.txt",
#          package="garnett")

# require(org.Hs.eg.db)

# load_all("~/bin/garnett")
# load_all("~/bin/org.Hs.eg.db")
marker_file_path <- "/net/trapnell/vol1/home/readdf/trapLabDir/fibrosis/data/garnettModels/ALI_Custom_Markers.txt"
marker_check <- check_markers(myCDS, marker_file_path,
                              db=org.Hs.eg.db,
                              # cds_gene_id_type = "SYMBOL",
                              cds_gene_id_type = "ENSEMBL",

                              marker_file_gene_id_type = "SYMBOL")

png(paste0("./plots/garnett_marker_check", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="sample", label_cell_groups=FALSE)
plot_markers(marker_check)
dev.off()

# Now get a classifier
set.seed(7)
airwayClassifier <- train_cell_classifier(cds = myCDS,
                                        marker_file = marker_file_path,
                                        db=org.Hs.eg.db,
                                        cds_gene_id_type="ENSEMBL",
                                        marker_file_gene_id_type="SYMBOL")

# Save the classifier:
classifierFileOutputPath = paste0("/net/trapnell/vol1/home/readdf/trapLabDir/fibrosis/data/garnettModels/",
            processingNotes, "ALI_5_model.rds")


saveRDS(airwayClassifier, classifierFileOutputPath)

# Load classifier
airwayClassifier <- readRDS(classifierFileOutputPath)

# Classify
myCDS <- classify_cells(myCDS, airwayClassifier,
                        db = org.Hs.eg.db,
                        cluster_extend=TRUE,
                        cds_gene_id_type="ENSEMBL")

# Plot

png(paste0("./plots/LungUMAP_with_Garnett_types", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="sample", label_cell_groups=FALSE)
plot_cells(myCDS, color_cells_by="cluster_ext_type", label_cell_groups=FALSE, label_leaves=FALSE,
                label_branch_points=FALSE)
dev.off()


library(dplyr)

counts = (as.data.frame(colData(myCDS)) %>% dplyr::group_by(sample_split_by_hash, cluster_ext_type) %>%
             tally())
clustProps = counts %>% dplyr::mutate(prop=prop.table(n))
clustProps

# str(clustProps %>% dplyr::ungroup() %>% dplyr::complete(sample_split_by_hash, cluster_ext_type, 
#                 fill = list(N=0, freq=0)) )

write.csv(as.data.frame(clustProps), 
    paste0("./plots/ALI5_", processingNotes, "garnett_proportions.csv"))


# Get some statistics
clustPropDF = as.data.frame(clustProps)
clustPropDF$sample_split_by_hash = as.character(clustPropDF$sample_split_by_hash)

cellTypes = c("Basal", "Ciliated", "Goblet")

# clustPropDF_identified[("T175" %in% clustPropDF_identified$sample_split_by_hash),]

clustPropDF$patient_origin = with(clustPropDF, ifelse( 
                                (clustPropDF$sample_split_by_hash %in% c("T178", "Hash_ALI_T178", "Trip_T178")), "T178",
                                ifelse( (clustPropDF$sample_split_by_hash %in% c("T175", "Hash_ALI_T175", "Trip_T175")), "T175",
                                ifelse(clustPropDF$sample_split_by_hash %in% c("T390", "Hash_ALI_T390", "Trip_390"), "T390", "Undefined")
                                    )            ) )

# Get the subset w/ a known sample (not unassigned hash sample)
clustPropDF_identified = clustPropDF[!(clustPropDF$patient_origin == "Undefined"),]

for (eachType in cellTypes){
    print(paste0("Type ", eachType))
    print(clustPropDF_identified[clustPropDF_identified$cluster_ext_type == eachType,])

}
















# Save the file

# # Chop out the sub-cluster that seems like doublets
# topHalfCDS <- 

# DE for each subset
# First, label each cell with it's patient of origin
colData(myCDS)$patient_origin = with(colData(myCDS), ifelse( 
                                (colData(myCDS)$sample_split_by_hash %in% c("T178", "Hash_ALI_T178")), "T178",
                                ifelse( (colData(myCDS)$sample_split_by_hash %in% c("T175", "Hash_ALI_T175")), "T175",
                                ifelse(colData(myCDS)$sample_split_by_hash %in% c("T390", "Hash_ALI_T390"), "T390", "Undefined")
                                    )            ) )


# Get violin plots
genesToViolin = c("ABCA13", "MAPK10", "FGFR3", "PTPRK")
mySubset = myCDS[rowData(myCDS)$gene_short_name %in% genesToViolin,]
png(paste0("./plots/Violin_Plot_Demo_Genes", processingNotes, ".png"), width=1000, height=1000, res=200)
#plot_cells(myCDS, color_cells_by="sample", label_cell_groups=FALSE)
plot_genes_violin(mySubset, group_cells_by="patient_origin", ncol=2)
dev.off()


# Get for each subset of cells
library(dplyr)
# thisType <- "Basal"
# thisTypeCDS <- myCDS[,myCDS$cluster_ext_type == thisType]
# # DE test
# gene_patient_fits <- fit_models(thisTypeCDS, model_formula_str = "~patient_origin")
# fit_coefs <- coefficient_table(gene_patient_fits)
# coefsAsDF <- as.data.frame(fit_coefs[,c("id","gene_short_name","term","estimate","q_value","p_value")])

# sortedCoefs <- coefsAsDF[order(coefsAsDF$q_value),] %>% filter(term %in% c("patient_originT178", "patient_originT175",
#                                                                            "patient_originT390")) %>%
#                     filter(q_value < 0.05)
# write.csv(sortedCoefs, paste0("./extraOutputFiles/", thisType, "_under_q_05.csv") )
# sortedCoefs[1:100,]

 # cellTypeFits_by_patient_underQ05 <- coefsAsDF %>% filter("patient_origin" %in% term) %>% filter(q_value < 0.05) %>%
#                                 dplyr::select(gene_short_name, term, q_value, estimate)


# # Goblet fits
# gobletCDS <- myCDS[,myCDS$cluster_ext_type == "Goblet"]
# #
# goblet_fits <- fit_models(gobletCDS, model_formula_str = "~patient_origin")
# gobletCoefs <- coefficient_table()

for (cellType in c("Goblet", "Ciliated", "Basal")){
    thisTypeCDS <- myCDS[,myCDS$cluster_ext_type == cellType]

    for (eachOrigin in levels(as.factor(colData(thisTypeCDS)$patient_origin) )){
        thisOrigin = as.character(eachOrigin)
        colData(thisTypeCDS)[eachOrigin] = with(colData(thisTypeCDS), ifelse((colData(thisTypeCDS)$patient_origin == thisOrigin), 1, 0) )
    }

    # DE test
    gene_patient_fits <- fit_models(thisTypeCDS, model_formula_str = "~patient_origin")
    fit_coefs <- coefficient_table(gene_patient_fits)
    # Save the results
    testResFile = paste0("./extraOutputFiles/", cellType, "_fit_by_patient.rds")
    saveRDS(gene_patient_fits, testResFile)
    # Also get a subset that shows only patient coefficients with q_values < 0.05
    coefsAsDF <- as.data.frame(fit_coefs[,c("id","gene_short_name","term","estimate","q_value","p_value")])

    sortedCoefs <- coefsAsDF[order(coefsAsDF$q_value),] %>% filter(term != "(Intercept)") %>%
                    filter(q_value < 0.05)

    # sortedCoefs <- coefsAsDF[order(coefsAsDF$q_value),] %>% filter(term %in% c("patient_originT178", "patient_originT175",
    #                                                                         "patient_originT390")) %>%
    #                 filter(q_value < 0.05)
    # Save
    write.csv(sortedCoefs, paste0("./extraOutputFiles/", cellType, "_under_q_05.csv") )
}

# colData(myCDS)$sample_split_by_hash = with(colData(myCDS), ifelse( 
#                                 (colData(myCDS)$RT_Sample=="Hash_ALI" &  colData(myCDS)$qval < .05),
#                                         paste(colData(myCDS)$RT_Sample, colData(myCDS)$top_oligo, sep="_"),
#                                         colData(myCDS)$RT_Sample
#                                          ) )





usableTerms = c("patient_originT178", "patient_originT175",   "patient_originT390")
weakSortedCoefs <- coefsAsDF[order(coefsAsDF$q_value),] %>% filter(mapply('%in%', term, usableTerms)) %>%
                    filter(q_value < 0.15)


weakSortedCoefs <- coefsAsDF[order(coefsAsDF$q_value),] %>% filter(term %in% c("patient_originT178", "patient_originT175",
                                                                            "patient_originT390")) %>%
                    filter(q_value < 0.15)


############################################################################################

print(str(colData(myCDS)))

theme_set(theme_gray(base_size = 8))

# Get the partitions data
colData(myCDS)$partitionLabel <- partitions(myCDS) # alternately, clusters(myCDS)
colData(myCDS)$partitionLabel <- as.character(colData(myCDS)$partitionLabel)

# Test for markers. Need to figure out hte unknown partitions (1 and 4, as of July 8th 2019)
#marker_test_res = top_markers(myCDS, group_cells_by="partition", reference_cells=100, cores=8) # Using 100, not 1000 as default in the monocle manual
marker_res_300 = top_markers(myCDS, group_cells_by="partition", reference_cells=300, cores=8) 

markerDF = as.data.frame(marker_res_300 %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(20, pseudo_R2))

markerDF[markerDF$cell_group == 3,]

# Get top markers for the partitions
top_specific_markers = marker_res_300 %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(5, pseudo_R2)
top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
# Plot these markers
png("./plots/markers_for_partitions.png", width=1400, height=1400, res=300)
plot_genes_by_group(myCDS,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)
dev.off()

# Bar Plots


# Get the barplot of cells in each partition, by sample
barTable <- with(colData(myCDS), table(sample_split_by_hash, partitionGuess))

barTableAsDF = as.data.frame(barTable)

# Get confidence intervals
ciVecUpper = c()
ciVecLower = c()
library(MultinomialCI)
for (eachRow in 1:nrow(barTable)){
    # thisRow = as.vector(barTable[eachRow,1:ncol(barTable)])
    thisRow = as.vector(barTable[eachRow, 1:ncol(barTable)])
    thisCI = multinomialCI(thisRow, .95)
    #ciVec <- c(ciVec, thisCI)
    ciVecUpper <- c(ciVecUpper, as.vector(thisCI[,2]))
    ciVecLower <- c(ciVecLower, as.vector(thisCI[,1]))
}
# Reformat to go the same way as the barTable DF 
UpperCImat = matrix(ciVecUpper, nrow=nrow(barTable), byrow=TRUE)
LowerCImat = matrix(ciVecLower, nrow=nrow(barTable), byrow=TRUE)
# Re-cast to vector, will fill by column and match formatting for plotting
ciVecUpper = as.vector(UpperCImat)
ciVecLower = as.vector(LowerCImat)



# Barplot of counts
png("./plots/cellCountsInPartitionBySample.png", width=1400, height=1000, res=250 )
barPlot <- ggplot(as.data.frame(barTable), aes(factor(partitionGuess), Freq,  fill=sample_split_by_hash)) + 
    geom_col(position="dodge")
    #geom_bar(position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #ggtitle("Cell Counts Passing UMI Minimum")
print(barPlot)
dev.off()

# Get proportions
# propTable <- prop.table(barTable, margin=1)
# png("./plots/cellProportionsInPartitionBySample.png", width=1400, height=1000, res=250 )
# barPlot <- ggplot(as.data.frame(propTable), aes(factor(partitionGuess), Freq,  fill=sample_split_by_hash)) + 
#     geom_col(position="dodge") 
#     #geom_bar(position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#     #ggtitle("Cell Counts Passing UMI Minimum")
# print(barPlot)
# dev.off()

# Get proportions with error bars
propTable <- prop.table(barTable, margin=1)
png("./plots/cellProportionsInPartitionBySample.png", width=1400, height=1000, res=250 )
barPlot <- ggplot(as.data.frame(propTable), aes(factor(partitionGuess), Freq,  fill=sample_split_by_hash)) + 
    geom_col(position="dodge") +
    geom_errorbar(aes(ymin=ciVecLower, ymax=ciVecUpper), width=.1, position=position_dodge(.9))
    #geom_bar(position="dodge", stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #ggtitle("Cell Counts Passing UMI Minimum")
print(barPlot)
dev.off()






print("All Done")