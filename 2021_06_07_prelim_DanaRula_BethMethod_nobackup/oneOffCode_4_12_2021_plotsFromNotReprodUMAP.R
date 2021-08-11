

hardAssignClustLabels <- function(inputCDS, kVal){
	if (kVal == 40){
		colData(inputCDS)$cluster_label = as.character(clusters(inputCDS))
        colData(inputCDS)$cluster_cell_type = colData(inputCDS)$cluster_label


        colData(inputCDS)$cluster_cell_type = ifelse(colData(inputCDS)$cluster_cell_type %in% c("12", "4"),
                                            "Cardiomyocyte", colData(inputCDS)$cluster_cell_type)
        colData(inputCDS)$cluster_cell_type = ifelse(colData(inputCDS)$cluster_cell_type  %in% c("1", "21", "8", "7", "3"),
                                            "Fibroblasts", colData(inputCDS)$cluster_cell_type)
        colData(inputCDS)$cluster_cell_type = ifelse(colData(inputCDS)$cluster_cell_type %in% c("2", "9", "11"),
                                            "Macrophages", colData(inputCDS)$cluster_cell_type) #RBPJ
        colData(inputCDS)$cluster_cell_type = ifelse(colData(inputCDS)$cluster_cell_type %in% c("6", "14", "13", "15"),
                                            "Endothelium", colData(inputCDS)$cluster_cell_type) #LDB@
        colData(inputCDS)$cluster_cell_type = ifelse(colData(inputCDS)$cluster_cell_type %in% c("10"),
                                            "Vascular Smooth Muscle", colData(inputCDS)$cluster_cell_type) #
        colData(inputCDS)$cluster_cell_type = ifelse(colData(inputCDS)$cluster_cell_type %in% c("5"),
                                            "T Cells", colData(inputCDS)$cluster_cell_type) #
        colData(inputCDS)$cluster_cell_type = ifelse(colData(inputCDS)$cluster_cell_type %in% c("19", "22", "16"),
                                            "B Cells", colData(inputCDS)$cluster_cell_type) #
        colData(inputCDS)$cluster_cell_type = ifelse(colData(inputCDS)$cluster_cell_type %in% c("17"),
                                            "Neuronal", colData(inputCDS)$cluster_cell_type) #
        colData(inputCDS)$cluster_cell_type = ifelse(colData(inputCDS)$cluster_cell_type %in% c("18"),
                                            "Adipocytes", colData(inputCDS)$cluster_cell_type) #
	}

	names(colData(inputCDS)$cluster_cell_type) = rownames(colData(inputCDS))
	inputCDS@clusters[["UMAP"]]$clusters <- as.factor(colData(inputCDS)$cluster_cell_type)
	plotUMAP_Monocle(inputCDS, processingNote, "cluster_cell_type", show_labels=TRUE)

	return(inputCDS)
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


hardAssignAnatomicalSitesDF <- function(inputDF){

	inputDF$Anatomical_Site = "Not_Specified"
	inputDF$Anatomical_Site = ifelse(grepl("Apex", inputDF$sample),
										"Apex", inputDF$Anatomical_Site)
	inputDF$Anatomical_Site = ifelse(grepl("Septum", inputDF$sample),
										"Septum", inputDF$Anatomical_Site)
	inputDF$Anatomical_Site = ifelse(grepl("Left.Vent", inputDF$sample),
										"Left_Vent", inputDF$Anatomical_Site)
	inputDF$Anatomical_Site = ifelse(grepl("Right.Vent", inputDF$sample),
										"Right_Vent", inputDF$Anatomical_Site)

	return(inputDF)
}



source("../singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()
print(modStatus)
# sessionInfo()

source("../backgroundMethods.R")

# processingNote = "HM10"

# useScrublet = FALSE
# writeMatrixForScrub = TRUE
# scrubCutoff = .2
# umiCellMin = 100
# bgUMImax = 15

# processingNote = paste0(processingNote, "UMI=", as.character(umiCellMin), 
#     "BG=", as.character(bgUMImax))

# print("Starting Now")

# dataPath <- "./final-output/"
# #myCDS <- readRDS(paste0(dataPath, "cds.RDS"))
# # myCDS <- load.cds(paste0(dataPath, "UMI.count.matrix"), paste0(dataPath, "gene.annotations"),
# #          paste0(dataPath, "cell.annotations"))
# # print("CDS Loaded")

# samplesToRead <- c("W134.Apex", "W135.Left.Vent", "W136.Apex", "W136.Left.Vent",
#                    #"W137.Apex", 
#                    "W139.Apex", "W139.Left.Vent", "W139.Right.Vent",
#                    "W139.Septum", "W142.Left.Vent", "W144.Apex", "W145.Apex",
#                    "W145.Left.Vent", "W146.Apex", "W146.Left.Vent")

# # myCDS <- readRDS(paste0(dataPath, "cdsWithHash.RDS"))
# myCDS <- makeCDS_fromBBIoutputDirs(samplesToRead)

# colData(myCDS)$Cell = rownames(colData(myCDS))
# colData(myCDS)$cell = colData(myCDS)$Cell


# # Get the background and main cell subsets
# bgCDS = myCDS[,colData(myCDS)$n.umi <= bgUMImax]
# fullCDS = myCDS
# myCDS = myCDS[,colData(myCDS)$n.umi >= umiCellMin]
# processingNote = paste0(processingNote, "UMI=", as.character(umiCellMin), 
#     "BG=", as.character(bgUMImax))
# cellNames = colData(myCDS)$Cell

# # Color by the PCR well
# splitCellCol = (separate(as.data.frame(colData(myCDS)), cell, c("P7", "P5", "RT_Well",
#          "Lig_Well", "SampleHolderNum"),
#         sep="_"))
# colData(myCDS)$P7_ind = splitCellCol$P7
# colData(myCDS)$P5_ind = splitCellCol$P5
# colData(myCDS)$RT_Well = splitCellCol$RT_Well
# colData(myCDS)$Lig_Well = splitCellCol$Lig_Well

# colData(myCDS)$sample = as.character(colData(myCDS)$sampleName)


# print("CDS Loaded")

# # Log10 scale is easier to plot
# colData(myCDS)$log_10_umi <- log10(colData(myCDS)$n.umi)


# ################################################################################################################
# # Barnyard plot and species labeling
# print("Making species assignments and barnyard plot")
# human.genes = rownames(rowData(myCDS))[grepl("^ENSG", rowData(myCDS)$id)]
# mouse.genes = rownames(rowData(myCDS))[grepl("^ENSMUSG", rowData(myCDS)$id)]
# print(str(human.genes))
# colData(myCDS)$n_human_umi = Matrix::colSums(exprs(myCDS)[human.genes,])
# colData(myCDS)$n_mouse_umi = Matrix::colSums(exprs(myCDS)[mouse.genes,])
# print("Loaded species-specific UMIs")

# # ggplot(colData(myCDS))  + 
# #   geom_point(aes(x = n_human_umi, y = n_mouse_umi), size = .75, stroke = 0) +
# #   scale_color_manual(values = c("TRUE" = "dimgrey", "FALSE" = "black")) +
# #   xlim(0,40000) + ylim(0,40000) + 
# #   monocle3:::monocle_theme_opts()
# # ggsave("./plots/Lung_barnyard.png")

# colData(myCDS)$Species <- ifelse(colData(myCDS)$n_human_umi/colData(myCDS)$n_mouse_umi > 10.0, "Human",
#                             ifelse(colData(myCDS)$n_mouse_umi/colData(myCDS)$n_human_umi > 10.0, "Mouse", "Mixed"))

# colData(myCDS)$sample <- as.character(colData(myCDS)$sample)

# boxplot_stat_by_X(myCDS, processingNote, "log_10_umi", "sample")
# # # See UMI distributions
# # ggplot(as.data.frame(colData(myCDS)), aes(x=sample, y=log_10_umi)) + 
# #     geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# # ggsave("./plots/UMIs_by_Sample.png")
# # # ggsave("./plots/Scrambled_UMIs_by_Sample.png")

# # boxplot_stat_by_X(myCDS, processingNote, "log_10_umi", "sample")
# counts = (as.data.frame(colData(myCDS)) %>% dplyr::group_by(sample) %>%
#              dplyr::summarise(n=n()))
# png(paste0("./plots/sampleCellCounts", processingNote, ".png"))
# ggplot(as.data.frame(counts), aes_string(x="sample", y="n")) +
#     geom_bar(position="dodge", stat="identity")+ 
#     theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#     labs(title="Cells per sample")
# dev.off()

# preUMAPCDS = myCDS
# preUMAPNotes = processingNote

# tuckerHeartMarkers = c("RBPJ", "LDB2", "DLC1", "NRXN1", "NEGR1",
#                     "GSN", "FHL2", "TTN")

# # Make UMAPs w/ color-coding for each sample, alone
# # for (eachSample in samplesToRead){
# #     subsetCDS = myCDS[,colData(myCDS)$sample == eachSample]
# #     # PCA, UMAP, show some heart markers
# #     print(paste0("Processing ", eachSample,  " alone"))
# #     set.seed(7)
# #     subsetCDS <- preprocess_cds(subsetCDS)
# #     subsetCDS <- reduce_dimension(subsetCDS)
# #     plotUMAP_Monocle_genes(subsetCDS, paste0(eachSample, "only", processingNote),
# #         tuckerHeartMarkers, "Tucker_et_al_Markers")
# # }

# # Get some summary plots
# medianStats = (as.data.frame(colData(myCDS)) %>%  group_by(sample) 
# 				%>% summarize(the_medians = median(n.umi)))

# # # Pre-process and plot
# set.seed(7)
# myCDS = preprocess_cds(myCDS) #, num_dim = 50)
# print("CDS Pre-processed")
# # png(paste0("./plots/", processingNote, "_varianceExplained.png"))
# # plot_pc_variance_explained(myCDS)
# # dev.off()

# tryAlign = TRUE
# set.seed(7)
# if (tryAlign){
#     print("Aligning data by sample")
#     alignBy = "sample"
#     myCDS = monocle3::align_cds(myCDS, num_dim=50, alignment_group = alignBy)
#     processingNote <- paste0(processingNote, "MNN=", alignBy)
# } else{
#     processingNote <- paste0(processingNote, "_noMNN")
# }

# set.seed(7)
# myCDS = reduce_dimension(myCDS, umap.fast_sgd=FALSE)

# # See some plots
# plotUMAP_Monocle(myCDS, processingNote, "sample", show_labels=FALSE)

# plotUMAP_Monocle_genes(myCDS, processingNote, 
#         tuckerHeartMarkers, "Markers_TuckerEtAl")




# plotUMAP_Monocle_genes(myCDS, processingNote, 
#         c("NEGR1", "ELN", "GSN", "NRXN1", 
#             "NRXN3", "CADM2", "ATP1B3", "CRISPLD2"), "_Fibroblast_vs_neuronal_TuckerEtAl")
# plotUMAP_Monocle_genes(myCDS, processingNote, 
#         c("PROX1", "FLT4", "PDPN", "BMX", "NPR3"), "_endothelial_subtyping_TuckerEtAl")
# plotUMAP_Monocle_genes(myCDS, processingNote, 
#         c("PRKG1", "RBMS3", "PLA2G5", "GMDS", "RP11-767I20.1"), "_endothelial_2_aka_endocardium_TuckerEtAl")

# plotUMAP_Monocle_genes(myCDS, processingNote, 
#         c("SYNE1", "PTPRM", "LDB2", "LINC00486"), "_pandEndothelial_TuckerEtAl")

# plotUMAP_Monocle_genes(myCDS, processingNote, 
#         c("KCNAB1", "CARMN", "ITGA8", "MYH11"), "_vascularSmoothMus_TuckerEtAl")

# plotUMAP_Monocle_genes(myCDS, processingNote, 
#         c("ACACB", "GPAM", "PDE3B", "EBF1"), "_adipocytes_TuckerEtAl")




# # Search for markers
# deRes = runDEtestingToID_markers(myCDS, processingNote, "cluster_label")

# # Code for hard-code labeling of heart clusters

# # Cluster for partition-based coloring
# print("Clustering Cells Now")
# set.seed(7)
# kVal = 40
# myCDS = cluster_cells(myCDS, k = kVal)#, resolution=c(10^seq(-6,-1))) #, 
# # myCDS = cluster_cells(myCDS), k = 3) #, resolution=c(10^seq(-6,-1)) )#, 
#        #python_home="/net/gs/vol3/software/modules-sw-python/2.7.3/louvain")# /381b7db/Linux/RHEL6/x86_64/lib/python2.7/site-packages/python_louvain-0.13-py2.7.egg/")

#         #community/")

# # Get the partitions data
# colData(myCDS)$partitionLabel <- partitions(myCDS) # alternately, clusters(myCDS)
# colData(myCDS)$partitionLabel <- as.character(colData(myCDS)$partitionLabel)
# colData(myCDS)$partitionGuess <- as.character(colData(myCDS)$partitionLabel)
# colData(myCDS)$cluster_label = as.character(clusters(myCDS))
# plotUMAP_Monocle(myCDS, paste0("k=", as.character(kVal), "_", processingNote), "cluster_label", show_labels=TRUE)
# plotUMAP_Monocle(myCDS, paste0("k=", as.character(kVal), "_", processingNote), "partitionLabel", show_labels=TRUE)

# # testTable = with(colData(myCDS), table(sample, partitionLabel))

# myCDS <- hardAssignClustLabels(myCDS, kVal)
# names(colData(myCDS)$cluster_cell_type) = rownames(colData(myCDS))
# myCDS@clusters[["UMAP"]]$clusters <- as.factor(colData(myCDS)$cluster_cell_type)
# plotUMAP_Monocle(myCDS, processingNote, "cluster_cell_type", show_labels=TRUE)



# # Save an RDS
# # saveRDS(myCDS, paste0("./rdsOutput/basicCellTyped_", processingNote, "_cds.RDS"))

processingNote = "HM10UMI=100BG=15MNN=sample"
# Note anatomical sites
myCDS = readRDS(paste0("./rdsOutput/basicCellTyped_", processingNote, "_cds.RDS"))


myCDS = hardAssignAnatomicalSites(myCDS)
# Plot the distributions

# testSplit = (as.data.frame(colData(myCDS)) %>% 
# 			group_by(sample, Anatomical_Site, cluster_cell_type) 
# 			 )
# str(testSplit)

processingNote = paste0("ReloadUnreproUMAP", processingNote)

countTable = with(colData(myCDS), table(sample, cluster_cell_type))
proportionTable <- prop.table(countTable, margin=1)
propDF = as.data.frame(proportionTable)

propDF

propDF <- hardAssignAnatomicalSitesDF(propDF)
# Remove the Unassigned clusters
propDF$cluster_cell_type = as.character(propDF$cluster_cell_type)
propDF = propDF[!(propDF$cluster_cell_type %in% c("20")),]
propDF$cluster_cell_type = ifelse(propDF$cluster_cell_type == "Vascular Smooth Muscle",
				"Vascular Smooth\nMuscle", propDF$cluster_cell_type)
colnames(propDF) = c("sample", "cluster_cell_type", "Freq", "Anatomical_Site")

png(paste0("./plots/", processingNote, "propsByAnatomicalSite.png"))
myPlot <- (ggplot(propDF, 
		aes_string(x="cluster_cell_type", y="Freq", fill="Anatomical_Site",
					color="Anatomical_Site")) +
		geom_boxplot() + 
		ylab("Proportion of Total") + xlab("Cell Type") +
		#guide_legend(title="Anatomical Site") +
		#scale_fill_discrete(name = "Anatomical Site")+
		theme(text = element_text(size = 12))    +
		 theme(axis.text.x = element_text(angle = 90, hjust = 1)))
print(myPlot)
dev.off()


# Split by facet


png(paste0("./plots/", processingNote, "_FACETED_propsByAnatomicalSite.png"))
myPlot <- (ggplot(propDF, 
        aes_string(x="Anatomical_Site", y="Freq", color="Anatomical_Site")) +
        
        #facet_grid(cols=vars(cluster_cell_type)) +

         facet_wrap(vars(cluster_cell_type), scales="free") +
        geom_boxplot() + 
        ylab("Proportion of Total") + xlab("") +

        #guide_legend(title="Anatomical Site") +
        #scale_fill_discrete(name = "Anatomical Site")+
        theme(text = element_text(size = 12))    +
         theme(axis.text.x = element_text(angle = 90, hjust = 1))   +
         theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())


         )
print(myPlot)
dev.off()













#
colData(myCDS)$cluster_cell_type = ifelse(
    colData(myCDS)$cluster_cell_type == "20", " ", colData(myCDS)$cluster_cell_type)


plotUMAP_Monocle <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, textSize=4){ #, xVal, yVal){
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



