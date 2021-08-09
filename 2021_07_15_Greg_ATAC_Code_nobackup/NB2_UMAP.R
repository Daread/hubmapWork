
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
dir.create(paste0(basepath, "archr/results/"))
out_dir = paste0(basepath, "archr/results/Heart_ATAC/")
dir.create(out_dir)
setwd(out_dir)

# DFR add
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
# modStatus <- loadMonoclePackages()


suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(reshape2)
  library(ggridges)
  library(UpSetR)
  library(ggrastr)
})
DelayedArray:::set_verbose_block_processing(verbose = TRUE)
options(DelayedArray.block.size=1000e7)


## Load saved CDS

# samples = c("W134.heart.apex.s1", 
#               "W135.heart.LV.s1", 
#               "W136.heart.apex.s1", "W136.heart.LV.s1",
#               "W137.heart.apex.s1", # Failed in RNA
#               "W139.heart.apex.s1", "W139.heart.LV.s1", "W139.heart.RV.s1", "W139.heart.septum.s1", 

#               "W142.heart.LV.s1",
#                "W144.heart.apex.s1", 
#                "W145.heart.apex.s1", "W145.heart.LV.s1",
#               "W146.heart.apex.s1", "W146.heart.LV.s1")
# samples = c("allHeartATAC", "W144.heart.apex.s1", "W135.heart.LV.s1")
# samples = c("allHeartATAC", "W144.heart.apex.s1",
#       "W146.heart.LV.s1", "W134.heart.apex.s1")
samples = c("allHeartATAC")

processingNote = "FRIP=0.3_FRIT=0.1UMI=1000"
samples = paste0(samples, processingNote)
useMNN = TRUE

dropFirstComponent = TRUE

# create summary table for samples 
df = data.frame(sample = character(), n_filteredCells = numeric(), n_filteredFeatures = numeric(),
                medUMI = numeric(), medFRIP = numeric(), medFRIT = numeric())

for(s in samples){
  load(paste0(basepath, "cds_objects/cds_p_", s))
  #load("/Users/gregorybooth/UW/booth/Trapnell_lab/Projects/atac/200207_hubmap_NOVASEQ/cds_p_Trisomy_brain_sentinel_new")
  p.cutoff = ncol(cds_p)*0.002
  FRIP_upper = quantile(colData(cds_p)$FRIP, .995)
  fcells = dplyr::filter(as.data.frame(colData(cds_p)), doublet == "singlet", FRIP < FRIP_upper)
  cds_f = cds_p[,colData(cds_p)$cell %in% fcells$cell]
  cds_f <- detect_genes(cds_f)
  cds_f = cds_f[rowData(cds_f)$num_cells_expressed > p.cutoff,]

  if (!("sampleName" %in% colnames(colData(cds_f)))){
    colData(cds_f)$sampleName = s
  }
  

  # reduce dimensions
  set.seed(2017) # ensures reproducibility of previous random number generation
  print("Prelim Processing")
  cds_f <- estimate_size_factors(cds_f)
  cds_f = preprocess_cds(cds_f, method = "LSI", num_dimensions=50)
  #cds_f = align_cds(cds_f, preprocess_method = "LSI", residual_model_formula_str = ~log(umi))
  print("UMAP Reduction Now")
  if (dropFirstComponent){
    reducedDim(cds_f) <- reducedDim(cds_f)[,2:50] # removes 1st LSI component
    processingNote = paste0("DropFirstComponent_", processingNote)
    } else{
      processingNote = paste0("KeepComponents_", processingNote)
    }
  
  # Now try MNN this time
  # DFR
  if (useMNN){
    processingNote = paste0(processingNote, "MNN_offLSI_")
    cds_f = align_cds(cds_f, preprocess_method = "LSI", residual_model_formula_str = ~log(umi))
    cds_f = reduce_dimension(cds_f, reduction_method = 'UMAP', preprocess_method = "Aligned")
  } else{
    processingNote = paste0(processingNote, "noMNN")
    cds_f = reduce_dimension(cds_f, reduction_method = 'UMAP', preprocess_method = "LSI")
  }
  
  print("Clustering cells")
  cds_f = cluster_cells(cds_f)
  
  print("Now plotting UMAP")
  ##################################################################################
  # custom plot
  colData(cds_f)$UMAP_1 <- reducedDims(cds_f)$UMAP[,1]
  colData(cds_f)$UMAP_2 <- reducedDims(cds_f)$UMAP[,2]
  colData(cds_f)$cluster <- cds_f@clusters$UMAP$clusters
  TCdat_pr = data.frame(colData(cds_f))
  
  # png(paste0("monocle3_UMAP", processingNote, "_", s, ".png"), width = 600, height = 400, res=200)
  #   print(ggplot(TCdat_pr, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  #           geom_point_rast(data = TCdat_pr, size = 0.75, stroke = 0) +
  #           theme(legend.position = "none", text = element_text(size = 12),  
  #               legend.key.width = unit(0.5,"line"), legend.key.height = unit(0.5,"line")) + 
  #           labs(#title = s, 
  #             x = "Component 1", 
  #             y = "Component 2") +
  #           guides(guides(colour = guide_legend(override.aes = list(size=1.5)))) +
  #           monocle3:::monocle_theme_opts())
  # dev.off()

  # Shuffle order
  cds_f = cds_f[sample(1:nrow(cds_f)), sample(1:ncol(cds_f))]

  outPath = paste0("./", s, "/")
  dir.create(outPath)

  plotUMAP_Monocle(cds_f, paste0(processingNote, "_", s), "cluster", outputPath=outPath)
  plotUMAP_Monocle(cds_f, paste0(processingNote, "_", s), "sampleName", outputPath = outPath, show_labels=FALSE)
  plotUMAP_Monocle(cds_f, paste0(processingNote, "_", s), "umi", outputPath = outPath, show_labels=FALSE)
  plotUMAP_Monocle(cds_f, paste0(processingNote, "_", s), "FRIP", outputPath = outPath, show_labels=FALSE)
  plotUMAP_Monocle(cds_f, paste0(processingNote, "_", s), "FRIT", outputPath = outPath, show_labels=FALSE)

  colData(cds_f)$log10UMI = log10(colData(cds_f)$umi)
  plotUMAP_Monocle(cds_f, paste0(processingNote, "_", s), "log10UMI", outputPath = outPath, show_labels=FALSE)
  colData(cds_f)$firstLSI = reducedDim(cds_f)[,1]
  plotUMAP_Monocle(cds_f, paste0(processingNote, "_", s), "firstLSI", outputPath = outPath, show_labels=FALSE)

  ## aggregate sample data for sumamry table
  tissue = s
  ncells = ncol(cds_f)
  nfeat = nrow(cds_f)
  medumi = median(colData(cds_f)$umi)
  medfrip = median(colData(cds_f)$FRIP)
  medfrit = median(colData(cds_f)$FRIT)
  sDAT = data.frame(sample = tissue, n_filteredCells = ncells, n_filteredFeatures =  nfeat,
                    medUMI = medumi, medFRIP = medfrip, medFRIT = medfrit)
  
  
  
  df = rbind(df, sDAT)
}





# # DFR
# ##################################
# basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
# dir.create(paste0(basepath, "plots/"))
# out_dir = paste0(basepath, "plots/")
# dir.create(out_dir)
# setwd(out_dir)


# processingNote = "noMNN"
# # Shuffle the cds
# set.seed(7)
# cds_f = cds_f[sample(1:nrow(cds_f)), sample(1:ncol(cds_f))]
# plotUMAP_Monocle(cds_f, processingNote, "cluster", outputPath=out_dir)
# plotUMAP_Monocle(cds_f, processingNote, "sampleName", outputPath=out_dir, show_labels=FALSE)

# # 
# alignedCDS = preprocess_cds(cds_f, method = "PCA", num_dimensions=50)
# alignedCDS = align_cds(alignedCDS, alignment_group = "sampleName")
# processingNote = "mnnBySample"
# alignedCDS = reduce_dimension(alignedCDS, preprocess_method="Aligned")
# plotUMAP_Monocle(alignedCDS, processingNote, "cluster", outputPath=out_dir)
# plotUMAP_Monocle(alignedCDS, processingNote, "sampleName", outputPath=out_dir, show_labels=FALSE)

# processingNote = paste0(processingNote, "dropFirstPC")



# reducedDim(alignedCDS) = reducedDim(alignedCDS)[,2:29]


# alignedCDSTwoOn = reduce_dimension(alignedCDS, preprocess_method="Aligned")
# plotUMAP_Monocle(alignedCDSTwoOn, processingNote, "cluster", outputPath=out_dir)
# plotUMAP_Monocle(alignedCDSTwoOn, processingNote, "sampleName", outputPath=out_dir, show_labels=FALSE)

#   # reducedDim(cds_f) <- reducedDim(cds_f)[,2:50] 


# # Now also try aligning and plotting
# processingNote = "mnnBySampleUsingLSI"
# alignByLSICDS = align_cds(cds_f, alignment_group = "sampleName", preprocess_method = "LSI")
# alignByLSICDS = reduce_dimension(alignByLSICDS, preprocess_method="Aligned")




# plotUMAP_Monocle(alignByLSICDS, processingNote, "cluster", outputPath=out_dir)
# plotUMAP_Monocle(alignByLSICDS, processingNote, "sampleName", outputPath=out_dir, show_labels=FALSE)







####################################################














write.table(df, file = paste0(out_dir,"SummaryTable.txt"), sep = "\t", quote =FALSE, row.names = FALSE, col.names = TRUE)

####################
# combine Col data from each sample for summary plots

df = data.frame(
  sample = character(), 
  cell = character(), 
  umi = numeric(),
  FRIP = numeric(), 
  FRIT = numeric(), 
  doublet = logical())

for(s in samples){
  load(paste0(basepath, "cds_objects/cds_p_", s))
  
  cdat = data.frame(colData(cds_p)) %>% 
    mutate(sample = s) %>% 
    select(sample, cell, umi, FRIP, FRIT, doublet)

  df = rbind(df, cdat)
}

# number of cells by sample
ncells = group_by(df, sample) %>% 
  dplyr::summarise(ncells = n())

ggplot(ncells) +
  geom_col(aes(x = sample, y = ncells)) +
  scale_color_brewer(palette='Set1') +
  #xlab("sample") +
  ylab("cells") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "barplot_cells_by_sample.png",
       width = 10, height = 3)

# UMI by dose
ggplot(df) +
  geom_boxplot(aes(x =  sample, y = log10(umi))) +
  #xlab("sample") +
  ylab("log10(Frags per cell)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "boxplot_FragsPerCell_by_sample.png",
       width = 10, height = 3)

# FRIP by sample 
ggplot(df) +
  geom_boxplot(aes(x =  sample, y = FRIP)) +
  #xlab("sample") +
  ylab("FRIP") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "boxplot_FRIP_by_sample.png",
       width = 10, height = 3)

# FRIT by sample
ggplot(df) +
  geom_boxplot(aes(x =  sample, y = FRIT)) +
  #xlab("sample") +
  ylab("FRIT") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(filename = "boxplot_FRIT_by_sample.png",
       width = 10, height = 3)


