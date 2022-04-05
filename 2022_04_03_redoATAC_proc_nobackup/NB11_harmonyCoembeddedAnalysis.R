
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/"
dir.create(paste0(basepath, "archr/results/NB11/"))
out_dir = paste0(basepath, "archr/results/NB11/")
setwd(out_dir)
set.seed(7)

# load requirements
suppressPackageStartupMessages({
  # library(ArchR)
  library(monocle3)
  library(dplyr)
  library(ggrastr)
  # library(Seurat)
})
library(tidyr)
library(plyr)
library(ggplot2)


source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

library("optparse")

# Get the passed parameters
option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
  			default="PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5", 
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-s", "--sampleRNAname"), type="character",  default= "All_Cells", #  "All_Cells",  # default="W142.Left.Vent",  "W144.Apex"
              help="Sample to integrate", metavar="character"),
    make_option(c("-u", "--useSubset"), type="character", default="Subset", 
              help="Subset or All", metavar="character"),

    make_option(c("-m", "--useMNN"), type="character", default="useMNN", 
              help="Subset or All", metavar="character"),

    make_option(c("-c", "--fripMin"), type="numeric", default=0.2, 
              help="Min FRIP value to permit", metavar="numeric"),

    make_option(c("-w", "--cdsToSave"), type="character", default="Peak_CDS",   # Peak_CDS
              help="Features of cds to save", metavar="character"),

    make_option(c("-k", "--kVal"), type="numeric", default=20, 
              help="K value for clustering", metavar="character"),

    make_option(c("-r", "--regressProto"), type="logical", default=TRUE, 
              help="Regress out protocadherins and re-run UMAP + KNN assignment?", metavar="logical"),

    make_option(c("-d", "--downsample"), type="logical", default=FALSE, 
              help="Down-sample to 1/10th of cells to speed up testing/debugging?", metavar="logical"),

    make_option(c("-z", "--pcToUse"), type="numeric", default=50, 
              help="Principel components/LSI coords to use", metavar="numeric"),

    make_option(c("-f", "--featureSet"), type="character", default="peakMat", # "peakMat" for greg/riza matrix, "bMat" for archr bins
                                                                              # "gMat" for archr activity scores
              help="peakMat, bMat, or gMat", metavar="character"),

    make_option(c("-a", "--ATACprocNote"), type="character", default="FRIP=0.1_FRIT=0.1UMI=1000DL=0.5", 
              help="How ATAC Cells were filtered", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(opt)

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


hardAssignATAC_clusters <- function(inputCDS, processingNote, kValUsed){
  # 8-6-21: Using k = 15 I made a prelim assignment
  if ((processingNote == "FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50") & (kValUsed == 20)){
    colData(inputCDS)$Assigned_Cell_Type = "Unassigned"

    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(22), 
                    "Endocardium", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(6, 23), 
                    "Macrophage", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(17,25), 
                    "Perivascular", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(16,2,11), 
                    "Fibroblast", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(24,14,10,8,9,19,12), 
                    "Vascular_Endothelium", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(1,15,21,4,13,20,7,3,5,18), 
                    "Cardiomyocyte", colData(inputCDS)$Assigned_Cell_Type)
  }
  return(inputCDS)
}

plotUMAP_MonocleModded <- function(dataCDS, processingNote, catToColor,
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

plotUMAP_MonocleModded_tech <- function(dataCDS, processingNote, catToColor,
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







# Get the saved CDS
rdsPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/rdsOutput/"
# fileName = paste0(opt$sampleRNAname, "_", opt$ATACprocNote, "gMatrixCDS_postTransfer.RDS")
fileName = paste0("All_Cells", "_", opt$ATACprocNote, "gMatrixCDS_postTransfer.RDS")

# Read in the data
cds.a.g = readRDS(paste0(rdsPath, fileName))
# bMatFile = paste0(opt$sampleRNAname, "_", opt$ATACprocNote, "bMatrixCDS_postTransfer.RDS")
bMatFile = paste0("All_Cells", "_", opt$ATACprocNote, "bMatrixCDS_postTransfer.RDS")
cds.a.b = readRDS(paste0(rdsPath, bMatFile))

# Get data from Greg's filtering process
# Added 7-29-21: Filter by checking names in a cds filtered by greg/riza's method
# filterCDSName = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/cds_objects/"
# processingNote = opt$ATACprocNote

# oldProcNote = "FRIP=0.2_FRIT=0.05UMI=1000"
# # cds_p holds data
# load(paste0(filterCDSName, "cds_p_allHeartATAC", oldProcNote))
peakCDSPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/archr/results/NB9/All_Cells/"
peakCDS_name = "FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20"

cds_p = readRDS(paste0 (peakCDSPath, peakCDS_name, ".rds"))


# Added 10-26-21: Get the CDS with gene activity scores calculated by cicero links
ciceroActivPath = paste0("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_24_ciceroWork_nobackup/fileOutputs/", 
                       "CiceroActivityCDS_", "Cicero_Defaults_cds_p_allHeartATAC", opt$ATACprocNote, "_ciceroConnections.RDS")
ciceroActivCDS = readRDS(ciceroActivPath)


# Added 11-10-21: Get the CDS with harmony alignments
harmonyPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/archr/results/NB10/All_Cells/"
harmonyFile = paste0(opt$processingNote, "_ArchRActivharmonyAligned_cds.rds")

harmonyCDS = readRDS(paste0(harmonyPath, harmonyFile))




# Only keep the ATAC data
harmonyATACandRNA = harmonyCDS
harmonyCDS = harmonyCDS[,colData(harmonyCDS)$tech == "ATAC"]


##################################################################################### harmony cds formatting
# Reformat sampleName
colData(harmonyCDS)$cellName = rownames(colData(harmonyCDS))
colDataDF = separate(data=as.data.frame(colData(harmonyCDS)), col = "cellName", 
                        into = c("Part1", "Part2", "Part3", "Discard"), sep="_", remove=TRUE)
# colData(harmonyCDS)$sampleName = colDataDF$sampleName
colData(harmonyCDS)$cellName = paste0(colDataDF$Part1, "_", colDataDF$Part2, "_", colDataDF$Part3)

rownames(colData(harmonyCDS)) = colData(harmonyCDS)$cellName
colnames(harmonyCDS) = rownames(colData(harmonyCDS))

#####################################################################################


##################################################################################### cicero activity cds formatting
# Reformat sampleName
colDataDF = separate(data=as.data.frame(colData(ciceroActivCDS)), col = "sampleName", 
                        into = c("sampleName", "Processing"), sep="FRIP", remove=TRUE)
# colData(ciceroActivCDS)$sampleName = colDataDF$sampleName
ciceroActivCDS$sampleName = colDataDF$sampleName

# Reformat the colnames of the activity matrix to match the bin/arrowActivity cds's
colDataDF$rowname = rownames(colDataDF)
colDataDF = separate(data=colDataDF, col = "rowname", 
                        into = c("Part1", "Part2", "Part3"), sep="_", remove=TRUE)
rownames(colDataDF) = paste0(colDataDF$Part1, "_", colDataDF$Part2, "_", colDataDF$Part3)

rownames(colData(ciceroActivCDS)) = paste0(colData(ciceroActivCDS)$sampleName, "#", rownames((colDataDF)))
colnames(ciceroActivCDS) = rownames(colData(ciceroActivCDS))

#####################################################################################



harmonyCDS = harmonyCDS[,colnames(cds.a.g)]

# Now get the portion of cds_p we want, in order
ciceroActivCDS = ciceroActivCDS[,colnames(cds.a.g)]

cds_p = cds_p[,colnames(cds.a.g)]
colData(cds.a.g)$sampleName = colData(cds_p)$sampleName
colData(cds.a.b)$sampleName = colData(cds_p)$sampleName

colData(harmonyCDS)$seuratCCA_celltype = colData(cds_p)$predicted.id
# All the CDS's are organized with matching names and orders.


if (opt$sampleRNAname != "All_Cells"){
  harmonyCDS = harmonyCDS[,colData(harmonyCDS) == samplesATACnames[opt$sampleRNAname]]
}

if (opt$downsample){
	print("Downsampling now")
	harmonyCDS = harmonyCDS[,sample(1:ncol(harmonyCDS), (ncol(harmonyCDS)/10), replace=FALSE)]
	# Downsample RNA as well
	harmonyRNA = harmonyATACandRNA[,colData(harmonyATACandRNA)$tech == "RNA"]
	harmonyRNA = harmonyRNA[,sample(1:ncol(harmonyRNA), (ncol(harmonyRNA)/10), replace=FALSE)]

	# Recombine for later work
	harmonyATACandRNA = combine_cds(list(harmonyRNA, harmonyCDS))
	opt$processingNote = paste0(opt$processingNote, "downSampled")

	reducedDims(harmonyATACandRNA)$PCA = rbind(reducedDims(harmonyRNA)$PCA, reducedDims(harmonyCDS)$PCA)
	reducedDims(harmonyATACandRNA)$UMAP = rbind(reducedDims(harmonyRNA)$UMAP, reducedDims(harmonyCDS)$UMAP)
} else{
	harmonyRNA = harmonyATACandRNA[,colData(harmonyATACandRNA)$tech == "RNA"]
	harmonyATACandRNA = combine_cds(list(harmonyRNA, harmonyCDS))
	reducedDims(harmonyATACandRNA)$PCA = rbind(reducedDims(harmonyRNA)$PCA, reducedDims(harmonyCDS)$PCA)
	reducedDims(harmonyATACandRNA)$UMAP = rbind(reducedDims(harmonyRNA)$UMAP, reducedDims(harmonyCDS)$UMAP)
}

# prevHarmonyCDS = harmonyATACandRNA

# colData(harmonyCDS)$n.umi = ifelse(is.na(colData(harmonyCDS)$n.umi), colData(harmonyCDS)$umi, colData(harmonyCDS)$n.umi)
# Run extra correction step to regress out protocadherin signal, then re-run UMAP and KNN label transfer?
if (opt$regressProto){
  print("Regressing out protocadherin signal")
  updatedHarmonyCDS = regressProtocadherinAndReumap(harmonyATACandRNA, opt)
  opt$processingNote = paste0("Regress_Protocadherin_", opt$processingNote)
} else{
	updatedHarmonyCDS = harmonyATACandRNA
	if (opt$downsample){
		updatedHarmonyCDS = reduce_dimension(updatedHarmonyCDS)
	}
}


harmonyCDS = updatedHarmonyCDS
harmonyCDS = harmonyCDS[,colData(harmonyCDS)$tech == "ATAC"]




# Save an RDS at this point


# Plot some outputs

out_dir = paste0(basepath, "archr/results/NB11/", opt$sampleRNAname, "/")
dir.create(out_dir)
outputNote = paste0("Harmony_Aligned_Coords", opt$processingNote)




# Tech mixutre
plotUMAP_Monocle(updatedHarmonyCDS, outputNote, "tech",
                      outputPath = out_dir, show_labels=FALSE)









# Test for mixture
plotUMAP_Monocle(harmonyCDS, outputNote, "sampleName",
                      outputPath = out_dir, show_labels=FALSE)

# See how KNN/harmony assigned types match up to the peak-based clusering
plotUMAP_Monocle(harmonyCDS, outputNote, "predicted.id",
                      outputPath = out_dir, show_labels=FALSE)
plotUMAP_Monocle(harmonyCDS, outputNote, "harmonyKNN_type",
                      outputPath = out_dir, show_labels=FALSE)

# See some marker genes
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("SKAP1", "CD247", "ARHGAP15", "RIPOR2" ), 
                      "tcellTopRNA",
                      outputPath = out_dir)

plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("THEMIS", "CD3E", "CD3D", "CD3G" ), 
                      "moreTCell",
                      outputPath = out_dir)


plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("TTN", "RBPJ","GSN", "LDB2", "MYH11", "THEMIS" ),
                      "generalMarkers",
                      outputPath = out_dir)
# 
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("PECAM1", "CD106", "CD62E", "SELE", "KDR", "ENG" ),
                      "coleEndothMark",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("MYH11", "SMTN", "CALD1", "CNN1", "CNN2" ),
                      "coleSmoothMuscMark",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("CD68", "CD14", "FCGR1", "MS4A7", "FCER1A", "CD163", "LY6C1", "FCN1", "MERTK" ),
                      "coleMonoMacMark",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("CD3E", "CD3G", 'CD3D", "CD4", "CD8A', "CD8B", "BCL11B" ),
                      "coleTCellMark",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("CD90", "THY1", "VIM", "DES" ),
                      "coleFibroblastMark",
                      outputPath = out_dir)

plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("PCDHGA1", "PCDHGA2" ),
                      "protocadherins",
                      outputPath = out_dir)

plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("F13A1", "RBM47",  "CD74", "RBPJ" ), # Skip RBPJ
                      "macTopRNA",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("SKAP1", "CD247", "ARHGAP15", "RIPOR2" ), # Skip RBPJ
                      "tcellTopRNA",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("VWF", "LDB2", "ANO2", "FLT1" ), # Skip RBPJ
                      "vascEndTopRNA",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("DLC1", "RGS5", "NOTCH3", "EGFLAM" ), # Skip RBPJ
                      "perivasTopRNA",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("NRXN1", "XKR4", "NRXN3", "CADM2" ), # Skip RBPJ
                      "neuronTopRNA",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("DCN", "GSN", "MGP", "C7" ), # Skip RBPJ
                      "fibroTopRNA",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("PKHD1L1", "PCDH15", "LDB2", "NPR3" ), # Skip RBPJ
                      "endocardTopRNA",
                      outputPath = out_dir)

plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("EMCN", "POSTN", "VWF" ), # Skip RBPJ
                      "endocardTopRNA_2",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("AQP1", "ADAMTSL1", "SMOC1" ), # Skip RBPJ
                      "endocardTopRNA_3",
                      outputPath = out_dir)

plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("TTN", "MYBPC3", "DB3", "RBM20" ), # Skip RBPJ
                      "cardiomyocyteTopRNA",
                      outputPath = out_dir)
plotUMAP_Monocle_genes(harmonyCDS, outputNote, c("BANK1", "PAX5", "CD79A", "RALGPS2" ), # Skip RBPJ
                      "bCellTopRNA",
                      outputPath = out_dir)


cdsAndDEtest = runDEtestingToID_markers(harmonyCDS, outputNote, "harmonyKNN_type",
                  howManyGenesToTest = 50, outputPath=out_dir)


# Save the test result
DEtestResFile = paste0(out_dir, outputNote, "_HarmonyKNN_based_Celltypes_DE_Testing.csv")
write.csv(cdsAndDEtest$marker_test_res[order(cdsAndDEtest$marker_test_res$cell_group),], DEtestResFile)




kVal = 20
updatedHarmonyCDS = cluster_cells(updatedHarmonyCDS, k=kVal)

colData(updatedHarmonyCDS)$cluster_label = clusters(updatedHarmonyCDS)

plotUMAP_Monocle(updatedHarmonyCDS, paste0(outputNote, "_k_", as.character(kVal)), "cluster_label",
                      outputPath = out_dir, show_labels=TRUE)

plotUMAP_MonocleModded_tech(updatedHarmonyCDS, paste0( "faceted",outputNote, "_k_", as.character(kVal)), "cluster_label",
                    show_labels=TRUE, textSize=10, outputPath = out_dir)





hardAssignClusterLables <- function(inputCDS, opt, kVal){

	return(inputCDS)
}




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


updatedHarmonyCDS = knnRNAtoATAC(updatedHarmonyCDS)

# Save output
cdsFile = paste0(out_dir, paste0(outputNote, "_k_", as.character(kVal), "_cds.rds"))
print(cdsFile)
saveRDS(updatedHarmonyCDS[,colData(updatedHarmonyCDS)$tech == "ATAC"], file=cdsFile)

# 
cdsFileRnaandATAC = paste0(out_dir, paste0(outputNote, "_k_", as.character(kVal), "_RNA_and_ATAC_cds.rds"))
print(cdsFileRnaandATAC)
saveRDS(updatedHarmonyCDS, file=cdsFileRnaandATAC)



#######################################################


harmonyCDS = updatedHarmonyCDS
harmonyCDS = harmonyCDS[,colData(harmonyCDS)$tech == "ATAC"]


##################################################################################### harmony cds formatting
# Reformat sampleName
colData(harmonyCDS)$cellName = rownames(colData(harmonyCDS))
colDataDF = separate(data=as.data.frame(colData(harmonyCDS)), col = "cellName", 
                        into = c("Part1", "Part2", "Part3", "Discard"), sep="_", remove=TRUE)
# colData(harmonyCDS)$sampleName = colDataDF$sampleName
colData(harmonyCDS)$cellName = paste0(colDataDF$Part1, "_", colDataDF$Part2, "_", colDataDF$Part3)

rownames(colData(harmonyCDS)) = colData(harmonyCDS)$cellName
colnames(harmonyCDS) = rownames(colData(harmonyCDS))

#####################################################################################

cds_p = cds_p[,colnames(harmonyCDS)]

colData(cds_p)$harmonyKNN_type = colData(harmonyCDS)$harmonyKNN_type
colData(cds_p)$harmony_clust_label = colData(harmonyCDS)$cluster_label




plotUMAP_Monocle(cds_p, paste0("cds_p_Peak_UMAP"), "Assigned_Cell_Type",
                      outputPath = out_dir, show_labels=FALSE)



plotUMAP_Monocle(cds_p, paste0("cds_p_Peak_UMAP"), "harmonyKNN_type",
                      outputPath = out_dir, show_labels=FALSE)

plotUMAP_Monocle(cds_p, paste0("cds_p_Peak_UMAP"), "cluster_label",
                      outputPath = out_dir, show_labels=FALSE)

plotGroupedProportions(cds_p, paste0("peak_clusts_vs_Harmony_Calls_bars"), "Assigned_Cell_Type",
                        "harmonyKNN_type",
                      pathToPlot = out_dir )


= function(inputCDS, processingNote,
                        groupValue, colForProportions,
                        pathToPlot="./plots/", widthToUse=1200,
                        heightToUse=800){



# Save output
cdsFile = paste0(out_dir, paste0(outputNote, "_k_", as.character(kVal), "_peak_cds.rds"))
print(cdsFile)



saveRDS(cds_p, file=cdsFile)










# cdsList = list(harmonyCDS)
# names(cdsList) = c("HarmonyCDS")

# preprocessMethods = c("PCA")
# names(preprocessMethods) = c("HarmonyCDS")

# # cdsList = list(cds_p, cds.a.g, cds.a.b, ciceroActivCDS)
# # names(cdsList) = c("Peak_CDS", "Activity_CDS", "Bin_CDS", "CiceroActiv_CDS")

# # preprocessMethods = c("LSI", "PCA", "LSI", "PCA")
# # names(preprocessMethods) = c("Peak_CDS", "Activity_CDS", "Bin_CDS", "CiceroActiv_CDS")


# # opt$sampleRNAname = "W144.heart.apex.s1"
# # opt$useMNN = "NoMNN"

# out_dir = paste0(basepath, "archr/results/NB9/", opt$sampleRNAname, "/")
# dir.create(out_dir)
# setwd(out_dir)

# # Loop, make a UMAP for each of these and 
# for (eachCDSname in names(cdsList)){
#   thisCDS = cdsList[[eachCDSname]]
#   colData(thisCDS)$Sample = colData(thisCDS)$sampleName
#   colData(thisCDS)$FRIP = colData(cds_p)$FRIP 
#   colData(thisCDS)$umi = colData(cds_p)$umi
#   colData(thisCDS)$log10umi = log10(colData(cds_p)$umi)
#   colData(thisCDS)$doublet_likelihood = colData(cds_p)$doublet_likelihood
#   colData(thisCDS)$predicted.id = colData(cdsList[["Activity_CDS"]])$predicted.id

#   # Subset?
#   if (!(opt$sampleRNAname == "All_Cells")){
#     # thisCDS = thisCDS[,colData(thisCDS)$sampleName == opt$sampleRNAname]
#     # miniCDS_g = cds.a.g[,colData(cds.a.g)$sampleName == opt$sampleRNAname]
#     thisCDS = thisCDS[,colData(thisCDS)$sampleName == samplesATACnames[opt$sampleRNAname]]
#     miniCDS_g = cds.a.g[,colData(cds.a.g)$sampleName == samplesATACnames[opt$sampleRNAname]]
#   } else {
#     miniCDS_g = cds.a.g
#   }

#   # Subset by FRIP
#   miniCDS_g$FRIP = thisCDS$FRIP
#   miniCDS_g = miniCDS_g[,colData(miniCDS_g)$FRIP > opt$fripMin]
#   thisCDS = thisCDS[,colData(thisCDS)$FRIP > opt$fripMin]

#   print(paste0("Working on ", eachCDSname))
#   print(str(colData(thisCDS)))

#   thisCDS = estimate_size_factors(thisCDS)
#   thisCDS = preprocess_cds(thisCDS, method=preprocessMethods[eachCDSname], num_dim=opt$pcToUse)

#   # MNN?
#   if (opt$useMNN == "useMNN"){
#     thisCDS = align_cds(thisCDS, preprocess_method = preprocessMethods[eachCDSname], alignment_group = "Sample", 
#                                          residual_model_formula_str = ~log10umi + FRIP + doublet_likelihood)
#     thisCDS = reduce_dimension(thisCDS, reduction_method = 'UMAP', preprocess_method = "Aligned")
#   } else{
#     thisCDS = reduce_dimension(thisCDS, reduction_method = 'UMAP', preprocess_method = preprocessMethods[eachCDSname])
#   }
#   # Reduce dimension

#   # Cluster
#   thisCDS = cluster_cells(thisCDS, k=opt$kVal)
#   colData(thisCDS)$cluster_label = as.character(clusters(thisCDS))
#   plotUMAP_Monocle(thisCDS, paste0(eachCDSname, "comp", as.character(opt$pcToUse), "k", as.character(opt$kVal), opt$ATACprocNote), 
#   						"cluster_label", outputPath=out_dir)

#   # See if there's a reason to save this cds
#   if (eachCDSname == opt$cdsToSave){
#   	# Assign clusters?
#   	thisCDS = hardAssignATAC_clusters(thisCDS, paste0(opt$ATACprocNote, "_", opt$useMNN, "_", eachCDSname, "_", as.character(opt$pcToUse)),
#   								opt$kVal)
#   	# Save it
#   	outputFile = paste0(out_dir, paste0(opt$ATACprocNote, "_", opt$useMNN, "_", eachCDSname, "_", as.character(opt$pcToUse), "k_", as.character(opt$kVal)), ".rds")
#   	saveRDS(thisCDS, file=outputFile)
#   	print(outputFile)

#   	# Get DE testing by cell type
#   	plotLabel = paste0(opt$ATACprocNote, "_", opt$useMNN, "_", eachCDSname, "_", as.character(opt$pcToUse), "k_", as.character(opt$kVal))
#   	cdsAndDEtest = runDEtestingToID_markers(thisCDS, plotLabel, "Assigned_Cell_Type",
#                   howManyGenesToTest = 50, outputPath=out_dir)
#   	DEtestResFile = paste0(out_dir,  paste0(opt$ATACprocNote, "_", opt$useMNN, "_", eachCDSname, "_", as.character(opt$pcToUse), "k_", as.character(opt$kVal)), "_DE.csv")
# 	write.csv(cdsAndDEtest$marker_test_res[order(cdsAndDEtest$marker_test_res$cell_group),], DEtestResFile)
#   }


# # "FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_peakMat_LSI1_50"

#   # Give these over to the g matrix
#   thisProcNote = paste0("UMAP_by_", eachCDSname, "_", opt$useMNN, "FRIP_", as.character(opt$fripMin), "comp", as.character(opt$pcToUse))

#   reducedDims(miniCDS_g)$UMAP = reducedDims(thisCDS)$UMAP
#   plotUMAP_MonocleModded(miniCDS_g, paste0(thisProcNote, "_faceted"), "sampleName",
#                         outputPath = out_dir, show_labels=FALSE)

#   plotUMAP_Monocle(miniCDS_g, thisProcNote, "sampleName",
#                         outputPath = out_dir, show_labels=FALSE)

#   # Now color by some marker genes
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("TTN", "RBPJ","GSN", "LDB2", "MYH11", "THEMIS" ),
#                         "generalMarkers",
#                         outputPath = out_dir)
#   # 
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("PECAM1", "CD106", "CD62E", "SELE", "KDR", "ENG" ),
#                         "coleEndothMark",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("MYH11", "SMTN", "CALD1", "CNN1", "CNN2" ),
#                         "coleSmoothMuscMark",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("CD68", "CD14", "FCGR1", "MS4A7", "FCER1A", "CD163", "LY6C1", "FCN1", "MERTK" ),
#                         "coleMonoMacMark",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("CD3E", "CD3G", 'CD3D", "CD4", "CD8A', "CD8B", "BCL11B" ),
#                         "coleTCellMark",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("CD90", "THY1", "VIM", "DES" ),
#                         "coleFibroblastMark",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("CD90", "THY1", "VIM", "DES" ),
#                         "coleFibroblastMark",
#                         outputPath = out_dir)



#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("CD90", "THY1", "VIM", "DES" ),
#                         "coleFibroblastMark",
#                         outputPath = out_dir)

#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("PCDHGA1", "PCDHGA2" ),
#                         "protocadherins",
#                         outputPath = out_dir)

#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("F13A1", "RBM47",  "CD74", "RBPJ" ), # Skip RBPJ
#                         "macTopRNA",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("SKAP1", "CD247", "ARHGAP15", "RIPOR2" ), # Skip RBPJ
#                         "tcellTopRNA",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("VWF", "LDB2", "ANO2", "FLT1" ), # Skip RBPJ
#                         "vascEndTopRNA",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("DLC1", "RGS5", "NOTCH3", "EGFLAM" ), # Skip RBPJ
#                         "perivasTopRNA",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("NRXN1", "XKR4", "NRXN3", "CADM2" ), # Skip RBPJ
#                         "neuronTopRNA",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("DCN", "GSN", "MGP", "C7" ), # Skip RBPJ
#                         "fibroTopRNA",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("PKHD1L1", "PCDH15", "LDB2", "NPR3" ), # Skip RBPJ
#                         "endocardTopRNA",
#                         outputPath = out_dir)

#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("EMCN", "POSTN", "VWF" ), # Skip RBPJ
#                         "endocardTopRNA_2",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("AQP1", "ADAMTSL1", "SMOC1" ), # Skip RBPJ
#                         "endocardTopRNA_3",
#                         outputPath = out_dir)

#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("TTN", "MYBPC3", "DB3", "RBM20" ), # Skip RBPJ
#                         "cardiomyocyteTopRNA",
#                         outputPath = out_dir)
#   plotUMAP_Monocle_genes(miniCDS_g, thisProcNote, c("BANK1", "PAX5", "CD79A", "RALGPS2" ), # Skip RBPJ
#                         "bCellTopRNA",
#                         outputPath = out_dir)
# }



