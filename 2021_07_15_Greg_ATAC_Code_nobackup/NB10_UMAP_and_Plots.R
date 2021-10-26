
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
dir.create(paste0(basepath, "archr/results/NB10/"))
out_dir = paste0(basepath, "archr/results/NB10/")
setwd(out_dir)
set.seed(7)

# load requirements
suppressPackageStartupMessages({
  # library(ArchR)
  library(monocle3)
  library(dplyr)
  library(ggrastr)
  library(Seurat)
})
library(tidyr)
library(plyr)
library(ggplot2)


source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

library("optparse")

# Get the passed parameters
option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
  			default="HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes", 
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-s", "--sampleRNAname"), type="character",  default="All_Cells",  # default="W142.Left.Vent", 
              help="Sample to integrate", metavar="character"),
    make_option(c("-u", "--useSubset"), type="character", default="Subset", 
              help="Subset or All", metavar="character"),

    make_option(c("-m", "--useMNN"), type="character", default="useMNN", 
              help="Subset or All", metavar="character"),

    # make_option(c("-c", "--fripMin"), type="numeric", default=0.2, 
    #           help="Min FRIP value to permit", metavar="numeric"),

    make_option(c("-o", "--lowDimOne"), type="numeric", default=1, 
              help="Coordinate to start in LSI", metavar="numeric"),

    make_option(c("-t", "--lowDimTwo"), type="numeric", default=50, 
              help="Highest coord in LSI", metavar="numeric"),


    make_option(c("-a", "--ATACprocNote"), type="character", default="FRIP=0.2_FRIT=0.05UMI=1000DL=0.5", 
              help="How ATAC Cells were filtered", metavar="character")
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

# Get the ATAC cds
atacProcNote =  opt$ATACprocNote # "FRIP=0.2_FRIT=0.05UMI=1000"


# }
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


# Get the saved CDS
rdsPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/rdsOutput/"
# fileName = paste0(opt$sampleRNAname, "_", opt$ATACprocNote, "gMatrixCDS_postTransfer.RDS")
fileName = paste0("All_Cells", "_", opt$ATACprocNote, "gMatrixCDS_postTransfer.RDS")

# Read in the data
cds.a.g = readRDS(paste0(rdsPath, fileName))
# bMatFile = paste0(opt$sampleRNAname, "_", opt$ATACprocNote, "bMatrixCDS_postTransfer.RDS")
bMatFile = paste0("All_Cells", "_", opt$ATACprocNote, "bMatrixCDS_postTransfer.RDS")
cds.a.b = readRDS(paste0(rdsPath, bMatFile))


# Get data from Greg's filtering process

# Added 7-29-21: Filter by checking names in a cds filtered by greg/riza's method
filterCDSName = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/cds_objects/"
# processingNote = "FRIP=0.2_FRIT=0.05UMI=1000"
# processingNote = "FRIP=0.3_FRIT=0.1UMI=1000"
processingNote = opt$ATACprocNote


oldProcNote = "FRIP=0.2_FRIT=0.05UMI=1000"
# cds_p holds data
load(paste0(filterCDSName, "cds_p_allHeartATAC", oldProcNote))

# Get the sampleName reformated
colDataDF = separate(data=as.data.frame(colData(cds_p)), col = "sampleName", 
                        into = c("sampleName", "Processing"), sep="FRIP", remove=TRUE)
cds_p$sampleName = colDataDF$sampleName

# Reformat the colnames of cds_p (original/peak matrix) to match the others
colDataDF$rowname = rownames(colDataDF)
colDataDF = separate(data=colDataDF, col = "rowname", 
                        into = c("Part1", "Part2", "Part3"), sep="_", remove=TRUE)
rownames(colDataDF) = paste0(colDataDF$Part1, "_", colDataDF$Part2, "_", colDataDF$Part3)


rownames(colData(cds_p)) = paste0(colData(cds_p)$sampleName, "#", rownames((colDataDF)))
colnames(cds_p) = rownames(colData(cds_p))

# Now get the portion of cds_p we want, in order
cds_p = cds_p[,colnames(cds.a.g)]
colData(cds.a.g)$sampleName = colData(cds_p)$sampleName
colData(cds.a.b)$sampleName = colData(cds_p)$sampleName
# Get data from 

cdsList = list(cds_p, cds.a.g, cds.a.b)
names(cdsList) = c("Peak_CDS", "Activity_CDS", "Bin_CDS")

preprocessMethods = c("LSI", "PCA", "LSI")
names(preprocessMethods) = c("Peak_CDS", "Activity_CDS", "Bin_CDS")


# opt$sampleRNAname = "W144.heart.apex.s1"
# opt$useMNN = "NoMNN"

out_dir = paste0(basepath, "archr/results/NB10/", opt$sampleRNAname, "/")
dir.create(out_dir)
setwd(out_dir)


# Get a UMAP and call clusters off of the reduced dimension
thisCDS = cds_p
thisCDS = estimate_size_factors(thisCDS)
thisCDS = preprocess_cds(thisCDS, method="LSI")


colData(thisCDS)$log10umi = log10(colData(cds_p)$umi)

# MNN?
if (opt$useMNN == "useMNN"){
  thisCDS = align_cds(thisCDS, preprocess_method = "LSI", alignment_group = "sampleName", 
                                       residual_model_formula_str = ~log10umi + FRIP + doublet_likelihood)
  thisCDS = reduce_dimension(thisCDS, reduction_method = 'UMAP', preprocess_method = "Aligned")
} else{
  thisCDS = reduce_dimension(thisCDS, reduction_method = 'UMAP', preprocess_method = "LSI")
}





hardAssignATAC_clusters <- function(inputCDS, processingNote, kValUsed){

  #
  if ((processingNote == "FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_peakMat_LSI1_50") & (kValUsed == 20)){
    colData(inputCDS)$Assigned_Cell_Type = "Ambiguous"

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

  if ((processingNote == "FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_peakMat_LSI1_50") & (kValUsed == 15)){
    colData(inputCDS)$Assigned_Cell_Type = "Ambiguous"

    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(32, 33), 
                    "Endocardium", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(18, 7), 
                    "Macrophage", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(26,27), 
                    "Perivascular", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(34, 3, 2, 13), 
                    "Fibroblast", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(30,17,28,5,15,11,8,20,16), 
                    "Vascular_Endothelium", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(1,22,29,14,6,10,4,25,23,12,31,19,24,9,21), 
                    "Cardiomyocyte", colData(inputCDS)$Assigned_Cell_Type)
  }


  return(inputCDS)
}



# Cluster
kValUsed = 10
thisCDS = cluster_cells(thisCDS, k = kValUsed)


colData(thisCDS)$cluster_label = clusters(thisCDS)
colData(thisCDS)$cluster_label = as.character(colData(thisCDS)$cluster_label)
colData(thisCDS)$partition_label = as.character(partitions(thisCDS))




plotUMAP_Monocle(thisCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kValUsed), 
            "cluster_label", outputPath=out_dir, show_labels=TRUE)
plotUMAP_Monocle(thisCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kValUsed), 
            "partition_label", outputPath=out_dir, show_labels=TRUE)

processingNote = paste0(processingNote, "_", opt$useMNN, "_peakMat_LSI1_50")

thisCDS = hardAssignATAC_clusters(thisCDS, processingNote, kValUsed)




plotUMAP_Monocle_modded <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, textSize=12, outputPath = "./plots/",
                    returnPlotObj=FALSE){ #, xVal, yVal){
    png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored.png"),
             width=1400, height=1000, res=200)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels,# group_cells_by="partition",
          cell_stroke=.1 , group_label_size=(textSize * .5)        
                ))
    myPlot = (myPlot + theme(text=element_text(size=textSize))) +
            scale_color_discrete(name = "Cell Type", labels = c("Cardiomyocyte", "Endocardium", "Fibroblast", 
                                                    "Macrophage", "Perivascular", "Vascular Endothelium"))+
            scale_fill_discrete(name = "Cell Type", labels = c("Cardiomyocyte", "Endocardium", "Fibroblast", 
                                                    "Macrophage", "Perivascular", "Vascular Endothelium"))+ 
            guides(fill=guide_legend(title="Cell Type"))
    print(myPlot)
    dev.off()   

    if (returnPlotObj){
        return(myPlot)
    }

}



plotUMAP_Monocle_modded(thisCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kValUsed), 
            "Assigned_Cell_Type", outputPath=out_dir, show_labels=FALSE, textSize=24)




















