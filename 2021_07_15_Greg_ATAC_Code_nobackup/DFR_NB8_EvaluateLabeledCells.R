
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
dir.create(paste0(basepath, "archr/results/NB8/"))
out_dir = paste0(basepath, "archr/results/NB8/")
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
library(tidyr)


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


# Get the saved CDS
rdsPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/rdsOutput/"
fileName = paste0(opt$sampleRNAname, "_", opt$ATACprocNote, "gMatrixCDS_postTransfer.RDS")

# Read in the data
cds.a.g = readRDS(paste0(rdsPath, fileName))

# Get data from Greg's filtering process

# Added 7-29-21: Filter by checking names in a cds filtered by greg/riza's method
filterCDSName = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/cds_objects/"
# processingNote = "FRIP=0.2_FRIT=0.05UMI=1000"
# processingNote = "FRIP=0.3_FRIT=0.1UMI=1000"
processingNote = opt$ATACprocNote


oldProcNote = "FRIP=0.2_FRIT=0.05UMI=1000"
# cds_p holds data
load(paste0(filterCDSName, "cds_p_allHeartATAC", oldProcNote))

qcInOrigFiltering = c("FRIP", "FRIT", "umi", "doublet_score", "doublet", "doublet_likelihood")
# 
dfTransferredColdata = as.data.frame(colData(cds.a.g))
dfOrigFilteredColdata = as.data.frame(colData(cds_p))

if (opt$sampleRNAname != "All_Cells"){
  dfOrigFilteredColdata = dfOrigFilteredColdata[c(qcInOrigFiltering, "cell")]
  dfOrigFilteredColdata$cell = paste0(samplesATACnames[opt$sampleRNAname], "#", dfOrigFilteredColdata$cell)
} else {
  dfOrigFilteredColdata = separate(data=dfOrigFilteredColdata, col = "sampleName", 
                        into = c("sampleName", "Processing"), sep="FRIP", remove=TRUE)
  dfOrigFilteredColdata = dfOrigFilteredColdata[c(qcInOrigFiltering, "cell", "sampleName")]
  dfOrigFilteredColdata$cell = paste0(dfOrigFilteredColdata$sampleName, "#", dfOrigFilteredColdata$cell)
}


dfTransferredColdata$cell = rownames(dfTransferredColdata)

library(plyr)


# mergedDF = merge(dfTransferredColdata, dfOrigFilteredColdata, by="cell")
mergedDF = join(dfTransferredColdata, dfOrigFilteredColdata, by="cell")
rownames(mergedDF) = rownames(dfTransferredColdata)

# Put into a new cds
annotatedCDS = new_cell_data_set(expression_data = exprs(cds.a.g), 
                            cell_metadata = mergedDF, 
                            gene_metadata = as.data.frame(rowData(cds.a.g)))

reducedDims(annotatedCDS)$UMAP = reducedDims(cds.a.g)$UMAP

# Plot doublet scores since I wonder if this is important in blurring cluster borders
# 8-6-21 addition of this plot
plotDir = out_dir
thisOutput = paste0(out_dir, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),  "Doublet_Likelihood_Hist.png")
png(thisOutput, width=1000, height=1000, res = 200)
myPlot = ggplot(as.data.frame(colData(annotatedCDS)) , aes_string("doublet_likelihood")) + 
      geom_histogram()
print(myPlot)
dev.off()


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

plotUMAP_MonocleModded(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                  "predicted.id", outputPath="./", show_labels=FALSE)


#shuffle
# annotatedCDS = annotatedCDS[sample(1:nrow(annotatedCDS)), sample(1:ncol(annotatedCDS))]


colData(annotatedCDS)$log10UMI = log10(colData(annotatedCDS)$umi)

colData(annotatedCDS)$monoFragProp = colData(annotatedCDS)$nMonoFrags / colData(annotatedCDS)$nFrags 
colData(annotatedCDS)$diFragProp = colData(annotatedCDS)$nDiFrags / colData(annotatedCDS)$nFrags 


for (eachCol in c("sampleName", qcInOrigFiltering, "log10UMI", "NucleosomeRatio", "TSSEnrichment", "PromoterRatio", "BlacklistRatio", 
            "monoFragProp", "diFragProp")){

  plotUMAP_MonocleModded(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                  eachCol, outputPath="./", show_labels=FALSE)
}


plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("TTN", "FHL2"), "CardiomyoctyeMarks", outputPath="./")
plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("GSN", "CRISPLD2", "NEGR1"), "FibroblastMarks", outputPath="./")
plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("RBPJ", "LYVE1"), "MacrophageMarks", outputPath="./")
plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("LDB2", "SYNE1"), "EndothelialMarks", outputPath="./")
plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("MYH11"), "VSM_Pericyte_Marks", plotSetTotals=TRUE,  outputPath="./")

plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("THEMIS", "CD34", "CD247", "SKAP1", "CD3D", "CD3E"), "T_Cell", plotSetTotals=TRUE,  outputPath="./")

plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("DCN", "PDGFRA"), "FibroblastMarksLitvin", outputPath="./")
plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("RGS5", "ABCC9", "KCNJ8", "TAGLN", "MYH11", "ACTA2"), plotSetTotals=TRUE,  "VSM_Peri_MarksLitvin", outputPath="./")

plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("SMOC1", "NPR3" ), "EndocardMarks", outputPath="./")
plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("SEMA3G", "EFNB2", "DLL4", "NR2F2", "ACKR1", "RGCC", "CA4" ), "Vasc_Endothelial_Marks",
                      plotSetTotals=TRUE, outputPath="./")

plotUMAP_Monocle_genes(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    c("MYH11"), "VSM_Pericyte_Marks", plotSetTotals=TRUE,  outputPath="./")


kVal=10
annotatedCDS = cluster_cells(annotatedCDS, k=kVal)

colData(annotatedCDS)$cluster_label = clusters(annotatedCDS)
colData(annotatedCDS)$cluster_label = as.character(colData(annotatedCDS)$cluster_label)

plotUMAP_Monocle(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
            "cluster_label", outputPath="./", show_labels=TRUE)


# Single big plot showing sample source
# plotUMAP_MonocleModded(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
#                   eachCol, outputPath="./", show_labels=FALSE)


plotUMAP_Monocle(annotatedCDS, paste0("ATAC_With_Transfer_No_Facet_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
            "sampleName", outputPath="./", show_labels=FALSE)



hardAssignATAC_clusters <- function(inputCDS, processingNote, kValUsed){

  # 8-6-21: Using k = 15 I made a prelim assignment
  if ((processingNote == "ATAC_With_Transfer_All_Cells_FRIP=0.2_FRIT=0.05UMI=1000") & (kValUsed == 15)){
    colData(inputCDS)$Assigned_Cell_Type = "Ambiguous"

    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(13), 
                    "Endocardium", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(7, 30), 
                    "Macrophage", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(5), 
                    "VSM_and_Perictye", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(28,33,10,11), 
                    "Fibroblast", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(24,27,14,25,32,18), 
                    "Vascular_Endothelium", colData(inputCDS)$Assigned_Cell_Type)
    colData(inputCDS)$Assigned_Cell_Type = ifelse(colData(inputCDS)$cluster_label %in% c(1,20,21,3,31,2,26,9,15,19,4,34,23,6), 
                    "Cardiomyocyte", colData(inputCDS)$Assigned_Cell_Type)
  }

  return(inputCDS)
}


annotatedCDS = hardAssignATAC_clusters(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
                    kVal)


plotUMAP_Monocle(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
            "Assigned_Cell_Type", outputPath="./", show_labels=FALSE)


# Given this  set of annotations, do we see representation of all samples within the different clusters?


makeFacetedProportionPlot_withFill(annotatedCDS, 
    paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
      "sampleName" , "Assigned_Cell_Type", 
      subsetPropsToShow=c("Cardiomyocyte", "Endocardium", "Fibroblast", "Macrophage", "Vascular_Endothelium", "VSM_and_Perictye"),
             outputPath="./" )


# Assign site/donor info to the dataframe
sampleInfoDF = tidyr::separate(data=as.data.frame(colData(annotatedCDS)), col = "sampleName", 
                        into = c("Donor", "Organ", "Anatomical_Site", "Replicate"), sep="\\.", remove=FALSE)


colData(annotatedCDS)$Donor = sampleInfoDF$Donor
colData(annotatedCDS)$Anatomical_Site = sampleInfoDF$Anatomical_Site



makeFacetedProportionPlot_withFill(annotatedCDS, 
    paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote, "_k=", kVal), 
      "sampleName" , "Assigned_Cell_Type", fillCol =  "Anatomical_Site", colorCol =  "Anatomical_Site",
      subsetPropsToShow=c("Cardiomyocyte", "Endocardium", "Fibroblast", "Macrophage", "Vascular_Endothelium", "VSM_and_Perictye"),
             outputPath="./" )





plotUMAP_Monocle(annotatedCDS, paste0("ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote), 
            "doublet_likelihood", outputPath="./", show_labels=FALSE)





















# plotUMAP_MonocleModded(cds.a.g, paste0("Not_annotated_ATAC_With_Transfer_", opt$sampleRNAname, "_", opt$ATACprocNote),
#                   "predicted.id", outputPath="./", show_labels=FALSE)










