
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/"
dir.create(paste0(basepath, "archr/results/finalPlots/"))
out_dir = paste0(basepath, "archr/results/finalPlots/")
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
library(RColorBrewer)

library(stringr)


source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

library("optparse")

# Get the passed parameters
option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
  			# default="PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5",
        default="PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5",  
              help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-s", "--sampleRNAname"), type="character",  default= "All_Cells", #  "All_Cells",  # default="W142.Left.Vent",  "W144.Apex"
              help="Sample to integrate", metavar="character"),
    make_option(c("-u", "--useSubset"), type="character", default="Subset", 
              help="Subset or All", metavar="character"),

    make_option(c("-m", "--useMNN"), type="character", default="useMNN", 
              help="Subset or All", metavar="character"),

    # make_option(c("-t", "--TTS_Scores"), type="logical", default=FALSE,
    #                 help="Set to true to over-ride and read in a gene score matrix from TSS-based read counts, instead of gene activity scores", metavar="logical"),

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
              help="peakMat, bMat, or gMat", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)



# outputNote = paste0("Harmony_Aligned_CoordsRegress_Protocadherin_", opt$processingNote)
# inputDir = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/archr/results/NB11/All_Cells/"
outputNote = opt$processingNote
inputDir = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/cds_objects/"
kVal = opt$kVal

# rdsDir = 
dir.create("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/plots/")
outDir = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/plots/QC_Outputs_Post_Filter/"
dir.create(outDir)

copyBBIoutputs = function(outDir){

  bbiQCDir = paste0(outDir, "Post_Initial_Pipeline_QC/")
  dir.create(bbiQCDir)

  dirList = list.dirs("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out", recursive=FALSE)
  dirList = dirList[str_detect(dirList, "analyze_out/W1")]

  for (eachDir in dirList){
    thisSample = sub('.*\\/', '', eachDir)

    myCommand = paste0("cp ", eachDir, "/summarize_cell_calls/", thisSample, "-called_cells_summary.stats.txt ", bbiQCDir)
    system(myCommand)
    # print(myCommand)
  }

}


plotUMAP_Monocle_paper_format <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, textSize=10, outputPath = "./plots/",
                    returnPlotObj=FALSE){ #, xVal, yVal){
    png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored.png"),
             width=2100, height=1500, res=300)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels,
          cell_stroke=.1 , group_label_size=textSize        
                ))
    myPlot = (myPlot + theme(text=element_text(size=textSize)))
    print(myPlot)
    dev.off()   

    if (returnPlotObj){
        return(myPlot)
    }
}


boxplot_stat_by_X_paper_format <- function(inputCDS, processingNote, statToPlot, statToSplit,
                xLabToUse=statToSplit, yLabToUse=statToPlot, outputPath = "./plots/"){
    # See doublet numbers (Scrublet Estimates)
    ggplot(as.data.frame(colData(inputCDS)), aes_string(x=statToSplit, y=statToPlot)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        xlab(xLabToUse) + ylab(yLabToUse)
        # facet_grid(. ~ statToPlot)
    ggsave(paste0(outputPath, processingNote, statToPlot,  "_by_", statToSplit, ".png"))
}


makeATAC_QC = function(inputCDS, opt, outDir, metricsToShow){

  # Now get extra QC plots from the harmony cells
  for (eachMetric in metricsToShow){
    plotUMAP_Monocle_paper_format(inputCDS, "ATAC_CDS_", eachMetric,
          show_labels = FALSE, outputPath = outDir
          )

    boxplot_stat_by_X_paper_format(inputCDS, "ATAC_CDS_", eachMetric, "Sample",
                                  outputPath = outDir)

  }

}

copyBBIoutputs(outDir)

harmonyFile = paste0(inputDir, outputNote, "_ArchRActivharmonyAligned_protocad_reg_cds.rds")
print(harmonyFile)
harmonyRNA_and_ATAC = readRDS(harmonyFile)

atacCDS = harmonyRNA_and_ATAC[,colData(harmonyRNA_and_ATAC)$tech == "ATAC"]

colnames(atacCDS) = sub("_[^_]*$", "", colnames(atacCDS))



makeATAC_QC(atacCDS, opt, outDir, c( "TSSEnrichment", "nFrags"))




cdsPfile = paste0(inputDir, outputNote, "_ArchRActivharmonyLabels_regressProtocad_cds_p.rds")
cds_p = readRDS(cdsPfile)

cds_p = cds_p[,colnames(atacCDS)]
reducedDims(cds_p)$UMAP = reducedDims(atacCDS)$UMAP


colData(cds_p)$Sample = colData(cds_p)$sampleName

makeATAC_QC(cds_p, opt, outDir, c("doublet_likelihood", "FRIT", "FRIP"))



