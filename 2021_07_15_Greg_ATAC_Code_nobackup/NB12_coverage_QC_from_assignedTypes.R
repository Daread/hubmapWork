
# Trying:
# qlogin -l mfree=20G -q trapnell-login.q -pe serial 16 -l centos=7

# Trying this on 7-20-21. (Didn't work, going back to 20g)
# qlogin -l mfree=40G -q trapnell-login.q -pe serial 16 -l centos=7

basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
dir.create(paste0(basepath, "archr/results/NB12/"))
out_dir = paste0(basepath, "archr/results/NB12/")
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
addArchRThreads(threads = 1) 
# addArchRThreads(threads = 16) 


# Parameters to get the co-embedded/annotated ATAC profiles
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

library("optparse")

library(rtracklayer)

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

    make_option(c("-a", "--ATACprocNote"), type="character", default="FRIP=0.2_FRIT=0.05UMI=1000DL=0.5", 
              help="How ATAC Cells were filtered", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
print(opt)

outputNote = paste0("Harmony_Aligned_CoordsRegress_Protocadherin_", opt$processingNote)
inputDir = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/archr/results/NB11/All_Cells/"
kVal = opt$kVal

# rdsDir = 
harmonyFile = paste0(inputDir, paste0(outputNote, "_k_", as.character(kVal), "_RNA_and_ATAC_cds.rds"))
print(harmonyFile)
harmonyRNA_and_ATAC = readRDS(harmonyFile)

atacCDS = harmonyRNA_and_ATAC[,harmonyRNA_and_ATAC$tech == "ATAC"]

procNoteFull = paste0(outputNote, "_k_", as.character(kVal))

# Format the cell names of 
archrFormatNames = sub('_[^_]*$', '', rownames(colData(atacCDS)))



# load filtered ArchR project
prj = loadArchRProject(path = paste0(basepath, "archr/Heart_filtered"),
                       showLogo = FALSE)



backupPrj = prj

cellsToKeepIDs = BiocGenerics::which(prj$cellNames %in% archrFormatNames)
archrCellsToKeep = prj$cellNames[cellsToKeepIDs]
prj = prj[archrCellsToKeep,]

# Format the ATAC data
colnames(atacCDS) = archrFormatNames
atacCDS = atacCDS[,prj$cellNames]

# Transfer labels over to the archr project
prj = addCellColData(ArchRProj = prj, data=as.character(colData(atacCDS)$harmonyKNN_type),
                    cells=prj$cellNames, name="harmonyKNN_type")






# load filtered ArchR project
copiedPrj = loadArchRProject(path = paste0(basepath, "archr"),
                       showLogo = FALSE)


copiedPrj = addCellColData(ArchRProj = copiedPrj, data= rep("1", length(copiedPrj$cellNames)),
                    cells=copiedPrj$cellNames, name="QC_Success")


png(paste0(out_dir, summaryToUse, "_EarlyArrow_UnfilteredByCells_TTN_single_groupobrowser.png"), res=200, height=1200, width=1000)
myPlot = plotBrowserTrack(ArchRProj = copiedPrj, groupBy="QC_Success",
                  plotSummary=c("bulkTrack", "scTrack", "geneTrack"),
                  scCellsMax=10000,
                      geneSymbol= c("TTN", "RYR2"), upstream=400000, downstream=400000)
print(grid::grid.draw(myPlot$TTN))
dev.off()





# Plot based on set coordiantes

# cdD_and_G_IVG_Range = GRanges(seqnames = "chr11", strand=c("+"), 
#                       ranges = IRanges(start=c(118212000), width=5000))

# png(paste0(out_dir, "CD3D_and_G_IGV_Coords_browser.png"), res=200, height=1200, width=1000)
# rangePlot = plotBrowserTrack(ArchRProj = prj, groupBy="harmonyKNN_type",
#                   plotSummary=c("bulkTrack", "scTrack", "geneTrack"),
#                   scCellsMax=10000,
#                       #geneSymbol= c("CD3D", "CD3E", "CD3G"), upstream=20000, downstream=20000)
#                       region= cdD_and_G_IVG_Range)
# print(grid::grid.draw(rangePlot))
# dev.off()


# TTN_IVG_Range = GRanges(seqnames = "chr2", strand=c("+"), 
#                       ranges = IRanges(start=c(179670000), width=10000))

# png(paste0(out_dir, "TTN_IGV_Coords_browser.png"), res=200, height=1200, width=1000)
# rangePlot = plotBrowserTrack(ArchRProj = prj, groupBy="harmonyKNN_type",
#                   plotSummary=c("bulkTrack", "scTrack", "geneTrack"),
#                   scCellsMax=10000,
#                       #geneSymbol= c("CD3D", "CD3E", "CD3G"), upstream=20000, downstream=20000)
#                       region= TTN_IVG_Range)
# grid::grid.newpage()
# print(grid::grid.draw(rangePlot))
# dev.off()


# png(paste0(out_dir, "CD3D_and_G_IGV_Coords_browser.png"), res=200, height=1200, width=1000)
# print(grid::grid.draw(rangePlot))
# dev.off()










# Plot browser tracks
png(paste0(out_dir, "CD3D_browser.png"), res=200, height=1200, width=1000)
myPlot = plotBrowserTrack(ArchRProj = prj, groupBy="harmonyKNN_type",
                  plotSummary=c("bulkTrack", "scTrack", "geneTrack"),
                  scCellsMax=10000,
                      geneSymbol= c("CD3D", "CD3E", "CD3G"), upstream=20000, downstream=20000)

print(grid::grid.draw(myPlot$CD3D))
dev.off()

png(paste0(out_dir, "CD3E_browser.png"), res=200, height=1200, width=1000)
print(grid::grid.draw(myPlot$CD3E))
dev.off()


png(paste0(out_dir, "CD3G_browser.png"), res=200, height=1200, width=1000)
print(grid::grid.draw(myPlot$CD3G))
dev.off()


summaryToUse = "bulkTrack" #  # "scTrack" not working

png(paste0(out_dir, summaryToUse, "_TTN_browser.png"), res=200, height=1200, width=1000)
myPlot = plotBrowserTrack(ArchRProj = prj, groupBy="harmonyKNN_type",
                  plotSummary=c("bulkTrack", "scTrack", "geneTrack"),
                  scCellsMax=10000,
                    #  scCellsMax = 1000,
                      geneSymbol= c("TTN", "RYR2"), upstream=500000, downstream=500000)
print(grid::grid.draw(myPlot$TTN))
dev.off()

png(paste0(out_dir, summaryToUse, "_RYR2_browser.png"), res=200, height=1200, width=1000)
print(grid::grid.draw(myPlot$RYR2))
dev.off()



# Test: Run grouping into a single "group"
colData(atacCDS)$QC_Success = as.character(colData(atacCDS)$PassQC)
prj = addCellColData(ArchRProj = prj, data=as.character(colData(atacCDS)$QC_Success),
                    cells=prj$cellNames, name="QC_Success")


png(paste0(out_dir, summaryToUse, "_TTN_single_groupobrowser.png"), res=200, height=1200, width=1000)
myPlot = plotBrowserTrack(ArchRProj = prj, groupBy="QC_Success",
                  plotSummary=c("bulkTrack", "scTrack", "geneTrack"),
                  scCellsMax=100000,
                      geneSymbol= c("TTN", "RYR2"), upstream=400000, downstream=400000)
print(grid::grid.draw(myPlot$TTN))
dev.off()










backupPrj = addCellColData(ArchRProj = backupPrj, data= rep("1", length(backupPrj$cellNames)),
                    cells=backupPrj$cellNames, name="QC_Success")


png(paste0(out_dir, summaryToUse, "_UnfilteredByCells_TTN_single_groupobrowser.png"), res=200, height=1200, width=1000)
myPlot = plotBrowserTrack(ArchRProj = backupPrj, groupBy="QC_Success",
                  plotSummary=c("bulkTrack", "scTrack", "geneTrack"),
                  scCellsMax=10000,
                      geneSymbol= c("TTN", "RYR2"), upstream=400000, downstream=400000)
print(grid::grid.draw(myPlot$TTN))
dev.off()





# Make bulk tracks






# Get some test information

# load PeakMatrix into memory
bMat = getMatrixFromProject(
  ArchRProj = prj, 
  useMatrix = "TileMatrix", 
  binarize = TRUE, 
  threads = getArchRThreads())


# Make a fragment distribution
png(paste0(out_dir, "bmatFragDistribution.png"))
myPlot = ggplot(as.data.frame(colData(bMat)), aes_string(x="nFrags")) + 
        geom_histogram()

print(myPlot)
dev.off()



TTN_RANGE = GRanges(seqnames = "chr2", strand=c("+"), 
                      ranges = IRanges(start=c(178760000), width=70000))
filteredFragments = getFragmentsFromProject(ArchRProj=prj, subsetBy=TTN_RANGE)
fragsList = unlist(as(filteredFragments, "GRangesList"))

# fragCoverage = coverage(fragsList)


export.bed(fragsList, con=paste0(out_dir, "TTN_filtered_Cell_reads.bed"))

bedgraphFrags = new("GraphTrackLine", type="bedGraph", name="MySample")
export(fragsList, trackLine=bedgraphFrags, format="bedGraph", con=paste0(out_dir, "TTN_filtered_Cell_reads.bedGraph"))

# Loop, get names by different cell types
typesToMakeBed = c("Cardiomyocyte", "Macrophage", "T_Cell", "Vascular_Endothelium", "Fibroblast")
for (eachType in typesToMakeBed){
  cellnamesToSave = colnames(atacCDS[,colData(atacCDS)$harmonyKNN_type == eachType])
  filteredFragments = getFragmentsFromProject(ArchRProj=prj, subsetBy=TTN_RANGE,
                    cellNames = cellnamesToSave)
  fragsList = unlist(as(filteredFragments, "GRangesList"))
  export.bed(fragsList, con=paste0(out_dir, eachType, "_TTN_filtered_Cell_reads.bed"))

  # Run the script to generate genome coverage bedGraphs
  bedFilePrefix=paste0(out_dir, eachType, "_TTN_filtered_Cell_reads")
  commandToRun = paste0("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/archr/results/NB12/convertBedfiles.sh ", bedFilePrefix)
  print(commandToRun)
  system(commandToRun)
}






CD3_RANGE = GRanges(seqnames = "chr11", strand=c("+"), 
                      ranges = IRanges(start=c(118330000), width=50000))
filteredFragments = getFragmentsFromProject(ArchRProj=prj, subsetBy=CD3_RANGE)
fragsList = unlist(as(filteredFragments, "GRangesList"))

# fragCoverage = coverage(fragsList)


export.bed(fragsList, con=paste0(out_dir, "CD3_filtered_Cell_reads.bed"))

bedFilePrefix=paste0(out_dir, "CD3_filtered_Cell_reads")
commandToRun = paste0("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/archr/results/NB12/convertBedfiles.sh ", bedFilePrefix)
print(commandToRun)
system(commandToRun)


# Loop, get names by different cell types
typesToMakeBed = c("Cardiomyocyte", "Macrophage", "T_Cell", "Vascular_Endothelium", "Fibroblast")
for (eachType in typesToMakeBed){
  cellnamesToSave = colnames(atacCDS[,colData(atacCDS)$harmonyKNN_type == eachType])
  filteredFragments = getFragmentsFromProject(ArchRProj=prj, subsetBy=CD3_RANGE,
                    cellNames = cellnamesToSave)
  fragsList = unlist(as(filteredFragments, "GRangesList"))
  export.bed(fragsList, con=paste0(out_dir, eachType, "_CD3_filtered_Cell_reads.bed"))

  # Run the script to generate genome coverage bedGraphs
  bedFilePrefix=paste0(out_dir, eachType, "_CD3_filtered_Cell_reads")
  commandToRun = paste0("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/archr/results/NB12/convertBedfiles.sh ", bedFilePrefix)
  print(commandToRun)
  system(commandToRun)
}





