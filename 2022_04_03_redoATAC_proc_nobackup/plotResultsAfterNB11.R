
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


# outputNote = paste0("Harmony_Aligned_CoordsRegress_Protocadherin_", opt$processingNote)
# inputDir = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/archr/results/NB11/All_Cells/"
outputNote = opt$processingNote
inputDir = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/cds_objects/"
kVal = opt$kVal

# rdsDir = 

# PC_20_All_CellsFRIP=0.1_FRIT=0.08UMI=1000DL=0.5_ArchRActivharmonyAligned_protocad_reg_cds.rds
harmonyFile = paste0(inputDir, outputNote, "_ArchRActivharmonyAligned_protocad_reg_cds.rds")
print(harmonyFile)
harmonyRNA_and_ATAC = readRDS(harmonyFile)



# replaceATAC_with_tssActivityScores <- function(originalCDS, opt){
#   # Get the original ATAC data and the new cds with tss scores
#   atacCDS = originalCDS[,colData(originalCDS)$tech == "ATAC"]
#   colnames(atacCDS) = sub('_[^_]*$', '', colnames(atacCDS))
#   tssCDS = readRDS("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results//2022_03_10_rerun_ArchR_nobackup/archr/results/NB4/geneActivityByPromoterRegionCDS.rds")

#   # Subset and match orders
#   tssCDS = tssCDS[,colnames(tssCDS) %in% colnames(atacCDS)]
#   tssCDS = tssCDS[,colnames(atacCDS)]

#   # Make a new CDS
#   mergedTSScds = monocle3::new_cell_data_set(
#   exprs(tssCDS), 
#   cell_metadata = colData(atacCDS),
#   gene_metadata = rowData(tssCDS))

#   mergedTSScds = estimate_size_factors(mergedTSScds)

#   # Merge back in with RNA data
#   cdsList = list("ATAC" = mergedTSScds, "RNA" = originalCDS[,colData(originalCDS)$tech == "RNA"])
#   newATAC_and_RNA_cds = combine_cds(cdsList, sample_col_name="InputCDS")

#   return(newATAC_and_RNA_cds)
# }




# procNoteFull = paste0(outputNote, "_k_", as.character(kVal))
procNoteFull = outputNote

# if (opt$TTS_Scores){
#   backupHarmony = harmonyRNA_and_ATAC
#   print("Replacing gene scores with TSS based read intersections")
#   harmonyRNA_and_ATAC <- replaceATAC_with_tssActivityScores(harmonyRNA_and_ATAC, opt)
#   out_dir = paste0(basepath, "archr/results/finalPlots/TSS_Scores/")
#   dir.create(out_dir)
# }







plotUMAP_Monocle_formatted <- function(dataCDS, processingNote, catToColor,
                    show_labels=FALSE, textSize=10, outputPath = "./plots/",
                    returnPlotObj=FALSE){ #, xVal, yVal){

  colData(dataCDS)[[catToColor]] = gsub("[.]", " ", colData(dataCDS)[[catToColor]])
  # Shuffle data to not overplot with the last subset
  dataCDS = dataCDS[,sample(ncol(dataCDS))]

    png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored_formatted.png"),
             width=1400, height=1000, res=200)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels,
          cell_stroke=.1 , group_label_size=textSize        
                ))
    myPlot = (myPlot + theme(text=element_text(size=textSize)))
    if (catToColor == "Assay"){
      myPlot = myPlot + scale_color_brewer(palette="Dark2")
    }
    print(myPlot)
    dev.off()   

    if (returnPlotObj){
        return(myPlot)
    }

}



formatCellType <- function(inputColumn){

  inputColumn = ifelse(inputColumn == "T_Cell", "T Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "VSM_and_Pericyte", "Perivascular Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "Vascular_Endothelium", "Vascular Endothelium", inputColumn)
  inputColumn = ifelse(inputColumn == "Lymphatic_Endothelium", "Lymphatic Endothelium", inputColumn)
  inputColumn = ifelse(inputColumn == "Mast_Cell", "Mast Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "B_Cell", "B Cell", inputColumn)

  return(inputColumn)
}


# Add a "tech" column ready to plot
colData(harmonyRNA_and_ATAC)$Assay = ifelse(colData(harmonyRNA_and_ATAC)$tech == "RNA", "RNA", "ATAC")
harmonyRNA_and_ATAC = harmonyRNA_and_ATAC[,sample(ncol(harmonyRNA_and_ATAC))]

# Plot by tech
plotUMAP_Monocle_formatted(harmonyRNA_and_ATAC, procNoteFull, "Assay", outputPath=out_dir, textSize=24, show_labels=FALSE)


plotUMAP_Monocle_formatted_knn <- function(dataCDS, processingNote, catToColor,
                    show_labels=FALSE, textSize=10, outputPath = "./plots/",
                    returnPlotObj=FALSE){ #, xVal, yVal){

  colData(dataCDS)[[catToColor]] = gsub("[.]", " ", colData(dataCDS)[[catToColor]])
  # Shuffle data to not overplot with the last subset
  dataCDS = dataCDS[,sample(ncol(dataCDS))]

    png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored_formatted.png"),
             width=1600, height=1000, res=200)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels,
          cell_stroke=.1 , group_label_size=textSize        
                ))
    myPlot = (myPlot + theme(text=element_text(size=textSize)) + 
            guides(col=guide_legend(title="KNN Assignment"))) +
             guides(colour = guide_legend(override.aes = list(size=10)))
    print(myPlot)
    dev.off()   

    if (returnPlotObj){
        return(myPlot)
    }

}

colData(harmonyRNA_and_ATAC)$harmonyKNN_type = as.character(colData(harmonyRNA_and_ATAC)$harmonyKNN_type)
colData(harmonyRNA_and_ATAC)$harmonyKNN_type = formatCellType(colData(harmonyRNA_and_ATAC)$harmonyKNN_type)


# Plot by tech
plotUMAP_Monocle_formatted_knn(harmonyRNA_and_ATAC, procNoteFull, "harmonyKNN_type", outputPath=out_dir, textSize=24, show_labels=FALSE)








# Also, make plots showing marker gene expression/activity based on classification by a cell type assignment

makeMarkerDotplot <- function(inputCDS, techToPlot, markersToPlot, markerSetName, colForAssignments, opt, outputPath = "./plots/",
						typePlotOrder="Default") {

	#Keep the tech desired
	plotCDS = inputCDS[,colData(inputCDS)$tech == techToPlot]


	if (!(typePlotOrder == "Default")) {
		colData(plotCDS)[[colForAssignments]] <- factor(colData(plotCDS)[[colForAssignments]], levels = typePlotOrder)
	}

	# Make the dotplot
	thisDotplot = plot_genes_by_group(plotCDS, markersToPlot, 
										group_cells_by = colForAssignments, ordering_type="none")
	thisDotplot = thisDotplot + 
				xlab("Cell Type") +
    				theme(text = element_text(size=18))

	# Output it
	png(paste0(outputPath, "Markers_for_", techToPlot, "_", markerSetName, ".png"), res=200, height = 1200, width=1500) 
	print(thisDotplot)
	dev.off()
}


makeDotplotsFromDE <- function(inputDEres, outputDir, myCDS, opt, qvalMin = .05, sortPolicy="marker_score",
									 nToPlot=15, assayToUse="ATAC", cellTypeCol="harmonyKNN_type", outputNote=""){

	cellTypes = levels(as.factor(inputDEres$cell_group))
	# Loop, get DE results for top markers per cell type by the policy selected
	for (eachType in cellTypes){
		print(paste0("Working on ", eachType))
		miniRes = inputDEres[(inputDEres$cell_group == eachType) & (inputDEres$marker_test_q_value < qvalMin),]
		if (nrow(miniRes) > 0){
			# Sort
			miniRes = miniRes[order(-miniRes[[sortPolicy]]),]
			# Plot 
			theseMarkers = miniRes[["gene_id"]][1:nToPlot]
			makeMarkerDotplot(myCDS, assayToUse, theseMarkers, paste0(outputNote, "Markers_for_", eachType, "_by_", sortPolicy), 
								cellTypeCol, opt, outputPath = outputDir)
		}
	}

}

library(dplyr)

calculateAndOutputATAC_vs_RNA_Props <- function(inputCDS, out_dir, opt, atacCol="harmonyKNN_type", 
                                                rnaCol="highLevelCellType"){
  # library(dplyr)
  # First get a dataframe with types and props
  atacCDS = inputCDS[,colData(inputCDS)$tech == "ATAC"]
  atacProps = as.data.frame(colData(atacCDS)) %>% dplyr::group_by(get(atacCol)) %>% 
                    dplyr::summarise(count=n()) %>% dplyr::mutate(atacProp=count / sum(count))
  atacProps = as.data.frame(atacProps)
  colnames(atacProps) = c("Cell_Type", "ATAC_Count", "ATAC_Prop")
  atacProps = atacProps[c("Cell_Type", "ATAC_Prop")]
  # browser()             
  rnaCDS = inputCDS[,colData(inputCDS)$tech == "RNA"]
  rnaProp = as.data.frame(colData(rnaCDS)) %>% dplyr::group_by(get(rnaCol)) %>% 
                    dplyr::summarise(count=n()) %>% dplyr::mutate(rnaProp=count / sum(count))
  rnaProp = as.data.frame(rnaProp)
  colnames(rnaProp) = c("Cell_Type", "RNA_Count", "RNA_Prop")
  rnaProp = rnaProp[c("Cell_Type", "RNA_Prop")]

  cellProps = merge(atacProps, rnaProp, by="Cell_Type")

  # First, make a scatter plot of the comparison
  propCorr = cor(cellProps$RNA_Prop, cellProps$ATAC_Prop)
  print(paste0("Correlation is ", propCorr))
  propScatter = paste0(out_dir, "RNA_vs_ATAC_CellType_Proportions.png")
  png(propScatter, res=200, height=800, width=1000)
  myPlot = ggplot(cellProps, aes_string(x="RNA_Prop", y="ATAC_Prop")) + 
                geom_point() + 
                theme(text=element_text(size=18)) + 
                xlab("Cell Prop in RNA") + ylab("Cell Prop in ATAC") 
  print( myPlot)
  dev.off()

  # print(cellProps)
  # Also output a csv of these proportions
  outputCSV = paste0(out_dir, "RNA_vs_ATAC_CellType_Proportions.csv")
  write.csv(cellProps, file=outputCSV)

}


colData(harmonyRNA_and_ATAC)$highLevelCellType = as.character(colData(harmonyRNA_and_ATAC)$highLevelCellType)
colData(harmonyRNA_and_ATAC)$highLevelCellType = formatCellType(colData(harmonyRNA_and_ATAC)$highLevelCellType)


calculateAndOutputATAC_vs_RNA_Props(harmonyRNA_and_ATAC, out_dir, opt)

cellTypePlotOrder = c("Cardiomyocyte", "Macrophage", "Fibroblast", "Vascular Endothelium", "Lymphatic Endothelium", "Endocardium", "Perivascular Cell",
						"T Cell", "B Cell", "Mast Cell", "Neuronal", "Adipocytes")


rnaGeneralMarkers = c("TTN", "RBPJ", "GSN", "KCNAB1", "LDB2", "THEMIS")
makeMarkerDotplot(harmonyRNA_and_ATAC, "RNA", rev(rnaGeneralMarkers), "General_RNA_Markers", 
						"highLevelCellType", outputPath=out_dir)


# PLCB1 empirically up in vascular. High in capillary ECs (https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.120.052318)
# NRG3 high in endocardial but not other ednothelial in hocker et al
# GJA4 high in pericytes in hocker et al
# RGS5 as measure of pericytes (https://pubmed.ncbi.nlm.nih.gov/18038251/)

rnaFigureMarkers = c("TTN", "MYH7", "RBPJ", "GSN", "DCN", "PLCB1", "CCL21", "VWF", "NRG3", "RGS5", "CD247", "JCHAIN", "KIT", "CPA3", "NRXN3", "PDE3B" )
# "KCNAB1", "LDB2",
makeMarkerDotplot(harmonyRNA_and_ATAC, "RNA", rev(rnaFigureMarkers), "Figure_Panel_RNA_Markers", 
						"highLevelCellType", outputPath=out_dir, typePlotOrder=cellTypePlotOrder)



makeMarkerDotplot(harmonyRNA_and_ATAC, "ATAC", rnaGeneralMarkers, "General_RNA_Markers", 
						"harmonyKNN_type", outputPath=out_dir)


#T markers
tMarkers = c("CD3D", "CD3G", "CD3E")
makeMarkerDotplot(harmonyRNA_and_ATAC, "ATAC", tMarkers	, "CD3_Markers", 
						"harmonyKNN_type", outputPath=out_dir)




atacDeRes1 = c("TTN", "MROH2A", "VWF", "CBS", "U2AF1", "CRYAA", "LRMBDA", "SULT1A4", "NRXN1", "LPP1", "FLT1")
makeMarkerDotplot(harmonyRNA_and_ATAC, "ATAC", atacDeRes1, "ATAC_DE_Marker_Tests_1", 
						"harmonyKNN_type", outputPath=out_dir)


hockerATACmarkers = c("DCN", "MYH7", "NPPA", "EGFL7", "GJA4", "MS4A6A", "IL7R", "ADIPOQ", "GPM6B")
makeMarkerDotplot(harmonyRNA_and_ATAC, "ATAC", hockerATACmarkers	, "Hocker_Markers", 
						"harmonyKNN_type", outputPath=out_dir)

empricalRNAmarkers = c("TTN", "MYH7", "GPAM" )
makeMarkerDotplot(harmonyRNA_and_ATAC, "ATAC", rnaGeneralMarkers, "General_RNA_Markers", 
						"harmonyKNN_type", outputPath=out_dir)




# # See if a marker file has already been found
# DE_ATAC_File = paste0(inputDir, outputNote, "_based_Celltype_DE_testing.csv")
# if (file.exists(DE_ATAC_File)){
# 	print("Reading calculated and stored DE results")
# 	deTestRes = read.csv(DE_ATAC_File)
# } else {
# 	print("Finding DE Results Now")
# 	deTestRes = runDEtestingToID_markers(harmonyRNA_and_ATAC[,colData(harmonyRNA_and_ATAC)$tech == "ATAC"], outputNote, "harmonyKNN_type",
#                   howManyGenesToTest = 50, outputPath=out_dir)
# 	deTestRes = deTestRes$marker_test_res[order(deTestRes$marker_test_res$cell_group),]
# 	write.csv(deTestRes, DE_ATAC_File)
# }

# DE_RNA_File = paste0(inputDir, outputNote, "_based_Celltype_RNA_DE_testing.csv")
# if (file.exists(DE_RNA_File)){
# 	print("Reading calculated and stored DE results")
# 	deTestResRNA = read.csv(DE_RNA_File)
# } else {
# 	print("Finding DE Results Now")
# 	deTestResRNA = runDEtestingToID_markers(harmonyRNA_and_ATAC[,colData(harmonyRNA_and_ATAC)$tech == "RNA"], outputNote, "highLevelCellType",
#                   howManyGenesToTest = 50, outputPath=out_dir)
# 	deTestResRNA = deTestResRNA$marker_test_res[order(deTestResRNA$marker_test_res$cell_group),]
# 	write.csv(deTestResRNA, DE_RNA_File)
# }







# markerTestDir = paste0(out_dir, "testingATAC_gene_markers/")
# dir.create(markerTestDir)

# makeDotplotsFromDE(deTestRes, markerTestDir, harmonyRNA_and_ATAC, opt)
# makeDotplotsFromDE(deTestRes, markerTestDir, harmonyRNA_and_ATAC, opt, sortPolicy="specificity")
# makeDotplotsFromDE(deTestRes, markerTestDir, harmonyRNA_and_ATAC, opt, sortPolicy="pseudo_R2")
# makeDotplotsFromDE(deTestRes, markerTestDir, harmonyRNA_and_ATAC, opt, sortPolicy="fraction_expressing")


# # Show markers from RNA DE testing

# markerTestDirRNA = paste0(out_dir, "testingRNA_gene_markers/")
# dir.create(markerTestDirRNA)

# makeDotplotsFromDE(deTestResRNA, markerTestDirRNA, harmonyRNA_and_ATAC, opt, assayToUse="RNA", cellTypeCol = "highLevelCellType")
# makeDotplotsFromDE(deTestResRNA, markerTestDirRNA, harmonyRNA_and_ATAC, opt, sortPolicy="specificity", assayToUse="RNA", cellTypeCol = "highLevelCellType")
# makeDotplotsFromDE(deTestResRNA, markerTestDirRNA, harmonyRNA_and_ATAC, opt, sortPolicy="pseudo_R2", assayToUse="RNA", cellTypeCol = "highLevelCellType")

# # Show these markers on the ATAC side
# makeDotplotsFromDE(deTestResRNA, markerTestDir, harmonyRNA_and_ATAC, opt, outputNote="RNA_DE_")

# makeDotplotsFromDE(deTestResRNA, markerTestDir, harmonyRNA_and_ATAC, opt, outputNote="RNA_DE_", sortPolicy="pseudo_R2", nToPlot=20)
# makeDotplotsFromDE(deTestResRNA, markerTestDir, harmonyRNA_and_ATAC, opt, outputNote="RNA_DE_", sortPolicy="specificity", nToPlot=20)

# makeDotplotsFromDE(deTestResRNA, markerTestDir, harmonyRNA_and_ATAC, opt, outputNote="RNA_DE_", nToPlot=25)






# setToCoPlot = c("TTN", "RBPJ", "BANK1", "PDGFRB", "VWF")
# makeMarkerDotplot(harmonyRNA_and_ATAC, "RNA", setToCoPlot, "Matching_Markers", 
# 						"highLevelCellType", outputPath=out_dir)
# makeMarkerDotplot(harmonyRNA_and_ATAC, "ATAC", setToCoPlot, "Matching_Markers", 
# 						"harmonyKNN_type", outputPath=out_dir)




# setToCoPlotMismatch = c("CD247", "FLT1", "PKDL1", "CCL21", "PROX1", "DCN")
# makeMarkerDotplot(harmonyRNA_and_ATAC, "RNA", setToCoPlotMismatch, "Mismatched_Markers", 
# 						"highLevelCellType", outputPath=out_dir)
# makeMarkerDotplot(harmonyRNA_and_ATAC, "ATAC", setToCoPlotMismatch, "Mismatched_Markers", 
# 						"harmonyKNN_type", outputPath=out_dir)



