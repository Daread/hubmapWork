
library(tidyr)
library(ggplot2)
library(glmnet)
library(ROCR)
library(data.table)
library(monocle3)
library(gridExtra)
library(RColorBrewer)

# Get shared functions between this and runRegressionModel.R
source("./utilityFunctions.R")

source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-y", "--predictionFraming"), type="character", 
        default="Regression",# "Classification",   # "Classification" or "Regression"
              help="Option of which RNA data to predict", metavar="character"),
   make_option(c("-c", "--cdsRNA"), type="character", 
        default="allCells_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes.rds", 
              help="Name of RNA cds to use", metavar="character"),
  make_option(c("-x", "--rnaPath"), type="character", 
        default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/", 
              help="Path to RNA cds", metavar="character"),

  make_option(c("-c", "--cdsATAC"), type="character", 
        default="Harmony_Aligned_CoordsRegress_Protocadherin_PC_20_All_CellsFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_k_20_peak_cds_PeakMotifCDS.rds", 
              help="Name of RNA cds to use", metavar="character"),
  make_option(c("-x", "--atacPath"), type="character", 
        default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/rdsOutput/", 
              help="Path to RNA cds", metavar="character")

)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

source("./utilityFunctions.R")

# outDir = paste0("./fileOutputs/", opt$predictionFraming, "/")

# Swapped out on 1-10-22
# outputDir = paste0("./plots/", opt$predictionFraming, "/violinPlots/")
outputDir = paste0("./plots/",  "expressionPlots/")
dir.create(outputDir)


# Get the RNA CDS 
rnaData = readRDS(paste0(opt$rnaPath, opt$cdsRNA))

rnaData = hardAssignDonorAges(rnaData)
rnaData = hardAssignDonorSexes(rnaData)

# Get the ATAC cds
atacData = readRDS(paste0(opt$atacPath, opt$cdsATAC))
rowData(atacData)$gene_short_name = rowData(atacData)$Motif


plotPaneledPositive <- function(miniCDS, genesToPlot, cellTypeSubset, setName, 
									outputDir = paste0("./plots/",  "expressionPlots/"), cellTypeCol, colToGroupCells, colForFill){

	groupsToPlot = as.character(levels(as.factor(colData(miniCDS)[[colToGroupCells]])))
	if (colToGroupCells == cellTypeCol){
		groupsToPlot = groupsToPlot[groupsToPlot %in% cellTypeSubset]
	}

	# Loop. For each gene, get a plot for each subset of the cells
	panelList = vector(mode="list", length = length(genesToPlot) * length(groupsToPlot))
	# browser()
	thisPalette = "Set1"
	panelInd = 1
	for (eachGene in genesToPlot){
		for (eachGroup in groupsToPlot){
			cdsPlotSubset = miniCDS[rowData(miniCDS)$gene_short_name == eachGene, colData(miniCDS)[[colToGroupCells]] == eachGroup]
			# browser()
			thisPlot = plot_percent_cells_positive(cdsPlotSubset, group_cells_by=colToGroupCells,
          ) + scale_fill_brewer(palette=thisPalette) + scale_color_brewer(palette=thisPalette)
			# browser()

          #+ theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_bar(stat="identity", width=0.03) + 
				#ggtitle(paste0(eachGroup, " by ", colForFill))
			panelList[[panelInd]] = thisPlot
			panelInd = panelInd + 1
		}
	}

	# Now we have the full list of plots. Make this into a multi-panel plot and display
	# browser()
	resToUse = 200
	png(paste0(outputDir, "percent_Positive_for_", setName, ".png"), res=resToUse, 
		width = resToUse * 2.5 * length(groupsToPlot), height=  resToUse * 1.75 * length(genesToPlot))
	panelPlot = do.call("grid.arrange", c(panelList, ncol=length(groupsToPlot)))
	print(panelPlot)
	dev.off()
}



plotViolins <- function(genesToPlot, setName, rnaData, outputDir = paste0("./plots/",  "expressionPlots/"),  #outputDir = paste0("./plots/", opt$predictionFraming, "/violinPlots/"),
                  cellTypeSubset=FALSE, cellTypeCol = "highLevelCellType",
                     colToGroupCells="highLevelCellType", colForFill = NULL, plotPercentPos = TRUE) {

  miniCDS = rnaData[rowData(rnaData)$gene_short_name %in% genesToPlot,]
  # Subset by a celltype?
  if (!(cellTypeSubset == FALSE)){
    miniCDS = miniCDS[,colData(miniCDS)[[cellTypeCol]] %in% c(cellTypeSubset)]
  } else {
  	cellTypeSubset = as.character(levels(as.factor(colData(miniCDS)[[cellTypeCol]])))
  }
  # browser()

  # Now make a violin plot for this set
  png(paste0(outputDir, "violins_for_", setName, ".png"), res=200, height=2400, width=2400)
  myPlot = plot_genes_violin_DFR_customPseudo(miniCDS, group_cells_by=colToGroupCells,
         ncol=1, log_scale=TRUE, pseudocount=.01, leaveFillBlank=TRUE
          ) + theme(axis.text.x=element_text(angle=45, hjust=1))
  # Add a split by 
  if (!(is.null(colForFill))){
    myPlot = myPlot + aes_string(fill=colForFill)
  }
  print(myPlot)
  dev.off()

  # Now make a violin plot for this set using linear scales
  png(paste0(outputDir, "linear_violins_for_", setName, ".png"), res=200, height=2400, width=2400)
  myPlot = plot_genes_violin_DFR_customPseudo(miniCDS, group_cells_by=colToGroupCells,
         ncol=1, log_scale=FALSE, leaveFillBlank=TRUE
          ) + theme(axis.text.x=element_text(angle=45, hjust=1))
  # myPlot = plot_genes_violin(miniCDS, group_cells_by="highLevelCellType",
  #         ncol=1, log_scale=FALSE) + theme(axis.text.x=element_text(angle=45, hjust=1))
  if (!(is.null(colForFill))){
    myPlot = myPlot + aes_string(fill=colForFill)
  }
  print(myPlot)
  dev.off()

  print("Plotting Percent Positive")
  # 1-11-22: Percent positive plot as well. For these, re-run the sampling procedure splitting by 
  if (plotPercentPos & (length(cellTypeSubset) > 1)){
  	plotPaneledPositive(miniCDS, genesToPlot, cellTypeSubset, setName, outputDir=outputDir, cellTypeCol,
  						colToGroupCells, colForFill)
  } else {
  	thisPalette = "Dark2"
    png(paste0(outputDir, "percent_Positive_for_", setName, ".png"), res=200, height=2400, width=2400)
    myPlot = plot_percent_cells_positive(miniCDS, group_cells_by=colToGroupCells,
           ncol=1, #log_scale=TRUE, pseudocount=.01
            ) + theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_bar(stat="identity", width=0.03) +
            theme(text = element_text(size = 40))   + ylab("Percent Positive Cells") +
             scale_fill_brewer(palette=thisPalette) + scale_color_brewer(palette=thisPalette)
    # Add a split by 
    # if (!(is.null(colForFill))){
    #   myPlot = myPlot + aes_string(colForFill)
    # }
    print(myPlot)
    dev.off()
  }
}



plotViolins(c("CCN2", "TGFB3"), "TGFB_In_Fibroblasts", rnaData, cellTypeSubset=c("Fibroblast"), 
            colToGroupCells="Sex", outputDir = outputDir)


plotViolins(c("PROM1", "IL1R1", "ALPL", "ACVRL1"), "Top_Sex_VascEnd", rnaData, cellTypeSubset=c("Vascular_Endothelium"), 
            colToGroupCells="Sex", outputDir = outputDir)


# plotViolins(c("SPI1", "CEBPA"), "Marker_test_MultiTypes", atacData, cellTypeCol = "harmonyKNN_type",
#            #cellTypeSubset=c("Fibroblast", "Macrophage", "VSM_and_Pericyte", "Vascular_Endothelium"), 
#             colToGroupCells= "harmonyKNN_type", colForFill = "harmonyKNN_type" )

# plotViolins(c("Smad4"), "Smad4_In_MultiTypes", atacData, cellTypeCol = "harmonyKNN_type",
#            cellTypeSubset=c("Fibroblast", "Macrophage", "VSM_and_Pericyte", "Vascular_Endothelium"), 
#             colToGroupCells= "harmonyKNN_type", colForFill = "Sex" )

# plotViolins(c("FOS", "FOS::JUNB"), "FOS_and_JUN_In_MultiTypes", atacData, cellTypeCol = "harmonyKNN_type",
#            cellTypeSubset=c("Fibroblast", "Macrophage", "VSM_and_Pericyte", "Vascular_Endothelium", "Cardiomyocyte", "T_Cell"), 
#             colToGroupCells= "harmonyKNN_type", colForFill = "Sex" )

# atacTypeCol = "harmonyKNN_type"
# testCDS = atacData

# for (eachType in levels(as.factor(colData(testCDS)[[atacTypeCol]]))){

#   miniCDS = testCDS[,colData(testCDS)[[atacTypeCol]] == eachType]
#   print(eachType)

#   thisAv = sum(exprs(miniCDS[rowData(miniCDS)$Motif == "SPI1",])) / ncol(miniCDS)
#   print(paste0("SPI1 expression is ", thisAv))

#   totalReadsAv = sum(exprs(miniCDS)) / ncol(miniCDS)
#   print(paste0("Average motifs in peaks is ",  as.character(totalReadsAv)))
  
#   umiAv = sum(colData(miniCDS)$umi) / ncol(miniCDS)
#   print(paste0("UMI Average ", as.character(umiAv)))
#   print(" ")
# }




plotViolins(c("MYBL1", "MYBL2"), "MYBL_Across_Types", rnaData, #cellTypeSubset=c("Fibroblast"), 
             outputDir = outputDir)


plotViolins(c("MEF2B", "MEF2A", "MEF2C", "MEF2D"), "MEF_Across_Types", rnaData, #cellTypeSubset=c("Fibroblast"), 
             outputDir = outputDir)





plotViolins(c("RFX1", "RFX2", "RFX3", "RFX4", "RFX5", "RFX6", "RFX7", "RFX8"), "RFX_Across_Types", rnaData, #cellTypeSubset=c("Fibroblast"), 
             outputDir = outputDir)




plotCoefBars <- function(genesToPlot, setName, rnaData, outputDir = paste0("./plots/",  "expressionPlots/"),  #outputDir = paste0("./plots/", opt$predictionFraming, "/violinPlots/"),
                  cellTypeSubset=FALSE, cellTypeCol = "highLevelCellType",
                     colToGroupCells="highLevelCellType", colForFill = NULL) {

  miniCDS = rnaData[rowData(rnaData)$gene_short_name %in% genesToPlot,]
  # Subset by a celltype?
  if (!(cellTypeSubset == FALSE)){
    miniCDS = miniCDS[,colData(miniCDS)[[cellTypeCol]] == cellTypeSubset]
  }
  # browser()

  # Now make a violin plot for this set
  png(paste0(outputDir, "coefficients_for_", setName, ".png"), res=200, height=2400, width=2400)
  myPlot = plot_genes_violin_DFR_customPseudo(miniCDS, group_cells_by=colToGroupCells,
         ncol=1, log_scale=TRUE, pseudocount=.01
          ) + theme(axis.text.x=element_text(angle=45, hjust=1))
  # Add a split by 
  if (!(is.null(colForFill))){
    myPlot = myPlot + aes_string(colForFill)
  }

  # myPlot = plot_genes_violin(miniCDS, group_cells_by="highLevelCellType",
  #         ncol=1, log_scale=FALSE) + theme(axis.text.x=element_text(angle=45, hjust=1))
  print(myPlot)
  dev.off()

}



plotViolins(c("IL1R1"), "IL1R1_In_Vasculature", rnaData, cellTypeSubset="Vascular_Endothelium", 
            colToGroupCells="Sex", outputDir = outputDir)




plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("SMAD6", "SMAD7", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD1"), 
          "SMADS", outputPath = outputDir)

plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("ACVRL1", "TGFBR1", "SMURF1", "SMURF2"), 
          "TGF_Access", outputPath = outputDir)



plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("TGFBR3", "ENG", "EP300", "CBP", "CSKI", "SNON"), 
          "TGF_Access_2", outputPath = outputDir)

plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("MSTN", "NFATC2", "CALC"),
 "Myostatin_Sig", outputPath = outputDir)




plotViolins(c("THEMIS"), "Just_THEMIS", rnaData)



plotViolins(c("IL1R1", "ALPL", "PROM1", "LCN6"), "Vasculature_Cardiac_Hits", rnaData, cellTypeSubset="Vascular_Endothelium", 
            colToGroupCells="Sex")



plotViolins(c("BACH1"), "Just_BACH1", rnaData)


plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("THEMIS"), "Just_THEMIS", outputPath = outputDir)

plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("BACH1"), "Just_BACH1", outputPath = outputDir)





plotViolins(c("SOX15", "MEF2B", "CEBPA", "SPIB", "SPIC"), "Known_Plus_SPIC", rnaData)

plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("SOX15", "MEF2B", "CEBPA", "SPIB", "SPIC", "BACH1"),
                      "Known_Plus_SPIC", outputPath = outputDir)


plotViolins(c("TFAP2C", "TFAP2B"), "TFAP2_Set", rnaData)

plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("TFAP2C", "TFAP2B"),
                      "TFAP2_Set", outputPath = outputDir)



plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("MTAP", "CDKN2A", "CDKN2B", "DMRTA1"),
                      "Locus_9p21", outputPath = outputDir)
# "MTAP", "CDKN2A", "CDKN2B",


plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("DLG2", "CAMTA1", "CTNNA3", "NPIPA1", "CCSER1", "MAGI2"),
                      "Harmony_ATAC_T_cell_hits", outputPath = outputDir)


plotUMAP_Monocle_genes(rnaData, opt$predictionFraming, c("PCDHGA1", "PCDHGA2"),
                      "Protocadherins", outputPath = outputDir)

# # Find and save marker genes by cell type
# cdsAndDEtest = runDEtestingToID_markers(rnaData, "all_RNA_cell_markers", "highLevelCellType",
#                   howManyGenesToTest = 50, outputPath=outputDir)

# write.csv(cdsAndDEtest$marker_test_res[order(cdsAndDEtest$marker_test_res$cell_group),], paste0(outputDir, "all_RNA_cell_markers.csv"))





# Try some ATAC motif-in-peak count plotting

# plotViolins(c("EHF", "HMBOX1"), "VascMotifHits", atacData, cellTypeCol = "harmonyKNN_type",
#            cellTypeSubset="Vascular_Endothelium", 
#             colToGroupCells="Sex" )

# plotViolins(c("Rarg(var.2)", "RUNX1", "DBP"), "Macrophage_Hits_Motifs_Sex", atacData, cellTypeCol = "harmonyKNN_type",
#            cellTypeSubset="Macrophage", 
#             colToGroupCells="Sex" )

plotViolins(c("Smad4"), "Smad4_In_MultiTypes", atacData, cellTypeCol = "harmonyKNN_type",
           cellTypeSubset=c("Fibroblast", "Macrophage", "VSM_and_Pericyte", "Vascular_Endothelium"), 
            colToGroupCells="harmonyKNN_type", colForFill = "Sex" )



plotViolins(c("Rarg(var.2)", "RUNX1", "DBP"), "Macrophage_Hits_Motifs_Sex", atacData, cellTypeCol = "harmonyKNN_type",
           cellTypeSubset="Macrophage", 
            colToGroupCells="Sex" )




plotViolins(c("TGFB3", "CCN2"), "TGFB_In_Fibroblasts", rnaData, cellTypeSubset="Fibroblast", 
            colToGroupCells="Sex", outputDir = outputDir)




plotPercentPosScatter <-function(genesToPlot, setName, rnaData, outputDir = paste0("./plots/",  "expressionPlots/"),  #outputDir = paste0("./plots/", opt$predictionFraming, "/violinPlots/"),
                  cellTypeSubset=FALSE, cellTypeCol = "highLevelCellType",
                     colToGroupCells="highLevelCellType", colForFill = NULL) {

  miniCDS = rnaData[rowData(rnaData)$gene_short_name %in% genesToPlot,]
  # Subset by a celltype?
  if (!(cellTypeSubset == FALSE)){
    miniCDS = miniCDS[,colData(miniCDS)[[cellTypeCol]] == cellTypeSubset]
  }

  allDonors = levels(as.factor(colData(rnaData)$Donor))
  # Loop for each gene. For each gene, get a df with each cell & the expression of the gene of interest
  for (eachGene in genesToPlot){
    plotDF = as.data.frame(colData(miniCDS))
    plotDF$Expression = as.vector(exprs(miniCDS[rowData(miniCDS)$gene_short_name == eachGene,]))

    # Get the nonzero proportion per donor
    percentPosDF = data.frame()
    for (eachDonor in allDonors){
      subsetDF = plotDF[plotDF$Donor == eachDonor,]
      subsetDF$nonzero = ifelse(subsetDF$Expression > 0, 1, 0)
      thisProp = sum(subsetDF$nonzero) * 1.0 / nrow(subsetDF)
      # Can shortcut to get the age of the donor from the subset df
      thisAge = median(subsetDF$Age)
      donorDF = data.frame("Donor" = eachDonor, "Cells_Expressing" = thisProp, "Age"=thisAge)
      percentPosDF = rbind(percentPosDF, donorDF)
    }

    # Now plot 
    thisPlotTitle = paste0(outputDir, "Perc_Pos_Scatter_", eachGene, "_", setName, ".png")
    png(thisPlotTitle, res=200, width=1200, height=1200)
    myPlot = ggplot(percentPosDF, aes_string(x="Age", y="Cells_Expressing")) + 
            geom_point() + ylab(paste0(cellTypeSubset, " Expressing ", eachGene)) +
            theme(text = element_text(size = 24))+
  theme(plot.margin = margin(20,10,10,10, "pt"))
    print(myPlot)
    dev.off()
  }
}





plotPercentPosScatter(c("CFH", "CMKLR1"), "Inflam_In_Aged_Fibroblasts", rnaData, cellTypeSubset="Fibroblast",
              colToGroupCells="Age", outputDir=outputDir)


















