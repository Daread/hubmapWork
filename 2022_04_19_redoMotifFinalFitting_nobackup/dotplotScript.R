
library(tidyr)
library(ggplot2)
library(glmnet)
library(ROCR)
library(data.table)
# library(monocle3)

# Get shared functions between this and runRegressionModel.R
source("./utilityFunctions.R")

library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-p", "--pValFIMOcutoff"), type="numeric", 
  			default=1e-4, 
              help="Max FIMO p value to retain match", metavar="character"),

  make_option(c("-f", "--featureSelection"), type="character", 
        default="Binary_Combined_Motif_Counts", # "Binary_Combined_Motif_Counts", "Binary_PromOnly_Motif_Counts"
              help="Option of which features to use for input", metavar="character"),

  make_option(c("-t", "--predictionTask"), type="character", 
        default="_log2_CPM",   # "_log2_ratio_Vs_AllTypeMean", #  "_log2_CPM"
              help="Option of which RNA data to predict", metavar="character"),

  make_option(c("-y", "--predictionFraming"), type="character", 
        default="Regression",# "Classification",   # "Classification" or "Regression"
              help="Option of which RNA data to predict", metavar="character"),

  make_option(c("-z", "--highCutoff"), type="numeric", 
        default=0, 
              help="Minimum value, above which label = 1", metavar="numeric"),

  make_option(c("-l", "--lowCutoff"), type="numeric", 
        default=0, 
              help="Max value, below which label = 0", metavar="numeric"),

  make_option(c("-a", "--alphaToUse"), type="numeric", 
        default=0.5, 
              help="Alpha value to use in glment", metavar="numeric"),

  make_option(c("-u", "--promoterUpstream"), type="numeric", 
        default=2000,
              help="Bases upstream of TSS to use for input features", metavar="numeric"),
  make_option(c("-d", "--promoterDownstream"), type="numeric", 
        default=1000,
              help="Bases downstream of TSS to use for input features", metavar="numeric"),

  make_option(c("-k", "--coaccessCutoff"), type="numeric", 
        default=.015,
              help="Cutoff for keeping coaccessible peaks", metavar="numeric"),

  make_option(c("-n", "--maxNdistalSites"), type="numeric", 
        default=20,
              help="Max number of sites to link to a gene's promoter", metavar="numeric"),

  make_option(c("-q", "--peakSize"), type="numeric", 
        default=1000,
              help="-1 -> Keep as is, or enter a positive integer to make all peaks the same size", metavar="numeric"),

  make_option(c("-g", "--gtfFile"), type="character", 
        default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/fileOutputs/protCodingOnlyHumanGTF.gtf",
              help="Path and file for human gtf", metavar="character"),
   make_option(c("-c", "--cdsRNA"), type="character", 
        default="allCells_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes.rds", 
              help="Name of RNA cds to use", metavar="character"),
  make_option(c("-x", "--rnaPath"), type="character", 
        default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/", 
              help="Path to RNA cds", metavar="character")


)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
opt$variableParams = paste0(opt$promoterUpstream, "_", opt$promoterDownstream, "_",
                           opt$coaccessCutoff, "_", opt$maxNdistalSites, "_", opt$peakSize )

outDirName = paste0("./plots/", opt$predictionFraming, "/",
            "Prot_Only_Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize,"pVal", as.character(opt$pValFIMOcutoff), "/" )

# Read in from the csv where outputs were made at the end of final fitting
outDir = paste0("./fileOutputs/", opt$predictionFraming, "/")
inputFile = paste0(outDir,  "All_Celltype_Results_", opt$predictionFraming, "_Final_FullFit_Coefficients.csv")
allCoefDF = read.csv(inputFile)
# Drop row numbers
allCoefDF = allCoefDF[,-which(colnames(allCoefDF) == "X")]


cellTypes = c("Adipocytes", "B_Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic_Endothelium", "Macrophage", "Mast_Cell", "Neuronal", "T_Cell",
                "Vascular_Endothelium", "VSM_and_Pericyte")




formatCellType <- function(inputColumn){

  inputColumn = ifelse(inputColumn == "T_Cell", "T Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "VSM_and_Pericyte", "Perviascular Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "Vascular_Endothelium", "Vascular Endothelium", inputColumn)
  inputColumn = ifelse(inputColumn == "B_Cell", "B Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "Mast_Cell", "Mast Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "Lymphatic_Endothelium", "Lymphatic Endothelium", inputColumn)

  return(inputColumn)
}


# cellTypes = c( "Macrophage", 

# # Make a sub-directory for violin plots
# violinDir = paste0(outDir, "violinPlots/")
# dir.create(violinDir)


# plotViolins <- function(genesToPlot, setName, opt, outputDir = paste0("./fileOutputs/", opt$predictionFraming, "/violinPlots/")){
#   # Get the RNA data
  

#   # Now make a violin plot for this set
#   png(outputDir, res=200, height=1200, width=1200)
#   myPlot = plot_genes_violin(rnaData[rnaData$gene_short_name %in% genesToPlot,], group_cells_by="Cell_Type",
#           ncol=2) + theme(axis.text.x=element_text(angle=45, hjust=1))
#   print(myPlot)
#   dev.off()

# }




plotNonzeroCoefficients <- function(plotDF, cellTypes, opt){
  plotDF = plotDF[,!(names(plotDF) %in% c("Motif")) ]
  colnames(plotDF) = formatCellType(colnames(plotDF))

  # For each, get the proportion that is nonzero
  binaryMatrix = (as.matrix(plotDF) != 0) + 0
  nonzeroProps = colSums(binaryMatrix) / nrow(plotDF)
  cellTypeNames = colnames(plotDF)

  # Make into a df for ggplot
  propDF = data.frame("Cell_Type"=cellTypeNames, 
                      "Nonzero_Coefficients" = nonzeroProps)
  # browser()
  print(str(propDF))
  # Plot
  plotFile = paste0("./plots/", opt$predictionFraming, "/Nonzero_Coefficients.png")
  png(plotFile, res=200, height=1000, width=1200)
  myPlot = ggplot(propDF, aes_string(x="Cell_Type", y="Nonzero_Coefficients")) + 
           geom_col() + xlab("Cell Type") + ylab("Motifs Used") +
           theme(text=element_text(size=20)) + 
       theme(axis.text.x = element_text(angle = 45, hjust=1)) + ylim(0,1.0) +
       theme(plot.margin = unit(c(10,10,10,10), "pt"))
  print(myPlot)
  dev.off()

}

plotNonzeroCoefficients(allCoefDF, cellTypes, opt)




# Make a plot showing motifs x cell types for some high-abs motifs
makeMotifByTypePlot <- function(allCoefDF, cellTypes, definedMotifs=FALSE, motifsToUse=NULL, nPerCelltype=3, plotNote="",
                minZerosToPlot=0, sizeScaleMax=40, widthToUse=2800, heightToUse=2400, textSizeToUse=40){
  # For each celltype, get n top hits to show
  if (definedMotifs == FALSE){
    print("Finding top motifs by cell type")
    hitMotifs = c()
    # Filter down to rows with enough zeros
    allCoefDF = allCoefDF[rowSums(allCoefDF == 0) >= minZerosToPlot,]

    # Get the top nPerCellType TF motifs by coefficient magnitude
    for (eachCelltype in cellTypes){
      miniDF = data.frame(theseCoefs = allCoefDF[[eachCelltype]],
                motifs = allCoefDF$Motif)
      miniDF$absCoef = abs(miniDF$theseCoefs)
      # Reorder
      miniDF = miniDF[order(miniDF$absCoef, decreasing=TRUE),]
      # Get the first few coefs
      sortedMotifs = miniDF$motifs 
      hitMotifs = c(hitMotifs, sortedMotifs[1:nPerCelltype])
    }
    # Now get the unique ones
    uniqueHitMotifs = unique(hitMotifs)
    fileName = paste0("./plots/", opt$predictionFraming, "/", plotNote, "_", opt$predictionFraming, "_", minZerosToPlot, "zeros_",
             "_n", nPerCelltype, "_", opt$variableParams, "_", plotNote, "_MotifCoefDotplot.png")
    # If also provided, add defined motifs to use as well
    uniqueHitMotifs = uniqueHitMotifs[!(uniqueHitMotifs %in% motifsToUse)]
    uniqueHitMotifs = rev( c(uniqueHitMotifs, motifsToUse) )
  } else {
    fileName =paste0("./plots/", opt$predictionFraming, "/", plotNote, "_", opt$predictionFraming,  opt$variableParams, "_MotifCoefDotplot.png")
    uniqueHitMotifs = rev( motifsToUse)
  }
  
  # Get a formatted dataframe
  meltedDF = data.table::melt(allCoefDF, id.vars=c("Motif"), variable.name="Cell_Type", value.name="Coefficient")
  meltedDF$absCoefficient = abs(meltedDF$Coefficient)
  meltedDF$Coef_Sign = ifelse(meltedDF$Coefficient > 0, "Positive", "Negative")

  # Format the melted DF to not have the AllSeq part of names
  meltedDF = separate(meltedDF, col="Motif", into=c("Motif", "Junk"), sep="_", remove=FALSE)
  # browser()
  # uniqueHitMotifs = meltedDF$Motif
  uniqueHitMotifs = as.character(as.vector(as.data.frame(strsplit(uniqueHitMotifs, "_"))[1,]))

  # Now make a plot of coefficients
  png(fileName,
      res=200, width=widthToUse, height=heightToUse)
  myPlot = ggplot(meltedDF[meltedDF$Motif %in% uniqueHitMotifs,],
             aes(x=Cell_Type, y=factor(Motif, level=uniqueHitMotifs), 
              col=Coef_Sign, size=ifelse(absCoefficient==0, NA, absCoefficient))) +
            geom_point()+ 
       theme(axis.text.x = element_text(angle = 45, hjust=1)) +
      scale_x_discrete(breaks=cellTypes, labels=c("Adipocyte", "B Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic Endothelium", "Macrophage", "Mast Cell", "Neuronal", "T Cell",
                "Vascular Endothelium", "Perivascular Cell")) +
      theme(text=element_text(size=textSizeToUse)) + xlab("Cell Type") + 
      # theme(legend.title=element_text(color="Black",
      #                                    face="bold",size=0))+ 
            guides(size=guide_legend(title="Coefficient"))+
            scale_size(range = c(0, sizeScaleMax)) +
            ylab("Motif") + 
            guides(col=guide_legend(title="Coefficient Sign", override.aes = list(size=20))) 
  print(myPlot)
  dev.off()
}




# makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, 
#                   motifsToUse=paste0(c("SOX15", "MEF2B", "TEAD3", "SPIC"), "_AllSeq"),
#                      plotNote= "Expr_Model_Selective_Hits")





makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, 
                  motifsToUse=paste0(c("MEF2B", "MEF2A", "SOX15", "SPIC", "SPI1",
                                        "CEBPA", "YY1", "ZNF384", "OTX2",
                                        "Smad4"
                                         ), "_AllSeq"),
                     plotNote= "Paper_Select_Hits")











makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=FALSE, nPerCelltype = 6, minZerosToPlot=1,
                  motifsToUse=paste0(c("SOX15", "MEF2B", "SPIC"), "_AllSeq"), sizeScaleMax = 10,
                     plotNote= "Presentation_Examples", widthToUse=3000, heightToUse=3800, textSizeToUse=30)


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, 
                  motifsToUse=paste0(c("SOX15", "MEF2B", "SPIC"), "_AllSeq"),
                     plotNote= "Presentation_Examples", textSizeToUse=40)


makeMotifByTypePlot(allCoefDF, cellTypes, nPerCelltype=10, minZerosToPlot=2)
makeMotifByTypePlot(allCoefDF, cellTypes, nPerCelltype=10, minZerosToPlot=5)

makeMotifByTypePlot(allCoefDF, cellTypes, nPerCelltype=3, minZerosToPlot=5)


makeMotifByTypePlot(allCoefDF, cellTypes, nPerCelltype=3, minZerosToPlot=8)

makeMotifByTypePlot(allCoefDF, cellTypes)


makeMotifByTypePlot(allCoefDF, cellTypes, nPerCelltype=5)


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("Smad2..Smad3_AllSeq", "SMAD2..SMAD3..SMAD4_AllSeq",
                                        "SMAD3_AllSeq", "Smad4_AllSeq", "FOXH1_AllSeq",
                                        "FOS_AllSeq", "FOS..JUN_AllSeq", "FOS..JUN.var.2._AllSeq",
                                        "JUN_AllSeq", "SNAI1_AllSeq"), plotNote="TGFBeta_Set")

makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("GATA4_AllSeq", "GATA5_AllSeq", "GATA6_AllSeq", "MEF2_AllSeq",
                                "FOXO_AllSeq", "NKX2.5_AllSeq", "YY1_AllSeq", "HEY2_AllSeq",
                                "MITF_AllSeq"), plotNote="Hypertrophy_Set")


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("CEBPA_AllSeq", "MEF2A_AllSeq",
                                        "MEIS1_AllSeq", "ERG_AllSeq", "EBF2_AllSeq", "SPIB_AllSeq", 
                                          "RUNX1_AllSeq", "CEBPA_AllSeq", "SOX10_AllSeq"), plotNote="HockerEtAl_MatchMotif_Set")


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("SPI1_AllSeq", "SPIB_AllSeq", "GATA6_AllSeq",
                                                    "GATA4_AllSeq", "GATA2_AllSeq", "GATA1_AllSeq", "ESRRG_AllSeq",
                                                      "ESRRB_AllSeq", "ESRRA_AllSeq"), plotNote="HockerEtAl_Integrative_Set")


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("TEAD3_AllSeq", "MYF6_AllSeq", "NKX2.3_AllSeq", 
                                                                          "AP4_AllSeq", "FOS..JUN_AllSeq", "SMAD3_AllSeq",
                                                                            "EWS..ERG_AllSeq", "EBF_AllSeq", "SPIC_AllSeq",
                                                                            "MEF2A_AllSeq", "GRE_AllSeq", "NF1_AllSeq", "CEBPA_AllSeq",
                                                                            "AP1_AllSeq", "PGR_AllSeq", "ATF3_AllSeq", "ELF3_AllSeq"), 
                                                                  plotNote="HockerEtAl_HeartFail_Set")


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("HEY2_AllSeq", # Should be encodardium and cardiomyocytes https://www.nature.com/articles/s41598-018-20917-w
                                                                            "RORA_AllSeq", # A repressor https://www.pnas.org/content/116/42/21140
                                                                            "RREB1_AllSeq", # Repressor or activator, y context https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7085234/
                                                                            "BACH2_AllSeq", # Known repressor, matches hits in fetal atlas
                                                                            "ERG_AllSeq", # Vascular endothelial hit in fetal atlas
                                                                            "SRF_AllSeq", # hit in endocardium in fetal atlas
                                                                            "ZBTB18_AllSeq" # hit in lymphatic endothelium in fetal atlas
                                                                            ), plotNote="Fetal_and_Characterized_1_Set")


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("JUN_AllSeq", "SOX15_AllSeq", "BLHLHE23_AllSeq", "OLIG1_AllSeq", "OLIG2_AllSeq"  
                                                                            ), plotNote="Fetal_Heart_wide_Set")


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("NFATC3_AllSeq", # Putative repressor in fetal atlas, in T cells
                                                                            "FOXO3_AllSeq", # Ambiguous if activator or repressor
                                                                            "REST_AllSeq", # Repress through non-chromatin mech https://www-science-org.offcampus.lib.washington.edu/doi/10.1126/science.aba7612
                                                                            "MEF2B_AllSeq" # Cardiomyocyte hit
                                                                            ), plotNote="Fetal_and_Characterized_2_Set")


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("MEF2B_AllSeq",# Supp file 3
                                                                             "MEF2A_AllSeq",
                                                                             "MEF2B_AllSeq",
                                                                             "ESRRB_AllSeq",
                                                                             "NFIC..TLX1_AllSeq",
                                                                             "NR3C1_AllSeq", 
                                                                             "ESRRG_AllSeq", 
                                                                             "DUX4_AllSeq",
                                                                             "NKX2.3_AllSeq", 
                                                                             "HAND1..TCF3_AllSeq",
                                                                             "PKNOX1_AllSeq", "TGIF2_AllSeq", "NR3C2_AllSeq", "MECOM_AllSeq", "RARB_AllSeq", "NFIC_AllSeq", "RORA_AllSeq",
                                                                              "MAFB_AllSeq", "ESRRA_AllSeq", "GATA5_AllSeq", "GATA3_AllSeq", "NEUROG2_AllSeq", "NKX2.8_AllSeq", "GATA4_AllSeq"
                                                                            ), plotNote="Fetal_atlas_cardiomyocytes_1_Set")


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, motifsToUse=c("FOXO1_AllSeq", "IRF9_AllSeq", "SOX15_AllSeq", "TWIST1_AllSeq",
                                                                            "RARA_AllSeq", # Var2?
                                                                            "ZBTB18_AllSeq", "RELA_AllSeq", "SOX9_AllSeq", "STAT1_AllSeq", "IRF2_AllSeq", "SRF_AllSeq", 
                                                                            "IRF3_AllSeq", "MEF2B_AllSeq", "CREB1_AllSeq", "SOX6_AllSeq", "SOX13_AllSeq", "IRF7_AllSeq",
                                                                              "ETS1_AllSeq", "FLI1_AllSeq", "ERG_AllSeq", "GLI2_AllSeq", "MEF2A_AllSeq", "MEF2D_AllSeq", 
                                                                              "JUND_AllSeq", "STAT4_AllSeq", "ETV6_AllSeq"
                                                                            ), plotNote="Fetal_atlas_vascEndoth_1_Set")



makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, 
                  motifsToUse=paste0(c("SOX15", "MEF2B", "Ddit3..Cebpa", "SPIB", "SPIC", "BACH1"), "_AllSeq"),
                     plotNote= "Known_Plus_SPIC")



makeMotifByTypePlot(allCoefDF, cellTypes, nPerCelltype=30)


makeMotifByTypePlot(allCoefDF, cellTypes, definedMotifs=TRUE, 
                  motifsToUse=paste0(c("SOX15", "MEF2B", "TEAD3", "SPIC"), "_AllSeq"),
                     plotNote= "Knwon_TEAD3_SPIC")

