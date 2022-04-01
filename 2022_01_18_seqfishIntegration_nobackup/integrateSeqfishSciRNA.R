

# load requirements
suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
})
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

library("optparse")
library("harmony")
library('class')
library(ggplot2)
# Get the passed parameters

option_list = list(
  # make_option(c("-p", "--processingNote"), type="character", 
  #       default="HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes", 
  #             help="Processing note from upstream CDS", metavar="character"),
    make_option(c("-s", "--sampleRNAname"), type="character",  default= "All_Cells", # "W144.Apex",  # default="W142.Left.Vent",  # "All_Cells"
              help="Sample to integrate", metavar="character"),

    make_option(c("-r", "--saveHarmonyCDS"), type="logical", default=TRUE,  # "ArchR", "Cicero"
              help="Save a cds of RNA+ATAC data after running harmony alignment", metavar="character"),

   make_option(c("-p", "--pcToUse"), type="integer", default=20,  # "ArchR", "Cicero"
              help="PCs to use during alignment/co-embedding", metavar="integer"),

    make_option(c("-q", "--seqfishPath"), type="character", default="./fileOutputs/", 
              help="Path to seqfish files", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



getSeqfishCDS <- function(opt, samplesToRead){
  # Read in all samples
  sampleList = vector(mode="list", length=length(samplesToRead))
  for (eachInd in 1:length(samplesToRead)){
    sampleList[[eachInd]] = read.csv(paste0(opt$seqfishPath, samplesToRead[eachInd]))
    rownames(sampleList[[eachInd]]) = sampleList[[eachInd]]$X 
    sampleList[[eachInd]] = sampleList[[eachInd]][,!(names(sampleList[[eachInd]]) %in% c("X"))]

  }
  # Combine into a single df
  fullExprs = bind_rows(sampleList)

  # This is cell x gene. Flip to gene x cell, and set up cell and gene metadata
  fullExprs = t(as.matrix(fullExprs))
  # browser()

  colData = data.frame(n.umi = colSums(fullExprs))
  rowData = data.frame(gene_short_name = rownames(fullExprs))

  rownames(colData) = colnames(fullExprs)
  rownames(rowData) = rownames(fullExprs)

  # Make a CDS
  seqfishCDS = new_cell_data_set(fullExprs,
                                cell_metadata = colData,
                                gene_metadata = rowData)
  return(seqfishCDS)
}

getSeqfishParams <- function(){
  seqfishParmas = list()
  seqfishParmas["umiMin"] = 100
  return(seqfishParmas)
}

makeQCplots <- function(inputCDS, procNote, outdir = "./plots/"){
  # Histogram of umis
  png(paste0(outdir, procNote, "_UMI_Hist.png"), res=200, height=800, width=1000)
  myPlot = ggplot(as.data.frame(colData(inputCDS)), aes_string("n.umi")) + 
          geom_histogram()
  print(myPlot)
  dev.off()
}

filterSeqfishCDS <- function(inputCDS, seqfishProcParams){
  # By umi
  inputCDS = inputCDS[, colData(inputCDS)$n.umi >= seqfishProcParams$umiMin]

  return(inputCDS)
}



seqfishSamples = c( #"W159_Heart_RV_table1.csv",  # Seeing repeats between this and table2, which appears to be a superset
              "W159_Heart_RV_table2.csv")

seqfishCDS = getSeqfishCDS(opt, seqfishSamples)
seqfishProcParams = getSeqfishParams()

dir.create('./plots/')
procNote = "preFilter"
makeQCplots(seqfishCDS, procNote)

# Filter
seqfishCDS = filterSeqfishCDS(seqfishCDS, seqfishProcParams)
procNote = "postFilter"
makeQCplots(seqfishCDS, procNote)

######################################################

# Read in the RNAseq cds
rnaPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/"
inputRNA_CDSname = "allCells_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes.rds"
rnaCDS = readRDS(paste0(rnaPath, inputRNA_CDSname))


###############################################
#Integrating seqfish with sci RNA

# First, get the names that match
commonGeneNames = rowData(seqfishCDS)$gene_short_name %in% rowData(rnaCDS)$gene_short_name
rna_genesInSeqfish = rnaCDS[rowData(rnaCDS)$gene_short_name %in% rowData(seqfishCDS)$gene_short_name,]
seqfish_genesInRNA = seqfishCDS[rowData(seqfishCDS)$gene_short_name %in% rowData(rnaCDS)$gene_short_name,]

geneShortnameList = rownames(rowData(rnaCDS))
names(geneShortnameList) = rowData(rnaCDS)$gene_short_name

rownames(seqfish_genesInRNA) = unname(geneShortnameList[rowData(seqfish_genesInRNA)$gene_short_name])

# See the UMI distribution of what remains in the RNA data
rna_genesInSeqfish$n.umi = colSums(exprs(rna_genesInSeqfish))
makeQCplots(rna_genesInSeqfish, "RNA_genes_in_seqfish")
# Keep above some nonzero cutoff
rnaSeqfishMatch_umiMin = 15
rna_genesInSeqfish = rna_genesInSeqfish[,colData(rna_genesInSeqfish)$n.umi >= rnaSeqfishMatch_umiMin]

# Note the tech
colData(rna_genesInSeqfish)$tech = "sciRNA"
colData(seqfish_genesInRNA)$tech = "seqfish"




# Try to simply merge the two cds objects, run through PCA, and align
combinedCDS = combine_cds(list(rna_genesInSeqfish, seqfish_genesInRNA))

set.seed(7)
combinedCDS = estimate_size_factors(combinedCDS)
combinedCDS = preprocess_cds(combinedCDS, num_dim = opt$pcToUse)

procNote = "coembedSharedGenes"
# Check: Does simply running UMAP look credible?
combinedCDS = reduce_dimension(combinedCDS)


plotUMAP_Monocle(combinedCDS, procNote, "tech",
                    show_labels=FALSE, textSize=10, outputPath = "./plots/",  returnPlotObj=FALSE)

plotUMAP_Monocle(combinedCDS, procNote, "highLevelCellType",
                    show_labels=FALSE, textSize=10, outputPath = "./plots/",  returnPlotObj=FALSE)


###################################################################
# Try running harmony and see the resulting UMAP from aligned PCA embedding
harmonyPCs = harmony::HarmonyMatrix(reducedDims(combinedCDS)$PCA, 
                    as.data.frame(colData(combinedCDS)), c("tech"), do_pca = FALSE, verbose=TRUE)

harmonyCDS = combinedCDS
reducedDims(harmonyCDS)$PCA = harmonyPCs

harmonyCDS = reduce_dimension(harmonyCDS, preprocess_method="PCA", reduction_method="UMAP")

plotUMAP_Monocle(harmonyCDS, paste0(procNote,  "harmonyAlign"),
                             "highLevelCellType", show_labels=FALSE)

plotUMAP_Monocle(harmonyCDS, paste0(procNote,  "harmonyAlign"),
                             "tech", show_labels=FALSE)










