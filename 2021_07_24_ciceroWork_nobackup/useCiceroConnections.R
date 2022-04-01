
# Get functions
library(monocle3)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
# modStatus <- loadMonoclePackages()
# library(tidyr)
library(cicero)
library(ggplot2)


library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-c", "--cdsName"), type="character", 
  			# default="cds_p_allHeartATAC", 
        # default="cds_p_W144.heart.apex.s1", 
        default="cds_p_allHeartATACFRIP=0.2_FRIT=0.05UMI=1000DL=0.5",
              help="Name of cds to read in and process", metavar="character"),

  make_option(c("-p", "--cdsPath"), type="character", 
  			default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/cds_objects/", 
              help="Path to cds to process", metavar="character"),

  make_option(c("-g", "--genomeFile"), type="character", 
        default="/net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished", 
        # default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/data/genomeData/bbiHg38.fasta",
              help="Path to genome fasta to use", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



processingNote = "Cicero_Defaults"

rdsFile = paste0("./rdsOutputs/", processingNote, "_", opt$cdsName, "_ciceroConnections.RDS")
connectionData =  readRDS(rdsFile)

# Get an idea of the distribution of co-accessibility scores
png(paste0("./plots/CoaccessHist", processingNote, "_", opt$cdsName, ".png"), 
        width=1000, height=1000, res = 200)
thisPlot = ggplot(connectionData, aes(coaccess)) + 
        geom_histogram()
print(thisPlot)
dev.off()

# Get an idea of distribution of positive, non-zero ones
png(paste0("./plots/PosNonzero_CoaccessHist", processingNote, "_", opt$cdsName, ".png"), 
        width=1000, height=1000, res = 200)
thisPlot = ggplot(connectionData[connectionData$coaccess > .01,], aes(coaccess)) + 
        geom_histogram()
print(thisPlot)
dev.off()







































