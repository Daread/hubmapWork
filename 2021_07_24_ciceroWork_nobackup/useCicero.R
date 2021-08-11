
# Get functions
library(monocle3)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
# modStatus <- loadMonoclePackages()
# library(tidyr)
library(cicero)


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





# Read in the CDS
pathToCDS = opt$cdsPath
loadedName = load(paste0(pathToCDS, opt$cdsName))

myCDS = cds_p

processingNote = "Cicero_Defaults"

# Cicero needs UMAP coordinates initially. Get these
set.seed(7)
myCDS = detect_genes(myCDS)
myCDS = estimate_size_factors(myCDS)
myCDS = preprocess_cds(myCDS, method='LSI')
myCDS = reduce_dimension(myCDS, reduction_method = "UMAP",
                          preprocess_method = "LSI")

# Set up a cicero cds object
ciceroCDS = make_cicero_cds(myCDS, reduced_coordinates = reducedDims(myCDS)$UMAP)

# Save this CDS. Will be used later
cdsOutput = paste0("./rdsOutputs/", processingNote, "_", opt$cdsName, "_ciceroCDS.RDS" )
saveRDS(ciceroCDS, file=cdsOutput)

inputCDSoutput = paste0("./rdsOutputs/", processingNote, "_", opt$cdsName, "_ciceroInput_CDS.RDS" )
saveRDS(myCDS, file=inputCDSoutput)

# Get genome coordinates for Hg38 (the reference used in the BBI pipeline)
# data("human.hg19.genome")

# Read in the chromosome lengths
genomeLengthDF = as.data.frame(read.table(file = "./hg38chromSizesFormattedForCicero.txt", sep='\t',
                                          header=FALSE))
colnames(genomeLengthDF) = c("V1", "V2")

# Run cicero
connections = run_cicero(ciceroCDS, genomeLengthDF)
# Save the output, return later

outFile = paste0("./rdsOutputs/", processingNote, "_", opt$cdsName, "_ciceroConnections.RDS")
saveRDS(connections, file=outFile)




