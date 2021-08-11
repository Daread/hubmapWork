
# Get functions
library(monocle3)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
# modStatus <- loadMonoclePackages()
library(tidyr)
# Note: Need to load bedtools prior to running this script!
# library("bedr")
# module load bedtools/2.29.2


library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-c", "--cdsName"), type="character", 
  			# default="cds_p_allHeartATAC", 
        default="cds_p_allHeartATACFRIP=0.2_FRIT=0.05UMI=1000DL=0.5",
              help="Name of cds to read in and process", metavar="character"),

  make_option(c("-p", "--cdsPath"), type="character", 
  			default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/cds_objects/", 
              help="Path to cds to process", metavar="character"),

  make_option(c("-g", "--genomeFile"), type="character", 
  			# default="/net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished", 
        # The default is a soft link to "/net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished"
  			default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/data/genomeData/bbiHg38.fasta",
              help="Path to genome fasta to use", metavar="character"),

  make_option(c("-s", "--seqLength"), type="numeric", 
        default=600, 
              help="Length of sequence to extract for each peak, centered at peak midpoint", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)





# Read in the CDS
pathToCDS = opt$cdsPath
loadedName = load(paste0(pathToCDS, opt$cdsName))

# From this, get the sequences. the "load" function will store the data as a variable named "cds_p"
peakLocDF = data.frame("peakPos" = rownames(as.data.frame(rowData(cds_p))) )

peakLocDF = separate(data=peakLocDF, col = "peakPos", into = c("chr", "start", "end"), sep="_", remove=FALSE)
peakLocDF = separate(data=peakLocDF, col = "chr", into = c("Null", "chr"), sep="chr", remove=TRUE)


# Re-order like a bed file
peakLocDF = peakLocDF[,c(3,4,5,1)]

# 8-10-21 add: get sequence that's centered at peak midpoint but with a defined size
peakLocDF$midpoint = as.integer((as.numeric(peakLocDF$start) + as.numeric(peakLocDF$end))/2.0)

# Replace start and end values
peakLocDF$start = peakLocDF$midpoint - (opt$seqLength / 2.0)
peakLocDF$end = peakLocDF$midpoint + (opt$seqLength / 2.0)

# Update the peakPos to use these new coordinates or merging will fail later
# peakLocDF$peakPos = paste0("chr", peakLocDF$chr, "_", as.character(peakLocDF$start), "_", as.character(peakLocDF$end))

peakLocDF = peakLocDF[,c("chr", "start", "end", "peakPos")]
peakLocDF$origPeak = peakLocDF$peakPos
peakLocDF$peakPos = paste0("chr", peakLocDF$chr, "_", as.character(peakLocDF$start), "_", as.character(peakLocDF$end))

# Write as a tsv (makes this a bed file)
bedFileName = paste0("./rdsOutput/", opt$cdsName, "_peaks.bed")
write.table(peakLocDF, file=bedFileName, 
		quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

seqFile = paste0("./rdsOutput/", opt$cdsName, "_sequences.tsv" )
# Use this newly generated bed file to get sequence
system(paste0("bedtools getfasta -fi ", opt$genomeFile, " -bed ", bedFileName, 
				" -tab > ", seqFile ))

# Now read in the table
peaksWithSeq = read.table(seqFile, sep='\t', header=FALSE,  col.names=c("peakPos", "seq"))
# Fix formatting to match the ATAC CDS
peaksWithSeq$peakPos = paste0("chr", peaksWithSeq$peakPos)
peaksWithSeq$peakPos = gsub("-", "_", peaksWithSeq$peakPos)
peaksWithSeq$peakPos = gsub(":", "_", peaksWithSeq$peakPos)




# Merge the sequence back in
combinedDF = merge(peakLocDF, peaksWithSeq, by="peakPos")



# Save the sequence-loaded df
write.csv(combinedDF, paste0("./rdsOutput/", opt$cdsName, "_peaksWithSeqDF", as.character(opt$seqLength), ".csv"))

# Also we need this info as a fasta format
fastaSeqFile = paste0("./rdsOutput/", opt$cdsName, "_peakSeqsAsFasta", as.character(opt$seqLength), ".fa")
system(paste0("bedtools getfasta -fi ", opt$genomeFile, " -bed ", bedFileName, 
				" > ", fastaSeqFile ))























