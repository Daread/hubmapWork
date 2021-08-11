
# Get functions
library(monocle3)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
# modStatus <- loadMonoclePackages()
library(tidyr)
library(ggplot2)
# Note: Need bedtools loaded for this! bedtools/2.29.2
library(cicero)

library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-c", "--ciceroFile"), type="character", 
  			default="Cicero_Defaults_cds_p_W144.heart.apex.s1_ciceroConnections.RDS", 
              help="Name of cds to read in and process", metavar="character"),

  make_option(c("-p", "--ciceroPath"), type="character", 
  			default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_24_ciceroWork_nobackup/rdsOutputs/", 
              help="Path to cds to process", metavar="character"),

  make_option(c("-f", "--ciceroCDS"), type="character", 
        default="Cicero_Defaults_cds_p_W144.heart.apex.s1_ciceroCDS.RDS", 
              help="Path to cds to process", metavar="character"),  

  make_option(c("-i", "--ciceroInputCDS"), type="character", 
        default="Cicero_Defaults_cds_p_W144.heart.apex.s1_ciceroInput_CDS.RDS", 
              help="CDS before binning", metavar="character"),  
 # "_ciceroInput_CDS.RDS" 

  make_option(c("-g", "--genomeFile"), type="character", 
  			# default="/net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished", 
        # The default is a soft link to "/net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished"
  			default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/data/genomeData/bbiHg38.fasta",
              help="Path to genome fasta to use", metavar="character"),

  make_option(c("-t", "--geneAnnotationsForTSS"), type="character", 
        default="/net/bbi/vol1/data/genomes_stage/human/human_atac/gene_bodies.bed.gz",
              help="Path to genome fasta to use", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Read the human gene_bodies.bed file into a dataframe
geneBodiesDF = as.data.frame(read.table(gzfile(opt$geneAnnotationsForTSS)))
# Process to get on TSS per gene. For now, use this as the first position in a gene_body entry (or last, if antisense)
geneBodiesDF$TSSpoint = ifelse(geneBodiesDF$V6 == "+", geneBodiesDF$V2, geneBodiesDF$V3)
geneBodiesDF$start = (geneBodiesDF$TSSpoint - 1)
geneBodiesDF$end =  (geneBodiesDF$TSSpoint + 1)

# Get a new DF
tssBed = data.frame("chr" = geneBodiesDF$V1, "start" = geneBodiesDF$start, "end" = geneBodiesDF$end, "geneName" = geneBodiesDF$V4,
                    "V5"=geneBodiesDF$V5, "orientation"=geneBodiesDF$V6)

# Write this output
outFileTSS = paste0("./fileOutputs/tssFromGeneBodyEnd.bed")
# saveRDS(tssBed, file=outFileTSS)
write.table(tssBed, file=outFileTSS, 
    quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

# Get the Cicero connection data
ciceroDF = readRDS(paste0(opt$ciceroPath, opt$ciceroFile))

ciceroBaseSites = data.frame("Peaks" = levels(as.factor(ciceroDF$Peak1)))

# Get peaks into bed format
ciceroPeaksAsBed = data.frame(separate(data=ciceroBaseSites, col = "Peaks", into = c("chr", "start", "end"), sep="_", remove=TRUE))
ciceroPeaksAsBed = separate(data=ciceroPeaksAsBed, col = "chr", into = c("Null", "chr"), sep="chr", remove=TRUE)
ciceroPeaksAsBed$name = "Peak"
ciceroPeaksAsBed = ciceroPeaksAsBed[c("chr", "start", "end", "name")]

# Output this
outFileATAC_sites = paste0("./fileOutputs/ATAC_Peaks_", opt$ciceroFile, ".bed")
write.table(ciceroPeaksAsBed, file=outFileATAC_sites, 
    quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

# Sort them both
sortedATAC = paste0("./fileOutputs/ATAC_Peaks_", opt$ciceroFile, "_sorted.bed")
system(paste0("bedtools sort -i ", outFileATAC_sites, " > ", sortedATAC))

sortedTSS = paste0("./fileOutputs/tssFromGeneBodyEndSorted.bed")
system(paste0("bedtools sort -i ", outFileTSS, " > ", sortedTSS))

# Intersect
intersectBed = paste0("./fileOutputs/TSS_Intersected_Peaks_", opt$ciceroFile, ".bed")
system(paste0("bedtools intersect -b ", sortedTSS, " -a ", sortedATAC, " -wa -wb > ", intersectBed))

# Get the intersected sites
intersectedDF = as.data.frame(read.table(intersectBed))
intersectedDF$formattedPeak = paste0("chr", intersectedDF$V5, "_", intersectedDF$V2, "_", intersectedDF$V3)

# Get the subset of cicero links going from a TSS
ciceroTSS_DF = ciceroDF[ciceroDF$Peak1 %in% intersectedDF$formattedPeak,]

# Take a look at the coaccessibility scores from TSS's.
png(paste0("./plots/CoaccessHist",  "_TSS_Only_", opt$ciceroFile, ".png"), 
        width=1000, height=1000, res = 200)
thisPlot = ggplot(ciceroTSS_DF, aes(coaccess)) + 
        geom_histogram()
print(thisPlot)
dev.off()

# Get an idea of distribution of positive, non-zero ones
png(paste0("./plots/PosNonzero_CoaccessHist",  "_TSS_Only_", opt$ciceroFile, ".png"), 
        width=1000, height=1000, res = 200)
thisPlot = ggplot(ciceroTSS_DF[ciceroTSS_DF$coaccess > .01,], aes(coaccess)) + 
        geom_histogram()
print(thisPlot)
dev.off()






# Use this intersection into within the cicero cds
ciceroCDS = readRDS(paste0(opt$ciceroPath, opt$ciceroInputCDS))
rowData(ciceroCDS)$Peak = rownames(rowData(ciceroCDS))

tssIntersectedPeakVec = intersectedDF$V8
names(tssIntersectedPeakVec) = intersectedDF$formattedPeak

# Update coldata to flag if a site intersects a TSS
rowData(ciceroCDS)$gene = ifelse(rowData(ciceroCDS)$Peak %in% names(tssIntersectedPeakVec),
                                     tssIntersectedPeakVec[rowData(ciceroCDS)$Peak],
                                     NA)

# geneVec = rowData(ciceroCDS)$gene
# geneVec[1:100]

# Try getting activity scores
unnorm_ga = build_gene_activity_matrix(ciceroCDS, ciceroDF)

unnorm_ga <- unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]

# make a list of num_genes_expressed
num_genes <- colData(ciceroCDS)$num_genes_expressed
names(num_genes) <- row.names(colData(ciceroCDS))

# normalize
cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)

# Make a cds out of this
activityDF = data.frame("GeneName" = rownames(cicero_gene_activities))
rownames(activityDF) = rownames(cicero_gene_activities)
activityCDS = new_cell_data_set(cicero_gene_activities,
                            cell_metadata=colData(ciceroCDS),
                            gene_metadata=activityDF)

set.seed(7)
activityCDS = estimate_size_factors(activityCDS)
activityCDS = preprocess_cds(activityCDS)
activityCDS = reduce_dimension(activityCDS)

# Plot
plotUMAP_Monocle(activityCDS, "Cicero_Activity_Based", "umi", show_labels=FALSE)

statsToPlot = c("RIP", "FRIP", "doublet_score", "num_genes_expressed")
for (eachStat in statsToPlot){
  plotUMAP_Monocle(activityCDS, "Cicero_Activity_Based", eachStat, show_labels=FALSE)
}




































