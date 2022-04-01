
# Get functions
library(monocle3)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
# modStatus <- loadMonoclePackages()
library(tidyr)
library(ggplot2)
# Note: Need bedtools loaded for this! bedtools/2.29.2
library(cicero)

# Add this to stop R from formatting into scientific notation, which causes errors in bedtools (which expects explciit genomic coordinates)
options(scipen=999)


library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-c", "--ciceroFile"), type="character", 
  			# default="Cicero_Defaults_cds_p_W144.heart.apex.s1_ciceroConnections.RDS", 
        default="Cicero_Defaults_cds_p_allHeartATACFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_ciceroConnections.RDS", 
              help="Name of cds to read in and process", metavar="character"),

  make_option(c("-p", "--ciceroPath"), type="character", 
  			default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_24_ciceroWork_nobackup/rdsOutputs/", 
              help="Path to cds to process", metavar="character"),

  make_option(c("-f", "--ciceroCDS"), type="character", 
        default="Cicero_Defaults_cds_p_allHeartATACFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_ciceroCDS.RDS", 
              help="Path to cds to process", metavar="character"),  

  make_option(c("-i", "--ciceroInputCDS"), type="character", 
        # default="Cicero_Defaults_cds_p_W144.heart.apex.s1_ciceroInput_CDS.RDS", 
        default="Cicero_Defaults_cds_p_allHeartATACFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_ciceroInput_CDS.RDS", 
              help="CDS before binning", metavar="character"),

  make_option(c("-k", "--coaccessCutoff"), type="numeric", 
        default=.05,
              help="Cutoff for keeping coaccessible peaks", metavar="numeric")
#   opt$cpmMin
# log2RatioVsMeanCutoff = opt$ratMin
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)





# Get the cicero connection file to use here
ciceroDF = readRDS(paste0(opt$ciceroPath, opt$ciceroFile))

factorCiceroDF = ciceroDF
factorCiceroDF$Peak1 = as.factor(factorCiceroDF$Peak1)
factorCiceroDF$Peak2 = as.factor(factorCiceroDF$Peak2)


# download and unzip
temp <- tempfile()
# download.file("ftp://ftp.ensembl.org/pub/release-65/gtf/mus_musculus/Mus_musculus.NCBIM37.65.gtf.gz", temp)
# ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz
download.file("ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz", temp)
gene_anno <- rtracklayer::readGFF(temp)
unlink(temp)

# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name


minSNPblock = 22072041
maxSNPblock = 22130390


getLinksInRange = function(inputDF, minPos, maxPos){
	inputDF$Peak1 = as.character(inputDF$Peak1)
	inputDF = separate(inputDF, col="Peak1", into=c("Peak_Chr", "MinPosPeak1", "MaxPosPeak1"), sep="_", remove=FALSE)

	# Keep only those in the proper range
	inputDF$MinPosPeak1 = as.numeric(inputDF$MinPosPeak1)
	inputDF$MaxPosPeak1 = as.numeric(inputDF$MaxPosPeak1)
	inputDF = inputDF[((inputDF$MinPosPeak1 > minPos) & (inputDF$MinPosPeak1 < maxPos)) | ((inputDF$MaxPosPeak1 > minPos) & (inputDF$MaxPosPeak1 < maxPos)),]

	# Only keep chr 9
	inputDF = inputDF[inputDF$Peak_Chr == "chr9",]

	return(inputDF[c("Peak1", "Peak2", "coaccess")])

}

linksFrom9p21 = getLinksInRange(ciceroDF, minSNPblock, maxSNPblock)

# set output
outDir = paste0("./plots/", "9p21_locus/")
dir.create(outDir)

png(paste0(outDir, "locus_nearby_Links_coaccess", as.character(opt$coaccessCutoff), ".png"),
			res=200, width=2000, height=1400)
# myPlot = plot_connections(factorCiceroDF, "chr9", 22072041, 22130390,
myPlot = plot_connections(factorCiceroDF, "chr9", 21750000, 22000000,
                 gene_model = gene_anno, 
                 coaccess_cutoff = opt$coaccessCutoff, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )

print(myPlot)
dev.off()



png(paste0(outDir, "locus_larger_Links_coaccess", as.character(opt$coaccessCutoff), ".png"),
			res=200, width=2000, height=1400)
# myPlot = plot_connections(factorCiceroDF, "chr9", 10965244, 12555547,
# myPlot = plot_connections(factorCiceroDF, "chr9", 22072041, 22130390,
myPlot = plot_connections(factorCiceroDF, "chr9", 21300000, 22500000,
                 gene_model = gene_anno, 
                 coaccess_cutoff = opt$coaccessCutoff, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )

print(myPlot)
dev.off()



find_overlapping_coordinates(factorCiceroDF$Peak1, "chr:22,072,041-22,130,390")



png(paste0(outDir, "outboundFrom_9p21_Snp_Block", as.character(opt$coaccessCutoff), ".png"),
			res=200, width=2000, height=1400)
# myPlot = plot_connections(factorCiceroDF, "chr9", 10965244, 12555547,
# myPlot = plot_connections(factorCiceroDF, "chr9", 22072041, 22130390,
myPlot = plot_connections(linksFrom9p21, "chr9",21300000, 22500000,
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.015, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )

print(myPlot)
dev.off()

png(paste0(outDir, "large_outboundFrom_9p21_Snp_Block", as.character(opt$coaccessCutoff), ".png"),
			res=200, width=2000, height=1400)
# myPlot = plot_connections(factorCiceroDF, "chr9", 10965244, 12555547,
# myPlot = plot_connections(factorCiceroDF, "chr9", 22072041, 22130390,
myPlot = plot_connections(linksFrom9p21, "chr9",21072041, 23130390,
                 gene_model = gene_anno, 
                 coaccess_cutoff = 0.015, 
                 connection_width = .5, 
                 collapseTranscripts = "longest" )

print(myPlot)
dev.off()





str(gene_anno[gene_anno$symbol == "MTAP",])
# MTAP from 21802543 to 21937651 (? Check if this is really longest)


# find_overlapping_coordinates(as.factor(linksFrom9p21$Peak1), "chr9:21,800,543-21,804,543")
overCutoff = linksFrom9p21[linksFrom9p21$coaccess > 0.015,]
overCutoff = overCutoff[!(is.na(overCutoff$Peak1)),]

# MTAP Promoter
find_overlapping_coordinates(as.factor(overCutoff$Peak2), "chr9:21,700,543-21,904,543")

overCutoff$Peak2 = as.character(overCutoff$Peak2)
# chr9_21801854_21803641 peak overlaps the MTAP promoter
overCutoff[overCutoff$Peak2 == "chr9_21801854_21803641",]
# > overCutoff[overCutoff$Peak2 == "chr9_21801854_21803641",]
#                           Peak1                  Peak2   coaccess
# 53811914 chr9_22079828_22080724 chr9_21801854_21803641 0.02058058
# 53812850 chr9_22103040_22103678 chr9_21801854_21803641 0.01562317
# 53812958 chr9_22103753_22104030 chr9_21801854_21803641 0.02581916

MTAP_Intersects =  linksFrom9p21[linksFrom9p21$Peak2 == "chr9_21801854_21803641",]
MTAP_Intersects[order(MTAP_Intersects$Peak2),]



# CDKN2A
find_overlapping_coordinates(as.factor(overCutoff$Peak2), "chr9:21,950,000-22,050,000")
# chr9_21994066_21996168 is a peak spanning TSS at 21,994,411
overCutoff[overCutoff$Peak2 == "chr9_21994066_21996168",]
# > overCutoff[overCutoff$Peak2 == "chr9_21994066_21996168",]
#                           Peak1                  Peak2   coaccess
# 53813305 chr9_22108729_22108933 chr9_21994066_21996168 0.02086954

cdkn2aIntersects =  linksFrom9p21[linksFrom9p21$Peak2 == "chr9_21994066_21996168",]
cdkn2aIntersects[order(cdkn2aIntersects$Peak2),]
# Interestingly, this finds a peak with a coaccessibility of -.054, and another with -.024. Negative regulators?
# at chr9_22103040_22103678 and chr9_22103753_22104030, respectively



# CDKN2B
find_overlapping_coordinates(as.factor(overCutoff$Peak2), "chr9:22,000,000-22,020,000")
# chr9_22008594_22009818 is a peak spanning TSS at 22,009,272
overCutoff[overCutoff$Peak2 == "chr9_22008594_22009818",]

# > overCutoff[overCutoff$Peak2 == "chr9_22008594_22009818",]
#                           Peak1                  Peak2   coaccess
# 53812770 chr9_22098749_22099337 chr9_22008594_22009818 0.02582328



cdkn2bIntersects =  linksFrom9p21[linksFrom9p21$Peak2 == "chr9_22008594_22009818",]
cdkn2bIntersects[order(cdkn2bIntersects$Peak2),]
# Also has a peak with -.0425 coaccessibility, at chr9_22124566_22125163 









