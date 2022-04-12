
# Greg Booth 2021
# This script takes the output from the BBI pipeline for sciATAC 
# for each sample
# performs some filtering (by QC and doublets)
# returns a cds object of peak by cell data. 

basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/"
out_path = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/"
blacklist = read.table("~/../gtb7/genomes/GRCh38/annotations/EncodeBlacklist/ENCFF356LFX.bed")
dir.create(paste0(out_path, "cds_objects/"))
out_dir = paste0(out_path, "cds_objects/")

# load requirements
# suppressPackageStartupMessages({
#   library(monocle3)
# })


.libPaths("2")


# DFR Mod: Load my packages
source("../../../sharedProjectCode/utility/singleCellUtilFuncs.R")

# modStatus <- loadMonoclePackages()
library(monocle3)


library(Seurat)



source("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/misc/atac_helper_functions.R")


# sample_df = read.csv(file = paste0(basepath, "merged_plots/merged.called_cells_summary.stats.csv"), 
#                      header = TRUE, sep = "\t") %>% 
#   select(samples = sample, med_FRIP = median_per_cell_frip, med_FRIT = median_per_cell_frit) %>% 
#   mutate(FRIP_cutoff = ifelse(med_FRIP - 0.05 < 0.1, 0.1, med_FRIP - 0.05), 
#          FRIT_cutoff = ifelse(med_FRIT - 0.025 < 0.05, 0.05, med_FRIT - 0.025))


FRIP_cutoffToUse = .15
FRIT_cutoffToUse = .05
umiCutoff = 1000
# DLcutoff=.7
processingNote = paste0("FRIP=", as.character(FRIP_cutoffToUse), "_FRIT=", as.character(FRIT_cutoffToUse),
                        "UMI=", as.character(umiCutoff) ) #, "DL=", as.character(DLcutoff))
print("Processing Note")


sample_df = data.frame(
  samples = c("W134.heart.apex.s1", 
              "W135.heart.LV.s1", 
              "W136.heart.apex.s1", "W136.heart.LV.s1",
              "W137.heart.apex.s1", # Failed in RNA
              "W139.heart.apex.s1", "W139.heart.LV.s1", "W139.heart.RV.s1", "W139.heart.septum.s1", 

              "W142.heart.LV.s1",
               "W144.heart.apex.s1", 
               "W145.heart.apex.s1", "W145.heart.LV.s1",
              "W146.heart.apex.s1", "W146.heart.LV.s1"),
  
  FRIP_cutoff = c(0.1, 0.1, 0.1, 0.1, 0.1, 
                  0.1, 0.1, 0.1, 0.1, 0.1,
                  0.1, 0.1, 0.1, 0.1, 0.1),
  
  FRIT_cutoff = c(0.05, 0.05, 0.05, 0.05, 0.05,
                  0.05, 0.05, 0.05, 0.05, 0.05,
                  0.05, 0.05, 0.05, 0.05, 0.05))


# Take a look at QC output
plotDir = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/plots/QC_Plots/"
dir.create(plotDir)



i = 0
for (row in 1:nrow(sample_df)){
  print(paste0("Working on row ", as.character(row)))

  s = sample_df[row, "samples"]
  FRIP_co  = FRIP_cutoffToUse # sample_df[row, "FRIP_cutoff"]
  FRIT_co =  FRIT_cutoffToUse #sample_df[row, "FRIT_cutoff"]
  i = i + 1
  # load peak matrix for sample. 
  pMat = load_mtx_file(paste0(basepath, s, "/make_matrices/",  s,"-peak_matrix.mtx.gz"))
  
  # add "chr" to beginning of chromosome names 
  rn = row.names(pMat)
  rn_fix = paste0(prefix = "chr", rn)
  row.names(pMat) <- rn_fix
  
  #binarize peak matrix
  pMat = binarize_matrix(pMat)
  
  cat("starting number of cells in sample ", s, " = ", dim(pMat)[2], "\n")
  cat("starting number of features in sample ", s, " = ", dim(pMat)[1], "\n")
  
  ######################################################################################
  # filter cells 
  ######################################################################################
  # load summary stats for adding to col Data and reorder to match Cell order in matrices
  cDat = read.table(paste0(basepath, s, "/count_report/",  s,"-count_report.txt"), head = T)
  cDat_f = cDat[match(colnames(pMat), cDat$cell),]
  row.names(cDat_f) <- cDat_f$cell
  colnames(cDat_f) <- c("cell", "total", "umi", "RIP", "RIT")
  cDat_f$FRIP = cDat_f$RIP/cDat_f$umi
  cDat_f$FRIT= cDat_f$RIT/cDat_f$umi

# # for (eachRow)
#   traitsToPlot = c("umi", "FRIP", "FRIT")
#   for (eachTrait in traitsToPlot){
#     thisOutput = paste0(plotDir, s, "_", eachTrait, ".png")
#     png(thisOutput, width=1000, height=1000, res = 200)
#     myPlot = ggplot(cDat_f, aes_string(eachTrait)) + 
#           geom_histogram()
#     print(myPlot)
#     dev.off()
#   }
  
  # filter clls based on Unique reads, FRIP and FRIT  (cutoffs set above)
  qc_cells = filter(cDat_f, umi > umiCutoff, FRIP > FRIP_co, FRIT > FRIT_co) %>% select(cell)
  pMat = pMat[,colnames(pMat) %in% qc_cells$cell]
  
  cat("number of cells after UMI, FRIP and FRIT cutoffs ", s, " = ", dim(pMat)[2], "\n")
  ######################################################################################
  # filter features
  ######################################################################################
  # remove outlier features (z-score based) .
  pMat = filter_features_z(bmat = pMat, lower_bound=-2, upper_bound=4, downsample=NULL)
  
  # the filter_features function from the atac_helper script removes features that overlap blackout regions of the genome. 
  # load blacklist
  colnames(blacklist) <- c('chrom', 'start', 'end') 
  features_f = filter_regions(features = row.names(pMat), blacklist_df = blacklist)
  pMat = pMat[row.names(pMat) %in% features_f,]
  
  cat("number of features after blacklist and outlier removal = ", dim(pMat)[1], "\n")
  ######################################################################################
  # run scrublet
  ######################################################################################
  scrub = atac_scrublet(bmat = pMat, k=NULL, fraction_sim_doublets=1, estimated_doublet_rate=0.1, dims=2:50)

  scrub_res = filter(scrub, simulated_doublet == FALSE)
  threshold = quantile(scrub_res$doublet_score, .9)
  scrub_res$doublet = sapply(scrub_res$doublet_score, function(x){
    ifelse(x < threshold, "singlet", "doublet")})
  
  cDat_f = inner_join(cDat_f, scrub_res, by = "cell")

  # # 8-6-21 addition of this plot
  # thisOutput = paste0(plotDir, s, "_", "doublet_likelihood", ".png")
  # png(thisOutput, width=1000, height=1000, res = 200)
  # myPlot = ggplot(cDat_f, aes_string("doublet_likelihood")) + 
  #       geom_histogram()
  # print(myPlot)
  # dev.off()

  cDat_f = cDat_f[match(colnames(pMat), cDat_f$cell),]
  row.names(cDat_f) <- cDat_f$cell
  # Filter by doublet likelihood
  # cDat_f = cDat_f[cDat_f$doublet_likelihood < DLcutoff,]
  #####################################################################################
  # Make CDS object for filtered cells and features. 
  cds_p = monocle3::new_cell_data_set(pMat, cell_metadata = cDat_f)
  cds_p <- monocle3::detect_genes(cds_p, min_expr = 0)
  save(cds_p, file = paste0(out_dir, "cds_p_", s, processingNote))
}



