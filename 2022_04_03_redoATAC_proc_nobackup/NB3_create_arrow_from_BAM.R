# on cluster, initiate qlogin session with 16 cores. 
#  qlogin -q trapnell-login.q -l mfree=4G -pe serial 16
# start R (4.0.0)

basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/"
out_path = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/"
out_dir = paste0(out_path, "archr/")
dir.create(out_dir)
setwd(out_dir)

suppressPackageStartupMessages({
  library(ArchR)
})

set.seed(1)

inputFiles = c()
for(s in c("W134.heart.apex.s1", 
              "W135.heart.LV.s1", 
              "W136.heart.apex.s1", "W136.heart.LV.s1",
              "W137.heart.apex.s1", # Failed in RNA
              "W139.heart.apex.s1", "W139.heart.LV.s1", "W139.heart.RV.s1", "W139.heart.septum.s1", 

              "W142.heart.LV.s1",
               "W144.heart.apex.s1", 
               "W145.heart.apex.s1", "W145.heart.LV.s1",
              "W146.heart.apex.s1", "W146.heart.LV.s1") ){
  inputFiles[s] = paste0(basepath, s, "/merge_bams/",  s, "-merged.bam")
}

inputFiles

addArchRGenome("hg38")
# create arrow file 

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 2, #Dont set this too high because you can always increase later
  filterFrags = 500, 
  bcTag = "qname",
  gsubExpression = ":.*", #If BamFile : Cell Barcode is everything after ":"
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE
)
# note that the output gives qc_metrix for each input sample 
