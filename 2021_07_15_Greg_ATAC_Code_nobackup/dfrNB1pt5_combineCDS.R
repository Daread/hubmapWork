
basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
# dir.create(paste0(basepath, "archr/"))
# dir.create(paste0(basepath, "archr/results/"))
out_path = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
blacklist = read.table("~/../gtb7/genomes/GRCh38/annotations/EncodeBlacklist/ENCFF356LFX.bed")
dir.create(paste0(out_path, "cds_objects/"))
out_dir = paste0(out_path, "cds_objects/")
setwd(out_dir)

source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
modStatus <- loadMonoclePackages()

source("~/../gtb7/sciatac_pipeline/src/sylvia_scripts/scripts/atac_helper_functions.R")


## Load saved CDS

# samples = c("W134.heart.apex.s1", 
#               "W135.heart.LV.s1",
#               "W136.heart.apex.s1", "W136.heart.LV.s1")

samples = c("W134.heart.apex.s1", 
              "W135.heart.LV.s1", 
              "W136.heart.apex.s1", "W136.heart.LV.s1",
              "W137.heart.apex.s1", # Failed in RNA
              "W139.heart.apex.s1", "W139.heart.LV.s1", "W139.heart.RV.s1", "W139.heart.septum.s1", 

              "W142.heart.LV.s1",
               "W144.heart.apex.s1", 
               "W145.heart.apex.s1", "W145.heart.LV.s1",
              "W146.heart.apex.s1", "W146.heart.LV.s1")

# Add a processing note?
processingNote = "FRIP=0.3_FRIT=0.1UMI=1000"
samples = paste0(samples, processingNote)

# Read in the CDS list, combine into a single one
sampleCDSlist = vector(mode="list", length = length(samples))

for (eachInd in 1:length(samples)){
  s = samples[eachInd]
  print(paste0("working on ", s))
  # Get the cds
  loadOutput = load(paste0(basepath, "cds_objects/cds_p_", s))
  cds_p$sampleName = s
  sampleCDSlist[[eachInd]] = cds_p
}

# Combine
cds_p = combine_cds(sampleCDSlist)

cds_p = detect_genes(cds_p)

# Write
save(cds_p, file = paste0(out_dir, "cds_p_allHeartATAC", processingNote))


