# Greg Booth 2021
# This script creates an ArchR project by combining all arrow files from HEART 
# It then filters cells based on the filtered cds objects previously created 
# for each sample

basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/"
dir.create(paste0(basepath, "archr/results/NB4/"))
out_dir = paste0(basepath, "archr/results/NB4/")
setwd(paste0(basepath, "archr"))

# load requirements
suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(dplyr)
})

# set genome
# addArchRGenome("hg19")
addArchRGenome("hg38")

addArchRThreads(threads = 32) 

# pointer to arrow files
# HeartArrowFiles = c(
#   "W134.heart.apex.arrow.s1","W135.heart.LV.arrow.s1", "W136.heart.apex.arrow.s1", 
#   "W136.heart.LV.arrow.s1",  "W137.heart.apex.arrow.s1", "W139.heart.apex.arrow.s1", 
#   "W139.heart.LV.arrow.s1", "W139.heart.RV.arrow.s1",  "W139.heart.septum.arrow.s1", 
#   "W142.heart.LV.arrow.s1", "W144.heart.apex.arrow.s1", "W145.heart.apex.arrow.s1", 
#   "W145.heart.LV.arrow.s1", "W146.heart.apex.arrow.s1", "W146.heart.LV.arrow.s1")

HeartArrowFiles =  paste0( c(
  "W134.heart.apex","W135.heart.LV", "W136.heart.apex", 
  "W136.heart.LV",  "W137.heart.apex", "W139.heart.apex", 
  "W139.heart.LV", "W139.heart.RV",  "W139.heart.septum", 
  "W142.heart.LV", "W144.heart.apex", "W145.heart.apex", 
  "W145.heart.LV", "W146.heart.apex", "W146.heart.LV"), ".s1.arrow")

prjHeart <- ArchRProject(
  ArrowFiles = HeartArrowFiles, 
  outputDirectory = "Heart_filtered",
  copyArrows = FALSE 
)

# create a filtered cell whitelist from filtered cds objects
samples = c(
  "W134.heart.apex.s1","W135.heart.LV.s1", "W136.heart.apex.s1", 
  "W136.heart.LV.s1", "W137.heart.apex.s1", "W139.heart.apex.s1", 
  "W139.heart.LV.s1", "W139.heart.RV.s1", "W139.heart.septum.s1", 
  "W142.heart.LV.s1", "W144.heart.apex.s1", "W145.heart.apex.s1", 
  "W145.heart.LV.s1", "W146.heart.apex.s1", "W146.heart.LV.s1")

cell_whitelist= c()

for(s in samples){
  load(paste0(basepath, "cds_objects/cds_p_", s))
  cdat = data.frame(colData(cds_p))
  cells = cdat %>% 
    mutate(cells = paste0(s, "#", cell)) %>% 
    select(cells)
  cell_whitelist = c(cell_whitelist, cells$cells)
}

# filter ArchR project based on whitelist
prjHeart_f = prjHeart[prjHeart$cellNames %in% cell_whitelist,]

saveArchRProject(
  ArchRProj = prjHeart_f,
  outputDirectory = getOutputDirectory(prjHeart_f),
  overwrite = TRUE,
  load = TRUE,
  dropCells = TRUE,
  logFile = createLogFile("saveArchRProject"),
  threads = getArchRThreads()
)
















  
  