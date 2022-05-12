# Note: This only runs if you get the right resources upon qlogin/qsub. Run:
#  qlogin -q trapnell-login.q -l mfree=4G -pe serial 16
# to start up

# Trying:
# qlogin -l mfree=20G -q trapnell-login.q -pe serial 16 -l centos=7

# Trying this on 7-20-21. (Didn't work, going back to 20g)
# qlogin -l mfree=40G -q trapnell-login.q -pe serial 16 -l centos=7

basepath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/"
dir.create(paste0(basepath, "archr/results/NB6/"))
out_dir = paste0(basepath, "archr/results/NB6/")
setwd(out_dir)

# load requirements
suppressPackageStartupMessages({
  library(ArchR)
  library(monocle3)
  library(dplyr)
  library(ggrastr)
  # library(Seurat)
})

# set genome
# addArchRGenome("hg19")
addArchRGenome("hg38")

# addArchRThreads(threads = 32) 
addArchRThreads(threads = 16) 

# load filtered ArchR project
prj = loadArchRProject(path = paste0(basepath, "archr/Heart_filtered"),
                       showLogo = FALSE)

########################
# Convert ATAC data to cds object

# load PeakMatrix into memory
bMat = getMatrixFromProject(
  ArchRProj = prj, 
  useMatrix = "TileMatrix", 
  binarize = TRUE, 
  threads = getArchRThreads())

# format rowData 
rd = data.frame(rowData(bMat)) %>% 
  dplyr::mutate(stop = start + 499,
                bin = paste(seqnames, start, stop, sep = "_")) %>% 
  dplyr::select(bin)
row.names(rd) <- rd$bin
row.names(bMat) <- rd$bin

# Create CDS from tile Matrix (SummarizedExperiment)
cds_b = monocle3::new_cell_data_set(
  assays(bMat)$TileMatrix, 
  cell_metadata = colData(bMat),
  gene_metadata = rd)


# load gene score Matrix into memory
gMat = getMatrixFromProject(
  ArchRProj = prj, 
  useMatrix = "GeneScoreMatrix", 
  binarize = FALSE, 
  threads = getArchRThreads())

# format rowData 
rd = data.frame(rowData(gMat)) %>% 
  dplyr::mutate(bin = paste(seqnames, start, end, sep = "_"), 
                gene_short_name = name) %>% 
  dplyr::select(bin, name, gene_short_name)
row.names(rd) <- rd$name
row.names(gMat) <- rd$name

# Create CDS from gene score Matrix (SummarizedExperiment)
cds_g = monocle3::new_cell_data_set(
  assays(gMat)$GeneScoreMatrix, 
  cell_metadata = colData(gMat),
  gene_metadata = rd)


# Added 7-29-21: Filter by checking names in a cds filtered by greg/riza's method
filterByGregCDS = TRUE
filterCDSName = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/cds_objects/"
processingNote = "FRIP=0.1_FRIT=0.08UMI=1000DL=0.5" #"MNN_sample"
# processingNote = "FRIP=0.3_FRIT=0.1UMI=1000"
# filterByDL = TRUE
# DLcutoff = .5


backupG = cds_g 
backupB = cds_b

if (filterByGregCDS){

  library(tidyr)
  # thisCDS object holds data
  load(paste0(filterCDSName, "cds_p_allHeartATACafterUMAP", processingNote, "MNN_sample"))

  sampleNameDF = separate(data=as.data.frame(colData(thisCDS)), col = "sampleName", 
            into = c("sampleName", "procNote"), sep="FRIP", remove=TRUE)
  colData(thisCDS)$sampleName = sampleNameDF$sampleName 
  sampleNameDF = NULL

  namesToSave = paste0(colData(thisCDS)$sampleName,"#", (colData(thisCDS)$cell))
  cds_b = cds_b[,(rownames(colData(cds_b)) %in% namesToSave) ]
  cds_g = cds_g[,(rownames(colData(cds_g)) %in% namesToSave) ]
  colnames(thisCDS) = namesToSave

  # Need to go the other way as well
  thisCDS = thisCDS[,colnames(thisCDS) %in% colnames(cds_g)]
  namesToSave = paste0(colData(thisCDS)$sampleName,"#", (colData(thisCDS)$cell))
} else{
  processingNote = ""
}

# Format so that the cells are in the same order and use the same syntax
# colnames(thisCDS) = namesToSave

cds_b = cds_b[,namesToSave]
cds_g = cds_g[,namesToSave]


# Keep bins with a basline level of expression
cds.atac.b = cds_b
cds.atac.b <- detect_genes(cds.atac.b)
ncells = ncol(cds.atac.b)*0.002
cds.atac.b = cds.atac.b[rowData(cds.atac.b)$num_cells_expressed > ncells,]

cds.atac.g = cds_g


# Output to the cds_
filteredCDSpath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_03_redoATAC_proc_nobackup/cds_objects/"
setwd(filteredCDSpath)

# save these cds objects for future use (saves time)
saveRDS(thisCDS,    file=paste0("HM10_all_heart_fullFilter_", processingNote, "_ATAC_cds_p.RDS"))
saveRDS(cds.atac.g, file=paste0("HM10_all_heart_fullFilter_", processingNote, "_ATAC_cds_g.RDS"))
saveRDS(cds.atac.b, file=paste0("HM10_all_heart_fullFilter_", processingNote, "_ATAC_cds_b.RDS"))





