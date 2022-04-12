

library(tidyr)
library(ggplot2)

library("optparse")
# Get the passed parameters
option_list = list(

  make_option(c("-p", "--pValFIMOcutoff"), type="numeric", 
  			default=1e-6, 
              help="Max FIMO p value to retain match", metavar="character"),
  make_option(c("-b", "--binarizeMotifs"), type="logical", 
        default=TRUE, 
              help="Binarize motif presence/absence, rather than track counts", metavar="logical"),
  make_option(c("-c", "--combineAllSeq"), type="logical", 
        default=TRUE, 
              help="Merge distal and promoter motif counts", metavar="logical"),
  make_option(c("-u", "--promoterUpstream"), type="numeric", 
        default=1000,
              help="Bases upstream of TSS to use for input features", metavar="numeric"),
  make_option(c("-d", "--promoterDownstream"), type="numeric", 
        default=200,
              help="Bases downstream of TSS to use for input features", metavar="numeric"),

  make_option(c("-k", "--coaccessCutoff"), type="numeric", 
        default=.05,
              help="Cutoff for keeping coaccessible peaks", metavar="numeric"),

  make_option(c("-n", "--maxNdistalSites"), type="numeric", 
        default=5,
              help="Max number of sites to link to a gene's promoter", metavar="numeric"),

  make_option(c("-q", "--peakSize"), type="numeric", 
        default=600,
              help="-1 -> Keep as is, or enter a positive integer to make all peaks the same size", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)




# Get previous output
dfName = paste0("Gene_Prom_Plus_Distal_WithSequence_Sites_Max", opt$maxNdistalSites, "_Upstream", opt$promoterUpstream,
                "_Downstream", opt$promoterDownstream, "_cicCuf", opt$coaccessCutoff,
                    "peakSize", opt$peakSize)
promoterDistalDFWithSeqs = readRDS(paste0("./rdsOutputs/", dfName, ".RDS"))

# Need to load the FIMO results
FIMOpath = paste0("../2022_04_07_redoFIMO_nobackup/", dfName, ".fa/")
fimoFile = paste0(FIMOpath, "fimo.tsv")
if (file.exists(fimoFile)){
  fimoDF = read.csv(fimoFile, sep="\t")
} else {
  FIMOpath = paste0("../2022_04_07_redoFIMO_nobackup/", dfName, ".RDS.fa/")
  fimoDF = read.csv(paste0(FIMOpath, "fimo.tsv"), sep="\t")
}


fimoDF = fimoDF[fimoDF$p.value < opt$pValFIMOcutoff,]

# Gene_Prom_Plus_Distal_WithSequence_Sites_Max5_Upstream1000_Downstream200_cicCuf0.05peakSize600.RDS.fa/



getFullMotifMatrix <- function(fimoDF, promDF, opt){

  # Get the names of all motifs to work with
  motifNames = as.character(levels(as.factor(fimoDF$motif_alt_id)))

  # Set up a new DF
  motifDF = promDF[,c("GeneID", "chr", "orientation", "linkedSites", "Model_Set")]
  rownames(motifDF) = motifDF$GeneID
  # Promoter sequences
  for (eachMotif in motifNames){
    motifDF[[paste0(eachMotif, "_Prom")]] = 0
  }
  for (eachMotif in motifNames){
    motifDF[[paste0(eachMotif, "_Dist")]] = 0
  }

  print("Getting motifs match matrix now")
  # Loop through the fimoDF. Update each appropriate entry
  for (eachInd in 1:nrow(fimoDF)){
    if (!(grepl(":", fimoDF[eachInd, "sequence_name"],fixed=TRUE) )) {
      # browser()
      next
    }

    # Count every 100k
    if ((eachInd %% 100000) == 0){
      print(eachInd)
    }

    # Get the name
    seqName = strsplit(fimoDF[eachInd, "sequence_name"], ":")[[1]]
    geneHit = seqName[1]
    siteHit = seqName[2] # Tells if a distal hit or promoter
    motifHit = fimoDF[eachInd, "motif_alt_id"] # Which motif found a match

    if (siteHit == "PromoterPos"){
      siteSuffix = "_Prom"
    } else {
      siteSuffix = "_Dist"
    }

    # Now we have the position exactly. Increment.
    motifDF[geneHit, paste0(motifHit, siteSuffix)] = motifDF[geneHit, paste0(motifHit, siteSuffix)] + 1

  }

  # Return the df
  return(motifDF)

}

fimoResMatrix = getFullMotifMatrix(fimoDF, promoterDistalDFWithSeqs, opt)



# Set up outputs
outDir = paste0("./fileOutputs/", dfName, "/")
dir.create(outDir)





getBinarizedData <- function(fimoResMatrix, motifNames){
  for (eachMotif in motifNames){
    fimoResMatrix[[paste0(eachMotif, "_Prom")]] = ifelse(fimoResMatrix[[paste0(eachMotif, "_Prom")]] > 0, 1, 0)
    fimoResMatrix[[paste0(eachMotif, "_Dist")]] = ifelse(fimoResMatrix[[paste0(eachMotif, "_Dist")]] > 0, 1, 0)
  }

  return(fimoResMatrix)
}

getCombinedData <- function(fimoResMatrix, motifNames){
  for (eachMotif in motifNames){
    fimoResMatrix[[paste0(eachMotif, "_AllSeq")]] = (fimoResMatrix[[paste0(eachMotif, "_Prom")]] +
                                                       fimoResMatrix[[paste0(eachMotif, "_Dist")]]     )
  }
  # Drop the others
  fimoResMatrix = fimoResMatrix[, !(colnames(fimoResMatrix) %in% paste0(motifNames, "_Prom"))]
  fimoResMatrix = fimoResMatrix[, !(colnames(fimoResMatrix) %in% paste0(motifNames, "_Dist"))]

  return(fimoResMatrix)
}
 
 binarizePooledMat <- function(unbinarizedPooled, motifNames){
  # Binarize all
  for (eachMotif in motifNames){
    unbinarizedPooled[[paste0(eachMotif, "_AllSeq")]] = ifelse(unbinarizedPooled[[paste0(eachMotif, "_AllSeq")]] > 0, 1, 0)
  }

  return (unbinarizedPooled)
 }


outputFormattedMatrices <- function(fimoResMatrix, outDir, dfName, motifNames, opt){


  # First, output as-is
  write.csv(fimoResMatrix, paste0(outDir, dfName, "pVal", 
                  as.character(opt$pValFIMOcutoff), "_Nonbinary_Uncombined_Motif_Counts.csv"))

  # Output the promoter only
  unbinarizedPromOnly = fimoResMatrix[, !(colnames(fimoResMatrix) %in% 
                                        paste0(motifNames, "_Dist"))]
  write.csv(unbinarizedPromOnly, paste0(outDir, dfName, "pVal", 
                  as.character(opt$pValFIMOcutoff), "_Nonbinary_PromOnly_Motif_Counts.csv"))

  # Now, get the version that's binarized
  binarizedUnpooled = getBinarizedData(fimoResMatrix, motifNames)

  write.csv(binarizedUnpooled, paste0(outDir, dfName, "pVal", 
                  as.character(opt$pValFIMOcutoff), "_Binary_Uncombined_Motif_Counts.csv"))

  binarizedPromOnly = binarizedUnpooled[, !(colnames(binarizedUnpooled) %in% 
                                        paste0(motifNames, "_Dist"))]
  write.csv(binarizedPromOnly, paste0(outDir, dfName, "pVal", 
                  as.character(opt$pValFIMOcutoff), "_Binary_PromOnly_Motif_Counts.csv"))

  # Now, get the version that's pooled (not binarized)
  unbinarizedPooled = getCombinedData(fimoResMatrix, motifNames)
  write.csv(unbinarizedPooled, paste0(outDir, dfName, "pVal", 
          as.character(opt$pValFIMOcutoff), "_Nonbinary_Combined_Motif_Counts.csv"))

  # Now, get the version that's pooled and then 
  binarizedPooled = binarizePooledMat(unbinarizedPooled, motifNames)
  write.csv(binarizedPooled, paste0(outDir, dfName, "pVal", 
          as.character(opt$pValFIMOcutoff), "_Binary_Combined_Motif_Counts.csv"))

}




motifNames = as.character(levels(as.factor(fimoDF$motif_alt_id)))

outputFormattedMatrices(fimoResMatrix, outDir, dfName, motifNames, opt)
















