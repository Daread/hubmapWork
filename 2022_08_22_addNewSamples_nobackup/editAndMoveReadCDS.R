
library(monocle3)


hardAssignDonorAges <- function(inputCDS){
  colData(inputCDS)$Age = 0
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W134", 43, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W135", 60, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W136", 43, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W137", 49, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W139", 45, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W142", 55, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W144", 53, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W145", 51, colData(inputCDS)$Age)
  colData(inputCDS)$Age =ifelse(colData(inputCDS)$Donor == "W146", 25, colData(inputCDS)$Age)

  colData(inputCDS)$Log10Age = log10(colData(inputCDS)$Age)

  return(inputCDS)
}


hardAssignDonorSexes <- function(inputCDS){
  colData(inputCDS)$Sex = "Not_Set"
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W134", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W135", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W136", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W137", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W139", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W142", "F", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W144", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W145", "M", colData(inputCDS)$Sex)
  colData(inputCDS)$Sex =ifelse(colData(inputCDS)$Donor == "W146", "F", colData(inputCDS)$Sex)

  return(inputCDS)
}



source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
source("./formattingFunctions.R")

library("optparse")
# Get the passed parameters
option_list = list(

  make_option(c("-c", "--cdsRNA"), type="character", 
        default="allCells_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes.rds", 
              help="Name of RNA cds to use", metavar="character"),
  make_option(c("-x", "--rnaPath"), type="character", 
        default="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/", 
              help="Path to RNA cds", metavar="character"),

  make_option(c("-v", "--variables"), type="character", 
        default="Age,Sex,Diabetes,Hypertension,Anatomical_Site,log10umi", 
              help="Variables to get for each dataset", metavar="character")

)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

readCDS = readRDS(paste0(opt$rnaPath, opt$cdsRNA))

# 
readCDS = hardAssignDonorSexes(readCDS)
readCDS = hardAssignDonorAges(readCDS)

# save the output
outfile = paste0("./read_et_al_data/", opt$cdsRNA)
saveRDS(readCDS, file=outfile)

