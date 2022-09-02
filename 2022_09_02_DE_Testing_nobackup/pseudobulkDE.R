

library(monocle3)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

source("./deHelperFunctions.R")

library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-c", "--cdsRNA"), type="character", 
        default="../2022_08_22_addNewSamples_nobackup/formattedData/allDatasetCDS.rds", 
              help="Name of RNA cds to use", metavar="character"), 
  make_option(c("-t", "--cellType"), type="character", 
        default="Myeloid", 
              help="Name of cell type to run DE tests", metavar="character"), 
  make_option(c("-v", "--variables"), type="character", 
        default="Donor,BMI,Age,Sex,Diabetes,Hypertension,Anatomical_Site", 
              help="Variables to get for each dataset", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Read in data
cdsPath = opt$cdsRNA
fullCDS = readRDS(cdsPath)

# Keep the subset desired
fullCDS = fullCDS[,colData(fullCDS)$Cell_Shared_Label == opt$cellType]

# Get the pseudobulked data
pseudobulkList = getPseudobulkData(fullCDS, opt)

str(pseudobulkList)






