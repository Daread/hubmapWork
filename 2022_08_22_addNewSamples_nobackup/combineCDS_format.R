
library(monocle3)

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
        default="Donor,BMI,Age,Sex,Diabetes,Hypertension,Anatomical_Site,log10_umi,Cell_Shared_Label", 
              help="Variables to get for each dataset", metavar="character")

)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

covarToGet = strsplit(opt$variables, ",")[[1]]

# datasetVec = c("Read", "Litvinukova", "Tucker", "Koenig", "Reichart")
datasetVec = c("Reichart", "Tucker") # c("Koenig", "Litvinukova", "Read")  #c("Read", "Litvinukova")

formatterList = getFormatterList(datasetVec, covarToGet)

cdsList = vector(mode="list", length = length(datasetVec))
counterVal = 1
for (eachDataset in datasetVec){
	print(paste0("Working on ", eachDataset))
	# Read in the CDS for this dataset
	thisDataset = getDataset(eachDataset, opt)

	# Now we need to format one column for each entry in covarToGet
	for (eachCovar in covarToGet){
		print(paste0("Working on ", eachCovar))
		# Apply the function named for this covariate, in entry in the list per dataset
		colData(thisDataset)[[eachCovar]] = apply(as.data.frame(colData(thisDataset)), 1,
												 formatterList[[eachDataset]][[eachCovar]])
		names(colData(thisDataset)[[eachCovar]]) = NULL
	}

	# Save this to the cds list
	cdsList[[counterVal]] = thisDataset
	counterVal = counterVal + 1
}









