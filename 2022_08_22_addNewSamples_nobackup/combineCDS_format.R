
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
        default="Donor,BMI,Age,Sex,Diabetes,Hypertension,Anatomical_Site,log10_umi,Cell_Shared_Label,DataType", 
              help="Variables to get for each dataset", metavar="character")

)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

covarToGet = strsplit(opt$variables, ",")[[1]]

# datasetVec = c("Read", "Litvinukova", "Tucker", "Koenig", "Reichart", "Chaffin")
datasetVec =  c("Read", "Litvinukova", "Tucker", "Koenig", "Reichart", "Chaffin") #c("Koenig", "Litvinukova") # c("Koenig", "Litvinukova", "Read")  #c("Read", "Litvinukova")

formattedGeneNames = getGeneNameVec(opt) # rownames(getDataset("Read", opt))

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

	colData(thisDataset) = colData(thisDataset)[,covarToGet]
	colData(thisDataset)$DataSource = eachDataset
	colData(thisDataset)$Instrument = ifelse(eachDataset == "Read", "sci", "10x")

	# Now format the rowdata and gene orders
	thisDataset = formatRowdata(thisDataset, eachDataset, formattedGeneNames)

	# Save this to the cds list
	cdsList[[counterVal]] = thisDataset
	counterVal = counterVal + 1
}

# Get a combined CDS
combinedCDS = combine_cds(cdsList)

# Save the combined output
combinedOutfile = "./formattedData/allDatasetCDS.rds"
saveRDS(combinedCDS, file=combinedOutfile)






# Looking at the number of genes that are widely shared across datasets vs. only used in a small number:
for (eachInd in 1:6){
	print(nrow(cdsList[[eachInd]][rownames(cdsList[[eachInd]]) %in% rownames(cdsList[[1]]), ]))
}

# See how many times various genes occur
geneUseCount = rep(1, length=nrow(cdsList[[1]]))
names(geneUseCount) = rownames(cdsList[[1]])

for (eachInd in 2:6){
	print(paste0("Dataset ", as.character(eachInd)))
	for (eachGeneInd in 1:nrow(cdsList[[eachInd]])){
		thisGene = rownames(cdsList[[eachInd]])[eachGeneInd]
		geneUseCount[thisGene] = 1 + geneUseCount[thisGene]
	}
}

# Save this info on gene use
geneUseFile = "./formattedData/geneUseVec.rds"
saveRDS(geneUseCount, file=geneUseFile)








# # Check overlaps of naming
# print(nrow(cdsList[[4]]))
# # Read match?
# print("Read:")
# print(sum(rownames(cdsList[[4]]) %in% rowData(cdsList[[1]])$gene_short_name))
# # Lit match
# print("Lit")
# print(sum(rownames(cdsList[[4]]) %in% rownames(cdsList[[2]])))
# # Tucker
# print("Tucker")
# print(sum(rownames(cdsList[[4]]) %in% rownames(cdsList[[3]])))

# # Reichart
# print("Reichart")
# print(sum(rownames(cdsList[[4]]) %in% as.character(rowData(cdsList[[5]])$feature_name)))
# # Chaffin
# print("Chaffin")
# print(sum(rownames(cdsList[[4]]) %in% rownames(cdsList[[6]])))



# rownames(cdsList[[4]])[!(rownames(cdsList[[4]]) %in% rownames(cdsList[[3]])) ]
