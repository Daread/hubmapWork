

library(monocle3)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

source("./deHelperFunctions.R")

library(DESeq2)
library(ggplot2)

library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-c", "--cdsRNA"), type="character", 
        default= "../2022_08_22_addNewSamples_nobackup/formattedData/NucleiOnlySharedGenesCDS.rds",   #"../2022_08_22_addNewSamples_nobackup/formattedData/allDatasetCDS.rds", 
              help="Name of RNA cds to use", metavar="character"), 
  make_option(c("-t", "--cellType"), type="character", 
        default="Endothelium", 
              help="Name of cell type to run DE tests", metavar="character"), 
  make_option(c("-v", "--variables"), type="character", 
        default= "Age,Sex,DataSource", #"Donor,BMI,Age,Sex,Diabetes,Hypertension", 
              help="Variables to get for each dataset", metavar="character"),
  make_option(c("-l", "--lvOnly"), type="logical", 
        default=TRUE, 
              help="Only use LV samples", metavar="logical"),
  make_option(c("-m","--minTPM"), type="numeric",
  		default = 5,
  			help ="Min reads per million to keep + test a gene", metavar = "numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

fixedVars = strsplit(opt$variables, ",", fixed=TRUE)[[1]]

# Read in data
cdsPath = opt$cdsRNA
fullCDS = readRDS(cdsPath)

processingNote = paste0("Pseudobulk_", opt$cellType)

# Keep the subset desired
fullCDS = fullCDS[,colData(fullCDS)$Cell_Shared_Label == opt$cellType]

# Only LV samples?
if (opt$lvOnly){
	processingNote = paste0(processingNote, "_LV_Only")
	fullCDS = fullCDS[,colData(fullCDS)$Anatomical_Site == "LV"]
} else{
	processingNote = paste0(processingNote, "_All_Sites")
}

# Get the pseudobulked data
pseudobulkList = getPseudobulkData(fullCDS, opt)

if (!(opt$lvOnly & (opt$variables == "Donor,BMI,Age,Sex,Diabetes,Hypertension,DataSource"))){
  # Write an output of organized metadata
  writeMetadata(pseudobulkList[[2]], opt)
}

str(pseudobulkList)

myFormula =  as.formula(paste0("~ ", paste(fixedVars, collapse = " + ")))

# Run the test
deInput = DESeqDataSetFromMatrix(countData = pseudobulkList[[1]], colData = as.data.frame(pseudobulkList[[2]]),
								design = myFormula)

deRes = DESeq(deInput)

resultsNames(deRes)

ageRes = results(deRes, name = "Age")
sexRes = results(deRes, name = "Sex_M_vs_F")

ageDF = as.data.frame(ageRes@listData)

sexRes = as.data.frame(sexRes@listData)

sexRes = sexRes[order(sexRes$pval),]






# Make some plots of pseudobulked values
genesToPlot = c("ACVRL1", "ALPL", "IL1R1", "PROM1", "LCN6", "CMKLR1", "IL1R1", "NRG3")
for (eachGene in genesToPlot){
  makeGeneAgeSexPlots(pseudobulkList, genesToPlot, processingNote)
}

















