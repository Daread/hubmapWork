
library(seurat)
library(monocle3)
library(SeuratObject)
#remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)


# Read in the seurat object for koenig et al:
koenigPath = "./koenig_et_al_data/koenig_nuclei.Robj"
load(koenigPath)
# object "nuclei" is the seurat object holding all nuclei
koenigNucleiCDS = as.cell_data_set(nuclei)

# Add information on donor attributes using donor metadata
donorDataPath = "./koenig_et_al_data/Koenig_Patient_Metadata_Nuclei.csv"
nucleiDonorMeta = read.csv(donorDataPath)

# Remove unwanted data from colData and add in metadata related to donor characteristics
allCol = colnames(colData(koenigNucleiCDS))
allCol = allCol[!grepl("SCT", allCol)]
colData(koenigNucleiCDS) = colData(koenigNucleiCDS)[colnames(colData(koenigNucleiCDS)) %in% allCol]

# Add the desired donor metadata
colToAdd = c("Age", "Sex")
for (eachCol in colToAdd){
	dataVec = nucleiDonorMeta[[eachCol]]
	names(dataVec) = nucleiDonorMeta[["Samples"]]

	colData(koenigNucleiCDS)[eachCol] = dataVec[colData(koenigNucleiCDS)$orig.ident]
	names(colData(koenigNucleiCDS)[[eachCol]]) = NULL
}


# rdsOutput = "./formattedData/koenigNucleiCDS.rds"

rdsOutput = "./koenig_et_al_data/koenigNucleiCDS.rds"
saveRDS(koenigNucleiCDS, file=rdsOutput)


# Get the seurat rds object
reichartPath = "./reichart_et_al_data/reichartSeurat.rds"
reichartSeurat = readRDS(reichartPath)

reichartCDS = as.cell_data_set(reichartSeurat)

# Get metadata
reichartMetadataFile = "./reichart_et_al_data/S1_Clinical_Metadata_Patient_Information.csv"
reichartMeta = read.csv(reichartMetadataFile)

# Keep only data from unique donors (not repeat data from litvinukova et al)
reichartUniqueDonors = c("H46", "H49", "H51", "H53", "H55", "H67")
reichartCDS = reichartCDS[, as.character(colData(reichartCDS)$Patient) %in% reichartUniqueDonors]

# Reset levels
colData(reichartCDS)$Patient = as.factor(as.character(colData(reichartCDS)$Patient))

# Add age
colToAdd = c("Age")
for (eachCol in colToAdd){
	dataVec = reichartMeta[[eachCol]]
	names(dataVec) = reichartMeta[["Donor"]]

	colData(reichartCDS)[eachCol] = dataVec[as.character(colData(reichartCDS)$Patient)]
	names(colData(reichartCDS)[[eachCol]]) = NULL
}


# Save output as a cds
cdsOutput = "./reichart_et_al_data/reichartCDS.rds"
saveRDS(reichartCDS, file=cdsOutput)












