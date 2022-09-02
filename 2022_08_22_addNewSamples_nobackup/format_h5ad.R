

# Need to use cluster installations of seruatDisk
.libPaths('2')
library(Seurat)
library(SeuratDisk)
library(monocle3)
library(SeuratWrappers)

getLitvAge = function(x){
	ages = as.numeric(strsplit(x["age_group"], "-", fixed=TRUE)[[1]] )
	ageMidpoint = (ages[1] + ages[2]) / 2.0
	return(ageMidpoint)
}


# Get litvinukova data
###################################################
litvinukovaPath = "./litvinukova_et_al_data/litvinukova_data.h5ad"
litvinukovaSeuratPath = "./litvinukova_et_al_data/litvinukova_data.h5seurat"
Convert(litvinukovaPath, dest = litvinukovaSeuratPath, overwrite=FALSE)

litvSeurat =  LoadH5Seurat(litvinukovaSeuratPath)

# Get a monocle CDS
litvCDS = as.cell_data_set(litvSeurat)

# Get the donor metadata
litvMetadata = "./litvinukova_et_al_data/Litv_Supp_Table_1_Donor_Metadata.csv"
litvDonorData = as.data.frame(t(read.csv(litvMetadata)))
litvDonorData$Donor = rownames(litvDonorData) #[2:length(rownames(litvDonorData))]

formattedLitvDonorData = litvDonorData
colnames(formattedLitvDonorData) = litvDonorData[1,]
formattedLitvDonorData = formattedLitvDonorData[2:nrow(formattedLitvDonorData),]

# Save formatted metadata
metadataOutput = "./litvinukova_et_al_data/Litv_Donor_Data_Formatted_DF.csv"
write.csv(formattedLitvDonorData, file=metadataOutput)
# saveRDS(formattedLitvDonorData, file=metadataOutput)

# Format to have "Age" and "Sex" columns (Check into this, is this actually gender or bio sex?)
colData(litvCDS)$Sex = colData(litvCDS)$gender
colData(litvCDS)$Age = apply(as.data.frame(colData(litvCDS)), 1, getLitvAge)

# Get gene IDs added
rowData(litvCDS)$gene_id = litvSeurat$RNA@meta.features[["gene_ids-Harvard-Nuclei"]]


# Save
litvCDSout = "./litvinukova_et_al_data/Litv_CDS.rds"
saveRDS(litvCDS, file=litvCDSout)

###########################################

# Get Tucker et al data
##########################################

tuckerPath = "./tuckerEtAlData/healthy_human_4chamber_map_unnormalized_V4.h5ad"
tuckerSeuratPath = "./tuckerEtAlData/tucker_data_V4.h5seurat"
Convert(tuckerPath, dest = tuckerSeuratPath, overwrite=FALSE)

tuckerSeurat = LoadH5Seurat(tuckerSeuratPath)

# Get a monocle CDS
tuckerCDS = as.cell_data_set(tuckerSeurat)

# Get metadata and format 
tuckerMetadataPath = "./tuckerEtAlData/Tucker_Donor_Metadata.csv"
tuckerMeta = read.csv(tuckerMetadataPath)

# Add Sex and Age (and other fields, if desired)
colData(tuckerCDS)$Donor = paste0("P", as.character(colData(tuckerCDS)$biological.individual))


# Add the desired donor metadata
colToAdd = c("Age", "Sex")
for (eachCol in colToAdd){
	dataVec = tuckerMeta[[eachCol]]
	names(dataVec) = tuckerMeta[["Donor"]]

	colData(tuckerCDS)[eachCol] = dataVec[colData(tuckerCDS)$Donor]
	names(colData(tuckerCDS)[[eachCol]]) = NULL
}

# Get Gene names
geneNameFile = "./tuckerEtAlData/genes_v2.tsv"
geneNameDF = as.data.frame(read.table(geneNameFile, sep='\t', header=FALSE))

colnames(geneNameDF) = c("Ensembl_ID", "Short_Name")
nameToID = geneNameDF$Ensembl_ID
names(nameToID) = geneNameDF$Short_Name

rowData(tuckerCDS)$gene_id = nameToID[rownames(tuckerCDS)]
names(rowData(tuckerCDS)$gene_id) = NULL

formattedTuckerPath = "./tuckerEtAlData/tuckerFormattedCDS.rds"
saveRDS(tuckerCDS, file = formattedTuckerPath)

############################################



# Get Chaffin et al data
##########################################

chaffinPath = "./chaffin_et_al_data/human_dcm_hcm_scportal_03.17.2022.h5ad"
chaffinSeuratPath = "./chaffin_et_al_data/chaffinData.h5seurat"
Convert(chaffinPath, dest = chaffinSeuratPath, overwrite = FALSE)

chaffinSeurat = LoadH5Seurat(chaffinSeuratPath)

# Get into a CDS
chaffinCDS = as.cell_data_set(chaffinSeurat)


# Only keep the non-failing hearts
colData(chaffinCDS)$disease = as.character(colData(chaffinCDS)$disease)
chaffinCDS = chaffinCDS[,colData(chaffinCDS)$disease == "NF"]


# Get Gene names
geneNameFile = "./chaffin_et_al_data/DCM_HCM_Expression_Matrix_genes_V1.tsv"
geneNameDF = as.data.frame(read.table(geneNameFile, sep='\t', header=FALSE))

colnames(geneNameDF) = c("Ensembl_ID", "Short_Name")
nameToID = geneNameDF$Ensembl_ID
names(nameToID) = geneNameDF$Short_Name

rowData(chaffinCDS)$gene_id = nameToID[rownames(chaffinCDS)]
names(rowData(chaffinCDS)$gene_id) = NULL

# Save the output
chaffinOutPath = "./chaffin_et_al_data/chaffinDataCDS.rds"
saveRDS(chaffinCDS, file = chaffinOutPath)




##########################################





