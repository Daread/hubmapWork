

getDataset <- function(datasetToGet, opt){
	# Read et al
	if (datasetToGet == "Read"){
		return (readRDS(paste0("./read_et_al_data/", opt$cdsRNA)))
	} 
	# Litvinukova
	if (datasetToGet == "Litvinukova"){
		return (readRDS("./litvinukova_et_al_data/Litv_CDS.rds"))
	}
	# Tucker
	if (datasetToGet == "Tucker"){
		return (readRDS("./tuckerEtAlData/tuckerFormattedCDS.rds"))
	}
	# Reichart
	if (datasetToGet == "Reichart"){
		return(readRDS("./reichart_et_al_data/reichartCDS.rds"))
	}
	# Koenig
	if (datasetToGet == "Koenig"){
		data = readRDS("./koenig_et_al_data/koenigNucleiCDS.rds")
		return(data[,colData(data)$Condition == "Donor"])
	}
	# Chaffin
	if (datasetToGet == "Chaffin"){
		return(readRDS("./chaffin_et_al_data/chaffinDataCDS.rds"))
	}
}

getGeneNameVec <- function(opt){
	cds = getDataset("Read", opt)
	geneVec = rowData(cds)$gene_short_name
	names(geneVec) = rownames(cds)
	return(geneVec)
}

############################################################################### Read Begin
# Get the metadata
readMetaFile = "./read_et_al_data/Read_Donor_Metadata_DFformat.csv"
readMeta = read.csv(readMetaFile)

readDonor = function(x){
	return(x["Donor"])
}

readAge = function(x){
	return(x["Age"])
}

readSex = function(x){
	return(x["Sex"])
}

readLog10umi = function(x){
	return( as.numeric(x[["log10_umi"]]) )
}

readAnatomical_Site = function(x){
	if (as.character(x["Anatomical_Site"]) == "Left_Vent"){
		return("LV")
	} else if (as.character(x["Anatomical_Site"]) == "Right_Vent"){
		return("RV")
	} else {
		return(x["Anatomical_Site"])
	}
}

diabetesList = as.character(readMeta$History.of.diabetes. )
names(diabetesList) = as.character(readMeta$id)
readDiabetes = function(x, donorStatus = diabetesList){
	return(ifelse(donorStatus[as.character(x["Donor"])] =="NO", "N", "Y" ))
}

hypertensionList = as.character(readMeta$History.of.hypertension.)
names(hypertensionList) = as.character(readMeta$id)
readHypertension = function(x, donorStatus = hypertensionList){
	return(ifelse(donorStatus[as.character(x["Donor"])] =="NO", "N", "Y" ))
}

# # [1]                           
#  [4]                 
#  [7]                          
# [10]                   
readTypeVec =          c("Adipocytes", "Endothelium",           "Fibroblast", "T_Cell", "B_Cell", "Myeloid",    "Neuron",  "Perivascular",      "Atrial_Cardiomyocytes", "Ventricular_Cardiomyocytes", "Mast_Cells", "Endothelium",           "Endothelium", "Lymphatic_Endothelium")
names(readTypeVec) =   c( "Adipocytes", "Vascular_Endothelium", "Fibroblast", "T_Cell", "B_Cell", "Macrophage", "Neuronal", "VSM_and_Pericyte", "Atrial_Cardiomyocytes", "Cardiomyocyte",              "Mast_Cell",   "Vascular_Endothelium", "Endocardium", "Lymphatic_Endothelium" )

readCell_Shared_Label = function(x, cellLabels = readTypeVec){
	return(readTypeVec[as.character(x["highLevelCellType"])])
}

bmiList = as.character(readMeta$Body.Mass.Index..BMI..)
names(bmiList) = as.character(readMeta$id)
readBMI = function(x, bmiLabels = bmiList){
	return(as.numeric(strsplit(bmiLabels[as.character(x["Donor"])], " ", fixed=TRUE)[[1]][1]))
}

readDatatype = function(x){
	return("Nuclei")
}

getReadFormatting <- function(covarToSetup){

	thisList = vector(mode="list", length = 0)

	# Get all the setup with a list of functions
	thisList = c(thisList, "Donor" = readDonor)
	thisList = c(thisList, "Age" = readAge)
	thisList = c(thisList, "Sex" = readSex)
	thisList = c(thisList, "log10_umi" = readLog10umi)
	thisList = c(thisList, "Anatomical_Site" = readAnatomical_Site)
	thisList = c(thisList, "Diabetes" = readDiabetes)
	thisList = c(thisList, "Hypertension" = readHypertension)
	thisList = c(thisList, "Cell_Shared_Label" = readCell_Shared_Label)
	thisList = c(thisList, "BMI" = readBMI)
	thisList = c(thisList, "DataType" = readDatatype)

	# Return this list
	return(thisList)
}
############################################################################### Read end

############################################################################### Litvinukova begin
litvMetaFile = "./litvinukova_et_al_data/Litv_Donor_Data_Formatted_DF.csv"
litvMeta = read.csv(litvMetaFile)

litvDonor = function(x){
	return(x["donor"])
}

litvAge = function(x){
	return(x["Age"])
}

litvSex = function(x){
	return( substr(as.character(x["Sex"]),1,1) )
}

litvLog10umi = function(x){
	return(log10(as.numeric(x["n_counts"])) )
}

litvAnatomical_Site = function(x){
	charSite = as.character(x["region"])
	if (charSite %in% c("LV", "RV", "LA", "RA")){
		return(charSite)
	} else if (charSite == "AX"){
		return("Apex")
	} else{
		return ("Septum")
	}
}

diabetesList = as.character(litvMeta$Diabetes)
names(diabetesList) = litvMeta$Donor
litvDiabetes = function(x, donorStatus = diabetesList){
	return(ifelse(donorStatus[x["Donor"]] =="N", "N", "Y" ))
}

hypertensionList = as.character(litvMeta$Hypertension)
names(hypertensionList) = litvMeta$Donor
litvHypertension = function(x, donorStatus = hypertensionList){
	return(ifelse(donorStatus[x["Donor"]] =="N", "N", "Y" ))
}

# No separate vascular or endocardium assignments here
# coreTypeVec = = c("Adipocytes", "Endothelium", "Fibroblast", "T_Cell","B_Cell", "Myeloid", "Neuron", "Perivascular", "Atrial_Cardiomyocytes", "Ventricular_Cardiomyocytes", "Mast_Cells", "Unassigned", "Doublets")
litTypeVec =         c("Adipocytes", "Endothelium", "Fibroblast", "T_Cell",   "B_Cell", "Myeloid", "Neuron", "Perivascular", "Perivascular",     "Atrial_Cardiomyocytes", "Ventricular_Cardiomyocytes", "Mast_Cells", "Endothelium", "Unassigned", "Doublets", "Lymphatic_Endothelium")
# Note: No B Cell cluster in Litv or lymphatic endothelium. No separate call for endocardium vs. vascular endothelium
names(litTypeVec) = c( "Adipocytes", "Endothelial", "Fibroblast", "Lymphoid", "B_Cell", "Myeloid", "Neuronal", "Pericytes", "Smooth_muscle_cells", "Atrial_Cardiomyocyte", "Ventricular_Cardiomyocyte", "Mast_Cells", "Mesothelial", "NotAssigned", "doublets", "Lymphatic_Endothelium")
# > levels(colData(thisDataset)$cell_type)
#  [1] "Adipocytes"                "Atrial_Cardiomyocyte"
#  [3] "Endothelial"               "Fibroblast"
#  [5] "Lymphoid"                  "Mesothelial"
#  [7] "Myeloid"                   "Neuronal"
#  [9] "NotAssigned"               "Pericytes"
# [11] "Smooth_muscle_cells"       "Ventricular_Cardiomyocyte"
# [13] "doublets"
litvCell_Shared_Label = function(x, litTypeMap=litTypeVec){
	return(litTypeMap[as.character(x["cell_type"])])
}

bmiList = as.character(litvMeta$BMI)
names(bmiList) = litvMeta$Donor
litvBMI = function(x, donorBMI = bmiList ){
	bmiHighLow = as.numeric(strsplit(donorBMI[x["Donor"]], "-", fixed=TRUE)[[1]])
	return((bmiHighLow[1] + bmiHighLow[2]) / 2.0)
}

litvDataType = function(x){
	if (as.character(x["cell_source"]) %in% c("Harvard-Nuclei", "Sanger-Nuclei")){
		return ("Nuclei")
	} else {
		return ("Cells")
	}
}

getLitvFormatting <- function(covarToSetup){
	thisList = vector(mode="list", length = 0)

	# Get all the setup with a list of functions
	thisList = c(thisList, "Donor" = litvDonor)
	thisList = c(thisList, "Age" = litvAge)
	thisList = c(thisList, "Sex" = litvSex)
	thisList = c(thisList, "log10_umi" = litvLog10umi)
	thisList = c(thisList, "Anatomical_Site" = litvAnatomical_Site)
	thisList = c(thisList, "Diabetes" = litvDiabetes)
	thisList = c(thisList, "Hypertension" = litvHypertension)
	thisList = c(thisList, "Cell_Shared_Label" = litvCell_Shared_Label)
	thisList = c(thisList, "BMI" = litvBMI)
	thisList = c(thisList, "DataType" = litvDataType)

	# Return this list
	return(thisList)
}

############################################################################### Litvinukova end

############################################################################### Koenig Start
koenigMetaFile = "./koenig_et_al_data/Koenig_Patient_Metadata_Nuclei.csv"
koenigMeta = read.csv(koenigMetaFile)

koenigDonor = function(x){
	return(x["orig.ident"])
}

koenigAge = function(x){
	return(x["Age"])
}

koenigSex = function(x){
	return( toupper(substr(as.character(x["Sex"]),1,1)) )
}

koenigLog10umi = function(x){
	return(log10(as.numeric(x["nCount_RNA"])))
}

koenigAnatomical_Site = function(x){
	return("LV")
}

koenigDiabetesList = koenigMeta$Diabetes 
names(koenigDiabetesList) = koenigMeta$Samples
koenigDiabetes = function(x, diabetesList=koenigDiabetesList){
	return(ifelse(as.numeric(diabetesList[x["orig.ident"]]) == 0, "N", "Y" ))
}

koenigHypertensionList = koenigMeta$HTN
names(koenigHypertensionList) = koenigMeta$Samples
koenigHypertension = function(x, hyperList = koenigHypertensionList){
	return(ifelse(as.numeric(hyperList[x["orig.ident"]]) == 0, "N", "Y" ))
}

koenigTypeVec =         c("Adipocytes", "Endothelium", "Fibroblast", "T_Cell",   "B_Cell",     "Myeloid", "Neuron", "Perivascular", "Perivascular",     "Atrial_Cardiomyocytes", "Ventricular_Cardiomyocytes", "Mast_Cells", "Endothelium", "Epicardium", "Myeloid", "Lymphatic_Endothelium")
# Note: No B Cell cluster in Litv or lymphatic endothelium. No separate call for endocardium vs. vascular endothelium
names(koenigTypeVec) = c( "Adipocytes", "Endothelium", "Fibroblasts", "T/NK_Cells", "B_Cells", "Macrophages", "Neurons", "Smooth_Muscle", "Pericytes", "Atrial_Cardiomyocytes",   "Cardiomyocytes",            "Mast_Cells", "Endocardium", "Epicardium", "Monocytes", "Lymphatic" )
# [1] "Fibroblasts"    "Cardiomyocytes" "Lymphatic"      "Neurons"
#  [5] "Endothelium"    "Epicardium"     "Mast_Cells"     "Adipocytes"
#  [9] "Macrophages"    "Monocytes"      "Pericytes"      "B_Cells"
# [13] "Endocardium"    "Smooth_Muscle"  "T/NK_Cells"

koenigCell_Shared_Label = function(x, celltypeVec = koenigTypeVec){
	return(celltypeVec[x["Names"]])
}

koenigBMIlist = koenigMeta$BMI 
names(koenigBMIlist) = koenigMeta$Samles 
koenigBMI = function(x, bmiList = koenigBMIlist){
	return(as.numeric(bmiList[x["orig.ident"]]))
}

koenigDataType = function(x){
	return("Nuclei")
}

getKoenigFormatting <- function(covarToSetup){
	thisList = vector(mode="list", length = 0)

	# Get all the setup with a list of functions
	thisList = c(thisList, "Donor" = koenigDonor)
	thisList = c(thisList, "Age" = koenigAge)
	thisList = c(thisList, "Sex" = koenigSex)
	thisList = c(thisList, "log10_umi" = koenigLog10umi)
	thisList = c(thisList, "Anatomical_Site" = koenigAnatomical_Site)
	thisList = c(thisList, "Diabetes" = koenigDiabetes)
	thisList = c(thisList, "Hypertension" = koenigHypertension)
	thisList = c(thisList, "Cell_Shared_Label" = koenigCell_Shared_Label)
	thisList = c(thisList, "BMI" = koenigBMI)
	thisList = c(thisList, "DataType" = koenigDataType)

	# Return this list
	return(thisList)
}

############################################################################### Koenig end

############################################################################### Tucker Start

tuckerMetaFile = "tuckerEtAlData/Tucker_Donor_Metadata.csv"
tuckerMeta = read.csv(tuckerMetaFile)

tuckerDonor = function(x){
	return(as.character(x["Donor"]))
}

tuckerAge = function(x){
	return(x["Age"])
}

tuckerSex = function(x){
	return(substr(as.character(x["Sex"]),1,1))
}

tuckerLog10umi = function(x){
	return(log10(as.numeric(x["nUMI"])))
}

tuckerAnatomical_Site = function(x){
	return(x["chamber"])
}

tuckerDiabetes = function(x){
	# No data
	return("UNKNOWN")
}

tuckerHypertension = function(x){
	# No data
	return("UNKNOWN")
}

# In the UMAP (Fig 1B of Tucker et al) it looks like cyto cardio II groups with ventricular cardiomyocytes, while cyto carido I groups with atrial
tuckerTypeVec =         c("Adipocytes",  "Endothelium",         "Endothelium",          "Fibroblast",       "Fibroblast",           "Fibroblast",          "T_Cell",   "B_Cell",     "Myeloid",         "Neuron",   "Perivascular",             "Perivascular",       "Atrial_Cardiomyocytes",  "Ventricular_Cardiomyocytes",         "Ventricular_Cardiomyocytes",        "Ventricular_Cardiomyocytes",          "Ventricular_Cardiomyocytes",      "Atrial_Cardiomyocytes")
# Note: No B Cell cluster in Litv or lymphatic endothelium. No separate call for endocardium vs. vascular endothelium
names(tuckerTypeVec) = c("11. Adipocyte", "09. Endothelium I", "10. Endothelium II", "01. Fibroblast I", "02. Fibroblast II", "14. Fibroblast III",  "17. Lymphocyte", "B_Cells", "08. Macrophage", "16. Neuronal", "07. Pericyte", "13. Vascular Smooth Muscle", "03. Atrial Cardiomyocyte", "04. Ventricular Cardiomyocyte I",  "06. Ventricular Cardiomyocyte II",  "15. Ventricular Cardiomyocyte III",  "12. Cytoplasmic Cardiomyocyte II", "05. Cytoplasmic Cardiomyocyte I" )
# [1] "01. Fibroblast I"                  "02. Fibroblast II"
#  [3] "03. Atrial Cardiomyocyte"          "04. Ventricular Cardiomyocyte I"
#  [5] "05. Cytoplasmic Cardiomyocyte I"   "06. Ventricular Cardiomyocyte II"
#  [7] "07. Pericyte"                      "08. Macrophage"
#  [9] "09. Endothelium I"                 "10. Endothelium II"
# [11] "11. Adipocyte"                     "12. Cytoplasmic Cardiomyocyte II"
# [13] "13. Vascular Smooth Muscle"        "14. Fibroblast III"
# [15] "15. Ventricular Cardiomyocyte III" "16. Neuronal"
# [17] "17. Lymphocyte"

tuckerCell_Shared_Label = function(x, celltypeVec = tuckerTypeVec){
	return(celltypeVec[as.character(x["Cluster"])])
}

tuckerBMIlist = (tuckerMeta$Weight_kg/ ((tuckerMeta$Height_cm/100.0)^2))
names(tuckerBMIlist) = tuckerMeta$Donor 

tuckerBMI = function(x, bmiList = tuckerBMIlist){
	return(bmiList[x["Donor"]])
}

tuckerDataType = function(x){
	return("Nuclei")
}

gettuckerFormatting <- function(covarToSetup){
	thisList = vector(mode="list", length = 0)

	# Get all the setup with a list of functions
	thisList = c(thisList, "Donor" = tuckerDonor)
	thisList = c(thisList, "Age" = tuckerAge)
	thisList = c(thisList, "Sex" = tuckerSex)
	thisList = c(thisList, "log10_umi" = tuckerLog10umi)
	thisList = c(thisList, "Anatomical_Site" = tuckerAnatomical_Site)
	thisList = c(thisList, "Diabetes" = tuckerDiabetes)
	thisList = c(thisList, "Hypertension" = tuckerHypertension)
	thisList = c(thisList, "Cell_Shared_Label" = tuckerCell_Shared_Label)
	thisList = c(thisList, "BMI" = tuckerBMI)
	thisList = c(thisList, "DataType" = tuckerDataType)

	# Return this list
	return(thisList)
}


############################################################################### Tucker end

############################################################################### Reichart Start

reichartMetaFile = "./reichart_et_al_data/S1_Clinical_Metadata_Patient_Information.csv"
reichartMeta = read.csv(reichartMetaFile)

reichartDonor = function(x){
	return(as.character(x["Patient"]))
}

reichartAge = function(x){
	return(x["Age"])
}

reichartSex = function(x){
	return( toupper(substr(as.character(x["sex"]), 1,1)) )
}

reichartLog10umi = function(x){
	return(log10(as.numeric(x["nCounts_RNA"])))
}

reichartAnatomical_Site = function(x){
	return(as.character(x["Region_x"]))
}

reichartDiabetsList = reichartMeta$Diabetes 
names(reichartDiabetsList) = reichartMeta$Donor
reichartDiabetes = function(x, diabetesList = reichartDiabetsList){
	#return(diabetesList[as.character(x["Patient"])])
	return("UNKNOWN")
}

reichartHypertensionList = reichartMeta$Hypertension
names(reichartHypertensionList) = reichartMeta$Donor 
reichartHypertension = function(x, hyperList = reichartHypertensionList){
	# return(hyperList[as.character(x["Patient"])])
	return("UNKNOWN")
}


reichartTypeVec =         c("Adipocytes",  "Endothelium",             "Fibroblast",              "Ventricular_Cardiomyocytes",  "T_Cell",    "Myeloid",         "Neuron",       "Perivascular", "Mast_Cells", "Native_Cell")
# Note: No B Cell cluster in Litv or lymphatic endothelium. No separate call for endocardium vs. vascular endothelium
names(reichartTypeVec) = c( "fat cell",    "endothelial cell",  "fibroblast of cardiac tissue", "cardiac muscle cell",        "lymphocyte", "myeloid cell",    "cardiac neuron", "mural cell",  "mast cell", "native cell" )
 # [1] "cardiac muscle cell"          "cardiac neuron"
 # [3] "endothelial cell"             "fat cell"
 # [5] "fibroblast of cardiac tissue" "lymphocyte"
 # [7] "mast cell"                    "mural cell"
 # [9] "myeloid cell"                 "native cell"

reichartCell_Shared_Label = function(x, celltypeVec = reichartTypeVec){
	return(celltypeVec[as.character(x["cell_type"])])
}

# reichartBMIlist = as.numeric(reichartMeta$Weight_Kg) / ((as.numeric(reichartMeta$Height_cm) / 100.0)^2) 
# names(reichartBMIlist) = reichartMeta$Donor 
# reichartBMI = function(x, bmiList = reichartBMIlist){
reichartBMI = function(x){
	return(x["BMI"])
}

reichartDataType = function(x){
	return("Nuclei")
}

getReichartFormatting <- function(covarToSetup){
	thisList = vector(mode="list", length = 0)

	# Get all the setup with a list of functions
	thisList = c(thisList, "Donor" = reichartDonor)
	thisList = c(thisList, "Age" = reichartAge)
	thisList = c(thisList, "Sex" = reichartSex)
	thisList = c(thisList, "log10_umi" = reichartLog10umi)
	thisList = c(thisList, "Anatomical_Site" = reichartAnatomical_Site)
	thisList = c(thisList, "Diabetes" = reichartDiabetes)
	thisList = c(thisList, "Hypertension" = reichartHypertension)
	thisList = c(thisList, "Cell_Shared_Label" = reichartCell_Shared_Label)
	thisList = c(thisList, "BMI" = reichartBMI)
	thisList = c(thisList, "DataType" = reichartDataType)

	# Return this list
	return(thisList)
}


############################################################################### Reichart End

############################################################################### Chaffin Start

chaffinMetaFile = "./chaffin_et_al_data/Formatted_Chaffin_Control_Metadata.csv"
chaffinMeta = read.csv(chaffinMetaFile)

chaffinDonor = function(x){
	return(as.character(x["donor_id"]))
}

chaffinAge = function(x){
	return(x["age"])
}

chaffinSex = function(x){
	return( toupper(substr(as.character(x["sex"]), 1, 1)) )
}

chaffinLog10umi = function(x){
	return(log10(as.numeric(x["cellbender_ncount"])))
}

# Not in h5ad file, but all are LV samples
chaffinAnatomical_Site = function(x){
	return("LV")
}

chaffinDiabetes = function(x){
	return("UNKNOWN")
}

chaffinHypertension = function(x){
	return("UNKNOWN")
}

chaffinCelltypeList =        c("Adipocytes", "Endothelium",     "Endothelium",  "Endothelium",   "Endothelium",      "Fibroblast",   "Fibroblast",      "Fibroblast",   "Ventricular_Cardiomyocytes","Ventricular_Cardiomyocytes","Ventricular_Cardiomyocytes",  "Myeloid",           "Myeloid",              "T_Cell",    "Neuron", "Perivascular", "Perivascular", "Perivascular", "Mast_Cells", "Lymphatic_Endothelium", "Epicardial")
names(chaffinCelltypeList) = c("Adipocyte",  "Endocardial", "Endothelial_I", "Endothelial_II", "Endothelial_III", "Fibroblast_I", "Fibroblast_II", "Activated_fibroblast", "Cardiomyocyte_I",           "Cardiomyocyte_II",           "Cardiomyocyte_III",       "Macrophage", "Proliferating_macrophage",  "Lymphocyte", "Neuronal", "Pericyte_I",  "Pericyte_II",      "VSMC",     "Mast_cell",  "Lymphatic_endothelial", "Epicardial"  )
# [1] "Activated_fibroblast"     "Adipocyte"
#  [3] "Cardiomyocyte_I"          "Cardiomyocyte_II"
#  [5] "Cardiomyocyte_III"        "Endocardial"
#  [7] "Endothelial_I"            "Endothelial_II"
#  [9] "Endothelial_III"          "Epicardial"
# [11] "Fibroblast_I"             "Fibroblast_II"
# [13] "Lymphatic_endothelial"    "Lymphocyte"
# [15] "Macrophage"               "Mast_cell"
# [17] "Neuronal"                 "Pericyte_I"
# [19] "Pericyte_II"              "Proliferating_macrophage"
# [21] "VSMC"

chaffinCell_Shared_Label = function(x, celltypeVec = chaffinCelltypeList){
	return(celltypeVec[as.character(x["cell_type_leiden0.6"])])
}

chaffinBMIlist = chaffinMeta$BMI 
names(chaffinBMIlist) = paste0("P", chaffinMeta$Individual)
chaffinBMI = function(x, bmiList = chaffinBMIlist){
	return(bmiList[as.character(x["donor_id"])])
}

chaffinDataType = function(x){
	return("Nuclei")
}


getChaffinFormatting <- function(covarToSetup){
	thisList = vector(mode="list", length = 0)

	# Get all the setup with a list of functions
	thisList = c(thisList, "Donor" = chaffinDonor)
	thisList = c(thisList, "Age" = chaffinAge)
	thisList = c(thisList, "Sex" = chaffinSex)
	thisList = c(thisList, "log10_umi" = chaffinLog10umi)
	thisList = c(thisList, "Anatomical_Site" = chaffinAnatomical_Site)
	thisList = c(thisList, "Diabetes" = chaffinDiabetes)
	thisList = c(thisList, "Hypertension" = chaffinHypertension)
	thisList = c(thisList, "Cell_Shared_Label" = chaffinCell_Shared_Label)
	thisList = c(thisList, "BMI" = chaffinBMI)
	thisList = c(thisList, "DataType" = chaffinDataType)

	# Return this list
	return(thisList)
}


############################################################################### Chaffin End

# Overall helper function to set up a list of lists, each of which contains named functions for helping
#    format the cds coldata entries for each input dataset
getFormatterList <- function(datasetToSetup, covarToSetup){

	# Set up an empty list
	formattingList = vector(mode="list", length = 0)

	# Add entries
	formattingList$Read = getReadFormatting(covarToSetup) #c(formattingList, "Read" = getReadFormatting(covarToSetup))
	formattingList$Litvinukova = getLitvFormatting(covarToSetup)
	formattingList$Koenig = getKoenigFormatting(covarToSetup)
	formattingList$Tucker = gettuckerFormatting(covarToSetup)
	formattingList$Reichart = getReichartFormatting(covarToSetup)
	formattingList$Chaffin = getChaffinFormatting(covarToSetup)

	# Return the ones we need for this use case
	# return(formattingList[names(formattingList) %in% datasetToSetup])
	return(formattingList)

}


formatRowdata <- function(inputCDS, datasetName, genesReadDataOrder){
	# Read et al
	if (datasetName == "Read"){
		# Remove extra rows
		rowData(inputCDS) = rowData(inputCDS)[c("gene_short_name")]
		return(inputCDS)
	} 
	# Litvinukova
	if (datasetName == "Litvinukova"){
		# rowData(inputCDS)$gene_short_name = rownames(inputCDS)
		rownames(inputCDS) = rowData(inputCDS)$gene_id

		genesInDataset = genesReadDataOrder[names(genesReadDataOrder) %in% rownames(inputCDS)]

		# Keep matching rows, in order matching Read et al
		inputCDS = inputCDS[names(genesInDataset),]

		# Get short names, matching the naming conventions used in the Read CDS
		rowData(inputCDS)$gene_short_name = genesInDataset
		
		# Remove extra rows
		rowData(inputCDS) = rowData(inputCDS)[c("gene_short_name")]
		return(inputCDS)
	}
	# Tucker
	if (datasetName == "Tucker"){
		# rowData(inputCDS)$gene_short_name = rownames(inputCDS)
		rownames(inputCDS) = rowData(inputCDS)$gene_id

		genesInDataset = genesReadDataOrder[names(genesReadDataOrder) %in% rownames(inputCDS)]

		# Keep matching rows, in order matching Read et al
		inputCDS = inputCDS[names(genesInDataset),]

		# Get short names, matching the naming conventions used in the Read CDS
		rowData(inputCDS)$gene_short_name = genesInDataset
		
		# Remove extra rows
		rowData(inputCDS) = rowData(inputCDS)[c("gene_short_name")]
		return(inputCDS)
	}
	# Reichart
	if (datasetName == "Reichart"){
		# Assign as would for the rest:
		# rownames(inputCDS) = rowData(inputCDS)$gene_id

		genesInDataset = genesReadDataOrder[names(genesReadDataOrder) %in% rownames(inputCDS)]

		# Keep matching rows, in order matching Read et al
		inputCDS = inputCDS[names(genesInDataset),]

		# Get short names, matching the naming conventions used in the Read CDS
		rowData(inputCDS)$gene_short_name = genesInDataset
		
		# Remove extra rows
		rowData(inputCDS) = rowData(inputCDS)[c("gene_short_name")]
		return(inputCDS)
	}
	# Koenig
	if (datasetName == "Koenig"){
		# browser()
		# First, need to use Tucker et al names to get a mapping of ensembl IDs for the short names in use here
		tuckerTempCDS = getDataset("Tucker", "EMPTY_OPT")
		tuckerShortNames = rowData(tuckerTempCDS)$gene_id
		names(tuckerShortNames) = rownames(tuckerTempCDS)
		tuckerTempCDS = NULL

		# Only keep reichart entries that match a gene short name in tucker
		inputCDS = inputCDS[rownames(inputCDS) %in% names(tuckerShortNames),]

		tuckerNameMatch = tuckerShortNames[names(tuckerShortNames) %in% rownames(inputCDS)]

		# Add a field with the ensembl IDs
		inputCDS = inputCDS[names(tuckerNameMatch),]
		rowData(inputCDS)$gene_id = tuckerNameMatch

		# Add as would for the rest
		rownames(inputCDS) = rowData(inputCDS)$gene_id
		genesInDataset = genesReadDataOrder[names(genesReadDataOrder) %in% rownames(inputCDS)]

		# Keep matching rows, in order matching Read et al
		inputCDS = inputCDS[names(genesInDataset),]

		# Get short names, matching the naming conventions used in the Read CDS
		rowData(inputCDS)$gene_short_name = genesInDataset
		
		# Remove extra rows
		rowData(inputCDS) = rowData(inputCDS)[c("gene_short_name")]
		return(inputCDS)
	}
	# Chaffin
	if (datasetName == "Chaffin"){
		# rowData(inputCDS)$gene_short_name = rownames(inputCDS)
		rownames(inputCDS) = rowData(inputCDS)$gene_id

		genesInDataset = genesReadDataOrder[names(genesReadDataOrder) %in% rownames(inputCDS)]

		# Keep matching rows, in order matching Read et al
		inputCDS = inputCDS[names(genesInDataset),]

		# Get short names, matching the naming conventions used in the Read CDS
		rowData(inputCDS)$gene_short_name = genesInDataset
		
		# Remove extra rows
		rowData(inputCDS) = rowData(inputCDS)[c("gene_short_name")]
		return(inputCDS)
	} else{
		print("No row data conversion method found")
		return(inputCDS)
	}
}



