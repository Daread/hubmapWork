

hardAssignAnatomicalSites <- function(inputCDS){
	colData(inputCDS)$Anatomical_Site = "Not_Specified"
	colData(inputCDS)$Anatomical_Site = ifelse(grepl("Apex", colData(inputCDS)$sample),
										"Apex", colData(inputCDS)$Anatomical_Site)
	colData(inputCDS)$Anatomical_Site = ifelse(grepl("Septum", colData(inputCDS)$sample),
										"Septum", colData(inputCDS)$Anatomical_Site)
	colData(inputCDS)$Anatomical_Site = ifelse(grepl("Left.Vent", colData(inputCDS)$sample),
										"Left_Vent", colData(inputCDS)$Anatomical_Site)
	colData(inputCDS)$Anatomical_Site = ifelse(grepl("Right.Vent", colData(inputCDS)$sample),
										"Right_Vent", colData(inputCDS)$Anatomical_Site)
	return(inputCDS)
}
hardAssignAnatomicalSitesDF <- function(inputDF){

	inputDF$Anatomical_Site = "Not_Specified"
	inputDF$Anatomical_Site = ifelse(grepl("Apex", inputDF$sample),
										"Apex", inputDF$Anatomical_Site)
	inputDF$Anatomical_Site = ifelse(grepl("Septum", inputDF$sample),
										"Septum", inputDF$Anatomical_Site)
	inputDF$Anatomical_Site = ifelse(grepl("Left.Vent", inputDF$sample),
										"Left_Vent", inputDF$Anatomical_Site)
	inputDF$Anatomical_Site = ifelse(grepl("Right.Vent", inputDF$sample),
										"Right_Vent", inputDF$Anatomical_Site)
	return(inputDF)
}
hardAssignDonors <- function(inputCDS){
	colData(inputCDS)$Donor = "Not_Specified"
	colData(inputCDS)$Donor = ifelse(grepl("W134", colData(inputCDS)$sample),
										"W134", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W135", colData(inputCDS)$sample),
										"W135", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W136", colData(inputCDS)$sample),
										"W136", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W139", colData(inputCDS)$sample),
										"W139", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W142", colData(inputCDS)$sample),
										"W142", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W144", colData(inputCDS)$sample),
										"W144", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W145", colData(inputCDS)$sample),
										"W145", colData(inputCDS)$Donor)
	colData(inputCDS)$Donor = ifelse(grepl("W146", colData(inputCDS)$sample),
										"W146", colData(inputCDS)$Donor)
	return(inputCDS)
}


hardAssignHighLevelCellTypes <- function(inputCDS, origProcessing){
	if (origProcessing == "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40"){
		# Assign
		colData(inputCDS)$highLevelCellType = "Unassigned"
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label == "1", 
						"Fibroblast", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("4"), 
						"Vascular_Endothelium", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("7"), 
						"Endocardium", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("8", "13"), 
						"Lymphatic_Endothelium", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("6"), 
						"VSM_and_Pericyte", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("2"), 
						"Macrophage", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("3"), 
						"Cardiomyocyte", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("5"), 
						"T_Cell", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("11"), 
						"Adipocytes", colData(inputCDS)$highLevelCellType)
		colData(inputCDS)$highLevelCellType = ifelse(colData(inputCDS)$partition_label %in% c("10"), 
						"Neuronal", colData(inputCDS)$highLevelCellType)
	}
	return(inputCDS)
}


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



hardAssignDonorAgesDF <- function(inputDF){
  inputDF$Age = 0
  inputDF$Age =ifelse(inputDF$Donor == "W134", 43, inputDF$Age)
  inputDF$Age =ifelse(inputDF$Donor == "W135", 60, inputDF$Age)
  inputDF$Age =ifelse(inputDF$Donor == "W136", 43, inputDF$Age)
  inputDF$Age =ifelse(inputDF$Donor == "W137", 49, inputDF$Age)
  inputDF$Age =ifelse(inputDF$Donor == "W139", 45, inputDF$Age)
  inputDF$Age =ifelse(inputDF$Donor == "W142", 55, inputDF$Age)
  inputDF$Age =ifelse(inputDF$Donor == "W144", 53, inputDF$Age)
  inputDF$Age =ifelse(inputDF$Donor == "W145", 51, inputDF$Age)
  inputDF$Age =ifelse(inputDF$Donor == "W146", 25, inputDF$Age)

  inputDF$Log10Age = log10(inputDF$Age)

  return(inputDF)
}


hardAssignDonorSexesDF <- function(inputDF){
  inputDF$Sex = "Not_Set"
  inputDF$Sex =ifelse(inputDF$Donor == "W134", "F", inputDF$Sex)
  inputDF$Sex =ifelse(inputDF$Donor == "W135", "M", inputDF$Sex)
  inputDF$Sex =ifelse(inputDF$Donor == "W136", "M", inputDF$Sex)
  inputDF$Sex =ifelse(inputDF$Donor == "W137", "F", inputDF$Sex)
  inputDF$Sex =ifelse(inputDF$Donor == "W139", "M", inputDF$Sex)
  inputDF$Sex =ifelse(inputDF$Donor == "W142", "F", inputDF$Sex)
  inputDF$Sex =ifelse(inputDF$Donor == "W144", "M", inputDF$Sex)
  inputDF$Sex =ifelse(inputDF$Donor == "W145", "M", inputDF$Sex)
  inputDF$Sex =ifelse(inputDF$Donor == "W146", "F", inputDF$Sex)

  return(inputDF)
}


# a help function to tidy the vgam model output - used in compare_abundance function
tidy.vglm = function(x, conf.int=FALSE, conf.level=0.95) {
  co <- as.data.frame(coef(summary(x)))
  names(co) <- c("estimate","std.error","statistic","p.value")
  if (conf.int) {
    qq <- qnorm((1+conf.level)/2)
    co <- transform(co,
                    conf.low=estimate-qq*std.error,
                    conf.high=estimate+qq*std.error)
  }
  co <- data.frame(term=rownames(co),co)
  rownames(co) <- NULL
  return(co)
}

