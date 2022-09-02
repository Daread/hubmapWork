










getPseudobulkData <- function(inputCDS, opt){

	# Track with covariates to put into a DF. Manually set in opt. 
	covarToGet = strsplit(opt$variables, ",")[[1]]

	colData(inputCDS)$biolSample = paste0(colData(inputCDS)$DataSource, "_", colData(inputCDS)$Donor, "_", colData(inputCDS)$Anatomical_Site)

	allSamples = levels(as.factor(colData(inputCDS)$biolSample))

	# Empty metadata DF
	metadataDF = data.frame(matrix(NA, nrow=length(allSamples), ncol = length(covarToGet)))
	colnames(metadataDF) = covarToGet
	rownames(metadataDF) = allSamples

	# Empty expression matrix
	bulkExprs = matrix(NA, nrow = nrow(inputCDS), ncol = length(allSamples))
	rownames(bulkExprs) = rownames(rowData(inputCDS)$gene_short_name)
	colnames(bulkExprs) = allSamples

	# Loop through and for each biological sample get:
	# 1: The pseudobulked counts
	# 2: Metadata added 
	for (eachBiolSample in allSamples){
		print(paste0("Getting pseudobulk data for ", eachBiolSample))
		# Get the cds subset for this source of cells
		miniCDS = inputCDS[,colData(inputCDS)$biolSample == eachBiolSample]

		# Manually I need to make sure I give inputs here that are shared across all cells in a biological sample.
		# Just grab the entries from the first row (cell) of coldata for all the covariates needed
		for (eachCovar in covarToGet){
			metadataDF[eachBiolSample, eachCovar] = colData(miniCDS)[[eachCovar]][1]
		}

		# Expression vector:
		exprsVec = rowSums(exprs(miniCDS))
		bulkExprs[,eachBiolSample] = exprsVec
	}

	# Return a list
	return(list(bulkExprs, metadataDF))
}






