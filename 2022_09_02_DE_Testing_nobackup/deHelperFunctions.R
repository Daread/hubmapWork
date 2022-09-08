










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
	rownames(bulkExprs) = paste0(rownames(inputCDS), "_", rowData(inputCDS)$gene_short_name)
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

	# Format the metadata to be numeric if needed
	numericVars = c("Age", "BMI")
	for (eachCol in colnames(metadataDF)){
		if (eachCol %in% numericVars){
			metadataDF[[eachCol]] = as.numeric(metadataDF[[eachCol]])
		}
	}

	# Now only keep the entries in the exprs matrix that are above the minimum TPM value
	tpmVec = rowSums(bulkExprs) * 1e6 / sum(bulkExprs)
	# browser()
	bulkExprs = bulkExprs[tpmVec > opt$minTPM,]

	# Return a list
	return(list(bulkExprs, metadataDF))
}


writeMetadata <- function(metadataDF, opt){
	dir.create("./fileOutputs/")
	outFile = paste0("./fileOutputs/", "AllSampleMetadata.csv" )
	write.csv(metadataDF, file = outFile)
}


makeGeneAgeSexPlots <- function(pseudobulkList, genesToPlot, processingNote){
	genePlotDir = "./plots/"
	dir.create(genePlotDir)

	# browser()
	exprsMat = pseudobulkList[[1]]
	rownames(exprsMat) = vapply(strsplit(rownames(exprsMat), "_"), "[", 2, FUN.VALUE=character(1)) #      strsplit(rownames(exprsMat), "_", fixed = TRUE)
	metaDF = pseudobulkList[[2]]
	# Loop and get a plot for each gene requested
	for (eachGene in genesToPlot){
		miniDF = metaDF
		miniDF$expressionFrac = exprsMat[eachGene,] / colSums(exprsMat)
		# Make a plot by age...

		png(paste0(genePlotDir, processingNote, eachGene, "_vs_Age.png"))
		myPlot = ggplot(miniDF, aes_string(x="Age", y="expressionFrac", col="DataSource")) +
				geom_point()
		print(myPlot)
		dev.off()

		# ... and by sex
		png(paste0(genePlotDir, processingNote, eachGene, "_vs_Sex.png"))
		myPlot = ggplot(miniDF, aes_string(x="Sex", y="expressionFrac")) +
				 geom_boxplot() + geom_point()
		print(myPlot)
		dev.off()

	}
}





































