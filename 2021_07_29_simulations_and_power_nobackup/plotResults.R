
library(stringr)
library(reshape2)
library(ggplot2)

print("Starting Now")

plotCurveFromRDSdataframes <- function(dfAsRDSfiles,groupTitle,
				colToPlot = c("posResultMM")){
	# read in the rds objects
	splitNames = str_split(dfAsRDSfiles, "_folds_")


	# Loop through. For each result make a separate plot
	for (eachResult in colToPlot){
		# Get the first DF:
		firstDF = readRDS(paste0("./rdsOutputs_nobackup/", dfAsRDSfiles[1]))
		paramName = str_split(splitNames[[1]][[2]], ".rds")[[1]][1]
		# Only keep the right column
		firstDF = firstDF[eachResult]
		colnames(firstDF) = c(paramName)
		firstDF["foldChange"] = rownames(firstDF)

		# browser()

		# Read each of the subsequence DF:
		if (length(dfAsRDSfiles) > 1){
			for (eachFileInd in 2:(length(dfAsRDSfiles))) {
				# Get the data and the set name
				thisDF = readRDS(paste0("./rdsOutputs_nobackup/", dfAsRDSfiles[eachFileInd]))
				thisSetName = str_split(splitNames[[eachFileInd]][[2]], ".rds")[[1]][1]

				# Keep the right column
				thisDF = thisDF[eachResult]
				colnames(thisDF) = c(thisSetName)
				thisDF["foldChange"] = rownames(thisDF)

				# Merge into the old one
				firstDF <- merge(firstDF, thisDF, by="foldChange")
			}
		}

		# Now plot:
		# browser()
		firstDF$foldChange = as.numeric(firstDF$foldChange)
		df.long <- melt(firstDF, id.vars="foldChange")
		
		png(paste0("./plots_nobackup/", eachResult, "_", groupTitle, ".png"))
		p <- ggplot(df.long, aes(foldChange, value, color=variable)) +
					geom_line() + ggtitle(paste0(eachResult, " by parameter")) +
					xlab("Fold Change") + ylab("Proportion") + ylim(0.0, 1.0)
		print(p)
		dev.off()
	}
}

plotScattersFromFileList <- function(rdsFileList, colForXaxis, colForYaxis){

}

plotPvalsFromFileList <- function(dfAsRDSfiles,groupTitle,
				pValLine = 2.5e-6, varToPlotBy = NULL){
	# read in the rds objects
	splitNames = str_split(dfAsRDSfiles, "_folds_")

	# Get the first DF:
	firstDF = readRDS(paste0("./rdsOutputs_nobackup/", dfAsRDSfiles[1]))
	paramName = str_split(splitNames[[1]][[2]], ".rds")[[1]][1]
	# Go through the pvals. Get the median pval
	pValMedians = rep(1, length(rownames(firstDF)))
	for (foldChangeInd in 1:(length(rownames(firstDF)))){
		pValMedians[foldChangeInd] = median(firstDF$pValsMM[[foldChangeInd]])
	}
	firstDF$pvalMedMM = pValMedians

	firstDF = firstDF["pvalMedMM"]
	colnames(firstDF) = c(paramName)
	firstDF["foldChange"] = rownames(firstDF)

	# Read each of the subsequence DF:
	if (length(dfAsRDSfiles) > 1){
		for (eachFileInd in 2:(length(dfAsRDSfiles))) {
			# Get the data and the set name
			thisDF = readRDS(paste0("./rdsOutputs_nobackup/", dfAsRDSfiles[eachFileInd]))
			thisSetName = str_split(splitNames[[eachFileInd]][[2]], ".rds")[[1]][1]
			# Go through the pvals. Get the median pval
			pValMedians = rep(1, length(rownames(thisDF)))
			for (foldChangeInd in 1:(length(rownames(thisDF)))){
				pValMedians[foldChangeInd] = median(thisDF$pValsMM[[foldChangeInd]])
			}
			thisDF$pvalMedMM = pValMedians

			# Keep the right column
			thisDF = thisDF["pvalMedMM"]
			colnames(thisDF) = c(thisSetName)
			thisDF["foldChange"] = rownames(thisDF)

			# Merge into the old one
			firstDF <- merge(firstDF, thisDF, by="foldChange")
		}
	}
	# Now plot:
	firstDF$foldChange = as.numeric(firstDF$foldChange)
	df.long <- melt(firstDF, id.vars="foldChange")

	print(df.long)
	
	png(paste0("./plots_nobackup/", "pValMedians_", groupTitle, ".png"))
	p <- ggplot(df.long, aes(foldChange, value, color=variable)) +
				geom_point() + ggtitle(paste0("Median MM Pval by param, line=", as.character(pValLine))) +
				xlab("Fold Change") + ylab("Median P Value") + scale_y_log10() +
				geom_hline(yintercept=pValLine) #ylim(0.0, 1.0) +
				#scale_x_log10()
	print(p)
	dev.off()
	# Make extra plots?
	if (!is.null(varToPlotBy)){
		if (varToPlotBy == "cells"){
		df.long$cells = as.numeric(sapply(str_split(df.long$variable, "_cells"), '[',1))
		df.long= df.long[c("foldChange", "cells", "value")]
		df.long$foldChange = as.character(df.long$foldChange)
		#Plot
		png(paste0("./plots_nobackup/", "pValMedians_Show_by_Cell_", groupTitle, ".png"))
		p <- ggplot(df.long, aes(cells, value, color=foldChange)) +
				geom_point() + ggtitle(paste0("Median MM Pval")) +
				xlab("Cell Count") + ylab("Median P Value") + scale_y_log10() 
				#geom_hline(yintercept=pValLine) #ylim(0.0, 1.0) +
				#scale_x_log10()
		print(p)
		dev.off()
		}
	}
}


powerCurveGivenPvalCutoff <- function(dfAsRDSfiles,groupTitle,
				pValCutoff = 2.5e-6){
	# read in the rds objects
	splitNames = str_split(dfAsRDSfiles, "_folds_")
	# Get the first DF:
	firstDF = readRDS(paste0("./rdsOutputs_nobackup/", dfAsRDSfiles[1]))
	paramName = str_split(splitNames[[1]][[2]], ".rds")[[1]][1]
	# Go through the pvals. Get the median pval
	pValProps = rep(1, length(rownames(firstDF)))
	for (foldChangeInd in 1:(length(rownames(firstDF)))){
		pValProps[foldChangeInd] = (sum(firstDF$pValsMM[[foldChangeInd]] < pValCutoff) /
							length(firstDF$pValsMM[[foldChangeInd]]))  #median(firstDF$pValsMM[[foldChangeInd]])
	}
	firstDF$power = pValProps
	firstDF = firstDF["power"]
	colnames(firstDF) = c(paramName)
	firstDF["foldChange"] = rownames(firstDF)

	# Read each of the subsequence DF:
	if (length(dfAsRDSfiles) > 1){
		for (eachFileInd in 2:(length(dfAsRDSfiles))) {
			# Get the data and the set name
			thisDF = readRDS(paste0("./rdsOutputs_nobackup/", dfAsRDSfiles[eachFileInd]))
			thisSetName = str_split(splitNames[[eachFileInd]][[2]], ".rds")[[1]][1]
			# Go through the pvals. Get the median pval
			pValProps = rep(1, length(rownames(thisDF)))
			for (foldChangeInd in 1:(length(rownames(thisDF)))){
				pValProps[foldChangeInd] = (sum(thisDF$pValsMM[[foldChangeInd]] < pValCutoff) /
									length(thisDF$pValsMM[[foldChangeInd]]))  #median(thisDF$pValsMM[[foldChangeInd]])
			}
			thisDF$power = pValProps

			thisDF = thisDF["power"]
			colnames(thisDF) = c(thisSetName)
			thisDF["foldChange"] = rownames(thisDF)

			# Merge into the old one
			firstDF <- merge(firstDF, thisDF, by="foldChange")
		}
	}
	# Now plot:
	firstDF$foldChange = as.numeric(firstDF$foldChange)
	df.long <- melt(firstDF, id.vars="foldChange")
	print(df.long)
	png(paste0("./plots_nobackup/", "Power_alpha=", as.character(pValCutoff),
			 "_", groupTitle, ".png"))
	p <- ggplot(df.long, aes(foldChange, value, color=variable)) +
				geom_point() + ggtitle(paste0("Power by alpha =", as.character(pValCutoff))) +
				xlab("Fold Change") + ylab("Power") + #scale_y_log10() +
				ylim(0.0, 1.0) 			
	print(p)
	dev.off()
}


powerCurveGivenPvalCutoffAndParamSet <- function(dfAsRDSfiles,groupTitle,
				pValCutoff = 2.5e-6, dataList = "pValsMM"){
	# read in the rds objects
	splitNames = str_split(dfAsRDSfiles, "_folds_")
	# Get the first DF:
	firstDF = readRDS(paste0("./rdsOutputs_nobackup/", dfAsRDSfiles[1]))
	paramName = str_split(splitNames[[1]][[2]], ".rds")[[1]][1]
	# Go through the pvals. Get the median pval
	pValProps = rep(1, length(rownames(firstDF)))
	for (foldChangeInd in 1:(length(rownames(firstDF)))){
		pValProps[foldChangeInd] = (sum(firstDF[[dataList]][[foldChangeInd]] < pValCutoff) /
							length(firstDF[[dataList]][[foldChangeInd]]))  #median(firstDF[[dataList]][[foldChangeInd]])
	}
	firstDF$power = pValProps
	firstDF = firstDF["power"]
	colnames(firstDF) = c(paramName)
	firstDF["foldChange"] = rownames(firstDF)

	# Read each of the subsequence DF:
	if (length(dfAsRDSfiles) > 1){
		for (eachFileInd in 2:(length(dfAsRDSfiles))) {
			# Get the data and the set name
			thisDF = readRDS(paste0("./rdsOutputs_nobackup/", dfAsRDSfiles[eachFileInd]))
			thisSetName = str_split(splitNames[[eachFileInd]][[2]], ".rds")[[1]][1]
			# Go through the pvals. Get the median pval
			pValProps = rep(1, length(rownames(thisDF)))
			for (foldChangeInd in 1:(length(rownames(thisDF)))){
				pValProps[foldChangeInd] = (sum(thisDF[[dataList]][[foldChangeInd]] < pValCutoff) /
									length(thisDF[[dataList]][[foldChangeInd]]))  #median(thisDF[[dataList]][[foldChangeInd]])
			}
			thisDF$power = pValProps

			thisDF = thisDF["power"]
			colnames(thisDF) = c(thisSetName)
			thisDF["foldChange"] = rownames(thisDF)

			# Merge into the old one
			firstDF <- merge(firstDF, thisDF, by="foldChange")
		}
	}
	# Now plot:
	firstDF$foldChange = as.numeric(firstDF$foldChange)
	df.long <- melt(firstDF, id.vars="foldChange")
	print(df.long)
	png(paste0("./plots_nobackup/", dataList, "_Power_alpha=", as.character(pValCutoff),
			 "_", groupTitle, ".png"))
	p <- ggplot(df.long, aes(foldChange, value, color=variable)) +
				geom_point() + ggtitle(paste0(dataList, " power by alpha =", as.character(pValCutoff))) +
				xlab("Fold Change") + ylab("Power") + #scale_y_log10() +
				ylim(0.0, 1.0) 			
	print(p)
	dev.off()
}


# rdsOutputFiles <- c("1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cellsPer_500_genes_3_indiv_per_group.rds",
# 						"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cellsPer_500_genes_5_indiv_per_group.rds",
# 						"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cellsPer_500_genes_8_indiv_per_group.rds")
# groupName = "FirstPowerRun"

# rdsOutputFiles <- c("1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cellsPer_500_genes_3_indiv_per_group.rds",
# 						"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cellsPer_500_genes_5_indiv_per_group.rds",
# 						"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_500_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 						"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_500_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds")
# groupName = "prelimRaiseNoise"

# rdsOutputFiles <- c("1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cellsPer_500_genes_3_indiv_per_group.rds",
# 						"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cellsPer_500_genes_5_indiv_per_group.rds",
# 						"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cellsPer_500_genes_8_indiv_per_group.rds",
# 						"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cellsPer_500_genes_10_indiv_per_group.rds"	)
# groupName = "prelimLowNoise"

# rdsOutputFiles <- c("1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.01sd_indiv_0.002sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.01sd_indiv_0.002sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.01sd_indiv_0.002sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.01sd_indiv_0.002sd_cell5_indiv.rds"
# 					# ,

# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.01sd_indiv_0.002sd_cell3_indiv.rds",
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.01sd_indiv_0.002sd_cell5_indiv.rds"
# 	)
# groupName = "miniVaryCellNum"

#############################################################################################################################
# rdsOutputFiles <- c(
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3090_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds"
# 	)
# groupName = "highNoiseVaryCellNum"

#############################################################################################################################
# rdsOutputFiles <- c(
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3090_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds"
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds"
# 	)
# groupName = "highNoiseThreeIndivVaryCellNum"


# rdsOutputFiles <- c(
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds" #,
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3090_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds"
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds"
# 	)
# groupName = "miniHighNoiseFiveIndivVaryCellNum"


# rdsOutputFiles <- c(
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",

# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds"#,
# 					# "1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds"
# 	)
# groupName = "miniHighNoiseVaryCellNum"

# rdsOutputFiles <- c(
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_500_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_500_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_500_genes_0.1sd_indiv_0.02sd_cell15_indiv.rds"
# 	)
# groupName = "miniHighNoise300cellVaryIndiv"

# rdsOutputFiles <- c(
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_500_genes_0.1sd_indiv_0.02sd_cell15_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_500_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds"
# 	)
# groupName = "15indivVs10000cell"

# rdsOutputFiles <- c(
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds"
# 	)
# groupName = "grid10indiv3000cell"

# rdsOutputFiles <- c(
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds"
# 	)
# groupName = "highNoise10indiv_VaryCell"

# rdsOutputFiles <- c(
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_300_cells_80_genes_0.01sd_indiv_0.002sd_cell10_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_1000_cells_80_genes_0.01sd_indiv_0.002sd_cell10_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_3000_cells_80_genes_0.01sd_indiv_0.002sd_cell10_indiv.rds",
# 					"1.1_1.25_1.5_2_3_4_5_7.5_10_folds_10000_cells_80_genes_0.01sd_indiv_0.002sd_cell10_indiv.rds"
# 	)
# groupName = "lowNoise10indiv_VaryCell"

# rdsOutputFiles <- c(
# 					"2_5_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"2_5_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					"2_5_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"2_5_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds"
# 	)
# groupName = "miniFoldFixedUMI"

# rdsOutputFiles <- c(
# 					"1_1.5_2_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					"1_1.5_2_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					"1_1.5_2_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					"1_1.5_2_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds"
# 	)
# groupName = "noUMIcovar_varyCells"

# rdsOutputFiles <- c(
# 					"1_1.5_2_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell3_indiv.rds",
# 					"1_1.5_2_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell5_indiv.rds",
# 					"1_1.5_2_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds",
# 					"1_1.5_2_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell15_indiv.rds"
# 	)
# groupName = "noUMIcovar_varyIndiv"

# rdsOutputFiles <- c(
# 					"1_1.5_2_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds",
# 					"1_1.5_2_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds",
# 					"1_1.5_2_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds"
# 	)
# groupName = "noUMIcovar10ind_varyCell"

rdsOutputFiles <- c(
			
			"1_1.5_2_5_7.5_10_folds_300_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds",
			"1_1.5_2_5_7.5_10_folds_1000_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds",
			"1_1.5_2_5_7.5_10_folds_3000_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds",
			"1_1.5_2_5_7.5_10_folds_10000_cells_80_genes_0.1sd_indiv_0.02sd_cell10_indiv.rds"
	)
groupName = "10indivVaryCells"

# plotCurveFromRDSdataframes(rdsOutputFiles, groupName)


# plotPvalsFromFileList(rdsOutputFiles, groupName, varToPlotBy="cells")

# powerCurveGivenPvalCutoff(rdsOutputFiles, groupName, .05)
# powerCurveGivenPvalCutoff(rdsOutputFiles, groupName, 2.5e-6)

powerCurveGivenPvalCutoffAndParamSet(rdsOutputFiles, groupName, .05, "pValsMM")
powerCurveGivenPvalCutoffAndParamSet(rdsOutputFiles, groupName, .05, "pValsFM")


powerCurveGivenPvalCutoffAndParamSet(rdsOutputFiles, groupName, 2.5e-6, "pValsMM")
powerCurveGivenPvalCutoffAndParamSet(rdsOutputFiles, groupName, 2.5e-6, "pValsFM")

print("All Done")




