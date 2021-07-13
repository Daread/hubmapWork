
library("ggplot2")
library("dplyr")

# Read in the csv files

geneGroupHere = "Markers_Only"

arterialMarks = c("SEMA3G", "EFNB2", "DLL4")
venousMarks = c("NR2F2", "ACKR1")
endocardMarks = c("SMOC1", "NPR3")
lymphatMarks = c("PROX1", "TBX1", "PDPN")
markerList = list(arterialMarks, venousMarks, endocardMarks, 
					lymphatMarks, c(lymphatMarks, venousMarks))
names(markerList) = c("Arterial_+_Capillary", "Venous_Endothelium",
					"Endocardium", "Lymphatic_Endothelium",
					"Lymphatic_+_Venous")

# Read in the appropriate csv
for (eachInd in 1:length(markerList)){
	thisSubtype = names(markerList)[eachInd]
	print(paste0("Working on ", thisSubtype))
	thisCSV = paste0("./rdsOutput/Endoth_DE_Tets_", geneGroupHere, "_for_",
				thisSubtype, ".csv")
	thisDF = read.csv(thisCSV)
	thisTypeMarkers = markerList[[eachInd]]
	# Color code by if each gene is one of the markers for this subtype
	thisDF$known_marker = ifelse(thisDF$gene_short_name %in% thisTypeMarkers,
								"Marker", "Not_Marker")
	#Plot
	png(paste0("./plots/endothelialSubtypeDE/", thisSubtype, "_DEtestingFor_",
			geneGroupHere, ".png"))
	thisDF$negLog10_pval = -1 * log10(thisDF$p_value)
	thisPlot = ggplot(thisDF, 
		aes_string(x="estimate", y="negLog10_pval", color="known_marker")) +
		geom_point() + 
		xlab("Log2 Fold Change vs. Other Endothelial Cells") +
		ylab("-log10(p value)") + labs(title=paste0(thisSubtype, " Expression Testing"))
	print(thisPlot)
	dev.off()
}



# --------------------------------------------------------------------
# DE testing by site
# -------------------------------------------------------------------

geneGroupHere = "Full_Transcr"
thisSubtype = "Endocardium"
thisCSV = paste0("./rdsOutput/DE_outputs/Endoth_Anatomic_Tests_", geneGroupHere, "_for_",
				thisSubtype, ".csv")
thisDF = read.csv(thisCSV)

# # Sort by p values
# thisDF = thisDF[order(thisDF$p_value),]

# From 7-7-21, ran without filtering by frac of cells expressed
# # Get all the results together
# overallDF = data.frame()
# for (eachInd in 1:length(markerList)){
# 	thisSubtype = names(markerList)[eachInd]
# 	print(paste0("Working on ", thisSubtype))
# 	thisCSV = paste0("./rdsOutput/DE_outputs/Endoth_Anatomic_Tests_", geneGroupHere, "_for_",
# 				thisSubtype, ".csv")
# 	thisDF = read.csv(thisCSV)
# 	# thisTypeMarkers = markerList[[eachInd]]
# 	overallDF = rbind(overallDF, thisDF)

# }

# 7-8-21
# Get all the results together
overallDF = data.frame()
for (eachInd in 1:length(markerList)){
	thisSubtype = names(markerList)[eachInd]
	print(paste0("Working on ", thisSubtype))
	thisCSV = paste0("./rdsOutput/DE_outputs/Endoth_Anatomic_Tests_", geneGroupHere, "_for_",
				thisSubtype, "min0.01cellsExpr", ".csv")
	thisDF = read.csv(thisCSV)
	# thisTypeMarkers = markerList[[eachInd]]
	overallDF = rbind(overallDF, thisDF)

}


# Filter
qCutoff = .2
filteredDF = (overallDF %>% filter(q_value < qCutoff) %>% 
			filter(term == "Anatomical_SiteLeft_Vent"))
# Which are from left ventriccle
notLymphaticFiltered = (filteredDF %>% filter(cellSubtype != "Lymphatic_Endothelium"))
# P filtered DF
pCutoff = .05
pValFilteredDF = (overallDF %>% filter(p_value < pCutoff) %>% 
			filter(term == "Anatomical_SiteLeft_Vent"))
filteredDF[c("gene_short_name", "estimate", "q_value", "cellSubtype")]


# Make a histogram of inter-donor effect estimates

donorDF = (overallDF %>% filter(q_value < qCutoff) %>% 
			filter(grepl("Donor", term)))



noFilterAllAnatomical = (overallDF %>% 
			filter(term == "Anatomical_SiteLeft_Vent"))
noFilterDonor = (overallDF  %>% 
			filter(grepl("Donor", term)))

plotHistHelper = function(inputDF, plotName){
	png(paste0("./plots/", plotName, "histogram.png"))
	myPlot = ggplot(inputDF, aes_string("estimate")) + 
			labs(title=paste0("Effect size distribution for ", plotName)) +
			geom_histogram( aes(y=..count../sum(..count..)))
	print(myPlot)
	dev.off()
}

plotHistHelper(noFilterAllAnatomical, "Anatomical_Coefficients")
plotHistHelper(noFilterDonor, "Donor_Coefficients")



