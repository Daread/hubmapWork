

###########################################################################################
# Load relevant packages
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
# modStatus <- loadMonoclePackages()
# print(modStatus)
library(monocle3)
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/scRNA_Seq_Background_Code/backgroundMethods.R")
source("./propUtilCode.R")
library(reshape2)
library("dplyr")
library(data.table)
library(tidyr)
library(VGAM)
library(ggplot2)
library(RColorBrewer)

library(stringr)

# Get the passed parameters
library("optparse")

option_list = list(
  make_option(c("-n", "--name"), type="character", default="TestRun", 
              help="Name of this subset analysis", metavar="character"),
    make_option(c("-c", "--colToUse"), type="character", default="highLevelCellType", 
              help="Column for which to run proportion tests", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

processingNote = paste0(opt$name, "_by_", opt$colToUse)

# With arguments, read in data
rdsPath = "../2021_05_10_HM10_Preprocessing_and_cell_annotations/rdsOutput/"
# oldProcNote = "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40"
oldProcNote = "HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes"
allCellCDS = readRDS(paste0(rdsPath, "allCells_", oldProcNote, ".rds"))

# Assign the subtypes
colData(allCellCDS)$sample = colData(allCellCDS)$sampleName
allCellCDS = hardAssignHighLevelCellTypes(allCellCDS, oldProcNote)
allCellCDS = hardAssignDonors(allCellCDS)
allCellCDS = hardAssignAnatomicalSites(allCellCDS)
allCellCDS = hardAssignDonorSexes(allCellCDS)
allCellCDS = hardAssignDonorAges(allCellCDS)

# Based off of this original CDS name, make a new directory for outputs
oldCDS_Path = paste0("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_14_cellTypeProportions_nobackup/plots/", oldProcNote, "/")
dir.create(file.path(oldCDS_Path), showWarnings=FALSE)
outputPath = paste0(oldCDS_Path, processingNote, "/")
dir.create(file.path(outputPath))


set.seed(7)
# Loop over each unique value in the colToUse
for (eachGroup in as.character(levels(as.factor(colData(allCellCDS)[[opt$colToUse]])))){
	# Plot proportions of this
	makeFacetedProportionPlot_withFill(allCellCDS, paste0(processingNote, "_", eachGroup), "sample", opt$colToUse, 
            fillCol = "Donor", subsetPropsToShow=c(eachGroup), colorCol="Anatomical_Site",
            outputPath=outputPath )
	
}

# Next, look and see if there are obvious trends in cell composition
# cellPropDF = (as.data.frame(colData(allCellCDS)) %>% )
countTable = with(as.data.frame(colData(allCellCDS)), table(get("sample"), get(opt$colToUse)))
propTable = prop.table(countTable, margin=1)
propDF = as.data.frame(propTable)



monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}


# Beta Binomial fitting


 getFormattedCountDF <- function(inputDF){
 	# Label the columns
 	colnames(inputDF) = c("Sample", "Cell_Type", "Count")
 	# Add info on anatomical site
 	inputDF = tidyr::separate(inputDF, col="Sample", into=c("Donor", "Anatomical_Site"), 
 					remove=FALSE, extra="merge", sep="[.]")
 	inputDF = hardAssignDonorSexesDF(inputDF)
 	inputDF = hardAssignDonorAgesDF(inputDF)

 	# I also need to add a column for the total cells per sample
 	sampleNames = as.character(unique(inputDF$Sample))
 	cellSampleTotals = rep(0, length(sampleNames))
 	names(cellSampleTotals) = sampleNames
 	# Loop and find the ttoal per sample
 	for (eachSample in sampleNames){
 		subsetTable = inputDF[inputDF$Sample == eachSample,]
 		cellSampleTotals[eachSample] = sum(subsetTable$Count)
 	}
 	# Add this to the df
 	inputDF$TotalCellsPerSample = cellSampleTotals[inputDF$Sample]

 	return(inputDF)
 }




# Format the dataframe to contain counts alongside information on anatomical site, age, and sex
countDF = melt(countTable)
cellTypeRegressionDF = getFormattedCountDF(countDF)




formatCellType <- function(inputColumn){

  inputColumn = ifelse(inputColumn == "T_Cell", "T Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "VSM_and_Pericyte", "Perivascular Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "Vascular_Endothelium", "Vascular Endothelium", inputColumn)
  inputColumn = ifelse(inputColumn == "Lymphatic_Endothelium", "Lymphatic Endothelium", inputColumn)
  inputColumn = ifelse(inputColumn == "Mast_Cell", "Mast Cell", inputColumn)
  inputColumn = ifelse(inputColumn == "B_Cell", "B Cell", inputColumn)

  return(inputColumn)
}


getBetaBinomialPropFits <- function(inputDF, modelFormula="countDF ~ Anatomical_Site + Sex + Age"){
	# Loop and fit a model per cell type
	cellTypesToTest = as.character(unique(inputDF$Cell_Type))
	inputDF$Cell_Type = as.character(inputDF$Cell_Type)
	fitResList = vector(mode="list", length(cellTypesToTest))
	origFitList = vector(mode="list", length(cellTypesToTest))
	counter = 1
	for (eachCelltype in cellTypesToTest){
		print(paste0("Fitting ", eachCelltype))
		subsetDF = inputDF[inputDF$Cell_Type == eachCelltype,]
		# Get the fit
		countDF = cbind(subsetDF$Count, subsetDF$TotalCellsPerSample - subsetDF$Count)
		# browser()
		fit =  vglm(as.formula(modelFormula), betabinomial, data = subsetDF, trace = TRUE)
		fit_df = tidy.vglm(fit)
		fit_df$Cell_Type = eachCelltype
		fitResList[[counter]] = fit_df
		origFitList[[counter]] = fit
		counter = counter + 1
	}

	names(fitResList) = cellTypesToTest
	names(origFitList) = cellTypesToTest

	# return(fitResList)
	return(list(fitResList, origFitList))
}


# Run the regression model on each

# betaBinomFits = getBetaBinomialPropFits(cellTypeRegressionDF)
fullAndTidyFits = getBetaBinomialPropFits(cellTypeRegressionDF)
betaBinomFits = fullAndTidyFits[[1]]
fullFits = fullAndTidyFits[[2]]


test_res = do.call(rbind, betaBinomFits)

nonCoefRes = test_res[!(test_res$term %in% c("(Intercept):1", "(Intercept):2")),]
nonCoefRes = nonCoefRes[order(nonCoefRes$p.value),]

nonCoefRes$q_value = p.adjust(nonCoefRes$p.value, method = "BH")

# Save the results
# nonCoefRes[1] = NULL
rownames(nonCoefRes) = NULL
testResultOutfile = paste0("./fileOutputs/betaBinomFitting.csv")
write.csv(nonCoefRes, testResultOutfile)

# countDF = as.data.frame(countTable)

# formattedCountDF = getFormattedCountDF(countDF)
















# Make simple jitter plots
jitterDir = paste0("./plots/boxplots/")
dir.create(jitterDir)

# Format the cell names for plotting
cellTypeRegressionDF$Cell_Type = as.character(cellTypeRegressionDF$Cell_Type)
cellTypeRegressionDF$Cell_Type = formatCellType(cellTypeRegressionDF$Cell_Type)

cellTypesToTest = as.character(unique(cellTypeRegressionDF$Cell_Type))

cellTypeRegressionDF$Proportion = (cellTypeRegressionDF$Count / cellTypeRegressionDF$TotalCellsPerSample)

# For every cell type, plot by age/sex/site
for (categoricalVar in c("Sex", "Anatomical_Site")){
	# for (eachCelltype in cellTypesToTest){
		# Make a box plot
		png(paste0(jitterDir, "Proportions_grouped_by_", categoricalVar, ".png"), 
				res=300, height=1800,width=2700)
		myPlot = ggplot(cellTypeRegressionDF, aes_string(x="Cell_Type", y="Proportion",
					# col=categoricalVar, 
					fill=categoricalVar)) + 
				geom_boxplot(outlier.shape=NA) + 
				monocle_theme_opts() + 
				geom_point(position=position_jitterdodge()) +
				theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
				theme(text=element_text(size=24)) + 
				xlab("Cell Type") + 
				scale_fill_brewer(palette= "Dark2")
				# geom_jitter()
		print(myPlot)
		dev.off()
	# }
}

# Scatter plot for each cell type showing age vs. 
for (eachCelltype in cellTypesToTest){

	subsetDF = cellTypeRegressionDF[cellTypeRegressionDF$Cell_Type == eachCelltype,]

	png(paste0(jitterDir, "Age_vs_prop_for_", eachCelltype, ".png"), 
			res=200, height=1000,width=1200)
	myPlot = ggplot(subsetDF, aes_string(x="Age", y="Proportion")) + 
			# geom_point()  +
			 geom_point(aes(col=Sex)) +
			 monocle_theme_opts() + 
			theme(text=element_text(size=20)) + 
			ylab(paste0(eachCelltype, " Proportion")) +
			geom_smooth(method = "lm", se = FALSE, col="black") + 
			 scale_color_brewer(palette="Dark2")
		  # stat_summary(fun.data= mean_cl_normal) + 
		  #geom_smooth(method='lm')
			# theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
			# geom_jitter()
	print(myPlot)
	dev.off()

	print(eachCelltype)
	print(cor(subsetDF$Age, subsetDF$Proportion))

}




getAdjustedProportions = function(inputPropDF, inputFits,
									regressOut){
	# Format the cell type names from beta binomial fit object
	inputFits$Cell_Type = formatCellType(inputFits$Cell_Type)

	cellTypes = as.character(levels(as.factor(inputPropDF$Cell_Type)))
	coefficientsFit = as.character(levels(as.factor(inputFits$term)))
	uncorrectedProps = inputPropDF$Proportion
	logitUncorrected = log(uncorrectedProps / (1 - uncorrectedProps))
	correctedProps = uncorrectedProps
	logitCorrected = logitUncorrected
	# Get the variables that appeared 
	for (eachVar in regressOut){
		coefMatchingVarInd = str_detect(coefficientsFit, eachVar)
		coefMatchingVar = coefficientsFit[coefMatchingVarInd]
		#Dummy variable this represents?
		dummyVarValues = str_replace(coefMatchingVar, eachVar, "")
		# browser()
		# Loop through the variables to regress out
		for (eachDummy in dummyVarValues){
			# browser()
			dummyVector = ifelse(inputPropDF[[eachVar]] == paste0(eachDummy), 1, 0)
			# Get the fit values for each cell type
			miniFit = inputFits[inputFits$term == paste0(eachVar, eachDummy),]
			coefsAsVec = miniFit$estimate
			names(coefsAsVec) = miniFit$Cell_Type
			cellTypeSpecCoefs = coefsAsVec[inputPropDF$Cell_Type]
			# Mulitply these together
			thisAdjustment = dummyVector * cellTypeSpecCoefs

			# browser()
			# Adjust
			logitCorrected = logitCorrected - thisAdjustment
		}
	}
	# Now that we have the corrected logit, get corrected proportion
	correctedProps = exp(logitCorrected) / (1 + exp(logitCorrected))

	return(correctedProps)
}






# newDFmale = data.frame("Sex" = c("M", "M", "M"), "Anatomical_Site" = c("Apex", "Apex", "Apex"), "Age" = c(25, 45, 60))
# newDFfemale = data.frame("Sex" = c("F", "F", "F"), "Anatomical_Site" = c("Apex", "Apex", "Apex"), "Age" = c(25, 45, 60))

newDFfemale = data.frame("Sex" = rep("F", 36), "Anatomical_Site" = rep("Apex", 36), "Age" = 25:60)
newDFmale = data.frame("Sex" = rep("M", 36), "Anatomical_Site" = rep("Apex", 36), "Age" = 25:60)


darkPalette = brewer.pal(8, "Dark2")


# Make age vs. adjusted proportion plots, regressing out sex + anatomical site

# cellTypeRegressionDF$adjustedPropSexSite = getAdjustedProportions(cellTypeRegressionDF, test_res, 
# 																regressOut=c("Sex", "Anatomical_Site"))

cellTypeRegressionDF$adjustedPropSexSite = getAdjustedProportions(cellTypeRegressionDF, test_res, 
																regressOut=c( "Anatomical_Site"))

# Format test_res to have cell type names match
test_res$Cell_Type = formatCellType(test_res$Cell_Type)
cellTypeRegressionDF$Cell_Type = formatCellType(cellTypeRegressionDF$Cell_Type)

# Scatter plot for adjusted props in each cell type vs. age
for (eachCelltype in cellTypesToTest){

	subsetDF = cellTypeRegressionDF[cellTypeRegressionDF$Cell_Type == eachCelltype,]

	malePred = predict(fullFits[[eachCelltype]], newdata = newDFmale, se=T)
	# browser()
	malePredDF = data.frame("logitlink_mu"=as.data.frame(malePred$fitted.values)[["logitlink(mu)"]], 
							"se" = as.data.frame(malePred$se.fit)[["logitlink(mu)"]],
							"adjustedPropSexSite" = exp(as.data.frame(malePred$fitted.values)[["logitlink(mu)"]]) / 
											(1 + exp(as.data.frame(malePred$fitted.values)[["logitlink(mu)"]]))	,
							   "Age" = newDFmale$Age)
	maleColor = darkPalette[2]

	femalePred = predict(fullFits[[eachCelltype]], newdata = newDFfemale, se=T)
	# browser()
	femalePredDF = data.frame("logitlink_mu"=as.data.frame(femalePred$fitted.values)[["logitlink(mu)"]], 
							"se" = as.data.frame(femalePred$se.fit)[["logitlink(mu)"]],
							"adjustedPropSexSite" = exp(as.data.frame(femalePred$fitted.values)[["logitlink(mu)"]]) / 
											(1 + exp(as.data.frame(femalePred$fitted.values)[["logitlink(mu)"]]))	,
							   "Age" = newDFfemale$Age)
	femaleColor = darkPalette[1]

	colnames(malePred$fitted.values)[1]

	# Get the predictions

	# Get the intercept and age coefficient for this cell type
	cellTypeFits = test_res[test_res$Cell_Type == eachCelltype,]
	cellTypeFitVec = cellTypeFits$estimate
	cellTypeSEvec = cellTypeFits$std.error
	names(cellTypeFitVec) = cellTypeFits$term
	names(cellTypeSEvec) =  cellTypeFits$term
	thisInt = cellTypeFitVec["(Intercept):1"]
	ageCoef = cellTypeFitVec["Age"]
	ageSE   = cellTypeSEvec["Age"]

	# browser()

	png(paste0(jitterDir, "Age_vs_adjustedProp_for_", eachCelltype, ".png"), 
			res=200, height=1000,width=1200)
	myPlot = ggplot(data = subsetDF, aes_string(x="Age", y="adjustedPropSexSite")) + 
			# geom_point()  +
			 geom_point(aes(col=Sex)) +
			 monocle_theme_opts() + 
			theme(text=element_text(size=20)) + 
			ylab(paste0(eachCelltype, " Adjusted Proportion")) +
			# geom_smooth(method = "lm", se = FALSE, col="black")
			# geom_function(fun = function(x) exp(thisInt + ageCoef*x)/(1 + exp(thisInt + ageCoef*x))) + 
			# geom_ribbon(aes(ymin =(  exp(thisInt + (ageCoef - 2*ageSE)*Age)/(1 + exp(thisInt + (ageCoef - 2*ageSE)*Age)) ) ,
			# 				ymax =(  exp(thisInt + (ageCoef + 2*ageSE)*Age)/(1 + exp(thisInt + (ageCoef + 2*ageSE)*Age)) ) ), fill="grey70", 
			# 				alpha=.15 ) + 
			# Male SE and predictions
			geom_line(data=malePredDF, color=maleColor) + 
			geom_ribbon(data=malePredDF,
				aes(ymin =(          exp((logitlink_mu - 2*se))/(1 + exp( (logitlink_mu - 2*se))) ) ,
							ymax =(  exp((logitlink_mu + 2*se))/(1 + exp( (logitlink_mu + 2*se))) ) ), fill=maleColor, 
							alpha=.15 ) +
			# ...and for female
			geom_line(data=femalePredDF, color=femaleColor) + 
			geom_ribbon(data=femalePredDF,
				aes(ymin =(          exp((logitlink_mu - 2*se))/(1 + exp( (logitlink_mu - 2*se))) ) ,
							ymax =(  exp((logitlink_mu + 2*se))/(1 + exp( (logitlink_mu + 2*se))) ) ), fill=femaleColor, 
							alpha=.15 ) +

			 scale_color_brewer(palette="Dark2")
	print(myPlot)
	dev.off()

	print(eachCelltype)

}



# Version getting CI based on predictions
































































################ End beta binomial fitting

##############################################################################################################
# Formatting for generating panel correlation plots
unmeltedPropDF = dcast(propDF, Var2 ~ Var1)
colnames(unmeltedPropDF)[1] = opt$colToUse
rownames(unmeltedPropDF) = unmeltedPropDF[[opt$colToUse]]
unmeltedPropDF = unmeltedPropDF[, -which(names(unmeltedPropDF) %in% c(opt$colToUse))]


# From https://r-coder.com/correlation-plot-r/

# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    Cor <- (cor(x, y)) # Remove abs function if desired
    txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
    if(missing(cex.cor)) {
        cex.cor <- 0.4 / strwidth(txt)
    }
    text(0.5, 0.5, txt,
         cex = 1 + cex.cor * abs(Cor)) # Resize the text by level of correlation
}
######


png(paste0(outputPath, "scatterPlotMatrix_", opt$colToUse, ".png"),
				width=2000, height=2000, res=200)
thisPlot = pairs(t(unmeltedPropDF),
				upper.panel = panel.cor)
print(thisPlot)
dev.off()


# Make a plot just showing fibroblast vs cardiomyocyte abundance
fibCardioDF = data.frame("Fibroblasts" = as.numeric(unmeltedPropDF["Fibroblast",]), 
						"Cardiomyocytes" = as.numeric(unmeltedPropDF["Cardiomyocyte",]))

png(paste0(outputPath, "_FibroblastVsCardiomyocytes",".png"),
				width=1400, height=1400, res=200)
thisPlot = ggplot(fibCardioDF, aes(x=Fibroblasts, y=Cardiomyocytes)) + 
			geom_point() + ggtitle("Proportion of Fibroblasts Vs Cardiomyocytes") + 
			theme(text=element_text(size=18)) 
print(thisPlot)
dev.off()





# Just make a plot to show the UMAP
plotUMAP_Monocle(allCellCDS, processingNote, "highLevelCellType", textSize=4,
				show_labels=TRUE, outputPath=outputPath)




# 8-17-21: Get a dataframe with proportions of each cell type for each sample, labeling donor and site.
#.         Then use this to look at inter-site vs inter-donor variation, as well as 

cellTypePropDF = transpose(unmeltedPropDF)
rownames(cellTypePropDF) = colnames(unmeltedPropDF)
colnames(cellTypePropDF) = rownames(unmeltedPropDF)

cellTypePropDF$sampleName = rownames(cellTypePropDF)

meltedPropDF = reshape::melt(cellTypePropDF, id="sampleName")
colnames(meltedPropDF) = c("sampleName", "Cell_Type", "Proportion")



meltedPropDF = tidyr::separate(meltedPropDF, col="sampleName", into=c("Donor", "Anatomical_Site"), 
				sep="[.]", remove=FALSE, extra="merge")






# Now, time to get 
varianceComparisonDF = data.frame("Mean_Value" = double(),
								"Comparison" = character(),
								  "Cell_Type" = character())

cellTypes = c("Adipocytes", "B_Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic_Endothelium", "Macrophage", "Mast_Cell", "Neuronal", "T_Cell",
                "Vascular_Endothelium", "VSM_and_Pericyte")



getInterDonorSetSite <- function(subsetDF, siteToUse){
	thisDF = subsetDF[subsetDF$Anatomical_Site == siteToUse,]
	# Loop over all combinations
	comboMeanVec = numeric()

	for (eachInd in 1:(nrow(thisDF) - 1)){
		for (eachComp in (eachInd + 1):nrow(thisDF)){
			thisAbsDiff = abs(thisDF[eachInd, "Proportion"] - thisDF[eachComp, "Proportion"] )

			comboMeanVec = c(comboMeanVec, thisAbsDiff)
		}
	}
	return ( mean(comboMeanVec))
}


getIntraDonorDiff <- function(subsetDF, siteOne, siteTwo){
	thisDF = subsetDF[subsetDF$Anatomical_Site %in% c(siteOne, siteTwo),]

	interSiteVec = numeric()
	# Get all the intra donor combos
	for (eachDonor in levels(as.factor(subsetDF$Donor))){
		miniDF = thisDF[thisDF$Donor == eachDonor,]
		# If only one, skip
		if (nrow(miniDF) == 1){
			next 
		}
		# Otherwise, get the comparison
		interSiteDiff = abs(miniDF[1, "Proportion"] - miniDF[2, "Proportion"])
		interSiteVec = c(interSiteVec, interSiteDiff)
	}

	return(mean(interSiteVec))
}



for (eachCellType in cellTypes){
	subsetDF = meltedPropDF[meltedPropDF$Cell_Type == eachCellType,]
	interLV = getInterDonorSetSite(subsetDF, "Left.Vent")
	interApex = getInterDonorSetSite(subsetDF, "Apex")

	interSiteIntraDonor = getIntraDonorDiff(subsetDF, "Left.Vent", "Apex")

	cellTypeDF = data.frame("Mean_Value" = c(interLV, interApex, interSiteIntraDonor))
	cellTypeDF$Comparison = c("Inter_Donor_LV", "Inter_Donor_Apex", "Inter_Site_Within_Donor")
	cellTypeDF$Cell_Type = eachCellType
	varianceComparisonDF = rbind(varianceComparisonDF, cellTypeDF)

}

# Make a plot


png(paste0(outputPath, "Compare_InterSite_vs_InterDonor_Differences", ".png"),
				width=2000, height=2000, res=200)
thisPlot = ggplot(varianceComparisonDF, aes_string(x="Cell_Type", y="Mean_Value", color="Comparison")) +
				geom_point()+ 
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(thisPlot)
dev.off()



# Calculate the t test for LV vs apex proportions
tTestDF = data.frame("Cell_Type" = cellTypes)
rownames(tTestDF) = tTestDF$Cell_Type

for (eachCellType in cellTypes){
	subsetDF = meltedPropDF[meltedPropDF$Cell_Type == eachCellType,]
	subsetDF = subsetDF[subsetDF$Anatomical_Site %in% c("Left.Vent", "Apex"),]

	# Get the t test
	xVec = subsetDF[subsetDF$Anatomical_Site == "Left.Vent",]$Proportion
	yVec = subsetDF[subsetDF$Anatomical_Site == "Apex",]$Proportion
	tRes = t.test(xVec, yVec)

	tTestDF[eachCellType, "P_Value_LV_vs_Apex"] = tRes$p.value

}

# Save this

write.csv(tTestDF, paste0(outputPath, "tTestsForProportions.csv"))

























