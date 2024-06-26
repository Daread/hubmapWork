

library("monocle3")
library(ggplot2)
library(stringr)

# Load helper functions
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")


monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    # theme(panel.grid.minor.x = element_blank(),
    #       panel.grid.minor.y = element_blank()) +
    # theme(panel.grid.major.x = element_blank(),
    #       panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank()) + 
    theme(
	  panel.grid.major = element_line(size = 0.25, linetype = 'solid',
	                                colour = "grey"), 
	  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                colour = "grey")
  )
}


getPropDF = function(opt){

	# Read in the file
	myDF = read.csv(opt$propData)

	return(myDF)
}


adjustProportions = function(thisPropDF, fitInput){
	# Need to remove the effect of datasource, sex, and anatomical site

	# Set up adjusted column
	# thisPropDF$Adjusted_Prop = thisPropDF$Proportion

	# Data source:
	# sourcesToRegress = c("Read", "Litvinukova", "Tucker", "Koenig", "Reichart") # Chaffin, alphabetically first, was used first arbitrarily
	# for (eachSource in sourcesToRegress){
	# 	# 
	# }

	logitUncorrected =  log( thisPropDF$Proportion / (1 - thisPropDF$Proportion))

	varToRegress = c("DataSource", "Anatomical_Site", "Sex")
	for (eachVar in varToRegress){
		# Loop for each of its levels
		varCoefs = fitInput$term[str_detect(fitInput$Term, eachVar)]
		# Get the levels
		dummyVars = unique(varCoefs)
		origLevels = gsub(eachVar, "", dummyVars)
		# Now loop through the proportion dataframe for each level
		for (eachLevel in origLevels){
			# Get index of which entries match this
			dummyMatch = ifelse(as.character(thisPropDF[[eachVar]]) == eachLevel, 1, 0)
			# Get the coefficient
			thisCoef = fitInput[fitInput$Term == paste0(eachVar, eachLevel),]$Estimate[1]

			# Adjust
			thisAdjustment = dummyMatch * thisCoef
			logitUncorrected = logitUncorrected - thisAdjustment
		}
	}

	# Back into proportions, not logit terms
	thisPropDF$Proportion = exp(logitUncorrected) / (1 + exp(logitUncorrected))

	return(thisPropDF)

}


plotProportionAgeFit = function(eachType, fitInput, propInput, opt, adjustProps = TRUE,
								pointColor = "DataSource"){



	# browser()
	# Get the data from this celltype
	thisPropDF = propInput[propInput$Cell_Type == eachType,]
	thisPropDF$Proportion = thisPropDF$Count * 1.0 / thisPropDF$TotalCellsPerSample
	thisFitDF = fitInput[fitInput$Cell_Type == eachType,]

	# If adjusting proportions, adjust here
	if (adjustProps){
		thisPropDF = adjustProportions(thisPropDF, thisFitDF)
		adjNote = "Adjusted_Props"
		y_lab_note = paste0(eachType, " Normalized Proportion")
	} else{
		adjNote = "Raw_Props"
		y_lab_note = paste0(eachType, " Proportion")
	}

	# Get the points for the curve 
	curveDF = data.frame("Age" = thisPropDF$Age )
	interceptVal = thisFitDF[thisFitDF$Term == "(Intercept)",]$Estimate[1]
	ageCoef = thisFitDF[thisFitDF$Term == "Age",]$Estimate[1]
	logitProp = interceptVal + (ageCoef * curveDF$Age)  #log(curveDF$Age / (1.0 - curveDF$Age))
	propCurve = exp(logitProp) / (1 + exp(logitProp))
	curveDF$LogitProp = logitProp
	curveDF$Proportion = propCurve

	# browser()

	# Get data needed for ribbons for +/- 2 SE
	curveDF$se = thisFitDF[thisFitDF$Term == "Age",]$StdErr[1]


	# Plot the result
	outDir = paste0("./plots/Proportions/")
	dir.create(outDir)

	# browser()

	if (pointColor == "DataSource"){
		pointColor = "Publication"
		thisPropDF$Publication = thisPropDF$DataSource
	}

	# Make the plot
	png(paste0(outDir, eachType, "_Age_Vs_", adjNote, "ColorBy_", pointColor, ".png"), res=200, height=1400, width=1200)
	myPlot = ggplot() + 
				geom_point(data=thisPropDF, aes_string(x="Age", y="Proportion", color=pointColor)) + 
				geom_line(data=curveDF, aes_string(x="Age", y="Proportion")) + 
				geom_ribbon(data=curveDF,
					aes(x=Age,
						 ymin = (exp((logitProp - 2*se*Age))/(1 + exp( (logitProp - 2*se*Age))) ),
						ymax = (  exp((logitProp + 2*se*Age))/(1 + exp( (logitProp + 2*se*Age))) ) ),
						 fill = "blue", alpha=.15) + 
				ylab(y_lab_note)+ 
				theme(text=element_text(size=24))  + 
				guides(colour = guide_legend(override.aes = list(size=10))) +
				monocle_theme_opts()

				# (data=malePredDF,
				# aes(ymin =(          exp((logitProp - 2*se))/(1 + exp( (logitProp - 2*se))) ) ,
				# 			ymax =(  exp((logitProp + 2*se))/(1 + exp( (logitProp + 2*se))) ) ), fill=, 
				# 			alpha=.15 )

	print(myPlot)
	dev.off()
}


set.seed(7)
library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-r", "--propResults"), type="character", 
        # default="../2022_08_22_addNewSamples_nobackup/formattedData/allDatasetCDS.rds", 
        default="./fileOutputs/No_Atrium_betaBinomFittingMixed.csv", 
              help="Path to BBmm outputs", metavar="character"),
  make_option(c("-p", "--propData"), type="character", 
        # default="../2022_08_22_addNewSamples_nobackup/formattedData/allDatasetCDS.rds", 
        default="./fileOutputs/CellProportions_No_Atrium_.csv", 
              help="Path to file with raw proportion data", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


processingNote = "MetaProp"

# Get the results of the fit models
fitRes = read.csv(opt$propResults)

# Get the dataframe that holds all raw proportions
allProps = getPropDF(opt)

# Plot proportions
ageCellTypes = c("Neuron")

for (eachType in ageCellTypes){
	plotProportionAgeFit(eachType, fitRes, allProps, opt)
}


















