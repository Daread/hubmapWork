
library(ggplot2)
library(ggrepel)

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


formatCellType = function(inputType){
	if (inputType == "Ventricular_Cardiomyocytes"){
		return("\nVentricular Cardiomyocytes")
	}
	if (inputType == "Neuron"){
		return("Neurons")
	}
	return(inputType)
}

plotCoefs = function(inputDF, genesToPlot, opt){

	# Subet to plotting DF
	plotDF = inputDF[inputDF$gene %in% genesToPlot,]
	# Add columns for upper/lower CI values
	plotDF$CI_Upper = plotDF$coefficientValue + 1.96 * plotDF$coef_std_error
	plotDF$CI_Lower = plotDF$coefficientValue - 1.96 * plotDF$coef_std_error

	# Add column for annotating q values
	# plotDF$Q_Annot = paste0("q=", as.character(format(round(plotDF$q_val, 2), nsmall = 2 )))
	plotDF$Q_Annot = paste0("q=", as.character(signif(plotDF$q_val, 3)))

	# Plot the desired coefficients with standard errors
	plotDir = "./plots/CoefficientPlots/"
	dir.create(plotDir)

	# Re-order genes to plot in order specified in input
	plotDF$gene = as.factor(plotDF$gene)
	plotDF$gene = factor(plotDF$gene, levels = rev(genesToPlot))

	# Get covariate label
	covarLabel = ifelse(opt$covar == "SexM", "Sex", opt$covar)
	thisCellType = formatCellType(opt$cellType)

	shift_in_set = c("SLC22A15", "KDM6A", "DOK5")
	plotDF$label_pos = ifelse(as.character(plotDF$gene) %in% shift_in_set, 
				plotDF$coefficientValue * .6, plotDF$coefficientValue)

	# Plot with ggplot
	png(paste0(plotDir, opt$cellType, "_Coef_", opt$genes, ".png"), res=300, height = 1500, width=1500)
	myPlot = ggplot(plotDF, aes_string( y = "gene", label="Q_Annot")) + 
				geom_point(aes_string(x="coefficientValue")) + 
				monocle_theme_opts() + 
				geom_text(aes_string(x="label_pos"), vjust=0, nudge_y=.1, size=6) +
				geom_errorbar(aes_string(xmin = "CI_Lower", xmax = "CI_Upper"), width = .1) +
				xlab(paste0("Coef in ", thisCellType, "\nby ", covarLabel)) +
				ylab("Gene")+ 
				theme(text=element_text(size=20)) + 
				xlim(min(plotDF$CI_Lower) * 1.2, max(plotDF$CI_Upper) * 1.2)

	print(myPlot)
	dev.off()

}

library("optparse")

# Get the passed parameters
option_list = list(
  make_option(c("-f", "--filePath"), type="character", 
  			# default="./plots/DE_Summaries/Combined_DE_by_",
  			default= "_fix_Anatomical_Site,Age,Sex,DataSource,log10_umi_rand_Donor_NucleiOnlySharedGenesCDS_noAtr_",
              help="Processing note from model fitting", metavar="character"),
  make_option(c("-q", "--qVal"), type="numeric", 
  			default=0.1,
              help="Q val cutoff for this ", metavar="character"),
  # make_option(c("-b", "--covariate"), type="character", 
  # 			default="SexM",   #"Age", # "SexM"
  #             help="Covariate to plot GSEA summary plot", metavar="character"),
  make_option(c("-c", "--cellType"), type="character", 
  			default="Neuron",
              help="Cell Type for which to plot coefficients", metavar="numeric"),
  make_option(c("-x", "--covar"), type="character", 
  			default="Age",
              help="Covariate to plot", metavar="numeric"),
  make_option(c("-g", "--genes"), type="character", 
  			default="LINGO1,KCNIP4,NXN",
              help="Covariate to plot", metavar="numeric")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

geneVec = strsplit(opt$genes, ",")[[1]]

# Read in the fit results
csvPath = paste0("./plots/", opt$cellType, opt$filePath, "/", opt$cellType, "_", opt$covar, "_AllTestsRun_Table.csv")
myData = read.csv(csvPath)

# Make the plot
plotCoefs(myData, geneVec, opt)




