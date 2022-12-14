
library(ggplot2)


plotCoefs = function(inputDF, genesToPlot, opt){

	# browser()
	# Subet to plotting DF
	plotDF = inputDF[inputDF$gene %in% genesToPlot,]
	# Add columns for upper/lower CI values
	plotDF$CI_Upper = plotDF$coefficientValue + 1.96 * plotDF$coef_std_error
	plotDF$CI_Lower = plotDF$coefficientValue - 1.96 * plotDF$coef_std_error

	# Plot the desired coefficients with standard errors
	plotDir = "./plots/CoefficientPlots/"
	dir.create(plotDir)

	# Plot with ggplot
	png(paste0(plotDir, opt$cellType, "_Coef_", opt$genes, ".png"), res=300, height = 1500, width=1500)
	myPlot = ggplot(plotDF, aes_string(x="coefficientValue", y = "gene")) + 
				geom_point() + 
				geom_errorbar(aes_string(xmin = "CI_Lower", xmax = "CI_Upper"))

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




