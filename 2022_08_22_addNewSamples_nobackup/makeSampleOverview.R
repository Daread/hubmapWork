

# Read in the cds holding all cells
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")
library(monocle3)
# print(modStatus)
library("optparse")
library(dplyr)
library(ggplot2)

library(RColorBrewer)


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


set.seed(7)
# Get the passed parameters
option_list = list(
  make_option(c("-p", "--processingNote"), type="character", 
  			default="NucleiOnlySharedGenesCDS", 
              help="Processing note from upstream CDS", metavar="character"),
    
    make_option(c("-a", "--noAtrium"), type="logical", default = TRUE,
              help="Exclude samples from atrium", metavar="logical")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


# Get a CDS with all data combined
rdsPath = "../2022_08_22_addNewSamples_nobackup/formattedData/"
oldProcNote = opt$processingNote
allCellCDS = readRDS(paste0(rdsPath, oldProcNote, ".rds"))

processingNote = ""


# First, group the cell totals by sample
cellDF = as.data.frame(colData(allCellCDS))

if (opt$noAtrium){
	cellDF = cellDF[!(cellDF$Anatomical_Site %in% c("RA", "LA")),]
	processingNote = paste0(processingNote, "No_Atrium")
}

cellDF$Sample = paste0(cellDF$DataSource, "_", cellDF$Donor, "_", cellDF$Anatomical_Site)

# # Get the grouped samples
colToRetain = c("DataSource", "Sample", "Donor", "BMI", "Age", "Sex", "Diabetes", "Hypertension", "Anatomical_Site", "Instrument")
sampleDF = as.data.frame(cellDF %>% group_by(across(all_of(colToRetain))) %>% summarise(Cell_Total = n()))

# Output this sample dataframe
outDir = "./formattedData/"
write.csv(sampleDF, paste0(outDir, processingNote, "_Sample_Info.csv"))

# Now also collapse down samples by Donor and get donor-by-donor information
colToUse = c("Donor")
donorDF = as.data.frame(sampleDF %>% group_by(across(all_of(colToUse))) %>%
				 mutate(Sites_Sampled = paste0(Anatomical_Site, collapse=",")) %>%
				 mutate(Sample = NULL, Instrument = NULL, Anatomical_Site=NULL, Cell_Total=NULL) %>%
				 distinct()
				 )

write.csv(donorDF, paste0(outDir, processingNote, "_Donor_Info.csv"))


# Make a plot that summarizes the donor information by sex and age
plotDir = "./plots/"
dir.create(plotDir)

donorDF$DataSource = as.character(donorDF$DataSource)
donorDF$Age = as.numeric(donorDF$Age)

png(paste0(plotDir, "Meta_Analysis_Donor_Vis.png"), res=300, height = 1200, width = 1800)
myPlot = ggplot(as.data.frame(donorDF), aes(x=DataSource, y=Age, #color=Sex,
								 fill=Sex)) + 
			geom_boxplot(outlier.shape=NA) +
			geom_point(position=position_jitterdodge(), color="black") + 
			xlab("Publication") + 
				scale_color_brewer(palette= "Dark2")+ 
				scale_fill_brewer(palette= "Dark2")+ 
				theme(text=element_text(size=24)) +
				theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
				monocle_theme_opts()
print(myPlot)
dev.off()










