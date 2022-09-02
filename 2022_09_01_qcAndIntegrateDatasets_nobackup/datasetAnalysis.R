
# library("devtools")
# library('spdep')
# load_all("/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/misc/monocle3mmVersion/monocle3")

library("monocle3")
library(ggplot2)

# Load helper functions
source("/net/trapnell/vol1/home/readdf/trapLabDir/sharedProjectCode/utility/singleCellUtilFuncs.R")

set.seed(7)
library("optparse")
# Get the passed parameters
option_list = list(
  make_option(c("-c", "--cdsRNA"), type="character", 
        default="../2022_08_22_addNewSamples_nobackup/formattedData/allDatasetCDS.rds", 
              help="Name of RNA cds to use", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Read in data
cdsPath = opt$cdsRNA
fullCDS = readRDS(cdsPath)

processingNote = "MNN_By_Dataset"


# Try the minimalist strategy of UMAPing all together, batch correct by source
fullCDS = estimate_size_factors(fullCDS)

fullCDS = preprocess_cds(fullCDS, num_dim=30)

fullCDS = align_cds(fullCDS, alignment_group="DataSource")

fullCDS = reduce_dimension(fullCDS, preprocess_method="Aligned")

# Plot outputs
dir.create('./plots/')

# Shuffle order for plotting
fullCDS = fullCDS[,sample(ncol(fullCDS))]

plotUMAP_Monocle(fullCDS, processingNote, "DataSource", 
	outputPath= "./plots/", show_labels=FALSE)


plotUMAP_Monocle(fullCDS, processingNote, "Cell_Shared_Label", 
	outputPath= "./plots/", show_labels=FALSE)


plotUMAP_Monocle(fullCDS, processingNote, "Anatomical_Site", 
	outputPath= "./plots/", show_labels=FALSE)


plotUMAP_Monocle(fullCDS, processingNote, "log10_umi", 
	outputPath= "./plots/", show_labels=FALSE)




dir.create("./rdsOutputs/")

# Save this cds if we want to come back later & not re-run the embedding steps
saveRDS(fullCDS, file=paste0("./rdsOutputs/", processingNote, "_CDS.rds"))



