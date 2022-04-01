
library(ggplot2)
library(dplyr)
library("optparse")
library(stringr)
library(msigdbr)
library(fgsea)
library(tidyr)
library(tidyverse)

print("Libraries loaded, starting now")

# Get the passed parameters
option_list = list(
  make_option(c("-m", "--modelNotes"), type="character", 
  			# default="_fix_Anatomical_Site_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult", 
  			default="_fix_Anatomical_Site,Age,Sex_rand_Donor_HM10UMI=100_mito=10Scrub=0.2noPackerMNN=sampleNameK=40addAllTypes_MMresult",
              help="Processing note from model fitting", metavar="character"),
  make_option(c("-c", "--cellType"), type="character", 
  			# default="Endocardium", 
  			default="Vascular_Endothelium",
              help="Cell type for which the model was fit", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


modifyResultDFwithTest <- function(thisResultDF, pathwaysToTest, pathwayStatMod, 
                              pathwayList){
  # Loop and modify pathways
  for (eachPath in pathwaysToTest){
    genesToModify = pathwayList[[eachPath]]
    thisResultDF$z_value = ifelse(thisResultDF$gene %in% genesToModify,
                                thisResultDF$z_value + pathwayStatMod, thisResultDF$z_value) 
  }

  return(thisResultDF)
}

cellType = opt$cellType
modelFitDescription = opt$modelNotes

# Read in the fit
rdsPath = "./rdsOutput/mixedModels/"
outputDir = paste0("./plots/", cellType, modelFitDescription, "/")

# fixedCoefsToGet = c("SexM", "Age", "Anatomical_SiteRight_Vent", 
# 					"Anatomical_SiteApex", "Anatomical_SiteSeptum" )
fixedCoefsToGet = c("SexM")

msigdbr_df = msigdbr(species = "Homo sapiens", category = "H")
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

pathwaysToTest = c("HALLMARK_ALLOGRAFT_REJECTION")
pathwayStatMod = -10.0



#fgseaResTidy[fgseaResTidy$pathway == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",]

# Run GSEA on each coefficient's output
for (eachCoef in fixedCoefsToGet){
  set.seed(7)
  print(paste0("Working on GSEA for ", eachCoef))
	thisResultDF = read.csv(paste0(outputDir, cellType, "_", eachCoef, "_AllTestsRun_Table.csv"))

  thisResultDF = modifyResultDFwithTest(thisResultDF, pathwaysToTest, pathwayStatMod, msigdbr_list)

  # 
  namesAndStatistics = thisResultDF[c("gene", "z_value")] # If rank by test statistic
  # namesAndStatistics = thisResultDF[c("gene", "coefficient")] # If rankingby coefficient (which seems more natural biologically...?)
  colnames(namesAndStatistics) = c("gene", "statistic")

  # namesAndStatistics$statistic = abs(namesAndStatistics$statistic)
  namesAndStatistics = namesAndStatistics[order(-(namesAndStatistics$statistic)),]
  
  # Reformat into a list for input into fgsea
  namesAndStatisticsList = deframe(namesAndStatistics)


  # Get the gsea result for this coefficient
	# gseaResult = fgsea(pathways = msigdbr_list, stats=namesAndStatisticsList, nperm=1000)
 #  fgseaResTidy <- (gseaResult %>% as_tibble() %>%  arrange(desc(NES)))

 #  browser()

  multiLevel = fgseaMultilevel(pathways = msigdbr_list, stats=namesAndStatisticsList)
  fgseaResTidy = (multiLevel %>% as_tibble() %>%  arrange(desc(NES)))

  # Save the output

  saveDF = as.data.frame(fgseaResTidy)
  # Change the leadingEdge column to a character column, for now, to save it
  # sapply(output$leadingEdge, paste, collapse=", ") 
  saveDF$leadingEdge <- vapply(saveDF$leadingEdge, paste, collapse = ", ", character(1L))
  outputDF = paste0(outputDir, "GSEA_", cellType, "_", eachCoef, ".csv")
  # browser()
  # write.csv(saveDF, file=outputDF)

  print("GSEA Results:")
  str(saveDF)
  print(saveDF)

}













# For a QC check, show relationship between z_value, coefficients, and p values
# qcDF = thisResultDF[c("z_value", "coefficientValue", "pval")]
# png(paste0(outputDir, "QC_Plot", cellType, "_", eachCoef, ".png"), res=200, height=1000, width=1200)
# myPlot = ggplot(qcDF, aes_string(x="coefficientValue", y="z_value", col="pval")) +
#           geom_point()
# print(myPlot)
# dev.off()




# # Repeat, no limitation to hallmark category
# msigdbr_df = msigdbr(species = "Homo sapiens")
# msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# # Run GSEA on each coefficient's output
# for (eachCoef in fixedCoefsToGet){
#   print(paste0("Working on GSEA for ", eachCoef))
#   thisResultDF = read.csv(paste0(outputDir, cellType, "_", eachCoef, "_AllTestsRun_Table.csv"))

#   # 
#   namesAndStatistics = thisResultDF[c("gene", "z_value")] # If rank by test statistic
#   # namesAndStatistics = thisResultDF[c("gene", "coefficient")] # If rankingby coefficient (which seems more natural biologically...?)
#   colnames(namesAndStatistics) = c("gene", "statistic")

#   # namesAndStatistics$statistic = abs(namesAndStatistics$statistic)
#   namesAndStatistics = namesAndStatistics[order(-(namesAndStatistics$statistic)),]
#   # Reformat into a list for input into fgsea
#   namesAndStatisticsList = deframe(namesAndStatistics)


#   # Get the gsea result for this coefficient
#   gseaResult = fgsea(pathways = msigdbr_list, stats=namesAndStatisticsList, nperm=1000)
#   fgseaResTidy <- (gseaResult %>% as_tibble() %>%  arrange(desc(NES)))

#   print("Plotting now")

#   # Plot the result
#   png(paste0(outputDir, "ALL_PATHWAY_GSEA_", cellType, "_", eachCoef, ".png"), res=200, height=2400, width=2600)
#   myPlot = ggplot(fgseaResTidy[fgseaResTidy$padj<0.1,], aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill=padj<0.05)) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title=paste0("All pathways GSEA ", cellType)) + 
#   theme_minimal()
#   print(myPlot)
#   dev.off()

#   # browser()

# }



