

centos6LoadMonoclePackages <- function(){

    modStatus <- tryCatch({
        suppressPackageStartupMessages({
            library(dplyr)
            library(SoupX)
            library(devtools)
            # library("uwot")
            library("irlba")
            library("Rtsne")
            library("slam")
            library("sp")
            library("sf")
            library("spdep")

            # Get custom uwot loaded
            load_all("~/bin/uwot")
            # Use only the libpath leading to the core monocle installation
            .libPaths('2')
            # load_all("~/bin/monocle3") # monocle3_0.1.0 
            load_all("~/bin/2020_03_10_monocle3_repull/monocle3")

            library(tidyr)
            # library(ggplot2)
            # library(reshape2)
            library(stringr)
            library(Matrix)
            
            modStatus = "Modules Loaded"
    }) 
    }, error=function(err){
        print(paste("MY_ERROR:  ",err))
        print("Make sure you have loaded the right modules")
        print("Run a command such as:")
        print("module load R/3.6.1 xerces-c/3.1.1 gdal/2.4.1 geos/3.7.1 proj/4.9.3")
        modStatus = ("Module Loading Error")
    }
    )
    print(modStatus)
    return(modStatus)
}




loadMonoclePackages <- function(){

    modStatus <- tryCatch({
        suppressPackageStartupMessages({
            library(dplyr)
            library(SoupX)
            library(devtools)
            # library("uwot")
            library("irlba")
            library("Rtsne")
            library("slam")
            library("sp")
            library("sf")
            library("spdep")

            # Get custom uwot loaded
            # load_all("~/bin/uwot")
            # Use only the libpath leading to the core monocle installation
            # .libPaths('2')
            # load_all("~/bin/monocle3") # monocle3_0.1.0 
            # load_all("~/bin/2020_03_10_monocle3_repull/monocle3")
            library(monocle3)
            # load_all("/net/trapnell/vol1/home/readdf/bin/2020_09_07_monocle_repull_master/monocle3")

            library(tidyr)
            library(ggplot2)
            # library(reshape2)
            library(stringr)
            library(irlba)
            # library(Matrix)
            
            modStatus = "Modules Loaded"
    }) 
    }, error=function(err){
        print(paste("MY_ERROR:  ",err))
        print("Make sure you have loaded the right modules")
        print("Run a command such as:")
        print("module load R/3.6.1 xerces-c/3.1.1 gdal/2.4.1 geos/3.7.1 proj/4.9.3")
        modStatus = ("Module Loading Error")
    }
    )
    print(modStatus)
    return(modStatus)
}

# Function for loading a CDS
load.cds = function(mat.path, gene.annotation.path, cell.annotation.path) {
    df = read.table(
        mat.path,
        col.names = c("gene.idx", "cell.idx", "count"),
        colClasses = c("integer", "integer", "integer"))
    
    gene.annotations = read.table(
        gene.annotation.path,
        col.names = c("id", "gene_short_name"),
        colClasses = c("character", "character"))
    
    cell.annotations = read.table(
        cell.annotation.path,
        col.names = c("cell", "sample"),
        colClasses = c("character", "factor"))
    
    rownames(gene.annotations) = gene.annotations$id
    rownames(cell.annotations) = cell.annotations$cell
    
    # add a dummy cell to ensure that all genes are included in the matrix
    # even if a gene isn't expressed in any cell
    df = rbind(df, data.frame(
        gene.idx = c(1, nrow(gene.annotations)),
        cell.idx = rep(nrow(cell.annotations)+1, 2),
        count = c(1, 1)))
    
    mat = sparseMatrix(i = df$gene.idx, j = df$cell.idx, x = df$count)
    mat = mat[, 1:(ncol(mat)-1)]
    
    rownames(mat) = gene.annotations$id
    colnames(mat) = cell.annotations$cell

    cds = new_cell_data_set(as(mat, "sparseMatrix"), cell_metadata = cell.annotations,
                                gene_metadata = gene.annotations)   

    colData(cds)$n.umi = apply(exprs(cds), 2, sum) 

    return(cds)
}

# Function for loading a from a matrix market input
load.cds.MM = function(mat.path, gene.annotation.path, cell.annotation.path) {

    cds <- load_mm_data(mat_path=mat.path,
                        feature_anno_path=gene.annotation.path,
                        cell_anno_path=cell.annotation.path)

    return(cds)
}

makeCDS_fromBBIoutputDirs = function(sampleList, pathToSample="./"){
    # Make an empty list, fill it with cds's in named directories
    CDSlist = vector(mode="list", length=length(sampleList))
    #Loop to fill
    for (eachIndex in 1:length(sampleList)){
        inputFile = paste0(pathToSample, sampleList[eachIndex], "/", sampleList[eachIndex], "_cds.RDS")
        CDSlist[eachIndex] = readRDS(inputFile)
        CDSlist[[eachIndex]]$sampleName = sampleList[eachIndex]
    }
    # Now combine and return
    newCDS = monocle3::combine_cds(CDSlist)
    # Add UMIs
    colData(newCDS)$n.umi = Matrix::colSums(exprs(newCDS))
    return(newCDS)
}

# makeCDS_fromBBIoutputDirs_definedPath = function(sampleList, pathTo){
#     # Make an empty list, fill it with cds's in named directories
#     CDSlist = vector(mode="list", length=length(sampleList))
#     #Loop to fill
#     for (eachIndex in 1:length(sampleList)){
#         inputFile = paste0("./", sampleList[eachIndex], "/", sampleList[eachIndex], "_cds.RDS")
#         CDSlist[eachIndex] = readRDS(inputFile)
#         CDSlist[[eachIndex]]$sampleName = sampleList[eachIndex]
#     }
#     # Now combine and return
#     newCDS = monocle3::combine_cds(CDSlist)
#     # Add UMIs
#     colData(newCDS)$n.umi = Matrix::colSums(exprs(newCDS))
#     return(newCDS)
# }

medianWithoutNA<-function(x) {
   median(x[which(!is.na(x))])
}


makeGarnettModelMouse <- function(inputCDS, inputMarkerFilePath, inputMarkerFileName,
                 processingNote){
    # Load Garnett    
    load_all("~/bin/garnett")
    library(org.Mm.eg.db)

    # Initial test of marker utility
    marker_file_path = paste0(inputMarkerFilePath, inputMarkerFileName, ".txt")
    marker_check <- check_markers(inputCDS, marker_file_path,
                              db=org.Mm.eg.db,
                              # cds_gene_id_type = "SYMBOL",
                              cds_gene_id_type = "ENSEMBL",
                              marker_file_gene_id_type = "SYMBOL")

    # Save the marker_check output to see marker specificity
    png(paste0("./plots/marker_check", processingNote,
                 inputMarkerFileName, ".png"), width=1000, height=1000, res=200)
    print(plot_markers(marker_check) )
    dev.off()

    # Train a classifier
    set.seed(7)
    thisClassifier <- train_cell_classifier(cds = inputCDS,
                                        marker_file = marker_file_path,
                                        db=org.Mm.eg.db,
                                        cds_gene_id_type="ENSEMBL",
                                        marker_file_gene_id_type="SYMBOL")

    # Save the classifier. Flag name with the processing notes for clarity on the cell
    #    set used in training. Unweildy names, but will save trouble later (I hope)
    classifierOutput = paste0(inputMarkerFilePath, inputMarkerFileName, "_trainedOn_",
            processingNote, ".rds")
    print(paste0("Saving model to ", classifierOutput))
    saveRDS(thisClassifier, classifierOutput)

    return(inputCDS)
}
applyGarnettModelMouse <- function(inputCDS, inputModelFile){
     # Load Garnett    
    load_all("~/bin/garnett")
    library(org.Mm.eg.db)

    # Get the model
    thisClassifier <- readRDS(inputModelFile)

    # Classify
    inputCDS <- classify_cells(inputCDS, thisClassifier,
            db= org.Mm.eg.db,
            cluster_extend=TRUE,
            cds_gene_id_type="ENSEMBL")

    return(inputCDS)
}


makeGarnettModelHuman <- function(inputCDS, inputMarkerFilePath, inputMarkerFileName,
                 processingNote, returnPath=TRUE){
	# If the file exists already, just use that
	classifierOutput = paste0(inputMarkerFilePath, inputMarkerFileName, "_trainedOn_",
            processingNote, ".rds")
	if (file.exists(classifierOutput)){
		print("Model already exists, returning path to file")
		# If this exists, just return the path to the already-set file
       	return(classifierOutput)
	}

    # Load Garnett    
    load_all("~/bin/garnett")
    library(org.Hs.eg.db)

    # Initial test of marker utility
    marker_file_path = paste0(inputMarkerFilePath, inputMarkerFileName, ".txt")
    marker_check <- check_markers(inputCDS, marker_file_path,
                              db=org.Hs.eg.db,
                              # cds_gene_id_type = "SYMBOL",
                              cds_gene_id_type = "ENSEMBL",
                              marker_file_gene_id_type = "SYMBOL")
    # Save the marker_check output to see marker specificity
    png(paste0("./plots/marker_check", processingNote,
                 inputMarkerFileName, ".png"), width=1000, height=1000, res=200)
    print(plot_markers(marker_check) )
    dev.off()

    # Train a classifier
    set.seed(7)
    thisClassifier <- train_cell_classifier(cds = inputCDS,
                                        marker_file = marker_file_path,
                                        db=org.Hs.eg.db,
                                        cds_gene_id_type="ENSEMBL",
                                        marker_file_gene_id_type="SYMBOL")

    # Save the classifier. Flag name with the processing notes for clarity on the cell
    #    set used in training. Unweildy names, but will save trouble later (I hope)
    print(paste0("Saving model to ", classifierOutput))
    saveRDS(thisClassifier, classifierOutput)

    if (!returnPath){
        return(inputCDS)
    } else{
        return(classifierOutput)
    }
}

applyGarnettModelHuman <- function(inputCDS, inputModelFile){
     # Load Garnett    
    load_all("~/bin/garnett")
    library(org.Hs.eg.db)

    # Get the model
    thisClassifier <- readRDS(inputModelFile)

    # Classify
    inputCDS <- classify_cells(inputCDS, thisClassifier,
            db= org.Hs.eg.db,
            cluster_extend=TRUE,
            cds_gene_id_type="ENSEMBL")

    return(inputCDS)
}


boxplot_stat_by_cluster <- function(inputCDS, processingNote, statToPlot){
    # See doublet numbers (Scrublet Estimates)
    ggplot(as.data.frame(colData(inputCDS)), aes_string(x="cluster_label", y=statToPlot)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))# +
        # facet_grid(. ~ statToPlot)
    ggsave(paste0("./plots/", processingNote, statToPlot,  "_by_cluster.png"))
}


boxplot_stat_by_X <- function(inputCDS, processingNote, statToPlot, statToSplit,
                xLabToUse=statToSplit, yLabToUse=statToPlot, outputPath = "./plots/"){
    # See doublet numbers (Scrublet Estimates)
    ggplot(as.data.frame(colData(inputCDS)), aes_string(x=statToSplit, y=statToPlot)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        xlab(xLabToUse) + ylab(yLabToUse)
        # facet_grid(. ~ statToPlot)
    ggsave(paste0(outputPath, processingNote, statToPlot,  "_by_", statToSplit, ".png"))
}

boxplot_stat_by_X_and_fill <- function(inputCDS, processingNote, statToPlot, statToSplit, 
            statToFill){
    # See doublet numbers (Scrublet Estimates)
    ggplot(as.data.frame(colData(inputCDS)), 
        aes_string(x=statToSplit, y=statToPlot, fill=statToFill)) + 
        geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))# +
        # facet_grid(. ~ statToPlot)
    ggsave(paste0("./plots/", processingNote, statToPlot, 
         "_by_", statToSplit, "and", statToFill, ".png"))
}


makeScatterPlot <- function(inputCDS, processingNote, inputXaxis, inputYaxis){
    ggplot(as.data.frame(colData(inputCDS)), aes_string(x=inputXaxis, y=inputYaxis)) + 
        geom_point(size=1, alpha=.4) #+ theme(axis.text.x = element_text(angle = 90, hjust = 1))# +
        # facet_grid(. ~ statToPlot)
    ggsave(paste0("./plots/", processingNote, inputXaxis, "_vs_", inputYaxis,  ".png"))
}


runDEtestingToID_markers <- function(inputCDS, processingNote, testGroup,
									howManyGenesToTest = 25, outputPath="./plots/"){
    #Get top markers for the clusters in the post-SoupX-correction UMAP
    marker_test_res = top_markers(inputCDS,
             group_cells_by=testGroup, 
             genes_to_test_per_group=howManyGenesToTest,
             # reference_cells=500, # Maybe modify later. 1000 is default
              cores=8) 

    # Get top markers for the clusters
    top_specific_markers = marker_test_res %>%
        filter(fraction_expressing >= 0.10) %>%
        group_by(cell_group) %>%
        top_n(3, pseudo_R2)
    top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
    # Plot these markers
    png(paste0(outputPath, processingNote,"_", testGroup,  "three_markers_DE.png"),
             width=1400, height=1400, res=300)
    print(plot_genes_by_group(inputCDS,
                        top_specific_marker_ids,
                        group_cells_by=testGroup,
                        ordering_type="maximal_on_diag",
                        max.size=3))
    dev.off()

    top_specific_markers = marker_test_res %>%
        filter(fraction_expressing >= 0.10) %>%
        group_by(cell_group) %>%
        top_n(2, pseudo_R2)
    top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
    # Plot these markers
    png(paste0(outputPath, processingNote, "_", testGroup, "two_markers_DE.png"),
             width=1700, height=1700, res=300)
    print(plot_genes_by_group(inputCDS,
                        top_specific_marker_ids,
                        group_cells_by=testGroup,
                        ordering_type="maximal_on_diag",
                        max.size=2) )
    dev.off()

    # Sorted
    return (list("CDS"=inputCDS, "marker_test_res"=marker_test_res))
}



# Requires:
#     inputCDS is a monocle CDS object
#     colToAddName is a string that'll name a new colData column
#     colForDef is the name of a column already in inputCDS
#     nameAssignList is a list object with named entries, where each named value
#        is a vector matching the datatype in the colForDef col in inputCDS
# Modifies:
#   inputCDS
# Effects:
#     Makes a new column named colToAddName in inputCDS.
#     Loop through every row in the colForDef col in colData of inputCDS. If
#     the value matches is in one of the named values (vectors) in nameAssignList,
#     make that row's colToAddName value be the corresponding name
#     Non-matched entries in colToAddName left as "Unknown"
assignLabelFromColData <- function(inputCDS, colToAddName, colForDef, nameAssignList){
    # Initialize the empty column
    colData(inputCDS)[colToAddName] = "Unknown"

    # Loop through the name assign list
    for (eachName in names(nameAssignList)){
        colData(inputCDS)[,colToAddName] = with(colData(inputCDS), ifelse(
            (colData(inputCDS)[,colForDef] %in% (nameAssignList[[eachName]])),
                eachName, colData(inputCDS)[,colToAddName]
            ))
    }

    # Return the data
    return(inputCDS)
}

# Requires: 
#    inputCDS is a monocle CDS object.
#    subsetName and catToTest are col names in colData(inputCDS)
#    subsetToUse is in names(levels(colData(inputCDS)$subsetName))
#    if included, covarToIncl must be a vector of names of columns in colData(inputCDS)
# Modifies: Nothing
# Effects:
#    Takes only the cells in inputCDS where colData(inputCDS)$subsetName == subsetToUse
#    (eg grab the epithelial cells when that's the entry in the cell type column).
#    Then, within that subset, run DE testing based on categorizing by "catToTest", whatever
#    the various entries may be there.  
runDEtestSubsetByCateg <- function(inputCDS, subsetName, subsetToUse, 
                                        catToTest, baselineName, covarToIncl = NULL,
                                        catToCompareVsRest=NULL, 
                                        errorDistribution = "quasipoisson"){
    # Get the subset matching the specified substeName
    subsetCDS = inputCDS[,colData(inputCDS)[[subsetName]] == subsetToUse]

    # In some contexts I just want to test one sub-group against all the rest.
    if (!(is.null(catToCompareVsRest))){
        # Binarize catToTest to be a test of "is this cell from catToCompareVsRest"
        colData(subsetCDS)[[catToTest]] = ifelse(
            colData(subsetCDS)[[catToTest]] == catToCompareVsRest, catToCompareVsRest,
                "allOthers")
        baselineName = "allOthers"
    }

    # Make the catToTest into a factor and assign the designated baseline category for the
    #   lowest factor level
    colData(subsetCDS)[,catToTest] = relevel(as.factor(colData(subsetCDS)[[catToTest]]), baselineName)
    # colData(subsetCDS)[[catToTest]] = relevel(colData(subsetCDS)[catToTest], baselineName)

    # print(str(colData(subsetCDS)))
    print("Running DE Test Now")
    # Now run DE testing
    if (length(covarToIncl) > 0){
        covarFormatted =paste0(" + ", paste(covarToIncl, collapse=" + "))
    } else {
        covarFormatted = ""
    }

    catFits <- fit_models(subsetCDS,
                     model_formula_str = paste0("~", catToTest, covarFormatted),
                     expression_family=errorDistribution )
    # fit_coefs <- coefficient_table(catFits)

    # Return the coefficients
    return(catFits)
}


plotGeneExprViolins <- function(inputCDS, geneVec, geneSetName, processingNote){
    subsetCDS = inputCDS[rowData(inputCDS)$gene_short_name %in% geneVec,]
   
    # Plot using monocle functions
    png(paste0("./plots/", processingNote, geneSetName, "violin_expression_log.png"),
         width=1000, height=(220*length(geneVec)), res=200)
    print( plot_genes_violin(subsetCDS, group_cells_by="scrublet_call", 
                panel_order=geneVec) )
    dev.off()   

# Plot using monocle functions, linear scacle
    png(paste0("./plots/", processingNote, geneSetName, "violin_expression_no_log.png"),
         width=1000, height=(220*length(geneVec)), res=200)
    print( plot_genes_violin(subsetCDS, group_cells_by="scrublet_call", 
                panel_order=geneVec, log_scale=FALSE, min_expr=1) )
    dev.off()   

}


plotGeneExprViolins_setSplit <- function(inputCDS, geneVec, geneSetName, processingNote, splitBy){
    subsetCDS = inputCDS[rowData(inputCDS)$gene_short_name %in% geneVec,]
   
    # Plot using monocle functions
    png(paste0("./plots/", processingNote, geneSetName, "splitBy", splitBy, "violin_expression_log.png"),
         width=1000, height=(220*length(geneVec)), res=200)
    print( plot_genes_violin(subsetCDS, group_cells_by=splitBy, 
                panel_order=geneVec, pseudocount=1) )
    dev.off()   

# Plot using monocle functions, linear scacle
    png(paste0("./plots/", processingNote, geneSetName, "splitBy", splitBy, "violin_expression_no_log.png"),
         width=1000, height=(220*length(geneVec)), res=200)
    print( plot_genes_violin(subsetCDS, group_cells_by=splitBy, 
                panel_order=geneVec, log_scale=FALSE #, min_expr=1
                ) )
    dev.off()   
}

plotGeneExprViolinsNonzero_setSplit <- function(inputCDS, geneVec, geneSetName, processingNote, splitBy){
    subsetCDS = inputCDS[rowData(inputCDS)$gene_short_name %in% geneVec,]
   
    # Plot using monocle functions
    png(paste0("./plots/", processingNote, geneSetName, "splitBy", splitBy, "violin_Nonzero_expression_log.png"),
         width=1000, height=(220*length(geneVec)), res=200)
    print( plot_genes_violin(subsetCDS, group_cells_by=splitBy, 
                panel_order=geneVec, min_expr=.1) )
    dev.off()   

# Plot using monocle functions, linear scacle
    png(paste0("./plots/", processingNote, geneSetName, "splitBy", splitBy, "violin_Nonzero_expression_no_log.png"),
         width=1000, height=(220*length(geneVec)), res=200)
    print( plot_genes_violin(subsetCDS, group_cells_by=splitBy, 
                panel_order=geneVec, log_scale=FALSE, min_expr=.1
                ) )
    dev.off()   
}

# Deprecated version in favor of one below where I can specify text size
# plotUMAP_Monocle <- function(dataCDS, processingNote, catToColor,
#                     show_labels=TRUE){ #, xVal, yVal){
#     png(paste0("./plots/", processingNote, "_UMAP_", catToColor, "colored.png"),
#              width=1400, height=1000, res=200)
#     print(plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
#         color_cells_by=catToColor, label_cell_groups=show_labels,
#           cell_stroke=.1 #, group_label_size=4        
#                 ))
#     dev.off()   
# }
plotUMAP_Monocle <- function(dataCDS, processingNote, catToColor,
                    show_labels=TRUE, textSize=10, outputPath = "./plots/",
                    returnPlotObj=FALSE){ #, xVal, yVal){
    png(paste0(outputPath, processingNote, "_UMAP_", catToColor, "colored.png"),
             width=1400, height=1000, res=200)
    myPlot <- (plot_cells(dataCDS, reduction_method="UMAP",    # x=xVal, y=yVal,
        color_cells_by=catToColor, label_cell_groups=show_labels,
          cell_stroke=.1 , group_label_size=textSize        
                ))
    myPlot = (myPlot + theme(text=element_text(size=textSize)))
    print(myPlot)
    dev.off()   

    if (returnPlotObj){
        return(myPlot)
    }

}

plotUMAP_Monocle_genes <- function(inputCDS, processingNote, genesToPlot,
                     geneSetName, outputPath = "./plots/", plotSetTotals=FALSE){
    # Now plot to show expression of common marker genes
    png(paste0(outputPath, processingNote, "_UMAP_", geneSetName, ".png"),
             width=1400, height=1000, res=200)
    print(plot_cells(inputCDS, genes=genesToPlot,
                cell_stroke=.1,  label_cell_groups=FALSE) )
    dev.off()

    # Also, plot a UMAP coloring by log2 expression of the sum of expression of these genes
    if (plotSetTotals){
        subsetCDS = inputCDS[rowData(inputCDS)$gene_short_name %in% genesToPlot,]
        colData(subsetCDS)$gene_set_log2_UMI = Matrix::colSums(exprs(subsetCDS))
        colData(subsetCDS)$gene_set_log2_UMI = log2(.1 + 
                (colData(subsetCDS)$gene_set_log2_UMI / colData(subsetCDS)$Size_Factor ) )
      # Plot
        png(paste0(outputPath, processingNote, "_UMAP_", geneSetName, "totals.png"),
                 width=1400, height=1000, res=200)
        print(plot_cells(subsetCDS, color_cells_by="gene_set_log2_UMI",
                    cell_stroke=.1,  label_cell_groups=FALSE) )
        dev.off()
    }
}

plotUMAP_grey_out_subset <- function(inputCDS, processingNote, colorBy=NULL, greyBy=NULL, subsetToShow=NULL){
    # browser()
    # Get the CDS of data to show, and the one to grey out
    showCDS = inputCDS[,colData(inputCDS)[[greyBy]] == subsetToShow]
    greyOutCDS = inputCDS[,!(colData(inputCDS)[[greyBy]] == subsetToShow)]

    colData(showCDS)$umap1 = (reducedDims(showCDS)$UMAP)[,1]
    colData(showCDS)$umap2 = (reducedDims(showCDS)$UMAP)[,2]
    colData(greyOutCDS)$umap1 = (reducedDims(greyOutCDS)$UMAP)[,1]
    colData(greyOutCDS)$umap2 = (reducedDims(greyOutCDS)$UMAP)[,2]

    # Now, make a plot by combining 2 ggplot2 calls, modd'ing the color on the grey'd out one.
    figScale = 2.0
    pointSize=.3
    png(paste0("./plots/", processingNote, "UMAP_show_", colorBy, "fadeAllbut", subsetToShow, ".png"),
            width=2.5, height=1.75/figScale, units="in", res=300)
    print(ggplot(as.data.frame(colData(showCDS)), aes_string(x="umap1", y="umap2", colour=colorBy)) +   
        # viridis::scale_fill_viridis(option= "plasma") +

        geom_point(size=pointSize,stroke=0) + 

        # theme(legend.title = "Log10(Viral UMI)" 
        #     + 
        theme(legend.key.size = unit(.05, "in")) + 
            labs(colour="Log10(Viral UMI)")+
        geom_point(data=as.data.frame(colData(greyOutCDS)), colour="grey", size=pointSize,stroke=0) +
        theme(text=element_text(size=6)) + 
        # viridis::scale_color_viridis(option= "plasma") +
        # scale_color_gradient(limits=c(0,3.71)) +
        scale_color_viridis(limits=c(0,3.71), option="plasma") +
        monocle3:::monocle_theme_opts()   + coord_fixed()
            ) 
    dev.off()
}



plotGroupedProportions = function(inputCDS, processingNote,
                        groupValue, colForProportions,
                        pathToPlot="./plots/", widthToUse=1200,
                        heightToUse=800){

    countTable = with(colData(inputCDS), 
                table(get(groupValue), get(colForProportions)))
    proportionTable <- prop.table(countTable, margin=1)
    propDF = as.data.frame(proportionTable)
    colnames(propDF) = c(groupValue, colForProportions, "Freq")

    # Plot
    # png(paste0(pathToPlot, processingNote, "props", colForProportions, "by", groupValue, ".png"),
    #          width=widthToUse, height=heightToUse, res=200 )
    # barPlot <- ggplot(propDF, aes_string((colForProportions), 
    #                 "Freq",  fill=groupValue)) + 
    #   geom_col(position="dodge") 
    png(paste0(pathToPlot, processingNote, "props", colForProportions, "by", groupValue, ".png"),
             width=widthToUse, height=heightToUse, res=200 )
    barPlot <- ggplot(propDF, aes_string((groupValue), 
                    "Freq",  fill=colForProportions)) + 
             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      geom_col(position="dodge") 
    print(barPlot)
    dev.off()

    return(propDF)
}




makeQCplots = function(inputCDS, processingNote, makeUMAP=FALSE, genesToShow = c("HBB"),
                        geneSetName = "HemoglobinB"){
    # If haven't made plots for this set of pre-processing already, make it here. dir.create doesn't error if already exists, just ignores
    dir.create(file.path("./plots", paste0("QC_", processingNote)), showWarnings=FALSE)
    outputPath = paste0("./plots/QC_", processingNote, "/")
    # Make plots of mitochondrial RNA content, UMIs-per-sample, and cells-per-sample
    boxplot_stat_by_X(inputCDS, processingNote, 
        "log10_umi", "sampleName", outputPath=outputPath)
    boxplot_stat_by_X(inputCDS, processingNote, 
        "perc_mitochondrial_umis", "sampleName", outputPath=outputPath)

    # Make a UMAP, coloring by some QC stats
    if (makeUMAP){
        set.seed(7)

        if (!("UMAP" %in% names(reducedDims(inputCDS)))){
            inputCDS = estimate_size_factors(inputCDS)
            inputCDS = preprocess_cds(inputCDS)
            inputCDS = reduce_dimension(inputCDS)
        }

        # Now make the plots
        plotUMAP_Monocle(inputCDS, processingNote,"log10_umi", outputPath=outputPath)
        plotUMAP_Monocle(inputCDS, processingNote,"scrublet_score", outputPath=outputPath)
        plotUMAP_Monocle(inputCDS, processingNote,"perc_mitochondrial_umis", outputPath=outputPath)
        # Plot genes
        plotUMAP_Monocle_genes(inputCDS, processingNote, genesToShow, geneSetName, outputPath=outputPath)
    }
}




makeFacetedProportionPlot_withFill <- function(inputCDS, processingNote, colForGroup, colForProps, 
            fillCol = colForGroup, colorCol = colForGroup,
            subsetPropsToShow=levels(as.factor(colData(inputCDS)[[colForProps]])),
            outputPath="./plots/" ){
    # Get the proportions
    countTable = with(colData(inputCDS), table(get(colForGroup), get(colForProps)))
    proportionTable = prop.table(countTable, margin=1)
    propDF = as.data.frame(proportionTable)

    # Get the alignment of colForGroup and fillCol
    countTable = with(colData(inputCDS), table(get(colForGroup), get(fillCol)))
    proportionTable = prop.table(countTable, margin=1)
    fillDF = as.data.frame(proportionTable)
    # Keep the ones that are synonymous
    fillDF = fillDF[fillDF$Freq == 1,]
    # Get the name vec
    groupToFillVec = as.character(fillDF$Var2)
    names(groupToFillVec) = as.character(fillDF$Var1)

########################
    # Get the alignment of colForGroup and colorCol
    countTable = with(colData(inputCDS), table(get(colForGroup), get(colorCol)))
    proportionTable = prop.table(countTable, margin=1)
    fillDF = as.data.frame(proportionTable)
    # Keep the ones that are synonymous
    fillDF = fillDF[fillDF$Freq == 1,]
    # Get the name vec
    groupToColorVec = as.character(fillDF$Var2)
    names(groupToColorVec) = as.character(fillDF$Var1)
##########################
#
    # Subset for only the subsets to show (by default this is all) and re-add the fill info
    colnames(propDF) = c(colForGroup, colForProps, "Proportion")
    propDF[[colForProps]] = as.character(propDF[[colForProps]])
    propDF[[colForGroup]] = as.character(propDF[[colForGroup]])
    propDF = propDF[propDF[[colForProps]] %in% subsetPropsToShow,]
    propDF[[fillCol]] = groupToFillVec[propDF[[colForGroup]]]
    propDF[[colorCol]] = groupToColorVec[propDF[[colForGroup]]]

    # Make the plot
    png(paste0(outputPath, processingNote, colForGroup, "propsBy", colForProps,
                    "groupedBy", fillCol, "colorBy", colorCol, ".png"),
                    width=1600, height=1200, res=200 )
    myPlot <- (ggplot(propDF, 
        aes_string(x=fillCol, y="Proportion", color=colorCol,
                                 fill=colForProps)) +
        geom_point() +
         facet_wrap(vars(get(colForProps)), scales="free") +
        # facet_wrap(vars(get(colForProps)), scales="free") +
      # geom_boxplot() + 
        # geom_jitter() + 
        # facet_wrap(vars(get(colForProps)), scales="free") +
        ylab("Proportion of Total") + #xlab("") +
        theme(text = element_text(size = 12))    +
         theme(axis.text.x = element_text(angle = 90, hjust = 1))  
         #  +
         # theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
         )
    print(myPlot)
    dev.off()

    return(myPlot)
}



plot_genes_violin_DFR_customPseudo <- function (cds_subset,
                               group_cells_by = NULL,
                               min_expr = 0,
                               nrow = NULL,
                               ncol = 1,
                               panel_order = NULL,
                               label_by_short_name = TRUE,
                               normalize = TRUE,
                               log_scale = TRUE,
                               pseudocount = 0) {

  assertthat::assert_that(methods::is(cds_subset, "cell_data_set"))

  if(!is.null(group_cells_by)) {
    assertthat::assert_that(group_cells_by %in% names(colData(cds_subset)),
                            msg = paste("group_cells_by must be a column in",
                                        "the colData table"))
  }

  assertthat::assert_that(assertthat::is.number(min_expr))

  if(!is.null(nrow)) {
    assertthat::assert_that(assertthat::is.count(nrow))
  }

  assertthat::assert_that(assertthat::is.count(ncol))
  assertthat::assert_that(assertthat::is.number(pseudocount))
  assertthat::assert_that(is.logical(label_by_short_name))
  if (label_by_short_name) {
    assertthat::assert_that("gene_short_name" %in% names(rowData(cds_subset)),
                            msg = paste("When label_by_short_name = TRUE,",
                                        "rowData must have a column of gene",
                                        "names called gene_short_name."))
  }
  if(!is.null(panel_order)) {
    if (label_by_short_name) {
      assertthat::assert_that(all(panel_order %in%
                                    rowData(cds_subset)$gene_short_name))
    } else {
      assertthat::assert_that(all(panel_order %in%
                                    row.names(rowData(cds_subset))))
    }
  }

  assertthat::assert_that(is.logical(normalize))
  assertthat::assert_that(is.logical(log_scale))


  assertthat::assert_that(nrow(rowData(cds_subset)) <= 100,
                          msg = paste("cds_subset has more than 100 genes -",
                                      "pass only the subset of the CDS to be",
                                      "plotted."))

  if (pseudocount > 0) {
    cds_exprs <- SingleCellExperiment::counts(cds_subset) + pseudocount
  } else {
    cds_exprs <- SingleCellExperiment::counts(cds_subset)
  }
  if (normalize) {
    cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  } else {
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  }

  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr

  cds_exprs <- merge(cds_exprs, rowData(cds_subset), by.x = "f_id",
                     by.y = "row.names")
  # cds_exprs <- merge(cds_exprs, colData(cds_subset), by.x = "Cell",  by.y = "row.names")
  # browser()
  # row.names(x)
  # cds_exprs$Cell = as.character(cds_exprs$Cell)
  cds_exprs <- merge(as.data.frame(cds_exprs), as.data.frame(colData(cds_subset)),
         by.x = "Cell",  by.y = "row.names")

  if (label_by_short_name) {
    if (!is.null(cds_exprs$gene_short_name)) {
      cds_exprs$feature_label <- cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    } else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  } else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }

  if (!is.null(panel_order)) {
    cds_exprs$feature_label = factor(cds_exprs$feature_label,
                                     levels = panel_order)
  }

  cds_exprs[,group_cells_by] <- as.factor(cds_exprs[,group_cells_by])

  q <- ggplot(aes_string(x = group_cells_by, y = "expression"),
              data = cds_exprs) +
    monocle_theme_opts()

  cds_exprs[,group_cells_by] <- as.factor(cds_exprs[,group_cells_by])
  q <- q + geom_violin(aes_string(fill = group_cells_by), scale="width") +
    guides(fill=FALSE)
  q <- q + stat_summary(fun=mean, geom="point", size=1, color="black")
  q <- q + facet_wrap(~feature_label, nrow = nrow,
                      ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }

  q <- q + ylab("Expression") + xlab(group_cells_by)

  if (log_scale){
    q <- q + scale_y_log10()
  }
  q
}



   
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




