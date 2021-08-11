#!/usr/bin/env Rscript

print("Hello world")

loadMonoclePackages <- function(){

    modStatus <- tryCatch({
        suppressPackageStartupMessages({
            # Adding this 5-14-21
            # .libPaths('2')

        	library(argparse)
        	library(Matrix)
            # library(dplyr)
            # # library(SoupX)
           # library(devtools)
            # # library("uwot")
            # library("irlba")
            # library("Rtsne")
            # library("slam")
            # library("sp")
            # library("sf")
            # # load_all("~/R/x86_64-pc-linux-gnu-library/4.0/sf")
            # # load_all("~/R/x86_64-pc-linux-gnu-library/4.0/spdep")
            # # load_all("~/R/x86_64-pc-linux-gnu-library/4.0/S4Vectors")
            # library("spdep")
            # # Get custom uwot loaded
            # # load_all("~/bin/uwot")
            # # Use only the libpath leading to the core monocle installation


            # 5-14-21: commenting out the .libPaths('2') line
            .libPaths('2')
            library(monocle3)
            # # load_all("~/bin/monocle3") # monocle3_0.1.0 
            # # load_all("~/bin/2020_03_10_monocle3_repull/monocle3")
            #load_all("~/R/x86_64-pc-linux-gnu-library/4.0/monocle3")
            # load_all("/net/trapnell/vol1/home/readdf/bin/2020_09_07_monocle_repull_master/monocle3")

            # library(tidyr)
            # library(ggplot2)
            # # library(reshape2)
            # library(stringr)
            # library(irlba)
            # # library(Matrix)
            modStatus = "Modules Loaded"
    }) 
    }, error=function(err){
        print(paste("MY_ERROR:  ",err))
        print(.libPaths())
        print("Make sure you have loaded the right modules")
        print("Run a command such as:")
        print("module load R/3.6.1 xerces-c/3.1.1 gdal/2.4.1 geos/3.7.1 proj/4.9.3")
        find.package("monocle3")
        find.package("Matrix")
        find.package("argparse")
        modStatus = ("Module Loading Error")
    }
    )
    print(modStatus)
    return(modStatus)
}

loadMonoclePackages()


print("All Done")