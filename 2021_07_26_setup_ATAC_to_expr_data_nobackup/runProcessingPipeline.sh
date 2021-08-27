#! /bin/bash

# Format of arguments:
# --promoterUpstream $1 --promoterDownstream $2 --coaccessCutoff $3 --maxNdistalSites $4 --peakSize $5
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/singleProcessingPipeline.sh 1500 500 0.035 5 600

 
# Submitting 8-25-21: Vary parameters now that I've got the full pipeline actually working

# 
# Allowing more distal sites for low-cutoff approaches
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/singleProcessingPipeline.sh 1500 500 0.035 10 600
# Lower cutoffs but otherwise match my initial runs
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/singleProcessingPipeline.sh 1000 200 0.05 5 600
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/singleProcessingPipeline.sh 1000 200 0.035 5 600
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/singleProcessingPipeline.sh 1000 200 0.015 5 600
# Vary cutoffs to match my first batch submission
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/singleProcessingPipeline.sh 1500 500 0.05 5 600
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/singleProcessingPipeline.sh 1500 500 0.035 5 600
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/singleProcessingPipeline.sh 1500 500 0.015 5 600
# Keep a more stringent cutoff but vary longer 'promoter' regions
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/singleProcessingPipeline.sh 3000 1000 0.05 5 600
qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_26_setup_ATAC_to_expr_data_nobackup/singleProcessingPipeline.sh 5000 2000 0.05 5 600




