#! /bin/bash

# # Format of arguments:
# # --promoterUpstream $1 --promoterDownstream $2 --coaccessCutoff $3 --maxNdistalSites $4 --peakSize $5
# qsub -P trapnelllab -l mfree=40G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/singleProcessingPipeline.sh 1500 500 0.035 5 600

 

#########################################################################

# 4-19-2022
# for maxSite in 5 10 20
# do
# 	for cutoff in 0.015 0.035 0.05
# 	do
# 		for peak in 600 1000
# 		do
# 			qsub -P trapnelllab -l mfree=32G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/singleProcessingPipeline.sh 1500 500 $cutoff $maxSite $peak
# 		done
# 	done
# done

# for maxSite in 5 10 20
# do
# 	for cutoff in 0.015 0.035 0.05
# 	do
# 		for peak in 600 1000
# 		do
# 			qsub -P trapnelllab -l mfree=32G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/singleProcessingPipeline.sh 2000 1000 $cutoff $maxSite $peak
# 		done
# 	done
# done




qsub -P trapnelllab -l mfree=64G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/singleProcessingPipeline.sh 1000 200 0.08 1 100


qsub -P trapnelllab -l mfree=64G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/singleProcessingPipeline.sh 5000 2000 0.08 1 100


###################################################################################################################


# qsub -P trapnelllab -l mfree=64G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/singleProcessingPipelinePostFIMO.sh 1000 200 .05 20 600














