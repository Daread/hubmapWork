#! /bin/bash



# 8-18-21 Runs
Rscript runGlmnetOnMotifs.R --featureSelection Binary_PromOnly_Motif_Counts
Rscript runGlmnetOnMotifs.R --featureSelection Binary_Combined_Motif_Counts




# 8-18-21 Runs
Rscript runGlmnetOnMotifs.R --featureSelection Binary_PromOnly_Motif_Counts --predictionTask _log2_CPM
Rscript runGlmnetOnMotifs.R --featureSelection Binary_Combined_Motif_Counts --predictionTask _log2_CPM


# More 8-18
Rscript runGlmnetOnMotifs.R --featureSelection Nonbinary_Combined_Motif_Counts --predictionTask _log2_CPM
Rscript runGlmnetOnMotifs.R --featureSelection Nonbinary_Uncombined_Motif_Counts --predictionTask _log2_CPM





# 8-18-21 Runs
Rscript runGlmnetOnMotifs.R --featureSelection Binary_Combined_Motif_Counts --predictionTask _log2_CPM --pValFIMOcutoff 1e-5
Rscript runGlmnetOnMotifs.R --featureSelection Binary_Combined_Motif_Counts --predictionTask _log2_ratio_Vs_AllTypeMean --pValFIMOcutoff 1e-5




# 8-19-21: Get organized and start submitting jobs
qsub -P trapnelllab -l mfree=30G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/singleFittingRun.sh Binary_Combined_Motif_Counts _log2_CPM  1e-5
qsub -P trapnelllab -l mfree=30G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/singleFittingRun.sh Binary_Combined_Motif_Counts _log2_ratio_Vs_AllTypeMean  1e-5
qsub -P trapnelllab -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/singleFittingRun.sh Binary_PromOnly_Motif_Counts _log2_CPM  1e-5
qsub -P trapnelllab -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/singleFittingRun.sh Binary_PromOnly_Motif_Counts _log2_ratio_Vs_AllTypeMean  1e-5


qsub -P trapnelllab -l mfree=30G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/singleFittingRun.sh Binary_Combined_Motif_Counts _log2_CPM  5e-6
qsub -P trapnelllab -l mfree=30G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/singleFittingRun.sh Binary_Combined_Motif_Counts _log2_ratio_Vs_AllTypeMean  5e-6
qsub -P trapnelllab -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/singleFittingRun.sh Binary_PromOnly_Motif_Counts _log2_CPM  5e-6
qsub -P trapnelllab -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_18_motifElasticnet_nobackup/singleFittingRun.sh Binary_PromOnly_Motif_Counts _log2_ratio_Vs_AllTypeMean  5e-6





