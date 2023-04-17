#! /bin/bash

module load bedtools/2.29.2

# echo "Starting peak matching and sequence extraction"
# # First, set wd
cd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup
# # For a defined promoter size and cicero linking hyperparams, get site info for each gene
Rscript definedPromoter_get_sites.R --promoterUpstream $1 --promoterDownstream $2 --coaccessCutoff $3 --maxNdistalSites $4 --peakSize $5 --linkSet $6

echo "Running FIMO now"
# Now run FIMO
cd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoFIMO_nobackup
./runFIMO.sh Gene_Prom_Plus_Distal_WithSequence_Sites_Max$4_Upstream$1_Downstream$2_cicCuf$3peakSize$5_Links_$6.fa

echo "FIMO outputs to matrix"
# Now run the processing of these FIMO outputs into matrices
cd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup
Rscript gatherFIMOresults.R  --promoterUpstream $1 --promoterDownstream $2 --coaccessCutoff $3 --maxNdistalSites $4 --peakSize $5 --pValFIMOcutoff 1e-6 --linkSet $6
Rscript gatherFIMOresults.R  --promoterUpstream $1 --promoterDownstream $2 --coaccessCutoff $3 --maxNdistalSites $4 --peakSize $5 --pValFIMOcutoff 1e-5 --linkSet $6
Rscript gatherFIMOresults.R  --promoterUpstream $1 --promoterDownstream $2 --coaccessCutoff $3 --maxNdistalSites $4 --peakSize $5 --pValFIMOcutoff 1e-4 --linkSet $6

echo "Now submitting jobs for glmnet"
# Now run the model fitting
cd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup

for PVAL in 1e-4 1e-5 1e-6
do 
	qsub -P trapnelllab -l mfree=30G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/singleFittingRun.sh $1 $2 $3 $4 $5 Binary_Combined_Motif_Counts _log2_CPM $PVAL $6
	qsub -P trapnelllab -l mfree=30G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/singleFittingRun.sh $1 $2 $3 $4 $5 Binary_PromOnly_Motif_Counts _log2_CPM $PVAL $6
	# qsub -P trapnelllab -l mfree=30G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/singleFittingRun.sh $1 $2 $3 $4 $5 Binary_Combined_Motif_Counts _log2_ratio_Vs_AllTypeMean $PVAL
	# qsub -P trapnelllab -l mfree=30G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/singleFittingRun.sh $1 $2 $3 $4 $5 Binary_PromOnly_Motif_Counts _log2_ratio_Vs_AllTypeMean $PVAL
done 

# # Classification models as well?
# for PVAL in 1e-4 1e-5 1e-6
# do 
# 	qsub -P trapnelllab -l mfree=30G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/singleClassificationRun.sh $1 $2 $3 $4 $5 Binary_Combined_Motif_Counts _log2_CPM $PVAL
# 	qsub -P trapnelllab -l mfree=30G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoMotifModel_nobackup/singleClassificationRun.sh $1 $2 $3 $4 $5 Binary_PromOnly_Motif_Counts _log2_CPM $PVAL
# done 




