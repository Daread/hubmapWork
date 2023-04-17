
module load bedtools/2.29.2

# echo "Starting peak matching and sequence extraction"
# # First, set wd
cd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup
# # For a defined promoter size and cicero linking hyperparams, get site info for each gene
Rscript definedPromoter_get_sites.R --promoterUpstream $1 --promoterDownstream $2 --coaccessCutoff $3 --maxNdistalSites $4 --peakSize $5 --linkSet $6 --seqInfOnly TRUE 
