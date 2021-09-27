

# Format in the runall.sh script
# coaccessCut=$1 # Cicero co-accessibility cutoff to count
# maxSites=$2 # Max number of sites to link
# peakSize=$3 # Size of distal sites
# upstream=$4
# downstream=$5



# 9-13-21 
# Test sending off batch submissions this way
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.05 0 600 0 0


qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.05 5 600 0 0


qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.035 20 600 0 0


qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.05 20 600 0 0




# 9-27-21: Sawa odd changes in my results from running with no linked sites...
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.95 0 600 0 0

qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.05 0 3000 0 0
