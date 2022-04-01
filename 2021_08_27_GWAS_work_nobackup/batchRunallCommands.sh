

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




# 9-27-21: Fixed 

# Verify the no-distal runs are looking identical
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.95 0 600 0 0 2 .1
# qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.05 0 600 0 0 2 .1

qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.05 0 3000 0 0 2 .1

# Add distal sites
# qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 2 .1

qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 4 0.1
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 6 0.1
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 8 0.1


qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 2 0.1

qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 2 0.1

qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 8 0.1



qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 4 0.1
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 6 0.1
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 8 0.1



# 9-29-21:
# I want to vary the log2(ratio) cutoff to take a look more at specific genes w/in a cell type, rather than generally expressed ones
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 4 2
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 4 2
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 2 3
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 2 3

# I want to vary the log2(ratio) cutoff to take a look more at specific genes w/in a cell type, rather than generally expressed ones
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 4 1
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 4 1

# Try lower expression cutoffs as well (2 -> barely expressed but clearly out of the non-expressed portion of bimodal distribution)
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 2 1
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 2 1

qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 2 2
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 2 2


# 9-30-21

# Try lower expression cutoffs as well (2 -> barely expressed but clearly out of the non-expressed portion of bimodal distribution)
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 3 1
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 3 1



# Try lower expression cutoffs as well (2 -> barely expressed but clearly out of the non-expressed portion of bimodal distribution)
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 3 1.585
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 3 1.585

qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 2 1.585
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 2 1.585


qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 4 1.585
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 4 1.585

# Try lower expression cutoffs as well (2 -> barely expressed but clearly out of the non-expressed portion of bimodal distribution)
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 5 1
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 5 1


# Try lower expression cutoffs as well (2 -> barely expressed but clearly out of the non-expressed portion of bimodal distribution)
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 20 600 0 0 6 1
qsub -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/runall.sh 0.015 0 600 0 0 6 1

