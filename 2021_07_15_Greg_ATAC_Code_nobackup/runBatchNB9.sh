 




qsub -l mfree=30G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh peakMat useMNN 2 50

qsub -l mfree=30G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh bMat useMNN 2 50




qsub -l mfree=30G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh peakMat useMNN 1 50
qsub -l mfree=30G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh peakMat useMNN 1 15
qsub -l mfree=30G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh peakMat useMNN 2 15



# 9-30-21: Run with varying FRIP and subsets of samples to use
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W144.heart.apex.s1 useMNN .2
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W144.heart.apex.s1 noMNN .2

qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W144.heart.apex.s1 useMNN .3

qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W146.heart.LV.s1 useMNN .2
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W144.heart.apex.s1 noMNN .3
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W146.heart.LV.s1 noMNN .2

qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W146.heart.LV.s1 useMNN .3
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W146.heart.LV.s1 noMNN .3



# Try to move the FRIP cutoff even further up?
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W146.heart.LV.s1 useMNN .4
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W146.heart.LV.s1 noMNN .4

qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W146.heart.LV.s1 useMNN .5
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W146.heart.LV.s1 noMNN .5



# Try to move the FRIP cutoff even further up?
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W144.heart.apex.s1 useMNN .4
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W144.heart.apex.s1 noMNN .4

qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W144.heart.apex.s1 useMNN .5
qsub -l mfree=40G -l centos=7 -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB9.sh W144.heart.apex.s1 noMNN .5





