
#! /bin/bash

# 7-21-21 run
for RNANAME in W134.Apex W135.Left.Vent W136.Apex W136.Left.Vent W139.Apex W139.Left.Vent W139.Right.Vent W139.Septum W142.Left.Vent W144.Apex W145.Apex W145.Left.Vent W146.Apex W146.Left.Vent; do
	qsub -P trapnelllab -pe serial 16 -l centos=7 -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB6.5.sh $RNANAME
done




# 7-30-21 run
for RNANAME in W134.Apex W135.Left.Vent W144.Apex W145.Apex W146.Apex W146.Left.Vent; do
	qsub -P trapnelllab -pe serial 16 -l centos=7 -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB6.5.sh $RNANAME
done





qsub -P trapnelllab -pe serial 16 -l centos=7 -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB6.5.sh All_Cells





# 8-2-21 run
for RNANAME in All_Cells W134.Apex W135.Left.Vent W136.Apex W136.Left.Vent W139.Apex W139.Left.Vent W139.Right.Vent W139.Septum W142.Left.Vent W144.Apex W145.Apex W145.Left.Vent W146.Apex W146.Left.Vent; do
	qsub -P trapnelllab -pe serial 16 -l centos=7 -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB6.5.sh $RNANAME
done





qsub -P trapnelllab -pe serial 16 -l centos=7 -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB6.5.sh All_Cells FRIP=0.2_FRIT=0.05UMI=1000DL=0.7
qsub -P trapnelllab -pe serial 16 -l centos=7 -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB6.5.sh All_Cells FRIP=0.2_FRIT=0.05UMI=1000DL=0.6
qsub -P trapnelllab -pe serial 16 -l centos=7 -l mfree=20G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/runSingleNB6.5.sh All_Cells FRIP=0.2_FRIT=0.05UMI=1000DL=0.5


