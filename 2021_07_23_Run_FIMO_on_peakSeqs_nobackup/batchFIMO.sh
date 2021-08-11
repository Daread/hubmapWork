#! /bin/bash



qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_23_Run_FIMO_on_peakSeqs_nobackup /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_23_Run_FIMO_on_peakSeqs_nobackup/runFIMO.sh cds_p_allHeartATACFRIP=0.2_FRIT=0.05UMI=1000DL=0.5_peakSeqsAsFasta600.fa

