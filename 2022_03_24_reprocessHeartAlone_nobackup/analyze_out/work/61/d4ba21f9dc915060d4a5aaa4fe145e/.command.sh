#!/bin/bash -ue
PROCESS_BLOCK='makeWindowedGenomeIntervalsProcess'
   SAMPLE_NAME="W142.heart.LV.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

outBed="W142.heart.LV.s1-genomic_windows.bed"

   bedtools makewindows         -g W142.heart.LV.s1-human.chromosome_sizes.sorted.txt         -w 5000         > ${outBed}

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'bedtools --version'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
