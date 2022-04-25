#!/bin/bash -ue
PROCESS_BLOCK='makeCountReportsProcess'
   SAMPLE_NAME="W137.heart.apex.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

outCountReport="W137.heart.apex.s1-count_report.txt"

   Rscript /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/make_count_report.R W137.heart.apex.s1-duplicate_report.txt W137.heart.apex.s1-peak_counts.txt W137.heart.apex.s1-tss_counts.txt ${outCountReport}

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'R --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
