#!/bin/bash -ue
PROCESS_BLOCK='combineReadCountsProcess'
    SAMPLE_NAME="W144.heart.apex.s1"
    START_TIME=`date '+%Y%m%d:%H%M%S'`

    outDuplicateReport="W144.heart.apex.s1-combined.duplicate_report.txt"

	python /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/combine_read_counts.py --input_duplicate_report W144.heart.apex.s1-duplicate_report.txt --input_mito_duplicate_report W144.heart.apex.s1-mito.duplicate_report.txt --output_combined_duplicate_report ${outDuplicateReport}
echo foo

    STOP_TIME=`date '+%Y%m%d:%H%M%S'`
    /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
