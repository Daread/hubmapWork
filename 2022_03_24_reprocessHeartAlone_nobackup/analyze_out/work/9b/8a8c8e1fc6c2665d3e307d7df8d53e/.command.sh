#!/bin/bash -ue
PROCESS_BLOCK='summarizeCellCallsProcess'
SAMPLE_NAME="W134.heart.apex.s1"
START_TIME=`date '+%Y%m%d:%H%M%S'`

outSummaryPlot="W134.heart.apex.s1-called_cells_summary.pdf"
outSummaryStats="W134.heart.apex.s1-called_cells_summary.stats.txt"

barnyardParams=""
if [ "0" = 1 ]
then
    barnyardParams="--window_matrices W134.heart.apex.s1-window_matrix.mtx.gz --barnyard"
fi

if [ "group_1" != "" ]
then
  PEAK_GROUP="group_1"
else
  PEAK_GROUP="-"
fi

if [ "" != "" ]
then
  PEAK_FILE="`basename `"
else
  PEAK_FILE="-"
fi

Rscript /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/summarize_cell_calls.R         --sample_name W134.heart.apex.s1         --read_count_tables W134.heart.apex.s1-count_report.txt         --stats_files W134.heart.apex.s1-called_cells_stats.json         --insert_size_tables W134.heart.apex.s1-insert_sizes.txt         --peak_groups ${PEAK_GROUP}         --peak_files ${PEAK_FILE}         --peak_call_files W134.heart.apex.s1-peaks.narrowPeak.gz         --merged_peaks W134.heart.apex.s1-merged_peaks.bed         --per_base_tss_region_coverage_files W134.heart.apex.s1-tss_region_coverage.txt.gz         --combined_duplicate_report W134.heart.apex.s1-combined.duplicate_report.txt         --plot ${outSummaryPlot}         --output_stats ${outSummaryStats} ${barnyardParams}

STOP_TIME=`date '+%Y%m%d:%H%M%S'`
/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'R --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
