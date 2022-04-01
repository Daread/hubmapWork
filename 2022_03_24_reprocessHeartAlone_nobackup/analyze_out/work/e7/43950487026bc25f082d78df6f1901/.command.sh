#!/bin/bash -ue
PROCESS_BLOCK='mergePeaksByFileProcess'
SAMPLE_NAME="W144.heart.apex.s1"
START_TIME=`date '+%Y%m%d:%H%M%S'`

outBed="W144.heart.apex.s1-merged_peaks.bed"
cat inBeds*         | cut -f1-3         | sort -k1,1V -k2,2n -k3,3n         | bedtools merge -i -         | sort -k1,1V -k2,2n -k3,3n > ${outBed}

outBedList="W144.heart.apex.s1-merge_by_file_beds.txt"
echo "W144.heart.apex.s1-group_merged_peaks.bed" > ${outBedList}

STOP_TIME=`date '+%Y%m%d:%H%M%S'`
/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'cut --version | head -1' 'sort --version | head -1' 'bedtools --version'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
