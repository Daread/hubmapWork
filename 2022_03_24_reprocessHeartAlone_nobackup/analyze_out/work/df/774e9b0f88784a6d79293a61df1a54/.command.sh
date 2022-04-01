#!/bin/bash -ue
PROCESS_BLOCK='sortTssBedProcess'
   SAMPLE_NAME="na"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

outBed="W142.heart.LV.s1-human.tss_file.sorted.bed.gz"

zcat /net/bbi/vol1/data/genomes_stage/human/human_atac/tss.bed.gz | sort -k1,1V -k2,2n -k3,3n | gzip > ${outBed}

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'zcat --version | head -1' 'sort --version | head -1' 'gzip --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
