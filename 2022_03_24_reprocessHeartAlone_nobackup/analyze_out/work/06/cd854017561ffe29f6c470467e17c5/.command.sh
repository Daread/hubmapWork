#!/bin/bash -ue
PROCESS_BLOCK='mergeMitoBamsProcess'
SAMPLE_NAME="W139.heart.LV.s1"
START_TIME=`date '+%Y%m%d:%H%M%S'`

outMitoBam="W139.heart.LV.s1-merged.mito.bam"

sambamba merge --nthreads 8 ${outMitoBam} W139.heart.LV.s1-RUN001_L003.mito.bam W139.heart.LV.s1-RUN001_L001.mito.bam W139.heart.LV.s1-RUN001_L004.mito.bam W139.heart.LV.s1-RUN001_L002.mito.bam
samtools index ${outMitoBam}

mkdir -p /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/W139.heart.LV.s1/genome_browser
pushd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/W139.heart.LV.s1/genome_browser
ln -sf ../merge_bams/${outMitoBam} .
ln -sf ../merge_bams/${outMitoBam}.bai .
popd

STOP_TIME=`date '+%Y%m%d:%H%M%S'`
/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'sambamba --version 2>&1 > /dev/null | head -2 | tail -1' 'samtools --version | head -2'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
