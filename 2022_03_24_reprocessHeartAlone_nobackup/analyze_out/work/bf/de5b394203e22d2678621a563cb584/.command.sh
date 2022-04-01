#!/bin/bash -ue
PROCESS_BLOCK='mergeBamsProcess'
   SAMPLE_NAME="W146.heart.apex.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

outBam="W146.heart.apex.s1-merged.bam"

   if [ 4 -gt 1 ]
   then
       sambamba merge --nthreads 8 ${outBam} W146.heart.apex.s1-RUN001_L003.bam W146.heart.apex.s1-RUN001_L001.bam W146.heart.apex.s1-RUN001_L002.bam W146.heart.apex.s1-RUN001_L004.bam
   else
       cp W146.heart.apex.s1-RUN001_L003.bam W146.heart.apex.s1-RUN001_L001.bam W146.heart.apex.s1-RUN001_L002.bam W146.heart.apex.s1-RUN001_L004.bam ${outBam}
   fi

samtools index ${outBam}

   mkdir -p /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/W146.heart.apex.s1/genome_browser
   pushd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/W146.heart.apex.s1/genome_browser
   ln -sf ../merge_bams/${outBam} .
   ln -sf ../merge_bams/${outBam}.bai .
   popd

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'sambamba --version 2>&1 > /dev/null | head -2 | tail -1' 'samtools --version | head -2'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
