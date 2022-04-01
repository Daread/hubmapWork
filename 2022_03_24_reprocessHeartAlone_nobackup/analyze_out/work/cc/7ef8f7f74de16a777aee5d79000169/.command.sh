#!/bin/bash -ue
PROCESS_BLOCK='getPerBaseCoverageTssProcess'
   SAMPLE_NAME="W137.heart.apex.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

outCoverage="W137.heart.apex.s1-tss_region_coverage.txt.gz"
tmpOut="W137.heart.apex.s1-temp_file.gz"

   # First get 2kb regions surrounding TSSs (not strand-specific here)
   # then calculate per-base coverage with bedtools
   # then write any non-zero entries to a file
   bedtools slop -i W137.heart.apex.s1-human.tss_file.sorted.bed.gz -g W137.heart.apex.s1-human.chromosome_sizes.sorted.txt -b 1000      | bedtools coverage -sorted -d -a stdin -b W137.heart.apex.s1-transposition_sites.bed.gz     | awk '{{ if ($8 > 0) print $0 }}'     | gzip > ${tmpOut}

   # Aggregate per-position coverage over all positions across genes, taking strand into account
   Rscript /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/aggregate_per_base_tss_region_counts.R ${tmpOut} ${outCoverage}

   rm ${tmpOut}

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'bedtools --version' 'awk --version | head -1' 'gzip --version | head -1' 'R --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
