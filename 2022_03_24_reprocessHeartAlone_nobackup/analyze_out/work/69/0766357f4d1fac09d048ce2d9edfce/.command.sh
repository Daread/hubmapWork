#!/bin/bash -ue
PROCESS_BLOCK='makeMergedPeakRegionCountsProcess'
   SAMPLE_NAME="W146.heart.apex.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/load_python_env_reqs.sh
   source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/python_env/bin/activate

outCounts="W146.heart.apex.s1-peak_counts.txt"
tmpRegions="W146.heart.apex.s1-temp_regions.gz"

   # TODO simplify this... maybe just have a stage that makes this file for TSS rather than complicating the stage itself
   bedtools slop -i W146.heart.apex.s1-merged_peaks.bed -g W146.heart.apex.s1-human.chromosome_sizes.sorted.txt -b 0     | bedtools merge -i stdin     | gzip > ${tmpRegions}

   python /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/get_region_counts.py         --transposition_sites_intersect <(bedtools intersect -sorted -a W146.heart.apex.s1-transposition_sites.bed.gz -b ${tmpRegions})         --output_file ${outCounts}

   rm ${tmpRegions}

   deactivate

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'bedtools --version' 'gzip --version | head -1' 'python3 --version'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
