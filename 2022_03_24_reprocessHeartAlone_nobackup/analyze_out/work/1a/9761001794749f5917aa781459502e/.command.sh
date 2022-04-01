#!/bin/bash -ue
PROCESS_BLOCK='mergePeaksByGroupProcess'
   SAMPLE_NAME="peaks.group_1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

outGroupBed="group_1-group_merged_peaks_set.bed"

   zcat inBeds*         | cut -f1-3         | sort -k1,1V -k2,2n -k3,3n         | bedtools merge -i -         | sort -k1,1V -k2,2n -k3,3n > ${outGroupBed}

   #
   # Make copies because Nextflow cannot pass on symbolic links.
   #
   for outSample in W134.heart.apex.s1 W135.heart.LV.s1 W136.heart.apex.s1 W136.heart.LV.s1 W137.heart.apex.s1 W139.heart.LV.s1 W139.heart.RV.s1 W139.heart.septum.s1 W139.heart.apex.s1 W142.heart.LV.s1 W144.heart.apex.s1 W145.heart.apex.s1 W145.heart.LV.s1 W146.heart.apex.s1 W146.heart.LV.s1
   do
     outBed="${outSample}-group_merged_peaks.bed"
     cp ${outGroupBed} ${outBed}

     outBedList="${outSample}-merge_by_group_beds.txt"
     echo "W146.heart.LV.s1-peaks.narrowPeak.gz W145.heart.LV.s1-peaks.narrowPeak.gz W144.heart.apex.s1-peaks.narrowPeak.gz W145.heart.apex.s1-peaks.narrowPeak.gz W142.heart.LV.s1-peaks.narrowPeak.gz W135.heart.LV.s1-peaks.narrowPeak.gz W134.heart.apex.s1-peaks.narrowPeak.gz W139.heart.LV.s1-peaks.narrowPeak.gz W139.heart.septum.s1-peaks.narrowPeak.gz W139.heart.RV.s1-peaks.narrowPeak.gz W139.heart.apex.s1-peaks.narrowPeak.gz W146.heart.apex.s1-peaks.narrowPeak.gz W137.heart.apex.s1-peaks.narrowPeak.gz W136.heart.LV.s1-peaks.narrowPeak.gz W136.heart.apex.s1-peaks.narrowPeak.gz" > ${outBedList}
   done

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'zcat --version | head -1' 'cut --version | head -1' 'bedtools --version' 'sort --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
