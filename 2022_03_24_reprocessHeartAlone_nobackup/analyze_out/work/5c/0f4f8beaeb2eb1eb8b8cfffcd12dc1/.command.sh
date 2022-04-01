#!/bin/bash -ue
PROCESS_BLOCK='callPeaksProcess'
   SAMPLE_NAME="W139.heart.septum.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

# MACS uses an underscore between the sample name and the file type
# so we need to rename output files.
outDir="call_peaks"
outMacsNarrowPeak="${outDir}/W139.heart.septum.s1_peaks.narrowPeak"
outNarrowPeak="W139.heart.septum.s1-peaks.narrowPeak.gz"
outMacsXls="${outDir}/W139.heart.septum.s1_peaks.xls"
outXls="W139.heart.septum.s1-peaks.xls"
outMacsSummits="${outDir}/W139.heart.septum.s1_summits.bed"
outSummits="W139.heart.septum.s1-summits.bed"

   # We used to add --shift -100 and --extsize, but the regions are now pre-shifted and extended
   # as output by other stages (ajh).
macs2 callpeak -t W139.heart.septum.s1-transposition_sites.bed.gz 		-f BED 		-g hs 		--nomodel 		--shift -100 		--extsize 200 		--keep-dup all 		--call-summits 		-n W139.heart.septum.s1 		--outdir ${outDir} 2> /dev/null

cat ${outMacsNarrowPeak} 		| sort -k1,1V -k2,2n -k3,3n 		| cut -f1-3 		| gzip > ${outNarrowPeak}

mv ${outMacsXls} ${outXls}
mv ${outMacsSummits} ${outSummits}

rm ${outMacsNarrowPeak}

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version' 'macs2 --version' 'sort --version | head -1' 'cut --version | head -1' 'gzip --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
