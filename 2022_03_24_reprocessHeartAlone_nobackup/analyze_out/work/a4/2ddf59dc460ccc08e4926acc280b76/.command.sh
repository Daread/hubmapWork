#!/bin/bash -ue
PROCESS_BLOCK='getUniqueFragmentsMitoProcess'
   SAMPLE_NAME="W136.heart.apex.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/load_python_env_reqs.sh
source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/python_env/bin/activate

outMitoFragments="W136.heart.apex.s1-mito.fragments.txt"
outMitoTranspositionSites="W136.heart.apex.s1-mito.transposition_sites.bed"
outMitoInsertSizes="W136.heart.apex.s1-mito.insert_sizes.txt"
outMitoDuplicateReport="W136.heart.apex.s1-mito.duplicate_report.txt"

python /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/get_unique_fragments.py 		W136.heart.apex.s1-merged.mito.bam 		--fragments ${outMitoFragments} 		--transposition_sites_bed ${outMitoTranspositionSites} 		--duplicate_read_counts ${outMitoDuplicateReport} 		--insert_sizes ${outMitoInsertSizes}

# Index BAM file / bgzip tabix index for fragments file and transposition_sites BED
bgzip -f ${outMitoFragments}
tabix -p bed ${outMitoFragments}.gz

bgzip -f ${outMitoTranspositionSites}
tabix -p bed ${outMitoTranspositionSites}.gz

   deactivate

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version' 'tabix 2>&1 > /dev/null | head -3 | tail -1' 'bedtools --version | head -1' 'awk --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
