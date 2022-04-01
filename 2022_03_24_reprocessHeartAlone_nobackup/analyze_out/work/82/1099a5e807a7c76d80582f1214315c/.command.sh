#!/bin/bash -ue
PROCESS_BLOCK='makeWindowMatrixProcess'
   SAMPLE_NAME="W146.heart.LV.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

   source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/load_python_env_reqs.sh
   source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/python_env/bin/activate

outWindowMatrix="W146.heart.LV.s1-window_matrix.mtx.gz"

   python /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/generate_sparse_matrix.py     --transposition_sites_intersect <(bedtools intersect -sorted -g W146.heart.LV.s1-human.chromosome_sizes.sorted.txt -a W146.heart.LV.s1-genomic_windows.bed -b W146.heart.LV.s1-transposition_sites.bed.gz -wa -wb)     --intervals W146.heart.LV.s1-genomic_windows.bed     --cell_whitelist W146.heart.LV.s1-called_cells_whitelist.txt     --matrix_output ${outWindowMatrix}

   deactivate

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version' 'bedtools --version'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
