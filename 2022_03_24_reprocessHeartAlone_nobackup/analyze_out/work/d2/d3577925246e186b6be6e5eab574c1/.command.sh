#!/bin/bash -ue
PROCESS_BLOCK='callCellsProcess'
   SAMPLE_NAME="W139.heart.septum.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

   source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/load_python_env_reqs.sh
   source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/python_env/bin/activate

outCalledCellsCounts="W139.heart.septum.s1-called_cells.txt"
outCellWhiteList="W139.heart.septum.s1-called_cells_whitelist.txt"
outCallCellsStats="W139.heart.septum.s1-called_cells_stats.json"

   python /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/call_cells.py ${outCalledCellsCounts}                                        ${outCellWhiteList}                                        --fit_metadata ${outCallCellsStats}                                        --count_report W139.heart.septum.s1-count_report.txt 

   mv ${outCellWhiteList} cell_whitelist.txt.tmp
   sort cell_whitelist.txt.tmp > ${outCellWhiteList}
   rm cell_whitelist.txt.tmp

   deactivate

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version' 'sort --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
