#!/bin/bash -ue
PROCESS_BLOCK='experimentDashboardProcess'
   SAMPLE_NAME="all"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

mkdir -p /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_dash/js
/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/make_run_data.py -i /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/args.json -o run_data.js

mkdir -p /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_dash/js /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_dash/img
cp /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/skeleton_dash/img/*.png /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_dash/img
cp /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/skeleton_dash/js/* /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_dash/js
cp -r /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/skeleton_dash/style /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_dash
cp /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/skeleton_dash/exp_dash.html /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_dash

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
