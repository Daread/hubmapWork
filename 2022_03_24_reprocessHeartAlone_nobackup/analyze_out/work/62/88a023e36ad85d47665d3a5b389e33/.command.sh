#!/bin/bash -ue
PROCESS_BLOCK='callMotifsProcess'
   SAMPLE_NAME="W139.heart.septum.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/load_python_env_reqs.sh
source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/python_env/bin/activate

gc_bin_padded=`echo 6 | awk '{printf("%02d",$1+1)}'`
outGcBinned="W139.heart.septum.s1-gc_${gc_bin_padded}-peak_calls.bb"

python /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/call_peak_motifs.py /net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished W139.heart.septum.s1-merged_peaks.bed /net/trapnell/vol1/ajh24/common_data/sciatac_pipeline_data/common_files/motifs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.pfm ${outGcBinned} --gc_bin 6 --pwm_threshold 1E-7

   deactivate

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'awk --version | head -1' 'python3 --version'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
