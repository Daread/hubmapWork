#!/bin/bash -ue
PROCESS_BLOCK='makeMotifMatrixProcess'
   SAMPLE_NAME="W146.heart.apex.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/load_python_env_reqs.sh
source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/python_env/bin/activate

outPeakTfMatrix="W146.heart.apex.s1-peak_motif_matrix.mtx.gz"

python /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/generate_motif_matrix.py 	--peak_motif_files W146.heart.apex.s1-gc_05-peak_calls.bb W146.heart.apex.s1-gc_13-peak_calls.bb W146.heart.apex.s1-gc_18-peak_calls.bb W146.heart.apex.s1-gc_11-peak_calls.bb W146.heart.apex.s1-gc_12-peak_calls.bb W146.heart.apex.s1-gc_06-peak_calls.bb W146.heart.apex.s1-gc_02-peak_calls.bb W146.heart.apex.s1-gc_16-peak_calls.bb W146.heart.apex.s1-gc_25-peak_calls.bb W146.heart.apex.s1-gc_08-peak_calls.bb W146.heart.apex.s1-gc_21-peak_calls.bb W146.heart.apex.s1-gc_15-peak_calls.bb W146.heart.apex.s1-gc_01-peak_calls.bb W146.heart.apex.s1-gc_09-peak_calls.bb W146.heart.apex.s1-gc_07-peak_calls.bb W146.heart.apex.s1-gc_22-peak_calls.bb W146.heart.apex.s1-gc_23-peak_calls.bb W146.heart.apex.s1-gc_24-peak_calls.bb W146.heart.apex.s1-gc_20-peak_calls.bb W146.heart.apex.s1-gc_04-peak_calls.bb W146.heart.apex.s1-gc_14-peak_calls.bb W146.heart.apex.s1-gc_17-peak_calls.bb W146.heart.apex.s1-gc_03-peak_calls.bb W146.heart.apex.s1-gc_10-peak_calls.bb W146.heart.apex.s1-gc_19-peak_calls.bb 	--fasta /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/data/genomeData/backupBBIHumanGenome_nobackup/Homo_sapiens.GRCh38.dna.toplevel.fa.finished 	--peaks W146.heart.apex.s1-merged_peaks.bed 	--motifs /net/trapnell/vol1/ajh24/common_data/sciatac_pipeline_data/common_files/motifs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.pfm 	--peak_tf_matrix ${outPeakTfMatrix}

   deactivate

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
