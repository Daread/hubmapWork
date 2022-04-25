#!/bin/bash -ue
PROCESS_BLOCK='makePromoterMatrixProcess'
   SAMPLE_NAME="W146.heart.LV.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

   source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/load_python_env_reqs.sh
   source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/python_env/bin/activate

outPromoterMatrix="W146.heart.LV.s1-promoter_matrix.mtx.gz"

   python /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/generate_sparse_matrix.py     --transposition_sites_intersect <(bedtools intersect -sorted -g W146.heart.LV.s1-human.chromosome_sizes.sorted.txt -a W146.heart.LV.s1-gene_regions.bed.gz -b W146.heart.LV.s1-transposition_sites.bed.gz -wa -wb)     --intervals W146.heart.LV.s1-gene_regions.bed.gz     --cell_whitelist W146.heart.LV.s1-called_cells_whitelist.txt     --matrix_output ${outPromoterMatrix}

   mv W146.heart.LV.s1-promoter_matrix.rows.txt promoter_matrix_rows.txt.no_metadata
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/add_gene_metadata.py --in_promoter_matrix_row_name_file promoter_matrix_rows.txt.no_metadata                                        --gene_metadata_file /net/bbi/vol1/data/genomes_stage/human/human_atac/gene_bodies_gene_map.txt                                        --out_promoter_matrix_row_name_file W146.heart.LV.s1-promoter_matrix.rows.txt

   deactivate

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
