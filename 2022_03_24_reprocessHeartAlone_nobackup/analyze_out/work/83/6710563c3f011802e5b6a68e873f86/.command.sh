#!/bin/bash -ue
PROCESS_BLOCK='makeMergedPlotFilesProcess'
SAMPLE_NAME="na"
START_TIME=`date '+%Y%m%d:%H%M%S'`

mkdir -p /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/merged_plots
/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/merge_summary_plots.py -i /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/args.json -o merged.called_cells_summary.pdf
/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/merge_umap_plots.py -i /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_out/args.json -o merged.umap_plots.pdf

header='sample cell_threshold fraction_hs fraction_tss median_per_cell_frip median_per_cell_frit tss_enrichment sample_peaks_called total_merged_peaks total_reads fraction_reads_in_cells total_barcodes number_of_cells median_reads_per_cell min_reads_per_cell max_reads_per_cell median_duplication_rate median_fraction_molecules_observed median_total_fragments total_deduplicated_reads fraction_mitochondrial_reads [bloom_collision_rate]'
stats_file='merged.called_cells_summary.stats.csv'
header_wtabs=`echo ${header} | sed 's/ /	/g'`
rm -f ${stats_file}
echo "${header_wtabs}" > ${stats_file}
lfil=`ls *-called_cells_summary.stats.txt`
for fil in ${lfil}
do
  tail -n +2 ${fil} >> ${stats_file}
done

STOP_TIME=`date '+%Y%m%d:%H%M%S'`
/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version' 'sed --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
