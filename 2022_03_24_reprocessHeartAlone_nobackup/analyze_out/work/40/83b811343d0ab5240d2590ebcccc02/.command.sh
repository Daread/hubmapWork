#!/bin/bash -ue
PROCESS_BLOCK='makeReducedDimensionMatrixProcess'
   SAMPLE_NAME="W142.heart.LV.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

   inPeakMatrix="W142.heart.LV.s1-peak_matrix.mtx.gz"
   inSampleName="W142.heart.LV.s1"

	outScrubletHistFile="W142.heart.LV.s1-scrublet_hist.png"
   outScrubletTableFile="W142.heart.LV.s1-scrublet_table.csv"
outLsiCoordsFile="W142.heart.LV.s1-lsi_coords.txt"
outUmapCoordsFile="W142.heart.LV.s1-umap_coords.txt"
outUmapPlotFile="W142.heart.LV.s1-umap_plot"
outMonocle3CdsFile="W142.heart.LV.s1-monocle3_cds.rds"
   outBlackListRegionsFile="W142.heart.LV.s1-blacklist_regions_file.log"
   outReduceDimensionsLogFile="W142.heart.LV.s1-reduce_dimensions.log"
   umi_cutoff=100
   frip_cutoff=0.1
   frit_cutoff=0.05
   num_lsi_dimensions=50
   cluster_resolution=1.0e-3
doublet_predict_top_ntile=0.1

   doublet_predict=""
   if [ "false" = "true" ]
   then
       doublet_predict=" --doublet_predict "
   fi

   black_list_file=""
   echo "blacklist_regions: /net/bbi/vol1/data/genomes_stage/human/human_encode/ENCFF356LFX.bed.sorted"

   if [[ "true" = "true" && "/net/bbi/vol1/data/genomes_stage/human/human_encode/ENCFF356LFX.bed.sorted" != "" ]]
   then
       black_list_file=" --black_list_file /net/bbi/vol1/data/genomes_stage/human/human_encode/ENCFF356LFX.bed.sorted "
       echo "blacklist_region file: /net/bbi/vol1/data/genomes_stage/human/human_encode/ENCFF356LFX.bed.sorted" > ${outBlackListRegionsFile}
   else
     echo "blacklist_region file: none" > ${outBlackListRegionsFile}
   fi

   if [ "false" = "true" ]
   then
       /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/run_scrublet.py --sample_name=${inSampleName} --mat_file=${inPeakMatrix} --umi_cutoff=$umi_cutoff
   fi

Rscript /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/reduce_dimensions.R     --sample_name ${inSampleName} 	--mat_file ${inPeakMatrix}     --count_file W142.heart.LV.s1-count_report.txt     --umi_cutoff ${umi_cutoff}     --frip_cutoff ${frip_cutoff}     --frit_cutoff ${frit_cutoff}     --doublet_predict_top_ntile ${doublet_predict_top_ntile}     --num_lsi_dimensions ${num_lsi_dimensions}     --cluster_resolution ${cluster_resolution}     --combined_read_count W142.heart.LV.s1-combined.duplicate_report.txt     --cds_file ${outMonocle3CdsFile}     --lsi_coords_file ${outLsiCoordsFile}     --umap_coords_file ${outUmapCoordsFile}     --umap_plot_file ${outUmapPlotFile} ${doublet_predict} ${black_list_file}

   if [ ! -e "${outScrubletHistFile}" ]
   then
     touch "${outScrubletHistFile}"
   fi

   if [ ! -e "${outScrubletTableFile}" ]
   then
     touch "${outScrubletTableFile}"
   fi

   if [ ! -e "${outLsiCoordsFile}" ]
   then
     touch "${outLsiCoordsFile}"
   fi

   if [ ! -e "${outUmapCoordsFile}" ]
   then
     touch "${outUmapCoordsFile}"
   fi

   if [ ! -e "${outUmapPlotFile}.png" ]
   then
     touch "${outUmapPlotFile}.null"
   fi

   if [ ! -e "${outMonocle3CdsFile}" ]
   then
     touch "${outMonocle3CdsFile}"
   fi

   if [ ! -e "${outBlackListRegionsFile}" ]
   then
     touch "${outBlackListRegionsFile}"
   fi

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version' 'R --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -f ${outReduceDimensionsLogFile} ${outBlackListRegionsFile}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
