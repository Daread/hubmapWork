#!/bin/bash -ue
PROCESS_BLOCK='runAlignProcess'
   SAMPLE_NAME="W135.heart.LV.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

   bowtieStderr="W135.heart.LV.s1-RUN001_L002.bowtie.stderr"
   samtoolsViewStderr="W135.heart.LV.s1-RUN001_L002.samtools_view.stderr"
   samtoolsSortStderr="W135.heart.LV.s1-RUN001_L002.samtools_sort.stderr"
outBam="W135.heart.LV.s1-RUN001_L002.bam"
   outMito="W135.heart.LV.s1-RUN001_L002.mito"

   if [ "true" == "true" ]
   then
     bowtie2 -3 1           -X 2000           -p 6           -x /net/bbi/vol1/data/genomes_stage/human/human_atac/human           -1 /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/demux_out/W135.heart.LV.s1/fastqs_trim/W135.heart.LV.s1-RUN001_L002_R1.trimmed.fastq.gz           -2 /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/demux_out/W135.heart.LV.s1/fastqs_trim/W135.heart.LV.s1-RUN001_L002_R2.trimmed.fastq.gz --seed 123456789 2> ${bowtieStderr}           | samtools view -h -L /net/bbi/vol1/data/genomes_stage/human/human_atac/whitelist_with_mt_regions.bed -f3 -F12 -q10 - 2> ${samtoolsViewStderr}           | /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/divert_mito_alignments.py -o ${outMito}.sam           | samtools sort -T ${outBam}.sorttemp --threads 4 -o ${outBam} - 2> ${samtoolsSortStderr}

     samtools sort -T ${outMito}.sorttemp --threads 4 ${outMito}.sam -o ${outMito}.bam
     rm ${outMito}.sam
   else
     bowtie2 -3 1           -X 2000           -p 6           -x /net/bbi/vol1/data/genomes_stage/human/human_atac/human           -1 /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/demux_out/W135.heart.LV.s1/fastqs_trim/W135.heart.LV.s1-RUN001_L002_R1.trimmed.fastq.gz           -2 /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/demux_out/W135.heart.LV.s1/fastqs_trim/W135.heart.LV.s1-RUN001_L002_R2.trimmed.fastq.gz --seed 123456789 2> ${bowtieStderr}           | samtools view -h -L /net/bbi/vol1/data/genomes_stage/human/human_atac/whitelist_regions.bed -f3 -F12 -q10 - 2> ${samtoolsViewStderr}           | samtools sort -T ${outBam}.sorttemp --threads 4 -o ${outBam} - 2> ${samtoolsSortStderr}

     samtools view -H -o ${outMito}.bam ${outBam}
   fi

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   MITO_READS="{\"mitochondrial_reads\": true}"
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'bowtie2 --version | head -1' 'samtools --version | head -2' 'bedtools --version' 'awk --version | head -1' 'sort --version | head -1' 'uniq --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -f ${bowtieStderr}     -j "${MITO_READS}"     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir