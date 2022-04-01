#!/bin/bash -ue
PROCESS_BLOCK='getUniqueFragmentsProcess'
   SAMPLE_NAME="W139.heart.septum.s1"
   START_TIME=`date '+%Y%m%d:%H%M%S'`

source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/load_python_env_reqs.sh
source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/python_env/bin/activate

outFragments="W139.heart.septum.s1-fragments.txt"
outTranspositionSites="W139.heart.septum.s1-transposition_sites.bed"
outInsertSizes="W139.heart.septum.s1-insert_sizes.txt"
outDuplicateReport="W139.heart.septum.s1-duplicate_report.txt"

python /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/get_unique_fragments.py 		W139.heart.septum.s1-merged.bam 		--fragments ${outFragments} 		--transposition_sites_bed ${outTranspositionSites} 		--duplicate_read_counts ${outDuplicateReport} 		--insert_sizes ${outInsertSizes}

# Index BAM file / bgzip tabix index for fragments file and transposition_sites BED
bgzip -f ${outFragments}
tabix -p bed ${outFragments}.gz

bgzip -f ${outTranspositionSites}
tabix -p bed ${outTranspositionSites}.gz

   # Make bedGraph file for genome browser.
   outBedGraph="W139.heart.septum.s1-read_alignments.bedgraph"
   outBigWig="W139.heart.septum.s1-read_alignments.bw"

   bedtools genomecov -ibam W139.heart.septum.s1-merged.bam -bga | awk 'BEGIN{OFS="	"}{$1="chr" $1}1' > ${outBedGraph}
   # bedGraphToBigWig ${outBedGraph} W139.heart.septum.s1-human.chromosome_sizes.sorted.txt ${outBigWig}

   deactivate

   STOP_TIME=`date '+%Y%m%d:%H%M%S'`
   /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-analyze/src/pipeline_logger.py     -r `cat /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/tmp/nextflow_run_name.txt`     -n ${SAMPLE_NAME}     -p ${PROCESS_BLOCK}     -v 'python3 --version' 'tabix 2>&1 > /dev/null | head -3 | tail -1' 'bedtools --version | head -1' 'awk --version | head -1'     -s ${START_TIME}     -e ${STOP_TIME}     -d /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/analyze_log_dir
