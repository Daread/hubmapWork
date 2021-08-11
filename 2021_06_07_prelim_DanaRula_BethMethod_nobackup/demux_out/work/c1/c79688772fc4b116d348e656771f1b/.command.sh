#!/bin/bash -ue
min_threads=$(((16/2)<4 ? (16/2):4))

bcl2fastq -R /net/shendure/vol9/seq/NEXTSEQ/210602_NB552332_0143_AHK7GWBGXH --output-dir ./lane_fastqs         --sample-sheet SampleSheet.csv         --loading-threads $min_threads         --processing-threads 16          --writing-threads $min_threads         --barcode-mismatches 1         --ignore-missing-positions         --ignore-missing-controls         --ignore-missing-filter         --ignore-missing-bcls         --minimum-trimmed-read-length 15         --mask-short-adapter-reads 15
