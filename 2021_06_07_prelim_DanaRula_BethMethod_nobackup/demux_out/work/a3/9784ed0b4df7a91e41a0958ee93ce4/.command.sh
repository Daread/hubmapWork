#!/bin/bash -ue
mkdir demux_out
make_sample_fastqs.py --run_directory /net/shendure/vol9/seq/NEXTSEQ/210602_NB552332_0143_AHK7GWBGXH         --read1 <(zcat Undetermined_S0_L001_R1_001.fastq.gz) --read2 <(zcat Undetermined_S0_L001_R2_001.fastq.gz)         --file_name Undetermined_S0_L001_R1_001.fastq.gz --sample_layout good_sample_sheet.csv         --p5_cols_used 1 --p7_rows_used C         --p5_wells_used 0 --p7_wells_used 0         --rt_barcode_file /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_06_07_prelim_DanaRula_BethMethod_nobackup/dfrBBItoTrapRTplates.txt         --multi_exp "0"         --output_dir ./demux_out --level 3
pigz -p 8 demux_out/*.fastq
