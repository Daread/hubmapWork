#!/bin/bash -ue
mkdir demux_out
make_sample_fastqs.py --run_directory /net/shendure/vol9/seq/NEXTSEQ/220110_NS500488_1276_AH23G5BGXL         --read1 <(zcat Undetermined_S0_L002_R1_001.fastq.gz) --read2 <(zcat Undetermined_S0_L002_R2_001.fastq.gz)         --file_name Undetermined_S0_L002_R1_001.fastq.gz --sample_layout good_sample_sheet.csv         --p5_cols_used 4 --p7_rows_used D         --p5_wells_used 0 --p7_wells_used 0         --rt_barcode_file /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/misc/RT_2level_indices_plate_format         --multi_exp "0"         --output_dir ./demux_out --level 2
pigz -p 8 demux_out/*.fastq
