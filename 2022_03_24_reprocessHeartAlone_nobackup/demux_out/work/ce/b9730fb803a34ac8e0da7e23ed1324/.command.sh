#!/bin/bash -ue
source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-demux/load_pypy_env_reqs.sh
PS1=${PS1:-}
source /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-demux/src/pypy_env/bin/activate

pypy /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-demux/src/barcode_correct_sciatac.py                        --samplesheet /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/dfr_samplesheet.json                        -1 <(zcat Undetermined_S0_L003_R1_001.fastq.gz)                        -2 <(zcat Undetermined_S0_L003_R2_001.fastq.gz)                        --filename Undetermined_S0_L003_R1_001.fastq.gz                        --out_dir .                        --stats_out 1                        --num_pigz_threads 6                        --write_buffer_blocks 8192                         --wells_384 --well_ids                        
deactivate
