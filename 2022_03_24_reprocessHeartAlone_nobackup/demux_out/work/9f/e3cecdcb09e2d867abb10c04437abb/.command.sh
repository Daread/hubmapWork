#!/bin/bash -ue
mkdir demux_dash
cp -R /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-demux/src/skeleton_dash/* demux_dash

/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-demux/src/make_dashboard_data.R  --input_folder="."                                      --samplesheet=/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/dfr_samplesheet.json                                      --project_name=2022_03_24_reprocessHeartAlone_nobackup                                      --image_output_folder="demux_dash/img"

/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-demux/src/make_run_data.py --input_folder="."                                --input_file="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup/demux_out/args.json"                                --output_file="demux_dash/js/run_data.js"
