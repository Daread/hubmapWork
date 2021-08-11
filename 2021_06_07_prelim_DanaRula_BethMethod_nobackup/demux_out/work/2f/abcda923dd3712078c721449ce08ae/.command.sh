#!/bin/bash -ue
mkdir demux_dash
cp -R /net/trapnell/vol1/home/readdf/bin/bbi-dmux/bin/skeleton_dash/* demux_dash/
generate_html.R         "." --p7_rows "C" --p5_cols "1" --p7_wells "0" --p5_wells "0" --level "3" --project_name "2021_06_07_prelim_DanaRula_BethMethod_nobackup" --sample_sheet "good_sample_sheet.csv"
