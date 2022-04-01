#!/bin/bash -ue
mkdir demux_dash
cp -R /net/trapnell/vol1/home/readdf/bin/bbi-dmux/bin/skeleton_dash/* demux_dash/
generate_html.R         "." --p7_rows "D" --p5_cols "4" --p7_wells "0" --p5_wells "0" --level "2" --project_name "2022_01_11_Dana_and_HM11_nobackup" --sample_sheet "good_sample_sheet.csv"
