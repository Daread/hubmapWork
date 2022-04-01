#!/bin/bash -ue
generate_dash_data.R all_sample_stats.csv /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup cell_counts.txt all_collision.txt false

mkdir exp_dash
cp -R /net/trapnell/vol1/home/readdf/bin/bbi-sci/bin/skeleton_dash/* exp_dash/
mv *.png exp_dash/img/

mv *.js exp_dash/js/
