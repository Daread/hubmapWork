#!/bin/bash -ue
generate_dash_data.R all_sample_stats.csv /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_06_07_prelim_DanaRula_BethMethod_nobackup cell_counts.txt all_collision.txt false

mkdir exp_dash
cp -R /net/trapnell/vol1/home/readdf/bin/bbi-sci/bin/skeleton_dash/* exp_dash/
mv *.png exp_dash/img/

mv *.js exp_dash/js/
