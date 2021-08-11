#!/bin/bash -ue
for f in Pancreas_cell_qc.csv Spleen_cell_qc.csv Liver_cell_qc.csv Lung_cell_qc.csv
do
  awk 'BEGIN {FS=","}; $2>100{c++} END{print FILENAME, "100", c-1}' $f >> cell_counts.txt
  awk 'BEGIN {FS=","}; $2>500{c++} END{print FILENAME, "500", c-1}' $f >> cell_counts.txt
  awk 'BEGIN {FS=","}; $2>1000{c++} END{print FILENAME, "1000", c-1}' $f >> cell_counts.txt
done
