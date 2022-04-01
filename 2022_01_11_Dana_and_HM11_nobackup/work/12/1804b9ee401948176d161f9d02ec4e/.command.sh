#!/bin/bash -ue
for f in HEK_cell_qc.csv
do
  awk 'BEGIN {FS=","}; $2>100{c++} END{print FILENAME, "100", c-1}' $f >> cell_counts.txt
  awk 'BEGIN {FS=","}; $2>500{c++} END{print FILENAME, "500", c-1}' $f >> cell_counts.txt
  awk 'BEGIN {FS=","}; $2>1000{c++} END{print FILENAME, "1000", c-1}' $f >> cell_counts.txt
done
