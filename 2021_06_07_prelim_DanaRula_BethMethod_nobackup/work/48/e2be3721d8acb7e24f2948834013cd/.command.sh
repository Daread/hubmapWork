#!/bin/bash -ue
sed -s 1d Pancreas_sample_stats.csv Liver_sample_stats.csv Spleen_sample_stats.csv Lung_sample_stats.csv > all_sample_stats.csv
