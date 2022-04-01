#!/bin/bash -ue
mkdir fastqc_lanes
fastqc *.fastq.gz -t 8 -o fastqc_lanes
