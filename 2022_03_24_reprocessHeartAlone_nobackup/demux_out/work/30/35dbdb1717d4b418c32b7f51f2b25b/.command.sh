#!/bin/bash -ue
mkdir fastqc_sample
fastqc *.fastq.gz -t 8 -o fastqc_sample
