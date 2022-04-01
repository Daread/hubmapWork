#!/bin/bash

# Load bedtools
module load bedtools/2.29.2

cd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/archr/results/NB12/

bedFileToUse=$1
# nb12Dir=/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/archr/results/NB12/


bedtools sort -i $bedFileToUse.bed > sortedTempBed.bed

bedtools genomecov -bg -i sortedTempBed.bed -g hg38.chrom.sizes  > $bedFileToUse.bedGraph







