#! /bin/bash


cd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup

# load bedtools
module load bedtools/2.29.2

# Make the mock bed file
# locusLine=chr9\t12555547\t10965244\n
# echo -e "chr9\t10965244\t12555547" > ./fileOutputs/gwasFollowUp/9p21_locus.bed
echo -e "chr9\t22072041\t22130390" > ./fileOutputs/gwasFollowUp/9p21_locus.bed
locusBed=./fileOutputs/gwasFollowUp/9p21_locus.bed

peaksBed=./fileOutputs/bedFilesForGWAS_0_0_0.015_20_600_2_0.1/0_0_0.015_20_600_2_0.1_masterGeneAndPeakFileSortedMerged.bed

head $locusBed

head $peaksBed


bedtools intersect -a $locusBed -b $peaksBed | head






