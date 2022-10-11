#! /bin/bash


for eachType in Adipocytes Endothelium Fibroblast Lymphocyte Neuron Perivascular Ventricular_Cardiomyocytes Mast_Cell Lymphatic_Endothelium Myeloid
do
	qsub -l mfree=80G -l centos=7 /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_09_02_DE_Testing_nobackup/runSingleBulkDE.sh eachType
done



