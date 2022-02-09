#! /bin/bash

Rscript /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_11_08_motifInPeaksRegression_nobackup/runMotifCelltypeRegression.R --cellType $1 --fixedEffects $2 --ATACprocNote FRIP=0.2_FRIT=0.05UMI=1000DL=0.5_useMNN_Peak_CDS_50k_20 --groupColumn Assigned_Cell_Type --cdsPath /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_07_15_Greg_ATAC_Code_nobackup/archr/results/NB9/All_Cells/

