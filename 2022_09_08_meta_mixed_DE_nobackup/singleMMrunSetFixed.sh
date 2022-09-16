#! /bin/bash

# Commands to set up R 4.1.2 to use 
module unload R/4.0.0
module unload pcre2/10.35
module load pcre2/10.39
module load R/4.1.2



Rscript /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_09_08_meta_mixed_DE_nobackup/runMixedModelOnHeartRNA.R --cellType $1 --fixedEffects $2 --nSplits $3 --split $4

