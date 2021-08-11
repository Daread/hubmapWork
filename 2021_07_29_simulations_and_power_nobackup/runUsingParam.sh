#! /bin/bash

# Rscript fileHere indiv cells ngenes indivSD cellSD
Rscript /net/trapnell/vol1/home/readdf/trapLabDir/fibrosis/results/2020_12_29_SymSim_Power_Work/runSimulations.R --indiv $1 --cells $2 --ngenes $3 --indSD $4 --cellSD $5
