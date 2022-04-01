#! /bin/bash


Rscript runGlmnetClassification.R   --promoterUpstream $1 --promoterDownstream $2 --coaccessCutoff $3 --maxNdistalSites $4 --peakSize $5 --featureSelection $6 --predictionTask $7 --pValFIMOcutoff $8

