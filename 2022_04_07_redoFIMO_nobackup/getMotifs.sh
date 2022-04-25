#! /bin/bash

# Get the FIMO motifs from JASPAR. Following the lead of the BBI at 
curl https://jaspar.genereg.net/download/data/2018/CORE.backup-20220216/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt -o fileOutputs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt




# Make a test file
cat fileOutputs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt | head -n 40 > fileOutputs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_meme.txt





# Copy over the motifs from the BBI pipeline
cp /net/trapnell/vol1/ajh24/common_data/sciatac_pipeline_data/common_files/motifs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.pfm fileOutputs/BBI_Pipeline_JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.pfm

