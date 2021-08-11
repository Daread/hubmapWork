#! /bin/bash

# Get the FIMO motifs from JASPAR. Following the lead of the BBI at 
curl http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt -o fileOutputs/JASPAR_2020_CORE_Vertebrates_non-redundant_meme.txt

# Make a test file
cat fileOutputs/JASPAR_2020_CORE_Vertebrates_non-redundant_meme.txt | head -n 40 > fileOutputs/JASPAR_2020_CORE_Vert_nonRedun_testTwoMotifs.txt


