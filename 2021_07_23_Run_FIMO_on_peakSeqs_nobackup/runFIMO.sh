#! /bin/bash

# FIMO needs a few things. It needs motifs in MEME format (already pulled from JASPAR)
MOTIF_FILE=./fileOutputs/JASPAR_2020_CORE_Vertebrates_non-redundant_meme.txt

inputSeqFile=$1

# cds_p_allHeartATAC_peakSeqsAsFasta.fa
# It also needs sequence to analyze. For this, I'm pulling a file (set here) of sequence within/centered at ATAC peaks
PEAK_SEQ_FILE=../2021_07_26_setup_ATAC_to_expr_data_nobackup/fileOutputs/$inputSeqFile #../2021_07_22_Get_sequenceFromATAC_CDS_nobackup/rdsOutput/$inputSeqFile

echo $PEAK_SEQ_FILE

# Load MEME to get FIMO
module load python/2.7.13
module load meme/5.0.5

# Run command
fimo --max-stored-scores 200000 --oc $inputSeqFile $MOTIF_FILE $PEAK_SEQ_FILE 







# Test


# TEST_MOTIF_FILE=./fileOutputs/JASPAR_2020_CORE_Vert_nonRedun_testTwoMotifs.txt
# fimo  --o testFIMO_output $TEST_MOTIF_FILE $PEAK_SEQ_FILE


