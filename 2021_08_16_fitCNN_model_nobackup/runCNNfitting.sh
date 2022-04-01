#! /bin/bash



# 8-25-21: Fixed failure to reverse complement sequence
python fitCNN_model.py --featureSet Promoter_Only


python fitCNN_model.py --featureSet Promoter_Only --cnnStruct 21filt_pool100_thenDense



python fitCNN_model.py --featureSet Promoter_and_Distal --cnnStruct 21filt_pool100_thenDense --convNum 500


