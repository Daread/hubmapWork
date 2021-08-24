# source activate exprKerasWork

import numpy as np 
import argparse
import pandas as pd
import math
import pdb

# from tensorflow import keras
# # from keras.models import Model 
# # from keras.layers import *
# from tensorflow.keras import Model
# from tensorflow.keras import layers
# import tensorflow as tf
# import tensorflow.keras.backend as K

parser = argparse.ArgumentParser(description="Parse args for processing notes and choices")
parser.add_argument("--promUpstream", type=int, 
										default=1000, help='Bases upstream of TSS used')
parser.add_argument("--promDownstream", type=int,
										default=200, help="Bases downstream of TSS")
parser.add_argument("--coaccessCutoff", type=float,
										default=.05, help="Coaccessibility min to link to promoters")
parser.add_argument("--maxNdistalSites", type=int,
										default=5, help = "Max distal sites to link to a promoter")
parser.add_argument("--peakSize", type=int, dest="peakSize",
										default=600, help="Size of distal peaks")

parser.add_argument("--predictionTask", type=str, dest="predictionTask",
										default="_log2_CPM", help="Measurement to predict in RNA")

parser.add_argument("--featureSet", type=str, dest="featureSet",
										default="Promoter_and_Distal", help="Promoter_and_Distal or Promoter_Only")

parser.add_argument("--cnnStruct", type=str, dest="cnnStruct", help="Structure of CNN",
						default="21filt_pool100_restFilt")
args = parser.parse_args()

cellTypes = ["Adipocytes", "B_Cell", "Cardiomyocyte", "Endocardium", "Fibroblast", 
              "Lymphatic_Endothelium", "Macrophage", "Mast_Cell", "Neuronal", "T_Cell",
                "Vascular_Endothelium", "VSM_and_Pericyte"]


# Get the one hot DF
oneHotName = ("OneHot_PromotersAndDist_Max" + str(args.maxNdistalSites)
				+ "_Up" + str(args.promUpstream) +
								"_Down" + str(args.promDownstream)
									 + "_cicCut" + str(args.coaccessCutoff) +
										"peakSize" + str(args.peakSize))
# oneHotDF = pd.read_csv(("../2021_07_26_setup_ATAC_to_expr_data_nobackup/fileOutputs/" + oneHotName + ".csv"))
oneHotDF = pd.read_pickle(("../2021_07_26_setup_ATAC_to_expr_data_nobackup/fileOutputs/" + oneHotName + ".csv"))

def getRNAdf(args, cellTypes):
	rnaDir = "../2021_07_26_setup_ATAC_to_expr_data_nobackup/rdsOutputs/RNA_Data/" 
	rnaDF = pd.read_csv(rnaDir + "Pseudobulked_RNA_countsPerMillionNotDeduplicated.csv")
	# Drop 
	rnaDF = rnaDF[["id", "gene_short_name"] + [x + args.predictionTask for x in cellTypes] ]
	return(rnaDF)


rnaDF = getRNAdf(args, cellTypes)

rnaDF = rnaDF.rename(columns={'id': "GeneID"})


def getCombinedDF(oneHotDF, rnaDF, args):
	# pdb.set_trace()
	# Get matching rows (genes in both)
	rnaDF = rnaDF.loc[rnaDF['GeneID'].isin(oneHotDF['GeneID'])]
	oneHotDF = oneHotDF.loc[oneHotDF["GeneID"].isin(rnaDF['GeneID'])]
	# Combine 
	comboDF = pd.merge(oneHotDF, rnaDF, on="GeneID")
	# comboDF = comboDF.drop(["X", "chr", "orientation", "linkedSites"], axis=1)
	comboDF = comboDF.drop(["Unnamed: 0", "chr", "orientation", "linkedSites"], axis=1)
	return(comboDF)


combinedDF = getCombinedDF(oneHotDF, rnaDF, args)

def getSubsetOfCombined(combinedDF, args, setToUse):
	combinedSubset = combinedDF.loc[combinedDF.Model_Set == setToUse]
	combinedSubset = combinedSubset.drop(["Model_Set", "GeneID", "gene_short_name"], axis=1)
	return (combinedSubset)


combinedTrain = getSubsetOfCombined(combinedDF, args, "Train")
combinedTest = getSubsetOfCombined(combinedDF, args, "Test")
combinedVal = getSubsetOfCombined(combinedDF, args, "Validation")

yCols = [x + args.predictionTask for x in cellTypes]

trainX = combinedTrain.drop(yCols, axis=1)
trainY = combinedTrain[yCols]



# def getCNNmodel(inputX, inputY, args):
# 	if (args.cnnStruct == "21filt_pool100_restFilt_poolRest"):
# 		promoterBlocks = int((args.promUpstream + args.promDownstream) / 100.0)
# 		# promoterPoolDownsize = ceil(21.0 / promoterBlocks )
# 		# Get a promoter model
# 		promoterFilters = layers.Conv2D(80, kernel_size=(4,21), activation="relu", padding="same")
# 		promMaxPoolOne = layers.MaxPooling2D(pool_size=(1,100))
# 		filtersFromPromPoolOne = layers.Conv2D(80, kernel_size=(1,promoterBlocks), activation="relu")
# 		# Inputs
# 		# promInput = np.dstack(inputX["PromoterPos"].tolist())
# 		promInput = (inputX["PromoterPos"].tolist())

# 		# formattedInput = ([np.array(x, np.int32) for x in promInput])

# 		formattedInput = np.dstack(promInput)
# 		formattedInput = np.rollaxis(formattedInput, -1)
# 		formattedInput = np.expand_dims(formattedInput, axis=3)


# 		inputTensor = tf.convert_to_tensor(formattedInput)
# 		inputTensor = tf.cast(inputTensor, tf.float32)

# 		# promInput = np.dstack(inputX["PromoterPos"].tolist())
# 		promoterConvPart = promoterFilters(inputTensor)
# 		promoterConvPart = promMaxPoolOne(promoterConvPart)
# 		promoterConvPart = filtersFromPromPoolOne(promoterConvPart)

		# 
		

		# # Enhancer processing
		# ########################################################
		# enhancerBlocks = (args.peakSize) / 100.0
		# # enhancerPoolDownsize = ceil(21.0 / enhancerBlocks )
		# # Get a enhancer model
		# enhancerFilters = Conv2D(80, kernel_size=(4,21), activation="relu", pading="same")
		# enhancerMaxPool = layers.MaxPooling2D(pool_size=(1,100))
		# filtersFromEnhancPool = Conv2D(80, kernel_size=(1,enhancerBlocks), activation="relu")

		# # Now assemble 




# Get the model structure
# trainCNNmodel = getCNNmodel(trainX, trainY, args)


import keras
from keras.models import Model 
import keras.layers as kl
import tensorflow as tf
from keras.utils import plot_model



# Take in a pandas column cast to a list, hodling 2d numpy arrays.
#   Run formatting and return a tensorflow tensor, appropriate for feeding into the model
def oneHotListToTensor(inputList):
	# Manipulate into a numpy array with the right dimensions
	formattedInput = np.dstack(inputList)
	formattedInput = np.rollaxis(formattedInput, -1)
	formattedInput = np.expand_dims(formattedInput, axis=3)
	# Now format to a tensor
	inputTensor = tf.convert_to_tensor(formattedInput)
	inputTensor = tf.cast(inputTensor, tf.float32)
	return(inputTensor)


def getDistalInputList(inputX, args):
	# Set up an empty list
	distalList = [None] * args.maxNdistalSites
	# Names of the columns
	distalNames = ["Distal_" + str(x + 1) for x in range(args.maxNdistalSites)]
	# Get the inputs
	for eachInd, eachName in enumerate(distalNames):
		distalList[eachInd] = oneHotListToTensor(inputX[eachName].tolist())
	return(distalList)







inputX = trainX
inputY = trainY

##################################################


promoterBlocks = int((args.promUpstream + args.promDownstream) / 100.0)
# Get a promoter model
firstConvFilt = 21
padSize = int((firstConvFilt - 1) / 2.0)
############################################################
promoterPadding = kl.ZeroPadding2D(padding=(0,padSize))
promoterFilters = kl.Conv2D(80, kernel_size=(4,21), activation="relu")
promMaxPoolOne = kl.MaxPooling2D(pool_size=(1,100))
filtersFromPromPoolOne = kl.Conv2D(80, kernel_size=(1,promoterBlocks), activation="relu")
# Inputs
promoterInputTensor = oneHotListToTensor(inputX["PromoterPos"].tolist())
# Assemble the layers
promKerasInput = kl.Input(shape=(promoterInputTensor.shape[1], promoterInputTensor.shape[2],
								 promoterInputTensor.shape[3],) )
paddedInput = promoterPadding(promKerasInput)
promoterConvPart = promoterFilters(paddedInput)
promoterConvPart = promMaxPoolOne(promoterConvPart)
promoterConvPart = filtersFromPromPoolOne(promoterConvPart)
promoterConvPart = kl.Flatten()(promoterConvPart)
###########################################################################

# See if this is running a promoter + distal site model. If so, set up the distal site layters
if (args.featureSet == "Promoter_and_Distal"):
	# Get the inputs, organized into a list
	distalTensorInputList = getDistalInputList(inputX, args)  # [None] * len(siteNames)
	distalBlocks = int((args.peakSize) / 100.0)
	# Layers I'll use
	distalPadding = kl.ZeroPadding2D(padding=(0,padSize))
	distalFilters = kl.Conv2D(80, kernel_size=(4,21), activation="relu")
	distalMaxPoolOne = kl.MaxPooling2D(pool_size=(1,100))
	filtersFromDistPoolOne = kl.Conv2D(80, kernel_size=(1,distalBlocks), activation="relu")
	# Now, for each input (aka each of args.maxNdistalSites I used) set up a network
	siteNames = ["Distal_" + str(x + 1) for x in range(args.maxNdistalSites)]
	distalInputList = [None] * len(siteNames)
	distalNetworkList = [None] * len(siteNames)
	# Get formatted input
	for eachInd, eachName in enumerate(siteNames):
		# Assemble layers for this
		distalInputList[eachInd] = kl.Input(shape=(distalTensorInputList[eachInd].shape[1], 
												distalTensorInputList[eachInd].shape[2],
											 distalTensorInputList[eachInd].shape[3],) )
		distalNetworkList[eachInd] = distalPadding(distalInputList[eachInd])
		distalNetworkList[eachInd] = distalFilters(distalNetworkList[eachInd])
		distalNetworkList[eachInd] = distalMaxPoolOne(distalNetworkList[eachInd])
		distalNetworkList[eachInd] = filtersFromDistPoolOne(distalNetworkList[eachInd])
		distalNetworkList[eachInd] = kl.Flatten()(distalNetworkList[eachInd])
	# Now add all these distal site intermediates together
	pooledDistal = kl.Add()(distalNetworkList)
	flattenedConv = kl.merge.concatenate([promoterConvPart, pooledDistal])
else:
	flattenedConv = promoterConvPart
# Into dense layers
denseOutput = kl.Dense(100, activation="relu")(flattenedConv)
finalOutput = kl.Dense(1, activation="linear")(denseOutput)

if (args.featureSet == "Promoter_and_Distal"):
	fullModel = keras.models.Model( ([promKerasInput] + distalInputList),
									 finalOutput)
else:
	fullModel = keras.models.Model(promKerasInput, finalOutput)

print(fullModel.summary())



if (args.featureSet == "Promoter_and_Distal"):
	plot_model(fullModel, to_file="./plots/promAndDistalArchitectureTest.png")
else:
	plot_model(fullModel, to_file="./plots/promOnlyModelTest.png")




# Run fitting














