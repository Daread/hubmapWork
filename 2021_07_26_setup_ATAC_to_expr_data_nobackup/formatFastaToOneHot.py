
import numpy as np 
import argparse
import pandas as pd
import math
import pdb


def oneHotFromSeq(inputStr, seqLength):
	 # See if this is empty
	# If so, return a zero'd vector
	# pdb.set_trace()
	if not ((isinstance(inputStr, str)) ):
		return (np.zeros((4, seqLength)))
	# Set up a numpy array for this
	oneHotVec = np.zeros((4,len(inputStr)))
	# 
	# Loop and fill in
	for eachInd, eachChar in enumerate(inputStr):
		# A
		if (eachChar == "A"):
			oneHotVec[0, eachInd] = 1
		elif(eachChar == "T"):
			oneHotVec[1, eachInd] = 1
		elif(eachChar == "G"):
			oneHotVec[2, eachInd] = 1
		elif(eachChar == "C"):
			oneHotVec[3, eachInd] = 1
		#
	# Return
	return (oneHotVec)
# End oneHotFromSeq

def getOneHotDF(inputDF, args):
	# Get the columns that need to be processed
	colToProcess = ["PromoterPos"]
	for eachDistalNum in range(1, (1 + args.maxNdistalSites)):
		colToProcess.append("Distal_" + str(eachDistalNum))
#
	# For each of these columns, get the one-hot version
	for eachCol in colToProcess:
		# Send a peak size to fill for NAs
		if (eachCol == "PromoterPos"):
			inputDF[eachCol] = (
				[oneHotFromSeq(x, (args.promUpstream + args.promDownstream)) for x in inputDF[eachCol]]
				)
		else:
			inputDF[eachCol] = (
				[oneHotFromSeq(x, args.peakSize) for x in inputDF[eachCol]]
				)
#
	# Return
	return(inputDF)
# End getOneHotDF


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

args = parser.parse_args()

# The sequences are output as a csv representing an R dataframe. Read this in,
#   accounting for optional arguments that describe processing choices
dfName = ("Gene_Prom_Plus_Distal_WithSequence_Sites_Max" + str(args.maxNdistalSites)
				+ "_Upstream" + str(args.promUpstream) +
								"_Downstream" + str(args.promDownstream)
									 + "_cicCuf" + str(args.coaccessCutoff) +
										"peakSize" + str(args.peakSize))

myDF = pd.read_csv("./rdsOutputs/" + dfName + ".csv")

oneHotDF = getOneHotDF(myDF, args)



myDF = pd.read_csv("./rdsOutputs/" + dfName + ".csv")


# Save the output:
oneHotName = ("OneHot_PromotersAndDist_Max" + str(args.maxNdistalSites)
				+ "_Up" + str(args.promUpstream) +
								"_Down" + str(args.promDownstream)
									 + "_cicCut" + str(args.coaccessCutoff) +
										"peakSize" + str(args.peakSize))


oneHotDF.to_csv(("./fileOutputs/" + oneHotName + ".csv"))




































