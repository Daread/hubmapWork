
library("optparse")
source("/net/trapnell/vol1/home/readdf/trapLabDir/fibrosis/results/2020_12_29_SymSim_Power_Work/simFuncsDavidMod.R")


print("Starting Now")

option_list = list(
	make_option(c("-i", "--indiv"), action="store", default=3, type="integer",
				help="number of individuals per group"),
	make_option(c("-c", "--cells"), action="store", default=300, type="integer",
				help="number of cells per individual")
	) # end opttion_list

opt = parse_args(OptionParser(option_list=option_list))


print(str(opt[["indiv"]]))

print(str(opt))
print(opt)


print("Now using a function")


optList = processParamArgs()

print(str(optList))










print("All Done")










