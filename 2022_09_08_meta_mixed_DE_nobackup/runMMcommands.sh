

splitNum=10
# for eachType in Endothelium Fibroblast Lymphocyte Myeloid Perivascular Ventricular_Cardiomyocytes Neuron Adipocytes
# for eachType in Adipocytes
for eachType in Endothelium Fibroblast Lymphocyte Myeloid Perivascular Ventricular_Cardiomyocytes Neuron Adipocytes
do
	# Loop through and submit a job for each sub-split
	for ((eachSplit=1; eachSplit<=$splitNum; eachSplit++))
	do
		qsub -P trapnelllab -l mfree=80G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_09_08_meta_mixed_DE_nobackup/logs/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_09_08_meta_mixed_DE_nobackup/singleMMrunSetFixed.sh $eachType Anatomical_Site,Age,Sex,DataSource,log10_umi $splitNum $eachSplit
	done
done



