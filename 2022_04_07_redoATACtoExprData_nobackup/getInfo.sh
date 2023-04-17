
for maxSite in 5 10 20
do
	for cutoff in 0.005 0.015 0.035 0.05 0.1
	do
		for peak in 600 1000
		do
			qsub -P trapnelllab -l mfree=32G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/singleInfo.sh 2000 1000 $cutoff $maxSite $peak LyonV1
		done
	done
done




for maxSite in 5
do
	for cutoff in 0.005
	do
		for peak in 600 1000
		do
			qsub -P trapnelllab -l mfree=32G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/singleInfo.sh 2000 1000 $cutoff $maxSite $peak LyonV1
		done
	done
done





for maxSite in 5
do
	for cutoff in 0.005
	do
		for peak in 600 1000
		do
			echo "qsub -P trapnelllab -l mfree=32G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/singleInfo.sh 2000 1000 $cutoff $maxSite $peak LyonV1"
		done
	done
done





for maxSite in 5
do
	for cutoff in 0.005
	do
		for peak in 600 1000
		do
			echo "qsub -P trapnelllab -l mfree=32G -wd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/ /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_04_07_redoATACtoExprData_nobackup/singleInfo.sh 2000 1000 $cutoff $maxSite $peak LyonV1"
		done
	done
done
