#! /bin/bash

# Run late Jan, 2021
echo "Starting batch sub"

for CELLCOUNT in 300 1000 3000 10000; do
	for INDIVCOUNT in 3 5 10 15; do
		qsub -P trapnelllab runUsingParam.sh $INDIVCOUNT $CELLCOUNT 80 .1 .02
	done
done

# Now with very low noise
for CELLCOUNT in 300 1000 3000 10000; do
	for INDIVCOUNT in 3 5 10 15; do
		qsub -P trapnelllab runUsingParam.sh $INDIVCOUNT $CELLCOUNT 80 .01 .002
	done
done

echo "Done with batch sub"


# Feb 3rd, 2021
qsub -P trapnelllab runUsingParam.sh 6 500 100 .35 .1
qsub -P trapnelllab runUsingParam.sh 6 500 100 .001 .0001

echo "Submitted"