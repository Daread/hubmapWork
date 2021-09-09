
import os
import argparse
import easygrid
import pdb # DFR add 8-31-21

testOutputPath = "/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/testEasyGridDir/"

pipeline = easygrid.JobManager(testOutputPath)

testOutFile = testOutputPath + "helloWorld.txt"
testCommand = 'echo "Hello World" ' 

pipeline.add(testCommand, name="TEST_PIPELINE", walltime="10:00:00", memory='10G', outputs=[testOutFile])


pipeline.run(infer_dependencies=False, dry=False)


