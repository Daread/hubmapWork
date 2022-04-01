#!/bin/bash

#
# Path to the Nextflow processing run configuration file.
#
WORKING_DIR="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup"
CONFIG_FILE="$WORKING_DIR/experiment.config" 

#
# Nextflow executable and pipeline script locations.
#
NEXTFLOW="/net/trapnell/vol1/home/readdf/nextflow"
NF_ANALYSIS="/net/trapnell/vol1/home/readdf/bin/bbi-sci/main.nf"

#
# Get the path to the analyze output directory from
# the configuration file and set the Nextflow work
# directory to be in the analyze output directory.
#
WORK_DIR="$WORKING_DIR/work"
#
# Use REPORT_FIL, TRACE_FIL, and TIMELINE_FIL so that
# Nextflow writes some informative processing information
# in the analyze output directory.
#

REPORT_FIL=$WORKING_DIR/analysis.report.html
TRACE_FIL=$WORKING_DIR/analysis.trace.tsv
TIMELINE_FIL=$WORKING_DIR/analysis.timeline.html
#
# Nextflow run parameters.
#
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"

date > ./run_start.txt

#
# Run Nextflow sci-analysis pipeline.
#
$NEXTFLOW run $NF_ANALYSIS $PARS

date > run_finish.txt

