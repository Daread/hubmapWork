#!/bin/bash
#
# Path to the Nextflow processing run configuration file.
#
WORKING_DIR="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_06_07_prelim_DanaRula_BethMethod_nobackup"
CONFIG_FILE="$WORKING_DIR/experiment.config" 
#

# Nextflow executable and pipeline script locations.
# 
# Note: the bbi-demux repo lives within the nextflow dir and in you bin (if you cloned it). 
# If you call the one from your bin you can make changes to the main.nf file and run it 
# If you make changes to the main.nf file within the nextflow path it will error out
#
NEXTFLOW="/net/trapnell/vol1/home/readdf/nextflow"
NF_DEMUX="/net/trapnell/vol1/home/readdf/bin/bbi-dmux/main.nf"
#
# Get the path to the analyze output directory from
# the configuration file and set the Nextflow work
# directory to be in the analyze output directory.
# Delete the work directory after running the pipeline. It holds copies of all fastqs.
#
DEMUX_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.demux_out"){print$2}}' | sed 's/"//g'`
WORK_DIR="$DEMUX_DIR/work"

#
# Use REPORT_FIL, TRACE_FIL, and TIMELINE_FIL so that
# Nextflow writes some informative processing information
# in the analyze output directory.
#
REPORT_FIL=$WORKING_DIR/demux.report.html
TRACE_FIL=$WORKING_DIR/demux.trace.tsv
TIMELINE_FIL=$WORKING_DIR/demux.timeline.html

#
# Nextflow run parameters.
#
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"

mkdir -p $DEMUX_DIR
pushd $DEMUX_DIR

date > ./run_start.txt

#
# Run Nextflow demux pipeline.
#
nextflow run $NF_DEMUX $PARS

date > run_finish.txt

popd

