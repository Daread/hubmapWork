#!/bin/bash

#
# Current date and time.
#
NOW=`date '+%Y%m%d.%H%M%S'`

#
# Path to the Nextflow processing run configuration file.
#
# PWD=`pwd`
PWD="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_03_24_reprocessHeartAlone_nobackup"
CONFIG_FILE="$PWD/experiment.config"

#
# Nextflow executable and pipeline script locations.
#
# NEXTFLOW="/net/trapnell/vol1/home/readdf/nextflow"
NEXTFLOW="/net/trapnell/vol1/home/readdf/nextflow"
# NEXTFLOW="/net/trapnell/vol1/home/readdf/bin/redownloadNextflow_2022_03_25/nextflow"
NF_DEMUX="/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-demux/main.nf"

# Greg locations
# NEXTFLOW="/net/gs/vol1/home/bge/bin/nextflow"
# NF_DEMUX="/net/trapnell/vol1/home/bge/git/bbi-sciatac-demux/main.nf"


#
# Get the path to the demux output directory from
# the configuration file and set the Nextflow work
# directory to be in the demux output directory.
# Notes:
#   o  bbi-sciatac-demux/main.nf and bbi-scripts-analyze/main.nf also define
#      DEMUX_DIR as ${params.output_dir}/demux_out so changing DEMUX_DIR here
#      requires changing those files as well.
#
OUTPUT_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.output_dir"){print$2}}' | sed 's/"//g'`
DEMUX_DIR="$OUTPUT_DIR/demux_out"
WORK_DIR="$DEMUX_DIR/work"

#
# Use REPORT_FIL, TRACE_FIL, and TIMELINE_FIL so that
# Nextflow writes some informative processing information
# in the analyze output directory.
#
REPORT_FIL=$DEMUX_DIR/run_reports/demux.report.html
TRACE_FIL=$DEMUX_DIR/run_reports/demux.trace.tsv
TIMELINE_FIL=$DEMUX_DIR/run_reports/demux.timeline.html

#
# Nextflow run parameters.
#
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"

mkdir -p $DEMUX_DIR/run_reports
pushd $DEMUX_DIR

date > run_reports/run_start.${NOW}.txt

#
# Run Nextflow sci-ATAC demux pipeline.
#
$NEXTFLOW run $NF_DEMUX $PARS

date > run_reports/run_finish.${NOW}.txt

popd
