#!/bin/bash -ue
head -n 2 run_scrublet.log > HEK_full.log
    printf "Nextflow version: 20.07.1
" >> HEK_full.log
    printf "Pipeline version: 2.2.0
" >> HEK_full.log
    printf "Git Repository, Version, Commit ID, Session ID: null, null, null, aa6515ab-3393-4d01-9e06-49f0f557f4e4

" >> HEK_full.log
    printf "Command:
nextflow run /net/trapnell/vol1/home/readdf/bin/bbi-sci/main.nf -c /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup/experiment.config -w /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup/work -with-report /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup/analysis.report.html -with-trace /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup/analysis.trace.tsv -with-timeline /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup/analysis.timeline.html -resume

" >> HEK_full.log
    printf "***** PARAMETERS *****: 

" >> HEK_full.log
    printf "    params.run_dir:               /net/shendure/vol9/seq/NEXTSEQ/220110_NS500488_1276_AH23G5BGXL
" >> HEK_full.log
    printf "    params.output_dir:            /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup
" >> HEK_full.log
    printf "    params.sample_sheet:          /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup/Hm11_dana_test_samplesheet.csv
" >> HEK_full.log
    printf "    params.demux_out:             /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup/demux_out
" >> HEK_full.log
    printf "    params.level:                 2
" >> HEK_full.log
    printf "    params.max_cores:             16
" >> HEK_full.log
    printf "    params.samples:               [HEK]
" >> HEK_full.log
    printf "    params.star_file:             /net/trapnell/vol1/home/readdf/bin/bbi-sci/bin/star_file.txt
" >> HEK_full.log
    printf "    params.gene_file:             /net/trapnell/vol1/home/readdf/bin/bbi-sci/bin/gene_file.txt
" >> HEK_full.log
    printf "    params.umi_cutoff:            50
" >> HEK_full.log
    printf "    params.rt_barcode_file:       /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/misc/RT_2level_indices_plate_format
" >> HEK_full.log
    printf "    params.hash_list:             false
" >> HEK_full.log
    printf "    params.max_wells_per_sample:  20

" >> HEK_full.log
    printf "    params.garnett_file:          false

" >> HEK_full.log
    printf "    params.skip_doublet_detect:   false

" >> HEK_full.log

    tail -n +2 run_scrublet.log >> HEK_full.log
    printf "
** End processes generate qc metrics and dashboard at: $(date)

" >> HEK_full.log
    printf "***** END PIPELINE *****: 

" >> HEK_full.log
    filename=HEK_full.log

    # Trimming:
    trim_start=`cat $filename | grep 'sequences processed in total' | awk -F ' ' '{sum += $1} END {print sum}'`
    trim_lost=`cat $filename | grep 'Sequences removed because they became shorter' | awk -F ' ' '{sum += $14} END {print sum}'`
    trim_end=$(($trim_start - $trim_lost))

    # Alignment:
    align_start=`cat $filename | grep 'Number of input reads' | awk -F '|' '{sum += $2} END {print sum}'`
    align_mapped=`cat $filename | grep 'Uniquely mapped reads number' | awk -F '|' '{sum += $2} END {print sum}'`
    align_totals=($(cat $filename | grep 'Number of input reads' | cut -d "|" -f 2 | awk '{print $1}'))
    align_multimapped=`cat $filename | grep 'Number of reads mapped to multiple loci' |  awk -F '|' '{sum += $2} END {print sum}'`
    align_too_short_arr=($(cat $filename | grep 'unmapped: too short' | cut -d "|" -f 2 | tr '%' ' ' | awk '{$1=$1/100;print}'))
    align_too_short=`a=0
    for i in ${align_too_short_arr[@]}
    do
        echo "${align_too_short_arr[$a]} * ${align_totals[$a]}" | bc
        a=$((a+1))
    done | awk '{sum += $1} END {printf "%1.0f", sum}'`

    # Sort and Filter:
    sf_start=`cat $filename | grep 'sort_and_filter starting reads' | awk -F ':' '{sum += $2} END {print sum}'`
    sf_end=`cat $filename | grep 'sort_and_filter ending reads' | awk -F ':' '{sum += $2} END {print sum}'`

    # Dups:
    dup_start=`cat $filename | grep 'remove_dups starting reads' | awk -F ':' '{sum += $2} END {print sum}'`
    dup_end=`cat $filename | grep 'remove_dups ending reads' | awk -F ':' '{sum += $2} END {print sum}'`

    # Assignment:
    assigned_exonic=`cat $filename | grep '    exonic     ' | awk -F ' ' '{sum += $2} END {print sum}'`
    assigned_intronic=`cat $filename | grep '    intronic     ' | awk -F ' ' '{sum += $2} END {print sum}'`
    assigned_end=$(($assigned_exonic + $assigned_intronic))

    # In real cells:
    reads_in_cells=`cat $filename | grep 'Total reads in cells with > 100 reads' | awk -F ':' '{sum += $2} END {print sum}'`

    printf "
            \"HEK\": {
            \"sample\": \"HEK\",
            \"alignment_start\" : \"$align_start\",
            \"alignment_mapped\" : \"$align_mapped\",
            \"align_multimapped\" : \"$align_multimapped\",
            \"align_too_short\" : \"$align_too_short\",
            \"sf_start\" : \"$sf_start\",
            \"sf_end\" : \"$sf_end\",
            \"dup_start\" : \"$dup_start\",
            \"dup_end\" : \"$dup_end\",
            \"assigned_exonic\" : \"$assigned_exonic\",
            \"assigned_intronic\" : \"$assigned_intronic\",
            \"reads_in_cells\" : \"$reads_in_cells\" }
      " > HEK_log_data.txt


    printf "***** PIPELINE READ STATS *****: 

" >> HEK_read_metrics.log

    printf "%20s %20s %20s %20s %20s
" "Process" "Starting reads" "Ending reads" "% lost" "% of total lost" >> HEK_read_metrics.log
    printf "========================================================================================================
" >> HEK_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f
" "Trimming" $trim_start $trim_end $(echo "($trim_start - $trim_end)/$trim_start * 100" | bc -l ) $(echo "($trim_start - $trim_end)/$trim_start * 100" | bc -l ) >> HEK_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f
" "Alignment" $align_start $sf_start $(echo "($align_start - $sf_start)/$align_start * 100" | bc -l ) $(echo "($align_start - $sf_start)/$trim_start * 100" | bc -l ) >> HEK_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f
" "Filtering" $sf_start $sf_end $(echo "($sf_start - $sf_end)/$sf_start * 100" | bc -l ) $(echo "($sf_start - $sf_end)/$trim_start * 100" | bc -l ) >> HEK_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f
" "Deduplication" $dup_start $dup_end $(echo "($dup_start - $dup_end)/$dup_start * 100" | bc -l ) $(echo "($dup_start - $dup_end)/$trim_start * 100" | bc -l ) >> HEK_read_metrics.log
    printf "%20s %20s %20s %20.2f %20.2f
" "Gene assignment" $dup_end $assigned_end $(echo "($dup_end - $assigned_end)/$dup_end * 100" | bc -l ) $(echo "($dup_end - $assigned_end)/$trim_start * 100" | bc -l ) >> HEK_read_metrics.log

    printf "
Alignment details: 
" >> HEK_read_metrics.log
    printf "%25s %20s %20s
" "" "Count" "Percent"  >> HEK_read_metrics.log
    printf "========================================================================================================
" >> HEK_read_metrics.log
    printf "%25s %20s %20s
" "Total reads processed:" $align_start "" >> HEK_read_metrics.log
    printf "%25s %20s %20.2f
" "Reads uniquely mapped:" $align_mapped $(echo "($align_mapped)/$align_start * 100" | bc -l ) >> HEK_read_metrics.log
    printf "%25s %20s %20.2f
" "Reads multi-mapped:" $align_multimapped $(echo "($align_multimapped)/$align_start * 100" | bc -l )  >> HEK_read_metrics.log
    printf "%25s %20s %20.2f
" "Reads too short:" $align_too_short $(echo "($align_too_short)/$align_start * 100" | bc -l ) >> HEK_read_metrics.log


    cat HEK_read_metrics.log >> HEK_full.log
