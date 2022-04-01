#!/bin/bash -ue
printf "BBI bbi-sci Pipeline Log

" > start.log
    printf "Run started at: $(date)

" >> start.log

    printf "***** BEGIN PIPELINE *****: 

" >> start.log
    printf "** Start process 'check_sample_sheet' at: $(date)

" >> start.log
    printf "    Process versions:
        $(python --version)

" >> start.log
    printf "    Process command:
        check_sample_sheet.py 
            --sample_sheet /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup/Hm11_dana_test_samplesheet.csv 
            --star_file /net/trapnell/vol1/home/readdf/bin/bbi-sci/bin/star_file.txt
            --level 2 --rt_barcode_file /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/misc/RT_2level_indices_plate_format
            --max_wells_per_samp 20

" >> start.log


    check_sample_sheet.py --sample_sheet /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_11_Dana_and_HM11_nobackup/Hm11_dana_test_samplesheet.csv --star_file /net/trapnell/vol1/home/readdf/bin/bbi-sci/bin/star_file.txt         --level 2 --rt_barcode_file /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/misc/RT_2level_indices_plate_format         --max_wells_per_samp 20


    printf "** End process 'check_sample_sheet' at: $(date)

" >> start.log
    cp start.log start.txt
