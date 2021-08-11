#!/bin/bash -ue
cat start.log > trim.log
    printf "** Start process 'trim_fastqs' for Spleen.fq.part2-L001.fastq.gz at: $(date)

" > piece.log
    printf "    Process versions:
        " >> piece.log
    python --version &>> piece.log
    printf "        trim_galore $(trim_galore -v | grep version | awk '{$1=$1;print}')
        cutadapt version $(cutadapt --version)

" >> piece.log

    printf "    Process command:
        trim_galore Spleen.fq.part2-L001.fastq.gz -a AAAAAAAA --three_prime_clip_R1 1
            --gzip -o ./trim_out/

    Process output:
" >> piece.log


    mkdir trim_out
    trim_galore Spleen.fq.part2-L001.fastq.gz         -a AAAAAAAA         --three_prime_clip_R1 1         --gzip         -o ./trim_out/


    cat trim_out/*trimming_report.txt | sed '/Overview of/,/RUN/{//!d}' | sed 's/Overview of removed sequences//' >> piece.log
    printf "** End process 'trim_fastqs' at: $(date)

" >> piece.log
    cp piece.log Spleen.fq.part2-L001_trim.txt
    cat piece.log >> trim.log
