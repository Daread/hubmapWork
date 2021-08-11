#!/bin/bash -ue
printf "** Start process 'sort_and_filter' for Liver.fq.part1-L001Aligned.out.bam at: $(date)

" > Liver.fq.part1-L001_piece.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> Liver.fq.part1-L001_piece.log
    printf "    Process command:
        samtools view -bh -q 30 -F 4 'Liver.fq.part1-L001Aligned.out.bam'
            | samtools sort -@ 10 - > 'Liver.fq.part1-L001.bam'

" >> Liver.fq.part1-L001_piece.log


    samtools view -bh -q 30 -F 4 "Liver.fq.part1-L001Aligned.out.bam"         | samtools sort -@ 10 -         > "Liver.fq.part1-L001.bam"


    printf "    Process stats:
        sort_and_filter starting reads: $(samtools view -c Liver.fq.part1-L001Aligned.out.bam)
        sort_and_filter ending reads  : $(samtools view -c Liver.fq.part1-L001.bam)

" >> Liver.fq.part1-L001_piece.log
    printf "** End process 'sort_and_filter' at: $(date)

" >> Liver.fq.part1-L001_piece.log

    cp Liver.fq.part1-L001_piece.log Liver.fq.part1-L001_sf.txt
