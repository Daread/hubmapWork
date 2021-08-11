#!/bin/bash -ue
printf "** Start process 'sort_and_filter' for Spleen.fq.part1-L002Aligned.out.bam at: $(date)

" > Spleen.fq.part1-L002_piece.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> Spleen.fq.part1-L002_piece.log
    printf "    Process command:
        samtools view -bh -q 30 -F 4 'Spleen.fq.part1-L002Aligned.out.bam'
            | samtools sort -@ 10 - > 'Spleen.fq.part1-L002.bam'

" >> Spleen.fq.part1-L002_piece.log


    samtools view -bh -q 30 -F 4 "Spleen.fq.part1-L002Aligned.out.bam"         | samtools sort -@ 10 -         > "Spleen.fq.part1-L002.bam"


    printf "    Process stats:
        sort_and_filter starting reads: $(samtools view -c Spleen.fq.part1-L002Aligned.out.bam)
        sort_and_filter ending reads  : $(samtools view -c Spleen.fq.part1-L002.bam)

" >> Spleen.fq.part1-L002_piece.log
    printf "** End process 'sort_and_filter' at: $(date)

" >> Spleen.fq.part1-L002_piece.log

    cp Spleen.fq.part1-L002_piece.log Spleen.fq.part1-L002_sf.txt
