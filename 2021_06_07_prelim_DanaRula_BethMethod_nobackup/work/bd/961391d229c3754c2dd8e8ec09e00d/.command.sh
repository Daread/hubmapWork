#!/bin/bash -ue
printf "** Start process 'sort_and_filter' for Lung.fq.part1-L003Aligned.out.bam at: $(date)

" > Lung.fq.part1-L003_piece.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> Lung.fq.part1-L003_piece.log
    printf "    Process command:
        samtools view -bh -q 30 -F 4 'Lung.fq.part1-L003Aligned.out.bam'
            | samtools sort -@ 10 - > 'Lung.fq.part1-L003.bam'

" >> Lung.fq.part1-L003_piece.log


    samtools view -bh -q 30 -F 4 "Lung.fq.part1-L003Aligned.out.bam"         | samtools sort -@ 10 -         > "Lung.fq.part1-L003.bam"


    printf "    Process stats:
        sort_and_filter starting reads: $(samtools view -c Lung.fq.part1-L003Aligned.out.bam)
        sort_and_filter ending reads  : $(samtools view -c Lung.fq.part1-L003.bam)

" >> Lung.fq.part1-L003_piece.log
    printf "** End process 'sort_and_filter' at: $(date)

" >> Lung.fq.part1-L003_piece.log

    cp Lung.fq.part1-L003_piece.log Lung.fq.part1-L003_sf.txt