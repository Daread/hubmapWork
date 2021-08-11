#!/bin/bash -ue
printf "** Start process 'sort_and_filter' for Lung.fq.part2-L002Aligned.out.bam at: $(date)

" > Lung.fq.part2-L002_piece.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> Lung.fq.part2-L002_piece.log
    printf "    Process command:
        samtools view -bh -q 30 -F 4 'Lung.fq.part2-L002Aligned.out.bam'
            | samtools sort -@ 10 - > 'Lung.fq.part2-L002.bam'

" >> Lung.fq.part2-L002_piece.log


    samtools view -bh -q 30 -F 4 "Lung.fq.part2-L002Aligned.out.bam"         | samtools sort -@ 10 -         > "Lung.fq.part2-L002.bam"


    printf "    Process stats:
        sort_and_filter starting reads: $(samtools view -c Lung.fq.part2-L002Aligned.out.bam)
        sort_and_filter ending reads  : $(samtools view -c Lung.fq.part2-L002.bam)

" >> Lung.fq.part2-L002_piece.log
    printf "** End process 'sort_and_filter' at: $(date)

" >> Lung.fq.part2-L002_piece.log

    cp Lung.fq.part2-L002_piece.log Lung.fq.part2-L002_sf.txt
