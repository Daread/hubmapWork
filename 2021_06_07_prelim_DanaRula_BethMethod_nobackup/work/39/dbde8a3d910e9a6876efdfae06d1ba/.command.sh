#!/bin/bash -ue
printf "** Start process 'sort_and_filter' for Pancreas.fq.part2-L004Aligned.out.bam at: $(date)

" > Pancreas.fq.part2-L004_piece.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> Pancreas.fq.part2-L004_piece.log
    printf "    Process command:
        samtools view -bh -q 30 -F 4 'Pancreas.fq.part2-L004Aligned.out.bam'
            | samtools sort -@ 10 - > 'Pancreas.fq.part2-L004.bam'

" >> Pancreas.fq.part2-L004_piece.log


    samtools view -bh -q 30 -F 4 "Pancreas.fq.part2-L004Aligned.out.bam"         | samtools sort -@ 10 -         > "Pancreas.fq.part2-L004.bam"


    printf "    Process stats:
        sort_and_filter starting reads: $(samtools view -c Pancreas.fq.part2-L004Aligned.out.bam)
        sort_and_filter ending reads  : $(samtools view -c Pancreas.fq.part2-L004.bam)

" >> Pancreas.fq.part2-L004_piece.log
    printf "** End process 'sort_and_filter' at: $(date)

" >> Pancreas.fq.part2-L004_piece.log

    cp Pancreas.fq.part2-L004_piece.log Pancreas.fq.part2-L004_sf.txt
