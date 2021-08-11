#!/bin/bash -ue
printf "** Start process 'sort_and_filter' for Spleen.fq.part2-L004Aligned.out.bam at: $(date)

" > Spleen.fq.part2-L004_piece.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> Spleen.fq.part2-L004_piece.log
    printf "    Process command:
        samtools view -bh -q 30 -F 4 'Spleen.fq.part2-L004Aligned.out.bam'
            | samtools sort -@ 10 - > 'Spleen.fq.part2-L004.bam'

" >> Spleen.fq.part2-L004_piece.log


    samtools view -bh -q 30 -F 4 "Spleen.fq.part2-L004Aligned.out.bam"         | samtools sort -@ 10 -         > "Spleen.fq.part2-L004.bam"


    printf "    Process stats:
        sort_and_filter starting reads: $(samtools view -c Spleen.fq.part2-L004Aligned.out.bam)
        sort_and_filter ending reads  : $(samtools view -c Spleen.fq.part2-L004.bam)

" >> Spleen.fq.part2-L004_piece.log
    printf "** End process 'sort_and_filter' at: $(date)

" >> Spleen.fq.part2-L004_piece.log

    cp Spleen.fq.part2-L004_piece.log Spleen.fq.part2-L004_sf.txt
