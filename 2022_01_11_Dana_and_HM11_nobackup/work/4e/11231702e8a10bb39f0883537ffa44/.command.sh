#!/bin/bash -ue
printf "** Start process 'sort_and_filter' for HEK-L004Aligned.out.bam at: $(date)

" > HEK-L004_piece.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> HEK-L004_piece.log
    printf "    Process command:
        samtools view -bh -q 30 -F 4 'HEK-L004Aligned.out.bam'
            | samtools sort -@ 10 - > 'HEK-L004.bam'

" >> HEK-L004_piece.log


    samtools view -bh -q 30 -F 4 "HEK-L004Aligned.out.bam"         | samtools sort -@ 10 -         > "HEK-L004.bam"


    printf "    Process stats:
        sort_and_filter starting reads: $(samtools view -c HEK-L004Aligned.out.bam)
        sort_and_filter ending reads  : $(samtools view -c HEK-L004.bam)

" >> HEK-L004_piece.log
    printf "** End process 'sort_and_filter' at: $(date)

" >> HEK-L004_piece.log

    cp HEK-L004_piece.log HEK-L004_sf.txt
