#!/bin/bash -ue
printf "** Start process 'sort_and_filter' for HEK-L002Aligned.out.bam at: $(date)

" > HEK-L002_piece.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> HEK-L002_piece.log
    printf "    Process command:
        samtools view -bh -q 30 -F 4 'HEK-L002Aligned.out.bam'
            | samtools sort -@ 10 - > 'HEK-L002.bam'

" >> HEK-L002_piece.log


    samtools view -bh -q 30 -F 4 "HEK-L002Aligned.out.bam"         | samtools sort -@ 10 -         > "HEK-L002.bam"


    printf "    Process stats:
        sort_and_filter starting reads: $(samtools view -c HEK-L002Aligned.out.bam)
        sort_and_filter ending reads  : $(samtools view -c HEK-L002.bam)

" >> HEK-L002_piece.log
    printf "** End process 'sort_and_filter' at: $(date)

" >> HEK-L002_piece.log

    cp HEK-L002_piece.log HEK-L002_sf.txt
