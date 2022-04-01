#!/bin/bash -ue
cat HEK_pre.log > merge_bams.log
    printf "** Start process 'merge_bams' at: $(date)

" >> merge_bams.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> merge_bams.log
    printf "    Process command:
        samtools merge HEK.bam HEK-L002.bam HEK-L001.bam HEK-L004.bam HEK-L003.bam

" >> merge_bams.log


    samtools merge -@ 8 HEK.bam HEK-L002.bam HEK-L001.bam HEK-L004.bam HEK-L003.bam


    printf "HEK	$(samtools view -c HEK.bam)" > HEK.read_count.txt

    printf "** End process 'merge_bams' at: $(date)

" >> merge_bams.log
