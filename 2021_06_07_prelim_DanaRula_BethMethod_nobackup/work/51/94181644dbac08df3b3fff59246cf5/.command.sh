#!/bin/bash -ue
cat Spleen_pre.log > merge_bams.log
    printf "** Start process 'merge_bams' at: $(date)

" >> merge_bams.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> merge_bams.log
    printf "    Process command:
        samtools merge Spleen.bam Spleen.fq.part1-L001.bam Spleen.fq.part2-L001.bam Spleen.fq.part2-L002.bam Spleen.fq.part2-L004.bam Spleen.fq.part1-L003.bam Spleen.fq.part2-L003.bam Spleen.fq.part1-L002.bam Spleen.fq.part1-L004.bam

" >> merge_bams.log


    samtools merge -@ 8 Spleen.bam Spleen.fq.part1-L001.bam Spleen.fq.part2-L001.bam Spleen.fq.part2-L002.bam Spleen.fq.part2-L004.bam Spleen.fq.part1-L003.bam Spleen.fq.part2-L003.bam Spleen.fq.part1-L002.bam Spleen.fq.part1-L004.bam


    printf "Spleen	$(samtools view -c Spleen.bam)" > Spleen.read_count.txt

    printf "** End process 'merge_bams' at: $(date)

" >> merge_bams.log
