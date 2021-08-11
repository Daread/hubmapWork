#!/bin/bash -ue
cat Liver_pre.log > merge_bams.log
    printf "** Start process 'merge_bams' at: $(date)

" >> merge_bams.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> merge_bams.log
    printf "    Process command:
        samtools merge Liver.bam Liver.fq.part2-L003.bam Liver.fq.part2-L004.bam Liver.fq.part1-L001.bam Liver.fq.part1-L002.bam Liver.fq.part2-L001.bam Liver.fq.part1-L004.bam Liver.fq.part1-L003.bam Liver.fq.part2-L002.bam

" >> merge_bams.log


    samtools merge -@ 8 Liver.bam Liver.fq.part2-L003.bam Liver.fq.part2-L004.bam Liver.fq.part1-L001.bam Liver.fq.part1-L002.bam Liver.fq.part2-L001.bam Liver.fq.part1-L004.bam Liver.fq.part1-L003.bam Liver.fq.part2-L002.bam


    printf "Liver	$(samtools view -c Liver.bam)" > Liver.read_count.txt

    printf "** End process 'merge_bams' at: $(date)

" >> merge_bams.log
