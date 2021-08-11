#!/bin/bash -ue
cat Lung_pre.log > merge_bams.log
    printf "** Start process 'merge_bams' at: $(date)

" >> merge_bams.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> merge_bams.log
    printf "    Process command:
        samtools merge Lung.bam Lung.fq.part2-L004.bam Lung.fq.part1-L001.bam Lung.fq.part1-L003.bam Lung.fq.part2-L001.bam Lung.fq.part2-L002.bam Lung.fq.part1-L002.bam Lung.fq.part1-L004.bam Lung.fq.part2-L003.bam

" >> merge_bams.log


    samtools merge -@ 8 Lung.bam Lung.fq.part2-L004.bam Lung.fq.part1-L001.bam Lung.fq.part1-L003.bam Lung.fq.part2-L001.bam Lung.fq.part2-L002.bam Lung.fq.part1-L002.bam Lung.fq.part1-L004.bam Lung.fq.part2-L003.bam


    printf "Lung	$(samtools view -c Lung.bam)" > Lung.read_count.txt

    printf "** End process 'merge_bams' at: $(date)

" >> merge_bams.log
