#!/bin/bash -ue
cat Pancreas_pre.log > merge_bams.log
    printf "** Start process 'merge_bams' at: $(date)

" >> merge_bams.log
    printf "    Process versions:
        $(samtools --version | tr '
' ' ')

" >> merge_bams.log
    printf "    Process command:
        samtools merge Pancreas.bam Pancreas.fq.part1-L004.bam Pancreas.fq.part2-L001.bam Pancreas.fq.part1-L001.bam Pancreas.fq.part1-L002.bam Pancreas.fq.part2-L003.bam Pancreas.fq.part2-L002.bam Pancreas.fq.part2-L004.bam Pancreas.fq.part1-L003.bam

" >> merge_bams.log


    samtools merge -@ 8 Pancreas.bam Pancreas.fq.part1-L004.bam Pancreas.fq.part2-L001.bam Pancreas.fq.part1-L001.bam Pancreas.fq.part1-L002.bam Pancreas.fq.part2-L003.bam Pancreas.fq.part2-L002.bam Pancreas.fq.part2-L004.bam Pancreas.fq.part1-L003.bam


    printf "Pancreas	$(samtools view -c Pancreas.bam)" > Pancreas.read_count.txt

    printf "** End process 'merge_bams' at: $(date)

" >> merge_bams.log
