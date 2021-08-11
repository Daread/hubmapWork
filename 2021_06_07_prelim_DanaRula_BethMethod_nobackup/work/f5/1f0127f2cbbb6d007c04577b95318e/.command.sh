#!/bin/bash -ue
rmdup.py --bam split.REF_9.bam --output_bam out.bam

samtools view -c out.bam > split.REF_9.bam_umi_count.txt

bedtools bamtobed -i out.bam -split             | sort -k1,1 -k2,2n -k3,3n -S 5G             > "split.REF_9.bam.bed"

bedtools map         -a "split.REF_9.bam.bed"         -b "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.exons.bed"         -nonamecheck -s -f 0.95 -c 7 -o distinct -delim '|'     | bedtools map         -a - -b "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.genes.bed"         -nonamecheck -s -f 0.95 -c 4 -o distinct -delim '|'     | sort -k4,4 -k2,2n -k3,3n -S 5G    | datamash         -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8     | assign-reads-to-genes.py "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.genes.bed"     | awk '$3 == "exonic" || $3 == "intronic" {{
        split($1, arr, "|")
        printf "%s_%s_%s	%s	%s\n", arr[3], arr[4], arr[5], $2, $3
}}'     | sort -k1,1 -k2,2 -S 5G > "split.REF_9.bam_ga.txt"
