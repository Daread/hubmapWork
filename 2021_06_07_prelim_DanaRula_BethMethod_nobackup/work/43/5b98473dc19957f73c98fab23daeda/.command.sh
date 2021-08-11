#!/bin/bash -ue
cat merge_bams.log > remove_dups.log
    printf "** Start processes 'remove duplicates, assign_genes, merge_assignment' at: $(date)

" >> remove_dups.log
    printf "    Process versions:
        $(bedtools --version)
        $(samtools --version | tr '
' ' ')
        $(bamtools --version | grep bamtools)
        $(python --version)

" >> remove_dups.log

    echo '    Process command:
        mkdir split_bams
        bamtools split -in Spleen.bam -reference -stub split_bams/split

        rmdup.py --bam in_bam --output_bam out.bam
    
        samtools view -c out.bam > split_bam_umi_count.txt

        bedtools bamtobed -i out.bam -split
                | sort -k1,1 -k2,2n -k3,3n -S 5G
                > "in_bam.bed"

        bedtools map
            -a in_bam.bed
            -b exon_index
            -nonamecheck -s -f 0.95 -c 7 -o distinct -delim "|"
        | bedtools map
            -a - -b gene_index
            -nonamecheck -s -f 0.95 -c 4 -o distinct -delim "|"
        | sort -k4,4 -k2,2n -k3,3n -S 5G
        | datamash
            -g 4 first 1 first 2 last 3 first 5 first 6 collapse 7 collapse 8
        | assign-reads-to-genes.py gene_index
        | awk $3 == "exonic" || $3 == "intronic" {{
                split($1, arr, "|")
                printf "%s_%s_%s	%s	%s\n", arr[3], arr[4], arr[5], $2, $3
        }}
        | sort -k2,2 -k1,1 -S 5G > in_bam.txt

        cat logfile > merge_assignment.log
        cat split_bed > key.bed
        sort -m -k1,1 -k2,2 split_gene_assign > key_ga.txt

        datamash -g 1,2 count 2 < key_ga.txt 
        | gzip > key.gz
        ' >> remove_dups.log

    printf "    Process stats:
        remove_dups starting reads: $(samtools view -c Spleen.bam)" >> remove_dups.log


    mkdir split_bams
    bamtools split -in Spleen.bam -reference -stub split_bams/split
    cd split_bams
    if [[ $(ls | grep "_[0-9A-Za-z\.]\{3,\}.bam$") ]]; then 
        ls | grep "_[0-9A-Za-z\.]\{3,\}.bam$" | samtools merge split.REFnonstand.bam -b -
        ls | grep "_[0-9A-Za-z\.]\{3,\}.bam$" | xargs -d"\n" rm
        mv split.REFnonstand.bam split.REF_nonstand.bam
    fi
