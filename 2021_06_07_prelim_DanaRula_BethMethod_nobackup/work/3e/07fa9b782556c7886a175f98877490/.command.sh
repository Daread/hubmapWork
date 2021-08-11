#!/bin/bash -ue
cat trim.log > align.log
    printf "** Start process 'align_reads' for Liver.fq.part2-L001_trimmed.fq.gz at: $(date)

" > piece.log
    printf "    Process versions:
        $(STAR --version)

" >> piece.log

    printf "    Process command:
        STAR --runThreadN 8 --genomeDir /net/bbi/vol1/data/genomes_stage/human/human_star
            --readFilesIn Liver.fq.part2-L001_trimmed.fq.gz --readFilesCommand zcat 
            --outFileNamePrefix ./align_out/Liver.fq.part2-L001 --outSAMtype BAM Unsorted 
            --outSAMmultNmax 1 --outSAMstrandField intronMotif


    Reference genome information:
      $(grep fastq_url /net/bbi/vol1/data/genomes_stage/human/human_star/../*gsrc/record.out | awk '{$1=$2=""; print $0}')
        FASTA download date: $(grep fastq_download_date /net/bbi/vol1/data/genomes_stage/human/human_star/../*gsrc/record.out | awk '{$1=$2=""; print $0}')
        Non REF sequences removed.

      $(grep gtf_url /net/bbi/vol1/data/genomes_stage/human/human_star/../*gsrc/record.out | awk '{$1=$2=""; print $0}')
        GTF download date: $(grep gtf_download_date /net/bbi/vol1/data/genomes_stage/human/human_star/../*gsrc/record.out | awk '{$1=$2=""; print $0}')
        $(grep gtf_include_biotypes /net/bbi/vol1/data/genomes_stage/human/human_star/record.out | awk '{$1=$2=""; print $0}')

    Process output:
" >> piece.log


    mkdir align_out
    STAR         --runThreadN 8         --genomeDir /net/bbi/vol1/data/genomes_stage/human/human_star         --readFilesIn Liver.fq.part2-L001_trimmed.fq.gz         --readFilesCommand zcat         --outFileNamePrefix ./align_out/Liver.fq.part2-L001          --outSAMtype BAM Unsorted         --outSAMmultNmax 1         --outSAMstrandField intronMotif


    cat align_out/*Log.final.out >> piece.log

    printf "
** End process 'align_reads' at: $(date)

" >> piece.log

    cp piece.log Liver.fq.part2-L001_align.txt
    cat piece.log >> align.log
