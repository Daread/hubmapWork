#!/bin/bash -ue
mkdir -p fastqs_trim
java -Xmx1G -jar /net/trapnell/vol1/home/readdf/bin/bbi-sciatac-demux/src/Trimmomatic-0.36/trimmomatic-0.36.jar        PE        -threads 4        W144.heart.apex.s1-RUN001_L004_R1.fastq.gz W144.heart.apex.s1-RUN001_L004_R2.fastq.gz        W144.heart.apex.s1-RUN001_L004_R1.trimmed.fastq.gz        W144.heart.apex.s1-RUN001_L004_R1.trimmed_unpaired.fastq.gz        W144.heart.apex.s1-RUN001_L004_R2.trimmed.fastq.gz        W144.heart.apex.s1-RUN001_L004_R2.trimmed_unpaired.fastq.gz        ILLUMINACLIP:/net/trapnell/vol1/home/readdf/bin/bbi-sciatac-demux/src/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:1:true        TRAILING:3        SLIDINGWINDOW:4:10        MINLEN:20
