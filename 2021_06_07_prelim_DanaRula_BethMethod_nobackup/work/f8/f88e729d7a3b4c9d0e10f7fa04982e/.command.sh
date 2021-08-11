#!/bin/bash -ue
cat remove_dups.log > merge_assignment.log
    cat split.REF_1.bam.bed split.REF_10.bam.bed split.REF_13.bam.bed split.REF_12.bam.bed split.REF_11.bam.bed split.REF_16.bam.bed split.REF_17.bam.bed split.REF_15.bam.bed split.REF_14.bam.bed split.REF_20.bam.bed split.REF_21.bam.bed split.REF_18.bam.bed split.REF_2.bam.bed split.REF_19.bam.bed split.REF_5.bam.bed split.REF_4.bam.bed split.REF_7.bam.bed split.REF_6.bam.bed split.REF_22.bam.bed split.REF_X.bam.bed split.REF_nonstand.bam.bed split.REF_Y.bam.bed split.REF_3.bam.bed split.REF_MT.bam.bed split.REF_9.bam.bed split.REF_8.bam.bed > "Spleen.bed"
    sort -m -k1,1 -k2,2 split.REF_1.bam_ga.txt split.REF_10.bam_ga.txt split.REF_13.bam_ga.txt split.REF_12.bam_ga.txt split.REF_11.bam_ga.txt split.REF_16.bam_ga.txt split.REF_17.bam_ga.txt split.REF_15.bam_ga.txt split.REF_14.bam_ga.txt split.REF_20.bam_ga.txt split.REF_21.bam_ga.txt split.REF_18.bam_ga.txt split.REF_2.bam_ga.txt split.REF_19.bam_ga.txt split.REF_5.bam_ga.txt split.REF_4.bam_ga.txt split.REF_7.bam_ga.txt split.REF_6.bam_ga.txt split.REF_22.bam_ga.txt split.REF_X.bam_ga.txt split.REF_nonstand.bam_ga.txt split.REF_Y.bam_ga.txt split.REF_3.bam_ga.txt split.REF_MT.bam_ga.txt split.REF_9.bam_ga.txt split.REF_8.bam_ga.txt > "Spleen_ga.txt"

    datamash -g 1,2 count 2 < "Spleen_ga.txt"     | gzip > "Spleen.gz"


    umi=`cat split.REF_1.bam_umi_count.txt split.REF_10.bam_umi_count.txt split.REF_13.bam_umi_count.txt split.REF_12.bam_umi_count.txt split.REF_11.bam_umi_count.txt split.REF_16.bam_umi_count.txt split.REF_17.bam_umi_count.txt split.REF_15.bam_umi_count.txt split.REF_14.bam_umi_count.txt split.REF_20.bam_umi_count.txt split.REF_21.bam_umi_count.txt split.REF_18.bam_umi_count.txt split.REF_2.bam_umi_count.txt split.REF_19.bam_umi_count.txt split.REF_5.bam_umi_count.txt split.REF_4.bam_umi_count.txt split.REF_7.bam_umi_count.txt split.REF_6.bam_umi_count.txt split.REF_22.bam_umi_count.txt split.REF_X.bam_umi_count.txt split.REF_nonstand.bam_umi_count.txt split.REF_Y.bam_umi_count.txt split.REF_3.bam_umi_count.txt split.REF_MT.bam_umi_count.txt split.REF_9.bam_umi_count.txt split.REF_8.bam_umi_count.txt | awk '{ sum += $1 } END { print sum }'`
    read=`cut -f2 Spleen.read_count.txt`
    perc=$(echo "100.0 * (1 - $umi/$read)" | bc -l)
    printf "%-18s   %10d    %10d    %7.1f\n" Spleen $read $umi $perc     >"Spleen.duplication_rate_stats.txt"

    printf "
        remove_dups ending reads  : $(wc -l Spleen.bed | awk '{print $1;}')


        Read assignments:
$(awk '{count[$3]++} END {for (word in count) { printf "            %-20s %10i\n", word, count[word]}}' Spleen_ga.txt)

" >> merge_assignment.log

    printf "** End processes 'remove duplicates, assign_genes, merge_assignment' at: $(date)

" >> merge_assignment.log
