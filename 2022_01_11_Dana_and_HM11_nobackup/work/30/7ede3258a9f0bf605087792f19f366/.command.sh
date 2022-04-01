#!/bin/bash -ue
cat merge_assignment.log > count_umis_by_sample.log
    printf "** Start process 'count_umis_by_sample' at: $(date)

" >> count_umis_by_sample.log
    printf "    Process versions:
        $(python --version)

" >> count_umis_by_sample.log
    printf "    Process command:
        tabulate_per_cell_counts.py
            --gene_assignment_files "HEK_ga.txt"
            --all_counts_file "HEK.UMIs.per.cell.barcode.txt"
            --intron_counts_file "HEK.UMIs.per.cell.barcode.intronic.txt"

"      >> count_umis_by_sample.log


    tabulate_per_cell_counts.py         --gene_assignment_files "HEK_ga.txt"         --all_counts_file "HEK.UMIs.per.cell.barcode.txt"         --intron_counts_file "HEK.UMIs.per.cell.barcode.intronic.txt"


    printf "    Process stats:
        Total cells                            : $(wc -l HEK.UMIs.per.cell.barcode.txt | awk '{print $1;}')
        Total cells > 100 reads                : $(awk '$2>100{c++} END{print c+0}' HEK.UMIs.per.cell.barcode.txt)
        Total cells > 1000 reads               : $(awk '$2>1000{c++} END{print c+0}' HEK.UMIs.per.cell.barcode.txt)
        Total reads in cells with > 100 reads  : $(awk '$2>100{c=c+$2} END{print c+0}' HEK.UMIs.per.cell.barcode.txt)

" >> count_umis_by_sample.log

    printf "** End process 'count_umis_by_sample' at: $(date)

" >> count_umis_by_sample.log
