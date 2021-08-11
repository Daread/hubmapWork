#!/bin/bash -ue
cat count_umis_by_sample.log > make_matrix.log
    printf "** Start process 'make_matrix' at: $(date)

" >> make_matrix.log

    echo '    Process command:
        make_matrix.py <(zcat Pancreas.gz) 
            --gene_annotation "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations" 
            --key "Pancreas"
        cat /net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations > "Pancreas.gene_annotations.txt"  ' >> make_matrix.log


    make_matrix.py <(zcat Pancreas.gz) --gene_annotation "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations" --key "Pancreas"
    cat "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations" > "Pancreas.gene_annotations.txt"


    printf "
** End process 'make_matrix' at: $(date)

" >> make_matrix.log
