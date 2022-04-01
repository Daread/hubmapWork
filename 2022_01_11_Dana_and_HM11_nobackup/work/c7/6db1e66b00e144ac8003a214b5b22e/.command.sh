#!/bin/bash -ue
cat count_umis_by_sample.log > make_matrix.log
    printf "** Start process 'make_matrix' at: $(date)

" >> make_matrix.log

    echo '    Process command:
        make_matrix.py <(zcat HEK.gz) 
            --gene_annotation "/net/bbi/vol1/data/genomes_stage/barnyard/barnyard_rna//latest.gene.annotations" 
            --key "HEK"
        cat /net/bbi/vol1/data/genomes_stage/barnyard/barnyard_rna//latest.gene.annotations > "HEK.gene_annotations.txt"  ' >> make_matrix.log


    make_matrix.py <(zcat HEK.gz) --gene_annotation "/net/bbi/vol1/data/genomes_stage/barnyard/barnyard_rna//latest.gene.annotations" --key "HEK"
    cat "/net/bbi/vol1/data/genomes_stage/barnyard/barnyard_rna//latest.gene.annotations" > "HEK.gene_annotations.txt"


    printf "
** End process 'make_matrix' at: $(date)

" >> make_matrix.log
