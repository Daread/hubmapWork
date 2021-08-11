#!/bin/bash -ue
cat count_umis_by_sample.log > make_matrix.log
    printf "** Start process 'make_matrix' at: $(date)

" >> make_matrix.log

    echo '    Process command:
        make_matrix.py <(zcat Liver.gz) 
            --gene_annotation "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations" 
            --key "Liver"
        cat /net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations > "Liver.gene_annotations.txt"  ' >> make_matrix.log


    make_matrix.py <(zcat Liver.gz) --gene_annotation "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations" --key "Liver"
    cat "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.gene.annotations" > "Liver.gene_annotations.txt"


    printf "
** End process 'make_matrix' at: $(date)

" >> make_matrix.log
