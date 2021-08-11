#!/bin/bash -ue
cat make_matrix.log > make_cds.log
    printf "** Start process 'make_cds' at: $(date)

" >> make_cds.log
    printf "    Process versions:
        $(R --version | grep 'R version')
            monocle3 version $(Rscript -e 'packageVersion("monocle3")')

" >> make_cds.log
    echo '    Process command:
        make_cds.R
            "Spleen.umi_counts.mtx"
            "Spleen.cell_annotations.txt"
            "Spleen.gene_annotations.txt"
            "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.genes.bed"
            "Spleen"
            "1"
' >> make_cds.log


    make_cds.R         "Spleen.umi_counts.mtx"        "Spleen.cell_annotations.txt"        "Spleen.gene_annotations.txt"        "/net/bbi/vol1/data/genomes_stage/human/human_rna//latest.genes.bed"        "Spleen"        "1"


    printf "** End process 'make_cds' at: $(date)

" >> make_cds.log
