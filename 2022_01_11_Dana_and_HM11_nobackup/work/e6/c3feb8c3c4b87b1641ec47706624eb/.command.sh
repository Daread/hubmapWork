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
            "HEK.umi_counts.mtx"
            "HEK.cell_annotations.txt"
            "HEK.gene_annotations.txt"
            "/net/bbi/vol1/data/genomes_stage/barnyard/barnyard_rna//latest.genes.bed"
            "HEK"
            "50"
' >> make_cds.log


    make_cds.R         "HEK.umi_counts.mtx"        "HEK.cell_annotations.txt"        "HEK.gene_annotations.txt"        "/net/bbi/vol1/data/genomes_stage/barnyard/barnyard_rna//latest.genes.bed"        "HEK"        "50"


    printf "** End process 'make_cds' at: $(date)

" >> make_cds.log
