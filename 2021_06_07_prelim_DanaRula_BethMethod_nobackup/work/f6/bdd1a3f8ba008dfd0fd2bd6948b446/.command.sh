#!/bin/bash -ue
cat make_cds.log > apply_garnett.log
    printf "** Start process 'apply_garnett' at: $(date)

" >> apply_garnett.log
    mkdir new_cds
    echo "No Garnett classifier provided for this sample" > garnett_error.txt
    if [ false == 'false' ]
    then
        cp Pancreas_cds.RDS new_cds/
    else
        apply_garnett.R Pancreas_cds.RDS false Pancreas
    fi

    cat garnett_error.txt >> apply_garnett.log
    printf "
** End process 'apply_garnett' at: $(date)

" >> apply_garnett.log
