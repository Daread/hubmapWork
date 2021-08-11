#!/bin/bash -ue
cat apply_garnett.log > run_scrublet.log
    printf "** Start process 'run_scrublet' at: $(date)

" >> run_scrublet.log
    printf "    Process versions:
        $(python --version)
            $(pip freeze | grep scrublet | tr '==' ' ')

" >> run_scrublet.log

    if [ false == 'false' ]
    then
        run_scrublet.py --key Spleen --mat Spleen_for_scrub.mtx
        echo '    Process command:
        run_scrublet.py --key Spleen --mat Spleen_for_scrub.mtx
'  >> run_scrublet.log
    else
        run_scrublet.py --key Spleen --mat Spleen_for_scrub.mtx --skip
        echo '    Process command:
        run_scrublet.py --key Spleen --mat Spleen_for_scrub.mtx --skip
'  >> run_scrublet.log
        printf "    Scrublet skipped by request

" >> run_scrublet.log
    fi

    printf "** End process 'run_scrublet' at: $(date)

" >> run_scrublet.log

    printf "** Start processes to generate qc metrics and dashboard at: $(date)

" >> run_scrublet.log
