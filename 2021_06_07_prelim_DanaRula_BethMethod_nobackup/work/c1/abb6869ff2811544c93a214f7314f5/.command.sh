#!/bin/bash -ue
spec=`awk 'BEGIN {FS=",";OFS=","}{split($2,a,"_fq_part");gsub("[_ /-]", ".", a[1]);print($1, a[1], $3)}' good_sample_sheet.csv | awk 'BEGIN {FS=","}; $2=="Pancreas" {print $3}' | uniq`
star_mem=`awk -v var="$spec" '$1==var {print $3}' /net/trapnell/vol1/home/readdf/bin/bbi-sci/bin/star_file.txt | uniq`
star_path=`awk -v var="$spec" '$1==var {print $2}' /net/trapnell/vol1/home/readdf/bin/bbi-sci/bin/star_file.txt | uniq`
gtf_path=`awk -v var="$spec" '$1==var {print $2}' /net/trapnell/vol1/home/readdf/bin/bbi-sci/bin/gene_file.txt | uniq`

# capture process environment
set +u
echo star_path=$star_path > .command.env
echo star_mem=$star_mem >> .command.env
echo gtf_path=$gtf_path >> .command.env
echo gtf_path=$gtf_path >> .command.env
