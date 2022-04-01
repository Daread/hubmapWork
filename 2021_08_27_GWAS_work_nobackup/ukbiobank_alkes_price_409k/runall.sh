# Get all the files
wget -r -np -nH --cut-dirs=3 -e robots=off -R "index.html*" https://data.broadinstitute.org/alkesgroup/UKBB/UKBB_409K/

mkdir -p ldsc_sumstats_files

module unload python
module load python/2.7.3
for fn in *.sumstats.gz; do
    output_prefix=ldsc_sumstats_files/${fn%.sumstats.gz}
    python ~ajh24/bin/ldsc/munge_sumstats.py --sumstats $fn --out $output_prefix
    rm $fn
done

python make_sumstats_metadata.py

