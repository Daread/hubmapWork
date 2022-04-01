# pull out RSIDs and p-values from some full studies that I downloaded (each of which has its own unique format)
# Note that I am trying to get N in some cases and for one or two studies it is a little tricky but fine for most part
module unload python
module load python/2.7.3

mkdir -p ldsc_sumstats_files

# N taken from https://www.cng.fr/gabriel/study_description.html
python ~ajh24/bin/ldsc/munge_sumstats.py --sumstats /net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/data/gwas_full_results/gabriel.asthma.txt \
--out ldsc_sumstats_files/gabriel.asthma.snps \
--merge-alleles ~ajh24/bin/ldsc/w_hm3.snplist \
--signed-sumstats OR_fix,1 \
--p P_fix \
--N-cas 10365 \
--N-con 16110

# This study already has N
python ~ajh24/bin/ldsc/munge_sumstats.py --sumstats /net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/data/gwas_full_results/ckdgen.chronic_kidney_disease.txt \
--out ldsc_sumstats_files/ckdgen.chronic_kidney_disease.snps \
--merge-alleles ~ajh24/bin/ldsc/w_hm3.snplist

python ~ajh24/bin/ldsc/munge_sumstats.py --sumstats /net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/data/gwas_full_results/bronson.iga_deficiency.txt \
--out ldsc_sumstats_files/bronson.iga_deficiency.snps \
--merge-alleles ~ajh24/bin/ldsc/w_hm3.snplist \
--N-cas 1635 \
--N-con 4852

# N from GWAS catalog for Male results Barban N
python ~ajh24/bin/ldsc/munge_sumstats.py --sumstats /net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/data/gwas_full_results/barban.number_children_born_male.txt \
--out ldsc_sumstats_files/barban.number_children_born_male.snps \
--merge-alleles ~ajh24/bin/ldsc/w_hm3.snplist \
--N 103909

python ~ajh24/bin/ldsc/munge_sumstats.py --sumstats /net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/data/gwas_full_results/barban.age_first_birth_male.txt \
--out ldsc_sumstats_files/barban.age_first_birth_male.snps \
--merge-alleles ~ajh24/bin/ldsc/w_hm3.snplist \
--N 48408 
