# Get fetal-specific peaks and shared fetal-adult peaks
bedtools intersect -wa -a ../../../data/reads/2018_12_21_bdrl_batch1_batch2_novaseq/nobackup/analysis/call_peaks/merged_peaks.bed -b /net/shendure/vol10/projects/scATAC/nobackup/human_atlas/Roadmap_Epi_Data/fetalorigin_onetoone_nobrain_combined_merged.bed > fetal_origin_peaks.bed
cat fetal_origin_peaks.bed | sed 's/\t/_/g' > fetal_origin_peaks.txt

bedtools intersect -wa -a ../../../data/reads/2018_12_21_bdrl_batch1_batch2_novaseq/nobackup/analysis/call_peaks/merged_peaks.bed -b /net/shendure/vol10/projects/scATAC/nobackup/human_atlas/Roadmap_Epi_Data/fetalspecific_onetoone_nobrain_combined_merged.bed > fetal_specific_peaks.bed
cat fetal_specific_peaks.bed | sed 's/\t/_/g' > fetal_specific_peaks.txt

cat fetal_origin_peaks.bed fetal_specific_peaks.bed | sort -k1,1V -k2,2n -k3,3n | uniq > fetal_peaks.master.bed

#samplesheet=samplesheet.final_annotations.txt
#samplesheet=samplesheet.brain_clusters.txt
samplesheet=samplesheet.fetal_specific.txt

# Generate a sample sheet with top N specific peaks per annotation from specificity scores
if [ ! -f $samplesheet ]; then
	Rscript generate_input_files.R
fi

merged_peak_set="fetal_peaks.master.bed"
#merged_peak_set="/net/trapnell/vol1/ajh24/proj/2018three_level_sciatac/data/reads/2018_12_21_bdrl_batch1_batch2_novaseq/nobackup/analysis/call_peaks/merged_peaks.bed"
gwas_set="/net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/data/gwas_full_results/ldsc_sumstats_files/sumstat_metadata.txt"

#output_directory="ldsc_results.final_annotations"
#output_directory="ldsc_results.brain_clusters"
output_directory="ldsc_results.fetal_specific"
mkdir -p $output_directory

python ~ajh24/proj/2017mouse_cell_atlas/bin/ldscore_regression/ldscore_peaks.py \
--sample_sheet $samplesheet \
--sumstats $gwas_set \
--output_directory $output_directory \
--master_peaks $merged_peak_set

gwas_set="/net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/data/ukbiobank_alkes_price_409k/ldsc_sumstats_files/sumstat_metadata.txt"
#output_directory=ldsc_results.ukbiobank_alkes_price_409k.final_annotations
#output_directory=ldsc_results.ukbiobank_alkes_price_409k.brain_clusters
output_directory=ldsc_results.ukbiobank_alkes_price_409k.fetal_specific
mkdir -p $output_directory

python ~ajh24/proj/2017mouse_cell_atlas/bin/ldscore_regression/ldscore_peaks.py --sample_sheet $samplesheet \
--sumstats $gwas_set \
--output_directory $output_directory \
--master_peaks $merged_peak_set \
--score_sumstats_options " --chisq-max 9999"
