# Get fetal-specific peaks and shared fetal-adult peaks
# bedtools intersect -wa -a ../../../data/reads/2018_12_21_bdrl_batch1_batch2_novaseq/nobackup/analysis/call_peaks/merged_peaks.bed -b /net/shendure/vol10/projects/scATAC/nobackup/human_atlas/Roadmap_Epi_Data/fetalorigin_onetoone_nobrain_combined_merged.bed > fetal_origin_peaks.bed
# cat fetal_origin_peaks.bed | sed 's/\t/_/g' > fetal_origin_peaks.txt

# bedtools intersect -wa -a ../../../data/reads/2018_12_21_bdrl_batch1_batch2_novaseq/nobackup/analysis/call_peaks/merged_peaks.bed -b /net/shendure/vol10/projects/scATAC/nobackup/human_atlas/Roadmap_Epi_Data/fetalspecific_onetoone_nobrain_combined_merged.bed > fetal_specific_peaks.bed
# cat fetal_specific_peaks.bed | sed 's/\t/_/g' > fetal_specific_peaks.txt

# cat fetal_origin_peaks.bed fetal_specific_peaks.bed | sort -k1,1V -k2,2n -k3,3n | uniq > fetal_peaks.master.bed




module load bedtools/2.29.2
#DFR: I've run gatherInputFiles.R, which has output a master bed file and bed files for each cell type's expressed genes + linked sites. All are as unsorted bedfiles.

# Note the hyperparameters used, detaililng site sizes and numbers
coaccessCut=0.05 # Cicero co-accessibility cutoff to count
maxSites=5 # Max number of sites to link
peakSize=600 # Size of distal sites
# Sort 
bedfilePathBase="./fileOutputs/bedFilesForGWAS_1000_200_"$coaccessCut"_"$maxSites"_"$peakSize"/1000_200_"$coaccessCut"_"$maxSites"_"$peakSize"_"
samplesheetPath="./1000_200_"$coaccessCut"_"$maxSites"_"$peakSize"_samplesheet.txt"

echo -e "sample_id\tsites" > $samplesheetPath

# Full set
########################
# # Loop and sort each one
# for cellFile in master VSM_and_Pericyte Vascular_Endothelium T_Cell Neuronal Macrophage Mast_Cell Lymphatic_Endothelium Fibroblast Endocardium Cardiomyocyte B_Cell Adipocytes
# do
# 	echo $bedfilePathBase$cellFile"GeneAndPeakFile.bed"
# 	bedtools sort -i $bedfilePathBase$cellFile"GeneAndPeakFile.bed" > $bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed"
# done
# # Don't add the master bedfile to the sample sheet
# for cellFile in VSM_and_Pericyte Vascular_Endothelium T_Cell Neuronal Macrophage Mast_Cell Lymphatic_Endothelium Fibroblast Endocardium Cardiomyocyte B_Cell Adipocytes
# do
# 	echo -e $cellFile"\t"$bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed" >> $samplesheetPath
# done
######################## # End of full set

# Testing
##############################
# Loop and sort each one
for cellFile in master Vascular_Endothelium Cardiomyocyte 
do
	echo $bedfilePathBase$cellFile"GeneAndPeakFile.bed"
	cat $bedfilePathBase$cellFile"GeneAndPeakFile.bed" | awk '$1 = "chr" $1' | awk '{print $1 "\t" $2 "\t" $3}' > $bedfilePathBase$cellFile"GeneAndPeakFileFormatted.bed"
	bedtools sort -i $bedfilePathBase$cellFile"GeneAndPeakFileFormatted.bed" > $bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed"
	# bedtools sort -i $bedfilePathBase$cellFile"GeneAndPeakFile.bed" > $bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed"
done
# Don't add the master bedfile to the sample sheet
for cellFile in Vascular_Endothelium Cardiomyocyte
do
	echo -e $cellFile"\t"$bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed" >> $samplesheetPath
done
############################## # End testing


# Format as in Andrew's setup. Format is
# chr1	start	end
# no other fields, chr in front of chrom number
# cat $bedfilePathBase$cellFile"GeneAndPeakFile.bed" | awk '$1 = "chr" $1' | awk '{print $1 "\t" $2 "\t" $3}' > $bedfilePathBase$cellFile"GeneAndPeakFileFormatted.bed"

#samplesheet=samplesheet.final_annotations.txt
#samplesheet=samplesheet.brain_clusters.txt
# samplesheet=samplesheet.fetal_specific.txt
samplesheet=$samplesheetPath # Paths to bedfiles


# Generate a sample sheet with top N specific peaks per annotation from specificity scores
# if [ ! -f $samplesheet ]; then
# 	Rscript generate_input_files.R
# fi

# merged_peak_set="fetal_peaks.master.bed"
merged_peak_set=$bedfilePathBase"masterGeneAndPeakFileSorted.bed"
#merged_peak_set="/net/trapnell/vol1/ajh24/proj/2018three_level_sciatac/data/reads/2018_12_21_bdrl_batch1_batch2_novaseq/nobackup/analysis/call_peaks/merged_peaks.bed"
# gwas_set="/net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/data/gwas_full_results/ldsc_sumstats_files/sumstat_metadata.txt"
gwas_set="./ukbiobank_alkes_price_409k/ldsc_sumstats_files/miniHeart_sumstat_metadata.txt"


#output_directory="ldsc_results.final_annotations"
#output_directory="ldsc_results.brain_clusters"
# output_directory="ldsc_results.fetal_specific"
output_directory="./fileOutputs/bedFilesForGWAS_1000_200_"$coaccessCut"_"$maxSites"_"$peakSize"/ldsc_results"

mkdir -p $output_directory

# python /net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/bin/ldscore_regression/ldscore_peaks.py \
# --sample_sheet $samplesheet \
# --sumstats $gwas_set \
# --output_directory $output_directory \
# --master_peaks $merged_peak_set



# gwas_set="/net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/data/ukbiobank_alkes_price_409k/ldsc_sumstats_files/sumstat_metadata.txt"
# #output_directory=ldsc_results.ukbiobank_alkes_price_409k.final_annotations
# #output_directory=ldsc_results.ukbiobank_alkes_price_409k.brain_clusters
# output_directory=ldsc_results.ukbiobank_alkes_price_409k.fetal_specific
# mkdir -p $output_directory

# 9-8-21: Not sure if I even need this? Set up .bashrc to load these
# export PYTHONPATH=""
# module add python/3.7.7
# module add drmaa/0.7.9
# module add easygrid/0.1



# python /net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/bin/ldscore_regression/ldscore_peaks.py --sample_sheet $samplesheet \
# --sumstats $gwas_set \
# --output_directory $output_directory \
# --master_peaks $merged_peak_set \
# --score_sumstats_options " --chisq-max 9999" 

# python /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/ldscore_peaks.py --sample_sheet $samplesheet \
# --sumstats $gwas_set \
# --output_directory $output_directory \
# --master_peaks $merged_peak_set \
# --score_sumstats_options " --chisq-max 9999" 




# # Debugging. Reroute output
# python /net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/bin/ldscore_regression/ldscore_peaks.py --sample_sheet $samplesheet \
# --sumstats $gwas_set \
# --output_directory $output_directory \
# --master_peaks $merged_peak_set \
# --score_sumstats_options " --chisq-max 9999" \
# &> $output_directory/ldscoreLog.txt



python /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/ldscore_regression_scripts/ldscore_peaks.py --sample_sheet $samplesheet \
--sumstats $gwas_set \
--output_directory $output_directory \
--master_peaks $merged_peak_set \
--score_sumstats_options " --chisq-max 9999" 


