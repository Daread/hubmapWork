

module load bedtools/2.29.2
#DFR: I've run gatherInputFiles.R, which has output a master bed file and bed files for each cell type's expressed genes + linked sites. All are as unsorted bedfiles.
cd /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup

# Note the hyperparameters used, detaililng site sizes and numbers
# coaccessCut=0.05 # Cicero co-accessibility cutoff to count
# maxSites=0 # Max number of sites to link
# peakSize=600 # Size of distal sites
# upstream=0
# downstream=0

coaccessCut=$1 # Cicero co-accessibility cutoff to count
maxSites=$2 # Max number of sites to link
peakSize=$3 # Size of distal sites
upstream=$4
downstream=$5

# coaccessCut=0.05 # Cicero co-accessibility cutoff to count
# maxSites=5 # Max number of sites to link
# peakSize=600 # Size of distal sites
# upstream=1000
# downstream=200


# Run the script that gets input bed files ready.
Rscript gatherInputFiles.R --promoterUpstream $upstream --promoterDownstream $downstream --coaccessCutoff $coaccessCut --maxNdistalSites $maxSites --peakSize $peakSize

# Sort 
# subPath="./fileOutputs/bedFilesForGWAS_"$upstream"_"$downstream"_"$coaccessCut"_"$maxSites"_"$peakSize"/"
# mkdir $subPath

bedfilePathBase="./fileOutputs/bedFilesForGWAS_"$upstream"_"$downstream"_"$coaccessCut"_"$maxSites"_"$peakSize"/"$upstream"_"$downstream"_"$coaccessCut"_"$maxSites"_"$peakSize"_"
samplesheetPath="./"$upstream"_"$downstream"_"$coaccessCut"_"$maxSites"_"$peakSize"_samplesheet.txt"

echo -e "sample_id\tsites" > $samplesheetPath

# Full set
########################
# # Loop and sort each one
for cellFile in master VSM_and_Pericyte Vascular_Endothelium T_Cell Neuronal Macrophage Mast_Cell Lymphatic_Endothelium Fibroblast Endocardium Cardiomyocyte B_Cell Adipocytes
do
	echo $bedfilePathBase$cellFile"GeneAndPeakFile.bed"
	cat $bedfilePathBase$cellFile"GeneAndPeakFile.bed" | awk '$1 = "chr" $1' | awk '{print $1 "\t" $2 "\t" $3}' > $bedfilePathBase$cellFile"GeneAndPeakFileFormatted.bed"
	bedtools sort -i $bedfilePathBase$cellFile"GeneAndPeakFileFormatted.bed" > $bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed"
	# bedtools sort -i $bedfilePathBase$cellFile"GeneAndPeakFile.bed" > $bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed"
done
# Don't add the master bedfile to the sample sheet
for cellFile in VSM_and_Pericyte Vascular_Endothelium T_Cell Neuronal Macrophage Mast_Cell Lymphatic_Endothelium Fibroblast Endocardium Cardiomyocyte B_Cell Adipocytes
do
	echo -e $cellFile"\t"$bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed" >> $samplesheetPath
done
######################## # End of full set

samplesheet=$samplesheetPath # Paths to bedfiles

merged_peak_set=$bedfilePathBase"masterGeneAndPeakFileSorted.bed"
#merged_peak_set="/net/trapnell/vol1/ajh24/proj/2018three_level_sciatac/data/reads/2018_12_21_bdrl_batch1_batch2_novaseq/nobackup/analysis/call_peaks/merged_peaks.bed"
# gwas_set="/net/trapnell/vol1/ajh24/proj/2017mouse_cell_atlas/data/gwas_full_results/ldsc_sumstats_files/sumstat_metadata.txt"
# gwas_set="./ukbiobank_alkes_price_409k/ldsc_sumstats_files/miniHeart_sumstat_metadata.txt" # For testing. Only includes 3 traits (high/low BP, and cardiac disease)
gwas_set="./ukbiobank_alkes_price_409k/ldsc_sumstats_files/heartOnly_sumstat_metadata.txt"


output_directory="./fileOutputs/bedFilesForGWAS_"$upstream"_"$downstream"_"$coaccessCut"_"$maxSites"_"$peakSize"/ldsc_results_full" #$$coaccessCut"_"$maxSites"_"$peakSize"/"$upstream"_"$downstream"_"$coaccessCut"_"$maxSites"_"$peakSize"_"

rm -r output_directory
mkdir -p $output_directory

# 9-8-21: Not sure if I even need this? Set up .bashrc to load these
# export PYTHONPATH=""
# module add python/3.7.7
# module add drmaa/0.7.9
# module add easygrid/0.1


# No liftover. Don't use this one
# python /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/ldscore_regression_scripts/ldscore_peaks.py --sample_sheet $samplesheet \
# --sumstats $gwas_set \
# --output_directory $output_directory \
# --master_peaks $merged_peak_set \
# --score_sumstats_options " --chisq-max 9999" 


# Liftover from hg19 (Andrew's SNP files) to hg38 (Current ATAC/RNA Annotation files)
python /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/ldscore_regression_scripts/ldscore_peaks.py --sample_sheet $samplesheet \
--sumstats $gwas_set \
--output_directory $output_directory \
--master_peaks $merged_peak_set \
--liftover_chain /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2021_08_27_GWAS_work_nobackup/liftoverFiles/hg19ToHg38.over.chain.gz \
--score_sumstats_options " --chisq-max 9999"  # At times Andrew ran without this flag. Look into this


Rscript andrew_generate_figures.R --promoterUpstream $upstream --promoterDownstream $downstream --coaccessCutoff $coaccessCut --maxNdistalSites $maxSites --peakSize $peakSize














# Testing
##############################
# Loop and sort each one
# for cellFile in master Vascular_Endothelium Cardiomyocyte 
# do
# 	echo $bedfilePathBase$cellFile"GeneAndPeakFile.bed"
# 	cat $bedfilePathBase$cellFile"GeneAndPeakFile.bed" | awk '$1 = "chr" $1' | awk '{print $1 "\t" $2 "\t" $3}' > $bedfilePathBase$cellFile"GeneAndPeakFileFormatted.bed"
# 	bedtools sort -i $bedfilePathBase$cellFile"GeneAndPeakFileFormatted.bed" > $bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed"
# 	# bedtools sort -i $bedfilePathBase$cellFile"GeneAndPeakFile.bed" > $bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed"
# done
# # Don't add the master bedfile to the sample sheet
# for cellFile in Vascular_Endothelium Cardiomyocyte
# do
# 	echo -e $cellFile"\t"$bedfilePathBase$cellFile"GeneAndPeakFileSorted.bed" >> $samplesheetPath
# done
############################## # End testing
