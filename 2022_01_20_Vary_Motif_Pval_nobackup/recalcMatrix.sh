#! /bin/bash



script_dir="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_20_Vary_Motif_Pval_nobackup/bbi-sciatac-analyze/src"
humanFASTA="/net/bbi/vol1/data/genomes_stage/human/human_atac/Homo_sapiens.GRCh38.dna.toplevel.fa.finished"
mergedPeaks="/net/trapnell/vol1/HuBMAP/novaseq/210111_Riza_sciATAC3_split_sample/analyze_out/W136.heart.LV.s1/call_peaks/W136.heart.LV.s1-merged_peaks.bed"

bbiMotifs="/net/trapnell/vol1/ajh24/common_data/sciatac_pipeline_data/common_files/motifs/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.pfm"
outDir="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_20_Vary_Motif_Pval_nobackup/fileOutputs"

# Re-run motif calling using a different p value

PVAL="1e-5" # "1e-7" is default in the pipeline

# GCBIN="01" # Use [1,25] 





# Setting upa  virtual environment to run the script(s) here
module load virtualenv/20.0.27
# virtualenv $script_dir/python_env # Ran initially to set up the virtualenv

source $script_dir/python_env/bin/activate



# python3 -m ensurepip

# pip install -r /net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_20_Vary_Motif_Pval_nobackup/bbi-sciatac-analyze/python_requirements.tx

# Hit the following error
: <<'END'
ERROR: pip's dependency resolver does not currently take into account all the packages that are installed. This behaviour is the source of the following dependency conflicts.
tox 3.24.0 requires six>=1.14.0, but you have six 1.12.0 which is incompatible.
END
################


# git clone https://github.com/andrewhill157/barcodeutils.git


# pushd barcodeutils
# python setup.py install
# popd



# python $script_dir/call_peak_motifs.py $humanFASTA $mergedPeaks $bbiMotifs $outDir/PEAK_CALLS_${GCBIN}pval$PVAL --gc_bin $GCBIN --pwm_threshold $PVAL


# Loop and run the peak calling for each of 25 GC blocks. 
# Within each loop, add the output file to a variable to pass into the next step, where we'll combine and make a peak x motif matrix
peakMotifList=""
for eachGCbin in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
do
	python $script_dir/call_peak_motifs.py $humanFASTA $mergedPeaks $bbiMotifs $outDir/PEAK_CALLS_${eachGCbin}pval$PVAL --gc_bin $eachGCbin --pwm_threshold $PVAL
	peakMotifList=${peakMotifList}" "$outDir"/PEAK_CALLS_"${eachGCbin}"pval"$PVAL
done


echo $peakMotifList



# Now make these output files into a single peak x motif matrix
outputMatrix="/net/trapnell/vol1/home/readdf/trapLabDir/hubmap/results/2022_01_20_Vary_Motif_Pval_nobackup/fileOutputs/peak_x_motif_matrix_"
# peakMotifList="$outDir/PEAK_CALLS_${GCBIN}pval$PVAL" # Testing only


python ${script_dir}/generate_motif_matrix.py \
	--peak_motif_files $peakMotifList \
	--fasta $humanFASTA \
	--peaks $mergedPeaks \
	--motifs $bbiMotifs \
	--peak_tf_matrix ${outputMatrix}pVal${PVAL}.mtx



















